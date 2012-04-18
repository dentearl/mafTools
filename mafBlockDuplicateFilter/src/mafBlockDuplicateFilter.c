/* 
 * Copyright (C) 2012 by 
 * Dent Earl (dearl@soe.ucsc.edu, dentearl@gmail.com)
 * ... and other members of the Reconstruction Team of David Haussler's 
 * lab (BME Dept. UCSC).
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE. 
 */
#include <getopt.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "common.h"

int g_verbose_flag = 0;
int g_debug_flag = 0;
const int kMaxSeqName = 1 << 8;

typedef struct mafLine {
    // a mafLine struct is a single line of 
    char *line; // the entire line
    char *species; // just the species name
    char *sequence; // just the sequence field
    double score;
    struct mafLine *next;
} mafLine_t;
typedef struct duplicate {
    // a duplicate is a species that shows up twice in a block
    char *species;
    mafLine_t *headMaf; // linked list of mafLine_t containing the duplicated lines
    bool reported;
    struct duplicate *next;
} duplicate_t;

void usage(void);
void speciesNameCopy(char *species, char *line, int n);
void sequenceCopy(char *seq, char *line);
void processBody(void);
void destroyBlock(mafLine_t *m);
void destroyDuplicates(duplicate_t *d);
void destroyStringArray(char **sArray, int n);
mafLine_t* newMafLine(void);
int cmp_by_score(const void *a, const void *b);

void parseOptions(int argc, char **argv) {
    extern int g_debug_flag;
    extern int g_verbose_flag;
    extern const int kMaxStringLength;
    int c;
    while (1) {
        static struct option long_options[] = {
            {"debug", no_argument, 0, 'd'},
            {"verbose", no_argument, 0, 'v'},
            {"help", no_argument, 0, 'h'},
            {0, 0, 0, 0}
        };
        int option_index = 0;
        c = getopt_long(argc, argv, "n:c:p:v",
                        long_options, &option_index);
        if (c == -1)
            break;
        switch (c) {
        case 0:
            break;
        case 'v':
            g_verbose_flag++;
            break;
        case 'd':
            g_debug_flag = 1;
            break;
        case 'h':
        case '?':
            usage();
            break;
        default:
            abort();
        }
    }
    // Check there's nothing left over on the command line 
    if (optind < argc) {
        char *errorString = de_malloc(kMaxStringLength);
        strcpy(errorString, "Unexpected arguments:");
        while (optind < argc) {
            strcat(errorString, " ");
            strcat(errorString, argv[optind++]);
        }
        fprintf(stderr, "%s\n", errorString);
        free(errorString);
        usage();
    }
}
void usage(void) {
    fprintf(stderr, "Usage: mafBlockDuplicateFilter < mafWithDuplicates.maf > pruned.maf \n\n"
            "mafBlockDuplicateFilter is a program that will \n\n");
    fprintf(stderr, "Options: \n"
            "  -h, --help     show this help message and exit.\n"
            "  -v, --verbose  turns on verbose output.\n");
    exit(EXIT_FAILURE);
}
mafLine_t* newMafLine(void) {
    mafLine_t *m = (mafLine_t *) de_malloc(sizeof(*m));
    m->species = NULL;
    m->sequence = NULL;
    m->line = NULL;
    m->next = NULL;
    m->score = 0.0;
    return m;
}
duplicate_t* newDuplicate(void) {
    duplicate_t *d = (duplicate_t *) de_malloc(sizeof(*d));
    d->species = NULL;
    d->headMaf = NULL;
    d->reported = false;
    d->next = NULL;
    return d;
}
int numberOfSequences(mafLine_t *m) {
    // count the number of actual sequence lines in the mafline linked list
    int s = 0;
    while (m != NULL) {
        if (m->line[0] == 's')
            ++s;
        m = m->next;
    }
    return s;
}
unsigned longestLine(mafLine_t *head) {
    // walk the mafline linked list and return the longest m->line value
    mafLine_t *m = head;
    unsigned max = 0;
    while (m != NULL) {
        if (max < strlen(m->line))
            max = strlen(m->line);
        m = m->next;
    }
    return max;
}
void printResidues(unsigned *r) {
    // debug function
    for (int i = 0; i < 5; ++i)
        printf("%d:%d, ", i, r[i]);
    printf("\n");
}
unsigned maxRes(unsigned residues[]) {
    // walk the residues array and return the largest value
    unsigned m = 0;
    for (int i = 0; i < 5; ++i) {
        if (m < residues[i])
            m = residues[i];
    }
    return m;
}
char consensusResidue(unsigned residues[]) {
    // given an unsigned array of counts of the 4 bases in the order A, C, G, T
    // return as a char the IUPAC for the consensus residue.
    unsigned m = maxRes(residues);
    bool maxA = false, maxC = false, maxG = false, maxT = false, allGap = false;
    if (residues[0] == m)
        maxA = true;
    if (residues[1] == m)
        maxC = true;
    if (residues[2] == m)
        maxG = true;
    if (residues[3] == m)
        maxT = true;
    if (residues[4] == m)
        allGap = true;
    if (maxA && maxC && maxG && maxT)
        return 'N';
    if (maxA && maxC && maxG)
        return 'V';
    if (maxA && maxC && maxT)
        return 'H';
    if (maxA && maxG && maxT)
        return 'D';
    if (maxC && maxG && maxT)
        return 'B';
    if (maxA && maxC)
        return 'M';
    if (maxA && maxG)
        return 'R';
    if (maxA && maxT)
        return 'W';
    if (maxC && maxG)
        return 'S';
    if (maxC && maxT)
        return 'Y';
    if (maxG && maxT)
        return 'K';
    if (maxA)
        return 'A';
    if (maxC)
        return 'C';
    if (maxG)
        return 'G';
    if (maxT)
        return 'T';
    if (allGap)
        return '-';
    return '?';
}
void buildConsensus(char *consensus, char **sequences, int numSeqs, unsigned lineno) {
    // given an empty string of the correct length, `consensus', a string array
    // containing all of the sequences, build the consensus sequence and store it
    // in the provided string.
    unsigned residues[5]; // order: {A, C, G, T, -}
    int n = (int) strlen(sequences[0]);
    for (int i = 0; i < n; ++i) {
        // columns
        memset(residues, 0, sizeof(residues));
        // printResidues(residues);
        for (int j = 1; j < numSeqs; ++j) {
            // rows
            switch (sequences[j][i]) {
            case 'A':
            case 'a':
                residues[0]++;
                break;
            case 'C':
            case 'c':
                residues[1]++;
                break;
            case 'G':
            case 'g':
                residues[2]++;
                break;
            case 'T':
            case 't':
                residues[3]++;
                break;
            case '-':
                break;
            default:
                fprintf(stderr, "Error, unanticipated character within sequence (%d,%d) "
                        "contained in block near line number %u: %c\n", 
                        j, i, lineno, sequences[j][i]);
                exit(EXIT_FAILURE);
            }
        }
        consensus[i] = consensusResidue(residues);
    }
    consensus[n] = '\0';
}
bool checkForDupes(char **species, int index, mafLine_t *m) {
    // walk through the species string array and check to see if m->species is contained
    // anywhere within.
    for (int i = 0; i < index; ++i) {
        if (!strcmp(species[i], m->species)) {
            return true;
        }
    }
    return false;
}
void reportBlock(mafLine_t *b) {
    // print out a maf block in the form of the mafline linked list
    while (b != NULL) {
        printf("%s\n", b->line);
        b = b->next;
    }
    printf("\n");
}
void reportBlockWithDuplicates(mafLine_t *mafHead, duplicate_t *dupHead) {
    // report the block represented by the linked list mafHead. If a given line
    // is a member of the duplicate linked list, report only the top scoring duplicate
    // which will be the one stored at the head of the mafline linkeded list (dup->headMaf).
    mafLine_t *m = mafHead;
    duplicate_t *d = dupHead;
    while (m != NULL) {
        d = dupHead;
        bool isDup = false;
        debug("  examining line for species \"%s\"\n", m->species);
        while (d != NULL && m->species != NULL && d->species != NULL) {
            debug("    cross checking against dup \"%s\"\n", d->species);
            if (!strcmp(m->species, d->species)) {
                if (numberOfSequences(d->headMaf) > 1) {
                    isDup = true;
                    if (!strcmp(m->line, d->headMaf->line) && !d->reported) {
                        printf("%s\n", d->headMaf->line);
                        d->reported = true;
                        break;
                    }
                }
            }
            d = d->next;
        }
        if (!isDup)
            printf("%s\n", m->line);
        m = m->next;
    }
    printf("\n");
}
void reportDuplicates(duplicate_t *dup) {
    // debugging function
    printf("Duplicates: ");
    while (dup->species != NULL) {
        if (numberOfSequences(dup->headMaf) < 2) {
            dup = dup->next;
            continue;
        }
        mafLine_t *m = dup->headMaf;
        printf("    %s\n", dup->species);
        while (m != NULL) {
            printf("        %.2f %s\n", m->score, m->sequence);
            m = m->next;
        }
        dup = dup->next;
    }
}
duplicate_t* findDuplicate(duplicate_t *dup, char *species) {
    // walk the dup linked list and search for dup->species equal to species.
    // return the pointer to the dup in question if found, NULL if not found.
    while(dup != NULL && dup->species != NULL) {
        if (!strcmp(dup->species, species)) {
            return dup;
        }
        dup = dup->next;
    }
    return NULL;
}
double bitScore(char a, char b) {
    // a is the truth, b is the prediction
    a = toupper(a);
    b = toupper(b);
    if (a == 'N')
        return 0.0;
    double f = 0.41503749927884381854626105605; // -log2(3/4)
    switch (a) {
    case 'A':
    case 'C':
    case 'G':
    case 'T':
        return a == b ? 2.0 : 0.0;
    case 'W':
        return (b == 'A' || b == 'T') ? 1.0 : 0.0;
    case 'S':
        return (b == 'C' || b == 'G') ? 1.0 : 0.0;
    case 'M':
        return (b == 'A' || b == 'C') ? 1.0 : 0.0;
    case 'K':
        return (b == 'G' || b == 'T') ? 1.0 : 0.0;
    case 'R':
        return (b == 'A' || b == 'G') ? 1.0 : 0.0;
    case 'Y':
        return (b == 'C' || b == 'T') ? 1.0 : 0.0;
    case 'B':
        return (b != 'A') ? f : 0.0;
    case 'D':
        return (b != 'C') ? f : 0.0;
    case 'H':
        return (b != 'G') ? f : 0.0;
    case 'V':
        return (b != 'T') ? f : 0.0;
    default:
        fprintf(stderr, "Unanticipated condition when calling bitScore(%c, %c)\n", a, b);
        exit(EXIT_FAILURE);
    }
}
double scoreSequence(char *consensus, char *seq) {
    // walk the sequence `seq' and calculate and return its score versus the sequence `consensus'
    double s = 0.0;
    for (unsigned i = 0; i < strlen(consensus); ++i)
        s += bitScore(consensus[i], seq[i]);
    return s;
}
void populateMafLineArray(mafLine_t *head, mafLine_t **array) {
    // given a mafLine linked list, create an array of the same. this array will be used
    // later to sort the mafline pointers based on their ->score value.
    unsigned i = 0;
    mafLine_t *m = head;
    while (m != NULL) {
        array[i++] = m;
        m = m->next;
    }
}
void findBestDupes(duplicate_t *head, char *consensus) {
    // For each duplicate, go through its mafline list and find the best line and move it
    // to the head of the list. 
    duplicate_t *d = head;
    mafLine_t *m = NULL;
    while (d != NULL) {
        m = d->headMaf;
        int n = numberOfSequences(m);
        if (n < 2) {
            d = d->next;
            continue;
        }
        // score all the dupes
        while (m != NULL) {
            m->score = scoreSequence(consensus, m->sequence);
            m = m->next;
        }
        // sort on scores
        mafLine_t *mafLineArray[n];
        populateMafLineArray(d->headMaf, mafLineArray);
        qsort(mafLineArray, n, sizeof(mafLine_t *), cmp_by_score);
        // move the top score to the head of the list
        d->headMaf = mafLineArray[0];
        m = d->headMaf;
        for (int i = 1; i < n; ++i) {
            m->next = mafLineArray[i];
            m = m->next;
        }
        m->next = NULL;
        d = d->next;
    }
}
int cmp_by_score(const void *a, const void *b) {
    // mafBlock_t * const *ia = a;
    // mafBlock_t * const *ib = b;
    mafLine_t **ia = (mafLine_t **) a;
    mafLine_t **ib = (mafLine_t **) b;
    // reverse sort
    return ((*ib)->score - (*ia)->score);
}
void checkBlock(mafLine_t *head, unsigned lineno) {
    // read through each line of a mafBlock and filter duplicates.
    // Report the top scoring duplication only.
    int n = numberOfSequences(head);
    char **species = (char **) de_malloc(sizeof(char *) * n);
    char **sequences = (char **) de_malloc(sizeof(char *) * n);
    int index = 0;
    bool containsDuplicates = false;
    mafLine_t *m = head;
    duplicate_t *d, *dupSpeciesHead = newDuplicate();
    d = dupSpeciesHead;
    while (m != NULL) {
        if (m->line[0] != 's') {
            // skip non-sequence lines
            m = m->next;
            continue;
        }
        species[index] = (char *) de_malloc(kMaxSeqName);
        sequences[index] = (char *) de_malloc(strlen(m->sequence) + 1);
        strcpy(species[index], m->species);
        strcpy(sequences[index], m->sequence);
        duplicate_t *thisDup = findDuplicate(dupSpeciesHead, m->species);
        if (thisDup == NULL) {
            // add new duplicate species
            debug("adding new species %s\n", m->species);
            d->species = (char *) de_malloc(kMaxSeqName);
            strcpy(d->species, m->species);
            // create the mafline linked list
            d->headMaf = newMafLine();
            d->headMaf->species = (char *) de_malloc(kMaxSeqName);
            d->headMaf->sequence = (char *) de_malloc(strlen(m->sequence) + 1);
            d->headMaf->line = (char *) de_malloc(strlen(m->line) + 1);
            strcpy(d->headMaf->species, m->species);
            strcpy(d->headMaf->sequence, m->sequence);
            strcpy(d->headMaf->line, m->line);
            // reset d to accept a new duplicate
            d->next = newDuplicate();
            d = d->next;
        } else {
            debug("extending duplicate on species %s\n", m->species);
            containsDuplicates = true;
            mafLine_t *ml = thisDup->headMaf;
            while (ml->next != NULL)
                ml = ml->next;
            ml->next = newMafLine();
            ml = ml->next;
            ml->species = (char *) de_malloc(kMaxSeqName);
            ml->sequence = (char *) de_malloc(strlen(m->sequence) + 1);
            ml->line = (char *) de_malloc(strlen(m->line) + 1);
            strcpy(ml->species, m->species);
            strcpy(ml->sequence, m->sequence);
            strcpy(ml->line, m->line);
        }
        ++index;
        m = m->next;
    }
    if (!containsDuplicates) {
        reportBlock(head);
        destroyStringArray(species, n);
        destroyStringArray(sequences, n);
        destroyDuplicates(dupSpeciesHead);
        return;
    }
    // this block contains duplicates
    char *consensus = (char *) de_malloc(longestLine(head) + 1);
    consensus[0] = '\0';
    buildConsensus(consensus, sequences, n, lineno);
    findBestDupes(dupSpeciesHead, consensus);
    reportBlockWithDuplicates(head, dupSpeciesHead);
    destroyStringArray(species, n);
    destroyStringArray(sequences, n);
    destroyDuplicates(dupSpeciesHead);
    free(consensus);
}
void destroyBlock(mafLine_t *m) {
    // free all memory associated with a mafline linked list
    mafLine_t *tmp = NULL;
    while (m != NULL) {
        tmp = m;
        m = m->next;
        free(tmp->line);
        free(tmp->species);
        free(tmp->sequence);
        free(tmp);
    }
}
void destroyDuplicates(duplicate_t *d) {
    // free all memory associated with a duplicate linked list
    duplicate_t *tmp = NULL;
    while (d != NULL) {
        tmp = d;
        d = d->next;
        free(tmp->species);
        destroyBlock(tmp->headMaf);
        free(tmp);
    }
}
void destroyStringArray(char **sArray, int n) {
    // free all memory associated with a string array
    for (int i = 0; i < n; ++i) {
        free(sArray[i]);
    }
    free(sArray);
}
void speciesNameCopy(char *species, char *line, int n) {
    // reads `line' and extracts the species name, copies it into
    // `species'
    // so the line:
    //s hg19.chr22 885203...
    // will copy out only the `hg19' part.
    if (line[0] != 's') {
        // if this is not a sequence line, species should be NULL
        species[0] = '\0';
        return;
    }
    line++; // skip the `^s'
    while (*line == ' ' || *line == '\t')
        line++;
    int i = 0;
    while (*line != ' ' && *line != '.' && *line != '\t' && i < n)
        species[i++] = *(line++);
    species[i] = '\0';
}
void sequenceCopy(char *seq, char *line) {
    // walk a maf line `line' until the sequence field, then copy that entire field into `seq'
    char *tline = (char *) de_malloc(strlen(line) + 1);
    strcpy(tline, line);
    char *tkn = NULL;
    tkn = strtok(tline, " \t"); // 's' field
    if (tkn == NULL) {
        failBadFormat();
    }
    if (tkn[0] != 's') {
        free(tline);
        seq[0] = '\0';
        return;
    }
    tkn = strtok(NULL, " \t"); // name field
    if (tkn == NULL) {
        failBadFormat();
    }
    tkn = strtok(NULL, " \t"); // start field
    if (tkn == NULL) {
        failBadFormat();
    }
    tkn = strtok(NULL, " \t"); // length field
    if (tkn == NULL) {
        failBadFormat();
    }
    tkn = strtok(NULL, " \t"); // strand field
    if (tkn == NULL) {
        failBadFormat();
    }
    tkn = strtok(NULL, " \t"); // sourceLength field
    if (tkn == NULL) {
        failBadFormat();
    }
    tkn = strtok(NULL, " \t"); // sequence field
    if (tkn == NULL) {
        failBadFormat();
    }
    strcpy(seq, tkn);
    free(tline);
}
void removeTrailingWhitespace(char *s) {
    // given a string `s', walk it backwards and cut off all trailing whitespace
    for (int i = strlen(s) - 1; 0 <= i; --i) {
        if (s[i] == ' ' || s[i] == '\t' || s[i] == '\n' || s[i] == '\r')
            s[i] = '\0';
        else
            return;
    }
}
void processBody(void) {
    // walk the body of the maf file and process it, block by block.
    extern const int kMaxStringLength;
    FILE *ifp = stdin;
    unsigned lineno = 0;
    int32_t n = kMaxStringLength;
    char *buf = (char *) de_malloc(n);
    char *cline = NULL; // freed in destroyBlock
    char *species = NULL; // freed in destroyBlock
    char *sequence = NULL; // freed in destroyBlock
    mafLine_t *head = NULL, *tail = NULL;
    while (de_getline(&buf, &n, ifp) != -1) {
        lineno++;
        if (*buf == 0 && head != NULL) {
            // empty line or end of file
            checkBlock(head, lineno);
            destroyBlock(head);
            head = NULL;
            tail = NULL;
            continue;
        }
        removeTrailingWhitespace(buf);
        if (head == NULL) {
            // new block
            head = newMafLine();
            cline = (char *) de_malloc(n);
            species = (char *) de_malloc(kMaxSeqName);
            sequence = (char *) de_malloc(n);
            strcpy(cline, buf);
            speciesNameCopy(species, buf, kMaxSeqName);
            sequenceCopy(sequence, buf);
            head->species = species;
            head->sequence = sequence;
            head->line = cline;
            tail = head;
        } else {
            // extend block
            tail->next = newMafLine();
            cline = (char *) de_malloc(n);
            species = (char *) de_malloc(kMaxSeqName);
            sequence = (char *) de_malloc(n);
            strcpy(cline, buf);
            speciesNameCopy(species, buf, kMaxSeqName);
            sequenceCopy(sequence, buf);
            tail = tail->next;
            tail->species = species;
            tail->sequence = sequence;
            tail->line = cline;
            tail->next = NULL;
        }
    }
    destroyBlock(head);
    free(buf);
}
int main(int argc, char **argv) {
    parseOptions(argc, argv);

    processHeader();
    processBody();
    
    return EXIT_SUCCESS;
}
