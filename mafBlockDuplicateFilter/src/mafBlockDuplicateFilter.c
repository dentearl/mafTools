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
#include <ctype.h>
#include <getopt.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "common.h"
#include "sharedMaf.h"
#include "buildVersion.h"

const char *g_version = "version 0.1 September 2012";

typedef struct scoredMafLine {
    // augmented data structure
    double score;
    mafLine_t *mafLine;
    struct scoredMafLine *next;
} scoredMafLine_t;
typedef struct duplicate {
    // a duplicate is a species that shows up twice in a block
    char *species;
    scoredMafLine_t *headScoredMaf; // linked list of scoredMafLine_t containing the duplicated lines
    bool reported;
    struct duplicate *next;
} duplicate_t;

void usage(void);
void version(void);
void processBody(mafFileApi_t *mfa);
void checkBlock(mafBlock_t *block);
// void destroyBlock(mafLine_t *m);
void destroyScoredMafLineList(scoredMafLine_t *sml);
void destroyDuplicates(duplicate_t *d);
void destroyStringArray(char **sArray, int n);
scoredMafLine_t* newScoredMafLine(void);
unsigned longestLine(mafBlock_t *mb);
unsigned numberOfSequencesScoredMafLineList(scoredMafLine_t *m);
void reportBlockWithDuplicates(mafBlock_t *mb, duplicate_t *dupHead);
int cmp_by_score(const void *a, const void *b);

void parseOptions(int argc, char **argv, char *filename) {
    int c;
    int setMName = 0;
    while (1) {
        static struct option longOptions[] = {
            {"debug", no_argument, 0, 'd'},
            {"verbose", no_argument, 0, 'v'},
            {"help", no_argument, 0, 'h'},
            {"version", no_argument, 0, 0},
            {"maf",  required_argument, 0, 'm'},            
            {0, 0, 0, 0}
        };
        int longIndex = 0;
        c = getopt_long(argc, argv, "d:m:h:v",
                        longOptions, &longIndex);
        if (c == -1)
            break;
        switch (c) {
        case 0:
            if (strcmp("version", longOptions[longIndex].name) == 0) {
                version();
                exit(EXIT_SUCCESS);
            }
            break;
        case 'm':
            setMName = 1;
            sscanf(optarg, "%s", filename);
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
    if (!setMName) {
        fprintf(stderr, "specify --maf\n");
        usage();
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
void version(void) {
    fprintf(stderr, "mafBlockDuplicateFilter, %s\nbuild: %s, %s, %s\n\n", g_version, g_build_date, 
            g_build_git_branch, g_build_git_sha);
}
void usage(void) {
    version();
    fprintf(stderr, "Usage: mafBlockDuplicateFilter --maf mafWithDuplicates.maf > pruned.maf \n\n"
            "mafBlockDuplicateFilter is a program to filter out duplications from a Multiple \n"
            "Alignment Format (maf) file. This program assumes the sequence name field is \n"
            "formatted as in \"speciesName.chromosomeName\" using the first period charater, \".\",\n"
            "as the delimiter between the species name and the chromosome name. For every \n"
            "block present in the alignment, mBDF looks for any duplicated species within the \n"
            "block. Instead of stripping out all copies of the duplication, the sequence with \n"
            "the highest similarity to the consensus of the block is left, all others are \n"
            "removed. Sequence similarity is computed as a bit score in comparison to the \n"
            "IUPAC-enabled consensus. Ties are resolved by picking the sequence that appears \n"
            "earliest in the file. \n\n");
    fprintf(stderr, "Options: \n");
    usageMessage('h', "help", "show this help message and exit.");
    usageMessage('m', "maf", "path to maf file.");
    usageMessage('v', "verbose", "turns on verbose output.");
    exit(EXIT_FAILURE);
}
scoredMafLine_t* newScoredMafLine(void) {
    scoredMafLine_t *m = (scoredMafLine_t *) de_malloc(sizeof(*m));
    m->mafLine = NULL;
    m->next = NULL;
    m->score = 0.0;
    return m;
}
duplicate_t* newDuplicate(void) {
    duplicate_t *d = (duplicate_t *) de_malloc(sizeof(*d));
    d->species = NULL;
    d->headScoredMaf = NULL;
    d->reported = false;
    d->next = NULL;
    return d;
}
unsigned longestLine(mafBlock_t *mb) {
    // walk the mafline linked list and return the longest m->line value
    mafLine_t *m = maf_mafBlock_getHeadLine(mb);
    unsigned max = 0;
    while (m != NULL) {
        if (max < strlen(maf_mafLine_getLine(m)))
            max = strlen(maf_mafLine_getLine(m));
        m = maf_mafLine_getNext(m);
    }
    return max;
}
unsigned numberOfSequencesScoredMafLineList(scoredMafLine_t *m) {
    // count the number of actual sequence lines in a mafLine_t list
    unsigned s = 0;
    while (m != NULL) {
        if (maf_mafLine_getType(m->mafLine) == 's')
            ++s;
        m = m->next;
    }
    return s;
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
    // given an unsigned array of counts of the 4 bases in the order A, C, G, T, N
    // return as a char the IUPAC for the consensus residue.
    unsigned m = maxRes(residues);
    bool maxA = false, maxC = false, maxG = false, maxT = false, maxN = false, allGap = false;
    if (residues[0] == m)
        maxA = true;
    if (residues[1] == m)
        maxC = true;
    if (residues[2] == m)
        maxG = true;
    if (residues[3] == m)
        maxT = true;
    if (residues[4] == m)
        maxN = true;
    if (residues[5] == m)
        allGap = true;
    if ((maxA && maxC && maxG && maxT) || (maxN))
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
    unsigned residues[6]; // order: {A, C, G, T, N, -}
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
            case 'N':
            case 'n':
                residues[4]++;
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
        if (!strcmp(species[i], maf_mafLine_getSpecies(m))) {
            return true;
        }
    }
    return false;
}
void reportBlock(mafBlock_t *b) {
    // print out a maf block in the form of the mafline linked list
    mafLine_t *ml = maf_mafBlock_getHeadLine(b);
    while (ml != NULL) {
        printf("%s\n", maf_mafLine_getLine(ml));
        ml = maf_mafLine_getNext(ml);
    }
    printf("\n");
}
void reportBlockWithDuplicates(mafBlock_t *mb, duplicate_t *dupHead) {
    // report the block represented by mb. If a given line
    // is a member of the duplicate linked list, report only the top scoring duplicate
    // which will be the one stored at the head of the mafline linkeded list (dup->headScoredMaf).
    mafLine_t *m = maf_mafBlock_getHeadLine(mb);
    duplicate_t *d = dupHead;
    while (m != NULL) {
        d = dupHead;
        bool isDup = false;
        while (d != NULL && maf_mafLine_getSpecies(m) != NULL && d->species != NULL) {
            if (!strcmp(maf_mafLine_getSpecies(m), d->species)) {
                if (numberOfSequencesScoredMafLineList(d->headScoredMaf) > 1) {
                    isDup = true;
                    if (!strcmp(maf_mafLine_getLine(m), maf_mafLine_getLine(d->headScoredMaf->mafLine))
                        && !d->reported) {
                        printf("%s\n", maf_mafLine_getLine(d->headScoredMaf->mafLine));
                        d->reported = true;
                        break;
                    }
                }
            }
            d = d->next;
        }
        if (!isDup)
            printf("%s\n", maf_mafLine_getLine(m));
        m = maf_mafLine_getNext(m);
    }
    printf("\n");
}
void reportDuplicates(duplicate_t *dup) {
    // debugging function
    printf("Duplicates: ");
    while (dup->species != NULL) {
        if (numberOfSequencesScoredMafLineList(dup->headScoredMaf) < 2) {
            dup = dup->next;
            continue;
        }
        scoredMafLine_t *m = dup->headScoredMaf;
        printf("    %s\n", dup->species);
        while (m != NULL) {
            printf("        %.2f %s\n", m->score, maf_mafLine_getSequence(m->mafLine));
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
    if ((a == 'N') || (b == 'N'))
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
void populateMafLineArray(scoredMafLine_t *head, scoredMafLine_t **array) {
    // given a mafLine linked list, create an array of the same. this array will be used
    // later to sort the mafline pointers based on their ->score value.
    unsigned i = 0;
    scoredMafLine_t *m = head;
    while (m != NULL) {
        array[i++] = m;
        m = m->next;
    }
}
void findBestDupes(duplicate_t *head, char *consensus) {
    // For each duplicate, go through its mafline list and find the best line and move it
    // to the head of the list. 
    duplicate_t *d = head;
    scoredMafLine_t *m = NULL;
    while (d != NULL) {
        m = d->headScoredMaf;
        unsigned n = numberOfSequencesScoredMafLineList(m);
        if (n < 2) {
            d = d->next;
            continue;
        }
        // score all the dupes
        while (m != NULL) {
            m->score = scoreSequence(consensus, maf_mafLine_getSequence(m->mafLine));
            m = m->next;
        }
        // sort on scores
        scoredMafLine_t *mafLineArray[n];
        populateMafLineArray(d->headScoredMaf, mafLineArray);
        qsort(mafLineArray, n, sizeof(scoredMafLine_t *), cmp_by_score);
        // move the top score to the head of the list
        d->headScoredMaf = mafLineArray[0];
        m = d->headScoredMaf;
        for (unsigned i = 1; i < n; ++i) {
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
    scoredMafLine_t **ia = (scoredMafLine_t **) a;
    scoredMafLine_t **ib = (scoredMafLine_t **) b;
    // reverse sort
    return ((*ib)->score - (*ia)->score);
}
void correctSpeciesNames(mafBlock_t *block) {
    // the sharedMaf.h block reading function takes the entire name field,
    // but we only want the name field up until the first '.' is observed.
    mafLine_t *m = maf_mafBlock_getHeadLine(block);
    while(m != NULL) {
        if (maf_mafLine_getType(m) != 's') {
            m = maf_mafLine_getNext(m);
            continue;
        }
        for (unsigned i = 0; i < strlen(maf_mafLine_getSpecies(m)); ++i) {
            if (maf_mafLine_getSpecies(m)[i] == '.') {
                maf_mafLine_getSpecies(m)[i] = '\0';
                break;
            }
        }
        m = maf_mafLine_getNext(m);
    }
}
void checkBlock(mafBlock_t *block) {
    // read through each line of a mafBlock and filter duplicates.
    // Report the top scoring duplication only.
    mafLine_t *m = maf_mafBlock_getHeadLine(block);
    unsigned n = maf_mafLine_getNumberOfSequences(m);
    char **species = (char **) de_malloc(sizeof(char *) * n);
    char **sequences = (char **) de_malloc(sizeof(char *) * n);
    int index = 0;
    bool containsDuplicates = false;
    duplicate_t *d = NULL, *dupSpeciesHead = NULL;
    while (m != NULL) {
        if (maf_mafLine_getType(m) != 's') {
            // skip non-sequence lines
            m = maf_mafLine_getNext(m);
            continue;
        }
        species[index] = (char *) de_malloc(kMaxSeqName);
        sequences[index] = (char *) de_malloc(strlen(maf_mafLine_getSequence(m)) + 1);
        strcpy(species[index], maf_mafLine_getSpecies(m));
        strcpy(sequences[index], maf_mafLine_getSequence(m));
        duplicate_t *thisDup = findDuplicate(dupSpeciesHead, maf_mafLine_getSpecies(m));
        if (thisDup == NULL) {
            // add new duplicate species
            if (dupSpeciesHead == NULL) {
                dupSpeciesHead = newDuplicate();
                d = dupSpeciesHead;
            } else {
                d->next = newDuplicate();
                d = d->next;
            }
            d->species = (char *) de_malloc(kMaxSeqName);
            strcpy(d->species, maf_mafLine_getSpecies(m));
            // create the mafline linked list
            d->headScoredMaf = newScoredMafLine();
            d->headScoredMaf->mafLine = m;
        } else {
            // this sequence is a duplicate, extend the duplicate list.
            containsDuplicates = true;
            scoredMafLine_t *ml = thisDup->headScoredMaf;
            while (ml->next != NULL)
                ml = ml->next;
            ml->next = newScoredMafLine();
            ml = ml->next;
            ml->mafLine = m;
        }
        ++index;
        m = maf_mafLine_getNext(m);
    }
    if (!containsDuplicates) {
        reportBlock(block);
        destroyStringArray(species, n);
        destroyStringArray(sequences, n);
        destroyDuplicates(dupSpeciesHead);
        return;
    }
    // this block contains duplicates
    char *consensus = (char *) de_malloc(longestLine(block) + 1);
    consensus[0] = '\0';
    buildConsensus(consensus, sequences, n, maf_mafLine_getLineNumber(maf_mafBlock_getHeadLine(block)));
    findBestDupes(dupSpeciesHead, consensus);
    reportBlockWithDuplicates(block, dupSpeciesHead);
    destroyStringArray(species, n);
    destroyStringArray(sequences, n);
    destroyDuplicates(dupSpeciesHead);
    free(consensus);
}
void destroyDuplicates(duplicate_t *d) {
    // free all memory associated with a duplicate linked list
    duplicate_t *tmp = NULL;
    while (d != NULL) {
        tmp = d;
        d = d->next;
        free(tmp->species);
        destroyScoredMafLineList(tmp->headScoredMaf);
        free(tmp);
    }
}
void destroyScoredMafLineList(scoredMafLine_t *sml) {
    scoredMafLine_t *tmp = NULL;
    while (sml != NULL) {
        tmp = sml;
        sml = sml->next;
        // note that sml->mafLine is NOT freed here
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
void processBody(mafFileApi_t *mfa) {
    // walk the body of the maf file and process it, block by block.
    mafBlock_t *thisBlock = NULL;
    thisBlock = maf_readBlock(mfa); // unused
    maf_destroyMafBlockList(thisBlock);
    while((thisBlock = maf_readBlock(mfa)) != NULL) {
        correctSpeciesNames(thisBlock);
        checkBlock(thisBlock);
        maf_destroyMafBlockList(thisBlock);
    }
}
int main(int argc, char **argv) {
    char filename[kMaxStringLength];
    parseOptions(argc, argv, filename);
    
    mafFileApi_t *mfa = maf_newMfa(filename, "r");
    processBody(mfa);

    maf_destroyMfa(mfa);

    return EXIT_SUCCESS;
}
