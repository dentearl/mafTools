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
    scoredMafLine_t *tailScoredMaf; // last element in ll
    bool reported; // whether or not this duplicate has been reported yet
    struct duplicate *next;
    uint32_t numSequences; // number of elements in the headScoredMaf ll
} duplicate_t;

void usage(void);
void version(void);
void printHeader(void);
void processBody(mafFileApi_t *mfa, char *seq, char strand);
void checkBlock(mafBlock_t *block, char *seq, char strand);
// void destroyBlock(mafLine_t *m);
void destroyScoredMafLineList(scoredMafLine_t *sml);
void destroyDuplicates(duplicate_t *d);
void destroyStringArray(char **sArray, int n);
scoredMafLine_t* newScoredMafLine(void);
unsigned longestLine(mafBlock_t *mb);
unsigned numberOfSequencesScoredMafLineList(scoredMafLine_t *m);
void reportBlockWithDuplicates(mafBlock_t *mb, duplicate_t *dupHead);
int cmp_by_score(const void *a, const void *b);

void parseOptions(int argc, char **argv, char *filename, char *seq, char *strand) {
    int c;
    bool setMaf = false;
    bool setSeq = false;
    while (1) {
        static struct option longOptions[] = {
            {"debug", no_argument, 0, 'd'},
            {"verbose", no_argument, 0, 'v'},
            {"help", no_argument, 0, 'h'},
            {"version", no_argument, 0, 0},
            {"maf",  required_argument, 0, 'm'},
            {"seq",  required_argument, 0, 0},
            {"strand",  required_argument, 0, 0},
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
            if (strcmp("seq", longOptions[longIndex].name) == 0) {
                setSeq = true;
                sscanf(optarg, "%s", seq);
                break;
            }
            if (strcmp("strand", longOptions[longIndex].name) == 0) {
                sscanf(optarg, "%c", strand);
                if (*strand != '+' && *strand != '-') {
                    fprintf(stderr, "--strand must be either + or -\n");
                    usage();
                }
                break;
            }
            break;
        case 'm':
            setMaf = true;
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
    if (!setMaf) {
        fprintf(stderr, "specify --maf\n");
        usage();
    }
    if (!setSeq) {
        fprintf(stderr, "specify --seq\n");
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
    fprintf(stderr, "mafBlockStrandCoercer, %s\nbuild: %s, %s, %s\n\n", g_version, g_build_date, 
            g_build_git_branch, g_build_git_sha);
}
void usage(void) {
    version();
    fprintf(stderr, "Usage: mafBlockStrandCoercer --maf alignment.maf --seq hg18 --strand + > positive.maf \n\n"
            "mafBlockStrandCoercer is a program to coerce a particular strandedness out for a block\n"
            "based the strandedness of a target sequence. If the block contains conflicing strands\n"
            "(i.e. both + and - strands are observed), then nothing is done.\n");
    fprintf(stderr, "Options: \n");
    usageMessage('h', "help", "show this help message and exit.");
    usageMessage('m', "maf", "input alignment maf file.");
    usageMessage('\0', "seq", "sequence to base block strandedness upon. (string comparison only done for length of input, i.e. --seq=hg18 will match hg18.chr1, hg18.chr2, etc etc)");
    usageMessage('\0', "strand", "strand to enforce, when possible. may be + or -, defaults to +.");
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
    d->tailScoredMaf = NULL;
    d->reported = false;
    d->next = NULL;
    d->numSequences = 1;
    return d;
}
void printHeader(void) {
    printf("##maf version=1\n\n");
}
void checkBlock(mafBlock_t *block, char *seq, char strand) {
    // read through each line of a mafBlock and check to see if a block needs to be reverse complemented.
    mafLine_t *ml = maf_mafBlock_getHeadLine(block);
    bool flipStrand = false;
    bool obsSeqPos = false; 
    bool obsSeqNeg = false; 
    size_t len = strlen(seq);
    while (ml != NULL) {
        if (maf_mafLine_getType(ml) != 's') {
            // skip non-sequence lines
            ml = maf_mafLine_getNext(ml);
            continue;
        }
        if (strncmp(maf_mafLine_getSpecies(ml), seq, len) == 0) {
            if (maf_mafLine_getStrand(ml) == '+') {
                obsSeqPos = true;
            } else {
                obsSeqNeg = true;
            }
            if (maf_mafLine_getStrand(ml) != strand) {
                flipStrand = true;
            }
        }
        ml = maf_mafLine_getNext(ml);
    }
    if (flipStrand && !(obsSeqPos && obsSeqNeg)) {
        maf_mafBlock_flipStrand(block);
    }
}
void processBody(mafFileApi_t *mfa, char *seq, char strand) {
    // walk the body of the maf file and process it, block by block.
    mafBlock_t *thisBlock = NULL;
    thisBlock = maf_readBlock(mfa); // header block, unused
    maf_destroyMafBlockList(thisBlock);
    printHeader();
    while((thisBlock = maf_readBlock(mfa)) != NULL) {
        checkBlock(thisBlock, seq, strand);
        maf_mafBlock_print(thisBlock);
        maf_destroyMafBlockList(thisBlock);
    }
}
int main(int argc, char **argv) {
    char filename[kMaxStringLength];
    char seq[kMaxStringLength];
    char strand = '+';
    parseOptions(argc, argv, filename, seq, &strand);
    mafFileApi_t *mfa = maf_newMfa(filename, "r");
    processBody(mfa, seq, strand);
    maf_destroyMfa(mfa);
    return EXIT_SUCCESS;
}
