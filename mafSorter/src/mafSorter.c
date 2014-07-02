/*
 * Copyright (C) 2011-2014 by
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
#include <assert.h>
#include <ctype.h>
#include <getopt.h>
#include <inttypes.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "common.h"
#include "sharedMaf.h"
#include "buildVersion.h"

const char *g_version = "version 0.1 September 2012";

typedef struct sortingMafBlock {
    // augmented data structure
    mafBlock_t *mafBlock; // pointer to actual mafBlock_t
    int64_t targetStart; // value to sort on, position in target sequence
} sortingMafBlock_t;

void version(void);
void usage(void);
void parseOptions(int argc, char **argv, char *filename, char *seqName);
int64_t max(int64_t a, int64_t b);
int64_t getTargetStartLine(mafLine_t *ml, char *targetSeq);
int64_t getTargetStartBlock(mafBlock_t *mb, char *targetSeq);
unsigned processBody(mafFileApi_t *mfa, mafBlock_t **head);
void populateArray(mafBlock_t *mb, sortingMafBlock_t **array, char *targetSequence);
int cmp_by_targetStart(const void *a, const void *b);
void reportBlock(sortingMafBlock_t *smb);
void reportBlocks(sortingMafBlock_t **array, unsigned numBlocks);
void destroyArray(sortingMafBlock_t **array, unsigned numBlocks);

void version(void) {
    fprintf(stderr, "mafSorter, %s\nbuild: %s, %s, %s\n\n", g_version, g_build_date,
            g_build_git_branch, g_build_git_sha);
}
void usage(void) {
    version();
    fprintf(stderr, "Usage: mafSorter --maf [maf file] --seq [sequence name (and possibly chr)] "
            "[options]\n\n"
            "mafSorter is a program that will sort the blocks of a maf in ascending\n"
            "order of the sequence start field of the specified sequence name. Blocks\n"
            "that do not contain the specified sequence will be output at the start of\n"
            "the maf in the order they appear in the input, followed by the sorted blocks.\n"
            "If a block has multiple instances of the target sequence, the largest starting\n"
            "position is used for that block.\n\n");
    fprintf(stderr, "Options: \n");
    usageMessage('h', "help", "show this help message and exit.");
    usageMessage('m', "maf", "path to the maf file.");
    usageMessage('s', "seq", "sequence name, e.g. `hg18.chr2'\n");
    usageMessage('v', "verbose", "turns on verbose output.");
    exit(EXIT_FAILURE);
}
void parseOptions(int argc, char **argv, char *filename, char *seqName) {
    extern int g_verbose_flag;
    extern int g_debug_flag;
    int c;
    bool setMName = false, setSName = false;
    while (1) {
        static struct option longOptions[] = {
            {"debug", no_argument, 0, 'd'},
            {"verbose", no_argument, 0, 'v'},
            {"help", no_argument, 0, 'h'},
            {"version", no_argument, 0, 0},
            {"maf",  required_argument, 0, 'm'},
            {"seq",  required_argument, 0, 's'},
            {0, 0, 0, 0}
        };
        int longIndex = 0;
        c = getopt_long(argc, argv, "d:v:h:m:s",
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
            setMName = true;
            strncpy(filename, optarg, kMaxSeqName);
            break;
        case 's':
            setSName = true;
            strncpy(seqName, optarg, kMaxSeqName);
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
    if (!(setMName && setSName)) {
        fprintf(stderr, "Error, specify --maf --seq\n");
        usage();
    }
    // Check there's nothing left over on the command line
    if (optind < argc) {
        char *errorString = de_malloc(kMaxSeqName);
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
int64_t max(int64_t a, int64_t b) {
    return (a > b ? a : b);
}
int64_t getTargetStartLine(mafLine_t *ml, char *targetSeq) {
    if (maf_mafLine_getType(ml) != 's')
        return INT64_MIN;
    if (strncmp(targetSeq, maf_mafLine_getSpecies(ml), strlen(targetSeq)) == 0) {
        return (int64_t) maf_mafLine_getPositiveLeftCoord(ml);
    }
    return INT64_MIN;
}
int64_t g_stableOrder = INT64_MIN;
int64_t getTargetStartBlock(mafBlock_t *mb, char *targetSeq) {
    mafLine_t *ml = maf_mafBlock_getHeadLine(mb);
    assert(ml != NULL);
    int64_t tStart = g_stableOrder++; // mac os x qsort is not a stable sort, impose a stable order
    while (ml != NULL) {
        tStart = max(tStart, getTargetStartLine(ml, targetSeq));
        ml = maf_mafLine_getNext(ml);
    }
    return tStart;
}
unsigned processBody(mafFileApi_t *mfa, mafBlock_t **head) {
    // process the body of the maf, block by block.
    *head = maf_readAll(mfa);
    return maf_mafBlock_getNumberOfBlocks(*head);
}
void populateArray(mafBlock_t *mb, sortingMafBlock_t **array, char *targetSequence) {
    // walk the linked list pointed to by head and stuff pointers to
    // the structs into the array.
    unsigned i = 0;
    while (mb != NULL) {
        de_debug("inserting block %2u: %s %" PRIi64 "\n", i,
                 maf_mafLine_getLine(maf_mafBlock_getHeadLine(mb)), getTargetStartBlock(mb, targetSequence));
        array[i] = (sortingMafBlock_t *) de_malloc(sizeof(sortingMafBlock_t));
        array[i]->mafBlock = mb;
        array[i++]->targetStart = getTargetStartBlock(mb, targetSequence);
        mb = maf_mafBlock_getNext(mb);
    }
}
int cmp_by_targetStart(const void *a, const void *b) {
    sortingMafBlock_t **ia = (sortingMafBlock_t **) a;
    sortingMafBlock_t **ib = (sortingMafBlock_t **) b;
    return ((*ia)->targetStart >= (*ib)->targetStart);
}
void reportBlock(sortingMafBlock_t *smb) {
    // print out the single block pointed to by mb
    mafLine_t *ml = maf_mafBlock_getHeadLine(smb->mafBlock);
    while(ml != NULL) {
        assert(maf_mafLine_getLine(ml) != NULL);
        printf("%s\n", maf_mafLine_getLine(ml));
        ml = maf_mafLine_getNext(ml);
    }
}
void reportBlocks(sortingMafBlock_t **array, unsigned numBlocks) {
    // look over the block array and print out all the blocks
    for (unsigned i = 0; i < numBlocks; ++i) {
        reportBlock(array[i]);
        printf("\n");
    }
}
void destroyArray(sortingMafBlock_t **array, unsigned numBlocks) {
    for (unsigned i = 0; i < numBlocks; ++i) {
        free(array[i]);
    }
}
int main(int argc, char **argv) {
    extern const int kMaxStringLength;
    char targetSequence[kMaxSeqName];
    char filename[kMaxStringLength];
    parseOptions(argc, argv, filename, targetSequence);

    mafFileApi_t *mfa = maf_newMfa(filename, "r");
    mafBlock_t *mb = NULL;
    unsigned numBlocks = processBody(mfa, &mb);
    sortingMafBlock_t *blockArray[numBlocks];
    populateArray(mb, blockArray, targetSequence);

    qsort(blockArray, numBlocks, sizeof(sortingMafBlock_t *), cmp_by_targetStart);
    reportBlocks(blockArray, numBlocks);
    destroyArray(blockArray, numBlocks);
    maf_destroyMfa(mfa);
    maf_destroyMafBlockList(mb);

    return(EXIT_SUCCESS);
}
