/* 
 * Copyright (C) 2011-2012 by 
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
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "common.h"

int g_verbose_flag = 0;
int g_debug_flag = 0;
const int kMaxSeqName = 1 << 8;

typedef struct mafLine {
    // a mafLine struct is a single line of a mafBlock
    char *line;
    struct mafLine *next;
} mafLine_t;

typedef struct mafBlock {
    // a mafBlock struct contains a maf block as a linked list
    // and itself can be part of a mafBlock linked list.
    mafLine_t *head;
    mafLine_t *tail;
    uint32_t targetStart;
    struct mafBlock *next;
} mafBlock_t;

void usage(void) {
    fprintf(stderr, "Usage: mafBlockSorter --seq [sequence name (and possibly chr)] "
            "[options] < myFile.maf\n\n"
            "mafBlockSorter is a program that will sort the blocks of a maf in ascending\n"
            "order of the sequence start field of the specified sequence name. Blocks\n"
            "that do not contain the specified sequence will be output at the start of\n"
            "the maf in the order they appear in the input, followed by the sorted blocks.\n"
            "If a block has multiple instances of the target sequence, the largest starting\n"
            "position is used for that block.\n\n");
    fprintf(stderr, "Options: \n"
            "  -h, --help     show this help message and exit.\n"
            "  -s, --seq      sequence name.chr e.g. `hg18.chr2.'\n"
            "  -v, --verbose  turns on verbose output.\n");
    exit(EXIT_FAILURE);
}
void parseOptions(int argc, char **argv, char *seqName) {
    extern int g_verbose_flag;
    extern int g_debug_flag;
    int c;
    bool setSName = false;
    while (1) {
        static struct option long_options[] = {
            {"debug", no_argument, 0, 'd'},
            {"verbose", no_argument, 0, 'v'},
            {"help", no_argument, 0, 'h'},
            {"seq",  required_argument, 0, 's'},
            {0, 0, 0, 0}
        };
        int option_index = 0;
        c = getopt_long(argc, argv, "n:c:p:v",
                        long_options, &option_index);
        if (c == -1)
            break;
        switch (c) {
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
    if (!(setSName)) {
        fprintf(stderr, "Error, specify --seq\n");
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
bool blankLine(char *s) {
    // return true of line is only whitespaces
    size_t n = strlen(s);
    for (size_t i = 0; i < n; ++i) {
        if (!isspace(*(s + i))) {
            return false;
        }
    }
    return true;
}
uint32_t getTargetStart(char *line, char *targetSeq) {
    // read a maf sequence line and if the line contains
    // the target sequence, return the value of the start field
    // otherwise return 0.
    char *cline = (char *) de_malloc(strlen(line) + 1);
    strcpy(cline, line);
    char *tkn = NULL;
    tkn = strtok(cline, " \t");
    long unsigned start;
    if (tkn == NULL) {
        failBadFormat();
    }
    if (tkn[0] != 's') {
        free (cline);
        return 0;
    }
    tkn = strtok(NULL, " \t"); // name field
    if (tkn == NULL) {
        printf("cline: %s\n", cline);
        failBadFormat();
    }
    if (strncmp(targetSeq, tkn, strlen(targetSeq)) == 0) {
        tkn = strtok(NULL, " \t"); // start position
        if (tkn == NULL)
            failBadFormat();
        start = strtoul(tkn, NULL, 10);
        free(cline);
        return start;
    }
    free(cline);
    return 0;
}
uint32_t max(uint32_t a, uint32_t b) {
    return (a > b ? a : b);
}
mafBlock_t* newMafBlock(void) {
    mafBlock_t *mb = (mafBlock_t *) de_malloc(sizeof(*mb));
    mb->next = NULL;
    mb->head = NULL;
    mb->tail = NULL;
    mb->targetStart = 0;
    return mb;
}
mafLine_t* newMafLine(void) {
    mafLine_t* ml = (mafLine_t *) de_malloc(sizeof(*ml));
    ml->line = NULL;
    ml->next = NULL;
    return ml;
}
unsigned processBody(mafBlock_t *head, char *targetSeq, char *lastLine) {
    // process the body of the maf, block by block.
    extern const int kMaxStringLength;
    FILE *ifp = stdin;
    int32_t n = kMaxStringLength;
    char *line = NULL;
    if (lastLine == NULL) {
        line = (char*) de_malloc(n);
        if (de_getline(&line, &n, ifp) == -1) {
            fprintf(stderr, "Error reading first body line of alignment!\n");
            exit(EXIT_FAILURE);
        }
    } else {
        if (n < (int32_t) strlen(lastLine) + 1)
            n = (int32_t) strlen(lastLine) + 1;
        line = (char*) de_malloc(n);
        strcpy(line, lastLine);
    }
    mafBlock_t *thisBlock = head;
    unsigned numBlocks = 1;
    bool prevLineBlank = false;
    do {
        if (blankLine(line)) {
            if (prevLineBlank) {
                continue;
            }
            prevLineBlank = true;
        } else {
            if (prevLineBlank) {
                // if this line is not blank and the previous line was blank
                // then this is the start of a new maf block.
                ++numBlocks;
                mafBlock_t *nextBlock = newMafBlock();
                thisBlock->next = nextBlock;
                thisBlock = nextBlock;
            }
            prevLineBlank = false;
            char *copy = (char *) de_malloc(n + 1); // freed in destroy lines
            strcpy(copy, line);
            thisBlock->targetStart = max(thisBlock->targetStart, getTargetStart(copy, targetSeq));
            if (thisBlock->head == NULL) {
                // if thisBlock->head is null then this is a brand new mafBlock and it needs
                // a new mafLine to be created.
                thisBlock->head = newMafLine();
                thisBlock->head->line = copy;
                thisBlock->tail = thisBlock->head;
            } else {
                thisBlock->tail->next = newMafLine();
                thisBlock->tail->next->line = copy;
                thisBlock->tail = thisBlock->tail->next;
            }
        }
    } while (de_getline(&line, &n, ifp) != -1);
    free(line);
    return numBlocks;
}
void populateArray(mafBlock_t *head, mafBlock_t **array) {
    // walk the linked list pointed to by head and stuff pointers to
    // the structs into the array.
    unsigned i = 0;
    mafBlock_t *thisBlock = head;
    while (thisBlock != NULL) {
        array[i++] = thisBlock;
        thisBlock = thisBlock->next;
    }
}
int cmp_by_targetStart(const void *a, const void *b) {
    // mafBlock_t * const *ia = a;
    // mafBlock_t * const *ib = b;
    mafBlock_t **ia = (mafBlock_t **) a;
    mafBlock_t **ib = (mafBlock_t **) b;
    return ((*ia)->targetStart - (*ib)->targetStart);
}
void reportBlock(mafBlock_t *mb) {
    // print out the single block pointed to by mb
    mafLine_t *ml = mb->head;
    while(ml != NULL) {
        assert(ml->line != NULL);
        printf("%s\n", ml->line);
        ml = ml->next;
    }
}
void reportBlocks(mafBlock_t **array, unsigned numBlocks) {
    // look over the block array and print out all the blocks
    for (unsigned i = 0; i < numBlocks; ++i) {
        reportBlock(array[i]);
        printf("\n");
    }
}
void destroyLines(mafLine_t *ml) {
    mafLine_t *tmp;
    while(ml != NULL) {
        tmp = ml;
        ml = ml->next;
        free(tmp->line);
        free(tmp);
    }
}
void destroyBlocks(mafBlock_t *mb) {
    mafBlock_t *tmp;
    while(mb != NULL) {
        tmp = mb;
        mb = mb->next;
        destroyLines(tmp->head);
        free(tmp);
    }
}
void printblockarrayvalues(mafBlock_t **blockArray, unsigned numBlocks) {
    // debug
    for (unsigned i = 0; i < numBlocks; ++i) {
        printf("%d%s", blockArray[i]->targetStart, i == numBlocks - 1 ? "\n" : ", ");
    }
}
int main(int argc, char **argv) {
    char targetSequence[kMaxSeqName];
    mafBlock_t *mafObj = newMafBlock();
    unsigned numBlocks = 0;
    parseOptions(argc, argv, targetSequence);
    // initialize
    char *lastLine = processHeader();
    // read input
    numBlocks = processBody(mafObj, targetSequence, lastLine);
    mafBlock_t *blockArray[numBlocks];
    // sort
    populateArray(mafObj, blockArray);
    qsort(blockArray, numBlocks, sizeof(mafBlock_t *), cmp_by_targetStart);
    // write output
    reportBlocks(blockArray, numBlocks);
    // cleanup
    destroyBlocks(mafObj);
    free(lastLine);
    return EXIT_SUCCESS;
}
