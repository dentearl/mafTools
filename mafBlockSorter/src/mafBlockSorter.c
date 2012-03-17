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
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "common.h"

int g_verbose_flag = 0;
int g_debug_flag = 0;
const int kMaxSeqName = 1<< 8;

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
    int targetStart;
    struct mafBlock *next;
} mafBlock_t;

void parseOptions(int argc, char **argv, char *seqName) {
    extern int g_verbose_flag;
    extern int g_debug_flag;
    strcpy(seqName, "hg18.chr7");
    argc++;
    argv[0][0] = seqName[0];
    g_verbose_flag = g_debug_flag;
}
void processHeader(void) {
    FILE *ifp = stdin;
    int32_t n = 1 << 10;
    char *line = (char*) de_malloc(n);
    int status = de_getline(&line, &n, ifp);
    if (status == -1) {
        fprintf(stderr, "Error, empty file\n");
        exit(EXIT_FAILURE);
    }
    if (strncmp(line, "##maf", 5) != 0) {
        fprintf(stderr, "Error, bad maf format. File should start with ##maf\n");
        exit(EXIT_FAILURE);
    }
    printf("%s\n", line);
    while (*line != 0 && status != -1) {
        status = de_getline(&line, &n, ifp);
        printf("%s\n", line);
    }
}
bool blankLine(char *s) {
    size_t n = strlen(s);
    for (size_t i = 0; i < n; ++i) {
        if (!isspace(*(s + i))) {
            return false;
        }
    }
    return true;
}
void badFormat(void) {
    fprintf(stderr, "The maf sequence lines are incorrectly formatted, exiting\n");
    exit(EXIT_FAILURE);
}
int getTargetStart(char *line, char *targetSeq) {
    char *cline = (char*) de_malloc(strlen(line) + 1);
    strcpy(cline, line);
    char *tkn = NULL;
    tkn = strtok(cline, " \t");
    if (tkn == NULL) {
        badFormat();
    }
    if (tkn[0] != 's') {
        free (cline);
        return 0;
    }
    tkn = strtok(NULL, " \t"); // name field
    if (tkn == NULL) {
        printf("cline: %s\n", cline);
        badFormat();
    }
    if (strncmp(targetSeq, tkn, strlen(targetSeq)) == 0) {
        tkn = strtok(NULL, " \t"); // start position
        if (tkn == NULL)
            badFormat();
        free(cline);
        return strtoul(tkn, NULL, 10);
    }
    free(cline);
    return 0;
}
int max(int a, int b) {
    return (a > b ? a : b);
}
unsigned processBody(mafBlock_t *head, char *targetSeq) {
    FILE *ifp = stdin;
    int32_t n = 1 << 10;
    char *line = (char*) de_malloc(n);
    mafBlock_t *thisBlock = head;
    unsigned numBlocks = 1;
    bool prevLineBlank = false;
    while (de_getline(&line, &n, ifp) != -1) {
        if (blankLine(line)) {
            if (prevLineBlank) {
                continue;
            }
            prevLineBlank = true;
        } else {
            if (prevLineBlank) {
                ++numBlocks;
                mafBlock_t *nextBlock = (mafBlock_t*) de_malloc(sizeof(mafBlock_t));
                nextBlock->next = NULL;
                nextBlock->head = NULL;
                nextBlock->tail = NULL;
                nextBlock->targetStart = 0;
                thisBlock->next = nextBlock;
                thisBlock = nextBlock;
            }
            prevLineBlank = false;
            char *copy = (char*) de_malloc(n + 1);
            strcpy(copy, line);
            thisBlock->targetStart = max(thisBlock->targetStart, getTargetStart(copy, targetSeq));
            if (thisBlock->head == NULL) {
                thisBlock->head = (mafLine_t*) de_malloc(sizeof(mafLine_t));
                thisBlock->head->next = NULL;
                thisBlock->head->line = copy;
                thisBlock->tail = thisBlock->head;
            } else {
                thisBlock->tail->next = (mafLine_t*) de_malloc(sizeof(mafLine_t));
                thisBlock->tail->next->next = NULL;
                thisBlock->tail->next->line = copy;
                thisBlock->tail = thisBlock->tail->next;
            }
        }
    }
    return numBlocks;
}
void populateArray(mafBlock_t *head, mafBlock_t **array) {
    // walk the linked list pointed to by head and stuff pointers to
    // the structs into the array.
    int i = 0;
    mafBlock_t *thisBlock = head;
    while(thisBlock != NULL) {
        array[i++] = thisBlock;
        thisBlock = thisBlock->next;
    }
}
int cmp_by_targetStart(const void *a, const void *b) {
    // mafBlock_t * const *ia = a;
    // mafBlock_t * const *ib = b;
    mafBlock_t **ia = (mafBlock_t**) a;
    mafBlock_t **ib = (mafBlock_t**) b;
    return ((*ia)->targetStart - (*ib)->targetStart);
}
void reportBlock(mafBlock_t *mb) {
    mafLine_t *ml = mb->head;
    while(ml != NULL) {
        assert(ml->line != NULL);
        printf("%s\n", ml->line);
        ml = ml->next;
    }
}
void reportBlocks(mafBlock_t **array, unsigned numBlocks) {
    for (unsigned i = 0; i < numBlocks; ++i) {
        reportBlock(array[i]);
        printf("\n");
    }
}
void destroyLines(mafLine_t *head) {
    mafLine_t *ml = head;
    while(ml != NULL) {
        free(ml->line);
        ml = ml->next;
    }
}
void destroyBlocks(mafBlock_t *head) {
    mafBlock_t *mb = head;
    while(mb != NULL) {
        destroyLines(mb->head);
        mb = mb->next;
    }
}
void printblockarrayvalues(mafBlock_t **blockArray, unsigned numBlocks) {
    for (unsigned i = 0; i < numBlocks; ++i) {
        printf("%d%s", blockArray[i]->targetStart, i == numBlocks - 1 ? "\n" : ", ");
    }
}
int main(int argc, char **argv) {
    char targetSequence[kMaxSeqName];
    mafBlock_t mafObj;
    unsigned numBlocks = 0;
    parseOptions(argc, argv, targetSequence);
    processHeader();
    
    mafObj.next = NULL;
    mafObj.head = NULL;
    mafObj.tail = NULL;
    mafObj.targetStart = 0;
    numBlocks = processBody(&mafObj, targetSequence);
    mafBlock_t *blockArray[numBlocks];
    populateArray(&mafObj, blockArray);
    qsort(blockArray, numBlocks, sizeof(mafBlock_t*), cmp_by_targetStart);
    reportBlocks(blockArray, numBlocks);
    destroyBlocks(&mafObj);
    
    return EXIT_SUCCESS;
}
