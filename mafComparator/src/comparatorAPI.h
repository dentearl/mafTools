/* 
 * Copyright (C) 2009-2012 by 
 * Dent Earl (dearl@soe.ucsc.edu, dentearl@gmail.com)
 * Benedict Paten (benedict@soe.ucsc.edu, benedictpaten@gmail.com)
 * Mark Diekhans (markd@soe.ucsc.edu)
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
#ifndef _COMPARATOR_API_H_
#define _COMPARATOR_API_H_

#include <assert.h>
#include <limits.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <getopt.h>
#include "avl.h"
#include "bioioC.h"
#include "commonC.h"
#include "hashTableC.h"
#include "hashTableC_itr.h"
#include "sonLib.h"
#include "sharedMaf.h"
#include "disjointset.h"

typedef struct _solo {
    char *name;
    int32_t pos;
} ASolo;
typedef struct _pair {
    char *seq1;
    char *seq2;
    uint32_t pos1;
    uint32_t pos2;
} APair;
typedef struct _resultPair {
    APair aPair;
    int32_t inAll;
    int32_t inBoth;
    int32_t inA;
    int32_t inB;
    int32_t inNeither;
    int32_t total;
    int32_t totalBoth;
    int32_t totalA;
    int32_t totalB;
    int32_t totalNeither;
} ResultPair;

void populateNames(const char *mAFFile, stHash *hash);
void intersectHashes(struct hashtable *h1, struct hashtable *h2, struct hashtable *h3);
void printNameHash(struct hashtable *h);
struct avl_table *compareMAFs_AB(const char *mAFFileA, const char *mAFFileB, uint32_t numberOfSamples, 
                                 stHash *legitimateSequences, stHash *intervalsHash, 
                                 int32_t verbose, uint32_t near);
void reportResults(struct avl_table *results_AB, const char *mAFFileA, const char *mAFFileB, 
                   FILE *fileHandle, uint32_t near, stHash *legitimateSequences, const char *bedFiles);
void aPair_destruct(APair *pair, void *extraArgument);
void writeXMLHeader( FILE *fileHandle );
bool* getLegitRows(char **names, uint32_t numSeqs, stHash *legitPairs);
uint64_t walkBlockCountingPairs(mafBlock_t *mb, stHash *legitPairs, uint64_t *chooseTwoArray);
uint64_t chooseTwo(uint32_t n);
uint64_t* buildChooseTwoArray(void);
uint64_t countPairsInMaf(const char *filename, stHash *legitPairs);
void pairIndicesToArrayIndex(uint64_t r, uint64_t c, uint64_t n, uint64_t *i);
void arrayIndexToPairIndices(uint64_t i, uint64_t n, uint64_t *r, uint64_t *c);
stHash* stHash_getIntersection(stHash *seqNames1, stHash *seqNames2);
void samplePairsFromColumn(char **mat, uint32_t c, bool *legitRows, 
                           char **names, double acceptProbability, struct avl_table *pairs, 
                           uint32_t numSeqs, uint64_t *chooseTwoArray,
                           mafLine_t **mlArray, uint32_t *positions);
void samplePairsFromColumnBruteForce(char **mat, uint32_t c, bool *legitRows, 
                                     char **names, double acceptProbability, struct avl_table *pairs, 
                                     uint32_t numSeqs, uint64_t *chooseTwoArray,
                                     mafLine_t **mlArray, uint32_t *positions, uint64_t numPairs);
void samplePairsFromColumnAnalytic(char **mat, uint32_t c, bool *legitRows, 
                                   char **names, double acceptProbability, struct avl_table *pairs, 
                                   uint32_t numSeqs, uint64_t *chooseTwoArray,
                                   mafLine_t **mlArray, uint32_t *positions, uint64_t numPairs);

#endif /* _COMPARATOR_API_H_ */
