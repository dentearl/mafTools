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
#include "bioioC.h"
#include "sonLib.h"
#include "sharedMaf.h"

typedef struct _solo {
    char *name;
    uint32_t pos;
} ASolo;
typedef struct _pair {
    char *seq1;
    char *seq2;
    uint32_t pos1;
    uint32_t pos2;
} APair;
typedef struct _position {
    char *name;
    uint32_t pos;
} APosition;
typedef struct _resultPair {
    APair aPair;
    uint32_t inAll;
    uint32_t inBoth;
    uint32_t inA;
    uint32_t inB;
    uint32_t inNeither;
    uint32_t total;
    uint32_t totalBoth;
    uint32_t totalA;
    uint32_t totalB;
    uint32_t totalNeither;
} ResultPair;
bool g_isVerboseFailures;

void populateNames(const char *mAFFile, stHash *hash);
stSortedSet* compareMAFs_AB(const char *mAFFileA, const char *mAFFileB, uint32_t numberOfSamples, 
                            uint32_t *numberOfPairsInFile,
                            stHash *legitimateSequences, stHash *intervalsHash, 
                            uint32_t near);
void findentprintf(FILE *fp, unsigned indent, char const *fmt, ...);
void reportResults(stSortedSet *results_AB, const char *mAFFileA, const char *mAFFileB, 
                   FILE *fileHandle, uint32_t near, stHash *legitimateSequences, 
                   const char *bedFiles);
APair* aPair_init(void);
APair* aPair_construct(const char *seq1, const char *seq2, uint32_t pos1, uint32_t pos2);
APair* aPair_copyConstruct(APair *pair);
void aPair_destruct(APair *pair);
void writeXMLHeader( FILE *fileHandle );
bool* getLegitRows(char **names, uint32_t numSeqs, stHash *legitPairs);
uint64_t walkBlockCountingPairs(mafBlock_t *mb, stHash *legitPairs, uint64_t *chooseTwoArray);
uint64_t chooseTwo(uint32_t n);
uint64_t* buildChooseTwoArray(void);
uint64_t countPairsInMaf(const char *filename, stHash *legitPairs);
void pairIndicesToArrayIndex(uint64_t r, uint64_t c, uint64_t n, uint64_t *i);
void arrayIndexToPairIndices(uint64_t i, uint64_t n, uint64_t *r, uint64_t *c);
stHash* stHash_getIntersection(stHash *seqNames1, stHash *seqNames2);
void samplePairsFromColumn(double acceptProbability, stSortedSet *pairs, 
                           uint32_t numPairs, uint64_t *chooseTwoArray,
                           char **nameArray, uint32_t *columnPositions);
void samplePairsFromColumnBruteForce(double acceptProbability, stSortedSet *pairs, 
                                     uint64_t *chooseTwoArray,
                                     char **nameArray, uint32_t *positions, uint32_t numSeqs, 
                                     uint32_t numPairs);
void samplePairsFromColumnAnalytic(double acceptProbability, stSortedSet *pairs, 
                                   uint64_t *chooseTwoArray,
                                   char **nameArray, uint32_t *positions, uint32_t numSeqs, 
                                   uint32_t numPairs);
void samplePairsFromColumnNaive(char **mat, uint32_t c, bool *legitRows, double acceptProbability, 
                                stSortedSet *pairs, uint64_t *chooseTwoArray, 
                                char **nameArray, uint32_t *positions, uint32_t numSeqs, 
                                uint32_t numPairs);
uint32_t countLegitPositions(char **mat, uint32_t c, uint32_t numRows);
mafLine_t** cullMlArrayByColumn(char **mat, uint32_t c, mafLine_t **mlArray, bool *legitRows, uint32_t numRows, uint32_t numLegitGaplessPositions);
uint32_t* cullPositionsByColumn(char **mat, uint32_t c, uint32_t *positions, bool *legitRows, uint32_t numRows, uint32_t numLegitGaplessPositions);
void walkBlockSamplingPairs(mafBlock_t *mb, stSortedSet *pairs, double acceptProbability, stHash *legitSequences, uint64_t *chooseTwoArray);
int aPair_cmpFunction(APair *aPair1, APair *aPair2);
uint32_t sumBoolArray(bool *legitRows, uint32_t numSeqs);
mafLine_t** createMafLineArray(mafBlock_t *mb, uint32_t numLegit, bool *legitRows);
void updatePositions(char **mat, uint32_t c, uint32_t *positions, int *strandInts, uint32_t numSeqs);
void printSortedSet(stSortedSet *pairs);
#endif /* _COMPARATOR_API_H_ */
