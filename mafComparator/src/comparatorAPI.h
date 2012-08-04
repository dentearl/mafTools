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
#include <unistd.h> // getpid
#include "bioioC.h"
#include "sonLib.h"
#include "sharedMaf.h"

typedef struct _options {
    // used to hold all the command line options
    char *logLevelString;
    char *mafFile1;
    char *mafFile2;
    char *outputFile;
    char *bedFiles;
    char *wigglePairs;
    char *legitSequences; // the intersection of sequence names between inputs
    char *numPairsString;
    uint32_t numberOfSamples;
    uint32_t randomSeed;
    uint32_t near;
    uint64_t numPairs1;
    uint64_t numPairs2;
    uint64_t wiggleBinLength;
} Options;
typedef struct _pair {
    // used for sampling pairs of aligned positions
    char *seq1;
    char *seq2;
    uint32_t pos1;
    uint32_t pos2;
} APair;
typedef struct _position {
    // used in homology testing on columns
    char *name;
    uint32_t pos;
} APosition;
typedef struct _resultPair {
    // result pairs are used to store the disposition of sampled
    // pairs at the sequence name level (i.e. they have no position,
    // they are amalgams of all positions in a given comparison)
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
typedef struct _wiggleContainer {
    // contains arrays of counts, used by the --wigglePair option
    char *ref;
    char *partner;
    uint64_t refLength;
    uint64_t numBins;
    uint64_t binLength;
    uint64_t *presentAtoB;
    uint64_t *presentBtoA;
    uint64_t *absentAtoB;
    uint64_t *absentBtoA;
} WiggleContainer;
bool g_isVerboseFailures;

Options* options_construct(void);
void options_destruct(Options* o);
void populateNames(const char *mAFFile, stSet *set, stHash *seqLengthHash);
stSortedSet* compareMAFs_AB(const char *mAFFileA, const char *mAFFileB, uint64_t *numberOfPairsInFile, 
                            stSet *legitimateSequences, stHash *intervalsHash, stHash *wigHash, bool isAtoB, 
                            Options *options);
void findentprintf(FILE *fp, unsigned indent, char const *fmt, ...);
void reportResults(stSortedSet *results_AB, const char *mAFFileA, const char *mAFFileB, 
                   FILE *fileHandle, uint32_t near, stSet *legitimateSequences, 
                   const char *bedFiles);
APair* aPair_init(void);
APair* aPair_construct(const char *seq1, const char *seq2, uint32_t pos1, uint32_t pos2);
APair* aPair_copyConstruct(APair *pair);
void aPair_destruct(APair *pair);
WiggleContainer* wiggleContainer_init(void);
void wiggleContainer_destruct(WiggleContainer *wc);
void writeXMLHeader( FILE *fileHandle );
bool* getLegitRows(char **names, uint32_t numSeqs, stSet *legitPairs);
uint64_t walkBlockCountingPairs(mafBlock_t *mb, stSet *legitPairs, uint64_t *chooseTwoArray);
int32_t* buildInt(int32_t n);
int64_t* buildInt64(int64_t n);
uint64_t chooseTwo(uint64_t n);
uint64_t* buildChooseTwoArray(void);
uint64_t countPairsInMaf(const char *filename, stSet *legitPairs);
void pairIndicesToArrayIndex(uint64_t r, uint64_t c, uint64_t n, uint64_t *i);
void arrayIndexToPairIndices(uint64_t i, uint64_t n, uint64_t *r, uint64_t *c);
void samplePairsFromColumn(double acceptProbability, stSortedSet *sampledPairs, 
                           uint32_t numSeqs, uint64_t *chooseTwoArray,
                           char **nameArray, uint32_t *columnPositions);
void samplePairsFromColumnBruteForce(double acceptProbability, stSortedSet *sampledPairs, 
                                     uint64_t *chooseTwoArray,
                                     char **nameArray, uint32_t *positions, uint32_t numSeqs, 
                                     uint64_t numPairs);
void samplePairsFromColumnAnalytic(double acceptProbability, stSortedSet *sampledPairs, 
                                   uint64_t *chooseTwoArray,
                                   char **nameArray, uint32_t *positions, uint32_t numSeqs, 
                                   uint64_t numPairs);
void samplePairsFromColumnNaive(char **mat, uint32_t c, bool *legitRows, double acceptProbability, 
                                stSortedSet *sampledPairs, uint64_t *chooseTwoArray, 
                                char **nameArray, uint32_t *positions, uint32_t numSeqs, 
                                uint64_t numPairs);
uint32_t countLegitPositions(char **mat, uint32_t c, uint32_t numRows);
mafLine_t** cullMlArrayByColumn(char **mat, uint32_t c, mafLine_t **mlArray, bool *legitRows, uint32_t numRows, uint32_t numLegitGaplessPositions);
uint32_t* cullPositionsByColumn(char **mat, uint32_t c, uint32_t *positions, bool *legitRows, uint32_t numRows, uint32_t numLegitGaplessPositions);
void walkBlockSamplingPairs(mafBlock_t *mb, stSortedSet *sampledPairs, double acceptProbability, stSet *legitSequences, uint64_t *chooseTwoArray);
int aPair_cmpFunction(APair *aPair1, APair *aPair2);
uint32_t sumBoolArray(bool *legitRows, uint32_t numSeqs);
mafLine_t** createMafLineArray(mafBlock_t *mb, uint32_t numLegit, bool *legitRows);
void updatePositions(char **mat, uint32_t c, uint32_t *positions, int *strandInts, uint32_t numSeqs);
void printSortedSet(stSortedSet *pairs);
unsigned countCommas(char *s);
bool patternMatches(char *a, char *b);
void buildWigglePairHash(stHash *sequenceLengthHash, stList *wigglePairPatternList, 
                         stHash *wigglePairHash, uint64_t wiggleBinLength);
void reportResultsForWiggles(stHash *wigglePairHash, FILE *fileHandle);
void buildSeqNamesSet(Options *options, stSet *seqNamesSet, stHash *sequenceLengthHash);
#endif /* _COMPARATOR_API_H_ */
