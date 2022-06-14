/*
 * Copyright (C) 2009-2013 by
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
    uint64_t wiggleRegionStart;
    uint64_t wiggleRegionStop;
    char *legitSequences; // the intersection of sequence names between inputs
    char *numPairsString;
    uint64_t numberOfSamples;
    uint64_t randomSeed;
    uint64_t near;
    uint64_t numPairs1;
    uint64_t numPairs2;
    uint64_t wiggleBinLength;
} Options;
typedef struct _pair {
    // used for sampling pairs of aligned positions
    char *seq1;
    char *seq2;
    uint64_t pos1;
    uint64_t pos2;
} APair;
typedef struct _position {
    // used in homology testing on columns
    char *name;
    uint64_t pos;
} APosition;
typedef struct _resultPair {
    // result pairs are used to store the disposition of sampled
    // pairs at the sequence name level (i.e. they have no position,
    // they are amalgams of all positions in a given comparison)
    char *seq1;
    char *seq2;
    uint64_t inAll;
    uint64_t inBoth;
    uint64_t inA;
    uint64_t inB;
    uint64_t inNeither;
    uint64_t total;
    uint64_t totalBoth;
    uint64_t totalA;
    uint64_t totalB;
    uint64_t totalNeither;
} ResultPair;
typedef struct _wiggleContainer {
    // contains arrays of counts, used by the --wigglePair option
    char *ref;
    char *partner;
    uint64_t refStart; // starting base. i.e. --wiggleRegionStart
    uint64_t refLength;
    uint64_t numBins;
    uint64_t binLength;
    uint64_t *presentAtoB;
    uint64_t *presentBtoA;
    uint64_t *absentAtoB;
    uint64_t *absentBtoA;
} WiggleContainer;
extern bool g_isVerboseFailures;

Options* options_construct(void);
void options_destruct(Options* o);
void populateNames(const char *mAFFile, stSet *set, stHash *seqLengthHash);
stSortedSet* compareMAFs_AB(const char *mAFFileA, const char *mAFFileB, uint64_t *numberOfPairsInFile,
                            stSet *legitimateSequences, stHash *intervalsHash, stHash *wigHash, bool isAtoB,
                            Options *options, stHash *sequenceLengthHash);
void findentprintf(FILE *fp, unsigned indent, char const *fmt, ...);
void reportResults(stSortedSet *results_AB, const char *mAFFileA, const char *mAFFileB,
                   FILE *fileHandle, uint64_t near, stSet *legitimateSequences,
                   const char *bedFiles);
APair* aPair_init(void);
void aPair_fillOut(APair *aPair, char *seq1, char *seq2, uint64_t pos1, uint64_t pos2);
APair* aPair_construct(const char *seq1, const char *seq2, uint64_t pos1, uint64_t pos2);
APair* aPair_copyConstruct(APair *pair);
void aPair_destruct(APair *pair);
WiggleContainer* wiggleContainer_init(void);
WiggleContainer* wiggleContainer_construct(char *ref, char *partner, uint64_t refStart,
                                           uint64_t refLength, uint64_t wiggleBinLength);
ResultPair *resultPair_construct(const char *seq1, const char *seq2);
void aPosition_fillOut(APosition *aPosition, char *name, uint64_t pos);
APosition* aPosition_init(void);
APosition* aPosition_construct(const char *name, uint64_t pos);
stHash* constructPositionHash(char **mat, uint64_t c, char **names, uint64_t numSeqs,
                              uint64_t *allPositions, bool *legitRows);
void aPosition_destruct(void *p);
void resultPair_destruct(ResultPair *rp);
int aPair_cmpFunction_seqsOnly(APair *aPair1, APair *aPair2);
int strcmpnull(const char *s1, const char *s2);
uint64_t aPositionKey(const void *p);
int aPositionEqualKey(const void *k1, const void *k2);
bool closeEnough(uint64_t p1, uint64_t p2, uint64_t near);
void wiggleContainer_destruct(WiggleContainer *wc);
void writeXMLHeader( FILE *fileHandle );
bool* getLegitRows(char **names, uint64_t numSeqs, stSet *legitPairs);
uint64_t walkBlockCountingPairs(mafBlock_t *mb, stSet *legitPairs, uint64_t *chooseTwoArray);
int64_t* buildInt(int64_t n);
int64_t* buildInt64(int64_t n);
uint64_t* buildUInt64(uint64_t n);
uint64_t* uint64Copy(uint64_t *i);
uint64_t chooseTwo(uint64_t n);
uint64_t* buildChooseTwoArray(void);
uint64_t countPairsInMaf(const char *filename, stSet *legitPairs);
uint64_t countPairsInColumn(char **mat, uint64_t c, uint64_t numSeqs, bool *legitRows, uint64_t *chooseTwoArray);
uint64_t countLegitGaplessPositions(char **mat, uint64_t c, uint64_t numRows, bool *legitRows);
void countPairs(APair *pair, stHash *intervalsHash, int64_t *counter,
                stSortedSet *legitPairs, void *a, uint64_t near);

void pairIndicesToArrayIndex(uint64_t r, uint64_t c, uint64_t n, uint64_t *i);
void arrayIndexToPairIndices(uint64_t i, uint64_t n, uint64_t *r, uint64_t *c);
void samplePairs(APair *thisPair, stHash *intervalsHash, stSortedSet *pairs,
                 double *acceptProbability, stHash *legitPairs, uint64_t near);
bool inInterval(stHash *intervalsHash, char *seq, uint64_t pos);
uint64_t findLowerBound(uint64_t pos, uint64_t near);
void recordNearPair(APair *thisPair, stSortedSet *sampledPairs, uint64_t near, stSet *positivePairs);
void samplePairsFromMaf(const char *filename, stSortedSet *pairs, double acceptProbability,
                        stSet *legitSequences, uint64_t *numPairs, stHash *sequenceLengthHash);
void samplePairsFromColumn(double acceptProbability, stSortedSet *sampledPairs,
                           uint64_t numSeqs, uint64_t *chooseTwoArray,
                           char **nameArray, uint64_t *columnPositions);
void samplePairsFromColumnBruteForce(double acceptProbability, stSortedSet *sampledPairs,
                                     uint64_t *chooseTwoArray,
                                     char **nameArray, uint64_t *positions, uint64_t numSeqs,
                                     uint64_t numPairs);
void samplePairsFromColumnAnalytic(double acceptProbability, stSortedSet *sampledPairs,
                                   uint64_t *chooseTwoArray,
                                   char **nameArray, uint64_t *positions, uint64_t numSeqs,
                                   uint64_t numPairs);
void samplePairsFromColumnNaive(char **mat, uint64_t c, bool *legitRows, double acceptProbability,
                                stSortedSet *sampledPairs, uint64_t *chooseTwoArray,
                                char **nameArray, uint64_t *positions, uint64_t numSeqs,
                                uint64_t numPairs);
void walkBlockTestingHomology(mafBlock_t *mb, stSortedSet *sampledPairs, stSet *positivePairs,
                              stSet *legitSequences, stHash *intervalsHash, uint64_t near);
void testHomologyOnColumn(char **mat, uint64_t c, uint64_t numSeqs, bool *legitRows, char **names,
                          stSortedSet *sampledPairs, stSet *positivePairs, mafLine_t **mlArray,
                          uint64_t *allPositions, stHash *intervalsHash, uint64_t near);
void performHomologyTests(const char *filename, stSortedSet *sampledPairs, stSet *positivePairs,
                          stSet *legitSequences, stHash *intervalsHash, uint64_t near);
void homologyTests1(APair *thisPair, stHash *intervalsHash, stSortedSet *pairs,
                    stSet *positivePairs, stSet *legitPairs, int64_t near);
void enumerateHomologyResults(stSortedSet *sampledPairs, stSortedSet *resultPairs, stHash *intervalsHash,
                              stSet *positivePairs, stHash *wigglePairHash, bool isAtoB,
                              uint64_t wiggleBinLength);
ResultPair *aggregateResult(void *(*getNextPair)(void *, void *), stSortedSet *set, void *seqName,
                            const char *name1, const char *name2);
void* addReferencesAndDups_getDups(void *iterator, void *seqName);
void* addReferencesAndDups_getReferences(void *iterator, void *seqName);
void addReferencesAndDups(stSortedSet *results_AB, stSet *legitSequences);
void* reportResults_fn(void *iterator, void *seqName);
void reportResultsForWigglesArrays(FILE *fileHandle, WiggleContainer *wc);
void reportResult(const char *tagName, double total, double totalTrue, FILE *fileHandle, unsigned tabLevel);
uint64_t* cullPositions(uint64_t *allPositions, uint64_t numSeqs, bool *legitRows, uint64_t numLegit);
int* cullStrandInts(int *allStrandInts, uint64_t numSeqs, bool *legitRows, uint64_t numLegit);
char** extractLegitGaplessNamesFromMlArrayByColumn(char **mat, uint64_t c, mafLine_t **mlArray, bool *legitRows,
                                                   uint64_t numRows, uint64_t numLegitGaplessPositions);
void validateMafBlockSourceLengths(const char *filename, mafBlock_t *mb, stHash *sequenceLengthHash);
bool pairMemberInPositionHash(stHash *positionHash, APair *thisPair);
void printmlarray(mafLine_t **mlArray, uint64_t n);
void printHash(stHash *hash);
void printAllPositions(uint64_t *allPositions, mafBlock_t *mb);
void printAllStrandInts(int *allStrandInts, mafBlock_t *mb);
uint64_t countLegitPositions(char **mat, uint64_t c, uint64_t numRows);
mafLine_t** cullMlArrayByColumn(char **mat, uint64_t c, mafLine_t **mlArray, bool *legitRows, uint64_t numRows, uint64_t numLegitGaplessPositions);
uint64_t* cullPositionsByColumn(char **mat, uint64_t c, uint64_t *positions, bool *legitRows, uint64_t numRows, uint64_t numLegitGaplessPositions);
void walkBlockSamplingPairs(const char *filename, mafBlock_t *mb, stSortedSet *sampledPairs, double acceptProbability, stSet *legitSequences, uint64_t *chooseTwoArray, uint64_t *numPairs, stHash *sequenceLengthHash);
int aPair_cmpFunction(APair *aPair1, APair *aPair2);
uint64_t sumBoolArray(bool *legitRows, uint64_t numSeqs);
mafLine_t** createMafLineArray(mafBlock_t *mb, uint64_t numLegit, bool *legitRows);
void updatePositions(char **mat, uint64_t c, uint64_t *positions, int *strandInts, uint64_t numSeqs);
void printSortedSet(stSortedSet *pairs);
unsigned countChars(char *s, char c);
bool patternMatches(char *a, char *b);
void buildWigglePairHash(stHash *sequenceLengthHash, stList *wigglePairPatternList,
                         stHash *wigglePairHash, uint64_t wiggleBinLength, uint64_t wiggleRegionStart,
                         uint64_t wiggleRegionStop);
void reportResultsForWiggles(stHash *wigglePairHash, FILE *fileHandle);
void buildSeqNamesSet(Options *options, stSet *seqNamesSet, stHash *sequenceLengthHash);
bool positionIsInWiggleRegion(WiggleContainer *wc, uint64_t *refPos);
#endif /* _COMPARATOR_API_H_ */
