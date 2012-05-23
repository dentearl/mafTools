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
#ifndef MAFTRANSITIVECLOSURE_H_
#define MAFTRANSITIVECLOSURE_H_
#include "sonLib.h"
#include "stPinchGraphs.h"
#include "common.h"
#include "CuTest.h"
#include "sharedMaf.h"

typedef struct mafTcSeq {
    // maf tc (trasitive closure) sequence
    char *name;
    char *sequence;
    uint32_t length;
} mafTcSeq_t;
typedef struct mafTcRegion {
    // region or interval
    uint32_t start;
    uint32_t end;
    struct mafTcRegion *next;
} mafTcRegion_t;
typedef struct mafTcComparisonOrder {
    // This struct is used in adding alignment information into the thread set.
    /* for this sequence matrix:
         0123456789
       0 AC---ACG-G
       1 ACG--ACGGC
       2 A-G-TACGGC
       3 ACGTTACGGC
       a comparison order is {3: [3, 3]}, {2: [4, 4]},
       {1: [8, 8]}, {1: [2, 2]}, {0: [9, 9]}, {0: [5, 7]}, {0: [0, 1]}
       [note that the ordering of these values is arbitrary, though consistent with the
       output of the algo in the code which appends new structs to the head of the list]
       e.g. the first block to process uses 0 as its reference and it starts at column 1 
       and ends at column 1. The second block to process uses 1 as its reference and it
       starts at column 2 and ends at column 2 (it is only one column), etc.
     */
    uint32_t ref; // 
    mafTcRegion_t *region;
    struct mafTcComparisonOrder *next;
} mafTcComparisonOrder_t;

void usage(void);
mafTcSeq_t* newMafTcSeq(char *name, unsigned length);
mafTcComparisonOrder_t* newMafTcComparisonOrder(void);
mafTcRegion_t* newMafTcRegion(uint32_t start, uint32_t end);
void destroyMafTcSeq(void *p);
void destroyMafTcRegionList(mafTcRegion_t *r);
void destroyMafTcRegion(mafTcRegion_t *r);
void destroyMafTcComparisonOrder(mafTcComparisonOrder_t *c);
uint32_t hashMafTcSeq(const mafTcSeq_t *mtcs);
int hashCompareMafTcSeq(const mafTcSeq_t *m1, const mafTcSeq_t *m2);
char* createNSequence(unsigned length);
void reverseComplementSequence(char *s);
void complementSequence(char *s);
char complementChar(char c);
void addSequenceValuesToMtcSeq(mafLine_t *ml, mafTcSeq_t *mtcs);
void parseOptions(int argc, char **argv, char *filename);
stPinchThreadSet* buildThreadSet(stHash *hash);
void walkBlockAddingAlignments(mafBlock_t *mb, stPinchThreadSet *threadSet);
void addAlignmentsToThreadSet(mafFileApi_t *mfa, stPinchThreadSet *threadSet);
void createSequenceHash(mafFileApi_t *mfa, stHash **hash, stHash **nameHash);
mafTcRegion_t* getComparisonOrderFromRow(char **mat, uint32_t row, mafTcComparisonOrder_t **done, mafTcRegion_t *todo);
mafTcComparisonOrder_t *getComparisonOrderFromMatrix(char **mat, uint32_t rowLength, uint32_t colLength);
void processPairForPinching(stPinchThreadSet *threadSet, stPinchThread *a, uint32_t aGlobalStart, 
                            uint32_t aGlobalLength, int aStrand, 
                            char *aSeq, stPinchThread *b, uint32_t bGlobalStart, uint32_t bGlobalLength,
                            int bStrand, char *bSeq, uint32_t regionStart, uint32_t regionEnd);
int64_t localSeqCoords(uint32_t p, char *s);
int64_t localSeqCoordsToGlobalPositiveCoords(int64_t c, uint32_t start, uint32_t sourceLength, char strand);
int64_t localSeqCoordsToGlobalPositiveStartCoords(int64_t c, uint32_t start, uint32_t sourceLength, 
                                                  char strand, uint32_t length);
// test suite
CuSuite* mafTransitiveClosure_TestSuite(void);
int mafTransitiveClosure_RunAllTests(void);

#endif // MAFTRANSITIVECLOSURE_H_
