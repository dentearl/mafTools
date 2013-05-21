/* 
 * Copyright (C) 2013 by 
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
#ifndef TEST_MAFTRANSITIVECLOSURE_H_
#define TEST_MAFTRANSITIVECLOSURE_H_
#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include "CuTest.h"
#include "common.h"
#include "sonLib.h"
#include "stPinchGraphs.h"
#include "mafTransitiveClosure.h"

void printRegionList(mafTcRegion_t *reg, FILE *ofp);
bool regionListsAreEqual(mafTcRegion_t *expected, mafTcRegion_t *obs, bool verbose);
void printTestComparisonOrder(mafTcComparisonOrder_t *co);
bool comparisonOrdersAreEqual(mafTcComparisonOrder_t *eo, mafTcComparisonOrder_t *oo, bool verbose);
bool mafBlocksAreEqual(mafBlock_t *input, mafBlock_t *expected, bool verbose);
void test_reverseComplement(CuTest *testCase);
void test_rowAlignmentBlockComparisonOrdering_0(CuTest *testCase);
void test_rowAlignmentBlockComparisonOrdering_1(CuTest *testCase);
void test_rowAlignmentBlockComparisonOrdering_2(CuTest *testCase);
void test_rowAlignmentBlockComparisonOrdering_3(CuTest *testCase);
void test_matrixAlignmentBlockComparisonOrdering_0(CuTest *testCase);
void test_matrixAlignmentBlockComparisonOrdering_1(CuTest *testCase);
void test_matrixAlignmentBlockComparisonOrdering_2(CuTest *testCase);
void test_matrixAlignmentBlockComparisonOrdering_3(CuTest *testCase);
void test_matrixAlignmentBlockComparisonOrdering_4(CuTest *testCase);
void test_addSequenceValuesToMtcSeq_0(CuTest *testCase);
void test_localSeqCoords_0(CuTest *testCase);
void test_localSeqCoordsToGlobalPositiveCoords_0(CuTest *testCase);
void test_localSeqCoordsToGlobalPositiveStartCoords_0(CuTest *testCase);
void test_coordinateTransforms_0(CuTest *testCase);
void test_coordinateTransforms_1(CuTest *testCase);
void test_mafBlockGapSorting_0(CuTest *testCase);
CuSuite* mafTransitiveClosure_TestSuite(void);

#endif // TEST_MAFTRANSITIVECLOSURE_H_
