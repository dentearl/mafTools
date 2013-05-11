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
#include <assert.h>
#include <inttypes.h>
#include <stdarg.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include "CuTest.h"
#include "common.h"
#include "sharedMaf.h"
#include "mafPairCoverageAPI.h"

static void test_is_wild_0(CuTest *testCase) {
    CuAssertTrue(testCase, is_wild("hg19*"));
    CuAssertTrue(testCase, is_wild("hg19.chr19*"));
    CuAssertTrue(testCase, !is_wild("hg19.chr19"));
    CuAssertTrue(testCase, !is_wild("aoeuaoeunstaoeunshtonuts.chrcrhrc.huaoeunsatohunt."));
}
static void test_searchMatched_0(CuTest *testCase) {
    mafLine_t *ml = maf_newMafLineFromString("s hg19.chr19        123480 13 + 1234870098734 ACGTACGTACGTA", 1);
    CuAssertTrue(testCase, searchMatched(ml, "hg19.chr19"));
    CuAssertTrue(testCase, searchMatched(ml, "hg19*"));
    CuAssertTrue(testCase, searchMatched(ml, "h*"));
    CuAssertTrue(testCase, searchMatched(ml, "*"));
    CuAssertTrue(testCase, !searchMatched(ml, "mm9"));
    maf_destroyMafLineList(ml);
}
static void test_compareLines_0(CuTest *testCase) {
    // just make sure that its counting the number of aligned positions correctly
    // test case 0
    mafLine_t *ml1 = maf_newMafLineFromString("s hg19.chr19      123480 13 + 1234870098734 ACGTACGTACGTA", 1);
    mafLine_t *ml2 = maf_newMafLineFromString("s mm9.chr1        123480 13 + 1234870098735 ACGTACGTACGTA", 1);
    stHash *seq1Hash = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, free, free);
    stHash *seq2Hash = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, free, free);
    uint64_t alignedPositions = 0;
    mafCoverageCount_t *mcct1 = createMafCoverageCount();
    mafCoverageCount_t *mcct2 = createMafCoverageCount();
    mafCoverageCount_setSourceLength(mcct1, 1234870098734);
    mafCoverageCount_setSourceLength(mcct2, 1234870098735);
    stHash_insert(seq1Hash, stString_copy("hg19.chr19"), mcct1);
    stHash_insert(seq2Hash, stString_copy("mm9.chr1"), mcct2);
    compareLines(ml1, ml2, seq1Hash, seq2Hash, &alignedPositions);
    CuAssertTrue(testCase, alignedPositions == 13);
    compareLines(ml1, ml2, seq1Hash, seq2Hash, &alignedPositions);
    CuAssertTrue(testCase, alignedPositions == 26);
    maf_destroyMafLineList(ml1);
    maf_destroyMafLineList(ml2);
    stHash_destruct(seq1Hash);
    stHash_destruct(seq2Hash);
    
    // test case 1
    alignedPositions = 0;
    ml1 = maf_newMafLineFromString("s hg19.chr19      123480 13 + 1234870098734 ACGTACGTACGTA", 1);
    ml2 = maf_newMafLineFromString("s mm9.chr2        123480  5 + 1234870098735 AC--------GTA", 1);
    seq1Hash = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, free, free);
    seq2Hash = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, free, free);
    mcct1 = createMafCoverageCount();
    mcct2 = createMafCoverageCount();
    mafCoverageCount_setSourceLength(mcct1, 1234870098734);
    mafCoverageCount_setSourceLength(mcct2, 1234870098735);
    stHash_insert(seq1Hash, stString_copy("hg19.chr19"), mcct1);
    stHash_insert(seq2Hash, stString_copy("mm9.chr2"), mcct2);
    compareLines(ml1, ml2, seq1Hash, seq2Hash, &alignedPositions);
    CuAssertTrue(testCase, alignedPositions == 5);
    maf_destroyMafLineList(ml1);
    maf_destroyMafLineList(ml2);
    stHash_destruct(seq1Hash);
    stHash_destruct(seq2Hash);
}
static void test_compareLines_1(CuTest *testCase) {
    // make sure that the hash is being correctly populated
    // test case 0
    mafLine_t *ml1 = maf_newMafLineFromString("s hg19.chr19      123480 13 + 1234870098734 ACGTACGTACGTA", 1);
    mafLine_t *ml2 = maf_newMafLineFromString("s mm9.chr1        123480 13 + 1234870098734 ACGTACGTACGTA", 1);
    stHash *seq1Hash = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, free, free);
    stHash *seq2Hash = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, free, free);
    uint64_t alignedPositions = 0;
    mafCoverageCount_t *mcct1 = createMafCoverageCount();
    mafCoverageCount_t *mcct2 = createMafCoverageCount();
    mafCoverageCount_setSourceLength(mcct1, 1234870098734);
    mafCoverageCount_setSourceLength(mcct2, 1234870098734);
    stHash_insert(seq1Hash, stString_copy("hg19.chr19"), mcct1);
    stHash_insert(seq2Hash, stString_copy("mm9.chr1"), mcct2);
    compareLines(ml1, ml2, seq1Hash, seq2Hash, &alignedPositions);
    CuAssertTrue(testCase, stHash_search(seq1Hash, "hg19.chr19") != NULL);
    CuAssertTrue(testCase, stHash_search(seq2Hash, "mm9.chr1") != NULL);
    CuAssertTrue(testCase, stHash_search(seq2Hash, "bannana") == NULL);
    maf_destroyMafLineList(ml1);
    maf_destroyMafLineList(ml2);
    stHash_destruct(seq1Hash);
    stHash_destruct(seq2Hash);
    // test case 1
}
CuSuite* pairCoverage_TestSuite(void) {
    CuSuite* suite = CuSuiteNew();
    (void) test_is_wild_0;
    (void) test_searchMatched_0;
    (void) test_compareLines_0;
    (void) test_compareLines_1;
    SUITE_ADD_TEST(suite, test_is_wild_0);
    SUITE_ADD_TEST(suite, test_searchMatched_0);
    SUITE_ADD_TEST(suite, test_compareLines_0);
    SUITE_ADD_TEST(suite, test_compareLines_1);
    return suite;
}
