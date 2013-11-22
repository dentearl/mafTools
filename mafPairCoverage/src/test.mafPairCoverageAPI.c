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
  CuAssertTrue(testCase, !is_wild("hg19.chr1*9"));
  CuAssertTrue(testCase, !is_wild("aoeuaoeunstaoeunshtonuts.chrcrhrc.huaoeunsatohunt."));
  CuAssertTrue(testCase, is_wild("aoeuaoeunstaoeunshtonuts.chrcrhrc.huaoeunsatohunt.*"));
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
  fprintf(stderr, "comparing lines\n");
  mafLine_t *ml1 = maf_newMafLineFromString("s hg19.chr19      123480 13 + 1234870098734 ACGTACGTACGTA", 1);
  mafLine_t *ml2 = maf_newMafLineFromString("s mm9.chr1        123480 13 + 1234870098735 ACGTACGTACGTA", 1);
  stHash *seq1Hash = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, free, free);
  stHash *seq2Hash = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, free, free);
  stHash *empty = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, free, free);
  uint64_t alignedPositions = 0;
  mafCoverageCount_t *mcct1 = createMafCoverageCount();
  mafCoverageCount_t *mcct2 = createMafCoverageCount();
  BinContainer *bc = binContainer_init();
  mafCoverageCount_setSourceLength(mcct1, 1234870098734);
  mafCoverageCount_setSourceLength(mcct2, 1234870098735);
  stHash_insert(seq1Hash, stString_copy("hg19.chr19"), mcct1);
  stHash_insert(seq2Hash, stString_copy("mm9.chr1"), mcct2);
  fprintf(stderr, "about to compare lines...\n");
  compareLines(ml1, ml2, seq1Hash, seq2Hash, &alignedPositions, empty, bc);
  CuAssertTrue(testCase, alignedPositions == 13);
  compareLines(ml1, ml2, seq1Hash, seq2Hash, &alignedPositions, empty, bc);
  CuAssertTrue(testCase, alignedPositions == 26);
  fprintf(stderr, "cleanning up.\n");
  maf_destroyMafLineList(ml1);
  maf_destroyMafLineList(ml2);
  stHash_destruct(seq1Hash);
  stHash_destruct(seq2Hash);
  stHash_destruct(empty);
  binContainer_destruct(bc);

  // test case 1
  alignedPositions = 0;
  ml1 = maf_newMafLineFromString("s hg19.chr19      123480 13 + 1234870098734 ACGTACGTACGTA", 1);
  ml2 = maf_newMafLineFromString("s mm9.chr2        123480  5 + 1234870098735 AC--------GTA", 1);
  seq1Hash = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, free, free);
  seq2Hash = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, free, free);
  empty = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, free, free);
  mcct1 = createMafCoverageCount();
  mcct2 = createMafCoverageCount();
  bc = binContainer_init();
  mafCoverageCount_setSourceLength(mcct1, 1234870098734);
  mafCoverageCount_setSourceLength(mcct2, 1234870098735);
  stHash_insert(seq1Hash, stString_copy("hg19.chr19"), mcct1);
  stHash_insert(seq2Hash, stString_copy("mm9.chr2"), mcct2);
  compareLines(ml1, ml2, seq1Hash, seq2Hash, &alignedPositions, empty, bc);
  CuAssertTrue(testCase, alignedPositions == 5);
  maf_destroyMafLineList(ml1);
  maf_destroyMafLineList(ml2);
  stHash_destruct(seq1Hash);
  stHash_destruct(seq2Hash);
  stHash_destruct(empty);
  binContainer_destruct(bc);
}
static void test_compareLines_1(CuTest *testCase) {
  // make sure that the hash is being correctly populated
  // test case 0
  mafLine_t *ml1 = maf_newMafLineFromString("s hg19.chr19      123480 13 + 1234870098734 ACGTACGTACGTA", 1);
  mafLine_t *ml2 = maf_newMafLineFromString("s mm9.chr1        123480 13 + 1234870098734 ACGTACGTACGTA", 1);
  stHash *seq1Hash = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, free, free);
  stHash *seq2Hash = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, free, free);
  stHash *empty = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, free, free);
  uint64_t alignedPositions = 0;
  mafCoverageCount_t *mcct1 = createMafCoverageCount();
  mafCoverageCount_t *mcct2 = createMafCoverageCount();
  BinContainer *bc = binContainer_init();
  mafCoverageCount_setSourceLength(mcct1, 1234870098734);
  mafCoverageCount_setSourceLength(mcct2, 1234870098734);
  stHash_insert(seq1Hash, stString_copy("hg19.chr19"), mcct1);
  stHash_insert(seq2Hash, stString_copy("mm9.chr1"), mcct2);
  compareLines(ml1, ml2, seq1Hash, seq2Hash, &alignedPositions, empty, bc);
  CuAssertTrue(testCase, stHash_search(seq1Hash, "hg19.chr19") != NULL);
  CuAssertTrue(testCase, stHash_search(seq2Hash, "mm9.chr1") != NULL);
  CuAssertTrue(testCase, stHash_search(seq2Hash, "bannana") == NULL);
  maf_destroyMafLineList(ml1);
  maf_destroyMafLineList(ml2);
  stHash_destruct(seq1Hash);
  stHash_destruct(seq2Hash);
  stHash_destruct(empty);
  binContainer_destruct(bc);
  // test case 1
}
static void test_intervalCheck_0(CuTest *testCase) {
  // make sure that the hash is being correctly populated
  // test case 0
  stHash *intervalsHash = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, free, free);
  stSortedSet *intervals = stSortedSet_construct3((int(*)(const void *, const void*)) stIntTuple_cmpFn,
                                                  (void(*)(void *)) stIntTuple_destruct);
  stHash_insert(intervalsHash, stString_copy("hg19.chr19"), intervals);
  stIntTuple *t = stIntTuple_construct2(123480, 123485);
  stSortedSet_insert(intervals, t);
  for (int i = 123480; i < 123485; ++i) {
    CuAssertTrue(testCase, inInterval(intervalsHash, "hg19.chr19", i) == true);
  }
  CuAssertTrue(testCase, inInterval(intervalsHash, "hg19.chr19", 123479) == false);
  CuAssertTrue(testCase, inInterval(intervalsHash, "hg19.chr19", 123486) == false);
  CuAssertTrue(testCase, inInterval(intervalsHash, "hg19.chr1", 123482) == false);
  CuAssertTrue(testCase, inInterval(intervalsHash, "hg19", 123482) == false);
  stSortedSet_destruct(intervals);
}
static void displayIntTuple(stIntTuple *t) {
  if (t != NULL) {
    printf("%" PRIi64 " %" PRIi64 "\n", stIntTuple_get(t, 0), stIntTuple_get(t, 1));
  } else {
    printf("NULL\n");
  }
}
static void displaySortedSet(stSortedSet *s) {
  stSortedSetIterator *sit = stSortedSet_getIterator(s);
  stIntTuple *t = NULL;
  while ((t = stSortedSet_getNext(sit)) != NULL) {
    displayIntTuple(t);
  }
}
static void test_compareLines_region_0(CuTest *testCase) {
  // make sure that the hash is being correctly populated
  // test case 0
  mafLine_t *ml1 = maf_newMafLineFromString("s hg19.chr19      "
                                            "123480 13 + 1234870098734 "
                                            "ACGTACGTACGTA", 1);
  mafLine_t *ml2 = maf_newMafLineFromString("s mm9.chr1        123480 13 + "
                                            "1234870098734 ACGTACGTACGTA", 1);
  stHash *seq1Hash = stHash_construct3(stHash_stringKey,
                                       stHash_stringEqualKey, free, free);
  stHash *seq2Hash = stHash_construct3(stHash_stringKey,
                                       stHash_stringEqualKey, free, free);
  stHash *intervalsHash = stHash_construct3(stHash_stringKey,
                                            stHash_stringEqualKey,
                                            free, free);
  stSortedSet *intervals = stSortedSet_construct3((int(*)(const void *, const void*)) stIntTuple_cmpFn,
                                                  (void(*)(void *)) stIntTuple_destruct);
  stHash_insert(intervalsHash, stString_copy("hg19.chr19"), intervals);
  stIntTuple *t = stIntTuple_construct2(123480, 123485);
  stIntTuple *q = stIntTuple_construct2(123482, 123484);
  stSortedSet_insert(intervals, t);
  uint64_t alignedPositions = 0;
  mafCoverageCount_t *mcct1 = createMafCoverageCount();
  mafCoverageCount_t *mcct2 = createMafCoverageCount();
  BinContainer *bc = binContainer_init();
  mafCoverageCount_setSourceLength(mcct1, 1234870098734);
  mafCoverageCount_setSourceLength(mcct2, 1234870098734);
  stHash_insert(seq1Hash, stString_copy("hg19.chr19"), mcct1);
  stHash_insert(seq2Hash, stString_copy("mm9.chr1"), mcct2);
  compareLines(ml1, ml2, seq1Hash, seq2Hash, &alignedPositions,
               intervalsHash, bc);
  CuAssertTrue(testCase, stHash_search(seq1Hash, "hg19.chr19") != NULL);
  CuAssertTrue(testCase, stHash_search(seq2Hash, "mm9.chr1") != NULL);
  CuAssertTrue(testCase, stHash_search(seq2Hash, "bannana") == NULL);
  stSortedSet *intSet = stHash_search(intervalsHash, "hg19.chr19");
  CuAssertTrue(testCase, intSet != NULL);
  stIntTuple *u = stSortedSet_search(intSet, q);
  CuAssertTrue(testCase, u == NULL);
  free(q);
  CuAssertTrue(testCase, mafCoverageCount_getInRegion(mcct1) == 4);
  CuAssertTrue(testCase, mafCoverageCount_getOutRegion(mcct1) == 9);
  CuAssertTrue(testCase, mafCoverageCount_getInRegion(mcct2) == 0);
  CuAssertTrue(testCase, mafCoverageCount_getOutRegion(mcct2) == 13);
  maf_destroyMafLineList(ml1);
  maf_destroyMafLineList(ml2);
  stHash_destruct(seq1Hash);
  stHash_destruct(seq2Hash);
  stSortedSet_destruct(intervals);
  binContainer_destruct(bc);
  // test case 1
}
static void test_binning_0(CuTest *testCase) {
  // make sure that the binning is working corretly
  // test case 0
  fprintf(stderr, "testing Binning\n");
  mafLine_t *ml1 = maf_newMafLineFromString("s hg19.chr19      "
                                            "123480 13 + 1234870098734 "
                                            "ACGTACGTACGTA", 1);
  mafLine_t *ml2 = maf_newMafLineFromString("s mm9.chr1        123480 13 + "
                                            "1234870098734 ACGTACGTACGTA", 1);
  stHash *seq1Hash = stHash_construct3(stHash_stringKey,
                                       stHash_stringEqualKey, free, free);
  stHash *seq2Hash = stHash_construct3(stHash_stringKey,
                                       stHash_stringEqualKey, free, free);
  stHash *intervalsHash = stHash_construct3(stHash_stringKey,
                                            stHash_stringEqualKey,
                                            free, free);
  stSortedSet *intervals = stSortedSet_construct3((int(*)(const void *, const void*)) stIntTuple_cmpFn,
                                                  (void(*)(void *)) stIntTuple_destruct);
  stHash_insert(intervalsHash, stString_copy("hg19.chr19"), intervals);
  stIntTuple *t = stIntTuple_construct2(123480, 123485);
  stSortedSet_insert(intervals, t);
  uint64_t alignedPositions = 0;
  mafCoverageCount_t *mcct1 = createMafCoverageCount();
  mafCoverageCount_t *mcct2 = createMafCoverageCount();
  BinContainer *bc = binContainer_construct(123480, 123500, 1);
  mafCoverageCount_setSourceLength(mcct1, 1234870098734);
  mafCoverageCount_setSourceLength(mcct2, 1234870098734);
  stHash_insert(seq1Hash, stString_copy("hg19.chr19"), mcct1);
  stHash_insert(seq2Hash, stString_copy("mm9.chr1"), mcct2);
  fprintf(stderr, "about to compare lines...\n");
  compareLines(ml1, ml2, seq1Hash, seq2Hash, &alignedPositions,
               intervalsHash, bc);
  fprintf(stderr, "lines compared.\n");
  CuAssertTrue(testCase, binContainer_getBinStart(bc) == 123480);
  CuAssertTrue(testCase, binContainer_getBinEnd(bc) == 123500);
  CuAssertTrue(testCase, binContainer_getBins(bc) != NULL);
  fprintf(stderr, "num bins: %" PRIi64 "\n", binContainer_getNumBins(bc));
  for (int i = 0; i < binContainer_getNumBins(bc); ++i) {
    if ((i != 0) && !(i % 5)) {
      fprintf(stderr, " ");
    }
    fprintf(stderr, "%" PRIu64 " ", binContainer_accessBin(bc, i));
  }
  fprintf(stderr, "\n");
  CuAssertTrue(testCase, binContainer_getNumBins(bc) == 21);
  for (int i = 0; i < 13; ++i) {
    CuAssertTrue(testCase, binContainer_accessBin(bc, i) == 1);
  }
  for (int i = 13; i < binContainer_getNumBins(bc); ++i) {
    CuAssertTrue(testCase, binContainer_accessBin(bc, i) == 0);
  }
  fprintf(stderr, "cleanning up.\n");
  maf_destroyMafLineList(ml1);
  maf_destroyMafLineList(ml2);
  stHash_destruct(seq1Hash);
  stHash_destruct(seq2Hash);
  stSortedSet_destruct(intervals);
  binContainer_destruct(bc);
  // test case 1
}
CuSuite* pairCoverage_TestSuite(void) {
  CuSuite* suite = CuSuiteNew();
  (void) displaySortedSet;
  (void) test_is_wild_0;
  (void) test_searchMatched_0;
  (void) test_compareLines_0;
  (void) test_compareLines_1;
  (void) test_compareLines_region_0;
  (void) test_intervalCheck_0;
  (void) test_binning_0;
  SUITE_ADD_TEST(suite, test_is_wild_0);
  SUITE_ADD_TEST(suite, test_searchMatched_0);
  SUITE_ADD_TEST(suite, test_compareLines_0);
  SUITE_ADD_TEST(suite, test_compareLines_1);
  SUITE_ADD_TEST(suite, test_intervalCheck_0);
  SUITE_ADD_TEST(suite, test_compareLines_region_0);
  SUITE_ADD_TEST(suite, test_binning_0);
  return suite;
}
