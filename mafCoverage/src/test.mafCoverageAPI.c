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
#include <limits.h>
#include "CuTest.h"
#include "common.h"
#include "sharedMaf.h"
#include "mafCoverageAPI.h"

static stHash *sequenceNamesToSequenceSizes = NULL;

static void setup() {
    sequenceNamesToSequenceSizes = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, free, (void(*)(void *)) stIntTuple_destruct);
    stHash_insert(sequenceNamesToSequenceSizes, stString_copy("bat.man"), stIntTuple_construct1(50));
    stHash_insert(sequenceNamesToSequenceSizes, stString_copy("spider.man"), stIntTuple_construct1(1));
    stHash_insert(sequenceNamesToSequenceSizes, stString_copy("bat.fink"), stIntTuple_construct1(7));
    stHash_insert(sequenceNamesToSequenceSizes, stString_copy("danger.mouse"), stIntTuple_construct1(12));
    stHash_insert(sequenceNamesToSequenceSizes, stString_copy("penfold"), stIntTuple_construct1(12));
}

static void teardown() {
    if (sequenceNamesToSequenceSizes != NULL) {
        stHash_destruct(sequenceNamesToSequenceSizes);
    }
}

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

static void test_getSpeciesNames(CuTest *testCase) {
    setup();
    stList *sequenceNames = stHash_getKeys(sequenceNamesToSequenceSizes);
    stSet *speciesNames = getSpeciesNames(sequenceNames, 0);
    CuAssertIntEquals(testCase, 4, stSet_size(speciesNames));
    CuAssertTrue(testCase, stSet_search(speciesNames, "bat") != NULL);
    CuAssertTrue(testCase, stSet_search(speciesNames, "spider") != NULL);
    CuAssertTrue(testCase, stSet_search(speciesNames, "danger") != NULL);
    CuAssertTrue(testCase, stSet_search(speciesNames, "penfold") != NULL);
    stSet_destruct(speciesNames);
    stList_destruct(sequenceNames);
    teardown();
}

static void test_getMapOfSequenceNamesToSequenceSizesForGivenSpecies(CuTest *testCase) {
    setup();
    stHash *sequenceNamesToSequenceSizeForGivenSpecies = getMapOfSequenceNamesToSequenceSizesForGivenSpeciesOrChr(sequenceNamesToSequenceSizes,
            "bat", 0);
    CuAssertIntEquals(testCase, 2, stHash_size(sequenceNamesToSequenceSizeForGivenSpecies));
    CuAssertIntEquals(testCase, 50, stIntTuple_get(stHash_search(sequenceNamesToSequenceSizeForGivenSpecies, "bat.man"), 0));
    CuAssertIntEquals(testCase, 7, stIntTuple_get(stHash_search(sequenceNamesToSequenceSizeForGivenSpecies, "bat.fink"), 0));
    stHash_destruct(sequenceNamesToSequenceSizeForGivenSpecies);
    teardown();
}

static void test_getTotalLengthOfSequences(CuTest *testCase) {
    setup();
    CuAssertIntEquals(testCase, 82, getTotalLengthOfSequences(sequenceNamesToSequenceSizes));
    teardown();
}

static void test_pairwiseCoverage(CuTest *testCase) {
    setup();
    PairwiseCoverage *pC = pairwiseCoverage_construct(sequenceNamesToSequenceSizes);
    //Check coverage is 0 when we start
    CuAssertDblEquals(testCase, 0.0, pairwiseCoverage_calculateCoverage(pC), 0.0);
    double *nCoverages = pairwiseCoverage_calculateNCoverages(pC);
    CuAssertDblEquals(testCase, 1.0, nCoverages[0], 0.0);
    for (int64_t i = 1; i <= SCHAR_MAX; i++) {
        CuAssertDblEquals(testCase, 0.0, nCoverages[i], 0.0);
    }
    free(nCoverages);

    //Add some coverage
    char *coverageArray = pairwiseCoverage_getCoverageArrayForSequence(pC, "spider.man");
    CuAssertTrue(testCase, coverageArray != NULL);
    pairwiseCoverageArray_increase(coverageArray, 0);
    coverageArray = pairwiseCoverage_getCoverageArrayForSequence(pC, "penfold");
    CuAssertTrue(testCase, coverageArray != NULL);
    pairwiseCoverageArray_increase(coverageArray, 2);
    pairwiseCoverageArray_increase(coverageArray, 2);

    //Now recalculate the coverages
    CuAssertDblEquals(testCase, 2.0/82.0, pairwiseCoverage_calculateCoverage(pC), 0.0);
    nCoverages = pairwiseCoverage_calculateNCoverages(pC);
    CuAssertDblEquals(testCase, 1.0, nCoverages[0], 0.0);
    CuAssertDblEquals(testCase, 2.0/82.0, nCoverages[1], 0.0);
    CuAssertDblEquals(testCase, 1.0/82.0, nCoverages[2], 0.0);
    free(nCoverages);
    nCoverages = pairwiseCoverage_calculateNCoverages(pC);
    for (int64_t i = 3; i <= SCHAR_MAX; i++) {
        CuAssertDblEquals(testCase, 0.0, nCoverages[i], 0.0);
    }
    free(nCoverages);

    pairwiseCoverage_destruct(pC);
    teardown();
}

static void test_nGenomeCoverage(CuTest *testCase) {
    setup();
    //Just build a single nGenomeCoverage and check the report functions work as expected.
    NGenomeCoverage *nGC = nGenomeCoverage_construct(sequenceNamesToSequenceSizes, "bat", 0);
    // nGenomeCoverage_reportHeader(stderr, 1);
    // nGenomeCoverage_report(nGC, stderr, 1);
    nGenomeCoverage_destruct(nGC);
    teardown();
}

CuSuite* coverage_TestSuite(void) {
    CuSuite* suite = CuSuiteNew();
    (void) test_is_wild_0;
    (void) test_searchMatched_0;
    SUITE_ADD_TEST(suite, test_is_wild_0);
    SUITE_ADD_TEST(suite, test_searchMatched_0);
    SUITE_ADD_TEST(suite, test_getSpeciesNames);
    SUITE_ADD_TEST(suite, test_getMapOfSequenceNamesToSequenceSizesForGivenSpecies);
    SUITE_ADD_TEST(suite, test_getTotalLengthOfSequences);
    SUITE_ADD_TEST(suite, test_pairwiseCoverage);
    SUITE_ADD_TEST(suite, test_nGenomeCoverage);
    return suite;
}
