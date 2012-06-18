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
#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include "CuTest.h"
#include "common.h"
#include "sonLib.h"
#include "ComparatorAPI.h"

static void printMat(uint64_t **mat, uint64_t n) {
    printf("printMat(mat, %" PRIi64 ")\n", n);
    for (uint64_t r = 0; r < n; ++r) {
        for (uint64_t c = 0; c < n; ++c) {
            if (r == c) {
                printf("  - ");
            } else if (c < r) {
                printf("    ");
            } else {
                printf(" %2" PRIi64 " ", mat[r][c]);
            }
        }
        printf("\n");
    }
    printf("\n");
}
static void test_mappingMatrixToArray_0(CuTest *testCase) {
    uint64_t **mat = NULL;
    uint64_t index, result;
    for (uint64_t n = 2; n < 11; ++n) {
        // build out mat
        mat = (uint64_t **) st_malloc(sizeof(*mat) * n);
        for (uint64_t r = 0; r < n; ++r) {
            mat[r] = st_malloc(sizeof(uint64_t *) * n);
        }
        index = 0;
        for (uint64_t r = 0; r < n; ++r) {
            for (uint64_t c = r + 1; c < n; ++c) {
                mat[r][c] = index++;
            }
        }
        // printMat(mat, n);
        // test
        for (uint64_t r = 0; r < n - 1; ++r) {
            for (uint64_t c = r + 1; c < n; ++c) {
                pairIndicesToArrayIndex(r, c, n, &result);
                CuAssertTrue(testCase, mat[r][c] == result);
            }
        }
        // clean up
        for (uint64_t r = 0; r < n; ++r) {
            free(mat[r]);
        }
        free(mat);
    }
}
static void test_mappingArrayToMatrix_0(CuTest *testCase) {
    uint64_t **mat = NULL;
    uint64_t index, resultRow, resultCol;
    for (uint64_t n = 2; n < 11; ++n) {
        // build out mat
        mat = (uint64_t **) st_malloc(sizeof(*mat) * n);
        for (uint64_t r = 0; r < n; ++r) {
            mat[r] = st_malloc(sizeof(uint64_t *) * n);
        }
        index = 0;
        for (uint64_t r = 0; r < n; ++r) {
            for (uint64_t c = r + 1; c < n; ++c) {
                mat[r][c] = index++;
            }
        }
        // printMat(mat, n);
        // test
        for (uint64_t r = 0; r < n - 1; ++r) {
            for (uint64_t c = r + 1; c < n; ++c) {
                arrayIndexToPairIndices(mat[r][c], n, &resultRow, &resultCol);
                CuAssertTrue(testCase, r == resultRow);
                CuAssertTrue(testCase, c == resultCol);
            }
        }
        // clean up
        for (uint64_t r = 0; r < n; ++r) {
            free(mat[r]);
        }
        free(mat);
    }
}
static void test_pairCounting_0(CuTest *testCase) {
    uint64_t *chooseTwoArray = buildChooseTwoArray();
    uint64_t expected, observed;
    uint32_t lineNumber, numSeqs, numLines;
    stHash *legitPairs = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, free, free);
    mafLine_t *ml = NULL;
    mafBlock_t *mb = NULL;
    numSeqs = 4;
    numLines = 5;
    expected = chooseTwo(4) * 13;
    char **input = (char**) st_malloc(sizeof(*input) * numLines);
    input[0] = stString_copy("a score=0.0");
    input[1] = stString_copy("s hg16.chr7    27707221 13 + 158545518 gcagctgaaaaca");
    input[2] = stString_copy("s panTro1.chr6 28869787 13 + 161576975 gcagctgaaaaca");
    input[3] = stString_copy("s baboon         249182 13 +   4622798 gcagctgaaaaca");
    input[4] = stString_copy("s mm4.chr6     53310102 13 + 151104725 ACAGCTGAAAATA");
    stHash_insert(legitPairs, stString_copy("hg16.chr7"), stString_copy(""));
    stHash_insert(legitPairs, stString_copy("panTro1.chr6"), stString_copy(""));
    stHash_insert(legitPairs, stString_copy("baboon"), stString_copy(""));
    stHash_insert(legitPairs, stString_copy("mm4.chr6"), stString_copy(""));
    mb = maf_newMafBlock();
    lineNumber = 3;
    maf_mafBlock_setLineNumber(mb, lineNumber);
    maf_mafBlock_setNumberOfSequences(mb, numSeqs);
    maf_mafBlock_setNumberOfLines(mb, numLines);
    maf_mafBlock_setHeadLine(mb, maf_newMafLineFromString(input[0], lineNumber));
    ml = maf_mafBlock_getHeadLine(mb);
    for (uint32_t i = 0; i < maf_mafBlock_getNumberOfLines(mb); ++i, ++lineNumber) {
        maf_mafLine_setNext(ml, maf_newMafLineFromString(input[i], lineNumber));
        ml = maf_mafLine_getNext(ml);
        maf_mafBlock_setTailLine(mb, ml);
    }
    observed = walkBlockCountingPairs(mb, legitPairs, chooseTwoArray);
    CuAssertTrue(testCase, observed == expected);
    // clean up
    for (unsigned i = 0; i < numLines; ++i) { 
        free(input[i]);
    }
    free(input);
    maf_destroyMafBlockList(mb);
    stHash_destruct(legitPairs);

    // new test
    expected = chooseTwo(3) * 13;
    legitPairs = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, free, free);
    numSeqs = 4;
    numLines = 5;
    input = (char**) st_malloc(sizeof(*input) * numLines);
    input[0] = stString_copy("a score=0.0");
    input[1] = stString_copy("s hg16.chr7    27707221 13 + 158545518 gcagctgaaaaca");
    input[2] = stString_copy("s panTro1.chr6 28869787 13 + 161576975 gcagctgaaaaca");
    input[3] = stString_copy("s baboon         249182 13 +   4622798 gcagctgaaaaca");
    input[4] = stString_copy("s mm4.chr6     53310102 13 + 151104725 ACAGCTGAAAATA");
    stHash_insert(legitPairs, stString_copy("hg16.chr7"), stString_copy(""));
    stHash_insert(legitPairs, stString_copy("panTro1.chr6"), stString_copy(""));
    stHash_insert(legitPairs, stString_copy("mm4.chr6"), stString_copy(""));
    mb = maf_newMafBlock();
    lineNumber = 3;
    maf_mafBlock_setLineNumber(mb, lineNumber);
    maf_mafBlock_setNumberOfSequences(mb, numSeqs);
    maf_mafBlock_setNumberOfLines(mb, numLines);
    maf_mafBlock_setHeadLine(mb, maf_newMafLineFromString(input[0], lineNumber));
    ml = maf_mafBlock_getHeadLine(mb);
    for (uint32_t i = 0; i < maf_mafBlock_getNumberOfLines(mb); ++i, ++lineNumber) {
        maf_mafLine_setNext(ml, maf_newMafLineFromString(input[i], lineNumber));
        ml = maf_mafLine_getNext(ml);
        maf_mafBlock_setTailLine(mb, ml);
    }
    observed = walkBlockCountingPairs(mb, legitPairs, chooseTwoArray);
    CuAssertTrue(testCase, observed == expected);
    // clean up
    for (unsigned i = 0; i < numLines; ++i) { 
        free(input[i]);
    }
    free(input);
    maf_destroyMafBlockList(mb);
    stHash_destruct(legitPairs);
    // new test
    expected = chooseTwo(3) + chooseTwo(4) * 37;
    legitPairs = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, free, free);
    numSeqs = 5;
    numLines = 6;
    input = (char**) st_malloc(sizeof(*input) * numLines);
    input[0] = stString_copy("a score=23262.0   ");
    input[1] = stString_copy("s hg18.chr7    27578828 38 + 158545518 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG");
    input[2] = stString_copy("s panTro1.chr6 28741140 38 + 161576975 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG");
    input[3] = stString_copy("s baboon         116834 38 +   4622798 AAA-GGGAATGTTAACCAAATGA---GTTGTCTCTTATGGTG");
    input[4] = stString_copy("s mm4.chr6     53215344 38 + 151104725 -AATGGGAATGTTAAGCAAACGA---ATTGTCTCTCAGTGTG");
    input[5] = stString_copy("s rn3.chr4     81344243 40 + 187371129 -AA-GGGGATGCTAAGCCAATGAGTTGTTGTCTCTCAATGTG");
    stHash_insert(legitPairs, stString_copy("hg18.chr7"), stString_copy(""));
    stHash_insert(legitPairs, stString_copy("panTro1.chr6"), stString_copy(""));
    stHash_insert(legitPairs, stString_copy("baboon"), stString_copy(""));
    stHash_insert(legitPairs, stString_copy("rn3.chr4"), stString_copy(""));
    mb = maf_newMafBlock();
    lineNumber = 3;
    maf_mafBlock_setLineNumber(mb, lineNumber);
    maf_mafBlock_setNumberOfSequences(mb, numSeqs);
    maf_mafBlock_setNumberOfLines(mb, numLines);
    maf_mafBlock_setHeadLine(mb, maf_newMafLineFromString(input[0], lineNumber));
    ml = maf_mafBlock_getHeadLine(mb);
    for (uint32_t i = 0; i < maf_mafBlock_getNumberOfLines(mb); ++i, ++lineNumber) {
        maf_mafLine_setNext(ml, maf_newMafLineFromString(input[i], lineNumber));
        ml = maf_mafLine_getNext(ml);
        maf_mafBlock_setTailLine(mb, ml);
    }
    observed = walkBlockCountingPairs(mb, legitPairs, chooseTwoArray);
    CuAssertTrue(testCase, observed == expected);
    // clean up
    for (unsigned i = 0; i < numLines; ++i) { 
        free(input[i]);
    }
    free(input);
    free(chooseTwoArray);
    maf_destroyMafBlockList(mb);
    stHash_destruct(legitPairs);
}
CuSuite* comparatorAPI_TestSuite(void) {
    (void) (printMat);
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, test_mappingMatrixToArray_0);
    SUITE_ADD_TEST(suite, test_mappingArrayToMatrix_0);
    SUITE_ADD_TEST(suite, test_pairCounting_0);
    return suite;
}
