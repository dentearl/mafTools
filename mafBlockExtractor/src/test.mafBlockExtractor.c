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
#include <inttypes.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include "CuTest.h"
#include "common.h"
#include "sharedMaf.h"
#include "mafBlockExtractorAPI.h"

static bool boolArraysAreEqual(bool *b1, bool *b2, uint32_t n) {
    for (uint32_t i = 0; i < n; ++i) {
        if (b1[i] != b2[i]) {
            return false;
        }
    }
    return true;
}
static void printBoolArray(bool *b, uint32_t n) {
    for (uint32_t i = 0; i < n; ++i) {
        if (b[i]) {
            printf("1");
        } else {
            printf("0");
        }
        if (i + 1 % 5 == 0) {
            printf(" ");
        }
    }
    printf("\n");
}
static bool mafLinesAreEqual(mafLine_t* ml1, mafLine_t *ml2) {
    if ((ml1 == NULL) || (ml2 == NULL)) {
        if ((ml1 == NULL) && (ml2 == NULL)) {
            return true;
        } else {
            fprintf(stderr, "mafLines differ: one mafLine is null, ml1:%p ml2:%p\n", (void*) ml1, (void*) ml2);
            return false;
        }
    }
    if (maf_mafLine_getLineNumber(ml1) != maf_mafLine_getLineNumber(ml2)) {
        fprintf(stderr, "mafLines differ in lineNumber:\n %3"PRIu32" %s\n %3"PRIu32" %s\n",
                maf_mafLine_getLineNumber(ml1),
                maf_mafLine_getLine(ml1), 
                maf_mafLine_getLineNumber(ml2),
                maf_mafLine_getLine(ml2));
        return false;
    }
    if (maf_mafLine_getType(ml1) != maf_mafLine_getType(ml2)) {
        fprintf(stderr, "mafLines differ in type\n");
        return false;
    }
    if (maf_mafLine_getStart(ml1) != maf_mafLine_getStart(ml2)) {
        fprintf(stderr, "mafLines differ in start: %3" PRIu32 " vs %3" PRIu32 "\n",
                maf_mafLine_getStart(ml1), maf_mafLine_getStart(ml2));
        return false;
    }
    if (maf_mafLine_getLength(ml1) != maf_mafLine_getLength(ml2)) {
        fprintf(stderr, "mafLines differ in length:\n  %3" PRIu32 ":%s\n  %3" PRIu32 ":%s\n",
                maf_mafLine_getLength(ml1),
                maf_mafLine_getSequence(ml1),
                maf_mafLine_getLength(ml2),
                maf_mafLine_getSequence(ml2));
        return false;
    }
    if (maf_mafLine_getStrand(ml1) != maf_mafLine_getStrand(ml2)) {
        fprintf(stderr, "mafLines differ in strand\n");
        return false;
    }
    if (maf_mafLine_getSourceLength(ml1) != maf_mafLine_getSourceLength(ml2)) {
        fprintf(stderr, "mafLines differ in source Length\n");
        return false;
    }
    if ((maf_mafLine_getSpecies(ml1) != NULL) && maf_mafLine_getSpecies(ml2) != NULL) {
        if (strcmp(maf_mafLine_getSpecies(ml1), maf_mafLine_getSpecies(ml2)) != 0) {
            fprintf(stderr, "mafLines differ in species name: ml1: [%s] ml2: [%s]\n",
                    maf_mafLine_getSpecies(ml1), maf_mafLine_getSpecies(ml2));
            return false;
        }
    } else {
        if (maf_mafLine_getSpecies(ml1) != maf_mafLine_getSpecies(ml2)) {
            // if true, one is NULL the other is not
            fprintf(stderr, "mafLines differ in species name, one is null: \n");
            if (maf_mafLine_getSpecies(ml1) == NULL) {
                fprintf(stderr, "  ml1: [NULL] ");
            } else {
                fprintf(stderr, "  ml1: [%s] ", maf_mafLine_getSpecies(ml1));
            }
            if (maf_mafLine_getSpecies(ml2) == NULL) {
                fprintf(stderr, "ml2: [NULL]\n");
            } else {
                fprintf(stderr, "ml2: [%s]\n", maf_mafLine_getSpecies(ml2));
            }
            return false;
        }
    }
    if ((maf_mafLine_getSequence(ml1) != NULL) && maf_mafLine_getSequence(ml2) != NULL) {
        if (strcmp(maf_mafLine_getSequence(ml1), maf_mafLine_getSequence(ml2)) != 0) {
            fprintf(stderr, "mafLines differ in sequence:\n  ml1: [%s]\n  ml2: [%s]\n",
                    maf_mafLine_getSequence(ml1), maf_mafLine_getSequence(ml2));
            return false;
        }
    } else {
        if (maf_mafLine_getSequence(ml1) != maf_mafLine_getSequence(ml2)) {
            // if true, one is NULL the other is not
            fprintf(stderr, "ml1 seq:%p ml2 seq:%p\n", 
                    (void*) maf_mafLine_getSequence(ml1), 
                    (void*) maf_mafLine_getSequence(ml2));
            return false;
        }
    }
    return true;
}
static bool mafBlocksAreEqual(mafBlock_t *mb1, mafBlock_t *mb2) {
    if ((mb1 == NULL) || (mb2 == NULL)) {
        if ((mb1 == NULL) && (mb2 == NULL)) {
            return true;
        } else {
            fprintf(stderr, "one mafBlock is null, mb1:%p mb2:%p\n", (void*) mb1, (void*) mb2);
            return false;
        }
    }
    if (maf_mafBlock_getLineNumber(mb1) != maf_mafBlock_getLineNumber(mb2)) {
        fprintf(stderr, "mafBlocks differ in lineNumber, %3" PRIu32 " vs %3" PRIu32 "\n",
                maf_mafBlock_getLineNumber(mb1), maf_mafBlock_getLineNumber(mb2));
        return false;
    }
    if (maf_mafBlock_getNumberOfLines(mb1) != maf_mafBlock_getNumberOfLines(mb2)) {
        fprintf(stderr, "mafBlocks differ in number of lines, %3"PRIu32" vs %3"PRIu32"\n",
                maf_mafBlock_getNumberOfLines(mb1), maf_mafBlock_getNumberOfLines(mb2));
        return false;
    }
    if (maf_mafBlock_getNumberOfSequences(mb1) != maf_mafBlock_getNumberOfSequences(mb2)) {
        fprintf(stderr, "mafBlocks differ in number of sequences\n");
        return false;
    }
    if (maf_mafBlock_getSequenceFieldLength(mb1) != maf_mafBlock_getSequenceFieldLength(mb2)) {
        fprintf(stderr, "mafBlocks differ in sequence field lengths, %3"PRIu32" vs %3"PRIu32"\n",
                maf_mafBlock_getSequenceFieldLength(mb1), maf_mafBlock_getSequenceFieldLength(mb2));
        fprintf(stderr, "mb1:\n");
        maf_mafBlock_print(mb1);
        fprintf(stderr, "mb2:\n");
        maf_mafBlock_print(mb2);
        return false;
    }
    mafLine_t *ml1, *ml2;
    ml1 = maf_mafBlock_getHeadLine(mb1);
    ml2 = maf_mafBlock_getHeadLine(mb2);
    while (ml1 != NULL) {
        if (! mafLinesAreEqual(ml1, ml2)) {
            fprintf(stderr, "mafLines differ\n");
            return false;
        }
        ml1 = maf_mafLine_getNext(ml1);
        ml2 = maf_mafLine_getNext(ml2);
    }
    return true;
}
static void targetColumnTest(CuTest *testCase, const char *mafString, uint32_t start, 
                             uint32_t stop, uint32_t expectedLen, bool expected[]) {
    mafBlock_t *ib = maf_newMafBlockFromString(mafString, 3);
    bool *targetColumns = NULL;
    uint32_t len = 0;
    getTargetColumns(&targetColumns, &len, ib, "theTarget.chr0", start, stop);
    CuAssertTrue(testCase, len == expectedLen);
    CuAssertTrue(testCase, boolArraysAreEqual(targetColumns, expected, len));
    maf_destroyMafBlockList(ib);
    free(targetColumns);
}
static void test_getTargetColumn_0(CuTest *testCase) {
    // test 0
    bool test0[] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    targetColumnTest(testCase, 
                     "a score=0\n"
                     "s theTarget.chr0 0 13 + 158545518 gcagctgaaaaca\n"
                     "s name.chr1      0 10 +       100 ATGT---ATGCCG\n"
                     "s name2.chr1     0 10 +       100 ATGT---ATGCCG\n",
                     0, 20, 13, test0);
    // test 1
    bool test1[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    targetColumnTest(testCase, 
                     "a score=0\n"
                     "s theTarget.chr0 0 13 + 158545518 gcagctgaaaaca\n"
                     "s name.chr1      0 10 +       100 ATGT---ATGCCG\n"
                     "s name2.chr1     0 10 +       100 ATGT---ATGCCG\n",
                     100, 200, 13, test1);
    // test 2
    bool test2[] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    targetColumnTest(testCase, 
                     "a score=0\n"
                     "s theTarget.chr0 0 13 - 158545518 gcagctgaaaaca\n"
                     "s name.chr1      0 10 +       100 ATGT---ATGCCG\n"
                     "s name2.chr1     0 10 +       100 ATGT---ATGCCG\n",
                     158545500, 158545518, 13, test2);
    // test 3
    bool test3[] = {1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1};
    targetColumnTest(testCase,
                     "a score=0\n"
                     "s theTarget.chr0 0 13 + 158545518 g--cagctgaaa--aca\n"
                     "s name.chr1      0 14 +       100 ATGT---ATGCCGAGGT\n"
                     "s name2.chr1     0 14 +       100 ATGT---ATGCCGAGGT\n",
                     0, 20, 17, test3);
    // test 4
    bool test4[] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1};
    targetColumnTest(testCase,
                     "a score=0\n"
                     "s theTarget.chr0 0 13 + 158545518 g--cagctgaaa--aca\n"
                     "s name.chr1      0 14 +       100 ATGT---ATGCCGAGGT\n"
                     "s name2.chr1     0 14 +       100 ATGT---ATGCCGAGGT\n"
                     "s theTarget.chr0 0 12 +       100 ATGT---ATGCC--GGT\n"
                     "s name3.chr1     0 14 +       100 ATGTAGCATGCCGAGGT\n",
                     0, 20, 17, test4);
    // test 5 
    bool test5[] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0};
    targetColumnTest(testCase, 
                     "a score=0\n"
                     "s theTarget.chr0 0 13 - 158545518 gcagctgaaaaca\n"
                     "s name.chr1      0 10 +       100 ATGT---ATGCCG\n"
                     "s name2.chr1     0 10 -       100 ATGT---ATGCCG\n",
                     158545507, 158545517, 13, test5);
}
static void trimTest(CuTest *testCase, const char *input, const char *expected, uint32_t n, bool isLeft) {
    mafBlock_t *ib = maf_newMafBlockFromString(input, 3);
    mafBlock_t *ob = maf_newMafBlockFromString(expected, 3);
    mafBlock_t *trimmed = trimBlock(ib, n, isLeft);
    CuAssertTrue(testCase, mafBlocksAreEqual(ob, trimmed));
    if (n == 0) {
        CuAssertTrue(testCase, ib == trimmed);
    }
    if (ib != trimmed)
        maf_destroyMafBlockList(trimmed);
    maf_destroyMafBlockList(ib);
    maf_destroyMafBlockList(ob);
}
static void test_leftTrim_0(CuTest *testCase) {
    // test 0
    trimTest(testCase, 
             "a score=0\n"
             "s theTarget.chr0 0 13 + 158545518 gcagctgaaaaca\n"
             "s name.chr1      0 10 +       100 ATGT---ATGCCG\n"
             "s name2.chr1     0 10 +       100 ATGT---ATGCCG\n",
             "a score=0\n"
             "s theTarget.chr0 2 11 + 158545518 agctgaaaaca\n"
             "s name.chr1      2  8 +       100 GT---ATGCCG\n"
             "s name2.chr1     2  8 +       100 GT---ATGCCG\n",
             2, true);
    // test 1
    trimTest(testCase,
             "a score=0\n"
             "s theTarget.chr0 0 13 + 158545518 gcagctgaaaaca\n"
             "s name.chr1      0 10 +       100 ATGT---ATGCCG\n"
             "s name2.chr1     0 10 +       100 ATGT---ATGCCG\n",
             "a score=0\n"
             "s theTarget.chr0 4 9 + 158545518 ctgaaaaca\n"
             "s name.chr1      4 6 +       100 ---ATGCCG\n"
             "s name2.chr1     4 6 +       100 ---ATGCCG\n",
             4, true);
    // test 2
    trimTest(testCase, 
             "a score=0\n"
             "s theTarget.chr0 0 13 + 158545518 gcagctgaaaaca\n"
             "s name.chr1      0  3 +       100 ATG----------\n"
             "s name2.chr1     0 10 +       100 ATGT---ATGCCG\n",
             "a score=0\n"
             "s theTarget.chr0 4 9 + 158545518 ctgaaaaca\n"
             "s name2.chr1     4 6 +       100 ---ATGCCG\n",
             4, true);
    // test 3
    trimTest(testCase,
             "a score=0\n"
             "s theTarget.chr0 0 13 + 158545518 gcagctgaaaaca\n"
             "s name.chr1      0 10 +       100 ATGG---ATGCCG\n"
             "s name2.chr1     0 10 +       100 ATGT---ATGCCG\n",
             "a score=0\n"
             "s theTarget.chr0 6 7 + 158545518 gaaaaca\n"
             "s name.chr1      4 6 +       100 -ATGCCG\n"
             "s name2.chr1     4 6 +       100 -ATGCCG\n",
             6, true);
    // test 4
    trimTest(testCase,
             "a score=0\n"
             "s theTarget.chr0 0 13 + 158545518 gcagctgaaaaca\n"
             "s name.chr1      0 10 +       100 ATGG---ATGCCG\n"
             "s name2.chr1     0 10 +       100 ATGT---ATGCCG\n",
             "a score=0\n"
             "s theTarget.chr0 12 1 + 158545518 a\n"
             "s name.chr1       9 1 +       100 G\n"
             "s name2.chr1      9 1 +       100 G\n",
             12, true);
    // test 5
    trimTest(testCase,
             "a score=0\n"
             "s theTarget.chr0 0 13 + 158545518 gcagctgaaaaca\n"
             "s name.chr1      0 10 +       100 ATGT---ATGCCG\n"
             "s name2.chr1     0 10 +       100 ATGT---ATGCCG\n",
             "a score=0\n"
             "s theTarget.chr0 0 13 + 158545518 gcagctgaaaaca\n"
             "s name.chr1      0 10 +       100 ATGT---ATGCCG\n"
             "s name2.chr1     0 10 +       100 ATGT---ATGCCG\n",
             0, true);
}
static void test_rightTrim_0(CuTest *testCase) {
    // test 0
    trimTest(testCase, 
             "a score=0\n"
             "s theTarget.chr0 0 13 + 158545518 gcagctgaaaaca\n"
             "s name.chr1      0 10 +       100 ATGT---ATGCCG\n"
             "s name2.chr1     0 10 +       100 ATGT---ATGCCG\n",
             "a score=0\n"
             "s theTarget.chr0 0 11 + 158545518 gcagctgaaaa\n"
             "s name.chr1      0  8 +       100 ATGT---ATGC\n"
             "s name2.chr1     0  8 +       100 ATGT---ATGC\n",
             2, false);
    // test 1
    trimTest(testCase,
             "a score=0\n"
             "s theTarget.chr0 0 13 + 158545518 gcagctgaaaaca\n"
             "s name.chr1      0 10 +       100 ATGT---ATGCCG\n"
             "s name2.chr1     0 10 +       100 ATGT---ATGCCG\n",
             "a score=0\n"
             "s theTarget.chr0 0  5 + 158545518 gcagc\n"
             "s name.chr1      0  4 +       100 ATGT-\n"
             "s name2.chr1     0  4 +       100 ATGT-\n",
             8, false);
    // test 2
    trimTest(testCase,
             "a score=0\n"
             "s theTarget.chr0 0 13 + 158545518 gcagctgaaaaca\n"
             "s name.chr1      0 10 +       100 ATGT---ATGCCG\n"
             "s name2.chr1     0 10 +       100 ATGT---ATGCCG\n",
             "a score=0\n"
             "s theTarget.chr0 0  1 + 158545518 g\n"
             "s name.chr1      0  1 +       100 A\n"
             "s name2.chr1     0  1 +       100 A\n",
             12, false);
}
static void splitTest(CuTest *testCase, const char *input, const char *expLeftStr, const char *expRightStr,
                      uint32_t n) {
    mafBlock_t *ib = maf_newMafBlockFromString(input, 3);
    mafBlock_t *expLeft = maf_newMafBlockFromString(expLeftStr, 3);
    mafBlock_t *expRight = NULL;
    if (expRightStr != NULL) {
        expRight = maf_newMafBlockFromString(expRightStr, 3);
    }
    mafBlock_t *obsLeft = NULL, *obsRight = NULL;
    splitBlock(ib, n, &obsLeft, &obsRight);
    CuAssertTrue(testCase, mafBlocksAreEqual(expLeft, obsLeft));
    CuAssertTrue(testCase, mafBlocksAreEqual(expRight, obsRight));
    maf_destroyMafBlockList(ib);
    maf_destroyMafBlockList(obsLeft);
    maf_destroyMafBlockList(obsRight);
    maf_destroyMafBlockList(expLeft);
    maf_destroyMafBlockList(expRight);
}
static void test_split_0(CuTest *testCase) {
    // test 0
    splitTest(testCase,
              "a score=0\n"
              "s theTarget.chr0 0 13 + 158545518 gcag--ctgaaaaca\n"
              "s name.chr1      0 15 +       100 ATGTATATTATGCCG\n"
              "s name2.chr1     0 15 +       100 ATGTATATAATGCCG\n",
              "a score=0\n"
              "s theTarget.chr0 0 4 + 158545518 gcag\n"
              "s name.chr1      0 4 +       100 ATGT\n"
              "s name2.chr1     0 4 +       100 ATGT\n",
              "a score=0\n"
              "s theTarget.chr0 4  9 + 158545518 --ctgaaaaca\n"
              "s name.chr1      4 11 +       100 ATATTATGCCG\n"
              "s name2.chr1     4 11 +       100 ATATAATGCCG\n",
              4);
    // test 1
    splitTest(testCase, 
              "a score=0\n"
              "s theTarget.chr0 0 13 + 158545518 gcagctgaaaaca\n"
              "s name.chr1      0 13 +       100 ATGTATTATGCCG\n"
              "s name2.chr1     0 13 +       100 ATGTATAATGCCG\n",
              "a score=0\n"
              "s theTarget.chr0 0 13 + 158545518 gcagctgaaaaca\n"
              "s name.chr1      0 13 +       100 ATGTATTATGCCG\n"
              "s name2.chr1     0 13 +       100 ATGTATAATGCCG\n",
              NULL,
              13);
}
static void processTrimTest(CuTest *testCase, const char *input, const char *expected, const char *seq,
                            uint32_t start, uint32_t stop, bool isLeft) {
    mafBlock_t *ib = maf_newMafBlockFromString(input, 3);
    mafBlock_t *ob = NULL;
    if (expected != NULL) {
        ob = maf_newMafBlockFromString(expected, 3);
    }
    mafBlock_t *trimmed = processBlockForTrim(ib, seq, start, stop, isLeft);
    CuAssertTrue(testCase, mafBlocksAreEqual(ob, trimmed));
    if (ib != trimmed)
        maf_destroyMafBlockList(trimmed);
    maf_destroyMafBlockList(ib);
    maf_destroyMafBlockList(ob);
}
                            
static void test_processLeftTrim_0(CuTest *testCase) {
    // test 0
    processTrimTest(testCase,
                    "a score=0\n"
                    "s theTarget.chr0 0 13 + 158545518 gcagctgaaaaca\n"
                    "s name.chr1      0 10 +       100 ATGT---ATGCCG\n"
                    "s name2.chr1     0 10 +       100 ATGT---ATGCCG\n",
                    "a score=0\n"
                    "s theTarget.chr0 2 11 + 158545518 agctgaaaaca\n"
                    "s name.chr1      2  8 +       100 GT---ATGCCG\n"
                    "s name2.chr1     2  8 +       100 GT---ATGCCG\n",
                    "theTarget.chr0", 2, 20, true);
    // test 1
    processTrimTest(testCase,
                    "a score=0\n"
                    "s theTarget.chr0 0 13 + 158545518 gcagctgaaaaca\n"
                    "s name.chr1      0 10 +       100 ATGT---ATGCCG\n"
                    "s name2.chr1     0 10 +       100 ATGT---ATGCCG\n",
                    "a score=0\n"
                    "s theTarget.chr0 4 9 + 158545518 ctgaaaaca\n"
                    "s name.chr1      4 6 +       100 ---ATGCCG\n"
                    "s name2.chr1     4 6 +       100 ---ATGCCG\n",
                    "theTarget.chr0", 4, 20, true);
    // test 2
    processTrimTest(testCase,
                    "a score=0\n"
                    "s theTarget.chr0 0 13 + 158545518 gcagctgaaaaca\n"
                    "s name.chr1      0  3 +       100 ATG----------\n"
                    "s name2.chr1     0 10 +       100 ATGT---ATGCCG\n",
                    "a score=0\n"
                    "s theTarget.chr0 4 9 + 158545518 ctgaaaaca\n"
                    "s name2.chr1     4 6 +       100 ---ATGCCG\n",
                    "theTarget.chr0", 4, 20, true);
    // test 3
    processTrimTest(testCase,
                    "a score=0\n"
                    "s theTarget.chr0 0 13 + 158545518 gcagctgaaaaca\n"
                    "s name.chr1      0 10 +       100 ATGG---ATGCCG\n"
                    "s name2.chr1     0 10 +       100 ATGT---ATGCCG\n",
                    "a score=0\n"
                    "s theTarget.chr0 6 7 + 158545518 gaaaaca\n"
                    "s name.chr1      4 6 +       100 -ATGCCG\n"
                    "s name2.chr1     4 6 +       100 -ATGCCG\n",
                    "theTarget.chr0", 6, 20, true);
    // test 4
    processTrimTest(testCase,
                    "a score=0\n"
                    "s theTarget.chr0 0 13 + 158545518 gcagctgaaaaca\n"
                    "s name.chr1      0 10 +       100 ATGG---ATGCCG\n"
                    "s name2.chr1     0 10 +       100 ATGT---ATGCCG\n",
                    "a score=0\n"
                    "s theTarget.chr0 12 1 + 158545518 a\n"
                    "s name.chr1       9 1 +       100 G\n"
                    "s name2.chr1      9 1 +       100 G\n",
                    "theTarget.chr0", 12, 20, true);
    // test 5
    processTrimTest(testCase,
                    "a score=0\n"
                    "s theTarget.chr0 0 13 + 158545518 gcagctgaaaaca\n"
                    "s name.chr1      0 10 +       100 ATGG---ATGCCG\n"
                    "s name2.chr1     0 10 +       100 ATGT---ATGCCG\n",
                    "a score=0\n"
                    "s theTarget.chr0 0 13 + 158545518 gcagctgaaaaca\n"
                    "s name.chr1      0 10 +       100 ATGG---ATGCCG\n"
                    "s name2.chr1     0 10 +       100 ATGT---ATGCCG\n",
                    "theTarget.chr0", 0, 20, true);
    // test 5
    processTrimTest(testCase,
                    "a score=0\n"
                    "s theTarget.chr0 0 13 + 158545518 gcagctgaaaaca\n"
                    "s name.chr1      0 10 +       100 ATGG---ATGCCG\n"
                    "s name2.chr1     0 10 +       100 ATGT---ATGCCG\n",
                    NULL,
                    "theTarget.chr0", 20, 30, true);
}
static void test_processLeftTrim_1(CuTest *testCase) {
    // test 0
    processTrimTest(testCase,
                    "a score=0\n"
                    "s theTarget.chr0 0 13 - 100 gcagctgaaaaca\n"
                    "s name.chr1      0 10 + 100 ATGT---ATGCCG\n"
                    "s name2.chr1     0 10 - 100 ATGT---ATGCCG\n",
                    "a score=0\n"
                    "s theTarget.chr0 2 11 - 100 agctgaaaaca\n"
                    "s name.chr1      2  8 + 100 GT---ATGCCG\n"
                    "s name2.chr1     2  8 - 100 GT---ATGCCG\n",
                    "theTarget.chr0", 80, 97, true);
    // test 1
    processTrimTest(testCase,
                    "a score=0\n"
                    "s theTarget.chr0 0 13 - 100 gcagctgaaaaca\n"
                    "s name.chr1      0 10 + 100 ATGT---ATGCCG\n"
                    "s name2.chr1     0 10 - 100 ATGT---ATGCCG\n",
                    "a score=0\n"
                    "s theTarget.chr0 4 9 - 100 ctgaaaaca\n"
                    "s name.chr1      4 6 + 100 ---ATGCCG\n"
                    "s name2.chr1     4 6 - 100 ---ATGCCG\n",
                    "theTarget.chr0", 80, 95, true);
}
static void test_processRightTrim_0(CuTest *testCase) {
    // test 0
    processTrimTest(testCase,
                    "a score=0\n"
                    "s theTarget.chr0 0 13 + 158545518 gcagctgaaaaca\n"
                    "s name.chr1      0 10 +       100 ATGT---ATGCCG\n"
                    "s name2.chr1     0 10 +       100 ATGT---ATGCCG\n",
                    "a score=0\n"
                    "s theTarget.chr0 0 11 + 158545518 gcagctgaaaa\n"
                    "s name.chr1      0  8 +       100 ATGT---ATGC\n"
                    "s name2.chr1     0  8 +       100 ATGT---ATGC\n",
                    "theTarget.chr0", 0, 10, false);
    // test 1
    processTrimTest(testCase,
                    "a score=0\n"
                    "s theTarget.chr0 0 13 + 158545518 gcagctgaaaaca\n"
                    "s name.chr1      0 10 +       100 ATGT---ATGCCG\n"
                    "s name2.chr1     0 10 +       100 ATGT---ATGCCG\n",
                    "a score=0\n"
                    "s theTarget.chr0 0  5 + 158545518 gcagc\n"
                    "s name.chr1      0  4 +       100 ATGT-\n"
                    "s name2.chr1     0  4 +       100 ATGT-\n",
                    "theTarget.chr0", 0, 4, false);
    // test 2
    processTrimTest(testCase,
                    "a score=0\n"
                    "s theTarget.chr0 0 13 + 158545518 gcagctgaaaaca\n"
                    "s name.chr1      0 10 +       100 ATGT---ATGCCG\n"
                    "s name2.chr1     0 10 +       100 ATGT---ATGCCG\n",
                    "a score=0\n"
                    "s theTarget.chr0 0  1 + 158545518 g\n"
                    "s name.chr1      0  1 +       100 A\n"
                    "s name2.chr1     0  1 +       100 A\n",
                    "theTarget.chr0", 0, 0, false);
    // test 3
    processTrimTest(testCase,
                    "a score=0\n"
                    "s theTarget.chr0 0 13 + 158545518 gcagctgaaaaca\n"
                    "s name.chr1      0 10 +       100 ATGG---ATGCCG\n"
                    "s name2.chr1     0 10 +       100 ATGT---ATGCCG\n",
                    "a score=0\n"
                    "s theTarget.chr0 0 13 + 158545518 gcagctgaaaaca\n"
                    "s name.chr1      0 10 +       100 ATGG---ATGCCG\n"
                    "s name2.chr1     0 10 +       100 ATGT---ATGCCG\n",
                    "theTarget.chr0", 0, 20, false);
}
static void test_processRightTrim_1(CuTest *testCase) {
    // test 0
    processTrimTest(testCase,
                    "a score=0\n"
                    "s theTarget.chr0 0 13 - 158545518 gcagctgaaaaca\n"
                    "s name.chr1      0 10 +       100 ATGT---ATGCCG\n"
                    "s name2.chr1     0 10 -       100 ATGT---ATGCCG\n",
                    "a score=0\n"
                    "s theTarget.chr0 0 11 - 158545518 gcagctgaaaa\n"
                    "s name.chr1      0  8 +       100 ATGT---ATGC\n"
                    "s name2.chr1     0  8 -       100 ATGT---ATGC\n",
                    "theTarget.chr0", 158545507, 158545517, false);
    // test 1
    processTrimTest(testCase,
                    "a score=0\n"
                    "s theTarget.chr0 0 13 - 100 gcagctgaaaaca\n"
                    "s name.chr1      0 10 + 100 ATGT---ATGCCG\n"
                    "s name2.chr1     0 10 - 100 ATGT---ATGCCG\n",
                    "a score=0\n"
                    "s theTarget.chr0 0  5 - 100 gcagc\n"
                    "s name.chr1      0  4 + 100 ATGT-\n"
                    "s name2.chr1     0  4 - 100 ATGT-\n",
                    "theTarget.chr0", 95, 99, false);
    // test 2
    processTrimTest(testCase,
                    "a score=0\n"
                    "s theTarget.chr0 0 13 - 158545518 gcagctgaaaaca\n"
                    "s name.chr1      0 10 +       100 ATGT---ATGCCG\n"
                    "s name2.chr1     0 10 -       100 ATGT---ATGCCG\n",
                    "a score=0\n"
                    "s theTarget.chr0 0  1 - 158545518 g\n"
                    "s name.chr1      0  1 +       100 A\n"
                    "s name2.chr1     0  1 -       100 A\n",
                    "theTarget.chr0", 158545517, 158545517, false);
    // test 3
    processTrimTest(testCase,
                    "a score=0\n"
                    "s theTarget.chr0 0 13 - 158545518 gcagctgaaaaca\n"
                    "s name.chr1      0 10 +       100 ATGG---ATGCCG\n"
                    "s name2.chr1     0 10 -       100 ATGT---ATGCCG\n",
                    "a score=0\n"
                    "s theTarget.chr0 0 13 - 158545518 gcagctgaaaaca\n"
                    "s name.chr1      0 10 +       100 ATGG---ATGCCG\n"
                    "s name2.chr1     0 10 -       100 ATGT---ATGCCG\n",
                    "theTarget.chr0", 158545500, 158545517, false);
}
static void processSplitTest(CuTest *testCase, const char *input, const char *expLeftStr, const char *expRightStr,
                             const char *seq, uint32_t start, uint32_t stop) {
    mafBlock_t *ib = maf_newMafBlockFromString(input, 3);
    mafBlock_t *expLeft = maf_newMafBlockFromString(expLeftStr, 3);
    mafBlock_t *expRight = NULL;
    if (expRightStr != NULL) {
        expRight = maf_newMafBlockFromString(expRightStr, 3);
    }
    mafBlock_t *obsLeft = NULL, *obsRight = NULL;
    processBlockForSplit(ib, seq, start, stop, &obsLeft, &obsRight);
    CuAssertTrue(testCase, mafBlocksAreEqual(expLeft, obsLeft));
    CuAssertTrue(testCase, mafBlocksAreEqual(expRight, obsRight));
    if (expRight == NULL) {
        CuAssertTrue(testCase, ib == obsLeft); // pointers should be the same
    }
    if (ib != obsLeft)
        maf_destroyMafBlockList(obsLeft);
    maf_destroyMafBlockList(ib);
    maf_destroyMafBlockList(obsRight);
    maf_destroyMafBlockList(expLeft);
    maf_destroyMafBlockList(expRight);
}
static void test_processSplit_0(CuTest *testCase) {
    // test 0
    processSplitTest(testCase, 
                     "a score=0\n"
                     "s theTarget.chr0 0 13 + 158545518 gcag--ctgaaaaca\n"
                     "s name.chr1      0 15 +       100 ATGTATATTATGCCG\n"
                     "s name2.chr1     0 15 +       100 ATGTATATAATGCCG\n",
                     "a score=0\n"
                     "s theTarget.chr0 0 4 + 158545518 gcag\n"
                     "s name.chr1      0 4 +       100 ATGT\n"
                     "s name2.chr1     0 4 +       100 ATGT\n",
                     "a score=0\n"
                     "s theTarget.chr0 4  9 + 158545518 --ctgaaaaca\n"
                     "s name.chr1      4 11 +       100 ATATTATGCCG\n"
                     "s name2.chr1     4 11 +       100 ATATAATGCCG\n",
                     "theTarget.chr0", 0, 20);
    // test 1
    processSplitTest(testCase, 
                     "a score=0\n"
                     "s theTarget.chr0 0 13 + 158545518 gcagctgaaaaca\n"
                     "s name.chr1      0 13 +       100 ATGTATTATGCCG\n"
                     "s name2.chr1     0 13 +       100 ATGTATAATGCCG\n",
                     "a score=0\n"
                     "s theTarget.chr0 0 13 + 158545518 gcagctgaaaaca\n"
                     "s name.chr1      0 13 +       100 ATGTATTATGCCG\n"
                     "s name2.chr1     0 13 +       100 ATGTATAATGCCG\n",
                     NULL,
                     "theTarget.chr0", 0, 20);
}
CuSuite* extractor_TestSuite(void) {
    CuSuite* suite = CuSuiteNew();
    (void) printBoolArray;
    SUITE_ADD_TEST(suite, test_getTargetColumn_0);
    SUITE_ADD_TEST(suite, test_leftTrim_0);
    SUITE_ADD_TEST(suite, test_rightTrim_0);
    SUITE_ADD_TEST(suite, test_split_0);
    SUITE_ADD_TEST(suite, test_processLeftTrim_0);
    SUITE_ADD_TEST(suite, test_processLeftTrim_1);
    SUITE_ADD_TEST(suite, test_processRightTrim_0);
    SUITE_ADD_TEST(suite, test_processRightTrim_1);
    SUITE_ADD_TEST(suite, test_processSplit_0);
    return suite;
}
