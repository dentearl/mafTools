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
#include "mafExtractorAPI.h"

static bool boolArraysAreEqual(bool *b1, bool *b2, uint64_t n) {
    for (uint64_t i = 0; i < n; ++i) {
        if (b1[i] != b2[i]) {
            return false;
        }
    }
    return true;
}
static void printBoolArray(bool *b, uint64_t n) {
    for (uint64_t i = 0; i < n; ++i) {
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
        fprintf(stderr, "mafLines differ in lineNumber:\n %3"PRIu64" %s\n %3"PRIu64" %s\n",
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
        fprintf(stderr, "mafLines differ in start:\n %3" PRIu64 ":%s\n %3" PRIu64 ":%s\n",
                maf_mafLine_getStart(ml1), 
                maf_mafLine_getLine(ml1), 
                maf_mafLine_getStart(ml2),
                maf_mafLine_getLine(ml2));
        return false;
    }
    if (maf_mafLine_getLength(ml1) != maf_mafLine_getLength(ml2)) {
        fprintf(stderr, "mafLines differ in length:\n  %3" PRIu64 ":%s\n  %3" PRIu64 ":%s\n",
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
        fprintf(stderr, "mafBlocks differ in lineNumber, %3" PRIu64 " vs %3" PRIu64 ".\n",
                maf_mafBlock_getLineNumber(mb1), maf_mafBlock_getLineNumber(mb2));
        printf("mb1:\n");
        maf_mafBlock_print(mb1);
        printf("mb2:\n");
        maf_mafBlock_print(mb2);
        return false;
    }
    if (maf_mafBlock_getNumberOfLines(mb1) != maf_mafBlock_getNumberOfLines(mb2)) {
        fprintf(stderr, "mafBlocks differ in number of lines, %3"PRIu64" vs %3"PRIu64"\n",
                maf_mafBlock_getNumberOfLines(mb1), maf_mafBlock_getNumberOfLines(mb2));
        return false;
    }
    if (maf_mafBlock_getNumberOfSequences(mb1) != maf_mafBlock_getNumberOfSequences(mb2)) {
        fprintf(stderr, "mafBlocks differ in number of sequences\n");
        return false;
    }
    if (maf_mafBlock_getSequenceFieldLength(mb1) != maf_mafBlock_getSequenceFieldLength(mb2)) {
        fprintf(stderr, "mafBlocks differ in sequence field lengths, %3"PRIu64" vs %3"PRIu64"\n",
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
static bool mafBlockListsAreEqual(mafBlock_t *head1, mafBlock_t *head2) {
    mafBlock_t *mb1 = head1, *mb2 = head2;
    if ((mb1 == NULL) || (mb2 == NULL)) {
        if ((mb1 == NULL) && (mb2 == NULL)) {
            return true;
        } else {
            fprintf(stderr, "one mafBlock list is null, mb1:%p mb2:%p\n", (void*) mb1, (void*) mb2);
            printf("block list1:\n");
            maf_mafBlock_printList(head1);
            printf("block list2:\n");
            maf_mafBlock_printList(head2);
            return false;
        }
    }
    uint64_t count1 = 0;
    while (mb1 != NULL) {
        mb1 = maf_mafBlock_getNext(mb1);
        ++count1;
    }
    uint64_t count2 = 0;
    while (mb2 != NULL) {
        mb2 = maf_mafBlock_getNext(mb2);
        ++count2;
    }
    if (count1 != count2) {
        fprintf(stderr, "mafBlock lists have different lengths!, mb1:%"PRIu64" mb2:%"PRIu64"\n", 
                count1, count2);
        printf("block list1:\n");
        maf_mafBlock_printList(head1);
        printf("block list2:\n");
        maf_mafBlock_printList(head2);
        return false;
    }
    mb1 = head1;
    mb2 = head2;
    while (mb1 != NULL) {
        if (!mafBlocksAreEqual(mb1, mb2)) {
            return false;
        }
        mb1 = maf_mafBlock_getNext(mb1);
        mb2 = maf_mafBlock_getNext(mb2);
    }
    return true;
}
static void targetColumnTest(CuTest *testCase, const char *mafString, uint64_t start, 
                             uint64_t stop, uint64_t expectedLen, bool expected[]) {
    mafBlock_t *ib = maf_newMafBlockFromString(mafString, 3);
    bool *targetColumns = NULL;
    uint64_t len = 0;
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
static void spliceTest(CuTest *testCase, const char *input, const char *expected, uint64_t l, 
                       uint64_t r, int64_t **offs) {
    mafBlock_t *ib = maf_newMafBlockFromString(input, 3);
    mafBlock_t *eb = maf_newMafBlockFromString(expected, 3);
    bool cleanOffs = false;
    if (offs == NULL) {
        cleanOffs = true;
        offs = createOffsets(maf_mafBlock_getNumberOfSequences(ib));
    }
    mafBlock_t *ob = spliceBlock(ib, l, r, offs);
    CuAssertTrue(testCase, mafBlocksAreEqual(eb, ob));
    if ((l == 0) && (r == maf_mafBlock_getSequenceFieldLength(ib) - 1)) {
        CuAssertTrue(testCase, ib == ob);
    }
    // clean up
    if (cleanOffs) {
        destroyOffsets(offs, maf_mafBlock_getNumberOfSequences(ib));
    }
    if (ib != ob)
        maf_destroyMafBlockList(ob);
    maf_destroyMafBlockList(ib);
    maf_destroyMafBlockList(eb);
}
static void test_splice_0(CuTest *testCase) {
    int64_t **offs = NULL;
    // int testcount = 0;
    // test 0
    // printf("test %d\n", testcount++);
    offs = createOffsets(3);
    /* for (uint64_t i = 0; i < 3; ++i) { */
    /*     printf("offs[%"PRIu64"][0] = %"PRIu64"\n", i, offs[i][0]); */
    /*     printf("offs[%"PRIu64"][1] = %"PRIu64"\n", i, offs[i][1]); */
    /* } */
    spliceTest(testCase, 
               "a score=0\n"
               "s theTarget.chr0 0 13 + 158545518 gcagctgaaaaca\n"
               "s name.chr1      0 10 +       100 ATGT---ATGCCG\n"
               "s name2.chr1     0 10 +       100 ATGT---ATGCCG\n",
               "a score=0\n"
               "s theTarget.chr0 2 11 + 158545518 agctgaaaaca\n"
               "s name.chr1      2  8 +       100 GT---ATGCCG\n"
               "s name2.chr1     2  8 +       100 GT---ATGCCG\n",
               2, 12, offs);
    /* printf("##########\n"); */
    /* for (uint64_t i = 0; i < 3; ++i) { */
    /*     printf("offs[%"PRIu64"][0] = %"PRIu64"\n", i, offs[i][0]); */
    /*     printf("offs[%"PRIu64"][1] = %"PRIu64"\n", i, offs[i][1]); */
    /* } */
    CuAssertTrue(testCase, offs[0][0] == 12); // seq field coord
    CuAssertTrue(testCase, offs[0][1] == 12); // non-gap offset
    for (uint64_t i = 1; i < 3; ++i) {
        CuAssertTrue(testCase, offs[i][0] == 12);
        CuAssertTrue(testCase, offs[i][1] == 9);
    }
    destroyOffsets(offs, 3);
    // test 1
    // printf("test %d\n", testcount++);
    spliceTest(testCase,
               "a score=0\n"
               "s theTarget.chr0 0 13 + 158545518 gcagctgaaaaca\n"
               "s name.chr1      0 10 +       100 ATGT---ATGCCG\n"
               "s name2.chr1     0 10 +       100 ATGT---ATGCCG\n",
               "a score=0\n"
               "s theTarget.chr0 4 9 + 158545518 ctgaaaaca\n"
               "s name.chr1      4 6 +       100 ---ATGCCG\n"
               "s name2.chr1     4 6 +       100 ---ATGCCG\n",
               4, 12, NULL);
    // test 2
    // printf("test %d\n", testcount++);
    spliceTest(testCase, 
               "a score=0\n"
               "s theTarget.chr0 0 13 + 158545518 gcagctgaaaaca\n"
               "s name.chr1      0  3 +       100 ATG----------\n"
               "s name2.chr1     0 10 +       100 ATGT---ATGCCG\n",
               "a score=0\n"
               "s theTarget.chr0 4 9 + 158545518 ctgaaaaca\n"
               "s name2.chr1     4 6 +       100 ---ATGCCG\n",
               4, 12, NULL);
    // test 3
    // printf("test %d\n", testcount++);
    spliceTest(testCase,
               "a score=0\n"
               "s theTarget.chr0 0 13 + 158545518 gcagctgaaaaca\n"
               "s name.chr1      0 10 +       100 ATGG---ATGCCG\n"
               "s name2.chr1     0 10 +       100 ATGT---ATGCCG\n",
               "a score=0\n"
               "s theTarget.chr0 6 7 + 158545518 gaaaaca\n"
               "s name.chr1      4 6 +       100 -ATGCCG\n"
               "s name2.chr1     4 6 +       100 -ATGCCG\n",
               6, 12, NULL);
    // test 4
    // printf("test %d\n", testcount++);
    spliceTest(testCase,
               "a score=0\n"
               "s theTarget.chr0 0 13 + 158545518 gcagctgaaaaca\n"
               "s name.chr1      0 10 +       100 ATGG---ATGCCG\n"
               "s name2.chr1     0 10 +       100 ATGT---ATGCCG\n",
               "a score=0\n"
               "s theTarget.chr0 12 1 + 158545518 a\n"
               "s name.chr1       9 1 +       100 G\n"
               "s name2.chr1      9 1 +       100 G\n",
               12, 12, NULL);
    // test 5
    // printf("test %d\n", testcount++);
    spliceTest(testCase,
               "a score=0\n"
               "s theTarget.chr0 0 13 + 158545518 gcagctgaaaaca\n"
               "s name.chr1      0 10 +       100 ATGT---ATGCCG\n"
               "s name2.chr1     0 10 +       100 ATGT---ATGCCG\n",
               "a score=0\n"
               "s theTarget.chr0 0 13 + 158545518 gcagctgaaaaca\n"
               "s name.chr1      0 10 +       100 ATGT---ATGCCG\n"
               "s name2.chr1     0 10 +       100 ATGT---ATGCCG\n",
               0, 12, NULL);
    // test 6
    // printf("test %d\n", testcount++);
    spliceTest(testCase,
               "a score=0\n"
               "s theTarget.chr0 0 13 + 158545518 gcag--ctgaaaaca\n"
               "s name.chr1      0 15 +       100 ATGTATATTATGCCG\n"
               "s name2.chr1     0 15 +       100 ATGTATATAATGCCG\n",
               "a score=0\n"
               "s theTarget.chr0 0 4 + 158545518 gcag\n"
               "s name.chr1      0 4 +       100 ATGT\n"
               "s name2.chr1     0 4 +       100 ATGT\n",
               0, 3, NULL);
    // test 7
    // printf("test %d\n", testcount++);
    offs = createOffsets(3);
    spliceTest(testCase,
               "a score=0\n"
               "s theTarget.chr0 0 13 + 158545518 gcag--ctgaaaaca\n"
               "s name.chr1      0 15 +       100 ATGTATATTATGCCG\n"
               "s name2.chr1     0 15 +       100 ATGTATATAATGCCG\n",
               "a score=0\n"
               "s theTarget.chr0 4  9 + 158545518 --ctgaaaaca\n"
               "s name.chr1      4 11 +       100 ATATTATGCCG\n"
               "s name2.chr1     4 11 +       100 ATATAATGCCG\n",
               4, 14, offs);
    CuAssertTrue(testCase, offs[0][0] == 14); // seq field coord
    CuAssertTrue(testCase, offs[0][1] == 12); // non-gap offset
    for (uint64_t i = 1; i < 3; ++i) {
        CuAssertTrue(testCase, offs[i][0] == 14);
        CuAssertTrue(testCase, offs[i][1] == 14);
    }
    destroyOffsets(offs, 3);
    // test 8
    // printf("test %d\n", testcount++);
    spliceTest(testCase, 
               "a score=0\n"
               "s theTarget.chr0 0 13 + 158545518 gcagctgaaaaca\n"
               "s name.chr1      0 13 +       100 ATGTATTATGCCG\n"
               "s name2.chr1     0 13 +       100 ATGTATAATGCCG\n",
               "a score=0\n"
               "s theTarget.chr0 5 6 + 158545518 tgaaaa\n"
               "s name.chr1      5 6 +       100 TTATGC\n"
               "s name2.chr1     5 6 +       100 TAATGC\n",
               5, 10, NULL);
    // test 9
    // printf("test %d\n", testcount++);
    offs = createOffsets(3);
    spliceTest(testCase, 
               "a score=0\n"
               "s theTarget.chr0 0 13 + 158545518 gcagctgaaaaca\n"
               "s name.chr1      0 13 +       100 ATGTATTATGCCG\n"
               "s name2.chr1     0 13 +       100 ATGTATAATGCCG\n",
               "a score=0\n"
               "s theTarget.chr0 0 13 + 158545518 gcagctgaaaaca\n"
               "s name.chr1      0 13 +       100 ATGTATTATGCCG\n"
               "s name2.chr1     0 13 +       100 ATGTATAATGCCG\n",
               0, 12, offs);
    for (uint64_t i = 0; i < 3; ++i) {
        // SPECIAL CASE! 
        // since the block is passed back directly, since no splice is performed,
        // there is no update made to the offset array. This makes sense because 
        // the array will not be used again.
        CuAssertTrue(testCase, offs[i][0] == 0);
        CuAssertTrue(testCase, offs[i][1] == -1);
    }
    destroyOffsets(offs, 3);

}
struct blockRecord {
    // ugly solution to task of keeping track of all blocks.
    mafBlock_t *b;
    struct blockRecord *next;
};
typedef struct blockRecord blockRecord_t;
static void processSpliceTest(CuTest *testCase, const char *seq, uint64_t start, uint64_t stop, 
                              const char *input, int numExpectedBlocks, ...) {
    // varargs that come in are a variable number of const char * items that should be built
    // into a mafBlock_t list.
    mafBlock_t *ibhead = maf_newMafBlockListFromString(input, 3);
    // printf("ibhead 0:\n");
    // maf_mafBlock_printList(ibhead);
    mafBlock_t *ib = ibhead, *eb = NULL, *ep = NULL, *obhead = NULL, *ob = NULL, *tmp = NULL;
    blockRecord_t *headRecord = NULL, *br = NULL;
    while (ib != NULL) {
        // since processBlockForSplice can actually alter the input it we need to make a record
        // of the input blocks so that they can all be free'd later
        if (headRecord == NULL) {
            headRecord = (blockRecord_t *) de_malloc(sizeof(blockRecord_t));
            headRecord->next = NULL;
            br = headRecord;
        } else {
            br->next = (blockRecord_t *) de_malloc(sizeof(blockRecord_t));
            br = br->next;
            br->next = NULL;
        }
        br->b = ib;
        ib = maf_mafBlock_getNext(ib);
    }
    ib = ibhead;
    va_list argp;
    const char *s;
    va_start(argp, numExpectedBlocks);
    for (int i = 0; i < numExpectedBlocks; ++i) {
        s = va_arg(argp, const char *);
        if (eb == NULL) {
            eb = maf_newMafBlockFromString(s, 3);
            ep = eb;
        } else {
            maf_mafBlock_setNext(ep, maf_newMafBlockFromString(s, 3));
            ep = maf_mafBlock_getNext(ep);
        }
    }
    va_end(argp);
    while (ib != NULL) {
        // process each member of the maf block linked list individiually
        tmp = processBlockForSplice(ib, 1, seq, start, stop, true);
        if (obhead == NULL) {
            ob = tmp;
            obhead = ob;
        } else {
            maf_mafBlock_setNext(ob, tmp);
            ob = maf_mafBlock_getNext(ob);
        }
        ib = maf_mafBlock_getNext(ib);
    }
    CuAssertTrue(testCase, mafBlockListsAreEqual(eb, obhead));
    // printf("ibhead 1:\n");
    // maf_mafBlock_printList(ibhead);
    // printf("obhead:\n");
    // maf_mafBlock_printList(obhead);
    if (ibhead != obhead) {
        // printf("  freeing obhead!\n");
        maf_destroyMafBlockList(obhead);
    } else {
        // printf("  leaving obhead!\n");
    }
    // maf_destroyMafBlockList(ibhead);
    while (headRecord != NULL) {
        // free all of the input maf blocks
        maf_destroyMafBlockList(headRecord->b);
        br = headRecord;
        headRecord = headRecord->next;
        free(br);
    }
    maf_destroyMafBlockList(eb);
}
static void test_processSplice_0(CuTest *testCase) {
    // int testcount = 0;
    // test 0
    // printf("test %d\n", testcount++);
    processSpliceTest(testCase, "theTarget.chr0", 2, 20,
                      "a score=0 test=0\n"
                      "s theTarget.chr0 0 13 + 158545518 gcagctgaaaaca\n"
                      "s name.chr1      0 10 +       100 ATGT---ATGCCG\n"
                      "s name2.chr1     0 10 +       100 ATGT---ATGCCG\n",
                      1, 
                      "a score=0\n"
                      "s theTarget.chr0 2 11 + 158545518 agctgaaaaca\n"
                      "s name.chr1      2  8 +       100 GT---ATGCCG\n"
                      "s name2.chr1     2  8 +       100 GT---ATGCCG\n");
    // test 1
    // printf("test %d\n", testcount++);
    processSpliceTest(testCase, "theTarget.chr0", 4, 20,
                      "a score=0 test=1\n"
                      "s theTarget.chr0 0 13 + 158545518 gcagctgaaaaca\n"
                      "s name.chr1      0 10 +       100 ATGT---ATGCCG\n"
                      "s name2.chr1     0 10 +       100 ATGT---ATGCCG\n",
                      1, 
                      "a score=0\n"
                      "s theTarget.chr0 4 9 + 158545518 ctgaaaaca\n"
                      "s name.chr1      4 6 +       100 ---ATGCCG\n"
                      "s name2.chr1     4 6 +       100 ---ATGCCG\n");
    // test 2
    // printf("test %d\n", testcount++);
    processSpliceTest(testCase, "theTarget.chr0", 4, 20,
                      "a score=0 test=2\n"
                      "s theTarget.chr0 0 13 + 158545518 gcagctgaaaaca\n"
                      "s name.chr1      0  3 +       100 ATG----------\n"
                      "s name2.chr1     0 10 +       100 ATGT---ATGCCG\n",
                      1, 
                      "a score=0\n"
                      "s theTarget.chr0 4 9 + 158545518 ctgaaaaca\n"
                      "s name2.chr1     4 6 +       100 ---ATGCCG\n");
    // test 3
    // printf("test %d\n", testcount++);
    processSpliceTest(testCase, "theTarget.chr0", 6, 20,
                      "a score=0 test=3\n"
                      "s theTarget.chr0 0 13 + 158545518 gcagctgaaaaca\n"
                      "s name.chr1      0 10 +       100 ATGG---ATGCCG\n"
                      "s name2.chr1     0 10 +       100 ATGT---ATGCCG\n",
                      1, 
                      "a score=0\n"
                      "s theTarget.chr0 6 7 + 158545518 gaaaaca\n"
                      "s name.chr1      4 6 +       100 -ATGCCG\n"
                      "s name2.chr1     4 6 +       100 -ATGCCG\n");
    // test 4
    // printf("test %d\n", testcount++);
    processSpliceTest(testCase, "theTarget.chr0", 12, 20,
                      "a score=0 test=4\n"
                      "s theTarget.chr0 0 13 + 158545518 gcagctgaaaaca\n"
                      "s name.chr1      0 10 +       100 ATGG---ATGCCG\n"
                      "s name2.chr1     0 10 +       100 ATGT---ATGCCG\n",
                      1, 
                      "a score=0\n"
                      "s theTarget.chr0 12 1 + 158545518 a\n"
                      "s name.chr1       9 1 +       100 G\n"
                      "s name2.chr1      9 1 +       100 G\n");
    // test 5
    // printf("test %d\n", testcount++);
    processSpliceTest(testCase, "theTarget.chr0", 0, 20,
                      "a score=0 test=5\n"
                      "s theTarget.chr0 0 13 + 158545518 gcagctgaaaaca\n"
                      "s name.chr1      0 10 +       100 ATGG---ATGCCG\n"
                      "s name2.chr1     0 10 +       100 ATGT---ATGCCG\n",
                      1, 
                      "a score=0\n"
                      "s theTarget.chr0 0 13 + 158545518 gcagctgaaaaca\n"
                      "s name.chr1      0 10 +       100 ATGG---ATGCCG\n"
                      "s name2.chr1     0 10 +       100 ATGT---ATGCCG\n");
    // test 6
    // printf("test %d\n", testcount++);
    processSpliceTest(testCase, "theTarget.chr0", 20, 30,
                      "a score=0 test=6\n"
                      "s theTarget.chr0 0 13 + 158545518 gcagctgaaaaca\n"
                      "s name.chr1      0 10 +       100 ATGG---ATGCCG\n"
                      "s name2.chr1     0 10 +       100 ATGT---ATGCCG\n",
                      0,
                      NULL);
    // test 7
    // printf("process test %d\n", testcount++);
    processSpliceTest(testCase, "theTarget.chr0", 0, 20,
                      "a score=0\n"
                      "s theTarget.chr0 0 13 + 158545518 --gcagctgaaaaca\n"
                      "s name.chr1      0 12 +       100 TTATGG---ATGCCG\n"
                      "s name2.chr1     0 12 +       100 TTATGT---ATGCCG\n",
                      1,
                      "a score=0\n"
                      "s theTarget.chr0 0 13 + 158545518 gcagctgaaaaca\n"
                      "s name.chr1      2 10 +       100 ATGG---ATGCCG\n"
                      "s name2.chr1     2 10 +       100 ATGT---ATGCCG\n"
                      );
    // test 8
    // printf("test %d\n", testcount++);
    processSpliceTest(testCase, "theTarget.chr0", 0, 20,
                      "a score=0\n"
                      "s theTarget.chr0 0 13 + 158545518 gcagctgaaaaca--\n"
                      "s name.chr1      0 12 +       100 ATGG---ATGCCGCC\n"
                      "s name2.chr1     0 12 +       100 ATGT---ATGCCGCC\n",
                      1,
                      "a score=0\n"
                      "s theTarget.chr0 0 13 + 158545518 gcagctgaaaaca\n"
                      "s name.chr1      0 10 +       100 ATGG---ATGCCG\n"
                      "s name2.chr1     0 10 +       100 ATGT---ATGCCG\n"
                      );
    // test 9
    // printf("test %d\n", testcount++);
    processSpliceTest(testCase, "theTarget.chr0", 0, 20,
                      "a score=0\n"
                      "s theTarget.chr0 0 13 + 158545518 gcagct---gaaaaca\n"
                      "s name.chr1      0 13 +       100 ATGG---TTTATGCCG\n"
                      "s name2.chr1     0 13 +       100 ATGT---TTTATGCCG\n",
                      2, 
                      "a score=0\n"
                      "s theTarget.chr0 0 6 + 158545518 gcagct\n"
                      "s name.chr1      0 4 +       100 ATGG--\n"
                      "s name2.chr1     0 4 +       100 ATGT--\n",
                      "a score=0\n"
                      "s theTarget.chr0 6 7 + 158545518 gaaaaca\n"
                      "s name.chr1      6 7 +       100 TATGCCG\n"
                      "s name2.chr1     6 7 +       100 TATGCCG\n"
                      );
    // test 10
    // printf("process test %d\n", testcount++);
    processSpliceTest(testCase, "theTarget.chr0", 2, 20,
                      "a score=0\n"
                      "s theTarget.chr0 0 13 + 158545518 gcagct---gaaa---aca\n"
                      "s name.chr1      0 16 +       100 ATGG---TTTATGCCGGGG\n"
                      "s name2.chr1     0 16 +       100 ATGT---TTTATGCCGGGG\n",
                      3, 
                      "a score=0\n"
                      "s theTarget.chr0 2 4 + 158545518 agct\n"
                      "s name.chr1      2 2 +       100 GG--\n"
                      "s name2.chr1     2 2 +       100 GT--\n",
                      "a score=0\n"
                      "s theTarget.chr0 6 4 + 158545518 gaaa\n"
                      "s name.chr1      6 4 +       100 TATG\n"
                      "s name2.chr1     6 4 +       100 TATG\n",
                      "a score=0\n"
                      "s theTarget.chr0 10 3 + 158545518 aca\n"
                      "s name.chr1      13 3 +       100 GGG\n"
                      "s name2.chr1     13 3 +       100 GGG\n"
                      );
    // test 11
    // printf("test %d\n", testcount++);
    processSpliceTest(testCase, "theTarget.chr0", 0, 20,
                      "a score=0\n"
                      "s theTarget.chr0 0  2 + 158545518 g--------------a\n"
                      "s name.chr1      0 13 +       100 ATGG---TTTATGCCG\n"
                      "s name2.chr1     0 13 +       100 ATGT---TTTATGCCG\n",
                      2, 
                      "a score=0\n"
                      "s theTarget.chr0 0 1 + 158545518 g\n"
                      "s name.chr1      0 1 +       100 A\n"
                      "s name2.chr1     0 1 +       100 A\n",
                      "a score=0\n"
                      "s theTarget.chr0  1 1 + 158545518 a\n"
                      "s name.chr1      12 1 +       100 G\n"
                      "s name2.chr1     12 1 +       100 G\n"
                      );
    // test 12
    // printf("test %d\n", testcount++);
    processSpliceTest(testCase, "theTarget.chr0", 8, 16,
                      "a score=0\n"
                      "s theTarget.chr0 0  0 -  20 gcagctgaaaaca\n"
                      "s name.chr1      0 10 + 100 ATTGT---AAGTG\n"
                      "s name2.chr1     0 10 + 100 ATTGT---AAGTG\n"
                      "s name3.chr1     0 10 + 100 ATTGT----AGTG\n",
                      1,
                      "a score=0\n"
                      "s theTarget.chr0 3 9 -  20 gctgaaaac\n"
                      "s name.chr1      3 6 + 100 GT---AAGT\n"
                      "s name2.chr1     3 6 + 100 GT---AAGT\n"
                      "s name3.chr1     3 5 + 100 GT----AGT\n"
                      );
    // test 13
    // printf("test %d\n", testcount++);
    processSpliceTest(testCase, "theTarget.chr0", 0, 20,
                      "a score=0\n"
                      "s theTarget.chr0 0 4 +  20 g-c-a-g\n"
                      "s name.chr1      0 6 + 100 ATTGTGG\n"
                      "s name2.chr1     0 5 + 100 ATTGT--\n"
                      "s name3.chr1     0 5 + 100 ATTGT--\n",
                      4,
                      "a score=0\n"
                      "s theTarget.chr0 0 1 +  20 g\n"
                      "s name.chr1      0 1 + 100 A\n"
                      "s name2.chr1     0 1 + 100 A\n"
                      "s name3.chr1     0 1 + 100 A\n",
                      "a score=0\n"
                      "s theTarget.chr0 1 1 +  20 c\n"
                      "s name.chr1      2 1 + 100 T\n"
                      "s name2.chr1     2 1 + 100 T\n"
                      "s name3.chr1     2 1 + 100 T\n",
                      "a score=0\n"
                      "s theTarget.chr0 2 1 +  20 a\n"
                      "s name.chr1      4 1 + 100 T\n"
                      "s name2.chr1     4 1 + 100 T\n"
                      "s name3.chr1     4 1 + 100 T\n",
                      "a score=0\n"
                      "s theTarget.chr0 3 1 +  20 g\n"
                      "s name.chr1      6 1 + 100 G\n"
                      );
    // test 14
    // printf("test %d\n", testcount++);
    processSpliceTest(testCase, "theTarget.chr0", 0, 20,
                      "a score=0\n"
                      "s theTarget.chr0 16 4 -  20 gc-ag\n"
                      "s name.chr1       0 5 + 100 ATTGG\n"
                      "s name2.chr1      0 5 + 100 ATTGG\n"
                      "s name3.chr1      0 5 + 100 ATTGG\n",
                      2,
                      "a score=0\n"
                      "s theTarget.chr0 16 2 -  20 gc\n"
                      "s name.chr1       0 2 + 100 AT\n"
                      "s name2.chr1      0 2 + 100 AT\n"
                      "s name3.chr1      0 2 + 100 AT\n",
                      "a score=0\n"
                      "s theTarget.chr0 18 2 -  20 ag\n"
                      "s name.chr1       3 2 + 100 GG\n"
                      "s name2.chr1      3 2 + 100 GG\n"
                      "s name3.chr1      3 2 + 100 GG\n"
                      );
    // test 15
    // printf("test %d\n", testcount++);
    processSpliceTest(testCase, "theTarget.chr0", 0, 20,
                      "a score=0\n"
                      "s theTarget.chr0 16 4 -  20 gc-ag\n"
                      "s name.chr1       0 5 + 100 ATTGG\n"
                      "s name2.chr1      0 5 + 100 ATTGG\n"
                      "s name3.chr1      0 5 + 100 ATTGG\n",
                      2,
                      "a score=0\n"
                      "s theTarget.chr0 16 2 -  20 gc\n"
                      "s name.chr1       0 2 + 100 AT\n"
                      "s name2.chr1      0 2 + 100 AT\n"
                      "s name3.chr1      0 2 + 100 AT\n",
                      "a score=0\n"
                      "s theTarget.chr0 18 2 -  20 ag\n"
                      "s name.chr1       3 2 + 100 GG\n"
                      "s name2.chr1      3 2 + 100 GG\n"
                      "s name3.chr1      3 2 + 100 GG\n"
                      );
    // test 15
    // printf("test %d\n", testcount++);
    processSpliceTest(testCase, "theTarget.chr0", 0, 17,
                      "a score=0\n"
                      "s theTarget.chr0  0 4 -  20 gc-ag\n"
                      "s name.chr1       0 5 + 100 ATTGG\n"
                      "s name2.chr1      0 5 + 100 ATTGG\n"
                      "s name3.chr1     99 1 - 100 ---G-\n",
                      1,
                      "a score=0\n"
                      "s theTarget.chr0  2 2 -  20 ag\n"
                      "s name.chr1       3 2 + 100 GG\n"
                      "s name2.chr1      3 2 + 100 GG\n"
                      "s name3.chr1     99 1 - 100 G-\n"
                      );
    // test 15
    // printf("test %d\n", testcount++);
    processSpliceTest(testCase, "theTarget.chr0", 0, 20,
                      "a score=0\n"
                      "s theTarget.chr0  0 2 +  20 ag---\n"
                      "s theTarget.chr0  0 2 -  20 ---ag\n"
                      "s name.chr1       0 5 + 100 ATTGG\n"
                      "s name2.chr1      0 5 + 100 ATTGG\n"
                      "s name3.chr1     99 1 - 100 ---G-\n",
                      2,
                      "a score=0\n"
                      "s theTarget.chr0  0 2 +  20 ag\n"
                      "s name.chr1       0 2 + 100 AT\n"
                      "s name2.chr1      0 2 + 100 AT\n",
                      "a score=0\n"
                      "s theTarget.chr0  0 2 -  20 ag\n"
                      "s name.chr1       3 2 + 100 GG\n"
                      "s name2.chr1      3 2 + 100 GG\n"
                      "s name3.chr1     99 1 - 100 G-\n"
                      );
    // test 16
    // printf("test %d\n", testcount++);
    processSpliceTest(testCase, "theTarget.chr0", 0, 20,
                      "a score=0\n"
                      "s theTarget.chr0  0 3 +  20 agg--\n"
                      "s theTarget.chr0  0 2 -  20 ---ag\n"
                      "s name.chr1       0 5 + 100 ATTGG\n"
                      "s name2.chr1      0 5 + 100 ATTGG\n"
                      "s name3.chr1     99 1 - 100 ---G-\n",
                      1,
                      "a score=0\n"
                      "s theTarget.chr0  0 3 +  20 agg--\n"
                      "s theTarget.chr0  0 2 -  20 ---ag\n"
                      "s name.chr1       0 5 + 100 ATTGG\n"
                      "s name2.chr1      0 5 + 100 ATTGG\n"
                      "s name3.chr1     99 1 - 100 ---G-\n"
                      );
    // test 16
    // printf("test %d\n", testcount++);
    processSpliceTest(testCase, "theTarget.chr0", 0, 20,
                      "a score=0\n"
                      "s theTarget.chr0  0 3 +  20 ag---\n"
                      "s theTarget.chr0  0 2 -  20 ---ag\n"
                      "s name.chr1       0 5 + 100 ATTGG\n"
                      "s name2.chr1      0 5 + 100 ATTGG\n"
                      "s name3.chr1     99 1 - 100 ---G-\n",
                      2,
                      "a score=0\n"
                      "s theTarget.chr0  0 2 +  20 ag\n"
                      "s name.chr1       0 2 + 100 AT\n"
                      "s name2.chr1      0 2 + 100 AT\n",
                      "a score=0\n"
                      "s theTarget.chr0  0 2 -  20 ag\n"
                      "s name.chr1       3 2 + 100 GG\n"
                      "s name2.chr1      3 2 + 100 GG\n"
                      "s name3.chr1     99 1 - 100 G-\n"
                      );
    // test 17
    // printf("test %d\n", testcount++);
    processSpliceTest(testCase, "dm3.chr3R", 12450223, 12950222,
                      "a score=0\n"
                      "s dm3.chr3R                  15282289    2 -  27905053 TG-----------------------------------------\n"
                      "s dm3.chr3R                  15282292   41 -  27905053 --CCAAACGAAATGAAATTTTCAGTTGAGTTGCCACAAGCCCT\n"
                      "s droRho.ctg7180000755597        5772    2 -      5776 TG-----------------------------------------\n"
                      "s droRho.ctg7180000755597        5775    1 -      5776 --C----------------------------------------\n",
                      1,
                      "a score=0\n"
                      "s dm3.chr3R                  15282289    2 -  27905053 TG-----------------------------------------\n"
                      "s dm3.chr3R                  15282292   41 -  27905053 --CCAAACGAAATGAAATTTTCAGTTGAGTTGCCACAAGCCCT\n"
                      "s droRho.ctg7180000755597        5772    2 -      5776 TG-----------------------------------------\n"
                      "s droRho.ctg7180000755597        5775    1 -      5776 --C----------------------------------------\n"
                      );
    // test 18
    // printf("test %d\n", testcount++);
    processSpliceTest(testCase, "dm3.chr3R", 12450223, 12950222,
                      "a score=0\n"
                      "s dm3.chr3R                  15282289    2 -  27905053 -TG-----------------------------------------\n"
                      "s dm3.chr3R                  15282292   41 -  27905053 ---CCAAACGAAATGAAATTTTCAGTTGAGTTGCCACAAGCCCT\n"
                      "s droRho.ctg7180000755597        5772    2 -      5776 -TG-----------------------------------------\n"
                      "s droRho.ctg7180000755597        5775    1 -      5776 ---C----------------------------------------\n",
                      1,
                      "a score=0\n"
                      "s dm3.chr3R                  15282289    2 -  27905053 TG-----------------------------------------\n"
                      "s dm3.chr3R                  15282292   41 -  27905053 --CCAAACGAAATGAAATTTTCAGTTGAGTTGCCACAAGCCCT\n"
                      "s droRho.ctg7180000755597        5772    2 -      5776 TG-----------------------------------------\n"
                      "s droRho.ctg7180000755597        5775    1 -      5776 --C----------------------------------------\n"
                      );
    // test 19
    // printf("test %d\n", testcount++);
    processSpliceTest(testCase, "theTarget.chr0", 2, 20,
                      "a score=0\n"
                      "s theTarget.chr0 0 13 + 158545518 -gcagctgaaaaca\n"
                      "s name.chr1      0 10 +       100 -ATGT---ATGCCG\n"
                      "s name2.chr1     0 10 +       100 -ATGT---ATGCCG\n",
                      1, 
                      "a score=0\n"
                      "s theTarget.chr0 2 11 + 158545518 agctgaaaaca\n"
                      "s name.chr1      2  8 +       100 GT---ATGCCG\n"
                      "s name2.chr1     2  8 +       100 GT---ATGCCG\n");
    // test 20
    // printf("test %d\n", testcount++);
    processSpliceTest(testCase, "theTarget.chr0", 2, 20,
                      "a score=0\n"
                      "s theTarget.chr0 0 13 + 158545518 -gcagctgaa--aaca\n"
                      "s name.chr1      0 11 +       100 -ATGT---ATG-TCCG\n"
                      "s name2.chr1     0 12 +       100 -ATGT---ATGTTCCA\n",
                      2, 
                      "a score=0\n"
                      "s theTarget.chr0 2 7 + 158545518 agctgaa\n"
                      "s name.chr1      2 4 +       100 GT---AT\n"
                      "s name2.chr1     2 4 +       100 GT---AT\n",
                      "a score=0\n"
                      "s theTarget.chr0 9 4 + 158545518 aaca\n"
                      "s name.chr1      7 4 +       100 TCCG\n"
                      "s name2.chr1     8 4 +       100 TCCA\n"
                      );
    // test 21
    // printf("test %d\n", testcount++);
    processSpliceTest(testCase, "simHuman.chrJ", 81139727, 81639726,
                      "a score=0\n"
                      "s simCow.chrB      60357748   16 +  86443571 -------------G------------TGGGGACAAGGTTTA--\n"
                      "s simHuman.chrJ     6786835   15 -  88398963 -GAGTAATGTTCAGTG---------------------------\n"
                      "s simHuman.chrJ     6786872   21 -  88398963 ----------------------TGTACAGCAGCCCTGCTTAGT\n",
                      2, 
                      "a score=0\n"
                      "s simCow.chrB      60357748    1 +  86443571 ------------G--\n"
                      "s simHuman.chrJ     6786835   15 -  88398963 GAGTAATGTTCAGTG\n",
                      "a score=0\n"
                      "s simCow.chrB      60357749   15 +  86443571 ----TGGGGACAAGGTTTA--\n"
                      "s simHuman.chrJ     6786872   21 -  88398963 TGTACAGCAGCCCTGCTTAGT\n"
                      );
    // test 21
    // printf("test %d\n", testcount++);
    processSpliceTest(testCase, "simHuman.chrJ", 69805407, 70305406,
                      "a score=0.0000000\n"
                      "s simRat.chrR      72133413   33 +  88137694 G--------TTGTTTTTATTGATGCTAGTAGTTTGACAACT\n"
                      "s simHuman.chrJ    18234935   38 -  88398963 GTT----GATAGTTGAGAATACACCAAGCTTGT-GCCTGTT\n",
                      3,
                      "a score=0\n"
                      "s simRat.chrR      72133413   1 +  88137694 G--\n"
                      "s simHuman.chrJ    18234935   3 -  88398963 GTT\n",
                      "a score=0\n"
                      "s simRat.chrR      72133414   24 +  88137694 --TTGTTTTTATTGATGCTAGTAGTT\n"
                      "s simHuman.chrJ    18234938   26 -  88398963 GATAGTTGAGAATACACCAAGCTTGT\n",
                      "a score=0\n"
                      "s simRat.chrR      72133439   7 +  88137694 GACAACT\n"
                      "s simHuman.chrJ    18234964   7 -  88398963 GCCTGTT\n"
                      );
    // test 22
    // printf("test %d\n", testcount++);
    processSpliceTest(testCase, "simHuman.chrJ", 69805407, 70305406,
                      "a score=375.0\n"
                      "s simHuman.chrJ    69889402  29 +  88398963 CAATTTAACTTGGATGATTGGTTGTCATA\n"
                      "s simDog.chrF      43936434  16 +  64906724 CAATTTAACTAGGGTG-------------\n\n"
                      "a score=6637.0\n"
                      "s simHuman.chrJ    69889669   19 +  88398963 -------------TCATATATCATGCCCT-TAA\n"
                      "s simDog.chrF      20970287   20 -  64906724 T-------------TGTATATCATGCCCTGTTC\n",
                      3,
                      "a score=0\n"
                      "s simHuman.chrJ    69889402  29 +  88398963 CAATTTAACTTGGATGATTGGTTGTCATA\n"
                      "s simDog.chrF      43936434  16 +  64906724 CAATTTAACTAGGGTG-------------\n",
                      "a score=0\n"
                      "s simHuman.chrJ    69889669   16 +  88398963 TCATATATCATGCCCT\n"
                      "s simDog.chrF      20970288   15 -  64906724 -TGTATATCATGCCCT\n",
                      "a score=0\n"
                      "s simHuman.chrJ    69889685    3 +  88398963 TAA\n"
                      "s simDog.chrF      20970304    3 -  64906724 TTC\n"
                      );

}
CuSuite* extractor_TestSuite(void) {
    CuSuite* suite = CuSuiteNew();
    (void) printBoolArray;
    (void) test_getTargetColumn_0;
    (void) test_splice_0;
    (void) test_processSplice_0;
    SUITE_ADD_TEST(suite, test_getTargetColumn_0);
    SUITE_ADD_TEST(suite, test_splice_0);
    SUITE_ADD_TEST(suite, test_processSplice_0);
    return suite;
}
