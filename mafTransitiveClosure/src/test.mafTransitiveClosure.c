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
#include "stPinchGraphs.h"
#include "mafTransitiveClosure.h"

bool regionListsAreEqual(mafTcRegion_t *expected, mafTcRegion_t *todo) {
    while (expected != NULL && todo != NULL) {
        // printf("checking expected against todo\n");
        if (expected->start != todo->start) {
            // fprintf(stderr, "expected->start %d != todo->start %d\n", expected->start, todo->start);
            return false;
        }
        if (expected->end != todo->end) {
            // fprintf(stderr, "expected->end %d != todo->end %d\n", expected->end, todo->end);
            return false;
        }
        // printf("okay pass\n");
        expected = expected->next;
        todo = todo->next;
    }
    if (expected != NULL) {
        return false;
    }
    if (todo != NULL) {
        return false;
    }
    return true;
}
bool comparisonOrdersAreEqual(mafTcComparisonOrder_t *eo, mafTcComparisonOrder_t *oo) {
    // printf("comparisonOrdersAreEqual()\n");
    while (eo != NULL && oo != NULL) {
        // printf("checking eo against oo\n");
        if (eo->ref != oo->ref) {
            // fprintf(stderr, "expected->ref %d != todo->ref %d\n", eo->ref, oo->ref);
            return false;
        }
        if (!regionListsAreEqual(eo->region, oo->region)){
            return false;
        }
        eo = eo->next;
        oo = oo->next;
    }
    if (eo != NULL) {
        return false;
    }
    if (oo != NULL) {
        return false;
    }
    return true;
}
void printRegion(mafTcRegion_t *reg) {
    printf("printRegion()\n");
    while (reg != NULL) {
        printf("  s: %d e: %d\n", reg->start, reg->end);
        reg = reg->next;
    }
}
void printComparisonOrder(mafTcComparisonOrder_t *co) {
    printf("printComparisonOrder()\n");
    while (co != NULL) {
        printf("  ref: %d ", co->ref);
        printRegion(co->region);
        co = co->next;
    }
}
static void test_reverseComplement(CuTest *testCase) {
    // test that the reverseComplement function works properly
    char *input = de_strdup("GGGGaaaaaaaatttatatat");
    char *output = de_strdup("ATATATAAATTTTTTTTCCCC");
    reverseComplementSequence(input);
    CuAssertTrue(testCase, strlen(input) == strlen(output));
    CuAssertTrue(testCase, strncmp(input, output, strlen(output)) == 0);
    free(input);
    free(output);
}
static void test_rowAlignmentBlockComparisonOrdering_0(CuTest *testCase) {
    // test that with known input that known output is generated.
    char **input = (char**) de_malloc(sizeof(char*));
    input[0] = de_strdup("AT");
    mafTcComparisonOrder_t *obsOrder = NULL;
    mafTcRegion_t *todo = newMafTcRegion(0, 1);
    todo = getComparisonOrderFromRow(input, 0, &obsOrder, todo);
    CuAssertTrue(testCase, todo == NULL);
    CuAssertTrue(testCase, obsOrder != NULL);
    mafTcComparisonOrder_t *expectedOrder = newMafTcComparisonOrder();
    mafTcComparisonOrder_t *eo = expectedOrder;
    eo->ref = 0;
    eo->region = newMafTcRegion(0, 1);
    CuAssertTrue(testCase, comparisonOrdersAreEqual(expectedOrder, obsOrder));
    // cleanup
    free(input[0]);
    free(input);
    destroyMafTcComparisonOrder(expectedOrder);
    destroyMafTcComparisonOrder(obsOrder);
    destroyMafTcRegionList(todo);
}
static void test_rowAlignmentBlockComparisonOrdering_1(CuTest *testCase) {
    // test that with known input that known output is generated.
    char **input = (char**) de_malloc(sizeof(char*));
    input[0] = de_strdup("AC---ACG-G");
    mafTcComparisonOrder_t *obsOrder = NULL;
    mafTcRegion_t *todo = newMafTcRegion(0, 9);
    todo = getComparisonOrderFromRow(input, 0, &obsOrder, todo);
    CuAssertTrue(testCase, todo != NULL);
    CuAssertTrue(testCase, obsOrder != NULL);
    mafTcRegion_t *expectedTodo = newMafTcRegion(2, 4);
    expectedTodo->next = newMafTcRegion(8, 8);
    CuAssertTrue(testCase, expectedTodo != NULL);
    CuAssertTrue(testCase, regionListsAreEqual(expectedTodo, todo));
    mafTcComparisonOrder_t *expectedOrder = newMafTcComparisonOrder();
    mafTcComparisonOrder_t *eo = expectedOrder;
    eo->ref = 0;
    eo->region = newMafTcRegion(9, 9);
    eo->next = newMafTcComparisonOrder();
    eo = eo->next;
    eo->ref = 0;
    eo->region = newMafTcRegion(5, 7);
    eo->next = newMafTcComparisonOrder();
    eo = eo->next;
    eo->ref = 0;
    eo->region = newMafTcRegion(0, 1);
    eo = eo->next;
    CuAssertTrue(testCase, eo == NULL);
    CuAssertTrue(testCase, comparisonOrdersAreEqual(expectedOrder, obsOrder));
    // cleanup
    free(input[0]);
    free(input);
    destroyMafTcComparisonOrder(expectedOrder);
    destroyMafTcComparisonOrder(obsOrder);
    destroyMafTcRegionList(expectedTodo);
    destroyMafTcRegionList(todo);
}
static void test_rowAlignmentBlockComparisonOrdering_2(CuTest *testCase) {
    // test that with known input that known output is generated.
    char **input = (char**) de_malloc(2 * sizeof(char*));
    input[0] = de_strdup("AC---ACG-G");
    input[1] = de_strdup("ACTG--CGGG");
    mafTcComparisonOrder_t *obsOrder = NULL;
    mafTcRegion_t *todo = newMafTcRegion(0, 9);
    todo = getComparisonOrderFromRow(input, 0, &obsOrder, todo);
    CuAssertTrue(testCase, todo != NULL);
    CuAssertTrue(testCase, obsOrder != NULL);
    todo = getComparisonOrderFromRow(input, 1, &obsOrder, todo);
    mafTcRegion_t *expectedTodo = newMafTcRegion(4, 4);
    CuAssertTrue(testCase, expectedTodo != NULL);
    CuAssertTrue(testCase, regionListsAreEqual(expectedTodo, todo));
    // printf("obsOrder: ");
    // printComparisonOrder(obsOrder);
    mafTcComparisonOrder_t *expectedOrder = newMafTcComparisonOrder();
    mafTcComparisonOrder_t *eo = expectedOrder;
    eo->ref = 1;
    eo->region = newMafTcRegion(8, 8);
    eo->next = newMafTcComparisonOrder();
    eo = eo->next;
    eo->ref = 1;
    eo->region = newMafTcRegion(2, 3);
    eo->next = newMafTcComparisonOrder();
    eo = eo->next;
    eo->ref = 0;
    eo->region = newMafTcRegion(9, 9);
    eo->next = newMafTcComparisonOrder();
    eo = eo->next;
    eo->ref = 0;
    eo->region = newMafTcRegion(5, 7);
    eo->next = newMafTcComparisonOrder();
    eo = eo->next;
    eo->ref = 0;
    eo->region = newMafTcRegion(0, 1);
    eo = eo->next;
    // printf("expected: ");
    // printComparisonOrder(expectedOrder);
    CuAssertTrue(testCase, comparisonOrdersAreEqual(expectedOrder, obsOrder));
    // cleanup
    free(input[0]);
    free(input[1]);
    free(input);
    destroyMafTcComparisonOrder(expectedOrder);
    destroyMafTcComparisonOrder(obsOrder);
    destroyMafTcRegionList(expectedTodo);
    destroyMafTcRegionList(todo);
}
static void test_matrixAlignmentBlockComparisonOrdering_1(CuTest *testCase) {
    // test that with known input that known output is generated.
    char **input = (char**) de_malloc(2 * sizeof(char*));
    input[0] = de_strdup("AC---ACG-G");
    input[1] = de_strdup("ACTG--CGGG");
    mafTcComparisonOrder_t *obsOrder = NULL;
    obsOrder = getComparisonOrderFromMatrix(input, 2, 10);
    CuAssertTrue(testCase, obsOrder != NULL);
    // printf("obsOrder: ");
    // printComparisonOrder(obsOrder);
    mafTcComparisonOrder_t *expectedOrder = newMafTcComparisonOrder();
    mafTcComparisonOrder_t *eo = expectedOrder;
    eo->ref = 1;
    eo->region = newMafTcRegion(8, 8);
    eo->next = newMafTcComparisonOrder();
    eo = eo->next;
    eo->ref = 1;
    eo->region = newMafTcRegion(2, 3);
    eo->next = newMafTcComparisonOrder();
    eo = eo->next;
    eo->ref = 0;
    eo->region = newMafTcRegion(9, 9);
    eo->next = newMafTcComparisonOrder();
    eo = eo->next;
    eo->ref = 0;
    eo->region = newMafTcRegion(5, 7);
    eo->next = newMafTcComparisonOrder();
    eo = eo->next;
    eo->ref = 0;
    eo->region = newMafTcRegion(0, 1);
    eo = eo->next;
    // printf("expected: ");
    // printComparisonOrder(expectedOrder);
    CuAssertTrue(testCase, comparisonOrdersAreEqual(expectedOrder, obsOrder));
    // cleanup
    free(input[0]);
    free(input[1]);
    free(input);
    destroyMafTcComparisonOrder(expectedOrder);
    destroyMafTcComparisonOrder(obsOrder);
}
static void test_addSequenceValuesToMtcSeq_0(CuTest *testCase) {
    mafLine_t *ml = maf_newMafLine();
    maf_mafLine_setType(ml, 's');
    maf_mafLine_setStrand(ml, '+');
    maf_mafLine_setStart(ml, 1);
    maf_mafLine_setLength(ml, 8);
    maf_mafLine_setSourceLength(ml, 20);
    maf_mafLine_setSequence(ml, de_strdup("ACGT----ACGTTT"));
    mafTcSeq_t *mtcs = newMafTcSeq(de_strdup("test.chr2"), maf_mafLine_getSourceLength(ml));
    CuAssertTrue(testCase, strcmp(mtcs->sequence, "NNNNNNNNNNNNNNNNNNNN") == 0);
    addSequenceValuesToMtcSeq(ml, mtcs);
    CuAssertTrue(testCase, strcmp(mtcs->sequence, "NACGTACGTTTNNNNNNNNN") == 0);
    maf_destroyMafLineList(ml);
    ml = maf_newMafLine();
    maf_mafLine_setType(ml, 's');
    maf_mafLine_setStrand(ml, '-');
    maf_mafLine_setStart(ml, 1);
    maf_mafLine_setLength(ml, 5);
    maf_mafLine_setSourceLength(ml, 20);
    maf_mafLine_setSequence(ml, de_strdup("AATT------G"));
    addSequenceValuesToMtcSeq(ml, mtcs);
    CuAssertTrue(testCase, strcmp(mtcs->sequence, "NACGTACGTTTNNNCAATTN") == 0);
    maf_destroyMafLineList(ml);
    ml = maf_newMafLine();
    maf_mafLine_setType(ml, 's');
    maf_mafLine_setStrand(ml, '+');
    maf_mafLine_setStart(ml, 10);
    maf_mafLine_setLength(ml, 5);
    maf_mafLine_setSourceLength(ml, 20);
    maf_mafLine_setSequence(ml, de_strdup("--T----G-G-G-C-"));
    addSequenceValuesToMtcSeq(ml, mtcs);
    CuAssertTrue(testCase, strcmp(mtcs->sequence, "NACGTACGTTTGGGCAATTN") == 0);
    maf_destroyMafLineList(ml);
    ml = maf_newMafLine();
    maf_mafLine_setType(ml, 's');
    maf_mafLine_setStrand(ml, '+');
    maf_mafLine_setStart(ml, 0);
    maf_mafLine_setLength(ml, 20);
    maf_mafLine_setSourceLength(ml, 20);
    maf_mafLine_setSequence(ml, de_strdup("AACGTACGTTTGGGCAATTG"));
    addSequenceValuesToMtcSeq(ml, mtcs);
    CuAssertTrue(testCase, strcmp(mtcs->sequence, "AACGTACGTTTGGGCAATTG") == 0);
    maf_destroyMafLineList(ml);
    ml = maf_newMafLine();
    maf_mafLine_setType(ml, 's');
    maf_mafLine_setStrand(ml, '-');
    maf_mafLine_setStart(ml, 0);
    maf_mafLine_setLength(ml, 20);
    maf_mafLine_setSourceLength(ml, 20);
    maf_mafLine_setSequence(ml, de_strdup("---CAATTGCCC---AAACGTAC----GTT"));
    addSequenceValuesToMtcSeq(ml, mtcs);
    CuAssertTrue(testCase, strcmp(mtcs->sequence, "AACGTACGTTTGGGCAATTG") == 0);
    maf_destroyMafLineList(ml);
    destroyMafTcSeq(mtcs);
}
static void test_localSeqCoords_0(CuTest *testCase) {
    char *s = de_strdup("ACGT");
    CuAssertTrue(testCase, localSeqCoords(3, s) == 3);
    free(s);
    s = de_strdup("-------ACGT");
    CuAssertTrue(testCase, localSeqCoords(10, s) == 3);
    free(s);
    s = de_strdup("-A-C-G-T");
    CuAssertTrue(testCase, localSeqCoords(7, s) == 3);
    free(s);
    s = de_strdup("AA-A-C-G-T");
    CuAssertTrue(testCase, localSeqCoords(9, s) == 5);
    free(s);
    s = de_strdup("AA-A-C-G-T");
    CuAssertTrue(testCase, localSeqCoords(3, s) == 2);
    free(s);
}
static void test_localSeqCoordsToGlobalPositiveCoords_0(CuTest *testCase) {
    // int64_t localSeqCoordsToGlobalPositiveCoords(int64_t c, uint32_t start, uint32_t sourceLength, 
    //                                              char strand);
    CuAssertTrue(testCase, localSeqCoordsToGlobalPositiveCoords(3, 0, 20, '+') == 3);
    CuAssertTrue(testCase, localSeqCoordsToGlobalPositiveCoords(3, 5, 20, '+') == 8);
    CuAssertTrue(testCase, localSeqCoordsToGlobalPositiveCoords(0, 0, 20, '-') == 19);
    CuAssertTrue(testCase, localSeqCoordsToGlobalPositiveCoords(3, 0, 20, '-') == 16);
    CuAssertTrue(testCase, localSeqCoordsToGlobalPositiveCoords(3, 5, 20, '-') == 11);
}
static void test_localSeqCoordsToGlobalPositiveStartCoords_0(CuTest *testCase) {
    // int64_t localSeqCoordsToGlobalPositiveStartCoords(int64_t c, uint32_t start, uint32_t sourceLength, 
    //                                                   char strand, uint32_t length);
    CuAssertTrue(testCase, localSeqCoordsToGlobalPositiveStartCoords(3, 0, 20, '+', 5) == 3);
    CuAssertTrue(testCase, localSeqCoordsToGlobalPositiveStartCoords(3, 5, 20, '+', 5) == 8);
    CuAssertTrue(testCase, localSeqCoordsToGlobalPositiveStartCoords(0, 0, 20, '-', 1) == 19);
    CuAssertTrue(testCase, localSeqCoordsToGlobalPositiveStartCoords(0, 0, 20, '-', 5) == 15);
    CuAssertTrue(testCase, localSeqCoordsToGlobalPositiveStartCoords(3, 0, 20, '-', 1) == 16);
    CuAssertTrue(testCase, localSeqCoordsToGlobalPositiveStartCoords(3, 0, 20, '-', 5) == 12);
    CuAssertTrue(testCase, localSeqCoordsToGlobalPositiveStartCoords(3, 5, 20, '-', 1) == 11);
    CuAssertTrue(testCase, localSeqCoordsToGlobalPositiveStartCoords(3, 5, 20, '-', 5) == 7);
    CuAssertTrue(testCase, localSeqCoordsToGlobalPositiveStartCoords(0, 0, 20, '-', 10) == 10);
}
static void test_coordinateTransforms_0(CuTest *testCase) {
    char *input = de_strdup("-AA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG");
    CuAssertTrue(testCase, localSeqCoordsToGlobalPositiveStartCoords(localSeqCoords(1, input), 
                                                                     0, 100, '+', 2) == 0);
    free(input);
}
CuSuite* mafTransitiveClosure_TestSuite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, test_reverseComplement);
    SUITE_ADD_TEST(suite, test_rowAlignmentBlockComparisonOrdering_0);
    SUITE_ADD_TEST(suite, test_rowAlignmentBlockComparisonOrdering_1);
    SUITE_ADD_TEST(suite, test_rowAlignmentBlockComparisonOrdering_2);
    SUITE_ADD_TEST(suite, test_matrixAlignmentBlockComparisonOrdering_1);
    SUITE_ADD_TEST(suite, test_addSequenceValuesToMtcSeq_0);
    SUITE_ADD_TEST(suite, test_localSeqCoords_0);
    SUITE_ADD_TEST(suite, test_localSeqCoordsToGlobalPositiveCoords_0);
    SUITE_ADD_TEST(suite, test_localSeqCoordsToGlobalPositiveStartCoords_0);
    SUITE_ADD_TEST(suite, test_coordinateTransforms_0);
    return suite;
}
