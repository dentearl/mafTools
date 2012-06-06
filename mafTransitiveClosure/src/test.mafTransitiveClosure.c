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
static void test_reverseComplement(CuTest *testCase) {
    // test that the reverseComplement function works properly
    char *input = de_strdup("GGGGaaaaaaaatttatatat");
    char *output = de_strdup("ATATATAAATTTTTTTTCCCC");
    reverseComplementSequence(input, strlen(input));
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
static void test_matrixAlignmentBlockComparisonOrdering_0(CuTest *testCase) {
    // test that with known input that known output is generated.
    char **input = (char**) de_malloc(2 * sizeof(char*));
    input[0] = de_strdup("AC---ACG-G");
    input[1] = de_strdup("ACTG--CGGG");
    mafTcComparisonOrder_t *obsOrder = NULL;
    obsOrder = getComparisonOrderFromMatrix(input, 2, 10);
    CuAssertTrue(testCase, obsOrder != NULL);
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
    CuAssertTrue(testCase, comparisonOrdersAreEqual(expectedOrder, obsOrder));
    // cleanup
    free(input[0]);
    free(input[1]);
    free(input);
    destroyMafTcComparisonOrder(expectedOrder);
    destroyMafTcComparisonOrder(obsOrder);
}
static void test_matrixAlignmentBlockComparisonOrdering_1(CuTest *testCase) {
    // test that with known input that known output is generated.
    char **input = (char**) de_malloc(4 * sizeof(char*));
    input[0] = de_strdup("acgcccag--------------ctcgcgaaatcgttc-------------------------ccggcattccgtccaggccgaagcgccaactgcagctccatctgtggcgtctctgttgcagcatcggagtgtgcaaatcatcccgacggccatcgtggtactcgtggtacacaggaaccaacaaataggagacgggtgccctgatcgacccgtgctccccggtgggcaccatcgatagattgctcgctgcggcgttccgtctcccc-----------------------ggtctgttcggcgactataaggtcacggacaggtaccttcaaaatcgacatggtgcttaaggtcag---------------------ccctaacttctccatgccttagctgacgatgttcgcgctaggtttaacgatatccgtcttgctgatgaacagtttcatcacccggcgacgatatccctggtgttgggctctgacgtttatcgtgatgtcatccaacccgggttcttaaatttggatga-gggctgcctgtcgcacaaagcac--tgtttggctggattatctcgggatcatgTAGACAATCTTAAGGCGCTAGACGTGTGGAAAATGTTGTTGCTTCGTTTGTCGTTGCA");
    input[1] = de_strdup("ccgctcggctttcaaggtcagtctcccggtttccttcaccatcggcaaggcaagcgcgtacaccggtatgacgtcc-------agcttcagtggcaaccccatc-------tcccgacggc-------------------------------catcgtgatgctggagt----ctggggccaccaccttcgagacggccgcgttgatcgacccgtgcactcccgtcagcaccatcgacagctccctggcaactgcattcaagttacccacgacgacagtgagaggtgaagaagtctgctcgtcgacgatccggtcaagaacgggtgatttccAGATCGACGTGCTCCTAAAGATAAA----------------------------------GCGCGCTTAGTGAAACTATGCGAGCCCAATTCAACGACATCCGTCTTGCCGATGAGCAGTTCCATCGCCCATATACAGTCTCGCTGGTGTTGGGCTCGGACGTATACCCCGACGTAATCCGGCCCGGGTTCCTAAATATACGTGATGGGCTGTCCGTCGCACAGGGCATG-TATTTGGTTGGGTAGTGTCTGGAGCATGCAGACACGCCTAAGG-GTTAATC-------------CGTTGTTACGCTCGCATTCGCA");
    input[2] = de_strdup("CCGCCCACATACCAGCAGTAGTCTCCCGGTCTCCTTCGCCAA------GGCAAGCGCTTACGCTGGTACGAcgtcc-------agcttcagcagcaaccccataagtagcttcgctactgc-------------------------------catcgtgttgctagagt----ctaggaccaccaactttgagacggcggcgttgatcgacccgtgcacttccgtcagcaccatagacagctccctggcagctgcgttcaagttacccacgacgacagtgagatgcgaagaagtctgttcgacgaccatccagtcaagaagaggcgatttccaaatcgacgtgctcctgaatattagccgaagtctacgcatccggaccc------------------------------------------------------------------------------------------------------------------------------gatccggcCCGGGCTCCTAAATATACGTGATACAT--------ATACAGGGCACGGTATTTGGATGAATCGTGTGTGGAGCACGCAGACACGCCTAAGG-ATTAGTC-------------CTTTGCTACGCTCGCCCTCGCA");
    input[3] = de_strdup("ccactcggctgtcaagggcagtctcgcggtctccatcgccagcggcacgggaagcacggacaccggtccgacgtcc-------agcgtcagcagcaactccatcagtggcggctctcctgcagcatcagagcctccacattcttccgacggctatcgtgttgctggagt----ctggggccaccaccttcgagacggctgcgttaatcgacccttgcacgcccgtcagcaccatcgacagctccctggcaactgcgttcaagttgcccacgacgacagtgagaggcgaagaagtctgctcgacgacaatccggtcgaggacgggcgatttccagatcgacgtgctgctgaagatcagtcagagtctacgaatccgcaccccta-----tccgcgcgctaagcgattctatgcgagcccaatttgatgatatccgtctggccgatgagcagttccatcgcccagcgacagtctcgctggtgttgggctcggacgtataccccgacgtaatcaggcccgggttcctaaatatacgagatgggctgcccgtcgcacagggcacggtctttggatgggtcgtgtctggagcatgcAGACACGCCTAAGG-ATTAATC-------------CTTTGCTACGTTCGCCCTCGCA");
    mafTcComparisonOrder_t *obsOrder = NULL;
    obsOrder = getComparisonOrderFromMatrix(input, 4, 648);
    CuAssertTrue(testCase, obsOrder != NULL);
    mafTcComparisonOrder_t *expectedOrder = newMafTcComparisonOrder();
    mafTcComparisonOrder_t *eo = expectedOrder;
    eo->ref = 2;
    eo->region = newMafTcRegion(561, 561);
    
    eo->next = newMafTcComparisonOrder();
    eo = eo->next;
    eo->ref = 2;
    eo->region = newMafTcRegion(357, 377);
    
    eo->next = newMafTcComparisonOrder();
    eo = eo->next;
    eo->ref = 1;
    eo->region = newMafTcRegion(560, 560);
    
    eo->next = newMafTcComparisonOrder();
    eo = eo->next;
    eo->ref = 1;
    eo->region = newMafTcRegion(536, 536);
    
    eo->next = newMafTcComparisonOrder();
    eo = eo->next;
    eo->ref = 1;
    eo->region = newMafTcRegion(268, 290);
    
    eo->next = newMafTcComparisonOrder();
    eo = eo->next;
    eo->ref = 1;
    eo->region = newMafTcRegion(37, 61);
    
    eo->next = newMafTcComparisonOrder();
    eo = eo->next;
    eo->ref = 1;
    eo->region = newMafTcRegion(8, 21);
    
    eo->next = newMafTcComparisonOrder();
    eo = eo->next;
    eo->ref = 0;
    eo->region = newMafTcRegion(562, 647);
    
    eo->next = newMafTcComparisonOrder();
    eo = eo->next;
    eo->ref = 0;
    eo->region = newMafTcRegion(537, 559);
    
    eo->next = newMafTcComparisonOrder();
    eo = eo->next;
    eo->ref = 0;
    eo->region = newMafTcRegion(378, 535);
    
    eo->next = newMafTcComparisonOrder();
    eo = eo->next;
    eo->ref = 0;
    eo->region = newMafTcRegion(291, 356);
    
    eo->next = newMafTcComparisonOrder();
    eo = eo->next;
    eo->ref = 0;
    eo->region = newMafTcRegion(62, 267);
    
    eo->next = newMafTcComparisonOrder();
    eo = eo->next;
    eo->ref = 0;
    eo->region = newMafTcRegion(22, 36);
    
    eo->next = newMafTcComparisonOrder();
    eo = eo->next;
    eo->ref = 0;
    eo->region = newMafTcRegion(0, 7);
    
    CuAssertTrue(testCase, comparisonOrdersAreEqual(expectedOrder, obsOrder));
    // cleanup
    for (int i = 0; i < 4; ++i)
        free(input[i]);
    free(input);
    destroyMafTcComparisonOrder(expectedOrder);
    destroyMafTcComparisonOrder(obsOrder);
}
static void test_matrixAlignmentBlockComparisonOrdering_2(CuTest *testCase) {
    // test that with known input that known output is generated.
    char **input = (char**) de_malloc(4 * sizeof(char*));
    input[0] = de_strdup("acgcccag--------------ctcgc---aatcgtt");
    input[1] = de_strdup("ccgctc--------------ctttcaag-tcagtctc");
    input[2] = de_strdup("CCGC--------CAGCAGTAGTCTCCC---CTCCTTC");
    input[3] = de_strdup("ccactcggctgtca---------tcgcggtctccatc");
    mafTcComparisonOrder_t *obsOrder = NULL;
    obsOrder = getComparisonOrderFromMatrix(input, 4, 37);
    CuAssertTrue(testCase, obsOrder != NULL);
    mafTcComparisonOrder_t *expectedOrder = newMafTcComparisonOrder();
    mafTcComparisonOrder_t *eo = expectedOrder;
    eo->ref = 3;
    eo->region = newMafTcRegion(28, 28);
    
    eo->next = newMafTcComparisonOrder();
    eo = eo->next;
    eo->ref = 3;
    eo->region = newMafTcRegion(8, 11);
    
    eo->next = newMafTcComparisonOrder();
    eo = eo->next;
    eo->ref = 2;
    eo->region = newMafTcRegion(12, 19);
    
    eo->next = newMafTcComparisonOrder();
    eo = eo->next;
    eo->ref = 1;
    eo->region = newMafTcRegion(29, 29);
    
    eo->next = newMafTcComparisonOrder();
    eo = eo->next;
    eo->ref = 1;
    eo->region = newMafTcRegion(27, 27);
    
    eo->next = newMafTcComparisonOrder();
    eo = eo->next;
    eo->ref = 1;
    eo->region = newMafTcRegion(20, 21);
    
    eo->next = newMafTcComparisonOrder();
    eo = eo->next;
    eo->ref = 0;
    eo->region = newMafTcRegion(30, 36);
    
    eo->next = newMafTcComparisonOrder();
    eo = eo->next;
    eo->ref = 0;
    eo->region = newMafTcRegion(22, 26);
    
    eo->next = newMafTcComparisonOrder();
    eo = eo->next;
    eo->ref = 0;
    eo->region = newMafTcRegion(0, 7);
    
    CuAssertTrue(testCase, comparisonOrdersAreEqual(expectedOrder, obsOrder));
    // cleanup
    for (int i = 0; i < 4; ++i)
        free(input[i]);
    free(input);
    destroyMafTcComparisonOrder(expectedOrder);
    destroyMafTcComparisonOrder(obsOrder);
}
static void test_matrixAlignmentBlockComparisonOrdering_3(CuTest *testCase) {
    // test that with known input that known output is generated.
    char **input = (char**) de_malloc(5 * sizeof(char*));
    input[0] = de_strdup("AATTG-----TCTCTCCCC--CTTTTT");
    input[1] = de_strdup("AATTGTC-----TCTGGCC--TTAATT");
    input[2] = de_strdup("CCCGGAGAG-----ACAAC--CTAATT");
    input[3] = de_strdup("ATTTAAATTTA-----GAG--ACAATC");
    input[4] = de_strdup("CCGGC-------GGTTGGGTTGTCTCT");
    mafTcComparisonOrder_t *obsOrder = NULL;
    obsOrder = getComparisonOrderFromMatrix(input, 5, 27);
    CuAssertTrue(testCase, obsOrder != NULL);
    mafTcComparisonOrder_t *expectedOrder = newMafTcComparisonOrder();
    mafTcComparisonOrder_t *eo = expectedOrder;
    eo->ref = 4;
    eo->region = newMafTcRegion(19, 20);
    
    eo->next = newMafTcComparisonOrder();
    eo = eo->next;
    eo->ref = 3;
    eo->region = newMafTcRegion(9, 9);
    
    eo->next = newMafTcComparisonOrder();
    eo = eo->next;
    eo->ref = 2;
    eo->region = newMafTcRegion(7, 8);
    
    eo->next = newMafTcComparisonOrder();
    eo = eo->next;
    eo->ref = 1;
    eo->region = newMafTcRegion(5, 6);
    
    eo->next = newMafTcComparisonOrder();
    eo = eo->next;
    eo->ref = 0;
    eo->region = newMafTcRegion(21, 26);
    
    eo->next = newMafTcComparisonOrder();
    eo = eo->next;
    eo->ref = 0;
    eo->region = newMafTcRegion(10, 18);
    
    eo->next = newMafTcComparisonOrder();
    eo = eo->next;
    eo->ref = 0;
    eo->region = newMafTcRegion(0, 4);
    
    CuAssertTrue(testCase, comparisonOrdersAreEqual(expectedOrder, obsOrder));
    // cleanup
    for (int i = 0; i < 5; ++i)
        free(input[i]);
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
    mafCoordinatePair_t mcp;
    mcp.a = 0;
    mcp.b = 0;
    char *s = de_strdup("ACGT");
    CuAssertTrue(testCase, localSeqCoords(3, s, &mcp) == 3);
    CuAssertTrue(testCase, mcp.a == 3);
    CuAssertTrue(testCase, mcp.b == 3);
    free(s);
    mcp.a = 0;
    mcp.b = 0;
    s = de_strdup("-------ACGT");
    CuAssertTrue(testCase, localSeqCoords(10, s, &mcp) == 3);
    CuAssertTrue(testCase, mcp.a == 10);
    CuAssertTrue(testCase, mcp.b == 3);
    free(s);
    mcp.a = 8;
    mcp.b = 1;
    s = de_strdup("-------ACGT");
    CuAssertTrue(testCase, localSeqCoords(10, s, &mcp) == 3);
    CuAssertTrue(testCase, mcp.a == 10);
    CuAssertTrue(testCase, mcp.b == 3);
    free(s);
    mcp.a = 10;
    mcp.b = 3;
    s = de_strdup("-------ACGT");
    CuAssertTrue(testCase, localSeqCoords(10, s, &mcp) == 3);
    CuAssertTrue(testCase, mcp.a == 10);
    CuAssertTrue(testCase, mcp.b == 3);
    free(s);
    mcp.a = 0;
    mcp.b = 0;
    s = de_strdup("-A-C-G-T");
    CuAssertTrue(testCase, localSeqCoords(7, s, &mcp) == 3);
    CuAssertTrue(testCase, mcp.a == 7);
    CuAssertTrue(testCase, mcp.b == 3);
    free(s);
    mcp.a = 0;
    mcp.b = 0;
    s = de_strdup("AA-A-C-G-T");
    CuAssertTrue(testCase, localSeqCoords(9, s, &mcp) == 5);
    CuAssertTrue(testCase, mcp.a == 9);
    CuAssertTrue(testCase, mcp.b == 5);
    free(s);
    mcp.a = 0;
    mcp.b = 0;
    s = de_strdup("AA-A-C-G-T");
    CuAssertTrue(testCase, localSeqCoords(3, s, &mcp) == 2);
    CuAssertTrue(testCase, mcp.a == 3);
    CuAssertTrue(testCase, mcp.b == 2);
    free(s);
    mcp.a = 1;
    mcp.b = 1;
    s = de_strdup("AA-A-C-G-T");
    CuAssertTrue(testCase, localSeqCoords(3, s, &mcp) == 2);
    CuAssertTrue(testCase, mcp.a == 3);
    CuAssertTrue(testCase, mcp.b == 2);
    free(s);
    mcp.a = 0;
    mcp.b = 0;
    s = de_strdup("GTTGTCTCTCAATGTG");
    CuAssertTrue(testCase, localSeqCoords(6, s, &mcp) == 6);
    CuAssertTrue(testCase, mcp.a == 6);
    CuAssertTrue(testCase, mcp.b == 6);
    free(s);
    s = de_strdup("GTTGTCTCTCAATGTG");
    CuAssertTrue(testCase, localSeqCoords(15, s, &mcp) == 15);
    CuAssertTrue(testCase, mcp.a == 15);
    CuAssertTrue(testCase, mcp.b == 15);
    free(s);
}
static void test_localSeqCoordsToGlobalPositiveCoords_0(CuTest *testCase) {
    // int64_t localSeqCoordsToGlobalPositiveCoords(localPosition, startField, sourceLength, strand);
    CuAssertTrue(testCase, localSeqCoordsToGlobalPositiveCoords(3, 0, 20, '+') == 3);
    CuAssertTrue(testCase, localSeqCoordsToGlobalPositiveCoords(3, 5, 20, '+') == 8);
    CuAssertTrue(testCase, localSeqCoordsToGlobalPositiveCoords(0, 0, 20, '-') == 19);
    CuAssertTrue(testCase, localSeqCoordsToGlobalPositiveCoords(3, 0, 20, '-') == 16);
    CuAssertTrue(testCase, localSeqCoordsToGlobalPositiveCoords(3, 5, 20, '-') == 11);
}
static void test_localSeqCoordsToGlobalPositiveStartCoords_0(CuTest *testCase) {
    // int64_t localSeqCoordsToGlobalPositiveStartCoords(localPosition, startField, sourceLength, 
    //                                                   strand, lengthField);
    CuAssertTrue(testCase, localSeqCoordsToGlobalPositiveStartCoords(3, 0, 20, '+', 5) == 3);
    CuAssertTrue(testCase, localSeqCoordsToGlobalPositiveStartCoords(3, 5, 20, '+', 5) == 8);
    CuAssertTrue(testCase, localSeqCoordsToGlobalPositiveStartCoords(3, 5, 20, '+', 10) == 8);
    CuAssertTrue(testCase, localSeqCoordsToGlobalPositiveStartCoords(0, 1, 100, '+', 16) == 1);
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
    mafCoordinatePair_t mcp;
    uint32_t start = 0, sourceLength = 100, seqLength = 37;
    mcp.a = 0;
    mcp.b = 0;
    CuAssertTrue(testCase, localSeqCoordsToGlobalPositiveStartCoords(localSeqCoords(1, input, &mcp), 
                                                                     start, sourceLength, '+', seqLength) == 0);
    CuAssertTrue(testCase, mcp.a == 1);
    CuAssertTrue(testCase, mcp.b == 0);
    CuAssertTrue(testCase, localSeqCoordsToGlobalPositiveStartCoords(localSeqCoords(4, input, &mcp), 
                                                                     start, sourceLength, '+', seqLength) == 2);
    CuAssertTrue(testCase, mcp.a == 4);
    CuAssertTrue(testCase, mcp.b == 2);
    CuAssertTrue(testCase, localSeqCoordsToGlobalPositiveStartCoords(localSeqCoords(26, input, &mcp), 
                                                                     start, sourceLength, '+', seqLength) == 21);
    CuAssertTrue(testCase, mcp.a == 26);
    CuAssertTrue(testCase, mcp.b == 21);
    CuAssertTrue(testCase, localSeqCoordsToGlobalPositiveStartCoords(localSeqCoords(41, input, &mcp), 
                                                                     start, sourceLength, '+', seqLength) == 36);
    CuAssertTrue(testCase, mcp.a == 41);
    CuAssertTrue(testCase, mcp.b == 36);
    // reverse strand.
    mcp.a = 0;
    mcp.b = 0;
    seqLength = 2;
    CuAssertTrue(testCase, localSeqCoordsToGlobalPositiveStartCoords(localSeqCoords(1, input, &mcp), 
                                                                     start, sourceLength, '-', seqLength) == 98);
    CuAssertTrue(testCase, mcp.a == 1);
    CuAssertTrue(testCase, mcp.b == 0);    
    seqLength = 19;
    CuAssertTrue(testCase, localSeqCoordsToGlobalPositiveStartCoords(localSeqCoords(4, input, &mcp), 
                                                                     start, sourceLength, '-', seqLength) == 79);
    CuAssertTrue(testCase, mcp.a == 4);
    CuAssertTrue(testCase, mcp.b == 2);
    seqLength = 16;
    CuAssertTrue(testCase, localSeqCoordsToGlobalPositiveStartCoords(localSeqCoords(26, input, &mcp), 
                                                                     start, sourceLength, '-', seqLength) == 63);
    CuAssertTrue(testCase, mcp.a == 26);
    CuAssertTrue(testCase, mcp.b == 21);
    free(input);
    input = de_strdup("GTTGTCTCTCAATGTG");
    mcp.a = 0;
    mcp.b = 0;
    start = 1;
    seqLength = 16;
    CuAssertTrue(testCase, localSeqCoordsToGlobalPositiveStartCoords(localSeqCoords(0, input, &mcp), 
                                                                     start, sourceLength, '+', seqLength) == 1);
    free(input);
}
CuSuite* mafTransitiveClosure_TestSuite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, test_reverseComplement);
    SUITE_ADD_TEST(suite, test_rowAlignmentBlockComparisonOrdering_0);
    SUITE_ADD_TEST(suite, test_rowAlignmentBlockComparisonOrdering_1);
    SUITE_ADD_TEST(suite, test_rowAlignmentBlockComparisonOrdering_2);
    SUITE_ADD_TEST(suite, test_matrixAlignmentBlockComparisonOrdering_0);
    SUITE_ADD_TEST(suite, test_matrixAlignmentBlockComparisonOrdering_1);
    SUITE_ADD_TEST(suite, test_matrixAlignmentBlockComparisonOrdering_2);
    SUITE_ADD_TEST(suite, test_matrixAlignmentBlockComparisonOrdering_3);
    SUITE_ADD_TEST(suite, test_addSequenceValuesToMtcSeq_0);
    SUITE_ADD_TEST(suite, test_localSeqCoords_0);
    SUITE_ADD_TEST(suite, test_localSeqCoordsToGlobalPositiveCoords_0);
    SUITE_ADD_TEST(suite, test_localSeqCoordsToGlobalPositiveStartCoords_0);
    SUITE_ADD_TEST(suite, test_coordinateTransforms_0);
    return suite;
}
