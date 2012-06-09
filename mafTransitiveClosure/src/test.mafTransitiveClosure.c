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

void printRegionList(mafTcRegion_t *reg, FILE *ofp) {
    while (reg != NULL) {
        fprintf(ofp, "[%" PRIu32 ", %" PRIu32 "], ", reg->start, reg->end);
        reg = reg->next;
    }
    fprintf(ofp, "\n");
}
bool regionListsAreEqual(mafTcRegion_t *expected, mafTcRegion_t *obs, bool verbose) {
    if (verbose) {
        printf("  regionListsAreEqual()\n");
    }
    while (expected != NULL && obs != NULL) {
        if (expected->start != obs->start) {
            if (verbose) {
                fprintf(stderr, "    expected->start %" PRIu32" != obs->start %" PRIu32"\n", 
                        expected->start, obs->start);
                fprintf(stderr, "    expected: ");
                printRegionList(expected, stderr);
                fprintf(stderr, "    observed: ");
                printRegionList(obs, stderr);
            }
            return false;
        }
        if (expected->end != obs->end) {
            if (verbose) {
                fprintf(stderr, "    expected->end %" PRIu32 " != obs->end %" PRIu32 "\n", 
                        expected->end, obs->end);
                fprintf(stderr, "    expected: ");
                printRegionList(expected, stderr);
                fprintf(stderr, "    observed: ");
                printRegionList(obs, stderr);
            }
            return false;
        }
        if (verbose) {
            printf("    okay pass\n");
        }
        expected = expected->next;
        obs = obs->next;
    }
    if (expected != NULL) {
        if (verbose) {
            printf("    premature end of observed\n");
        }
        return false;
    }
    if (obs != NULL) {
        if (verbose) {
            printf("    premature end of expected\n");
        }
        return false;
    }
    return true;
}
static void printComparisonOrder(mafTcComparisonOrder_t *co) {
    fprintf(stderr, "printComparisonOrder()\n");
    while (co != NULL) {
        fprintf(stderr, " ref: %2" PRIu32 " \n", co->ref);
        printRegionList(co->region, stderr);
        co = co->next;
    }
}
bool comparisonOrdersAreEqual(mafTcComparisonOrder_t *eo, mafTcComparisonOrder_t *oo, bool verbose) {
    if (verbose) {
        printf("comparisonOrdersAreEqual()\n");
        printf("expected:\n");
        printComparisonOrder(eo);
        printf("observed:\n");
        printComparisonOrder(oo);
    }
    while (eo != NULL && oo != NULL) {
        if (eo->ref != oo->ref) {
            if (verbose) {
                fprintf(stderr, "expected->ref %d != obs->ref %d\n", eo->ref, oo->ref);
            }
            return false;
        }
        if (!regionListsAreEqual(eo->region, oo->region, verbose)){
            if (verbose) {
                fprintf(stderr, "regionListsAreEqual() = %s", 
                        regionListsAreEqual(eo->region, oo->region, verbose) ? "true" : "false");
            }
            return false;
        }
        eo = eo->next;
        oo = oo->next;
    }
    if (eo != NULL) {
        if (verbose) {
            printf("premature end of observed comparison order\n");
            printf("remaining expected:\n");
            printComparisonOrder(eo);
        }
        return false;
    }
    if (oo != NULL) {
        if (verbose) {
            printf("premature end of expected comparison order\n");
            printf("remaining observed:\n");
            printComparisonOrder(oo);
        }
        return false;
    }
    return true;
}
static bool mafBlocksAreEqual(mafBlock_t *input, mafBlock_t *expected, bool verbose) {
    mafLine_t *m1 = NULL, *m2 = NULL;
    m1 = maf_mafBlock_getHeadLine(input);
    m2 = maf_mafBlock_getHeadLine(expected);
    if (verbose) {
        printf("mafBlocksAreEqual():\n");
    }
    while (m1 != NULL) {
        if (m2 == NULL) { 
            return false;
        }
        if (verbose) {
            printf("Comparing lines:\n");
            printf("input   : %s\n", maf_mafLine_getLine(m1));
            printf("expected: %s\n", maf_mafLine_getLine(m2));
        }
        if (maf_mafLine_getType(m1) != maf_mafLine_getType(m2)) {
            if (verbose) {
                printf("Types differ, input:%c expected:%c\n", 
                       maf_mafLine_getType(m1), maf_mafLine_getType(m2));
            }
            return false;
        }
        if (maf_mafBlock_getNumberOfLines(input) != maf_mafBlock_getNumberOfLines(expected)) {
            if (verbose) {
                printf("Number of lines differ:\n  input:%" PRIu32 "\n  expected:%" PRIu32 "\n", 
                       maf_mafBlock_getNumberOfLines(input), maf_mafBlock_getNumberOfLines(expected));
            }
            return false;
        }
        if (strcmp(maf_mafLine_getLine(m1), maf_mafLine_getLine(m2)) != 0) {
            if (verbose) {
                printf("Lines differ:\n  input:%s\n  expected:%s\n", 
                       maf_mafLine_getLine(m1), maf_mafLine_getLine(m2));
            }
            return false;
        }
        if (maf_mafLine_getType(m1) != 's') {
            m1 = maf_mafLine_getNext(m1);
            m2 = maf_mafLine_getNext(m2);
            continue;
        }
        if (strcmp(maf_mafLine_getSpecies(m1), maf_mafLine_getSpecies(m2)) != 0) {
            if (verbose) {
                printf("Species differ:\n  input:%s\n  expected:%s\n", 
                       maf_mafLine_getSpecies(m1), maf_mafLine_getSpecies(m2));
            }
            return false;
        }
        if (maf_mafLine_getStart(m1) != maf_mafLine_getStart(m2)) {
            if (verbose) {
                printf("Starts differ:\n  input:%" PRIu32 "\n  expected:%" PRIu32"\n", 
                       maf_mafLine_getStart(m1), maf_mafLine_getStart(m2));
            }
            return false;
        }
        if (maf_mafLine_getLength(m1) != maf_mafLine_getLength(m2)) {
            if (verbose) {
                printf("Lengths differ:\n  input:%" PRIu32 "\n  expected:%" PRIu32"\n", 
                       maf_mafLine_getLength(m1), maf_mafLine_getLength(m2));
            }
            return false;
        }
        if (maf_mafLine_getStrand(m1) != maf_mafLine_getStrand(m2)) {
            if (verbose) {
                printf("Strands differ:\n  input:%c\n  expected:%c\n", 
                       maf_mafLine_getStrand(m1), maf_mafLine_getStrand(m2));
            }
            return false;
        }
        if (maf_mafLine_getSourceLength(m1) != maf_mafLine_getSourceLength(m2)) {
            if (verbose) {
                printf("Source lengths differ:\n  input:%" PRIu32 "\n  expected:%" PRIu32"\n", 
                       maf_mafLine_getSourceLength(m1), maf_mafLine_getSourceLength(m2));
            }
            return false;
        }
        m1 = maf_mafLine_getNext(m1);
        m2 = maf_mafLine_getNext(m2);
    }
    if (m2 != NULL) {
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
    todo = getComparisonOrderFromRow(input, 0, &obsOrder, todo, 0);
    CuAssertTrue(testCase, todo == NULL);
    CuAssertTrue(testCase, obsOrder != NULL);
    mafTcComparisonOrder_t *expectedOrder = newMafTcComparisonOrder();
    mafTcComparisonOrder_t *eo = expectedOrder;
    eo->ref = 0;
    eo->region = newMafTcRegion(0, 1);
    CuAssertTrue(testCase, comparisonOrdersAreEqual(expectedOrder, obsOrder, false));
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
    todo = getComparisonOrderFromRow(input, 0, &obsOrder, todo, 1);
    CuAssertTrue(testCase, todo != NULL);
    CuAssertTrue(testCase, obsOrder != NULL);
    mafTcRegion_t *expectedTodo = newMafTcRegion(2, 4);
    expectedTodo->next = newMafTcRegion(8, 8);
    CuAssertTrue(testCase, expectedTodo != NULL);
    CuAssertTrue(testCase, regionListsAreEqual(expectedTodo, todo, false));
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
    CuAssertTrue(testCase, comparisonOrdersAreEqual(expectedOrder, obsOrder, false));
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
    char **input = (char**) de_malloc(sizeof(char*) * 2);
    input[0] = de_strdup("AC---ACG-G");
    input[1] = de_strdup("ACTG--CGGG");
    mafTcComparisonOrder_t *obsOrder = NULL;
    mafTcRegion_t *todo = newMafTcRegion(0, 9);
    todo = getComparisonOrderFromRow(input, 0, &obsOrder, todo, 1);
    CuAssertTrue(testCase, todo != NULL);
    CuAssertTrue(testCase, obsOrder != NULL);
    todo = getComparisonOrderFromRow(input, 1, &obsOrder, todo, 1);
    mafTcRegion_t *expectedTodo = newMafTcRegion(4, 4);
    CuAssertTrue(testCase, expectedTodo != NULL);
    CuAssertTrue(testCase, regionListsAreEqual(expectedTodo, todo, false));
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
    CuAssertTrue(testCase, comparisonOrdersAreEqual(expectedOrder, obsOrder, false));
    // cleanup
    free(input[0]);
    free(input[1]);
    free(input);
    destroyMafTcComparisonOrder(expectedOrder);
    destroyMafTcComparisonOrder(obsOrder);
    destroyMafTcRegionList(expectedTodo);
    destroyMafTcRegionList(todo);
}
static void test_rowAlignmentBlockComparisonOrdering_3(CuTest *testCase) {
    // test that with known input that known output is generated.
    char **input = (char**) de_malloc(sizeof(char*));
    input[0] = de_strdup("AATTGTC-----TCTGGCC--TTAATT");
    mafTcComparisonOrder_t *obsOrder = NULL;
    mafTcRegion_t *todo = newMafTcRegion(5, 9);
    todo->next = newMafTcRegion(19, 20);
    todo = getComparisonOrderFromRow(input, 0, &obsOrder, todo, 1);
    CuAssertTrue(testCase, todo != NULL);
    CuAssertTrue(testCase, obsOrder != NULL);
    mafTcRegion_t *expectedTodo = newMafTcRegion(7, 9);
    expectedTodo->next = newMafTcRegion(19, 20);
    CuAssertTrue(testCase, expectedTodo != NULL);
    CuAssertTrue(testCase, regionListsAreEqual(expectedTodo, todo, false));
    mafTcComparisonOrder_t *expectedOrder = newMafTcComparisonOrder();
    mafTcComparisonOrder_t *eo = expectedOrder;
    eo->ref = 0;
    eo->region = newMafTcRegion(5, 6);
    eo = eo->next;
    CuAssertTrue(testCase, comparisonOrdersAreEqual(expectedOrder, obsOrder, false));
    // cleanup
    free(input[0]);
    free(input);
    destroyMafTcComparisonOrder(expectedOrder);
    destroyMafTcComparisonOrder(obsOrder);
    destroyMafTcRegionList(expectedTodo);
    destroyMafTcRegionList(todo);
}
static void test_matrixAlignmentBlockComparisonOrdering_0(CuTest *testCase) {
    // test that with known input that known output is generated.
    char **input = (char**) de_malloc(sizeof(char*) * 2);
    uint32_t *lengths = de_malloc(sizeof(uint32_t) * 2);
    input[0] = de_strdup("AC---ACG-G");
    input[1] = de_strdup("ACTG--CGGG");
    lengths[0] = 6;
    lengths[1] = 8;
    mafTcComparisonOrder_t *obsOrder = NULL;
    int **dummyVizMatrix = NULL;
    obsOrder = getComparisonOrderFromMatrix(input, 2, 10, lengths, dummyVizMatrix);
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
    CuAssertTrue(testCase, comparisonOrdersAreEqual(expectedOrder, obsOrder, false));
    // cleanup
    free(input[0]);
    free(input[1]);
    free(input);
    free(lengths);
    destroyMafTcComparisonOrder(expectedOrder);
    destroyMafTcComparisonOrder(obsOrder);
}
static void test_matrixAlignmentBlockComparisonOrdering_1(CuTest *testCase) {
    // test that with known input that known output is generated.
    char **input = (char**) de_malloc(sizeof(char*) * 4);
    uint32_t *lengths = de_malloc(sizeof(uint32_t) * 4);
    input[0] = de_strdup("acgcccag--------------ctcgcgaaatcgttc-------------------------ccggcattccgtccaggccgaagcgccaactgcagctccatctgtggcgtctctgttgcagcatcggagtgtgcaaatcatcccgacggccatcgtggtactcgtggtacacaggaaccaacaaataggagacgggtgccctgatcgacccgtgctccccggtgggcaccatcgatagattgctcgctgcggcgttccgtctcccc-----------------------ggtctgttcggcgactataaggtcacggacaggtaccttcaaaatcgacatggtgcttaaggtcag---------------------ccctaacttctccatgccttagctgacgatgttcgcgctaggtttaacgatatccgtcttgctgatgaacagtttcatcacccggcgacgatatccctggtgttgggctctgacgtttatcgtgatgtcatccaacccgggttcttaaatttggatga-gggctgcctgtcgcacaaagcac--tgtttggctggattatctcgggatcatgTAGACAATCTTAAGGCGCTAGACGTGTGGAAAATGTTGTTGCTTCGTTTGTCGTTGCA");
    input[1] = de_strdup("ccgctcggctttcaaggtcagtctcccggtttccttcaccatcggcaaggcaagcgcgtacaccggtatgacgtcc-------agcttcagtggcaaccccatc-------tcccgacggc-------------------------------catcgtgatgctggagt----ctggggccaccaccttcgagacggccgcgttgatcgacccgtgcactcccgtcagcaccatcgacagctccctggcaactgcattcaagttacccacgacgacagtgagaggtgaagaagtctgctcgtcgacgatccggtcaagaacgggtgatttccAGATCGACGTGCTCCTAAAGATAAA----------------------------------GCGCGCTTAGTGAAACTATGCGAGCCCAATTCAACGACATCCGTCTTGCCGATGAGCAGTTCCATCGCCCATATACAGTCTCGCTGGTGTTGGGCTCGGACGTATACCCCGACGTAATCCGGCCCGGGTTCCTAAATATACGTGATGGGCTGTCCGTCGCACAGGGCATG-TATTTGGTTGGGTAGTGTCTGGAGCATGCAGACACGCCTAAGG-GTTAATC-------------CGTTGTTACGCTCGCATTCGCA");
    input[2] = de_strdup("CCGCCCACATACCAGCAGTAGTCTCCCGGTCTCCTTCGCCAA------GGCAAGCGCTTACGCTGGTACGAcgtcc-------agcttcagcagcaaccccataagtagcttcgctactgc-------------------------------catcgtgttgctagagt----ctaggaccaccaactttgagacggcggcgttgatcgacccgtgcacttccgtcagcaccatagacagctccctggcagctgcgttcaagttacccacgacgacagtgagatgcgaagaagtctgttcgacgaccatccagtcaagaagaggcgatttccaaatcgacgtgctcctgaatattagccgaagtctacgcatccggaccc------------------------------------------------------------------------------------------------------------------------------gatccggcCCGGGCTCCTAAATATACGTGATACAT--------ATACAGGGCACGGTATTTGGATGAATCGTGTGTGGAGCACGCAGACACGCCTAAGG-ATTAGTC-------------CTTTGCTACGCTCGCCCTCGCA");
    input[3] = de_strdup("ccactcggctgtcaagggcagtctcgcggtctccatcgccagcggcacgggaagcacggacaccggtccgacgtcc-------agcgtcagcagcaactccatcagtggcggctctcctgcagcatcagagcctccacattcttccgacggctatcgtgttgctggagt----ctggggccaccaccttcgagacggctgcgttaatcgacccttgcacgcccgtcagcaccatcgacagctccctggcaactgcgttcaagttgcccacgacgacagtgagaggcgaagaagtctgctcgacgacaatccggtcgaggacgggcgatttccagatcgacgtgctgctgaagatcagtcagagtctacgaatccgcaccccta-----tccgcgcgctaagcgattctatgcgagcccaatttgatgatatccgtctggccgatgagcagttccatcgcccagcgacagtctcgctggtgttgggctcggacgtataccccgacgtaatcaggcccgggttcctaaatatacgagatgggctgcccgtcgcacagggcacggtctttggatgggtcgtgtctggagcatgcAGACACGCCTAAGG-ATTAATC-------------CTTTGCTACGTTCGCCCTCGCA");
    lengths[0] = 562;
    lengths[1] = 550;
    lengths[2] = 452;
    lengths[3] = 618;
    mafTcComparisonOrder_t *obsOrder = NULL;
    int **dummyVizMatrix = NULL;
    obsOrder = getComparisonOrderFromMatrix(input, 4, 648, lengths, dummyVizMatrix);
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
    
    CuAssertTrue(testCase, comparisonOrdersAreEqual(expectedOrder, obsOrder, false));
    // cleanup
    for (int i = 0; i < 4; ++i)
        free(input[i]);
    free(input);
    free(lengths);
    destroyMafTcComparisonOrder(expectedOrder);
    destroyMafTcComparisonOrder(obsOrder);
}
static void test_matrixAlignmentBlockComparisonOrdering_2(CuTest *testCase) {
    // test that with known input that known output is generated.
    char **input = (char**) de_malloc(sizeof(char*) * 4);
    uint32_t *lengths = de_malloc(sizeof(uint32_t) * 4);
    input[0] = de_strdup("acgcccag--------------ctcgc---aatcgtt");
    input[1] = de_strdup("ccgctc--------------ctttcaag-tcagtctc");
    input[2] = de_strdup("CCGC--------CAGCAGTAGTCTCCC---CTCCTTC");
    input[3] = de_strdup("ccactcggctgtca---------tcgcggtctccatc");
    lengths[0] = 20;
    lengths[1] = 22;
    lengths[2] = 26;
    lengths[3] = 28;
    mafTcComparisonOrder_t *obsOrder = NULL;
    int **dummyVizMatrix = NULL;
    obsOrder = getComparisonOrderFromMatrix(input, 4, 37, lengths, dummyVizMatrix);
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
    
    CuAssertTrue(testCase, comparisonOrdersAreEqual(expectedOrder, obsOrder, false));
    // cleanup
    for (int i = 0; i < 4; ++i)
        free(input[i]);
    free(input);
    free(lengths);
    destroyMafTcComparisonOrder(expectedOrder);
    destroyMafTcComparisonOrder(obsOrder);
}
static void test_matrixAlignmentBlockComparisonOrdering_3(CuTest *testCase) {
    // test that with known input that known output is generated.
    char **input = (char**) de_malloc(sizeof(char*) * 5);
    uint32_t *lengths = de_malloc(sizeof(uint32_t) * 5);
    input[0] = de_strdup("AATTG-----TCTCTCCCC--CTTTTT");
    input[1] = de_strdup("AATTGTC-----TCTGGCC--TTAATT");
    input[2] = de_strdup("CCCGGAGAG-----ACAAC--CTAATT");
    input[3] = de_strdup("ATTTAAATTTA-----GAG--ACAATC");
    input[4] = de_strdup("CCGGC-------GGTTGGGTTGTCTCT");
    lengths[0] = 20;
    lengths[1] = 20;
    lengths[2] = 20;
    lengths[3] = 20;
    lengths[4] = 20;
    mafTcComparisonOrder_t *obsOrder = NULL;
    int **dummyVizMatrix = NULL;
    obsOrder = getComparisonOrderFromMatrix(input, 5, 27, lengths, dummyVizMatrix);
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
    
    CuAssertTrue(testCase, comparisonOrdersAreEqual(expectedOrder, obsOrder, false));
    // cleanup
    for (int i = 0; i < 5; ++i)
        free(input[i]);
    free(input);
    free(lengths);
    destroyMafTcComparisonOrder(expectedOrder);
    destroyMafTcComparisonOrder(obsOrder);
}
static void test_matrixAlignmentBlockComparisonOrdering_4(CuTest *testCase) {
    // test that with known input that known output is generated.
    char **input = (char**) de_malloc(sizeof(char*) * 5);
    uint32_t *lengths = de_malloc(sizeof(uint32_t) * 5);
    input[0] = de_strdup("AATTG-----TCTCTCC-CC--CT---T-TTT");
    input[1] = de_strdup("AATTGTC-----TCTG-GCC--TTA-AT---T");
    input[2] = de_strdup("CCCGGAGAG-----ACA-AC--CTA--A--TT");
    input[3] = de_strdup("ATTTAAATTTA-----G-AG--ACAA--T--C");
    input[4] = de_strdup("CCGGC-------GGTTG-GGTTGTC-T-C--T");
    lengths[0] = 20;
    lengths[1] = 20;
    lengths[2] = 20;
    lengths[3] = 20;
    lengths[4] = 20;
    mafTcComparisonOrder_t *obsOrder = NULL;
    int **dummyVizMatrix = NULL;
    obsOrder = getComparisonOrderFromMatrix(input, 5, 32, lengths, dummyVizMatrix);
    CuAssertTrue(testCase, obsOrder != NULL);
    mafTcComparisonOrder_t *expectedOrder = newMafTcComparisonOrder();
    mafTcComparisonOrder_t *eo = expectedOrder;
    eo->ref = 4;
    eo->region = newMafTcRegion(20, 21);
    
    eo->next = newMafTcComparisonOrder();
    eo = eo->next;
    eo->ref = 3;
    eo->region = newMafTcRegion(28, 28);
    
    eo->next = newMafTcComparisonOrder();
    eo = eo->next;
    eo->ref = 3;
    eo->region = newMafTcRegion(25, 25);
    
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
    eo->region = newMafTcRegion(26, 26);
    
    eo->next = newMafTcComparisonOrder();
    eo = eo->next;
    eo->ref = 1;
    eo->region = newMafTcRegion(24, 24);
    
    eo->next = newMafTcComparisonOrder();
    eo = eo->next;
    eo->ref = 1;
    eo->region = newMafTcRegion(17, 17);
    
    eo->next = newMafTcComparisonOrder();
    eo = eo->next;
    eo->ref = 1;
    eo->region = newMafTcRegion(5, 6);
    
    eo->next = newMafTcComparisonOrder();
    eo = eo->next;
    eo->ref = 0;
    eo->region = newMafTcRegion(29, 31);

    eo->next = newMafTcComparisonOrder();
    eo = eo->next;
    eo->ref = 0;
    eo->region = newMafTcRegion(27, 27);
    
    eo->next = newMafTcComparisonOrder();
    eo = eo->next;
    eo->ref = 0;
    eo->region = newMafTcRegion(22, 23);
    
    eo->next = newMafTcComparisonOrder();
    eo = eo->next;
    eo->ref = 0;
    eo->region = newMafTcRegion(18, 19);
    
    eo->next = newMafTcComparisonOrder();
    eo = eo->next;
    eo->ref = 0;
    eo->region = newMafTcRegion(10, 16);
    
    eo->next = newMafTcComparisonOrder();
    eo = eo->next;
    eo->ref = 0;
    eo->region = newMafTcRegion(0, 4);
    
    CuAssertTrue(testCase, comparisonOrdersAreEqual(expectedOrder, obsOrder, false));
    // cleanup
    for (int i = 0; i < 5; ++i)
        free(input[i]);
    free(input);
    free(lengths);
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
    mcp.b = -1;
    char *s = de_strdup("ACGT");
    CuAssertTrue(testCase, localSeqCoords(3, s, &mcp, 0) == 3);
    CuAssertTrue(testCase, mcp.a == 3);
    CuAssertTrue(testCase, mcp.b == 3);
    free(s);
    mcp.a = 0;
    mcp.b = -1;
    s = de_strdup("-------ACGT");
    // printf("expecting 3 observing %" PRIi64 "\n", localSeqCoords(10, s, &mcp, 1));
    mcp.a = 0;
    mcp.b = -1;
    CuAssertTrue(testCase, localSeqCoords(10, s, &mcp, 1) == 3);
    CuAssertTrue(testCase, mcp.a == 10);
    CuAssertTrue(testCase, mcp.b == 3);
    free(s);
    mcp.a = 0;
    mcp.b = 0;
    s = de_strdup("A-------ACGT");    
    // printf("expecting 4 observing %" PRIi64 "\n", localSeqCoords(11, s, &mcp, 1));
    mcp.a = 0;
    mcp.b = 0;
    CuAssertTrue(testCase, localSeqCoords(11, s, &mcp, 1) == 4);
    CuAssertTrue(testCase, mcp.a == 11);
    CuAssertTrue(testCase, mcp.b == 4);
    free(s);
    mcp.a = 8;
    mcp.b = 1;
    s = de_strdup("-------ACGT");
    CuAssertTrue(testCase, localSeqCoords(10, s, &mcp, 1) == 3);
    CuAssertTrue(testCase, mcp.a == 10);
    CuAssertTrue(testCase, mcp.b == 3);
    free(s);
    mcp.a = 10;
    mcp.b = 3;
    s = de_strdup("-------ACGT");
    CuAssertTrue(testCase, localSeqCoords(10, s, &mcp, 1) == 3);
    CuAssertTrue(testCase, mcp.a == 10);
    CuAssertTrue(testCase, mcp.b == 3);
    free(s);
    mcp.a = 0;
    mcp.b = -1;
    s = de_strdup("-A-C-G-T");
    CuAssertTrue(testCase, localSeqCoords(7, s, &mcp, 1) == 3);
    CuAssertTrue(testCase, mcp.a == 7);
    CuAssertTrue(testCase, mcp.b == 3);
    free(s);
    mcp.a = 0;
    mcp.b = -1;
    s = de_strdup("AA-A-C-G-T");
    // printf("expecting 5 observing %" PRIi64 "\n", localSeqCoords(9, s, &mcp, 1));
    CuAssertTrue(testCase, localSeqCoords(9, s, &mcp, 1) == 5);
    CuAssertTrue(testCase, mcp.a == 9);
    CuAssertTrue(testCase, mcp.b == 5);
    free(s);
    mcp.a = 0;
    mcp.b = -1;
    s = de_strdup("AA-A-C-G-T");
    CuAssertTrue(testCase, localSeqCoords(3, s, &mcp, 1) == 2);
    CuAssertTrue(testCase, mcp.a == 3);
    CuAssertTrue(testCase, mcp.b == 2);
    free(s);
    mcp.a = 1;
    mcp.b = 1;
    s = de_strdup("AA-A-C-G-T");
    CuAssertTrue(testCase, localSeqCoords(3, s, &mcp, 1) == 2);
    CuAssertTrue(testCase, mcp.a == 3);
    CuAssertTrue(testCase, mcp.b == 2);
    free(s);
    mcp.a = 0;
    mcp.b = -1;
    s = de_strdup("GTTGTCTCTCAATGTG");
    CuAssertTrue(testCase, localSeqCoords(6, s, &mcp, 0) == 6);
    CuAssertTrue(testCase, mcp.a == 6);
    CuAssertTrue(testCase, mcp.b == 6);
    free(s);
    s = de_strdup("GTTGTCTCTCAATGTG");
    CuAssertTrue(testCase, localSeqCoords(15, s, &mcp, 0) == 15);
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
    // printf("19 ?= %" PRIi64 "\n", localSeqCoordsToGlobalPositiveStartCoords(0, 0, 20, '-', 1));
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
    uint32_t start = 0, sourceLength = 100, seqLength = 37, localStart = 1;
    int64_t expectation = 0;
    char strand = '+';
    bool containsGaps = true;
    mcp.a = 0;
    mcp.b = -1;
    CuAssertTrue(testCase, localSeqCoordsToGlobalPositiveStartCoords(localSeqCoords(localStart, input, &mcp, containsGaps), start, sourceLength, strand, seqLength) == expectation);
    CuAssertTrue(testCase, mcp.a == 1);
    CuAssertTrue(testCase, mcp.b == 0);
    CuAssertTrue(testCase, localSeqCoordsToGlobalPositiveStartCoords(localSeqCoords(localStart, input, &mcp, containsGaps), start, sourceLength, strand, seqLength) == expectation);
    CuAssertTrue(testCase, mcp.a == 1);
    CuAssertTrue(testCase, mcp.b == 0);
    localStart = 4;
    expectation = 2;
    CuAssertTrue(testCase, localSeqCoordsToGlobalPositiveStartCoords(localSeqCoords(localStart, input, &mcp, containsGaps), start, sourceLength, strand, seqLength) == expectation);
    CuAssertTrue(testCase, mcp.a == 4);
    CuAssertTrue(testCase, mcp.b == 2);
    localStart = 26;
    expectation = 21;
    CuAssertTrue(testCase, localSeqCoordsToGlobalPositiveStartCoords(localSeqCoords(localStart, input, &mcp, containsGaps), start, sourceLength, strand, seqLength) == expectation);
    CuAssertTrue(testCase, mcp.a == 26);
    CuAssertTrue(testCase, mcp.b == 21);
    localStart = 41;
    expectation = 36;
    CuAssertTrue(testCase, localSeqCoordsToGlobalPositiveStartCoords(localSeqCoords(localStart, input, &mcp, containsGaps), start, sourceLength, strand, seqLength) == expectation);
    CuAssertTrue(testCase, mcp.a == 41);
    CuAssertTrue(testCase, mcp.b == 36);
    // reverse strand.
    mcp.a = 0;
    mcp.b = -1;
    seqLength = 2;
    localStart = 1;
    expectation = 98;
    strand = '-';
    CuAssertTrue(testCase, localSeqCoordsToGlobalPositiveStartCoords(localSeqCoords(localStart, input, &mcp, containsGaps), start, sourceLength, strand, seqLength) == expectation);
    CuAssertTrue(testCase, mcp.a == 1);
    CuAssertTrue(testCase, mcp.b == 0);    
    seqLength = 19;
    localStart = 4;
    expectation = 79;
    CuAssertTrue(testCase, localSeqCoordsToGlobalPositiveStartCoords(localSeqCoords(localStart, input, &mcp, containsGaps), start, sourceLength, strand, seqLength) == expectation);
    CuAssertTrue(testCase, mcp.a == 4);
    CuAssertTrue(testCase, mcp.b == 2);
    seqLength = 16;
    localStart = 26;
    expectation = 63;
    CuAssertTrue(testCase, localSeqCoordsToGlobalPositiveStartCoords(localSeqCoords(localStart, input, &mcp, containsGaps), start, sourceLength, strand, seqLength) == expectation);
    CuAssertTrue(testCase, mcp.a == 26);
    CuAssertTrue(testCase, mcp.b == 21);
    free(input);
    input = de_strdup("GTTGTCTCTCAATGTG");
    mcp.a = 0;
    mcp.b = 0;
    start = 1;
    seqLength = 16;
    localStart = 0;
    expectation = 1;
    strand = '+';
    CuAssertTrue(testCase, localSeqCoordsToGlobalPositiveStartCoords(localSeqCoords(localStart, input, &mcp, 0), start, sourceLength, strand, seqLength) == expectation);
    free(input);
}
static void test_coordinateTransforms_1(CuTest *testCase) {
    char *input = de_strdup("CGTCCGCAGATCGTTAACTTAATTGTTCCGCTTGAAATCCGAAAACT---------GCCAA-CGCTATTGTTGCGACTGATAGTCGTGAATGGCCGTGACCACGCCCCCAATCCC----CTAAGCCCCCCTTT--------------TGGCAAA---------------------------------------CGGC---GCCTATGG--CTGG-----g-aa---------acag-------------------------------------------------------------------------------ga------------------------------------------acaga-----------------------------------------------------------------------------------------aacaggaacagGAATGCAATAAAA---TTGGCGTGACTAACTCAGCACTGGGATGCGAT---------GCGAGCATTG--CAG---------------------ATGAGCTGAG------GATTGGAG--CTTGAAAGTGGAGGAGGA---------T------TGGGGGGG-----GATGAAGGGGT----------------------------------------TC-----TGGGATTGGATGCC------CCAA---TGTGGCAGC---------CACAGAAGGGC--------------GCAAGTCGTGCGTGC---------CTCGGCGAAAC------------GTTGACGC------T-----------------GTCATGCAATCAGCAAATAGGCGACCGCAGCAAAAGTCGCGTAATTAACGC");
    mafCoordinatePair_t mcp;
    uint32_t start = 47275, localStart = 566, sourceLength = 47773, seqLength = 1, localCoord = 234;
    int64_t expectation = 263;
    char strand = '-';
    bool containsGaps = true;
    mcp.a = 0;
    mcp.b = 0;
    CuAssertTrue(testCase, localSeqCoordsToGlobalPositiveStartCoords(localSeqCoords(localStart, input, &mcp, containsGaps), start, sourceLength, strand, seqLength) == expectation);
    CuAssertTrue(testCase, mcp.a == localStart);
    CuAssertTrue(testCase, mcp.b == localCoord);

    localStart = 629;    
    seqLength = 11;
    expectation = 218;
    localCoord = 269;
    // printf("localSeqCoords: %" PRIi64 "\n", localSeqCoords(localStart, input, &mcp, containsGaps));
    // printf("expecting %" PRIi64 " observing %" PRIi64 "\n", expectation, localSeqCoordsToGlobalPositiveStartCoords(localSeqCoords(localStart, input, &mcp, containsGaps), start, sourceLength, strand, seqLength));
    CuAssertTrue(testCase, localSeqCoordsToGlobalPositiveStartCoords(localSeqCoords(localStart, input, &mcp, containsGaps), start, sourceLength, strand, seqLength) == expectation);
    CuAssertTrue(testCase, mcp.a == localStart);
    CuAssertTrue(testCase, mcp.b == localCoord);
    
    free(input);
}
static void test_mafBlockGapSorting_0(CuTest *testCase) {
    // test that the mafBlock_sortBlockByIncreasingGap() function works as expected
    // HEY YOU! qsort is NOT stable on all systems, do NOT write a test case that assumes stability.
    // 
    // build input
    mafBlock_t *input = maf_newMafBlock();
    mafBlock_t *expected = maf_newMafBlock();
    maf_mafBlock_setHeadLine(input, maf_newMafLineFromString("a score=0.0", 1));
    mafLine_t *ml = maf_mafBlock_getHeadLine(input);
    maf_mafLine_setNext(ml, maf_newMafLineFromString("s hg16.chr7    27707221 13 + 158545518 gcagctgaaaaca", 1));
    ml = maf_mafLine_getNext(ml);
    maf_mafLine_setNext(ml, maf_newMafLineFromString("s panTro1.chr6 28869787 13 + 161576975 gcagctgaaaaca", 1));
    ml = maf_mafLine_getNext(ml);
    maf_mafLine_setNext(ml, maf_newMafLineFromString("s baboon         249182 13 +   4622798 gcagctgaaaaca", 1));
    ml = maf_mafLine_getNext(ml);
    maf_mafLine_setNext(ml, maf_newMafLineFromString("s mm4.chr6     53310102 13 + 151104725 ACAGCTGAAAATA", 1));
    maf_mafBlock_setTailLine(input, ml);
    maf_mafBlock_setNumberOfLines(input, 5);
    // build expected
    maf_mafBlock_setHeadLine(expected, maf_newMafLineFromString("a score=0.0", 1));
    ml = maf_mafBlock_getHeadLine(expected);
    maf_mafLine_setNext(ml, maf_newMafLineFromString("s hg16.chr7    27707221 13 + 158545518 gcagctgaaaaca", 1));
    ml = maf_mafLine_getNext(ml);
    maf_mafLine_setNext(ml, maf_newMafLineFromString("s panTro1.chr6 28869787 13 + 161576975 gcagctgaaaaca", 1));
    ml = maf_mafLine_getNext(ml);
    maf_mafLine_setNext(ml, maf_newMafLineFromString("s baboon         249182 13 +   4622798 gcagctgaaaaca", 1));
    ml = maf_mafLine_getNext(ml);
    maf_mafLine_setNext(ml, maf_newMafLineFromString("s mm4.chr6     53310102 13 + 151104725 ACAGCTGAAAATA", 1));
    maf_mafBlock_setTailLine(expected, ml);
    maf_mafBlock_setNumberOfLines(expected, 5);
    // test
    mafBlock_sortBlockByIncreasingGap(input);
    CuAssertTrue(testCase, mafBlocksAreEqual(input, expected, false));
    // cleanup
    maf_destroyMafBlockList(input);
    maf_destroyMafBlockList(expected);
    
    // build input
    input = maf_newMafBlock();
    expected = maf_newMafBlock();
    maf_mafBlock_setHeadLine(input, maf_newMafLineFromString("a score=0.0", 1));
    ml = maf_mafBlock_getHeadLine(input);
    maf_mafLine_setNext(ml, maf_newMafLineFromString("s hg16.chr7    27707221 14 + 158545518 gcagctgaaa--Taca", 1));
    ml = maf_mafLine_getNext(ml);
    maf_mafLine_setNext(ml, maf_newMafLineFromString("s panTro1.chr6 28869787 16 + 161576975 gcagctgaaaTTTaca", 1));
    ml = maf_mafLine_getNext(ml);
    maf_mafLine_setNext(ml, maf_newMafLineFromString("s baboon         249182 15 +   4622798 gcagctgaaa-TTaca", 1));
    ml = maf_mafLine_getNext(ml);
    maf_mafLine_setNext(ml, maf_newMafLineFromString("s mm4.chr6     53310102 13 + 151104725 ACAGCTGAAA---ATA", 1));
    maf_mafBlock_setTailLine(input, ml);
    maf_mafBlock_setNumberOfLines(input, 5);
    // build expected
    maf_mafBlock_setHeadLine(expected, maf_newMafLineFromString("a score=0.0", 1));
    ml = maf_mafBlock_getHeadLine(expected);
    maf_mafLine_setNext(ml, maf_newMafLineFromString("s panTro1.chr6 28869787 16 + 161576975 gcagctgaaaTTTaca", 1));
    ml = maf_mafLine_getNext(ml);
    maf_mafLine_setNext(ml, maf_newMafLineFromString("s baboon         249182 15 +   4622798 gcagctgaaa-TTaca", 1));
    ml = maf_mafLine_getNext(ml);
    maf_mafLine_setNext(ml, maf_newMafLineFromString("s hg16.chr7    27707221 14 + 158545518 gcagctgaaa--Taca", 1)); 
    ml = maf_mafLine_getNext(ml);
    maf_mafLine_setNext(ml, maf_newMafLineFromString("s mm4.chr6     53310102 13 + 151104725 ACAGCTGAAA---ATA", 1));
    maf_mafBlock_setTailLine(expected, ml);
    maf_mafBlock_setNumberOfLines(expected, 5);
    // test
    mafBlock_sortBlockByIncreasingGap(input);
    CuAssertTrue(testCase, mafBlocksAreEqual(input, expected, false));
    // cleanup
    maf_destroyMafBlockList(input);
    maf_destroyMafBlockList(expected);
}
CuSuite* mafTransitiveClosure_TestSuite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, test_reverseComplement);
    SUITE_ADD_TEST(suite, test_rowAlignmentBlockComparisonOrdering_0);
    SUITE_ADD_TEST(suite, test_rowAlignmentBlockComparisonOrdering_1);
    SUITE_ADD_TEST(suite, test_rowAlignmentBlockComparisonOrdering_2);
    SUITE_ADD_TEST(suite, test_rowAlignmentBlockComparisonOrdering_3);
    SUITE_ADD_TEST(suite, test_matrixAlignmentBlockComparisonOrdering_0);
    SUITE_ADD_TEST(suite, test_matrixAlignmentBlockComparisonOrdering_1);
    SUITE_ADD_TEST(suite, test_matrixAlignmentBlockComparisonOrdering_2);
    SUITE_ADD_TEST(suite, test_matrixAlignmentBlockComparisonOrdering_3);
    SUITE_ADD_TEST(suite, test_matrixAlignmentBlockComparisonOrdering_4);
    SUITE_ADD_TEST(suite, test_addSequenceValuesToMtcSeq_0);
    SUITE_ADD_TEST(suite, test_localSeqCoords_0);
    SUITE_ADD_TEST(suite, test_localSeqCoordsToGlobalPositiveCoords_0);
    SUITE_ADD_TEST(suite, test_localSeqCoordsToGlobalPositiveStartCoords_0);
    SUITE_ADD_TEST(suite, test_coordinateTransforms_0);
    SUITE_ADD_TEST(suite, test_coordinateTransforms_1);
    SUITE_ADD_TEST(suite, test_mafBlockGapSorting_0);
    return suite;
}
