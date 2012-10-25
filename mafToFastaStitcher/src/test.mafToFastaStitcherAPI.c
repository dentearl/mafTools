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
#include <string.h>
#include <time.h>
#include "CuTest.h"
#include "common.h"
#include "sonLib.h"
#include "mafToFastaStitcher.h"
#include "mafToFastaStitcherAPI.h"
#include "test.mafToFastaStitcherAPI.h"

static void printBlockHash(stHash *hash, const char *title) {
    stHashIterator *hit = stHash_getIterator(hash);
    char *key = NULL;
    row_t *r = NULL;
    printf("%s:\n", title);
    while ((key = stHash_getNext(hit)) != NULL) {
        r = stHash_search(hash, key);
        printf("%20s %6"PRIu32" %6"PRIu32" %c %9"PRIu32" %s\n", r->name ,r->start, r->length, 
               r->strand, r->sourceLength, r->sequence);
    }
    stHash_destructIterator(hit);
}
static void printList(stList *list, const char *title) {
    char *key = NULL;
    printf("%s:\n", title);
    for (int32_t i = 0; i < stList_length(list); ++i) {
        key = (char*)stList_get(list, i);
        printf("%s ", key);
    }
    printf("\n");
}
static bool rowsAreEqual(row_t *a, row_t *b) {
    if (a->name != NULL && b->name != NULL) {
        if (strcmp(a->name, b->name) != 0) {
            fprintf(stderr, "rows differ: names a:%s, b:%s\n", a->name, b->name);
            return false;
        }
    }
    if (a->prevName != NULL && b->prevName != NULL) {
        if (strcmp(a->prevName, b->prevName) != 0) {
            fprintf(stderr, "%s rows differ: prevNames a:%s, b:%s\n", a->name, a->prevName, b->prevName);
            return false;
        }
    }
    if (a->sequence != NULL && b->sequence != NULL) {
        if (strcmp(a->sequence, b->sequence) != 0) {
            fprintf(stderr, "%s rows differ: sequences:\n    %s\n    %s\n", a->name, a->sequence, b->sequence);
            return false;
        }
    }
    if (a->multipleNames != b->multipleNames) {
        fprintf(stderr, "%s rows differ: multipleNames: %d %d\n", a->name, a->multipleNames, b->multipleNames);
        return false;
    }
    if (a->start != b->start) {
        fprintf(stderr, "%s rows differ: start: %"PRIu32" %"PRIu32"\n", a->name, a->start, b->start);
        return false;
    }
    if (a->length != b->length) {
        fprintf(stderr, "%s rows differ: length: %"PRIu32" %"PRIu32"\n"
                "    %s\n    %s\n", a->name, a->length, b->length, a->sequence, b->sequence);
        return false;
    }
    if (a->prevRightPos != b->prevRightPos) {
        fprintf(stderr, "%s rows differ: prevRightPos: %"PRIu32" %"PRIu32"\n", 
                a->name, a->prevRightPos, b->prevRightPos);
        return false;
    }
    if (a->strand != b->strand) {
        fprintf(stderr, "%s rows differ: strand: %c %c\n", a->name, a->strand, b->strand);
        return false;
    }
    if (a->prevStrand != b->prevStrand) {
        fprintf(stderr, "%s rows differ: prevStrand: %c %c\n", a->name, a->prevStrand, b->prevStrand);
        return false;
    }
    if (a->sourceLength != b->sourceLength) {
        fprintf(stderr, "%s rows differ: sourceLength: %"PRIu32" %"PRIu32"\n", a->name, 
                a->sourceLength, b->sourceLength);
        return false;
    }
    if (a->memLength != b->memLength) {
        fprintf(stderr, "%s rows differ: memLength: %"PRIu32" %"PRIu32"\n", a->name, 
                a->memLength, b->memLength);
        return false;
    }
    if (a->index != b->index) {
        fprintf(stderr, "%s rows differ: index: %"PRIu32" %"PRIu32"\n", a->name, 
                a->index, b->index);
        return false;
    }
    return true;
}
static bool hashesAreEqual(stHash *observedHash, stHash *expectedHash) {
    stHashIterator *hit = stHash_getIterator(observedHash);
    char *key;
    while ((key = stHash_getNext(hit)) != NULL) {
        if (stHash_search(expectedHash, key) == NULL) {
            printBlockHash(observedHash, "observed");
            printBlockHash(expectedHash, "expected");
            return false;
        }
        if (!rowsAreEqual(stHash_search(observedHash, key), stHash_search(expectedHash, key))) {
            printBlockHash(observedHash, "observed");
            printBlockHash(expectedHash, "expected");
            return false;
        }
    }
    stHash_destructIterator(hit);
    hit = stHash_getIterator(expectedHash);
    while ((key = stHash_getNext(hit)) != NULL) {
        if (stHash_search(observedHash, key) == NULL) {
            printBlockHash(observedHash, "observed");
            printBlockHash(expectedHash, "expected");
            return false;
        }
        if (!rowsAreEqual(stHash_search(observedHash, key), stHash_search(expectedHash, key))) {
            printBlockHash(observedHash, "observed");
            printBlockHash(expectedHash, "expected");
            return false;
        }
    }
    stHash_destructIterator(hit);
    return true;
}
static bool listsAreEqual(stList *observedList, stList *expectedList) {
    if (stList_length(observedList) != stList_length(expectedList)) {
        fprintf(stderr, "stList lengths are not equal: %"PRIu32" %"PRIu32"\n", 
                stList_length(observedList), stList_length(expectedList));
        printList(observedList, "observed");
        printList(expectedList, "expected");
        return false;
    }
    for (int32_t i = 0; i < stList_length(observedList); ++i) {
        if (strcmp(stList_get(observedList, i), stList_get(expectedList, i)) != 0) {
            fprintf(stderr, "stList elements are not equal at index %"PRIu32": %s %s\n", 
                    i, (char *)stList_get(observedList, i), (char *)stList_get(expectedList, i));
            return false;
        }
    }
    return true;
}
static mtfseq_t* newMtfseqFromString(char *s) {
    mtfseq_t *mtfs = newMtfseq(strlen(s) + 1);
    seq_copyIn(mtfs, s);
    return mtfs;
}
static stHash *createBlockHashFromString(char *input, stList *orderList) {
    mafBlock_t *ibhead = maf_newMafBlockListFromString(input, 3);
    stHash *hash = mafBlockToBlockHash(ibhead, orderList);
    maf_destroyMafBlockList(ibhead);
    return hash;
}
static stHash *createSeqHashFromString(char *name, char *input) {
    mtfseq_t *mtfs = newMtfseq(strlen(input));
    stHash *hash = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, free, destroyMtfseq);
    seq_copyIn(mtfs, input);
    stHash_insert(hash, stString_copy(name), mtfs);
    return hash;
}
static void test_newBlockHashFromBlock_0(CuTest *testCase) {
    stList *orderList = stList_construct3(0, free);
    stHash *observedHash = createBlockHashFromString("a score=0 test=0\n"
                                                     "s reference.chr0 0 13 + 158545518 gcagctgaaaaca\n"
                                                     "s name.chr1      0 10 +       100 ATGT---ATGCCG\n"
                                                     "s name2.chr1     0 10 +       100 ATGT---ATGCCG\n"
                                                     "s name3.chr1     5 13 -       100 ATGTgggATGCCG\n",
                                                     orderList);
    CuAssertTrue(testCase, stHash_size(observedHash) == 4);
    row_t *key = NULL;
    key = stHash_search(observedHash, "reference");
    CuAssertTrue(testCase, key != NULL);
    CuAssertTrue(testCase, strcmp(key->name, "reference.chr0") == 0);
    CuAssertTrue(testCase, strcmp(key->prevName, "reference.chr0") == 0);
    CuAssertTrue(testCase, key->multipleNames == false);
    CuAssertTrue(testCase, key->start == 0);
    CuAssertTrue(testCase, key->length == 13);
    CuAssertTrue(testCase, key->sourceLength == 158545518);
    CuAssertTrue(testCase, key->prevRightPos == 12);
    CuAssertTrue(testCase, key->strand == '+');
    CuAssertTrue(testCase, key->prevStrand == '+');
    CuAssertTrue(testCase, strcmp(key->sequence, "gcagctgaaaaca") == 0);
    // row 2
    key = stHash_search(observedHash, "name");
    CuAssertTrue(testCase, key != NULL);
    CuAssertTrue(testCase, strcmp(key->name, "name.chr1") == 0);
    CuAssertTrue(testCase, strcmp(key->prevName, "name.chr1") == 0);
    CuAssertTrue(testCase, key->multipleNames == false);
    CuAssertTrue(testCase, key->start == 0);
    CuAssertTrue(testCase, key->length == 10);
    CuAssertTrue(testCase, key->sourceLength == 100);
    CuAssertTrue(testCase, key->prevRightPos == 9);
    CuAssertTrue(testCase, key->strand == '+');
    CuAssertTrue(testCase, key->prevStrand == '+');
    CuAssertTrue(testCase, strcmp(key->sequence, "ATGT---ATGCCG") == 0);
    // row 3
    key = stHash_search(observedHash, "name2");
    CuAssertTrue(testCase, key != NULL);
    CuAssertTrue(testCase, strcmp(key->name, "name2.chr1") == 0);
    CuAssertTrue(testCase, strcmp(key->prevName, "name2.chr1") == 0);
    CuAssertTrue(testCase, key->multipleNames == false);
    CuAssertTrue(testCase, key->start == 0);
    CuAssertTrue(testCase, key->length == 10);
    CuAssertTrue(testCase, key->sourceLength == 100);
    CuAssertTrue(testCase, key->prevRightPos == 9);
    CuAssertTrue(testCase, key->strand == '+');
    CuAssertTrue(testCase, key->prevStrand == '+');
    CuAssertTrue(testCase, strcmp(key->sequence, "ATGT---ATGCCG") == 0);
    // row 4
    key = stHash_search(observedHash, "name3");
    CuAssertTrue(testCase, key != NULL);
    CuAssertTrue(testCase, strcmp(key->name, "name3.chr1") == 0);
    CuAssertTrue(testCase, strcmp(key->prevName, "name3.chr1") == 0);
    CuAssertTrue(testCase, key->multipleNames == false);
    CuAssertTrue(testCase, key->start == 5);
    CuAssertTrue(testCase, key->length == 13);
    CuAssertTrue(testCase, key->sourceLength == 100);
    CuAssertTrue(testCase, key->prevRightPos == 17);
    CuAssertTrue(testCase, key->strand == '-');
    CuAssertTrue(testCase, key->prevStrand == '-');
    CuAssertTrue(testCase, strcmp(key->sequence, "ATGTgggATGCCG") == 0);
    stList_destruct(orderList);
    stHash_destruct(observedHash);
}
static void test_addMafLineToRow_0(CuTest *testCase) {
    row_t *obs = newRow(20);
    obs->name = stString_copy("seq1.chr0");
    obs->prevName = stString_copy("seq1.chr0");
    obs->multipleNames = false;
    obs->start = 3;
    obs->length = 10;
    obs->prevRightPos = 20;
    obs->strand = '+';
    obs->prevStrand = '+';
    obs->sourceLength = 100;
    row_copyIn(obs, "acgtacgtac");
    mafLine_t *ml = maf_newMafLineFromString("s seq1.chr0 13 5 + 100 ACGTA", 10);
    addMafLineToRow(obs, ml);
    row_t *exp = newRow(20);
    exp->name = stString_copy("seq1.chr0");
    exp->prevName = stString_copy("seq1.chr0");
    exp->multipleNames = false;
    exp->start = 3;
    exp->length = 15;
    exp->prevRightPos = 17;
    exp->strand = '+';
    exp->prevStrand = '+';
    exp->sourceLength = 100;
    row_copyIn(exp, "acgtacgtacACGTA");
    CuAssertTrue(testCase, rowsAreEqual(obs, exp));
    destroyRow(obs);
    destroyRow(exp);
    maf_destroyMafLineList(ml);
}
static void test_addMafLineToRow_1(CuTest *testCase) {
    row_t *obs = newRow(20);
    obs->name = stString_copy("seq1.chr0");
    obs->prevName = stString_copy("seq1.amazing");
    obs->multipleNames = true;
    obs->start = 3;
    obs->length = 10;
    obs->prevRightPos = 20;
    obs->strand = '+';
    obs->prevStrand = '+';
    obs->sourceLength = 100;
    row_copyIn(obs, "acgtacgtac");
    mafLine_t *ml = maf_newMafLineFromString("s seq1.chr_bleh 13 5 + 100 ACGTA", 10);
    addMafLineToRow(obs, ml);
    row_t *exp = newRow(20);
    exp->name = stString_copy("seq1.chr0");
    exp->prevName = stString_copy("seq1.chr_bleh");
    exp->multipleNames = true;
    exp->start = 3;
    exp->length = 15;
    exp->prevRightPos = 17;
    exp->strand = '+';
    exp->prevStrand = '+';
    exp->sourceLength = 15;
    row_copyIn(exp, "acgtacgtacACGTA");
    CuAssertTrue(testCase, rowsAreEqual(obs, exp));
    destroyRow(obs);
    destroyRow(exp);
    maf_destroyMafLineList(ml);
}

static void test_penalize_0(CuTest *testCase) {
    stList *observedList = stList_construct3(0, free);
    stList *expectedList = stList_construct3(0, free);
    stHash *observedHash = createBlockHashFromString("a score=0\n"
                                                     "s reference.chr0 0 13 + 158545518 gcagctgaaaaca\n"
                                                     "s name.chr1      0 10 +       100 ATGT---ATGCCG\n"
                                                     "s name2.chr1     0 10 +       100 ATGT---ATGCCG\n",
                                                     observedList);
    stHash *expectedHash = createBlockHashFromString("a score=0\n"
                                                     "s reference.chr0 0 13 + 158545518 gcagctgaaaaca-----\n"
                                                     "s name.chr1      0 15 +       100 ATGT---ATGCCGNNNNN\n"
                                                     "s name2.chr1     0 10 +       100 ATGT---ATGCCG-----\n",
                                                     expectedList);
    penalize(observedHash, "name.chr1", 5);
    CuAssertTrue(testCase, hashesAreEqual(observedHash, expectedHash));
    CuAssertTrue(testCase, listsAreEqual(observedList, expectedList));
    // clean up
    stHash_destruct(observedHash);
    stHash_destruct(expectedHash);
    stList_destruct(observedList);
    stList_destruct(expectedList);
}
static void test_interstitial_0(CuTest *testCase) {
    stList *observedList = stList_construct3(0, free);
    stList *expectedList = stList_construct3(0, free);
    stHash *observedHash = createBlockHashFromString("a score=0\n"
                                                     "s reference.chr0 0 13 + 158545518 gcagctgaaaaca\n"
                                                     "s name.chr1      0 10 +       100 ATGT---ATGCCG\n"
                                                     "s name2.chr1     0 10 +       100 ATGT---ATGCCG\n",
                                                     observedList);
    stHash *expectedHash = createBlockHashFromString("a score=0\n"
                                                     "s reference.chr0 0 13 + 158545518 gcagctgaaaaca-----\n"
                                                     "s name.chr1      0 15 +       100 ATGT---ATGCCGaaaTa\n"
                                                     "s name2.chr1     0 10 +       100 ATGT---ATGCCG-----\n",
                                                     expectedList);
    stHash *seqHash = createSeqHashFromString("name.chr1", "ATGTATGCCGaaaTaTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
                                              "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT");
    interstitialInsert(observedHash, seqHash, "name.chr1", 10, '+', 5);
    CuAssertTrue(testCase, hashesAreEqual(observedHash, expectedHash));
    CuAssertTrue(testCase, listsAreEqual(observedList, expectedList));
    // clean up
    stHash_destruct(observedHash);
    stHash_destruct(expectedHash);
    stHash_destruct(seqHash);
    stList_destruct(observedList);
    stList_destruct(expectedList);
    observedList = stList_construct3(0, free);
    expectedList = stList_construct3(0, free);
    observedHash = createBlockHashFromString("a score=0\n"
                                             "s reference.chr0 0 13 + 158545518 gcagctgaaaaca\n"
                                             "s name.chr1      0 10 -       100 ATGT---ATGCCG\n"
                                             "s name2.chr1     0 10 +       100 ATGT---ATGCCG\n",
                                             observedList);
    expectedHash = createBlockHashFromString("a score=0\n"
                                             "s reference.chr0 0 13 + 158545518 gcagctgaaaaca-----\n"
                                             "s name.chr1      0 15 -       100 ATGT---ATGCCGaaaTa\n"
                                             "s name2.chr1     0 10 +       100 ATGT---ATGCCG-----\n",
                                             expectedList);
    seqHash = createSeqHashFromString("name.chr1", "ggggggggggggTTgggggggggggggggggggggggggggggggggggg" // 50
                                      "gggaagggGgggCgggTgggAgggcgggtgggagg" // 35
                                      "tAtttCGGCATACAT");
    interstitialInsert(observedHash, seqHash, "name.chr1", 10, '-', 5);
    CuAssertTrue(testCase, hashesAreEqual(observedHash, expectedHash));
    CuAssertTrue(testCase, listsAreEqual(observedList, expectedList));
    // clean up
    stHash_destruct(observedHash);
    stHash_destruct(expectedHash);
    stHash_destruct(seqHash);
    stList_destruct(observedList);
    stList_destruct(expectedList);
}
static void test_addBlockToHash_0(CuTest *testCase) {
    // basic concatenation
    options_t *options = options_construct();
    options->breakpointPenalty = 10;
    options->interstitialSequence = 5;
    stList *observedList = stList_construct3(0, free);
    stList *expectedList = stList_construct3(0, free);
    stHash *observedHash = createBlockHashFromString("a score=0\n"
                                                     "s reference.chr0 0 13 + 158545518 gcagctgaaaaca\n"
                                                     "s name.chr1      0 10 +       100 ATGT---ATGCCG\n"
                                                     "s name2.chr1     0 10 +       100 ATGT---ATGCCG\n",
                                                     observedList);
    mafBlock_t *mb = maf_newMafBlockListFromString("a score=0 test\n"
                                                   "s reference.chr0 13  5 + 158545518 ACGTA\n"
                                                   "s name.chr1      10  5 +       100 acgtc\n"
                                                   "s name2.chr1     10  5 +       100 ATGTg\n"
                                                   , 3);
    stHash *expectedHash = createBlockHashFromString("a score=0\n"
                                                     "s reference.chr0 0 18 + 158545518 gcagctgaaaacaACGTA\n"
                                                     "s name.chr1      0 15 +       100 ATGT---ATGCCGacgtc\n"
                                                     "s name2.chr1     0 15 +       100 ATGT---ATGCCGATGTg\n",
                                                     expectedList);
    stHash *seqHash = createSeqHashFromString("name.chr1", "ATGTATGCCGacgtc"
                                              "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG"
                                              "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG");
    mtfseq_t *mtfs = newMtfseqFromString("gcagctgaaaacaACGTA"
                                         "tttttttttttttttttttttttttttttttt"
                                         "tttttttttttttttttttttttttttttttttttttttttttttttttt");
    stHash_insert(seqHash, stString_copy("reference.chr0"), mtfs);
    mtfs = newMtfseqFromString("ATGTATGCCGATGTg"
                               "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"
                               "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC");
    stHash_insert(seqHash, stString_copy("name2.chr1"), mtfs);
    addMafBlockToRowHash(observedHash, seqHash, observedList, mb, options);
    CuAssertTrue(testCase, hashesAreEqual(observedHash, expectedHash));
    CuAssertTrue(testCase, listsAreEqual(observedList, expectedList));
    // clean up
    stHash_destruct(observedHash);
    stHash_destruct(expectedHash);
    stHash_destruct(seqHash);
    stList_destruct(observedList);
    stList_destruct(expectedList);
    maf_destroyMafBlockList(mb);
    destroyOptions(options);
}
static void test_addBlockToHash_1(CuTest *testCase) {
    // concatenation with 2 bases of interstitial
    options_t *options = options_construct();
    options->breakpointPenalty = 10;
    options->interstitialSequence = 5;
    stList *observedList = stList_construct3(0, free);
    stList *expectedList = stList_construct3(0, free);
    stHash *observedHash = createBlockHashFromString("a score=0\n"
                                                     "s reference.chr0 0 13 + 158545518 gcagctgaaaaca\n"
                                                     "s name.chr1      0 10 +       100 ATGT---ATGCCG\n"
                                                     "s name2.chr1     0 10 +       100 ATGT---ATGCCG\n",
                                                     observedList);
    mafBlock_t *mb = maf_newMafBlockListFromString("a score=0 test\n"
                                                   "s reference.chr0 13  5 + 158545518 ACGTA\n"
                                                   "s name.chr1      12  5 +       100 gtcGG\n"
                                                   "s name2.chr1     10  5 +       100 ATGTg\n"
                                                   , 3);
    stHash *expectedHash = createBlockHashFromString("a score=0\n"
                                                     "s reference.chr0 0 18 + 158545518 gcagctgaaaaca--ACGTA\n"
                                                     "s name.chr1      0 17 +       100 ATGT---ATGCCGacgtcGG\n"
                                                     "s name2.chr1     0 15 +       100 ATGT---ATGCCG--ATGTg\n",
                                                     expectedList);
    stHash *seqHash = createSeqHashFromString("name.chr1", "ATGTATGCCGacgtc"
                                              "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG"
                                              "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG");
    mtfseq_t *mtfs = newMtfseqFromString("gcagctgaaaacaACGTA"
                                         "tttttttttttttttttttttttttttttttt"
                                         "tttttttttttttttttttttttttttttttttttttttttttttttttt");
    stHash_insert(seqHash, stString_copy("reference.chr0"), mtfs);
    mtfs = newMtfseqFromString("ATGTATGCCGATGTg"
                               "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"
                               "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC");
    stHash_insert(seqHash, stString_copy("name2.chr1"), mtfs);
    addMafBlockToRowHash(observedHash, seqHash, observedList, mb, options);
    CuAssertTrue(testCase, hashesAreEqual(observedHash, expectedHash));
    CuAssertTrue(testCase, listsAreEqual(observedList, expectedList));
    // clean up
    stHash_destruct(observedHash);
    stHash_destruct(expectedHash);
    stHash_destruct(seqHash);
    stList_destruct(observedList);
    stList_destruct(expectedList);
    maf_destroyMafBlockList(mb);
    destroyOptions(options);
}
static void test_addBlockToHash_2(CuTest *testCase) {
    // concatenation with 2 bases of interstitial AND a previously unobserved sequence
    options_t *options = options_construct();
    options->breakpointPenalty = 10;
    options->interstitialSequence = 5;
    stList *observedList = stList_construct3(0, free);
    stList *expectedList = stList_construct3(0, free);
    stHash *observedHash = createBlockHashFromString("a score=0\n"
                                                     "s reference.chr0 0 13 + 158545518 gcagctgaaaaca\n"
                                                     "s name.chr1      0 10 +       100 ATGT---ATGCCG\n"
                                                     "s name2.chr1     0 10 +       100 ATGT---ATGCCG\n",
                                                     observedList);
    mafBlock_t *mb = maf_newMafBlockListFromString("a score=0 test\n"
                                                   "s reference.chr0 13  5 + 158545518 ACGTA\n"
                                                   "s name.chr1      12  5 +       100 gTcGG\n"
                                                   "s name2.chr1     10  5 +       100 ATGTg\n"
                                                   "s name3.chr@      0  5 +        20 aaccg\n"
                                                   , 3);
    stHash *expectedHash = createBlockHashFromString("a score=0\n"
                                                     "s reference.chr0 0 18 + 158545518 gcagctgaaaaca--ACGTA\n"
                                                     "s name.chr1      0 17 +       100 ATGT---ATGCCGacgTcGG\n"
                                                     "s name2.chr1     0 15 +       100 ATGT---ATGCCG--ATGTg\n"
                                                     "s name3.chr@     0  5 +        20 ---------------aaccg\n",
                                                     expectedList
                                                     );
    stHash *seqHash = createSeqHashFromString("name.chr1", "ATGTATGCCGacgTc"
                                              "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG"
                                              "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG");
    mtfseq_t *mtfs = newMtfseqFromString("gcagctgaaaacaACGTA"
                                         "tttttttttttttttttttttttttttttttt"
                                         "tttttttttttttttttttttttttttttttttttttttttttttttttt");
    stHash_insert(seqHash, stString_copy("reference.chr0"), mtfs);
    mtfs = newMtfseqFromString("ATGTATGCCGATGTg"
                               "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"
                               "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC");
    stHash_insert(seqHash, stString_copy("name2.chr1"), mtfs);
    mtfs = newMtfseqFromString("aaccgTTTTTTTTTTTTTTT");
    stHash_insert(seqHash, stString_copy("name3.chr@"), mtfs);
    addMafBlockToRowHash(observedHash, seqHash, observedList, mb, options);
    CuAssertTrue(testCase, hashesAreEqual(observedHash, expectedHash));
    CuAssertTrue(testCase, listsAreEqual(observedList, expectedList));
    // clean up
    stHash_destruct(observedHash);
    stHash_destruct(expectedHash);
    stHash_destruct(seqHash);
    stList_destruct(observedList);
    stList_destruct(expectedList);
    maf_destroyMafBlockList(mb);
    destroyOptions(options);
}
static void test_addBlockToHash_3(CuTest *testCase) {
    // concatenation with 2 bases of interstitial and a sequence length breakpoint
    options_t *options = options_construct();
    options->breakpointPenalty = 10;
    options->interstitialSequence = 5;
    stList *observedList = stList_construct3(0, free);
    stList *expectedList = stList_construct3(0, free);
    stHash *observedHash = createBlockHashFromString("a score=0\n"
                                                     "s reference.chr0 0 13 + 158545518 gcagctgaaaaca\n"
                                                     "s name.chr1      0 10 +       100 ATGT---ATGCCG\n"
                                                     "s name2.chr1     0 10 +       100 ATGT---ATGCCG\n"
                                                     "s name3.chr1     0 13 +       100 GCAGCTGAAAACA\n",
                                                     observedList
                                                     );
    mafBlock_t *mb = maf_newMafBlockListFromString("a score=0 test\n"
                                                   "s reference.chr0 13  5 + 158545518 ACGTA\n"
                                                   "s name.chr1      12  5 +       100 gtcGG\n"
                                                   "s name2.chr1     10  5 +       100 ATGTg\n"
                                                   "s name3.chr1     50  5 +       100 CCCCC\n"
                                                   , 3);
    stHash *expectedHash = NULL;
    expectedHash = createBlockHashFromString("a score=0\n"
                                             "s reference.chr0 0 18 + 158545518 gcagctgaaaaca------------ACGTA\n"
                                             "s name.chr1      0 17 +       100 ATGT---ATGCCGac----------gtcGG\n"
                                             "s name2.chr1     0 15 +       100 ATGT---ATGCCG------------ATGTg\n"
                                             "s name3          0 28 +        28 GCAGCTGAAAACA--NNNNNNNNNNCCCCC\n",
                                             expectedList
                                             );
    row_t *r = stHash_search(expectedHash, "name3");
    r->prevRightPos = 54;
    free(r->prevName);
    r->prevName = stString_copy("name3.chr1");
    r->multipleNames = true;
    stHash *seqHash = createSeqHashFromString("name.chr1", "ATGTATGCCGacgtc"
                                              "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG"
                                              "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG");
    mtfseq_t *mtfs = newMtfseqFromString("gcagctgaaaacaACGTA"
                                         "tttttttttttttttttttttttttttttttt"
                                         "tttttttttttttttttttttttttttttttttttttttttttttttttt");
    stHash_insert(seqHash, stString_copy("reference.chr0"), mtfs);
    mtfs = newMtfseqFromString("ATGTATGCCGATGTg"
                               "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"
                               "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC");
    stHash_insert(seqHash, stString_copy("name2.chr1"), mtfs);
    mtfs = newMtfseqFromString("GCAGCTGAAAACAggggggggggggggggggggggggggggggggggggg"
                               "CCCCCaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"
                               );
    stHash_insert(seqHash, stString_copy("name3.chr1"), mtfs);
    addMafBlockToRowHash(observedHash, seqHash, observedList, mb, options);
    CuAssertTrue(testCase, hashesAreEqual(observedHash, expectedHash));
    CuAssertTrue(testCase, listsAreEqual(observedList, expectedList));
    // clean up
    stHash_destruct(observedHash);
    stHash_destruct(expectedHash);
    stHash_destruct(seqHash);
    stList_destruct(observedList);
    stList_destruct(expectedList);
    maf_destroyMafBlockList(mb);
    destroyOptions(options);
}
static void test_addBlockToHash_4(CuTest *testCase) {
    // concatenation with sequnece breakpoint due to *strand* alone
    // note that name3 is well within the interstitial boundary, the two blocks
    // essentially looking like >>>>>>>>>>>>> <<<<< (strand diffs)
    options_t *options = options_construct();
    options->breakpointPenalty = 10;
    options->interstitialSequence = 5;
    stList *observedList = stList_construct3(0, free);
    stList *expectedList = stList_construct3(0, free);
    stHash *observedHash = createBlockHashFromString("a score=0\n"
                                                     "s reference.chr0 0 13 + 158545518 gcagctgaaaaca\n"
                                                     "s name.chr1      0 10 +       100 ATGT---ATGCCG\n"
                                                     "s name2.chr1     0 10 +       100 ATGT---ATGCCG\n"
                                                     "s name3.chr1     0 13 +       100 GCAGCTGAAAACA\n",
                                                     observedList
                                                     );
    mafBlock_t *mb = maf_newMafBlockListFromString("a score=0 test\n"
                                                   "s reference.chr0 13  5 + 158545518 ACGTA\n"
                                                   "s name.chr1      12  5 +       100 gtcGG\n"
                                                   "s name2.chr1     10  5 +       100 ATGTg\n"
                                                   "s name3.chr1     82  5 -       100 GGGGG\n"
                                                   , 3);
    stHash *expectedHash = NULL;
    expectedHash = createBlockHashFromString("a score=0\n"
                                             "s reference.chr0 0 18 + 158545518 gcagctgaaaaca------------ACGTA\n"
                                             "s name.chr1      0 17 +       100 ATGT---ATGCCGac----------gtcGG\n"
                                             "s name2.chr1     0 15 +       100 ATGT---ATGCCG------------ATGTg\n"
                                             "s name3.chr1     0 28 +       100 GCAGCTGAAAACA--NNNNNNNNNNGGGGG\n",
                                             expectedList
                                             );
    row_t *r = stHash_search(expectedHash, "name3");
    r->prevRightPos = 86;
    r->strand = '*';
    r->prevStrand = '-';
    stHash *seqHash = createSeqHashFromString("name.chr1", "ATGTATGCCGacgtc"
                                              "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG"
                                              "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG");
    mtfseq_t *mtfs = newMtfseqFromString("gcagctgaaaacaACGTA"
                                         "tttttttttttttttttttttttttttttttt"
                                         "tttttttttttttttttttttttttttttttttttttttttttttttttt");
    stHash_insert(seqHash, stString_copy("reference.chr0"), mtfs);
    mtfs = newMtfseqFromString("ATGTATGCCGATGTg"
                               "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"
                               "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC");
    stHash_insert(seqHash, stString_copy("name2.chr1"), mtfs);
    mtfs = newMtfseqFromString("GCAGCTGAAAACACCCCCgggggggggggggggggggggggggggggggg"
                               "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"
                               );
    stHash_insert(seqHash, stString_copy("name3.chr1"), mtfs);
    addMafBlockToRowHash(observedHash, seqHash, observedList, mb, options);
    CuAssertTrue(testCase, hashesAreEqual(observedHash, expectedHash));
    CuAssertTrue(testCase, listsAreEqual(observedList, expectedList));
    // clean up
    stHash_destruct(observedHash);
    stHash_destruct(expectedHash);
    stHash_destruct(seqHash);
    stList_destruct(observedList);
    stList_destruct(expectedList);
    maf_destroyMafBlockList(mb);
    destroyOptions(options);
}
static void test_addBlockToHash_5(CuTest *testCase) {
    // concatenation with sequnece breakpoint due to distance and strand
    options_t *options = options_construct();
    options->breakpointPenalty = 10;
    options->interstitialSequence = 5;
    stList *observedList = stList_construct3(0, free);
    stList *expectedList = stList_construct3(0, free);
    stHash *observedHash = createBlockHashFromString("a score=0\n"
                                                     "s reference.chr0 0 13 + 158545518 gcagctgaaaaca\n"
                                                     "s name.chr1      0 10 +       100 ATGT---ATGCCG\n"
                                                     "s name2.chr1     0 10 +       100 ATGT---ATGCCG\n"
                                                     "s name3.chr1     0 13 +       100 GCAGCTGAAAACA\n",
                                                     observedList
                                                     );
    mafBlock_t *mb = maf_newMafBlockListFromString("a score=0 test\n"
                                                   "s reference.chr0 13  5 + 158545518 ACGTA\n"
                                                   "s name.chr1      12  5 +       100 gtcGG\n"
                                                   "s name2.chr1     10  5 +       100 ATGTg\n"
                                                   "s name3.chr1     50  5 -       100 GGGGG\n"
                                                   , 3);
    stHash *expectedHash = NULL;
    expectedHash = createBlockHashFromString("a score=0\n"
                                             "s reference.chr0 0 18 + 158545518 gcagctgaaaaca------------ACGTA\n"
                                             "s name.chr1      0 17 +       100 ATGT---ATGCCGac----------gtcGG\n"
                                             "s name2.chr1     0 15 +       100 ATGT---ATGCCG------------ATGTg\n"
                                             "s name3.chr1     0 28 +       100 GCAGCTGAAAACA--NNNNNNNNNNGGGGG\n",
                                             expectedList
                                             );
    row_t *r = stHash_search(expectedHash, "name3");
    r->prevRightPos = 54;
    r->strand = '*';
    r->prevStrand = '-';
    stHash *seqHash = createSeqHashFromString("name.chr1", "ATGTATGCCGacgtc"
                                              "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG"
                                              "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG");
    mtfseq_t *mtfs = newMtfseqFromString("gcagctgaaaacaACGTA"
                                         "tttttttttttttttttttttttttttttttt"
                                         "tttttttttttttttttttttttttttttttttttttttttttttttttt");
    stHash_insert(seqHash, stString_copy("reference.chr0"), mtfs);
    mtfs = newMtfseqFromString("ATGTATGCCGATGTg"
                               "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"
                               "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC");
    stHash_insert(seqHash, stString_copy("name2.chr1"), mtfs);
    mtfs = newMtfseqFromString("GCAGCTGAAAACAggggggggggggggggggggggggggggggggggggg"
                               "CCCCCaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"
                               );
    stHash_insert(seqHash, stString_copy("name3.chr1"), mtfs);
    addMafBlockToRowHash(observedHash, seqHash, observedList, mb, options);
    CuAssertTrue(testCase, hashesAreEqual(observedHash, expectedHash));
    CuAssertTrue(testCase, listsAreEqual(observedList, expectedList));
    // clean up
    stHash_destruct(observedHash);
    stHash_destruct(expectedHash);
    stHash_destruct(seqHash);
    stList_destruct(observedList);
    stList_destruct(expectedList);
    maf_destroyMafBlockList(mb);
    destroyOptions(options);
}
static void test_addBlockToHash_6(CuTest *testCase) {
    // concatenation with sequence breakpoint due to sequence name
    options_t *options = options_construct();
    options->breakpointPenalty = 10;
    options->interstitialSequence = 5;
    stList *observedList = stList_construct3(0, free);
    stList *expectedList = stList_construct3(0, free);
    stHash *observedHash = createBlockHashFromString("a score=0\n"
                                                     "s reference.chr0 0 13 + 158545518 gcagctgaaaaca\n"
                                                     "s name.chr1      0 10 +       100 ATGT---ATGCCG\n"
                                                     "s name2.chr1     0 10 +       100 ATGT---ATGCCG\n"
                                                     "s name3.chr1     0 13 +       100 GCAGCTGAAAACA\n",
                                                     observedList
                                                     );
    mafBlock_t *mb = maf_newMafBlockListFromString("a score=0 test\n"
                                                   "s reference.chr0 13  5 + 158545518 ACGTA\n"
                                                   "s name.chr1      12  5 +       100 gtcGG\n"
                                                   "s name2.chr1     10  5 +       100 ATGTg\n"
                                                   "s name3.chr2     13  5 +       100 aaaaa\n"
                                                   , 3);
    stHash *expectedHash = NULL;
    expectedHash = createBlockHashFromString("a score=0\n"
                                             "s reference.chr0 0 18 + 158545518 gcagctgaaaaca------------ACGTA\n"
                                             "s name.chr1      0 17 +       100 ATGT---ATGCCGac----------gtcGG\n"
                                             "s name2.chr1     0 15 +       100 ATGT---ATGCCG------------ATGTg\n"
                                             "s name3          0 28 +        28 GCAGCTGAAAACA--NNNNNNNNNNaaaaa\n",
                                             expectedList
                                             );
    row_t *r = stHash_search(expectedHash, "name3");
    r->prevRightPos = 17;
    r->strand = '*';
    r->prevStrand = '+';
    r->multipleNames = true;
    free(r->prevName);
    r->prevName = stString_copy("name3.chr2");
    stHash *seqHash = createSeqHashFromString("name.chr1", "ATGTATGCCGacgtc"
                                              "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG"
                                              "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG");
    mtfseq_t *mtfs = newMtfseqFromString("gcagctgaaaacaACGTA"
                                         "tttttttttttttttttttttttttttttttt"
                                         "tttttttttttttttttttttttttttttttttttttttttttttttttt");
    stHash_insert(seqHash, stString_copy("reference.chr0"), mtfs);
    mtfs = newMtfseqFromString("ATGTATGCCGATGTg"
                               "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"
                               "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC");
    stHash_insert(seqHash, stString_copy("name2.chr1"), mtfs);
    mtfs = newMtfseqFromString("GCAGCTGAAAACAggggggggggggggggggggggggggggggggggggg"
                               "CCCCCaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"
                               );
    stHash_insert(seqHash, stString_copy("name3.chr1"), mtfs);
    mtfs = newMtfseqFromString("TTTTTTTTTTTTTaaaaaGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG"
                               "CCCCCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
                               );
    stHash_insert(seqHash, stString_copy("name3.chr2"), mtfs);
    addMafBlockToRowHash(observedHash, seqHash, observedList, mb, options);
    CuAssertTrue(testCase, hashesAreEqual(observedHash, expectedHash));
    CuAssertTrue(testCase, listsAreEqual(observedList, expectedList));
    // clean up
    stHash_destruct(observedHash);
    stHash_destruct(expectedHash);
    stHash_destruct(seqHash);
    stList_destruct(observedList);
    stList_destruct(expectedList);
    maf_destroyMafBlockList(mb);
    destroyOptions(options);
}
CuSuite* mafToFastaStitcher_TestSuite(void) {
    // listing the tests as void allows us to quickly comment out certain tests
    // when trying to isolate bugs highlighted by one particular test
    (void) test_newBlockHashFromBlock_0;
    (void) test_addMafLineToRow_0;
    (void) test_addMafLineToRow_1;
    (void) test_penalize_0;
    (void) test_interstitial_0;
    (void) test_addBlockToHash_0;
    (void) test_addBlockToHash_1;
    (void) test_addBlockToHash_2;
    (void) test_addBlockToHash_3;
    (void) test_addBlockToHash_4;
    (void) test_addBlockToHash_5;
    (void) test_addBlockToHash_6;
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, test_newBlockHashFromBlock_0);
    SUITE_ADD_TEST(suite, test_addMafLineToRow_0);
    SUITE_ADD_TEST(suite, test_addMafLineToRow_1);
    SUITE_ADD_TEST(suite, test_penalize_0);
    SUITE_ADD_TEST(suite, test_interstitial_0);
    SUITE_ADD_TEST(suite, test_addBlockToHash_0);
    SUITE_ADD_TEST(suite, test_addBlockToHash_1);
    SUITE_ADD_TEST(suite, test_addBlockToHash_2);
    SUITE_ADD_TEST(suite, test_addBlockToHash_3);
    SUITE_ADD_TEST(suite, test_addBlockToHash_4);
    SUITE_ADD_TEST(suite, test_addBlockToHash_5);
    SUITE_ADD_TEST(suite, test_addBlockToHash_6);
    return suite;
}
