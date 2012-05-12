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
#ifndef TEST_SHAREDMAF_H_
#define TEST_SHAREDMAF_H_
#include <assert.h>
#include <stdbool.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <unistd.h>
#include "CuTest.h"
#include "common.h"
#include "sharedMaf.h"

int createTmpFolder(void) {
    // create a temporary directory with permissions only for the owner
    return mkdir("test_tmp", S_IRWXU | S_IRUSR | S_IXUSR | S_IWUSR);
}
void writeStringToTmpFile(char *s) {
    FILE *f = de_open("test_tmp/test.maf", "w+");
    fprintf(f, "%s", s);
    fclose(f);
}
bool filesAreIdentical(char *fileA, char *fileB) {
    extern const int kMaxStringLength;
    int32_t n = kMaxStringLength;
    FILE *ifpA = fopen(fileA, "r");
    FILE *ifpB = fopen(fileB, "r");
    char *lineA = (char*) de_malloc(n);
    char *lineB = (char*) de_malloc(n);
    while(de_getline(&lineA, &n, ifpA) != -1) {
        if (de_getline(&lineB, &n, ifpB) == -1) {
            free(lineA);
            free(lineB);
            fclose(ifpA);
            fclose(ifpB);
            return false;
        }
        if (strlen(lineA) != strlen(lineB)) {
            free(lineA);
            free(lineB);
            fclose(ifpA);
            fclose(ifpB);
            return false;
        }
        for (unsigned i = 0; i < strlen(lineA); ++i) {
            if (lineA[i] != lineB[i]) {
                free(lineA);
                free(lineB);
                fclose(ifpA);
                fclose(ifpB);
                return false;
            }
        }
            
    }
    if (de_getline(&lineB, &n, ifpB) != -1) {
        free(lineA);
        free(lineB);
        fclose(ifpA);
        fclose(ifpB);
        return false;
    }
    free(lineA);
    free(lineB);
    fclose(ifpA);
    fclose(ifpB);
    return true;
}
static void test_newMafLineFromString(CuTest *testCase) {
    // verify that a maf line string is correctly parsed
    assert(testCase != NULL);
    // case 1
    char *input = de_strdup("s hg16.chr7    27707221 13 + 158545518 gcagctgaaaaca");
    mafLine_t *ml = maf_newMafLineFromString(input, 1);
    CuAssertUInt32Equals(testCase, ml->type, 's');
    CuAssertStrEquals(testCase, ml->line, input);
    CuAssertStrEquals(testCase, ml->species, "hg16.chr7");
    CuAssertUInt32Equals(testCase, (int) ml->start, 27707221);
    CuAssertUInt32Equals(testCase, (int) ml->length, 13);
    CuAssertUInt32Equals(testCase, (int) ml->sourceLength, 158545518);
    CuAssertIntEquals(testCase, ml->strand, '+');
    CuAssertStrEquals(testCase, ml->sequence, "gcagctgaaaaca");
    CuAssertUInt32Equals(testCase, (int) ml->lineNumber, 1);
    CuAssertTrue(testCase, ml->next == NULL);
    free(input);
    maf_destroyMafLineList(ml);
    // case 2
    input = de_strdup("i panTro1.chr6 N 0 C 0");
    ml = maf_newMafLineFromString(input, 2);
    CuAssertUInt32Equals(testCase, ml->type, 'i');
    CuAssertStrEquals(testCase, ml->line, input);
    CuAssertUInt32Equals(testCase, (int) ml->lineNumber, 2);
    CuAssertTrue(testCase, ml->species == NULL);
    CuAssertTrue(testCase, ml->start == 0);
    CuAssertTrue(testCase, ml->length == 0);
    CuAssertTrue(testCase, ml->sourceLength == 0);
    CuAssertTrue(testCase, ml->strand == 0);
    CuAssertTrue(testCase, ml->sequence == NULL);
    free(input);
    maf_destroyMafLineList(ml);
    // case 3
    input = de_strdup("e mm4.chr6     53310102 13 + 151104725 I");
    ml = maf_newMafLineFromString(input, 3);
    CuAssertUInt32Equals(testCase, ml->type, 'e');
    CuAssertStrEquals(testCase, ml->line, input);
    CuAssertUInt32Equals(testCase, (int) ml->lineNumber, 3);
    CuAssertTrue(testCase, ml->species == NULL);
    CuAssertTrue(testCase, ml->start == 0);
    CuAssertTrue(testCase, ml->length == 0);
    CuAssertTrue(testCase, ml->sourceLength == 0);
    CuAssertTrue(testCase, ml->strand == 0);
    CuAssertTrue(testCase, ml->sequence == NULL);
    free(input);
    maf_destroyMafLineList(ml);
    // case 4
    input = de_strdup("s hg18.chr7    27578828 38 + 158545518 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG");
    ml = maf_newMafLineFromString(input, 4);
    CuAssertUInt32Equals(testCase, ml->type, 's');
    CuAssertStrEquals(testCase, ml->line, input);
    CuAssertUInt32Equals(testCase, (int) ml->lineNumber, 4);
    CuAssertStrEquals(testCase, ml->species, "hg18.chr7");
    CuAssertUInt32Equals(testCase, (int) ml->start, 27578828);
    CuAssertUInt32Equals(testCase, (int) ml->length, 38);
    CuAssertUInt32Equals(testCase, (int) ml->sourceLength, 158545518);
    CuAssertIntEquals(testCase, ml->strand, '+');
    CuAssertStrEquals(testCase, ml->sequence, "AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG");
    CuAssertTrue(testCase, ml->next == NULL);
    free(input);
    maf_destroyMafLineList(ml);
    // case 5
    input = de_strdup("a score=23262.0 ");
    ml = maf_newMafLineFromString(input, 5);
    CuAssertUInt32Equals(testCase, ml->type, 'a');
    CuAssertStrEquals(testCase, ml->line, input);
    CuAssertUInt32Equals(testCase, ml->lineNumber, 5);
    CuAssertTrue(testCase, ml->species == NULL);
    CuAssertTrue(testCase, ml->start == 0);
    CuAssertTrue(testCase, ml->length == 0);
    CuAssertTrue(testCase, ml->sourceLength == 0);
    CuAssertTrue(testCase, ml->strand == 0);
    CuAssertTrue(testCase, ml->sequence == NULL);
    free(input);
    maf_destroyMafLineList(ml);
    // case 6
    input = de_strdup("s rn3     81444246 6 - 187371129 taagga");
    ml = maf_newMafLineFromString(input, 6);
    CuAssertUInt32Equals(testCase, ml->type, 's');
    CuAssertStrEquals(testCase, ml->line, input);
    CuAssertUInt32Equals(testCase, ml->lineNumber, 6);
    CuAssertStrEquals(testCase, ml->species, "rn3");
    CuAssertUInt32Equals(testCase, ml->start, 81444246);
    CuAssertUInt32Equals(testCase, ml->length, 6);
    CuAssertUInt32Equals(testCase, ml->sourceLength, 187371129);
    CuAssertIntEquals(testCase, ml->strand, '-');
    CuAssertStrEquals(testCase, ml->sequence, "taagga");
    CuAssertTrue(testCase, ml->next == NULL);
    free(input);
    maf_destroyMafLineList(ml);
}
static void test_readBlock(CuTest *testCase) {
    // verify we read a header and a block correctly
    assert(testCase != NULL);
    createTmpFolder();
    char *input = de_strdup("track name=euArc visibility=pack \n\
##maf version=1 scoring=tba.v8 \n\
# tba.v8 (((human chimp) baboon) (mouse rat)) \n\
                   \n\
\n\
a score=23262.0     \n\
s hg18.chr7    27578828 38 + 158545518 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG\n\
s panTro1.chr6 28741140 38 + 161576975 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG\n\
s baboon         116834 38 +   4622798 AAA-GGGAATGTTAACCAAATGA---GTTGTCTCTTATGGTG\n\
s mm4.chr6     53215344 38 + 151104725 -AATGGGAATGTTAAGCAAACGA---ATTGTCTCTCAGTGTG\n\
s rn3.chr4     81344243 40 + 187371129 -AA-GGGGATGCTAAGCCAATGAGTTGTTGTCTCTCAATGTG\n\
                   \n\
a score=5062.0                    \n\
s hg18.chr7    27699739 6 + 158545518 TAAAGA\n\
s panTro1.chr6 28862317 6 + 161576975 TAAAGA\n\
s baboon         241163 6 +   4622798 TAAAGA \n\
s mm4.chr6     53303881 6 + 151104725 TAAAGA\n\
s rn3.chr4     81444246 6 + 187371129 taagga\n\
\n\
# non block comment line \n\
\n\
\n\
a score=6636.0\n\
s hg18.chr7    27707221 13 + 158545518 gcagctgaaaaca\n\
s panTro1.chr6 28869787 13 + 161576975 gcagctgaaaaca\n\
s baboon         249182 13 +   4622798 gcagctgaaaaca\n\
s mm4.chr6     53310102 13 + 151104725 ACAGCTGAAAATA\n\
\n\n");
    writeStringToTmpFile(input);
    mafFileApi_t *mapi = maf_newMfa("test_tmp/test.maf", "r");
    mafBlock_t *mb = maf_readBlock(mapi);
    CuAssertTrue(testCase, mb != NULL);
    mafLine_t *ml = mb->headLine;
    CuAssertTrue(testCase, maf_numberOfSequencesMafLineList(ml) == 0);
    CuAssertStrEquals(testCase, ml->line, "track name=euArc visibility=pack ");
    CuAssertTrue(testCase, ml->next != NULL);
    ml = ml->next;
    CuAssertStrEquals(testCase, ml->line, "##maf version=1 scoring=tba.v8 ");
    CuAssertTrue(testCase, ml->next != NULL);
    ml = ml->next;
    CuAssertStrEquals(testCase, ml->line, "# tba.v8 (((human chimp) baboon) (mouse rat)) ");
    CuAssertTrue(testCase, ml->next == NULL);
    CuAssertTrue(testCase, mb->next == NULL);
    maf_destroyMafBlockList(mb);
    mb = maf_readBlock(mapi);
    CuAssertTrue(testCase, mb != NULL);
    ml = mb->headLine;
    CuAssertStrEquals(testCase, ml->line, "a score=23262.0     ");
    CuAssertTrue(testCase, ml->next != NULL);
    ml = ml->next;
    CuAssertTrue(testCase, maf_numberOfSequencesMafLineList(ml) == 5);
    CuAssertStrEquals(testCase, ml->line, "s hg18.chr7    27578828 38 + 158545518 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG");
    CuAssertStrEquals(testCase, ml->species, "hg18.chr7");
    CuAssertUInt32Equals(testCase, ml->start, 27578828);
    CuAssertUInt32Equals(testCase, ml->length, 38);
    CuAssertIntEquals(testCase, ml->strand, '+');
    CuAssertUInt32Equals(testCase, ml->sourceLength, 158545518);
    CuAssertStrEquals(testCase, ml->sequence, "AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG");
    CuAssertTrue(testCase, ml->next != NULL);
    ml = ml->next;
    // clean up
    maf_destroyMafBlockList(mb);
    unlink("test_tmp/test.maf");
    rmdir("test_tmp");
    maf_destroyMfa(mapi);
    free(input);
}
static void test_readBlock2(CuTest *testCase) {
    // This variation looks at a maf with no blank line after the header
    assert(testCase != NULL);
    createTmpFolder();
    // case 1
    char *input = de_strdup("track name=euArc visibility=pack \n\
##maf version=1 scoring=tba.v8 \n\
# tba.v8 (((human chimp) baboon) (mouse rat)) \n\
a score=23262.0     \n\
s hg18.chr7    27578828 38 + 158545518 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG\n\
s panTro1.chr6 28741140 38 + 161576975 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG\n\
s baboon         116834 38 +   4622798 AAA-GGGAATGTTAACCAAATGA---GTTGTCTCTTATGGTG\n\
s mm4.chr6     53215344 38 + 151104725 -AATGGGAATGTTAAGCAAACGA---ATTGTCTCTCAGTGTG\n\
s rn3.chr4     81344243 40 + 187371129 -AA-GGGGATGCTAAGCCAATGAGTTGTTGTCTCTCAATGTG\n\
                   \n\
a score=5062.0                    \n\
s hg18.chr7    27699739 6 + 158545518 TAAAGA\n\
s panTro1.chr6 28862317 6 + 161576975 TAAAGA\n\
s baboon         241163 6 +   4622798 TAAAGA \n\
s mm4.chr6     53303881 6 + 151104725 TAAAGA\n\
s rn3.chr4     81444246 6 + 187371129 taagga\n\
\n\
a score=6636.0\n\
s hg18.chr7    27707221 13 + 158545518 gcagctgaaaaca\n\
s panTro1.chr6 28869787 13 + 161576975 gcagctgaaaaca\n\
s baboon         249182 13 +   4622798 gcagctgaaaaca\n\
s mm4.chr6     53310102 13 + 151104725 ACAGCTGAAAATA\n\
\n\n");
    writeStringToTmpFile(input);
    mafFileApi_t *mapi = maf_newMfa("test_tmp/test.maf", "r");
    mafBlock_t *mb = maf_readBlock(mapi);
    CuAssertTrue(testCase, mb != NULL);
    mafLine_t *ml = mb->headLine;
    CuAssertTrue(testCase, maf_numberOfSequencesMafLineList(ml) == 0);
    CuAssertStrEquals(testCase, ml->line, "track name=euArc visibility=pack ");
    CuAssertTrue(testCase, ml->next != NULL);
    ml = ml->next;
    CuAssertStrEquals(testCase, ml->line, "##maf version=1 scoring=tba.v8 ");
    CuAssertTrue(testCase, ml->next != NULL);
    ml = ml->next;
    CuAssertStrEquals(testCase, ml->line, "# tba.v8 (((human chimp) baboon) (mouse rat)) ");
    CuAssertTrue(testCase, ml->next == NULL);
    CuAssertTrue(testCase, mb->next == NULL);
    maf_destroyMafBlockList(mb);
    mb = maf_readBlock(mapi);
    CuAssertTrue(testCase, mb != NULL);
    ml = mb->headLine;
    CuAssertStrEquals(testCase, ml->line, "a score=23262.0     ");
    CuAssertTrue(testCase, ml->next != NULL);
    ml = ml->next;
    CuAssertTrue(testCase, maf_numberOfSequencesMafLineList(ml) == 5);
    CuAssertStrEquals(testCase, ml->line, "s hg18.chr7    27578828 38 + 158545518 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG");
    CuAssertStrEquals(testCase, ml->species, "hg18.chr7");
    CuAssertUInt32Equals(testCase, ml->start, 27578828);
    CuAssertUInt32Equals(testCase, ml->length, 38);
    CuAssertIntEquals(testCase, ml->strand, '+');
    CuAssertUInt32Equals(testCase, ml->sourceLength, 158545518);
    CuAssertStrEquals(testCase, ml->sequence, "AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG");
    CuAssertTrue(testCase, ml->next != NULL);
    ml = ml->next;
    // clean up
    maf_destroyMafBlockList(mb);
    unlink("test_tmp/test.maf");
    rmdir("test_tmp");
    maf_destroyMfa(mapi);
    free(input);
}
static void test_lineNumbers(CuTest *testCase) {
    // This test makes sure that the line numbers are being correctly recorded in mafLine_t
    assert(testCase != NULL);
    createTmpFolder();
    // case 1
    char *input = de_strdup("track name=euArc visibility=pack \n\
##maf version=1 scoring=tba.v8 \n\
# tba.v8 (((human chimp) baboon) (mouse rat)) \n\
\n\
\n\
a score=23262.0     \n\
s hg18.chr7    27578828 38 + 158545518 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG\n\
# generic comment\n\
s panTro1.chr6 28741140 38 + 161576975 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG\n\
s baboon         116834 38 +   4622798 AAA-GGGAATGTTAACCAAATGA---GTTGTCTCTTATGGTG\n\
s mm4.chr6     53215344 38 + 151104725 -AATGGGAATGTTAAGCAAACGA---ATTGTCTCTCAGTGTG\n\
s rn3.chr4     81344243 40 + 187371129 -AA-GGGGATGCTAAGCCAATGAGTTGTTGTCTCTCAATGTG\n\
\n\n");
    writeStringToTmpFile(input);
    mafFileApi_t *mapi = maf_newMfa("test_tmp/test.maf", "r");
    mafBlock_t *mb = maf_readBlock(mapi);
    CuAssertTrue(testCase, mb != NULL);
    mafLine_t *ml = mb->headLine;
    CuAssertTrue(testCase, maf_numberOfSequencesMafLineList(ml) == 0);
    CuAssertStrEquals(testCase, ml->line, "track name=euArc visibility=pack ");
    CuAssertTrue(testCase, ml->lineNumber == 1);
    CuAssertTrue(testCase, ml->next != NULL);
    ml = ml->next;
    CuAssertStrEquals(testCase, ml->line, "##maf version=1 scoring=tba.v8 ");
    CuAssertTrue(testCase, ml->lineNumber == 2);
    CuAssertTrue(testCase, ml->next != NULL);
    ml = ml->next;
    CuAssertStrEquals(testCase, ml->line, "# tba.v8 (((human chimp) baboon) (mouse rat)) ");
    CuAssertTrue(testCase, ml->lineNumber == 3);
    CuAssertTrue(testCase, ml->next == NULL);
    CuAssertTrue(testCase, mb->next == NULL);
    maf_destroyMafBlockList(mb);
    mb = maf_readBlock(mapi);
    CuAssertTrue(testCase, mb != NULL);
    ml = mb->headLine;
    CuAssertStrEquals(testCase, ml->line, "a score=23262.0     ");
    CuAssertTrue(testCase, ml->lineNumber == 6);
    CuAssertTrue(testCase, ml->next != NULL);
    ml = ml->next;
    CuAssertTrue(testCase, maf_numberOfSequencesMafLineList(ml) == 5);
    CuAssertStrEquals(testCase, ml->line, "s hg18.chr7    27578828 38 + 158545518 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG");
    CuAssertTrue(testCase, ml->lineNumber == 7);
    CuAssertStrEquals(testCase, ml->species, "hg18.chr7");
    CuAssertUInt32Equals(testCase, ml->start, 27578828);
    CuAssertUInt32Equals(testCase, ml->length, 38);
    CuAssertIntEquals(testCase, ml->strand, '+');
    CuAssertUInt32Equals(testCase, ml->sourceLength, 158545518);
    CuAssertStrEquals(testCase, ml->sequence, "AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG");
    CuAssertTrue(testCase, ml->next != NULL);
    ml = ml->next;
    CuAssertStrEquals(testCase, ml->line, "# generic comment");
    CuAssertTrue(testCase, ml->lineNumber == 8);
    CuAssertTrue(testCase, ml->next != NULL);
    ml = ml->next;
    CuAssertStrEquals(testCase, ml->line, "s panTro1.chr6 28741140 38 + 161576975 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG");
    CuAssertTrue(testCase, ml->lineNumber == 9);
    // clean up
    maf_destroyMafBlockList(mb);
    unlink("test_tmp/test.maf");
    rmdir("test_tmp");
    maf_destroyMfa(mapi);
    free(input);
}
static void test_readWriteMaf(CuTest *testCase) {
    // make sure that we can do a complete read and complete write and that
    // the output is identical to the input.
    assert(testCase != NULL);
    createTmpFolder();
    // case 1
    char *input = de_strdup("track name=euArc visibility=pack \n\
##maf version=1 scoring=tba.v8 \n\
# tba.v8 (((human chimp) baboon) (mouse rat)) \n\
\n\
a score=23262.0     \n\
s hg18.chr7    27578828 38 + 158545518 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG\n\
# generic comment\n\
s panTro1.chr6 28741140 38 + 161576975 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG\n\
s baboon         116834 38 +   4622798 AAA-GGGAATGTTAACCAAATGA---GTTGTCTCTTATGGTG\n\
s mm4.chr6     53215344 38 + 151104725 -AATGGGAATGTTAAGCAAACGA---ATTGTCTCTCAGTGTG\n\
s rn3.chr4     81344243 40 + 187371129 -AA-GGGGATGCTAAGCCAATGAGTTGTTGTCTCTCAATGTG\n\
\n\n");
    writeStringToTmpFile(input);
    mafFileApi_t *mfaRead = maf_newMfa("test_tmp/test.maf", "r");
    CuAssertTrue(testCase, mfaRead != NULL);
    mafFileApi_t *mfaWrite = maf_newMfa("test_tmp/testWrite.maf", "w");
    CuAssertTrue(testCase, mfaWrite != NULL);
    mafBlock_t *mb = maf_readAll(mfaRead);
    CuAssertTrue(testCase, mb != NULL);
    maf_writeAll(mfaWrite, mb);
    CuAssertTrue(testCase, mfaRead->lineNumber == mfaWrite->lineNumber);
    CuAssertTrue(testCase, filesAreIdentical(mfaRead->filename, mfaWrite->filename) == true);
    // clean up
    maf_destroyMafBlockList(mb);
    unlink("test_tmp/test.maf");
    unlink("test_tmp/testWrite.maf");
    rmdir("test_tmp");
    maf_destroyMfa(mfaRead);
    maf_destroyMfa(mfaWrite);
    free(input);
}
CuSuite* mafShared_TestSuite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, test_newMafLineFromString);
    SUITE_ADD_TEST(suite, test_readBlock);
    SUITE_ADD_TEST(suite, test_readBlock2);
    SUITE_ADD_TEST(suite, test_lineNumbers);
    SUITE_ADD_TEST(suite, test_readWriteMaf);
    return suite;
}
#endif // TEST_SHAREDMAF_H_
