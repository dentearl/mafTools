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
#include "comparatorRandom.h"

static void test_rbinom_zero_0(CuTest *testCase) {
    // small n means using one of the small algorithms
    uint64_t a;
    uint64_t N = 100000;
    uint64_t n = 10;
    for (uint32_t i = 0; i < N; ++i) {
        a = rbinom(n, 0.0);
        CuAssertTrue(testCase, a == 0);
    }
}
static void test_rbinom_zero_1(CuTest *testCase) {
    // large n means using the btpe algorithm
    uint64_t a;
    uint64_t N = 100000;
    uint64_t n = 100;
    for (uint32_t i = 0; i < N; ++i) {
        a = rbinom(n, 0.0);
        CuAssertTrue(testCase, a == 0);
    }
}
static void test_rbinom_one_0(CuTest *testCase) {
    // small n means using one of the small algorithms
    uint64_t a;
    uint64_t N = 100000;
    uint64_t n = 10;
    for (uint64_t i = 0; i < N; ++i) {
        a = rbinom(n, 1.0);
        CuAssertTrue(testCase, a == n);
    }
}
static void test_rbinom_one_1(CuTest *testCase) {
    // large n means using the btpe algorithm
    uint64_t a;
    uint64_t N = 100000;
    uint64_t n = 100;
    for (uint64_t i = 0; i < N; ++i) {
        a = rbinom(n, 1.0);
        CuAssertTrue(testCase, a == n);
    }
}
static void test_rbinom_one_2(CuTest *testCase) {
    // small n means using one of the small algorithms
    uint64_t a;
    uint64_t N = 100000;
    uint64_t n = 10;
    for (uint64_t i = 0; i < N; ++i) {
        a = rbinom(n, 100.0);
        CuAssertTrue(testCase, a == n);
    }
}
static void test_rbinom_one_3(CuTest *testCase) {
    // large n means using the btpe algorithm
    uint64_t a;
    uint64_t N = 100000;
    uint64_t n = 100;
    for (uint64_t i = 0; i < N; ++i) {
        a = rbinom(n, 1.0);
        CuAssertTrue(testCase, a == n);
    }
}
static double runningAverage(double m, uint64_t x, uint64_t i) {
    m += (x - m) / (i + 1);
    return m;
}
static void test_rbinom_clt_0(CuTest *testCase) {
    // central limit theorm for small n
    double m = 0.0;
    uint64_t N = 1000000;
    uint64_t n = 10;
    for (uint64_t i = 0; i < N; ++i) {
        m = runningAverage(m, rbinom(n, 0.5), i);
    }
    CuAssertTrue(testCase, m > n * 0.49);
    CuAssertTrue(testCase, m < n * 0.51);
}
static void test_rbinom_clt_1(CuTest *testCase) {
    // central limit theorm for large n
    double m = 0.0;
    uint64_t N = 1000000;
    uint64_t n = 100;
    for (uint64_t i = 0; i < N; ++i) {
        m = runningAverage(m, rbinom(n, 0.5), i);
    }
    CuAssertTrue(testCase, m > n * 0.49);
    CuAssertTrue(testCase, m < n * 0.51);
}
CuSuite* comparatorRandom_TestSuite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, test_rbinom_zero_0);
    SUITE_ADD_TEST(suite, test_rbinom_zero_1);
    SUITE_ADD_TEST(suite, test_rbinom_one_0);
    SUITE_ADD_TEST(suite, test_rbinom_one_1);
    SUITE_ADD_TEST(suite, test_rbinom_one_2);
    SUITE_ADD_TEST(suite, test_rbinom_one_3);
    SUITE_ADD_TEST(suite, test_rbinom_clt_0);
    SUITE_ADD_TEST(suite, test_rbinom_clt_1);

    return suite;
}
