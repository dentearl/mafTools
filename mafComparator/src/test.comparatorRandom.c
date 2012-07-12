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
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include "CuTest.h"
#include "common.h"
#include "sonLib.h"
#include "comparatorRandom.h"

static int cmpu64(const void *a, const void *b);
static void test_rbinom_zero_0(CuTest *testCase) {
    // small n*p means using one of the small algorithms
    uint64_t a;
    uint64_t N = 100000;
    uint64_t n = 10;
    for (uint32_t i = 0; i < N; ++i) {
        a = rbinom(n, 0.0);
        CuAssertTrue(testCase, a == 0);
    }
}
static void test_rbinom_zero_1(CuTest *testCase) {
    // large n*p means using the btpe algorithm
    uint64_t a;
    uint64_t N = 100000;
    uint64_t n = 100;
    for (uint32_t i = 0; i < N; ++i) {
        a = rbinom(n, 0.0);
        CuAssertTrue(testCase, a == 0);
    }
}
static void test_rbinom_one_0(CuTest *testCase) {
    // small n*p means using one of the small algorithms
    uint64_t a;
    uint64_t N = 100000;
    uint64_t n = 10;
    for (uint64_t i = 0; i < N; ++i) {
        a = rbinom(n, 1.0);
        CuAssertTrue(testCase, a == n);
    }
}
static void test_rbinom_one_1(CuTest *testCase) {
    // large n*p means using the btpe algorithm
    uint64_t a;
    uint64_t N = 100000;
    uint64_t n = 100;
    for (uint64_t i = 0; i < N; ++i) {
        a = rbinom(n, 1.0);
        CuAssertTrue(testCase, a == n);
    }
}
static void test_rbinom_one_2(CuTest *testCase) {
    // small n*p means using one of the small algorithms
    uint64_t a;
    uint64_t N = 100000;
    uint64_t n = 10;
    for (uint64_t i = 0; i < N; ++i) {
        a = rbinom(n, 100.0);
        CuAssertTrue(testCase, a == n);
    }
}
static void test_rbinom_one_3(CuTest *testCase) {
    // large n*p means using the btpe algorithm
    uint64_t a;
    uint64_t N = 100000;
    uint64_t n = 100;
    for (uint64_t i = 0; i < N; ++i) {
        a = rbinom(n, 1.0);
        CuAssertTrue(testCase, a == n);
    }
}
static void test_rbinom_one_4(CuTest *testCase) {
    // large n*p means using the btpe algorithm
    uint64_t a;
    uint64_t N = 100000;
    uint64_t n = INT32_MAX;
    n *= 2;
    for (uint64_t i = 0; i < N; ++i) {
        a = rbinom(n, 1.0);
        CuAssertTrue(testCase, a == n);
    }
}
static void test_rbinom_ranges_0(CuTest *testCase) {
    // ranges for small n*p
    uint64_t N = 1000000;
    uint64_t n = 10;
    uint64_t *array = (uint64_t*) st_malloc(sizeof(*array) * N);
    for (uint64_t i = 0; i < N; ++i) {
        array[i] = rbinom(n, 0.01);
        CuAssertTrue(testCase, array[i] <= n);
    }
    free(array);
}
static void test_rbinom_ranges_1(CuTest *testCase) {
    // ranges for large n*p
    uint64_t N = 1000000;
    uint64_t n = 10000;
    uint64_t *array = (uint64_t*) st_malloc(sizeof(*array) * N);
    for (uint64_t i = 0; i < N; ++i) {
        array[i] = rbinom(n, 0.01);
        CuAssertTrue(testCase, array[i] <= n);
    }
    free(array);
}
static void headarray(uint64_t *a, uint64_t n) {
    int c = 0;
    for (uint64_t i = 0; i < n; ++i) {
        ++c;
        printf("%"PRIu64",%s", a[i], c % 20 ? " " : "\n");
    }
    printf("\n");
}
static void tailarray(uint64_t *a, uint64_t N, uint64_t n) {
    int c = 0;
    for (uint64_t i = N - 1 - n; i < N; ++i) {
        ++c;
        printf("%"PRIu64",%s", a[i], c % 20 ? " " : "\n");
    }
    printf("\n");
}
static void test_rbinom_ranges_2(CuTest *testCase) {
    // ranges for large n*p
    (void) (headarray);
    (void) (tailarray);
    uint64_t N = 1000000;
    uint64_t n = 10;
    uint64_t *array = (uint64_t*) st_malloc(sizeof(*array) * N);
    for (uint64_t i = 0; i < N; ++i) {
        array[i] = rbinom(n, 0.5);
        CuAssertTrue(testCase, array[i] <= n);
    }
    qsort(array, N, sizeof(uint64_t), cmpu64);
    CuAssertTrue(testCase, array[0] == 0);
    CuAssertTrue(testCase, array[N - 1] == n);
    free(array);
}
static void test_rbinom_ranges_3(CuTest *testCase) {
    // ranges for large n*p
    uint64_t N = 1000000;
    uint64_t n = INT32_MAX;
    n *= 3;
    uint64_t *array = (uint64_t*) st_malloc(sizeof(*array) * N);
    for (uint64_t i = 0; i < N; ++i) {
        array[i] = rbinom(n, 0.01);
        CuAssertTrue(testCase, array[i] <= n);
    }
    free(array);
}
static void test_rbinom_ranges_4(CuTest *testCase) {
    // ranges for large n*p
    uint64_t N = 1000000;
    uint64_t n = INT32_MAX;
    uint64_t m = n * 2.5; // we must be careful of overflowing constant expressions here
    n *=  3;
    uint64_t *array = (uint64_t*) st_malloc(sizeof(*array) * N);
    for (uint64_t i = 0; i < N; ++i) {
        array[i] = rbinom(n, 0.99);
        CuAssertTrue(testCase, array[i] > m);
    }
    free(array);
}
static double runningAverage(double m, uint64_t x, uint64_t i) {
    m += (x - m) / (i + 1);
    return m;
}
static void test_rbinom_clt_0(CuTest *testCase) {
    // central limit theorm for small n*p
    double m = 0.0, p = 0.5;
    uint64_t N = 1000000;
    uint64_t n = 10;
    for (uint64_t i = 0; i < N; ++i) {
        m = runningAverage(m, rbinom(n, p), i);
    }
    CuAssertTrue(testCase, m < n * (p + 0.001));
    CuAssertTrue(testCase, m > n * (p - 0.001));
}
static void test_rbinom_clt_1(CuTest *testCase) {
    // central limit theorm for large n*p
    double m = 0.0, p = 0.5;
    uint64_t N = 1000000;
    uint64_t n = 100;
    for (uint64_t i = 0; i < N; ++i) {
        m = runningAverage(m, rbinom(n, p), i);
    }
    CuAssertTrue(testCase, m < n * (p + 0.001));
    CuAssertTrue(testCase, m > n * (p - 0.001));
}
static int cmpu64(const void *a, const void *b) {
    const uint64_t *ua = (const uint64_t *)a;
    const uint64_t *ub = (const uint64_t *)b;
    if (*ua < *ub) {
        return -1;
    } else {
        return 1;
    }
    return 0;
}
static double median(uint64_t *array, uint64_t n) {
    qsort(array, n, sizeof(uint64_t), cmpu64);
    if (n % 2) {
        return (double) array[(n + 1) / 2 - 1];
    } else {
        return (double) ((array[n / 2 - 1] + array[(n / 2)]) / 2.0);
    }
}
static double mean(uint64_t *array, uint64_t n) {
    double mu = 0.0;
    for (uint64_t i = 0; i < n; ++i) {
        mu += array[i];
    }
    return mu / n;
}
static double sv(uint64_t *array, uint64_t n, double mu) {
    // sample variance.
    double sv = 0.0;
    assert(n > 1);
    for (uint64_t i = 0; i < n; ++i) {
        sv += (array[i] - mu) * (array[i] - mu);
    }
    sv /= (n - 1);
    return sv;
}
static void test_rbinom_distribution_0(CuTest *testCase) {
    // central limit theorm for small n*p
    double mu, var, med, p = 0.5;
    uint64_t N = 100000;
    uint64_t n = 10;
    uint64_t *array = NULL;
    for (uint64_t j = 0; j < 100; ++j) {
        mu = 0.0;
        array = (uint64_t*) st_malloc(sizeof(*array) * N);
        for (uint64_t i = 0; i < N; ++i) {
            array[i] = rbinom(n, p);
            mu = runningAverage(mu, array[i], i);
        }
        var = sv(array, N, mu);
        med = median(array, N);
        CuAssertTrue(testCase, array[0] <= array[N-1]);
        CuAssertTrue(testCase, array[0] <= n); // unsigned wraparound    
        CuAssertTrue(testCase, array[N-1] <= n);
        CuAssertTrue(testCase, med == 5.0);
        CuAssertTrue(testCase, mu > n * (p - 0.01));
        CuAssertTrue(testCase, mu < n * (p + 0.01));
        CuAssertTrue(testCase, var < 1.1 * n * p * (1.0 - p));
        CuAssertTrue(testCase, var > 0.9 * n * p * (1.0 - p));
        free(array);
    }
}
static void parray(uint64_t *array, uint64_t n) {
    int c = 0;
    for (uint64_t i = 0; i < n; ++i) {
        ++c;
        printf("%"PRIu64",%s", array[i], c % 20 ? " " : "\n");
    }
    printf("\n");
}
static void test_rbinom_distribution_1(CuTest *testCase) {
    // central limit theorm for large n*p
    (void) (parray);
    double mu, med, var, p = 0.5;
    uint64_t N = 100000;
    uint64_t n = 100;
    uint64_t *array = NULL;
    for (uint64_t j = 0; j < 100; ++j) {
        mu = 0.0;
        array = (uint64_t*) st_malloc(sizeof(*array) * N);
        for (uint64_t i = 0; i < N; ++i) {
            array[i] = rbinom(n, p);
            mu = runningAverage(mu, array[i], i);
        }
        var = sv(array, N, mu);
        med = median(array, N);
        CuAssertTrue(testCase, array[0] <= array[N-1]);
        CuAssertTrue(testCase, array[0] <= n); // unsigned wraparound    
        CuAssertTrue(testCase, array[N-1] <= n);
        CuAssertTrue(testCase, med == 50.0);
        CuAssertTrue(testCase, mu > n * (p - 0.01));
        CuAssertTrue(testCase, mu < n * (p + 0.01));
        CuAssertTrue(testCase, var < 1.1 * n * p * (1.0 - p));
        CuAssertTrue(testCase, var > 0.9 * n * p * (1.0 - p));
        free(array);
    }
}
static void test_rbinom_distribution_2(CuTest *testCase) {
    // central limit theorm for small-ish n*p
    (void) (parray);
    double mu, var, p, med;
    uint64_t N = 100000;
    uint64_t n = 10;
    uint64_t *array = NULL;
    for (uint64_t j = 0; j < 100; ++j) {
        p = st_random();
        mu = 0.0;
        array = (uint64_t*) st_malloc(sizeof(*array) * N);
        for (uint64_t i = 0; i < N; ++i) {
            array[i] = rbinom(n, p);
            mu = runningAverage(mu, array[i], i);
        }
        var = sv(array, N, mu);
        med = median(array, N);
        
        CuAssertTrue(testCase, array[0] <= array[N-1]);
        CuAssertTrue(testCase, array[0] <= n); // unsigned wraparound    
        CuAssertTrue(testCase, array[N-1] <= n);
        CuAssertTrue(testCase, ((med == floor(n * p)) || 
                                (med == ceil(n * p))));
        CuAssertTrue(testCase, mu > n * (p - 0.01));
        CuAssertTrue(testCase, mu < n * (p + 0.01));
        CuAssertTrue(testCase, var < 1.1 * n * p * (1.0 - p));
        CuAssertTrue(testCase, var > 0.9 * n * p * (1.0 - p));
        free(array);
    }
}
static void test_rbinom_distribution_3(CuTest *testCase) {
    // central limit theorm for large-ish n*p
    (void) (parray);
    double mu, med, var, p;
    uint64_t N = 100000;
    uint64_t n = 1000;
    uint64_t *array = NULL;
    for (uint64_t j = 0; j < 100; ++j) {
        p = st_random() * 0.9 + + 0.05;
        mu = 0.0;
        array = (uint64_t*) st_malloc(sizeof(*array) * N);
        for (uint64_t i = 0; i < N; ++i) {
            array[i] = rbinom(n, p);
            mu = runningAverage(mu, array[i], i);
        }
        var = sv(array, N, mu);
        med = median(array, N);
        CuAssertTrue(testCase, array[0] <= array[N-1]);
        CuAssertTrue(testCase, array[0] <= n); // unsigned wraparound    
        CuAssertTrue(testCase, array[N-1] <= n);
        CuAssertTrue(testCase, ((med == floor(n * p)) || 
                                (med == ceil(n * p))));
        CuAssertTrue(testCase, mu > n * (p - 0.01));
        CuAssertTrue(testCase, mu < n * (p + 0.01));
        CuAssertTrue(testCase, var < 1.1 * n * p * (1.0 - p));
        CuAssertTrue(testCase, var > 0.9 * n * p * (1.0 - p));
        free(array);
    }
}
static uint64_t* createSequence(uint64_t n) {
    uint64_t *a = (uint64_t*) st_malloc(sizeof(*a) * n);
    for (uint64_t i = 0; i < n; ++i) {
        a[i] = i;
    }
    return a;
}
static void test_median_0(CuTest *testCase) {
    // tested against values from R
    uint64_t n = 40;
    uint64_t *a = createSequence(n);
    CuAssertTrue(testCase, median(a, n) == 19.5);
    free(a);
    n = 41;
    a = createSequence(n);
    CuAssertTrue(testCase, median(a, n) == 20.0);
    free(a);
    n = 51;
    a = createSequence(n);
    CuAssertTrue(testCase, median(a, n) == 25.0);
    free(a);
    n = 6;
    a = createSequence(n);
    CuAssertTrue(testCase, median(a, n) == 2.5);
    free(a);
    n = 5;
    a = createSequence(n);
    CuAssertTrue(testCase, median(a, n) == 2.0);
    free(a);
}
static void test_sampleVariance_0(CuTest *testCase) {
    // tested against values from R
    uint64_t n = 39;
    uint64_t *a = createSequence(n);
    CuAssertTrue(testCase, sv(a, n, mean(a, n)) == 130.0);
    free(a);
    n = 41;
    a = createSequence(n);
    CuAssertTrue(testCase, sv(a, n, mean(a, n)) == 143.5);
    free(a);
    n = 51;
    a = createSequence(n);
    CuAssertTrue(testCase, sv(a, n, mean(a, n)) == 221.0);
    free(a);
    n = 6;
    a = createSequence(n);
    CuAssertTrue(testCase, sv(a, n, mean(a, n)) == 3.5);
    free(a);
    n = 5;
    a = createSequence(n);
    CuAssertTrue(testCase, sv(a, n, mean(a, n)) == 2.5);
    free(a);
}
CuSuite* comparatorRandom_TestSuite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, test_median_0);
    SUITE_ADD_TEST(suite, test_sampleVariance_0);
    SUITE_ADD_TEST(suite, test_rbinom_zero_0);
    SUITE_ADD_TEST(suite, test_rbinom_zero_1);
    SUITE_ADD_TEST(suite, test_rbinom_one_0);
    SUITE_ADD_TEST(suite, test_rbinom_one_1);
    SUITE_ADD_TEST(suite, test_rbinom_one_2);
    SUITE_ADD_TEST(suite, test_rbinom_one_3);
    SUITE_ADD_TEST(suite, test_rbinom_one_4);
    SUITE_ADD_TEST(suite, test_rbinom_ranges_0);
    SUITE_ADD_TEST(suite, test_rbinom_ranges_1);
    SUITE_ADD_TEST(suite, test_rbinom_ranges_2);
    SUITE_ADD_TEST(suite, test_rbinom_ranges_3);
    SUITE_ADD_TEST(suite, test_rbinom_ranges_4);
    SUITE_ADD_TEST(suite, test_rbinom_clt_0);
    SUITE_ADD_TEST(suite, test_rbinom_clt_1);
    SUITE_ADD_TEST(suite, test_rbinom_distribution_0);
    SUITE_ADD_TEST(suite, test_rbinom_distribution_1);
    SUITE_ADD_TEST(suite, test_rbinom_distribution_2);
    SUITE_ADD_TEST(suite, test_rbinom_distribution_3);
    return suite;
}
