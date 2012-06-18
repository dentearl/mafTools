/* 
 * Copyright (C) 2009-2012 by 
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
#include <stdbool.h>
#include <stdlib.h>
#include <math.h>
#include "ComparatorRandom.h"
#include "sonLib.h"

enum STATE {
    done,
    selectRegion,
    parallelograms,
    leftExponentialTail,
    rightExponentialTail,
    acceptRejectTest,
    acceptRejectRecursive,
    acceptRejectSqueeze,
    acceptRejectFinal,
    correctPvalue
};

typedef struct btpeCalc {
    uint64_t n;
    double p, r, q, f_M, M, p1, x_M, x_L, x_R, c, l_L, l_R, p2, p3, p4, u, v, y, m, A, k;
    enum STATE nextStep; // used to traverse the inner subroutines
    uint64_t result;   // the ultimate result, a binomal draw from binom(n, p);
} btpeCalc_t;

static uint64_t rbinom_smallNaive(const uint64_t n, const double p);
static uint64_t rbinom_smallLog(const uint64_t n, const double p);
static uint64_t rbinom_smallCdf(const uint64_t n, const double p);
static btpeCalc_t* newBtpeCalc(const uint64_t n, const double p);
static void rbinom_btpe_selectRegion_1(btpeCalc_t *b);
static void rbinom_btpe_parallelograms_2(btpeCalc_t *b);
static void rbinom_btpe_leftExponentialTail_3(btpeCalc_t *b);
static void rbinom_btpe_rightExponentialTail_4(btpeCalc_t *b);
static void rbinom_btpe_acceptRejectTest_5_0(btpeCalc_t *b);
static void rbinom_btpe_acceptRejectRecursive_5_1(btpeCalc_t *b);
static void rbinom_btpe_acceptRejectSqueeze_5_2(btpeCalc_t *b);
static void rbinom_btpe_acceptRejectFinal_5_3(btpeCalc_t *b);
static void rbinom_btpe_correctPvalue_6(btpeCalc_t *b);
static uint64_t rbinom_btpe(uint64_t n, double p);
    // BTPE (Binomial, Trinagle, Parallelogram, Exponential)
    // Kachitvichyanukul, Voratas and Schmeiser, Bruce W. (1988)
    // Binomial Random Variate Generation, Communications of the ACM, 31(2): 216-222

uint64_t rbinom(const uint64_t n, const double p) {
    // make a draw from a binomial distribution with parameters n and p
    (void) (rbinom_smallNaive);
    (void) (rbinom_smallLog);
    (void) (rbinom_smallCdf);
    if (n < 20)
        return rbinom_smallCdf(n, p);
    else
        return rbinom_btpe(n, p);
}
static uint64_t rbinom_smallNaive(const uint64_t n, const double p) {
    // speed proportional to n
    uint64_t x = 0;
    for (uint64_t i = 0; i < n; ++i) {
        if (st_random() <= p) {
            ++x;
        }
    }
    return x;
}
static uint64_t rbinom_smallLog(const uint64_t n, const double p) {
    // Binomial via Geometric, a la Devroye
    // speed proportional to n * p
    uint64_t x = 0, y = 0;
    double c = log(1.0 - p), u = 0.0;
    if (c == 0.0) {
        return x;
    }
    while (true) {
        u = st_random();
        while (u == 0.0) {
            u = st_random();
        }
        y += floor(log(u) / c) + 1;
        if (y <= n) {
            // note that this is incorrectly (y < n) in K&S 1988.
            ++x;
        } else {
            break;
        }
    }
    return x;
}
static uint64_t rbinom_smallCdf(const uint64_t n, const double p) {
    // Binomial via inverse
    // speed proportional to n * p
    uint64_t x = 0;
    double q, s, a, r, u;
    q = 1.0 - p;
    s = p / q;
    a = (n + 1.0) * s;
    r = pow(q, n);
    u = st_random();
    while(u > r) {
        u -= r;
        ++x;
        r *= (a / (double)x) - s;
    }
    return x;
}
static btpeCalc_t* newBtpeCalc(uint64_t n, double p) {
    /* step 0. 
       Set-up constants as functions of n and p. Execute whenever the value of n or p change
     */
    btpeCalc_t *b = (btpeCalc_t*) st_malloc(sizeof(*b));
    double a;
    b->n = n;
    b->p = p;
    b->result = 0;
    b->nextStep = selectRegion;
    if (b->p < 1.0 - b->p) {
        b->r = b->p;
    } else {
        b->r = 1.0 - b->p;
    }
    b->q = 1.0 - b->r;
    b->f_M = n * b->r + b->r;
    b->M = floor(b->f_M);
    b->p1 = floor(2.195 * sqrt(b->n * b->r * b->q) - 4.6 * b->q);
    b->x_M = b->M + 0.5;
    b->x_L = b->x_M - b->p1;
    b->x_R = b->x_M + b->p1;
    b->c = 0.134 + 20.5 / (15.3 + b->M);
    a = (b->f_M - b->x_L) / (b->f_M - b->x_L * b->r);
    b->l_L = a * (a + a / 2.0);
    a = (b->x_R - b->f_M) / (b->x_R * b->q);
    b->l_R = a * (1.0 + a / 2.0);
    b->p2 = b->p1 * (1.0 + 2.0 * b->c);
    b->p3 = b->p2 + b->c / b->l_L;
    b->p4 = b->p3 + b->c / b->l_R;
    b->m = b->n * b->p + b->p;
    return b;
}
static void rbinom_btpe_selectRegion_1(btpeCalc_t *b) {
    /* step 1.
       Generate u ~ U(0, p4) for selecting the region. If region 1 is selected, generate 
       a triangularly distributed variate.
     */
    b->u = st_random() * b->p4;
    b->v = st_random();
    if (b->u > b->p1) {
        b->nextStep = parallelograms;
    } else {
        b->y = floor(b->x_M - b->p1 * b->v + b->u);
        b->nextStep = correctPvalue;
    }
}
static void rbinom_btpe_parallelograms_2(btpeCalc_t *b) {
    /* step 2.
       Region 2, parallelograms. Check if region 2 is used.
       If so, generate y ~ U(x_L - 0.5, x_R - 0.5)
     */
    double x;
    if (b->u > b->p2) {
        b->nextStep = leftExponentialTail;
    } else {
        x = b->x_L + (b->u - b->p1) / b->c;
        b->v = b->v * b->c + 1.0 - fabs(b->M - x + 0.5) / b->p1;
        if (b->v > 1.0) {
            b->nextStep = selectRegion;
        } else {
            b->y = floor(x);
            b->nextStep = acceptRejectTest;
        }
    }
}
static void rbinom_btpe_leftExponentialTail_3(btpeCalc_t *b) {
    /* step 3.
       Region 3, left exponential tail.
     */
    if (b->u > b->p3) {
        b->nextStep = rightExponentialTail;
    } else {
        b->y = floor(b->x_L + log(b->v) / b->l_L);
        if (b->y < 0) {
            b->nextStep = selectRegion;
        } else {
            b->v = b->v * (b->u - b->p2) * b->l_L;
            b->nextStep = acceptRejectTest;
        }
    }
}
static void rbinom_btpe_rightExponentialTail_4(btpeCalc_t *b) {
    /* step 4.
       Region 4, right exponential tail.
     */
    b->y = floor(b->x_R - log(b->v) / b->l_R);
    if (b->y > b->n) {
        b->nextStep = selectRegion;
    } else {
        b->v = b->v * (b->u - b->p3) * b->l_R;
        b->nextStep = acceptRejectTest;
    }
}
static void rbinom_btpe_acceptRejectTest_5_0(btpeCalc_t *b) {
    /* step 5.0
       Acceptance / Rejection comparison. Test for appropriate method of evaluating f(y).
     */
    b->k = fabs(b->y - b->M);
    if ((b->k > 20) && (b->k < ((b->n * b->r * b->q) / 2.0 - 1.0))) {
        b->nextStep = acceptRejectSqueeze;
    } else {
        b->nextStep = acceptRejectRecursive;
    }
}
static void rbinom_btpe_acceptRejectRecursive_5_1(btpeCalc_t *b) {
    /* step 5.1
       Evaluate f(y) via the recursive relationship f(y) = f(y - 1)(a/x - s).
       Start the search from the mode.
     */
    double F, s, i, a;
    s = b->r / b->q;
    a = s * (b->n + 1);
    F = 1.0;
    if (b->m < b->y) {
        for (i = b->M + 1; i <= b-> y; ++i) {
            F = F * (a / i - s);
        }
    } else {
        if (b->M > b->y) {
            for (i = b->y + 1; i <= b->M; ++i) {
                F = F / ( a / i - s);
            }
        }
    }
    if (b->v > F) {
        b->nextStep = selectRegion;
    } else {
        b->nextStep = correctPvalue;
    }
}
static void rbinom_btpe_acceptRejectSqueeze_5_2(btpeCalc_t *b) {
    /* step 5.2
       Sequeezing. Check the value of ln(v) against upper and lower bound of ln(f(y)).
     */
    double rho, t;
    bool setStep = false;
    rho = (b->k / (b->n * b->r * b->q)) * ((b->k * (b->k / 3.0 + 0.625) + 0.16666666666666666) / 
                                           (b->n * b->r * b->q) + 0.5);
    t = -(b->k * b->k) / (2.0 * b->n * b->r * b->q);
    b->A = log(b->v);
    if (b->A < t - rho) {
        b->nextStep = correctPvalue;
        setStep = true;
    }
    if (b->A > t + rho) {
        b->nextStep = selectRegion;
        setStep = true;
    }
    if (!setStep) {
        b->nextStep = acceptRejectFinal;
    }
}
static void rbinom_btpe_acceptRejectFinal_5_3(btpeCalc_t *b) {
    /* step 5.3
       Final Acceptance / Rejection Test.
     */
    double x1, x2, f1, f2, z1, z2, w1, w2;
    x1 = b->y + 1.0;
    f1 = b->M + 1.0;
    z1 = b->n + 1.0 - b->M;
    w1 = b->n - b->y + 1.0;
    x2 = x1 * x1;
    f2 = f1 * f1;
    z2 = z1 * z1;
    w2 = w1 * w1;
    if (b->A > (b->x_M * log(f1 / x1) + (b->n - b->M + 0.5) * log(z1 / w1) 
                + (b->y - b->M) * log(w1 * b->r / x1 * b->q)
                + (13860.0 - (462.0 - (132.0 - (99.0 - 140.0 / f2) / f2) / f2) / f2) / f1 / 166320.0
                + (13860.0 - (462.0 - (132.0 - (99.0 - 140.0 / z2) / z2) / z2) / z2) / z1 / 166320.0
                + (13860.0 - (462.0 - (132.0 - (99.0 - 140.0 / x2) / x2) / x2) / x2) / x1 / 166320.0
                + (13860.0 - (462.0 - (132.0 - (99.0 - 140.0 / w2) / w2) / w2) / w2) / w1 / 166320.0)) {
        b->nextStep = selectRegion;
    } else {
        b->nextStep = correctPvalue;
    }
}
static void rbinom_btpe_correctPvalue_6(btpeCalc_t *b) {
    /* step 6.
       Check if original probability p is greater than 1/2.
       If so, return n - y;
     */
    if (b->p > 0.5) {
        b->result = b->n - b->y;
    } else {
        b->result = b->y;
    }
    b->nextStep = done;
}
static uint64_t rbinom_btpe(uint64_t n, double p) {
    // BTPE (Binomial, Trinagle, Parallelogram, Exponential)
    // Kachitvichyanukul, Voratas and Schmeiser, Bruce W. (1988)
    // Binomial Random Variate Generation, Communications of the ACM, 31(2): 216-222
    btpeCalc_t *b = newBtpeCalc(n, p);
    while (b->nextStep != done) {
        switch (b->nextStep) {
        case selectRegion:
            rbinom_btpe_selectRegion_1(b);
            break;
        case parallelograms:
            rbinom_btpe_parallelograms_2(b);
            break;
        case leftExponentialTail:
            rbinom_btpe_leftExponentialTail_3(b);
            break;
        case rightExponentialTail:
            rbinom_btpe_rightExponentialTail_4(b);
            break;
        case acceptRejectTest:
            rbinom_btpe_acceptRejectTest_5_0(b);
            break;
        case acceptRejectRecursive:
            rbinom_btpe_acceptRejectRecursive_5_1(b);
            break;
        case acceptRejectSqueeze:
            rbinom_btpe_acceptRejectSqueeze_5_2(b);
            break;
        case acceptRejectFinal:
            rbinom_btpe_acceptRejectFinal_5_3(b);
            break;
        case correctPvalue:
            rbinom_btpe_correctPvalue_6(b);
            break;
        case done:
            break;
        }
    }
    uint64_t result = b->result;
    free(b);
    return result;
}
