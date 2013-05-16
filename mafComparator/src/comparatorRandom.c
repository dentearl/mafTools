/* 
 * Copyright (C) 2009-2013 by 
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
#include "comparatorRandom.h"
#include "sonLib.h"

enum STATE {
    done,
    initialization,
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
    uint64_t n, nsave;
    double p, psave, r, q, fm, p1, xm, xl, xr, c, laml, lamr, p2, p3, p4, u, v, A;
    enum STATE nextStep; // used to traverse the inner subroutines
    enum STATE prevStep; // debugging
    uint64_t result; 
    int64_t y, k, m;
} btpeCalc_t;
static char*stateStrings[] = {"done", "initilaization", "selectRegion", "parallelograms", "leftExponentialTail", "rightExponentialTail", "acceptRejectTest", "acceptRejectRecursive", "acceptRejectSqueeze", "acceptRejectFinal", "correctPvalue"};

static double dmin(double a, double b);
static uint64_t rbinom_smallNaive(const uint64_t n, const double p);
static uint64_t rbinom_smallamlog(const uint64_t n, const double p);
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
    (void) (rbinom_smallamlog);
    (void) (rbinom_btpe);
    (void) (rbinom_smallCdf);
    (void) (stateStrings);
    uint64_t result;
    double q = dmin(p, 1.0 - p);
    if (n > 9223372036854775807) {
        // since we cast down to int64_t we must verify that n will fit
        // 2^63
        fprintf(stderr, "Error in rbinom(), n greater than 9223372036854775807: %"PRIu64"\n", n);
        exit(EXIT_FAILURE);
    }
    if (p == 0.0) {
        return 0;
    }
    if (p < 0.0) {
        fprintf(stderr, "Error in rbinom(), p value less than 0: %f\n", p);
        exit(EXIT_FAILURE);
    }
    if (p >= 1.0) {
        return n;
    }
    if (n * q <= 30) {
        result = rbinom_smallCdf(n, p);
    } else {
        result = rbinom_btpe(n, p);
    }
    assert(result <= n);
    return result;
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
static uint64_t rbinom_smallamlog(const uint64_t n, const double p) {
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
    r = pow(q, (double) n);
    u = st_random();
    while(u > r) {
        u -= r;
        ++x;
        r *= (a / (double)x) - s;
    }
    return x;
}
static btpeCalc_t* emptyBtpeCalc(void) {
    btpeCalc_t *b = (btpeCalc_t*) st_malloc(sizeof(*b));
    b->nsave = 0;
    b->psave = -1.0;
    return b;
}
static double dmin(double a, double b) {
    if (a < b)
        return a;
    else
        return b;
}
static btpeCalc_t* newBtpeCalc(uint64_t n, double p) {
    /* step 0. 
       Set-up constants as functions of n and p. Execute whenever the value of n or p change
     */
    //////////////////////////////////////////////////
    // not multithread safe due to use of static pointer here
    //////////////////////////////////////////////////
    static btpeCalc_t *b;
    if (b == NULL) {
         b = emptyBtpeCalc();
    }
    if ((b->nsave != n) | (b->psave != p)) {
        double a;
        b->n = n;
        b->p = p;
        b->r = dmin(b->p, 1.0 - b->p);
        b->q = 1.0 - b->r;
        b->fm = b->n * b->r + b->r;
        b->m = (int64_t) floor(b->fm);
        b->p1 = floor(2.195 * sqrt(b->n * b->r * b->q) - 4.6 * b->q) + 0.5;
        b->xm = b->m + 0.5;
        b->xl = b->xm - b->p1;
        b->xr = b->xm + b->p1;
        b->c = 0.134 + 20.5 / (15.3 + b->m);
        a = (b->fm - b->xl) / (b->fm - b->xl * b->r);
        b->laml = a * (1.0 + a / 2.0);
        a = (b->xr - b->fm) / (b->xr * b->q);
        b->lamr = a * (1.0 + a / 2.0);
        b->p2 = b->p1 * (1.0 + 2.0 * b->c);
        b->p3 = b->p2 + b->c / b->laml;
        b->p4 = b->p3 + b->c / b->lamr;
    }
    b->prevStep = initialization;
    b->nextStep = selectRegion;
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
        b->y = (int64_t) floor(b->xm - b->p1 * b->v + b->u);
        b->nextStep = correctPvalue;
    }
    b->prevStep = selectRegion;
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
        x = b->xl + (b->u - b->p1) / b->c;
        b->v = b->v * b->c + 1.0 - fabs(b->m - x + 0.5) / b->p1;
        if (b->v > 1.0) {
            b->nextStep = selectRegion;
        } else {
            b->y = (int64_t) floor(x);
            b->nextStep = acceptRejectTest;
        }
    }
    b->prevStep = parallelograms;
}
static void rbinom_btpe_leftExponentialTail_3(btpeCalc_t *b) {
    /* step 3.
       Region 3, left exponential tail.
     */
    if (b->u > b->p3) {
        b->nextStep = rightExponentialTail;
    } else {
        b->y = (int64_t) floor(b->xl + log(b->v) / b->laml);
        if (b->y < 0) {
            b->nextStep = selectRegion;
        } else {
            b->v = b->v * (b->u - b->p2) * b->laml;
            b->nextStep = acceptRejectTest;
        }
    }
    b->prevStep = leftExponentialTail;
}
static void rbinom_btpe_rightExponentialTail_4(btpeCalc_t *b) {
    /* step 4.
       Region 4, right exponential tail.
     */
    b->y = (int64_t) floor(b->xr - log(b->v) / b->lamr);
    if (b->y > b->n) {
        b->nextStep = selectRegion;
    } else {
        b->v = b->v * (b->u - b->p3) * b->lamr;
        b->nextStep = acceptRejectTest;
    }
    b->prevStep = rightExponentialTail;
}
static void rbinom_btpe_acceptRejectTest_5_0(btpeCalc_t *b) {
    /* step 5.0
       Acceptance / Rejection comparison. Test for appropriate method of evaluating f(y).
     */
    b->k = (int64_t) fabs((double) (b->y - b->m));
    if ((b->k > 20) && (b->k < ((b->n * b->r * b->q) / 2.0 - 1.0))) {
        b->nextStep = acceptRejectSqueeze;
    } else {
        b->nextStep = acceptRejectRecursive;
    }
    b->prevStep = acceptRejectTest;
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
        for (i = b->m; i <= b->y; ++i) {
            F *= (a / i - s);
        }
        b->nextStep = acceptRejectSqueeze;
    } else if (b->m > b->y) {
        for (i = b->y; i <= b->m; ++i) {
            F /= (a / i - s);
        }
        b->nextStep = acceptRejectSqueeze;
    } else {
        if (b->v > F) {
            b->nextStep = selectRegion;
        } else {
            b->nextStep = correctPvalue;
        }
    }
    b->prevStep = acceptRejectRecursive;
}
static void rbinom_btpe_acceptRejectSqueeze_5_2(btpeCalc_t *b) {
    /* step 5.2
       Sequeezing. Check the value of ln(v) against upper and lower bound of ln(f(y)).
     */
    double rho, t;
    bool setStep = false;
    rho = (b->k / (b->n * b->r * b->q)) * ((b->k * (b->k / 3.0 + 0.625) + 0.16666666666666666) / 
                                           (b->n * b->r * b->q) + 0.5);
    t = -b->k * b->k / (2.0 * b->n * b->r * b->q);
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
    b->prevStep = acceptRejectSqueeze;
}
static void rbinom_btpe_acceptRejectFinal_5_3(btpeCalc_t *b) {
    /* step 5.3
       Final Acceptance / Rejection Test.
     */
    double x1, x2, f1, f2, z1, z2, w1, w2;
    x1 = b->y + 1.0;
    f1 = b->m + 1.0;
    z1 = b->n + 1.0 - b->m;
    w1 = b->n - b->y + 1.0;
    x2 = x1 * x1;
    f2 = f1 * f1;
    z2 = z1 * z1;
    w2 = w1 * w1;
    if (b->A > (b->xm * log(f1 / x1) 
                + (b->n - b->m + 0.5) * log(z1 / w1) 
                + (b->y - b->m) * log((w1 * b->r) / (x1 * b->q))
                + (13860. - (462. - (132. - (99. - 140. / f2) / f2) / f2) / f2) / f1 / 166320.
                + (13860. - (462. - (132. - (99. - 140. / z2) / z2) / z2) / z2) / z1 / 166320.
                + (13860. - (462. - (132. - (99. - 140. / x2) / x2) / x2) / x2) / x1 / 166320.
                + (13860. - (462. - (132. - (99. - 140. / w2) / w2) / w2) / w2) / w1 / 166320.)) {
        b->nextStep = selectRegion;
    } else {
        b->nextStep = correctPvalue;
    }
    b->prevStep = acceptRejectFinal;
}
static void rbinom_btpe_correctPvalue_6(btpeCalc_t *b) {
    /* step 6.
       Check if original probability p is greater than 1/2.
       If so, return n - y;
     */
    if (b->p > 0.5) {
        b->y = b->n - b->y;
    }
    b->nextStep = done;
    b->prevStep = correctPvalue;
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
        case initialization:
            break;
        }
    }
    assert(b->y <= b->n); // unsigned wrap around
    uint64_t result = (uint64_t) b->y;
    // free(b); // don't do this now that b is a local static pointer.
    return result;
}
