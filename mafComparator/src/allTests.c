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

#include <stdio.h>
#include <stdlib.h>
#include "CuTest.h"
#include "comparatorAPI.h"
#include "test.comparatorAPI.h"
#include "test.comparatorRandom.h"

CuSuite* comparatorAPI_TestSuite(void);
CuSuite* comparatorRandom_TestSuite(void);

int comparator_RunAllTests(void) {
    CuString *output = CuStringNew();
    CuSuite *suite = CuSuiteNew();
    CuSuite *comparatorAPI_s = comparatorAPI_TestSuite();
    CuSuite *comparatorRandom_s = comparatorRandom_TestSuite();
    CuSuiteAddSuite(suite, comparatorAPI_s);
    CuSuiteAddSuite(suite, comparatorRandom_s);
    CuSuiteRun(suite);
    CuSuiteSummary(suite, output);
    CuSuiteDetails(suite, output);
    printf("%s\n", output->buffer);
    CuStringDelete(output);
    int status = (suite->failCount > 0);
    free(comparatorAPI_s);
    free(comparatorRandom_s);
    CuSuiteDelete(suite);
    return status;
}
int main(void) {
    return comparator_RunAllTests();
}
