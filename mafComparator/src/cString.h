/* 
 * Copyright (C) 2009-2013 by 
 * Dent Earl (dearl@soe.ucsc.edu, dentearl@gmail.com)
 * Benedict Paten (benedict@soe.ucsc.edu, benedictpaten@gmail.com)
 * Mark Diekhans (markd@soe.ucsc.edu)
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



#ifndef CSTRING_H_
#define CSTRING_H_

#include "commonC.h"

#include <assert.h>
#include <ctype.h>
#include <string.h>
#include <stdint.h>

/*
 * Comparison function to sort strings alphabetically
 */
int cStr_compare(const void *a, const void *b);

/*
 * Comparison function to sort strings in descending order
 */
int cStr_compareDesc(const void *a, const void *b);

/*
 * In-place substitution to lower-case string
 */
void cStr_lowerCase(char *string);

/*
 * In-place substitution to upper-case string
 */
void cStr_upperCase(char *string);

/* 
 * Check if "string" starts with "query" and ignores case
 *   if "ignorecase" == 1
 */
int cStr_startsWith(char *string, char *query, int ignorecase);

int64_t cStr_getIntLength(int64_t n);

void cStr_reverse(char *s);

void cStr_itoa(int n, char *s);

void cStr_appendChar(char *s, char c);

char *cStr_getStringFromIntArray(int64_t *array, int64_t size, const char sep);

#endif /* CSTRING_H_ */
