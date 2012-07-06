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
#ifndef COMMON_H_
#define COMMON_H_
#include <stdio.h>
#include <stdint.h>

int g_verbose_flag;
int g_debug_flag;
const int kMaxStringLength;
const int kMaxMessageLength;
const int kMaxSeqName;

void de_verbose(char const *fmt, ...);
void de_debug(char const *fmt, ...);
void* de_malloc(int32_t n);
int32_t de_getline(char **s, int32_t *n, FILE *f);
FILE* de_fopen(const char *s, char const *mode);
char* de_strdup(const char *s);
void failBadFormat(void);
void usageMessage(char shortopt, const char *name, const char *description);
char* stringCommasToSpaces(const char *string);

#endif // COMMON_H_
