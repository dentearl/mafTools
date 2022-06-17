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

extern int g_verbose_flag;
extern int g_debug_flag;
extern const int kMaxStringLength;
extern const int kMaxMessageLength;
extern const int kMaxSeqName;

void de_verbose(char const *fmt, ...);
void de_debug(char const *fmt, ...);
void* de_malloc(size_t n);
int64_t de_getline(char **s, int64_t *n, FILE *f);
FILE* de_fopen(const char *s, char const *mode);
char* de_strdup(const char *s);
char* de_strndup(const char *s, size_t n);
void failBadFormat(void);
void usageMessage(char shortopt, const char *name, const char *description);
char* stringReplace(const char *string, const char a, const char b);
int minint(int a, int b);
char* de_strtok(char **s, char t);
unsigned countChar(char *s, const char c);
char** extractSubStrings(char *nameList, unsigned n, const char delineator);

#endif // COMMON_H_
