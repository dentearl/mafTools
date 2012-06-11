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
#include <ctype.h>
#include <stdarg.h>
#include <errno.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include "CuTest.h"
#include "common.h"

int g_verbose_flag = 0;
int g_debug_flag = 0;
const int kMaxStringLength = 2048;
const int kMaxMessageLength = 1024;
const int kMaxSeqName = 1 << 9;

void* de_malloc(size_t n) {
    void *i;
    i = malloc(n);
    if (i == NULL) {
        fprintf(stderr, "malloc failed on a request for %zu bytes\n", n);
        exit(EXIT_FAILURE);
    }
    return i;
}
int32_t de_getline(char **s, int32_t *n, FILE *f) {
    register int32_t nMinus1 = ((*n) - 1), i = 0;
    char *s2 = *s;
    while (1) {
        register int32_t ch = (char) getc(f);
        if (ch == '\r')
            ch = getc(f);
        if (i == nMinus1) {
            *n = 2 * (*n) + 1;
            *s = realloc(*s, (*n + 1));
            assert(*s != NULL);
            s2 = *s + i;
            nMinus1 = ((*n) - 1);
        }
        if ((ch == '\n') || (ch == EOF)) {
            *s2 = '\0';
            return(feof(f) ? -1 : i);
        } else {
            *s2 = ch;
            s2++;
        }
        ++i;
    }
}
FILE* de_open(const char *filename, char const *mode) {
    FILE *f = fopen(filename, mode);
    if (f == NULL) {
        if (errno == ENOENT) {
            fprintf(stderr, "ERROR, file %s does not exist.\n", filename);
        } else {
            fprintf(stderr, "ERROR, unable to open file %s for mode \"%s\"\n", filename, mode);
        }
        exit(EXIT_FAILURE);
    }
    return f;
}
char* de_strdup(const char *s) {
    size_t n = strlen(s) + 1;
    char *copy = de_malloc(n);
    strcpy(copy, s);
    return copy;
}
static void de_message(char const *type, char const *fmt, ...) {
    va_list args;
    va_start(args, fmt);
    fprintf(stderr, "%s: ", type);
    vfprintf(stderr, fmt, args);
    va_end(args);
}
void de_verbose(char const *fmt, ...) {
    if (!g_verbose_flag) {
        return;
    }
    char str[kMaxMessageLength];
    va_list args;
    va_start(args, fmt);
    int n = vsprintf(str, fmt, args);
    if (n >= kMaxMessageLength) {
        fprintf(stderr, "Error, failure in verbose(), (n = %d) > "
                "(kMaxMessageLength %d)\n", n, kMaxMessageLength);
        exit(EXIT_FAILURE);
    }
    de_message("Verbose", str, args);
    va_end(args);
}
void de_debug(char const *fmt, ...) {
    if (!g_debug_flag) {
        return;
    }
    char str[kMaxMessageLength];
    va_list args;
    va_start(args, fmt);
    int n = vsprintf(str, fmt, args);
    if (n >= kMaxMessageLength) {
        fprintf(stderr, "Error, failure in debug(), (n = %d) > "
                "(kMaxMessageLength %d)\n", n, kMaxMessageLength);
        exit(EXIT_FAILURE);
    }
    de_message("Debug", str, args);
    va_end(args);
}
void failBadFormat(void) {
    fprintf(stderr, "The maf sequence lines are incorrectly formatted, exiting\n");
    exit(EXIT_FAILURE);
}
