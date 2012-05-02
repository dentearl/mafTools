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
#include <assert.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

void verbose(char const *fmt, ...);
void debug(char const *fmt, ...);
void message(char const *type, char const *fmt, ...);
char* processHeader(void);
void failBadFormat(void);

void* de_malloc(size_t n);
int32_t de_getline(char **s, int32_t *n, FILE *f);

const int kMaxStringLength = 2048;
const int kMaxMessageLength = 1024;

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
void verbose(char const *fmt, ...) {
    extern int g_verbose_flag;
    char str[kMaxMessageLength];
    va_list args;
    va_start(args, fmt);
    if (g_verbose_flag) {
        int n = vsprintf(str, fmt, args);
        if (n >= kMaxMessageLength) {
            fprintf(stderr, "Error, failure in verbose(), (n = %d) > "
                    "(kMaxMessageLength %d)\n", n, kMaxMessageLength);
            exit(EXIT_FAILURE);
        }
        message("Verbose", str, args);
    }
    va_end(args);
}
void debug(char const *fmt, ...) {
    extern int g_debug_flag;
    char str[kMaxMessageLength];
    va_list args;
    va_start(args, fmt);
    if (g_debug_flag) {
        int n = vsprintf(str, fmt, args);
        if (n >= kMaxMessageLength) {
            fprintf(stderr, "Error, failure in debug(), (n = %d) > "
                    "(kMaxMessageLength %d)\n", n, kMaxMessageLength);
            exit(EXIT_FAILURE);
        }
        message("Debug", str, args);
    }
    va_end(args);
}
void message(char const *type, char const *fmt, ...) {
    va_list args;
    va_start(args, fmt);
    fprintf(stderr, "%s: ", type);
    vfprintf(stderr, fmt, args);
    va_end(args);
}
char* processHeader(void) {
    // Read in the header and spit it back out
    // IF the header does not end in an empty line,
    // but instead transitions directly to an alignment line, 
    // the return value is a pointer to that alignment line. 
    // ELSE the return value is null.
    FILE *ifp = stdin;
    int32_t n = kMaxStringLength;
    char *line = (char*) de_malloc(n);
    int status = de_getline(&line, &n, ifp);
    if (status == -1) {
        fprintf(stderr, "Error, empty file\n");
        free(line);
        exit(EXIT_FAILURE);
    }
    if (strncmp(line, "##maf", 5) != 0) {
        fprintf(stderr, "line: %s\n", line);
        fprintf(stderr, "Error, bad maf format. File should start with ##maf\n");
        free(line);
        exit(EXIT_FAILURE);
    }
    printf("%s\n", line);
    while (*line != 0 && status != -1) {
        status = de_getline(&line, &n, ifp);
        if (line[0] == 'a')
            return line;
        printf("%s\n", line);
    }
    free(line);
    return 0;
}
void failBadFormat(void) {
    fprintf(stderr, "The maf sequence lines are incorrectly formatted, exiting\n");
    exit(EXIT_FAILURE);
}

#endif // COMMON_H_
