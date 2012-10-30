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
#include <inttypes.h>
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
        fprintf(stderr, "(de_) malloc failed on a request for %zu bytes\n", n);
        exit(EXIT_FAILURE);
    }
    return i;
}
int32_t de_getline(char **s, int32_t *n, FILE *f) {
    register int32_t nMinus1 = ((*n) - 1), i = 0;
    char *s2 = *s;
    while (1) {
        register int32_t ch = (char) getc(f);
        if (ch == '\r') {
            ch = getc(f);
        }
        if (i == nMinus1) {
            *n = 2 * (*n) + 1;
            *s = realloc(*s, (*n + 1) * sizeof(char));
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
FILE* de_fopen(const char *filename, char const *mode) {
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
char* de_strndup(const char *s, size_t n) {
    char *copy = de_malloc(n + 1);
    strncpy(copy, s, n);
    copy[n] = '\0';
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
    de_message("Verbose", str);
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
    de_message("Debug", str);
    va_end(args);
}
void failBadFormat(void) {
    fprintf(stderr, "The maf sequence lines are incorrectly formatted, exiting\n");
    exit(EXIT_FAILURE);
}
int minint(int a, int b) {
    if (a < b) {
        return a;
    } else {
        return b;
    }
}
void usageMessage(char shortopt, const char *name, const char *description) {
    // pretty print a usage() message a la:
    //     -o, --option     blah blah blah, explanation blah blah.
    //     --option2        blah blah blah, explanation blah blah.
    // can accept a third argument, a CHAR for the short option name.
    int lineLength = 70;
    int indent = 2;
    int shortlen = 0;
    int length = strlen(name);
    int linePos = length + indent + 5;
    int strPos = 0;
    for (int i = 0; i < indent; ++i) {
        fprintf(stderr, " ");
    }
    if (shortopt != '\0') {
        fprintf(stderr, "-%c, ", shortopt);
        shortlen = 4;
    }
    fprintf(stderr, "--%s", name);
    if (indent + length + 2 + shortlen < 25) {
        for (int i = indent + length + 2 + shortlen; i < 25; ++i) {
            fprintf(stderr, " ");
        }
    } else {
        fprintf(stderr, "\n");
        for (int i = 0; i < 25; ++i) {
            fprintf(stderr, " ");
        }
    }
    while (strPos < (int)strlen(description)) {
        fprintf(stderr, "%c", description[strPos++]);
        ++linePos;
        if (linePos >= lineLength) {
            if (description[strPos] == ' ') {
                linePos = 0;
                fprintf(stderr, "\n");
                for (int i = 0; i < 25; ++i) {
                    fprintf(stderr, " ");
                }
                while (description[strPos] != '\0' && description[strPos] == ' ') {
                    ++strPos;
                }
            }
        }
    }
    fprintf(stderr, "\n");
}
char* stringReplace(const char *string, const char a, const char b) {
    /* swap all chars equal to a, for b.
     */
    char *s = de_strdup(string);
    int i = 0;
    for (i = 0; i < (int)strlen(s); i++) {
        if (s[i] == a)
            s[i] = b;
    }
    return s;
}
char* de_strtok(char **s, char t) {
    // takes a string and returns a COPY that excludes all content that comes
    // BEFORE the token. Meant to be iterated. Destructive to input.
    if (*s == NULL) {
        return NULL;
    }
    if ((**s) == t) {
        (*s)++;
    }
    char *start = *s;
    while (**s != '\0' && (**s) != t) {
        (*s)++;
    }
    if ((*s) - start > 0) {
        char *result = memcpy(de_malloc(((*s) - start + 1) * sizeof(char)), start, ((*s) - start) * sizeof(char));
        result[(*s) - start] = '\0';
        return result;
    }
    return NULL;
}
unsigned countChar(char *s, const char c) {
    // counts the number of a specified character in a string.
    unsigned i, m, n = 0;
    m = strlen(s);
    if (m == 0) {
        return 0;
    }
    for (i = 0; i < m; ++i) {
        if (s[i] == c) {
            ++n;
        }
    }
    return n;
}
char** extractSubStrings(char *nameList, unsigned n, const char delineator) {
    // n is the number of names in the name list
    // create a char array to hold each element in the *namelist comma separated list
    if (n == 0) {
        return NULL;
    }
    char t[2];
    t[0] = delineator;
    t[1] = '\0';
    char **mat = (char**) de_malloc(sizeof(char*) * n);
    unsigned index = 0;
    char *tkn = NULL;
    char *copy = de_strdup(nameList);
    tkn = strtok(copy, t);
    while (tkn != NULL) {
        mat[index] = (char*) de_malloc(sizeof(char) * (strlen(tkn) + 1));
        strcpy(mat[index++], tkn);
        tkn = strtok(NULL, t);
    }
    free(copy);
    copy = NULL;
    return mat;
}
