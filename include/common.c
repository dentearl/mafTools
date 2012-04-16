#include <assert.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "common.h"

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
void processHeader(void) {
    // Read in the header and spit it back out
    FILE *ifp = stdin;
    int32_t n = kMaxStringLength;
    char *line = (char*) de_malloc(n);
    int status = de_getline(&line, &n, ifp);
    if (status == -1) {
        fprintf(stderr, "Error, empty file\n");
        exit(EXIT_FAILURE);
    }
    if (strncmp(line, "##maf", 5) != 0) {
        fprintf(stderr, "line: %s\n", line);
        fprintf(stderr, "Error, bad maf format. File should start with ##maf\n");
        exit(EXIT_FAILURE);
    }
    printf("%s\n", line);
    while (*line != 0 && status != -1) {
        status = de_getline(&line, &n, ifp);
        printf("%s\n", line);
    }
}
void failBadFormat(void) {
    fprintf(stderr, "The maf sequence lines are incorrectly formatted, exiting\n");
    exit(EXIT_FAILURE);
}
