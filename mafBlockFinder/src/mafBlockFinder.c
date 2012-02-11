/* 
 * Copyright (C) 2011-2012 by 
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
#include <getopt.h>
#include <limits.h>
#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#define d_MAX_STRING_LENGTH 2048
#define d_MAX_MESSAGE_LENGTH 1024

static int verbose_flag;
static int debug_flag;

void usage(void);
void verbose(char const *fmt, ...);
void debug(char const *fmt, ...);
void message(char const *type, char const *fmt, ...);
void parseOptions(int argc, char **argv, char *seqName, uint32_t *position);
void checkRegion(unsigned lineno, char *fullname, uint32_t pos, uint32_t start, 
                 uint32_t length, uint32_t sourceLength, char strand);
void searchInput(FILE *ifp, char *fullname, unsigned long pos);
void* de_malloc(size_t n);
int32_t de_getline(char **s, int32_t *n, FILE *f);

void* de_malloc(size_t n) {
    void *i;
    i = malloc(n);
    if (i == NULL) {
        fprintf(stderr, "malloc failed on a request for %zu bytes\n", n);
        exit(1);
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
    extern int verbose_flag;
    char str[d_MAX_MESSAGE_LENGTH];
    va_list args;
    va_start(args, fmt);
    if (verbose_flag) {
        int n = vsprintf(str, fmt, args);
        if (n >= d_MAX_MESSAGE_LENGTH) {
            fprintf(stderr, "Error, failure in verbose(), (n = %d) > "
                    "(d_MAX_MESSAGE_LENGTH %d)\n", n, d_MAX_MESSAGE_LENGTH);
            exit(1);
        }
        message("Verbose", str, args);
    }
    va_end(args);
}

void debug(char const *fmt, ...) {
    extern int debug_flag;
    char str[d_MAX_MESSAGE_LENGTH];
    va_list args;
    va_start(args, fmt);
    if (debug_flag) {
        int n = vsprintf(str, fmt, args);
        if (n >= d_MAX_MESSAGE_LENGTH) {
            fprintf(stderr, "Error, failure in debug(), (n = %d) > "
                    "(d_MAX_MESSAGE_LENGTH %d)\n", n, d_MAX_MESSAGE_LENGTH);
            exit(1);
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

void usage(void) {
    fprintf(stderr, "Usage: mafBlockFinder --seq [sequence name (and possibly chr)"
            "] --pos [position to search for] [options] < myFile.maf\n\n");
    fprintf(stderr, "Options: \n"
            "  -h, --help     show this help message and exit.\n"
            "  -s, --seq      sequence name.chr e.g. `hg18.chr2.'\n"
            "  -p, --pos      position along the chromosome you are searching for.\n"
            "                 Must be a positive number.\n"
            "  -v, --verbose  turns on verbose output.\n");
    exit(1);
}

void parseOptions(int argc, char **argv, char *seqName, uint32_t *position) {
    int c;
    int setSName = 0, setPos = 0;
    int32_t tempPos = 0;
    while (1) {
        static struct option long_options[] = {
            {"debug", no_argument, &debug_flag, 1},
            {"verbose", no_argument, 0, 'v'},
            {"help", no_argument, 0, 'h'},
            {"seq",  required_argument, 0, 's'},
            {"pos",  required_argument, 0, 'p'},
            {0, 0, 0, 0}
        };
        int option_index = 0;
        c = getopt_long(argc, argv, "n:c:p:v",
                        long_options, &option_index);
        if (c == -1) {
            break;
        }
        switch (c) {
        case 0:
            break;
        case 's':
            setSName = 1;
            sscanf(optarg, "%s", seqName);
            break;
        case 'p':
            setPos = 1;
            tempPos = strtoll(optarg, NULL, 10);
            if (tempPos < 1) {
                fprintf(stderr, "Error, --pos %d must be positive\n", tempPos);
                usage();
            }
            *position = tempPos;
            break;
        case 'v':
            verbose_flag++;
            break;
        case 'h':
        case '?':
            usage();
            break;
        default:
            abort();
        }
    }
    if (!(setSName && setPos)) {
        fprintf(stderr, "specify --seq --position\n");
        usage();
    }
    // Check there's nothing left over on the command line 
    if (optind < argc) {
        char errorString[30] = "Unexpected arguments:";
        while (optind < argc) {
            strcat(errorString, " ");
            strcat(errorString, argv[optind++]);
        }
        fprintf(stderr, "%s\n", errorString);
        usage();
    }
}

void checkRegion(unsigned lineno, char *fullname, uint32_t pos, uint32_t start, 
                 uint32_t length, uint32_t sourceLength, char strand) {
    // check to see if pos is in this block
    static unsigned long absStart, absEnd;
    if (strand == '-') {
        absStart =  sourceLength - (start + length);
        absEnd = sourceLength - start - 1;
    } else {
        absStart = start;
        absEnd = start + length - 1;
    }
    if ((absStart <= pos) && (absEnd >= pos))
        printf("%u: s %s %u %u %c %u ...\n", lineno, fullname, start, 
               length, strand, sourceLength);
}

void searchInput(FILE *ifp, char *fullname, unsigned long pos) {
    int32_t n = d_MAX_STRING_LENGTH;
    char *buffer = (char *) de_malloc(n);
    const char search[] = " \t";
    char *tkn = NULL;
    char *endPtr = NULL;
    char strand = '\0';
    unsigned lineno = 0;
    uint32_t start, length, sourceLength;
    verbose("searching for %s %u\n", fullname, pos);
    while (de_getline(&buffer, &n, ifp) != -1) {
        lineno++;
        tkn = strtok(buffer, search); // ^a or ^s
        if (tkn == NULL)
            continue;
        if (tkn[0] == '#')
            continue;
        if (tkn[0] != 's')
            continue;
        tkn = strtok(NULL, search); // name field
        if (strncmp(fullname, tkn, strlen(fullname)) == 0) {
            tkn = strtok(NULL, search); // start position
            start = strtoul(tkn, &endPtr, 10);
            tkn = strtok(NULL, search); // length position
            length = strtoul(tkn, &endPtr, 10);
            if (length == UINT32_MAX) {
                fprintf(stderr, "Error on line %u, length (%u) greater than "
                        "or equal to UINT32_MAX.\n",
                        lineno, UINT32_MAX);
                exit(1);
            }
            tkn = strtok(NULL, search); // strand
            if ((strcmp(tkn, "+") == 0) || (strcmp(tkn, "-") == 0))
                strcpy(&strand, tkn);
            tkn = strtok(NULL, search); // source length
            sourceLength = strtoul(tkn, &endPtr, 10);
            if (sourceLength == UINT32_MAX) {
                fprintf(stderr, "Error on line %u, source length (%u) greater "
                        "than or equal to UINT32_MAX.\n",
                        lineno, UINT32_MAX);
                exit(1);
            }
            tkn = strtok(NULL, search); // sequence field
            checkRegion(lineno, fullname, pos, start, length, sourceLength, strand);
        }
    }
}

int main(int argc, char **argv) {
    char targetName[d_MAX_STRING_LENGTH];
    uint32_t targetPos;

    parseOptions(argc, argv,  targetName, &targetPos);
    searchInput(stdin, targetName, targetPos);
   
    return 0;
}
