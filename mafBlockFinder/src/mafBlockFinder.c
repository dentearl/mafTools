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
#include "common.h"

int g_verbose_flag = 0;
int g_debug_flag = 0;

void parseOptions(int argc, char **argv, char *seqName, uint32_t *position);
void checkRegion(unsigned lineno, char *fullname, uint32_t pos, uint32_t start, 
                 uint32_t length, uint32_t sourceLength, char strand);
void searchInput(FILE *ifp, char *fullname, unsigned long pos);

void usage(void) {
    fprintf(stderr, "Usage: mafBlockFinder --seq [sequence name (and possibly chr)"
            "] --pos [position to search for] [options] < myFile.maf\n\n"
            "mafBlockFinder is a program that will look through a maf file for a\n"
            "particular sequence name and location. If a match is found the line\n"
            "number and first few fields are returned. If no match is found\n"
            "nothing is returned.\n\n");
    fprintf(stderr, "Options: \n"
            "  -h, --help     show this help message and exit.\n"
            "  -s, --seq      sequence name.chr e.g. `hg18.chr2.'\n"
            "  -p, --pos      position along the chromosome you are searching for.\n"
            "                 Must be a positive number.\n"
            "  -v, --verbose  turns on verbose output.\n");
    exit(EXIT_FAILURE);
}

void parseOptions(int argc, char **argv, char *seqName, uint32_t *position) {
    extern int g_debug_flag;
    extern int g_verbose_flag;
    int c;
    int setSName = 0, setPos = 0;
    int32_t tempPos = 0;
    while (1) {
        static struct option long_options[] = {
            {"debug", no_argument, &g_debug_flag, 1},
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
            if (tempPos < 0) {
                fprintf(stderr, "Error, --pos %d must be nonnegative.\n", tempPos);
                usage();
            }
            *position = tempPos;
            break;
        case 'v':
            g_verbose_flag++;
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
    extern const int kMaxStringLength;
    int32_t n = kMaxStringLength;
    char *buffer = (char *) de_malloc(n);
    char *tkn = NULL;
    char *endPtr = NULL;
    char strand = '\0';
    unsigned lineno = 0;
    uint32_t start, length, sourceLength;
    verbose("searching for %s %u\n", fullname, pos);
    while (de_getline(&buffer, &n, ifp) != -1) {
        lineno++;
        tkn = strtok(buffer, " \t"); // ^a or ^s
        if (tkn == NULL)
            continue;
        if (tkn[0] == '#')
            continue;
        if (tkn[0] != 's')
            continue;
        tkn = strtok(NULL, " \t"); // name field
        if (strncmp(fullname, tkn, strlen(fullname)) == 0) {
            tkn = strtok(NULL, " \t"); // start position
            start = strtoul(tkn, &endPtr, 10);
            tkn = strtok(NULL, " \t"); // length position
            length = strtoul(tkn, &endPtr, 10);
            if (length == UINT32_MAX) {
                fprintf(stderr, "Error on line %u, length (%u) greater than "
                        "or equal to UINT32_MAX (%u).\n",
                        lineno, length, UINT32_MAX);
                exit(EXIT_FAILURE);
            }
            tkn = strtok(NULL, " \t"); // strand
            if ((strcmp(tkn, "+") == 0) || (strcmp(tkn, "-") == 0))
                strcpy(&strand, tkn);
            tkn = strtok(NULL, " \t"); // source length
            sourceLength = strtoul(tkn, &endPtr, 10);
            if (sourceLength == UINT32_MAX) {
                fprintf(stderr, "Error on line %u, source length (%u) greater "
                        "than or equal to UINT32_MAX (%u).\n",
                        lineno, length, UINT32_MAX);
                exit(EXIT_FAILURE);
            }
            tkn = strtok(NULL, " \t"); // sequence field
            checkRegion(lineno, fullname, pos, start, length, sourceLength, strand);
        }
    }
}

int main(int argc, char **argv) {
    extern const int kMaxStringLength;
    char targetName[kMaxStringLength];
    uint32_t targetPos;

    parseOptions(argc, argv,  targetName, &targetPos);
    searchInput(stdin, targetName, targetPos);
   
    return EXIT_SUCCESS;
}
