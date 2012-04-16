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
#include <getopt.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "common.h"

const int kMaxSeqName = 1 << 9;

int g_verbose_flag = 0;
int g_debug_flag = 0;

typedef struct mafBlock {
    struct mafBlock *next;
    char *line;
} mafBlock_t;

void usage(void);

void parseOptions(int argc, char **argv, char *seqName, uint32_t *start, uint32_t *stop) {
    extern int g_debug_flag;
    extern int g_verbose_flag;
    int c;
    bool setSName = false, setStart = false, setStop = false;
    int32_t value = 0;
    while (1) {
        static struct option long_options[] = {
            {"debug", no_argument, 0, 'd'},
            {"verbose", no_argument, 0, 'v'},
            {"help", no_argument, 0, 'h'},
            {"seq",  required_argument, 0, 's'},
            {"start",  required_argument, 0, 0},
            {"stop",  required_argument, 0, 0},
            {0, 0, 0, 0}
        };
        int option_index = 0;
        c = getopt_long(argc, argv, "n:c:p:v",
                        long_options, &option_index);
        if (c == -1)
            break;
        switch (c) {
        case 0:
            if (strcmp("start", long_options[option_index].name) == 0) {
                value = strtoll(optarg, NULL, 10);
                if (value < 0) {
                    fprintf(stderr, "Error, --start %d must be nonnegative.\n", value);
                    usage();
                }
                *start = value;
                setStart = true;
            } else if (strcmp("stop", long_options[option_index].name) == 0) {
                value = strtoll(optarg, NULL, 10);
                if (value < 0) {
                    fprintf(stderr, "Error, --stop %d must be nonnegative.\n", value);
                    usage();
                }
                *stop = value;
                setStop = true;
            }
            break;
        case 's':
            setSName = true;
            strncpy(seqName, optarg, kMaxSeqName);
            break;
        case 'v':
            g_verbose_flag++;
            break;
        case 'd':
            g_debug_flag = 1;
            break;
        case 'h':
        case '?':
            usage();
            break;
        default:
            abort();
        }
    }
    if (!(setSName && setStart && setStop)) {
        fprintf(stderr, "Error, specify --seq --start --stop\n");
        usage();
    }
    if (*start > *stop) {
        uint32_t t = *start;
        *start = *stop;
        *stop = t;
    }
    // Check there's nothing left over on the command line 
    if (optind < argc) {
        char *errorString = de_malloc(kMaxSeqName);
        strcpy(errorString, "Unexpected arguments:");
        while (optind < argc) {
            strcat(errorString, " ");
            strcat(errorString, argv[optind++]);
        }
        fprintf(stderr, "%s\n", errorString);
        usage();
    }
}

void usage(void) {
    fprintf(stderr, "Usage: mafBlockExtractor --seq [sequence name (and possibly chr)] "
            "--start [start of region, inclusive, 0 based] --stop [end of region, inclusive] "
            "[options] < myFile.maf\n\n"
            "mafBlockExtractor is a program that will look through a maf file for a\n"
            "particular sequence name and region. If a match is found then the block\n"
            "containing the querry will be printed to standard out.\n\n");
    fprintf(stderr, "Options: \n"
            "  -h, --help     show this help message and exit.\n"
            "  -s, --seq      sequence name.chr e.g. `hg18.chr2.'\n"
            "  --start        start of region, inclusive, 0 based.\n"
            "  --stop         end of region, inclusive, 0 based.\n"
            "  -v, --verbose  turns on verbose output.\n");
    exit(EXIT_FAILURE);
}

void processHeader(void) {
    FILE *ifp = stdin;
    int32_t n = kMaxSeqName;
    char *line = (char*) de_malloc(n);
    int status = de_getline(&line, &n, ifp);
    if (status == -1) {
        fprintf(stderr, "Error, empty file\n");
        exit(EXIT_FAILURE);
    }
    if (strncmp(line, "##maf", 5) != 0) {
        printf("line: %s\n", line);
        fprintf(stderr, "Error, bad maf format. File should start with ##maf\n");
        exit(EXIT_FAILURE);
    }
    printf("%s\n", line);
    while (*line != 0 && status != -1) {
        status = de_getline(&line, &n, ifp);
        printf("%s\n", line);
    }
}

bool checkRegion(uint32_t start, uint32_t stop, uint32_t str, 
                 uint32_t lng, uint32_t src, char strand) {
    // check to see if pos is in this block
    static unsigned long absStart, absEnd;
    if (strand == '-') {
        absStart = src - (str + lng);
        absEnd = lng - str - 1;
    } else {
        absStart = str;
        absEnd = str + lng - 1;
    }
    if (absEnd < start)
        return false;
    if (stop < absStart)
        return false;
    if ((absStart <= start) && (start <= absEnd))
        return true;
    if ((absStart <= stop) && (stop <= absEnd))
        return true;
    if ((start <= absStart) && (absEnd <= stop))
        return true;
    return false;
}

void reportBlock(mafBlock_t *b) {
    while (b != NULL) {
        printf("%s\n", b->line);
        b = b->next;
    }
    printf("\n");
}

void badFormat(void) {
    fprintf(stderr, "The maf sequence lines are incorrectly formatted, exiting\n");
    exit(EXIT_FAILURE);
}

void hardPrint(char *s, size_t n) {
    printf("hardPrint %s(%zu): ", s, n);
    for (unsigned i = 0; i < n; ++i){
        if (s[i] == 0)
            printf("_");
        else
            printf("%c", s[i]);
    }
    printf("\n");
}

int searchMatched(char *origLine, char *seq, uint32_t start, uint32_t stop) {
    // report 0 if search did not match, 1 if it did
    char *tkn = NULL;
    char strand;
    uint32_t str, lng, src;
    size_t n;
    char *cline = (char*) de_malloc(strlen(origLine) + 1);
    strcpy(cline, origLine);
    n = strlen(cline);
    debug("line: %s\n", cline);
    tkn = strtok(cline, " \t");
    if (tkn[0] != 's') {
        free(cline);
        return 0;
    }
    debug("tkn: %s\n", tkn);
    tkn = strtok(NULL, " \t"); // name field
    if (tkn == NULL)
        badFormat();
    debug("banana: %s\n", tkn);
    if (strcmp(seq, tkn) == 0) {
        tkn = strtok(NULL, " \t"); // start position
        if (tkn == NULL)
            badFormat();
        debug("apple: %s\n", tkn);
        str = strtoul(tkn, NULL, 10);
        tkn = strtok(NULL, " \t"); // length position
        if (tkn == NULL)
            badFormat();
        debug("carrot: %s\n", tkn);
        lng = strtoul(tkn, NULL, 10);
        if (lng == UINT32_MAX) {
            fprintf(stderr, "Error, length (%u) greater than "
                    "or equal to UINT32_MAX (%u).\n", lng, UINT32_MAX);
            exit(EXIT_FAILURE);
        }
        tkn = strtok(NULL, " \t"); // strand
        if (tkn == NULL)
            badFormat();
        debug("plum: %s\n", tkn);
        if ((strcmp(tkn, "+") == 0) || (strcmp(tkn, "-") == 0)) {
            strcpy(&strand, tkn);
        } else {
            fprintf(stderr, "Error, bad strand character: %s\n", tkn);
            exit(EXIT_FAILURE);
        }
        tkn = strtok(NULL, " \t"); // source length
        if (tkn == NULL)
            badFormat();
        debug("strawberry: %s\n", tkn);
        src = strtoul(tkn, NULL, 10);
        if (src == UINT32_MAX) {
            fprintf(stderr, "Error, source length (%u) greater "
                    "than or equal to UINT32_MAX (%u).\n", src, UINT32_MAX);
            exit(EXIT_FAILURE);
        }
        tkn = strtok(NULL, " \t"); // sequence field
        if (tkn == NULL)
            badFormat();
        debug("cherry: %s\n", tkn);
        if (checkRegion(start, stop, str, lng, src, strand)) {
            free(cline);
            return 1;
        }
    }
    free(cline);
    return 0;
}
void checkBlock(mafBlock_t *b, char *seq, uint32_t start, uint32_t stop) {
    // read through each line of a mafBlock and if the sequence matches the region
    // we're looking for, report the block.
    mafBlock_t *head = b;
    while (b != NULL) {
        if (searchMatched(b->line, seq, start, stop)) {
            reportBlock(head);
            break;
        } 
        b = b->next;
    }
}
void destroyBlock(mafBlock_t *b) {
    mafBlock_t *tmp = NULL;
    while (b != NULL) {
        tmp = b;
        b = b->next;
        free(tmp);
    }
}
void processBody(char *seq, uint32_t start, uint32_t stop) {
    FILE *ifp = stdin;
    int32_t n = kMaxSeqName;
    char *line = (char*) de_malloc(n);
    char *cline = NULL;
    mafBlock_t *head = NULL, *tail = NULL;
    while (de_getline(&line, &n, ifp) != -1) {
        if (*line == 0 && head != NULL) {
            // empty line or end of file
            checkBlock(head, seq, start, stop);
            destroyBlock(head);
            head = NULL;
            tail = NULL;
            continue;
        }
        if (head == NULL) {
            // new block
            head = (mafBlock_t *) de_malloc(sizeof(mafBlock_t *));
            cline = (char *) de_malloc(n);
            strcpy(cline, line);
            head->line = cline;
            head->next = NULL;
            tail = head;
        } else {
            // extend block
            tail->next = (mafBlock_t *) de_malloc(sizeof(mafBlock_t *));
            cline = (char *) de_malloc(n);
            strcpy(cline, line);
            tail = tail->next;
            tail->line = cline;
            tail->next = NULL;
        }
    }
}

int main(int argc, char **argv) {
    char seq[kMaxSeqName];
    uint32_t start, stop;
    parseOptions(argc, argv, seq, &start, &stop);

    processHeader();
    processBody(seq, start, stop);
    
    return EXIT_SUCCESS;
}
