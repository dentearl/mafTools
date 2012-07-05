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
#include "sharedMaf.h"

void usage(void);
void processBody(mafFileApi_t *mfa, char *seq, uint32_t start, uint32_t stop);
bool searchMatched(mafLine_t *ml, char *seq, uint32_t start, uint32_t stop);

void parseOptions(int argc, char **argv, char *filename, char *seqName, uint32_t *start, uint32_t *stop) {
    extern int g_debug_flag;
    extern int g_verbose_flag;
    int c;
    bool setSName = false, setStart = false, setStop = false, setMName = false;
    int32_t value = 0;
    while (1) {
        static struct option long_options[] = {
            {"debug", no_argument, 0, 'd'},
            {"verbose", no_argument, 0, 'v'},
            {"help", no_argument, 0, 'h'},
            {"maf",  required_argument, 0, 'm'},
            {"seq",  required_argument, 0, 's'},
            {"start",  required_argument, 0, 0},
            {"stop",  required_argument, 0, 0},
            {0, 0, 0, 0}
        };
        int option_index = 0;
        c = getopt_long(argc, argv, "m:s:h:v:d",
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
        case 'm':
            setMName = true;
            strncpy(filename, optarg, kMaxSeqName);
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
    if (!(setMName && setSName && setStart && setStop)) {
        fprintf(stderr, "Error, specify --maf --seq --start --stop\n");
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
    fprintf(stderr, "Usage: mafBlockExtractor --maf [maf file] --seq [sequence name (and possibly chr)] "
            "--start [start of region, inclusive, 0 based] --stop [end of region, inclusive] "
            "[options]\n\n"
            "mafBlockExtractor is a program that will look through a maf file for a\n"
            "particular sequence name and region. If a match is found then the block\n"
            "containing the querry will be printed to standard out.\n\n");
    fprintf(stderr, "Options: \n");
    usageMessage('h', "help", "show this help message and exit.");
    usageMessage('m', "maf", "path to maf file.");
    usageMessage('s', "seq", "sequence name, e.g. `hg18.chr2'.");
    usageMessage('\0', "start", "start of region, inclusive, 0 based.");
    usageMessage('\0', "stop", "end of region, inclusive, 0 based.");
    usageMessage('v', "verbose", "turns on verbose output.");
    exit(EXIT_FAILURE);
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
    mafLine_t *ml = maf_mafBlock_getHeadLine(b);
    while (ml != NULL) {
        printf("%s\n", maf_mafLine_getLine(ml));
        ml = maf_mafLine_getNext(ml);
    }
    printf("\n");
}
bool searchMatched(mafLine_t *ml, char *seq, uint32_t start, uint32_t stop) {
    // report false if search did not match, true if it did
    if (maf_mafLine_getType(ml) != 's')
        return false;
    if (!(strcmp(maf_mafLine_getSpecies(ml), seq) == 0))
        return false;
    if (checkRegion(start, stop, maf_mafLine_getStart(ml), maf_mafLine_getLength(ml), 
                    maf_mafLine_getSourceLength(ml), maf_mafLine_getStrand(ml)))
        return true;
    return false;
}
void printHeader(void) {
    printf("##maf version=1\n\n");
}
void checkBlock(mafBlock_t *b, char *seq, uint32_t start, uint32_t stop, bool *printedHeader) {
    // read through each line of a mafBlock and if the sequence matches the region
    // we're looking for, report the block.
    mafLine_t *ml = maf_mafBlock_getHeadLine(b);
    while (ml != NULL) {
        if (searchMatched(ml, seq, start, stop)) {
            if (!*printedHeader) {
                printHeader();
                *printedHeader = true;
            }
            reportBlock(b);
            break;
        } 
        ml = maf_mafLine_getNext(ml);
    }
}
void processBody(mafFileApi_t *mfa, char *seq, uint32_t start, uint32_t stop) {
    mafBlock_t *thisBlock = NULL;
    bool printedHeader = false;
    while ((thisBlock = maf_readBlock(mfa)) != NULL) {
        checkBlock(thisBlock, seq, start, stop, &printedHeader);
        maf_destroyMafBlockList(thisBlock);
    }
}

int main(int argc, char **argv) {
    extern const int kMaxStringLength;
    char seq[kMaxSeqName];
    char filename[kMaxStringLength];
    uint32_t start, stop;
    parseOptions(argc, argv, filename, seq, &start, &stop);
    mafFileApi_t *mfa = maf_newMfa(filename, "r");

    processBody(mfa, seq, start, stop);
    maf_destroyMfa(mfa);
    
    return EXIT_SUCCESS;
}
