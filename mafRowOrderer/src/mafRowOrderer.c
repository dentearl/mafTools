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
#include <errno.h> // file existence via ENOENT
#include <getopt.h>
#include <inttypes.h>
#include <limits.h>
#include <stdarg.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>
#include "common.h"
#include "sharedMaf.h"
#include "buildVersion.h"

const char *g_version = "version 0.1 October 2012";

void version(void);
void usage(void);
void parseOptions(int argc, char **argv, char *filename, char *orderlist);
void checkRegion(unsigned lineno, char *fullname, uint32_t pos, uint32_t start, 
                 uint32_t length, uint32_t sourceLength, char strand);
void printHeader(void);

void version(void) {
    fprintf(stderr, "mafRowOrderer, %s\nbuild: %s, %s, %s\n\n", g_version, g_build_date, 
            g_build_git_branch, g_build_git_sha);
}
void usage(void) {
    version();
    fprintf(stderr, "Usage: mafBlockFilter --maf [path to maf] --order [comma separated list of species] "
            "[options]\n\n"
            "mafRowOrderer is a program that will look through a maf file block by block and order\n"
            "the maf lines within a block according to the order provided. Species not in the\n"
            "established ordered are excised. Comments are excised. Non sequnece lines ('^s') are excised.\n"
            );
    fprintf(stderr, "Options: \n");
    usageMessage('h', "help", "show this help message and exit.");
    usageMessage('m', "maf", "path to maf file.");
    usageMessage('\0', "order", "comma separated list of sequence names.");
    usageMessage('v', "verbose", "turns on verbose output.");
    exit(EXIT_FAILURE);
}
void parseOptions(int argc, char **argv, char *filename, char *orderlist) {
    extern int g_debug_flag;
    extern int g_verbose_flag;
    int c;
    bool setMafName = false, setOrder = false;
    while (1) {
        static struct option longOptions[] = {
            {"debug", no_argument, &g_debug_flag, 1},
            {"verbose", no_argument, 0, 'v'},
            {"help", no_argument, 0, 'h'},
            {"version", no_argument, 0, 0},
            {"maf",  required_argument, 0, 'm'},
            {"order",  required_argument, 0, 0},
            {0, 0, 0, 0}
        };
        int longIndex = 0;
        c = getopt_long(argc, argv, "m:i:e:g:l:v:h",
                        longOptions, &longIndex);
        if (c == -1) {
            break;
        }
        switch (c) {
        case 0:
            if (strcmp("version", longOptions[longIndex].name) == 0) {
                version();
                exit(EXIT_SUCCESS);
            }
            if (strcmp("order", longOptions[longIndex].name) == 0) {
                setOrder = true;
                sscanf(optarg, "%s", orderlist);
                break;
            }
            break;
        case 'm':
            setMafName = true;
            sscanf(optarg, "%s", filename);
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
    if (!setMafName) {
        fprintf(stderr, "specify --maf\n");
        usage();
    }
    if (!setOrder) {
        fprintf(stderr, "specify --order\n");
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
void printHeader(void) {
    printf("##maf version=1\n\n");
}
void checkBlock(mafBlock_t *mb, char **order, unsigned n) {
    // the plan:
    // create an array of mafLine_t linked lists, of length n
    // walk the block, *copying* mafLines into the linked list at the coresponding array element
    // then attach all existing array elements head to tail to form a single, ordered, 
    // linked list. Report this. Free everything.
    mafLine_t **lineArrayHeads = (mafLine_t**) de_malloc(sizeof(mafLine_t*) * n);
    mafLine_t **lineArrayTails = (mafLine_t**) de_malloc(sizeof(mafLine_t*) * n);
    unsigned i;
    for (i = 0; i < n; ++i) {
        lineArrayHeads[i] = NULL;
        lineArrayTails[i] = NULL;
    }
    mafLine_t *ml = maf_mafBlock_getHeadLine(mb);
    // build array
    while (ml != NULL) {
        if (maf_mafLine_getType(ml) != 's') {
            ml = maf_mafLine_getNext(ml);
            continue;
        }
        for (i = 0; i < n; ++i) {
            if (strncmp(order[i], maf_mafLine_getSpecies(ml), strlen(order[i])) == 0) {
                if (lineArrayHeads[i] == NULL) {
                    lineArrayHeads[i] = maf_copyMafLine(ml);
                    lineArrayTails[i] = lineArrayHeads[i];
                } else {
                    maf_mafLine_setNext(lineArrayTails[i], maf_copyMafLine(ml));
                    lineArrayTails[i] = maf_mafLine_getNext(lineArrayTails[i]);
                }
                break;
            }
        }
        ml = maf_mafLine_getNext(ml);
    }
    // connect array heads / tails
    mafLine_t *head = NULL;
    int prev = -1;
    bool reportBlock = false;
    for (i = 0; i < n; ++i) {
        if (lineArrayHeads[i] == NULL) {
            continue;
        }
        reportBlock = true;
        if (head == NULL) {
            head = lineArrayHeads[i];
        }
        if (prev > -1) {
            maf_mafLine_setNext(lineArrayTails[prev], lineArrayHeads[i]);
        }
        prev = i;
    }
    // put into a dummy block
    mafBlock_t *orderedBlock = maf_newMafBlock();
    maf_mafBlock_setHeadLine(orderedBlock, maf_newMafLineFromString("a ordered=true", 0));
    maf_mafLine_setNext(maf_mafBlock_getHeadLine(orderedBlock), head);
    // report block
    if (reportBlock) {
        maf_mafBlock_print(orderedBlock);
    }
    maf_destroyMafBlockList(orderedBlock);
    free(lineArrayHeads);
    free(lineArrayTails);
}
void orderInput(mafFileApi_t *mfa, char **order, unsigned n) {
    mafBlock_t *thisBlock = NULL;
    bool headBlock = true;
    printHeader();
    while ((thisBlock = maf_readBlock(mfa)) != NULL) {
        if (headBlock) {
            headBlock = false;
            maf_destroyMafBlockList(thisBlock);
            continue;
        }
        checkBlock(thisBlock, order, n);
        maf_destroyMafBlockList(thisBlock);
    }
}
void destroyNameList(char **names, unsigned n) {
    for (unsigned i = 0; i < n; ++i) {
        free(names[i]);
    }
    free(names);
}
int main(int argc, char **argv) {
    char filename[kMaxStringLength];
    char orderlist[kMaxStringLength];
    orderlist[0] = '\0';
    parseOptions(argc, argv,  filename, orderlist);
    unsigned n = 1 + countChar(orderlist, ',');
    char **order = extractSubStrings(orderlist, n, ',');
    mafFileApi_t *mfa = maf_newMfa(filename, "r");
    orderInput(mfa, order, n);
    maf_destroyMfa(mfa);
    destroyNameList(order, n);
    return EXIT_SUCCESS;
}
