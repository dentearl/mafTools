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

void usage(void);
void parseOptions(int argc, char **argv, char *filename, char *nameList, bool *isInclude);
void checkRegion(unsigned lineno, char *fullname, uint32_t pos, uint32_t start, 
                 uint32_t length, uint32_t sourceLength, char strand);

void usage(void) {
    fprintf(stderr, "Usage: mafBlockFilter --maf [path to maf] "
            "[options]\n\n"
            "mafBlockFilter is a program that will look through a\n"
            "maf file block by block and excise out sequence lines\n"
            "that match criteria established by the user on the\n"
            "command line. For example one can filter out all\n"
            "sequence lines that start with 'hg18' using --exclude\n"
            "or filter for sequence lines starting with only 'hg19',\n"
            "'mm9' and 'rn4' using --include.\n\n");
    fprintf(stderr, "Options: \n");
    usageMessage('h', "help", "show this help message and exit.");
    usageMessage('m', "maf", "path to maf file.");
    usageMessage('i', "include", "comma separated list of sequence names to include.");
    usageMessage('e', "exclude", "comma separated list of sequence names to exclude.");
    usageMessage('v', "verbose", "turns on verbose output.");
    exit(EXIT_FAILURE);
}
void parseOptions(int argc, char **argv, char *filename, char *nameList, bool *isInclude) {
    extern int g_debug_flag;
    extern int g_verbose_flag;
    int c;
    bool setMName = false, setName = false;
    while (1) {
        static struct option long_options[] = {
            {"debug", no_argument, &g_debug_flag, 1},
            {"verbose", no_argument, 0, 'v'},
            {"help", no_argument, 0, 'h'},
            {"maf",  required_argument, 0, 'm'},
            {"include",  required_argument, 0, 'i'},
            {"exclude",  required_argument, 0, 'e'},
            {0, 0, 0, 0}
        };
        int option_index = 0;
        c = getopt_long(argc, argv, "m:i:e:v:h",
                        long_options, &option_index);
        if (c == -1) {
            break;
        }
        switch (c) {
        case 0:
            break;
        case 'm':
            setMName = true;
            sscanf(optarg, "%s", filename);
            break;
        case 'i':
            setName = true;
            *isInclude = true;
            sscanf(optarg, "%s", nameList);
            break;
        case 'e':
            setName = true;
            *isInclude = false;
            sscanf(optarg, "%s", nameList);
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
    if (!(setMName && setName)) {
        fprintf(stderr, "specify --maf and [--include or --exclude]\n");
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
bool nameOnList(char *name, char **namelist, unsigned n) {
    for (unsigned i = 0; i < n; ++i) {
        if ((strncmp(name, namelist[i], strlen(namelist[i])) == 0)) {
            return true;
        }
    }
    return false;
}
void reportBlock(mafBlock_t *mb, char **names, unsigned n, bool isInclude) {
    // report the block being mindful of only including or excluding.
    mafLine_t *ml = maf_mafBlock_getHeadLine(mb);
    while (ml != NULL) {
        if (maf_mafLine_getType(ml) != 's') {
            // report all sequence lines
            printf("%s\n", maf_mafLine_getLine(ml));
            ml = maf_mafLine_getNext(ml);
            continue;
        }
        if (isInclude) {
            if (nameOnList(maf_mafLine_getSpecies(ml), names, n)) {
                printf("%s\n", maf_mafLine_getLine(ml));
                ml = maf_mafLine_getNext(ml);
                continue;
            }
        } else {
            if (!nameOnList(maf_mafLine_getSpecies(ml), names, n)) {
                printf("%s\n", maf_mafLine_getLine(ml));
                ml = maf_mafLine_getNext(ml);
                continue;
            }
        }
        ml = maf_mafLine_getNext(ml);
    }
    printf("\n");
}
void checkBlock(mafBlock_t *mb, char **names, unsigned n, bool isInclude) {
    // walk through the maf lines and see if this block should be reported
    mafLine_t *ml = maf_mafBlock_getHeadLine(mb);
    while (ml != NULL) {
        if (maf_mafLine_getType(ml) != 's') {
            ml = maf_mafLine_getNext(ml);
            ml = maf_mafLine_getNext(ml);
            continue;
        }
        if (isInclude) {
            if (nameOnList(maf_mafLine_getSpecies(ml), names, n)) {
                reportBlock(mb, names, n, isInclude);
                return;
            }
        } else {
            if (!nameOnList(maf_mafLine_getSpecies(ml), names, n)) {
                reportBlock(mb, names, n, isInclude);
                return;
            }
        }
        ml = maf_mafLine_getNext(ml);
    }
}
void filterInput(mafFileApi_t *mfa, char **names, unsigned n, bool isInclude) {
    mafBlock_t *thisBlock = NULL;
    bool headBlock = true;
    while ((thisBlock = maf_readBlock(mfa)) != NULL) {
        if (headBlock) {
            reportBlock(thisBlock, names, n, isInclude);
            headBlock = false;
            continue;
        }
        checkBlock(thisBlock, names, n, isInclude);
        maf_destroyMafBlockList(thisBlock);
    }
}
unsigned countNames(char *s) {
    unsigned i, n = 1;
    for (i = 0; i < strlen(s); ++i) {
        if (s[i] == ',') {
            ++n;
        }
    }
    return n;
}
char** extractNames(char *nameList, unsigned n) {
    // n is the number of names in the name list
    char **mat = (char**) de_malloc(sizeof(char*) * n);
    unsigned index = 0;
    char *tkn = NULL;
    char *copy = de_strdup(nameList);
    tkn = strtok(copy, ",");
    while (tkn != NULL) {
        mat[index] = (char*) de_malloc(sizeof(char) * (strlen(tkn) + 1));
        strcpy(mat[index++], tkn);
        tkn = strtok(NULL, ",");
    }
    free(copy);
    return mat;
}
void destroyNameList(char **names, unsigned n) {
    for (unsigned i = 0; i < n; ++i) {
        free(names[i]);
    }
    free(names);
}
void reportNames(char **names, unsigned n) {
    for (unsigned i = 0; i < n; ++i) {
        printf("name: %s\n", names[i]);
    }
}
int main(int argc, char **argv) {
    (void) (reportNames);
    char filename[kMaxStringLength];
    char nameList[kMaxStringLength];
    bool isInclude = true; // if 0 then we are in exclude mode. 1 is include mode.
    parseOptions(argc, argv,  filename, nameList, &isInclude);
    unsigned n = countNames(nameList);
    char **names = extractNames(nameList, n);
    mafFileApi_t *mfa = maf_newMfa(filename, "r");
    
    filterInput(mfa, names, n, isInclude);
    
    maf_destroyMfa(mfa);
    destroyNameList(names, n);
    
    return EXIT_SUCCESS;
}
