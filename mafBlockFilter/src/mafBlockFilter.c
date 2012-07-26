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
void parseOptions(int argc, char **argv, char *filename, char *nameList, bool *isInclude, int64_t *blockDegLT, int64_t *blockDegGT);
void checkRegion(unsigned lineno, char *fullname, uint32_t pos, uint32_t start, 
                 uint32_t length, uint32_t sourceLength, char strand);

void usage(void) {
    fprintf(stderr, "Usage: mafBlockFilter --maf [path to maf] "
            "[options]\n\n"
            "mafBlockFilter is a program that will look through a\n"
            "maf file block by block and excise out sequence lines\n"
            "that match criteria established by the user on the\n"
            "command line. For example one can filter out all\n"
            "sequence lines that start with 'hg18' using --excludeSeq\n"
            "or filter for sequence lines starting with only 'hg19',\n"
            "'mm9' and 'rn4' using --includeSeq.\n\n");
    fprintf(stderr, "Options: \n");
    usageMessage('h', "help", "show this help message and exit.");
    usageMessage('m', "maf", "path to maf file.");
    usageMessage('i', "includeSeq", "comma separated list of sequence names to include.");
    usageMessage('e', "excludeSeq", "comma separated list of sequence names to exclude.");
    usageMessage('g', "noDegreeGT", "filter out all blocks with degree greater than this value.");
    usageMessage('l', "noDegreeLT", "filter out all blocks with degree less than this value.");
    usageMessage('v', "verbose", "turns on verbose output.");
    exit(EXIT_FAILURE);
}
void parseOptions(int argc, char **argv, char *filename, char *nameList, bool *isInclude, int64_t *blockDegGt, int64_t *blockDegLt) {
    extern int g_debug_flag;
    extern int g_verbose_flag;
    int c;
    bool setMafName = false, setNames = false, setBlockLimits = false;
    while (1) {
        static struct option long_options[] = {
            {"debug", no_argument, &g_debug_flag, 1},
            {"verbose", no_argument, 0, 'v'},
            {"help", no_argument, 0, 'h'},
            {"maf",  required_argument, 0, 'm'},
            {"includeSeq",  required_argument, 0, 'i'},
            {"excludeSeq",  required_argument, 0, 'e'},
            {"noDegreeGT", required_argument, 0, 'g'},
            {"noDegreeLT", required_argument, 0, 'l'},
            {0, 0, 0, 0}
        };
        int option_index = 0;
        c = getopt_long(argc, argv, "m:i:e:g:l:v:h",
                        long_options, &option_index);
        if (c == -1) {
            break;
        }
        switch (c) {
        case 0:
            break;
        case 'm':
            setMafName = true;
            sscanf(optarg, "%s", filename);
            break;
        case 'i':
            setNames = true;
            *isInclude = true;
            sscanf(optarg, "%s", nameList);
            break;
        case 'e':
            setNames = true;
            *isInclude = false;
            sscanf(optarg, "%s", nameList);
            break;
        case 'g':
            setNames = true;
            sscanf(optarg, "%" PRIi64, blockDegGt);
            break;
        case 'l':
            setNames = true;
            sscanf(optarg, "%" PRIi64, blockDegLt);
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
    if ((setNames && setBlockLimits) || !(setNames || setBlockLimits)) {
        fprintf(stderr, "specify *one* from [--includeSeq or --excludeSeq] or [--noDegreeGT or --noDegreeLT]\n");
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
        if (n > 0) {
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
        } else {
            // report entire block, this came from one of the blockDegree options
            printf("%s\n", maf_mafLine_getLine(ml));
        }
        ml = maf_mafLine_getNext(ml);
    }
    printf("\n");
}
void checkBlock(mafBlock_t *mb, char **names, unsigned n, bool isInclude, int64_t excludeBlockDegreeGT,
                 int64_t excludeBlockDegreeLT) {
    // walk through the maf lines and see if this block should be reported
    mafLine_t *ml = maf_mafBlock_getHeadLine(mb);
    while (ml != NULL) {
        if (maf_mafLine_getType(ml) != 's') {
            ml = maf_mafLine_getNext(ml);
            continue;
        }
        if (n > 0) {
            // filtering on names
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
        } else {
            // filtering on block degrees
            uint32_t m = maf_mafBlock_getNumberOfSequences(mb);
            if (excludeBlockDegreeGT != -1 && excludeBlockDegreeLT != -1) {
                if (m >= excludeBlockDegreeLT && m <= excludeBlockDegreeGT) {
                    reportBlock(mb, names, n, isInclude);
                    return;
                }
            } else if (excludeBlockDegreeGT != -1) {
                if (m <= excludeBlockDegreeGT) {
                    reportBlock(mb, names, n, isInclude);
                    return;
                }
            } else {
                if (m >= excludeBlockDegreeLT) {
                    reportBlock(mb, names, n, isInclude);
                    return;
                }
            }
        }
        ml = maf_mafLine_getNext(ml);
    }
}
void filterInput(mafFileApi_t *mfa, char **names, unsigned n,
                 bool isInclude, int64_t excludeBlockDegreeGT,
                 int64_t excludeBlockDegreeLT) {
    mafBlock_t *thisBlock = NULL;
    bool headBlock = true;
    while ((thisBlock = maf_readBlock(mfa)) != NULL) {
        if (headBlock) {
            reportBlock(thisBlock, names, n, isInclude);
            headBlock = false;
            maf_destroyMafBlockList(thisBlock);
            continue;
        }
        checkBlock(thisBlock, names, n, isInclude, excludeBlockDegreeGT, excludeBlockDegreeLT);
        maf_destroyMafBlockList(thisBlock);
    }
}
unsigned countNames(char *s) {
    unsigned i, n;
    if (strlen(s) == 0) {
        return 0;
    }
    n = 1;
    for (i = 0; i < strlen(s); ++i) {
        if (s[i] == ',') {
            ++n;
        }
    }
    return n;
}
char** extractNames(char *nameList, unsigned n) {
    // n is the number of names in the name list
    if (n == 0) {
        return NULL;
    }
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
    nameList[0] = '\0';
    int64_t excludeBlockDegreeGT = -1;
    int64_t excludeBlockDegreeLT = -1;
    bool isInclude = true; // if 0 then we are in exclude mode. 1 is include mode.
    parseOptions(argc, argv,  filename, nameList, &isInclude, &excludeBlockDegreeGT, &excludeBlockDegreeLT);
    unsigned n = countNames(nameList);
    char **names = extractNames(nameList, n);
    mafFileApi_t *mfa = maf_newMfa(filename, "r");
    
    filterInput(mfa, names, n, isInclude, excludeBlockDegreeGT, excludeBlockDegreeLT);
    
    maf_destroyMfa(mfa);
    destroyNameList(names, n);
    
    return EXIT_SUCCESS;
}
