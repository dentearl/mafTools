/* 
 * Copyright (C) 2009-2012 by 
 * Dent Earl (dearl@soe.ucsc.edu, dentearl@gmail.com)
 * Benedict Paten (benedict@soe.ucsc.edu, benedictpaten@gmail.com)
 * Mark Diekhans (markd@soe.ucsc.edu)
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
#include <limits.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <getopt.h>
#include <stdarg.h>
#include <unistd.h> // getpid

#include "sonLib.h"
#include "comparatorAPI.h"
#include "common.h"

const char *g_version = "version 0.6 July 2012";
bool g_isVerboseFailures = false;

/*
 * The script takes two MAF files and for each ordered pair of sequences 
 * in the MAFS calculates a predefined number of sample homology tests 
 * (see below), then reports the statistics in an XML formatted file.
 * It is suitable for running over very large numbers of alignments, 
 * because it does not attempt to hold everything in memory, and instead 
 * takes a sampling approach.
 *
 * For two sets of pairwise alignments, A and B, a homology test is 
 * defined as follows. Pick a pair of aligned positions in A, called a 
 * homology pair - the AB homology test returns true if the pair is in B, 
 * otherwise it returns false. The set of possible homology tests for the 
 * ordered pair (A, B) is not necessarily equivalent to the set of 
 * possible (B, A) homology tests. We call the proportion of true tests 
 * as a percentage of the total of a set of homology tests C from 
 * (A, B)  A~B.
 *
 * If A is the set of true pairwise alignments and B the predicted set of 
 * alignments then A~B (over large enough  C), is a proxy to sensitivity 
 * of B in predicted the set of correctly aligned pairs in A. Conversely 
 * B~A (over large enough C) is a proxy to the specificity of the 
 * aligned pairs in B with respect to the set of correctly aligned pairs 
 * in A.
 */

void parseBedFile(const char *filepath, stHash *intervalsHash) {
    /* 
     * takes a filepath and the intervalsHash, opens and reads the file,
     * adding intervals taken from each bed line to the intervalsHash.
     */
    FILE *fileHandle = fopen(filepath, "r");
    st_logDebug("Parsing the bed file: %s\n", filepath);
    if (fileHandle == NULL) {
        if (errno == ENOENT)
            fprintf(stderr, "ERROR, file %s does not exist.\n", filepath);
        else
            fprintf(stderr, "ERROR, unable to open %s\n", filepath);
        exit(EXIT_FAILURE);
    }
    int nBytes = 100;
    char *cA2 = st_malloc(nBytes + 1);
    int32_t bytesRead = benLine(&cA2, &nBytes, fileHandle);
    while (bytesRead != -1) {
        if (bytesRead > 0) {
            int32_t start, stop;
            char *cA3 = stString_copy(cA2);
            int32_t i = sscanf(cA2, "%s %" PRIi32 " %" PRIi32 "", cA3, &start, &stop);
            assert(i == 3);
            stSortedSet *intervals = stHash_search(intervalsHash, cA3);
            if (intervals == NULL) {
                intervals = stSortedSet_construct3(
                        (int(*)(const void *, const void *)) stIntTuple_cmpFn,
                        (void(*)(void *)) stIntTuple_destruct);
                stHash_insert(intervalsHash, stString_copy(cA3), intervals);
            }
            stIntTuple *j = stIntTuple_construct(2, start, stop);
            stIntTuple *k = stSortedSet_searchLessThanOrEqual(intervals,
                    j);
            if (k != NULL) {
                if (stIntTuple_getPosition(k, 1) > start) {
                    st_errAbort(
                            "Found an overlapping interval in the bed file: %s %" PRIi32 " "
                            "%" PRIi32 " overlaps %s %" PRIi32 " %" PRIi32 "",
                            cA3, start, stop, cA2,
                            stIntTuple_getPosition(k, 0),
                            stIntTuple_getPosition(k, 1));
                }
            }
            k = stSortedSet_searchGreaterThanOrEqual(intervals, j);
            if (k != NULL) {
                if (stIntTuple_getPosition(k, 1) < stop) {
                    st_errAbort(
                            "Found an overlapping interval in the bed file: %s %" PRIi32 " "
                            "%" PRIi32 " overlaps %s %" PRIi32 " %" PRIi32 "",
                            cA3, start, stop, cA2,
                            stIntTuple_getPosition(k, 0),
                            stIntTuple_getPosition(k, 1));
                }
            }
            assert(stSortedSet_search(intervals, j) == NULL);
            st_logDebug("Adding in an interval: %s %" PRIi32 " %" PRIi32 "\n", cA3, start, stop);
            stSortedSet_insert(intervals, j);
            free(cA3);
        }
        bytesRead = benLine(&cA2, &nBytes, fileHandle);
    }
    free(cA2);
    fclose(fileHandle);
    st_logDebug("Finished parsing the bed file: %s\n", filepath);
}
void parseBedFiles(const char *commaSepFiles, stHash *bedFileHash) {
    /*
     * takes input from the command line, swaps spaces, ' ', for commas
     * and passes each bed file to parseBedFile().
     */
    st_logDebug("Starting to parse bed files\n");
    char *spaceSepFiles = stringCommasToSpaces(commaSepFiles);
    char *currentLocation = spaceSepFiles;
    char *currentWord;
    while ((currentWord = stString_getNextWord(&currentLocation)) != NULL) {
        parseBedFile(currentWord, bedFileHash);
        free(currentWord);
    }
    free(spaceSepFiles);
    st_logDebug("Done parsing bed files\n");
}
void usage(void) {
    fprintf(stderr, "mafComparator, %s\n\n", g_version);
    fprintf(stderr, "Usage: $ mafComparator --maf1=FILE1 --maf2=FILE2 --out=OUT.xml [options]\n\n");
    fprintf(stderr, "This program takes two MAF files and compares them to one another.\n"
            "Specifically, for each ordered pair of sequences in the first MAF it \n"
            "samples a predefined number of sample homology tests (see below), then \n"
            "reads the second MAF checking to see which, if any, of the sampled pairs, \n"
            "is present. The comparison is then reversed and repeated. Statistics are\n"
            "then reported in an XML formatted file. \n\n");
    fprintf(stderr, "Options:\n");
    usageMessage('h', "help", "Print this help screen.");
    usageMessage('\0', "maf1", "The location of the first MAF file. "
                 "If comparing true to predicted "
                 "alignments, this is the truth.");
    usageMessage('\0', "maf2", "The location of the second MAF file. "
                 "If comparing true to predicted "
                 "alignments, this is the prediction.");
    usageMessage('\0', "out", "The output XML formatted results file.");
    usageMessage('\0', "samples", "The ideal number of sample homology tests to perform for the "
                 "two comparisons (i.e. file1 -> file and file2 -> file1). This "
                 "number is an ideal because pairs are sampled and thus the "
                 "actual number may be slightly higher or slightly lower than "
                 "this value. If this value is equal to or greater than the "
                 "total number of pairs in a file, then all pairs will be "
                 "tested. [default 1000000]");
    usageMessage('\0', "bedFiles", "The location of bed file(s) used to filter the "
                 "pairwise comparisons, separated by commas.");
    usageMessage('\0', "near", "The number of bases in either sequence to allow a match "
                 "to slip by. I.e. --near=n (where _n_ is a non-negative integer) "
                 "will consider a homology test for a given pair (S1:_x_, S2:_y_) "
                 "where S1 and S2 are sequences and _x_ and _y_ are positions in "
                 "the respective sequences, to be a true homology test so long as "
                 "there is a pair within the other alignment (S1:_w_, S2:_z_) where "
                 "EITHER (_w_ is equal to _x_ and _y_ - _n_ <= _z_ <= _y_ + _n_) "
                 "OR (_x_ - _n_ <= _w_ <= _x_ + _n_ and _y_ is equal to _z_).");
    usageMessage('\0', "logLevel", "Set the log level. [off, critical, info, debug] "
                 "in ascending order.");
    usageMessage('\0', "printFailed", "Print tab-delimited details about failed "
                 "tests to stderr.");
    usageMessage('\0', "seed", "an integer used to seed the random number generator "
                 "used to perform sampling. If omitted a seed is pseudorandomly "
                 "generated. The seed value is always stored in the output xml.");
    usageMessage('v', "version", "Print current version number.");
}
void version(void) {
    fprintf(stderr, "mafComparator, %s\n", g_version);
}
int parseArgs(int argc, char **argv, char **mafFile1, char **mafFile2, char **outputFile, 
              char **bedFiles, uint32_t *sampleNumber, uint32_t *randomSeed, uint32_t *near, 
              char **logLevelString) {
    static const char *optString = "a:b:c:d:e:p:v:h:f:g:s:";
    static const struct option longOpts[] = {
        {"logLevel", required_argument, 0, 'a'},
        {"mafFile1", required_argument, 0, 'b'},
        {"maf1", required_argument, 0, 0},
        {"mafFile2", required_argument, 0, 'c'},
        {"maf2", required_argument, 0, 0},
        {"outputFile", required_argument, 0, 'd'},
        {"out", required_argument, 0, 0},
        {"sampleNumber", required_argument, 0, 'e'},
        {"samples", required_argument, 0, 0},
        {"printFailed", no_argument, 0, 'p'},
        {"version", no_argument, 0, 'v'},
        {"help", no_argument, 0, 'h'},
        {"bedFiles", required_argument, 0, 'f'},
        {"near", required_argument, 0, 'g'},
        {"seed", required_argument, 0, 's'},
        {0, 0, 0, 0 }};
    int longIndex = 0;
    size_t i;
    int key = getopt_long(argc, argv, optString, longOpts, &longIndex);
    while (key != -1) {
        switch (key) {
        case 0:
            if (strcmp("maf1", longOpts[longIndex].name) == 0) {
                *mafFile1 = stString_copy(optarg);
                break;
            }
            if (strcmp("maf2", longOpts[longIndex].name) == 0) {
                *mafFile2 = stString_copy(optarg);
                break;
            }
            if (strcmp("out", longOpts[longIndex].name) == 0) {
                *outputFile = stString_copy(optarg);
                break;
            }
            if (strcmp("samples", longOpts[longIndex].name) == 0) {
                i = sscanf(optarg, "%" PRIu32, sampleNumber);
                assert(i == 1);
                break;
            }
        case 'a':
            *logLevelString = stString_copy(optarg);
            break;
        case 'b':
            *mafFile1 = stString_copy(optarg);
            break;
        case 'c':
            *mafFile2 = stString_copy(optarg);
            break;
        case 'd':
            *outputFile = stString_copy(optarg);
            break;
        case 'v':
            version();
            exit(EXIT_SUCCESS);
        case 'p':
            g_isVerboseFailures = true;
            break;
        case 'e':
            i = sscanf(optarg, "%" PRIu32, sampleNumber);
            assert(i == 1);
            break;
        case 's':
            i = sscanf(optarg, "%" PRIu32, randomSeed);
            assert(i == 1);
            break;
        case 'h':
            usage();
            exit(EXIT_SUCCESS);
        case 'f':
            *bedFiles = stString_copy(optarg);
            break;
        case 'g':
            i = sscanf(optarg, "%" PRIu32, near);
            assert(i == 1);
            break;
        default:
            usage();
            exit(EXIT_SUCCESS);
        }
        key = getopt_long(argc, argv, optString, longOpts, &longIndex);
    }
    if (*mafFile1 == NULL) {
        usage();
        fprintf(stderr, "\nError, specify --mafFile1\n");
        exit(2);
    }
    if (*mafFile2 == NULL) {
        usage();
        fprintf(stderr, "\nError, specify --mafFile2\n");
        exit(2);
    }
    if (*outputFile == NULL) {
        usage();
        fprintf(stderr, "\nError, specify --outputFile\n");
        exit(2);
    }
    FILE *fileHandle = de_fopen(*mafFile1, "r");
    fclose(fileHandle);
    fileHandle = de_fopen(*mafFile2, "r");
    fclose(fileHandle);
    return optind;
}
int main(int argc, char **argv) {
    char *logLevelString = NULL;
    char *mafFile1 = NULL;
    char *mafFile2 = NULL;
    char *outputFile = NULL;
    FILE *fileHandle = NULL;
    stHash *intervalsHash = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, free,
                                              (void(*)(void *)) stSortedSet_destruct);
    char *bedFiles = NULL;
    uint32_t sampleNumber = 1000000; // by default do a million samples per pair.
    uint32_t randomSeed = (time(NULL) << 16) | (getpid() & 65535); // Likely to be unique
    uint32_t near = 0;
    // (0) Parse the inputs handed by genomeCactus.py / setup stuff.
    parseArgs(argc, argv, &mafFile1, &mafFile2, &outputFile, &bedFiles, &sampleNumber, 
              &randomSeed, &near, &logLevelString);
    // Set up logging
    st_setLogLevelFromString(logLevelString);
    st_logDebug("Seeding the random number generator with the value %lo\n", randomSeed);
    st_randomSeed(randomSeed);
    // (0) Check the inputs.
    //Parse the bed file hashes
    if(bedFiles != NULL) {
        parseBedFiles(bedFiles, intervalsHash);
    }
    else {
        st_logDebug("No bed files specified\n");
    }
    //Log (some of) the inputs
    st_logInfo("MAF file 1 name : %s\n", mafFile1);
    st_logInfo("MAF file 2 name : %s\n", mafFile2);
    st_logInfo("Output stats file : %s\n", outputFile);
    st_logInfo("Bed files parsed : %" PRIi32 "\n", stHash_size(intervalsHash));
    st_logInfo("Number of samples %" PRIu32 "\n", sampleNumber);
    // note that random seed has already been logged.
    // Create sequence name hashtable from the first MAF file.
    stHash *seqNames1 = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, free, free);
    populateNames(mafFile1, seqNames1);
    stHash *seqNames2 = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, free, free);
    populateNames(mafFile2, seqNames2);
    stHash *seqNames = stHash_getIntersection(seqNames1, seqNames2);
    // Do comparisons.
    if (g_isVerboseFailures) {
        fprintf(stderr, "# Sampling from %s, comparing to %s\n", mafFile1, mafFile2);
        fprintf(stderr, "# seq1\tabsPos1\torigPos1\tseq2\tabsPos2\torigPos2\n");
    }
    uint32_t numPairs1 = 0;
    stSortedSet *results_12 = compareMAFs_AB(mafFile1, mafFile2, sampleNumber, &numPairs1, seqNames, 
                                             intervalsHash, near);
    if (g_isVerboseFailures) {
        fprintf(stderr, "# Sampling from %s, comparing to %s\n", mafFile2, mafFile1);
        fprintf(stderr, "# seq1\tabsPos1\torigPos1\tseq2\tabsPos2\torigPos2\n");
    }
    uint32_t numPairs2 = 0;
    stSortedSet *results_21 = compareMAFs_AB(mafFile2, mafFile1, sampleNumber, &numPairs2, seqNames, 
                                             intervalsHash, near);
    fileHandle = de_fopen(outputFile, "w");
    // Report results.
    writeXMLHeader(fileHandle);
    char bedString[kMaxStringLength];
    if (bedFiles != NULL) {
        sprintf(bedString, " bedFiles=\"%s\"", bedFiles);
    } else {
        bedString[0] = '\0';
    }
    fprintf(fileHandle, "<alignmentComparisons sampleNumber=\"%" PRIu32 "\" "
            "near=\"%" PRIu32 "\" seed=\"%" PRIu32 "\" maf1=\"%s\" maf2=\"%s\" "
            "numberOfPairsInMaf1=\"%" PRIu32 "\" "
            "numberOfPairsInMaf2=\"%" PRIu32 "\"%s>\n",
            sampleNumber, near, randomSeed, mafFile1, mafFile2, numPairs1, numPairs2, bedString);
        
    reportResults(results_12, mafFile1, mafFile2, fileHandle, near, seqNames, bedFiles);
    reportResults(results_21, mafFile2, mafFile1, fileHandle, near, seqNames, bedFiles);
    fprintf(fileHandle, "</alignmentComparisons>\n");
    fclose(fileHandle);
    // Clean up.
    stSortedSet_destruct(results_12);
    stSortedSet_destruct(results_21);
    free(mafFile1);
    free(mafFile2);
    free(bedFiles);
    free(outputFile);
    free(logLevelString);
    stHash_destruct(seqNames);
    stHash_destruct(seqNames1);
    stHash_destruct(seqNames2);
    stHash_destruct(intervalsHash);
    return(EXIT_SUCCESS);
}
