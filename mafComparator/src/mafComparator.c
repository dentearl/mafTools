/* 
 * Copyright (C) 2009-2011 by 
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
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <getopt.h>
#include <stdarg.h>
#include <unistd.h> // getpid

//#include "cactus.h"
#include "avl.h"
#include "commonC.h"
#include "hashTableC.h"
#include "hashTableC_itr.h"
#include "bioioC.h"
#include "sonLibRandom.h" // seeds and sampling

#include "ComparatorAPI.h"

float VERSION = 0.3;
int32_t VERBOSEFAILURES = 0;

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
        exit(1);
    }
    int nBytes = 100;
    char *cA2 = st_malloc(nBytes + 1);
    int32_t bytesRead = benLine(&cA2, &nBytes, fileHandle);

    while (bytesRead != -1) {
        if (bytesRead > 0) {
            int32_t start, stop;
            char *cA3 = stString_copy(cA2);
            int32_t i = sscanf(cA2, "%s %i %i", cA3, &start, &stop);
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
                            "Found an overlapping interval in the bed file: %s %i %i overlaps %s %i %i",
                            cA3, start, stop, cA2,
                            stIntTuple_getPosition(k, 0),
                            stIntTuple_getPosition(k, 1));
                }
            }

            k = stSortedSet_searchGreaterThanOrEqual(intervals, j);
            if (k != NULL) {
                if (stIntTuple_getPosition(k, 1) < stop) {
                    st_errAbort(
                            "Found an overlapping interval in the bed file: %s %i %i overlaps %s %i %i",
                            cA3, start, stop, cA2,
                            stIntTuple_getPosition(k, 0),
                            stIntTuple_getPosition(k, 1));
                }
            }

            assert(stSortedSet_search(intervals, j) == NULL);
            st_logDebug("Adding in an interval: %s %i %i\n", cA3, start, stop);
            stSortedSet_insert(intervals, j);
            free(cA3);
        }
        bytesRead = benLine(&cA2, &nBytes, fileHandle);
    }
    free(cA2);
    fclose(fileHandle);
    st_logDebug("Finished parsing the bed file: %s\n", filepath);
}

char *stringCommasToSpaces(const char *string) {
    /* swap all commas, ',', for spaces ' '.
     */
    char *cA = stString_copy(string);
    int i = 0;
    for (i = 0; i < strlen(cA); i++) {
        if (cA[i] == ',')
            cA[i] = ' ';
    }
    return cA;
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

void usage() {
    fprintf(stderr, "mafComparator, version %.1f\n", VERSION);
    fprintf(stderr,
            "-a --logLevel : Set the log level. [off, critical, info, debug] "
            "in ascending order.\n");
    fprintf(stderr,
            "-b --mafFile1 : The location of the first MAF file (used to "
            "create sequence name hash.)\n");
    fprintf(stderr, "-c --mafFile2 : The location of the second MAF file\n");
    fprintf(stderr,
            "-d --outputFile : The output XML formatted results file.\n");
    fprintf(stderr,
            "-e --sampleNumber : The number of sample homology tests to perform "
            "(total) [default 1000000].\n");
    fprintf(stderr,
            "-p --printFailures : Print tab-delimited details about failed "
            "tests to stderr.\n");
    fprintf(stderr,
            "-f --bedFiles : The location of bed file(s) used to filter the "
            "pairwise comparisons, separated by commas.\n");
    fprintf(stderr,
            "-g --near : The number of bases in either sequence to allow a "
            "match to slip by.\n");
    fprintf(stderr,
            "-s --seed : an integer used to seed the random number generator "
            "used to perform sampling. If omitted a seed is pseudorandomly "
            "generated.\n");
    fprintf(stderr, "-v --version : Print current version number\n");
    fprintf(stderr, "-h --help : Print this help screen\n");
}

void version() {
    fprintf(stderr, "mafComparator, version %.1f\n", VERSION);
}

int main(int argc, char *argv[]) {
    /*
     * Arguments/options
     */
    char * logLevelString = NULL;
    char * mAFFile1 = NULL;
    char * mAFFile2 = NULL;
    char * outputFile = NULL;
    stHash* intervalsHash =
            stHash_construct3(stHash_stringKey, stHash_stringEqualKey, free,
                    (void(*)(void *)) stSortedSet_destruct);
    char * bedFiles = NULL;
    int32_t sampleNumber = 1000000; // by default do a million samples per pair.
    int32_t randomSeed = (time(NULL) << 16) | (getpid() & 65535); // Likely to be unique
    int32_t near = 0;
    int32_t i;

    ///////////////////////////////////////////////////////////////////////////
    // (0) Parse the inputs handed by genomeCactus.py / setup stuff.
    ///////////////////////////////////////////////////////////////////////////

    while (1) {
        static struct option long_options[] = { 
            { "logLevel", required_argument, 0, 'a' }, 
            { "mafFile1", required_argument, 0, 'b' }, 
            { "mafFile2", required_argument, 0, 'c' }, 
            { "outputFile", required_argument, 0, 'd' }, 
            { "sampleNumber", required_argument, 0, 'e' }, 
            { "printFailed", no_argument, 0, 'p' }, 
            { "version", no_argument, 0, 'v' }, 
            { "help", no_argument, 0, 'h' },
            { "bedFiles", required_argument, 0, 'f' }, 
            { "near", required_argument, 0, 'g' },
            { "seed", required_argument, 0, 's' },
            { 0, 0, 0, 0 } };

        int option_index = 0;

        int key = getopt_long(argc, argv, "a:b:c:d:e:p:v:h:f:g:s:", long_options,
                              &option_index);

        if (key == -1)
            break;

        switch (key) {
            case 'a':
                logLevelString = stString_copy(optarg);
                break;
            case 'b':
                mAFFile1 = stString_copy(optarg);
                break;
            case 'c':
                mAFFile2 = stString_copy(optarg);
                break;
            case 'd':
                outputFile = stString_copy(optarg);
                break;
            case 'v':
                version();
                return 0;
            case 'p':
                VERBOSEFAILURES = 1;
                break;
            case 'e':
                i = sscanf(optarg, "%i", &sampleNumber);
                assert(i == 1);
                break;
            case 's':
                i = sscanf(optarg, "%i", &randomSeed);
                assert(i == 1);
                break;
            case 'h':
                usage();
                return 0;
            case 'f':
                bedFiles = stString_copy(optarg);
                break;
            case 'g':
                i = sscanf(optarg, "%i", &near);
                assert(i == 1);
                break;
            default:
                usage();
                return 1;
        }
    }

    //////////////////////////////////////////////
    //Set up logging
    //////////////////////////////////////////////

    st_setLogLevelFromString(logLevelString);
    st_logDebug("Seeding the random number generator with the value %lo\n", randomSeed);
    st_randomSeed(randomSeed);
    
    ///////////////////////////////////////////////////////////////////////////
    // (0) Check the inputs.
    ///////////////////////////////////////////////////////////////////////////
    if (mAFFile1 == NULL) {
        usage();
        fprintf(stderr, "\nSet --mAFFile1\n");
        exit(2);
    }
    if (mAFFile2 == NULL) {
        usage();
        fprintf(stderr, "\nSet --mAFFile2\n");
        exit(2);
    }
    if (outputFile == NULL) {
        usage();
        fprintf(stderr, "\nSet --outputFile\n");
        exit(2);
    }
    FILE *fileHandle = fopen(mAFFile1, "r");
    if (fileHandle == NULL) {
        if (errno == ENOENT)
            fprintf(stderr, "ERROR, file %s does not exist.\n", mAFFile1);
        else
            fprintf(stderr, "ERROR, unable to open %s\n", mAFFile1);
        exit(1);
    }
    fclose(fileHandle);
    fileHandle = fopen(mAFFile2, "r");
    if (fileHandle == NULL) {
        if (errno == ENOENT)
            fprintf(stderr, "ERROR, file %s does not exist.\n", mAFFile2);
        else
            fprintf(stderr, "ERROR, unable to open %s\n", mAFFile2);
        exit(1);
    }
    fclose(fileHandle);

    //Parse the bed file hashes
    if(bedFiles != NULL) {
        parseBedFiles(bedFiles, intervalsHash);
    }
    else {
        st_logDebug("No bed files specified\n");
    }

    //////////////////////////////////////////////
    //Log (some of) the inputs
    //////////////////////////////////////////////

    st_logInfo("MAF file 1 name : %s\n", mAFFile1);
    st_logInfo("MAF file 2 name : %s\n", mAFFile2);
    st_logInfo("Output stats file : %s\n", outputFile);
    st_logInfo("Bed files parsed : %i\n", stHash_size(intervalsHash));
    st_logInfo("Number of samples %i\n", sampleNumber);

    //////////////////////////////////////////////
    // Create sequence name hashtable from the first MAF file.
    //////////////////////////////////////////////

    stSortedSet *seqNames1 = stSortedSet_construct3((int(*)(const void *, const void *)) strcmp, free);
    populateNames(mAFFile1, seqNames1);
    stSortedSet *seqNames2 = stSortedSet_construct3((int(*)(const void *, const void *)) strcmp, free);
    populateNames(mAFFile2, seqNames2);
    stSortedSet *seqNames = stSortedSet_getIntersection(seqNames1, seqNames2);

    //////////////////////////////////////////////
    //Do comparisons.
    //////////////////////////////////////////////

    if (VERBOSEFAILURES) {
        fprintf(stderr, "# Comparing %s to %s\n", mAFFile1, mAFFile2);
        fprintf(stderr, "# seq1\tabsPos1\torigPos1\tseq2\tabsPos2\torigPos2\n");
    }
    struct avl_table *results_12 = compareMAFs_AB(mAFFile1, mAFFile2, sampleNumber, seqNames, 
                                                  intervalsHash, VERBOSEFAILURES, near);
    if (VERBOSEFAILURES) {
        fprintf(stderr, "# Comparing %s to %s\n", mAFFile2, mAFFile1);
        fprintf(stderr, "# seq1\tabsPos1\torigPos1\tseq2\tabsPos2\torigPos2\n");
    }
    struct avl_table *results_21 = compareMAFs_AB(mAFFile2, mAFFile1, sampleNumber, seqNames, 
                                                  intervalsHash, VERBOSEFAILURES, near);
    fileHandle = fopen(outputFile, "w");
    if (fileHandle == NULL) {
        fprintf(stderr, "ERROR, unable to open %s for writing.\n", outputFile);
        exit(1);
    }

    //////////////////////////////////////////////
    //Report results.
    //////////////////////////////////////////////

    writeXMLHeader(fileHandle);
    if (bedFiles == NULL)
        fprintf(fileHandle, "<alignment_comparisons sampleNumber=\"%i\" "
                "near=\"%i\" seed=\"%u\">\n",
                sampleNumber, near, randomSeed);
    else
        fprintf(fileHandle,
                "<alignment_comparisons sampleNumber=\"%i\" near=\"%i\" seed=\"%u\" "
                "bedFiles=\"%s\">\n",
                sampleNumber, near, randomSeed, bedFiles);
    reportResults(results_12, mAFFile1, mAFFile2, fileHandle, near, seqNames, bedFiles);
    reportResults(results_21, mAFFile2, mAFFile1, fileHandle, near, seqNames, bedFiles);
    fprintf(fileHandle, "</alignment_comparisons>\n");
    fclose(fileHandle);

    ///////////////////////////////////////////////////////////////////////////
    // Clean up.
    ///////////////////////////////////////////////////////////////////////////

    avl_destroy(results_12, (void(*)(void *, void *)) aPair_destruct);
    avl_destroy(results_21, (void(*)(void *, void *)) aPair_destruct);

    free(mAFFile1);
    free(mAFFile2);
    free(outputFile);
    free(logLevelString);

    stSortedSet_destruct(seqNames);
    stSortedSet_destruct(seqNames1);
    stSortedSet_destruct(seqNames2);

    return 0;
}
