#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <getopt.h>
#include <stdarg.h>

//#include "cactus.h"
#include "avl.h"
#include "commonC.h"
#include "hashTableC.h"
#include "hashTableC_itr.h"
#include "bioioC.h"

#include "ComparatorAPI.h"

float VERSION = 0.1;
int32_t ULTRAVERBOSE = 0;

/*
 * The script takes two MAF files and for each ordered pair of sequences in the MAFS calculates
 * a predefined number of sample homology tests (see below), then reports the statistics in an XML formatted file.
 * It is suitable for running over very large numbers of alignments, because it does not attempt to hold
 * everything in memory, and instead takes a sampling approach.
 *
 * For two sets of pairwise alignments, A and B, a homology test is defined as follows.
 * Pick a pair of aligned positions in A, called a homology pair - the AB homology test returns true if the
 * pair is in B, otherwise it returns false. The set of possible homology tests for the ordered pair (A, B)
 * is not necessarily equivalent to the set of possible (B, A) homology tests.
 * We call the proportion of true tests as a percentage of the total of a set of homology tests C from (A, B)  A~B.
 *
 * If A is the set of true pairwise alignments and B the predicted set of alignments then A~B (over large enough
 * C), is a proxy to sensitivity of B in predicted the set of correctly aligned pairs in A. Conversely B~A (over large enough C)
 * is a proxy to the specificity of the aligned pairs in B with respect to the set of correctly aligned pairs in A.
 */

void parseBedFile(const char *cA, stHash *intervalsHash) {
    FILE *fileHandle = fopen(cA, "r");

    int nBytes = 100;
    char *cA2 = st_malloc(nBytes + 1);
    int32_t bytesRead = benLine (&cA2, &nBytes, fileHandle);

    //read through lines until reaching a line starting with an 's':
    while(bytesRead != -1) {
        if (bytesRead > 0) {
            int32_t start, stop;
            char *cA3 = stString_copy(cA2);
            int32_t i = sscanf(cA2, "%s %i %i", cA3, &start, &stop);
            assert(i == 3);
            stSortedSet *intervals = stHash_search(intervalsHash, cA3);
            if(intervals == NULL) {
                intervals = stSortedSet_construct3((int (*)(const void *, const void *))stIntTuple_cmpFn, (void (*)(void *))stIntTuple_destruct);
                stHash_insert(intervalsHash, stString_copy(cA3), intervals);
            }
            free(cA3);
            stIntTuple *j = stIntTuple_construct(2, start, stop);
            assert(stSortedSet_search(intervals, j) == NULL);
            stSortedSet_insert(intervals, j);
        }
        bytesRead = benLine (&cA2, &nBytes, fileHandle);
    }

    free(cA2);
    fclose(fileHandle);
}

void parseBedFiles(const char *cA, stHash *bedFileHash) {
    char *cA2 = stString_copy(cA);
    char *cA3 = cA2;
    char *cA4;
    while((cA4 = stString_getNextWord(&cA3)) != NULL) {
        parseBedFile(cA4, bedFileHash);
        free(cA4);
    }
    free(cA2);
}

void usage() {
   fprintf(stderr, "MAFComparator, version %.1f\n", VERSION);
   fprintf(stderr, "-a --logLevel : Set the log level\n");
   fprintf(stderr, "-b --mAFFile1 : The location of the first MAF file (used to create sequence name hash.)\n");
   fprintf(stderr, "-c --mAFFile2 : The location of the second MAF file\n");
   fprintf(stderr, "-d --outputFile : The output XML formatted results file.\n");
   fprintf(stderr, "-e --sampleNumber : The number of sample homology tests to perform (total) [default 1000000].\n");
   fprintf(stderr, "-u --ultraVerbose : Print tab-delimited details about failed tests to stderr.\n");
   fprintf(stderr, "-v --version : Print current version number\n");
   fprintf(stderr, "-h --help : Print this help screen\n");
   fprintf(stderr, "-f --bedFiles : The location of bed file used to filter the pairwise comparisons.\n");
   fprintf(stderr, "-g --near : The number of bases in either sequence to allow a match to slip by.\n");
}

void version() {
   fprintf(stderr, "MAFComparator, version %.1f\n", VERSION);
}

int main(int argc, char *argv[]) {
    /*
     * Arguments/options
     */
    char * logLevelString = NULL;
    char * mAFFile1 = NULL;
    char * mAFFile2 = NULL;
    char * outputFile = NULL;
    stHash* intervalsHash = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, free, (void (*)(void *))stSortedSet_destruct);
    int32_t sampleNumber = 1000000; // by default do a million samples per pair.
    int32_t near = 0;
    int32_t i;

    ///////////////////////////////////////////////////////////////////////////
    // (0) Parse the inputs handed by genomeCactus.py / setup stuff.
    ///////////////////////////////////////////////////////////////////////////

    while(1) {
        static struct option long_options[] = {
            { "logLevel", required_argument, 0, 'a' },
            { "mAFFile1", required_argument, 0, 'b' },
            { "mAFFile2", required_argument, 0, 'c' },
            { "outputFile", required_argument, 0, 'd' },
            { "sampleNumber", required_argument, 0, 'e' },
            { "ultraVerbose", no_argument, 0, 'u'},
            { "version", no_argument, 0, 'v'},
            { "help", no_argument, 0, 'h' },
            { "bedFiles", required_argument, 0, 'f' },
            { "near", required_argument, 0, 'g' },
            { 0, 0, 0, 0 }
        };

        int option_index = 0;

        int key = getopt_long(argc, argv, "a:b:c:d:e:u:v:hf:g:", long_options, &option_index);

        if(key == -1) {
            break;
        }

        switch(key) {
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
            case 'u':
                ULTRAVERBOSE = 1;
                break;
            case 'e':
                i = sscanf(optarg, "%i", &sampleNumber);
                assert(i == 1);
                break;
            case 'h':
                usage();
                return 0;
            case 'f':
                //Parse the bed file hashes
                parseBedFiles(optarg, intervalsHash);
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

    ///////////////////////////////////////////////////////////////////////////
    // (0) Check the inputs.
    ///////////////////////////////////////////////////////////////////////////
    if (mAFFile1 == NULL){
       usage();
       fprintf(stderr, "\nSet --mAFFile1 .\n");
       exit(2);
    }
    if (mAFFile2 == NULL){
       usage();
       fprintf(stderr, "\nSet --mAFFile2 .\n");
       exit(2);
    }
    if (outputFile == NULL){
       usage();
       fprintf(stderr, "\nSet --outputFile .\n");
       exit(2);
    }
    FILE *fileHandle = fopen(mAFFile1, "r");
    if(fileHandle == NULL){
       fprintf(stderr, "ERROR, unable to open `%s', is path correct?\n", mAFFile1);
       exit(1);
    }
    fclose(fileHandle);
    fileHandle = fopen(mAFFile2, "r");
    if(fileHandle == NULL){
       fprintf(stderr, "ERROR, unable to open `%s', is path correct?\n", mAFFile2);
       exit(1);
    }
    fclose(fileHandle);

    //////////////////////////////////////////////
    //Set up logging
    //////////////////////////////////////////////

    if(logLevelString != NULL && strcmp(logLevelString, "INFO") == 0) {
        st_setLogLevel(ST_LOGGING_INFO);
    }
    if(logLevelString != NULL && strcmp(logLevelString, "DEBUG") == 0) {
        st_setLogLevel(ST_LOGGING_DEBUG);
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

    stSortedSet *seqNames1 = stSortedSet_construct3((int (*)(const void *, const void *))strcmp, free);
    populateNames(mAFFile1, seqNames1);
    stSortedSet *seqNames2 = stSortedSet_construct3((int (*)(const void *, const void *))strcmp, free);
    populateNames(mAFFile2, seqNames2);
    stSortedSet *seqNames = stSortedSet_getIntersection(seqNames1, seqNames2);

    //////////////////////////////////////////////
    //Do comparisons.
    //////////////////////////////////////////////

    if (ULTRAVERBOSE){
       fprintf(stderr, "# Comparing %s to %s\n", mAFFile1, mAFFile2);
       fprintf(stderr, "# seq1\tabsPos1\torigPos1\tseq2\tabsPos2\torigPos2\n");
    }
    struct avl_table *results_12 = compareMAFs_AB(mAFFile1, mAFFile2, sampleNumber, seqNames, intervalsHash, ULTRAVERBOSE, near);
    if (ULTRAVERBOSE){
       fprintf(stderr, "# Comparing %s to %s\n", mAFFile2, mAFFile1);
       fprintf(stderr, "# seq1\tabsPos1\torigPos1\tseq2\tabsPos2\torigPos2\n");
    }
    struct avl_table *results_21 = compareMAFs_AB(mAFFile2, mAFFile1, sampleNumber, seqNames, intervalsHash, ULTRAVERBOSE, near);
    fileHandle = fopen(outputFile, "w");
    if(fileHandle == NULL){
       fprintf(stderr, "ERROR, unable to open `%s' for writing.\n", outputFile);
       exit(1);
    }

    //////////////////////////////////////////////
    //Report results.
    //////////////////////////////////////////////
    
    writeXMLHeader( fileHandle );
    fprintf(fileHandle, "<alignment_comparisons sampleNumber=\"%i\">\n", sampleNumber);
    reportResults(results_12, mAFFile1, mAFFile2, fileHandle, near, seqNames);
    reportResults(results_21, mAFFile2, mAFFile1, fileHandle, near, seqNames);
    fprintf(fileHandle, "</alignment_comparisons>\n");
    fclose(fileHandle);

    ///////////////////////////////////////////////////////////////////////////
    // Clean up.
    ///////////////////////////////////////////////////////////////////////////

    avl_destroy(results_12, (void (*)(void *, void *))aPair_destruct);
    avl_destroy(results_21, (void (*)(void *, void *))aPair_destruct);

    free(mAFFile1);
    free(mAFFile2);
    free(outputFile);
    free(logLevelString);

    stSortedSet_destruct(seqNames);
    stSortedSet_destruct(seqNames1);
    stSortedSet_destruct(seqNames2);

    return 0;
}
