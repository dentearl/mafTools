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
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <getopt.h>

#include "cactus.h"
#include "avl.h"
#include "commonC.h"
#include "hashTableC.h"
#include "hashTableC_itr.h"
#include "bioioC.h"

#include "ComparatorAPI.h"
#include "cString.h"

/*
 * The script takes two MAF files and for each ordered trio in the MAFS calculates
 * a predefined number of sample phylogeny trio tests, then reports the statistics in an XML formatted file.
 * It is suitable for running over very large numbers of alignments, because it does not attempt to hold
 * everything in memory, and instead takes a sampling approach.
 */

void freeTrioNames(void *trio) {
    free(((TrioNames *)trio)->speciesA);
    free(((TrioNames *)trio)->speciesB);
    free(((TrioNames *)trio)->speciesC);
    free(trio);
}

void usage() {
    fprintf(stderr, "eval_PhyloComparator, version 0.1\n");
    fprintf(stderr, "\t-a --logLevel : Set the log level\n");
    fprintf(stderr, "\t-b --mAFFile1 : The location of the first MAF file\n");
    fprintf(stderr, "\t-c --mAFFile2 : The location of the second MAF file\n");
    fprintf(stderr, "\t-d --outputFile : The output XML formatted results file.\n");
    fprintf(stderr, "\t-e --sampleNumber : The number of sample phylotrio tests to perform (total).\n");
    fprintf(stderr, "\t-f --trioFile: The location of the trio species file\n");
    fprintf(stderr, "\t-h --help : Print this help screen\n");
}

struct List *parseTrioFile(char *trioFile) {
    FILE *fileHandle = fopen(trioFile, "r");
    int bytesRead;
    int nBytes = 100;
    char *cA;
    int j;

    char *species[3];

    struct List *speciesList = NULL;
    speciesList = constructEmptyList(0, freeTrioNames);

    cA = st_malloc(nBytes + 1);
    bytesRead = benLine(&cA, &nBytes, fileHandle);

    while(bytesRead != -1) {
        if (bytesRead > 0) {
            species[0] = st_malloc(sizeof(char) * (1 + (bytesRead)));
            species[1] = st_malloc(sizeof(char) * (1 + (bytesRead)));
            species[2] = st_malloc(sizeof(char) * (1 + (bytesRead)));
            j = sscanf(cA, "%s\t%s\t%s", species[0], species[1], species[2]);
            if (j != 3) {
                fprintf(stderr, "Invalid triple line '%s' in '%s'\n", cA, trioFile);
                exit(1);
            }

            cStr_lowerCase(species[0]);
            cStr_lowerCase(species[1]);
            cStr_lowerCase(species[2]);
            qsort(species, 3, sizeof(char *), cStr_compare);

            TrioNames *trio = st_malloc(sizeof(TrioNames));
            trio->speciesA = species[0];
            trio->speciesB = species[1];
            trio->speciesC = species[2];
            
            listAppend(speciesList, trio);
        }
        bytesRead = benLine(&cA, &nBytes, fileHandle);
    }
    fclose(fileHandle);

    free(cA);

    return speciesList;
}

int main(int argc, char *argv[]) {
    /*
     * Arguments/options
     */
    char * logLevelString = NULL;
    char * mAFFile1 = NULL;
    char * mAFFile2 = NULL;
    char * outputFile = NULL;
    int32_t sampleNumber = 1000000; // by default do a million samples per pair.
    char * trioFile = NULL;

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
            { "trioFile", required_argument, 0, 'f' },
            { "help", no_argument, 0, 'h' },
            { 0, 0, 0, 0 }
        };

        int option_index = 0;

        int key = getopt_long(argc, argv, "a:b:c:d:e:f:h", long_options, &option_index);

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
            case 'e':
                assert(sscanf(optarg, "%i", &sampleNumber) == 1);
                break;
            case 'f':
                trioFile = stString_copy(optarg);
                break;
            case 'h':
                usage();
                return 0;
            default:
                usage();
                return 1;
        }
    }

    ///////////////////////////////////////////////////////////////////////////
    // (0) Check the inputs.
    ///////////////////////////////////////////////////////////////////////////

    if (argc == 1) {
        usage();
        return 1;
    }

    assert(mAFFile1 != NULL);
    assert(mAFFile2 != NULL);
    assert(outputFile != NULL);
    assert(trioFile != NULL);

    FILE *fileHandle = fopen(mAFFile1, "r");
    if (fileHandle == NULL) {
        fprintf(stderr, "ERROR, unable to open `%s', is path correct?\n", mAFFile1);
        exit(1);
    }
    fclose(fileHandle);
    fileHandle = fopen(mAFFile2, "r");
    if (fileHandle == NULL) {
        fprintf(stderr, "ERROR, unable to open `%s', is path correct?\n", mAFFile2);
        exit(1);
    }
    fclose(fileHandle);
    fileHandle = fopen(trioFile, "r");
    if (fileHandle == NULL) {
        fprintf(stderr, "ERROR, unable to open `%s', is path correct?\n", trioFile);
        exit(1);
    }
    fclose(fileHandle);

    //////////////////////////////////////////////
    //Set up logging
    //////////////////////////////////////////////

    st_setLogLevelFromString(logLevelString);

    //////////////////////////////////////////////
    //Log (some of) the inputs
    //////////////////////////////////////////////

    st_logInfo("MAF file 1 name : %s\n", mAFFile1);
    st_logInfo("MAF file 2 name : %s\n", mAFFile2);
    st_logInfo("Trio species name : %s\n", trioFile);
    st_logInfo("Output stats file : %s\n", outputFile);

    /* Parse the trioFile triples into a list */
    struct List *speciesList = parseTrioFile(trioFile);

    //////////////////////////////////////////////
    // Create hashtable for the first MAF file.
    //////////////////////////////////////////////

    struct hashtable *seqNameHash;
    seqNameHash = create_hashtable(256, hashtable_stringHashKey, hashtable_stringEqualKey, free, free);
    populateNameHash(mAFFile1, seqNameHash);

    // TODO: Check if query species are in maf file

    //////////////////////////////////////////////
    //Do comparisons.
    //////////////////////////////////////////////

    struct avl_table *results_12 = compareMAFs_AB_Trio(mAFFile1, mAFFile2, sampleNumber, seqNameHash, speciesList);
    struct avl_table *results_21 = compareMAFs_AB_Trio(mAFFile2, mAFFile1, sampleNumber, seqNameHash, speciesList);

    fileHandle = fopen(outputFile, "w");
    if (fileHandle == NULL) {
        fprintf(stderr, "ERROR, unable to open `%s' for writing.\n", outputFile);
        exit(1);
    }

    fprintf(fileHandle, "<trio_comparisons sampleNumber=\"%i\">\n", sampleNumber);
    reportResultsTrio(results_12, mAFFile1, mAFFile2, fileHandle);
    reportResultsTrio(results_21, mAFFile2, mAFFile1, fileHandle);
    fprintf(fileHandle, "</trio_comparisons>\n");
    fclose(fileHandle);

    ///////////////////////////////////////////////////////////////////////////
    // Clean up.
    ///////////////////////////////////////////////////////////////////////////

    free(mAFFile1);
    free(mAFFile2);
    free(outputFile);
    free(trioFile);
    free(logLevelString);

    hashtable_destroy(seqNameHash, 1, 1);

    avl_destroy(results_12, (void (*)(void *, void *))aTrio_destruct);
    avl_destroy(results_21, (void (*)(void *, void *))aTrio_destruct);

    destructList(speciesList);

    return 0;
}
