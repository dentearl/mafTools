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

#include <getopt.h>

#include "sonLib.h"
#include "common.h"
#include "comparatorAPI.h"

const char *g_version = "version 0.1 July 2012";

void usage(void) {
    fprintf(stderr, "mafPairCounter, %s\n\n", g_version);
    fprintf(stderr, "Usage: $ mafPairCounter --maf=FILE\n\n");
    fprintf(stderr, "Options:\n");
    usageMessage('\0', "maf", "The location of the MAF file. "
                 "The number of pairs contained in the file will be counted and "
                 "reported in stdout.");
    usageMessage('\0', "sequences", "Comma separated list of sequences allowed to be in pairs. " 
                 "To allow all sequences, either specify *every* sequence or don't invoke "
                 "this option. Leaving --sequences off results in all sequences being used.");
    usageMessage('\0', "maf2", "IF Specificied, this is the  location of the second MAF file. "
                 "Using this option causes --sequences option to be ignored. Sequences will "
                 "be discovered by intersection of sequences present in both maf files, pairs "
                 "reported will be from the --maf option.");
    usageMessage('v', "version", "Print current version number.");
    usageMessage('h', "help", "Print this help screen.");
}
void version(void) {
    fprintf(stderr, "mafPairCounter, %s\n", g_version);
}
int parseArgs(int argc, char **argv, char **maf, char **maf2, char **seqList) {
    static const char *optString = "v:h:";
    static const struct option longOpts[] = {
        {"maf", required_argument, 0, 0},
        {"maf2", required_argument, 0, 0},
        {"sequences", required_argument, 0, 0},
        {"version", no_argument, 0, 'v'},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0 }};
    int longIndex = 0;
    int key = getopt_long(argc, argv, optString, longOpts, &longIndex);
    while (key != -1) {
        switch (key) {
        case 0:
            if (strcmp("maf", longOpts[longIndex].name) == 0) {
                *maf = stString_copy(optarg);
                break;
            }
            if (strcmp("maf2", longOpts[longIndex].name) == 0) {
                *maf2 = stString_copy(optarg);
                break;
            }
            if (strcmp("sequences", longOpts[longIndex].name) == 0) {
                *seqList = stString_copy(optarg);
                break;
            }
        case 'v':
            version();
            exit(EXIT_SUCCESS);
            break;
        case 'h':
            usage();
            exit(EXIT_SUCCESS);
            break;
        default:
            usage();
            exit(EXIT_SUCCESS);
            break;
        }
        key = getopt_long(argc, argv, optString, longOpts, &longIndex);
    }
    if (*maf == NULL) {
        usage();
        fprintf(stderr, "\nError, specify --maf\n");
        exit(2);
    }
    FILE *fileHandle = de_fopen(*maf, "r");
    fclose(fileHandle);
    if (*maf2 != NULL) {
        fileHandle = de_fopen(*maf2, "r");
        fclose(fileHandle);
        if (*seqList != NULL) {
            free(seqList);
            seqList = NULL;
        }
    }
    return optind;
}
stHash* buildHash(char *listOfLegitSequences) {
    char *spaceSepFiles = stringCommasToSpaces(listOfLegitSequences);
    char *currentLocation = spaceSepFiles;
    char *currentWord;
    stHash *legitSeqsHash = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, free, free);
    while ((currentWord = stString_getNextWord(&currentLocation)) != NULL) {
        stHash_insert(legitSeqsHash, stString_copy(currentWord), stString_copy(""));
        free(currentWord);
    }
    free(spaceSepFiles);
    return legitSeqsHash;
}
int main(int argc, char **argv) {
    char *maf = NULL;
    char *maf2 = NULL;
    char *listOfLegitSequences = NULL;
    stHash *legitSeqsHash = NULL;
    stHash *maf1Hash = NULL;
    stHash *maf2Hash = NULL;
    parseArgs(argc, argv, &maf, &maf2, &listOfLegitSequences);
    if (listOfLegitSequences != NULL) {
        legitSeqsHash = buildHash(listOfLegitSequences);
    }
    if (maf2 != NULL) {
        // build legitHash by intersection
        maf1Hash = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, free, free);
        maf2Hash = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, free, free);
        populateNames(maf, maf1Hash);
        populateNames(maf2, maf2Hash);
        legitSeqsHash = stHash_getIntersection(maf1Hash, maf2Hash);
    }
    uint64_t numberOfPairs = countPairsInMaf(maf, legitSeqsHash);
    printf("%"PRIu64"\n", numberOfPairs);
    // clean up
    if (legitSeqsHash != NULL) {
        stHash_destruct(legitSeqsHash);
    }
    if (maf1Hash != NULL) {
        stHash_destruct(maf1Hash);
        stHash_destruct(maf2Hash);
    }
    free(maf);
    free(maf2);
    free(listOfLegitSequences);
    return(EXIT_SUCCESS);
}
