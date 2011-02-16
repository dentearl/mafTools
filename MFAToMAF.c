#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <getopt.h>

#include "bioioC.h"
#include "commonC.h"
#include "sonLib.h"

#include "eTreeExtras.h"

/*
 * Takes an MFA file and creates a MAF file containing a single MAF block,
 * representing the alignment.
 */

void usage() {
    fprintf(stderr, "eval_MFAtoMAF, version 0.1\n");
    fprintf(stderr, "\t-a --logLevel : Set the log level\n");
    fprintf(stderr, "\t-b --mFAFile : The location of the MFA file\n");
    fprintf(stderr, "\t-d --outputFile : The output file MAF file.\n");
    fprintf(stderr, "\t-t --treeFile: The location of the tree file\n");
    fprintf(stderr, "\t-h --help : Print this help screen\n");
}

int main(int argc, char *argv[]) {
    /*
     * Arguments/options
     */
    char *logLevelString = NULL;
    char *mFAFile = NULL;
    char *outputFile = NULL;
    char *treeFile = NULL;

    ///////////////////////////////////////////////////////////////////////////
    // (0) Parse the inputs handed by genomeCactus.py / setup stuff.
    ///////////////////////////////////////////////////////////////////////////

    while (1) {
        static struct option long_options[] = {
            { "logLevel", required_argument, 0, 'a' },
            { "mFAFile", required_argument, 0, 'b' },
            { "outputFile", required_argument, 0, 'd' },
            { "treeFile", optional_argument, 0, 't' },
            { "help", no_argument, 0, 'h' },
            { 0, 0, 0, 0 }
        };

        int option_index = 0;

        int key = getopt_long(argc, argv, "a:b:d:t:h", long_options, &option_index);

        if (key == -1) {
            break;
        }

        switch (key) {
            case 'a':
                logLevelString = stString_copy(optarg);
                break;
            case 'b':
                mFAFile = stString_copy(optarg);
                break;
            case 'd':
                outputFile = stString_copy(optarg);
                break;
            case 't':
                treeFile = stString_copy(optarg);
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
        exit(1);
    }

    assert(mFAFile != NULL);
    assert(outputFile != NULL);

    //////////////////////////////////////////////
    //Set up logging
    //////////////////////////////////////////////

    if (logLevelString != NULL && strcmp(logLevelString, "INFO") == 0) {
        st_setLogLevel(ST_LOGGING_INFO);
    }
    if (logLevelString != NULL && strcmp(logLevelString, "DEBUG") == 0) {
        st_setLogLevel(ST_LOGGING_DEBUG);
    }

    //////////////////////////////////////////////
    //Log (some of) the inputs
    //////////////////////////////////////////////

    st_logInfo("MFA file  name : %s\n", mFAFile);
    st_logInfo("Output MAF file : %s\n", outputFile);
    st_logInfo("Tree file name: %s\n", treeFile == NULL ? "null" : treeFile);

    //////////////////////////////////////////////
    //Get the MFA alignment
    //////////////////////////////////////////////

    //get the alignment
    struct List *sequences = constructEmptyList(0, free);
    struct List *seqLengths = constructEmptyList(0, (void (*)(void *))destructInt);
    struct List *fastaNames = constructEmptyList(0, free);
    FILE *fileHandle = fopen(mFAFile, "r");
    if (fileHandle == NULL) {
        usage();
        exit(1);
    }
    fastaRead(fileHandle, sequences, seqLengths, fastaNames);
    fclose(fileHandle);

    //////////////////////////////////////////////
    //Get the tree alignment
    //////////////////////////////////////////////

    ETree *tree = NULL;
    LeafPtrArray *leafArray = NULL;

    int32_t leafCount = 0;
    if (treeFile != NULL) {
        tree = eTreeX_getTreeFromFile(treeFile);

        eTreeX_postOrderTraversal(tree, eTreeX_countLeaves, &leafCount);

        leafArray = eTreeX_constructLeafPtrArray(leafCount);
        eTreeX_postOrderTraversal(tree, eTreeX_getLeafArray, (void *) leafArray);
    }

    //////////////////////////////////////////////
    //Write the MFA alignment.
    //////////////////////////////////////////////

    fileHandle = fopen(outputFile, "w");
    //write the header.
    fprintf(fileHandle, "##maf version=1 scoring=NULL\n");
    fprintf(fileHandle, "# converted_from_MFA\n\n");

    //write the score line
    char *treeString = NULL;
    if (treeFile != NULL) {
        treeString = eTree_getNewickTreeString(tree);
        fprintf(fileHandle, "a score=0 tree=\"%s\"\n", treeString);
    }
    else {
        fprintf(fileHandle, "a score=0\n");
        leafCount = sequences->length;
    }

    //write the alignment
    int32_t i, j;
    int32_t ii;
    const char *label;
    for (ii=0; ii<leafCount; ii++) {
        if (treeFile != NULL) {
            label = eTree_getLabel((ETree *) leafArray->ptrArray[ii]);

            /* Do a brute force search to find the appropriate sequence that matches "label" */
            for (i=0; i<sequences->length; i++) {
                char *fastaHeader = fastaNames->list[i];
                char *sequenceName = st_malloc(sizeof(char) *(1 + strlen(fastaHeader)));
                sscanf(fastaHeader, "%s", sequenceName); //take the sequence name to be the first word of the sequence.
                if (strcmp(label, sequenceName) == 0) {
                    free(sequenceName);
                    break;
                }
                free(sequenceName);
            }
        }
        else {
            i = ii;
        }

        char *sequence = sequences->list[i];
        int32_t seqLength = *((int32_t *)seqLengths->list[i]);
        assert(seqLength == (int32_t)strlen(sequence));
        char *fastaHeader = fastaNames->list[i];
        char *sequenceName = st_malloc(sizeof(char) *(1 + strlen(fastaHeader)));
        sscanf(fastaHeader, "%s", sequenceName); //take the sequence name to be the first word of the sequence.
        int32_t length = 0;
        for (j=0; j<(int32_t)strlen(sequence); j++) {
            if (sequence[j] != '-') {
                length++;
            }
        }
        fprintf(fileHandle, "s\t%s\t%i\t%i\t%s\t%i\t%s\n", sequenceName, 0, length, "+", length, sequence);
        free(sequenceName);
    }

    fclose(fileHandle);

    //////////////////////////////////////////////
    //Clean up.
    //////////////////////////////////////////////

    free(mFAFile);
    free(outputFile);
    free(treeFile);

    if (treeFile != NULL) {
        eTree_destruct(tree);
        free(treeString);
        eTreeX_destructLeafPtrArray(leafArray);
    }

    destructList(sequences);
    destructList(seqLengths);
    destructList(fastaNames);

    return 0;
}
