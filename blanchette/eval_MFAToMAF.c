#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <getopt.h>

#include "commonC.h"
#include "bioioC.h"

/*
 * Takes an MFA file and creates a MAF file containing a single MAF block,
 * representing the alignment.
 */

void usage() {
	fprintf(stderr, "eval_MFAtoMAF, version 0.1\n");
	fprintf(stderr, "-a --logLevel : Set the log level\n");
	fprintf(stderr, "-b --mFAFile : The location of the MFA file\n");
	fprintf(stderr, "-d --outputFile : The output file MAF file.\n");
	fprintf(stderr, "-h --help : Print this help screen\n");
}

int main(int argc, char *argv[]) {
	/*
	 * Arguments/options
	 */
	char * logLevelString = NULL;
	char * mFAFile = NULL;
	char * outputFile = NULL;

	///////////////////////////////////////////////////////////////////////////
	// (0) Parse the inputs handed by genomeCactus.py / setup stuff.
	///////////////////////////////////////////////////////////////////////////

	while(1) {
		static struct option long_options[] = {
			{ "logLevel", required_argument, 0, 'a' },
			{ "mFAFile", required_argument, 0, 'b' },
			{ "outputFile", required_argument, 0, 'd' },
			{ "help", no_argument, 0, 'h' },
			{ 0, 0, 0, 0 }
		};

		int option_index = 0;

		int key = getopt_long(argc, argv, "a:b:d:h", long_options, &option_index);

		if(key == -1) {
			break;
		}

		switch(key) {
			case 'a':
				logLevelString = stringCopy(optarg);
				break;
			case 'b':
				mFAFile = stringCopy(optarg);
				break;
			case 'd':
				outputFile = stringCopy(optarg);
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

	assert(mFAFile != NULL);
	assert(outputFile != NULL);

	//////////////////////////////////////////////
	//Set up logging
	//////////////////////////////////////////////

	if(logLevelString != NULL && strcmp(logLevelString, "INFO") == 0) {
		setLogLevel(LOGGING_INFO);
	}
	if(logLevelString != NULL && strcmp(logLevelString, "DEBUG") == 0) {
		setLogLevel(LOGGING_DEBUG);
	}

	//////////////////////////////////////////////
	//Log (some of) the inputs
	//////////////////////////////////////////////

	logInfo("MFA file  name : %s\n", mFAFile);
	logInfo("Output MAF file : %s\n", outputFile);

	//////////////////////////////////////////////
	//Get the MFA alignment
	//////////////////////////////////////////////

	//get the alignment
	struct List *sequences = constructEmptyList(0, free);
	struct List *seqLengths = constructEmptyList(0, (void (*)(void *))destructInt);
	struct List *fastaNames = constructEmptyList(0, free);
	FILE *fileHandle = fopen(mFAFile, "r");
	fastaRead(fileHandle, sequences, seqLengths, fastaNames);
	fclose(fileHandle);

	//////////////////////////////////////////////
	//Write the MFA alignment.
	//////////////////////////////////////////////

	fileHandle = fopen(outputFile, "w");
	//write the header.
	fprintf(fileHandle, "##maf version=1 scoring=NULL\n");
	fprintf(fileHandle, "# converted_from_MFA\n\n");

	//write the score line
	fprintf(fileHandle, "a score=0\n");

	//write the alignment
	int32_t i, j;
	for(i=0; i<sequences->length; i++) {
		char *sequence = sequences->list[i];
		int32_t seqLength = *((int32_t *)seqLengths->list[i]);
		assert(seqLength == (int32_t)strlen(sequence));
		char *fastaHeader = fastaNames->list[i];
		char *sequenceName = malloc(sizeof(char) *(1 + strlen(fastaHeader)));
		sscanf(fastaHeader, "%s", sequenceName); //take the sequence name to be the first word of the sequence.
		int32_t length = 0;
		for(j=0; j<(int32_t)strlen(sequence); j++) {
			if(sequence[j] != '-') {
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

	destructList(sequences);
	destructList(seqLengths);
	destructList(fastaNames);

	return 0;
}

