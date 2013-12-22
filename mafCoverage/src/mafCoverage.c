/*
 * Copyright (C) 2011-2013 by
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
#include <getopt.h>
#include <math.h>  // ceil()
#include <stdbool.h>
#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>
#include <stdlib.h>
#include <string.h>
#include "common.h"
#include "sharedMaf.h"
#include "mafCoverage.h"
#include "mafCoverageAPI.h"
#include "buildVersion.h"
#include "sonLib.h"

static char *mafFileName = NULL;
static stSet *speciesOrChromosomeNames = NULL;
static bool nCoverage = 0, identity = 0, ignoreSpecies = 0;

const char *g_version = "version 0.1 May 2013";
uint64_t getRegionSize(char *seq1, stHash *intervalsHash);

void version(void) {
    fprintf(stderr, "mafCoverage, %s\nbuild: %s, %s, %s\n\n", g_version, g_build_date, g_build_git_branch, g_build_git_sha);
}

void usage(void) {
    version();
    fprintf(stderr, "Usage: mafCoverage [maf file] \n\n"
        "Reports the pairwise (n-)coverage between a specified genome and all other genomes in the given maf, using a tab delimited format.\n"
        "Output table format has fields: querySpecies\ttargetSpecies\tlengthOfQueryGenome\tcoverage\tn-coverages (if specified)\n"
        "For a pair of genomes A and B, the coverage of B on A is the proportion of sites in A that align to a base in B.\n"
        "The n-coverage of B on A is the proportion of sites in A that align to n or more sites in B.\n");
    fprintf(stderr, "Options: \n");
    usageMessage('h', "help", "show this help message and exit.");
    usageMessage('m', "maf", "path to maf file.");
    usageMessage('s', "speciesOrChr",
            "species or species.chromosome name, e.g. `hg19' or 'hg19.chr1', if not specified reports results for every possible species."
                "wildcard at the end.");
    usageMessage('n', "nCoverage", "report all n-coverages, for 1 <= n <= 128 instead of just for n=1 (the default).");
    usageMessage('i', "identity", "report coverage of identical bases.");
    usageMessage('l', "logLevel", "Set logging level, either 'CRITICAL'/'INFO'/'DEBUG'.");
    usageMessage('a', "ignoreSpecies", "Do all chromosomes-against-all-chromosomes coverage.");
    exit(EXIT_FAILURE);
}

static void parseOptions(int argc, char **argv) {
    int c;
    speciesOrChromosomeNames = stSet_construct3(stHash_stringKey, stHash_stringEqualKey, free);
    while (1) {
        static struct option longOptions[] = { { "help", no_argument, 0, 'h' }, { "maf", required_argument, 0, 'm' }, { "speciesOrChr",
                required_argument, 0, 's' }, { "nCoverage", no_argument, 0, 'n' }, { "identity", no_argument, 0, 'i' }, { "logLevel",
                required_argument, 0, 'l' }, { "ignoreSpecies", no_argument, 0, 'a' }, { 0, 0, 0, 0 } };
        int longIndex = 0;
        c = getopt_long(argc, argv, "m:s:hnl:a", longOptions, &longIndex);
        if (c == -1)
            break;
        switch (c) {
            case 's':
                stSet_insert(speciesOrChromosomeNames, stString_copy(optarg));
                break;
            case 'm':
                mafFileName = stString_copy(optarg);
                break;
            case 'n':
                nCoverage = 1;
                break;
            case 'i':
                identity = 1;
                break;
            case 'l':
                st_setLogLevelFromString(optarg);
                break;
            case 'h':
                usage();
                break;
            case 'a':
                ignoreSpecies = 1;
                break;
            default:
                abort();
        }
    }
    //Check we have the essentials.
    if (mafFileName == NULL) {
        fprintf(stderr, "Error, specify --maf\n");
        usage();
    }
    // Check there's nothing left over on the command line
    if (optind < argc) {
        fprintf(stderr, "Unexpected input arguments\n");
        usage();
    }
}

int main(int argc, char **argv) {
    parseOptions(argc, argv);
    //Work out the structure of the chromosomes of the query sequence
    stHash *sequenceNamesToSequenceSizes = getMapOfSequenceNamesToSizesFromMaf(mafFileName);
    stHashIterator *sequenceNameIt = stHash_getIterator(sequenceNamesToSequenceSizes);
    char *sequenceName;
    while ((sequenceName = stHash_getNext(sequenceNameIt)) != NULL) {
        st_logDebug("Got a sequence name: %s with length %" PRId64 "\n", sequenceName, stIntTuple_get(stHash_search(sequenceNamesToSequenceSizes, sequenceName), 0));
    }
    stHash_destructIterator(sequenceNameIt);
    stList *sequenceNames = stHash_getKeys(sequenceNamesToSequenceSizes);
    //If the species/chr name is not specified then replace with all possible species name.
    if (stSet_size(speciesOrChromosomeNames) == 0) {
        st_logInfo("As no species name was specified, using all possible species names\n");
        stSet_destruct(speciesOrChromosomeNames);
        speciesOrChromosomeNames = getSpeciesNames(sequenceNames, ignoreSpecies);
    } else { //Sanity checks on the input species/chr
        stList *names = stSet_getList(speciesOrChromosomeNames);
        assert(stList_length(names) == 1);
        char *speciesOrChrName = stList_get(names, 0);
        stList_destruct(names);
        if (ignoreSpecies) {
            if (stHash_search(sequenceNamesToSequenceSizes, speciesOrChrName) == NULL) {
                st_errAbort("Chromosome name not recognised (perhaps you gave a species name but have specified --ignoreSpecies?): %s\n",
                        speciesOrChrName);
            }
        } else {
            stSet *speciesNames = getSpeciesNames(sequenceNames, ignoreSpecies);
            if (stSet_search(speciesNames, speciesOrChrName) == NULL && stHash_search(sequenceNamesToSequenceSizes, speciesOrChrName)
                    == NULL) {
                st_errAbort("Species or chr name name not recognised: %s\n", speciesOrChrName);
            }
            stSet_destruct(speciesNames);
        }
    }
    //Print header
    nGenomeCoverage_reportHeader(stdout, nCoverage);
    //For each of the chosen species calculate species
    stSetIterator *speciesOrChrNamesIt = stSet_getIterator(speciesOrChromosomeNames);
    char *speciesOrChrName;
    while ((speciesOrChrName = stSet_getNext(speciesOrChrNamesIt)) != NULL) {
        st_logInfo("Computing the coverages for species/chr: %s\n", speciesOrChrName);
        //Build the coverage data structure
        NGenomeCoverage *nGC = nGenomeCoverage_construct(sequenceNamesToSequenceSizes, speciesOrChrName, ignoreSpecies);
        nGenomeCoverage_populate(nGC, mafFileName, identity);
        //Report
        nGenomeCoverage_report(nGC, stdout, nCoverage);
        //cleanup loop
        nGenomeCoverage_destruct(nGC);
    }
    //Cleanup
    stList_destruct(sequenceNames);
    stSet_destructIterator(speciesOrChrNamesIt);
    stHash_destruct(sequenceNamesToSequenceSizes);
    stSet_destruct(speciesOrChromosomeNames);
    free(mafFileName);
    // while(1);
    return EXIT_SUCCESS;
}
