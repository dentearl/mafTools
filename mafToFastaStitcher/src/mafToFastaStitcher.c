/* 
 * Copyright (C) 2012 by 
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

#include <ctype.h> // mac os x toupper()
#include <getopt.h>
#include <inttypes.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include "common.h"
#include "CuTest.h"
#include "sharedMaf.h"
#include "sonLib.h"
#include "mafToFastaStitcher.h"
#include "mafToFastaStitcherAPI.h"
#include "buildVersion.h"

const char *g_version = "v0.1 Oct 2012";

void version(void);
void usage(void);

void parseOptions(int argc, char **argv, options_t *options) {
    int c;
    bool setMafName = false, setSeqNames = false, setOutName = false, 
        setBreakpointPenalty = false, setInterstitialSequence = false,
        setReference = false;
    size_t i;
    while (1) {
        static struct option long_options[] = {
            {"debug", no_argument, 0, 'd'},
            {"verbose", no_argument, 0, 'v'},
            {"help", no_argument, 0, 'h'},
            {"version", no_argument, 0, 0},
            {"maf",  required_argument, 0, 0},
            {"seqs",  required_argument, 0, 0},
            {"outMfa",  required_argument, 0, 0},
            {"outMaf",  required_argument, 0, 0},
            {"breakpointPenalty",  required_argument, 0, 0},
            {"interstitialSequence",  required_argument, 0, 0},
            {"referenceSequence",  required_argument, 0, 0},
            {0, 0, 0, 0}
        };
        int option_index = 0;
        c = getopt_long(argc, argv, "d:m:s:h:v:t", long_options, &option_index);
        if (c == -1)
            break;
        switch (c) {
        case 0:
            if (strcmp("version", long_options[option_index].name) == 0) {
                version();
                exit(EXIT_SUCCESS);
            }
            if (strcmp("maf", long_options[option_index].name) == 0) {
                setMafName = true;
                options->maf = stString_copy(optarg);
                break;
            }
            if (strcmp("seqs", long_options[option_index].name) == 0) {
                setSeqNames = true;
                options->seqs = stString_copy(optarg);
                break;
            }
            if (strcmp("outMaf", long_options[option_index].name) == 0) {
                setOutName = true;
                options->outMaf = stString_copy(optarg);
                break;
            }
            if (strcmp("outMfa", long_options[option_index].name) == 0) {
                setOutName = true;
                options->outMfa = stString_copy(optarg);
                break;
            }
            if (strcmp("breakpointPenalty", long_options[option_index].name) == 0) {
                setBreakpointPenalty = true;
                i = sscanf(optarg, "%" PRIu32, &(options->breakpointPenalty));
                assert(i == 1);
                break;
            }
            if (strcmp("interstitialSequence", long_options[option_index].name) == 0) {
                setInterstitialSequence = true;
                i = sscanf(optarg, "%" PRIu32, &(options->interstitialSequence));
                assert(i == 1);
                break;
            }
            if (strcmp("referenceSequence", long_options[option_index].name) == 0) {
                setReference = true;
                options->reference = stString_copy(optarg);
                break;
            }
            break;
        case 'v':
            g_verbose_flag++;
            break;
        case 'd':
            g_debug_flag = 1;
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
    if (!setOutName) {
        fprintf(stderr, "specify --outMaf or --outMfa\n");
        usage();
    }
    if (!setSeqNames) {
        fprintf(stderr, "specify --seqs\n");
        usage();
    }
    if (!setBreakpointPenalty) {
        fprintf(stderr, "specify --breakpointPenalty\n");
        usage();
    }
    if (!setInterstitialSequence) {
        fprintf(stderr, "specify --interstitialSequence\n");
        usage();
    }
    // Check there's nothing left over on the command line 
    if (optind < argc) {
        char *errorString = st_malloc(kMaxStringLength);
        strcpy(errorString, "Unexpected arguments:");
        while (optind < argc) {
            strcat(errorString, " ");
            strcat(errorString, argv[optind++]);
        }
        fprintf(stderr, "%s\n", errorString);
        free(errorString);
        usage();
    }
}
void version(void) {
    fprintf(stderr, "mafToFastaStitcher, %s\nbuild: %s, %s, %s\n\n", g_version, g_build_date, 
            g_build_git_branch, g_build_git_sha);
}
void usage(void) {
    version();
    fprintf(stderr, "Usage: mafToFastaStitcher --maf mafFile.maf --seqs seq1.fa,seq2.fa[,...] --breakpointPenalty 5 --interstitialSequence 20 --outMfa output.mfa \n\n"
            "\n\n");
    fprintf(stderr, "Options: \n");
    usageMessage('h', "help", "show this message and exit.");
    usageMessage('m', "maf", "path to the maf file.");
    usageMessage('\0', "seqs", "comma separated list of fasta sequences. each fasta may contain multiple entries. all sequences in the input alignment must be accounted for with an element in a fasta.");
    usageMessage('\0', "outMfa", "multiple sequence fasta output file.");
    usageMessage('\0', "breakpointPenalty", "number of `N' characters to insert into a sequence when a breakpoint is detected.");
    usageMessage('\0', "interstitialSequence", "maximum length of interstitial sequence to be added (from a fasta) into the fasta before a breakpoint is declared and the <code>--breakpointPenalty</code> number of <code>N</code>'s is added instead.");
    usageMessage('\0', "outMaf", "multiple alignment format output file.");
    usageMessage('\0', "reference", "optional. The name of the reference sequence. All intervening reference sequence between the first and last block of the input --maf will be read out in the output.");
    usageMessage('v', "verbose", "turns on verbose output.");
    exit(EXIT_FAILURE);
}
int main(int argc, char **argv) {
    options_t *options = options_construct();
    stHash *sequenceHash = NULL; // keyed on fasta headers, valued with mtfseq_t pointers
    stHash *alignmentHash = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, free, destroyRow); // keyed on species names, valued with row_t pointers
    stList *rowOrder = stList_construct3(0, free); // when adding keys to alignmentHash, append to this list
    parseOptions(argc, argv, options);
    // read fastas, populate sequenceHash
    de_verbose("Creating sequence hash.\n");
    sequenceHash = createSequenceHash(options->seqs);
    mafFileApi_t *mfapi = maf_newMfa(options->maf, "r");
    de_verbose("Creating alignment hash.\n");
    buildAlignmentHash(mfapi, alignmentHash, sequenceHash, rowOrder, options);
    if (options->outMfa != NULL) {
        // fasta output
        de_verbose("Writing fasta output.\n");
        writeFastaOut(alignmentHash, rowOrder, options);
    }
    if (options->outMaf != NULL) {
        // maf output
        de_verbose("Writing maf output.\n");
        writeMafOut(alignmentHash, rowOrder, options);
    }
    // cleanup
    maf_destroyMfa(mfapi);
    stHash_destruct(alignmentHash);
    stHash_destruct(sequenceHash);
    stList_destruct(rowOrder);
    destroyOptions(options);
    return(EXIT_SUCCESS);
}
