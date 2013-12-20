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
#include <assert.h>
#include <errno.h> // file existence via ENOENT, errno
#include <getopt.h>
#include <inttypes.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <ctype.h>
#include "common.h"
#include "sharedMaf.h"
#include "mafCoverageAPI.h"
#include "bioioC.h" // benLine()
#include "sonLib.h"

bool is_wild(const char *s) {
    // return true if char array ends in *, false otherwise
    if (s[strlen(s) - 1] == '*')
        return true;
    return false;
}
bool searchMatched(mafLine_t *ml, const char *seq) {
    // report false if search did not match, true if it did
    if (maf_mafLine_getType(ml) != 's')
        return false;
    return searchMatched_(maf_mafLine_getSpecies(ml), seq);
}
bool searchMatched_(const char *target, const char *seq) {
    if (is_wild(seq)) {
        // only compare up to the wildcard character
        if (!(strncmp(target, seq, strlen(seq) - 1) == 0)) {
            return false;
        }
    } else {
        if (!(strcmp(target, seq) == 0)) {
            return false;
        }
    }
    return true;
}

stHash *getSequenceSizes(char *mafFileName) {
    stHash *sequenceSizes = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, free, (void(*)(void *)) stIntTuple_destruct);

    //Walk through the MAF
    mafFileApi_t *mfa = maf_newMfa(mafFileName, "r");
    mafBlock_t *thisBlock = NULL;
    while ((thisBlock = maf_readBlock(mfa)) != NULL) {
        mafLine_t *ml = maf_mafBlock_getHeadLine(thisBlock);
        while (ml != NULL) {
            char *sequenceName = copySpeciesName(maf_mafLine_getSpecies(ml));
            if (stHash_search(sequenceSizes, sequenceName) == NULL) {
                stHash_insert(sequenceSizes, sequenceName, stIntTuple_construct1(maf_mafLine_getSourceLength(ml)));
            } else {
                assert(stIntTuple_get(stHash_search(sequenceSizes, sequenceName), 0) == maf_mafLine_getSourceLength(ml));
                free(sequenceName);
            }
            ml = maf_mafLine_getNext(ml);
        }
        maf_destroyMafBlockList(thisBlock);
    }
    maf_destroyMfa(mfa);

    return sequenceSizes;
}

stSet *getSpecies(stHash *sequenceSizes) {
    stSet *speciesNames = stSet_construct3(stHash_stringKey, stHash_stringEqualKey, free);
    stHashIterator *it = stHash_getIterator(sequenceSizes);
    char *sequenceName;
    while ((sequenceName = stHash_getNext(it)) != NULL) {
        char *speciesName = copySpeciesName(sequenceName);
        if (stSet_search(speciesNames, speciesName) == NULL) {
            stSet_insert(speciesNames, speciesName);
        } else {
            free(speciesName);
        }
    }
    stHash_destructIterator(it);
    return speciesNames;
}

stHash *getSequenceSizesForGivenSpecies(stHash *allSequenceSizes, char *speciesName) {
    stHash *sequenceSizes = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, NULL, NULL);
    stHashIterator *it = stHash_getIterator(allSequenceSizes);
    char *sequenceName;
    while ((sequenceName = stHash_getNext(it)) != NULL) {
        if (stString_eq(speciesName, copySpeciesName(sequenceName))) {
            stHash_insert(sequenceSizes, sequenceName, stHash_search(allSequenceSizes, sequenceName));
        }
    }
    stHash_destructIterator(it);
    return sequenceSizes;
}

int64_t getTotalLength(stHash *sequenceSizes) {
    int64_t length = 0;
    stHashIterator *it = stHash_getIterator(sequenceSizes);
    char *sequenceName;
    while ((sequenceName = stHash_getNext(it)) != NULL) {
        length += stIntTuple_get(stHash_search(sequenceSizes, sequenceName), 0);
    }
    stHash_destructIterator(it);
    return length;
}

/*
 * The Pairwise Coverage structure.
 */

struct _pairwiseCoverage {
    stHash *sequenceCoverages; //A hash of sequence names to char arrays describing the coverage of each base of the given species.
    char *speciesName; //The name of species.
    stHash *sequenceSizes; //A hash of sequence names to their lengths.
};

PairwiseCoverage *pairwiseCoverage_construct(stHash *allSequenceSizes, char *speciesName) {
    PairwiseCoverage *pC = st_malloc(sizeof(PairwiseCoverage));
    pC->sequenceCoverages = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, free, free);
    pC->speciesName = stString_copy(speciesName);
    pC->sequenceSizes = getSequenceSizesForGivenSpecies(allSequenceSizes, speciesName);
    stHashIterator *it = stHash_getIterator(pC->sequenceSizes);
    char *sequenceName;
    while ((sequenceName = stHash_getNext(it)) != NULL) {
        stHash_insert(pC->sequenceCoverages, sequenceName,
                st_calloc(stIntTuple_get(stHash_search(pC->sequenceSizes, sequenceName), 0), sizeof(char)));
    }
    stHash_destructIterator(it);
    return pC;
}

void pairwiseCoverage_destruct(PairwiseCoverage *pC) {
    stHash_destruct(pC->sequenceCoverages);
    free(pC->speciesName);
    stHash_destruct(pC->sequenceSizes);
    free(pC);
}

char *pairwiseCoverage_getCoverageArrayForSequence(PairwiseCoverage *pC, char *sequenceName) {
    char *sequenceCoverageArray = stHash_search(pC->sequenceCoverages, sequenceName);
    assert(sequenceCoverageArray != NULL);
    return sequenceCoverageArray;
}

void pairwiseCoverageArray_increase(char *sequenceCoverageArray, int64_t position) {
    assert(position >= 0);
    if ((int64_t) sequenceCoverageArray[position] < SCHAR_MAX) {
        sequenceCoverageArray[position]++;
    }
}

double *pairwiseCoverage_calculateNCoverages(PairwiseCoverage *pC) {
    int64_t length = getTotalLength(pC->sequenceSizes);
    stHashIterator *it = stHash_getIterator(pC->sequenceSizes);
    double *nCoverages = st_calloc(SCHAR_MAX + 1, sizeof(double));
    char *sequenceName;
    while ((sequenceName = stHash_getNext(it)) != NULL) {
        int64_t chromosomeLength = stIntTuple_get(stHash_search(pC->sequenceSizes, sequenceName), 0);
        char *chromosomeCoverage = stHash_search(pC->sequenceCoverages, sequenceName);
        for (int64_t i = 0; i < chromosomeLength; i++) {
            for (int64_t j = chromosomeCoverage[i]; j > 0; j--) {
                nCoverages[j]++;
            }
        }
    }
    for (int64_t i = 0; i <= SCHAR_MAX; i++) {
        nCoverages[i] /= length;
    }
    return nCoverages;
}

double pairwiseCoverage_calculateCoverage(PairwiseCoverage *pC) {
    double *nCoverages = pairwiseCoverage_calculateNCoverages(pC);
    double coverage = nCoverages[1];
    free(nCoverages);
    return coverage;
}

struct _nGenomeCoverage {
    //Top level hash is hash of other species names to pairwise coverage structs.
    char *speciesName;
    int64_t speciesLength;
    stHash *pairwiseCoverages;
};

NGenomeCoverage *nGenomeCoverage_construct(stHash *sequenceSizes, char *speciesName, bool nCoverage) {
    NGenomeCoverage *nGC = st_malloc(sizeof(NGenomeCoverage));
    nGC->speciesName = stString_copy(speciesName);
    stHash *sequencesForSpecies = getSequenceSizesForGivenSpecies(sequenceSizes, speciesName);
    nGC->speciesLength = getTotalLength(sequencesForSpecies);
    stHash_destruct(sequencesForSpecies);
    nGC->pairwiseCoverages = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, free, (void(*)(void *)) pairwiseCoverage_destruct);
    //Build the N different pairwise coverage structures.
    stSet *speciesNames = getSpecies(sequenceSizes);
    stSetIterator *it = stSet_getIterator(speciesNames);
    char *speciesName2;
    while ((speciesName2 = stSet_getNext(it)) != NULL) {
        stHash_insert(nGC->pairwiseCoverages, speciesName2, pairwiseCoverage_construct(sequenceSizes, speciesName));
    }
    stSet_destructIterator(it);
    //Cleanup
    stSet_destruct(speciesNames);
    return nGC;
}

void nGenomeCoverage_populate(NGenomeCoverage *nGC, char *mafFileName, bool requireIdentityForMatch) {
    mafFileApi_t *mfa = maf_newMfa(mafFileName, "r");
    mafBlock_t *thisBlock = NULL;
    while ((thisBlock = maf_readBlock(mfa)) != NULL) {
        mafLine_t *ml = maf_mafBlock_getHeadLine(thisBlock);
        //Get any lines for out target species
        stList *querySpeciesLines = stList_construct();
        stList *targetSpeciesLines = stList_construct();
        stList *targetSpeciesPairwiseCoverages = stList_construct();
        while (ml != NULL) {
            if (maf_mafLine_getType(ml) == 's') {
                char *lineSpeciesName = copySpeciesName(maf_mafLine_getSpecies(ml));
                if (stString_eq(nGC->speciesName, lineSpeciesName)) {
                    stList_append(querySpeciesLines, ml);
                } else {
                    stList_append(targetSpeciesLines, ml);
                    stList_append(targetSpeciesPairwiseCoverages, stHash_search(nGC->pairwiseCoverages, lineSpeciesName));
                }
                free(lineSpeciesName);
            }
            ml = maf_mafLine_getNext(ml);
        }
        for (int64_t i = 0; i < stList_length(querySpeciesLines); i++) {
            mafLine_t *qML = stList_get(querySpeciesLines, i);
            char *querySequence = maf_mafLine_getSequence(qML);
            char *querySequenceName = maf_mafLine_getSpecies(qML);
            for (int64_t j = 0; j < stList_length(targetSpeciesLines); j++) {
                PairwiseCoverage *pC = stList_get(targetSpeciesPairwiseCoverages, j);
                assert(pC != NULL);
                char *coverageArray = pairwiseCoverage_getCoverageArrayForSequence(pC, querySequenceName);
                mafLine_t *tML = stList_get(targetSpeciesLines, j);
                char *targetSequence = maf_mafLine_getSequence(tML);
                int64_t position = maf_mafLine_getPositiveCoord(qML);
                assert(maf_mafLine_getLength(qML) == maf_mafLine_getLength(tML));
                assert(strlen(querySequence) == strlen(targetSequence));
                assert(strlen(querySequence) == maf_mafLine_getLength(qML));
                for (int64_t k = 0; k < maf_mafLine_getLength(qML); k++) {
                    if (querySequence[k] != '-' && targetSequence[k] != '-') {
                        if(!requireIdentityForMatch || (toupper(querySequence[k]) != 'N' && toupper(querySequence[k]) == toupper(targetSequence[k]))) {
                            pairwiseCoverageArray_increase(coverageArray, position);
                            position += maf_mafLine_getStrand(qML) ? 1 : -1;
                        }
                    }
                }
            }
        }
        //Cleanup
        maf_destroyMafBlockList(thisBlock);
        stList_destruct(querySpeciesLines);
        stList_destruct(targetSpeciesLines);
        stList_destruct(targetSpeciesPairwiseCoverages);
    }
    maf_destroyMfa(mfa);
}

void nGenomeCoverage_destruct(NGenomeCoverage *nGC) {
    free(nGC->speciesName);
    stHash_destruct(nGC->pairwiseCoverages);
    free(nGC);
}

void nGenomeCoverage_reportHeader(FILE *out, bool includeNCoverage) {
    fprintf(out, "querySpecies\ttargetSpecies\tlength\tcoverage\t");
    if (includeNCoverage) {
        for (int i = 2; i < SCHAR_MAX; i++) {
            fprintf(out, "\t%i-coverage", i);
        }
    }
    fprintf(out, "\n");
}

void nGenomeCoverage_report(NGenomeCoverage *nGC, FILE *out, bool includeNCoverage) {
    stHashIterator *it = stHash_getIterator(nGC->pairwiseCoverages);
    char *speciesName;
    while ((speciesName = stHash_getNext(it)) != NULL) {
        PairwiseCoverage *pC = stHash_search(nGC->pairwiseCoverages, speciesName);
        fprintf(out, "%s\t%s\t%" PRId64 "\t%f", nGC->speciesName, speciesName, getTotalLength(pC->sequenceSizes), pairwiseCoverage_calculateCoverage(pC));
        if (includeNCoverage) {
            double *nCoverages = pairwiseCoverage_calculateNCoverages(pC);
            for (int i = 2; i < SCHAR_MAX; i++) {
                fprintf(out, "\t%f", nCoverages[i]);
            }
        }
        fprintf(out, "\n");
    }
    stHash_destructIterator(it);
}
