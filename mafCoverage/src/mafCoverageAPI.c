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

/*
 * The following functions are not currently used.
 */

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

/*
 * These functions are used to gather and sift the sequence names and the lengths of the sequences in mafs.
 */

stHash *getMapOfSequenceNamesToSizesFromMaf(char *mafFileName) {
    stHash *sequenceNamesToSequenceSizes = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, free, (void(*)(void *)) stIntTuple_destruct);

    //Walk through the MAF
    mafFileApi_t *mfa = maf_newMfa(mafFileName, "r");
    mafBlock_t *thisBlock = NULL;
    while ((thisBlock = maf_readBlock(mfa)) != NULL) {
        mafLine_t *ml = maf_mafBlock_getHeadLine(thisBlock);
        while (ml != NULL) {
            if (maf_mafLine_getType(ml) == 's') {
                char *sequenceName = maf_mafLine_getSpecies(ml);
                if (stHash_search(sequenceNamesToSequenceSizes, sequenceName) == NULL) {
                    stHash_insert(sequenceNamesToSequenceSizes, stString_copy(sequenceName), stIntTuple_construct1(maf_mafLine_getSourceLength(ml)));
                } else {
                    assert(stIntTuple_get(stHash_search(sequenceNamesToSequenceSizes, sequenceName), 0) == maf_mafLine_getSourceLength(ml));
                }
            }
            ml = maf_mafLine_getNext(ml);
        }
        maf_destroyMafBlockList(thisBlock);
    }
    maf_destroyMfa(mfa);

    return sequenceNamesToSequenceSizes;
}

static char *copySpeciesName2(const char *sequenceName, bool ignoreSpeciesNames) {
    return ignoreSpeciesNames ? stString_copy(sequenceName) : copySpeciesName(sequenceName);
}

stSet *getSpeciesNames(stList *sequenceNames, bool ignoreSpeciesNames) {
    stSet *speciesNames = stSet_construct3(stHash_stringKey, stHash_stringEqualKey, free);
    for(int64_t i=0; i<stList_length(sequenceNames); i++) {
        char *speciesName = copySpeciesName2(stList_get(sequenceNames, i), ignoreSpeciesNames);
        if (stSet_search(speciesNames, speciesName) == NULL) {
            stSet_insert(speciesNames, speciesName);
        } else {
            free(speciesName);
        }
    }
    return speciesNames;
}

stHash *getMapOfSequenceNamesToSequenceSizesForGivenSpeciesOrChr(stHash *sequenceNamesToSequenceSizes, char *speciesOrChrName, bool ignoreSpeciesNames) {
    stHash *sequenceNamesToSequenceSizeForGivenSpecies = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, NULL, NULL);
    stHashIterator *it = stHash_getIterator(sequenceNamesToSequenceSizes);
    char *sequenceName;
    while ((sequenceName = stHash_getNext(it)) != NULL) {
        char *speciesNameOfSequence = copySpeciesName2(sequenceName, ignoreSpeciesNames);
        if (stString_eq(speciesOrChrName, speciesNameOfSequence) || stString_eq(speciesOrChrName, sequenceName)) {
            stHash_insert(sequenceNamesToSequenceSizeForGivenSpecies, sequenceName, stHash_search(sequenceNamesToSequenceSizes, sequenceName));
        }
        free(speciesNameOfSequence);
    }
    stHash_destructIterator(it);
    return sequenceNamesToSequenceSizeForGivenSpecies;
}

int64_t getTotalLengthOfSequences(stHash *sequenceNamesToSequenceSizes) {
    int64_t length = 0;
    stHashIterator *it = stHash_getIterator(sequenceNamesToSequenceSizes);
    char *sequenceName;
    while ((sequenceName = stHash_getNext(it)) != NULL) {
        length += stIntTuple_get(stHash_search(sequenceNamesToSequenceSizes, sequenceName), 0);
    }
    stHash_destructIterator(it);
    return length;
}

/*
 * The Pairwise Coverage structure.
 */

struct _pairwiseCoverage {
    stHash *sequenceCoverages; //A hash of sequence names to char arrays describing the coverage of each base of the given species.
    const stHash *sequenceNamesToSequenceSizeForGivenSpecies; //A hash of sequence names to their lengths, memory not owned by the object.
};

PairwiseCoverage *pairwiseCoverage_construct(const stHash *sequenceNamesToSequenceSizeForGivenSpecies) {
    PairwiseCoverage *pC = st_malloc(sizeof(PairwiseCoverage));
    pC->sequenceCoverages = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, free, free);
    pC->sequenceNamesToSequenceSizeForGivenSpecies = sequenceNamesToSequenceSizeForGivenSpecies;
    //Now build the sequence coverage arrays.
    stHashIterator *it = stHash_getIterator((stHash *)pC->sequenceNamesToSequenceSizeForGivenSpecies);
    char *sequenceName;
    while ((sequenceName = stHash_getNext(it)) != NULL) {
        stHash_insert(pC->sequenceCoverages, stString_copy(sequenceName),
                st_calloc(stIntTuple_get(stHash_search((stHash *)pC->sequenceNamesToSequenceSizeForGivenSpecies, sequenceName), 0), sizeof(char)));
    }
    stHash_destructIterator(it);
    return pC;
}

void pairwiseCoverage_destruct(PairwiseCoverage *pC) {
    stHash_destruct(pC->sequenceCoverages);
    free(pC);
}

char *pairwiseCoverage_getCoverageArrayForSequence(PairwiseCoverage *pC, char *sequenceName) {
    char *sequenceCoverageArray = stHash_search(pC->sequenceCoverages, sequenceName);
    assert(sequenceCoverageArray != NULL);
    return sequenceCoverageArray;
}

bool pairwiseCoverageArray_increase(char *sequenceCoverageArray, int64_t position) {
    assert(position >= 0);
    if ((int64_t) sequenceCoverageArray[position] < SCHAR_MAX) {
        sequenceCoverageArray[position]++;
        return 0;
    }
    return 1;
}

double *pairwiseCoverage_calculateNCoverages(PairwiseCoverage *pC) {
    stHashIterator *it = stHash_getIterator((stHash *)pC->sequenceNamesToSequenceSizeForGivenSpecies);
    double *nCoverages = st_calloc(SCHAR_MAX + 1, sizeof(double));
    char *sequenceName;
    while ((sequenceName = stHash_getNext(it)) != NULL) {
        int64_t chromosomeLength = stIntTuple_get(stHash_search((stHash *)pC->sequenceNamesToSequenceSizeForGivenSpecies, sequenceName), 0);
        char *chromosomeCoverage = stHash_search(pC->sequenceCoverages, sequenceName);
        for (int64_t i = 0; i < chromosomeLength; i++) {
            for (int64_t j = chromosomeCoverage[i]; j >= 0; j--) {
                nCoverages[j]++;
            }
        }
    }
    stHash_destructIterator(it);
    int64_t genomeLength = getTotalLengthOfSequences((stHash *)pC->sequenceNamesToSequenceSizeForGivenSpecies);
    for (int64_t i = 0; i <= SCHAR_MAX; i++) {
        nCoverages[i] /= genomeLength;
    }
    return nCoverages;
}

double pairwiseCoverage_calculateCoverage(PairwiseCoverage *pC) {
    double *nCoverages = pairwiseCoverage_calculateNCoverages(pC);
    double coverage = nCoverages[1];
    free(nCoverages);
    return coverage;
}

/*
 * The one-species-to-all-other-species structure.
 */

struct _nGenomeCoverage {
    //Top level hash is hash of other species names to pairwise coverage structs.
    char *speciesOrChrName;
    stHash *pairwiseCoverages;
    stHash *sequenceNamesToSequenceSizeForGivenSpeciesOrChr;
    bool ignoreSpeciesNames;
};

NGenomeCoverage *nGenomeCoverage_construct(stHash *sequenceNamesToSequenceSizes, char *speciesOrChrName, bool ignoreSpeciesNames) {
    NGenomeCoverage *nGC = st_malloc(sizeof(NGenomeCoverage));
    nGC->speciesOrChrName = stString_copy(speciesOrChrName);
    //The subset of sequence names and species we care about.
    nGC->sequenceNamesToSequenceSizeForGivenSpeciesOrChr = getMapOfSequenceNamesToSequenceSizesForGivenSpeciesOrChr(sequenceNamesToSequenceSizes, speciesOrChrName, ignoreSpeciesNames);
    //Build the N different pairwise coverage structures.
    nGC->pairwiseCoverages = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, free, (void(*)(void *)) pairwiseCoverage_destruct);
    stList *sequenceNames = stHash_getKeys(sequenceNamesToSequenceSizes);
    stSet *speciesNames = getSpeciesNames(sequenceNames, ignoreSpeciesNames);
    stSetIterator *it = stSet_getIterator(speciesNames);
    char *speciesName2;
    while ((speciesName2 = stSet_getNext(it)) != NULL) {
        stHash_insert(nGC->pairwiseCoverages, stString_copy(speciesName2), pairwiseCoverage_construct(nGC->sequenceNamesToSequenceSizeForGivenSpeciesOrChr));
    }
    nGC->ignoreSpeciesNames = ignoreSpeciesNames;
    //Cleanup
    stSet_destructIterator(it);
    stList_destruct(sequenceNames);
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
                char *lineSpeciesName = copySpeciesName2(maf_mafLine_getSpecies(ml), nGC->ignoreSpeciesNames);
                if (stString_eq(nGC->speciesOrChrName, lineSpeciesName) || stString_eq(nGC->speciesOrChrName, maf_mafLine_getSpecies(ml))) {
                    stList_append(querySpeciesLines, ml);
                }
                stList_append(targetSpeciesLines, ml);
                stList_append(targetSpeciesPairwiseCoverages, stHash_search(nGC->pairwiseCoverages, lineSpeciesName));
                free(lineSpeciesName);
            }
            ml = maf_mafLine_getNext(ml);
        }
        for (int64_t i = 0; i < stList_length(querySpeciesLines); i++) {
            mafLine_t *qML = stList_get(querySpeciesLines, i);
            char *querySequenceFragment = maf_mafLine_getSequence(qML);
            char *querySequenceName = maf_mafLine_getSpecies(qML);
            //Asserts on coordinates
            int64_t querySequenceLength = stIntTuple_get(stHash_search(nGC->sequenceNamesToSequenceSizeForGivenSpeciesOrChr, querySequenceName), 0);
            (void)querySequenceLength;
            assert(querySequenceLength == maf_mafLine_getSourceLength(qML));
            if(maf_mafLine_getLength(qML) > 0) { //This if is because mafs produced by hal2maf sometimes break these assumptions for all gap lines.
                assert(maf_mafLine_getPositiveCoord(qML) >= 0);
                assert(maf_mafLine_getPositiveCoord(qML) <= querySequenceLength);
                if(maf_mafLine_getStrand(qML) == '+') {
                    assert(maf_mafLine_getPositiveCoord(qML) + maf_mafLine_getLength(qML) <= maf_mafLine_getSourceLength(qML));
                }
                else {
                    assert(maf_mafLine_getStrand(qML) == '-');
                    assert((int64_t)(maf_mafLine_getPositiveCoord(qML) - maf_mafLine_getLength(qML)) >= -1);
                }
            }
            assert(strlen(querySequenceFragment) == maf_mafLine_getSequenceFieldLength(qML));
            for (int64_t j = 0; j < stList_length(targetSpeciesLines); j++) {
                mafLine_t *tML = stList_get(targetSpeciesLines, j);
                if(qML != tML) { //To allow self alignments we must ignore identity alignments
                    PairwiseCoverage *pC = stList_get(targetSpeciesPairwiseCoverages, j);
                    assert(pC != NULL);
                    char *coverageArray = pairwiseCoverage_getCoverageArrayForSequence(pC, querySequenceName);
                    assert(coverageArray != NULL);
                    char *targetSequenceFragment = maf_mafLine_getSequence(tML);
                    int64_t position = maf_mafLine_getPositiveCoord(qML);
                    assert(maf_mafLine_getSequenceFieldLength(qML) == maf_mafLine_getSequenceFieldLength(tML));
                    assert(strlen(targetSequenceFragment) == maf_mafLine_getSequenceFieldLength(tML));
                    bool saturated = 1;
                    for (int64_t k = 0; k < maf_mafLine_getSequenceFieldLength(qML); k++) {
                        assert(querySequenceFragment[k] != '\0');
                        assert(querySequenceFragment[k] != ' ');
                        assert(querySequenceFragment[k] != '\n');
                        if (querySequenceFragment[k] != '-') {
                            assert(position >= 0);
                            assert(position < querySequenceLength);
                            if(targetSequenceFragment[k] != '-' && (!requireIdentityForMatch || (toupper(querySequenceFragment[k]) != 'N' && toupper(querySequenceFragment[k]) == toupper(targetSequenceFragment[k])))) {
                                if(pairwiseCoverageArray_increase(coverageArray, position)) {
                                    saturated = 0;
                                }
                            }
                            else {
                                saturated = 0;
                            }
                            position += maf_mafLine_getStrand(qML) == '+' ? 1 : -1;
                        }
                    }
                    if(saturated) {
                        break;
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
    stHash_destruct(nGC->pairwiseCoverages);
    stHash_destruct(nGC->sequenceNamesToSequenceSizeForGivenSpeciesOrChr);
    free(nGC->speciesOrChrName);
    free(nGC);
}

void nGenomeCoverage_reportHeader(FILE *out, bool includeNCoverage) {
    fprintf(out, "referenceSpecies/Chr\tquerySpecies/Chr\tlengthOfReferenceGenome\tcoverage\t");
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
    int64_t speciesGenomeLength = getTotalLengthOfSequences(nGC->sequenceNamesToSequenceSizeForGivenSpeciesOrChr);
    while ((speciesName = stHash_getNext(it)) != NULL) {
        PairwiseCoverage *pC = stHash_search(nGC->pairwiseCoverages, speciesName);
        assert(pC != NULL);
        fprintf(out, "%s\t%s\t%" PRId64 "\t%f", nGC->speciesOrChrName, speciesName, speciesGenomeLength, pairwiseCoverage_calculateCoverage(pC));
        if (includeNCoverage) {
            double *nCoverages = pairwiseCoverage_calculateNCoverages(pC);
            for (int i = 2; i < SCHAR_MAX; i++) {
                fprintf(out, "\t%f", nCoverages[i]);
            }
            free(nCoverages);
        }
        fprintf(out, "\n");
    }
    stHash_destructIterator(it);
}
