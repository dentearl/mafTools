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

#ifndef _MAF_COVERAGE_API_H_
#define _MAF_COVERAGE_API_H_

#include <stdbool.h>
#include <inttypes.h>
#include "common.h"
#include "sharedMaf.h"
#include "sonLib.h"
#include "mafCoverage.h"

bool is_wild(const char *s);
bool searchMatched(mafLine_t *ml, const char *seq);
bool searchMatched_(const char *target, const char *seq);

/*
 * Iterates through the maf and builds a hash of sequence names to coordinates.
 * Lengths are specified by an stIntTuple.
 */
stHash *getSequenceSizes(char *mafFileName);

/*
 * Each sequence name is comprised of two fields separated by a period. The first is the species field, the second is the
 * chromosome field. This function returns the set of distinct species names.
 */
stSet *getSpecies(stHash *sequenceSizes);

/*
 * Gets the subset of the hash for all sequences involving the given species.
 */
stHash *getSequenceSizesForGivenSpecies(stHash *allSequenceSizes, char *speciesName);

/*
 * Returns the combined length of all the sequences in the set.
 */
int64_t getTotalLength(stHash *sequenceSizes);

/*
 * The pairwise coverage object.
 */

typedef struct _pairwiseCoverage PairwiseCoverage;

PairwiseCoverage *pairwiseCoverage_construct(stHash *allSequenceSizes, char *speciesName);

void pairwiseCoverage_destruct(PairwiseCoverage *pC);



/*
 * Returns the coverage of the target genome on query species, that is the proportion of bases in the query aligned to one
 * or more positions in the target.
 */
double pairwiseCoverage_calculateCoverage(PairwiseCoverage *pC);

/*
 * Returns an array of the n-coverages upto but excluding 128, with the index corresponding to n.
 */
double *pairwiseCoverage_calculateNCoverages(PairwiseCoverage *pC);

/*
 * Increases the coverage count of a given sequence position.
 */
char *pairwiseCoverage_getCoverageArrayForSequence(PairwiseCoverage *pC, char *sequenceName);

void pairwiseCoverageArray_increase(char *sequenceCoverageArray, int64_t position);

/*
 * An all-against-a-given-species object.
 */

typedef struct _nGenomeCoverage NGenomeCoverage;

void nGenomeCoverage_destruct(NGenomeCoverage *nGC);

NGenomeCoverage *nGenomeCoverage_construct(stHash *sequenceSizes, char *speciesName, bool nCoverage);

/*
 * Iterate through a maf file and populate the species coverages.
 */
void nGenomeCoverage_populate(NGenomeCoverage *nGC, char *mafFileName, bool requireIdentityForMatch);

/*
 * Reports stats in tab delimited format.
 */
void nGenomeCoverage_reportHeader(FILE *out, bool includeNCoverage);
void nGenomeCoverage_report(NGenomeCoverage *nGC, FILE *out, bool includeNCoverage);

#endif // _MAF_COVERAGE_API_H_
