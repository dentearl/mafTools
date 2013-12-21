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

#ifndef _PAIR_COVERAGE_API_H_
#define _PAIR_COVERAGE_API_H_

#include <stdbool.h>
#include <inttypes.h>
#include "common.h"
#include "sharedMaf.h"
#include "sonLib.h"
#include "mafPairCoverage.h"

typedef struct mafCoverageCount mafCoverageCount_t;
typedef struct _BinContainer BinContainer;

mafCoverageCount_t* createMafCoverageCount(void);
uint64_t mafCoverageCount_getSourceLength(mafCoverageCount_t *mcct);
uint64_t mafCoverageCount_getObservedLength(mafCoverageCount_t *mcct);
uint64_t mafCoverageCount_getCount(mafCoverageCount_t *mcct);
uint64_t mafCoverageCount_getInRegion(mafCoverageCount_t *mcct);
uint64_t mafCoverageCount_getOutRegion(mafCoverageCount_t *mcct);
void mafCoverageCount_setSourceLength(mafCoverageCount_t *mcct, uint64_t n);
void mafCoverageCount_setCount(mafCoverageCount_t *mcct, uint64_t n);
void mafCoverageCount_setInRegion(mafCoverageCount_t *mcct, uint64_t n);
void mafCoverageCount_setOutRegion(mafCoverageCount_t *mcct, uint64_t n);
int64_t binContainer_getBinStart(BinContainer *bc);
int64_t binContainer_getBinEnd(BinContainer *bc);
int64_t binContainer_getBinLength(BinContainer *bc);
int64_t binContainer_getNumBins(BinContainer *bc);
uint64_t* binContainer_getBins(BinContainer *bc);
uint64_t binContainer_accessBin(BinContainer *bc, int64_t i);
void binContainer_setBinStart(BinContainer *bc, int64_t i);
void binContainer_setBinEnd(BinContainer *bc, int64_t i);
void binContainer_setBinLength(BinContainer *bc, int64_t);
void binContainer_incrementPosition(BinContainer *bc, int64_t i);
void binContainer_incrementBin(BinContainer *bc, int64_t i);
void binContainer_setBinValue(BinContainer *bc, int64_t i, int64_t v);
bool is_wild(const char *s);
bool inInterval(stHash *intervalsHash, char *seq, uint64_t pos);
bool searchMatched(mafLine_t *ml, const char *seq);
bool searchMatched_(const char *target, const char *seq);
void compareLines(mafLine_t *ml1, mafLine_t *ml2, stHash *seq1Hash,
                  stHash *seq2Hash, uint64_t *alignedPositions,
                  stHash *intervalsHash, BinContainer *bc);
void wrapDestroyMafLine(void *p);
void checkBlock(mafBlock_t *b, const char *seq1, const char *seq2,
                stHash *seq1Hash, stHash *seq2Hash, uint64_t *alignedPositions,
                stHash *intervalsHash, BinContainer *bc);
void processBody(mafFileApi_t *mfa, char *seq1, char *seq2, stHash *seq1Hash,
                 stHash *seq2Hash,
                 uint64_t *alignedPositions, stHash *intervalsHash,
                 BinContainer *bc);
void parseBedFile(const char *filepath, stHash *intervalsHash);
void reportResultsBins(char *seq1, char *seq2, BinContainer *bin_container);
BinContainer* binContainer_init(void);
BinContainer* binContainer_construct(int64_t bin_start, int64_t bin_end,
                                     int64_t bin_length);
void binContainer_destruct(BinContainer *bc);

#endif // _PAIR_COVERAGE_API_H_
