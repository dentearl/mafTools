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

typedef struct mafCoverageCount mafCoverageCount_t;

mafCoverageCount_t* createMafCoverageCount(void);
uint64_t mafCoverageCount_getSourceLength(mafCoverageCount_t *mcct);
uint64_t mafCoverageCount_getCount(mafCoverageCount_t *mcct);
void mafCoverageCount_setSourceLength(mafCoverageCount_t *mcct, uint64_t n);
void mafCoverageCount_setCount(mafCoverageCount_t *mcct, uint64_t n);
bool is_wild(const char *s);
bool searchMatched(mafLine_t *ml, const char *seq);
void compareLines(mafLine_t *ml1, mafLine_t *ml2, stHash *seq1Hash, stHash *seq2Hash, uint64_t *alignedPositions);
void wrapDestroyMafLine(void *p);
void checkBlock(mafBlock_t *b, const char *seq1, const char *seq2, 
                stHash *seq1Hash, stHash *seq2Hash, uint64_t *alignedPositions);
void processBody(mafFileApi_t *mfa, char *seq1, char *seq2, stHash *seq1hash, stHash *seq2Hash,
                 uint64_t *alignedPositions);

#endif // _PAIR_COVERAGE_API_H_
