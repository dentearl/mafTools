/*
 * Copyright (C) 2011-2014 by
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

#ifndef _BLOCK_EXTRACTOR_API_H_
#define _BLOCK_EXTRACTOR_API_H_

#include <stdbool.h>
#include <inttypes.h>
#include "common.h"
#include "sharedMaf.h"

bool checkRegion(uint64_t targetStart, uint64_t targetStop, uint64_t lineStart,
                 uint64_t length, uint64_t sourceLength, char strand);
bool searchMatched(mafLine_t *ml, const char *seq, uint64_t start, uint64_t stop);
void printHeader(void);
uint64_t getTargetColumns(bool **targetColumns, uint64_t *n, mafBlock_t *b, const char *seq,
                          uint64_t start, uint64_t stop);
void printTargetColumns(bool *targetColumns, uint64_t n);
int64_t **createOffsets(uint64_t n);
void destroyOffsets(int64_t **offs, uint64_t n);
mafBlock_t *processBlockForSplice(mafBlock_t *b, uint64_t blockNumber, const char *seq,
                                  uint64_t start, uint64_t stop, bool store);
mafBlock_t *spliceBlock(mafBlock_t *mb, uint64_t l, uint64_t r, int64_t **offsetArray);
void checkBlock(mafBlock_t *b, uint64_t blockNumber, const char *seq, uint64_t start,
                uint64_t stop, bool *printedHeader, bool isSoft);
void processBody(mafFileApi_t *mfa, char *seq, uint64_t start, uint64_t stop, bool isSoft);
uint64_t sumBool(bool *array, uint64_t n);
void printOffsetArray(int64_t **offsetArray, uint64_t n);

#endif // _BLOCK_EXTRACTOR_API_H_
