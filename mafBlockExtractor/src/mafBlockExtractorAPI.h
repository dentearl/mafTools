/* 
 * Copyright (C) 2011-2012 by 
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

bool checkRegion(uint32_t targetStart, uint32_t targetStop, uint32_t lineStart, 
                 uint32_t length, uint32_t sourceLength, char strand);
bool searchMatched(mafLine_t *ml, const char *seq, uint32_t start, uint32_t stop);
void printHeader(void);
uint32_t getTargetColumns(bool **targetColumns, uint32_t *n, mafBlock_t *b, const char *seq, 
                          uint32_t start, uint32_t stop);
void printTargetColumns(bool *targetColumns, uint32_t n);
mafBlock_t *processBlockForTrim(mafBlock_t *b, const char *seq, uint32_t start, uint32_t stop, bool isLeft);
mafBlock_t *trimBlock(mafBlock_t *b, uint32_t n, bool isLeft);
void processBlockForSplit(mafBlock_t *b, const char *seq, uint32_t start, uint32_t stop, 
                          mafBlock_t **left, mafBlock_t **right);
void splitBlock(mafBlock_t *b, uint32_t n, mafBlock_t **mb_left, mafBlock_t **mb_right);
int32_t **createOffsets(uint32_t n);
void destroyOffsets(int32_t **offs, uint32_t n);
mafBlock_t *processBlockForSplice(mafBlock_t *b, uint32_t blockNumber, const char *seq, 
                                  uint32_t start, uint32_t stop, bool store);
mafBlock_t *spliceBlock(mafBlock_t *mb, uint32_t l, uint32_t r, int32_t **offsetArray);
void trimAndReportBlock(mafBlock_t *orig, const char *seq, uint32_t start, uint32_t stop);
void checkBlock(mafBlock_t *b, const char *seq, uint32_t start, uint32_t stop, bool *printedHeader, bool isSoft);
void processBody(mafFileApi_t *mfa, char *seq, uint32_t start, uint32_t stop, bool isSoft);
uint32_t sumBool(bool *array, uint32_t n);

#endif // _BLOCK_EXTRACTOR_API_H_
