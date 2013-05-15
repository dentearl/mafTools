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
#ifndef _PAIR_COVERAGE_H_
#define _PAIR_COVERAGE_H_

#include <stdbool.h>
#include <stdint.h>
#include "common.h"
#include "sharedMaf.h"
#include "sonLib.h"

void version(void);
void usage(void);
void parseOptions(int argc, char **argv, char *filename, char *seq1Name, char *seq2Name, stHash *intervalsHashn);
void reportResults(char *seq1, char *seq2, stHash *seq1Hash, stHash *seq2Hash, uint64_t *alignedPositions);
void reportResultsRegion(char *seq1, char *seq2, stHash *seq1Hash, stHash *seq2Hash, uint64_t *alignedPositions,
                   stHash *intervalsHash);

#endif // _PAIR_COVERAGE_H_
