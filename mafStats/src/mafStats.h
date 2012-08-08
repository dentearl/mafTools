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
#ifndef _MAFSTATS_H_
#define _MAFSTATS_H_

#include <stdio.h>
#include <stdlib.h>

typedef struct stats {
    char *filename;
    uint64_t numLines;
    uint64_t numHeaderLines;
    uint64_t numSeqLines;
    uint64_t numBlocks;
    uint64_t numELines;
    uint64_t numILines;
    uint64_t numQLines;
    uint64_t numCommentLines;
    uint64_t numGapCharacters;
    uint64_t numSeqCharacters;
    uint64_t numColumns;
    uint64_t sumSeqField;
    uint64_t maxSeqField;
    uint64_t sumNumSpeciesInBlock;
    uint64_t maxNumSpeciesInBlock;
    uint64_t sumBlockArea;
    uint64_t maxBlockArea;
    stHash *seqHash; // keyed with names, valued with uint64_t count of bases present
} stats_t;
typedef struct seq {
    char *name;
    uint64_t count;
} seq_t;

void usage(void);
void parseArgs(int argc, char **argv, char **filename);
stats_t* stats_create(char *filename);
void stats_destroy(stats_t *stats);
void countCharacters(char *seq, stats_t *stats);
void processBlock(mafBlock_t *mb, stats_t *stats);
void recordStats(mafFileApi_t *mfa, stats_t *stats);
void readFilesize(struct stat *fileStat, char **filesizeString);
int cmp_seq(const void *a, const void *b);
void reportHash(stHash *hash);
void reportStats(stats_t *stats);

#endif // _MAFSTATS_H_
