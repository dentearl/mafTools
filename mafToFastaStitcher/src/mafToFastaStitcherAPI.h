/* 
 * Copyright (C) 2013 by 
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
#ifndef MAFTOFASTASTITCHER_API_H_
#define MAFTOFASTASTITCHER_API_H_
#include <stdint.h>
#include "common.h"
#include "CuTest.h"
#include "sharedMaf.h"
#include "sonLib.h"

typedef struct _options {
    // used to hold all the command line options
    char *maf;
    char *seqs;
    char *outMfa;
    char *outMaf;
    char *reference;
    uint64_t breakpointPenalty;
    uint64_t interstitialSequence;
} options_t;
typedef struct _sequence {
    // used to store fasta sequence elements
    char *seq; // DNA sequence
    uint64_t index; // first empty position in *seq
    uint64_t memLength; // size of the *seq buffer
} mtfseq_t;
typedef struct _row {
    // used to store the ultimate output of the utility,
    // a single element of either a multiple fasta alignment (mfa)
    // or a single row in a multiple alignment format (maf) file.
    char *name;
    char *prevName; // 
    char *sequence;
    bool multipleNames; // initalized false, if prevName is ever != name, then this should be set permanently true
    uint64_t start; 
    uint64_t length;
    uint64_t prevRightPos; // rightmost position in the sequence, 0 based
    char strand; // `+' `-' or `*' when both strands have been observed (multipleNames should be set true)
    char prevStrand; //
    uint64_t sourceLength;
    uint64_t index; // first empty position in *sequence
    uint64_t memLength; //size of the *sequence buffer
} row_t;

options_t* options_construct(void);
void destroyOptions(options_t *o);
mtfseq_t* newMtfseq(uint64_t length);
void resizeMtfseq(mtfseq_t *m);
void resizeRowSequence(row_t *r);
void destroyMtfseq(void *p);
row_t* newRow(uint64_t length);
void destroyRow(void *row);
row_t* mafLineToRow(mafLine_t *ml);
stHash* mafBlockToBlockHash(mafBlock_t *mb, stList *orderList);
stHash* createSequenceHash(char *fastas);
void seq_copyIn(mtfseq_t *mtfss, char *src);
void row_copyIn(row_t *row, char *src);
void addSequencesToHash(stHash *hash, char *filename);
void reportSequenceHash(stHash *hash);
void penalize(stHash *hash, char *name, uint64_t n);
void extendSequence(row_t *r, uint64_t n);
void interstitialInsert(stHash *alignHash, stHash *seqHash, char *name, uint64_t pos, char strand, uint64_t n);
char* extractSubSequence(mtfseq_t *mtfs, char strand, uint64_t pos, uint64_t n);
void addMafLineToRow(row_t *row, mafLine_t *ml);
void addMafBlockToRowHash(stHash *alignHash, stHash *seqHash, stList *order, mafBlock_t *mb, options_t *options);
void prependGaps(row_t *r, uint64_t n);
void buildAlignmentHash(mafFileApi_t *mfapi, stHash *alignmentHash, stHash *sequenceHash, 
                        stList *rowOrder, options_t *options);
void writeFastaOut(stHash *alignmentHash, stList *rowOrder, options_t *options);
void writeMafOut(stHash *alignmentHash, stList *rowOrder, options_t *options);
uint64_t nearestTwo(uint64_t n);
#endif // MAFTOFASTASTITCHER_API_H_
