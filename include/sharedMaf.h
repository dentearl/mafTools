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
#ifndef SHAREDMAF_H_
#define SHAREDMAF_H_
#include <stdbool.h>
#include <stdint.h>

typedef struct mafFileApi mafFileApi_t;
typedef struct mafBlock mafBlock_t;
typedef struct mafLine mafLine_t;

// creators, destroyers
mafFileApi_t* maf_newMfa(const char *filename, char const *mode);
mafBlock_t* maf_newMafBlock(void);
mafBlock_t* maf_newMafBlockFromString(const char *s, uint32_t lineNumber);
mafLine_t* maf_newMafLine(void);
mafLine_t* maf_newMafLineFromString(const char *s, uint32_t lineNumber);
mafBlock_t* maf_copyMafBlock(mafBlock_t *orig);
mafLine_t* maf_copyMafLine(mafLine_t *orig);
void maf_destroyMafLineList(mafLine_t *ml);
void maf_destroyMafBlockList(mafBlock_t *mb);
void maf_destroyMfa(mafFileApi_t *mfa);
void maf_mafBlock_destroySequenceMatrix(char **mat, unsigned n);
// read / write
mafBlock_t* maf_readAll(mafFileApi_t *mfa);
mafBlock_t* maf_readBlock(mafFileApi_t *mfa);
mafBlock_t* maf_readBlockHeader(mafFileApi_t *mfa);
mafBlock_t* maf_readBlockBody(mafFileApi_t *mfa);
void maf_writeAll(mafFileApi_t *mfa, mafBlock_t *mb);
void maf_writeBlock(mafFileApi_t *mfa, mafBlock_t *mb);
uint32_t maf_mafFileApi_getLineNumber(mafFileApi_t *mfa);
// getters
char* maf_mafFileApi_getFilename(mafFileApi_t *mfa);
uint32_t maf_mafFileApi_getLineNumber(mafFileApi_t *mfa);
mafLine_t* maf_mafBlock_getHeadLine(mafBlock_t *mb);
mafLine_t* maf_mafBlock_getTailLine(mafBlock_t *mb);
uint32_t maf_mafBlock_getLineNumber(mafBlock_t *mb);
uint32_t maf_mafBlock_getNumberOfLines(mafBlock_t *b);
uint32_t maf_mafBlock_getNumberOfSequences(mafBlock_t *b);
char* maf_mafBlock_getStrandArray(mafBlock_t *mb);
int* maf_mafBlock_getStrandIntArray(mafBlock_t *mb);
uint32_t* maf_mafBlock_getPosCoordStartArray(mafBlock_t *mb);
uint32_t* maf_mafBlock_getPosCoordLeftArray(mafBlock_t *mb);
uint32_t* maf_mafBlock_getStartArray(mafBlock_t *mb);
uint32_t* maf_mafBlock_getSourceLengthArray(mafBlock_t *mb);
uint32_t* maf_mafBlock_getSequenceLengthArray(mafBlock_t *mb);
char** maf_mafBlock_getSpeciesArray(mafBlock_t *mb);
mafBlock_t* maf_mafBlock_getNext(mafBlock_t *mb);
char** maf_mafBlock_getSequenceMatrix(mafBlock_t *mb, unsigned n, unsigned m);
mafLine_t** maf_mafBlock_getMafLineArray_seqOnly(mafBlock_t *mb);
uint32_t maf_mafBlock_getSequenceFieldLength(mafBlock_t *mb);
char* maf_mafLine_getLine(mafLine_t *ml);
uint32_t maf_mafLine_getLineNumber(mafLine_t *ml);
char maf_mafLine_getType(mafLine_t *ml);
char* maf_mafLine_getSpecies(mafLine_t *ml);
uint32_t maf_mafLine_getStart(mafLine_t *ml);
uint32_t maf_mafLine_getLength(mafLine_t *ml);
char maf_mafLine_getStrand(mafLine_t *ml);
uint32_t maf_mafLine_getSourceLength(mafLine_t *ml);
char* maf_mafLine_getSequence(mafLine_t *ml);
mafLine_t* maf_mafLine_getNext(mafLine_t *ml);
// setters
void maf_mafBlock_setHeadLine(mafBlock_t *mb, mafLine_t *ml);
void maf_mafBlock_setTailLine(mafBlock_t *mb, mafLine_t *ml);
void maf_mafBlock_setNumberOfLines(mafBlock_t *mb, uint32_t n);
void maf_mafBlock_incrementNumberOfLines(mafBlock_t *mb);
void maf_mafBlock_decrementNumberOfLines(mafBlock_t *mb);
void maf_mafBlock_setNumberOfSequences(mafBlock_t *mb, uint32_t n);
void maf_mafBlock_incrementNumberOfSequences(mafBlock_t *mb);
void maf_mafBlock_decrementNumberOfSequences(mafBlock_t *mb);
void maf_mafBlock_setLineNumber(mafBlock_t *mb, uint32_t n);
void maf_mafBlock_incrementLineNumber(mafBlock_t *mb);
void maf_mafBlock_decrementLineNumber(mafBlock_t *mb);
void maf_mafBlock_setSequenceFieldLength(mafBlock_t *mb, uint32_t sfl);
void maf_mafLine_setLine(mafLine_t *ml, char *line);
void maf_mafLine_setLineNumber(mafLine_t *ml, uint32_t n);
void maf_mafLine_setType(mafLine_t *ml, char c);
void maf_mafLine_setSpecies(mafLine_t *ml, char *s);
void maf_mafLine_setStrand(mafLine_t *ml, char c);
void maf_mafLine_setStart(mafLine_t *ml, uint32_t n);
void maf_mafLine_setLength(mafLine_t *ml, uint32_t n);
void maf_mafLine_setSourceLength(mafLine_t *ml, uint32_t n);
void maf_mafLine_setSequence(mafLine_t *ml, char *s);
void maf_mafLine_setNext(mafLine_t *ml, mafLine_t *next);
// utilities
char* maf_mafLine_imputeLine(mafLine_t* ml);
unsigned maf_mafBlock_getNumberOfBlocks(mafBlock_t *b);
bool maf_mafBlock_containsSequence(mafBlock_t *m);
uint32_t maf_mafLine_getNumberOfSequences(mafLine_t *m);
uint32_t maf_mafLine_getPositiveCoord(mafLine_t *ml);
uint32_t maf_mafLine_getPositiveLeftCoord(mafLine_t *ml);
unsigned umax(unsigned a, unsigned b);
// print
void maf_mafBlock_print(mafBlock_t *m);
#endif // SHAREDMAF_H_
