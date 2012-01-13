/* 
 * Copyright (C) 2009-2012 by 
 * Dent Earl (dearl@soe.ucsc.edu, dentearl@gmail.com)
 * Benedict Paten (benedict@soe.ucsc.edu, benedictpaten@gmail.com)
 * Mark Diekhans (markd@soe.ucsc.edu)
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


#ifndef _COMPARATOR_API_H_
#define _COMPARATOR_API_H_

#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <getopt.h>

#include "avl.h"
#include "bioioC.h"
//#include "cactus.h"
#include "commonC.h"
#include "hashTableC.h"
#include "hashTableC_itr.h"
#include "sonLib.h"

#include "disjointset.h"

typedef struct _solo {
    char *name;
    int32_t pos;
} ASolo;

typedef struct _pair {
    char *seq1;
    char *seq2;
    int32_t pos1;
    int32_t pos2;
    int32_t origPos1;
    int32_t origPos2;
} APair;

typedef struct _resultPair {
    APair aPair;
    int32_t inAll;
    int32_t inBoth;
    int32_t inA;
    int32_t inB;
    int32_t inNeither;
    int32_t total;
    int32_t totalBoth;
    int32_t totalA;
    int32_t totalB;
    int32_t totalNeither;
} ResultPair;

typedef struct _trio {
    char *seq1;
    char *seq2;
    char *seq3;
    int32_t pos1;
    int32_t pos2;
    int32_t pos3;
    int32_t top; // tree topology
    int32_t topMat[10];
} ATrio;

typedef struct _trioNames {
        char *speciesA;
        char *speciesB;
        char *speciesC;
} TrioNames;

typedef struct _trioDecoder {
    char **nodeLabelArray;
    char **leafLabelArray;
    struct hashtable *treeLabelHash;
    int32_t **lcaMatrix;
    int32_t nodeNum;
    int32_t leafNum;
} TrioDecoder;

void populateNames(const char *mAFFile, stSortedSet *htp);
void intersectHashes(struct hashtable *h1, struct hashtable *h2, struct hashtable *h3);
void printNameHash(struct hashtable *h);
struct avl_table *compareMAFs_AB(const char *mAFFileA, const char *mAFFileB, int32_t numberOfSamples, 
                                 stSortedSet *legitimateSequences, stHash *intervalsHash, 
                                 int32_t verbose, int32_t near);
struct avl_table *compareMAFs_AB_Trio(const char *mAFFileA, const char *mAFFileB, int32_t numberOfSamples, 
                                      struct hashtable *ht, struct List *speciesList);
void reportResults(struct avl_table *results_AB, const char *mAFFileA, const char *mAFFileB, 
                   FILE *fileHandle, int32_t near, stSortedSet *legitimateSequences, const char *bedFiles);
void reportResultsTrio(struct avl_table *results_AB, const char *mAFFileA, const char *mAFFileB, FILE *fileHandle);
void aPair_destruct(APair *pair, void *extraArgument);
void aTrio_destruct(ATrio *trio, void *extraArgument);

TrioDecoder *trioDecoder_construct(char *treestring);
int32_t calcTrioState(TrioDecoder *decoder, int32_t i, int32_t j, int32_t k);
void writeXMLHeader( FILE *fileHandle );

#endif /* _COMPARATOR_API_H_ */
