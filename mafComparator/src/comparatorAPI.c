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

#include "sonLib.h"
#include "comparatorAPI.h"
#include "comparatorRandom.h"

void aPair_fillOut(APair *aPair, char *seq1, char *seq2, uint32_t pos1, uint32_t pos2) {
    int i = strcmp(seq1, seq2);
    if (i > 0 || (i == 0 && pos1 > pos2)) { 
        // we construct the sequence so that seq1 < seq2 || seq1 == seq2 and pos1 < pos2
        aPair_fillOut(aPair, seq2, seq1, pos2, pos1);
    } else {
        aPair->seq1 = seq1;
        aPair->seq2 = seq2;
        aPair->pos1 = pos1;
        aPair->pos2 = pos2;
    }
}
APair* aPair_construct(const char *seq1, const char *seq2, uint32_t pos1, uint32_t pos2) {
    APair *aPair = st_malloc(sizeof(APair));
    aPair_fillOut(aPair, (char *) seq1, (char *) seq2, pos1, pos2);
    aPair->seq1 = stString_copy(aPair->seq1);
    aPair->seq2 = stString_copy(aPair->seq2);
    return aPair;
}
ResultPair *resultPair_construct(const char *seq1, const char *seq2) {
    APair *aPair = aPair_construct(seq1, seq2, 0, 0);
    ResultPair *resultPair = st_calloc(1, sizeof(ResultPair));
    resultPair->aPair = *aPair;
    free(aPair);
    return resultPair;
}
APair *aPair_copyConstruct(APair *pair) {
   return aPair_construct(pair->seq1, pair->seq2, pair->pos1, pair->pos2);
}
void aPair_destruct(APair *pair, void *extraArgument) {
    assert(extraArgument == NULL); //because destruction takes place in avl_tree
    free(pair->seq1);
    free(pair->seq2);
    free(pair);
}
void resultPair_destruct(ResultPair *pair, void *extraArgument) {
    assert(extraArgument == NULL); //because destruction takes place in avl_tree
    free(pair->aPair.seq1);
    free(pair->aPair.seq2);
    free(pair);
}
int32_t aPair_cmpFunction_seqsOnly(APair *aPair1, APair *aPair2, void *a) {
    /*
     * Compares only the sequence name arguments of the a pairs.
     */
    assert(a == NULL);
    int32_t i = strcmp(aPair1->seq1, aPair2->seq1);
    if (i == 0) {
        i = strcmp(aPair1->seq2, aPair2->seq2);
    }
    return i;
}
int32_t aPair_cmpFunction(APair *aPair1, APair *aPair2, void *a) {
    /*
     * Compares first the sequence then the position arguments of the a pairs.
     */
    int32_t i = aPair_cmpFunction_seqsOnly(aPair1, aPair2, a);
    if (i == 0) {
        i = aPair1->pos1 - aPair2->pos1;
        if (i == 0) {
            i = aPair1->pos2 - aPair2->pos2;
        }
    }
    return i;
}
void getPairsP(void(*passPairFn)(APair *pair, stHash *intervalsHash, void *extraArgument1, void *extraArgument2,
                                 void *extraArgument3, int32_t, int32_t), 
               stHash *intervalsHash, void *extraArgument1, void *extraArgument2,
               void *extraArgument3, int32_t isVerboseFailures, int32_t near, 
               int *bytesRead, int *nBytes, char **cA, FILE *fileHandle) {
    int32_t length, start, seqLength, i, j, k, pos1, pos2, inc1, inc2, origPos1, origPos2;
    struct List *ranges = constructEmptyList(0, free);
    char *seqName;
    char *sequence;
    char strand;
    static APair aPair;
    // process block function iterates through successive lines while we have not reached
    // the end of the file and the newline starts with 's'
    // for each line grep the sequence, start position and the length.
    // all the lengths must be equal.
    *bytesRead = benLine(cA, nBytes, fileHandle);
    length = INT32_MAX;
    while (*bytesRead > 0 && (*cA)[0] == 's') {
        seqName = st_malloc(sizeof(char) * (1 + (*bytesRead)));
        sequence = st_malloc(sizeof(char) * (1 + (*bytesRead)));
        //uglyf("Got the line :##%s#%i %i\n", *cA, (int)strlen(*cA), *bytesRead);
        //uglyf("I read %i \n", sscanf(*cA, "s %s %i %i %c %i %s", seqName, &start, &i /*ignore the length field*/, &strand, &seqLength, sequence));
        //uglyf("%s,  %i %i %c %i %s\n", seqName, start, i /*ignore the length field*/, strand, seqLength, sequence);
        j = sscanf(*cA, "s %s %i %i %c %i %s", seqName, &start, &i, /*ignore the length field*/
                   &strand, &seqLength, sequence);
        assert(j == 6 || (j == 5 && seqLength == 0));
        if (j == 5) {
            free(sequence);
            sequence = stString_print("");
        }
        length = strlen(sequence);
        assert(strand == '+' || strand == '-');
        listAppend(ranges, seqName);
        if (strand == '+') {
            listAppend(ranges, constructInt(start));
            listAppend(ranges, constructInt(1));
            listAppend(ranges, constructInt(start));
        } else {
            listAppend(ranges, constructInt(seqLength - 1 - start));
            listAppend(ranges, constructInt(-1));
            listAppend(ranges, constructInt(start));
        }
        listAppend(ranges, sequence);
        *bytesRead = benLine(cA, nBytes, fileHandle);
    }
    //Now call the pair function for every pair of aligned bases.
    for (i = 0; i < ranges->length; i += 5) {
        char *seq1 = ranges->list[i];
        inc1 = *((int32_t *) ranges->list[i + 2]);
        const char *sequence1 = ranges->list[i + 4];
        for (j = i + 5; j < ranges->length; j += 5) {
            char *seq2 = ranges->list[j];
            pos2 = *((int32_t *) ranges->list[j + 1]);
            inc2 = *((int32_t *) ranges->list[j + 2]);
            origPos2 = *((int32_t *) ranges->list[j + 3]);
            const char *sequence2 = ranges->list[j + 4];
            pos1 = *((int32_t *) ranges->list[i + 1]);
            origPos1 = *((int32_t *) ranges->list[i + 3]);
            // fprintf(stderr, "%s %d %d %d\n", seq1, pos1, inc1, origPos1);
            // fprintf(stderr, "%s %d %d %d\n", seq2, pos2, inc2, origPos2);
            assert((int32_t)strlen(sequence1) == length);
            assert((int32_t)strlen(sequence2) == length);
            for (k = 0; k < length; k++) {
                if (sequence1[k] != '-') {
                    if (sequence2[k] != '-') {
                        aPair_fillOut(&aPair, seq1, seq2, pos1, pos2);
                        passPairFn(&aPair, intervalsHash, extraArgument1, extraArgument2, 
                                   extraArgument3, isVerboseFailures, near);
                        pos2 += inc2;
                        origPos2 ++;
                    }
                    pos1 += inc1;
                    origPos1 ++;
                } else {
                    if (sequence2[k] != '-') {
                        pos2 += inc2;
                        origPos2 ++;
                    }
                }
            }
        }
    }
    //cleanup the list
    destructList(ranges);
}
struct stringSortIdx {
    char *name;
    int32_t index;
};
int stringSortIdx_cmp(const void *a, const void *b) {
    struct stringSortIdx *ia = (struct stringSortIdx *) a;
    struct stringSortIdx *ib = (struct stringSortIdx *) b;
    return strcmp(ia->name, ib->name);
}
void getPairs(const char *mafFile1, 
              void(*passPairFn)(APair *pair, stHash *intervalsHash, void *extraArgument1, void *extraArgument2, 
                                void *extraArgument3, int32_t isVerboseFailures, int32_t near), 
              stHash *intervalsHash,
              void *extraArgument1, void *extraArgument2, void *extraArgument3, 
              int32_t isVerboseFailures, int32_t near) {
    /*
     * Iterates through all pairs of bases in a set of MAFs, calling the given function for each one.
     */
    FILE *fileHandle = fopen(mafFile1, "r");
    int bytesRead;
    int nBytes = 100;
    char *cA;
    cA = st_malloc(nBytes + 1);
    bytesRead = benLine(&cA, &nBytes, fileHandle);
    // read through lines until reaching a line starting with an 'a':
    // then call process block function passing it the next line.
    while (bytesRead != -1) {
        if (bytesRead > 0 && cA[0] == 'a') {
            getPairsP(passPairFn, intervalsHash, extraArgument1, extraArgument2, 
                      extraArgument3, isVerboseFailures, near,
                      &bytesRead, &nBytes, &cA, fileHandle);
        } else { //deal with empty, white space lines.
            bytesRead = benLine(&cA, &nBytes, fileHandle);
        }
    }
    //Cleanup
    free(cA);
    fclose(fileHandle);
}
int32_t encodeTopoMatIdx(int32_t top1, int32_t top2) {
    int32_t tmp;
    if (top1 > top2) {
        tmp = top1;
        top1 = top2;
        top2 = tmp;
    }
    return 4 * top1 + top2 - ((top1 + 1) * top1 / 2);
}
void pairIndicesToArrayIndex(uint64_t r, uint64_t c, uint64_t n, uint64_t *i) {
    // r = row, c = column, n = number of sequences
    /* A proof by picture:
           c
           0 1 2 3 4
       r 0 - 0 1 2 3
         1   - 4 5 6
         2     - 7 8
         3       - 9
         4         -

       so (3, 4) can expand row by row:
       = (n - 1) + (n - 2) + (n - 3) + c - r - 1
       = 3 * n - (1 + 2 + 3) + c - r - 1
       = r * n - (sum_{i=1}^{r} i) + c - r - 1
       = r * n - (r (r + 1) / 2) + c - r - 1 
       = 3 * 5 - (3 * 4 / 2) + 4 - 3 - 1
       = 15 - 6 + 4 - 3 - 1
       = 9
       so (1, 3) yeilds 5 by 4 + 2 - 1, the 4 comes from 5 - sum_{i=1}^{r} i = r(r+1)/2
       and the 2 comes from c - r.
     */
    assert(n > 0);
    *i = (r * n - (r * (r + 1) / 2) + c - r - 1);
}
void arrayIndexToPairIndices(uint64_t i, uint64_t n, uint64_t *r, uint64_t *c) { 
    // i is the index, n is the length of one side of the square, 
    // r and c are pointers to record the row and colunm of the index.
    assert(n > 0);
    if (n == 1) {
        *r = 0;
        *c = 1;
        return;
    }
    assert(i < ((n * (n - 1)) << 1));
    uint64_t x = (n * (n - 1)) / 2 - 1 - i;
    uint64_t k = floor((sqrt(8 * x + 1) - 1) / 2);
    *r = n - 2 - k;
    *c = i - ((*r) * n - ((*r) * ((*r) + 1)) / 2 - (*r) -1);
    assert((*r) < n);
    assert((*c) < n);
}
bool* getLegitRows(char **names, uint32_t numSeqs, stHash *legitSequences) {
    bool *legitRows = (bool *) st_malloc(sizeof(bool) * numSeqs);
    for (uint32_t i = 0; i < numSeqs; ++i) {
        if (stHash_search(legitSequences, names[i]) != NULL) {
            legitRows[i] = true;
        } else {
            legitRows[i] = false;
        }
    }
    return legitRows;
}
uint64_t countPairsInColumn(char **mat, uint32_t c, uint32_t numSeqs, 
                            bool *legitRows, uint64_t *chooseTwoArray) {
    uint32_t possiblePartners = 0;
    for (uint32_t r = 0; r < numSeqs; ++r) {
        if (!legitRows[r]) {
            continue;
        }
        if (mat[r][c] != '-') {
            ++possiblePartners;
        }
    }
    if (possiblePartners < 101) {
        return chooseTwoArray[possiblePartners];
    } else {
        return chooseTwo(possiblePartners);
    }
}
uint64_t walkBlockCountingPairs(mafBlock_t *mb, stHash *legitSequences, uint64_t *chooseTwoArray) {
    // size is the MAXIMUM INDEX of the chooseTwo array
    uint64_t count = 0;
    uint32_t numSeqs = maf_mafBlock_getNumberOfSequences(mb);
    if (numSeqs < 1) {
        return 0;
    }
    uint32_t seqFieldLength = maf_mafBlock_longestSequenceField(mb);
    char **names = maf_mafBlock_getSpeciesArray(mb);
    char **mat = maf_mafBlock_getSequenceMatrix(mb, numSeqs, seqFieldLength);
    bool *legitRows = getLegitRows(names, numSeqs, legitSequences);
    for (uint32_t c = 0; c < seqFieldLength; ++c) {
        count += countPairsInColumn(mat, c, numSeqs, legitRows, chooseTwoArray);
    }
    // clean up
    for (uint32_t i = 0; i < numSeqs; ++i) {
        free(names[i]);
    }
    free(names);
    maf_mafBlock_destroySequenceMatrix(mat, numSeqs);
    free(legitRows);
    return count;
}
uint64_t chooseTwo(uint32_t n) {
    if ((n == 0) || (n == 1)) {
        return 0;
    } else {
        return ((n * (n - 1)) >> 1);
    }
}
uint64_t* buildChooseTwoArray(void) {
    // pre-calculate a bunch of smaller sizes
    uint64_t *cta = (uint64_t *) st_malloc(sizeof(int64_t) * 101);
    for (uint32_t i = 0; i < 101; ++i) {
        cta[i] = chooseTwo(i);
    }
    return cta;
}
uint64_t countPairsInMaf(const char *filename, stHash *legitSequences) {
    mafFileApi_t *mfa = maf_newMfa(filename, "r");
    mafBlock_t *mb = NULL;
    uint64_t counter = 0;
    uint64_t *chooseTwoArray = buildChooseTwoArray();
    while ((mb = maf_readBlock(mfa)) != NULL) {
        counter += walkBlockCountingPairs(mb, legitSequences, chooseTwoArray);
        maf_destroyMafBlockList(mb);
    }
    // clean up
    free(chooseTwoArray);
    maf_destroyMfa(mfa);
    return counter;
}
static int int32EqualKey(const void *key1, const void *key2) {
    return *((int32_t *) key1) == *((int32_t*) key2);
}
static uint32_t int32Key(const void *k) {
    return *((int32_t*)k);
}
int32_t* copyInt32(int32_t *i) {
    int32_t *j = st_malloc(sizeof(*j));
    *j = *i;
    return j;
}
void printmlarray(mafLine_t **mlArray, uint32_t n) {
    for (uint32_t i = 0; i < n; ++i) {
        printf("%" PRIu32 ":%s, ", i, maf_mafLine_getSpecies(mlArray[i]));
    }
    printf("\n");
}
uint32_t countLegitGaplessPositions(char **mat, uint32_t c, uint32_t numRows, bool *legitRows) {
    uint32_t a = 0;
    for (uint32_t r = 0; r < numRows; ++r) {
        if (!legitRows[r]) {
            continue;
        }
        if (mat[r][c] != '-') {
            ++a;
        }
    }
    return a;
}
void samplePairsFromColumn(char **mat, uint32_t c, bool *legitRows, double acceptProbability, 
                           struct avl_table *pairs, 
                           uint32_t numRows, uint32_t numLegit, uint64_t *chooseTwoArray, 
                           mafLine_t **mlArray, uint32_t *positions) {
    // mat is the matrix of characters representing the alignment, c is the current column,
    // legtRows is a boolean array with true for a legit sequence, acceptProbability is the per
    // base accept probability, pairs is where we store pairs, numRows is the number of total sequences
    // in the matrix, numLegit is the number of legit sequences -- also the length
    // of the mlArray and positions arrays.
    mafLine_t **columnMlArray = NULL;
    uint32_t *columnPositions = NULL;
    uint32_t numLegitGaplessPositions = countLegitGaplessPositions(mat, c, numRows, legitRows);
    // printmlarray(mlArray, numRows);
    columnMlArray = cullMlArrayByColumn(mat, c, numRows, legitRows, mlArray, numLegitGaplessPositions);
    // printmlarray(mlArray, sumBoolArray(legitRows, numRows));
    columnPositions = cullPositionsByColumn(mat, c, numRows, legitRows, positions, numLegitGaplessPositions);
    uint64_t numPairs = countPairsInColumn(mat, c, numRows, legitRows, chooseTwoArray);
    if (numPairs < 5) {
        // keep in mind the sequence of v = k choose 2 is 
        // k: 2 3 4  5 ...
        // v: 1 3 6 10 ...
        samplePairsFromColumnBruteForce(acceptProbability, pairs,
                                        numLegitGaplessPositions, chooseTwoArray, columnMlArray, 
                                        columnPositions, numPairs);
    } else {
        samplePairsFromColumnAnalytic(acceptProbability, pairs,
                                      numLegitGaplessPositions, chooseTwoArray, columnMlArray, 
                                      columnPositions, numPairs);
    }
    free(columnMlArray);
    free(columnPositions);
}
void samplePairsFromColumnBruteForce(double acceptProbability, struct avl_table *pairs, 
                                     uint32_t numSeqs, uint64_t *chooseTwoArray,
                                     mafLine_t **mlArray, uint32_t *positions, uint64_t numPairs) {
    uint64_t p1, p2;
    for (uint64_t i = 0; i < numPairs; ++i) {
        if (st_random() <= acceptProbability) {
            arrayIndexToPairIndices(i, numSeqs, &p1, &p2);
            APair *aPair = aPair_construct(maf_mafLine_getSpecies(mlArray[p1]),
                                           maf_mafLine_getSpecies(mlArray[p2]),
                                           positions[p1], positions[p2]);
            avl_insert(pairs, aPair);
        }
    }
}
void samplePairsFromColumnAnalytic(double acceptProbability, struct avl_table *pairs, 
                                   uint32_t numSeqs, uint64_t *chooseTwoArray,
                                   mafLine_t **mlArray, uint32_t *positions, uint64_t numPairs) {
    uint64_t n = rbinom(numPairs, acceptProbability);
    if (n == 0) {
        return;
    }
    stHash *hash = stHash_construct3(int32Key, int32EqualKey, free, free);
    int32_t *tmp = st_malloc(sizeof(*tmp));
    uint64_t m = 0;
    if (((double) n > numPairs / 2.0) && (numPairs > n)) {
        // sample (numPairs - n) many pairs
        m = numPairs - n;
    } else {
        // sample n many pairs
        m = n;
    }
    uint32_t i = 0;
    while (i < m) {
        *tmp = st_randomInt(0, numPairs);
        if (stHash_search(hash, tmp) == NULL) {
            stHash_insert(hash, copyInt32(tmp), copyInt32(tmp));
            ++i;
        }
    }
    stHashIterator *hit = stHash_getIterator(hash);
    int32_t *key = NULL;
    uint64_t p1, p2;
    if (m == n) {
        // items in hash have been sampled
        while ((key = stHash_getNext(hit)) != NULL) {
            // use mlArray
            arrayIndexToPairIndices((uint64_t)(*key), numSeqs, &p1, &p2);
            APair *aPair = aPair_construct(maf_mafLine_getSpecies(mlArray[p1]),
                                           maf_mafLine_getSpecies(mlArray[p2]),
                                           positions[p1], positions[p2]);
            avl_insert(pairs, aPair);
        }
    } else {
        // items not in hash have been sampled
        for (uint64_t i = 0; i < numSeqs; ++i) {
            *tmp = i;
            if (stHash_search(hash, tmp) == NULL) {
                arrayIndexToPairIndices(i, numSeqs, &p1, &p2);
                APair *aPair = aPair_construct(maf_mafLine_getSpecies(mlArray[p1]),
                                               maf_mafLine_getSpecies(mlArray[p2]),
                                               positions[p1], positions[p2]);
                avl_insert(pairs, aPair);
            }
        }
    }
    // clean up
    stHash_destructIterator(hit);
    stHash_destruct(hash);
    free(tmp);
}
void samplePairsFromColumnNaive(char **mat, uint32_t c, bool *legitRows, double acceptProbability, 
                                struct avl_table *pairs, 
                                uint32_t numRows, uint32_t numLegit, uint64_t *chooseTwoArray, 
                                mafLine_t **mlArray, uint32_t *positions, uint64_t numPairs) {
    // mat is the matrix of characters representing the alignment, c is the current column,
    // legtRows is a boolean array with true for a legit sequence, acceptProbability is the per
    // base accept probability, pairs is where we store pairs, numRows is the number of total sequences
    // in the matrix, numLegit is the number of legit sequences -- also the length
    // of the mlArray and positions arrays.
    uint64_t p1, p2;
    for (uint64_t i = 0; i < numPairs; ++i) {
        arrayIndexToPairIndices(i, numRows, &p1, &p2);
        if ((!legitRows[p1]) || (!legitRows[p2])) {
            continue;
        }
        if ((mat[p1][c] == '-') || (mat[p2][c]) == '-') {
            continue;
        }
        if (st_random() <= acceptProbability) {
            APair *aPair = aPair_construct(maf_mafLine_getSpecies(mlArray[p1]),
                                           maf_mafLine_getSpecies(mlArray[p2]),
                                           positions[p1], positions[p2]);
            avl_insert(pairs, aPair);
        }
    }
}
mafLine_t** createMafLineArray(mafBlock_t *mb, uint32_t numLegit, bool *legitRows) {
    if (numLegit == 0) {
        return NULL;
    }
    mafLine_t *ml = maf_mafBlock_getHeadLine(mb);
    uint32_t i, j;
    mafLine_t **mlArray = (mafLine_t**) st_malloc(sizeof(*mlArray) * numLegit);
    i = 0;
    j = 0;
    while (ml != NULL) {
        if (maf_mafLine_getType(ml) == 's') {
            if (legitRows[i]) {
                mlArray[j++] = ml;
            }
            ++i;
        }
        ml = maf_mafLine_getNext(ml);
    }
    return mlArray;
}
void updatePositions(char **mat, uint32_t c, uint32_t *positions, int *strandInts, uint32_t numSeqs) {
    for (uint32_t i = 0; i < numSeqs; ++i) {
        if (mat[i][c] != '-') {
            positions[i] += strandInts[i];
        }
    }
}
uint32_t* cullPositions(uint32_t *allPositions, uint32_t numSeqs, bool *legitRows, uint32_t numLegit) {
    uint32_t *positions = (uint32_t*) st_malloc(sizeof(*positions) * numLegit);
    uint32_t j = 0;
    for (uint32_t i = 0; i < numSeqs; ++i) {
        if (legitRows[i]) {
            positions[j++] = allPositions[i];
        }
    }
    return positions;
}
int* cullStrandInts(int *allStrandInts, uint32_t numSeqs, bool *legitRows, uint32_t numLegit) {
    int *strandInts = (int*) st_malloc(sizeof(*strandInts) * numLegit);
    uint32_t j = 0;
    for (uint32_t i = 0; i < numSeqs; ++i) {
        if (legitRows[i]) {
            strandInts[j++] = allStrandInts[i];
        }
    }
    return strandInts;
}
uint32_t sumBoolArray(bool *legitRows, uint32_t numSeqs) {
    uint32_t a = 0;
    for (uint32_t i = 0; i < numSeqs; ++i) {
        if (legitRows[i])
            ++a;
    }
    return a;
}
uint32_t countLegitPositions(char **mat, uint32_t c, uint32_t numRows) {
    uint32_t n = 0;
    for (uint32_t r = 0; r < numRows; ++r) {
        if (mat[r][c] != '-') {
            ++n;
        }
    }
    return n;
}
mafLine_t** cullMlArrayByColumn(char **mat, uint32_t c, uint32_t numRows, bool *legitRows, 
                                mafLine_t **mlArray, uint32_t numLegitGaplessPositions) {
    // create an array of mafLine_t for a given column, excluding all sequences that contain gaps 
    mafLine_t **colMlArray = (mafLine_t**) st_malloc(sizeof(*mlArray) * numLegitGaplessPositions);
    uint32_t j = 0;
    for (uint32_t r = 0; r < numRows; ++r) {
        if (!legitRows[r]) {
            continue;
        }
        if (mat[r][c] != '-') {
            colMlArray[j++] = mlArray[r];
        }
    }
    return colMlArray;
}
uint32_t* cullPositionsByColumn(char **mat, uint32_t c, uint32_t numRows, bool *legitRows,
                                uint32_t *positions, uint32_t numLegitGaplessPositions) {
    // create an array of positive coordinate position values that excludes all sequences that contain gaps
    uint32_t *colPositions = (uint32_t*) st_malloc(sizeof(*positions) * numLegitGaplessPositions);
    uint32_t j = 0;
    for (uint32_t r = 0; r < numRows; ++r) {
        if (!legitRows[r]) {
            continue;
        }
        if (mat[r][c] != '-') {
            colPositions[j++] = positions[r];
        }
    }
    return colPositions;
}
void walkBlockSamplingPairs(mafBlock_t *mb, struct avl_table *pairs,
                            double acceptProbability, stHash *legitSequences, 
                            uint64_t *chooseTwoArray) {
    uint32_t numSeqs = maf_mafBlock_getNumberOfSequences(mb);
    if (numSeqs < 2) {
        return;
    }
    uint32_t seqFieldLength = maf_mafBlock_longestSequenceField(mb);
    char **names = maf_mafBlock_getSpeciesArray(mb);
    char **mat = maf_mafBlock_getSequenceMatrix(mb, numSeqs, seqFieldLength);
    bool *legitRows = getLegitRows(names, numSeqs, legitSequences);
    uint32_t numLegit = sumBoolArray(legitRows, numSeqs);
    if (numLegit < 2) {
        return;
    }
    mafLine_t **mlArray = createMafLineArray(mb, numLegit, legitRows);
    uint32_t *allPositions = maf_mafBlock_getPosCoordStartArray(mb);
    int *allStrandInts = maf_mafBlock_getStrandIntArray(mb);
    int *strandInts = cullStrandInts(allStrandInts, numSeqs, legitRows, numLegit);
    for (uint32_t c = 0; c < seqFieldLength; ++c) {
        samplePairsFromColumn(mat, c, legitRows, acceptProbability, pairs, numSeqs, numLegit, chooseTwoArray, 
                              mlArray, allPositions);
        updatePositions(mat, c, allPositions, strandInts, numLegit);
    }
    // clean up
    free(mlArray);
    free(allPositions);
    free(allStrandInts);
    free(strandInts);
    for (uint32_t i = 0; i < numLegit; ++i) {
         free(names[i]);
    }
    free(names);
    maf_mafBlock_destroySequenceMatrix(mat, numSeqs);
    free(legitRows);
}
void samplePairsFromMaf(const char *filename, struct avl_table *pairs,
                        double acceptProbability, stHash *legitSequences) {
    mafFileApi_t *mfa = maf_newMfa(filename, "r");
    mafBlock_t *mb = NULL;
    uint64_t *chooseTwoArray = buildChooseTwoArray();
    while ((mb = maf_readBlock(mfa)) != NULL) {
        walkBlockSamplingPairs(mb, pairs, acceptProbability, legitSequences, chooseTwoArray);
        maf_destroyMafBlockList(mb);
    }
    // clean up
    free(chooseTwoArray);
    maf_destroyMfa(mfa);
}
void countPairs(APair *pair, stHash *intervalsHash, int64_t *counter, 
                stSortedSet *legitPairs, void *a, int32_t isVerboseFailures, int32_t near) {
    /*
     * Counts the number of pairs in the MAF file.
     */
    assert(pair != NULL);
    assert(a == NULL);
    if (stSortedSet_search(legitPairs, pair->seq1) != NULL) {
        if (stSortedSet_search(legitPairs, pair->seq2) != NULL) {
            (*counter)++;
        }
    }
}
void samplePairs(APair *thisPair, stHash *intervalsHash, struct avl_table *pairs, 
                 double *acceptProbability, stHash *legitPairs, int32_t isVerboseFailures, 
                 int32_t near) {
    /*
     * Adds *thisPair to *pairs with a given probability.
     */
    if (stHash_search(legitPairs, thisPair->seq1) != NULL)
        if (stHash_search(legitPairs, thisPair->seq2) != NULL)
            if (st_random() <= *acceptProbability)
                avl_insert(pairs, aPair_copyConstruct(thisPair));
}
bool inInterval(stHash *intervalsHash, char *seq, int32_t position) {
    /* 
     * 
     */
    stSortedSet *intervals = stHash_search(intervalsHash, seq);
    if (intervals == NULL) {
        return 0;
    }
    stIntTuple *i = stIntTuple_construct(2, position, INT32_MAX);
    stIntTuple *j = stSortedSet_searchLessThanOrEqual(intervals, i);
    stIntTuple_destruct(i);
    if (j == NULL) {
        return 0;
    }
    assert(stIntTuple_length(j) == 2);
    return stIntTuple_getPosition(j, 0) <= position && position < stIntTuple_getPosition(j, 1);
}
void getNearPairs(APair *thisPair, struct avl_table *pairs, int32_t near, stSortedSet *positivePairs) {
    /* given thisPair, if thisPair is in the table `pairs' then record it in the positivePairs set.
     * if the `near' option is set, do this not only for thisPair but for all pairs within +- `near'.
     * if near = 0 this will just look at thisPair->pos1 and thisPair->pos2 and record those values.
     */
    APair *aPair;
    int32_t i = thisPair->pos1;
    // Try modifying position 1
    for (thisPair->pos1 -= near; thisPair->pos1 < i + near + 1; thisPair->pos1++) {
        if ((aPair = avl_find(pairs, thisPair)) != NULL) {
            stSortedSet_insert(positivePairs, aPair);
        }
    }
    thisPair->pos1 = i;
    // Try modifying position 2
    i = thisPair->pos2;
    for (thisPair->pos2 -= near; thisPair->pos2 < i + near + 1; thisPair->pos2++) {
        if ((aPair = avl_find(pairs, thisPair)) != NULL) {
            stSortedSet_insert(positivePairs, aPair);
        }
    }
    thisPair->pos2 = i;
}
void homologyTests1(APair *thisPair, stHash *intervalsHash, struct avl_table *pairs, 
                    stSortedSet *positivePairs, stHash *legitPairs, int32_t isVerboseFailures, 
                    int32_t near) {
    /*
     * If both members of *thisPair are in the intersection of MAFFileA and MAFFileB,
     * and *thisPair is in the set *pairs then adds to the result pair a positive result.
     */
    if ((stHash_search(legitPairs, thisPair->seq1) != NULL)
        && (stHash_search(legitPairs, thisPair->seq2) != NULL)) {
        getNearPairs(thisPair, pairs, near, positivePairs);
    }
}
void homologyTests2(struct avl_table *pairs, struct avl_table *resultPairs, stHash *intervalsHash,
                    stSortedSet *positivePairs, int32_t isVerboseFailures) {
    /*
     * For every pair in 'pairs', add 1 to the total number of homology tests for the sequence-pair.
     * We don't bother with the seqNames hashtable here because the ht was used to build *pairs
     * in the first place.
     */
    static struct avl_traverser iterator;
    avl_t_init(&iterator, pairs);
    APair *pair;
    while ((pair = avl_t_next(&iterator)) != NULL) {
        ResultPair *resultPair = avl_find(resultPairs, pair);
        if (resultPair == NULL) {
            resultPair = resultPair_construct(pair->seq1, pair->seq2);
            avl_insert(resultPairs, resultPair);
        }
        bool b = stSortedSet_search(positivePairs, pair) != NULL;
        if (inInterval(intervalsHash, pair->seq1, pair->pos1)) {
            if (inInterval(intervalsHash, pair->seq2, pair->pos2)) {
                resultPair->totalBoth++;
                if (b) {
                    resultPair->inBoth++;
                }
            } else {
                resultPair->totalA++;
                if (b) {
                    resultPair->inA++;
                }
            }
        } else {
            if (inInterval(intervalsHash, pair->seq2, pair->pos2)) {
                resultPair->totalB++;
                if (b) {
                    resultPair->inB++;
                }
            } else {
                resultPair->totalNeither++;
                if (b) {
                    resultPair->inNeither++;
                }
            }
        }
        resultPair->total++;
        if (b) {
            resultPair->inAll++;
        } else {
           if (isVerboseFailures){
              fprintf(stderr, "%s\t%d\t\t%s\t%d\n", pair->seq1, pair->pos1, 
                      pair->seq2, pair->pos2);
           }
        }
    }
}
struct avl_table* compareMAFs_AB(const char *mafFileA, const char *mafFileB, uint32_t numberOfSamples,
                                 stHash *legitSequences, stHash *intervalsHash, 
                                 int32_t isVerboseFailures, uint32_t near) {
    uint64_t numberOfPairs = 0;
    // count the number of pairs in mafFileA
    numberOfPairs = countPairsInMaf(mafFileA, legitSequences);
    if (numberOfPairs == 0) {
        return avl_create((int32_t(*)(const void *, const void *, void *)) aPair_cmpFunction_seqsOnly, NULL, NULL);
    }
    double acceptProbability = ((double) numberOfSamples) / (double) numberOfPairs;
    struct avl_table *pairs = avl_create((int32_t(*)(const void *, const void *, void *)) aPair_cmpFunction, 
                                         NULL, NULL);
    // sample pairs from mafFileA
    samplePairsFromMaf(mafFileA, pairs, acceptProbability, legitSequences);
    // perform homology tests on mafFileB using sampled pairs from mafFileA
    stSortedSet *positivePairs = stSortedSet_construct();
    getPairs(mafFileB, (void(*)(APair *, stHash *, void *, void *, void *, int32_t, int32_t)) homologyTests1,
             intervalsHash, pairs, positivePairs, legitSequences, isVerboseFailures, near);
    struct avl_table *resultPairs = avl_create((int32_t(*)(const void *, const void *, void *)) aPair_cmpFunction_seqsOnly, NULL, NULL);
    homologyTests2(pairs, resultPairs, intervalsHash, positivePairs, isVerboseFailures);
    // clean up
    avl_destroy(pairs, (void(*)(void *, void *)) aPair_destruct);
    stSortedSet_destruct(positivePairs);
    return resultPairs;
}
void reportResult(const char *tagName, double total, double totalTrue, FILE *fileHandle, unsigned tabLevel) {
    assert(total >= totalTrue);
    findentprintf(fileHandle, tabLevel, "<%s totalTests=\"%" PRIi32 "\" totalTrue=\"%" PRIi32 "\" "
                  "totalFalse=\"%" PRIi32 "\" average=\"%f\"/>\n",
                  tagName, (int32_t) total, (int32_t) totalTrue, 
                  (int32_t) (total - totalTrue), total == 0 ? 0.0 : totalTrue / total);
}
ResultPair* aggregateResult(void *(*getNextPair)(void *, void *), void *arg1, void *arg2, 
                            const char *name1, const char *name2) {
    /* loop through all ResultPairs available via the getNextPair() iterator and aggregate their
     * results into a single ResultPair, return this struct.
     */
    ResultPair *resultPair;
    ResultPair *resultPair2 = resultPair_construct(name1, name2);
    while ((resultPair = getNextPair(arg1, arg2)) != NULL) {
        assert(resultPair->inAll == resultPair->inBoth + resultPair->inA + resultPair->inB + resultPair->inNeither);
        assert(resultPair->total == resultPair->totalBoth + resultPair->totalA + resultPair->totalB + resultPair->totalNeither);
        assert(resultPair->total >= resultPair->inAll);
        assert(resultPair->totalBoth >= resultPair->inBoth);
        assert(resultPair->totalA >= resultPair->inA);
        assert(resultPair->totalB >= resultPair->inB);
        assert(resultPair->totalNeither >= resultPair->inNeither);
        resultPair2->inAll += resultPair->inAll;
        resultPair2->inBoth += resultPair->inBoth;
        resultPair2->inA += resultPair->inA;
        resultPair2->inB += resultPair->inB;
        resultPair2->inNeither += resultPair->inNeither;
        resultPair2->total += resultPair->total;
        resultPair2->totalBoth += resultPair->totalBoth;
        resultPair2->totalA += resultPair->totalA;
        resultPair2->totalB += resultPair->totalB;
        resultPair2->totalNeither += resultPair->totalNeither;
    }
    return resultPair2;
}
void *addReferencesAndDups_getDups(void *arg, void *arg2) {
    ResultPair *resultPair;
    while ((resultPair = avl_t_next(arg)) != NULL) {
        if(strcmp(resultPair->aPair.seq1, resultPair->aPair.seq2) == 0) {
            break;
        }
    }
    return resultPair;
}
void *addReferencesAndDups_getReferences(void *arg, void *arg2) {
    ResultPair *resultPair;
    while ((resultPair = avl_t_next(arg)) != NULL) {
        if(strcmp(resultPair->aPair.seq1, arg2) == 0 || strcmp(resultPair->aPair.seq2, arg2) == 0) {
            break;
        }
    }
    return resultPair;
}
void addReferencesAndDups(struct avl_table *results_AB, stHash *legitSequences) {
    /*
     * Adds tags for all against each species.
     */
    stHashIterator *speciesIt = stHash_getIterator(legitSequences);
    static struct avl_traverser iterator;
    char *species;
    avl_t_init(&iterator, results_AB);
    stList *list = stList_construct();
    stList_append(list, aggregateResult(addReferencesAndDups_getDups, &iterator, NULL, "self", "self"));
    // Add references
    while((species = stHash_getNext(speciesIt)) != NULL) {
        avl_t_init(&iterator, results_AB);
        stList_append(list, aggregateResult(addReferencesAndDups_getReferences, &iterator, 
                                            species, species, "aggregate"));
    }
    stHash_destructIterator(speciesIt);
    for(int32_t i = 0; i < stList_length(list); i++) {
        avl_insert(results_AB, stList_get(list, i));
    }
    stList_destruct(list);
}
void *reportResults_fn(void *arg, void *arg2) {
   return avl_t_next(arg);
}
void findentprintf(FILE *fp, unsigned indent, char const *fmt, ...) {
    va_list args;
    va_start(args, fmt);
    char str[16384];
    int n = vsprintf(str, fmt, args);
        if (n >= 16384) {
            fprintf(stderr, "Error, failure in findentprintf, (n = %u) >= 16384\n", n);
            exit(EXIT_FAILURE);
        }
    for (unsigned i = 0; i < indent; ++i) {
        fprintf(fp, "\t");
    }
    fprintf(fp, str);
    va_end(args);
}
void reportResults(struct avl_table *results_AB, const char *mafFileA, const char *mafFileB, FILE *fileHandle,
                   uint32_t near, stHash *legitSequences, const char *bedFiles) {
    /*
     * Report results in an XML formatted document.
     */
    static struct avl_traverser iterator;
    avl_t_init(&iterator, results_AB);
    ResultPair *resultPair;
    ResultPair *aggregateResults = aggregateResult(reportResults_fn, &iterator, NULL, "", "");
    addReferencesAndDups(results_AB, legitSequences);
    findentprintf(fileHandle, 1, "<homologyTests fileA=\"%s\" fileB=\"%s\">\n", mafFileA, mafFileB);
    findentprintf(fileHandle, 2, "<aggregateResults>\n");
    reportResult("all", aggregateResults->total, aggregateResults->inAll, fileHandle, 3);
    if (bedFiles != NULL){
        reportResult("both", aggregateResults->totalBoth, aggregateResults->inBoth, fileHandle, 3);
        reportResult("A", aggregateResults->totalA, aggregateResults->inA, fileHandle, 3);
        reportResult("B", aggregateResults->totalB, aggregateResults->inB, fileHandle, 3);
        reportResult("neither", aggregateResults->totalNeither, aggregateResults->inNeither, fileHandle, 3);
    }
    findentprintf(fileHandle, 2, "</aggregateResults>\n");
    findentprintf(fileHandle, 2, "<homologyPairTests>\n");
    while ((resultPair = avl_t_prev(&iterator)) != NULL) {
        findentprintf(fileHandle, 3, "<homologyTest sequenceA=\"%s\" sequenceB=\"%s\">\n",
                      resultPair->aPair.seq1, resultPair->aPair.seq2);
        findentprintf(fileHandle, 4, "<aggregateResults>\n");
        reportResult("all", resultPair->total, resultPair->inAll, fileHandle, 5);
        if (bedFiles != NULL){
            reportResult("both", resultPair->totalBoth, resultPair->inBoth, fileHandle, 5);
            reportResult("A", resultPair->totalA, resultPair->inA, fileHandle, 5);
            reportResult("B", resultPair->totalB, resultPair->inB, fileHandle, 5);
            reportResult("neither", resultPair->totalNeither, resultPair->inNeither, fileHandle, 5);
        }
                      
        findentprintf(fileHandle, 4, "</aggregateResults>\n");
        findentprintf(fileHandle, 4, "<singleHomologyTests>\n");
        findentprintf(fileHandle, 5, "<singleHomologyTest sequenceA=\"%s\" sequenceB=\"%s\">\n", 
                      resultPair->aPair.seq1, resultPair->aPair.seq2);
        findentprintf(fileHandle, 6, "<aggregateResults>\n");
        reportResult("all", resultPair->total, resultPair->inAll, fileHandle, 6);
        if (bedFiles != NULL){
            reportResult("both", resultPair->totalBoth, resultPair->inBoth, fileHandle, 6);
            reportResult("A", resultPair->totalA, resultPair->inA, fileHandle, 6);
            reportResult("B", resultPair->totalB, resultPair->inB, fileHandle, 6);
            reportResult("neither", resultPair->totalNeither, resultPair->inNeither, fileHandle, 6);
        }
        findentprintf(fileHandle, 6, "</aggregateResults>\n");
        findentprintf(fileHandle, 5, "</singleHomologyTest>\n");
        findentprintf(fileHandle, 4, "</singleHomologyTests>\n");
        findentprintf(fileHandle, 3, "</homologyTest>\n");
    }
    findentprintf(fileHandle, 2, "</homologyPairTests>\n");
    findentprintf(fileHandle, 1, "</homologyTests>\n");
    resultPair_destruct(aggregateResults, NULL);
    return;
}
void populateNames(const char *filename, stHash *hash) {
    /*
     * populates a hash table with the names of sequences from a MAF file.
     */
    mafFileApi_t *mfa = maf_newMfa(filename, "r");
    mafBlock_t *mb = NULL;
    mafLine_t *ml = NULL;
    char *name = NULL;
    while ((mb = maf_readBlock(mfa)) != NULL) {
        ml = maf_mafBlock_getHeadLine(mb);
        while (ml != NULL) {
            if (maf_mafLine_getType(ml) == 's') {
                name = maf_mafLine_getSpecies(ml);
                if (stHash_search(hash, name) == NULL) {
                    stHash_insert(hash, stString_copy(name), stString_copy(""));
                }
            }
            ml = maf_mafLine_getNext(ml);
        }
        maf_destroyMafBlockList(mb);
    }
    // clean up
    maf_destroyMfa(mfa);
}
stHash* stHash_getIntersection(stHash *seqNames1, stHash *seqNames2) {
    stHash *hash = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, free, free);
    stHashIterator *hit = stHash_getIterator(seqNames1);
    char *key = NULL;
    while ((key = stHash_getNext(hit)) != NULL) {
        if (stHash_search(seqNames2, key) != NULL) {
            stHash_insert(hash, stString_copy(key), stString_copy(""));
        }
    }
    stHash_destructIterator(hit);
    return hash;
}
void printNameHash(struct hashtable *h) {
    // Debug function to iterate through a hash and print the contents.
    char *k;
    struct hashtable_itr *itr = NULL;
    if (hashtable_count(h) > 0) {
        itr = hashtable_iterator(h);
        do {
            k = hashtable_iterator_key(itr);
            fprintf(stderr, "%s\n", k);
        } while (hashtable_iterator_advance(itr));
    }
    if (itr != NULL) {
        free(itr);
    }
}
void intersectHashes(struct hashtable *h1, struct hashtable *h2, struct hashtable *h3) {
    // intersects two hashes. (hash 1 \cap  hash 2) = hash 3.
    char *k;
    struct hashtable_itr *itr = NULL;
    if (hashtable_count(h1) > 0) {
        itr = hashtable_iterator(h1);
        do {
            k = hashtable_iterator_key(itr);
            if (hashtable_search(h2, k) != NULL) {
                hashtable_insert(h3, k, constructInt(1));
            }
        } while (hashtable_iterator_advance(itr));
    }
    if (itr != NULL) {
        free(itr);
    }
}
int32_t countNodes(stTree *node) {
    int32_t i = 0;
    int32_t count = 0;
    if (node == NULL) {
        return 0;
    } else {
        count = 1;
        for (i = 0; i < stTree_getChildNumber(node); i++) {
            count += countNodes(stTree_getChild(node, i));
        }
        return count;
    }
}
int32_t countLeaves(stTree *node) {
    int32_t i = 0;
    int32_t count = 0;
    if (node == NULL) {
        return 0;
    } else {
        if (stTree_getChildNumber(node) == 0) {
            return 1;
        } else {
            for (i = 0; i < stTree_getChildNumber(node); i++) {
                count += countLeaves(stTree_getChild(node, i));
            }
            return count;
        }
    }
}
void postOrderLabelNodes(stTree *node, int32_t *index, char **labelArray) {
    int32_t i = 0;
    for (i = 0; i < stTree_getChildNumber(node); i++) {
        postOrderLabelNodes(stTree_getChild(node, i), index, labelArray);
    }
    labelArray[*index] = stString_copy(stTree_getLabel(node));
    *index += 1;
    return;
}
void postOrderLabelLeaves(stTree *node, int32_t *index, char **labelArray) {
    int32_t i = 0;
    for (i = 0; i < stTree_getChildNumber(node); i++) {
        postOrderLabelLeaves(stTree_getChild(node, i), index, labelArray);
    }
    if (stTree_getChildNumber(node) == 0) {
        labelArray[*index] = stString_copy(stTree_getLabel(node));
        *index += 1;
    }
    return;
}
void lcaP(stTree *node, struct djs *uf, int32_t *ancestor, int32_t *color, 
          struct hashtable *ht, int32_t size, int32_t **lcaMatrix) {
    /*
     * 
     */
    if (node == NULL) {
        return;
    }
    int32_t u = *((int *) hashtable_search(ht, (void *) stTree_getLabel(node)));
    djs_makeset(uf, u);
    ancestor[djs_findset(uf, u)] = u;
    int32_t v;
    int32_t i = 0;
    for (i = 0; i < stTree_getChildNumber(node); i++) {
        lcaP(stTree_getChild(node, i), uf, ancestor, color, ht, size, lcaMatrix);
        v = *((int *) hashtable_search(ht, (void *) stTree_getLabel(stTree_getChild(node, i))));
        djs_union(uf, u, v);
        ancestor[djs_findset(uf, u)] = u;
    }
    color[u] = 1;
    for (i = 0; i < size; i++) {
        if (color[i] == 1) {
            lcaMatrix[u][i] = ancestor[djs_findset(uf, i)];
            lcaMatrix[i][u] = ancestor[djs_findset(uf, i)];
        }
    }
    return;
}
int32_t **lca(stTree *root, struct hashtable *ht) {
    int32_t nodeNum = countNodes(root);
    struct djs *uf = NULL;
    uf = djs_new(nodeNum);
    int32_t i = 0;
    int32_t *ancestor = NULL;
    ancestor = st_malloc(sizeof(int32_t) * nodeNum);
    for (i = 0; i < nodeNum; i++) {
        ancestor[i] = 0;
    }
    int32_t *color = NULL;
    color = st_malloc(sizeof(int32_t) * nodeNum);
    for (i = 0; i < nodeNum; i++) {
        color[i] = 0;
    }
    int32_t **lcaMatrix = NULL;
    lcaMatrix = st_malloc(sizeof(int32_t *) * nodeNum);
    for (i = 0; i < nodeNum; i++) {
        lcaMatrix[i] = st_malloc(sizeof(int32_t) * nodeNum);
    }
    lcaP(root, uf, ancestor, color, ht, nodeNum, lcaMatrix);
    djs_free(uf);
    free(ancestor);
    free(color);
    return lcaMatrix;
}
void lcaMatrix_free(int32_t **lcaMatrix, int32_t nodeNum) {
    int32_t i = 0;
    for (i = 0; i < nodeNum; i++) {
        free(lcaMatrix[i]);
    }
    free(lcaMatrix);
    return;
}
char **createNodeLabelArray(stTree *tree, int32_t nodeNum) {
    char **labelArray = st_malloc(sizeof(char *) * nodeNum);
    int32_t po_index = 0;
    postOrderLabelNodes(tree, &po_index, labelArray);
    return labelArray;
}
char **createLeafLabelArray(stTree *tree, int32_t nodeNum) {
    char **leafLabelArray = st_malloc(sizeof(char *) * nodeNum);
    int32_t po_index = 0;
    postOrderLabelLeaves(tree, &po_index, leafLabelArray);
    return leafLabelArray;
}
void labelArray_destruct(char **labelArray, int32_t num) {
    int32_t i = 0;
    for (i = 0; i < num; i++) {
        free(labelArray[i]);
    }
    free(labelArray);
    return;
}
struct hashtable * getTreeLabelHash(char **labelArray, int32_t nodeNum) {
    struct hashtable *treeLabelHash = NULL;
    int32_t i = 0;
    treeLabelHash = create_hashtable(256, hashtable_stringHashKey, hashtable_stringEqualKey, free, free);
    for (i = 0; i < nodeNum; i++) {
        hashtable_insert(treeLabelHash, stString_copy(labelArray[i]), constructInt(i));
    }
    return treeLabelHash;
}
void writeXMLHeader(FILE *fileHandle){
   fprintf(fileHandle, "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\" ?>\n");
   return;
}
