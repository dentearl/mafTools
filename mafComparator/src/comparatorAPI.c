/*
 * Copyright (C) 2009-2013 by
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

#include <math.h>
#include <string.h>
#include "sonLib.h"
#include "common.h"
#include "comparatorAPI.h"
#include "comparatorRandom.h"

const unsigned kChooseTwoCacheLength = 101;
bool g_isVerboseFailures = false;

void aPair_fillOut(APair *aPair, char *seq1, char *seq2, uint64_t pos1, uint64_t pos2) {
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
APair* aPair_init(void) {
    APair *aPair = st_malloc(sizeof(*aPair));
    aPair->seq1 = NULL;
    aPair->seq2 = NULL;
    aPair->pos1 = 0;
    aPair->pos2 = 0;
    return aPair;
}
WiggleContainer* wiggleContainer_init(void) {
    WiggleContainer *wc = st_malloc(sizeof(*wc));
    wc->ref = NULL;
    wc->partner = NULL;
    wc->refStart = 0;
    wc->refLength = 0;
    wc->numBins = 0;
    wc->binLength = 0;
    wc->presentAtoB = NULL;
    wc->presentBtoA = NULL;
    wc->absentAtoB = NULL;
    wc->absentBtoA = NULL;
    return wc;
}
Options* options_construct(void) {
    Options *o = (Options*) st_malloc(sizeof(*o));
    o->logLevelString = NULL;
    o->mafFile1 = NULL;
    o->mafFile2 = NULL;
    o->outputFile = NULL;
    o->bedFiles = NULL;
    o->wigglePairs = NULL;
    o->wiggleRegionStart = 0;
    o->wiggleRegionStop = 0;
    o->numPairsString = NULL;
    o->legitSequences = NULL;
    o->numberOfSamples = 1000000; // by default do a million samples per file pair.
    o->randomSeed = (time(NULL) << 16) | (getpid() & 65535); // Likely to be unique
    o->near = 0;
    o->numPairs1 = 0;
    o->numPairs2 = 0;
    o->wiggleBinLength = 100000; // by default have bins of length 100,000
    return o;
}
APair* aPair_construct(const char *seq1, const char *seq2, uint64_t pos1, uint64_t pos2) {
    APair *aPair = aPair_init();
    aPair_fillOut(aPair, (char *) seq1, (char *) seq2, pos1, pos2);
    aPair->seq1 = stString_copy(aPair->seq1);
    aPair->seq2 = stString_copy(aPair->seq2);
    return aPair;
}
APair* aPair_copyConstruct(APair *pair) {
    return aPair_construct(stString_copy(pair->seq1), stString_copy(pair->seq2), pair->pos1, pair->pos2);
}
ResultPair *resultPair_construct(const char *seq1, const char *seq2) {
    ResultPair *resultPair = st_calloc(1, sizeof(ResultPair));
    resultPair->seq1 = stString_copy(seq1);
    resultPair->seq2 = stString_copy(seq2);
    return resultPair;
}
WiggleContainer* wiggleContainer_construct(char *ref, char *partner, uint64_t refStart,
                                           uint64_t refLength, uint64_t wiggleBinLength) {
    WiggleContainer *wc = wiggleContainer_init();
    wc->ref = stString_copy(ref);
    wc->partner = stString_copy(partner);
    wc->refStart = refStart;
    wc->refLength = refLength;
    wc->numBins = (uint64_t)ceil((double)refLength / wiggleBinLength);
    wc->binLength = wiggleBinLength;
    wc->presentAtoB = st_calloc(wc->numBins, sizeof(uint64_t));
    wc->presentBtoA = st_calloc(wc->numBins, sizeof(uint64_t));
    wc->absentAtoB = st_calloc(wc->numBins, sizeof(uint64_t));
    wc->absentBtoA = st_calloc(wc->numBins, sizeof(uint64_t));
    return wc;
}
void aPosition_fillOut(APosition *aPosition, char *name, uint64_t pos) {
    aPosition->name = name;
    aPosition->pos = pos;
}
APosition* aPosition_init(void) {
    APosition *aPosition = st_malloc(sizeof(*aPosition));
    aPosition->name = NULL;
    aPosition->pos = 0;
    return aPosition;
}
APosition* aPosition_construct(const char *name, uint64_t pos) {
    APosition *aPosition = st_malloc(sizeof(*aPosition));
    aPosition_fillOut(aPosition, (char *) name, pos);
    return aPosition;
}
void options_destruct(Options *o) {
    if (o == NULL) {
        return;
    }
    free(o->logLevelString);
    free(o->mafFile1);
    free(o->mafFile2);
    free(o->outputFile);
    free(o->bedFiles);
    free(o->wigglePairs);
    free(o->legitSequences);
    free(o->numPairsString);
    free(o);
    o = NULL;
}
void aPair_destruct(APair *pair) {
    if (pair == NULL) {
        return;
    }
    if (pair->seq1 != NULL) {
        free(pair->seq1);
    }
    if (pair->seq2 != NULL) {
        free(pair->seq2);
    }
    free(pair);
    pair = NULL;
}
void aPosition_destruct(void *p) {
    if (p == NULL) {
        return;
    }
    free(((APosition*)p)->name);
    free(p);
    p = NULL;
}
void resultPair_destruct(ResultPair *rp) {
    if (rp == NULL) {
        return;
    }
    free(rp->seq1);
    free(rp->seq2);
    free(rp);
    rp = NULL;
}
void wiggleContainer_destruct(WiggleContainer *wc) {
    free(wc->ref);
    free(wc->partner);
    free(wc->presentAtoB);
    free(wc->presentBtoA);
    free(wc->absentAtoB);
    free(wc->absentBtoA);
    free(wc);
    wc = NULL;
}
int aPair_cmpFunction_seqsOnly(APair *aPair1, APair *aPair2) {
    /*
     * Compares only the sequence name arguments of the a pairs.
     */
    int i;
    if (aPair1->seq1 == NULL && aPair2->seq1 == NULL) {
        i = 0;
    } else if (aPair2->seq1 == NULL) {
        i = 1;
    } else if (aPair1->seq1 == NULL) {
        i = -1;
    } else {
        i = strcmp(aPair1->seq1, aPair2->seq1);
    }
    if (i == 0) {
        // NOTE THAT WE ALLOW seq2 to be NULL! This is solely because of
        // the search that takes place in testHomolygOnColumn()
        if (aPair1->seq2 != NULL) {
            if (aPair2->seq2 != NULL) {
                i = strcmp(aPair1->seq2, aPair2->seq2);
            } else {
                // NULL entries will be considered equal, sequence name wise
                i = 0;
            }
        } else {
            i = 0;
        }
    }
    return i;
}
int strcmpnull(const char *s1, const char *s2) {
    // do strcmp with NULL strings. Totally problematic.
    int i;
    if ((s1 == NULL) && (s2 == NULL)) {
        i = 0;
    } else if (s2 == NULL) {
        i = 1;
    } else if (s1 == NULL) {
        i = -1;
    } else {
        i = strcmp(s1, s2);
    }
    return i;
}
int aPair_cmpFunction(APair *p1, APair *p2) {
    /*
     * Compares first the sequence then the first position, then the second sequence,
     * then the second position. Allows us to traverse the pairs sortedSet quickly
     */
    int i = strcmpnull(p1->seq1, p2->seq1);
    if (i == 0) {
        if (p1->pos1 < p2->pos1) {
            i = -1;
        } else if (p1->pos1 > p2->pos1) {
            i = 1;
        } else {
            i = strcmpnull(p1->seq2, p2->seq2);
            if (i == 0) {
                if (p1->pos2 < p2->pos2) {
                    i = -1;
                } else if (p1->pos2 > p2->pos2) {
                    i = 1;
                } else {
                    i = 0;
                }
            }
        }
    }
    return i;
}
uint64_t aPositionKey(const void *p) {
    // maybe not the best composite hash function,
    APosition *pos = (APosition*) p;
    return stHash_stringKey((const void*) pos->name) + pos->pos;
}
int aPositionEqualKey(const void *k1, const void *k2) {
    // return 1 if equal, 0 if not equal.
    APosition *ik1 = (APosition*) k1;
    APosition *ik2 = (APosition*) k2;
    if (!strcmp(ik1->name, ik2->name)) {
        if (ik1->pos == ik2->pos) {
            return 1;
        }
    }
    return 0;
}
bool closeEnough(uint64_t p1, uint64_t p2, uint64_t near) {
    if (p1 == p2) {
        return true;
    } else if (p1 < p2 ) {
        if (p1 + near >= p2) {
            return true;
        } else {
            return false;
        }
    } else {
        // (p1 > p2)
        if (p1 <= p2 + near) {
            return true;
        } else {
            return false;
        }
    }
}
int64_t* buildInt(int64_t n) {
    int64_t *p = (int64_t*) st_malloc(sizeof(*p));
    *p = n;
    return p;
}
int64_t* buildInt64(int64_t n) {
    int64_t *p = (int64_t*) st_malloc(sizeof(*p));
    *p = n;
    return p;
}
uint64_t* buildUInt64(uint64_t n) {
    uint64_t *p = (uint64_t*) st_malloc(sizeof(*p));
    *p = n;
    return p;
}
void pairIndicesToArrayIndex(uint64_t r, uint64_t c, uint64_t n, uint64_t *i) {
    // r = row, c = column, n = number of sequences
    /* A proof by picture:
             c
             0 1 2 3 4
            __________
       r 0 | - 0 1 2 3
         1 |   - 4 5 6
         2 |     - 7 8
         3 |       - 9
         4 |         -

       so (3, 4) can expand row by row:
       = (n - 1) + (n - 2) + (n - 3) + c - r - 1
       = 3 * n - (1 + 2 + 3) + c - r - 1
       = r * n - (sum_{i=1}^{r} i) + c - r - 1
       = r * n - (r (r + 1) / 2) + c - r - 1
       = r * (n - 1) - (r (r + 1) / 2) + c - 1
       = 3 * 5 - (3 * 4 / 2) + 4 - 3 - 1
       = 15 - 6 + 4 - 3 - 1
       = 9
       so (1, 3) yeilds 5 by 4 + 2 - 1, the 4 comes from 5 - sum_{i=1}^{r} i = r(r+1)/2
       and the 2 comes from c - r.
       Essentially, this function is the inverse of arrayIndexToPairIndices().
     */
    assert(n > 1);
    if (n == 2) {
        *i = 0;
        return;
    }
    *i = (r * (n - 1) - (r * (r + 1) / 2) + c - 1);
}
void arrayIndexToPairIndices(uint64_t i, uint64_t n, uint64_t *r, uint64_t *c) {
    // i is the index, n is the length of one side of the square (i.e. the number
    // of sequences),  r and c are pointers to record the row and column (respectively)
    // of the index. Essentially, this function is the inverse of  pairIndicesToArrayIndex().
    /*
                   c
             0 1 2 3 4
            __________
       r 0 | - 0 1 2 3
         1 |   - 4 5 6
         2 |     - 7 8
         3 |       - 9
         4 |         -

      Our solution for finding the row is based on the fact that by counting from
      lower rows to upper rows the number of elements is 1 + 2 + 3 + 4... with the
      general form sum_{j=1}^{m} j, which simplifies to m*(m-1)/2.
      We find the row, r, first by way of x, where
      x is the number of positions greater than i
      x = n * (n - 1) / 2 - i - 1
      the term n*(n-1)/2 is the total number of elements, the -1 transforms
      from 1 based to 0 based coordinates.
      Now we wish to find the largest k such that k(k-1)/2 <= x
      k * (k - 1) / 2 <= x
            k * (k-1) <= 2 * x
    (k - 1/2)^2 - 1/4 <= 2x
          (k - 1/2)^2 <= 2x + 1/4
    2^2 * (k - 1/2)^2 <= 8 * x + 1
        2 * (k - 1/2) <= sqrt(8 * x + 1)
                2 * k <= sqrt(8 * x + 1) + 1
                    k <= (sqrt(8 * x + 1) + 1) / 2
                    we seek the largest k,
                    k = floor((sqrt(8 * x + 1) + 1) / 2)
      from here we find r as the difference from the number of rows and k, with a 1 to 0 transform
      r = n - k - 1
      Finally we find the column index, c, given i, n and r.
      c = i + 1 - (r * (n - 1) - r * (r + 1) / 2)
      where the +1 term shifts the column index from 0 to 1 (see picture above, c can never be 0).
      the term (r * (n - 1) - r * (r + 1) / 2) provides the column offset as a function of the row.
     */
    assert(n > 1);
    if (n == 2) {
        *r = 0;
        *c = 1;
        return;
    }
    if (i >= ((n * (n - 1)) / 2)) {
        fprintf(stderr, "Bad ju-ju, i:%" PRIu64 ", n:%" PRIu64", ((n*(n-1))/2:%" PRIu64 "\n", i, n, ((n*(n-1))/2));
    }
    assert(i < ((n * (n - 1)) / 2));
    uint64_t x = (n * (n - 1)) / 2 - i - 1;
    uint64_t k = floor((sqrt(8 * x + 1) + 1) / 2);
    *r = n - k - 1;
    *c = i - ((*r) * (n - 1) - ((*r) * ((*r) + 1)) / 2) + 1;
    assert((*r) < n);
    assert((*c) < n);
}
bool* getLegitRows(char **names, uint64_t numSeqs, stSet *legitSequences) {
    bool *legitRows = (bool *) st_malloc(sizeof(*legitRows) * numSeqs);
    for (uint64_t i = 0; i < numSeqs; ++i) {
        if (legitSequences != NULL) {
            if (stSet_search(legitSequences, names[i]) != NULL) {
                legitRows[i] = true;
            } else {
                legitRows[i] = false;
            }
        } else {
            legitRows[i] = true;
        }
    }
    return legitRows;
}
uint64_t countPairsInColumn(char **mat, uint64_t c, uint64_t numSeqs,
                            bool *legitRows, uint64_t *chooseTwoArray) {
    uint64_t possiblePartners = 0;
    for (uint64_t r = 0; r < numSeqs; ++r) {
        if (!legitRows[r]) {
            continue;
        }
        if (mat[r][c] != '-') {
            ++possiblePartners;
        }
    }
    if (possiblePartners < kChooseTwoCacheLength) {
        return chooseTwoArray[possiblePartners];
    } else {
        return chooseTwo(possiblePartners);
    }
}
uint64_t walkBlockCountingPairs(mafBlock_t *mb, stSet *legitSequences, uint64_t *chooseTwoArray) {
    // size is the MAXIMUM INDEX of the chooseTwo array
    uint64_t count = 0;
    uint64_t numSeqs = maf_mafBlock_getNumberOfSequences(mb);
    if (numSeqs < 1) {
        return 0;
    }
    uint64_t seqFieldLength = maf_mafBlock_getSequenceFieldLength(mb);
    char **names = maf_mafBlock_getSpeciesArray(mb);
    char **mat = maf_mafBlock_getSequenceMatrix(mb, numSeqs, seqFieldLength);
    bool *legitRows = getLegitRows(names, numSeqs, legitSequences);
    for (uint64_t c = 0; c < seqFieldLength; ++c) {
        count += countPairsInColumn(mat, c, numSeqs, legitRows, chooseTwoArray);
    }
    // clean up
    for (uint64_t i = 0; i < numSeqs; ++i) {
        free(names[i]);
    }
    free(names);
    maf_mafBlock_destroySequenceMatrix(mat, numSeqs);
    free(legitRows);
    return count;
}
uint64_t chooseTwo(uint64_t n) {
    if ((n == 0) || (n == 1)) {
        return 0;
    } else {
        return ((n * (n - 1)) / 2);
    }
}
uint64_t* buildChooseTwoArray(void) {
    // pre-calculate a bunch of smaller sizes
    uint64_t *cta = (uint64_t *) st_malloc(sizeof(*cta) * kChooseTwoCacheLength);
    for (uint64_t i = 0; i < kChooseTwoCacheLength; ++i) {
        cta[i] = chooseTwo(i);
    }
    return cta;
}
uint64_t countPairsInMaf(const char *filename, stSet *legitSequences) {
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
static uint64_t uint64Key(const void *k) {
    uint64_t p = *(uint64_t*)k;
    if (p > INT64_MAX) {
        while (p > INT64_MAX) {
            p -= INT64_MAX;
        }
    }
    return p;
}
uint64_t* uint64Copy(uint64_t *i) {
    uint64_t *j = st_malloc(sizeof(*j));
    *j = *i;
    return j;
}
static int uint64EqualKey(const void *key1, const void *key2) {
    return *((uint64_t *) key1) == *((uint64_t*) key2);
}
void printmlarray(mafLine_t **mlArray, uint64_t n) {
    for (uint64_t i = 0; i < n; ++i) {
        printf("%" PRIu64 ":%s, ", i, maf_mafLine_getSpecies(mlArray[i]));
    }
    printf("\n");
}
uint64_t countLegitGaplessPositions(char **mat, uint64_t c, uint64_t numRows, bool *legitRows) {
    uint64_t a = 0;
    for (uint64_t r = 0; r < numRows; ++r) {
        if (!legitRows[r]) {
            continue;
        }
        if (mat[r][c] != '-') {
            ++a;
        }
    }
    return a;
}
void samplePairsFromColumn(double acceptProbability,
                           stSortedSet *pairs, uint64_t numSeqs, uint64_t *chooseTwoArray,
                           char **nameArray, uint64_t *columnPositions) {
    // acceptProbability is the per base accept probability, pairs is where we store pairs,
    // numSeqs is the number of sequences in either array, columnMlArray is an array that contains
    // pointers to mafLine_t's and columnPositions is an array that contains the current position
    // of the sequence
    uint64_t numPairs;
    if (numSeqs < kChooseTwoCacheLength) {
        numPairs = chooseTwoArray[numSeqs];
    } else {
        numPairs = chooseTwo(numSeqs);
    }
    if (numPairs < 5) {
        // keep in mind the sequence of values for v = k choose 2 is
        // k: 2 3 4  5 ...
        // v: 1 3 6 10 ...
        samplePairsFromColumnBruteForce(acceptProbability, pairs,
                                        chooseTwoArray, nameArray,
                                        columnPositions, numSeqs, numPairs);
    } else {
        samplePairsFromColumnAnalytic(acceptProbability, pairs,
                                      chooseTwoArray, nameArray,
                                      columnPositions, numSeqs, numPairs);
    }
    free(nameArray);
    free(columnPositions);
}
void samplePairsFromColumnBruteForce(double acceptProbability, stSortedSet *pairs,
                                     uint64_t *chooseTwoArray,
                                     char **nameArray, uint64_t *positions, uint64_t numSeqs,
                                     uint64_t numPairs) {
    uint64_t p1, p2;
    for (uint64_t i = 0; i < numPairs; ++i) {
        if (st_random() <= acceptProbability) {
            arrayIndexToPairIndices(i, numSeqs, &p1, &p2);
            // printf("1. adding pair (%s %u):(%s %u)\n", maf_mafLine_getSpecies(mlArray[p1]),
            //        positions[p1], maf_mafLine_getSpecies(mlArray[p2]), positions[p2]);
            APair *aPair = aPair_construct(nameArray[p1], nameArray[p2],
                                           positions[p1], positions[p2]);
            stSortedSet_insert(pairs, aPair);
        }
    }
}
void samplePairsFromColumnAnalytic(double acceptProbability, stSortedSet *pairs,
                                   uint64_t *chooseTwoArray,
                                   char **nameArray, uint64_t *positions, uint64_t numSeqs,
                                   uint64_t numPairs) {
    uint64_t n = rbinom(numPairs, acceptProbability);
    if (n == 0) {
        return;
    }
    stSet *set = stSet_construct3(uint64Key, uint64EqualKey, free);
    uint64_t *randPair = st_malloc(sizeof(*randPair));
    uint64_t numPairsToSample = 0;
    int offset = 0; // used when numPairs > INT64_MAX, offset is the number of times to multiply by INT64_MAX
    (void) offset;
    if (((double) n > numPairs / 2.0) && (numPairs > n)) {
        // sample (numSeqs - n) many pairs
        numPairsToSample = numPairs - n;
    } else {
        // sample n many pairs
        numPairsToSample = n;
    }
    uint64_t i = 0;
    if (numPairs > UINT64_MAX) {
        // panic!
        fprintf(stderr, "Error in samplePairsFromColumnAnalytic(), numPairs (%" PRIu64
                ") > UINT64_MAX (%" PRIu64 "), this will cause st_randomInt64() to fail.\n", numPairs, UINT64_MAX);
        exit(EXIT_FAILURE);
    } else if (numPairs > INT64_MAX) {
        // will have to sample from a shifted range using st_randomInt64()
        int64_t imin = INT64_MAX - numPairs; // shift range
        assert(imin < 0);
        int64_t randi;
        while (i < numPairsToSample) {
            randi = st_randomInt64(imin, INT64_MAX); // sample from shifted range
            assert((randi - imin) >= 0);
            *randPair = randi - imin; // shift range back
            assert(*randPair <= numPairs);
            if (stSet_search(set, randPair) == NULL) {
                stSet_insert(set, uint64Copy(randPair));
                ++i;
            }
        }
    } else if (numPairs > INT64_MAX){
        // sample straight away using st_randomInt64()
        while (i < numPairsToSample) {
            *randPair = st_randomInt64(0, numPairs);
            assert(*randPair <= numPairs);
            if (stSet_search(set, randPair) == NULL) {
                stSet_insert(set, uint64Copy(randPair));
                ++i;
            }
        }
    } else {
        // sample straight away using st_randomInt()
        while (i < numPairsToSample) {
            *randPair = st_randomInt(0, numPairs);
            assert(*randPair <= numPairs);
            if (stSet_search(set, randPair) == NULL) {
                stSet_insert(set, uint64Copy(randPair));
                ++i;
            }
        }
    }
    stSetIterator *sit = stSet_getIterator(set);
    uint64_t *key = NULL;
    uint64_t p1, p2;
    if (numPairsToSample == n) {
        // items in set *have* been sampled
        while ((key = stSet_getNext(sit)) != NULL) {
            // use nameArray
            arrayIndexToPairIndices(*key, numSeqs, &p1, &p2);
            APair *aPair = aPair_construct(nameArray[p1], nameArray[p2],
                                           positions[p1], positions[p2]);
            stSortedSet_insert(pairs, aPair);
        }
    } else {
        // items in set *have not* been sampled
        for (i = 0; i < numPairs; ++i) {
            *randPair = i;
            if (stSet_search(set, randPair) == NULL) {
                arrayIndexToPairIndices(i, numSeqs, &p1, &p2);
                APair *aPair = aPair_construct(nameArray[p1], nameArray[p2],
                                               positions[p1], positions[p2]);
                stSortedSet_insert(pairs, aPair);
            }
        }
    }
    // clean up
    stSet_destructIterator(sit);
    stSet_destruct(set);
    free(randPair);
}
void samplePairsFromColumnNaive(char **mat, uint64_t c, bool *legitRows, double acceptProbability,
                                stSortedSet *pairs,
                                uint64_t *chooseTwoArray,
                                char **nameArray, uint64_t *positions, uint64_t numSeqs,
                                uint64_t numPairs) {
    // mat is the matrix of characters representing the alignment, c is the current column,
    // legtRows is a boolean array with true for a legit sequence, acceptProbability is the per
    // base accept probability, pairs is where we store pairs, numRows is the number of total sequences
    // in the matrix, numLegit is the number of legit sequences -- also the length
    // of the nameArray and positions arrays.
    uint64_t p1, p2;
    for (uint64_t i = 0; i < numPairs; ++i) {
        arrayIndexToPairIndices(i, numSeqs, &p1, &p2);
        if ((!legitRows[p1]) || (!legitRows[p2])) {
            continue;
        }
        if ((mat[p1][c] == '-') || (mat[p2][c]) == '-') {
            continue;
        }
        if (st_random() <= acceptProbability) {
            APair *aPair = aPair_construct(nameArray[p1], nameArray[p2],
                                           positions[p1], positions[p2]);
            stSortedSet_insert(pairs, aPair);
        }
    }
}
mafLine_t** createMafLineArray(mafBlock_t *mb, uint64_t numLegit, bool *legitRows) {
    if (numLegit == 0) {
        return NULL;
    }
    mafLine_t *ml = maf_mafBlock_getHeadLine(mb);
    uint64_t i = 0, j = 0;
    mafLine_t **mlArray = (mafLine_t**) st_malloc(sizeof(*mlArray) * numLegit);
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
void updatePositions(char **mat, uint64_t c, uint64_t *allPositions, int *allStrandInts, uint64_t numSeqs) {
    for (uint64_t i = 0; i < numSeqs; ++i) {
        if (mat[i][c] != '-') {
            allPositions[i] += allStrandInts[i];
        }
    }
}
uint64_t* cullPositions(uint64_t *allPositions, uint64_t numSeqs, bool *legitRows, uint64_t numLegit) {
    uint64_t *positions = (uint64_t*) st_malloc(sizeof(*positions) * numLegit);
    uint64_t j = 0;
    for (uint64_t i = 0; i < numSeqs; ++i) {
        if (legitRows[i]) {
            positions[j++] = allPositions[i];
        }
    }
    return positions;
}
int* cullStrandInts(int *allStrandInts, uint64_t numSeqs, bool *legitRows, uint64_t numLegit) {
    int *strandInts = (int*) st_malloc(sizeof(*strandInts) * numLegit);
    uint64_t j = 0;
    for (uint64_t i = 0; i < numSeqs; ++i) {
        if (legitRows[i]) {
            strandInts[j++] = allStrandInts[i];
        }
    }
    return strandInts;
}
uint64_t sumBoolArray(bool *legitRows, uint64_t numSeqs) {
    uint64_t a = 0;
    for (uint64_t i = 0; i < numSeqs; ++i) {
        if (legitRows[i])
            ++a;
    }
    return a;
}
uint64_t countLegitPositions(char **mat, uint64_t c, uint64_t numRows) {
    uint64_t n = 0;
    for (uint64_t r = 0; r < numRows; ++r) {
        if (mat[r][c] != '-') {
            ++n;
        }
    }
    return n;
}
mafLine_t** cullMlArrayByColumn(char **mat, uint64_t c, mafLine_t **mlArray, bool *legitRows,
                                uint64_t numRows, uint64_t numLegitGaplessPositions) {
    // create an array of mafLine_t for a given column, excluding all sequences that contain gaps
    mafLine_t **colMlArray = (mafLine_t**) st_malloc(sizeof(*mlArray) * numLegitGaplessPositions);
    uint64_t j = 0;
    for (uint64_t r = 0; r < numRows; ++r) {
        if (legitRows[r] && mat[r][c] != '-') {
            colMlArray[j++] = mlArray[r];
        }
    }
    return colMlArray;
}
char** extractLegitGaplessNamesFromMlArrayByColumn(char **mat, uint64_t c, mafLine_t **mlArray, bool *legitRows,
                                                   uint64_t numRows, uint64_t numLegitGaplessPositions) {
    // winner of longest function name award
    // create an array of mafLine_t for a given column, excluding all sequences that contain gaps
    char **nameArray = (char **) st_malloc(sizeof(*nameArray) * numLegitGaplessPositions);
    uint64_t j = 0;
    for (uint64_t r = 0; r < numRows; ++r) {
        if (legitRows[r] && mat[r][c] != '-') {
            // NOTE THAT THIS IS NOT MAKING A COPY,
            // ELEMENTS OF nameArray SHOULD NOT BE MODIFIED.
            nameArray[j++] = maf_mafLine_getSpecies(mlArray[r]);
        }
    }
    return nameArray;
}
uint64_t* cullPositionsByColumn(char **mat, uint64_t c, uint64_t *positions, bool *legitRows,
                                uint64_t numRows, uint64_t numLegitGaplessPositions) {
    // create an array of positive coordinate position values that excludes all sequences that contain gaps
    uint64_t *colPositions = (uint64_t*) st_malloc(sizeof(*positions) * numLegitGaplessPositions);
    uint64_t j = 0;
    for (uint64_t r = 0; r < numRows; ++r) {
        if (!legitRows[r]) {
            continue;
        }
        if (mat[r][c] != '-') {
            colPositions[j++] = positions[r];
        }
    }
    return colPositions;
}
void validateMafBlockSourceLengths(const char *filename, mafBlock_t *mb, stHash *sequenceLengthHash) {
    // make sure that the information contained in the maf matches the information
    // contained in the sequenceLengthHash data.
    mafLine_t *ml = maf_mafBlock_getHeadLine(mb);
    while (ml != NULL) {
        if (maf_mafLine_getType(ml) == 's') {
            uint64_t *len = (uint64_t*)stHash_search(sequenceLengthHash, maf_mafLine_getSpecies(ml));
            if (len != NULL) {
                if ((uint64_t)maf_mafLine_getSourceLength(ml) != *len) {
                    fprintf(stderr, "Error, conflicting source length information for sequence in maf. "
                            "%s source length was first %" PRIu64 " but on line %" PRIu64 " of %s the value "
                            "is %" PRIu64 ".\n",
                            maf_mafLine_getSpecies(ml), *len, maf_mafLine_getLineNumber(ml), filename,
                            (uint64_t)maf_mafLine_getSourceLength(ml));
                    exit(EXIT_FAILURE);
                }
            }
        }
        ml = maf_mafLine_getNext(ml);
    }
}
void walkBlockSamplingPairs(const char *filename, mafBlock_t *mb, stSortedSet *sampledPairs,
                            double acceptProbability, stSet *legitSequences,
                            uint64_t *chooseTwoArray, uint64_t *numPairs, stHash *sequenceLengthHash) {
    uint64_t numSeqs = maf_mafBlock_getNumberOfSequences(mb);
    uint64_t numLegitGaplessPositions; // number of legit gapless sequences in the given column
    if (numSeqs < 2) {
        return;
    }
    validateMafBlockSourceLengths(filename, mb, sequenceLengthHash);

    uint64_t seqFieldLength = maf_mafBlock_getSequenceFieldLength(mb);
    char **names = maf_mafBlock_getSpeciesArray(mb);
    char **mat = maf_mafBlock_getSequenceMatrix(mb, numSeqs, seqFieldLength);
    bool *legitRows = getLegitRows(names, numSeqs, legitSequences);
    uint64_t numLegit = sumBoolArray(legitRows, numSeqs);
    if (numLegit < 2) {
        return;
    }
    mafLine_t **mlArray = maf_mafBlock_getMafLineArray_seqOnly(mb);
    uint64_t *allPositions = maf_mafBlock_getPosCoordStartArray(mb);
    int *allStrandInts = maf_mafBlock_getStrandIntArray(mb);
    char **gaplessNameArray = NULL;
    uint64_t *gaplessPositions = NULL;
    // walk over each column in the block
    for (uint64_t c = 0; c < seqFieldLength; ++c) {
        numLegitGaplessPositions = countLegitGaplessPositions(mat, c, numSeqs, legitRows);
        // create arrays that contain *only* the valid (legit and non gap) sequences for this column
        gaplessNameArray = extractLegitGaplessNamesFromMlArrayByColumn(mat, c, mlArray, legitRows,
                                                                       numSeqs, numLegitGaplessPositions);
        gaplessPositions = cullPositionsByColumn(mat, c, allPositions, legitRows,
                                                 numSeqs, numLegitGaplessPositions);
        samplePairsFromColumn(acceptProbability, sampledPairs, numLegitGaplessPositions, chooseTwoArray,
                              gaplessNameArray, gaplessPositions);
        updatePositions(mat, c, allPositions, allStrandInts, numSeqs);
        // double check:
        if (numLegitGaplessPositions < kChooseTwoCacheLength) {
            *numPairs += chooseTwoArray[numLegitGaplessPositions];
        } else {
            *numPairs += chooseTwo(numLegitGaplessPositions);
        }
    }
    // clean up
    free(mlArray);
    free(allPositions);
    free(allStrandInts);
    for (uint64_t i = 0; i < numSeqs; ++i) {
         free(names[i]);
    }
    free(names);
    maf_mafBlock_destroySequenceMatrix(mat, numSeqs);
    free(legitRows);
}
void samplePairsFromMaf(const char *filename, stSortedSet *pairs, double acceptProbability,
                        stSet *legitSequences, uint64_t *numPairs, stHash *sequenceLengthHash) {
    mafFileApi_t *mfa = maf_newMfa(filename, "r");
    mafBlock_t *mb = NULL;
    uint64_t *chooseTwoArray = buildChooseTwoArray();
    while ((mb = maf_readBlock(mfa)) != NULL) {
        walkBlockSamplingPairs(filename, mb, pairs, acceptProbability, legitSequences, chooseTwoArray,
                               numPairs, sequenceLengthHash);
        maf_destroyMafBlockList(mb);
    }
    // clean up
    free(chooseTwoArray);
    maf_destroyMfa(mfa);
}
void countPairs(APair *pair, stHash *intervalsHash, int64_t *counter,
                stSortedSet *legitPairs, void *a, uint64_t near) {
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
void samplePairs(APair *thisPair, stHash *intervalsHash, stSortedSet *pairs,
                 double *acceptProbability, stHash *legitPairs, uint64_t near) {
    /*
     * Adds *thisPair to *pairs with a given probability.
     */
    if (stHash_search(legitPairs, thisPair->seq1) != NULL)
        if (stHash_search(legitPairs, thisPair->seq2) != NULL)
            if (st_random() <= *acceptProbability)
                stSortedSet_insert(pairs, aPair_copyConstruct(thisPair));
}
bool inInterval(stHash *intervalsHash, char *seq, uint64_t pos) {
    /*
     * check to see if the sequence and position are within the intervals hash.
     */
    stSortedSet *intervals = stHash_search(intervalsHash, seq);
    if (intervals == NULL) {
        return false;
    }
    stIntTuple *i = stIntTuple_construct2( pos, INT64_MAX);
    stIntTuple *j = stSortedSet_searchLessThanOrEqual(intervals, i);
    stIntTuple_destruct(i);
    if (j == NULL) {
        return false;
    }
    assert(stIntTuple_length(j) == 2);
    return stIntTuple_get(j, 0) <= pos && pos < stIntTuple_get(j, 1);
}
uint64_t findLowerBound(uint64_t pos, uint64_t near) {
    // since we have unsigned values we must be careful about subtracting
    // the "near" value willy-nilly.
    if (near >= pos) {
        return 0;
    } else {
        return pos - near;
    }
}
void recordNearPair(APair *thisPair, stSortedSet *sampledPairs, uint64_t near, stSet *positivePairs) {
    /* given thisPair, if thisPair is in the table `sampledPairs' then record it in the positivePairs set.
     * if the `near' option is set, do this not only for thisPair but for all pairs within +- `near'.
     * if near = 0 this will just look at thisPair->pos1 and thisPair->pos2 and record those values.
     */
    APair *aPair;
    uint64_t i = thisPair->pos1;
    // Try modifying position 1
    for (thisPair->pos1 = findLowerBound(thisPair->pos1, near); thisPair->pos1 < i + near + 1; thisPair->pos1++) {
        if ((aPair = stSortedSet_search(sampledPairs, thisPair)) != NULL) {
            if (strcmp(thisPair->seq1, aPair->seq1) == 0 && strcmp(thisPair->seq2, aPair->seq2) == 0 &&
                thisPair->pos1 == aPair->pos1 && thisPair->pos2 == aPair->pos2) {
                // printf("positivePair1: (%s %u, %s %u):(%s %u, %s %u)\n",
                //        thisPair->seq1, thisPair->pos1, thisPair->seq2, thisPair->pos2,
                //        aPair->seq1, aPair->pos1, aPair->seq2, aPair->pos2);
                stSet_insert(positivePairs, aPair);
            }
        }
    }
    thisPair->pos1 = i; // reset pos 1
    // Try modifying position 2
    i = thisPair->pos2;
    for (thisPair->pos2 = findLowerBound(thisPair->pos2, near); thisPair->pos2 < i + near + 1; thisPair->pos2++) {
        if ((aPair = stSortedSet_search(sampledPairs, thisPair)) != NULL) {
            if (strcmp(thisPair->seq1, aPair->seq1) == 0 && strcmp(thisPair->seq2, aPair->seq2) == 0 &&
                thisPair->pos1 == aPair->pos1 && thisPair->pos2 == aPair->pos2) {
                // printf("positivePair2: (%s %u, %s %u):(%s %u, %s %u)\n",
                //        thisPair->seq1, thisPair->pos1, thisPair->seq2, thisPair->pos2,
                //        aPair->seq1, aPair->pos1, aPair->seq2, aPair->pos2);
                stSet_insert(positivePairs, aPair);
            }
        }
    }
    thisPair->pos2 = i; // reset pos 2
}
stHash* constructPositionHash(char **mat, uint64_t c, char **names, uint64_t numSeqs,
                              uint64_t *allPositions, bool *legitRows) {
    stHash *posHash = stHash_construct3(aPositionKey, aPositionEqualKey, aPosition_destruct, free);
    APosition *pos = NULL;
    // printf("constructing position hash on col %"PRIu64", numseqs: %"PRIu64"\n", c, numSeqs);
    for (uint64_t r = 0; r < numSeqs; ++r) {
        if (!legitRows[r]) {
            // printf("row %"PRIu64" not legit.\n", r);
            continue;
        }
        if (mat[r][c] == '-') {
            // printf("row %"PRIu64" col %"PRIu64" is gap.\n", r, c);
            continue;
        }
        pos = aPosition_construct(stString_copy(names[r]), allPositions[r]);
        if (stHash_search(posHash, pos) == NULL) {
            // printf("adding position to posHash (%s %u %c)\n", names[r], allPositions[r], mat[r][c]);
            stHash_insert(posHash, pos, stString_copy(""));
        } else {
            aPosition_destruct(pos);
        }
    }
    return posHash;
}
bool pairMemberInPositionHash(stHash *positionHash, APair *thisPair) {
    // we know seq1 and pos1 are in the position hash. what matters is if seq2 and pos2 are in.
    APosition *p = aPosition_init();
    p->name = stString_copy(thisPair->seq2);
    p->pos = thisPair->pos2;
    // printf("checking to see if pairmemberinpositionhash, (%s %u)... ", p->name, p->pos);
    if (stHash_search(positionHash, p) != NULL) {
        // printf("YES :D\n");
        aPosition_destruct(p);
        return true;
    } else {
        // printf("NO :(\n");
    }
    aPosition_destruct(p);
    return false;
}
void printSortedSet(stSortedSet *pairs) {
    stSortedSetIterator *sit = NULL;
    sit = stSortedSet_getIterator(pairs);
    APair *pair;
    printf("printSortedSet()\n");
    while ((pair = stSortedSet_getNext(sit)) != NULL) {
        printf("pair: (%s %" PRIu64", %s %" PRIu64 ")\n",
               pair->seq1, pair->pos1, pair->seq2, pair->pos2);
    }
}
void printHash(stHash *hash) {
    stHashIterator *hit = NULL;
    hit = stHash_getIterator(hash);
    APosition *key = NULL;
    printf("Position hash: ");
    while ((key = stHash_getNext(hit)) != NULL) {
        printf("(%s %" PRIi64 "), ", key->name, key->pos);
    }
    printf("\n");
    stHash_destructIterator(hit);
}
void testHomologyOnColumn(char **mat, uint64_t c, uint64_t numSeqs, bool *legitRows, char **names,
                          stSortedSet *sampledPairs, stSet *positivePairs, mafLine_t **mlArray,
                          uint64_t *allPositions, stHash *intervalsHash, uint64_t near) {
    /* For a given column,
       1) hash all the positions in the column
       2) For each position in the hash:
       ..a) Iterate over pairs involving the position in the sortedSet,
       .....i) check if other aligned positions are in the hash of step 1.
     */
    APair *thisPair = NULL;
    APair *thatPair = NULL;
    APair *otherPair = NULL;
    APosition *key = NULL;
    stHash *positionHash = NULL;
    stHashIterator *hit = NULL;
    stSortedSetIterator *sit = NULL;
    // 1.
    positionHash = constructPositionHash(mat, c, names, numSeqs, allPositions, legitRows);
    hit = stHash_getIterator(positionHash);
    // 2.
    while ((key = stHash_getNext(hit)) != NULL) {
        assert(thisPair == NULL);
        thisPair = aPair_init();
        thisPair->seq1 = stString_copy(key->name);
        thisPair->pos1 = key->pos;
        if ((thatPair = stSortedSet_searchGreaterThanOrEqual(sampledPairs, thisPair)) != NULL) {
            if (strcmp(thisPair->seq1, thatPair->seq1) != 0){
                aPair_destruct(thisPair);
                thisPair = NULL;
                continue;
            }
            if (thisPair->pos1 != thatPair->pos1) {
                aPair_destruct(thisPair);
                thisPair = NULL;
                continue;
            }
            sit = stSortedSet_getIteratorFrom(sampledPairs, thatPair);
            // 2a.
            while ((otherPair = stSortedSet_getNext(sit)) != NULL) {
                if ((strcmp(thatPair->seq1, thatPair->seq1) != 0) ||
                    (!closeEnough(thatPair->pos1, otherPair->pos1, near))) {
                    // bail out on iteration once we've overstepped the range of interest
                    aPair_destruct(thisPair);
                    thisPair = NULL;
                    break;
                }
                thisPair->seq2 = stString_copy(otherPair->seq2);
                thisPair->pos2 = otherPair->pos2;
                // 2ai.
                if (pairMemberInPositionHash(positionHash, thisPair)) {
                    recordNearPair(thisPair, sampledPairs, near, positivePairs);
                }
                free(thisPair->seq2);
                thisPair->seq2 = NULL;
            }
            stSortedSet_destructIterator(sit);
        }
        aPair_destruct(thisPair);
        thisPair = NULL;
    }
    assert(thisPair == NULL);
    stHash_destructIterator(hit);
    stHash_destruct(positionHash);
}
void printAllPositions(uint64_t *allPositions, mafBlock_t *mb) {
    printf("allPositions: [");
    for (uint64_t i = 0; i < maf_mafBlock_getNumberOfSequences(mb); ++i) {
        printf("%"PRIu64", ", allPositions[i]);
    }
    printf("]\n");
}
void printAllStrandInts(int *allStrandInts, mafBlock_t *mb) {
    printf("allStrandInts: [");
    for (uint64_t i = 0; i < maf_mafBlock_getNumberOfSequences(mb); ++i) {
        printf("%+i, ", allStrandInts[i]);
    }
    printf("]\n");
}
void walkBlockTestingHomology(mafBlock_t *mb, stSortedSet *sampledPairs, stSet *positivePairs,
                              stSet *legitSequences, stHash *intervalsHash, uint64_t near) {
    uint64_t numSeqs = maf_mafBlock_getNumberOfSequences(mb);
    if (numSeqs < 2) {
        return;
    }
    uint64_t seqFieldLength = maf_mafBlock_getSequenceFieldLength(mb);
    char **names = maf_mafBlock_getSpeciesArray(mb);
    char **mat = maf_mafBlock_getSequenceMatrix(mb, numSeqs, seqFieldLength);
    bool *legitRows = getLegitRows(names, numSeqs, legitSequences);
    uint64_t numLegit = sumBoolArray(legitRows, numSeqs);
    if (numLegit < 2) {
        return;
    }
    mafLine_t **mlArray = createMafLineArray(mb, numLegit, legitRows);
    uint64_t *allPositions = maf_mafBlock_getPosCoordStartArray(mb);
    int *allStrandInts = maf_mafBlock_getStrandIntArray(mb);
    for (uint64_t c = 0; c < seqFieldLength; ++c) {
        testHomologyOnColumn(mat, c, numSeqs, legitRows, names, sampledPairs, positivePairs,
                             mlArray, allPositions, intervalsHash, near);
        updatePositions(mat, c, allPositions, allStrandInts, numSeqs);
    }
    // clean up
    free(mlArray);
    free(allPositions);
    free(allStrandInts);
    for (uint64_t i = 0; i < numSeqs; ++i) {
         free(names[i]);
    }
    free(names);
    maf_mafBlock_destroySequenceMatrix(mat, numSeqs);
    free(legitRows);
}
void performHomologyTests(const char *filename, stSortedSet *sampledPairs, stSet *positivePairs,
                          stSet *legitSequences, stHash *intervalsHash, uint64_t near) {
    mafFileApi_t *mfa = maf_newMfa(filename, "r");
    mafBlock_t *mb = NULL;
    while ((mb = maf_readBlock(mfa)) != NULL) {
        walkBlockTestingHomology(mb, sampledPairs, positivePairs, legitSequences, intervalsHash, near);
        maf_destroyMafBlockList(mb);
    }
    // clean up
    maf_destroyMfa(mfa);
}
void homologyTests1(APair *thisPair, stHash *intervalsHash, stSortedSet *pairs,
                    stSet *positivePairs, stSet *legitPairs, int64_t near) {
    /*
     * If both members of *thisPair are in the intersection of maf1 and maf2,
     * and *thisPair is in the set *pairs then adds to the result pair a positive result.
     */
    if ((stSet_search(legitPairs, thisPair->seq1) != NULL)
        && (stSet_search(legitPairs, thisPair->seq2) != NULL)) {
        recordNearPair(thisPair, pairs, near, positivePairs);
    }
}
bool positionIsInWiggleRegion(WiggleContainer *wc, uint64_t *refPos) {
    // verify that refPos, a positive coordinate, is within the wiggle
    // container's region of interest
    if (wc == NULL) {
        return false;
    }
    if (refPos == NULL) {
        return false;
    }
    if (*refPos < wc->refStart) {
        return false;
    }
    if (*refPos > wc->refStart + wc->refLength) {
        return false;
    }
    return true;
}
void enumerateHomologyResults(stSortedSet *sampledPairs, stSortedSet *resultPairs, stHash *intervalsHash,
                              stSet *positivePairs, stHash *wigglePairHash, bool isAtoB,
                              uint64_t wiggleBinLength) {
    /*
     * For every pair in 'sampledPairs', add 1 to the total number of homology tests for the sequence-pair
     * (the ResultPair).
     */
    static stSortedSetIterator *sit;
    sit = stSortedSet_getIterator(sampledPairs);
    APair *pair = NULL;
    ResultPair *thisResultPair = NULL;
    WiggleContainer *wc = NULL;
    uint64_t *refPos = NULL;
    uint64_t localPos = 0; // local offset within the region of interest (0 is wc->refStart)
    char wigKey[kMaxStringLength];
    wigKey[0] = '\0';
    while ((pair = stSortedSet_getNext(sit)) != NULL) {
        if ((thisResultPair = stSortedSet_search(resultPairs, pair)) == NULL) {
            // the stSortedSet resultPairs is searched only based on sequence names.
            thisResultPair = resultPair_construct(pair->seq1, pair->seq2);
            stSortedSet_insert(resultPairs, thisResultPair);
        }
        sprintf(wigKey, "%s-%s", pair->seq1, pair->seq2);
        if ((wc = stHash_search(wigglePairHash, wigKey)) == NULL) {
            // seq1 is not the ref
            sprintf(wigKey, "%s-%s", pair->seq2, pair->seq1);
            if ((wc = stHash_search(wigglePairHash, wigKey)) == NULL) {
                // seq2 is also not the ref
                wigKey[0] = '\0';
            } else {
                refPos = buildUInt64(pair->pos2);
            }
        } else {
            refPos = buildUInt64(pair->pos1);
        }
        bool foundPair = stSet_search(positivePairs, pair) != NULL;
        if (inInterval(intervalsHash, pair->seq1, pair->pos1)) {
            if (inInterval(intervalsHash, pair->seq2, pair->pos2)) {
                ++(thisResultPair->totalBoth);
                if (foundPair) {
                    ++(thisResultPair->inBoth);
                }
            } else {
                ++(thisResultPair->totalA);
                if (foundPair) {
                    ++(thisResultPair->inA);
                }
            }
        } else {
            if (inInterval(intervalsHash, pair->seq2, pair->pos2)) {
                ++(thisResultPair->totalB);
                if (foundPair) {
                    ++(thisResultPair->inB);
                }
            } else {
                ++(thisResultPair->totalNeither);
                if (foundPair) {
                    ++(thisResultPair->inNeither);
                }
            }
        }

        ++(thisResultPair->total);
        // put results in wiggle pairs
        if (positionIsInWiggleRegion(wc, refPos)) {
            localPos = *refPos - wc->refStart;
            if (isAtoB) {
                ++(wc->absentAtoB[(int)floor(localPos / wiggleBinLength)]);
            } else {
                ++(wc->absentBtoA[(int)floor(localPos / wiggleBinLength)]);
            }
        }
        if (foundPair) {
            ++(thisResultPair->inAll);
            if (positionIsInWiggleRegion(wc, refPos)) {
                localPos = *refPos - wc->refStart;
                if (isAtoB) {
                    --(wc->absentAtoB[(int)floor(localPos / wiggleBinLength)]);
                    ++(wc->presentAtoB[(int)floor(localPos / wiggleBinLength)]);
                } else {
                    --(wc->absentBtoA[(int)floor(localPos / wiggleBinLength)]);
                    ++(wc->presentBtoA[(int)floor(localPos / wiggleBinLength)]);
                }
            }
        } else {
           if (g_isVerboseFailures){
              fprintf(stderr, "sampled pair not present in comparison: (%s, %" PRIu64 "):(%s, %" PRIu64 ")\n",
                      pair->seq1, pair->pos1, pair->seq2, pair->pos2);
           }
        }
        if (refPos != NULL) {
            free(refPos);
            refPos = NULL;
        }
    }
    stSortedSet_destructIterator(sit);
}
stSortedSet *compareMAFs_AB(const char *mafFileA, const char *mafFileB, uint64_t *numberOfPairs,
                            stSet *legitSequences, stHash *intervalsHash, stHash *wigglePairHash,
                            bool isAtoB, Options *options, stHash *sequenceLengthHash) {
    // count the number of pairs in mafFileA
    if (*numberOfPairs == 0) {
        // can be manually set via the command line
        *numberOfPairs = countPairsInMaf(mafFileA, legitSequences);
    }
    if (*numberOfPairs == 0) {
        return stSortedSet_construct3((int(*)(const void *, const void *)) aPair_cmpFunction_seqsOnly, (void(*)(void *)) aPair_destruct);
    }
    double acceptProbability = ((double) options->numberOfSamples) / (double) *numberOfPairs;
    stSortedSet *pairs = stSortedSet_construct3((int(*)(const void *, const void *)) aPair_cmpFunction, (void(*)(void *)) aPair_destruct);
    // sample pairs from mafFileA
    uint64_t verifiedNumberOfPairs = 0;
    samplePairsFromMaf(mafFileA, pairs, acceptProbability, legitSequences, &verifiedNumberOfPairs,
                       sequenceLengthHash);
    if (verifiedNumberOfPairs != *numberOfPairs) {
        fprintf(stderr, "Error, differing numberOfPairs values, %"PRIu64" != %"PRIu64"\n",
                verifiedNumberOfPairs, *numberOfPairs);
        exit(EXIT_FAILURE);
    }
    // perform homology tests on mafFileB using sampled pairs from mafFileA
    stSet *positivePairs = stSet_construct(); // comparison by pointer
    performHomologyTests(mafFileB, pairs, positivePairs, legitSequences, intervalsHash, options->near);
    stSortedSet *resultPairs = stSortedSet_construct3((int(*)(const void *, const void *)) aPair_cmpFunction_seqsOnly, (void(*)(void *)) aPair_destruct);
    enumerateHomologyResults(pairs, resultPairs, intervalsHash, positivePairs, wigglePairHash, isAtoB,
                             options->wiggleBinLength);
    // clean up
    stSortedSet_destruct(pairs);
    stSet_destruct(positivePairs);
    return resultPairs;
}
ResultPair *aggregateResult(void *(*getNextPair)(void *, void *), stSortedSet *set, void *seqName,
                            const char *name1, const char *name2) {
    /* loop through all ResultPairs available via the getNextPair() iterator and aggregate their
     * results into a single ResultPair, return this struct.
     */
    ResultPair *resultPair;
    ResultPair *resultPair2 = resultPair_construct(name1, name2);
    stSortedSetIterator *iterator = stSortedSet_getIterator(set);
    while ((resultPair = getNextPair(iterator, seqName)) != NULL) {
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
    stSortedSet_destructIterator(iterator);
    return resultPair2;
}
void* addReferencesAndDups_getDups(void *iterator, void *seqName) {
    ResultPair *resultPair;
    while ((resultPair = stSortedSet_getNext(iterator)) != NULL) {
        if (strcmp(resultPair->seq1, resultPair->seq2) == 0) {
            break;
        }
    }
    return resultPair;
}
void* addReferencesAndDups_getReferences(void *iterator, void *seqName) {
    ResultPair *resultPair;
    while ((resultPair = stSortedSet_getNext(iterator)) != NULL) {
        if (strcmp(resultPair->seq1, seqName) == 0 || strcmp(resultPair->seq2, seqName) == 0) {
            break;
        }
    }
    return resultPair;
}
void addReferencesAndDups(stSortedSet *results_AB, stSet *legitSequences) {
    /*
     * Adds tags for all against each species.
     */
    stSetIterator *speciesIt = stSet_getIterator(legitSequences);
    char *species;
    stList *list = stList_construct();
    // add duplicates
    stList_append(list, aggregateResult(addReferencesAndDups_getDups, results_AB, NULL, "self", "self"));
    // add references
    while ((species = stSet_getNext(speciesIt)) != NULL) {
        stList_append(list, aggregateResult(addReferencesAndDups_getReferences, results_AB,
                                            species, species, "aggregate"));
    }
    stSet_destructIterator(speciesIt);
    for(int64_t i = 0; i < stList_length(list); i++) {
        stSortedSet_insert(results_AB, stList_get(list, i));
    }
    stList_destruct(list);
}
void* reportResults_fn(void *iterator, void *seqName) {
   return stSortedSet_getNext(iterator);
}
void findentprintf(FILE *fp, unsigned indent, char const *fmt, ...) {
    va_list args;
    va_start(args, fmt);
    char str[kMaxStringLength];
    int n = vsprintf(str, fmt, args);
        if (n >= kMaxStringLength) {
            fprintf(stderr, "Error, failure in findentprintf, (n = %d) >= %d\n", n, kMaxStringLength);
            exit(EXIT_FAILURE);
        }
    for (unsigned i = 0; i < indent; ++i) {
        fprintf(fp, "\t");
    }
    fprintf(fp, "%s", str);
    va_end(args);
}
void reportResultsForWigglesArrays(FILE *fileHandle, WiggleContainer *wc) {
    // present A to B
    findentprintf(fileHandle, 2, "<presentMaf1ToMaf2 description=\"Binned data, counts. "
                  "Values are the number of samples drawn from maf1 that are present in maf2 and "
                  "whose position on the reference sequence coincides with the given bin.\">");
    for (uint64_t i = 0; i < wc->numBins; ++i) {
        fprintf(fileHandle, "%" PRIu64 "%s", wc->presentAtoB[i], (i == wc->numBins - 1) ? "" : ",");
    }
    fprintf(fileHandle, "</presentMaf1ToMaf2>\n");
    // present B to A
    findentprintf(fileHandle, 2, "<presentMaf2ToMaf1 description=\"Binned data, counts. "
                  "Values are the number of samples drawn from maf2 that are present in maf1 and "
                  "whose position on the reference sequence coincides with the given bin.\">");
    for (uint64_t i = 0; i < wc->numBins; ++i) {
        fprintf(fileHandle, "%" PRIu64 "%s", wc->presentBtoA[i], (i == wc->numBins - 1) ? "" : ",");
    }
    fprintf(fileHandle, "</presentMaf2ToMaf1>\n");
    // absent A to B
    findentprintf(fileHandle, 2, "<absentMaf1ToMaf2 description=\"Binned data, counts. "
                  "Values are the number of samples drawn from maf1 that are absent in maf2 and "
                  "whose position on the reference sequence coincides with the given bin.\">");
    for (uint64_t i = 0; i < wc->numBins; ++i) {
        fprintf(fileHandle, "%" PRIu64 "%s", wc->absentAtoB[i], (i == wc->numBins - 1) ? "" : ",");
    }
    fprintf(fileHandle, "</absentMaf1ToMaf2>\n");
    // absent B to A
    findentprintf(fileHandle, 2, "<absentMaf2ToMaf1 description=\"Binned data, counts. "
                  "Values are the number of samples drawn from maf2 that are absent in maf1 and "
                  "whose position on the reference sequence coincides with the given bin.\">");
    for (uint64_t i = 0; i < wc->numBins; ++i) {
        fprintf(fileHandle, "%" PRIu64 "%s", wc->absentBtoA[i], (i == wc->numBins - 1) ? "" : ",");
    }
    fprintf(fileHandle, "</absentMaf2ToMaf1>\n");
}
void reportResultsForWiggles(stHash *wigglePairHash, FILE *fileHandle) {
    //
    stHashIterator *hit = stHash_getIterator(wigglePairHash);
    WiggleContainer *wc = NULL;
    char *key = NULL;
    findentprintf(fileHandle, 1, "<wigglePairs>\n");
    while ((key = stHash_getNext(hit)) != NULL) {
        wc = stHash_search(wigglePairHash, key);
        findentprintf(fileHandle, 2, "<wigglePair reference=\"%s\" partner=\"%s\" "
                      "referenceStart=\"%" PRIu64 "\" referenceLength=\"%" PRIu64 "\" "
                      "numberOfBins=\"%" PRIu64 "\" "
                      "binLength=\"%" PRIu64 "\" >\n",
                      wc->ref, wc->partner, wc->refStart, wc->refLength, wc->numBins, wc->binLength);

        reportResultsForWigglesArrays(fileHandle, wc);
        findentprintf(fileHandle, 2, "</wigglePair>\n");
    }
    findentprintf(fileHandle, 1, "</wigglePairs>\n");
    stHash_destructIterator(hit);
}
void reportResult(const char *tagName, double total, double totalTrue, FILE *fileHandle, unsigned tabLevel) {
    assert(total >= totalTrue);
    findentprintf(fileHandle, tabLevel, "<%s totalTests=\"%" PRIu64 "\" totalTrue=\"%" PRIu64 "\" "
                  "totalFalse=\"%" PRIu64 "\" average=\"%f\"/>\n",
                  tagName, (uint64_t) total, (uint64_t) totalTrue,
                  (uint64_t) (total - totalTrue), total == 0 ? 0.0 : totalTrue / total);
}
void reportResults(stSortedSet *results_AB, const char *mafFileA, const char *mafFileB, FILE *fileHandle,
                   uint64_t near, stSet *legitSequences, const char *bedFiles) {
    /*
     * Report results in an XML formatted document.
     */
    stSortedSetIterator *iterator = NULL;
    ResultPair *resultPair;
    unsigned tabLevel = 1;
    ResultPair *aggregateResults = aggregateResult(reportResults_fn, results_AB, NULL, "", "");
    addReferencesAndDups(results_AB, legitSequences);
    findentprintf(fileHandle, tabLevel++, "<homologyTests fileA=\"%s\" fileB=\"%s\">\n", mafFileA, mafFileB);
    findentprintf(fileHandle, tabLevel++, "<aggregateResults>\n");
    reportResult("all", aggregateResults->total, aggregateResults->inAll, fileHandle, tabLevel);
    if (bedFiles != NULL){
        reportResult("both", aggregateResults->totalBoth, aggregateResults->inBoth, fileHandle, tabLevel);
        reportResult("A", aggregateResults->totalA, aggregateResults->inA, fileHandle, tabLevel);
        reportResult("B", aggregateResults->totalB, aggregateResults->inB, fileHandle, tabLevel);
        reportResult("neither", aggregateResults->totalNeither, aggregateResults->inNeither, fileHandle, tabLevel);
    }
    findentprintf(fileHandle, --tabLevel, "</aggregateResults>\n");
    findentprintf(fileHandle, tabLevel++, "<homologyPairTests>\n");
    iterator = stSortedSet_getIterator(results_AB);
    while ((resultPair = stSortedSet_getPrevious(iterator)) != NULL) {
        findentprintf(fileHandle, tabLevel++, "<homologyTest sequenceA=\"%s\" sequenceB=\"%s\">\n",
                      resultPair->seq1, resultPair->seq2);
        findentprintf(fileHandle, tabLevel++, "<aggregateResults>\n");
        reportResult("all", resultPair->total, resultPair->inAll, fileHandle, tabLevel);
        if (bedFiles != NULL){
            reportResult("both", resultPair->totalBoth, resultPair->inBoth, fileHandle, tabLevel);
            reportResult("A", resultPair->totalA, resultPair->inA, fileHandle, tabLevel);
            reportResult("B", resultPair->totalB, resultPair->inB, fileHandle, tabLevel);
            reportResult("neither", resultPair->totalNeither, resultPair->inNeither, fileHandle, tabLevel);
        }

        findentprintf(fileHandle, --tabLevel, "</aggregateResults>\n");
        findentprintf(fileHandle, tabLevel++, "<singleHomologyTests>\n");
        findentprintf(fileHandle, tabLevel++, "<singleHomologyTest sequenceA=\"%s\" sequenceB=\"%s\">\n",
                      resultPair->seq1, resultPair->seq2);
        findentprintf(fileHandle, tabLevel++, "<aggregateResults>\n");
        reportResult("all", resultPair->total, resultPair->inAll, fileHandle, tabLevel);
        if (bedFiles != NULL){
            reportResult("both", resultPair->totalBoth, resultPair->inBoth, fileHandle, tabLevel);
            reportResult("A", resultPair->totalA, resultPair->inA, fileHandle, tabLevel);
            reportResult("B", resultPair->totalB, resultPair->inB, fileHandle, tabLevel);
            reportResult("neither", resultPair->totalNeither, resultPair->inNeither, fileHandle, tabLevel);
        }
        findentprintf(fileHandle, --tabLevel, "</aggregateResults>\n");
        findentprintf(fileHandle, --tabLevel, "</singleHomologyTest>\n");
        findentprintf(fileHandle, --tabLevel, "</singleHomologyTests>\n");
        findentprintf(fileHandle, --tabLevel, "</homologyTest>\n");
    }
    findentprintf(fileHandle, --tabLevel, "</homologyPairTests>\n");
    findentprintf(fileHandle, --tabLevel, "</homologyTests>\n");
    resultPair_destruct(aggregateResults);
    stSortedSet_destructIterator(iterator);
    assert(tabLevel == 1);
    return;
}
void populateNames(const char *filename, stSet *set, stHash *sequenceLengthHash) {
    /*
     * populates a set with the names of sequences from a MAF file.
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
                if (stHash_search(sequenceLengthHash, name) == NULL) {
                    stHash_insert(sequenceLengthHash, stString_copy(name),
                                  buildInt64(maf_mafLine_getSourceLength(ml)));
                } else {
                    if (*(int64_t*)stHash_search(sequenceLengthHash, name) != maf_mafLine_getSourceLength(ml)) {
                        fprintf(stderr, "Inconsistency detected in a maf. Previous source length for sequence "
                                "%s was %" PRIu64 " but changed to %" PRIu64 " on line %" PRIu64 "\n",
                                name, maf_mafLine_getSourceLength(ml),
                                *(int64_t*)stHash_search(sequenceLengthHash, name),
                                maf_mafLine_getLineNumber(ml));
                        exit(EXIT_FAILURE);
                    }
                }
                if (stSet_search(set, name) == NULL) {
                    stSet_insert(set, stString_copy(name));
                }
            }
            ml = maf_mafLine_getNext(ml);
        }
        maf_destroyMafBlockList(mb);
    }
    // clean up
    maf_destroyMfa(mfa);
}
void writeXMLHeader(FILE *fileHandle){
    fprintf(fileHandle, "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\" ?>\n");
    return;
}
unsigned countChars(char *s, char c) {
    unsigned v = 0;
    if (s == 0) {
        return v;
    }
    unsigned l = strlen(s);
    for (unsigned i = 0; i < l; ++i) {
        if (s[i] == c) {
            ++v;
        }
    }
    return v;
}
bool patternMatches(char *a, char *b) {
    // special matching function, recognizes * as a wildcard character if at the end
    // of the string, or not at all.
    if (a[strlen(a) - 1] == '*') {
        return strncmp(a, b, strlen(a) - 1) == 0;
    } else {
        return strcmp(a, b) == 0;
    }
}
void buildWigglePairHash(stHash *sequenceLengthHash, stList *wigglePairPatternList, stHash *wigglePairHash,
                         uint64_t wiggleBinLength, uint64_t wiggleRegionStart, uint64_t wiggleRegionStop) {
    stHashIterator *hit1 = NULL;
    stHashIterator *hit2 = NULL;
    char *key1 = NULL;
    char *key2 = NULL;
    char *name = NULL;
    WiggleContainer *wc = NULL;
    uint64_t regionLength = 0;
    for (int64_t i = 0; i < stList_length(wigglePairPatternList) - 1; i = i + 2) {
        // for every wiggle pair
        hit1 = stHash_getIterator(sequenceLengthHash);
        while ((key1 = stHash_getNext(hit1)) != NULL) {
            if (patternMatches(stList_get(wigglePairPatternList, i), key1)) {
                hit2 = stHash_getIterator(sequenceLengthHash);
                while ((key2 = stHash_getNext(hit2)) != NULL) {
                    if (patternMatches(stList_get(wigglePairPatternList, i + 1), key2)) {
                        if (wiggleRegionStop == 0) {
                            // if the option is not set, use the source length field
                            regionLength = *(uint64_t*)stHash_search(sequenceLengthHash, key1);
                        } else {
                            regionLength = 1 + wiggleRegionStop - wiggleRegionStart;
                        }
                        name = st_malloc(kMaxStringLength);
                        sprintf(name, "%s-%s", key1, key2);
                        wc = wiggleContainer_construct(key1, key2,
                                                       wiggleRegionStart,
                                                       regionLength,
                                                       wiggleBinLength);
                        stHash_insert(wigglePairHash, stString_copy(name), wc);
                        free(name);
                    }
                }
                stHash_destructIterator(hit2);
            }
        }
        stHash_destructIterator(hit1);
    }
}
void buildSeqNamesSet(Options *options, stSet *seqNamesSet, stHash *sequenceLengthHash) {
    uint64_t length = 0;
    if (options->legitSequences == NULL) {
        // read the input maf files and construct the set and hash from them
        stSet *seqNamesSet1 = stSet_construct3(stHash_stringKey, stHash_stringEqualKey, free);
        stSet *seqNamesSet2 = stSet_construct3(stHash_stringKey, stHash_stringEqualKey, free);
        populateNames(options->mafFile1, seqNamesSet1, sequenceLengthHash);
        populateNames(options->mafFile2, seqNamesSet2, sequenceLengthHash);
        stSet *seqNamesSetTmp = stSet_getIntersection(seqNamesSet1, seqNamesSet2);
        stSetIterator *sit = stSet_getIterator(seqNamesSetTmp);
        char *key = NULL;
        while ((key = stSet_getNext(sit)) != NULL) {
            // the intersection does not make copies and seqNamesSet1 and 2 need to be free'd
            stSet_insert(seqNamesSet, stString_copy(key));
        }
        stSet_destructIterator(sit);
        stSet_destruct(seqNamesSet1);
        stSet_destruct(seqNamesSet2);
        stSet_destruct(seqNamesSetTmp);
    } else {
        // trust the command line input from the user, and use those to build the set and hash
        char *spaceSep = stringReplace(options->legitSequences, ',', ' ');
        char *currentLocation = spaceSep;
        char *currentWord = NULL;
        while ((currentWord = stString_getNextWord(&currentLocation)) != NULL) {
            // separate out the sequence:length pairs by comma
            char *colonSep = stringReplace(currentWord, ':', ' ');
            char *currentLocation2 = colonSep;
            char *currentWord2 = NULL;
            int i = -1;
            char *seqName = NULL;
            while ((currentWord2 = stString_getNextWord(&currentLocation2)) != NULL) {
                // separate out the sequence string from the length value
                ++i;
                if (i % 2) {
                    // odd, sequence length
                    assert(1 == sscanf(currentWord2, "%" PRIu64, &length));
                    stHash_insert(sequenceLengthHash, stString_copy(seqName), buildInt64(length));
                    free(seqName);
                    seqName = NULL;
                } else {
                    // even, sequence string
                    assert(seqName == NULL);
                    seqName = stString_copy(currentWord2);
                    stSet_insert(seqNamesSet, stString_copy(seqName));
                }
                free(currentWord2);
            }
            free(colonSep);
            free(currentWord);
        }
        free(spaceSep);
    }
}
