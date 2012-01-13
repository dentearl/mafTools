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

#include "ComparatorAPI.h"
#include "cString.h"

int aSolo_cmpFunction(const void *a, const void *b);
void trioDecoder_destruct(TrioDecoder *decoder);
int32_t calcTrioState(TrioDecoder *decoder, int32_t spA, int32_t spB, int32_t spC);

void aPair_fillOut(APair *aPair, char *seq1, char *seq2, int32_t pos1, 
                   int32_t pos2, int32_t origPos1, int32_t origPos2) {
    int32_t i = strcmp(seq1, seq2);
    if (i > 0 || (i == 0 && pos1 > pos2)) { 
        // we construct the sequence so that seq1 < seq2 || seq1 == seq2 and pos1 < pos2
        aPair_fillOut(aPair, seq2, seq1, pos2, pos1, origPos2, origPos1);
    } else {
        aPair->seq1 = seq1;
        aPair->seq2 = seq2;
        aPair->pos1 = pos1;
        aPair->pos2 = pos2;
        aPair->origPos1 = origPos1;
        aPair->origPos2 = origPos2;
    }
}

void aTrio_fillOut(ATrio *aTrio, char *seq1, char *seq2, char *seq3, 
                   int32_t pos1, int32_t pos2, int32_t pos3, int32_t top) {
    ASolo structs[3];

    structs[0].name = seq1;
    structs[0].pos = pos1;
    structs[1].name = seq2;
    structs[1].pos = pos2;
    structs[2].name = seq3;
    structs[2].pos = pos3;

    qsort(structs, 3, sizeof(ASolo), aSolo_cmpFunction);

    aTrio->seq1 = structs[0].name;
    aTrio->pos1 = structs[0].pos;
    aTrio->seq2 = structs[1].name;
    aTrio->pos2 = structs[1].pos;
    aTrio->seq3 = structs[2].name;
    aTrio->pos3 = structs[2].pos;

    aTrio->top = top;
}

APair *aPair_construct(const char *seq1, const char *seq2, int32_t pos1, 
                       int32_t pos2, int32_t origPos1, int32_t origPos2) {
    APair *aPair = st_malloc(sizeof(APair));
    aPair_fillOut(aPair, (char *) seq1, (char *) seq2, pos1, pos2, origPos1, origPos2);
    aPair->seq1 = stString_copy(aPair->seq1);
    aPair->seq2 = stString_copy(aPair->seq2);
    return aPair;
}

ResultPair *resultPair_construct(const char *seq1, const char *seq2) {
    APair *aPair = aPair_construct(seq1, seq2, 0, 0, 0, 0);
    ResultPair *resultPair = st_calloc(1, sizeof(ResultPair));
    resultPair->aPair = *aPair;
    free(aPair);
    return resultPair;
}

ATrio *aTrio_construct(const char *seq1, const char *seq2, const char *seq3, 
                       int32_t pos1, int32_t pos2, int32_t pos3, int32_t top) {
    int32_t i;
    ATrio *aTrio = st_malloc(sizeof(ATrio));
    aTrio_fillOut(aTrio, (char *) seq1, (char *) seq2, (char *) seq3, pos1, pos2, pos3, top);
    aTrio->seq1 = stString_copy(aTrio->seq1);
    aTrio->seq2 = stString_copy(aTrio->seq2);
    aTrio->seq3 = stString_copy(aTrio->seq3);
    for (i = 0; i < 10; i++) {
        aTrio->topMat[i] = 0;
    }
    return aTrio;
}

APair *aPair_copyConstruct(APair *pair) {
   return aPair_construct(pair->seq1, pair->seq2, pair->pos1, pair->pos2, pair->origPos1, pair->origPos2);
}

ATrio *aTrio_copyConstruct(ATrio *trio) {
    return aTrio_construct(trio->seq1, trio->seq2, trio->seq3, trio->pos1, trio->pos2, trio->pos3, trio->top);
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

void aTrio_destruct(ATrio *trio, void *extraArgument) {
    assert(extraArgument == NULL); //because destruction takes place in avl_tree
    free(trio->seq1);
    free(trio->seq2);
    free(trio->seq3);
    free(trio);
}

int aSolo_cmpFunction(const void *a, const void *b) {
    ASolo *ia = (ASolo *) a;
    ASolo *ib = (ASolo *) b;

    int cmp = strcmp(ia->name, ib->name);

    if (cmp == 0) {
        return ia->pos - ib->pos;
    }
    return cmp;
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

int32_t aTrio_cmpFunction_seqsOnly(ATrio *aTrio1, ATrio *aTrio2, void *a) {
    /*
     * Compares only the sequence name arguments of the a trios.
     */
    assert(a == NULL);
    int32_t i = strcmp(aTrio1->seq1, aTrio2->seq1);
    if (i == 0) {
        i = strcmp(aTrio1->seq2, aTrio2->seq2);
        if (i == 0) {
            i = strcmp(aTrio1->seq3, aTrio2->seq3);
        }
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

int32_t aTrio_cmpFunction_seqAndPos(ATrio *aTrio1, ATrio *aTrio2, void *a) {
    /*
     * Compares first the sequence then the position arguments of the a trios.
     */
    int32_t i = aTrio_cmpFunction_seqsOnly(aTrio1, aTrio2, a);
    if (i == 0) {
        i = aTrio1->pos1 - aTrio2->pos1;
        if (i == 0) {
            i = aTrio1->pos2 - aTrio2->pos2;
            if (i == 0) {
                i = aTrio1->pos3 - aTrio2->pos3;
            }
        }
    }
    return i;
}

int32_t aTrio_cmpFunction(ATrio *aTrio1, ATrio *aTrio2, void *a) {
    /*
     * Compares first the sequence, then the position, then the top arguments of the a trios.
     */
    int32_t i = aTrio_cmpFunction_seqAndPos(aTrio1, aTrio2, a);
    //  if (i == 0) {
    //      i = aTrio1->top - aTrio2->top;
    //  }
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
                       aPair_fillOut(&aPair, seq1, seq2, pos1, pos2, origPos1, origPos2);
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

void getTriosP(void(*passTrioFn)(ATrio *trio, void *extraArgument1, void *extraArgument2, void *extraArgument3),
               void *extraArgument1, void *extraArgument2, void *extraArgument3, 
               int *bytesRead, int *nBytes, char **cA, FILE *fileHandle, 
               struct List *speciesList, char *treeString) {
    int32_t length, start, seqLength, i, j, k, l, pos1, pos2, inc1, inc2;
    int32_t pos3, inc3, top;
    struct List *ranges = constructEmptyList(0, free);
    char *seqName;
    char *sequence, *sequence2, *sequence3;
    char strand;
    static ATrio aTrio;

    TrioDecoder *decoder = trioDecoder_construct(treeString);

    //process block function iterates through successive lines while we have not reached
    //the end of the file and the newline starts with 's'
    //for each line grep the sequence, start position and the length.
    //all the lengths must be equal.
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

        listAppend(ranges, seqName);

        assert(strand == '+' || strand == '-');
        if (strand == '+') {
            listAppend(ranges, constructInt(start));
            listAppend(ranges, constructInt(1));
        } else {
            listAppend(ranges, constructInt(seqLength - 1 - start));
            listAppend(ranges, constructInt(-1));
        }
        listAppend(ranges, sequence);

        *bytesRead = benLine(cA, nBytes, fileHandle);
    }

    struct stringSortIdx tmpNameArray[3];
    TrioNames *trio = NULL;
    int32_t index = 0;
    int32_t flag = 0;

    int32_t seq1_rowIndex = 0;
    int32_t seq2_rowIndex = 0;
    int32_t seq3_rowIndex = 0;

    //Now call the pair function for every pair of aligned bases.
    for (i = 0; i < ranges->length; i += 4) {
        char *seq1 = ranges->list[i];
        pos1 = *((int32_t *) ranges->list[i + 1]);
        inc1 = *((int32_t *) ranges->list[i + 2]);
        sequence = ranges->list[i + 3];

        for (j = i + 4; j < ranges->length; j += 4) {
            char *seq2 = ranges->list[j];
            pos2 = *((int32_t *) ranges->list[j + 1]);
            inc2 = *((int32_t *) ranges->list[j + 2]);
            sequence2 = ranges->list[j + 3];

            for (k = j + 4; k < ranges->length; k += 4) {
                char *seq3 = ranges->list[k];
                pos3 = *((int32_t *) ranges->list[k + 1]);
                inc3 = *((int32_t *) ranges->list[k + 2]);
                sequence3 = ranges->list[k + 3];

                assert((int32_t)strlen(sequence) == length);
                assert((int32_t)strlen(sequence2) == length);
                assert((int32_t)strlen(sequence3) == length);

                flag = 0;

                seq1_rowIndex = i / 4;
                seq2_rowIndex = j / 4;
                seq3_rowIndex = k / 4;

                tmpNameArray[0].name = seq1;
                tmpNameArray[1].name = seq2;
                tmpNameArray[2].name = seq3;
                tmpNameArray[0].index = seq1_rowIndex;
                tmpNameArray[1].index = seq2_rowIndex;
                tmpNameArray[2].index = seq3_rowIndex;

                qsort(tmpNameArray, 3, sizeof(struct stringSortIdx), stringSortIdx_cmp);
                for (index = 0; index < speciesList->length; index++) {
                    trio = speciesList->list[index];
                    if (cStr_startsWith(tmpNameArray[0].name, trio->speciesA, 1) == 1 && 
                        cStr_startsWith(tmpNameArray[1].name, trio->speciesB, 1) == 1 && 
                        cStr_startsWith(tmpNameArray[2].name, trio->speciesC, 1) == 1) {
                        flag = 1;
                        break;
                    }
                }

                if (!flag) {
                    continue;
                }

                top = calcTrioState(decoder, tmpNameArray[0].index, tmpNameArray[1].index, tmpNameArray[2].index);
                assert(top != -1);

                for (l = 0; l < length; l++) {
                    if (sequence[l] != '-' && sequence2[l] != '-' && sequence3[l] != '-') {
                        aTrio_fillOut(&aTrio, seq1, seq2, seq3, pos1, pos2, pos3, top);
                        passTrioFn(&aTrio, extraArgument1, extraArgument2, extraArgument3);
                    }
                    if (sequence[l] != '-') {
                        pos1 += inc1;
                    }
                    if (sequence2[l] != '-') {
                        pos2 += inc2;
                    }
                    if (sequence3[l] != '-') {
                        pos3 += inc3;
                    }
                }
            }
        }
    }

    //cleanup the list
    destructList(ranges);
    trioDecoder_destruct(decoder);
}

void getPairs(const char *mAFFile1, 
              void(*passPairFn)(APair *pair, stHash *intervalsHash, void *extraArgument1, void *extraArgument2, 
                                void *extraArgument3, int32_t isVerboseFailures, int32_t near), 
              stHash *intervalsHash,
              void *extraArgument1, void *extraArgument2, void *extraArgument3, 
              int32_t isVerboseFailures, int32_t near) {
    /*
     * Iterates through all pairs of bases in a set of MAFs, calling the given function for each one.
     */
    FILE *fileHandle = fopen(mAFFile1, "r");
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

void getTrios(const char *mAFFile1, 
              void(*passTrioFn)(ATrio *trio, void *extraArgument1, void *extraArgument2, void *extraArgument3), 
              void *extraArgument1, void *extraArgument2, void *extraArgument3, struct List *speciesList) {
    /*
     * Iterates through all trios of bases in a set of MAFs, calling the given function for each one.
     */
    FILE *fileHandle = fopen(mAFFile1, "r");
    int bytesRead;
    int nBytes = 100;
    char *cA;

    char *treeString = NULL;
    char *match = NULL;
    int sscanfResultCnt = 0;

    cA = st_malloc(nBytes + 1);
    bytesRead = benLine(&cA, &nBytes, fileHandle);

    //read through lines until reaching a line starting with an 'a':
    //then call process block function passing it the next line.
    while (bytesRead != -1) {
        if (bytesRead > 0 && cA[0] == 'a') {
            // Parse out the newick tree from the alignment block line
            // Sample line:
            //   a score=5 tree="(A,B)R;"
            match = strstr(cA, "tree");
            assert(match != NULL);
            treeString = st_malloc(sizeof(char) * (bytesRead + 1));
            sscanfResultCnt = sscanf(match, "tree=%*[\'\"]%[^\'\"]", treeString);
            assert(sscanfResultCnt != 0);
            getTriosP(passTrioFn, extraArgument1, extraArgument2, extraArgument3, &bytesRead, &nBytes, &cA, fileHandle,
                    speciesList, treeString);
        } else { //deal with empty, white space lines.
            bytesRead = benLine(&cA, &nBytes, fileHandle);
        }
    }

    //Cleanup
    free(cA);
    free(treeString);
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

void countTrios(ATrio *trio, int64_t *counter, struct hashtable *legitTrios, void *a) {
    /*
     * Counts the number of pairs in the MAF file.
     */
    assert(trio != NULL);
    assert(a == NULL);
    if (hashtable_search(legitTrios, trio->seq1) != NULL)
        if (hashtable_search(legitTrios, trio->seq2) != NULL)
            if (hashtable_search(legitTrios, trio->seq3) != NULL)
                (*counter)++;
}

void samplePairs(APair *thisPair, stHash *intervalsHash, struct avl_table *pairs, 
                 double *acceptProbability, stSortedSet *legitPairs, int32_t isVerboseFailures, 
                 int32_t near) {
    /*
     * Adds *thisPair to *pairs with a given probability.
     */
    if (stSortedSet_search(legitPairs, thisPair->seq1) != NULL)
        if (stSortedSet_search(legitPairs, thisPair->seq2) != NULL)
            if (st_random() <= *acceptProbability)
                avl_insert(pairs, aPair_copyConstruct(thisPair));
}

void sampleTrios(ATrio *trio, struct avl_table *trios, double *acceptProbability, struct hashtable *legitTrios) {
    /*
     * Adds *trio to *trios with a given probability.
     */
    if (hashtable_search(legitTrios, trio->seq1) != NULL)
        if (hashtable_search(legitTrios, trio->seq2) != NULL)
            if (hashtable_search(legitTrios, trio->seq3) != NULL)
                if (st_random() >= *acceptProbability)
                    avl_insert(trios, aTrio_copyConstruct(trio));
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
                    stSortedSet *positivePairs, stSortedSet *legitPairs, int32_t isVerboseFailures, 
                    int32_t near) {
    /*
     * If both members of *thisPair are in the intersection of MAFFileA and MAFFileB,
     * and *thisPair is in the set *pairs then adds to the result pair a positive result.
     */
    if ((stSortedSet_search(legitPairs, thisPair->seq1) != NULL)
        && (stSortedSet_search(legitPairs, thisPair->seq2) != NULL)) {
        getNearPairs(thisPair, pairs, near, positivePairs);
    }
}

void homologyTestsTrio1(ATrio *trio, struct avl_table *trios, struct avl_table *resultTrios,
                        struct hashtable *legitTrios) {
    /*
     * If all 3 members of *trio are in the intersection of MAFFileA and MAFFileB,
     * and *trio is in the set *trios and tops are equal then adds to the result trio a positive result.
     */

    int32_t idx = -1;
    if ((hashtable_search(legitTrios, trio->seq1) != NULL) && 
        (hashtable_search(legitTrios, trio->seq2) != NULL) &&
        (hashtable_search(legitTrios, trio->seq3) != NULL)) {
        ATrio *resultTrio = avl_find(resultTrios, trio);
        if (resultTrio == NULL) {
            resultTrio = aTrio_construct(trio->seq1, trio->seq2, trio->seq3, 0, 0, 0, -1);
            avl_insert(resultTrios, resultTrio);
        }

        ATrio *triosTrio = avl_find(trios, trio);
        if (triosTrio != NULL) { //we use the positions as indices as to the outcome of the trio test.
            idx = encodeTopoMatIdx(trio->top, triosTrio->top);
            resultTrio->topMat[idx]++;
            resultTrio->pos1++;
        }
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
              fprintf(stderr, "%s\t%d\t%d\t%s\t%d\t%d\n", pair->seq1, pair->pos1, 
                      pair->origPos1, pair->seq2, pair->pos2, pair->origPos2);
           }
        }
    }
}

void homologyTestsTrio2(struct avl_table *trios, struct avl_table *resultTrios) {
    /*
     * For every trio in 'trios', add 1 to the total number of trios tests for the sequence-trio.
     * We don't bother with the  seqNames hashtable here because the ht was used to build *pairs
     * in the first place.
     */
    static struct avl_traverser iterator;
    avl_t_init(&iterator, trios);
    ATrio *trio;
    while ((trio = avl_t_next(&iterator)) != NULL) {
        ATrio *resultTrio = avl_find(resultTrios, trio);
        if (resultTrio == NULL) {
            resultTrio = aTrio_construct(trio->seq1, trio->seq2, trio->seq3, 0, 0, 0, -1);
            avl_insert(resultTrios, resultTrio);
        }
        resultTrio->pos2++;
    }
}

struct avl_table *compareMAFs_AB(const char *mAFFileA, const char *mAFFileB, int32_t numberOfSamples,
                                 stSortedSet *legitimateSequences, stHash *intervalsHash, 
                                 int32_t isVerboseFailures, int32_t near) {
    /*
     * Gets samples.
     */
    int64_t pairNumber = 0;
    getPairs(mAFFileA, (void(*)(APair *, stHash *, void *, void *, void *, int32_t, int32_t)) countPairs,
             intervalsHash, &pairNumber, legitimateSequences, NULL, isVerboseFailures, near);

    double acceptProbability = ((double) numberOfSamples) / pairNumber;
    struct avl_table *pairs = avl_create((int32_t(*)(const void *, const void *, void *)) aPair_cmpFunction, 
                                         NULL, NULL);
    getPairs(mAFFileA, (void(*)(APair *, stHash *, void *, void *, void *, int32_t, int32_t)) samplePairs,
             intervalsHash, pairs, &acceptProbability, legitimateSequences, isVerboseFailures, near);

    stSortedSet *positivePairs = stSortedSet_construct();
    getPairs(mAFFileB, (void(*)(APair *, stHash *, void *, void *, void *, int32_t, int32_t)) homologyTests1,
             intervalsHash, pairs, positivePairs, legitimateSequences, isVerboseFailures, near);

    struct avl_table *resultPairs = avl_create((int32_t(*)(const void *, const void *, void *)) aPair_cmpFunction_seqsOnly, NULL, NULL);

    homologyTests2(pairs, resultPairs, intervalsHash, positivePairs, isVerboseFailures);
    avl_destroy(pairs, (void(*)(void *, void *)) aPair_destruct);
    stSortedSet_destruct(positivePairs);
    return resultPairs;
}

struct avl_table *compareMAFs_AB_Trio(const char *mAFFileA, const char *mAFFileB, int32_t numberOfSamples,
                                      struct hashtable *ht, struct List *speciesList) {
    /*
     * Gets samples.
     */
    int64_t trioNumber = 0;
    getTrios(mAFFileA, (void(*)(ATrio *, void *, void *, void *)) countTrios, &trioNumber, ht, NULL, speciesList);

    double acceptProbability = 1.0 - ((double) numberOfSamples) / trioNumber;
    struct avl_table *trios =
            avl_create((int32_t(*)(const void *, const void *, void *)) aTrio_cmpFunction, NULL, NULL);
    getTrios(mAFFileA, (void(*)(ATrio *, void *, void *, void *)) sampleTrios, trios, &acceptProbability, ht,
            speciesList);

    struct avl_table *resultTrios = avl_create(
            (int32_t(*)(const void *, const void *, void *)) aTrio_cmpFunction_seqsOnly, NULL, NULL);
    getTrios(mAFFileB, (void(*)(ATrio *, void *, void *, void *)) homologyTestsTrio1, trios, resultTrios, ht,
            speciesList);

    homologyTestsTrio2(trios, resultTrios);
    avl_destroy(trios, (void(*)(void *, void *)) aTrio_destruct);
    return resultTrios;
}

void reportResult(const char *tagName, double total, double totalTrue, FILE *fileHandle, int tabLevel) {
    assert(total >= totalTrue);
    int i;
    for (i = 0; i < tabLevel; i++)
        fprintf(fileHandle, "\t");
    fprintf(fileHandle, "<%s totalTests=\"%i\" totalTrue=\"%i\" totalFalse=\"%i\" average=\"%f\"/>\n",
            tagName, (int32_t) total, (int32_t) totalTrue, 
            (int32_t) (total - totalTrue), total == 0 ? 0.0 : totalTrue / total);
}

ResultPair *aggregateResult(void *(*getNextPair)(void *, void *), void *arg1, void *arg2, 
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

void addReferencesAndDups(struct avl_table *results_AB, stSortedSet *legitimateSequences) {
    /*
     * Adds tags for all against each species.
     */
    stSortedSetIterator *speciesIt = stSortedSet_getIterator(legitimateSequences);
    static struct avl_traverser iterator;
    char *species;
    avl_t_init(&iterator, results_AB);
    stList *list = stList_construct();
    stList_append(list, aggregateResult(addReferencesAndDups_getDups, &iterator, NULL, "self", "self"));
    //Add references
    while((species = stSortedSet_getNext(speciesIt)) != NULL) {
        avl_t_init(&iterator, results_AB);
        stList_append(list, aggregateResult(addReferencesAndDups_getReferences, &iterator, 
                                            species, species, "aggregate"));
    }
    stSortedSet_destructIterator(speciesIt);
    for(int32_t i=0; i<stList_length(list); i++) {
        avl_insert(results_AB, stList_get(list, i));
    }
    stList_destruct(list);
}

void *reportResults_fn(void *arg, void *arg2) {
   return avl_t_next(arg);
}

void reportResults(struct avl_table *results_AB, const char *mAFFileA, const char *mAFFileB, FILE *fileHandle,
                   int32_t near, stSortedSet *legitimateSequences, const char *bedFiles) {
    /*
     * Report results in an XML formatted document.
     */
    static struct avl_traverser iterator;
    avl_t_init(&iterator, results_AB);
    ResultPair *resultPair;

    ResultPair *aggregateResults = aggregateResult(reportResults_fn, &iterator, NULL, "", "");

    addReferencesAndDups(results_AB, legitimateSequences);

    fprintf(fileHandle, "\t<homologyTests fileA=\"%s\" fileB=\"%s\">\n"
            "\t\t<aggregateResults>\n", 
            mAFFileA, mAFFileB);
    reportResult("all", aggregateResults->total, aggregateResults->inAll, fileHandle, 3);
    if (bedFiles != NULL){
        reportResult("both", aggregateResults->totalBoth, aggregateResults->inBoth, fileHandle, 3);
        reportResult("A", aggregateResults->totalA, aggregateResults->inA, fileHandle, 3);
        reportResult("B", aggregateResults->totalB, aggregateResults->inB, fileHandle, 3);
        reportResult("neither", aggregateResults->totalNeither, aggregateResults->inNeither, fileHandle, 3);
    }
    fprintf(fileHandle, "\t\t</aggregateResults>\n\t\t<homologyPairTests>\n");

    while ((resultPair = avl_t_prev(&iterator)) != NULL) {
        fprintf(fileHandle, "\t\t\t<homologyTest sequenceA=\"%s\" sequenceB=\"%s\">\n"
                "\t\t\t\t<aggregateResults>\n",
                resultPair->aPair.seq1, resultPair->aPair.seq2);
        reportResult("all", resultPair->total, resultPair->inAll, fileHandle, 5);
        if (bedFiles != NULL){
            reportResult("both", resultPair->totalBoth, resultPair->inBoth, fileHandle, 5);
            reportResult("A", resultPair->totalA, resultPair->inA, fileHandle, 5);
            reportResult("B", resultPair->totalB, resultPair->inB, fileHandle, 5);
            reportResult("neither", resultPair->totalNeither, resultPair->inNeither, fileHandle, 5);
        }
        fprintf(fileHandle, "\t\t\t\t</aggregateResults>\n\t\t\t\t<singleHomologyTests>\n");
        fprintf(fileHandle, "\t\t\t\t\t<singleHomologyTest sequenceA=\"%s\" sequenceB=\"%s\">\n"
                "\t\t\t\t\t\t<aggregateResults>\n",
                resultPair->aPair.seq1, resultPair->aPair.seq2);
        reportResult("all", resultPair->total, resultPair->inAll, fileHandle, 6);
        if (bedFiles != NULL){
            reportResult("both", resultPair->totalBoth, resultPair->inBoth, fileHandle, 6);
            reportResult("A", resultPair->totalA, resultPair->inA, fileHandle, 6);
            reportResult("B", resultPair->totalB, resultPair->inB, fileHandle, 6);
            reportResult("neither", resultPair->totalNeither, resultPair->inNeither, fileHandle, 6);
        }
        fprintf(fileHandle,"\t\t\t\t\t\t</aggregateResults>\n"
                "\t\t\t\t\t</singleHomologyTest>\n"
                "\t\t\t\t</singleHomologyTests>\n"
                "\t\t\t</homologyTest>\n");
    }

    fprintf(fileHandle, "\t\t</homologyPairTests>\n\t</homologyTests>\n");
    return;
}

void reportResultsTrio(struct avl_table *results_AB, const char *mAFFileA, 
                       const char *mAFFileB, FILE *fileHandle) {
    /*
     * Report results in an XML formatted document.
     */
    int32_t i, j;
    static struct avl_traverser iterator;
    avl_t_init(&iterator, results_AB);
    ATrio *resultTrio;
    double positive = 0.0;
    double total = 0.0;
    while ((resultTrio = avl_t_next(&iterator)) != NULL) {
        assert(resultTrio->pos2 >= resultTrio->pos1);
        positive += resultTrio->pos1;
        total += resultTrio->pos2;
    }

    //  char *topString = NULL;
    double goodDiagonal = 0.0;
    double goodTriangle = 0.0;
    double totalTriangle = 0.0;

    fprintf(
            fileHandle,
            "\t<trioTests fileA=\"%s\" fileB=\"%s\" totalTests=\"%.0f\" "
            "totalTrue=\"%.0f\" totalFalse=\"%.0f\" average=\"%f\">\n",
            mAFFileA, mAFFileB, total, positive, total - positive, positive / total);
    while ((resultTrio = avl_t_prev(&iterator)) != NULL) {
        assert(resultTrio->pos2 >= resultTrio->pos1);
        fprintf(
                fileHandle,
                "\t\t<trioTest sequenceA=\"%s\" sequenceB=\"%s\" sequenceC=\"%s\" "
                "totalTests=\"%i\" totalTrue=\"%i\" totalFalse=\"%i\" average=\"%f\">\n",
                resultTrio->seq1, resultTrio->seq2, resultTrio->seq3, resultTrio->pos2, resultTrio->pos1,
                resultTrio->pos2 - resultTrio->pos1, ((double) resultTrio->pos1) / resultTrio->pos2);

        //      topString = arrayToString(resultTrio->topMat, 10, ',');
        goodDiagonal = 0.0;
        goodTriangle = 0.0;
        totalTriangle = 0.0;
        for (i = 0; i < 3; i++) {
            goodDiagonal += resultTrio->topMat[encodeTopoMatIdx(i, i)];
        }
        for (i = 0; i < 3; i++) {
            for (j = i; j < 3; j++) {
                goodTriangle += resultTrio->topMat[encodeTopoMatIdx(i, j)];
            }
        }
        for (i = 0; i < 4; i++) {
            for (j = i; j < 4; j++) {
                totalTriangle += resultTrio->topMat[encodeTopoMatIdx(i, j)];
            }
        }

        fprintf(
                fileHandle,
                "\t\t\t<singleTrioTest sequenceA=\"%s\" sequenceB=\"%s\" "
                "sequenceC=\"%s\" totalTests=\"%i\" totalTrue=\"%i\" totalFalse=\"%i\" "
                "average=\"%f\" topGoodAvg=\"%f\" topAllAvg=\"%f\">\n",
                resultTrio->seq1, resultTrio->seq2, resultTrio->seq3, resultTrio->pos2, resultTrio->pos1,
                resultTrio->pos2 - resultTrio->pos1, ((double) resultTrio->pos1) / resultTrio->pos2, goodDiagonal
                        / goodTriangle, goodDiagonal / totalTriangle);

        for (i = 0; i < 4; i++) {
            for (j = i; j < 4; j++) {
                fprintf(fileHandle,
                        "\t\t\t\t<trioTestInstance topStateA=\"%d\" "
                        "topStateB = \"%d\" count = \"%d\"/>\n", i, j,
                        resultTrio->topMat[encodeTopoMatIdx(i, j)]);
            }
        }

        fprintf(fileHandle, "\t\t\t</singleTrioTest>\n");
        //      free(topString);
        //      topString = NULL;

        fprintf(fileHandle, "\t\t</trioTest>\n");
    }
    fprintf(fileHandle, "\t</trioTests>\n");

    return;
}

void populateNames(const char *mAFFile, stSortedSet *htp) {
    /*
     * populates a hash table with the names of sequences from a MAF file.
     */
    FILE *fileHandle = fopen(mAFFile, "r");
    int bytesRead;
    int nBytes = 100;
    char *cA;
    int32_t start, seqLength, i, j;
    char *seqName;
    char *sequence;
    char strand;
    cA = st_malloc(nBytes + 1);
    bytesRead = benLine(&cA, &nBytes, fileHandle);

    //read through lines until reaching a line starting with an 's':
    while (bytesRead != -1) {
        if (bytesRead > 0 && cA[0] == 's') {
            seqName = st_malloc(sizeof(char) * (1 + (bytesRead)));
            sequence = st_malloc(sizeof(char) * (1 + (bytesRead)));
            //assert(sscanf(cA, "s %s %i %i %c %i %s", seqName, &start, &i /*ignore the length field*/, &strand, &seqLength, sequence) == 6);
            j = sscanf(cA, "s %s %i %i %c %i %s", seqName, &start, &i, /*ignore the length field*/
                       &strand, &seqLength, sequence);
            assert(j == 6 || (j == 5 && seqLength == 0));

            if (stSortedSet_search(htp, seqName) == NULL) {
                stSortedSet_insert(htp, seqName);
            }
            free(sequence);
        }
        bytesRead = benLine(&cA, &nBytes, fileHandle);
    }

    //Cleanup
    fclose(fileHandle);
    free(cA);
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

int32_t calcTrioState(TrioDecoder *decoder, int32_t spAIdx, int32_t spBIdx, int32_t spCIdx) {
    int32_t lca_AB = -1, lca_AC = -1, lca_BC = -1;

    char *spA = decoder->leafLabelArray[spAIdx];
    char *spB = decoder->leafLabelArray[spBIdx];
    char *spC = decoder->leafLabelArray[spCIdx];

    int32_t aIdx = *(int *) hashtable_search(decoder->treeLabelHash, (void *) spA);
    int32_t bIdx = *(int *) hashtable_search(decoder->treeLabelHash, (void *) spB);
    int32_t cIdx = *(int *) hashtable_search(decoder->treeLabelHash, (void *) spC);

    lca_AB = decoder->lcaMatrix[aIdx][bIdx];
    lca_AC = decoder->lcaMatrix[aIdx][cIdx];
    lca_BC = decoder->lcaMatrix[bIdx][cIdx];

    if (lca_AB == lca_AC && lca_AC == lca_BC) {
        /* Handle the case with multi-furcations */
        //printf("3: (%d,%d,%d);\n", spAIdx, spBIdx, spCIdx);
        return 3;
    } else if (lca_AC == lca_BC) {
        //printf("0: ((%d,%d),%d);\n", spAIdx, spBIdx, spCIdx);
        return 0;
    } else if (lca_AB == lca_AC) {
        //printf("1: (%d,(%d,%d));\n", spAIdx, spBIdx, spCIdx);
        return 1;
    } else if (lca_AB == lca_BC) {
        //printf("2: ((%d,%d),%d);\n", spAIdx, spCIdx, spBIdx);
        return 2;
    } else {
        //printf("?\n");
        return -1;
    }
}

TrioDecoder *trioDecoder_construct(char *treestring) {
    stTree *tree = NULL;
    tree = stTree_parseNewickString(treestring);

    int32_t nodeNum = countNodes(tree);
    int32_t leafNum = countLeaves(tree);

    char **nodeLabelArray = createNodeLabelArray(tree, nodeNum);
    char **leafLabelArray = createLeafLabelArray(tree, leafNum);
    struct hashtable *treeLabelHash = getTreeLabelHash(nodeLabelArray, nodeNum);

    int32_t **lcaMatrix = NULL;
    lcaMatrix = lca(tree, treeLabelHash);

    TrioDecoder *decoder = st_malloc(sizeof(TrioDecoder));
    decoder->nodeLabelArray = nodeLabelArray;
    decoder->leafLabelArray = leafLabelArray;
    decoder->treeLabelHash = treeLabelHash;
    decoder->lcaMatrix = lcaMatrix;
    decoder->nodeNum = nodeNum;
    decoder->leafNum = leafNum;

    stTree_destruct(tree);
    return decoder;
}

void trioDecoder_destruct(TrioDecoder *decoder) {
    labelArray_destruct(decoder->nodeLabelArray, decoder->nodeNum);
    labelArray_destruct(decoder->leafLabelArray, decoder->leafNum);
    hashtable_destroy(decoder->treeLabelHash, 1, 1);
    lcaMatrix_free(decoder->lcaMatrix, decoder->nodeNum);
    free(decoder);

    return;
}

void writeXMLHeader(FILE *fileHandle){
   fprintf(fileHandle, "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\" ?>\n");
   return;
}
