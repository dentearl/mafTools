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
#include <assert.h>
#include <getopt.h>
#include <inttypes.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "common.h"
#include "sharedMaf.h"
#include "mafBlockExtractorAPI.h"

bool checkRegion(uint32_t targetStart, uint32_t targetStop, uint32_t lineStart, 
                 uint32_t length, uint32_t sourceLength, char strand) {
    // check to see if pos is in this block
    uint32_t absStart, absEnd;
    if (strand == '-') {
        absStart = sourceLength - (lineStart + length);
        absEnd = sourceLength - 1 - lineStart;
    } else {
        absStart = lineStart;
        absEnd = lineStart + length - 1;
    }
    if (absEnd < targetStart)
        return false;
    if (targetStop < absStart)
        return false;
    if ((absStart <= targetStart) && (targetStart <= absEnd))
        return true;
    if ((absStart <= targetStop) && (targetStop <= absEnd))
        return true;
    if ((targetStart <= absStart) && (absEnd <= targetStop))
        return true;
    return false;
}
bool searchMatched(mafLine_t *ml, const char *seq, uint32_t start, uint32_t stop) {
    // report false if search did not match, true if it did
    if (maf_mafLine_getType(ml) != 's')
        return false;
    if (!(strcmp(maf_mafLine_getSpecies(ml), seq) == 0))
        return false;
    if (checkRegion(start, stop, maf_mafLine_getStart(ml), maf_mafLine_getLength(ml), 
                    maf_mafLine_getSourceLength(ml), maf_mafLine_getStrand(ml)))
        return true;
    return false;
}
void printHeader(void) {
    printf("##maf version=1\n\n");
}
void printTargetColumns(bool *targetColumns, uint32_t n) {
    for (uint32_t i = 0; i < n; ++i) {
        if (targetColumns[i])
            printf("1");
        else 
            printf("0");
        // if (!((i + 1) % 5))
        //     printf(" ");
    }
    printf("\n");
}
uint32_t getTargetColumns(bool **targetColumns, uint32_t *len, mafBlock_t *b, const char *seqName,
                          uint32_t start, uint32_t stop) {
    // given a block and a target, create an array of bools where true means the target
    // is present in that column and false means it is absent
    mafLine_t *ml = maf_mafBlock_getHeadLine(b);
    uint32_t sum = 0;
    while(maf_mafLine_getType(ml) != 's') {
        ml = maf_mafLine_getNext(ml);
    }
    char *seq = NULL;
    *len = maf_mafBlock_getSequenceFieldLength(b);
    // printf("target columns len: %" PRIu32 "\n", *len);
    if (*len == 0) {
        *targetColumns = NULL;
        return sum;
    }
    *targetColumns = (bool*) de_malloc(sizeof(bool*) * (*len));
    memset(*targetColumns, false, sizeof(bool*) * (*len));
    // printf("target columns len: %" PRIu32 "\n", *len);
    int it = 0;
    uint64_t pos = 0;
    while (ml != NULL) {
        if (maf_mafLine_getType(ml) != 's') {
            ml = maf_mafLine_getNext(ml);
            continue;
        }
        if (maf_mafLine_getSpecies(ml) == NULL) {
            ml = maf_mafLine_getNext(ml);
            continue;
        }
        // make sure name is the same
        if (strcmp(maf_mafLine_getSpecies(ml), seqName) != 0) {
            ml = maf_mafLine_getNext(ml);
            continue;
        }
        if (maf_mafLine_getStrand(ml) == '+') {
            it = 1;
        } else {
            it = -1;
        }
        seq = maf_mafLine_getSequence(ml);
        pos = maf_mafLine_getPositiveCoord(ml);
        for (uint32_t i = 0; i < (*len); ++i) {
            if (seq[i] != '-') {
                if ((start <= pos) && (pos <= stop)) {
                    if ((*targetColumns)[i] == 0) {
                        ++sum;
                        (*targetColumns)[i] = 1;
                    }
                }
                pos += it;
            }
        }
        ml = maf_mafLine_getNext(ml);
    }
    return sum;
}
uint32_t sumBool(bool *array, uint32_t n) {
    uint32_t a = 0;
    for (uint32_t i = 0; i < n; ++i) {
        if (array[i]) {
            // printf("1");
            ++a;
        } else {
            // printf("0");
        }
        /* if ((i + 1) % 5 == 0) { */
        /*     printf(" "); */
        /* } */
    }
    // printf("\n");
    return a;
}
void printBoolArray(bool *targetColumns, uint32_t len) {
    printf("printBoolArray: len=%" PRIu32 "\n", len);
    for (uint32_t i = 0; i < len; ++i) {
        if (targetColumns[i]) {
            printf("1");
        } else {
            printf("0");
        }
        if (i + 1 % 5 == 0) {
            printf(" ");
        }
    }
    printf("\n");
}
mafBlock_t *processBlockForTrim(mafBlock_t *b, const char *seq, uint32_t start, uint32_t stop, bool isLeft) {
    // take in a prospective block, and some target information
    // create a targetColumns array and then check to see if this block needs to be trimmed
    bool *targetColumns = NULL;
    uint32_t len = 0, n = 0, sum = 0;
    sum = getTargetColumns(&targetColumns, &len, b, seq, start, stop);
    if (sum == 0) {
        free(targetColumns);
        return NULL;
    }
    if (sum == len) {
        // bail out, nothing to trim here.
        free(targetColumns);
        return b;
    }
    if (isLeft) {
        for (uint32_t i = 0; i < len; ++i) {
            if (!targetColumns[i]) {
                ++n;
            } else {
                break;
            }
        }
    } else {
        for (uint32_t i = len - 1; i < UINT32_MAX; --i) {
            if (!targetColumns[i]) {
                ++n;
            } else {
                break;
            }
        }
    }
    free(targetColumns);
    if (n == 0) {
        // nothing to trim, return a copy of the original
        return maf_copyMafBlock(b);
    } else if (n == maf_mafBlock_getSequenceFieldLength(b)) {
        // this would be the case where every single element in targetColumns is 0,
        // which would mean this block contains nothing relevant.
        return NULL;
    } else {
        // need to trim off `n' bases from the left / right
        // where n is greater than 0 and less than the total length of the sequence
        return trimBlock(b, n, isLeft);
    }
}
mafBlock_t *trimBlock(mafBlock_t *b, uint32_t n, bool isLeft) {
    // create a new mafBlock_t* containing a left- or right-trimmed block,
    // Trim off `n' bases from the left/right of all lines, updating the start, 
    // legth, sequence fields:
    // printf("trimBlock()\n");
    if (((isLeft) && (n == 0)) || 
        ((!isLeft) && (n == maf_mafBlock_getSequenceFieldLength(b)))) {
        return b;
    }
    mafBlock_t *mb = maf_newMafBlock();
    mafLine_t *ml1 = NULL, *ml2 = NULL;
    ml1 = maf_mafBlock_getHeadLine(b);
    uint32_t lineNumber = maf_mafLine_getLineNumber(ml1);
    ml1 = maf_mafLine_getNext(ml1);
    maf_mafBlock_setHeadLine(mb, maf_newMafLineFromString("a score=0", lineNumber));
    ml2 = maf_mafBlock_getHeadLine(mb);
    maf_mafBlock_setLineNumber(mb, lineNumber);
    maf_mafBlock_incrementNumberOfLines(mb);
    uint32_t m, len;
    char *seq = NULL;
    bool prevLineUsed = true; // used when a mafline is dropped because it's length becomes 0
    bool emptyBlock = true;
    while (ml1 != NULL) {
        // loop through all maf lines in a block
        if (prevLineUsed) {
            // advance ml2
            maf_mafBlock_setTailLine(mb, ml2); // assign this first, in case last line is not used.
            maf_mafLine_setNext(ml2, maf_newMafLine());
            ml2 = maf_mafLine_getNext(ml2);
        }
        if (maf_mafLine_getLine(ml1) == NULL) {
            ml1 = maf_mafLine_getNext(ml1);
            prevLineUsed = false;
            continue;
        }
        maf_mafLine_setType(ml2, maf_mafLine_getType(ml1));
        maf_mafLine_setLineNumber(ml2, ++lineNumber);
        maf_mafBlock_incrementNumberOfLines(mb);
        maf_mafBlock_incrementLineNumber(mb);
        if (maf_mafLine_getType(ml2) != 's') {
            maf_mafLine_setLine(ml2, de_strdup(maf_mafLine_getLine(ml1)));
            ml1 = maf_mafLine_getNext(ml1);
            prevLineUsed = true;
            continue;
        }
        maf_mafBlock_incrementNumberOfSequences(mb);
        seq = maf_mafLine_getSequence(ml1);
        m = 0;
        len = maf_mafBlock_getSequenceFieldLength(b);
        if (maf_mafBlock_getSequenceFieldLength(mb) == 0) {
            maf_mafBlock_setSequenceFieldLength(mb, len);
        }
        if (isLeft) {
            for (uint32_t i = 0; i < n; ++i) {
                if (seq[i] != '-') {
                    ++m;
                }
            }
        } else {
            for (uint32_t i = (len - 1); len - 1 - n < i  ; --i) {
                if (seq[i] != '-') {
                    ++m;
                }
            }
        }
        if (maf_mafLine_getLength(ml1) <= m) {
            ml1 = maf_mafLine_getNext(ml1);
            prevLineUsed = false;
            maf_mafBlock_decrementLineNumber(mb);
            maf_mafBlock_decrementNumberOfLines(mb);
            maf_mafBlock_decrementNumberOfSequences(mb);
            --lineNumber;
            continue;
        }
        if (isLeft) {
            assert(maf_mafLine_getStart(ml1) + m < maf_mafLine_getSourceLength(ml1));
            maf_mafLine_setStart(ml2, maf_mafLine_getStart(ml1) + m);
            maf_mafLine_setSequence(ml2, de_strdup(seq + n));
        } else {
            maf_mafLine_setStart(ml2, maf_mafLine_getStart(ml1));
            if (len <= n) {
                fprintf(stderr, "Well, this is bullshit.\n");
                fprintf(stderr, "Trying to trim %" PRIu32 " off of a sequence of length %"PRIu32"\n", n, len);
                fprintf(stderr, "line: %s\n", seq);
                exit(EXIT_FAILURE);
            }
            maf_mafLine_setSequence(ml2, de_strndup(seq, len - n));
        }
        maf_mafLine_setLength(ml2, maf_mafLine_getLength(ml1) - m);
        maf_mafBlock_setSequenceFieldLength(mb, strlen(maf_mafLine_getSequence(ml2)));
        maf_mafLine_setStrand(ml2, maf_mafLine_getStrand(ml1));
        maf_mafLine_setSourceLength(ml2, maf_mafLine_getSourceLength(ml1));
        maf_mafLine_setSpecies(ml2, de_strdup(maf_mafLine_getSpecies(ml1)));
        maf_mafLine_setLine(ml2, maf_mafLine_imputeLine(ml2));
        prevLineUsed = true;
        ml1 = maf_mafLine_getNext(ml1);
        emptyBlock = false;
    }
    maf_mafBlock_incrementLineNumber(mb); // extra \n at the end of a block
    if (prevLineUsed) {
        maf_mafBlock_setTailLine(mb, ml2);
    }
    if (!emptyBlock) {
        return mb;
    } else {
        // this condition should never be tripped, we have a condition set up
        // in the "process" wrapper above this function.
        printf("block was empty, this should never happen in trimBlock()\n");
        assert(2 + 2 == 5);
        maf_destroyMafBlockList(mb);
        return NULL;
    }
}
void processBlockForSplit(mafBlock_t *b, const char *seq, uint32_t start, uint32_t stop, 
                          mafBlock_t **left, mafBlock_t **right) {
    bool *targetColumns = NULL;
    uint32_t len = 0, n = 0, sum = 0;
    sum = getTargetColumns(&targetColumns, &len, b, seq, start, stop);
    for (uint32_t i = 0; i < len; ++i) {
        if (targetColumns[i]) {
            ++n;
        } else {
            break;
        }
    }
    free(targetColumns);
    if (sum == len) {
        // bail out, nothing to trim here.
        *left = b;
        *right = NULL;
    } else {
        // split the block in twain
        splitBlock(b, n, left, right);
    }
}
void splitBlock(mafBlock_t *b, uint32_t n, mafBlock_t **mbLeft, mafBlock_t **mbRight) {
    /* n is the number of bases to be included in the left block
     * 
     */
    assert(*mbLeft == NULL);
    assert(*mbRight == NULL);
    uint32_t m = maf_mafBlock_getSequenceFieldLength(b);
    if ((n > 0) && (n < m - 1)) {
        // Yes, this looks confusing, but think of it this way:
        // the left block must be right-trimmed and the right block must be left-trimmed.
        *mbLeft = trimBlock(b, m - n, false);
        *mbRight = trimBlock(b, n, true);
    } else {
        // either n == 0 or n == m - 1, i.e. do nothing.
        // printf("closing in\n");
        *mbLeft = maf_copyMafBlock(b);
    }
}
int32_t **createOffsets(uint32_t n) {
    // create a matrix to store relative offsets for walking through a maf block
    int32_t **offs = (int32_t**) de_malloc(sizeof(int32_t*) * n);
    for (uint32_t i = 0; i < n; ++i) {
        offs[i] = (int32_t*) de_malloc(sizeof(int32_t) * 2); // 0: index 1: offset
        offs[i][0] = 0;
        offs[i][1] = 0;
    }
    return offs;
}
void destroyOffsets(int32_t **offs, uint32_t n) {
    for (uint32_t i = 0; i < n; ++i) {
        free(offs[i]);
        offs[i] = NULL;
    }
    free(offs);
    offs = NULL;
}
mafBlock_t *processBlockForSplice(mafBlock_t *b, const char *seq, uint32_t start, uint32_t stop, bool store) {
    // walks mafBlock_t b, returns a mafBlock_t (using the linked list feature) of all spliced out bits.
    // if store is true, will return a mafBlock_t linked list of all sub-blocks. If store is false,
    // will report each sub-block (maf_mafBlock_print()) as it comes in and immediatly destroy that block.
    bool *targetColumns = NULL;
    uint32_t len = 0, sum = 0;
    mafBlock_t *head = NULL, *mb = NULL;
    sum = getTargetColumns(&targetColumns, &len, b, seq, start, stop);
    // printTargetColumns(targetColumns, len);
    int32_t **offsets = createOffsets(maf_mafBlock_getNumberOfSequences(b));
    uint32_t l = 0, r = 0, ri = 0;
    while (l < len) {
        if (!targetColumns[l]) {
            // find the left most element
            ++l;
            r = l;
            continue;
        }
        while (targetColumns[r] && r < len) {
            // find the end of the right 
            ++r;
        }
        // set ri equal to the index of the last element
        ri = r - 1;
        if (store) {
            // used in unit tests
            if (head == NULL) {
                head = spliceBlock(b, l, ri, offsets);
                mb = head;
            } else {
                maf_mafBlock_setNext(mb, spliceBlock(b, l, ri, offsets));
                mb = maf_mafBlock_getNext(mb);
            }
        } else {
            // used in production
            mb = spliceBlock(b, l, ri, offsets);
            maf_mafBlock_print(mb);
            if (mb != b) {
                maf_destroyMafBlockList(mb);
            }
        }
        l = r;
        if (l == len - 1) { 
            break;
        }
    }
    // clean up
    free(targetColumns);
    destroyOffsets(offsets, maf_mafBlock_getNumberOfSequences(b));
    return head;
}
void printOffsetArray(int32_t **offsetArray, uint32_t n) {
    for (uint32_t i = 0; i < n; ++i) {
        printf("(%" PRIi32 ", %" PRIi32 "), ", offsetArray[i][0], offsetArray[i][1]);
    }
}
mafBlock_t *spliceBlock(mafBlock_t *b, uint32_t l, uint32_t r, int32_t **offsetArray) {
    // b is the input maf block
    // l is the left index in the sequence field, the start of inclusion
    // r is the right index in the sequence field, the stop of inclusion
    if ((l == 0) && (r == maf_mafBlock_getSequenceFieldLength(b) - 1)) {
        return b;
    }
    // printf("spliceBlock(b, l=%"PRIu32" r%"PRIu32" offsetArray=[", l, r);
    // printOffsetArray(offsetArray, maf_mafBlock_getNumberOfSequences(b));
    // printf("])\n");
    mafBlock_t *mb = maf_newMafBlock();
    mafLine_t *ml1 = NULL, *ml2 = NULL;
    ml1 = maf_mafBlock_getHeadLine(b);
    uint32_t lineNumber = maf_mafLine_getLineNumber(ml1);
    uint32_t numberOfSequences = maf_mafBlock_getNumberOfSequences(b);
    for (uint32_t i = 0; i < numberOfSequences; ++i) {
        assert(offsetArray[i][0] <= (int32_t)l);
    }
    assert(r < maf_mafBlock_getSequenceFieldLength(b));
    ml1 = maf_mafLine_getNext(ml1);
    maf_mafBlock_setHeadLine(mb, maf_newMafLineFromString("a score=0", lineNumber));
    ml2 = maf_mafBlock_getHeadLine(mb);
    maf_mafBlock_setLineNumber(mb, lineNumber);
    maf_mafBlock_incrementNumberOfLines(mb);
    uint32_t len;
    char *seq = NULL;
    bool prevLineUsed = true; // used when a mafline is dropped because it's length becomes 0
    bool emptyBlock = true;
    uint32_t si = 0; // sequence index, for addressing into offsetArray
    while (ml1 != NULL) {
        // loop through all maf lines in a block
        if (prevLineUsed) {
            // advance ml2
            maf_mafBlock_setTailLine(mb, ml2); // assign this first, in case last line is not used.
            maf_mafLine_setNext(ml2, maf_newMafLine());
            ml2 = maf_mafLine_getNext(ml2);
        }
        if (maf_mafLine_getLine(ml1) == NULL) {
            ml1 = maf_mafLine_getNext(ml1);
            prevLineUsed = false;
            continue;
        }
        maf_mafLine_setType(ml2, maf_mafLine_getType(ml1));
        maf_mafLine_setLineNumber(ml2, ++lineNumber);
        maf_mafBlock_incrementNumberOfLines(mb);
        maf_mafBlock_incrementLineNumber(mb);
        if (maf_mafLine_getType(ml2) != 's') {
            // copy the line over, move on to next ml
            maf_mafLine_setLine(ml2, de_strdup(maf_mafLine_getLine(ml1)));
            ml1 = maf_mafLine_getNext(ml1);
            prevLineUsed = true;
            continue;
        }
        maf_mafBlock_incrementNumberOfSequences(mb);
        seq = maf_mafLine_getSequence(ml1);
        len = maf_mafBlock_getSequenceFieldLength(b);
        if (maf_mafBlock_getSequenceFieldLength(mb) == 0) {
            maf_mafBlock_setSequenceFieldLength(mb, len);
        }
        // find first non-gap position, update offset
        while (seq[offsetArray[si][0]] == '-' && 
               offsetArray[si][0] <= (int32_t)l) {
            ++offsetArray[si][0];
            if (seq[offsetArray[si][0]] != '-' && offsetArray[si][1] > 0) {
                // if we've advanced to a non-gap char, advance the offset
                ++offsetArray[si][1]; // advance offset
            }
        }
        // offsets
        for (int32_t i = offsetArray[si][0] + 1; i <= (int32_t)l; ++i) {
            // figure out the non-gap offset for the splice-in point, `l'
            offsetArray[si][0] = i; // local pos
            if (seq[i] != '-') {
                ++offsetArray[si][1]; // advance offset
            }
        }
        bool allGaps = true;
        for (uint32_t i = l; i <= r; ++i) {
            if (seq[i] != '-') {
                allGaps = false;
                break;
            }
        }
        if (allGaps) {
            // this sequence is all gaps in this region, exclude it
            ml1 = maf_mafLine_getNext(ml1);
            prevLineUsed = false;
            maf_mafBlock_decrementLineNumber(mb);
            maf_mafBlock_decrementNumberOfLines(mb);
            maf_mafBlock_decrementNumberOfSequences(mb);
            --lineNumber;
            ++si;
            continue;
        }
        // offsets
        while (seq[offsetArray[si][0]] == '-') {
            ++offsetArray[si][0];
            if (seq[offsetArray[si][0]] != '-') {
                ++offsetArray[si][1];
            }
        }
        if (offsetArray[si][1] == -1) {
            offsetArray[si][1] = 0;
        }
        maf_mafLine_setStart(ml2, maf_mafLine_getStart(ml1) + offsetArray[si][1]);
        // update offsetArray:
        for (uint32_t i = offsetArray[si][0] + 1; i <= r; ++i) {
            offsetArray[si][0] = i;
            if (seq[i] != '-') {
                ++offsetArray[si][1];
            }
        }
        maf_mafLine_setSequence(ml2, de_strndup(seq + l, 1 + r - l));
        maf_mafLine_setLength(ml2, countNonGaps(maf_mafLine_getSequence(ml2)));
        maf_mafBlock_setSequenceFieldLength(mb, strlen(maf_mafLine_getSequence(ml2)));
        maf_mafLine_setStrand(ml2, maf_mafLine_getStrand(ml1));
        maf_mafLine_setSourceLength(ml2, maf_mafLine_getSourceLength(ml1));
        maf_mafLine_setSpecies(ml2, de_strdup(maf_mafLine_getSpecies(ml1)));
        maf_mafLine_setLine(ml2, maf_mafLine_imputeLine(ml2));
        prevLineUsed = true;
        ml1 = maf_mafLine_getNext(ml1);
        ++si;
        emptyBlock = false;
    }
    maf_mafBlock_incrementLineNumber(mb); // extra \n at the end of a block
    if (prevLineUsed) {
        maf_mafBlock_setTailLine(mb, ml2);
    }
    if (!emptyBlock) {
        return mb;
    } else {
        // this condition should never be tripped, we have a condition set up
        // in the "process" wrapper above this function.
        printf("block was empty, this should never happen in spliceBlock()\n");
        assert(2 + 2 == 5);
        maf_destroyMafBlockList(mb);
        return NULL;
    }
}
void trimAndReportBlock(mafBlock_t *orig, const char *seq, uint32_t start, uint32_t stop) {
    // we can not modify orig, it is cleaned up in processBody()
    mafBlock_t *leftTrimmed = NULL, *rightTrimmed = NULL, *copy = NULL;
    copy = maf_copyMafBlock(orig);
    leftTrimmed = processBlockForTrim(copy, seq, start, stop, true);
    if (leftTrimmed == NULL) {
        maf_destroyMafBlockList(copy);
        return;
    }
    if (leftTrimmed != copy) {
        // i.e. we actually trimmed off some stuff
        maf_destroyMafBlockList(copy);
        copy = leftTrimmed;
    }
    mafBlock_t *leftPart = NULL, *rightPart = NULL;
    processBlockForSplit(copy, seq, start, stop, &leftPart, &rightPart);
    if (leftPart != copy) {
        maf_destroyMafBlockList(copy);
        copy = leftPart;
    }
    rightTrimmed = processBlockForTrim(copy, seq, start, stop, false);
    if (rightTrimmed != NULL)
        maf_mafBlock_print(rightTrimmed);
    if (rightTrimmed != copy) {
        maf_destroyMafBlockList(copy);
        copy = rightTrimmed;
    }
    maf_destroyMafBlockList(copy);
    if (rightPart != NULL) {
        // recurse for a spell if we're doing block splitting
        trimAndReportBlock(rightPart, seq, start, stop);
        maf_destroyMafBlockList(rightPart);
    }
}
void checkBlock(mafBlock_t *b, const char *seq, uint32_t start, uint32_t stop, bool *printedHeader, bool isSoft) {
    // read through each line of a mafBlock and if the sequence matches the region
    // we're looking for, report the block.
    mafLine_t *ml = maf_mafBlock_getHeadLine(b);
    mafBlock_t *dummy = NULL;
    while (ml != NULL) {
        if (searchMatched(ml, seq, start, stop)) {
            if (!*printedHeader) {
                printHeader();
                *printedHeader = true;
            }
            if (isSoft) {
                maf_mafBlock_print(b);
                break;
            } else {
                dummy = processBlockForSplice(b, seq, start, stop, false);
                assert(dummy == NULL);
                // trimAndReportBlock(b, seq, start, stop);
                break;
            }
        } 
        ml = maf_mafLine_getNext(ml);
    }
}
void processBody(mafFileApi_t *mfa, char *seq, uint32_t start, uint32_t stop, bool isSoft) {
    mafBlock_t *thisBlock = NULL;
    bool printedHeader = false;
    while ((thisBlock = maf_readBlock(mfa)) != NULL) {
        checkBlock(thisBlock, seq, start, stop, &printedHeader, isSoft);
        maf_destroyMafBlockList(thisBlock);
    }
}
