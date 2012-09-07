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
void getTargetColumns(bool **targetColumns, uint32_t *len, mafBlock_t *b, const char *seqName,
                      uint32_t start, uint32_t stop) {
    // given a block and a target, create an array of bools where true means the target
    // is present in that column and false means it is absent
    mafLine_t *ml = maf_mafBlock_getHeadLine(b);
    while(maf_mafLine_getType(ml) != 's') {
        ml = maf_mafLine_getNext(ml);
    }
    char *seq = NULL;
    *len = maf_mafBlock_getSequenceFieldLength(b);
    // printf("target columns len: %" PRIu32 "\n", *len);
    if (*len == 0) {
        *targetColumns = NULL;
        return;
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
                (*targetColumns)[i] = (*targetColumns)[i] | ((start <= pos) && (pos <= stop));
                pos += it;
            }
        }
        ml = maf_mafLine_getNext(ml);
    }
    // printf("I'm out, son: len=%" PRIu32 "\n", *len);
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
    uint32_t len = 0, n = 0;
    // printf("processBlockForTrim(), blocklen: %" PRIu32 "\n", maf_mafBlock_getSequenceFieldLength(b));
    getTargetColumns(&targetColumns, &len, b, seq, start, stop);
    // printf("back from getTargetColumns(): len=%" PRIu32 "\n", len);
    // printBoolArray(targetColumns, len);
    if (sumBool(targetColumns, len) == 0) {
        free(targetColumns);
        return NULL;
    }
    if (sumBool(targetColumns, len) == len) {
        // bail out, nothing to trim here.
        // printf("bailing out, sumBool:%"PRIu32" == len %"PRIu32"\n", sumBool(targetColumns, len), len);
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
    // printf("processBlockForTrim(), n:%" PRIu32 "\n", n);
    if (n == 0) {
        // nothing to trim, return a copy of the original
        free(targetColumns);
        return maf_copyMafBlock(b);
    } else if (n == maf_mafBlock_getSequenceFieldLength(b)) {
        // this would be the case where every single element in targetColumns is 0,
        // which would mean this block contains nothing relevant.
        return NULL;
    } else {
        // need to trim off `n' bases from the left / right
        // where n is greater than 0 and less than the total length of the sequence
        free(targetColumns);
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
    uint32_t linenumber = maf_mafLine_getLineNumber(ml1);
    ml1 = maf_mafLine_getNext(ml1);
    maf_mafBlock_setHeadLine(mb, maf_newMafLineFromString("a score=0", linenumber));
    ml2 = maf_mafBlock_getHeadLine(mb);
    maf_mafBlock_setLineNumber(mb, linenumber);
    maf_mafBlock_incrementNumberOfLines(mb);
    uint32_t m, len;
    char *seq = NULL;
    bool prevLineUsed = true; // used when a mafline is dropped because it's length becomes 0
    bool emptyBlock = true;
    while (ml1 != NULL) {
        // loop through all maf lines in a block
        if (prevLineUsed) {
            maf_mafLine_setNext(ml2, maf_newMafLine());
            ml2 = maf_mafLine_getNext(ml2);
        }
        if (maf_mafLine_getLine(ml1) == NULL) {
            ml1 = maf_mafLine_getNext(ml1);
            prevLineUsed = false;
            continue;
        }
        maf_mafLine_setType(ml2, maf_mafLine_getType(ml1));
        maf_mafLine_setLineNumber(ml2, ++linenumber);
        maf_mafBlock_incrementNumberOfLines(mb);
        maf_mafBlock_incrementLineNumber(mb);
        if (maf_mafLine_getType(ml2) != 's') {
            maf_mafLine_setLine(ml2, de_malloc(strlen(maf_mafLine_getLine(ml1)) + 1));
            strcpy(maf_mafLine_getLine(ml2), maf_mafLine_getLine(ml1));
            ml1 = maf_mafLine_getNext(ml1);
            prevLineUsed = true;
            continue;
        }
        maf_mafBlock_incrementNumberOfSequences(mb);
        seq = maf_mafLine_getSequence(ml1);
        
        m = 0;
        len = strlen(seq);
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
            --linenumber;
            continue;
        }
        if (isLeft) {
            assert(maf_mafLine_getStart(ml1) + m < maf_mafLine_getSourceLength(ml1));
            maf_mafLine_setStart(ml2, maf_mafLine_getStart(ml1) + m);
            maf_mafLine_setSequence(ml2, de_strdup(seq + n));
        } else {
            maf_mafLine_setStart(ml2, maf_mafLine_getStart(ml1));
            size_t len = strlen(seq);
            if (len <= n) {
                fprintf(stderr, "Well, this is bullshit.\n");
                fprintf(stderr, "Trying to trim %" PRIu32 " off of a sequence of length %zu\n", n, len);
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
    maf_mafBlock_setTailLine(mb, ml2);
    if (!emptyBlock) {
        return mb;
    } else {
        // this condition should never be tripped, we have a condition set up
        // in the "process" wrapper above this function.
        printf("block was empty, this should not happen\n");
        assert(2 + 2 == 5);
        maf_destroyMafBlockList(mb);
        return NULL;
    }
}
void processBlockForSplit(mafBlock_t *b, const char *seq, uint32_t start, uint32_t stop, 
                          mafBlock_t **left, mafBlock_t **right) {
    bool *targetColumns = NULL;
    uint32_t len = 0, n = 0;
    // printf("processBlockForSplit()\n");
    getTargetColumns(&targetColumns, &len, b, seq, start, stop);
    // printf("a different place for target columns...\n");
    // printBoolArray(targetColumns, len);
    for (uint32_t i = 0; i < len; ++i) {
        if (targetColumns[i]) {
            ++n;
        } else {
            break;
        }
    }
    if (sumBool(targetColumns, len) == len) {
        // bail out, nothing to trim here.
        *left = b;
        *right = NULL;
    } else {
        // split the block in twain
        splitBlock(b, n, left, right);
    }
    free(targetColumns);
}
void splitBlock(mafBlock_t *b, uint32_t n, mafBlock_t **mbLeft, mafBlock_t **mbRight) {
    /* n is the number of bases to be included in the left block
     * 
     */
    // printf("splitBlock()\n");
    assert(*mbLeft == NULL);
    assert(*mbRight == NULL);
    uint32_t m = maf_mafBlock_getSequenceFieldLength(b);
    if ((n > 0) && (n < m - 1)) {
        // printf("trimming stuff? n:%"PRIu32" m:%" PRIu32"\n", n, m);
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
        maf_destroyMafBlockList(copy);
        copy = leftTrimmed;
    }
    mafBlock_t *left = NULL, *right = NULL;
    processBlockForSplit(copy, seq, start, stop, &left, &right);
    if (left != copy) {
        maf_destroyMafBlockList(copy);
        copy = left;
    }
    if (right != NULL) {
        // recurse for a spell if we're doing block splitting
        trimAndReportBlock(right, seq, start, stop);
        maf_destroyMafBlockList(right);
    }
    rightTrimmed = processBlockForTrim(copy, seq, start, stop, false);
    if (rightTrimmed != NULL)
        maf_mafBlock_print(rightTrimmed);
    if (rightTrimmed != copy) {
        maf_destroyMafBlockList(copy);
        copy = rightTrimmed;
    }
    maf_destroyMafBlockList(copy);
}
void quickPrintBlock(mafBlock_t *b) {
    mafLine_t *ml = maf_mafBlock_getHeadLine(b);
    while (ml != NULL) {
        printf("%s\n", maf_mafLine_getLine(ml));
        ml = maf_mafLine_getNext(ml);
    }
    printf("\n");
}
void checkBlock(mafBlock_t *b, const char *seq, uint32_t start, uint32_t stop, bool *printedHeader, bool isSoft) {
    // read through each line of a mafBlock and if the sequence matches the region
    // we're looking for, report the block.
    mafLine_t *ml = maf_mafBlock_getHeadLine(b);
    while (ml != NULL) {
        if (searchMatched(ml, seq, start, stop)) {
            if (!*printedHeader) {
                printHeader();
                *printedHeader = true;
            }
            if (isSoft) {
                quickPrintBlock(b);
                break;
            } else {
                trimAndReportBlock(b, seq, start, stop);
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
