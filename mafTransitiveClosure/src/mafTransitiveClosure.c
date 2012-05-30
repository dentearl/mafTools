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
#include <ctype.h> // mac os x toupper()
#include <getopt.h>
#include <inttypes.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include "common.h"
#include "CuTest.h"
#include "sharedMaf.h"
#include "sonLib.h"
#include "stPinchGraphs.h"
#include "mafTransitiveClosure.h"
#include "test.mafTransitiveClosure.h"

const uint32_t kPinchThreshold = 50000000;
const char *kVersion = "v0.1 June 2012";

void parseOptions(int argc, char **argv, char *filename) {
    int c;
    int setMName = 0;
    while (1) {
        static struct option long_options[] = {
            {"debug", no_argument, 0, 'd'},
            {"verbose", no_argument, 0, 'v'},
            {"help", no_argument, 0, 'h'},
            {"test", no_argument, 0, 't'},
            {"maf",  required_argument, 0, 'm'},            
            {0, 0, 0, 0}
        };
        int option_index = 0;
        c = getopt_long(argc, argv, "d:m:h:v:t",
                        long_options, &option_index);
        if (c == -1)
            break;
        switch (c) {
        case 0:
            break;
        case 'm':
            setMName = 1;
            sscanf(optarg, "%s", filename);
            break;
        case 'v':
            g_verbose_flag++;
            break;
        case 'd':
            g_debug_flag = 1;
            break;
        case 't':
            exit(mafTransitiveClosure_RunAllTests());
            break;
        case 'h':
        case '?':
            usage();
            break;
        default:
            abort();
        }
    }
    if (!setMName) {
        fprintf(stderr, "specify --maf\n");
        usage();
    }
    // Check there's nothing left over on the command line 
    if (optind < argc) {
        char *errorString = de_malloc(kMaxStringLength);
        strcpy(errorString, "Unexpected arguments:");
        while (optind < argc) {
            strcat(errorString, " ");
            strcat(errorString, argv[optind++]);
        }
        fprintf(stderr, "%s\n", errorString);
        free(errorString);
        usage();
    }
}
void usage(void) {
    fprintf(stderr, "Usage: mafTransitiveClosure --maf mafFile.maf > transitivelyClosed.maf \n\n"
            "mafTransitiveClosure is a program to perform the transitive closure on\n"
            "an alignment. That is it checks every column of the alignment and looks\n"
            "for situations where a position A is aligned to B in one part of a file\n"
            "and B is aligned to C in another part of the file. The transitive closure\n"
            "of this relationship would be a single column with A, B and C all present.\n"
            "Useful for when you have pairwise alignments and you wish to turn them into\n"
            "something more resembling a multiple alignment."
            "\n\n");
    fprintf(stderr, "Options: \n"
            "  -h, --help     show this help message and exit.\n"
            "  -m, --maf      path to maf file.\n"
            "  -v, --verbose  turns on verbose output.\n");
    exit(EXIT_FAILURE);
}
uint32_t hashMafTcSeq(const mafTcSeq_t *mtcs) {
    return stHash_stringKey(mtcs->name);
}
int hashCompareMafTcSeq(const mafTcSeq_t *m1, const mafTcSeq_t *m2) {
    return strcmp(m1->name, m2->name) == 0;
}
char* createNSequence(uint32_t length) {
    char *seq = (char*) de_malloc(length + 1);
    memset(seq, 'N', length);
    seq[length] = '\0';
    return seq;
}
mafTcSeq_t* newMafTcSeq(char *name, uint32_t length) {
    mafTcSeq_t *mtcs = (mafTcSeq_t *) de_malloc(sizeof(*mtcs));
    mtcs->name = name;
    mtcs->sequence = createNSequence(length);
    mtcs->length = length;
    return mtcs;
}
int g_mafRegions = 0;
mafTcRegion_t* newMafTcRegion(uint32_t start, uint32_t end) {
    // de_debug("create, g_mafRegions = %d, [%u, %u]\n", ++g_mafRegions, start, end);
    mafTcRegion_t *reg = (mafTcRegion_t *) de_malloc(sizeof(*reg));
    reg->start = start;
    reg->end = end;
    reg->next = NULL;
    return reg;
}
mafTcComparisonOrder_t* newMafTcComparisonOrder(void) {
    mafTcComparisonOrder_t *co = (mafTcComparisonOrder_t *) de_malloc(sizeof(*co));
    co->ref = 0;
    co->region = NULL;
    co->next = NULL;
    return co;
}
void destroyMafTcRegionList(mafTcRegion_t *r) {
    mafTcRegion_t *tmp = NULL;
    while (r != NULL) {
        tmp = r;
        r = r->next;
        destroyMafTcRegion(tmp);
    }
}
void destroyMafTcRegion(mafTcRegion_t *r) {
    // de_debug("destroy, g_mafRegions = %d [%u, %u]\n", --g_mafRegions, r->start, r->end);
    free(r);
}
void destroyMafTcSeq(void *p) {
    free(((mafTcSeq_t *)p)->name);
    free(((mafTcSeq_t *)p)->sequence);
    free(p);
}
void destroyMafTcComparisonOrder(mafTcComparisonOrder_t *co) {
    mafTcComparisonOrder_t *tmp = NULL;
    while (co != NULL) {
        tmp = co;
        co = co->next;
        destroyMafTcRegion(tmp->region);
        free(tmp);
    }
}
void reverseComplementSequence(char *s) {
    // accepts upper and lower case, returns upper case 
   int c, i, j;
    for (i = 0, j = strlen(s) - 1; i < j; ++i, j--) {
        c = s[i];
        s[i] = s[j];
        s[j] = c;
    }
    complementSequence(s);
}
void complementSequence(char *s) {
    // accepts upper and lower case, returns upper case
    for (unsigned i = 0; i < strlen(s); ++i)
        s[i] = complementChar(s[i]);
}
char complementChar(char c) {
    switch (toupper(c)) {
    case 'A': 
        return 'T';
    case 'C':
        return 'G';
    case 'G':
        return 'C';
    case 'T':
        return 'A';
    case 'M':
        return 'K';
    case 'R':
        return 'Y';
    case 'W':
        return 'W';
    case 'S':
        return 'S';
    case 'Y':
        return 'R';
    case 'K':
        return 'M';
    case 'V':
        return 'B';
    case 'H':
        return 'D';
    case 'D':
        return 'H';
    case 'B':
        return 'V';
    case 'N':
    case '-':
    case 'X':
        return c;
    default:
        fprintf(stderr, "Error, unanticipated character in sequence: %c\n", c);
        exit(EXIT_FAILURE);
    }
}
void addSequenceValuesToMtcSeq(mafLine_t *ml, mafTcSeq_t *mtcs) {
    // add sequence values to maf transitive closure sequence
    int s; // transformed pos coordinate start (zero based)
    int e; // transformed pos coordinate stop (one based)
    if (maf_mafLine_getStrand(ml) == '+') {
        s = maf_mafLine_getStart(ml);
    } else {
        // THIS IS A DESTRUCTIVE OPERATION ON THE MAF LINE ml:
        reverseComplementSequence(maf_mafLine_getSequence(ml));
        s = maf_mafLine_getSourceLength(ml) - (maf_mafLine_getStart(ml) + maf_mafLine_getLength(ml));
    }
    e = s + maf_mafLine_getLength(ml);
    char *seq = maf_mafLine_getSequence(ml);
    for (unsigned i = 0, p = 0; i < strlen(maf_mafLine_getSequence(ml)); ++i) {
        // p is the current position coordinate within the sequence (zero based)
        if (seq[i] != '-') {
            if (mtcs->sequence[s + p] != 'N') {
                // sanity check
                if (toupper(mtcs->sequence[s + p]) != toupper(seq[i])) {
                    fprintf(stderr, "Error, maf file is inconsistent with regard to sequence. "
                            "On line number %" PRIu32 " sequence %s position %" PRIu32" is %c, but previously "
                            "observed value is %c.\n", maf_mafLine_getLineNumber(ml), maf_mafLine_getSpecies(ml), 
                            s + p, seq[i], mtcs->sequence[s + p]);
                    exit(EXIT_FAILURE);
                }
            }
            mtcs->sequence[s + p] = toupper(seq[i]);
            ++p;
        }
    }
}
void walkBlockAddingSequence(mafBlock_t *mb, stHash *hash, stHash *nameHash) {
    de_debug("walkBlockAddingSequence()\n");
    mafLine_t *ml = maf_mafBlock_getHeadLine(mb);
    mafTcSeq_t *mtcs = NULL;
    while ((ml = maf_mafLine_getNext(ml)) != NULL) {
        if (maf_mafLine_getType(ml) == 's') {
            if (stHash_search(hash, maf_mafLine_getSpecies(ml)) == NULL) {
                mtcs = newMafTcSeq(stString_copy(maf_mafLine_getSpecies(ml)), maf_mafLine_getSourceLength(ml));
                de_debug("Adding new sequence to hash: %s\n", maf_mafLine_getSpecies(ml));
                addSequenceValuesToMtcSeq(ml, mtcs);
                stHash_insert(hash, stString_copy(maf_mafLine_getSpecies(ml)), mtcs);
                if (stHash_search(nameHash, (void*)(int64_t)stHash_stringKey(maf_mafLine_getSpecies(ml))) == NULL) {
                    stHash_insert(nameHash, (void*)(int64_t)stHash_stringKey(maf_mafLine_getSpecies(ml)), 
                                  stString_copy(maf_mafLine_getSpecies(ml)));
                }
            } else {
                addSequenceValuesToMtcSeq(ml, stHash_search(hash, (mafTcSeq_t *)maf_mafLine_getSpecies(ml)));
            }
        }
    }
}
void createSequenceHash(mafFileApi_t *mfa, stHash **hash, stHash **nameHash) {
    de_debug("createSequenceHash()\n");
    *hash = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, free, destroyMafTcSeq);
    *nameHash = stHash_construct2(NULL, free);
    mafBlock_t *mb = NULL;
    while ((mb = maf_readBlock(mfa)) != NULL) {
        walkBlockAddingSequence(mb, *hash, *nameHash);
        maf_destroyMafBlockList(mb);
    }
}
void reportSequenceHash(stHash *hash, stHash *nameHash) {
    stHashIterator *hit = stHash_getIterator(hash);
    char *key = NULL;
    while ((key = stHash_getNext(hit)) != NULL) {
        printf("found key: %s: ", key);
        printf("value: %" PRIu32 "\n", ((mafTcSeq_t *)stHash_search(hash, key))->length);
        printf("   %s\n", ((mafTcSeq_t *)stHash_search(hash, key))->sequence);
    }
    stHash_destructIterator(hit);
    hit = stHash_getIterator(nameHash);
    key = NULL;
    while ((key = stHash_getNext(hit)) != NULL) {
        printf("found key: %" PRIi64 ": ", (int64_t) key);
        printf("value: %s\n", ((char *)stHash_search(nameHash, key)));
    }
    stHash_destructIterator(hit);   
}
stPinchThreadSet* buildThreadSet(stHash *hash) {
    stPinchThreadSet *ts = stPinchThreadSet_construct();
    stHashIterator *hit = stHash_getIterator(hash);
    char *key = NULL;
    while ((key = stHash_getNext(hit)) != NULL) {
        stPinchThreadSet_addThread(ts, stHash_stringKey(key), 0, ((mafTcSeq_t *)stHash_search(hash, key))->length);
    }
    stHash_destructIterator(hit);
    return ts;
}
mafTcComparisonOrder_t *getComparisonOrderFromMatrix(char **mat, uint32_t rowLength, uint32_t colLength) {
    mafTcRegion_t *todo = newMafTcRegion(0, colLength - 1); // the entire region needs to be done
    mafTcComparisonOrder_t *co = NULL;
    uint32_t r = 0;
    while (todo != NULL && r < rowLength) {
        todo = getComparisonOrderFromRow(mat, r++, &co, todo);
    }
    destroyMafTcRegionList(todo);
    return co;
}
mafTcRegion_t* getComparisonOrderFromRow(char **mat, uint32_t row, mafTcComparisonOrder_t **done, mafTcRegion_t *todo) {
    // 
    mafTcRegion_t *reg = NULL;
    mafTcRegion_t *newTodo = NULL;
    mafTcRegion_t *headTodo = todo;
    mafTcComparisonOrder_t *co = NULL;
    bool inGap = false;
    while (todo != NULL) {
        for (uint32_t i = todo->start; i <= todo->end; ++i) {
            de_debug("position %" PRIu32 "\n", i);
            if (mat[row][i] == '-') {
                // inside a gap region
                de_debug("inside a gap region\n");
                if (!inGap) {
                    de_debug("!inGap\n");
                    // transition from sequence into gap
                    inGap = true;
                    if (reg != NULL) {
                        de_debug("reg != NULL\n");
                        reg->next = newMafTcRegion(i, i);
                        reg = reg->next;
                    } else {
                        de_debug("reg == NULL\n");
                        reg = newMafTcRegion(i, i);
                    }
                    if (newTodo == NULL) {
                        de_debug("newTodo == NULL\n");
                        newTodo = reg;
                    }
                } else {
                    de_debug("inGap\n");
                    // remain in gap
                    reg->end = i;
                }
            } else {
                de_debug("inside a sequence region\n");
                // inside a sequence region
                if (inGap) {
                    de_debug("inGap\n");
                    // transition from gap to sequence
                    inGap = false;
                    co = newMafTcComparisonOrder();
                    co->ref = row;
                    co->region = newMafTcRegion(i, i);
                    de_debug("creating new region row:%" PRIu32 " at %" PRIu32 "\n", row, i);
                    // insert into the start of the list
                    if (done != NULL) {
                        de_debug("done != NULL\n");
                        co->next = *done;
                    }
                    *done = co;
                } else {
                    de_debug("remain in sequence\n");
                    // remain in sequence
                    if (co != NULL) {
                        de_debug("co != NULL, moving ref:%" PRIu32 " end to %" PRIu32 "\n", co->ref, i);
                        co->region->end = i;
                    } else {
                        de_debug("co == NULL, ref: %" PRIu32 "\n", row);
                        co = newMafTcComparisonOrder();
                        co->ref = row;
                        co->region = newMafTcRegion(i, i);
                        // insert into the start of the `done' list
                        if (done != NULL) {
                            de_debug("done != NULL\n");
                            co->next = *done;
                        } else {
                            de_debug("done == NULL\n");
                            done = (mafTcComparisonOrder_t **) de_malloc(sizeof(*done));
                        }
                        *done = co;
                    }
                }
            }
        }
        // de_debug("advance todo:%p to %p\n", (void*)todo, (void*)todo->next);
        todo = todo->next;
    }
    destroyMafTcRegionList(headTodo);
    // de_debug("done: %p, ->ref: %d ->region: [%d, %d]\n", 
    //          (void*)(*done), (*done)->ref, (*done)->region->start, (*done)->region->end);
    // de_debug("newTodo: %p\n", (void*)newTodo);
    return newTodo;
}
int64_t localSeqCoordsToGlobalPositiveStartCoords(int64_t c, uint32_t start, uint32_t sourceLength, 
                                                  char strand, uint32_t length) {
    // given a coordinate inside of a sequence with respect to the start (which is 0), 
    // transform the START coordinate of a BLOCK to positive full source coordinates.
    if (strand == '+') {
        return localSeqCoordsToGlobalPositiveCoords(c, start, sourceLength, strand);
    } else {
        return localSeqCoordsToGlobalPositiveCoords(c, start + (length - 1), sourceLength, strand);
    }
}
int64_t localSeqCoordsToGlobalPositiveCoords(int64_t c, uint32_t start, uint32_t sourceLength, char strand) {
    // given a coordinate inside of a sequence with respect to the start (which is 0), 
    // transform the coordinate to positive full source coordinates.
    if (strand == '+') {
        return (int64_t) (start + c);
    } else {
        return (int64_t) (sourceLength - 1 - start - c);
    }
}
int64_t localSeqCoords(uint32_t p, char *s) {
    // given a coordinate inside of a sequence with respect to the start (which is 0), 
    // walk backwards and anytime you see a gap drop one from the coordinate to return
    int64_t value = (int64_t) p; // coordinate to return
    for (int64_t i = p; i >= 0; --i) {
        if (s[i] == '-') {
            --value;
        }
    }
    return value;
}
uint32_t g_numPinches = 0;
void processPairForPinching(stPinchThreadSet *threadSet, stPinchThread *a, uint32_t aGlobalStart, 
                            uint32_t aGlobalLength, int aStrand, 
                            char *aSeq, stPinchThread *b, uint32_t bGlobalStart, uint32_t bGlobalLength,
                            int bStrand, char *bSeq, uint32_t regionStart, uint32_t regionEnd) {
    // perform a pinch operation for regions of bSeq that are not gaps, i.e. `-'
    // GlobalStart is the positive strand position (zero based) coordinate of the start of this block
    de_debug("processPairForPinching() g_numPinches: %u\n", g_numPinches);
    uint32_t l = 0, pos = 0, regionLength = regionEnd - regionStart + 1;
    uint32_t blockStart = pos;
    de_debug("aGlobalStart: %" PRIu32 ", bGlobalStart: %" PRIu32 ", pos & regionStart: %" PRIu32 ", regionLength: %" PRIu32 "\n", 
             aGlobalStart, bGlobalStart, regionStart, regionLength);
    de_debug("aSeq: %s\n", aSeq);
    de_debug("bSeq: %s\n", bSeq);
    bool inBlock = false;
    for (pos = regionStart; pos < regionStart + regionLength; ++pos) {
        de_debug("pos: %" PRIu32 " \n", pos);
        if (bSeq[pos] == '-') {
            if (inBlock) {
                de_debug("bSeq[pos] == - && inBlock\n");
                inBlock = false;
                de_debug("maybe little pinch (1)? p:%" PRIu32 ", l:%" PRIu32 "\n", blockStart, l);
                de_debug("a coords local: %" PRIi64 " global: %" PRIi64 
                         ", b coords local: %" PRIi64 ", global: %" PRIi64 "\n", 
                         localSeqCoords(blockStart, aSeq),
                         localSeqCoordsToGlobalPositiveStartCoords(localSeqCoords(blockStart, aSeq),
                                                              aGlobalStart, aGlobalLength, aStrand, l),
                         localSeqCoords(blockStart, bSeq),
                         localSeqCoordsToGlobalPositiveStartCoords(localSeqCoords(blockStart, bSeq),
                                                                   bGlobalStart, bGlobalLength, bStrand, l));
                stPinchThread_pinch(a, b,
                                    localSeqCoordsToGlobalPositiveStartCoords(localSeqCoords(blockStart, aSeq), 
                                                                              aGlobalStart,
                                                                              aGlobalLength,
                                                                              aStrand, l),
                                    localSeqCoordsToGlobalPositiveStartCoords(localSeqCoords(blockStart, bSeq),
                                                                              bGlobalStart,
                                                                              bGlobalLength,
                                                                              bStrand,l),
                                    l, (aStrand == bStrand));
                ++g_numPinches;
                if (g_numPinches > kPinchThreshold) {
                    stPinchThreadSet_joinTrivialBoundaries(threadSet);
                    g_numPinches = 0;
                }
                l = 0;
            }
        } else {
            if (!inBlock) {
                de_debug("bSeq[pos] != - && !inBlock\n");
                inBlock = true;
                blockStart = pos;
            }
            ++l; // length of this alignment
        }
    }
    de_debug("done walking comparison region\n");
    if (inBlock) {
        inBlock = false;
        de_debug("maybe little pinch (2)? p:%" PRIu32 ", l:%" PRIu32 "\n", blockStart, l);
        de_debug("a coords local: %" PRIi64 " global: %" PRIi64 
                         ", b coords local: %" PRIi64 ", global: %" PRIi64 "\n", 
                         localSeqCoords(blockStart, aSeq),
                         localSeqCoordsToGlobalPositiveStartCoords(localSeqCoords(blockStart, aSeq),
                                                                   aGlobalStart, aGlobalLength, aStrand, l),
                         localSeqCoords(blockStart, bSeq),
                         localSeqCoordsToGlobalPositiveStartCoords(localSeqCoords(blockStart, bSeq),
                                                                   bGlobalStart, bGlobalLength, bStrand, l));
        stPinchThread_pinch(a, b,
                            localSeqCoordsToGlobalPositiveStartCoords(localSeqCoords(blockStart, aSeq), 
                                                                      aGlobalStart,
                                                                      aGlobalLength,
                                                                      aStrand, l),
                            localSeqCoordsToGlobalPositiveStartCoords(localSeqCoords(blockStart, bSeq),
                                                                      bGlobalStart,
                                                                      bGlobalLength,
                                                                      bStrand, l),
                            l, (aStrand == bStrand));
        ++g_numPinches;
    }
    if (g_numPinches > kPinchThreshold) {
        stPinchThreadSet_joinTrivialBoundaries(threadSet);
        g_numPinches = 0;
    }
    de_debug("exiting processPairForPinching()...\n");
}
static void printMatrix(char **mat, uint32_t n) {
    fprintf(stderr, "printMatrix()\n");
    for (uint32_t i = 0; i < n; ++i)
        fprintf(stderr, "%" PRIu32 ": %s\n", i, mat[i]);
    fprintf(stderr, "\n");
}
static void printu32Array(uint32_t *a, uint32_t n) {
    for (uint32_t i = 0; i < n; ++i)
        fprintf(stderr, "%" PRIu32 "%s", a[i], (i == n - 1) ? "\n" : ", ");
}
static void printComparisonOrder(mafTcComparisonOrder_t *co) {
    de_debug("printComparisonOrder()\n");
    while (co != NULL) {
        de_debug("{%" PRIu32 ": [%" PRIu32 ", %" PRIu32 "]}%s", co->ref, co->region->start, co->region->end,
                 (co->next == NULL) ? "\n" : ", ");
        co = co->next;
    }
}
void walkBlockAddingAlignments(mafBlock_t *mb, stPinchThreadSet *threadSet) {
    // for a given block, add the alignment information to the threadset.
    de_debug("walkBlockAddingAlignments()\n");
    if (!maf_mafBlock_containsSequence(mb))
        return;
    uint32_t numSeqs = maf_mafBlock_getNumberOfSequences(mb);
    uint32_t seqFieldLength = maf_mafBlock_longestSequenceField(mb);
    char **mat = maf_mafBlock_getSequenceMatrix(mb, numSeqs, seqFieldLength);
    char *strands = maf_mafBlock_getStrandArray(mb);
    char **names = maf_mafBlock_getSpeciesMatrix(mb);
    uint32_t *starts = maf_mafBlock_getStartArray(mb);
    uint32_t *lengths = maf_mafBlock_getSourceLengthArray(mb);
    // comparison order coordinates are relative to the block
    mafTcComparisonOrder_t *co = getComparisonOrderFromMatrix(mat, numSeqs, seqFieldLength);
    printComparisonOrder(co);
    mafTcComparisonOrder_t *c = co;
    stPinchThread *a = NULL, *b = NULL;
    while (c != NULL) {
        de_debug("c != NULL, c->ref=%" PRIu32 ", c->region->start=%" PRIu32 " c->region->end=%" PRIu32 "\n", 
                 c->ref, c->region->start, c->region->end);
        a = stPinchThreadSet_getThread(threadSet, stHash_stringKey(names[c->ref]));
        assert(a != NULL);
        for (uint32_t r = c->ref + 1; r < numSeqs; ++r) {
            b = stPinchThreadSet_getThread(threadSet, stHash_stringKey(names[r]));
            assert(b != NULL);
            de_debug("going to pinch ref %" PRIu32 ":%s to %" PRIu32 ":%s\n", c->ref, names[c->ref], r, names[r]);
            de_debug("a %" PRIu32 ": %s\n", c->ref, mat[c->ref]);
            de_debug("b %" PRIu32 ": %s\n", r, mat[r]);
            processPairForPinching(threadSet, a, starts[c->ref], lengths[c->ref], strands[c->ref],
                                   mat[c->ref], b, starts[r], lengths[r], strands[r], mat[r], c->region->start,
                                   c->region->end);
        }
        c = c->next;
    }
    
    // cleanup
    destroyMafTcComparisonOrder(co);
    maf_mafBlock_destroySequenceMatrix(mat, maf_mafBlock_getNumberOfSequences(mb));
    free(strands);
    free(starts);
    free(lengths);
    for (uint32_t i = 0; i < numSeqs; ++i) {
        free(names[i]);
    }
    free(names);
}
void addAlignmentsToThreadSet(mafFileApi_t *mfa, stPinchThreadSet *threadSet) {
    mafBlock_t *mb = NULL;
    while ((mb = maf_readBlock(mfa)) != NULL) {
        de_debug("addAlignmentsToThreadSet(), read a block\n");
        walkBlockAddingAlignments(mb, threadSet);
        maf_destroyMafBlockList(mb);
    }
    stPinchThreadSet_joinTrivialBoundaries(threadSet);
}
uint32_t getMaxNameLength(stHash *hash) {
    stHashIterator *hit = stHash_getIterator(hash);
    char *key = NULL;
    uint32_t max = 0;
    while ((key = stHash_getNext(hit)) != NULL) {
        if (max < strlen(key) + 2)
            max = strlen(key) + 2;
    }
    stHash_destructIterator(hit);
    return max;
}
void getMaxFieldLengths(stHash *hash, stHash *nameHash, stPinchBlock *block, uint32_t *maxStart, 
                        uint32_t *maxLength, uint32_t *maxSource) {
    stPinchBlockIt thisSegIt = stPinchBlock_getSegmentIterator(block);
    stPinchSegment *thisSeg = NULL;
    *maxStart = 0;
    *maxLength = 0;
    *maxSource = 0;
    char *temp = NULL;
    char *key = NULL;
    char strand = '\0';
    while ((thisSeg = stPinchBlockIt_getNext(&thisSegIt)) != NULL) {
        strand = stPinchSegment_getBlockOrientation(thisSeg) == 1 ? '+' : '-';
        key = (char*)stHash_search(nameHash, (void *)stPinchSegment_getName(thisSeg));
        temp = (char*) de_malloc(kMaxStringLength);
        if (strand == '+') {
            sprintf(temp, "%" PRIi64, stPinchSegment_getStart(thisSeg));
        } else {
            sprintf(temp, "%" PRIi64,
                    (((int64_t)((mafTcSeq_t*)stHash_search(hash, key))->length) - 
                     stPinchSegment_getStart(thisSeg) - stPinchSegment_getLength(thisSeg)));
        }
        if (*maxStart < strlen(temp))
            *maxStart = strlen(temp);

        free(temp);
        temp = (char*) de_malloc(kMaxStringLength);
        sprintf(temp, "%" PRIi64, stPinchSegment_getLength(thisSeg));
        if (*maxLength < strlen(temp))
            *maxLength = strlen(temp);
        free(temp);
        temp = (char*) de_malloc(kMaxStringLength);
        sprintf(temp, "%" PRIu32, ((mafTcSeq_t*)stHash_search(hash, key))->length);
        if (*maxSource < strlen(temp))
            *maxSource = strlen(temp);
        free(temp);
    }
}
char* getSequenceSubset(char *seq, int64_t start, char strand, int64_t length) {
    char *out = (char*) de_malloc(sizeof(*out) * length + 1);
    int64_t i, j;
    for (i = start, j = 0; i < start + length; ++i, ++j) {
        out[j] = seq[i];
    }
    out[j] = '\0';
    if (strand == '-')
        reverseComplementSequence(out);
    return out;
}
void reportTransitiveClosure(stPinchThreadSet *threadSet, stHash *hash, stHash *nameHash) {
    stPinchThreadSetBlockIt thisBlockIt = stPinchThreadSet_getBlockIt(threadSet);
    stPinchBlock *thisBlock = NULL;
    stPinchBlockIt thisSegIt;
    stPinchSegment *thisSeg = NULL;
    char *key = NULL;
    char *seq = NULL;
    char strand = '\0';
    (void) hash;
    (void) nameHash;
    printf("##maf version=1\n");
    printf("# mafTransitiveClosure %s\n\n", kVersion);
    uint32_t maxNameLength, maxStartLength, maxLengthLength, maxSourceLengthLength;
    int64_t xformedStart;
    maxNameLength = getMaxNameLength(hash);
    while ((thisBlock = stPinchThreadSetBlockIt_getNext(&thisBlockIt)) != NULL) {
        getMaxFieldLengths(hash, nameHash, thisBlock, &maxStartLength, 
                           &maxLengthLength, &maxSourceLengthLength);
        printf("a degree=%" PRIu32 "\n", stPinchBlock_getDegree(thisBlock));
        thisSegIt = stPinchBlock_getSegmentIterator(thisBlock);
        while ((thisSeg = stPinchBlockIt_getNext(&thisSegIt)) != NULL) {
            strand = stPinchSegment_getBlockOrientation(thisSeg) == 1 ? '+' : '-';
            key = (char*)stHash_search(nameHash, (void *)stPinchSegment_getName(thisSeg));
            seq = getSequenceSubset(((mafTcSeq_t*)stHash_search(hash, key))->sequence,
                                    stPinchSegment_getStart(thisSeg),
                                    strand,
                                    stPinchSegment_getLength(thisSeg));
            if (strand == '+') {
                xformedStart = stPinchSegment_getStart(thisSeg);
            } else {
                xformedStart = (((int64_t)((mafTcSeq_t*)stHash_search(hash, key))->length) - 
                                stPinchSegment_getStart(thisSeg) - stPinchSegment_getLength(thisSeg));
            }
            printf("s %-*s %*" PRIi64 " %*" PRIi64 " %c %*" PRIu32 " %s\n",
                   maxNameLength, key,
                   maxStartLength, xformedStart,
                   maxLengthLength, stPinchSegment_getLength(thisSeg),
                   strand,
                   maxSourceLengthLength, ((mafTcSeq_t*)stHash_search(hash, key))->length,
                   seq
                   );
            free(seq);
        }
        printf("\n");
    }
    printf("\n");
} 
int main(int argc, char **argv) {
    (void) (printMatrix);
    (void) (printu32Array);
    (void) (reportSequenceHash);
    char filename[kMaxStringLength];
    stHash *sequenceHash, *nameHash;
    parseOptions(argc, argv, filename);
    
    // first pass, build sequence hash
    mafFileApi_t *mfa = maf_newMfa(filename, "r");
    createSequenceHash(mfa, &sequenceHash, &nameHash);
    stPinchThreadSet *threadSet = buildThreadSet(sequenceHash);
    maf_destroyMfa(mfa);
    
    // second pass, build pinch graph
    mfa = maf_newMfa(filename, "r");
    addAlignmentsToThreadSet(mfa, threadSet);
    maf_destroyMfa(mfa);
    
    // consolidate and report
    reportTransitiveClosure(threadSet, sequenceHash, nameHash);
    
    // cleanup
    stHash_destruct(sequenceHash);
    stHash_destruct(nameHash);
    stPinchThreadSet_destruct(threadSet);

    return EXIT_SUCCESS;
}
