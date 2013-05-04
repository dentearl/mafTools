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
#include "buildVersion.h"

const uint64_t kPinchThreshold = 50000000;
const char *g_version = "v0.2 May 2013";
bool g_isSort = false;

void version(void);
void usage(void);

void parseOptions(int argc, char **argv, char *filename) {
    int c;
    int setMName = 0;
    while (1) {
        static struct option long_options[] = {
            {"debug", no_argument, 0, 'd'},
            {"verbose", no_argument, 0, 'v'},
            {"help", no_argument, 0, 'h'},
            {"version", no_argument, 0, 0},
            {"test", no_argument, 0, 't'},
            {"maf",  required_argument, 0, 'm'},
            {"sort", no_argument, 0, 's'},
            {0, 0, 0, 0}
        };
        int option_index = 0;
        c = getopt_long(argc, argv, "d:m:s:h:v:t", long_options, &option_index);
        if (c == -1)
            break;
        switch (c) {
        case 0:
            if (strcmp("version", long_options[option_index].name) == 0) {
                version();
                exit(EXIT_SUCCESS);
            }
            break;
        case 'm':
            setMName = 1;
            sscanf(optarg, "%s", filename);
            break;
        case 's':
            g_isSort = true;
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
void version(void) {
    fprintf(stderr, "mafTransitiveClosure, %s\nbuild: %s, %s, %s\n\n", g_version, g_build_date, 
            g_build_git_branch, g_build_git_sha);
}
void usage(void) {
    version();
    fprintf(stderr, "Usage: mafTransitiveClosure --maf mafFile.maf > transitivelyClosed.maf \n\n"
            "mafTransitiveClosure is a program to perform the transitive closure on\n"
            "an alignment. That is it checks every column of the alignment and looks\n"
            "for situations where a position A is aligned to B in one part of a file\n"
            "and B is aligned to C in another part of the file. The transitive closure\n"
            "of this relationship would be a single column with A, B and C all present.\n"
            "Useful for when you have pairwise alignments and you wish to turn them into\n"
            "something more resembling a multiple alignment."
            "\n\n");
    fprintf(stderr, "Options: \n");
    usageMessage('h', "help", "show this message and exit.");
    usageMessage('m', "maf", "path to the maf file.");
    usageMessage('v', "verbose", "turns on verbose output..");
    exit(EXIT_FAILURE);
}
uint64_t hashMafTcSeq(const mafTcSeq_t *mtcs) {
    return stHash_stringKey(mtcs->name);
}
int hashCompareMafTcSeq(const mafTcSeq_t *m1, const mafTcSeq_t *m2) {
    return strcmp(m1->name, m2->name) == 0;
}
char* createNSequence(uint64_t length) {
    char *seq = (char*) de_malloc(length + 1);
    memset(seq, 'N', length);
    seq[length] = '\0';
    return seq;
}
mafTcSeq_t* newMafTcSeq(char *name, uint64_t length) {
    mafTcSeq_t *mtcs = (mafTcSeq_t *) de_malloc(sizeof(*mtcs));
    mtcs->name = name;
    mtcs->sequence = createNSequence(length);
    mtcs->length = length;
    return mtcs;
}
int g_mafRegions = 0;
mafTcRegion_t* newMafTcRegion(uint64_t start, uint64_t end) {
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
static void printRegion(mafTcRegion_t *reg) {
    while (reg != NULL) {
        de_debug("  s: %3" PRIu64 " e: %3" PRIu64 "\n", reg->start, reg->end);
        reg = reg->next;
    }
}
static void printComparisonOrder(mafTcComparisonOrder_t *co) {
    de_debug("printComparisonOrder()\n");
    while (co != NULL) {
        de_debug(" ref: %2" PRIu64 " \n", co->ref);
        printRegion(co->region);
        co = co->next;
    }
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
    free(r);
}
void destroyMafTcSeq(void *p) {
    // the extra casting here is due to the fact that this is called by the stHash destructor
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
void addSequenceValuesToMtcSeq(mafLine_t *ml, mafTcSeq_t *mtcs) {
    // add sequence values to maf transitive closure sequence
    int64_t s; // transformed pos coordinate start (zero based)
    if (maf_mafLine_getStrand(ml) == '+') {
        s = maf_mafLine_getStart(ml);
    } else {
        // THIS IS A DESTRUCTIVE OPERATION ON THE MAF LINE ml:
        reverseComplementSequence(maf_mafLine_getSequence(ml), maf_mafLine_getSequenceFieldLength(ml));
        s = maf_mafLine_getSourceLength(ml) - (maf_mafLine_getStart(ml) + maf_mafLine_getLength(ml));
    }
    char *seq = maf_mafLine_getSequence(ml);
    uint64_t n = maf_mafLine_getSequenceFieldLength(ml);
    for (uint64_t i = 0, p = 0; i < n; ++i) {
        // p is the current position coordinate within the sequence (zero based)
        if (seq[i] != '-') {
            if (mtcs->sequence[s + p] != 'N') {
                // sanity check
                if (toupper(mtcs->sequence[s + p]) != toupper(seq[i])) {
                    fprintf(stderr, "Error, maf file is inconsistent with regard to sequence. "
                            "On line number %" PRIu64 " sequence %s position %" PRIu64" is %c, but previously "
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
    mafLine_t *ml = maf_mafBlock_getHeadLine(mb);
    mafTcSeq_t *mtcs = NULL;
    char *name = NULL;
    while ((ml = maf_mafLine_getNext(ml)) != NULL) {
        if (maf_mafLine_getType(ml) == 's') {
            name = maf_mafLine_getSpecies(ml);
            if (stHash_search(hash, name) == NULL) {
                mtcs = newMafTcSeq(stString_copy(name), maf_mafLine_getSourceLength(ml));
                addSequenceValuesToMtcSeq(ml, mtcs);
                stHash_insert(hash, stString_copy(name), mtcs);
                int64_t *key = (int64_t *)st_malloc(sizeof(*key));
                *key = (int64_t) stHash_stringKey(name);
                if (stHash_search(nameHash, key) == NULL) {
                    stHash_insert(nameHash, key, stString_copy(name));
                } else {
                    free(key);
                }
            } else {
                addSequenceValuesToMtcSeq(ml, stHash_search(hash, name));
            }
        }
    }
}
static uint64_t uint64Return(const void *key) {
    return *((int64_t*) key);
}
static int int64EqualKey(const void *key1, const void *key2) {
    return *((int64_t *) key1) == *((int64_t*) key2);
}
void createSequenceHash(mafFileApi_t *mfa, stHash **hash, stHash **nameHash) {
    *hash = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, free, destroyMafTcSeq);
    *nameHash = stHash_construct3(uint64Return, int64EqualKey, free, free);
    mafBlock_t *mb = NULL;
    while ((mb = maf_readBlock(mfa)) != NULL) {
        walkBlockAddingSequence(mb, *hash, *nameHash);
        maf_destroyMafBlockList(mb);
    }
}
void reportSequenceHash(stHash *hash, stHash *nameHash) {
    stHashIterator *hit = stHash_getIterator(hash);
    char *key = NULL;
    printf("Sequence Hash:\n");
    while ((key = stHash_getNext(hit)) != NULL) {
        printf("found key: %s: ", key);
        printf("length: %" PRIu64 "\n", ((mafTcSeq_t *)stHash_search(hash, key))->length);
        printf("   %s\n", ((mafTcSeq_t *)stHash_search(hash, key))->sequence);
    }
    stHash_destructIterator(hit);
    hit = stHash_getIterator(nameHash);
    key = NULL;
    printf("Name Hash:\n");
    while ((key = stHash_getNext(hit)) != NULL) {
        printf("found key: %" PRIi64 ": ", *((int64_t *)key));
        printf("value: %s\n", ((char *)stHash_search(nameHash, key)));
    }
    stHash_destructIterator(hit);
}
stPinchThreadSet* buildThreadSet(stHash *hash) {
    stPinchThreadSet *ts = stPinchThreadSet_construct();
    stHashIterator *hit = stHash_getIterator(hash);
    char *key = NULL;
    while ((key = stHash_getNext(hit)) != NULL) {
        stPinchThreadSet_addThread(ts, stHash_stringKey(key), 0, 
                                   ((mafTcSeq_t *)stHash_search(hash, key))->length);
    }
    stHash_destructIterator(hit);
    return ts;
}
mafTcComparisonOrder_t *getComparisonOrderFromMatrix(char **mat, uint64_t numRows, uint64_t numCols, 
                                                     uint64_t *lengths, int **vizMat) {
    /* given a char matrix and its dimensions (and a debugging int visualization matrix) generate
       a linked list of mafTcComparisonOrder_t that contains gapless sequences that produce coverage
       of the entire matrix. Example:
       ATGT---AGT--A
       AG-GTG-AGT-TA
       ATGTGGTAGTTTA
       becomes, essentially:
       ****---***--*
       ..-.**-...-*.
       ......*...*..
       where * represent the comparison "reference" and . represent some non-gap character.
    */
    mafTcRegion_t *todo = newMafTcRegion(0, numCols - 1); // the entire region needs to be done
    mafTcComparisonOrder_t *co = NULL;
    uint64_t r = 0;
    if (g_debug_flag) {
        printTodoArray(todo, numCols);
        printVizMatrix(vizMat, numRows, numCols);
    }
    while (todo != NULL && r < numRows) {
        todo = getComparisonOrderFromRow(mat, r, &co, todo, (lengths[r] != numCols));
        r++;
        if (g_debug_flag) {            
            updateVizMatrix(vizMat, co);
            printTodoArray(todo, numCols);
            printVizMatrix(vizMat, numRows, numCols);
        }
    }
    destroyMafTcRegionList(todo);
    return co;
}
mafTcRegion_t* getComparisonOrderFromRow(char **mat, uint64_t row, mafTcComparisonOrder_t **done, 
                                         mafTcRegion_t *todo, int containsGaps) {
    // proudce a comparison order given a sequence matrix (mat), a row index (row), a list of already
    // completed regions (done) and a list of regions that still need to be done (todo)
    mafTcRegion_t *newTodo = NULL, *newTodoTail = NULL;
    mafTcRegion_t *headTodo = todo;
    mafTcComparisonOrder_t *co = NULL;
    bool inGap;
    de_debug("getComparisonOrderFromRow(mat, %u, done, todo)\n", row);
    while (todo != NULL) {
        // walk the todo linked list and see if we can fill in any regions with the current row.
        de_debug("starting to walk todo [%" PRIu64 ", %" PRIu64 "]\n", todo->start, todo->end);
        inGap = false;
        if (containsGaps == 0) {
            // if this row does not contain gaps we can shortcut this by accepting all todo regions.
            co = newMafTcComparisonOrder();
            co->ref = row;
            co->region = newMafTcRegion(todo->start, todo->end);
            de_debug("creating new region for row:%" PRIu64 ", position %" PRIu64 "\n", row, todo->start);
            // insert into the start of the list
            if (done != NULL) {
                co->next = *done;
            }
            *done = co;
            co = NULL;
            todo = todo->next;
            continue;
        }
        for (uint64_t i = todo->start; i <= todo->end; ++i) {
            // walk the current region.
            de_debug("position %" PRIu64 "\n", i);
            if (mat[row][i] == '-') {
                // inside a gap region
                de_debug("now inside a gap region\n");
                if (!inGap) {
                    de_debug("previously was not inGap\n");
                    // transition from sequence into gap
                    // no matter what we create a new todoRegion:
                    inGap = true;
                    if (newTodo == NULL) {
                        de_debug("newTodo == NULL, creating newTodo & tail at %" PRIu64 "\n", i);
                        newTodo = newMafTcRegion(i, i);
                        newTodoTail = newTodo;
                    } else {
                        de_debug("newTodo != NULL, creating new newTodoTail at %" PRIu64 "\n", i);
                        newTodoTail->next = newMafTcRegion(i, i);
                        newTodoTail = newTodoTail->next;
                    }
                } else {
                    de_debug("remaining inGap, extending gap end to %" PRIu64 "\n", i);
                    // remain in gap
                    de_debug("newTodoTail->start: %" PRIu64 ", ->end: %" PRIu64 "\n", 
                             newTodoTail->start, newTodoTail->end);
                    newTodoTail->end = i;
                    de_debug("newTodoTail->start: %" PRIu64 ", ->end: %" PRIu64 "\n", 
                             newTodoTail->start, newTodoTail->end);
                }
                de_debug("current newTodo:\n");
                printRegion(newTodo);
            } else {
                de_debug("now inside a sequence region\n");
                // inside a sequence region
                if (inGap) {
                    de_debug("previously was inGap\n");
                    // transition from gap to sequence
                    inGap = false;
                    co = newMafTcComparisonOrder();
                    co->ref = row;
                    co->region = newMafTcRegion(i, i);
                    de_debug("creating new region for row:%" PRIu64 ", position %" PRIu64 "\n", row, i);
                    // insert into the start of the list
                    if (done != NULL) {
                        // de_debug("done != NULL, inserting this comparison into the start of the list.\n");
                        co->next = *done;
                    }
                    *done = co;
                } else {
                    de_debug("remaining in sequence\n");
                    // remain in sequence
                    if (co != NULL) {
                        de_debug("co != NULL, moving ref:%" PRIu64 " with "
                                 "start %" PRIu64 ", end to %" PRIu64 "\n", co->ref, co->region->start, i);
                        co->region->end = i;
                    } else {
                        de_debug("co == NULL, ref: %" PRIu64 ", must create new region at %" PRIu64 "\n", 
                                 row, i);
                        co = newMafTcComparisonOrder();
                        co->ref = row;
                        co->region = newMafTcRegion(i, i);
                        // insert into the start of the `done' list
                        if (done != NULL) {
                            de_debug("done != NULL, inserting this comparison into head of list\n");
                            co->next = *done;
                        } else {
                            de_debug("done == NULL, inserting a new comparison into head of list\n");
                            done = (mafTcComparisonOrder_t **) de_malloc(sizeof(*done));
                        }
                        *done = co;
                    }
                }
            }
        }
        co = NULL;
        todo = todo->next;
    }
    destroyMafTcRegionList(headTodo);
    return newTodo;
}
int64_t localSeqCoordsToGlobalPositiveStartCoords(int64_t c, uint64_t start, uint64_t sourceLength, 
                                                  char strand, uint64_t length) {
    // given a coordinate inside of a sequence with respect to the start (which is 0), 
    // transform the START coordinate of a BLOCK to positive full source coordinates.
    if (strand == '+') {
        return localSeqCoordsToGlobalPositiveCoords(c, start, sourceLength, strand);
    } else {
        return localSeqCoordsToGlobalPositiveCoords(c, start + length - 1, sourceLength, strand);
    }
}
int64_t localSeqCoordsToGlobalPositiveCoords(int64_t c, uint64_t start, uint64_t sourceLength, char strand) {
    // given a coordinate inside of a sequence with respect to the start (which is 0), 
    // transform the coordinate to positive full source coordinates.
    if (strand == '+') {
        return (int64_t) (start + c);
    } else {
        return (int64_t) (sourceLength - 1 - (start + c));
    }
}
int64_t localSeqCoords(uint64_t p, char *s, mafCoordinatePair_t *b, int containsGaps) {
    // given a coordinate inside of a sequence with respect to the start (which is 0), 
    // walk backwards and anytime you see a gap drop one from the coordinate to return
    if (containsGaps == 0) {
        // store bookmark
        b->a = p;
        b->b = p;
        return b->b;
    }
    int64_t bases = 0;
    for (int64_t i = p; i >= 0; --i) {
        if (b->a == i) {
            if (s[i] != '-' && b->b == -1) {
                bases += 1;
            }
            b->a = p;
            b->b += bases;
            return b->b;
        }
        if (s[i] != '-') {
            ++bases;
        }
    }
    b->a = p;
    b->b = bases;
    return b->b;
}
uint64_t g_numPinches = 0;
void processPairForPinching(stPinchThreadSet *threadSet, stPinchThread *a, uint64_t aGlobalStart, 
                            uint64_t aGlobalLength, int aStrand, 
                            char *aSeq, stPinchThread *b, uint64_t bGlobalStart, uint64_t bGlobalLength,
                            int bStrand, char *bSeq, uint64_t regionStart, uint64_t regionEnd,
                            mafCoordinatePair_t aBookmark, mafCoordinatePair_t bBookmark, 
                            int aContainsGaps, int bContainsGaps, 
                            void (*pinchFunction)(stPinchThread *, stPinchThread *, int64_t, int64_t, int64_t, bool)) {
    // perform a pinch operation for regions of bSeq that are not gaps, i.e. `-'
    // GlobalStart is the positive strand position (zero based) coordinate of the start of this block
    // the nasty construction of passing in the pinch function is so that we can more easily unit test 
    // this code. Sorry about that.
    (void) (threadSet);
    uint64_t length = 0, localPos = 0;
    uint64_t localBlockStart = localPos;
    int64_t aLocalPosCoords, bLocalPosCoords, aGlobalPosCoords, bGlobalPosCoords;
    de_debug("aGlobalStart: %" PRIu64 ", bGlobalStart: %" PRIu64 ", pos & regionStart: %" PRIu64 ", regionEnd: %" PRIu64 "\n", 
             aGlobalStart, bGlobalStart, regionStart, regionEnd);
    de_debug("aSeq: %s\n", aSeq);
    de_debug("bSeq: %s\n", bSeq);
    bool inBlock = false;
    for (localPos = regionStart; localPos < (regionEnd + 1); ++localPos) {
        de_debug("pos: %" PRIu64 " \n", localPos);
        if (bSeq[localPos] == '-') {
            if (inBlock) {
                inBlock = false;
                aLocalPosCoords = localSeqCoords(localBlockStart, aSeq, &aBookmark, aContainsGaps);
                aGlobalPosCoords = localSeqCoordsToGlobalPositiveStartCoords(aLocalPosCoords, aGlobalStart, 
                                                                             aGlobalLength, aStrand, length);
                bLocalPosCoords = localSeqCoords(localBlockStart, bSeq, &bBookmark, bContainsGaps);
                bGlobalPosCoords = localSeqCoordsToGlobalPositiveStartCoords(bLocalPosCoords, bGlobalStart, 
                                                                             bGlobalLength, bStrand, length);
                de_debug("bSeq[pos] is a gap && inBlock\n");
                de_debug("maybe little pinch (#1)? p:%" PRIu64 ", l:%" PRIu64 "\n", localBlockStart, length);
                de_debug("a coords local: %" PRIi64 " global: %" PRIi64
                         ", b coords local: %" PRIi64 ", global: %" PRIi64 "\n",
                         aLocalPosCoords, aGlobalPosCoords, bLocalPosCoords, bGlobalPosCoords);
                pinchFunction(a, b, aGlobalPosCoords, bGlobalPosCoords, length, (aStrand == bStrand));
                ++g_numPinches;
                length = 0;
            }
        } else {
            if (!inBlock) {
                de_debug("bSeq[pos] is not a gap char && !inBlock, starting "
                         "new block at %" PRIu64 "\n", localPos);
                inBlock = true;
                localBlockStart = localPos;
            }
            ++length; // length of this segment
        }
    }
    de_debug("done walking comparison region\n");
    if (inBlock) {
        inBlock = false;
        aLocalPosCoords = localSeqCoords(localBlockStart, aSeq, &aBookmark, aContainsGaps);
        aGlobalPosCoords = localSeqCoordsToGlobalPositiveStartCoords(aLocalPosCoords, aGlobalStart, 
                                                                     aGlobalLength, aStrand, length);
        bLocalPosCoords = localSeqCoords(localBlockStart, bSeq, &bBookmark, bContainsGaps);
        bGlobalPosCoords = localSeqCoordsToGlobalPositiveStartCoords(bLocalPosCoords, bGlobalStart, 
                                                                     bGlobalLength, bStrand, length);
        de_debug("maybe little pinch (#2)? p:%" PRIu64 ", l:%" PRIu64 "\n", localBlockStart, length);
        de_debug("a coords local: %" PRIi64 " global: %" PRIi64
                         ", b coords local: %" PRIi64 ", global: %" PRIi64 "\n",
                         aLocalPosCoords, aGlobalPosCoords, bLocalPosCoords, bGlobalPosCoords);
        pinchFunction(a, b, aGlobalPosCoords, bGlobalPosCoords, length, (aStrand == bStrand));
        length = 0;
        ++g_numPinches;
    }
}
static void printMatrix(char **mat, uint64_t n) {
    fprintf(stderr, "printMatrix()\n");
    for (uint64_t i = 0; i < n; ++i)
        fprintf(stderr, "%" PRIu64 ": %s\n", i, mat[i]);
    fprintf(stderr, "\n");
}
static void printu32Array(uint64_t *a, uint64_t n) {
    for (uint64_t i = 0; i < n; ++i)
        fprintf(stderr, "%" PRIu64 "%s", a[i], (i == n - 1) ? "\n" : ", ");
}
mafCoordinatePair_t* newCoordinatePairArray(uint64_t numSeqs, char **seqs) {
    // create a new array of mafCoordinatePair_t* and then initialize it using the 
    // decomposed information from the mafBlock_t (strands, starts, sourceLengths, seqs, etc)
    mafCoordinatePair_t *a = (mafCoordinatePair_t*) de_malloc(sizeof(*a) * numSeqs);
    uint64_t pos;
    for (uint64_t i = 0; i < numSeqs; ++i) {
        pos = 0;
        while ((seqs[i][pos] == '-') && (pos < strlen(seqs[i]))) {
            // find the first non gap position in the sequence and start the bookmark there.
            ++pos;
        }
        a[i].a = (int64_t) pos;
        a[i].b = -1;
    }
    return a;
}
void destroyCoordinatePairArray(mafCoordinatePair_t *cp) {
    free(cp);
}
int** getVizMatrix(mafBlock_t *mb, unsigned n, unsigned m) {
    // currently this is not stored and must be built
    // should return a matrix containing the alignment, one row per sequence
    int** matrix = NULL;
    matrix = (int**) de_malloc(sizeof(int*) * n);
    unsigned i;
    unsigned len;
    char *s = NULL;
    for (i = 0; i < n; ++i)
        matrix[i] = (int*) de_malloc(sizeof(int) * m);
    mafLine_t *ml = maf_mafBlock_getHeadLine(mb);
    i = 0;
    while (ml != NULL) {
        while (maf_mafLine_getType(ml) != 's') {
            ml = maf_mafLine_getNext(ml);
        }
        s = maf_mafLine_getSequence(ml);
        len = maf_mafLine_getSequenceFieldLength(ml);
        for (unsigned j = 0; j < len; ++j) {
            if (s[j] == '-') {
                matrix[i][j] = 0;
            } else {
                matrix[i][j] = 1;
            }
        }
        ++i;
        ml = maf_mafLine_getNext(ml);
    }
    return matrix;
}
void updateVizMatrix(int **mat, mafTcComparisonOrder_t *co) {
    // given an int matrix that contains debug information about a block and 
    // a comparison order list that contains information about which sequence
    // segments to use as reference for pinches, update the matrix.
    while (co != NULL) {
        for (uint64_t i = co->region->start; i <= co->region->end; ++i) {
            mat[co->ref][i] = 2;
        }
        co = co->next;
    }
}
void printVizMatrix(int **mat, uint64_t n, uint64_t m) {
    for (uint64_t i = 0; i < n; ++i) {
        for (uint64_t j = 0; j < m; ++j) {
            if (mat[i][j] == 0) {
                printf("-");
            } else if (mat[i][j] == 1) {
                printf(".");
            } else if (mat[i][j] == 2) {
                printf("*");
            }
        }
        printf("\n");
    }
    printf("\n");
}
void printTodoArray(mafTcRegion_t *reg, unsigned max) {
    unsigned i = 0;
    while (reg != NULL) {
        while (i < reg->start) {
            printf("_");
            ++i;
        }
        while (i <= reg->end) {
            printf("@");
            ++i;
        }
        reg = reg->next;
    }
    while (i < max) {
        printf("_");
        ++i;
    }
    printf("\n");
}
void destroyVizMatrix(int **mat, unsigned n) {
    if (mat == NULL) {
        return;
    }
    for (unsigned i = 0; i < n; ++i) {
        free(mat[i]);
    }
    free(mat);
}
int cmp_by_gaps(const void *a, const void *b) {
    mafBlockSort_t **ia = (mafBlockSort_t **) a;
    mafBlockSort_t **ib = (mafBlockSort_t **) b;
    return ((*ia)->value >= (*ib)->value);
}
int64_t g_stableOrder = INT64_MIN;
void mafBlock_sortBlockByIncreasingGap(mafBlock_t *mb) {
    // take a pointer to mafblock, sort the maflines by 
    // the number of gaps they contain smallest to largest. All non-sequence
    // lines are stored in order of appearance and in the output will appear
    // before the sequence lines. THIS MEANS THAT THIS SORT WILL BREAK UP
    // RELATIONSHIPS BETWEEN COMMENT LINES AND 'e', 'q' and 'i' LINES AND
    // THEIR PARTNER 's' LINES. BEWARE.
    mafLine_t *next = NULL, *ml = maf_mafBlock_getHeadLine(mb);
    assert(ml != NULL);
    int64_t i, n = maf_mafBlock_getNumberOfLines(mb);
    int64_t value;
    mafBlockSort_t **array = (mafBlockSort_t **) de_malloc(sizeof(mafBlockSort_t *) * n);
    for (i = 0; i < n; ++i) {
        array[i] = (mafBlockSort_t *) de_malloc(sizeof(mafBlockSort_t));
        array[i]->ml = ml;
        if (maf_mafLine_getType(ml) == 's') {
            value = maf_mafLine_getSequenceFieldLength(ml) - maf_mafLine_getLength(ml);
            if (value == 0) {
                // force stability for blocks without gaps
                array[i]->value = ++g_stableOrder;
            } else {
                array[i]->value = value;
            }
        } else {
            array[i]->value = ++g_stableOrder;
        }
        ml = maf_mafLine_getNext(ml);
    }
    // sort
    qsort(array, n, sizeof(mafBlockSort_t *), cmp_by_gaps);
    // rebuild mafBlockList
    maf_mafBlock_setHeadLine(mb, array[0]->ml);
    for (i = 0; i < (n - 1); ++i) {
        ml = array[i]->ml;
        maf_mafLine_setNext(ml, array[i + 1]->ml);
    }
    if (ml != NULL) {
        next = maf_mafLine_getNext(ml);
        if (next != NULL) {
            maf_mafLine_setNext(next, NULL);
        }
    }
    if (n > 1) {
        i = n - 1;
    } else {
        i = 0;
    }
    maf_mafBlock_setTailLine(mb, array[i]->ml);
    // clean up
    for (i = 0; i < n; ++i)
        free(array[i]);
    free(array);
}
void walkBlockAddingAlignments(mafBlock_t *mb, stPinchThreadSet *threadSet) {
    // for a given block, add the alignment information to the threadset.
    de_debug("walkBlockAddingAlignments():\n");
    if (g_isSort)
        mafBlock_sortBlockByIncreasingGap(mb);
    uint64_t numSeqs = maf_mafBlock_getNumberOfSequences(mb);
    if (numSeqs < 1) 
        return;
    uint64_t seqFieldLength = maf_mafBlock_getSequenceFieldLength(mb);
    char **mat = maf_mafBlock_getSequenceMatrix(mb, numSeqs, seqFieldLength);
    int **vizMat = NULL;
    if (g_debug_flag) {
        vizMat = getVizMatrix(mb, numSeqs, seqFieldLength);
    }
    char *strands = maf_mafBlock_getStrandArray(mb);
    char **names = maf_mafBlock_getSpeciesArray(mb);
    uint64_t *starts = maf_mafBlock_getStartArray(mb);
    uint64_t *sourceLengths = maf_mafBlock_getSourceLengthArray(mb);
    uint64_t *lengths = maf_mafBlock_getSequenceLengthArray(mb);
    // coordinate bookmarks are used to store the mapping between local block position
    // and local sequence coordinate positions, ie local block position minus gap positions.
    mafCoordinatePair_t *bookmarks = newCoordinatePairArray(numSeqs, mat);
    // comparison order coordinates are relative to the block
    mafTcComparisonOrder_t *c = getComparisonOrderFromMatrix(mat, numSeqs, seqFieldLength, lengths, vizMat);
    de_debug("comparisonOrder_t obtained\n");
    mafTcComparisonOrder_t *tmp = NULL;
    stPinchThread *a = NULL, *b = NULL;
    while (c != NULL) {
        a = stPinchThreadSet_getThread(threadSet, stHash_stringKey(names[c->ref]));
        assert(a != NULL);
        for (uint64_t r = c->ref + 1; r < numSeqs; ++r) {
            b = stPinchThreadSet_getThread(threadSet, stHash_stringKey(names[r]));
            assert(b != NULL);
            processPairForPinching(threadSet, a, starts[c->ref], sourceLengths[c->ref], strands[c->ref],
                                   mat[c->ref], b, starts[r], sourceLengths[r], strands[r], mat[r], 
                                   c->region->start, c->region->end, bookmarks[c->ref], bookmarks[r],
                                   (lengths[c->ref] != seqFieldLength), (lengths[r] != seqFieldLength), 
                                   stPinchThread_pinch);
        }
        tmp = c;
        c = c->next;
        destroyMafTcRegion(tmp->region);
        free(tmp);
    }
    // cleanup
    maf_mafBlock_destroySequenceMatrix(mat, numSeqs);
    destroyVizMatrix(vizMat, numSeqs);
    free(strands);
    free(starts);
    free(sourceLengths);
    free(lengths);
    destroyCoordinatePairArray(bookmarks);
    for (uint64_t i = 0; i < numSeqs; ++i) {
        free(names[i]);
    }
    free(names);
}
void addAlignmentsToThreadSet(mafFileApi_t *mfa, stPinchThreadSet *threadSet) {
    mafBlock_t *mb = NULL;
    while ((mb = maf_readBlock(mfa)) != NULL) {
        walkBlockAddingAlignments(mb, threadSet);
        maf_destroyMafBlockList(mb);
    }
    stPinchThreadSet_joinTrivialBoundaries(threadSet);
}
uint64_t getMaxNameLength(stHash *hash) {
    // utility function to find out the length of the longest sequence name in the hash.
    stHashIterator *hit = stHash_getIterator(hash);
    char *key = NULL;
    uint64_t max = 0;
    while ((key = stHash_getNext(hit)) != NULL) {
        if (max < strlen(key) + 2)
            max = strlen(key) + 2;
    }
    stHash_destructIterator(hit);
    return max;
}
void getMaxFieldLengths(stHash *hash, stHash *nameHash, stPinchBlock *block, uint64_t *maxStart, 
                        uint64_t *maxLength, uint64_t *maxSource) {
    // utility function to find out the length of the longest field members.
    stPinchBlockIt thisSegIt = stPinchBlock_getSegmentIterator(block);
    stPinchSegment *thisSeg = NULL;
    *maxStart = 0;
    *maxLength = 0;
    *maxSource = 0;
    char *temp = NULL;
    char *key = NULL;
    char strand = '\0';
    int64_t *intKey = NULL;
    while ((thisSeg = stPinchBlockIt_getNext(&thisSegIt)) != NULL) {
        strand = stPinchSegment_getBlockOrientation(thisSeg) == 1 ? '+' : '-';
        intKey = (int64_t*) st_malloc(sizeof(*intKey));
        *intKey = stPinchSegment_getName(thisSeg);
        key = (char*) stHash_search(nameHash, (void*) intKey);
        temp = (char*) de_malloc(kMaxStringLength);
        if (strand == '+') {
            sprintf(temp, "%" PRIi64, stPinchSegment_getStart(thisSeg));
        } else {
            sprintf(temp, "%" PRIi64,
                    (((int64_t)((mafTcSeq_t*)stHash_search(hash, key))->length) - 
                     stPinchSegment_getStart(thisSeg) - stPinchSegment_getLength(thisSeg)));
        }
        if (*maxStart < strlen(temp)) {
            *maxStart = strlen(temp);
        }
        free(temp);
        temp = (char*) de_malloc(kMaxStringLength);
        sprintf(temp, "%" PRIi64, stPinchSegment_getLength(thisSeg));
        if (*maxLength < strlen(temp))
            *maxLength = strlen(temp);
        free(temp);
        temp = (char*) de_malloc(kMaxStringLength);
        sprintf(temp, "%" PRIu64, ((mafTcSeq_t*)stHash_search(hash, key))->length);
        if (*maxSource < strlen(temp))
            *maxSource = strlen(temp);
        free(temp);
        free(intKey);
    }
}
char* getSequenceSubset(char *seq, int64_t start, char strand, int64_t length) {
    // used to extract a copy of a small subset of a sequnce, for use when printing
    // out an alignment.
    char *out = (char*) de_malloc(sizeof(*out) * length + 1);
    int64_t i, j;
    for (i = start, j = 0; i < start + length; ++i, ++j) {
        out[j] = seq[i];
    }
    out[j] = '\0';
    if (strand == '-')
        reverseComplementSequence(out, strlen(out));
    return out;
}
void reportTransitiveClosure(stPinchThreadSet *threadSet, stHash *hash, stHash *nameHash) {
    // walk the completed threadSet and report back the blocks that form the transitive closure
    // of the alignment.
    stPinchThreadSetBlockIt thisBlockIt = stPinchThreadSet_getBlockIt(threadSet);
    stPinchBlock *thisBlock = NULL;
    stPinchBlockIt thisSegIt;
    stPinchSegment *thisSeg = NULL;
    char *key = NULL;
    char *seq = NULL;
    char strand = '\0';
    printf("##maf version=1\n");
    printf("# mafTransitiveClosure %s, build: %s, %s, %s\n\n", g_version, g_build_date, 
           g_build_git_branch, g_build_git_sha);
    uint64_t maxNameLength, maxStartLength, maxLengthLength, maxSourceLengthLength;
    int64_t xformedStart;
    int64_t *intKey = NULL;
    maxNameLength = getMaxNameLength(hash);
    while ((thisBlock = stPinchThreadSetBlockIt_getNext(&thisBlockIt)) != NULL) {
        getMaxFieldLengths(hash, nameHash, thisBlock, &maxStartLength,
                           &maxLengthLength, &maxSourceLengthLength);
        printf("a degree=%" PRIu64 "\n", stPinchBlock_getDegree(thisBlock));
        thisSegIt = stPinchBlock_getSegmentIterator(thisBlock);
        while ((thisSeg = stPinchBlockIt_getNext(&thisSegIt)) != NULL) {
            intKey = (int64_t *) st_malloc(sizeof(*intKey));
            *intKey = stPinchSegment_getName(thisSeg);
            strand = stPinchSegment_getBlockOrientation(thisSeg) == 1 ? '+' : '-';
            key = (char*)stHash_search(nameHash, (void *)intKey);
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
            printf("s %-*s %*" PRIi32 " %*" PRIi32 " %c %*" PRIu32 " %s\n",
                   (uint32_t)maxNameLength, key,
                   (uint32_t)maxStartLength, (int32_t)xformedStart,
                   (uint32_t)maxLengthLength, (uint32_t)stPinchSegment_getLength(thisSeg),
                   strand,
                   (uint32_t)maxSourceLengthLength, (uint32_t)((mafTcSeq_t*)stHash_search(hash, key))->length,
                   seq
                   );
            free(seq);
            free(intKey);
        }
        printf("\n");
    }
    printf("\n");
} 
int main(int argc, char **argv) {
    (void) (printMatrix);
    (void) (printu32Array);
    (void) (reportSequenceHash);
    (void) (printComparisonOrder);
    (void) (printRegion);
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
