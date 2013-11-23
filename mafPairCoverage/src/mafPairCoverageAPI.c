/*
 * Copyright (C) 2011-2013 by
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
#include <errno.h> // file existence via ENOENT, errno
#include <getopt.h>
#include <inttypes.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "common.h"
#include "sharedMaf.h"
#include "mafPairCoverageAPI.h"
#include "bioioC.h" // benLine()

struct mafCoverageCount {
  // used to slice the coverage data down to the sequence level
  uint64_t sourceLength; // length of the genome / chromosme / sequence in question
  uint64_t observedLength; // Number of residues actually observed in the alignment file
  uint64_t count;
  uint64_t inRegion;
  uint64_t outRegion;
};

struct _BinContainer {
  // contains arrays of counts, used by the --bin* options
  int64_t bin_start; // starting base. i.e. --bin_start
  int64_t bin_end; // starting base. i.e. --bin_start
  int64_t bin_length;
  int64_t num_bins;
  uint64_t *bins;
};

mafCoverageCount_t* createMafCoverageCount(void) {
  mafCoverageCount_t *mcct = (mafCoverageCount_t *)st_malloc(sizeof(*mcct));
  mcct->sourceLength = 0; // sequence source length
  mcct->observedLength = 0; // sequence source length
  mcct->count = 0; // number of aligned positions
  mcct->inRegion = 0; // number of aligned positions inside of bed specified region,
  mcct->outRegion = 0; // number of aligned positions outside of bed specified region,
  return mcct;
}
uint64_t mafCoverageCount_getSourceLength(mafCoverageCount_t *mcct) {
  return mcct->sourceLength;
}
uint64_t mafCoverageCount_getObservedLength(mafCoverageCount_t *mcct) {
  return mcct->observedLength;
}
uint64_t mafCoverageCount_getCount(mafCoverageCount_t *mcct) {
  return mcct->count;
}
uint64_t mafCoverageCount_getInRegion(mafCoverageCount_t *mcct) {
  return mcct->inRegion;
}
uint64_t mafCoverageCount_getOutRegion(mafCoverageCount_t *mcct) {
  return mcct->outRegion;
}
void mafCoverageCount_setSourceLength(mafCoverageCount_t *mcct, uint64_t n) {
  mcct->sourceLength = n;
}
void mafCoverageCount_setCount(mafCoverageCount_t *mcct, uint64_t n) {
  mcct->count = n;
}
void mafCoverageCount_setInRegion(mafCoverageCount_t *mcct, uint64_t n) {
  mcct->inRegion = n;
}
void mafCoverageCount_setOutRegion(mafCoverageCount_t *mcct, uint64_t n) {
  mcct->outRegion = n;
}
int64_t binContainer_getBinStart(BinContainer *bc) {
  return bc->bin_start;
}
int64_t binContainer_getBinEnd(BinContainer *bc) {
  return bc->bin_end;
}
int64_t binContainer_getBinLength(BinContainer *bc) {
  return bc->bin_length;
}
int64_t binContainer_getNumBins(BinContainer *bc) {
  return bc->num_bins;
}
uint64_t* binContainer_getBins(BinContainer *bc) {
  return bc->bins;
}
uint64_t binContainer_accessBin(BinContainer *bc, int64_t i) {
  assert (bc->num_bins > i);
  return bc->bins[i];
}
void binContainer_setBinStart(BinContainer *bc, int64_t i) {
  bc->bin_start = i;
}
void binContainer_setBinEnd(BinContainer *bc, int64_t i) {
  bc->bin_end = i;
}
void binContainer_setBinLength(BinContainer *bc, int64_t i) {
  bc->bin_length = i;
}
void binContainer_incrementPosition(BinContainer *bc, int64_t pos) {
  // given a genomic position, figure out the local position relative
  // to the bin_start value, then figure out the proper bin index,
  // then request that bin index to increment.
  if (bc == NULL) {
    // this container has not even been initialized
    return;
  }
  if (binContainer_getBins(bc) == NULL) {
    // this container has only been initialized, not created.
    return;
  }
  if ((pos < binContainer_getBinStart(bc)) ||
      (pos > binContainer_getBinEnd(bc))) {
    // this position is outside of the range of the binning.
    return;
  }
  int64_t local_pos = pos - binContainer_getBinStart(bc);
  int64_t i = (int)floor(local_pos / binContainer_getBinLength(bc));
  // fprintf(stderr, "incrementing pos %" PRIi64", local_pos %" PRIi64 ", "
  //         "bin %" PRIi64 "\n", pos, local_pos, i);
  binContainer_incrementBin(bc, i);
}
void binContainer_incrementBin(BinContainer *bc, int64_t i) {
  // increment a bin at an index
  assert (bc->num_bins > i);
  ++(bc->bins[i]);
}
void binContainer_setBinValue(BinContainer *bc, int64_t i, int64_t v) {
  // set the value of a bin at an index.
  assert (bc->num_bins > i);
  bc->bins[i] = v;
}
bool is_wild(const char *s) {
  // return true if char array ends in *, false otherwise
  if (s[strlen(s) - 1] == '*')
    return true;
  return false;
}
bool searchMatched(mafLine_t *ml, const char *seq) {
  // report false if search did not match, true if it did
  if (maf_mafLine_getType(ml) != 's')
    return false;
  return searchMatched_(maf_mafLine_getSpecies(ml), seq);
}
bool searchMatched_(const char *target, const char *seq) {
  if (is_wild(seq)) {
    // only compare up to the wildcard character
    if (!(strncmp(target, seq, strlen(seq) - 1) == 0)) {
      return false;
    }
  } else {
    if (!(strcmp(target, seq) == 0)) {
      return false;
    }
  }
  return true;
}
static void quickSetup(mafLine_t *ml1, mafLine_t *ml2, uint64_t *pos1,
                       uint64_t *pos2, int *strand1, int *strand2) {
  // pulling this cruft into its own function
  if (maf_mafLine_getStrand(ml1) == '+') {
    *strand1 = 1;
  } else {
    *strand1 = -1;
  }
  if (maf_mafLine_getStrand(ml2) == '+') {
    *strand2 = 1;
  } else {
    *strand2 = -1;
  }
  *pos1 = maf_mafLine_getPositiveCoord(ml1);
  *pos2 = maf_mafLine_getPositiveCoord(ml2);
}
void compareLines(mafLine_t *ml1, mafLine_t *ml2, stHash *seq1Hash,
                  stHash *seq2Hash, uint64_t *alignedPositions,
                  stHash *intervalsHash, BinContainer *bin_container) {
  // look through the sequences position by position and count the number of
  // places where the two sequences contained aligned residues, i.e. neither
  // position contains a gap character.
  // ml1 and seq1Hash are both from the --seq1 command line, treat them as the
  // reference when evaluating bins / shards.
  char *s1 = maf_mafLine_getSequence(ml1);
  char *s2 = maf_mafLine_getSequence(ml2);
  char *seqName1 = maf_mafLine_getSpecies(ml1);
  char *seqName2 = maf_mafLine_getSpecies(ml2);
  assert(s1 != NULL);
  assert(s2 != NULL);
  uint64_t n = strlen(s1);
  assert(n == strlen(s2));
  mafCoverageCount_t *mcct1 = stHash_search(seq1Hash,
                                            maf_mafLine_getSpecies(ml1));
  mafCoverageCount_t *mcct2 = stHash_search(seq2Hash,
                                            maf_mafLine_getSpecies(ml2));
  if (mcct2 == NULL) {
    printf("hash2 is devoid of %s\n", maf_mafLine_getSpecies(ml2));
  }
  assert(mcct1 != NULL);
  assert(mcct2 != NULL);
  uint64_t s1_start = maf_mafLine_getStart(ml1);
  if (stHash_size(intervalsHash) == 0) {
    // no intervals: yay, life is simple! :D
    for (uint64_t i = 0; i < n; ++i) {
      if ((s1[i] != '-') && (s2[i] != '-')) {
        ++(*alignedPositions);
        ++(mcct1->count);
        ++(mcct2->count);
        binContainer_incrementPosition(bin_container, s1_start + i);
      }
    }
  } else {
    // intervals: boo, this shit is complicated! >:(
    uint64_t pos1, pos2;
    int strand1, strand2;
    quickSetup(ml1, ml2, &pos1, &pos2, &strand1, &strand2);
    for (uint64_t i = 0; i < n; ++i) {
      pos1 += strand1;
      pos2 += strand2;
      if ((s1[i] != '-') && (s2[i] != '-')) {
        ++(*alignedPositions);
        ++(mcct1->count);
        ++(mcct2->count);
        binContainer_incrementPosition(bin_container, s1_start + i);
        if (inInterval(intervalsHash, seqName1, pos1)) {
          // seq 1 is in the interval
          ++(mcct1->inRegion);
        } else {
          // seq 1 is not in the interval
          ++(mcct1->outRegion);
        }
        if (inInterval(intervalsHash, seqName2, pos2)) {
          // seq 2 is in the interval
          ++(mcct2->inRegion);
        } else {
          // seq2 is not in the interval
          ++(mcct2->outRegion);
        }
      }
    }
  }
}
void wrapDestroyMafLine(void *p) {
  maf_destroyMafLineList((mafLine_t *) p);
}
void checkBlock(mafBlock_t *b, const char *seq1, const char *seq2,
                stHash *seq1Hash, stHash *seq2Hash, uint64_t *alignedPositions,
                stHash *intervalsHash, BinContainer *bin_container) {
  // read through each line of a mafBlock and if the sequence matches the
  // region we're looking for, report the block.
  mafLine_t *ml1 = maf_mafBlock_getHeadLine(b);
  // do a quick scan for either seq before doing the full n^2 comparison.
  // i'm doing this because i know there are some transitively closed mafs
  // that contain upwards of tens of millions of rows... :/
  //
  bool has1 = false, has2 = false;
  stList *seq1List = stList_construct3(0, wrapDestroyMafLine);
  stList *seq2List = stList_construct3(0, wrapDestroyMafLine);
  mafCoverageCount_t *mcct1 = NULL, *mcct2 = NULL;
  while (ml1 != NULL) {
    if (searchMatched(ml1, seq1)) {
      has1 = true;
      stList_append(seq1List, maf_copyMafLine(ml1));
      // create an item in the hash for this sequence
      mcct1 = NULL;
      if ((mcct1 = stHash_search(seq1Hash, maf_mafLine_getSpecies(ml1)))
          == NULL) {
        // new sequence, add to the hash
        mcct1 = createMafCoverageCount();
        mcct1->sourceLength = maf_mafLine_getSourceLength(ml1);
        mcct1->observedLength = maf_mafLine_getLength(ml1);
        stHash_insert(seq1Hash,
                      stString_copy(maf_mafLine_getSpecies(ml1)),
                      mcct1);
      } else {
        assert(mcct1->sourceLength == maf_mafLine_getSourceLength(ml1));
        mcct1->observedLength += maf_mafLine_getLength(ml1);
      }
    }
    if (searchMatched(ml1, seq2)) {
      has2 = true;
      stList_append(seq2List, maf_copyMafLine(ml1));
      // create an item in the hash for this sequence
      mcct2 = NULL;
      if ((mcct2 = stHash_search(seq2Hash, maf_mafLine_getSpecies(ml1)))
          == NULL) {
        // new sequence, add to the hash
        mcct2 = createMafCoverageCount();
        mcct2->sourceLength = maf_mafLine_getSourceLength(ml1);
        mcct2->observedLength = maf_mafLine_getLength(ml1);
        stHash_insert(seq2Hash,
                      stString_copy(maf_mafLine_getSpecies(ml1)),
                      mcct2);
      } else {
        assert(mcct2->sourceLength == maf_mafLine_getSourceLength(ml1));
        mcct2->observedLength += maf_mafLine_getLength(ml1);
      }
    }
    ml1 = maf_mafLine_getNext(ml1);
  }
  if (!(has1 && has2)) {
    // if this block does not contain both seq1 and seq2, do nothing
    stList_destruct(seq1List);
    stList_destruct(seq2List);
    return;
  }
  // perform the full n^2 scan on the instances of seq1 and seq2 matches
  ml1 = NULL;
  mafLine_t *ml2 = NULL;
  stListIterator *sl_it1 = stList_getIterator(seq1List);
  stListIterator *sl_it2 = NULL;
  while ((ml1 = stList_getNext(sl_it1)) != NULL) {
    sl_it2 = stList_getIterator(seq2List);
    while ((ml2 = stList_getNext(sl_it2)) != NULL) {
      compareLines(ml1, ml2, seq1Hash, seq2Hash,
                   alignedPositions, intervalsHash, bin_container);
    }
    stList_destructIterator(sl_it2);
  }
  stList_destructIterator(sl_it1);
  stList_destruct(seq1List);
  stList_destruct(seq2List);
}


void processBody(mafFileApi_t *mfa, char *seq1, char *seq2, stHash *seq1Hash,
                 stHash *seq2Hash,
                 uint64_t *alignedPositions, stHash *intervalsHash,
                 BinContainer *bin_container) {
  mafBlock_t *thisBlock = NULL;
  *alignedPositions = 0;
  while ((thisBlock = maf_readBlock(mfa)) != NULL) {
    checkBlock(thisBlock, seq1, seq2, seq1Hash, seq2Hash,
               alignedPositions, intervalsHash, bin_container);
    maf_destroyMafBlockList(thisBlock);
  }
}


void parseBedFile(const char *filepath, stHash *intervalsHash) {
  /*
   * WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING
   * This function is near copy-pasta from mafComparator.c
   * WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING
   * takes a filepath and the intervalsHash, opens and reads the file,
   * adding intervals taken from each bed line to the intervalsHash.
   * intervals hash is keyed with the name of the sequence and valued with
   * a st_sortedset
   */
  FILE *fileHandle = fopen(filepath, "r");
  st_logDebug("Parsing the bed file: %s\n", filepath);
  if (fileHandle == NULL) {
    if (errno == ENOENT)
      fprintf(stderr, "ERROR, file %s does not exist.\n", filepath);
    else
      fprintf(stderr, "ERROR, unable to open %s\n", filepath);
    exit(EXIT_FAILURE);
  }
  int64_t nBytes = 100;
  char *cA2 = st_malloc(nBytes + 1);
  int64_t bytesRead = benLine(&cA2, &nBytes, fileHandle);
  while (bytesRead != -1) {
    if (bytesRead > 0) {
      int64_t start, stop;
      char *cA3 = stString_copy(cA2);
      int64_t i = sscanf(cA2, "%s %" PRIi64 " %" PRIi64 "", cA3,
                         &start, &stop);
      if (i != 3) {
        fprintf(stderr, "Error while parsing bed file %s, expected 3 "
                "column bed format:\n"
                "sequence_name\tstart\tend\n", filepath);
        exit(EXIT_FAILURE);
      }
      stSortedSet *intervals = stHash_search(intervalsHash, cA3);
      if (intervals == NULL) {
        intervals = stSortedSet_construct3(
                                           (int(*)(const void *, const void *))
                                           stIntTuple_cmpFn,
                                           (void(*)(void *))
                                           stIntTuple_destruct);
        stHash_insert(intervalsHash, stString_copy(cA3), intervals);
      }
      stIntTuple *j = stIntTuple_construct2(start, stop);
      stIntTuple *k = stSortedSet_searchLessThanOrEqual(intervals, j);
      if (k != NULL) {
        if (stIntTuple_get(k, 1) > start) {
          st_errAbort(
                      "Found an overlapping interval in the bed file: %s "
                      "%" PRIi64 " "
                      "%" PRIi64 " overlaps %s %" PRIi64 " %" PRIi64 "",
                      cA3, start, stop, cA2,
                      stIntTuple_get(k, 0),
                      stIntTuple_get(k, 1));
        }
      }
      k = stSortedSet_searchGreaterThanOrEqual(intervals, j);
      if (k != NULL) {
        if (stIntTuple_get(k, 1) < stop) {
          st_errAbort(
                      "Found an overlapping interval in the bed file: %s "
                      "%" PRIi64 " "
                      "%" PRIi64 " overlaps %s %" PRIi64 " %" PRIi64 "",
                      cA3, start, stop, cA2,
                      stIntTuple_get(k, 0),
                      stIntTuple_get(k, 1));
        }
      }
      assert(stSortedSet_search(intervals, j) == NULL);
      st_logDebug("Adding in an interval: %s %" PRIi64 " %" PRIi64 "\n",
                  cA3, start, stop);
      stSortedSet_insert(intervals, j);
      free(cA3);
    }
    bytesRead = benLine(&cA2, &nBytes, fileHandle);
  }
  free(cA2);
  fclose(fileHandle);
  st_logDebug("Finished parsing the bed file: %s\n", filepath);
}


bool inInterval(stHash *intervalsHash, char *seq, uint64_t pos) {
  /*
   * WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING
   * This function is copy pasta from mafComparatator/src/comparatorAPI.c
   * WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING
   * check to see if the sequence and position are within the intervals hash.
   *
   */
  stSortedSet *intervals = stHash_search(intervalsHash, seq);
  if (intervals == NULL) {
    // printf("%s %" PRIu64 " is not in the region, bad seq.\n", seq, pos);
    return false;
  }
  stIntTuple *i = stIntTuple_construct2( pos, INT64_MAX);
  stIntTuple *j = stSortedSet_searchLessThanOrEqual(intervals, i);
  stIntTuple_destruct(i);
  if (j == NULL) {
    // printf("%s %" PRIu64 " is not in the region, bad pos.\n", seq, pos);
    return false;
  }
  assert(stIntTuple_length(j) == 2);
  /* if (stIntTuple_get(j, 0) <= pos && pos < stIntTuple_get(j, 1)) { */
  /*      printf("%s %" PRIu64 " IS in the region.\n", seq, pos); */
  /* } else { */
  /*     printf("%s %" PRIu64 " IS NOT in the region.\n", seq, pos); */
  /* } */
  return stIntTuple_get(j, 0) <= pos && pos < stIntTuple_get(j, 1);
}

BinContainer* binContainer_init(void) {
  BinContainer *bc = st_malloc(sizeof(*bc));
  bc->bin_start = 0;
  bc->bin_end = 0;
  bc->bin_length = 0;
  bc->num_bins = 0;
  bc->bins = NULL;
  return bc;
}

BinContainer* binContainer_construct(int64_t bin_start, int64_t bin_end,
                                     int64_t bin_length) {
  BinContainer *bc = binContainer_init();
  bc->bin_start = bin_start;
  bc->bin_end = bin_end;
  bc->bin_length = bin_length;
  bc->num_bins = (int64_t)ceil((double) ((bin_end + 1) - bin_start) /
                               bin_length);
  bc->bins = st_calloc(bc->num_bins, sizeof(uint64_t));
  return bc;
}

void binContainer_destruct(BinContainer *bc) {
  if (bc == NULL) {
    return;
  }
  if (bc->bins != NULL) {
    free(bc->bins);
    bc->bins = NULL;
  }
  free(bc);
  bc = NULL;
}

void reportResultsBins(char *seq1, char *seq2, BinContainer *bc) {
  if (bc == NULL) {
    return;
  }
  if (binContainer_getBins(bc) == NULL) {
    return;
  }
  printf("# Bins\n# Ref_Bin_Starting_Pos\tCoverage\n");
  int64_t cur_pos = bc->bin_start;
  assert(bc->bin_length > 0);
  for (int64_t i = 0; i < bc->num_bins; ++i) {
    printf("%" PRIi64 "\t%15e\n", cur_pos,
           bc->bins[i] / (double)bc->bin_length);
  }
}
