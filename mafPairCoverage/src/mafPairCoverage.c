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
#include <getopt.h>
#include <math.h>  // ceil()
#include <stdbool.h>
#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>
#include <stdlib.h>
#include <string.h>
#include "common.h"
#include "sharedMaf.h"
#include "mafPairCoverage.h"
#include "mafPairCoverageAPI.h"
#include "buildVersion.h"

const char *g_version = "version 0.1 May 2013";
uint64_t getRegionSize(char *seq1, stHash *intervalsHash);

void version(void) {
  fprintf(stderr, "mafPairCoverage, %s\nbuild: %s, %s, %s\n\n", g_version,
          g_build_date, g_build_git_branch, g_build_git_sha);
}
void usage(void) {
  version();
  fprintf(stderr, "Usage: mafPairExtractor --maf [maf file] --seq1 [sequence "
          "name (and possibly chr),\n"
          "may end in * wildcard] --seq2 [sequence name (and possibly chr),\n"
          "may end in * wildcard [options]\n\n"
          "mafPairCoverage is a program that will look through a maf file "
          "block \n"
          "by block and check for a particular pair of sequences (allowing "
          "input \n"
          "sequence to end in wildcard *) and count the number of columns "
          "where \n"
          "the two sequences have residues aligned. Coverage of genome A "
          "onto \n"
          "genome B is then symmetrically calculated as the number of "
          "shared \n"
          "columns divided by the total size of genome B.\n\n");
  fprintf(stderr, "Options: \n");
  usageMessage('h', "help", "show this help message and exit.");
  usageMessage('m', "maf", "path to maf file.");
  usageMessage('\0', "seq1", "sequence name, e.g. `hg19*'. Accepts * "
               "wildcard at the end.");
  usageMessage('\0', "seq2", "sequence name, e.g. `mm9.chr9'. Accepts * "
               "wildcard at the end.");
  usageMessage('\0', "bed", "path to 3 column bedfile that will define "
               "regions of interest in output.");
  usageMessage('\0', "bin_start", "starting position (inclusive) of the sub-"
               "region to analyze.");
  usageMessage('\0', "bin_end", "ending position (inclusive) of the sub-"
               "region to analyze.");
  usageMessage('\0', "bin_length", "the length of each bin within the "
               "region. default=1000");
  usageMessage('v', "verbose", "turns on verbose output.");
  exit(EXIT_FAILURE);
}


void parseOptions(int argc, char **argv, char *filename, char *seq1Name,
                  char *seq2Name, stHash *intervalsHash, int64_t *bin_start,
                  int64_t *bin_end, int64_t *bin_length) {
  extern int g_debug_flag;
  extern int g_verbose_flag;
  int c;
  bool setMName = false, setS1Name = false, setS2Name = false;
  while (1) {
    static struct option longOptions[] = {
      {"debug", no_argument, 0, 'd'},
      {"verbose", no_argument, 0, 'v'},
      {"help", no_argument, 0, 'h'},
      {"version", no_argument, 0, 0},
      {"maf", required_argument, 0, 'm'},
      {"seq1", required_argument, 0, 0},
      {"seq2", required_argument, 0, 0},
      {"bed", required_argument, 0, 0},
      {"bin_start", required_argument, 0, 0},
      {"bin_end", required_argument, 0, 0},
      {"bin_length", required_argument, 0, 0},
      {0, 0, 0, 0}
    };
    size_t i;
    int longIndex = 0;
    c = getopt_long(argc, argv, "m:s:h:v:d",
                    longOptions, &longIndex);
    if (c == -1)
      break;
    switch (c) {
    case 0:
      if (strcmp("seq1", longOptions[longIndex].name) == 0) {
        setS1Name = true;
        strncpy(seq1Name, optarg, kMaxSeqName);
      } else if (strcmp("seq2", longOptions[longIndex].name) == 0) {
        setS2Name = true;
        strncpy(seq2Name, optarg, kMaxSeqName);
      } else if (strcmp("version", longOptions[longIndex].name) == 0) {
        version();
        exit(EXIT_SUCCESS);
      } else if (strcmp("bed", longOptions[longIndex].name) == 0) {
        parseBedFile(optarg, intervalsHash);
      } else if (strcmp("bin_start", longOptions[longIndex].name) == 0) {
        i = sscanf(optarg, "%" PRIi64, bin_start);
        assert(i == 1);
      } else if (strcmp("bin_end", longOptions[longIndex].name) == 0) {
        i = sscanf(optarg, "%" PRIi64, bin_end);
        assert(i == 1);
      } else if (strcmp("bin_length", longOptions[longIndex].name) == 0) {
        i = sscanf(optarg, "%" PRIi64, bin_length);
        assert(i == 1);
      }
      break;
    case 'm':
      setMName = true;
      strncpy(filename, optarg, kMaxSeqName);
      break;
    case 'v':
      g_verbose_flag++;
      break;
    case 'd':
      g_debug_flag = 1;
      break;
    case 'h':
    case '?':
      usage();
    break;
    default:
      abort();
    }
  }
  if (!(setMName && setS1Name && setS2Name)) {
    fprintf(stderr, "Error, specify --maf --seq1 --seq2\n");
    usage();
  }
  // Check there's nothing left over on the command line
  if (optind < argc) {
    char *errorString = de_malloc(kMaxSeqName);
    strcpy(errorString, "Unexpected arguments:");
    while (optind < argc) {
      strcat(errorString, " ");
      strcat(errorString, argv[optind++]);
    }
    fprintf(stderr, "%s\n", errorString);
    usage();
  }
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
  bc->num_bins = (int64_t)ceil((double) (bin_end - bin_start) / bin_length);
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

uint64_t getRegionSize(char *seq1, stHash *intervalsHash) {
  // go through intervalsHash and see if seq1 matches any of the keys we pull
  // out. if so, add up the size
  uint64_t n = 0;
  stHashIterator *hit = stHash_getIterator(intervalsHash);
  char *key = NULL;
  while ((key = stHash_getNext(hit)) != NULL) {
    stSortedSetIterator *sit = NULL;
    if (!searchMatched_(key, seq1)) {
      continue;
    }
    sit = stSortedSet_getIterator(stHash_search(intervalsHash, key));
    stIntTuple *inttup = NULL;
    while ((inttup = stSortedSet_getNext(sit)) != NULL) {
      n += stIntTuple_get(inttup, 1)  - stIntTuple_get(inttup, 0);
    }
    stSortedSet_destructIterator(sit);
  }
  stHash_destructIterator(hit);
  return n;
}


void reportResultsRegion(char *seq1, char *seq2, stHash *seq1Hash,
                         stHash *seq2Hash, uint64_t *alignedPositions,
                         stHash *intervalsHash) {
  /*
   * If there are results within the intervals hash, report those
   */
  if (stHash_size(intervalsHash) == 0) {
    return;
  }
  uint64_t tot1 = 0, tot1in = 0, tot1out = 0;
  stHashIterator *hit = stHash_getIterator(seq1Hash);
  double cov1in = 0.;
  char *key;
  // get the total in-region coverage number via the hash
  if (stHash_size(seq1Hash) > 0) {
    hit = stHash_getIterator(seq1Hash);
    while ((key = stHash_getNext(hit)) != NULL) {
      tot1 += getRegionSize(key, intervalsHash);
      tot1in += mafCoverageCount_getInRegion(stHash_search(seq1Hash, key));
      tot1out += mafCoverageCount_getOutRegion(stHash_search(seq1Hash, key));
    }
    stHash_destructIterator(hit);
    if (tot1 == 0) {
      cov1in = 0.;
    } else {
      cov1in = (double) tot1in / (double) tot1;
    }
  }

  // print the totals
  printf("\n# BED Regions\n");
  printf("#%19s\t%15s\t%15s\t%15s\t%15s\n", "Sequence",
         "Region Length", "In Regions", "Out of Regions", "Coverage");
  printf("%20s\t%15" PRIu64 "\t%15" PRIu64 "\t%15" PRIu64 "\t%15e\n",
         seq1, tot1, tot1in, tot1out, cov1in);

  // print each seq1 on its own
  if (stHash_size(seq1Hash) > 0) {
    hit = stHash_getIterator(seq1Hash);
    while ((key = stHash_getNext(hit)) != NULL) {
      if (stHash_search(intervalsHash, key) == NULL) {
        // don't report sequences that do not show up in the bed region file
        continue;
      }
      if ((double) mafCoverageCount_getInRegion(stHash_search(seq1Hash, key)) == 0) {
        cov1in = 0;
      } else {
        cov1in = (double) mafCoverageCount_getInRegion(stHash_search(seq1Hash, key)) /
          (double) getRegionSize(key, intervalsHash);
      }
      printf("%20s\t%15" PRIu64 "\t%15" PRIu64 "\t%15" PRIu64 "\t%15e\n",
             key, getRegionSize(key, intervalsHash),
             mafCoverageCount_getInRegion(stHash_search(seq1Hash, key)),
             mafCoverageCount_getOutRegion(stHash_search(seq1Hash, key)),
             cov1in);
    }
    stHash_destructIterator(hit);
  }
  // print each seq2 on its own
  if (stHash_size(seq2Hash) > 0) {
    hit = stHash_getIterator(seq2Hash);
    while ((key = stHash_getNext(hit)) != NULL) {
      if (stHash_search(intervalsHash, key) == NULL) {
        // don't report sequences that do not show up in the bed region file
        continue;
      }
      if ((double) mafCoverageCount_getInRegion(stHash_search(seq2Hash, key)) == 0) {
        cov1in = 0;
      } else {
        cov1in = (double) mafCoverageCount_getInRegion(stHash_search(seq2Hash, key)) /
          (double) getRegionSize(key, intervalsHash);

      }
      printf("%20s\t%15" PRIu64 "\t%15" PRIu64 "\t%15" PRIu64 "\t%15e\n",
             key, getRegionSize(key, intervalsHash),
             mafCoverageCount_getInRegion(stHash_search(seq2Hash, key)),
             mafCoverageCount_getOutRegion(stHash_search(seq2Hash, key)),
             cov1in);
    }
    stHash_destructIterator(hit);
  }
}


void reportResults(char *seq1, char *seq2, stHash *seq1Hash, stHash *seq2Hash,
                   uint64_t *alignedPositions) {
  uint64_t tot1 = 0, tot2 = 0, obstot1 = 0, obstot2 = 0;
  stHashIterator *hit = stHash_getIterator(seq1Hash);
  mafCoverageCount_t *mcct = NULL;
  double cov1 = 0., cov2 = 0., obscov1 = 0., obscov2 = 0.;
  char *key;
  while ((key = stHash_getNext(hit)) != NULL) {
    mcct = stHash_search(seq1Hash, key);
    assert(mcct != NULL);
    tot1 += mafCoverageCount_getSourceLength(mcct);
    obstot1 += mafCoverageCount_getObservedLength(mcct);
  }
  if (tot1 == 0) {
    cov1 = 0.0;
  } else {
    cov1 = *alignedPositions / (double) tot1;
  }
  if (obstot1 == 0) {
    obscov1 = 0.0;
  } else {
    obscov1 = *alignedPositions / (double) obstot1;
  }
  stHash_destructIterator(hit);
  hit = stHash_getIterator(seq2Hash);
  key = NULL;
  while ((key = stHash_getNext(hit)) != NULL) {
    mcct = (mafCoverageCount_t *) stHash_search(seq2Hash, key);
    assert(mcct != NULL);
    tot2 += mafCoverageCount_getSourceLength(mcct);
    obstot2 += mafCoverageCount_getObservedLength(mcct);
  }
  if (tot2 == 0) {
    cov2 = 0.0;
  } else {
    cov2 = *alignedPositions / (double) tot2;
  }
  if (obstot2 == 0) {
    obscov2 = 0.0;
  } else {
    obscov2 = *alignedPositions / (double) obstot2;
  }
  stHash_destructIterator(hit);
  // print the totals
  printf("# Overall\n");
  printf("#%19s\t%15s\t%15s\t%15s\t%15s\t%15s\n", "Sequence", "Source Length",
         "Obs Length",
         "Aligned Pos.", "Coverage", "Obs Coverage");
  printf("%20s\t%15" PRIu64 "\t%15" PRIu64 "\t%15" PRIu64 "\t%15e\t%15e\n",
         seq1, tot1, obstot1, *alignedPositions, cov1, obscov1);
  printf("%20s\t%15" PRIu64 "\t%15" PRIu64 "\t%15" PRIu64 "\t%15e\t%15e\n\n",
         seq2, tot2, obstot2, *alignedPositions, cov2, obscov2);

  printf("# Individual Sequences\n");
  printf("#%19s\t%15s\t%15s\t%15s\t%15s\t%15s\n", "Sequence", "Source Length",
         "Obs Length",
         "Aligned Pos.", "Coverage", "Obs Coverage");
  // print each seq1 on its own
  if (stHash_size(seq1Hash) > 0) {
    hit = stHash_getIterator(seq1Hash);
    while ((key = stHash_getNext(hit)) != NULL) {
      if  (mafCoverageCount_getSourceLength(stHash_search(seq1Hash, key)) == 0) {
        cov1 = 0.;
      } else {
        cov1 = (double) mafCoverageCount_getCount(stHash_search(seq1Hash, key)) /
          (double) mafCoverageCount_getSourceLength(stHash_search(seq1Hash, key));
      }
      if (mafCoverageCount_getObservedLength(stHash_search(seq1Hash, key)) == 0) {
        printf("%s has an observed length of 0!\n", key);
        obscov1 = 0.;
      } else {
        obscov1 = (double) mafCoverageCount_getCount(stHash_search(seq1Hash, key)) /
          (double) mafCoverageCount_getObservedLength(stHash_search(seq1Hash, key));
      }
      printf("%20s\t%15" PRIu64 "\t%15" PRIu64 "\t%15" PRIu64 "\t%15e\t%15e\n",
             key, mafCoverageCount_getSourceLength(stHash_search(seq1Hash, key)),
             mafCoverageCount_getObservedLength(stHash_search(seq1Hash, key)),
             mafCoverageCount_getCount(stHash_search(seq1Hash, key)), cov1, obscov1);
    }
    stHash_destructIterator(hit);
  }
  // print each seq2 on its own
  if (stHash_size(seq2Hash) > 0) {
    hit = stHash_getIterator(seq2Hash);
    while ((key = stHash_getNext(hit)) != NULL) {
      if  (mafCoverageCount_getSourceLength(stHash_search(seq2Hash, key)) == 0){
        cov2 = 0.;
      } else {
        cov2 = (double) mafCoverageCount_getCount(stHash_search(seq2Hash, key)) /
          (double) mafCoverageCount_getSourceLength(stHash_search(seq2Hash, key));
      }
      if (mafCoverageCount_getObservedLength(stHash_search(seq2Hash, key)) == 0){
        obscov2 = 0.;
      } else {
        obscov2 = (double) mafCoverageCount_getCount(stHash_search(seq2Hash, key)) /
          (double) mafCoverageCount_getObservedLength(stHash_search(seq2Hash, key));
      }
      printf("%20s\t%15" PRIu64 "\t%15" PRIu64 "\t%15" PRIu64 "\t%15e\t%15e\n",
             key, mafCoverageCount_getSourceLength(stHash_search(seq2Hash, key)),
             mafCoverageCount_getObservedLength(stHash_search(seq2Hash, key)),
             mafCoverageCount_getCount(stHash_search(seq2Hash, key)), cov2, obscov2);
    }
    stHash_destructIterator(hit);
  }
}


void reportResultsBins(char *seq1, char *seq2, BinContainer *bc) {
  if (bc == NULL) {
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


int main(int argc, char **argv) {
  extern const int kMaxStringLength;
  char seq1[kMaxSeqName];
  char seq2[kMaxSeqName];
  char filename[kMaxStringLength];
  int64_t bin_start = -1;  // sentinel value. real values > 0
  int64_t bin_end = -1;  // sentinel value. real values > 0
  int64_t bin_length = 1000;
  BinContainer *bin_container = NULL;
  stHash *intervalsHash = stHash_construct3(stHash_stringKey,
                                            stHash_stringEqualKey, free,
                                            (void(*)(void *))
                                            stSortedSet_destruct);
  parseOptions(argc, argv, filename, seq1, seq2, intervalsHash,
               &bin_start, &bin_end, &bin_length);
  if ((bin_start > 0) && (bin_end > 0) && (bin_length > 0)) {
    bin_container = binContainer_construct(bin_start, bin_end,
                                                         bin_length);
    if (is_wild(seq1)) {
      fprintf(stderr, "--seq1 may not contain a wild card character, *, if "
              "the binning options are used.\n");
      usage();
    }
  } else {
    bin_container = NULL;
  }
  mafFileApi_t *mfa = maf_newMfa(filename, "r");
  stHash *seq1Hash = stHash_construct3(stHash_stringKey, stHash_stringEqualKey,
                                       free, free);
  stHash *seq2Hash = stHash_construct3(stHash_stringKey, stHash_stringEqualKey,
                                       free, free);
  uint64_t alignedPositions = 0;
  processBody(mfa, seq1, seq2, seq1Hash, seq2Hash, &alignedPositions,
              intervalsHash, bin_container);
  reportResults(seq1, seq2, seq1Hash, seq2Hash, &alignedPositions);
  reportResultsRegion(seq1, seq2, seq1Hash, seq2Hash, &alignedPositions,
                      intervalsHash);
  reportResultsBins(seq1, seq2, bin_container);
  maf_destroyMfa(mfa);
  stHash_destruct(seq1Hash);
  stHash_destruct(seq2Hash);
  stHash_destruct(intervalsHash);
  binContainer_destruct(bin_container);
  return EXIT_SUCCESS;
}
