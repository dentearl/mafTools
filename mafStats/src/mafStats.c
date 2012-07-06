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

#include <getopt.h>
#include <math.h>
#include <stdbool.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include "sonLib.h"
#include "common.h"
#include "sharedMaf.h"
#include "mafStats.h"

const char *kVersion = "v0.1 July 2012";

typedef struct stats {
    char *filename;
    uint32_t numLines;
    uint32_t numHeaderLines;
    uint32_t numSeqLines;
    uint32_t numBlocks;
    uint32_t numELines;
    uint32_t numILines;
    uint32_t numQLines;
    uint32_t numCommentLines;
    uint32_t numGapCharacters;
    uint32_t numSeqCharacters;
    uint32_t sumSeqField;
    uint32_t maxSeqField;
    uint32_t sumNumSpeciesInBlock;
    uint32_t maxNumSpeciesInBlock;
    uint32_t sumBlockArea;
    uint32_t maxBlockArea;
    stHash *seqHash; // keyed with names, valued with uint32_t count of bases present
} stats_t;
typedef struct seq {
    char *name;
    uint32_t count;
} seq_t;
void usage(void) {
    fprintf(stderr, "mafStats, %s.\n", kVersion);
    fprintf(stderr, "Usage:  --maf [maf file] [options]\n\n"
            "description goes here.\n\n");
    fprintf(stderr, "Options: \n");
    usageMessage('h', "help", "show this help message and exit.");
    usageMessage('m', "maf", "path to the maf file.");
    usageMessage('v', "verbose", "turns on verbose output.");
    exit(EXIT_FAILURE);
}
void parseArgs(int argc, char **argv, char **filename) {
    extern int g_verbose_flag;
    extern int g_debug_flag;
    int c;
    bool setMName = false;
    while (1) {
        static struct option long_options[] = {
            {"debug", no_argument, 0, 'd'},
            {"verbose", no_argument, 0, 'v'},
            {"help", no_argument, 0, 'h'},
            {"maf",  required_argument, 0, 'm'},
            {0, 0, 0, 0}
        };
        int option_index = 0;
        c = getopt_long(argc, argv, "d:v:h:m:s",
                        long_options, &option_index);
        if (c == -1)
            break;
        switch (c) {
        case 'm':
            setMName = true;
            *filename = stString_copy(optarg);
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
    if (!(setMName)) {
        fprintf(stderr, "Error, specify --maf\n");
        usage();
    }
    // Check there's nothing left over on the command line 
    if (optind < argc) {
        char *errorString = st_malloc(kMaxSeqName);
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
stats_t* stats_create(char *filename) {
    stats_t *stats = (stats_t*) st_malloc(sizeof(*stats));
    stats->filename = filename; // NOT to be free'd in _destroy
    stats->numLines = 0;
    stats->numHeaderLines = 1; // the ##maf line
    stats->numBlocks = 0;
    stats->numSeqLines = 0;
    stats->numELines = 0;
    stats->numILines = 0;
    stats->numQLines = 0;
    stats->numCommentLines = 0;
    stats->numGapCharacters = 0;
    stats->numSeqCharacters = 0;
    stats->sumSeqField = 0;
    stats->maxSeqField = 0;
    stats->sumNumSpeciesInBlock = 0;
    stats->maxNumSpeciesInBlock = 0;
    stats->sumBlockArea = 0;
    stats->maxBlockArea = 0;
    stats->seqHash = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, free, free);
    return stats;
}
void stats_destroy(stats_t *stats) {
    stHash_destruct(stats->seqHash);
    free(stats);
    stats = NULL;
}
void countCharacters(char *seq, stats_t *stats) {
    size_t len = strlen(seq);
    for (size_t i = 0; i < len; ++i) {
        if (seq[i] == '-') {
            ++(stats->numGapCharacters);
        } else {
            ++(stats->numSeqCharacters);
        }
    }
}
void processBlock(mafBlock_t *mb, stats_t *stats) {
    mafLine_t *ml = maf_mafBlock_getHeadLine(mb);
    char t = '\0';
    char *name = NULL;
    uint32_t *v = NULL;
    uint32_t blockSeqFieldLength = 0;
    t = maf_mafLine_getType(ml);
    if (t == '#') {
        ++(stats->numCommentLines);
    } else if (t == 'a') {
        ++(stats->numBlocks);
    }
    while ((ml = maf_mafLine_getNext(ml)) != NULL) {
        t = maf_mafLine_getType(ml);
        if (t == 's') {
            ++(stats->numSeqLines);
            if (blockSeqFieldLength == 0) {
                blockSeqFieldLength = strlen(maf_mafLine_getSequence(ml));
                if (stats->maxSeqField < blockSeqFieldLength) {
                    stats->maxSeqField = blockSeqFieldLength;
                }
            }
            name = maf_mafLine_getSpecies(ml);
            stats->sumSeqField += maf_mafLine_getLength(ml);
            if (stHash_search(stats->seqHash, name) == NULL) {
                v = (uint32_t *) st_malloc(sizeof(*v));
                *v = maf_mafLine_getLength(ml);
                stHash_insert(stats->seqHash, stString_copy(name), v);
            } else {
                v = stHash_search(stats->seqHash, name);
                *v += maf_mafLine_getLength(ml);
            }
            countCharacters(maf_mafLine_getSequence(ml), stats);
        } else if (t == '#') {
            ++(stats->numCommentLines);
        } else if (t == 'e') {
            ++(stats->numELines);
        } else if (t == 'i') {
            ++(stats->numILines);
        } else if (t == 'q') {
            ++(stats->numQLines);
        } else if (t == 'h') {
            ++(stats->numHeaderLines);
        }
    }
    if (stats->maxBlockArea < maf_mafBlock_getNumberOfSequences(mb) * blockSeqFieldLength) {
        stats->maxBlockArea = maf_mafBlock_getNumberOfSequences(mb) * blockSeqFieldLength;
    }
    stats->sumBlockArea += maf_mafBlock_getNumberOfSequences(mb) * blockSeqFieldLength;
    if (stats->maxNumSpeciesInBlock < maf_mafBlock_getNumberOfSequences(mb)) {
        stats->maxNumSpeciesInBlock = maf_mafBlock_getNumberOfSequences(mb);
    }
    stats->sumNumSpeciesInBlock += maf_mafBlock_getNumberOfSequences(mb);
}
void recordStats(mafFileApi_t *mfa, stats_t *stats) {
    mafBlock_t *mb = NULL;
    while ((mb = maf_readBlock(mfa)) != NULL) {
        processBlock(mb, stats);
        maf_destroyMafBlockList(mb);
    }
    stats->numLines = maf_mafFileApi_getLineNumber(mfa);
}
void readFilesize(struct stat *fileStat, char **filesizeString) {
    char *s = st_malloc(kMaxStringLength);
    *filesizeString = s;
    double fs = fileStat->st_size;
    if (fs > 1024.0) {
        // KB
        fs /= 1024.0;
    } else {
        sprintf(s, "%d Bytes", (int)fileStat->st_size);
        return;
    }
    if (fs > 1024.0) {
        // MB
        fs /= 1024.0;
    } else {
        sprintf(s, "%.2f KB", fs);
        return;
    }
    if (fs > 1024.0) {
        // GB
        fs /= 1024.0;
    } else {
        sprintf(s, "%.2f MB", fs);
        return;
    }
    if (fs > 1024.0) {
        // TB
        fs /= 1024.0;
        sprintf(s, "%.2f TB", fs);
        return;
    } else {
        sprintf(s, "%.2f GB", fs);
        return;
    }
}
int cmp_seq(const void *a, const void *b) {
    seq_t **ia = (seq_t **) a;
    seq_t **ib = (seq_t **) b;
    return ((*ia)->count < (*ib)->count);
}
void reportHash(stHash *hash) {
    int32_t n = stHash_size(hash);
    seq_t **order = (seq_t **) st_malloc(sizeof(*order) * n);
    stHashIterator *hit = stHash_getIterator(hash);
    char *key = NULL;
    int32_t i = 0;
    uint32_t total = 0;
    while ((key = stHash_getNext(hit)) != NULL) {
        order[i] = (seq_t*) st_malloc(sizeof(seq_t));
        order[i]->name = key;
        order[i]->count = *(uint32_t*)(stHash_search(hash, key));
        total += order[i++]->count;
    }
    qsort(order, n, sizeof(seq_t*), cmp_seq);
    for (i = 0; i < n; ++i) {
        printf("%25s: %12" PRIu32 " (%6.2f%%)\n", order[i]->name, order[i]->count, 
               100.0 * order[i]->count / total);
    }
    printf("%25s: %12" PRIu32 " (100.00%%)\n", "total", total);
}
void reportStats(stats_t *stats) {
    struct stat fileStat;
    stat(stats->filename, &fileStat);
    printf("%s\n", stats->filename);
    printf("------------------------------\n");
    char *filesizeString;
    readFilesize(&fileStat, &filesizeString);
    printf("File size:     %10s\n", filesizeString);
    printf("Lines:         %10" PRIu32 "\n", stats->numLines);
    printf("Header lines:  %10" PRIu32 "\n", stats->numHeaderLines);
    printf("s lines:       %10" PRIu32 "\n", stats->numSeqLines);
    printf("e lines:       %10" PRIu32 "\n", stats->numELines);
    printf("i lines:       %10" PRIu32 "\n", stats->numILines);
    printf("q lines:       %10" PRIu32 "\n", stats->numQLines);
    printf("Blank lines:   %10" PRIu32 "\n", stats->numLines - stats->numSeqLines - 
           stats->numELines - stats->numILines - stats->numQLines - 
           stats->numHeaderLines - stats->numCommentLines);
    printf("Comment lines: %10" PRIu32 "\n", stats->numCommentLines);
    printf("Sequence chars: %12" PRIu32 " (%6.2f%%)\n", stats->numSeqCharacters, 
           100.0 * (double) stats->numSeqCharacters / (stats->numSeqCharacters + stats->numGapCharacters));
    printf("Gap chars:      %12" PRIu32 " (%6.2f%%)\n", stats->numGapCharacters,
           100.0 * (double) stats->numGapCharacters / (stats->numSeqCharacters + stats->numGapCharacters));
    printf("Blocks:                 %10" PRIu32 "\n", stats->numBlocks);
    printf("Ave block area:         %10.2f\n", (double) stats->sumBlockArea / stats->numBlocks);
    printf("Max block area:         %10" PRIu32 "\n", stats->maxBlockArea);
    printf("Ave seq field length:   %10.2f\n", (double) stats->sumSeqField / stats->numSeqLines);
    printf("Max seq field length:   %10" PRIu32 "\n", stats->maxSeqField);
    printf("Ave seq count in block: %10.2f\n", (double) stats->sumNumSpeciesInBlock / stats->numBlocks);
    printf("Max seq count in block: %10" PRIu32 "\n", stats->maxNumSpeciesInBlock);
    printf("%" PRIi32 " unique sequences, ordered by # bases present:\n", stHash_size(stats->seqHash));
    reportHash(stats->seqHash);
    printf("\n");
}
int main(int argc, char **argv) {
    char *maf = NULL;
    parseArgs(argc, argv, &maf);
    mafFileApi_t *mfa = maf_newMfa(maf, "r");
    stats_t *stats = stats_create(maf);

    recordStats(mfa, stats);
    reportStats(stats);

    // clean up
    free(maf);
    maf_destroyMfa(mfa);
    stats_destroy(stats);
    return(EXIT_SUCCESS);
}
