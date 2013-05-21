/* 
 * Copyright (C) 2013 by 
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
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include "common.h"
#include "CuTest.h"
#include "sharedMaf.h"
#include "sonLib.h"
#include "bioioC.h" // benLine()
#include "mafToFastaStitcher.h"
#include "mafToFastaStitcherAPI.h"
#include "buildVersion.h"

options_t* options_construct(void) {
    options_t *o = (options_t*) st_malloc(sizeof(*o));
    o->maf = NULL;
    o->seqs = NULL;
    o->outMfa = NULL;
    o->outMaf = NULL;
    o->reference = NULL;
    o->breakpointPenalty = 0;
    o->interstitialSequence = 0;
    return o;
}
void destroyOptions(options_t *o) {
    if (o->maf != NULL) {
        free(o->maf);
        o->maf = NULL;
    }
    if (o->seqs != NULL) {
        free(o->seqs);
        o->seqs = NULL;
    }
    if (o->outMfa != NULL) {
        free(o->outMfa);
        o->outMfa = NULL;
    }
    if (o->outMaf != NULL) {
        free(o->outMaf);
        o->outMaf = NULL;
    }
    if (o->reference != NULL) {
        free(o->reference);
        o->reference = NULL;
    }
    free(o);
    o = NULL;
}
mtfseq_t* newMtfseq(uint64_t length) {
    assert(length > 0);
    mtfseq_t *mtfs = (mtfseq_t *) st_malloc(sizeof(*mtfs));
    mtfs->seq = (char *) st_malloc(length);
    mtfs->memLength = length;
    mtfs->index = 0;
    mtfs->seq[mtfs->index] = '\0';
    return mtfs;
}
void resizeMtfseq(mtfseq_t *m) {
    // double the size of the mtfseq_t, updating all members
    uint64_t n = 2 << 16;
    if (n < m->memLength) {
        n = m->memLength;
    }
    n *= 2;
    // de_debug("Doubling the size of mtfseq_t->seq, from %"PRIu64" to %"PRIu64"\n", m->memLength, n);
    m->memLength = n;
    char *new = (char*) st_malloc(n);
    new[0] = '\0';
    strcpy(new, m->seq);
    free(m->seq);
    m->seq = new;
}
void resizeRowSequence(row_t *r) {
    // double the size of the mtfseq_t, updating all members
    uint64_t n = 2 << 16;
    if (n < r->memLength) {
        n = r->memLength;
    }
    n *= 2;
    r->memLength = n;
    char *new = (char *) st_malloc(n);
    new[0] = '\0';
    strcpy(new, r->sequence);
    free(r->sequence);
    r->sequence = new;
}
void destroyMtfseq(void *p) {
    // extra casting due to function being called by stHash destructor
    free(((mtfseq_t *)p)->seq);
    free(p);
}
row_t* newRow(uint64_t n) {
    // n should be a power of two
    row_t *r = (row_t *) st_malloc(sizeof(*r));
    r->name = NULL;
    r->multipleNames = false;
    r->start = 0;
    r->length = 0;
    r->prevRightPos = 0;
    r->prevName = NULL;
    r->strand = '0';
    r->prevStrand = '0';
    r->sourceLength = 0;
    r->memLength = n;
    r->sequence = (char*) st_malloc(n);
    r->sequence[0] = '\0';
    r->index = 0;
    return r;
}
void destroyRow(void *row) {
    // extra casting due to function being called by stHash destructor
    free(((row_t *)row)->name);
    free(((row_t *)row)->prevName);
    free(((row_t *)row)->sequence);
    free(row);
}
row_t* mafLineToRow(mafLine_t *ml) {
    // take a mafLine_t pointer and turn it into a valid row_t pointer
    row_t *r = newRow(nearestTwo(maf_mafLine_getSequenceFieldLength(ml)));
    assert(r->name == NULL);
    assert(r->prevName == NULL);
    assert(r->sequence[0] == '\0');
    r->name = stString_copy(maf_mafLine_getSpecies(ml));
    r->prevName = stString_copy(maf_mafLine_getSpecies(ml));
    row_copyIn(r, maf_mafLine_getSequence(ml)); // copy in sequence
    r->start = maf_mafLine_getStart(ml);
    r->length = maf_mafLine_getLength(ml);
    r->prevRightPos = maf_mafLine_getStart(ml) + maf_mafLine_getLength(ml) - 1;
    r->strand = maf_mafLine_getStrand(ml);
    r->prevStrand = r->strand;
    r->sourceLength = maf_mafLine_getSourceLength(ml);
    return r;
}
stHash* createSequenceHash(char *fastas) {
    unsigned n = 1 + countChar(fastas, ',');
    char **fastaArray = extractSubStrings(fastas, n, ',');
    stHash *sequenceHash = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, free, destroyMtfseq);
    for (unsigned i = 0; i < n; ++i) {
        de_verbose("Reading fasta %s\n", fastaArray[i]);
        addSequencesToHash(sequenceHash, fastaArray[i]);
        free(fastaArray[i]);
    }
    free(fastaArray);
    return sequenceHash;
}
stHash* mafBlockToBlockHash(mafBlock_t *mb, stList *orderList) {
    // create an alignment hash, keyed by species names, valued by row_t pointers
    stHash *bh = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, free, destroyRow);
    mafLine_t *ml = maf_mafBlock_getHeadLine(mb);
    row_t *dummy = NULL;
    char *name = NULL;
    row_t *row = NULL;
    while (ml != NULL) {
        if (maf_mafLine_getType(ml) != 's') {
            ml = maf_mafLine_getNext(ml);
            continue;
        }
        name = copySpeciesName(maf_mafLine_getSpecies(ml));
        dummy = stHash_search(bh, name);
        if (dummy != NULL) {
            fprintf(stderr, "Error, multiple instances of the same species, %s, were detected "
                    "in a maf block near line %" PRIu64 "\n", name, maf_mafBlock_getLineNumber(mb));
            stHash_destruct(bh);
            free(name);
            exit(EXIT_FAILURE);
        }
        row = mafLineToRow(ml);
        stHash_insert(bh, stString_copy(name), row);
        stList_append(orderList, stString_copy(name));
        free(name);
        ml = maf_mafLine_getNext(ml);
    }
    return bh;
}
void seq_copyIn(mtfseq_t *mtfs, char *src) {
    // copy src into mtfseq_t->seq starting at ->index. 
    unsigned n = strlen(src);
    while (mtfs->index + n + 1 >= mtfs->memLength) {
        resizeMtfseq(mtfs);
    }
    for (unsigned i = 0; i < n; ++i) {
        // copy the new sequence into ->seq at the index;
        mtfs->seq[(mtfs->index)++] = src[i];
    }
    mtfs->seq[mtfs->index] = '\0';
    // printf("done inside copyIn, memLength: %"PRIu64" index: %"PRIu64"\n", 
    //        (*mtfs)->memLength, (*mtfs)->index);
}
void row_copyIn(row_t *row, char *src) {
    // copy src into row_t starting with ->index. 
    unsigned n = strlen(src);
    extendSequence(row, n + 1);
    for (unsigned i = 0; i < n; ++i) {
        // copy the new sequence into ->sequence at the index;
        row->sequence[(row->index)++] = src[i];
    }
    row->sequence[row->index] = '\0';
}
void addSequencesToHash(stHash *hash, char *filename) {
    // add sequences from a fasta file into the seqHash containing mtfseq_t values
    FILE *ifp = de_fopen(filename, "r");
    int64_t n = kMaxStringLength;
    char *line = (char*) st_malloc(n);
    char *copy = NULL, *copyHead = NULL;
    char *name = NULL;
    mtfseq_t *mtfs = NULL;
    // uint64_t lineNumber = 1;
    while (benLine(&line, &n, ifp) != -1) {
        // if ((lineNumber++ % 25) == 0) {
        //     de_debug("Reading %s line number %"PRIu64"...\n", filename, lineNumber);
        // }
        if (copyHead != NULL) {
            free(copyHead);
        }
        copy = stString_copy(line);
        copyHead = copy;
        if (copy[0] == '>') {
            // sequence header
            if (name != NULL) {
                // record previous sequence before moving on
                de_debug("    ....Done\n");
                stHash_insert(hash, stString_copy(name), mtfs);
                free(name);
                mtfs = NULL;
            }
            name = stString_getNextWord(&copy);
            if (strlen(name) < 2 && name[0] == '>') {
                // chuck the > as a name, we want an actual sequence name
                free(name);
                name = stString_getNextWord(&copy);
                de_debug("Reading sequence %s from %s\n", name, filename);
            }
            if (name[0] == '>') {
                // we don't want the > at the start of a name, get rid of it
                char *tmp = name;
                name = stString_copy(name + 1);
                de_debug("Reading sequence %s from %s\n", name, filename);
                free(tmp);
            }
            mtfs = newMtfseq(2 << 16);
        } else {
            // sequence line
            seq_copyIn(mtfs, line);
        }
    }
    if (name != NULL) {
        // record last sequence
        de_debug("    ....Done\n");
        stHash_insert(hash, stString_copy(name), mtfs);
        free(name);
    }
    if (copy != NULL) {
        free(copy);
    }
    fclose(ifp);
    free(line);
    de_verbose("Finished reading fasta %s\n", filename);
}
void reportSequenceHash(stHash *hash) {
    stHashIterator *hit = stHash_getIterator(hash);
    char *key = NULL;
    printf("Sequence Hash:\n");
    while ((key = stHash_getNext(hit)) != NULL) {
        printf("found key    : %s\n", key);
        printf("    memLength: %" PRIu64 "\n", ((mtfseq_t *)stHash_search(hash, key))->memLength);
        printf("        index: %" PRIu64 "\n", ((mtfseq_t *)stHash_search(hash, key))->index);
        printf("               %s\n", ((mtfseq_t *)stHash_search(hash, key))->seq);
    }
    stHash_destructIterator(hit);
}
void extendSequence(row_t *r, uint64_t n) {
    // ensure there is enough room to write an additional `n' chars to ->sequence
    while (r->index + n + 1 >= r->memLength) {
        resizeRowSequence(r);
    }
}
void penalize(stHash *hash, char *name, uint64_t n) {
    // walk the hash looking for a row_t with ->name equal to input *name,
    // penalize that sequence
    stHashIterator *hit = stHash_getIterator(hash);
    char *key = NULL;
    row_t *row = NULL;
    char fill = '-';
    char *sppName = copySpeciesName(name);
    char *rowSppName = NULL;
    while ((key = stHash_getNext(hit)) != NULL) {
        row = stHash_search(hash, key);
        extendSequence(row, n); // make space
        fill = '-';
        rowSppName = copySpeciesName(row->name);
        if (strcmp(rowSppName, sppName) == 0) {
            // printf("   laying the hurt down on %20s: ", rowSppName);
            // penalize this row
            fill = 'N';
            row->length += n;
            row->prevRightPos += n;
        } else {
            // printf("   just going to gap       %20s: ", rowSppName);
        }
        for (uint64_t i = 0; i < n; ++i) {
            row->sequence[row->index] = fill;
            ++(row->index);
        }
        row->sequence[row->index] = '\0';
        // printf("%s\n", row->sequence);
        free(rowSppName);
        rowSppName = NULL;
    }
    free(sppName);
    stHash_destructIterator(hit);
}
void interstitialInsert(stHash *alignHash, stHash *seqHash, char *name, uint64_t pos, 
                        char strand, uint64_t n) {
    // for row *name, insert the correct sequence on the end. pad other sequences with gap characters.
    stHashIterator *hit = stHash_getIterator(alignHash);
    char *key = NULL;
    row_t *row = NULL;
    mtfseq_t *mtfs = NULL;
    char *seq = NULL;
    while ((key = stHash_getNext(hit)) != NULL) {
        row = stHash_search(alignHash, key);
        extendSequence(row, n); // make space
        if (strcmp(row->name, name) == 0) {
            // printf("    going to interstitialize %20s: ", row->name);
            // insert into this row
            row->length += n;
            row->prevRightPos += n;
            mtfs = stHash_search(seqHash, name);
            if (mtfs == NULL) {
                fprintf(stderr, "Error, unable to locate sequnce %s in the sequence hash. "
                        "Check your input fasta files.\n", name);
                exit(EXIT_FAILURE);
            }
            seq = extractSubSequence(mtfs, strand, pos, n);
            for (uint64_t i = 0; i < n; ++i) {
                row->sequence[row->index] = seq[i];
                ++(row->index);
            }
            // printf("%s\n", row->sequence);
            free(seq);
        } else {
            // printf("    just going to gap        %20s: ", row->name);
            // these aren't the droids you're looking for, write some gaps instead
            for (uint64_t i = 0; i < n; ++i) {
                row->sequence[row->index] = '-';
                ++(row->index);
            }
            // printf("%s\n", row->sequence);
        }
        row->sequence[row->index] = '\0';
    }
    stHash_destructIterator(hit);
}
char* extractSubSequence(mtfseq_t *mtfs, char strand, uint64_t pos, uint64_t n) {
    // make a copy of a region of a mtfseq_t structure, performing coordinate transform and 
    // reverse complementation if the strand is -
    char *seq = (char*) st_malloc(n + 1);
    if (strand == '+') {
        for (uint64_t i = 0; i < n; ++i) {
            seq[i] = mtfs->seq[pos + i];
        }
        seq[n] = '\0';
    } else {
        uint64_t p = mtfs->index - pos - n;
        for (uint64_t i = 0; i < n; ++i) {
            seq[i] = mtfs->seq[p + i];
        }
        seq[n] = '\0';
        reverseComplementSequence(seq, n);
    }
    return seq;
}
void addMafLineToRow(row_t *row, mafLine_t *ml) {
    // given a row_t and a mafLine_t, add the information from the mafLine_t to the row_t
    char *seq = maf_mafLine_getSequence(ml);
    size_t n = strlen(seq);
    extendSequence(row, n + 1);
    row_copyIn(row, seq);
    free(row->prevName);
    row->prevName = stString_copy(maf_mafLine_getSpecies(ml));
    row->prevRightPos = maf_mafLine_getStart(ml) + maf_mafLine_getLength(ml) - 1;
    row->prevStrand = maf_mafLine_getStrand(ml);
    row->length += maf_mafLine_getLength(ml);
    if (row->multipleNames) {
        row->sourceLength = row->length;
    }
}
void prependGaps(row_t *r, uint64_t n) {
    // add `n' many gap characters, '-', to the begining of row_t *r
    extendSequence(r, n + 1);
    char *new = (char*) st_malloc(r->memLength);
    new[0] = '\0';
    for (uint64_t i = 0; i < n; ++i) {
        new[i] = '-';
    }
    new[n] = '\0';
    if (r->sequence != NULL) {
        strcat(new, r->sequence);
        free(r->sequence);
    }
    r->sequence = new;
    r->index += n;
}
uint64_t nearestTwo(uint64_t n) {
    // return the smallest power of two result that is greater than n
    uint64_t t = 64;
    while (t < n) {
        t *= 2;
    }
    return t;
}
void addMafBlockToRowHash(stHash *alignHash, stHash *seqHash, stList *orderList, 
                          mafBlock_t *mb, options_t *options) {
    mafLine_t *ml = maf_mafBlock_getHeadLine(mb);
    row_t *r = NULL;
    char *seqName = NULL; // full name
    char *sppName = NULL; // just the bit until the first '.' character
    // first loop, penalize and interstitialize the existing hash as necessary:
    while (ml != NULL) {
        if (maf_mafLine_getType(ml) != 's') {
            ml = maf_mafLine_getNext(ml);
            continue;
        }
        seqName = maf_mafLine_getSpecies(ml);
        sppName = copySpeciesName(seqName);
        if (alignHash == NULL) {
            alignHash = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, free, destroyMtfseq);
        }
        r = stHash_search(alignHash, sppName);
        // printf("observed sequence %s,\n", seqName);
        uint64_t n;
        if (r == NULL) {
            // printf("sequence %s is novel, adding to hash\n", seqName);
            // add this row to the hash
            stHashIterator *hit = stHash_getIterator(alignHash);
            char *key = stHash_getNext(hit);
            if (key != NULL) {
                // if key is not null then figure out how many gap chars to put in front of this sequence
                row_t *value = stHash_search(alignHash, key);
                n = value->index;
            } else {
                n = 0;
            }
            stHash_destructIterator(hit);
            if (n > 0) {
                r = newRow(nearestTwo(n));
                // printf("prepend some gaps (%" PRIu64 ") on %s\n", n, seqName);
                prependGaps(r, n);
            } else { 
                r = newRow(2 << 7); // 256 seems like an okay starting point
            }
            // empty row_t structure, populate it:
            assert(r->name == NULL);
            assert(r->prevName == NULL);
            r->name = stString_copy(seqName);
            r->prevName = stString_copy(seqName);
            r->start = maf_mafLine_getStart(ml);
            r->sourceLength = maf_mafLine_getSourceLength(ml);
            r->strand = maf_mafLine_getStrand(ml);
            r->prevStrand = r->strand;
            stHash_insert(alignHash, stString_copy(sppName), r);
            de_debug("inserted new row into hash: %s %s\n", sppName, r->name);
            stList_append(orderList, stString_copy(sppName));
        } else {
            // row already in hash
            if ((options->reference != NULL) && (strcmp(options->reference, r->name) == 0)) {
                // extend the reference if necessary
                if (r->prevRightPos + 1 < maf_mafLine_getStart(ml)) {
                    interstitialInsert(alignHash, seqHash, seqName, 
                                       r->prevRightPos + 1,
                                       maf_mafLine_getStrand(ml), 
                                       maf_mafLine_getStart(ml) - r->prevRightPos - 1);
                }
            }else if (r->prevStrand != maf_mafLine_getStrand(ml)) {
                // different strands is a breakpoint
                // printf("penalize 0 (%"PRIu64") %s\n", options->breakpointPenalty, seqName);
                penalize(alignHash, seqName, options->breakpointPenalty);
                r->strand = '*';
            } else if (strcmp(r->prevName, seqName) != 0) {
                // different names implies diff. chromosomes, is a breakpoint
                // printf("penalize 1 (%"PRIu64") %s\n", options->breakpointPenalty, seqName);
                penalize(alignHash, seqName, options->breakpointPenalty);
                r->strand = '*';
                r->multipleNames = true;
                de_debug("penalizing %s because name has changed to %s.\n", r->prevName, seqName);
                free(r->name);
                r->name = copySpeciesName(seqName);
                free(r->prevName);
                r->prevName = stString_copy(seqName);
                r->sourceLength = r->length;
                r->start = 0;
            } else if (r->prevRightPos + options->interstitialSequence < maf_mafLine_getStart(ml)) {
                // same chromosome but beyond the accepted interstitial range, breakpoint
                // printf("penalize 2 (%"PRIu64") %s\n", options->breakpointPenalty, seqName);
                penalize(alignHash, seqName, options->breakpointPenalty);
                r->multipleNames = true;
                de_debug("penalizing %s because of interstitial distance. prev: %"PRIu64" current:%"PRIu64", "
                         "line nubber: %"PRIu64"\n", 
                         r->prevName, r->prevRightPos, maf_mafLine_getStart(ml), maf_mafLine_getLineNumber(ml));
                free(r->name);
                r->name = copySpeciesName(seqName);
                free(r->prevName);
                r->prevName = stString_copy(seqName);
                r->start = 0;
                r->sourceLength = r->length;
            } else if ((r->prevRightPos + 1 < maf_mafLine_getStart(ml)) && 
                       (r->prevRightPos + 1 + options->interstitialSequence >= maf_mafLine_getStart(ml))) {
                // same chromosome and within the accepted interstitial range, insert sequence
                // printf("interstitialize %s\n", seqName);
                interstitialInsert(alignHash, seqHash, seqName, 
                                   r->prevRightPos + 1, 
                                   maf_mafLine_getStrand(ml), 
                                   maf_mafLine_getStart(ml) - r->prevRightPos - 1);
            }
        }
        assert(r->name != NULL);
        free(sppName);
        ml = maf_mafLine_getNext(ml);
    }
    // secord loop, append the block's sequence to the hash rows, update the prev* entries
    ml = maf_mafBlock_getHeadLine(mb);
    stSet *presentSet = stSet_construct3(stHash_stringKey, stHash_stringEqualKey, free);
    uint64_t currentIndex = 0;
    while (ml != NULL) {
        if (maf_mafLine_getType(ml) != 's') {
            ml = maf_mafLine_getNext(ml);
            continue;
        }
        seqName = maf_mafLine_getSpecies(ml);
        sppName = copySpeciesName(seqName);
        r = stHash_search(alignHash, sppName);
        stSet_insert(presentSet, stString_copy(sppName));
        assert(r != NULL);
        addMafLineToRow(r, ml);
        currentIndex = r->index;
        free(sppName);
        ml = maf_mafLine_getNext(ml);
    }
    stHashIterator *hit = stHash_getIterator(alignHash);
    char *key = NULL;
    // printf(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> RESULTS <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n");
    // third loop, gap out all sequences that are in the alignHash but were not in the mafBlock
    while ((key = stHash_getNext(hit)) != NULL) {
        if (stSet_search(presentSet, key) == NULL) {
            // this species was not presentSet in the current mafBlock, gap out the sequence.
            r = stHash_search(alignHash, key);
            extendSequence(r, 1 + currentIndex - r->index);
            for (uint64_t i = r->index; i < currentIndex; ++i) {
                r->sequence[i] = '-';
            }
            r->sequence[currentIndex] = '\0';
            r->index = currentIndex;
        }
        r = stHash_search(alignHash, key);
        // printf("       result              %20s: %s\n", r->name, r->sequence);
    }
    stSet_destruct(presentSet);
    stHash_destructIterator(hit);
}
void buildAlignmentHash(mafFileApi_t *mfapi, stHash *alignmentHash, stHash *sequenceHash, 
                        stList *rowOrder, options_t *options) {
    mafBlock_t *mb = NULL;
    assert(alignmentHash != NULL);
    assert(sequenceHash != NULL);
    assert(rowOrder != NULL);
    while ((mb = maf_readBlock(mfapi)) != NULL) {
        // printf("working on this block:\n");
        // maf_mafBlock_print(mb);
        addMafBlockToRowHash(alignmentHash, sequenceHash, rowOrder, mb, options);
        // printf("   ...block done.\n\n");
        maf_destroyMafBlockList(mb);
    }
}
void writeFastaOut(stHash *alignmentHash, stList *rowOrder, options_t *options) {
    row_t *r = NULL;
    FILE *fa = de_fopen(options->outMfa, "w");
    // printf("printing fasta out!\n");
    for (int64_t i = 0; i < stList_length(rowOrder); ++i) {
        r = stHash_search(alignmentHash, stList_get(rowOrder, i));
        assert(r != NULL);
        fprintf(fa, "> %s\n", r->name);
        // printf("> %s\n%s\n", r->name, r->sequence);
        for (uint64_t j = 0; j < r->index; ++j) {
            fprintf(fa, "%c", r->sequence[j]);
            if (!((j + 1) % 50) && j != r->index - 1) {
                fprintf(fa, "\n");
            }
        }
        fprintf(fa, "\n");
    }
    fclose(fa);
}
void writeMafOut(stHash *alignmentHash, stList *rowOrder, options_t *options) {
    row_t *r = NULL;
    FILE *maf = de_fopen(options->outMaf, "w");
    // fprintf(stderr, "printing Maf out!\n");
    uint64_t maxName = 1, maxStart = 1, maxLen = 1, maxSource = 1;
    char fmtName[10] = "\0", fmtStart[32] = "\0", fmtLen[32] = "\0", fmtSource[32] = "\0", *fmtLine = NULL;
    fprintf(maf, "##maf version=1\n\n");
    if (stList_length(rowOrder) == 0) {
        // There's nothing to write out.
        fclose(maf);
        return;
    }
    for (int64_t i = 0; i < stList_length(rowOrder); ++i) {
        // first loop, get formating correct
        r = stHash_search(alignmentHash, stList_get(rowOrder, i));
        // printf("collecting info on ");
        // printf("%s\n", (char*)stList_get(rowOrder, i));
        assert(r != NULL);
        if (maxName < strlen(r->name)) {
            maxName = strlen(r->name);
        }
        if (maxStart < r->start) {
            maxStart = r->start;
        }
        if (maxLen < r->length) {
            maxLen = r->length;
        }
        if (maxSource < r->sourceLength) {
            maxSource = r->sourceLength;
        }
    }
    fmtLine = (char*) st_malloc(4 + (int)log10(maxName) + 3 + (int)log10(maxStart) + 3 + (int)log10(maxLen) + 
                                3 + (int)log10(maxSource) + 4 + r->index);
    fmtLine[0] = '\0';
    sprintf(fmtName, " %%-%" PRIu64 "s", maxName + 2);
    sprintf(fmtStart, " %%%d" PRIu64, (int)log10(maxStart) + 2);
    sprintf(fmtLen, " %%%d" PRIu64, (int)log10(maxLen) + 2);
    sprintf(fmtSource, " %%%d" PRIu64, (int)log10(maxSource) + 2);
    strcat(fmtLine, "s");
    strcat(fmtLine, fmtName);
    strcat(fmtLine, fmtStart);
    strcat(fmtLine, fmtLen);
    strcat(fmtLine, " %c");
    strcat(fmtLine, fmtSource);
    strcat(fmtLine, " %s\n");
    char strand;
    fprintf(maf, "a stitched=true\n");
    // printf("\nPrinting actual block now, fmtline: %s", fmtLine);
    for (int64_t i = 0; i < stList_length(rowOrder); ++i) {
        // second loop, print!
        r = stHash_search(alignmentHash, stList_get(rowOrder, i));
        assert(r != NULL);
        if (r->strand == '*') {
            strand = '+';
        } else {
            strand = r->strand;
        }
        fprintf(maf, fmtLine, r->name, r->start, r->length, strand, r->sourceLength, r->sequence);
    }
    fprintf(maf, "\n");
    free(fmtLine);
    fclose(maf);
}
