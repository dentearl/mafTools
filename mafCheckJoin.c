#include "common.h"
#include "options.h"
#include "genome.h"
#include "sonLibHash.h"
#include "jkmaf.h"
#include "dnautil.h"
#include <stdbool.h>

// FIXME: this is too simplistic due to indels, really needed
// to check that a bases continues to be aligned to the
// same sequences in the output, but that is very expensive

/* command line option specifications */
static struct optionSpec optionSpecs[] = {
    {NULL, 0}
};

static char *usageMsg =
    "mafCheckJoin [options] mafIn1 mafIn2 mafJoined outBed\n"
    "\n"
    "Compare to MAFs that were joined with the resulting joined MAF.\n"
    "this checks to see if the number of bases a given bases is aligned\n"
    "with remains constant.\n"
    "\n";

/* usage msg and exit */
static void usage(char *msg) {
    errAbort("Error: %s\n%s", msg, usageMsg);
}

/* per-bases stats */
struct Stat {
    bool aligned;  // was this base aligned
};

struct SeqCounts {
    struct Seq *seq;
    struct Stat *stats;  // by base
};

/* constructor */
static struct SeqCounts *seqCounts_construct(struct Seq *seq) {
    struct SeqCounts *sc;
    AllocVar(sc);
    sc->seq = seq;
    sc->stats = needLargeZeroedMem(seq->size * sizeof(struct Stat));
    return sc;
}

/* flag a column as aligned */
static void seqCounts_markAligned(struct SeqCounts *sc, int iSeq) {
    assert((0 <= iSeq) && (iSeq < sc->seq->size));
    sc->stats[iSeq].aligned = true;
}

/* is a column aligned */
static bool seqCounts_isAligned(struct SeqCounts *sc, int iSeq) {
    assert((0 <= iSeq) && (iSeq < sc->seq->size));
    return sc->stats[iSeq].aligned;
}

/* counts over one or more alignments */
struct AlignCounts {
    struct Genomes *genomes;  // all genome objects
    stHash *bySeq;  // mapping of Seq object to SeqCounts
};

/* constructor */
static struct AlignCounts *alignCounts_construct(struct Genomes *genomes) {
    struct AlignCounts *ac;
    AllocVar(ac);
    ac->genomes = genomes;
    ac->bySeq = stHash_construct();
    return ac;
}

/* obtain the a SeqCounts object, constructing if it doesn't exist */
static struct SeqCounts *alignCounts_obtain(struct AlignCounts *ac, struct Seq *seq) {
    struct SeqCounts *sc = stHash_search(ac->bySeq, seq);
    if (sc == NULL) {
        sc = seqCounts_construct(seq);
        stHash_insert(ac->bySeq, seq, sc);
    }
    return sc;
}

/* get the number of sequences */
static int alignCounts_getNumSeqs(struct AlignCounts *ac) {
    return stHash_size(ac->bySeq);
}

/* count a mafAli component */
static void alignCounts_countComp(struct AlignCounts *ac, struct mafAli *blk, struct mafComp *comp) {
    struct Seq *seq = genomesObtainSeqForOrgSeqName(ac->genomes, comp->src, comp->srcSize);
    struct SeqCounts *sc = alignCounts_obtain(ac, seq);
    int chromStart = comp->start, chromEnd = comp->start + comp->size;
    if (comp->strand == '-') {
        reverseIntRange(&chromStart, &chromEnd, comp->srcSize);
    }
    int iSeq = chromStart;
    for (int iCol = 0; iCol < blk->textSize; iCol++) {
        if (isBase(comp->text[iCol])) {
            seqCounts_markAligned(sc, iSeq);
            iSeq++;
        }
    }
}

/* count a mafAli block */
static void alignCounts_countBlk(struct AlignCounts *ac, struct mafAli *blk) {
    for (struct mafComp *comp = blk->components; comp != NULL; comp = comp->next) {
        alignCounts_countComp(ac, blk, comp);
    }
}

/* count all blocks in a maf file */
static void alignCounts_countMaf(struct AlignCounts *ac, char *mafFile) {
    struct mafFile *fh = mafOpen(mafFile);
    struct mafAli *blk;
    while ((blk = mafNext(fh)) != NULL) {
        alignCounts_countBlk(ac, blk);
        mafAliFree(&blk);
    }
    mafFileFree(&fh);
    
}

/* advance to the next position where the aligned states is are not the same */
static int advanceToDiffState(struct SeqCounts *inSc, struct SeqCounts *joinedSc, int iSeq) {
    while ((iSeq < inSc->seq->size) && (seqCounts_isAligned(inSc, iSeq) == seqCounts_isAligned(joinedSc, iSeq))) {
        iSeq++;
    }
    return iSeq;
}

/* advance to the next position were the state changes */
static int advanceToChangedState(struct SeqCounts *inSc, struct SeqCounts *joinedSc, int iSeq) {
    bool inState = seqCounts_isAligned(inSc, iSeq);
    bool joinedState = seqCounts_isAligned(joinedSc, iSeq);
    iSeq++;
    while ((iSeq < inSc->seq->size) && (seqCounts_isAligned(inSc, iSeq) == inState) && (seqCounts_isAligned(joinedSc, iSeq) == joinedState)) {
        iSeq++;
    }
    return iSeq;
}

/* return t for aligned, f for not aligned */
static char stateToChar(bool state) {
    return (state) ? 't' : 'f';
}

/* compare a range of a sequence that has the same counts */
static bool reportSeqRangeDiffs(struct SeqCounts *inSc, struct SeqCounts *joinedSc, int iSeq, int iSeqNext, FILE *outBedFh) {
    bool isOk = true;
    fprintf(outBedFh, "%s\t%d\t%d\t%c/%c\n", inSc->seq->orgSeqName, iSeq, iSeqNext, stateToChar(seqCounts_isAligned(inSc, iSeq)), stateToChar(seqCounts_isAligned(joinedSc, iSeq)));
    return isOk;
}

/* compare counts for a sequence */
static bool compareSeqCounts(struct SeqCounts *inSc, struct SeqCounts *joinedSc, FILE *outBedFh) {
    bool isOk = true;
    int iSeq = 0;
    while (iSeq < inSc->seq->size) {
        iSeq = advanceToDiffState(inSc, joinedSc, iSeq);
        if (iSeq < inSc->seq->size) {
            int iSeqNext = advanceToChangedState(inSc, joinedSc, iSeq);
            reportSeqRangeDiffs(inSc, joinedSc, iSeq, iSeqNext, outBedFh);
            isOk = false;
            iSeq = iSeqNext;
        }
    }
    return isOk;
}

/* compare counts */
static bool compareCounts(struct AlignCounts *inAc, struct AlignCounts *joinedAc, FILE *outBedFh) {
    bool isOk = true;
    if (alignCounts_getNumSeqs(joinedAc) > alignCounts_getNumSeqs(inAc)) {
        errAbort("Error: number of joined sequences (%d) is great than number of input sequences %d; perhaps wrong MAFs were supplied", alignCounts_getNumSeqs(joinedAc), alignCounts_getNumSeqs(inAc));
    }

    stHashIterator *iter = stHash_getIterator(inAc->bySeq);
    struct Seq *seq;
    while ((seq = stHash_getNext(iter)) != NULL) {
        if (!compareSeqCounts(alignCounts_obtain(inAc, seq), alignCounts_obtain(joinedAc, seq), outBedFh)) {
            isOk = false;
        }
    }
    stHash_destructIterator(iter);
    return isOk;
}

/* check joins */
static void mafCheckJoin(char *mafIn1File, char *mafIn2File, char *mafJoinedFile, char *outBedFile) {
    struct Genomes *genomes = genomesNew();
    struct AlignCounts *inCnts = alignCounts_construct(genomes);
    alignCounts_countMaf(inCnts, mafIn1File);
    alignCounts_countMaf(inCnts, mafIn2File);
    struct AlignCounts *joinedCnts = alignCounts_construct(genomes);
    alignCounts_countMaf(joinedCnts, mafJoinedFile);

    FILE *outBedFh = mustOpen(outBedFile, "w");
    int isOk = compareCounts(inCnts, joinedCnts, outBedFh);
    carefulClose(&outBedFh);
    if (!isOk) {
        errAbort("Error: counts don't match");
    }
}

/* Process command line. */
int main(int argc, char *argv[]) {
    optionInit(&argc, argv, optionSpecs);
    if (optionExists("help")) {
        usage("Usage:");
    }
    if (argc != 5)  {
        usage("Error: wrong number of arguments");
    }

    mafCheckJoin(argv[1], argv[2], argv[3], argv[4]);
    return 0;
}
