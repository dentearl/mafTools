#include "common.h"
#include "options.h"
#include "basicBed.h"
#include "jkmaf.h"
#include "dnautil.h"
#include <stdbool.h>

/* command line option specifications */
static struct optionSpec optionSpecs[] = {
    {NULL, 0}
};

static char *usageMsg =
    "mafOverlap [options] inMaf selectBed outMaf\n"
    "\n"
    "Select records from inMaf with any overlap with records in selectbed.\n"
    "The chrom column in bed should be org.chrom.\n";

/* usage msg and exit */
static void usage(char *msg) {
    errAbort("Error: %s\n%s", msg, usageMsg);
}

/* see if any beds overlap this component */
static bool anyOverlaps(struct bed *selectBeds, struct mafComp *comp) {
    int start1 = comp->start, end1 = comp->start + comp->size;
    int start2 = start1, end2 = end1;
    reverseIntRange(&start1, &end1, comp->srcSize);
    
    for (struct bed *bed = selectBeds; bed != NULL; bed = bed->next) {
        if (sameString(bed->chrom, comp->src)
            && (positiveRangeIntersection(bed->chromStart, bed->chromEnd, start1, end1)
                || positiveRangeIntersection(bed->chromStart, bed->chromEnd, start2, end2))) {
            return true;
        }
    }
    return false;
}

/* determine if a maf block should be kept */
static bool keepMafBlk(struct bed *selectBeds, struct mafAli *blk) {
    for (struct mafComp *comp = blk->components; comp != NULL; comp = comp->next) {
        if (anyOverlaps(selectBeds, comp)) {
            return true;
        }
    }
    return false;
}

/* select overlapping MAF records */
static void mafOverlap(char *inMafFile, char *selectBedFile, char *outMafFile) {
    struct bed *selectBeds = bedLoadAll(selectBedFile);
    struct mafFile *inMafFh = mafOpen(inMafFile);
    FILE *outMafFh = mustOpen(outMafFile, "w");
    mafWriteStart(outMafFh, NULL);
    
    struct mafAli *blk;
    while ((blk = mafNext(inMafFh)) != NULL) {
        if (keepMafBlk(selectBeds, blk)) {
            mafWrite(outMafFh, blk);
        }
        mafAliFree(&blk);
    }

    mafWriteEnd(outMafFh);
    carefulClose(&outMafFh);
    mafFileFree(&inMafFh);
}

/* Process command line. */
int main(int argc, char *argv[]) {
    optionInit(&argc, argv, optionSpecs);
    if (optionExists("help")) {
        usage("Usage:");
    }
    if (argc != 4)  {
        usage("Error: wrong number of arguments");
    }

    mafOverlap(argv[1], argv[2], argv[3]);
    return 0;
}
