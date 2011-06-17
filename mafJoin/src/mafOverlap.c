#include "common.h"
#include "options.h"
#include "basicBed.h"
#include "jkmaf.h"
#include "dnautil.h"
#include <stdbool.h>

/* command line option specifications */
static struct optionSpec optionSpecs[] = {
    {"overlapComps", OPTION_BOOLEAN},
    {NULL, 0}
};

static char *usageMsg =
    "mafOverlap [options] inMaf selectBed outMaf\n"
    "\n"
    "Select records from inMaf with any overlap with records in selectbed.\n"
    "The chrom column in bed should be org.chrom.  If the BED record includes\n"
    "a strand, both the strand and it's reverse complement on the opposite\n"
    "strand are checked.   Without a strand it the BED, the BED coordinates will\n"
    "be tried as both positive and negative strand coordinates.\n"
    "\n"
    "  -overlapComps - drop components that don't overlap the regions\n";

static boolean overlapComps = FALSE;

/* usage msg and exit */
static void usage(char *msg) {
    errAbort("Error: %s\n%s", msg, usageMsg);
}

/* see if a bed overlaps a component */
static bool isOverlapped(struct bed *bed, struct mafComp *comp) {
    if (sameString(bed->chrom, comp->src)) {
        int compStart = comp->start;
        int compEnd = compStart + comp->size;
        if (bed->strand[0] != '\0') {
            if (bed->strand[0] != comp->strand) {
                reverseIntRange(&compStart, &compEnd, comp->srcSize);
            }
            return (positiveRangeIntersection(bed->chromStart, bed->chromEnd, compStart, compEnd) > 0);
        } else {
            // try both ways
            if (positiveRangeIntersection(bed->chromStart, bed->chromEnd, compStart, compEnd) > 0) {
                return true;
            }
            reverseIntRange(&compStart, &compEnd, comp->srcSize);
            if (positiveRangeIntersection(bed->chromStart, bed->chromEnd, compStart, compEnd) > 0) {
                return true;
            }
        }
    }
    return false;
}

/* see if any beds overlap this component */
static bool anyOverlaps(struct bed *selectBeds, struct mafComp *comp) {
    for (struct bed *bed = selectBeds; bed != NULL; bed = bed->next) {
        if (isOverlapped(bed, comp)) {
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

/* only write components that overlap */
static void writeOverlapComps(struct bed *selectBeds, struct mafAli *blk, FILE *outMafFh) {
    struct mafComp *comp, *holdComps = NULL, *useComps = NULL;
    while ((comp = slPopHead(&blk->components)) != NULL) {
        slAddHead((anyOverlaps(selectBeds, comp) ? &useComps : &holdComps), comp);
    }
    slReverse(&useComps);
    if (holdComps != NULL) {
        blk->components = useComps;
    }
    mafWrite(outMafFh, blk);
    if (holdComps != NULL) {
        blk->components = slCat(useComps, holdComps);
    }
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
            if (overlapComps) {
                writeOverlapComps(selectBeds, blk, outMafFh);
            } else {
                mafWrite(outMafFh, blk);
            }
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
    overlapComps = optionExists("overlapComps");

    mafOverlap(argv[1], argv[2], argv[3]);
    return 0;
}
