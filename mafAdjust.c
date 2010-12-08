#include "common.h"
#include "options.h"
#include "genome.h"
#include "mafIO.h"
#include "malnBlk.h"
#include "mafTree.h"
#include "jkmaf.h"
#include <limits.h>

/* command line option specifications */
static struct optionSpec optionSpecs[] = {
    {"maxBlkWidth", OPTION_INT},
    {"dumpFile", OPTION_STRING},
    {"help", OPTION_BOOLEAN},
    {NULL, 0}
};

static char *usageMsg =
    "mafAdjust [options] inMaf outMaf\n"
    "\n"
    "Options:\n"
    "  -help\n"
    "  -maxBlkWidth=n - Limit size of MAF blocks to this value.\n"
    "  -dumpDir=dir - dump info about MAFs at various points during the\n"
    "   process to files in this directory.\n"
    "\n";


/* usage msg and exit */
static void usage(char *msg) {
    errAbort("Error: %s\n%s", msg, usageMsg);
}

/* make adjustment to a block */
static void adjustBlk(struct malnBlk *inBlk, FILE *outMafFh, int maxBlkWidth) {
    int alnStart = 0;
    while (alnStart < inBlk->alnWidth) {
        int alnEnd = alnStart + maxBlkWidth;
        if (alnEnd > inBlk->alnWidth) {
            alnEnd = inBlk->alnWidth;
        }
        struct malnBlk *outBlk = malnBlk_constructSubrange(inBlk, alnStart, alnEnd);
        slReverse(&outBlk->comps);  // FIXME: hack for removing trees
        mafIO_malnBlkWrite(outBlk, outMafFh);
        malnBlk_destruct(outBlk);
        alnStart = alnEnd;
    }
}

/* process one maf block */
static void mafAliProcess(struct Genomes *genomes, int maxBlkWidth, struct mafAli *mafAli, FILE *outMafFh, FILE *dumpFh) {
    struct malnBlk *inBlk = mafIO_malnBlkRead(genomes, mafAli, 0.0, NULL, false);
    if (dumpFh != NULL) {
    }
    // FIXME: tmp hack to remove trees until we have accurate trees
    mafTree_destruct(inBlk->mTree);
    inBlk->mTree = NULL;
    adjustBlk(inBlk, outMafFh, maxBlkWidth);
    malnBlk_destruct(inBlk);
}

/* make adjustments to mash mafs */
static void mafAdjust(char *inMaf, char *outMaf, int maxBlkWidth, char *dumpFile) {
    struct Genomes *genomes = genomesNew();
    struct mafFile *inMafFh = mafOpen(inMaf);
    FILE *outMafFh = mustOpen(outMaf, "w");
    mafWriteStart(outMafFh, NULL);
    FILE *dumpFh = (dumpFile != NULL) ? mustOpen(dumpFile, "w") : NULL;
    struct mafAli *mafAli;
    while ((mafAli = mafNext(inMafFh)) != NULL) {
        mafAliProcess(genomes, maxBlkWidth, mafAli, outMafFh, dumpFh);
        mafAliFree(&mafAli);
    }

    carefulClose(&dumpFh);
    carefulClose(&outMafFh);
    mafFileFree(&inMafFh);
    genomesFree(genomes);
}

/* Process command line. */
int main(int argc, char *argv[]) {
    optionInit(&argc, argv, optionSpecs);
    if (optionExists("help")) {
        usage("Usage:");
    }
    if (argc != 3)  {
        usage("Error: wrong number of arguments");
    }

    mafAdjust(argv[1], argv[2],
              optionInt("maxBlkWidth", INT_MAX),
              optionVal("dumpFile", NULL));
    return 0;
}
