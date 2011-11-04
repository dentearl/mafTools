#include "common.h"
#include "options.h"
#include "genome.h"
#include "malnSet.h"
#include "malnBlkSet.h"
#include "malnAdjust.h"
#include <limits.h>

/* set to true to cleanup all memory for memory leak checks */
static bool memLeakDebugCleanup = true;

/* command line option specifications */
static struct optionSpec optionSpecs[] = {
    {"maxBlkWidth", OPTION_INT},
    {"dumpDir", OPTION_STRING},
    {"help", OPTION_BOOLEAN},
    {NULL, 0}
};

static char *usageMsg =
    "mafAdjust [options] inMaf outMaf\n"
    "\n"
    "Options:\n"
    "  -help\n"
    "  -maxBlkWidth=n - Limit size of MAF blocks to this value.\n"
    "   Due to current issues with tree, blocks will be cut to include\n"
    "   at least one column of the root.\n"
    "  -dumpDir=dir - dump info about MAFs at various points during the\n"
    "   process to files in this directory.\n"
    "\n";


/* usage msg and exit */
static void usage(char *msg) {
    errAbort("Error: %s\n%s", msg, usageMsg);
}

/* make adjustments to mash mafs */
static void mafAdjust(char *inMaf, char *outMaf, int maxBlkWidth, char *dumpDir) {
    struct Genomes *genomes = genomesNew();
    struct malnSet *malnSetIn = malnSet_constructFromMaf(genomes, inMaf, INT_MAX, 0.0, NULL);
    malnSet_dumpToDir(malnSetIn, dumpDir, "in", "1.input");
    
    struct malnSet *malnSetOut = malnSet_construct(genomes, NULL);

    malnAdjust_adjustSet(malnSetIn, malnSetOut, maxBlkWidth);

    malnSet_dumpToDir(malnSetOut, dumpDir, "out", "1.output");
    malnSet_writeMaf(malnSetOut, outMaf);

    if (memLeakDebugCleanup) {
        malnSet_destruct(malnSetIn);
        malnSet_destruct(malnSetOut);
        genomesFree(genomes);
    }
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
              optionVal("dumpDir", NULL));
    return 0;
}
