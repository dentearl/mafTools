#include "common.h"
#include "options.h"
#include "genome.h"
#include "malnSet.h"
#include "malnJoinWithinSet.h"
#include "malnJoinSets.h"
#include "malnMergeComps.h"
#include "sonLibETree.h"
#include "malnMultiParents.h"
#include <limits.h>

/* command line option specifications */
static struct optionSpec optionSpecs[] = {
    {"treelessRoot", OPTION_STRING},
    {"help", OPTION_BOOLEAN},
    {NULL, 0}
};

static char *usageMsg =
    "mafJoin [options] inMaf droppedOut\n"
    "\n"
    "Options:\n"
    "  -help\n"
    "  -multiParentDropped=file - write regions were one copy was drop due to\n"
    "   having multiple parents.\n"
    "  -inMafDump=file - dump info about input inMaf after merging dups\n"
    "\n";


/* usage msg and exit */
static void usage(char *msg) {
    errAbort("Error: %s\n%s", msg, usageMsg);
}

/* check a maf */
static void mafMultiParentCheck(char *inMaf, char *multiParentDroppedFile, char *treelessRootName, char *inMafDump) {
    struct Genomes *genomes = genomesNew();
    struct Genome *treelessRootGenome = (treelessRootName != NULL) ? genomesObtainGenome(genomes, treelessRootName) : NULL;
    FILE *dropLogFh = malnMultiParents_openResolveDropLog(multiParentDroppedFile);
    struct malnSet *malnSet = malnSet_constructFromMaf(genomes, inMaf, 10000, 0.1, treelessRootGenome);
    int numResolved = malnMultiParents_resolve(malnSet, dropLogFh);
    if (inMafDump != NULL) {
        malnSet_dumpFile(malnSet, inMafDump, "maf");
    }
    if (numResolved > 0) {
        errAbort("Error: found %d multi-parent instances", numResolved);
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

    mafMultiParentCheck(argv[1], argv[2], optionVal("treelessRoot", NULL), optionVal("inMafDump", NULL));
    return 0;
}
