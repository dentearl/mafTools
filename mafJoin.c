#include "common.h"
#include "options.h"
#include "genome.h"
#include "malnSet.h"
#include "malnJoinDups.h"
#include "malnJoinSets.h"
#include "sonLibETree.h"
#include "malnMultiParents.h"

/*
 * Notes: 
 *  - this would probably be simpler (and for sure more efficient) if
 *    it was done with block based rather than column based merges.
 */

/* set to true to cleanup all memory for memory leak checks */
static bool memLeakDebugCleanup = true;

/* command line option specifications */
static struct optionSpec optionSpecs[] = {
    {"branchLength", OPTION_DOUBLE},
    {"treelessRoot1", OPTION_STRING},
    {"treelessRoot2", OPTION_STRING},
    {"maf1Copy", OPTION_STRING},
    {"maf2Copy", OPTION_STRING},
    {"discardTwoParents", OPTION_BOOLEAN},
    {"help", OPTION_BOOLEAN},
    {NULL, 0}
};

static char *usageMsg =
    "mafJoin [options] refGenome inMaf1 inMaf2 outMaf\n"
    "\n"
    "Options:\n"
    "  -help\n"
    "  -branchLength=0.1 - branch length to use when generating\n"
    "   a trees for MAF. Defaults to 0.1.\n"
    "  -treelessRoot1=genome - root genome for inMaf1 blocks\n"
    "   that do not have trees (see below).\n"
    "  -treelessRoot2=genome - root genome for inMaf2 blocks\n"
    "   that do not have trees.\n"
    "  -maf1Out=maf1Copy - output maf1 for debugging after adding trees.\n"
    "  -maf2Out=maf2Copy - output maf2 for debugging after adding trees.\n"
    "  -discardTwoParents - if a block with two parents is detected when\n"
    "    joining dups, then discard an arbitrary one. Otherwise it's an error.\n"
    "\n"
    "If MAF blocks (mafAli) don't have a tree associated with them, one\n"
    "will be created.  The root genome for the tree is chosen based on\n"
    "the genome specified by the -treelessRoot1 or -treelessRoot2 options.\n"
    "One sequence from that genome becomes the root and the remainder\n"
    "become it's direct children.\n";

/* usage msg and exit */
static void usage(char *msg) {
    errAbort("Error: %s\n%s", msg, usageMsg);
}

/* join two mafs */
static void mafJoin(char *refGenomeName, char *inMaf1, char *inMaf2, char *outMaf, double defaultBranchLength,
                    char *treelessRoot1Name, char *treelessRoot2Name,
                    char *maf1Copy, char *maf2Copy, bool discardTwoParents) {
    struct Genomes *genomes = genomesNew();
    struct Genome *refGenome = genomesObtainGenome(genomes, refGenomeName);
    struct Genome *treelessRoot1Genome = (treelessRoot1Name != NULL) ? genomesObtainGenome(genomes, treelessRoot1Name) : NULL;
    struct Genome *treelessRoot2Genome = (treelessRoot2Name != NULL) ? genomesObtainGenome(genomes, treelessRoot2Name) : NULL;

    struct malnSet *malnSet1 = malnSet_constructFromMaf(genomes, inMaf1, defaultBranchLength, treelessRoot1Genome);
    malnMultiParents_check(malnSet1, discardTwoParents);
    malnJoin_joinSetDups(malnSet1);
    if (maf1Copy != NULL) {
        malnSet_writeMaf(malnSet1, maf1Copy);
    }

    struct malnSet *malnSet2 = malnSet_constructFromMaf(genomes, inMaf2, defaultBranchLength, treelessRoot2Genome);
    malnMultiParents_check(malnSet2, discardTwoParents);
    malnJoin_joinSetDups(malnSet2);
    if (maf2Copy != NULL) {
        malnSet_writeMaf(malnSet2, maf2Copy);
    }

    // join and then merge overlapping blocks that were created
    struct malnSet *malnSetJoined = malnJoinSets(refGenome, malnSet1, malnSet2);
    malnJoin_joinSetDups(malnSetJoined);
    malnSet_writeMaf(malnSetJoined, outMaf);

    if (memLeakDebugCleanup) {
        malnSet_destruct(malnSet1);
        malnSet_destruct(malnSet2);
        malnSet_destruct(malnSetJoined);
        genomesFree(genomes);
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
    mafJoin(argv[1], argv[2], argv[3], argv[4], optionDouble("branchLength", 0.1), 
            optionVal("treelessRoot1", NULL), optionVal("treelessRoot2", NULL),
            optionVal("maf1Copy", NULL), optionVal("maf2Copy", NULL),
            optionExists("discardTwoParents"));
    return 0;
}
