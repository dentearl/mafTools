#include "common.h"
#include "options.h"
#include "genome.h"
#include "malnSet.h"
#include "malnJoinDups.h"
#include "malnJoinSets.h"
#include "sonLibETree.h"

/*
 * Notes: 
 *  - this would probably be simpler (and for sure more efficient) if
 *    it was done with block based rather than column based merges.
 */

/* set to true to cleanup all memory for memory leak checks */
static bool memLeakDebugCleanup = false;

/* command line option specifications */
static struct optionSpec optionSpecs[] = {
    {"branchLength", OPTION_DOUBLE},
    {"maf1Copy", OPTION_STRING},
    {"maf2Copy", OPTION_STRING},
    {NULL, 0}
};

static char *usageMsg =
    "mafJoin [options] refGenome inMaf1 inMaf2 outMaf\n"
    "\n"
    "Options:\n"
    "  -branchLength=0.1 - branch length to use when generating\n"
    "   a tree from a pair-wise MAF. Defaults to 0.1.\n"
    "  -maf1Out=maf1Copy - output maf1 for debugging\n"
    "  -maf2Out=maf2Copy - output maf2 for debugging\n";

/* usage msg and exit */
static void usage(char *msg) {
    errAbort("Error: %s\n%s", msg, usageMsg);
}

/* join two mafs */
static void mafJoin(char *refGenomeName, char *inMaf1, char *inMaf2, char *outMaf, double defaultBranchLength, char *maf1Copy, char *maf2Copy) {
    struct Genomes *genomes = genomesNew();
    struct Genome *refGenome = genomesObtainGenome(genomes, refGenomeName);
    struct malnSet *malnSet1 = malnSet_constructFromMaf(genomes, inMaf1, defaultBranchLength);
    malnJoin_joinSetDups(malnSet1);
    if (maf1Copy != NULL) {
        malnSet_writeMaf(malnSet1, maf1Copy);
    }
    struct malnSet *malnSet2 = malnSet_constructFromMaf(genomes, inMaf2, defaultBranchLength);
    malnJoin_joinSetDups(malnSet2);
    if (maf2Copy != NULL) {
        malnSet_writeMaf(malnSet2, maf2Copy);
    }
    struct malnSet *malnSetJoined = malnJoinSets(refGenome, malnSet1, malnSet2);
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
    if (argc != 5)  {
        usage("Error: wrong number of arguments");
    }
    mafJoin(argv[1], argv[2], argv[3], argv[4], optionDouble("branchLength", 0.1), optionVal("maf1Copy", NULL), optionVal("maf2Copy", NULL));
    return 0;
}
