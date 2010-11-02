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

static bool debug = false;  // FIXME: tmp

/* set to true to cleanup all memory for memory leak checks */
static bool memLeakDebugCleanup = true;

/* command line option specifications */
static struct optionSpec optionSpecs[] = {
    {"branchLength", OPTION_DOUBLE},
    {"treelessRoot1", OPTION_STRING},
    {"treelessRoot2", OPTION_STRING},
    {"maxInputBlkWidth", OPTION_INT},
    {"maxBlkWidth", OPTION_INT},
    {"multiParentDropped", OPTION_STRING},
    {"inMaf1Dump", OPTION_STRING},
    {"inMaf2Dump", OPTION_STRING},
    {"outMafDump", OPTION_STRING},
    {"help", OPTION_BOOLEAN},
    {NULL, 0}
};

static char *usageMsg =
    "mafJoin [options] guideGenome inMaf1 inMaf2 outMaf\n"
    "\n"
    "Options:\n"
    "  -help\n"
    "  -branchLength=0.1 - branch length to use when generating\n"
    "   a trees for MAF. Defaults to 0.1.\n"
    "  -treelessRoot1=genome - root genome for inMaf1 blocks\n"
    "   that do not have trees (see below).\n"
    "  -treelessRoot2=genome - root genome for inMaf2 blocks\n"
    "   that do not have trees.\n"
    "  -maxInputBlkWidth=n - split input block to this size to limit CPU and\n"
    "   memory consumption during merge. NOT RECOMMENDED\n"
    "  -maxBlkWidth=n - set a cutoff on the maximum width of an alignment\n"
    "   block when doing the final join blocks in the output MAF.  This is\n"
    "   used to limit size of blocks, near the Evolover root were long,\n"
    "   contiguous overlapping regions can cause exponential growth in run\n"
    "   time and memory requirements\n"
    "  -multiParentDropped=file - write regions were one copy was drop due to\n"
    "   having multiple parents.\n"
    "  -inMaf1Dump=file - dump info about input inMaf1 after merging dups\n"
    "  -inMaf2Dump=file - dump info about input inMaf2 after merging dups\n"
    "  -outMafDump=file - dump info about final output maf\n"
    "\n"
    "If MAF blocks (mafAli) don't have a tree associated with them, one\n"
    "will be created.  The root genome for the tree is chosen based on\n"
    "the genome specified by the -treelessRoot1 or -treelessRoot2 options.\n"
    "One sequence from that genome becomes the root and the remainder\n"
    "become it's direct children.  If  -treelessRoot option is specified, it\n"
    "also triggers merging of the duplication blocks that Evolver outputs.\n";


/* usage msg and exit */
static void usage(char *msg) {
    errAbort("Error: %s\n%s", msg, usageMsg);
}

/* load a MAF and do internal joining.  */
static struct malnSet *loadMaf(struct Genomes *genomes, char *inMaf, int maxInputBlkWidth, double defaultBranchLength,
                               char *treelessRootName, char *setName, FILE *dropLogFh) {
    struct Genome *treelessRootGenome = (treelessRootName != NULL) ? genomesObtainGenome(genomes, treelessRootName) : NULL;
    struct malnSet *malnSet = malnSet_constructFromMaf(genomes, inMaf, maxInputBlkWidth, defaultBranchLength, treelessRootGenome);
    if (debug) {
        malnSet_dump(malnSet, stderr, "%s input", setName);
    }
    if (treelessRootGenome != NULL) {
        malnMultiParents_resolve(malnSet, dropLogFh);
        malnJoinWithinSet_joinDups(malnSet);
        if (debug) {
            malnSet_dump(malnSet, stderr, "%s: joined-dups", setName);
        }
    }
    malnMultiParents_check(malnSet); // FIXME: no longer needed due to above resolve step
    return malnSet;
}

/* join two mafs */
static void mafJoin(char *guideGenomeName, char *inMaf1, char *inMaf2, char *outMaf, double defaultBranchLength,
                    char *treelessRoot1Name, char *treelessRoot2Name, int maxInputBlkWidth, int maxBlkWidth,
                    char *multiParentDroppedFile, char *inMaf1Dump, char *inMaf2Dump, char *outMafDump) {
    struct Genomes *genomes = genomesNew();
    struct Genome *guideGenome = genomesObtainGenome(genomes, guideGenomeName);
    FILE *dropLogFh = (multiParentDroppedFile != NULL) ? malnMultiParents_openResolveDropLog(multiParentDroppedFile) : NULL;
    struct malnSet *malnSet1 = loadMaf(genomes, inMaf1, maxInputBlkWidth, defaultBranchLength, treelessRoot1Name, "set1", dropLogFh);
    if (inMaf1Dump != NULL) {
        malnSet_dumpFile(malnSet1, inMaf1Dump, "set1");
    }
    struct malnSet *malnSet2 = loadMaf(genomes, inMaf2, maxInputBlkWidth, defaultBranchLength, treelessRoot2Name, "set2", dropLogFh);
    if (inMaf2Dump != NULL) {
        malnSet_dumpFile(malnSet2, inMaf2Dump, "set2");
    }
    carefulClose(&dropLogFh);

    // join and then merge overlapping blocks that were created
    struct malnSet *malnSetJoined = malnJoinSets(guideGenome, malnSet1, malnSet2);
    if (debug) {
        malnSet_dump(malnSetJoined, stderr, "out: joined");
    }
    malnJoinWithinSet_joinOverlapAdjacent(malnSetJoined, maxBlkWidth);
    if (debug) {
        malnSet_dump(malnSetJoined, stderr, "out: joined-adjacent");
    }
    malnMergeComps_merge(malnSetJoined);
    if (debug) {
        malnSet_dump(malnSetJoined, stderr, "out: merged");
    }
    malnMultiParents_check(malnSetJoined);
    if (outMafDump != NULL) {
        malnSet_dumpFile(malnSetJoined, outMafDump, "out");
    }
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
            optionInt("maxInputBlkWidth", INT_MAX), optionInt("maxBlkWidth", INT_MAX),
            optionVal("multiParentDropped", NULL),
            optionVal("inMaf1Dump", NULL), optionVal("inMaf2Dump", NULL), optionVal("outMafDump", NULL));
    return 0;
}
