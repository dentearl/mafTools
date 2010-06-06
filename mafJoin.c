#include "common.h"
#include "options.h"
#include "genome.h"
#include "mafInvert.h"
#include "mafInvertJoin.h"
#include "mafInvertWrite.h"

/*
 * Notes: 
 *  - this would probably be simpler (and for sure more efficient) if
 *    it was done with block based rather than column based merges.
 */


/* command line option specifications */
static struct optionSpec optionSpecs[] = {
    {"maf1Copy", OPTION_STRING},
    {"maf2Copy", OPTION_STRING},
    {NULL, 0}
};

static char *usageMsg =
    "mafJoin [options] refGenome inMaf1 inMaf2 outMaf\n"
    "\n"
    "-maf1Out=maf1Copy - output maf1 for debugging\n"
    "-maf2Out=maf2Copy - output maf1 for debugging\n";

/* usage msg and exit */
static void usage(char *msg) {
    errAbort("Error: %s\n%s", msg, usageMsg);
}

/* join two mafs */
static void mafJoin(char *refGenomeName, char *inMaf1, char *inMaf2, char *outMaf, char *maf1Copy, char *maf2Copy) {
    struct Genomes *genomes = genomesNew();
    struct Genome *ref = genomesObtainGenome(genomes, refGenomeName);
    struct MafInvert *maf1Inv = mafInvertBuild(genomes, ref, inMaf1);
    if (maf1Copy != NULL) {
        mafInvertToMaf(maf1Inv, maf1Copy);
    }
    struct MafInvert *maf2Inv = mafInvertBuild(genomes, ref, inMaf2);
    if (maf2Copy != NULL) {
        mafInvertToMaf(maf2Inv, maf2Copy);
    }
    struct MafInvert *joinMafInv = mafInvertJoin(genomes, maf1Inv, maf2Inv);
    mafInvertToMaf(joinMafInv, outMaf);
}

/* Process command line. */
int main(int argc, char *argv[]) {
    optionInit(&argc, argv, optionSpecs);
    if (argc != 5)  {
        usage("Error: wrong number of arguments");
    }
    mafJoin(argv[1], argv[2], argv[3], argv[4], optionVal("maf1Copy", NULL), optionVal("maf1Copy", NULL));
    return 0;
}
