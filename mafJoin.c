#include "common.h"
#include "options.h"
#include "genome.h"
#include "mafInvert.h"
#include "mafInvertJoin.h"
#include "mafRevert.h"

/*
 * Notes: 
 *  - this would probably be simpler (and for sure more efficient) if
 *    it was done with block based rather than column based merges.
 */


/* command line option specifications */
static struct optionSpec optionSpecs[] = {
    {NULL, 0}
};

static char *usageMsg =
    "mafJoin refGenome inMaf1 inMaf2 outMaf\n";

/* usage msg and exit */
static void usage(char *msg) {
    errAbort("Error: %s\n%s", msg, usageMsg);
}

/* join two mafs */
static void mafJoin(char *refGenomeName, char *inMaf1, char *inMaf2, char *outMaf) {
    struct Genomes *genomes = genomesNew();
    struct Genome *ref = genomesObtainGenome(genomes, refGenomeName);
    struct MafInvert *maf1Inv = mafInvertBuild(genomes, ref, inMaf1);
    struct MafInvert *maf2Inv = mafInvertBuild(genomes, ref, inMaf2);
    struct MafInvert *joinMafInv = mafInvertJoin(genomes, maf1Inv, maf2Inv);
    mafInvertToMaf(joinMafInv, outMaf);
}

/* Process command line. */
int main(int argc, char *argv[]) {
    optionInit(&argc, argv, optionSpecs);
    if (argc != 5)  {
        usage("Error: wrong number of arguments");
    }
    mafJoin(argv[1], argv[2], argv[3], argv[4]);
    return 0;
}
