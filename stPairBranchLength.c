#include "stPairBranchLength.h"
#include "sonLibETree.h"
#include "sonLibString.h"
#include <stdbool.h>

/* FIXME: end up not using this and just having a fixed branch length as
 * input, but keeping it around just in case this is wrong.
 */
#if 0
// stuff in mafJoin.c
    {"initBranchLenTree", OPTION_STRING},
    "  -initBranchLenTree=tree.nh - newick tree use to obtain branch lengths when generating\n"
    "   a tree from a pair-wise MAF.\n"

    ETree *initBranchLenTree = NULL;
    if (initBranchLenTreeFile != NULL) {
        initBranchLenTree = eTree_parseNewickString(initBranchLenTreeFile);
    }

#endif

/* keep track of where we are in the search */
struct tracker {
    char *org1;
    char *org2;
    double distance;
};

/* DFS search measuring second org, returns true when time to stop  */
static bool dfsMeasureSecond(ETree *root, struct tracker *tracker, char *org, double distance) {
    if (stString_eq(eTree_getLabel(root), org)) {
        tracker->distance = distance;
        return true;
    } else {
        for (int i = 0; i < eTree_getChildNumber(root); i++) {
            if (dfsMeasureSecond(eTree_getChild(root, i), tracker, org, distance+eTree_getBranchLength(root))) {
                return true;
            }
        }
        return false;
    }
}

/* DFS search for first org, returns true when time to stop  */
static bool dfsFindFirst(ETree *root, struct tracker *tracker) {
    if (stString_eq(eTree_getLabel(root), tracker->org1)) {
        return dfsMeasureSecond(root, tracker, tracker->org2, 0.0);
    } else if (stString_eq(eTree_getLabel(root), tracker->org2)) {
        return dfsMeasureSecond(root, tracker, tracker->org1, 0.0);
    } else {
        for (int i = 0; i < eTree_getChildNumber(root); i++) {
            if (dfsFindFirst(eTree_getChild(root, i), tracker)) {
                return true;
            }
        }
        return false;
    }
}


/* get distance between two species in the tree as an estimate of length between
 * species. */
double stPairBranchLength(ETree *initBranchLenTree, char *org1, char *org2) {
    struct tracker tracker = {org1, org2, 0.0};

    return tracker.distance;
}

