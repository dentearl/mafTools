#ifndef stMafTree_h
#define stMafTree_h
#include "sonLibTypes.h"
#include "mafJoinTypes.h"
struct mafAli;

/* construct a MafTree object from an eTree and corresponding MafComp block */
stMafTree *stMafTree_constructFromMaf(ETree *eTree, struct mafAli *ali);

/* destroy a MafTree object */
void stMafTree_destruct(stMafTree *mTree);

/* validate tree order and nodes matches mafAli */
void stMafTree_validateWithMaf(stMafTree *mTree, struct mafAli *ali);

/* get the associated ETree */
ETree *stMafTree_getETree(stMafTree *mTree);

/* Find MAF order number for an organism/sequence, return -1 if not found */
int stMafTree_findMafOrder(stMafTree *mTree, const char *orgSeq);

/* Get MAF order number for an organism/sequence, error if not found */
int stMafTree_getMafOrder(stMafTree *mTree, const char *orgSeq);

/* Get an eTree node by MAF order number */
ETree *stMafTree_getByMafOrder(stMafTree *mTree, int mafOrder);

/* join two trees at the nodes identified by their maf order */
stMafTree *stMafTree_join(stMafTree *mTree1, int mafOrder1, stMafTree *mTree2, int mafOrder2);
#endif
