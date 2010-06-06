#ifndef stMafTreeTbl_h
#define stMafTreeTbl_h
#include "mafJoinTypes.h"
struct mafAli;

/* constructor */
stMafTreeTbl *stMafTreeTbl_construct();

/* destructor */
void stMafTreeTbl_destruct(stMafTreeTbl *tbl);

/* obtain a stMafTree objects corresponding to the tree in the mafAli, reusing a
 * existing tree if available. */
stMafTree *stMafTreeTbl_obtainForMafAli(stMafTreeTbl *tbl, struct mafAli *ali);

/* obtain the join of the specified trees, reusing an existing tree if available */
stMafTree *stMafTreeTbl_obtainJoin(stMafTreeTbl *tbl, stMafTree *mTree1, int mafOrder1, stMafTree *mTree2, int mafOrder2);

#endif
