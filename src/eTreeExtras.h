#ifndef ETREEEXTRAS_H_
#define ETREEEXTRAS_H_

#include <assert.h>

#include "bioioC.h"
#include "commonC.h"
#include "sonLib.h"

typedef struct leafPtrArray {
        void **ptrArray;
        int32_t index;
} LeafPtrArray;

void eTreeX_countLeaves(ETree *node, void *data);
void eTreeX_getLeafArray(ETree *node, void *data);
void eTreeX_postOrderTraversal(ETree *node, void (*func)(ETree *, void *), void *data);
ETree *eTreeX_getTreeFromFile(char *treeFile);

LeafPtrArray *eTreeX_constructLeafPtrArray(int32_t leafCount);
void eTreeX_destructLeafPtrArray(LeafPtrArray *leafArray);

#endif /* ETREEEXTRAS_H_ */
