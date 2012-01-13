/* Copyright (C) 2009-2012 by  */
/* Dent Earl (dearl@soe.ucsc.edu, dentearl@gmail.com) */
/* Benedict Paten (benedict@soe.ucsc.edu, benedictpaten@gmail.com) */
/* Mark Diekhans (markd@soe.ucsc.edu) */
/* ... and other members of the Reconstruction Team of David Haussler's  */
/* lab (BME Dept. UCSC). */

/* Permission is hereby granted, free of charge, to any person obtaining a copy */
/* of this software and associated documentation files (the "Software"), to deal */
/* in the Software without restriction, including without limitation the rights */
/* to use, copy, modify, merge, publish, distribute, sublicense, and/or sell */
/* copies of the Software, and to permit persons to whom the Software is */
/* furnished to do so, subject to the following conditions: */

/* The above copyright notice and this permission notice shall be included in */
/* all copies or substantial portions of the Software. */

/* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR */
/* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, */
/* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE */
/* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER */
/* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, */
/* OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN */
/* THE SOFTWARE. */


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

void eTreeX_countLeaves(stTree *node, void *data);
void eTreeX_getLeafArray(stTree *node, void *data);
void eTreeX_postOrderTraversal(stTree *node, void (*func)(stTree *, void *), void *data);
stTree *eTreeX_getTreeFromFile(char *treeFile);

LeafPtrArray *eTreeX_constructLeafPtrArray(int32_t leafCount);
void eTreeX_destructLeafPtrArray(LeafPtrArray *leafArray);

#endif /* ETREEEXTRAS_H_ */
