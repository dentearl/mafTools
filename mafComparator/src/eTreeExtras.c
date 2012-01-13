/* 
 * Copyright (C) 2009-2012 by 
 * Dent Earl (dearl@soe.ucsc.edu, dentearl@gmail.com)
 * Benedict Paten (benedict@soe.ucsc.edu, benedictpaten@gmail.com)
 * Mark Diekhans (markd@soe.ucsc.edu)
 * ... and other members of the Reconstruction Team of David Haussler's 
 * lab (BME Dept. UCSC).
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE. 
*/


#include "eTreeExtras.h"

void eTreeX_countLeaves(stTree *node, void *data) {
    if (stTree_getChildNumber(node) == 0) {
        *((int *) data) += 1;
    }
    return;
}

void eTreeX_getLeafArray(stTree *node, void *data) {
    LeafPtrArray *p = (LeafPtrArray *) data;

    if (stTree_getChildNumber(node) == 0) {
        p->ptrArray[p->index] = (void *) node;
        p->index++;
    }
    return;
}

void eTreeX_postOrderTraversal(stTree *node, void (*func)(stTree *, void *), void *data) {
    int32_t i = 0;
    stTree *child = NULL;
    for (i=0; i<stTree_getChildNumber(node); i++) {
        child = stTree_getChild(node, i);
        eTreeX_postOrderTraversal(child, func, data);
    }
    func(node, data);
    return;
}

stTree *eTreeX_getTreeFromFile(char *treeFile) {
    int bytesRead;
    int nBytes = 100;
    char *cA = st_malloc(nBytes + 1);

    FILE *fileHandle = fopen(treeFile, "r");
    if (fileHandle == NULL) {
        fprintf(stderr, "Error: Can't open file '%s' for reading", treeFile);
        exit(1);
    }
    bytesRead = benLine(&cA, &nBytes, fileHandle);
    assert(bytesRead == -1); // Only read a single line from file
    fclose(fileHandle);

    stTree *tree = NULL;
    tree = stTree_parseNewickString(cA);

    free(cA);

    return tree;
}

LeafPtrArray *eTreeX_constructLeafPtrArray(int32_t leafCount) {
    LeafPtrArray *leafArray = NULL;
    leafArray = st_malloc(sizeof(LeafPtrArray));
    leafArray->ptrArray = st_malloc(sizeof(void *) * leafCount);
    leafArray->index = 0;

    return leafArray;
}

void eTreeX_destructLeafPtrArray(LeafPtrArray *leafArray) {
    free(leafArray->ptrArray);
    free(leafArray);

    return;
}
