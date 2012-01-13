/* Disjoint-Set Data Structure
 * 
 * Implementation of the disjoint-set forest data structure.
 * Based off pseudocode in "Introduction to Algorithms" by Cormen, et al.
 * 
 * Copyright (C) 2009-2012 by 
 * Bernard Suh (bsuh@soe.ucsc.edu)
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


#ifndef DISJOINTSET_H_
#define DISJOINTSET_H_

#include "commonC.h"

/* Structure for a disjoint-set forest data structure */
struct djs 
{
    uint32_t *p;    /* Parent of node x is p[x] */
    uint32_t *rank; /* Rank of node x is rank[x] */
    uint32_t size;  /* Maximum number of elements in set */
};

/* 
 * Create a new disjoint-set with a max set size equal to "size"
 */
struct djs *djs_new(uint32_t size);

/* 
 * Free the memory associated with a disjoint-set
 */
void djs_free(struct djs *unionFind);

/* 
 * Adds the element "x" to the set
 */
void djs_makeset(struct djs *unionFind, uint32_t x);

/* 
 * Returns the set that element "x" belongs to.
 */
uint32_t djs_findset (struct djs *unionFind, uint32_t x);

/* 
 * Unites the sets that contain elements "x" and "y"
 */
void djs_union(struct djs *unionFind, uint32_t x, uint32_t y);

#endif /* DISJOINTSET_H_ */
