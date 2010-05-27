/*
 * Disjoint-Set Data Structure
 *
 * Implementation of the disjoint-set forest data structure.
 * Based off pseudocode in "Introduction to Algorithms" by Cormen, et al.
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
