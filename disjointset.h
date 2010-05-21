#ifndef _DISJOINTSET_H_
#define _DISJOINTSET_H_

#include <stdint.h>
#include <stdlib.h>
#include "commonC.h"

struct djs 
{
  int32_t *p;
  int32_t *rank;
};

struct djs * djs_new (int32_t size);
void djs_free (struct djs *unionFind);

void djs_makeset (struct djs *unionFind, int32_t x);
void djs_link (struct djs *unionFind, int32_t x, int32_t y);
int32_t djs_findset (struct djs *unionFind, int32_t x);
void djs_union (struct djs *unionFind, int32_t x, int32_t y);

#endif /* _DISJOINTSET_H_ */
