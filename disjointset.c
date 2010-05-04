#include "disjointset.h"

struct djs * djs_new (int32_t size)
{
  struct djs *unionFind = NULL;
  unionFind = malloc(sizeof(struct djs));

  unionFind->p = malloc(sizeof(int32_t) * size);
  unionFind->rank = malloc(sizeof(int32_t) * size);

  return unionFind;
}

void djs_free (struct djs *unionFind)
{
  free(unionFind->p);
  free(unionFind->rank);
  free(unionFind);
}

void djs_makeset (struct djs *unionFind, int32_t x)
{
  unionFind->p[x] = x;
  unionFind->rank[x] = 0;
}


void djs_link (struct djs *unionFind, int32_t x, int32_t y)
{
  if (unionFind->rank[x] > unionFind->rank[y]) {
    unionFind->p[y] = x;
  } else {
    unionFind->p[x] = y;
    if (unionFind->rank[x] == unionFind->rank[y]) {
      unionFind->rank[y] += 1;
    }
  }
}

int32_t djs_findset (struct djs *unionFind, int32_t x)
{
  if (x != unionFind->p[x]) {
    unionFind->p[x] = djs_findset(unionFind, unionFind->p[x]);
  }
  return unionFind->p[x];
}

void djs_union (struct djs *unionFind, int32_t x, int32_t y)
{
  djs_link(unionFind, djs_findset(unionFind, x), djs_findset(unionFind, y));
}
