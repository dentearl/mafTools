/* merging of adjacent components within a block */
#ifndef malnMergeComps_h
#define malnMergeComps_h
struct malnSet *malnSet;

/* Merge adjacent components within the blocks of a set. */
void malnMergeComps_merge(struct malnSet *malnSet);

#endif
