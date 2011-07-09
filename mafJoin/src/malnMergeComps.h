/* merging of adjacent components within a block */
#ifndef malnMergeComps_h
#define malnMergeComps_h
struct malnSet *malnSet;

/* Merge overlapping and adjacent components within blocks. */
void malnMergeComps_merge(struct malnSet *malnSet);

#endif
