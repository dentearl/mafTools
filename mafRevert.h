/* convert a mafInvert object back to a MAF file */
#ifndef mafRevert_h
#define mafRevert_h
struct MafInvert;

/* convert a MafInvert to a MAF */
void mafInvertToMaf(struct MafInvert *mi, char *mafFile);

#endif
