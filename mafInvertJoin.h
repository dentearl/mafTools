/* join two mafInvert objects */
#ifndef mafInvertJoin_h
#define mafInvertJoin_h
struct Genomes;
struct MafInvert;

/* Join two inverted mafs based on the reference sequence */
struct MafInvert *mafInvertJoin(struct Genomes *genomes, struct MafInvert *mi1, struct MafInvert *mi2);

#endif
