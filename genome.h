/* genome information */
#ifndef genome_h
#define genome_h
#include <stdbool.h>
#include <ctype.h>

/* a [0-n) range */
struct stRange {
    int start;
    int end;
};

/* range representing nothing {-1, -1} */
extern const struct stRange stNullRange;

/* is a stRange null */
static inline bool stRangeIsNull(struct stRange range) {
    return (range.start == stNullRange.start) && (range.end == stNullRange.end);
}

/* are two stRanges equal */
static inline bool stRangeEq(struct stRange range0, struct stRange range1) {
    return (range0.start == range1.start) && (range0.end == range1.end);
}

/* a sequence in the genome */
struct Seq {
    struct Seq *next;
    struct Genome *genome;
    char *name;
    char *orgSeqName;  // db.seq
    int size;  // -1 of not known
};

/* an organism's genome */
struct Genome {
    struct Genome *next;
    char *name;
    struct hash *seqMap;  /* map and list of sequences */
    struct Seq *seqs;
};

/* does a character represent a base */
static inline bool isBase(char base) {
    // n.b. isalpha doesn't return 0/1, might be out of bool range
    return isalpha(base) != 0;
}

/* are two bases the same, ignoring case */
static inline bool baseEq(char base1, char base2) {
    return toupper(base1) == toupper(base2);
}

/* Table of all genomes. Genomes are normally added in a lazy manner while
 * scanning a MAF. */
struct Genomes {
    struct hash *genomeMap;
};

/* compare two sequences for deterministic sorting */
int seqCmp(struct Seq *seq1, struct Seq *seq2);

/* obtain a new Seq object, creating if it doesn't exist */
struct Seq *genomeObtainSeq(struct Genome *genome, char *name, int size);

/* obtain get a Seq object, error doesn't exist */
struct Seq *genomeGetSeq(struct Genome *genome, char *name);

/* compare two genomes for deterministic sorting */
int genomeCmp(struct Genome *genome1, struct Genome *genome2);

/* constructor */
struct Genomes *genomesNew(void);

/* destructor */
void genomesFree(struct Genomes *genomes);

/* obtain a genome object, constructing a new one if it doesn't exist */
struct Genome *genomesObtainGenome(struct Genomes *genomes, char *name);

/* get a genome object, error if it doesn't exist */
struct Genome *genomesGetGenome(struct Genomes *genomes, char *name);

/* Obtain a new Seq object, creating the genome and seq objects it they don't
 * exist. If size is -1, then it will not be initialized until a request is
 * made with the size. */
struct Seq *genomesObtainSeq(struct Genomes *genomes, char *genomeName, char *seqName, int size);

/* Obtain a new Seq object given organism.seq, creating the genome and seq objects it they
 * don't exist. If size is -1, then it will not be initialized until a request is
 * made with the size. */
struct Seq *genomesObtainSeqForOrgSeqName(struct Genomes *genomes, char *orgSeqName, int size);

#endif
