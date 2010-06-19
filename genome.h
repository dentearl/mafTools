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

/* a sequence in the genome */
struct Seq {
    struct Seq *next;
    struct Genome *genome;
    char *name;
    char *orgSeqName;  // db.seq
    int size;
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

/* obtain a new Seq object, creating if it doesn't exist */
struct Seq *genomeObtainSeq(struct Genome *genome, char *name, int size);

/* obtain get a Seq object, error doesn't exist */
struct Seq *genomeGetSeq(struct Genome *genome, char *name);

/* constructor */
struct Genomes *genomesNew(void);

/* obtain a genome object, constructing a new one if it doesn't exist */
struct Genome *genomesObtainGenome(struct Genomes *genomes, char *name);

/* get a genome object, error if it doesn't exist */
struct Genome *genomesGetGenome(struct Genomes *genomes, char *name);

/* Obtain a new Seq object, creating the genome and seq objects it they
 * don't exist. */
struct Seq *genomesObtainSeq(struct Genomes *genomes, char *genomeName, char *seqName, int size);

#endif
