/* genome information */
#include "common.h"
#include "genome.h"
#include "hash.h"
#include "jkmaf.h"

/* make db.seq name */
char *seqMkName(struct Seq *seq) {
    int bufSize = strlen(seq->genome->name) + strlen(seq->name) + 2;
    char *buf = needMem(bufSize);
    safef(buf, bufSize, "%s.%s", seq->genome->name, seq->name);
    return buf;
}

/* constructor  */
static struct Seq *seqNew(struct Genome *genome, char *name, int size) {
    struct Seq *seq;
    AllocVar(seq);
    seq->genome = genome;
    seq->name = cloneString(name);
    seq->size = size;
    return seq;
};

/* constructor */
static struct Genome *genomeNew(char *name) {
    struct Genome *genome;
    AllocVar(genome);
    genome->name = cloneString(name);
    genome->seqMap = hashNew(8);
    return genome;
}

/* obtain a new Seq object, creating if it doesn't exist */
struct Seq *genomeObtainSeq(struct Genome *genome, char *name, int size) {
    struct Seq *seq = hashFindVal(genome->seqMap, name);
    if (seq == NULL) {
        seq = seqNew(genome, name, size);
        hashAdd(genome->seqMap, name, seq);
        slAddHead(&genome->seqs, seq);
    }
    return seq;
}

/* obtain get a Seq object, error doesn't exist */
struct Seq *genomeGetSeq(struct Genome *genome, char *name) {
    struct Seq *seq = hashFindVal(genome->seqMap, name);
    if (seq == NULL) {
        errAbort("can;t find seq \"%s\" in genome \"%s\"", name, genome->name);
    }
    return seq;
}

/* constructor */
struct Genomes *genomesNew() {
    struct Genomes *genomes;
    AllocVar(genomes);
    genomes->genomeMap = hashNew(8);
    return genomes;
}

/* obtain a genome object, constructing a new one if it doesn't exist. */
struct Genome *genomesObtainGenome(struct Genomes *genomes, char *name) {
    struct Genome *genome = hashFindVal(genomes->genomeMap, name);
    if (genome == NULL) {
        genome = genomeNew(name);
        hashAdd(genomes->genomeMap, name, genome);
    }
    return genome;
}

/* get a genome object, error if it doesn't exist */
struct Genome *genomesGetGenome(struct Genomes *genomes, char *name) {
    struct Genome *genome = hashFindVal(genomes->genomeMap, name);
    if (genome == NULL) {
        errAbort("can't find genome \"%s\"", name);
    }
    return genome;
}


/* Obtain a new Seq object, creating the genome and seq objects it they
 * don't exist.  If the genomeName is NULL, the seqName is used as the
 * genome name.  This behavior corresponds to having a maf with only
 * one name. */
struct Seq *genomesObtainSeq(struct Genomes *genomes, char *genomeName, char *seqName, int size) {
    if (genomeName == NULL) {
        genomeName = seqName;
    }
    struct Genome *genome = genomesObtainGenome(genomes, genomeName);
    return genomeObtainSeq(genome, seqName, size);
}
