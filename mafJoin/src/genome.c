/* genome information */
#include "common.h"
#include "genome.h"
#include "hash.h"
#include "jkmaf.h"
#include "sonLibString.h"

// FIXME: make naming consistent with other modules

/* range representing nothing */
const struct stRange stNullRange = {-1, -1};

/* make db.seq name */
static char *seqMkName(struct Seq *seq) {
    int bufSize = strlen(seq->genome->name) + strlen(seq->name) + 2;
    char *buf = needMem(bufSize);
    safef(buf, bufSize, "%s.%s", seq->genome->name, seq->name);
    return buf;
}

/* constructor  */
static struct Seq *seqNew(struct Genome *genome, const char *name, int size) {
    struct Seq *seq;
    AllocVar(seq);
    seq->genome = genome;
    seq->name = stString_copy(name);
    seq->size = size;
    seq->orgSeqName = seqMkName(seq);
    return seq;
}

/* destructor */
static void seqFree(struct Seq *seq) {
    freeMem(seq->name);
    freeMem(seq->orgSeqName);
    freeMem(seq);
}

/* compare two sequences for deterministic sorting */
int seqCmp(struct Seq *seq1, struct Seq *seq2) {
    if (seq1->genome == seq2->genome) {
        return strcmp(seq1->name, seq2->name);
    } else {
        return genomeCmp(seq1->genome, seq2->genome);
    }
}

/* constructor */
static struct Genome *genomeNew(const char *name) {
    struct Genome *genome;
    AllocVar(genome);
    genome->name = stString_copy(name);
    genome->seqMap = hashNew(8);
    return genome;
}

/* destructor */
static void genomeFree(struct Genome *genome) {
    struct Seq *seq;
    while ((seq = slPopHead(&genome->seqs)) != NULL) {
        seqFree(seq);
    }
    freeMem(genome->name);
    hashFree(&genome->seqMap);
    freeMem(genome);
}

/* obtain a new Seq object, creating if it doesn't exist */
struct Seq *genomeObtainSeq(struct Genome *genome, const char *name, int size) {
    struct Seq *seq = hashFindVal(genome->seqMap, (char*)name);
    if (seq == NULL) {
        seq = seqNew(genome, name, size);
        hashAdd(genome->seqMap, (char*)name, seq);
        slAddHead(&genome->seqs, seq);
    }
    return seq;
}

/* obtain get a Seq object, error doesn't exist */
struct Seq *genomeGetSeq(struct Genome *genome, const char *name) {
    struct Seq *seq = hashFindVal(genome->seqMap, (char*)name);
    if (seq == NULL) {
        errAbort("can;t find seq \"%s\" in genome \"%s\"", name, genome->name);
    }
    return seq;
}

/* compare two genomes for deterministic sorting */
int genomeCmp(struct Genome *genome1, struct Genome *genome2) {
    return strcmp(genome1->name, genome2->name);
}

/* constructor */
struct Genomes *genomesNew(void) {
    struct Genomes *genomes;
    AllocVar(genomes);
    genomes->genomeMap = hashNew(8);
    return genomes;
}

/* destructor */
void genomesFree(struct Genomes *genomes) {
    struct hashCookie cookie = hashFirst(genomes->genomeMap);
    struct hashEl *hel;
    while ((hel = hashNext(&cookie)) != NULL) {
        genomeFree(hel->val);
    }
    hashFree(&genomes->genomeMap);
    freeMem(genomes);
}

/* obtain a genome object, constructing a new one if it doesn't exist. */
struct Genome *genomesObtainGenome(struct Genomes *genomes, const char *name) {
    struct Genome *genome = hashFindVal(genomes->genomeMap, (char*)name);
    if (genome == NULL) {
        genome = genomeNew(name);
        hashAdd(genomes->genomeMap, (char*)name, genome);
    }
    return genome;
}

/* get a genome object, error if it doesn't exist */
struct Genome *genomesGetGenome(struct Genomes *genomes, const char *name) {
    struct Genome *genome = hashFindVal(genomes->genomeMap, (char*)name);
    if (genome == NULL) {
        errAbort("can't find genome \"%s\"", name);
    }
    return genome;
}


/* Obtain a new Seq object, creating the genome and seq objects it they don't
 * exist. If size is -1, then it will not be initialized until a request is
 * made with the size. */
struct Seq *genomesObtainSeq(struct Genomes *genomes, const char *genomeName, const char *seqName, int size) {
    struct Genome *genome = genomesObtainGenome(genomes, genomeName);
    struct Seq *seq = genomeObtainSeq(genome, seqName, size);
    if ((seq->size < 0) && (size >= 0)) {
        seq->size = size;
    }
    return seq;
}

/* Obtain a new Seq object given organism.seq, creating the genome and seq objects it they
 * don't exist. If size is -1, then it will not be initialized until a request is
 * made with the size. */
struct Seq *genomesObtainSeqForOrgSeqName(struct Genomes *genomes, const char *orgSeqName, int size) {
    char nameBuf[strlen(orgSeqName)+1];
    strcpy(nameBuf, orgSeqName);
    char *dot = strchr(nameBuf, '.');
    if (dot == NULL) {
        errAbort("sequence name not in the form org.seq: %s", orgSeqName);
    }
    *dot = '\0';
    char *genomeName = nameBuf;
    char *seqName = dot+1;
    return genomesObtainSeq(genomes, genomeName, seqName, size);
}

