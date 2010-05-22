#include "eval_ComparatorAPI.h"
#include "cstring.h"

int aSolo_cmpFunction(const void *a, const void *b);
void trioDecoder_destruct(TrioDecoder *decoder);
int32_t calcTrioState(TrioDecoder *decoder, int32_t spA, int32_t spB, int32_t spC);

void aPair_fillOut(APair *aPair, char *seq1, char *seq2, int32_t pos1, int32_t pos2) {
	int32_t i = strcmp(seq1, seq2);
	if(i > 0 || (i == 0 && pos1 > pos2)) { //we construct the sequence so that seq1 < seq2 || seq1 == seq2 and pos1 < pos2
		return aPair_fillOut(aPair, seq2, seq1, pos2, pos1);
	}
	aPair->seq1 = seq1;
	aPair->seq2 = seq2;
	aPair->pos1 = pos1;
	aPair->pos2 = pos2;
}

void aTrio_fillOut(ATrio *aTrio, char *seq1, char *seq2, char *seq3, int32_t pos1, int32_t pos2, int32_t pos3, int32_t top) {
	ASolo structs[3];

	structs[0].name = seq1;
	structs[0].pos = pos1;
	structs[1].name = seq2;
	structs[1].pos = pos2;
	structs[2].name = seq3;
	structs[2].pos = pos3;

	qsort(structs, 3, sizeof(ASolo), aSolo_cmpFunction);
	
	aTrio->seq1 = structs[0].name;
	aTrio->pos1 = structs[0].pos;
	aTrio->seq2 = structs[1].name;
	aTrio->pos2 = structs[1].pos;
	aTrio->seq3 = structs[2].name;
	aTrio->pos3 = structs[2].pos;

	aTrio->top = top;
}

APair *aPair_construct(const char *seq1, const char *seq2, int32_t pos1, int32_t pos2) {
	APair *aPair = mallocLocal(sizeof(APair));
	aPair_fillOut(aPair, (char *)seq1, (char *)seq2, pos1, pos2);
	aPair->seq1 = stringCopy(aPair->seq1);
	aPair->seq2 = stringCopy(aPair->seq2);
	return aPair;
}

ATrio *aTrio_construct(const char *seq1, const char *seq2, const char *seq3, int32_t pos1, int32_t pos2, int32_t pos3, int32_t top) {
	ATrio *aTrio = mallocLocal(sizeof(ATrio));
	aTrio_fillOut(aTrio, (char *)seq1, (char *)seq2, (char *)seq3, pos1, pos2, pos3, top);
	aTrio->seq1 = stringCopy(aTrio->seq1);
	aTrio->seq2 = stringCopy(aTrio->seq2);
	aTrio->seq3 = stringCopy(aTrio->seq3);
	return aTrio;
}

APair *aPair_copyConstruct(APair *pair) {
	return aPair_construct(pair->seq1, pair->seq2, pair->pos1, pair->pos2);
}

ATrio *aTrio_copyConstruct(ATrio *trio) {
	return aTrio_construct(trio->seq1, trio->seq2, trio->seq3, trio->pos1, trio->pos2, trio->pos3, trio->top);
}

void aPair_destruct(APair *pair, void *extraArgument) {
	assert(extraArgument == NULL); //because destruction takes place in avl_tree
	free(pair->seq1);
	free(pair->seq2);
	free(pair);
}

void aTrio_destruct(ATrio *trio, void *extraArgument) {
	assert(extraArgument == NULL); //because destruction takes place in avl_tree
	free(trio->seq1);
	free(trio->seq2);
	free(trio->seq3);
	free(trio);
}

int aSolo_cmpFunction(const void *a, const void *b) {
	ASolo *ia = (ASolo *) a;
	ASolo *ib = (ASolo *) b;

	int cmp = strcmp(ia->name, ib->name);

	if (cmp == 0) {
		return ia->pos - ib->pos;
	}
	return cmp;
}

int32_t aPair_cmpFunction_seqsOnly(APair *aPair1, APair *aPair2, void *a) {
	/*
	 * Compares only the sequence name arguments of the a pairs.
	 */
	assert(a == NULL);
	int32_t i = strcmp(aPair1->seq1, aPair2->seq1);
	if(i == 0) {
		i = strcmp(aPair1->seq2, aPair2->seq2);
	}
	return i;
}

int32_t aTrio_cmpFunction_seqsOnly(ATrio *aTrio1, ATrio *aTrio2, void *a) {
	/*
	 * Compares only the sequence name arguments of the a trios.
	 */
	assert(a == NULL);
	int32_t i = strcmp(aTrio1->seq1, aTrio2->seq1);
	if(i == 0) {
		i = strcmp(aTrio1->seq2, aTrio2->seq2);
		if (i == 0) {
			i = strcmp(aTrio1->seq3, aTrio2->seq3);
		}
	}
	return i;
}

int32_t aPair_cmpFunction(APair *aPair1, APair *aPair2, void *a) {
	/*
	 * Compares first the sequence then the position arguments of the a pairs.
	 */
	int32_t i = aPair_cmpFunction_seqsOnly(aPair1, aPair2, a);
	if(i == 0) {
		i = aPair1->pos1 - aPair2->pos1;
		if(i == 0) {
			i = aPair1->pos2 - aPair2->pos2;
		}
	}
	return i;
}

int32_t aTrio_cmpFunction(ATrio *aTrio1, ATrio *aTrio2, void *a) {
	/*
	 * Compares first the sequence then the position arguments of the a trios.
	 */
	int32_t i = aTrio_cmpFunction_seqsOnly(aTrio1, aTrio2, a);
	if(i == 0) {
		i = aTrio1->pos1 - aTrio2->pos1;
		if(i == 0) {
			i = aTrio1->pos2 - aTrio2->pos2;
			if (i == 0) {
				i = aTrio1->pos3 - aTrio2->pos3;
			}
		}
	}
	return i;
}

void getPairsP(void (*passPairFn)(APair *pair, void *extraArgument1, void *extraArgument2, void *extraArgument3),
		void *extraArgument1, void *extraArgument2, void *extraArgument3, int *bytesRead, int *nBytes, char **cA, FILE *fileHandle) {
	int32_t length, start, seqLength, i, j, k, pos1, pos2, inc1, inc2;
	struct List *ranges = constructEmptyList(0, free);
	char *seqName;
	char *sequence;
	char strand;
	static APair aPair;

	//process block function iterates through successive lines while we have not reached
	//the end of the file and the newline starts with 's'
	//for each line grep the sequence, start position and the length.
	//all the lengths must be equal.
	*bytesRead = benLine (cA, nBytes, fileHandle);

	length = INT32_MAX;
	while(*bytesRead > 0 && (*cA)[0] == 's') {
		seqName = mallocLocal(sizeof(char) * (1+ (*bytesRead)));
		sequence = mallocLocal(sizeof(char) * (1+ (*bytesRead)));
		//uglyf("Got the line :##%s#%i %i\n", *cA, (int)strlen(*cA), *bytesRead);
		//uglyf("I read %i \n", sscanf(*cA, "s %s %i %i %c %i %s", seqName, &start, &i /*ignore the length field*/, &strand, &seqLength, sequence));
		//uglyf("%s,  %i %i %c %i %s\n", seqName, start, i /*ignore the length field*/, strand, seqLength, sequence);
		j = sscanf(*cA, "s %s %i %i %c %i %s", seqName, &start, &i /*ignore the length field*/, &strand, &seqLength, sequence);
		assert(j == 6 || (j == 5 && seqLength == 0));
		if(j == 5) {
			free(sequence);
			sequence = stringPrint("");
		}
		length = strlen(sequence);
		assert(strand == '+' || strand == '-');
		listAppend(ranges, seqName);
		if(strand == '+') {
			listAppend(ranges, constructInt(start));
			listAppend(ranges, constructInt(1));
		}
		else {
			listAppend(ranges, constructInt(seqLength-1-start));
			listAppend(ranges, constructInt(-1));
		}
		listAppend(ranges, sequence);
		*bytesRead = benLine (cA, nBytes, fileHandle);
	}

	//Now call the pair function for every pair of aligned bases.
	for(i=0; i<ranges->length; i+=4) {
		char *seq1 = ranges->list[i];
		inc1 = *((int32_t *)ranges->list[i+2]);
		sequence = ranges->list[i+3];

		for(j=i+4; j<ranges->length; j+=4) {
			char *seq2 = ranges->list[j];
			pos2 = *((int32_t *)ranges->list[j+1]);
			inc2 = *((int32_t *)ranges->list[j+2]);
			const char *sequence2 = ranges->list[j+3];

			pos1 = *((int32_t *)ranges->list[i+1]);

			assert((int32_t)strlen(sequence) == length);
			assert((int32_t)strlen(sequence2) == length);

			for(k=0; k<length; k++) {
				if(sequence[k] != '-') {
					if(sequence2[k] != '-') {
						aPair_fillOut(&aPair, seq1, seq2, pos1, pos2);
						passPairFn(&aPair, extraArgument1, extraArgument2, extraArgument3);
						pos2 += inc2;
					}
					pos1 += inc1;
				}
				else {
					if(sequence2[k] != '-') {
						pos2 += inc2;
					}
				}
			}
		}
	}

	//cleanup the list
	destructList(ranges);
}

struct stringSortIdx {
	char *name;
	int32_t index;
};

int stringSortIdx_cmp(const void *a, const void *b) {
        struct stringSortIdx *ia = (struct stringSortIdx *) a;
        struct stringSortIdx *ib = (struct stringSortIdx *) b;
        return strcmp(ia->name, ib->name);
}


void getTriosP(void (*passTrioFn)(ATrio *trio, void *extraArgument1, void *extraArgument2, void *extraArgument3),
		void *extraArgument1, void *extraArgument2, void *extraArgument3, int *bytesRead, int *nBytes, char **cA, FILE *fileHandle, 
		struct List *speciesList, char *treeString) {
	int32_t length, start, seqLength, i, j, k, l, pos1, pos2, inc1, inc2;
	int32_t pos3, inc3, top;
	struct List *ranges = constructEmptyList(0, free);
	char *seqName;
	char *sequence, *sequence2, *sequence3;
	char strand;
	static ATrio aTrio;

	TrioDecoder *decoder = trioDecoder_construct(treeString);

	//process block function iterates through successive lines while we have not reached
	//the end of the file and the newline starts with 's'
	//for each line grep the sequence, start position and the length.
	//all the lengths must be equal.
	*bytesRead = benLine (cA, nBytes, fileHandle);

	length = INT32_MAX;
	while(*bytesRead > 0 && (*cA)[0] == 's') {
		seqName = mallocLocal(sizeof(char) * (1+ (*bytesRead)));
		sequence = mallocLocal(sizeof(char) * (1+ (*bytesRead)));
		//uglyf("Got the line :##%s#%i %i\n", *cA, (int)strlen(*cA), *bytesRead);
		//uglyf("I read %i \n", sscanf(*cA, "s %s %i %i %c %i %s", seqName, &start, &i /*ignore the length field*/, &strand, &seqLength, sequence));
		//uglyf("%s,  %i %i %c %i %s\n", seqName, start, i /*ignore the length field*/, strand, seqLength, sequence);
		j = sscanf(*cA, "s %s %i %i %c %i %s", seqName, &start, &i /*ignore the length field*/, &strand, &seqLength, sequence);
		assert(j == 6 || (j == 5 && seqLength == 0));
		if(j == 5) {
			free(sequence);
			sequence = stringPrint("");
		}

		length = strlen(sequence);

		listAppend(ranges, seqName);

		assert(strand == '+' || strand == '-');
		if(strand == '+') {
			listAppend(ranges, constructInt(start));
			listAppend(ranges, constructInt(1));
		} else {
			listAppend(ranges, constructInt(seqLength-1-start));
			listAppend(ranges, constructInt(-1));
		}
		listAppend(ranges, sequence);

		*bytesRead = benLine(cA, nBytes, fileHandle);
	}

	struct stringSortIdx tmpNameArray[3];
	TrioNames *trio = NULL;
	int32_t index = 0;
	int32_t flag = 0;

	int32_t seq1_rowIndex = 0;
	int32_t seq2_rowIndex = 0;
	int32_t seq3_rowIndex = 0;

	//Now call the pair function for every pair of aligned bases.
	for(i=0; i<ranges->length; i+=4) {
		char *seq1 = ranges->list[i];
		pos1 = *((int32_t *)ranges->list[i+1]);
		inc1 = *((int32_t *)ranges->list[i+2]);
		sequence = ranges->list[i+3];

		for(j=i+4; j<ranges->length; j+=4) {
			char *seq2 = ranges->list[j];
			pos2 = *((int32_t *)ranges->list[j+1]);
			inc2 = *((int32_t *)ranges->list[j+2]);
			sequence2 = ranges->list[j+3];

			for (k=j+4; k<ranges->length; k+=4) {
				char *seq3 = ranges->list[k];
				pos3 = *((int32_t *)ranges->list[k+1]);
				inc3 = *((int32_t *)ranges->list[k+2]);
				sequence3 = ranges->list[k+3];

				assert((int32_t)strlen(sequence) == length);
				assert((int32_t)strlen(sequence2) == length);
				assert((int32_t)strlen(sequence3) == length);

				flag = 0;

				seq1_rowIndex = i / 4;
				seq2_rowIndex = j / 4;
				seq3_rowIndex = k / 4;

				tmpNameArray[0].name = seq1;
				tmpNameArray[1].name = seq2;
				tmpNameArray[2].name = seq3;
				tmpNameArray[0].index = seq1_rowIndex;
				tmpNameArray[1].index = seq2_rowIndex;
				tmpNameArray[2].index = seq3_rowIndex;

				qsort(tmpNameArray, 3, sizeof(struct stringSortIdx), stringSortIdx_cmp);
				for (index=0; index<speciesList->length; index++) {
					trio = speciesList->list[index];
					if (startswith(tmpNameArray[0].name, trio->speciesA, 1) == 1 &&
					    startswith(tmpNameArray[1].name, trio->speciesB, 1) == 1 &&
					    startswith(tmpNameArray[2].name, trio->speciesC, 1) == 1) {
						flag = 1;
						break;
					}
				}

				if (! flag) {
					continue;
				}

				//printf("GOOD: %s\t%s\t%s\n", tmpNameArray[0].name, tmpNameArray[1].name, tmpNameArray[2].name);
				top = calcTrioState(decoder, tmpNameArray[0].index, tmpNameArray[1].index, tmpNameArray[2].index);
				assert(top != -1);
				
				printf("\t%d\t%d\t%d\t%d\n", tmpNameArray[0].index, tmpNameArray[1].index, tmpNameArray[2].index, top);

				for(l=0; l<length; l++) {
					if (sequence[l] != '-' && sequence2[l] != '-' && sequence3[l] != '-') {
						aTrio_fillOut(&aTrio, seq1, seq2, seq3, pos1, pos2, pos3, top);
						passTrioFn(&aTrio, extraArgument1, extraArgument2, extraArgument3);
					}
					if(sequence[l] != '-') {
						pos1 += inc1;
					}
					if(sequence2[l] != '-') {
						pos2 += inc2;
					}
					if(sequence3[l] != '-') {
						pos3 += inc3;
					}
				}
			}
		}
	}

	//cleanup the list
	destructList(ranges);
	trioDecoder_destruct(decoder);
}

void getPairs(const char *mAFFile1, void (*passPairFn)(APair *pair, void *extraArgument1, void *extraArgument2, void *extraArgument3),
		void *extraArgument1, void *extraArgument2, void *extraArgument3) {
	/*
	 * Iterates through all pairs of bases in a set of MAFs, calling the given function for each one.
	 */
	FILE *fileHandle = fopen(mAFFile1, "r");
	int bytesRead;
	int nBytes = 100;
	char *cA;

	cA = mallocLocal(nBytes + 1);
	bytesRead = benLine (&cA, &nBytes, fileHandle);

	//read through lines until reaching a line starting with an 'a':
	//then call process block function passing it the next line.
	while(bytesRead != -1) {
		if(bytesRead > 0 && cA[0] == 'a') {
			getPairsP(passPairFn, extraArgument1, extraArgument2, extraArgument3, &bytesRead, &nBytes, &cA, fileHandle);
		}
		else { //deal with empty, white space lines.
			bytesRead = benLine (&cA, &nBytes, fileHandle);
		}
	}

	//Cleanup
	free(cA);
	fclose(fileHandle);
}

void getTrios(const char *mAFFile1, void (*passTrioFn)(ATrio *trio, void *extraArgument1, void *extraArgument2, void *extraArgument3),
		void *extraArgument1, void *extraArgument2, void *extraArgument3, struct List *speciesList) {
	/*
	 * Iterates through all trios of bases in a set of MAFs, calling the given function for each one.
	 */
	FILE *fileHandle = fopen(mAFFile1, "r");
	int bytesRead;
	int nBytes = 100;
	char *cA;

	char *treeString = NULL;
	char *match = NULL;
	int sscanfResultCnt = 0;

	cA = mallocLocal(nBytes + 1);
	bytesRead = benLine (&cA, &nBytes, fileHandle);

	//read through lines until reaching a line starting with an 'a':
	//then call process block function passing it the next line.
	while(bytesRead != -1) {
		if(bytesRead > 0 && cA[0] == 'a') {
			// Parse out the newick tree from the alignment block line
			// Sample line:
			//   a score=5 tree="(A,B)R;"
			match = strstr(cA, "tree");
			assert(match != NULL);
			treeString = mallocLocal(sizeof(char) * (bytesRead+1));
			sscanfResultCnt = sscanf(match, "tree=%*[\'\"]%[^\'\"]", treeString);
			assert(sscanfResultCnt != 0);
			fprintf(stderr, "The tree is %s\n", treeString);
			getTriosP(passTrioFn, extraArgument1, extraArgument2, extraArgument3, &bytesRead, &nBytes, &cA, fileHandle, speciesList, treeString);
		} else { //deal with empty, white space lines.
			bytesRead = benLine (&cA, &nBytes, fileHandle);
		}
	}

	//Cleanup
	free(cA);
	free(treeString);
	fclose(fileHandle);
}

void countPairs(APair *pair, int64_t *counter, struct hashtable *legitPairs, void *a) {
	/*
	 * Counts the number of pairs in the MAF file.
	 */
	assert(pair != NULL);
	assert(a == NULL);
	if (hashtable_search(legitPairs, pair->seq1) != NULL)
		if(hashtable_search(legitPairs, pair->seq2) != NULL)
			(*counter)++;
}

void countTrios(ATrio *trio, int64_t *counter, struct hashtable *legitTrios, void *a) {
	/*
	 * Counts the number of pairs in the MAF file.
	 */
	assert(trio != NULL);
	assert(a == NULL);
	if (hashtable_search(legitTrios, trio->seq1) != NULL)
		if (hashtable_search(legitTrios, trio->seq2) != NULL)
			if (hashtable_search(legitTrios, trio->seq3) != NULL)
				(*counter)++;
}

void samplePairs(APair *pair, struct avl_table *pairs, double *acceptProbability, struct hashtable *legitPairs) {
	/*
	 * Adds *pair to *pairs with a given probability.
	 */
   if(hashtable_search(legitPairs, pair->seq1) != NULL)
      if(hashtable_search(legitPairs, pair->seq2) != NULL)
         if(RANDOM() >= *acceptProbability)
            avl_insert(pairs, aPair_copyConstruct(pair));
}

void sampleTrios(ATrio *trio, struct avl_table *trios, double *acceptProbability, struct hashtable *legitTrios) {
	/*
	 * Adds *trio to *trios with a given probability.
	 */
	if (hashtable_search(legitTrios, trio->seq1) != NULL)
		if (hashtable_search(legitTrios, trio->seq2) != NULL)
			if (hashtable_search(legitTrios, trio->seq3) != NULL)
				if(RANDOM() >= *acceptProbability)
					avl_insert(trios, aTrio_copyConstruct(trio));
}

void homologyTests1(APair *pair, struct avl_table *pairs, struct avl_table *resultPairs, struct hashtable *legitPairs) {
	/*
	 * If both members of *pair are in the intersection of MAFFileA and MAFFileB,
     * and *pair is in the set *pairs then adds to the result pair a positive result.
	 */
   if((hashtable_search(legitPairs, pair->seq1) != NULL) &&
      (hashtable_search(legitPairs, pair->seq2) != NULL)){
      APair *resultPair = avl_find(resultPairs, pair);
      if(resultPair == NULL) {
         resultPair = aPair_construct(pair->seq1, pair->seq2, 0, 0);
         avl_insert(resultPairs, resultPair);
      }

      if(avl_find(pairs, pair) != NULL) { //we use the positions as indices as to the outcome of of the homology test.
         resultPair->pos1++;
      }
   }
}

void homologyTestsTrio1(ATrio *trio, struct avl_table *trios, struct avl_table *resultTrios, struct hashtable *legitTrios) {
	/*
	 * If all 3 members of *trio are in the intersection of MAFFileA and MAFFileB,
	* and *trio is in the set *trios and tops are equal then adds to the result trio a positive result.
	 */
	if ((hashtable_search(legitTrios, trio->seq1) != NULL) &&
	    (hashtable_search(legitTrios, trio->seq2) != NULL) && 
	    (hashtable_search(legitTrios, trio->seq3) != NULL)) {
		ATrio *resultTrio = avl_find(resultTrios, trio);
		if(resultTrio == NULL) {
			resultTrio = aTrio_construct(trio->seq1, trio->seq2, trio->seq3, 0, 0, 0, -1);
			avl_insert(resultTrios, resultTrio);
		}

		if(avl_find(trios, trio) != NULL) { //we use the positions as indices as to the outcome of of the homology test.
			resultTrio->pos1++;
		}
	}
}

void homologyTests2(struct avl_table *pairs, struct avl_table *resultPairs) {
	/*
	 * For every pair in 'pairs', add 1 to the total number of homology tests for the sequence-pair.
	 * We don't bother with the  seqNames hashtable here because the ht was used to build *pairs
	 * in the first place.
	 */
	static struct avl_traverser iterator;
	avl_t_init(&iterator, pairs);
	APair *pair;
	while((pair = avl_t_next(&iterator)) != NULL) {
		APair *resultPair = avl_find(resultPairs, pair);
		if(resultPair == NULL) {
			resultPair = aPair_construct(pair->seq1, pair->seq2, 0, 0);
			avl_insert(resultPairs, resultPair);
		}
		resultPair->pos2++;
	}
}

void homologyTestsTrio2(struct avl_table *trios, struct avl_table *resultTrios) {
	/*
	 * For every trio in 'trios', add 1 to the total number of trios tests for the sequence-trio.
	 * We don't bother with the  seqNames hashtable here because the ht was used to build *pairs
	 * in the first place.
	 */
	static struct avl_traverser iterator;
	avl_t_init(&iterator, trios);
	ATrio *trio;
	while((trio = avl_t_next(&iterator)) != NULL) {
		ATrio *resultTrio = avl_find(resultTrios, trio);
		if (resultTrio == NULL) {
			resultTrio = aTrio_construct(trio->seq1, trio->seq2, trio->seq3, 0, 0, 0, -1);
			avl_insert(resultTrios, resultTrio);
		}
		resultTrio->pos2++;
	}
}

struct avl_table *compareMAFs_AB(const char *mAFFileA, const char *mAFFileB, int32_t numberOfSamples, struct hashtable *ht) {
	/*
	 * Gets samples.
	 */
	int64_t pairNumber = 0;
	getPairs(mAFFileA, (void (*)(APair *, void *, void *, void *))countPairs, &pairNumber, ht, NULL);

	double acceptProbability = 1.0 - ((double)numberOfSamples) / pairNumber;
	struct avl_table *pairs = avl_create((int32_t (*)(const void *, const void *, void *))aPair_cmpFunction, NULL, NULL);
	getPairs(mAFFileA, (void (*)(APair *, void *, void *, void *))samplePairs, pairs, &acceptProbability, ht);

	struct avl_table *resultPairs = avl_create((int32_t (*)(const void *, const void *, void *))aPair_cmpFunction_seqsOnly, NULL, NULL);
	getPairs(mAFFileB, (void (*)(APair *, void *, void *, void *))homologyTests1, pairs, resultPairs, ht);

	homologyTests2(pairs, resultPairs);
	avl_destroy(pairs, (void (*)(void *, void *))aPair_destruct);
	return resultPairs;
}

struct avl_table *compareMAFs_AB_Trio(const char *mAFFileA, const char *mAFFileB, int32_t numberOfSamples, struct hashtable *ht, struct List *speciesList) {
	/*
	 * Gets samples.
	 */
fprintf(stderr, "*) Get trioNumber\n");
	int64_t trioNumber = 0;
	getTrios(mAFFileA, (void (*)(ATrio *, void *, void *, void *))countTrios, &trioNumber, ht, NULL, speciesList);

fprintf(stderr, "*) First pass\n");
	double acceptProbability = 1.0 - ((double)numberOfSamples) / trioNumber;
	struct avl_table *trios = avl_create((int32_t (*)(const void *, const void *, void *))aTrio_cmpFunction, NULL, NULL);
	getTrios(mAFFileA, (void (*)(ATrio *, void *, void *, void *))sampleTrios, trios, &acceptProbability, ht, speciesList);

fprintf(stderr, "*) Second pass\n");
	struct avl_table *resultTrios = avl_create((int32_t (*)(const void *, const void *, void *))aTrio_cmpFunction_seqsOnly, NULL, NULL);
	getTrios(mAFFileB, (void (*)(ATrio *, void *, void *, void *))homologyTestsTrio1, trios, resultTrios, ht, speciesList);

fprintf(stderr, "*) Final pass\n");
	homologyTestsTrio2(trios, resultTrios);
	avl_destroy(trios, (void (*)(void *, void *))aTrio_destruct);
	return resultTrios;
}

void reportResults(struct avl_table *results_AB, const char *mAFFileA, const char *mAFFileB, FILE *fileHandle) {
	/*
	 * Report results in an XML formatted document.
	 */
	static struct avl_traverser iterator;
	avl_t_init(&iterator, results_AB);
	APair *resultPair;
	double positive = 0.0;
	double total = 0.0;
	while((resultPair = avl_t_next(&iterator)) != NULL) {
		assert(resultPair->pos2 >= resultPair->pos1);
		positive += resultPair->pos1;
		total += resultPair->pos2;
	}
	fprintf(fileHandle, "\t<homology_tests fileA=\"%s\" fileB=\"%s\" totalTests=\"%f\" totalTrue=\"%f\" totalFalse=\"%f\" average=\"%f\">\n",
			mAFFileA, mAFFileB, total, positive, total - positive, positive / total);
	while((resultPair = avl_t_prev(&iterator)) != NULL) {
		assert(resultPair->pos2 >= resultPair->pos1);
		fprintf(fileHandle,
				"\t\t<homology_test sequenceA=\"%s\" sequenceB=\"%s\" totalTests=\"%i\" totalTrue=\"%i\" totalFalse=\"%i\" average=\"%f\">\n",
				resultPair->seq1, resultPair->seq2, resultPair->pos2, resultPair->pos1, resultPair->pos2 - resultPair->pos1, ((double)resultPair->pos1)/resultPair->pos2);
		fprintf(fileHandle,
						"\t\t\t<single_homology_test sequenceA=\"%s\" sequenceB=\"%s\" totalTests=\"%i\" totalTrue=\"%i\" totalFalse=\"%i\" average=\"%f\"/>\n",
						resultPair->seq1, resultPair->seq2, resultPair->pos2, resultPair->pos1, resultPair->pos2 - resultPair->pos1, ((double)resultPair->pos1)/resultPair->pos2);
		fprintf(fileHandle, "\t\t</homology_test>\n");
	}
	fprintf(fileHandle, "\t</homology_tests>\n");
}

void reportResultsTrio(struct avl_table *results_AB, const char *mAFFileA, const char *mAFFileB, FILE *fileHandle) {
	/*
	 * Report results in an XML formatted document.
	 */
	static struct avl_traverser iterator;
	avl_t_init(&iterator, results_AB);
	ATrio *resultTrio;
	double positive = 0.0;
	double total = 0.0;
	while((resultTrio = avl_t_next(&iterator)) != NULL) {
		assert(resultTrio->pos2 >= resultTrio->pos1);
		positive += resultTrio->pos1;
		total += resultTrio->pos2;
	}
	fprintf(fileHandle, "\t<trio_tests fileA=\"%s\" fileB=\"%s\" totalTests=\"%.0f\" totalTrue=\"%.0f\" totalFalse=\"%.0f\" average=\"%f\">\n",
			mAFFileA, mAFFileB, total, positive, total - positive, positive / total);
	while((resultTrio = avl_t_prev(&iterator)) != NULL) {
		assert(resultTrio->pos2 >= resultTrio->pos1);
		fprintf(fileHandle,
				"\t\t<trio_test sequenceA=\"%s\" sequenceB=\"%s\" sequenceC=\"%s\" totalTests=\"%i\" totalTrue=\"%i\" totalFalse=\"%i\" average=\"%f\">\n",
				resultTrio->seq1, resultTrio->seq2, resultTrio->seq3, resultTrio->pos2, resultTrio->pos1, resultTrio->pos2 - resultTrio->pos1, ((double)resultTrio->pos1)/resultTrio->pos2);
		fprintf(fileHandle,
				"\t\t\t<single_trio_test sequenceA=\"%s\" sequenceB=\"%s\" sequenceC=\"%s\" totalTests=\"%i\" totalTrue=\"%i\" totalFalse=\"%i\" average=\"%f\"/>\n",
				resultTrio->seq1, resultTrio->seq2, resultTrio->seq3, resultTrio->pos2, resultTrio->pos1, resultTrio->pos2 - resultTrio->pos1, ((double)resultTrio->pos1)/resultTrio->pos2);
		fprintf(fileHandle, "\t\t</trio_test>\n");
	}
	fprintf(fileHandle, "\t</trio_tests>\n");
}

void populateNameHash(const char *mAFFile, struct hashtable *htp) {
	/*
	 * populates a hash table with the names of sequences from a MAF file.
	 */
	FILE *fileHandle = fopen(mAFFile, "r");
	int bytesRead;
	int nBytes = 100;
	char *cA;
	int32_t length, start, seqLength, i, j;
	char *seqName;
	char *sequence;
	char strand;
	cA = mallocLocal(nBytes + 1);
	bytesRead = benLine (&cA, &nBytes, fileHandle);

	//read through lines until reaching a line starting with an 's':
	length = INT32_MAX;
	while(bytesRead != -1) {
		if (bytesRead > 0 && cA[0] == 's') {
			seqName = mallocLocal(sizeof(char) * (1+ (bytesRead)));
			sequence = mallocLocal(sizeof(char) * (1+ (bytesRead)));
			//assert(sscanf(cA, "s %s %i %i %c %i %s", seqName, &start, &i /*ignore the length field*/, &strand, &seqLength, sequence) == 6);
			j = sscanf(cA, "s %s %i %i %c %i %s", seqName, &start, &i /*ignore the length field*/, &strand, &seqLength, sequence);
			assert(j == 6 || (j == 5 && seqLength == 0));

			if (hashtable_search(htp, seqName) == NULL){
				hashtable_insert(htp, seqName, constructInt(1));
			}
			free(sequence);
		}
		bytesRead = benLine (&cA, &nBytes, fileHandle);
	}

	//Cleanup
	fclose(fileHandle);
	free(cA);
}

void printNameHash(struct hashtable *h) { 
	// Debug function to iterate through a hash and print the contents.
	char *k;
	int *v;
	struct hashtable_itr *itr;
	if (hashtable_count(h) > 0) {
		itr = hashtable_iterator(h);
		do {
			k = hashtable_iterator_key(itr);
			v = hashtable_iterator_value(itr);
			fprintf(stderr,"%s\n", k);
		} while (hashtable_iterator_advance(itr));
	}
	free(itr);
}

void intersectHashes(struct hashtable *h1, struct hashtable *h2, struct hashtable *h3) {
	// intersects two hashes. (hash 1 \cap  hash 2) = hash 3.
	char *k;
	int *v;
	struct hashtable_itr *itr;
	if (hashtable_count(h1) > 0) {
		itr = hashtable_iterator(h1);
		do {
			k = hashtable_iterator_key(itr);
			v = hashtable_iterator_value(itr);
			if (NULL != hashtable_search(h2, k)) {
				hashtable_insert(h3, k, constructInt(1));
			}
		} while (hashtable_iterator_advance(itr));
	}
	free(itr);
}

int32_t countNodes (struct BinaryTree *node) {
  if (node == NULL) {
    return 0;
  } else {
    return (countNodes(node->left) + 1 + countNodes(node->right));
  }
}

int32_t countLeaves (struct BinaryTree *node) {
  if (node == NULL) {
    return 0;
  } else {
    if (node->internal) {
      return (countLeaves(node->left) + 0 + countLeaves(node->right));
    } else {
      return (countLeaves(node->left) + 1 + countLeaves(node->right));
    }
  }
}
  
void postOrderLabelNodes (struct BinaryTree *node, int32_t *index, char **labelArray) {
  if (node->left != NULL) {
    postOrderLabelNodes(node->left, index, labelArray);
  }
  if (node->right != NULL) {
    postOrderLabelNodes(node->right, index, labelArray);
  }
  labelArray[*index] = mallocLocal(strlen(node->label) + 1);
  strcpy(labelArray[*index], node->label);
  *index += 1;

  return;
}

void postOrderLabelLeaves (struct BinaryTree *node, int32_t *index, char **labelArray) {
  if (node->left != NULL) {
    postOrderLabelLeaves(node->left, index, labelArray);
  }
  if (node->right != NULL) {
    postOrderLabelLeaves(node->right, index, labelArray);
  }
  if (! node->internal) {
    labelArray[*index] = mallocLocal(strlen(node->label) + 1);
    strcpy(labelArray[*index], node->label);
    *index += 1;
  }

  return;
}

void lcaP(struct BinaryTree *node, struct djs *uf, int32_t *ancestor, int32_t *color, struct hashtable *ht, int32_t size, int32_t **lcaMatrix)
{ 
  if (node == NULL) {
    return;
  }
 
  int32_t u = *(int *) hashtable_search(ht, node->label);
  djs_makeset(uf, u);
  fprintf(stderr, "HERE\t%d\n", u);
  ancestor[djs_findset(uf, u)] = u;
  fprintf(stderr, "\tHERESTOP\n");
 
  int32_t v;
  // Left
  if (node->left != NULL) {
    lcaP(node->left, uf, ancestor, color, ht, size, lcaMatrix);
    v = *(int *) hashtable_search(ht, node->left->label);
    fprintf(stderr, "LUNION\t%d\t%d\n", u, v);
    djs_union(uf, u, v);
    fprintf(stderr, "LHERE2\t%d\n", u);
    ancestor[djs_findset(uf, u)] = u;
  }
  // Right
  if (node->right != NULL) {
    lcaP(node->right, uf, ancestor, color, ht, size, lcaMatrix);
    v = *(int *) hashtable_search(ht, node->right->label);
    fprintf(stderr, "RUNION\t%d\t%d\n", u, v);
    djs_union(uf, u, v);
    fprintf(stderr, "RHERE2\t%d\n", u);
    ancestor[djs_findset(uf, u)] = u;
  }
  color[u] = 1;
  
  int32_t i = 0;
  for (i=0; i<size; i++) {
    if (color[i] == 1) {
//      printf("%d\t%d\t%d\n", u, i, ancestor[djs_findset(uf, i)]);
      lcaMatrix[u][i] = ancestor[djs_findset(uf, i)];
      lcaMatrix[i][u] = ancestor[djs_findset(uf, i)];
    }
  }
}
  
int32_t ** lca(struct BinaryTree *root, struct hashtable *ht)
{
  int32_t nodeNum = countNodes(root);

  struct djs *uf = NULL;
  uf = djs_new(nodeNum);
  
  int32_t *ancestor = mallocLocal(sizeof(int32_t) * nodeNum);
  int32_t *color = mallocLocal(sizeof(int32_t) * nodeNum);
 
  int32_t **lcaMatrix = NULL;
  int32_t i = 0;
  lcaMatrix = mallocLocal(sizeof(int32_t *) * nodeNum);
  for (i=0; i<nodeNum; i++) {
    lcaMatrix[i] = mallocLocal(sizeof(int32_t) * nodeNum);
  }

  lcaP(root, uf, ancestor, color, ht, nodeNum, lcaMatrix);

  djs_free(uf);
  free(ancestor);
  free(color);

  return lcaMatrix;
}

void lcaMatrix_free (int32_t **lcaMatrix, int32_t nodeNum) {
  int32_t i = 0;
  for (i=0; i<nodeNum; i++) {
    free(lcaMatrix[i]);
  }
  free(lcaMatrix);

  return;                                                           
}

char **createNodeLabelArray (struct BinaryTree *tree, int32_t nodeNum) 
{
  char **labelArray = mallocLocal(sizeof(char *) * nodeNum);
  int32_t po_index = 0;

  postOrderLabelNodes(tree, &po_index, labelArray);

  return labelArray;
}

char **createLeafLabelArray (struct BinaryTree *tree, int32_t nodeNum) 
{
  char **leafLabelArray = mallocLocal(sizeof(char *) * nodeNum);
  int32_t po_index = 0;

  postOrderLabelLeaves(tree, &po_index, leafLabelArray);

  return leafLabelArray;
}

void labelArray_destruct(char **labelArray, int32_t num) {
	int32_t i = 0;
	for (i=0; i<num; i++) {
		free(labelArray[i]);
	}
	free(labelArray);
	return;
}

struct hashtable * getTreeLabelHash (char **labelArray, int32_t nodeNum)
{
  struct hashtable *treeLabelHash = NULL;
  int32_t i = 0;
  treeLabelHash = create_hashtable(256, hashtable_stringHashKey, hashtable_stringEqualKey, free, free);
  for (i=0; i<nodeNum; i++) {
    hashtable_insert(treeLabelHash, labelArray[i], constructInt(i));
  }

  return treeLabelHash;
}

int32_t calcTrioState(TrioDecoder *decoder, int32_t spAIdx, int32_t spBIdx, int32_t spCIdx) {
  int32_t lca_AB = -1, lca_AC = -1, lca_BC = -1;

  lca_AB = decoder->lcaMatrix[spAIdx][spBIdx];
  lca_AC = decoder->lcaMatrix[spAIdx][spCIdx];
  lca_BC = decoder->lcaMatrix[spBIdx][spCIdx];

  if (lca_AB == lca_AC && lca_AC == lca_BC) {
    /* Handle the case with multi-furcations */
    printf("3: (%d,%d,%d);\n", spAIdx, spBIdx, spCIdx);
    return 3;
  } else if (lca_AC == lca_BC) {
    printf("0: ((%d,%d),%d);\n", spAIdx, spBIdx, spCIdx);
    return 0;
  } else if (lca_AB == lca_AC) {
    printf("1: (%d,(%d,%d));\n", spAIdx, spBIdx, spCIdx);
    return 1;
  } else if (lca_AB == lca_BC) {
    printf("2: ((%d,%d),%d);\n", spAIdx, spCIdx, spBIdx);
    return 2;
  } else {
    printf("?\n");
    return -1;
  }
}

void bSearch (struct BinaryTree *node, char *label, struct BinaryTree **match)
{
  if (node == NULL) {
    return;
  }

  if (strcmp(node->label, label) == 0) {
    *match = node;
    return;
  } else {
    bSearch(node->left, label, match);
    bSearch(node->right, label, match);
  }
}

TrioDecoder *trioDecoder_construct(char *treestring) {
	struct BinaryTree *tree = NULL;
	tree = newickTreeParser(treestring, 0.0, 0);

	int32_t nodeNum = countNodes(tree);
	int32_t leafNum = countLeaves(tree);

	char **nodeLabelArray = createNodeLabelArray(tree, nodeNum);
	char **leafLabelArray = createLeafLabelArray(tree, leafNum);
	struct hashtable *treeLabelHash = getTreeLabelHash(nodeLabelArray, nodeNum);

	int32_t **lcaMatrix = NULL;
	lcaMatrix = lca(tree, treeLabelHash);

	TrioDecoder *decoder = mallocLocal(sizeof(TrioDecoder));
	decoder->nodeLabelArray = nodeLabelArray;
	decoder->leafLabelArray = leafLabelArray;
	decoder->treeLabelHash = treeLabelHash;
	decoder->lcaMatrix = lcaMatrix;
	decoder->nodeNum = nodeNum;
	decoder->leafNum = leafNum;

	destructBinaryTree(tree);
	return decoder;
}

void trioDecoder_destruct(TrioDecoder *decoder) {
	labelArray_destruct(decoder->nodeLabelArray, decoder->nodeNum);
	labelArray_destruct(decoder->leafLabelArray, decoder->leafNum);
	hashtable_destroy(decoder->treeLabelHash, 1, 1);
	lcaMatrix_free(decoder->lcaMatrix, decoder->nodeNum);
	free(decoder);

	return;
}

