#include "eTreeExtras.h"

void eTreeX_countLeaves(ETree *node, void *data) {
	if (eTree_getChildNumber(node) == 0) {
		*((int *) data) += 1;
	}
	return;
}

void eTreeX_getLeafArray(ETree *node, void *data) {
	LeafPtrArray *p = (LeafPtrArray *) data;

	if (eTree_getChildNumber(node) == 0) {
		p->ptrArray[p->index] = (void *) node;
		p->index++;
	}
	return;
}

void eTreeX_postOrderTraversal(ETree *node, void (*func)(ETree *, void *), void *data) {
	int32_t i = 0;
	ETree *child = NULL;
	for (i=0; i<eTree_getChildNumber(node); i++) {
		child = eTree_getChild(node, i);
		eTreeX_postOrderTraversal(child, func, data);
	}
	func(node, data);
	return;
}

ETree *eTreeX_getTreeFromFile(char *treeFile) {
	int bytesRead;
	int nBytes = 100;
	char *cA = st_malloc(nBytes + 1);

	FILE *fileHandle = fopen(treeFile, "r");
	if (fileHandle == NULL) {
		fprintf(stderr, "Error: Can't open file '%s' for reading", treeFile);
		exit(1);
	}
	bytesRead = benLine(&cA, &nBytes, fileHandle);
	assert(bytesRead == -1); // Only read a single line from file
	fclose(fileHandle);

	ETree *tree = NULL;
	tree = eTree_parseNewickString(cA);

	free(cA);

	return tree;
}

LeafPtrArray *eTreeX_constructLeafPtrArray(int32_t leafCount) {
	LeafPtrArray *leafArray = NULL;
	leafArray = st_malloc(sizeof(LeafPtrArray));
	leafArray->ptrArray = st_malloc(sizeof(void *) * leafCount);
	leafArray->index = 0;

	return leafArray;
}

void eTreeX_destructLeafPtrArray(LeafPtrArray *leafArray) {
	free(leafArray->ptrArray);
	free(leafArray);

	return;
}
