/*
 * Copyright (C) 2012 by
 * Dent Earl (dearl@soe.ucsc.edu, dentearl@gmail.com)
 * ... and other members of the Reconstruction Team of David Haussler's
 * lab (BME Dept. UCSC).
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */
#include <assert.h>
#include <ctype.h>
#include <inttypes.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "common.h"
#include "CuTest.h"
#include "sharedMaf.h"

struct mafFileApi {
  // a mafFileApi struct provides an interface into a maf file.
  // Allows for easy reading of files in entirety or block by block via
  // functions
  uint64_t lineNumber; // last read line / wrote
  FILE *mfp; // maf file pointer
  char *filename; // filename of the maf
  char *lastLine; /* a temporary cache in case the header fails to have a blank
                   * line before the first alignment block.
                   */
};
struct mafLine {
  // a mafLine struct is a single line of a mafBlock
  char *line; // the entire line, unparsed
  uint64_t lineNumber; // line number in the maf file
  char type; // either a, s, i, q, e, h, f where h is header (an internal code)
  char *species; // species name
  uint64_t start;
  uint64_t length;
  char strand;
  uint64_t sourceLength;
  char *sequence; // sequence field
  uint64_t sequenceFieldLength;
  struct mafLine *next;
};
struct mafBlock {
  // a mafBlock struct contains a maf block as a linked list
  // and itself can be part of a mafBlock linked list.
  mafLine_t *headLine;
  mafLine_t *tailLine;
  uint64_t lineNumber; // line number of start of the block (i.e. the `a' line)
  uint64_t numberOfLines; // number of mafLine_t structures in the *headLine list
  uint64_t numberOfSequences;
  uint64_t sequenceFieldLength;
  struct mafBlock *next;
};
static bool maf_isBlankLine(char *s) {
  // return true if line is only whitespaces
  size_t n = strlen(s);
  for (size_t i = 0; i < n; ++i) {
    if (!isspace(*(s + i))) {
      return false;
    }
  }
  return true;
}
static void maf_checkForPrematureMafEnd(char *filename, int status, char *line) {
  if (status == -1) {
    free(line);
    line = NULL;
    fprintf(stderr, "Error, premature end to maf file: %s\n", filename);
    exit(EXIT_FAILURE);
  }
}
static void maf_failBadFormat(uint64_t lineNumber, char *errorMessage) {
  fprintf(stderr, "The maf sequence at line %" PRIi64 " is incorrectly formatted: %s\n",
          lineNumber, errorMessage);
  exit(EXIT_FAILURE);
}
mafLine_t* maf_newMafLine(void) {
  mafLine_t *ml = (mafLine_t *) de_malloc(sizeof(*ml));
  ml->line = NULL;
  ml->lineNumber = 0;
  ml->type = '\0';
  ml->species = NULL;
  ml->start = 0;
  ml->length = 0;
  ml->strand = 0;
  ml->sourceLength = 0;
  ml->sequence = NULL;
  ml->sequenceFieldLength = 0;
  ml->next = NULL;
  return ml;
}
mafLine_t* maf_copyMafLineList(mafLine_t *orig) {
  // create and return a copy of orig, a mafLine_t linked list
  if (orig == NULL) {
    return NULL;
  }
  mafLine_t *head = NULL, *ml = NULL, *tmp = NULL;
  while (orig != NULL) {
    if (head == NULL) {
      tmp = maf_copyMafLine(orig);
      ml = tmp;
      head = ml;
    } else {
      tmp = maf_copyMafLine(orig);
      maf_mafLine_setNext(ml, tmp);
      ml = maf_mafLine_getNext(ml);
    }
    orig = maf_mafLine_getNext(orig);
  }
  return head;
}
mafLine_t* maf_copyMafLine(mafLine_t *orig) {
  // create and return a copy of a single mafLine_t structure
  if (orig == NULL) {
    return NULL;
  }
  mafLine_t *ml = maf_newMafLine();
  if (orig->line != NULL) {
    ml->line = de_strdup(orig->line);
  }
  ml->lineNumber = orig->lineNumber;
  ml->type = orig->type;
  if (orig->species != NULL) {
    ml->species = de_strdup(orig->species);
  }
  ml->start = orig->start;
  ml->length = orig->length;
  ml->strand = orig->strand;
  ml->sourceLength = orig->sourceLength;
  if (orig->sequence != NULL) {
    ml->sequence = de_strdup(orig->sequence);
  }
  ml->sequenceFieldLength = orig->sequenceFieldLength;
  return ml;
}
mafBlock_t* maf_newMafBlockListFromString(const char *s, uint64_t lineNumber) {
  // given a string, walk through and create a mafBlock_t linked list
  // for all maf blocks in the string
  mafBlock_t *head = NULL, *mb = NULL, *tmp = NULL;
  char *block_s = NULL;
  size_t len, start = 0, stop;
  len = strlen(s);
  for (stop = start + 1; stop < len; ++stop) {
    if (s[stop] == 'a' && s[stop - 1] == '\n') {
      block_s = (char *) de_malloc(stop - start + 2);
      strncpy(block_s, s + start, stop - start + 1);
      block_s[stop - start + 1] = '\0';
      tmp = maf_newMafBlockFromString(block_s, lineNumber);
      if (head == NULL) {
        mb = tmp;
        head = mb;
      } else {
        maf_mafBlock_setNext(mb, tmp);
        mb = maf_mafBlock_getNext(mb);
      }
      free(block_s);
      block_s = NULL;
      start = stop;
    }
  }
  block_s = (char *) de_malloc(stop - start + 2);
  strncpy(block_s, s + start, stop - start + 1);
  block_s[stop - start + 1] = '\0';
  tmp = maf_newMafBlockFromString(block_s, lineNumber);
  if (head == NULL) {
    mb = tmp;
    head = mb;
  } else {
    maf_mafBlock_setNext(mb, tmp);
    mb = maf_mafBlock_getNext(mb);
  }
  free(block_s);
  block_s = NULL;
  return head;
}
mafBlock_t* maf_newMafBlockFromString(const char *s, uint64_t lineNumber) {
  if (s[0] != 'a') {
    char *error = de_malloc(kMaxStringLength);
    sprintf(error,
            "Unable to create maf block from input, "
            "first line does not start with 'a': %s", s);
    maf_failBadFormat(lineNumber, error);
  }
  mafBlock_t* mb = maf_newMafBlock();
  mafLine_t* ml = NULL;
  maf_mafBlock_setLineNumber(mb, lineNumber);
  char *cline = (char *) de_malloc(strlen(s) + 1);
  strcpy(cline, s);
  char *cline_orig = cline;
  char **ptr = &cline;
  char *tkn = NULL;
  tkn = de_strtok(ptr, '\n');
  maf_mafBlock_incrementLineNumber(mb);
  maf_mafBlock_incrementNumberOfLines(mb);
  ml = maf_newMafLineFromString(tkn, lineNumber++);
  maf_mafBlock_setHeadLine(mb, ml);
  free(tkn);
  tkn = NULL;
  tkn = de_strtok(ptr, '\n');
  while (tkn != NULL) {
    maf_mafBlock_incrementLineNumber(mb);
    maf_mafBlock_incrementNumberOfLines(mb);
    maf_mafLine_setNext(ml, maf_newMafLineFromString(tkn, lineNumber++));
    ml = maf_mafLine_getNext(ml);
    if (maf_mafLine_getType(ml) == 's') {
      maf_mafBlock_incrementNumberOfSequences(mb);
      mb->sequenceFieldLength = maf_mafLine_getSequenceFieldLength(ml);
    }
    maf_mafBlock_setTailLine(mb, ml);
    free(tkn);
    tkn = NULL;
    tkn = de_strtok(ptr, '\n');
  }
  free(cline_orig);
  cline_orig = NULL;
  return mb;
}
mafLine_t* maf_newMafLineFromString(const char *s, uint64_t lineNumber) {
  extern const int kMaxStringLength;
  mafLine_t *ml = (mafLine_t *) de_malloc(sizeof(*ml));
  char *copy = (char *) de_malloc(strlen(s) + 1);
  char *cline = (char *) de_malloc(strlen(s) + 1);
  strcpy(copy, s);
  strcpy(cline, s);
  ml->next = NULL;
  ml->line = copy;
  ml->lineNumber = lineNumber;
  ml->species = NULL;
  ml->sequence = NULL;
  ml->start = 0;
  ml->length = 0;
  ml->sourceLength = 0;
  ml->strand = 0;
  ml->type = ml->line[0];
  if (ml->type != 's') {
    free(cline);
    cline = NULL;
    return ml;
  }
  char *tkn = NULL;
  tkn = strtok(cline, " \t");
  if (tkn == NULL) {
    free(cline);
    cline = NULL;
    char *error = de_malloc(kMaxStringLength);
    sprintf(error, "Unable to separate line on tabs and spaces at line definition field:\n%s", s);
    maf_failBadFormat(lineNumber, error);
  }
  tkn = strtok(NULL, " \t"); // name field
  if (tkn == NULL) {
    free(cline);
    cline = NULL;
    maf_failBadFormat(lineNumber, "Unable to separate line on tabs and spaces at name field.");
  }
  char *species = (char *) de_malloc(strlen(tkn) + 1);
  strcpy(species, tkn);
  ml->species = species;
  tkn = strtok(NULL, " \t"); // start position
  if (tkn == NULL) {
    free(cline);
    cline = NULL;
    maf_failBadFormat(lineNumber, "Unable to separate line on tabs and spaces at start position field.");
  }
  ml->start = strtoul(tkn, NULL, 10);
  tkn = strtok(NULL, " \t"); // length position
  if (tkn == NULL){
    free(cline);
    cline = NULL;
    maf_failBadFormat(lineNumber, "Unable to separate line on tabs and spaces at length position field.");
  }
  ml->length = strtoul(tkn, NULL, 10);
  tkn = strtok(NULL, " \t"); // strand
  if (tkn == NULL) {
    free(cline);
    cline = NULL;
    maf_failBadFormat(lineNumber, "Unable to separate line on tabs and spaces at strand field.");
  }
  if (tkn[0] != '-' && tkn[0] != '+') {
    char *error = (char*) de_malloc(kMaxStringLength);
    sprintf(error, "Strand must be either + or -, not %c.", tkn[0]);
    maf_failBadFormat(lineNumber, error);
  }
  ml->strand = tkn[0];
  tkn = strtok(NULL, " \t"); // source length position
  if (tkn == NULL) {
    free(cline);
    cline = NULL;
    maf_failBadFormat(lineNumber, "Unable to separate line on tabs and spaces at source length field.");
  }
  ml->sourceLength = strtoul(tkn, NULL, 10);
  tkn = strtok(NULL, " \t"); // sequence field
  if (tkn == NULL) {
    free(cline);
    cline = NULL;
    char *error = de_malloc(kMaxStringLength);
    sprintf(error, "Unable to separate line on tabs and spaces at sequence field:\n%s", s);
    maf_failBadFormat(lineNumber, error);
  }
  char *seq = (char *) de_malloc(strlen(tkn) + 1);
  strcpy(seq, tkn);
  ml->sequence = seq;
  ml->sequenceFieldLength = strlen(ml->sequence);
  free(cline);
  cline = NULL;
  return ml;
}
mafBlock_t* maf_newMafBlock(void) {
  mafBlock_t *mb = (mafBlock_t *) de_malloc(sizeof(*mb));
  mb->next = NULL;
  mb->headLine = NULL;
  mb->tailLine = NULL;
  mb->lineNumber = 0;
  mb->numberOfSequences = 0;
  mb->numberOfLines = 0;
  mb->sequenceFieldLength = 0;
  return mb;
}
mafBlock_t* maf_copyMafBlockList(mafBlock_t *orig) {
  if (orig == NULL) {
    return NULL;
  }
  mafBlock_t *head = NULL, *mb = NULL, *tmp = NULL;
  while (orig != NULL) {
    if (head == NULL) {
      tmp = maf_copyMafBlock(orig);
      mb = tmp;
      head = mb;
    } else {
      tmp = maf_copyMafBlock(orig);
      maf_mafBlock_setNext(mb, tmp);
      mb = maf_mafBlock_getNext(mb);
    }
    orig = maf_mafBlock_getNext(orig);
  }
  return head;
}
mafBlock_t* maf_copyMafBlock(mafBlock_t *orig) {
  // copies a SINGLE maf block. Does not copy down the linked list of blocks.
  if (orig == NULL) {
    return NULL;
  }
  mafBlock_t *mb = maf_newMafBlock();
  // copy mafLine_t linked list
  mb->headLine = maf_copyMafLineList(orig->headLine);
  // record the mafBlock tail line
  mafLine_t *ml = mb->headLine;
  while (ml != NULL) {
    mb->tailLine = ml;
    ml = ml->next;
  }
  // copy attributes
  mb->lineNumber = orig->lineNumber;
  mb->numberOfSequences = orig->numberOfSequences;
  mb->numberOfLines = orig->numberOfLines;
  mb->sequenceFieldLength = orig->sequenceFieldLength;
  return mb;
}
mafFileApi_t* maf_newMfa(const char *filename, char const *mode) {
  mafFileApi_t *mfa = (mafFileApi_t *) de_malloc(sizeof(*mfa));
  mfa->lineNumber = 0;
  mfa->lastLine = NULL;
  mfa->mfp = de_fopen(filename, mode);
  mfa->filename = de_strdup(filename);
  return mfa;
}
void maf_destroyMafLineList(mafLine_t *ml) {
  // walk down a mafLine_t following the ->next pointers, search and destroy
  if (ml == NULL) {
    return;
  }
  mafLine_t *tmp = NULL;
  while(ml != NULL) {
    tmp = ml;
    ml = ml->next;
    free(tmp->line);
    tmp->line = NULL;
    if (tmp->species != NULL) {
      // you can have a maf line without a species member
      free(tmp->species);
      tmp->species = NULL;
    }
    if (tmp->sequence != NULL) {
      // you can have a maf line without a sequence member
      free(tmp->sequence);
      tmp->sequence = NULL;
    }
    free(tmp);
    tmp = NULL;
  }
}
void maf_destroyMafBlockList(mafBlock_t *mb) {
  if (mb == NULL) {
    return;
  }
  mafBlock_t *tmp = NULL;
  while(mb != NULL) {
    tmp = mb;
    mb = mb->next;
    if (tmp->headLine != NULL)
      maf_destroyMafLineList(tmp->headLine);
    free(tmp);
    tmp = NULL;
  }
}
void maf_destroyMfa(mafFileApi_t *mfa) {
  if (mfa->mfp != NULL) {
    fclose(mfa->mfp);
    mfa->mfp = NULL;
  }
  free(mfa->lastLine);
  mfa->lastLine = NULL;
  free(mfa->filename);
  mfa->filename = NULL;
  free(mfa);
  mfa = NULL;
}
char* maf_mafFileApi_getFilename(mafFileApi_t *mfa) {
  return mfa->filename;
}
uint64_t maf_mafFileApi_getLineNumber(mafFileApi_t *mfa) {
  return mfa->lineNumber;
}
mafLine_t* maf_mafBlock_getHeadLine(mafBlock_t *mb) {
  return mb->headLine;
}
mafLine_t* maf_mafBlock_getTailLine(mafBlock_t *mb) {
  return mb->tailLine;
}
uint64_t maf_mafBlock_getLineNumber(mafBlock_t *mb) {
  return mb->lineNumber;
}
uint64_t maf_mafBlock_getNumberOfSequences(mafBlock_t *mb) {
  return mb->numberOfSequences;
}
mafBlock_t* maf_mafBlock_getNext(mafBlock_t *mb) {
  return mb->next;
}
char** maf_mafBlock_getSequenceMatrix(mafBlock_t *mb, unsigned n, unsigned m) {
  // currently this is not stored and must be built
  // should return a matrix containing the alignment, one row per sequence
  char** matrix = NULL;
  matrix = (char**) de_malloc(sizeof(char*) * n);
  unsigned i;
  for (i = 0; i < n; ++i) {
    matrix[i] = (char*) de_malloc(sizeof(char) * (m + 1));
  }
  mafLine_t *ml = maf_mafBlock_getHeadLine(mb);
  i = 0;
  while (ml != NULL) {
    while (ml != NULL && maf_mafLine_getType(ml) != 's') {
      ml = maf_mafLine_getNext(ml);
    }
    if (ml == NULL) {
      break;
    }
    strncpy(matrix[i], maf_mafLine_getSequence(ml), m);
    matrix[i++][m] = '\0';
    ml = maf_mafLine_getNext(ml);
  }
  return matrix;
}
void maf_mafBlock_destroySequenceMatrix(char **mat, unsigned n) {
  // currently this is not stored and must be built
  // should return a matrix containing the alignment, one row per sequence
  for (unsigned i = 0; i < n; ++i) {
    free(mat[i]);
    mat[i] = NULL;
  }
  free(mat);
  mat = NULL;
}
char* maf_mafBlock_getStrandArray(mafBlock_t *mb) {
  // currently this is not stored and must be built
  // should return a char array containing an in-order list of strandedness
  // for all sequence lines. Either + or - char permitted.
  if (maf_mafBlock_getNumberOfSequences(mb) == 0) {
    return NULL;
  }
  char *a = (char*) de_malloc(sizeof(*a) * (maf_mafBlock_getNumberOfSequences(mb) + 1));
  mafLine_t *ml = maf_mafBlock_getHeadLine(mb);
  unsigned i = 0;
  while (ml != NULL) {
    if (ml->type == 's')
      a[i++] = ml->strand;
    ml = ml->next;
  }
  a[i] = '\0';
  return a;
}
mafLine_t** maf_mafBlock_getMafLineArray_seqOnly(mafBlock_t *mb) {
  // currently this is not stored and must be built
  // should return an array of mafLine_t pointers to all sequences in mb
  if (maf_mafBlock_getNumberOfSequences(mb) == 0) {
    return NULL;
  }
  mafLine_t **a = (mafLine_t**) de_malloc(sizeof(*a) * (maf_mafBlock_getNumberOfSequences(mb)));
  mafLine_t *ml = maf_mafBlock_getHeadLine(mb);
  unsigned i = 0;
  while (ml != NULL) {
    if (ml->type == 's')
      a[i++] = ml;
    ml = ml->next;
  }
  return a;
}
int* maf_mafBlock_getStrandIntArray(mafBlock_t *mb) {
  // currently this is not stored and must be built
  // should return an int array containing an in-order list of strandedness
  // for all sequence lines. Either 1 or -1
  int *a = (int*) de_malloc(sizeof(*a) * maf_mafBlock_getNumberOfSequences(mb));
  mafLine_t *ml = maf_mafBlock_getHeadLine(mb);
  unsigned i = 0;
  while (ml != NULL) {
    if (ml->type == 's') {
      if (ml->strand == '+') {
        a[i++] = 1;
      } else {
        a[i++] = -1;
      }
    }
    ml = ml->next;
  }
  return a;
}
uint64_t* maf_mafBlock_getStartArray(mafBlock_t *mb) {
  // currently this is not stored and must be built
  // should return a uint64_t array containing an in-order list of source lengths
  // for all sequence lines.
  uint64_t *a = (uint64_t*) de_malloc(sizeof(*a) * maf_mafBlock_getNumberOfSequences(mb));
  mafLine_t *ml = maf_mafBlock_getHeadLine(mb);
  unsigned i = 0;
  while (ml != NULL) {
    if (ml->type == 's') {
      a[i++] = ml->start;
    }
    ml = ml->next;
  }
  return a;
}
uint64_t* maf_mafBlock_getPosCoordStartArray(mafBlock_t *mb) {
  // currently this is not stored and must be built
  // should return a uint64_t array containing an in-order list of the start position in
  // positive coordinates.
  uint64_t *a = (uint64_t*) de_malloc(sizeof(*a) * maf_mafBlock_getNumberOfSequences(mb));
  mafLine_t *ml = maf_mafBlock_getHeadLine(mb);
  unsigned i = 0;
  while (ml != NULL) {
    if (ml->type == 's') {
      if (ml->strand == '+')
        a[i++] = ml->start;
      else
        a[i++] = ml->sourceLength - ml->start - 1;
    }
    ml = ml->next;
  }
  return a;
}
uint64_t* maf_mafBlock_getPosCoordLeftArray(mafBlock_t *mb) {
  // currently this is not stored and must be built
  // should return a uint64_t array containing an in-order list of the left-most positon
  // of the block in positive coordinates.
  uint64_t *a = (uint64_t*) de_malloc(sizeof(*a) * maf_mafBlock_getNumberOfSequences(mb));
  mafLine_t *ml = maf_mafBlock_getHeadLine(mb);
  unsigned i = 0;
  while (ml != NULL) {
    if (ml->type == 's') {
      if (ml->strand == '+')
        a[i++] = ml->start;
      else
        a[i++] = ml->sourceLength - (ml->start + ml->length);
    }
    ml = ml->next;
  }
  return a;
}
uint64_t* maf_mafBlock_getSourceLengthArray(mafBlock_t *mb) {
  // currently this is not stored and must be built
  // should return a uint64_t array containing an in-order list of positive
  // coordinate start positions for all sequence lines.
  uint64_t *a = (uint64_t*) de_malloc(sizeof(*a) * maf_mafBlock_getNumberOfSequences(mb));
  mafLine_t *ml = maf_mafBlock_getHeadLine(mb);
  unsigned i = 0;
  while (ml != NULL) {
    if (ml->type == 's') {
      a[i++] = ml->sourceLength;
    }
    ml = ml->next;
  }
  return a;
}
uint64_t* maf_mafBlock_getSequenceLengthArray(mafBlock_t *mb) {
  // currently this is not stored and must be built
  // should return a uint64_t array containing an in-order list of
  // sequence length field values
  uint64_t *a = (uint64_t*) de_malloc(sizeof(*a) * maf_mafBlock_getNumberOfSequences(mb));
  mafLine_t *ml = maf_mafBlock_getHeadLine(mb);
  unsigned i = 0;
  while (ml != NULL) {
    if (ml->type == 's') {
      a[i++] = ml->length;
    }
    ml = ml->next;
  }
  return a;
}
char** maf_mafBlock_getSpeciesArray(mafBlock_t *mb) {
  // currently this is not stored and must be built
  // should return an array of char pointers containing an in-order list of
  // sequence name fields for all sequences.
  char** m = NULL;
  m = (char**) de_malloc(sizeof(char*) * maf_mafBlock_getNumberOfSequences(mb));
  mafLine_t *ml = maf_mafBlock_getHeadLine(mb);
  unsigned i = 0;
  while (ml != NULL) {
    if (ml->type == 's')
      m[i++] = de_strdup(ml->species);
    ml = ml->next;
  }
  return m;
}
char* maf_mafLine_getLine(mafLine_t *ml) {
  return ml->line;
}
uint64_t maf_mafLine_getLineNumber(mafLine_t *ml) {
  return ml->lineNumber;
}
char maf_mafLine_getType(mafLine_t *ml) {
  return ml->type;
}
char* maf_mafLine_getSpecies(mafLine_t *ml) {
  return ml->species;
}
uint64_t maf_mafLine_getStart(mafLine_t *ml) {
  return ml->start;
}
uint64_t maf_mafLine_getLength(mafLine_t *ml) {
  return ml->length;
}
char maf_mafLine_getStrand(mafLine_t *ml) {
  return ml->strand;
}
uint64_t maf_mafLine_getSourceLength(mafLine_t *ml) {
  return ml->sourceLength;
}
char* maf_mafLine_getSequence(mafLine_t *ml) {
  return ml->sequence;
}
uint64_t maf_mafLine_getSequenceFieldLength(mafLine_t *ml) {
  return ml->sequenceFieldLength;
}
mafLine_t* maf_mafLine_getNext(mafLine_t *ml) {
  return ml->next;
}
uint64_t maf_mafBlock_getSequenceFieldLength(mafBlock_t *mb) {
  return mb->sequenceFieldLength;
}
unsigned maf_mafBlock_getNumberOfBlocks(mafBlock_t *b) {
  unsigned n = 0;
  while (b != NULL) {
    ++n;
    b = b->next;
  }
  return n;
}
unsigned umax(unsigned a, unsigned b) {
  return (a > b ? a : b);
}
bool maf_mafBlock_containsSequence(mafBlock_t *mb) {
  if (mb->numberOfSequences > 0) {
    return true;
  } else {
    return false;
  }
}
uint64_t maf_mafBlock_getNumberOfLines(mafBlock_t *mb) {
  // count the number of mafLine_t lines in a headLine list
  return mb->numberOfLines;
}
uint64_t maf_mafLine_getNumberOfSequences(mafLine_t *ml) {
  // count the number of actual sequence lines in a mafLine_t list
  uint64_t s = 0;
  while (ml != NULL) {
    if (ml->type == 's')
      ++s;
    ml = ml->next;
  }
  return s;
}
uint64_t maf_mafLine_getPositiveCoord(mafLine_t *ml) {
  // return the start field coordinate in postive zero based coordinates.
  // NOTE THAT FOR - STRANDS, THIS COORDINATE WILL BE THE RIGHT-MOST (END POINT)
  // OF THE SEQUENCE. TO GET THE LEFT-MOST (START POINT) YOU WOULD NEED TO SUBTRACT
  //
  if (ml->strand == '+') {
    return ml->start;
  } else {
    return ml->sourceLength - (ml->start + 1);
  }
}
uint64_t maf_mafLine_getPositiveLeftCoord(mafLine_t *ml) {
  // return the left most coordinate in postive zero based coordinates.
  // for - strands this includes the length of the sequence.
  if (ml->strand == '+') {
    return ml->start;
  } else {
    return ml->sourceLength - (ml->start + ml->length);
  }
}
void maf_mafBlock_setHeadLine(mafBlock_t *mb, mafLine_t *ml) {
  mb->headLine = ml;
}
void maf_mafBlock_setTailLine(mafBlock_t *mb, mafLine_t *ml) {
  mb->tailLine = ml;
}
void maf_mafBlock_setNumberOfSequences(mafBlock_t *mb, uint64_t n) {
  mb->numberOfSequences = n;
}
void maf_mafBlock_incrementNumberOfSequences(mafBlock_t *mb) {
  ++(mb->numberOfSequences);
}
void maf_mafBlock_decrementNumberOfSequences(mafBlock_t *mb) {
  --(mb->numberOfSequences);
}
void maf_mafBlock_setNumberOfLines(mafBlock_t *mb, uint64_t n) {
  mb->numberOfLines = n;
}
void maf_mafBlock_incrementNumberOfLines(mafBlock_t *mb) {
  ++(mb->numberOfLines);
}
void maf_mafBlock_decrementNumberOfLines(mafBlock_t *mb) {
  --(mb->numberOfLines);
}
void maf_mafBlock_setLineNumber(mafBlock_t *mb, uint64_t n) {
  mb->lineNumber = n;
}
void maf_mafBlock_incrementLineNumber(mafBlock_t *mb) {
  ++(mb->lineNumber);
}
void maf_mafBlock_decrementLineNumber(mafBlock_t *mb) {
  --(mb->lineNumber);
}
void maf_mafBlock_setSequenceFieldLength(mafBlock_t *mb, uint64_t sfl) {
  mb->sequenceFieldLength = sfl;
}
void maf_mafBlock_setNext(mafBlock_t *mb, mafBlock_t *next) {
  mb->next = next;
}
void maf_mafLine_setLine(mafLine_t *ml, char *line) {
  ml->line = line;
}
void maf_mafLine_setLineNumber(mafLine_t *ml, uint64_t n) {
  ml->lineNumber = n;
}
void maf_mafLine_setType(mafLine_t *ml, char c) {
  ml->type = c;
}
void maf_mafLine_setSpecies(mafLine_t *ml, char *s) {
  ml->species = s;
}
void maf_mafLine_setStrand(mafLine_t *ml, char c) {
  ml->strand = c;
}
void maf_mafLine_setStart(mafLine_t *ml, uint64_t n) {
  ml->start = n;
}
void maf_mafLine_setLength(mafLine_t *ml, uint64_t n) {
  ml->length = n;
}
void maf_mafLine_setSourceLength(mafLine_t *ml, uint64_t n) {
  ml->sourceLength = n;
}
void maf_mafLine_setSequence(mafLine_t *ml, char *s) {
  ml->sequence = s;
  ml->sequenceFieldLength = strlen(ml->sequence);
}
void maf_mafLine_setNext(mafLine_t *ml, mafLine_t *next) {
  ml->next = next;
}
mafBlock_t* maf_readBlockHeader(mafFileApi_t *mfa) {
  extern const int kMaxStringLength;
  int64_t n = kMaxStringLength;
  char *line = (char*) de_malloc(n);
  mafBlock_t *header = maf_newMafBlock();
  int status = de_getline(&line, &n, mfa->mfp);
  bool validHeader = false;
  ++(mfa->lineNumber);
  maf_checkForPrematureMafEnd(maf_mafFileApi_getFilename(mfa), status, line);
  if (strncmp(line, "track", 5) == 0) {
    // possible first line of a maf
    validHeader = true;
    mafLine_t *ml = maf_newMafLine();
    char *copy = (char *) de_malloc(n + 1); // freed in destroy lines
    strcpy(copy, line);
    ml->line = copy;
    ml->type = 'h';
    ml->lineNumber = mfa->lineNumber;
    header->headLine = ml;
    header->tailLine = ml;
    status = de_getline(&line, &n, mfa->mfp);
    ++(mfa->lineNumber);
    header->lineNumber = mfa->lineNumber;
    ++(header->numberOfLines);
    maf_checkForPrematureMafEnd(maf_mafFileApi_getFilename(mfa), status, line);
  }
  if (strncmp(line, "##maf", 5) == 0) {
    // possible first or second line of maf
    validHeader = true;
    mafLine_t *ml = maf_newMafLine();
    char *copy = (char *) de_malloc(n + 1); // freed in destroy lines
    strcpy(copy, line);
    ml->line = copy;
    ml->type = 'h';
    ml->lineNumber = mfa->lineNumber;
    if (header->headLine == NULL) {
      header->headLine = ml;
      header->tailLine = ml;
    } else {
      header->headLine->next = ml;
      header->tailLine = ml;
    }
    status = de_getline(&line, &n, mfa->mfp);
    ++(mfa->lineNumber);
    header->lineNumber = mfa->lineNumber;
    ++(header->numberOfLines);
    maf_checkForPrematureMafEnd(maf_mafFileApi_getFilename(mfa), status, line);
  }
  if (!validHeader) {
    fprintf(stderr, "Error, maf file %s does not contain a valid header!\n", mfa->filename);
    exit(EXIT_FAILURE);
  }
  mafLine_t *thisMl = header->tailLine;
  while(line[0] != 'a' && !maf_isBlankLine(line)) {
    // eat up the file until we hit the first alignment block
    mafLine_t *ml = maf_newMafLine();
    char *copy = (char *) de_malloc(n + 1); // freed in destroy lines
    strcpy(copy, line);
    ml->line = copy;
    ml->type = 'h';
    ml->lineNumber = mfa->lineNumber;
    thisMl->next = ml;
    thisMl = ml;
    header->tailLine = thisMl;
    status = de_getline(&line, &n, mfa->mfp);
    ++(mfa->lineNumber);
    header->lineNumber = mfa->lineNumber;
    ++(header->numberOfLines);
    maf_checkForPrematureMafEnd(maf_mafFileApi_getFilename(mfa), status, line);
  }
  if (line[0] == 'a') {
    // stuff this line in ->lastLine for processesing
    char *copy = (char *) de_malloc(n + 1); // freed in destroy lines
    strcpy(copy, line);
    mfa->lastLine = copy;
  }
  free(line);
  line = NULL;
  return header;
}
mafBlock_t* maf_readBlockBody(mafFileApi_t *mfa) {
  extern const int kMaxStringLength;
  mafBlock_t *thisBlock = maf_newMafBlock();
  if (mfa->lastLine != NULL) {
    // this is only invoked when the header is not followed by a blank line
    mafLine_t *ml = maf_newMafLineFromString(mfa->lastLine, mfa->lineNumber);
    if (ml->type == 's') {
      ++(thisBlock->numberOfSequences);
      if (thisBlock->sequenceFieldLength == 0) {
        thisBlock->sequenceFieldLength = maf_mafLine_getSequenceFieldLength(ml);
      }
    }
    ++(thisBlock->numberOfLines);
    thisBlock->headLine = ml;
    thisBlock->tailLine = ml;
    free(mfa->lastLine);
    mfa->lastLine = NULL;
  }
  int64_t n = kMaxStringLength;
  char *line = (char*) de_malloc(n);
  thisBlock->lineNumber = mfa->lineNumber;
  while(de_getline(&line, &n, mfa->mfp) != -1) {
    ++(mfa->lineNumber);
    if (maf_isBlankLine(line)) {
      if (thisBlock->headLine == NULL) {
        // this handles multiple blank lines in a row
        continue;
      } else {
        break;
      }
    }
    mafLine_t *ml = maf_newMafLineFromString(line, mfa->lineNumber);
    if (thisBlock->headLine == NULL) {
      thisBlock->headLine = ml;
      thisBlock->tailLine = ml;
    } else {
      thisBlock->tailLine->next = ml;
      thisBlock->tailLine = ml;
    }
    if (ml->type == 's') {
      ++(thisBlock->numberOfSequences);
      if (thisBlock->sequenceFieldLength == 0) {
        thisBlock->sequenceFieldLength = maf_mafLine_getSequenceFieldLength(ml);
      }
    }
    ++(thisBlock->numberOfLines);
  }
  free(line);
  line = NULL;
  return thisBlock;
}
mafBlock_t* maf_readBlock(mafFileApi_t *mfa) {
  // either returns a pointer to the next mafBlock in the maf file,
  // or a NULL pointer if the end of the file has been reached.
  if (mfa->lineNumber == 0) {
    // header
    mafBlock_t *header = maf_readBlockHeader(mfa);
    if (header->headLine != NULL) {
      return header;
    } else {
      maf_destroyMafBlockList(header);
      return NULL;
    }
  } else {
    // body
    mafBlock_t *mb = maf_readBlockBody(mfa);
    if (mb->headLine != NULL) {
      return mb;
    } else {
      maf_destroyMafBlockList(mb);
      return NULL;
    }
  }
}
mafBlock_t* maf_readAll(mafFileApi_t *mfa) {
  // read an entire mfa, creating a linked list of mafBlock_t, returning the head.
  mafBlock_t *head = maf_readBlock(mfa);
  mafBlock_t *mb = head;
  mafBlock_t *tmp = NULL;
  while ((tmp = maf_readBlock(mfa)) != NULL) {
    mb->next = tmp;
    mb = tmp;
  }
  return head;
}
void maf_writeAll(mafFileApi_t *mfa, mafBlock_t *mb) {
  // write an entire mfa, creating a linked list of mafBlock_t, returning the head.
  while (mb != NULL) {
    maf_writeBlock(mfa, mb);
    mb = mb->next;
  }
  fprintf(mfa->mfp, "\n");
  ++(mfa->lineNumber);
  fclose(mfa->mfp);
  mfa->mfp = NULL;
}
void maf_writeBlock(mafFileApi_t *mfa, mafBlock_t *mb) {
  mafLine_t *ml = mb->headLine;
  while (ml != NULL) {
    fprintf(mfa->mfp, "%s\n", ml->line);
    ++(mfa->lineNumber);
    ml = ml->next;
  }
  fprintf(mfa->mfp, "\n");
  ++(mfa->lineNumber);
}
void maf_mafBlock_appendToAlignmentBlock(mafBlock_t *m, char *s) {
  mafLine_t *ml = maf_mafBlock_getHeadLine(m);
  char *line = maf_mafLine_getLine(ml);
  assert(line[0] == 'a');
  char *newline = (char*) de_malloc(strlen(line) + strlen(s) + 1);
  newline[0] = '\0';
  strcat(newline, line);
  strcat(newline, s);
  free(ml->line);
  ml->line = newline;
}
void maf_mafBlock_printList(mafBlock_t *m) {
  while (m != NULL) {
    maf_mafBlock_print(m);
    m = maf_mafBlock_getNext(m);
  }
}
void maf_mafBlock_print(mafBlock_t *m) {
  // pretty print a mafBlock.
  if (m == NULL) {
    printf("..block NULL\n");
    return;
  }
  mafLine_t* ml = maf_mafBlock_getHeadLine(m);
  char *line = NULL;
  uint64_t maxName = 1, maxStart = 1, maxLen = 1, maxSource = 1;
  char fmtName[32] = "\0", fmtStart[32] = "\0", fmtLen[32] = "\0", fmtSource[32] = "\0", fmtLine[256] = "\0";
  while (ml != NULL) {
    line = maf_mafLine_getLine(ml);
    if (line == NULL) {
      break;
    }
    if (maf_mafLine_getType(ml) != 's') {
      ml = maf_mafLine_getNext(ml);
      continue;
    }
    if (maxName < strlen(maf_mafLine_getSpecies(ml))) {
      maxName = strlen(maf_mafLine_getSpecies(ml));
    }
    if (maxStart < maf_mafLine_getStart(ml)) {
      maxStart = maf_mafLine_getStart(ml);
    }
    if (maxLen < maf_mafLine_getLength(ml)) {
      maxLen = maf_mafLine_getLength(ml);
    }
    if (maxSource < maf_mafLine_getSourceLength(ml)) {
      maxSource = maf_mafLine_getSourceLength(ml);
    }
    ml = maf_mafLine_getNext(ml);
  }
  ml = maf_mafBlock_getHeadLine(m);
  if (ml != NULL) {
    sprintf(fmtName, " %%-%"PRIu64"s", maxName + 2);
    sprintf(fmtStart, " %%%d"PRIu64, (int)log10(maxStart) + 2);
    sprintf(fmtLen, " %%%d"PRIu64, (int)log10(maxLen) + 2);
    sprintf(fmtSource, " %%%d"PRIu64, (int)log10(maxSource) + 2);
    strcat(fmtLine, "s");
    strcat(fmtLine, fmtName);
    strcat(fmtLine, fmtStart);
    strcat(fmtLine, fmtLen);
    strcat(fmtLine, " %c");
    strcat(fmtLine, fmtSource);
    strcat(fmtLine, " %s\n");
  }
  while (ml != NULL) {
    line = maf_mafLine_getLine(ml);
    if (line == NULL) {
      // printf("..line NULL\n");
      break;
    }
    if (maf_mafLine_getType(ml) != 's') {
      printf("%s\n", line);
    } else {
      printf(fmtLine, maf_mafLine_getSpecies(ml), maf_mafLine_getStart(ml), maf_mafLine_getLength(ml),
             maf_mafLine_getStrand(ml), maf_mafLine_getSourceLength(ml), maf_mafLine_getSequence(ml));
    }
    ml = maf_mafLine_getNext(ml);
  }
  printf("\n");
}
static int intmax(int a, int b) {
  if (a > b) {
    return a;
  } else {
    return b;
  }
}
char* maf_mafLine_imputeLine(mafLine_t* ml) {
  // HEY!
  // IF you want to print a maf block, use maf_mafBlock_print(mafBlock_t *m). That is all.
  //
  // given a mafLine_t *, return a fresh char *
  // that represents a pretty print of the maf line.
  // i.e. it could be printed directly into a .maf file
  // uint32max = 4294967296 which is ~ 5 * 10^10
  char *s = (char*) de_malloc(2 + intmax(strlen(maf_mafLine_getSpecies(ml)), 15) +
                              15 + 15 + 3 + 15 + maf_mafLine_getSequenceFieldLength(ml) +
                              1 + 32); // extra 32 for possbile overflow from formatting
  sprintf(s, "s %-15s %10" PRIu64 " %10" PRIu64 " %c %10" PRIu64 " %s",
          maf_mafLine_getSpecies(ml),
          maf_mafLine_getStart(ml),
          maf_mafLine_getLength(ml),
          maf_mafLine_getStrand(ml),
          maf_mafLine_getSourceLength(ml),
          maf_mafLine_getSequence(ml));
  return s;
}
uint64_t countNonGaps(char *seq) {
  uint64_t n = strlen(seq);
  uint64_t m = 0;
  for (uint64_t i = 0; i < n; ++i) {
    if (seq[i] != '-') {
      ++m;
    }
  }
  return m;
}
void maf_mafBlock_flipStrand(mafBlock_t *mb) {
  // take a maf block and perform an in-place strand flip (including reverse complementing the
  // sequence, transforming the start coords) on all maf lines in the block.
  mafLine_t *ml = maf_mafBlock_getHeadLine(mb);
  while (ml != NULL) {
    if (maf_mafLine_getType(ml) != 's') {
      ml = maf_mafLine_getNext(ml);
      continue;
    }
    // rc sequence
    reverseComplementSequence(maf_mafLine_getSequence(ml), maf_mafBlock_getSequenceFieldLength(mb));
    // coordinate transform
    maf_mafLine_setStart(ml, maf_mafLine_getSourceLength(ml) -
                         (maf_mafLine_getStart(ml) + maf_mafLine_getLength(ml)));
    // strand flip
    if (maf_mafLine_getStrand(ml) == '+') {
      maf_mafLine_setStrand(ml, '-');
    } else {
      maf_mafLine_setStrand(ml, '+');
    }
    ml = maf_mafLine_getNext(ml);
  }
}
void reverseComplementSequence(char *s, size_t n) {
  // accepts upper and lower case, full iupac
  int c, i, j;
  for (i = 0, j = n - 1; i < j; ++i, j--) {
    c = s[i];
    s[i] = s[j];
    s[j] = c;
  }
  complementSequence(s, n);
}
void complementSequence(char *s, size_t n) {
  // accepts upper and lower case, full iupac
  for (unsigned i = 0; i < n; ++i)
    s[i] = complementChar(s[i]);
}
char complementChar(char c) {
  // accepts upper and lower case, full iupac
  bool wasUpper = false;
  char a = '\0';
  if (toupper(c) == c) {
    wasUpper = true;
  }
  switch (toupper(c)) {
  case 'A':
    a = 't';
    break;
  case 'C':
    a = 'g';
    break;
  case 'G':
    a = 'c';
    break;
  case 'T':
    a = 'a';
    break;
  case 'M':
    a = 'k';
    break;
  case 'R':
    a = 'y';
    break;
  case 'W':
    a = 'w';
    break;
  case 'S':
    a = 's';
    break;
  case 'Y':
    a = 'r';
    break;
  case 'K':
    a = 'm';
    break;
  case 'V':
    a = 'b';
    break;
  case 'H':
    a = 'd';
    break;
  case 'D':
    a = 'h';
    break;
  case 'B':
    a = 'v';
    break;
  case 'N':
  case '-':
  case 'X':
    a = c;
  break;
  default:
    fprintf(stderr, "Error, unanticipated character in DNA sequence: %c\n", c);
    exit(EXIT_FAILURE);
  }
  if (wasUpper) {
    a = toupper(a);
  }
  return a;
}
char *copySpeciesName(const char *s) {
  // return a copy of the string, minus chromosome / contig information
  // hg18.chr1 -> hg18
  unsigned n, l = 0;
  n = strlen(s);
  for (l = 0; l < n; ++l) {
    // find the first instance of a `.' character
    if (s[l] == '.') {
      break;
    }
  }
  char *copy = (char *) de_malloc(l + 1);
  strncpy(copy, s, l);
  copy[l] = '\0';
  return copy;
}
char *copyChromosomeName(const char *s) {
  // return a copy of the string, minus species information
  // hg18.chr1 -> chr1
  unsigned n, l = 0;
  n = strlen(s);
  for (l = 0; l < n; ++l) {
    // find the first instance of a `.' character
    if (s[l] == '.') {
      l += 1;
      break;
    }
  }
  if (l < n) {
    char *copy = (char *) de_malloc(n - l + 1);
    strncpy(copy, s + l, n - l);
    copy[n - l] = '\0';
    return copy;
  } else {
    char *copy = (char *) de_malloc(1);
    copy[0] = '\0';
    return copy;
  }
}
