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
#ifndef SHAREDMAF_H_
#define SHAREDMAF_H_
#include <assert.h>
#include <stdbool.h>
#include <stdint.h>
#include "common.h"
#include "CuTest.h"

typedef struct mafFileApi {
    // a mafFileApi struct provides an interface into a maf file.
    // Allows for easy reading of files in entirety or block by block via
    // functions
    uint32_t lineNumber; // last read line / wrote
    FILE *mfp; // maf file pointer
    char *filename; // filename of the maf
    char *lastLine; /* a temporary cache in case the header fails to have a blank 
                     * line before the first alignment block. 
                     */
} mafFileApi_t;
typedef struct mafLine {
    // a mafLine struct is a single line of a mafBlock
    char *line; // the entire line, unparsed
    uint32_t lineNumber; // line number in the maf file
    char type; // either a, s, i, q, e, h, f where h is header (an internal code)
    char *species; // species name
    uint32_t start;
    uint32_t length;
    char strand;
    uint32_t sourceLength;
    char *sequence; // sequence field
    struct mafLine *next;
} mafLine_t;
typedef struct mafBlock {
    // a mafBlock struct contains a maf block as a linked list
    // and itself can be part of a mafBlock linked list.
    mafLine_t *headLine;
    mafLine_t *tailLine;
    uint32_t lineNumber; // line number of start of the block (i.e. the `a' line)
    struct mafBlock *next;
} mafBlock_t;

mafFileApi_t* maf_newMfa(char *filename, char const *mode);
mafBlock_t* maf_readAll(mafFileApi_t *mfa);
mafBlock_t* maf_readBlock(mafFileApi_t *mfa);
mafBlock_t* maf_readHeader(mafFileApi_t *mfa);
void maf_writeAll(mafFileApi_t *mfa, mafBlock_t *mb);
void maf_writeBlock(mafFileApi_t *mfa, mafBlock_t *mb);
mafBlock_t* maf_newMafBlock(void);
mafLine_t* maf_newMafLine(void);
mafLine_t* maf_newMafLineFromString(char *s, uint32_t lineNumber);
unsigned maf_numberOfBlocks(mafBlock_t *b);
unsigned maf_numberOfSequences(mafBlock_t *b);
unsigned maf_numberOfSequencesMafLineList(mafLine_t *m);
void maf_destroyMafLineList(mafLine_t *ml);
void maf_destroyMafBlockList(mafBlock_t *mb);
void maf_destroyMfa(mafFileApi_t *mfa);
bool maf_blankLine(char *s);
void maf_checkForPrematureMafEnd(int status, char *line);
void maf_failBadFormat(uint32_t lineNumber, char *errorMessage);

mafLine_t* maf_newMafLine(void) {
    mafLine_t *ml = (mafLine_t *) de_malloc(sizeof(*ml));
    ml->line = NULL;
    ml->lineNumber = 0;
    ml->type = 0;
    ml->species = NULL;
    ml->start = 0;
    ml->length = 0;
    ml->strand = 0;
    ml->sourceLength = 0;
    ml->sequence = NULL;
    ml->next = NULL;
    return ml;
}
mafLine_t* maf_newMafLineFromString(char *s, uint32_t lineNumber) {
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
        return ml;
    }
    char *tkn = NULL;
    tkn = strtok(cline, " \t");
    if (tkn == NULL) {
        free(cline);
        maf_failBadFormat(lineNumber, "Unable to separate line on tabs and spaces at line definition field.");
    }
    tkn = strtok(NULL, " \t"); // name field
    if (tkn == NULL) {
        free(cline);
        maf_failBadFormat(lineNumber, "Unable to separate line on tabs and spaces at name field.");
    }
    char *species = (char *) de_malloc(strlen(tkn) + 1);
    strcpy(species, tkn);
    ml->species = species;
    tkn = strtok(NULL, " \t"); // start position
    if (tkn == NULL) {
        free(cline);
        maf_failBadFormat(lineNumber, "Unable to separate line on tabs and spaces at start position field.");
    }
    ml->start = strtoul(tkn, NULL, 10);
    tkn = strtok(NULL, " \t"); // length position
    if (tkn == NULL){
        free(cline);
        maf_failBadFormat(lineNumber, "Unable to separate line on tabs and spaces at length position field.");
    }
    ml->length = strtoul(tkn, NULL, 10);
    tkn = strtok(NULL, " \t"); // strand
    if (tkn == NULL) {
        free(cline);
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
        maf_failBadFormat(lineNumber, "Unable to separate line on tabs and spaces at source length field.");
    }
    ml->sourceLength = strtoul(tkn, NULL, 10);
    tkn = strtok(NULL, " \t"); // sequence field
    if (tkn == NULL) {
        free(cline);
        maf_failBadFormat(lineNumber, "Unable to separate line on tabs and spaces at sequence field.");
    }
    char *seq = (char *) de_malloc(strlen(tkn) + 1);
    strcpy(seq, tkn);
    ml->sequence = seq;
    free(cline);
    return ml;
}
mafBlock_t* maf_newMafBlock(void) {
    mafBlock_t *mb = (mafBlock_t *) de_malloc(sizeof(*mb));
    mb->next = NULL;
    mb->headLine = NULL;
    mb->tailLine = NULL;
    return mb;
}
mafFileApi_t* maf_newMfa(char *filename, char const *mode) {
    mafFileApi_t *mfa = (mafFileApi_t *) de_malloc(sizeof(*mfa));
    mfa->lineNumber = 0;
    mfa->lastLine = NULL;
    mfa->mfp = de_open(filename, mode);
    mfa->filename = de_strdup(filename);
    return mfa;
}
void maf_destroyMafLineList(mafLine_t *ml) {
    // walk down a mafLine_t following the ->next pointers, search and destroy
    mafLine_t *tmp = NULL;
    while(ml != NULL) {
        tmp = ml;
        ml = ml->next;
        free(tmp->line);
        free(tmp->species);
        free(tmp->sequence);
        free(tmp);
    }
}
void maf_destroyMafBlockList(mafBlock_t *mb) {
    mafBlock_t *tmp = NULL;
    while(mb != NULL) {
        tmp = mb;
        mb = mb->next;
        if (tmp->headLine != NULL)
            maf_destroyMafLineList(tmp->headLine);
        free(tmp);
    }
}
void maf_destroyMfa(mafFileApi_t *mfa) {
    if (mfa->mfp != NULL) {
        fclose(mfa->mfp);
        mfa->mfp = NULL;
    }
    free(mfa->lastLine);
    free(mfa->filename);
    free(mfa);
}
bool maf_blankLine(char *s) {
    // return true if line is only whitespaces
    size_t n = strlen(s);
    for (size_t i = 0; i < n; ++i) {
        if (!isspace(*(s + i))) {
            return false;
        }
    }
    return true;
}
unsigned maf_numberOfBlocks(mafBlock_t *b) {
    unsigned n = 0;
    while (b != NULL) {
        ++n;
        b = b->next;
    }
    return n;
}
unsigned maf_numberOfSequences(mafBlock_t *b) {
    // count the number of actual sequence lines in the mafBlock.
    return maf_numberOfSequencesMafLineList(b->headLine);
}
unsigned maf_numberOfSequencesMafLineList(mafLine_t *m) {
    // count the number of actual sequence lines in a mafLine_t list
    unsigned s = 0;
    while (m != NULL) {
        if (m->type == 's')
            ++s;
        m = m->next;
    }
    return s;
}
void maf_checkForPrematureMafEnd(int status, char *line) {
    if (status == -1) {
        free(line);
        fprintf(stderr, "Error, premature end to maf file\n");
        exit(EXIT_FAILURE);
    }
}
void maf_failBadFormat(uint32_t lineNumber, char *errorMessage) {
    fprintf(stderr, "The maf sequence at line %u is incorrectly formatted: %s\n", 
            lineNumber, errorMessage);
    exit(EXIT_FAILURE);
}
mafBlock_t* maf_readBlockHeader(mafFileApi_t *mfa) {
    extern const int kMaxStringLength;
    int32_t n = kMaxStringLength;
    char *line = (char*) de_malloc(n);
    mafBlock_t *header = maf_newMafBlock();
    int status = de_getline(&line, &n, mfa->mfp);
    ++(mfa->lineNumber);
    maf_checkForPrematureMafEnd(status, line);
    if (strncmp(line, "track", 5) == 0) {
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
        maf_checkForPrematureMafEnd(status, line);
    }
    if (strncmp(line, "##maf", 5) == 0) {
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
        maf_checkForPrematureMafEnd(status, line);
    }
    mafLine_t *thisMl = header->tailLine;
    while(line[0] != 'a' && !maf_blankLine(line)) {
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
        maf_checkForPrematureMafEnd(status, line);

    }
    if (line[0] == 'a') {
        char *copy = (char *) de_malloc(n + 1); // freed in destroy lines
        strcpy(copy, line);
        mfa->lastLine = copy;
    }
    free(line);
    return header;
}
mafBlock_t* maf_readBlockBody(mafFileApi_t *mfa) {
    extern const int kMaxStringLength;
    mafBlock_t *thisBlock = maf_newMafBlock();
    if (mfa->lastLine != NULL) {
        // this is only invoked when the header is not followed by a blank line
        mafLine_t *ml = maf_newMafLineFromString(mfa->lastLine, mfa->lineNumber);
        thisBlock->headLine = ml;
        thisBlock->tailLine = ml;
        free(mfa->lastLine);
        mfa->lastLine = NULL;
    }
    int32_t n = kMaxStringLength;
    char *line = (char*) de_malloc(n);
    while(de_getline(&line, &n, mfa->mfp) != -1) {
        ++(mfa->lineNumber);
        if (maf_blankLine(line)) {
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
    }
    free(line);
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
#endif // SHAREDMAF_H_
