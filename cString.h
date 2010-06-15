#ifndef CSTRING_H_
#define CSTRING_H_

#include "commonC.h"

#include <assert.h>
#include <ctype.h>
#include <string.h>
#include <stdint.h>

/*
 * Comparison function to sort strings alphabetically
 */
int cStr_compare(const void *a, const void *b);

/*
 * Comparison function to sort strings in descending order
 */
int cStr_compareDesc(const void *a, const void *b);

/*
 * In-place substitution to lower-case string
 */
void cStr_lowerCase(char *string);

/*
 * In-place substitution to upper-case string
 */
void cStr_upperCase(char *string);

/* 
 * Check if "string" starts with "query" and ignores case
 *   if "ignorecase" == 1
 */
int cStr_startsWith(char *string, char *query, int ignorecase);

int32_t cStr_getIntLength(int32_t n);

void cStr_reverse(char *s);

void cStr_itoa(int n, char *s);

void cStr_appendChar(char *s, char c);

char *cStr_getStringFromIntArray(int32_t *array, int32_t size, const char sep);

#endif /* CSTRING_H_ */
