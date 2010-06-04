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
int cstring_cmp(const void *a, const void *b);

/*
 * Comparison function to sort strings in descending order
 */
int cstring_cmp_desc(const void *a, const void *b);

/*
 * In-place substitution to lower-case string
 */
void lowerCase(char *string);

/*
 * In-place substitution to upper-case string
 */
void upperCase(char *string);

/* 
 * Check if "string" starts with "query" and ignores case
 *   if "ignorecase" == 1
 */
int startswith(char *string, char *query, int ignorecase);

int32_t countDigits(int32_t n);

void reverse(char *s);

void itoa(int n, char *s);

void append(char *s, char c);

char *arrayToString(int32_t *array, int32_t size, const char sep);

#endif /* CSTRING_H_ */
