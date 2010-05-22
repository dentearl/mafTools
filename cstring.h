#ifndef CSTRING_H_
#define CSTRING_H_

#include <assert.h>
#include <ctype.h>
#include <string.h>

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

#endif /* CSTRING_H_ */
