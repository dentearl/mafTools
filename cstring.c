#include "cstring.h"

/*
 * Comparison function to sort strings alphabetically
 */
int cstring_cmp(const void *a, const void *b) {
        const char **ia = (const char **)a;
        const char **ib = (const char **)b;
        return strcmp(*ia, *ib);
}

/*
 * Comparison function to sort strings in descending order
 */
int cstring_cmp_desc(const void *a, const void *b) {
        const char **ia = (const char **)a;
        const char **ib = (const char **)b;
        return -1 * strcmp(*ia, *ib);
}

/*
 * In-place substitution to lower-case string
 */
void lowerCase(char *string) {
        char *p;
        for (p=string; *p != '\0'; p++) {
                *p = tolower(*p);
        }
}

/*
 * In-place substitution to upper-case string
 */
void upperCase(char *string) {
        char *p;
        for (p=string; *p != '\0'; p++) {
                *p = toupper(*p);
        }
}

