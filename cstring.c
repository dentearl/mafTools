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

/*
 * Check if "string" starts with "query" and ignores case
 *   if "ignorecase" == 1
 */
int startswith(char *string, char *query, int ignorecase) {
	assert(strlen(query) <= strlen(string));
	assert(strlen(string) > 0);
	assert(strlen(query) > 0);

	int i = 0;
	while(1) {
		if (query[i] == '\0') {
			return 1;
		}
		if (ignorecase) {
			if (tolower(string[i]) != tolower(query[i])) {
				return 0;
			}
		} else {
			if (string[i] != query[i]) {
				return 0;
			}
		}

		i++;
	}
}

