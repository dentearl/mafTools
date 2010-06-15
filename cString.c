#include "cString.h"

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

int32_t countDigits(int32_t n) {
	int32_t count = 0;
	do {
		count++;
	} while ((n /= 10) > 0);

	if (n < 0) {
		count++;
	}
	return count;
}

/* reverse:  reverse string s in place */
void reverse(char *s) {
	int i, j;
	char c;

	for (i = 0, j = strlen(s) - 1; i < j; i++, j--) {
		c = s[i];
		s[i] = s[j];
		s[j] = c;
	}
}

/* itoa:  convert n to characters in s */
void itoa(int n, char *s) {
	int i, sign;

	if ((sign = n) < 0)	/* record sign */
		n = -n;		/* make n positive */
	i = 0;
	do {			/* generate digits in reverse order */
		s[i++] = n % 10 + '0';	/* get next digit */
	} while ((n /= 10) > 0);/* delete it */
	if (sign < 0)
		s[i++] = '-';
	s[i] = '\0';
	reverse(s);
}

void append(char *s, char c) {
	int len = strlen(s);
	s[len] = c;
	s[len + 1] = '\0';
}

char *arrayToString(int32_t *array, int32_t size, const char sep) {
	int32_t i;
	int numChars = 0;
	char *string = NULL;
	char buffer[64];

	for (i = 0; i < size; i++) {
		numChars += countDigits(array[i]);
	}
	numChars += (size - 1);

	string = st_malloc(sizeof(char) * (numChars + 1));
	string[0] = '\0';

	itoa(array[0], buffer);
	strcat(string, buffer);
	for (i = 1; i < size; i++) {
		append(string, sep);
		itoa(array[i], buffer);
		strcat(string, buffer);
	}

	i = strlen(string);
	string[i+1] = '\0';

	return string;
}
