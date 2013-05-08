/* 
 * Copyright (C) 2009-2013 by 
 * Dent Earl (dearl@soe.ucsc.edu, dentearl@gmail.com)
 * Benedict Paten (benedict@soe.ucsc.edu, benedictpaten@gmail.com)
 * Mark Diekhans (markd@soe.ucsc.edu)
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


#include "cString.h"

/*
 * Comparison function to sort strings alphabetically
 */
int cStr_compare(const void *a, const void *b) {
        const char **ia = (const char **)a;
        const char **ib = (const char **)b;
        return strcmp(*ia, *ib);
}

/*
 * Comparison function to sort strings in descending order
 */
int cStr_compareDesc(const void *a, const void *b) {
        const char **ia = (const char **)a;
        const char **ib = (const char **)b;
        return -1 * strcmp(*ia, *ib);
}

/*
 * In-place substitution to lower-case string
 */
void cStr_lowerCase(char *string) {
        char *p;
        for (p=string; *p != '\0'; p++) {
                *p = tolower(*p);
        }
}

/*
 * In-place substitution to upper-case string
 */
void cStr_upperCase(char *string) {
        char *p;
        for (p=string; *p != '\0'; p++) {
                *p = toupper(*p);
        }
}

/*
 * Check if "string" starts with "query" and ignores case
 *   if "ignorecase" == 1
 */
int cStr_startsWith(char *string, char *query, int ignorecase) {
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

int64_t cStr_getIntLength(int64_t n) {
	int64_t count = 0;
	do {
		count++;
	} while ((n /= 10) > 0);

	if (n < 0) {
		count++;
	}
	return count;
}

/* reverse:  reverse string s in place */
void cStr_reverse(char *s) {
	int i, j;
	char c;

	for (i = 0, j = strlen(s) - 1; i < j; i++, j--) {
		c = s[i];
		s[i] = s[j];
		s[j] = c;
	}
}

/* itoa:  convert n to characters in s */
void cStr_itoa(int n, char *s) {
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
	cStr_reverse(s);
}

void cStr_appendChar(char *s, char c) {
	int len = strlen(s);
	s[len] = c;
	s[len + 1] = '\0';
}

char *cStr_getStringFromIntArray(int64_t *array, int64_t size, const char sep) {
	int64_t i;
	int numChars = 0;
	char *string = NULL;
	char buffer[64];

	for (i = 0; i < size; i++) {
		numChars += cStr_getIntLength(array[i]);
	}
	numChars += (size - 1);

	string = st_malloc(sizeof(char) * (numChars + 1));
	string[0] = '\0';

	cStr_itoa(array[0], buffer);
	strcat(string, buffer);
	for (i = 1; i < size; i++) {
		cStr_appendChar(string, sep);
		cStr_itoa(array[i], buffer);
		strcat(string, buffer);
	}

	i = strlen(string);
	string[i+1] = '\0';

	return string;
}
