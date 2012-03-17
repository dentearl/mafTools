#ifndef COMMON_H_
#define COMMON_H_
#include <stdint.h>

void verbose(char const *fmt, ...);
void debug(char const *fmt, ...);
void message(char const *type, char const *fmt, ...);

void* de_malloc(size_t n);
int32_t de_getline(char **s, int32_t *n, FILE *f);

#endif // COMMON_H_
