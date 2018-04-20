#ifndef _UTIL_H_

#define _UTIL_H_

#define min(a,b) ((a)>(b)?(b):(a))

#include <stdio.h>
#include <stdarg.h>

void exit_error(FILE *f_log, const char *fmt, ...);

void *malloc_or_die(size_t size, const char *func, const char *desc);

#endif
