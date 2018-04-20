#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <errno.h>
#include <string.h>
#include <mpi.h>

void
exit_error(FILE *f_log, const char *fmt, ...)
{
    int init_flag;
    va_list argptr;
    /* print formatted error message */
    va_start(argptr, fmt);
    vfprintf(f_log, fmt, argptr);
    va_end(argptr);
    /* exit, either via exit() or MPI_Abort() */
    MPI_Initialized(&init_flag);
    if (init_flag)
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    else
        exit(EXIT_FAILURE);
}

void *
malloc_or_die(size_t size, const char *func, const char *desc)
{
    void *p = malloc(size);
    if (p == NULL)
        exit_error(stdout, "Error [%s]: cannot allocate %li bytes for %s in %s"
                   " - %s\n", __func__, size, desc, func, strerror(errno));
    return p;
}
