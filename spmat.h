#ifndef _SPMAT_H
#define _SPMAT_H

typedef struct _spmat {
    /* Matrix size (n*n) */
    int n;
    int M;
    int *k;

    /* Adds row i the matrix. Called before any other call,
     * exactly n times in order (i = 0 to n-1) */
    void (*add_row)(struct _spmat *A, int *row, int i, int k);

    /* Frees all resources used by A */
    void (*free)(struct _spmat *A);

    /* Multiplies matrix A by vector v, into result (result is pre-allocated) */
    void (*mult)(const struct _spmat *A, const double *v, double *result, int group, const int *groupid);

    void (*printSprase)(struct _spmat *A);

    double (*matShifting)(struct _spmat *A, int group, const int *groupid, double *vecF);
    int (*isVal)(struct _spmat *A, int row, int col,int group, int *groupID);
    /* Private field for inner implementation.
     * Should not be read or modified externally */
    void *private;
} spmat;

spmat *readGraph(FILE *input);


#endif
