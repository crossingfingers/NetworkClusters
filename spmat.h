#ifndef _SPMAT_H
#define _SPMAT_H

typedef struct _spIterator {
    int colIdx;
    int rowIdx;
    int kRow;
    int kCol;
    void *private;
} spIterator;

typedef struct _spmat {
    /* Matrix size (n*n) */
    int n;
    int M;

    /* Adds row i the matrix. Called before any other call,
     * exactly n times in order (i = 0 to n-1) */
    void (*add_row)(struct _spmat *A, int *row, int i, int k);

    /* Frees all resources used by A */
    void (*free)(struct _spmat *A);

    /* Multiplies matrix A by vector v, into result (result is pre-allocated) */
    void (*mult)(const struct _spmat *A, const double *v, double *result);

    void (*print_list)(struct _spmat *A);

    int (*hasNext)(struct _spmat *A);

    void (*initIterator)(struct _spmat *A);

    void (*iterNext)(struct _spmat *A);

    /* Private field for inner implementation.
     * Should not be read or modified externally */
    void *private;

    spIterator *iter;
} spmat;

/* Allocates a new linked-lists sparse matrix of size n */
spmat *spmat_allocate_list(int n);

/* Allocates a new arrays sparse matrix of size n with nnz non-zero elements */
spmat *spmat_allocate_array(int n, int nnz);

#endif
