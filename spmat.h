#ifndef _SPMAT_H
#define _SPMAT_H

#include <stdio.h>


#define IS_POSITIVE(X) ((X) > 0.00001)
#define PIERROR 2
#define ALLOCERROR 1
#define READVALERROR 3
#define ARGSERROR 4
#define ZERODIV 5
#define FILECORR 6
#define FILEOUT 7

/**
@file spmat.h
**Author:** Ofek Bransky & Gal Cohen
**Date:**  18.9.2020
## This is the Sparse matrix H file, maintains the sparse matrix structure and networks structure
*/


/**Sparse matrix structure with various methods, maintained by a CSR array
 * @param n : size of matrix ( n x n)
 * @param M : total rank of full graph
 * @param k : an array containing ranks of all vertices in the sparse matrix
 */
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
    void (*mult)(const struct _spmat *A, const double *vec, double *result, int groupSize);

    /*prints matrix*/
    void (*printSparse)(struct _spmat *A);

    /*calculates matrix shifting value to get positive eigen values*/
    double (*matShifting)(struct _spmat *A, const int *group, int groupSize, const double *F);

    /*method to split sparse matrix into two subgroups (each sparse matrix contains only vertices in subgroup)*/
    struct _spmat **(*splitGraph)(struct _spmat *currSp, double *s, int *group, int g1Size, int g2Size);

    /*iterates over a specific row in the sparse matrix, returns 0 if theres a next value*/
    int (*hasNextARow)(struct _spmat *A, int i, const int *ptr);

    /*gets the next value in the row from the sparse matrix*/
    int *(*getARowIterator)(struct _spmat *A, int i);

    /* Private field for inner implementation.
     * Should not be read or modified externally, in our case a pointer to the array struct */
    void *private;
} spmat;

spmat *readGraphA(FILE *input);

#endif
