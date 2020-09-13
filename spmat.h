#ifndef _SPMAT_H
#define _SPMAT_H
#include <stdio.h>

/**
@file spmat.h
**Author:** Ofek Bransky & Gal Cohen
**Date:**  18.9.2020
## This is the Sparse matrix H file, maintains the sparse matrix structure
*/



struct  _networks;

/**Sparse matrix structure with various methods, maintained by a CSR array
 * @param n : size of matrix ( n x n)
 * @param M : total rank of full graph
 * @param *k : array containing ranks of all vertices in the sparse matrix
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
    void (*splitGraph)(struct _networks *graphs, int groupIdx, int newGroupIdx, double *s, int *group, int groupSize, int g1Size,
                       int g2Size);
    /*gets value of sparse matrix in coordinates (i,j) */
    int (*findAij)(struct _spmat *sp, int i, int j);

    /* Private field for inner implementation.
     * Should not be read or modified externally, in our case a pointer to the array struct */
    void *private;
} spmat;

/** a struct that contains all the sparse matrices in the program
 * @param **A : a array of sparse matrices
 * @paran n : number of vertices in the graph
 */
typedef struct _networks{
    spmat **A;
    int n;
    void (*free)(struct _networks *graphs, int numOfGroups);
}networks;

networks *readGraph(FILE *input);

#endif
