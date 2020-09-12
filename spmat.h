#ifndef _SPMAT_H
#define _SPMAT_H
#include <stdio.h>


struct  _networks;
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
    void (*mult)(const struct _spmat *A, const double *v, double *result, const int *group, int groupSize,const int *groupToVertice);

    void (*printSprase)(struct _spmat *A);

    double (*matShifting)(struct _spmat *A, const int *group, int groupSize, const int *vertexToGroup, int groupIdx, double *vecF);

    void (*splitGraph)(struct _networks *graphs, int groupIdx, int newGroupIdx, double *s, int *group, int groupSize, int g1Size,
                       int g2Size);

    /* Private field for inner implementation.
     * Should not be read or modified externally */
    void *private;
} spmat;

typedef struct _networks{
    spmat **A;
    int n;
    void (*free)(struct _networks *graphs, int numOfGroups);
}networks;

networks *readGraph(FILE *input, int type);

#endif
