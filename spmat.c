#include <stdlib.h>
#include <stdio.h>
#include "spmat.h"
#include <math.h>
#include <time.h>
/**
@file spmat.c
**Author:** Ofek Bransky & Gal Cohen
**Date:**  18.9.2020
## This is the Sparse matrix C file, maintains the sparse matrix and networks structures with methods
*/

/**
 * initializes the rank to vertice array to zero values
 * @param A : the sparse matrix
 */
void initk(spmat *A) {
    int i;
    for (i = 0; i < A->n; ++i) {
        A->k[i] = 0;
    }
}
/**
 * definition of CSR array to maintain the sparse matrices
 * @param *colind : pointer to column indexes of non zero values in sparse matrix
 * @param *rowptr : pointer to row pointers in sparse matrix
 * @param lastindex : used to add values to array when reading graph from input
 * @param lastRowPtr : used to add values to array when reading graph from input
 * @param nnz : the number of non zero values in the sparse matrix
 */
typedef struct _array {
    int *colind;
    int *rowptr;
    int lastindex;
    int lastRowPtr;
    int nnz;
} array;

//TODO remove values array, no need....
void add_row_array(struct _spmat *A, int *row, int i, int k) {
    array *sparray = (array *) A->private;
    int *rowInput = row;
    int ci;
    A->M += k;
    A->k[i] = k;
    sparray->rowptr[sparray->lastRowPtr + 1] = sparray->rowptr[sparray->lastRowPtr] + k;
    sparray->lastRowPtr++;
/*updates values array*/
    for (ci = 0; ci < k; ci++) {
        sparray->colind[sparray->lastindex] = *rowInput;
        rowInput++;
        sparray->lastindex++;
    }
}

int findAijArray(spmat *sp, int i, int j) {
    array *sparray = (array *) sp->private;
    int *rowPtr = sparray->rowptr;
    register int *cols;
    int colStart;
    register int size;
    register int counter = 0;
    int colEnd;
    colStart = *(rowPtr + i);
    cols = sparray->colind + colStart;
    colEnd = *(rowPtr + i + 1);
    size = colEnd - colStart;
    while (counter < size && *cols <= j) {
        if (*cols == j) {
            return 1;
        }
        counter++;
        cols++;
    }
    return 0;
}

void mult_array(const struct _spmat *A, const double *vec, double *result, int groupSize) {

    array *sparray = (array *) A->private;
    int *rowPtr = sparray->rowptr;
    register int *cols;
    int i;
    register int colStart;
    register int colEnd;
    register double res;
    for (i = 0; i < groupSize; ++i) {
        res = 0;
        colStart = *(rowPtr + i);
        cols = sparray->colind + colStart;
        colEnd = *(rowPtr + i + 1);
        while (colStart < colEnd) {
            res += vec[*(cols)];
            colStart++;
            if (colStart < colEnd) { cols++; }
        }
        result[i] = res;
    }
}

void free_array(struct _spmat *A) {
    array *sparray = (array *) A->private;
    free(sparray->rowptr);
    free(sparray->colind);
    free(A->k);
    free(A->private);
}

void printMatrix(spmat *A) {

    int i;
    array *sparr = (array *) A->private;


    printf("\ncollarr:\n");
    for (i = 0; i < sparr->nnz; i++) {
        printf("%d", sparr->colind[i]);
    }
    printf("\nrowptrarr:\n");
    for (i = 0; i <= A->n; i++) {
        printf("%d", sparr->rowptr[i]);
    }

}

void print_array(struct _spmat *A) {
    int i;
    int j;
    array *sparray = (array *) A->private;
    int *colindarr = (int *) sparray->colind;
    int counter = 0;
    for (i = 0; i < A->n; ++i) {
        printf("%d - \t", i);

        for (j = 0; j < A->k[i]; ++j) {
            printf("%d\t", colindarr[counter]);
            counter++;
        }
        printf("\n");
    }
    printMatrix(A);
    printf("\n");
}

double arrayShifting(spmat *A, const int *group, int groupSize, const double *F) {
    double max = 0;
    double sum;
    int val;
    array *sparray = A->private;
    int *rowPtr = sparray->rowptr;
    int *cols;
    int colStart;
    int colEnd;
    int i;
    int j;
    int ki;
    int kj;
    int M = A->M;
    int vertice1;
    double Fi;

    for (i = 0; i < groupSize; ++i) {
        sum = 0;
        vertice1 = group[i];
        ki = A->k[i];
        Fi = (double) F[i];
        colStart = *(rowPtr + i);
        cols = sparray->colind + colStart;
        colEnd = *(rowPtr + i + 1);
        for (j = 0; j < groupSize; ++j) {
            kj = A->k[j];
            if ((*(cols) == group[j]) && (colStart < colEnd)) {
                val = 1;
                colStart++;
                if (colStart < colEnd) { cols++; }
            } else {
                val = 0;
            }
            if (vertice1 != group[j]) {
                sum += fabs((double) val - ((double) (ki * kj) / M));
            } else {
                sum += fabs(
                        (double) val - ((double) (ki * kj) / M) - Fi);
            }
        }
        max = (max >= sum) ? max : sum;
    }
    return max;
}

void splitGraphArray(networks *graphs, int groupIdx, int newGroupIdx, double *s, int *group, int groupSize, int g1Size,
                     int g2Size);

spmat *spmat_allocate_array(int n, int nnz) {
    spmat *sp;
    array *sparray = malloc(sizeof(array));
    if (sparray == NULL) {
        printf("ERROR - memory allocation unsuccessful");
        exit(EXIT_FAILURE);
    }
    sparray->colind = malloc(nnz*sizeof(int));
    sparray->rowptr = malloc((n + 1)* sizeof(int));
    if ((sparray->rowptr==NULL)||(sparray->colind == NULL)) {
        printf("ERROR - memory allocation unsuccessful");
        exit(EXIT_FAILURE);
    }
    sparray->lastindex = 0;
    sparray->lastRowPtr = 0;
    sparray->rowptr[0] = 0;
    sparray->nnz = nnz;
    sp = malloc(sizeof(spmat));
    if (sp==NULL) {
        printf("ERROR - memory allocation unsuccessful");
        exit(EXIT_FAILURE);
    }
    sp->n = n;
    sp->add_row = add_row_array;
    sp->free = free_array;
    sp->mult = mult_array;
    sp->private = sparray;
    sp->printSparse = print_array;
    sp->splitGraph = splitGraphArray;
    sp->M = 0;
    sp->findAij = findAijArray;
    sp->k = malloc(sizeof(int) * n);
    if (sp->k==NULL) {
        printf("ERROR - memory allocation unsuccessful");
        exit(EXIT_FAILURE);
    }
    initk(sp);
    sp->matShifting = arrayShifting;
    return sp;
}

int find_nnz(FILE *input) {
    int res;
    int size;
    fread(&size, sizeof(int), 1, input);
    fseek(input, 0L, SEEK_END);
//    printf("size of file is %ld\n", ftell(input));
//    printf("size is %d\n", size);
    res = ftell(input) - (int) ((1 + size) * sizeof(int));
//    printf("res is %d\n", res);
    fseek(input, 0, SEEK_SET);
    return res / 4;
}

void freeNetworks(networks *graphs, int numOfGroups) {
    int i;
    spmat *sp, **mats = graphs->A;
    //TODO change 1 to numOfGroups!!
    for (i = 0; i < numOfGroups; ++i) {
        sp = *mats++;
        sp->free(sp);
        free(sp);
    }
    free(graphs->A);
}

networks *allocateNetworks(int n) {
    networks *graphs = malloc(sizeof(networks));
    if (graphs == NULL) {
        printf("ERROR - memory allocation unsuccessful");
        exit(EXIT_FAILURE);
    }
    graphs->n = n;
    graphs->A = malloc(sizeof(spmat *) * n);
    if (graphs->A == NULL) {
        printf("ERROR - memory allocation unsuccessful");
        exit(EXIT_FAILURE);
    }
    graphs->free = freeNetworks;
    return graphs;
}

networks *readArray(FILE *input) {
    spmat *graph;
    int i, size, elem, *row, j;
    networks *graphs;
    unsigned int n;
    int nnz;
    nnz = find_nnz(input);
//    printf("nnz is %d\n", nnz);
    n = fread(&elem, sizeof(int), 1, input);
    if (n != 1) {
        printf("ERROR - mismatch reading value");
        exit(EXIT_FAILURE);
    }
    size = elem;
//    printf("n: %d, nnz: %d\n", elem, nnz);
    graph = spmat_allocate_array(size, nnz);
    graphs = allocateNetworks(size);
    row = malloc(sizeof(int) * size);
    if (row == NULL) {
        printf("ERROR - memory allocation unsuccessful");
        exit(EXIT_FAILURE);
    }
    for (i = 0; i < size; ++i) {
        n = fread(&elem, sizeof(int), 1, input);
        if (n != 1) {
            printf("ERROR - mismatch reading value");
            exit(EXIT_FAILURE);
        }
        n = fread(row, sizeof(int), elem, input);
        if (n != elem) {
            printf("ERROR - mismatch reading value");
            exit(EXIT_FAILURE);
        }
        graph->add_row(graph, row, i, elem);
    }
    free(row);
    graphs->A[0] = graph;
//    graph->printSprase(graph);
    return graphs;
}


networks *readGraph(FILE *input) {
    return readArray(input);
}


void getNewNnz(spmat *sp, double *s, int *group, int groupSize, int *newNnz) {
    int *g2Nnz = &newNnz[1], i, counterNnz1 = 0, counterNnz2 = 0, t, *groupCopy = group;
    array *spArray = (array *) sp->private;
    int *rowPtr = spArray->rowptr, curr, j, valsInRow, idx, *colIdx, flag, *counter;
    for (i = 0; i < groupSize; ++i) {
        idx = group[i];
        if (s[i] == 1) {
            flag = 1;
            counter = &counterNnz1;
        } else {
            flag = -1;
            counter = &counterNnz2;
        }
        curr = *(rowPtr + i);
        valsInRow = *(rowPtr + i + 1) - curr;
        colIdx = spArray->colind + curr;
        j = 0;
        groupCopy = group;
        t = 0;
        while (j < valsInRow && t < groupSize) {
            if (*colIdx == t) {
                if (s[t] == flag) {
                    (*counter)++;
                }
                j++;
                t++;
                groupCopy++;
                colIdx++;
            } else if (*colIdx < t) {
                j++;
                colIdx++;
            } else {
                t++;
                groupCopy++;
            }
        }
    }
    *newNnz = counterNnz1;
    *g2Nnz = counterNnz2;
}

void splitArray(spmat *currSp, spmat *g1Sp, spmat *g2Sp, double *s, int *group, int groupSize) {
    int i, j, t, valsCounter, idxI, flag;
    int *groupCopy;
    array *currSpArray = (array *) currSp->private, *g1SpArray = (array *) g1Sp->private, *g2SpArray = (array *) g2Sp->private;
    array *currArray;
    int *rowPtr = currSpArray->rowptr, curr, valsInRow, *oldColIdx, *currColIdx, currGroupCounter;
    int *g1K = g1Sp->k, *g2K = g2Sp->k;
    double *sI = s, *sT;
    *(g1SpArray->rowptr) = 0;
    *(g2SpArray->rowptr) = 0;
    g1SpArray->lastindex = 0;
    g2SpArray->lastindex = 0;
    for (i = 0; i < groupSize; ++i) {
        idxI = group[i];
        if (*sI++ == 1) {
            currArray = g1SpArray;
            flag = 1;
            *g1K = currSp->k[i];
            g1K++;
        } else {
            currArray = g2SpArray;
            flag = -1;
            *g2K = currSp->k[i];
            g2K++;
        }
        curr = *(rowPtr + i);
        valsInRow = *(rowPtr + i + 1) - curr;
        j = 0;
        t = 0;
        oldColIdx = currSpArray->colind + curr;
        groupCopy = group;
        currColIdx = currArray->colind + currArray->lastindex;
        valsCounter = 0;
        currGroupCounter = 0;
        sT = s;
        while (j < valsInRow && t < groupSize) {
            if (*oldColIdx == t) {
                if (*sT == flag) {
//                    *currColIdx = *groupCopy;
                    *currColIdx = currGroupCounter;
                    currColIdx++;
                    valsCounter++;
                    currGroupCounter++;
                }
                j++;
                t++;
                sT++;
                groupCopy++;
                oldColIdx++;
            } else if (*oldColIdx < t) {
                j++;
                oldColIdx++;
            } else {
                if (*sT == flag)
                    currGroupCounter++;
                t++;
                sT++;
                groupCopy++;
            }
        }
        int test = *(currArray->rowptr + currArray->lastRowPtr);
        *(currArray->rowptr + currArray->lastRowPtr + 1) = *(currArray->rowptr + currArray->lastRowPtr) + valsCounter;
        currArray->lastindex += valsCounter;
        currArray->lastRowPtr++;
    }
}

void
splitGraphArray(networks *graphs, int groupIdx, int newGroupIdx, double *s, int *group, int groupSize, int g1Size,
                int g2Size) {
    spmat **Amats = graphs->A;
    spmat *currSp = Amats[groupIdx], *g1Sp, *g2Sp;
    int *newNnz = malloc(sizeof(int) * 2), size = graphs->n;
    if (newNnz == NULL) {
        printf("ERROR - memory allocation unsuccessful");
        exit(EXIT_FAILURE);
    }
    getNewNnz(currSp, s, group, groupSize, newNnz);
    g1Sp = spmat_allocate_array(g1Size, newNnz[0]);
    g2Sp = spmat_allocate_array(g2Size, newNnz[1]);
    splitArray(currSp, g1Sp, g2Sp, s, group, groupSize);
    g1Sp->M = currSp->M;
    g2Sp->M = currSp->M;
    currSp->free(currSp);
    free(currSp);
    free(newNnz);
    graphs->A[groupIdx] = g1Sp;
    graphs->A[newGroupIdx] = g2Sp;
}


