#include <stdlib.h>
#include <stdio.h>
#include "spmat.h"
#include <math.h>

/**
@file spmat.c
**Author:** Ofek Bransky & Gal Cohen
**Date:**  18.9.2020
**Summary:
 * This is the Sparse matrix C file, maintains the sparse matrix struct,  inner array implementation, and the errors which can be met during runtime.
 * Functions:
 * error - function that activates a runtime error
 * initk - initializes an array to store vertice ranks (with zero values)
 * add_row_array - adds a row of values to the sparse matrix array
 * hasNextARowArray - part of the row iterator, checks if theres a next value in the row
 * getARowIteratorArray - iterates on row values
 * mult_array - multiplication of a vector with the sparse matrix
 * free_array - frees all sparse matrix allocated memory
 * arrayShifting - calculates matrix shifting value
 * spmat_allocate_array - allocates an sparse matrix array struct
 * find_nnz - finds the non zero values in a file input
 * readArray - reads an graph from a file into a sparse matrix array
 * getNewNnz - gets the non zero values of the sparse matrix after a split of a group into two
 * splitArray - splits a sparse matrix array into two small ones (based on modularity split)
 * splitGraphArray - creates to new sparse matrix struct during group split
 * readGraphA - reads an graph from a file into a sparse matrix array
*/

/**
 * definition of CSR array to maintain the sparse matrices,
 * this is an inner implementation of the spmat struct declared in the header file
 * @param colind : pointer to an array containing column indexes of non zero values in sparse matrix
 * @param rowptr : pointer to an array containing row pointers in sparse matrix
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

/**
 * Method that receives a type of error, and prints the error event
 * @param errorCode : the type of error
 */
void error(int errorCode) {
    switch (errorCode) {
        case ALLOCERROR:
            printf("ERROR - memory allocation unsuccessful\n");
            break;
        case PIERROR:
            printf("Power Iteration can't converge\n");
            break;
        case READVALERROR:
            printf("ERROR - mismatch in data during reading input values (file cannot be found or is corrupt)\n");
            break;
        case ARGSERROR:
            printf("ERROR - there aren't 2 input arguments\n");
            break;
        case ZERODIV:
            printf("ERROR - dividing with value that equals to zero (no edges found!)\n");
            break;
        case FILECORR:
            printf("ERROR - input file not found or is corrupt\n");
            break;
        case FILEOUT:
            printf("ERROR - output file not found or is corrupt\n");
            break;
        default:
            printf("unexpected error\n");
    }
}

/**
 * initializes the rank to vertice array to zero values
 * @param A : the sparse matrix
 */
void initk(spmat *A) {
    int i, *k = A->k;
    for (i = 0; i < A->n; ++i) {
        *k = 0;
        k++;
    }
}

/**
 * Adds a row to the sparse matrix
 * @param A : the sparse matrix
 * @param row :the row to be added
 * @param i : the index of the added row
 * @param k : number of values in the row to be added
 */
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

/**
 * Part of a row value iterator, checks if there's a next value in a specific row
 * @param A : The sparse matrix to iterate upon
 * @param i : the row to check
 * @param ptr : the pointer to the current column value
 * @return : returns 1 if true, 0 if false
 */
int hasNextARowArray(spmat *A, int i, const int *ptr) {
    array *spArray = (array *) A->private;
    return ptr - (spArray->colind) < *(spArray->rowptr + i + 1) - 1;
}

/**
 * returns an iterator to itererate upon a specific row
 * @param A : the sparse matrix
 * @param i : the row to iterate upon
 * @return  : a pointer to the start of the row
 */
int *getARowIteratorArray(spmat *A, int i) {
    array *spArray = (array *) A->private;
    if (*(spArray->rowptr + i + 1) - *(spArray->rowptr + i) > 0) {
        return spArray->colind + *(spArray->rowptr + i);
    }
    return NULL;
}

/**
 * Multiplies a sparse matrix array (CSR) with a vector
 * @param A : the sparse matrix, stored in CSR format
 * @param vec : the vector to be multiplied
 * @param result : the result of the multiplication
 * @param groupSize : the size of the subgroup which is represented by the sparse matrix
 */
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

/**
 * Frees the allocated sparse matrix, including inner array implementation
 * @param A : the sparse matrix to be freed
 */
void free_array(struct _spmat *A) {
    array *sparray = (array *) A->private;
    free(sparray->rowptr);
    free(sparray->colind);
    free(A->k);
    free(A->private);
}

/**calculates matrix shifting value to get positive eigen values
 * @param A : the sparse matrix of the current subgroup
 * @param group : the current subgroup
 * @return the shifting value
 * */
double arrayShifting(spmat *A, int groupSize, const double *F) {
    double max = 0;
    register double sum;
    register int val;
    array *sparray = A->private;
    int *rowPtr = sparray->rowptr;
    int *cols;
    int colStart;
    int colEnd;
    int i;
    int j;
    int ki;
    int kj;
    register int M = A->M;
    int vertice1;
    double Fi;

    for (i = 0; i < groupSize; ++i) {
        sum = 0;
        vertice1 = i;
        ki = A->k[i];
        Fi = (double) F[i];
        colStart = *(rowPtr + i);
        cols = sparray->colind + colStart;
        colEnd = *(rowPtr + i + 1);
        for (j = 0; j < groupSize; ++j) {
            kj = A->k[j];
            if ((*(cols) == j) && (colStart < colEnd)) {
                val = 1;
                colStart++;
                if (colStart < colEnd) { cols++; }
            } else {
                val = 0;
            }
            if (vertice1 != j) {
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

/**
 * the function splits the sparse matrix into two sparse matrix, once a division is found, each representing a subgroup
 * @param currSp : the current sparse matrix to be split
 * @param s : the division vector
 * @param group : the original vertex group
 * @param g1Size : group 1's size
 * @param g2Size : group 2's size
 * @return an array with two pointers-> one to each new sparse matrix (after split into two different matrices)
 */
spmat **splitGraphArray(spmat *currSp, double *s, int *group, int g1Size, int g2Size);

/**
 * Allocates a sparse matrix array in CSR format
 * @param n : size of the graph (n x n vertices)
 * @param nnz : the number of non zero values in the sparse matrix
 * @return a pointer to the sparse matrix
 */
spmat *spmat_allocate_array(int n, int nnz) {
    spmat *sp;
    array *sparray = malloc(sizeof(array));
    if (sparray == NULL) {
        error(ALLOCERROR);
        exit(EXIT_FAILURE);
    }
    sparray->colind = malloc(nnz * sizeof(int));
    sparray->rowptr = malloc((n + 1) * sizeof(int));
    if ((sparray->rowptr == NULL) || (sparray->colind == NULL)) {
        error(ALLOCERROR);
        exit(EXIT_FAILURE);
    }
    sparray->lastindex = 0;
    sparray->lastRowPtr = 0;
    sparray->rowptr[0] = 0;
    sparray->nnz = nnz;
    sp = malloc(sizeof(spmat));
    if (sp == NULL) {
        error(ALLOCERROR);
        exit(EXIT_FAILURE);
    }
    sp->n = n;
    sp->add_row = add_row_array;
    sp->free = free_array;
    sp->mult = mult_array;
    sp->private = sparray;
    sp->splitGraph = splitGraphArray;
    sp->getARowIterator = getARowIteratorArray;
    sp->hasNextARow = hasNextARowArray;
    sp->M = 0;
    sp->k = malloc(sizeof(int) * n);
    if (sp->k == NULL) {
        error(ALLOCERROR);
        exit(EXIT_FAILURE);
    }
    initk(sp);
    sp->matShifting = arrayShifting;
    return sp;
}

/**
 * Counts the number of non zero values in the initial graph (to allocate a correct size of sparse matrix)
 * @param input : the input file
 * @return : the number NNZ
 */
int find_nnz(FILE *input) {
    int res;
    int size;
    fread(&size, sizeof(int), 1, input);
    fseek(input, 0L, SEEK_END);
    res = ftell(input) - (int) ((1 + size) * sizeof(int));
    fseek(input, 0, SEEK_SET);
    return res / 4;
}

/**
 * Reads the initial graph from a file input
 * @param input : the input file
 * @return a pointer to the struct
 */
spmat *readArray(FILE *input) {
    spmat *graph;
    int i, size, elem, *row;
    unsigned int n;
    int nnz;
    if (input == NULL) {
        error(FILECORR);
        exit(EXIT_FAILURE);
    }
    nnz = find_nnz(input);
    n = fread(&elem, sizeof(int), 1, input);
    if (n != 1) {
        error(READVALERROR);
        exit(EXIT_FAILURE);
    }
    size = elem;
    graph = spmat_allocate_array(size, nnz);
    row = malloc(sizeof(int) * size);
    if (row == NULL) {
        error(ALLOCERROR);
        exit(EXIT_FAILURE);
    }
    for (i = 0; i < size; ++i) {
        n = fread(&elem, sizeof(int), 1, input);
        if (n != 1) {
            error(READVALERROR);
            exit(EXIT_FAILURE);
        }
        n = fread(row, sizeof(int), elem, input);
        if (n != elem) {
            error(READVALERROR);
            exit(EXIT_FAILURE);
        }
        graph->add_row(graph, row, i, elem);
    }
    free(row);
    return graph;
}

/**
 *  this function is used to count the non zero values in the original sparse matrix, belonging to the divided subgroup,Once a division is found.
 * this method is used when splitting a sparse matrix, then allocating a correct new sparse matrix with the found size
 * @param sp : the original sparse matrix
 * @param s  : the division vector
 * @param group : the original group
 * @param groupSize : the original group size
 * @param newNnz : the new non zero value pointer to insert the found NNZ
 */
void getNewNnz(spmat *sp, const double *s, int *group, int groupSize, int *newNnz) {
    int *g2Nnz = &newNnz[1], i, counterNnz1 = 0, counterNnz2 = 0, t, *groupCopy;
    array *spArray = (array *) sp->private;
    int *rowPtr = spArray->rowptr, curr, j, valsInRow, *colIdx, flag, *counter;
    for (i = 0; i < groupSize; ++i) {
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

/**
 * this function copies the values of a sparse matrix belonging to a group into two subgroups (2 arrays), based on a division vector input
 * @param currSp : the original Sparse matrix array
 * @param g1Sp : the new sparse matrix for group 1
 * @param g2Sp : the new sparse matrix for group 2
 * @param s : the division vector
 * @param group : the original group
 * @param groupSize : the original group size
 */
void splitArray(spmat *currSp, spmat *g1Sp, spmat *g2Sp, double *s, int *group, int groupSize) {
    int i, j, t, valsCounter, flag;
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
        *(currArray->rowptr + currArray->lastRowPtr + 1) = *(currArray->rowptr + currArray->lastRowPtr) + valsCounter;
        currArray->lastindex += valsCounter;
        currArray->lastRowPtr++;
    }
}

/**
 * the function splits the sparse matrix into two sparse matrix, once a division is found, each representing a subgroup
 * @param currSp : the current sparse matrix to be split
 * @param s : the division vector
 * @param group : the original vertex group
 * @param g1Size : group 1's size
 * @param g2Size : group 2's size
 * @return an array with two pointers-> one to each new sparse matrix (after split into two different matrices)
 */
spmat **splitGraphArray(spmat *currSp, double *s, int *group, int g1Size, int g2Size) {
    int groupSize = g1Size + g2Size;
    spmat *g1Sp, *g2Sp, **newSpMats;
    int *newNnz = malloc(sizeof(int) * 2);
    newSpMats = malloc(sizeof(spmat *) * 2);
    if (newNnz == NULL || newSpMats == NULL) {
        error(ALLOCERROR);
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
    newSpMats[0] = g1Sp;
    newSpMats[1] = g2Sp;
    return newSpMats;
}
//TODO why do we need this?
/**
 * when called, reads the array from a file into a sparse matrix
 * @param input : a pointer to the input binary file
 * @return a pointer to the new sparse matrix
 */
spmat *readGraphA(FILE *input) {
    return readArray(input);
}


