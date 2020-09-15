//
// Created by gal21 on 05/09/2020.
//
#include "utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

/**
@file utils.c
**Author:** Ofek Bransky & Gal Cohen
**Date:**  18.9.2020
## This is the utility C file, contains various functions for different calculations in the program
*/

/**
 * Method that recieves a type of error, and prints the error event
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
            printf("ERROR - divide in zero\n");
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
 * Prints the vector
 * @param vec : the vector to be printed
 * @param n : the size of the vector
 */
void printVector(double *vec, int n) {
    int i, idx;
    for (i = 0; i < n; ++i) {
        printf("%f\t", vec[i]);
    }
    printf("\n");
}

/**
 * Prints a integer vector
 * @param vec : the vector to be printed
 * @param n : the size of the vector
 */
void printIntVector(int *vec, int n) {
    int i;
    for (i = 0; i < n; ++i) {
        printf("%d\t", vec[i]);
    }
    printf("\n");
}

/**
 * inserts random variables into an initialized vector
 * @param size : the size of the vector
 * @param vec : the initialized vector
 * @param groupSize : the number of values inserted, inserted into the first indexes
 * @param group : the subgroup array containing the vertices of the group
 */
void randomizeVec(int size, double *vec, int groupSize, int *group) {
    int i;
    for (i = 0; i < groupSize; i++) {
//        vec[i] = rand();
        vec[i] = group[i];
    }//TODO- Randomize vec before due date
}

/** gets to vectors and returns the sum of the vector with b0*shifting (shifting is the value from matrix shifting)*/
void vecSum(double *vec, const double *b0, double shifting, const int *group, int n) {
    int i;
    for (i = 0; i < n; ++i) {
        *vec += *b0++ * shifting;
        vec++;

    }
}

/**multiplies a scalar value with a vector*/
void scalarMult(double *vec, double x, const int *group, int n) {
    int i;
    for (i = 0; i < n; ++i) {
        *vec++ *= x;
    }
}

/**calculates the dot product of integer vector with double vector*/
double dotProd(const int *vec1, const double *vec2, const int *group, int n) {
    register double res = 0;
    int i;
    for (i = 0; i < n; ++i) {
        res += *vec1++ * *vec2++;
    }
    return res;
}

/**calculates the dot product of two double vectors*/
double dotDoubleProd(const double *vec1, const double *vec2, const int *group, int n) {
    register double res = 0;
    int i;
    for (i = 0; i < n; ++i) {

        res += *vec1++ * *vec2++;
    }
    return res;
}

/**normalizes a vector*/
void normalize(int size, double *vec, const int *group, int groupSize) {
    double res = dotDoubleProd(vec, vec, group, groupSize);
    res = sqrt(res);
    scalarMult(vec, 1 / res, group, groupSize);
}

/**the calculation of multiplying B matrix with a vector v by spliting B : to A (sparse matrix) and KiKj matrix (rank matrix)*/
void vecDecK(double *vec1, spmat *sp, int n, double dotM) {
    int i, *k = sp->k;
    for (i = 0; i < n; ++i) {
        *vec1 -= *k++ * dotM;
        vec1++;
    }
}

/**
 * multiplies a vector and the modularity matrix B
 * @param sp : the sparse matrix
 * @param vec : the vector to be multiplied by
 * @param res : the vector result of the multiplication
 * @param groupSize : the number of values inserted, inserted into the first indexes
 * @param group : the subgroup array containing the vertices of the group
 * @param debug
 * @return : the time taken to finish the method
 */
double multBv(spmat *sp, double *vec, const int *group, double *res, int groupSize) {
    double dot;
    double start, end;
    dot = dotProd(sp->k, vec, group, groupSize);
    start = clock();
    sp->mult(sp, vec, res, groupSize);
    end = clock();
    vecDecK(res, sp, groupSize, (double) dot / sp->M);
    return ((double) (end - start) / CLOCKS_PER_SEC);
}

/**
 * initializes a vector with all values identical
 * @param unitVec : the vector to be initialized
 * @param n : the size of the vector
 * @param group : the subgroup array containing the vertices of the group
 * @param val : the value to be inserted
 */
void initOneValVec(double *unitVec, int n, const int *group, int val) {
    int i;
    for (i = 0; i < n; ++i) {
        *unitVec = val;
        unitVec++;
    }
}

//TODO what does this do? for documentation
void vecDecFv(double *res, double *vecF, double *v, int n) {
    int i;
    for (i = 0; i < n; ++i) {
        *res -= *vecF++ * *v++;
        res++;
    }
}

/**
 * multiplies a vector and the modularity matrix B-roof, of a specific subgroup
 * @return
 * @param sp : the sparse matrix
 * @param vec : the vector to be multiplied by
 * @param res : the vector result of the multiplication
 * @param groupSize : the number of values inserted, inserted into the first indexes
 * @param group : the subgroup array containing the vertices of the group
 * @param vecF : an array containing the sum of values for each column
 * @param debug
 * @return : the time taken to finish the method
 */
double
multBRoof(spmat *sp, double *vec, const int *group, int groupSize, double *res, double *vecF) {
    int size = sp->n;
    double total;
    total = multBv(sp, vec, group, res, groupSize);
    vecDecFv(res, vecF, vec, size);
    return total;
}


/**
 * Power iteration method, multiplying a random vector with the matrix B, to find max eigenvector
 * @param sp : sparse matrix
 * @param b0  : initial random vector
 * @param shifting  : shifting value to shift the matrix & find max positive eigenvalue
 * @param result : the vector result of the power iteration
 * @param groupSize : the number of values inserted, inserted into the first indexes
 * @param group : the subgroup array containing the vertices of the group
 * @param vecF : an array containing the sum of values for each column
 * @param debug
 */
void powerIter(spmat *sp, double *b0, double shifting, int *group, int groupSize, double *result, double *vecF) {
    int flag = 1, i;
    int size = sp->n;
    int counter = 0;
    double MAX_ITERS = 0.5 * (size * size) + 5000 * size + 10000;
    double total = 0;
    while (flag == 1 && counter < MAX_ITERS) {
        flag = 0;
        total += multBRoof(sp, b0, group, groupSize, result, vecF);
        vecSum(result, b0, shifting, group, groupSize);
        normalize(size, result, group, groupSize);
        for (i = 0; i < groupSize; i++) {
            if (IS_POSITIVE(fabs(result[i] - b0[i])))
                flag = 1;
            b0[i] = result[i];
        }
        counter++;
    }
    if (counter == 10000) {
        error(PIERROR);
        exit(EXIT_FAILURE);
    }
}

/**
 * Calculates the eigenvalue of a vector
 * @param sp : sparse matrix
 * @param vec : the eigenvector
 * @param group : the subgroup array containing the vertices of the group
 * @param groupSize : the number of values inserted, inserted into the first indexes
 * @param vecF : an array containing the sum of values for each column
 * @return : the eigenvalue
 */
double eigenValue(spmat *sp, double *vec, const int *group, int groupSize, double *vecF) {
    int size = sp->n;
    double *tmp = malloc(sizeof(double) * size);
    if (tmp == NULL) {
        error(ALLOCERROR);
        exit(EXIT_FAILURE);
    }
    double res;
    multBRoof(sp, vec, group, groupSize, tmp, vecF);
    res = dotDoubleProd(tmp, vec, group, groupSize);
    free(tmp);
    return res;
}

/**
 * calculates the modularity of a division vector for a subgroup
 * @param sp : sparse matrix
 * @param vec : the eigenvector
 * @param group : the subgroup array containing the vertices of the group
 * @param groupSize : the number of values inserted, inserted into the first indexes
 * @param vecF : an array containing the sum of values for each column
 * @return the modularity of the division
 */
double modularityCalc(spmat *sp, double *vec, int *group, int groupSize, double *vecF) {
    double res = 0;
    int size = sp->n;
    double *tmp = malloc(sizeof(double) * size);
    if (tmp == NULL) {
        error(ALLOCERROR);
        exit(EXIT_FAILURE);
    }
    multBRoof(sp, vec, group, groupSize, tmp, vecF);
    res = dotDoubleProd(tmp, vec, group, groupSize);
    free(tmp);
    return res / 2;
}
