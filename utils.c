//
// Created by gal21 on 05/09/2020.
//
#include "utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "spmat.h"
#include <time.h>
/**
@file utils.c
**Author:** Ofek Bransky & Gal Cohen
**Date:**  18.9.2020
## This is the utility C file, contains various functions for different calculations in the program
*/

/**
 * Method that recieves a type of error, and prints the error event
 * @param errorCode
 */
void error(int errorCode) {
    switch (errorCode) {
        case ALLOCERROR:
            printf("ERROR - memory allocation unsuccessful\n");
            break;
        case PIERROR:
            printf("PI can't converge\n");
            break;
        case READVALERROR:
            printf("ERROR - mismatch reading value\n");
            break;
        case ARGSERROR:
            printf("ERROR - there is not 2 arguments\n");
            break;
        case ZERODIV:
            printf("ERROR - divide in zero\n");
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
//        vec[i] = i;
    }
//    for (; i < size; ++i) {
//        vec[i] = 0;
//    }
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

//TODO remove debug
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
double multBv(spmat *sp, double *vec, const int *group, double *res, int groupSize, int debug) {
    double dot;
    int size = sp->n;
//    double *res1 = malloc(size * sizeof(double));
//    if (res1 == NULL) {
//        printf("ERROR - memory allocation unsuccessful");
//        exit(EXIT_FAILURE);
//    }
    double start, end;
    dot = dotProd(sp->k, vec, group, groupSize);
//    if (debug == 1){
//        printf("dot is %f\n", dot);
//    }
    if (sp->M == 0) {
        error(ZERODIV);
        exit(EXIT_FAILURE);
    }
//    copyVec(sp->k, res1, group, groupSize);
//    scalarMult(res1, (double) dot / sp->M, group, groupSize);
//    if (debug == 1){
//        printf("res1 is : ");
//        printVector(res1, size);
//    }
    start = clock();
    sp->mult(sp, vec, res, groupSize);
    end = clock();
//    if (debug == 2) {
//        printf("dot is %f\n", dot);
//        printf("res after mult A is : ");
//        printVector(res, groupSize, group);
//        printf("vec is : ");
//        printVector(vec, groupSize,group);
//        printf("groupsize is %d\n and group is : ", groupSize);
//        printIntVector(group, groupSize);
//    }
//    vecDec(res, res1, group, groupSize);
//    free(res1);
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
void vecDecFv(double *res, double* vecF, double *v, int n){
    int i;
    for(i = 0; i < n; ++i){
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
//TODO remove debug
double
multBRoof(spmat *sp, double *vec, const int *group, int groupSize, double *res, double *vecF, int debug) {
    int size = sp->n;
//    double *unitVec = malloc(size * sizeof(double));
//    double *vecFCopy = malloc(size * sizeof(double));
    double total, start, end;
//    if (vecFCopy == NULL) {
//        printf("ERROR - memory allocation unsuccessful");
//        exit(EXIT_FAILURE);
//    }
    start = clock();
//    initOneValVec(unitVec, groupSize, group, 1);
//    multBv(sp, unitVec, group, vecF, groupSize, 0);
//    copyDoubleVec(vecF, vecFCopy, group, groupSize);
//    vecMult(vecFCopy, vec, group, groupSize);
    total = multBv(sp, vec, group, res, groupSize, debug);
//    if (debug == 1) {
//        printVector(vec, size, group);
//            printVector(result, size, group);
//        printf("PI multBv total time took %f\n", total);
//    }
//    vecDec(res, vecFCopy, group, groupSize);
//    free(vecFCopy);
    vecDecFv(res, vecF, vec, size);
    debug = 0;
//    free(unitVec);
    end = clock();
//    return ((double)(end-start)/ CLOCKS_PER_SEC);
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
void powerIter(spmat *sp, double *b0, double shifting, int *group, int groupSize, double *result, double *vecF, int debug) {
    int flag = 1, i, idx;
    int size = sp->n;
    int counter = 0;
    double MAX_ITERS = 0.5 * (size * size) + 5000 * size + 10000;
    double total = 0;
    while (flag == 1 && counter < MAX_ITERS) {
        flag = 0;
        total += multBRoof(sp, b0, group, groupSize, result, vecF, debug);
        vecSum(result, b0, shifting, group, groupSize);
        normalize(size, result, group, groupSize);
        for (i = 0; i < groupSize; i++) {
            idx = group[i];
            if (IS_POSITIVE(fabs(result[i] - b0[i])))
                flag = 1;
            b0[i] = result[i];
        }
        counter++;
        debug = 0;
    }
    if (counter == 10000) {
        error(PIERROR);
        exit(EXIT_FAILURE);
    }
//    printf("took %d iterations\n", counter);
//    printf("mult A in PI took %f\n", total);
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
    multBRoof(sp, vec, group, groupSize, tmp, vecF, 0);
    res = dotDoubleProd(tmp, vec, group, groupSize);
    free(tmp);
    //   printf("eigen value is %f\n", res);
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
    multBRoof(sp, vec, group, groupSize, tmp, vecF, 0);
    res = dotDoubleProd(tmp, vec, group, groupSize);
    free(tmp);
    return res / 2;
}
