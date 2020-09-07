//
// Created by gal21 on 05/09/2020.
//
#include "utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "spmat.h"
#include <assert.h>
/*get a vector and the size of it and init any value to a random value*/
void randomizeVec(int size, double *vec) {
    int i;
    assert(vec != NULL);
    for (i = 0; i < size; i++) {
        vec[i] = rand();
    }
}

/* used for the F vector as a Matrix to multiply it by the v vector*/
void vecMult(double *vec1, const double *vec2, const int *group, int size) {
    int i, idx;
    for (i = 0; i < size; ++i) {
        idx = group[i];
        vec1[idx] = vec1[idx] * vec2[idx];
    }
}

/* gets to vectors and return the sum of vec with b0*shifting (shifting is the value from matrix shifting*/
void vecSum(double *vec, const double *b0, double shifting, const int *group, int n) {
    int i, idx;
    for (i = 0; i < n; ++i) {
        idx = group[i];
        vec[idx] += b0[idx] * shifting;
    }
}

/*decrease vec1 by vec2*/
void vecDec(double *vec1, const double *vec2, const int *group, int n) {
    int i, idx;
    for (i = 0; i < n; ++i) {
        idx = group[i];
        vec1[idx] -= vec2[idx];
    }
}

void scalarMult(double *vec, double x, const int *group, int n) {
    int i, idx;
    for (i = 0; i < n; ++i) {
        idx = *group;
        vec[idx] *= x;
        group++;
    }
}

/*this is a dot product of int vector with double vector*/
double dotProd(const int *vec1, const double *vec2, const int *group, int n) {
    double res = 0;
    int i, idx;
    for (i = 0; i < n; ++i) {
        idx = *(group + i);
        res += vec1[idx] * vec2[idx];
    }
    return res;
}

/* dot product of two double vectors*/
double dotDoubleProd(const double *vec1, const double *vec2, const int *group, int n) {
    double res = 0;
    int i, idx;
    for (i = 0; i < n; ++i) {
        idx = group[i];
        res += vec1[idx] * vec2[idx];
    }
    return res;
}

void copyVec(const int *src, double *dst, const int *group, int n) {
    int i, idx;
    for (i = 0; i < n; ++i) {
        idx = group[i];
        dst[idx] = src[idx];
    }
}

void copyDoubleVec(const double *src, double *dst, const int *group, int n) {
    int i, idx;
    for (i = 0; i < n; ++i) {
        idx = group[i];
        dst[idx] = src[idx];
    }
}

void normalize(int size, double *vec, const int *group, int groupSize) {
    double res = dotDoubleProd(vec, vec, group, groupSize);
    res = sqrt(res);
    scalarMult(vec, 1 / res, group, groupSize);
}


/*the calculation of Bv by split B to A, and KiKj matrix*/
void multBv(spmat *sp, double *vec, const int *group, double *res, int groupSize, int debug,const int *groupToVertice) {
    double dot;
    int size = sp->n;
    double *res1 = malloc(size * sizeof(double));
    dot = dotProd(sp->k, vec, group, groupSize);
//    if (debug == 1){
//        printf("dot is %f\n", dot);
//    }
    if (sp->M == 0) {
        printf("ERROR - divide in zero");
        exit(EXIT_FAILURE);
    }
    copyVec(sp->k, res1, group, groupSize);
    scalarMult(res1, (double) dot / sp->M, group, groupSize);
//    if (debug == 1){
//        printf("res1 is : ");
//        printVector(res1, size);
//    }
    sp->mult(sp, vec, res, group, groupSize,groupToVertice);

//    if (debug == 1){
//        printf("res after mult A is : ");
//        printVector(res, size);
//    }
    vecDec(res, res1, group, groupSize);
    free(res1);
}

/*get a vector and initialize it values to 1*/
void initOneValVec(double *unitVec, int n, const int *group, int val) {
    int i;
    for (i = 0; i < n; ++i) {
        *(unitVec + *(group + i)) = val;
    }
}


/* calculate the vector B^*v to res, by split the B^ into B and F vector as values of diag matrix*/
void multBRoof(spmat *sp, double *vec, const int *group, int groupSize, double *res, double *vecF, const int *vertexToGroup) {
    vecMult(vecF, vec, group, groupSize);

    multBv(sp, vec, group, res, groupSize, 0,vertexToGroup);

    vecDec(res, vecF, group, groupSize);
}


/*power iteration on B^ to calculate the leading eigenvalue, using matrix shifting*/
void powerIter(spmat *sp, double *b0, double shifting, int *group, int groupSize,double *vecF, double *result,const int *vertexToGroup) {
    int flag = 1, i, idx;
    int size = sp->n;
    int counter = 0;
    double *vecFCopy = malloc(size * sizeof(double));
    if (vecFCopy == NULL) {
        printf("ERROR - memory allocation unsuccessful");
        exit(EXIT_FAILURE);
    }
    while (flag == 1 && counter < 10000) {
        flag = 0;
        copyDoubleVec(vecF, vecFCopy, group, groupSize);
        multBRoof(sp, b0, group, groupSize, result, vecFCopy,vertexToGroup);
        vecSum(result, b0, shifting, group, groupSize);
        normalize(size, result, group, groupSize);
        for (i = 0; i < groupSize; i++) {
            idx = group[i];
            if (IS_POSITIVE(fabs(result[idx] - b0[idx])))
                flag = 1;
            b0[idx] = result[idx];
        }
        counter++;
    }
//    printf("took %d iterations\n", counter);
    free(vecFCopy);
}

/*calculate the eigenvalue of the leading eigenVector found*/
double eigenValue(spmat *sp, double *vec, const int *group, int groupSize,const int *vertexToGroup) {
    int size = sp->n;
    double *tmp = malloc(sizeof(double) * size);
    double res;
    double *vecF = malloc(size * sizeof(double));
    if (tmp == NULL || vecF == NULL) {
        printf("ERROR - memory allocation unsuccessful");
        exit(EXIT_FAILURE);
    }
    initOneValVec(tmp, groupSize, group, 1);
    multBv(sp, tmp, group, vecF, groupSize, 0,vertexToGroup);
    multBRoof(sp, vec, group, groupSize, tmp, vecF,vertexToGroup);
    res = dotDoubleProd(tmp, vec, group, groupSize);
    free(tmp);
    free(vecF);
//    printf("eigen value is %f\n", res);
    return res;
}

/*modularity calculation by multiply +-1 vector with B^*/
double modularityCalc(spmat *sp, double *vec, int *group, int groupSize,const int *vertexToGroup) {
    double res = 0;
    int size = sp->n;
    double *tmp = malloc(sizeof(double) * size);
    double *vecF = malloc(size * sizeof(double));
    if (tmp == NULL || vecF == NULL) {
        printf("ERROR - memory allocation unsuccessful");
        exit(EXIT_FAILURE);
    }
    initOneValVec(tmp, groupSize, group, 1);
    multBv(sp, tmp, group, vecF, groupSize, 0,vertexToGroup);
    multBRoof(sp, vec, group, groupSize, tmp, vecF,vertexToGroup);
    res = dotDoubleProd(tmp, vec, group, groupSize);
    free(tmp);
    free(vecF);
    return res / 2;
}

