//
// Created by gal21 on 05/09/2020.
//
#include "utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "spmat.h"
#include <time.h>

void printVector(double *vec, int n, const int *group) {
    int i, idx;
    for (i = 0; i < n; ++i) {
//        idx = group[i];
        printf("%f\t", vec[i]);
    }
    printf("\n");
}

void printIntVector(int *vec, int n) {
    int i;
    for (i = 0; i < n; ++i) {
        printf("%d\t", vec[i]);
    }
    printf("\n");
}


//*get a vector and the size of it and init any value to a random value*/
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

/* used for the F vector as a Matrix to multiply it by the v vector*/
void vecMult(double *vec1, const double *vec2, const int *group, int size) {
    int i, delta, prev = *group;
//    vec1 += prev;
//    vec2 += prev;
    for (i = 0; i < size; ++i) {
//        idx = group[i];
//        vec1Ptr += (*groupPtr);
        *vec1 = (*vec1) * (*vec2);
        vec1++;
        vec2++;
//        vec1[idx] = vec1[idx] * vec2[idx];
//        delta = *++group - prev;
//        prev = *group;
//        vec1 += delta;
//        vec2 += delta;
    }
}

/* gets to vectors and return the sum of vec with b0*shifting (shifting is the value from matrix shifting*/
void vecSum(double *vec, const double *b0, double shifting, const int *group, int n) {
    int i, idx, prev = *group, delta;
//    b0 += prev;
//    vec += prev;
    for (i = 0; i < n; ++i) {
//        idx = group[i];
//        vec[idx] += b0[idx] * shifting;
//        *vec += *b0 * shifting;
//        delta = *++group - prev;
//        prev = *group;
//        vec += delta;
//        b0 += delta;
        *vec += *b0++ * shifting;
        vec++;

    }
}

/*decrease vec1 by vec2*/
void vecDec(double *vec1, const double *vec2, const int *group, int n) {
    int i, idx, delta, prev = *group;
//    vec1 += prev;
//    vec2 += prev;
    for (i = 0; i < n; ++i) {
//        idx = group[i];
//        vec1[idx] -= vec2[idx];
//        *vec1 -= *vec2;
//        delta = *++group - prev;
//        prev = *group;
//        vec1 += delta;
//        vec2 += delta;
        *vec1 -= *vec2++;
        vec1++;
    }
}

void scalarMult(double *vec, double x, const int *group, int n) {
    int i, idx, prev = *group, delta;
//    vec += prev;
    for (i = 0; i < n; ++i) {
//        idx = *group;
//        vec[idx] *= x;
//        group++;
//        *vec *= x;
//        delta = *++group - prev;
//        prev = *group;
//        vec += delta;
        *vec++ *= x;
    }
}


/*this is a dot product of int vector with double vector*/
double dotProd(const int *vec1, const double *vec2, const int *group, int n) {
    double res = 0;
    int i, idx, prev = *group, delta;
//    vec1 += prev;
//    vec2 += prev;
    for (i = 0; i < n; ++i) {
//        idx = group[i];
//        res += vec1[idx] * vec2[idx];
//        res += *vec1 * *vec2;
//        delta = *++group - prev;
//        prev = *group;
//        vec1 += delta;
//        vec2 += delta;
        res += *vec1++ * *vec2++;
    }
    return res;
}

/* dot product of two double vectors*/
double dotDoubleProd(const double *vec1, const double *vec2, const int *group, int n) {
    double res = 0;
    int i, idx, prev = *group, delta;
//    vec1 += prev;
//    vec2 += prev;
    for (i = 0; i < n; ++i) {
//        idx = group[i];
//        res += vec1[idx] * vec2[idx];
//        res += *vec1 * *vec2;
//        delta = *++group - prev;
//        prev = *group;
//        vec1 += delta;
//        vec2 += delta;
        res += *vec1++ * *vec2++;
    }
    return res;
}

void copyDoubleVec(const double *src, double *dst, const int *group, int n) {
    int i, idx, prev = *group, delta;
//    src += prev;
//    dst += prev;
    for (i = 0; i < n; ++i) {
//        idx = group[i];
//        dst[idx] = src[idx];
//        *dst = *src;
//        delta = *++group - prev;
//        prev = *group;
//        dst += delta;
//        src += delta;
        *dst = *src++;
        dst++;
    }
}

void copyVec(const int *src, double *dst, const int *group, int n) {
    int i, idx, prev = *group, delta;
//    src += prev;
//    dst += prev;
    for (i = 0; i < n; ++i) {
//        idx = group[i];
//        dst[idx] = src[idx];
//        *dst = *src;
//        delta = *++group - prev;
//        prev = *group;
//        dst += delta;
//        src += delta;

        *dst = *src++;
        dst++;
    }
}

void normalize(int size, double *vec, const int *group, int groupSize) {
    double res = dotDoubleProd(vec, vec, group, groupSize);
    res = sqrt(res);
    scalarMult(vec, 1 / res, group, groupSize);
}


/*the calculation of Bv by split B to A, and KiKj matrix*/
double multBv(spmat *sp, double *vec, const int *group, double *res, int groupSize, int debug, int *verticeToGroup) {
    double dot;
    int size = sp->n;
    double *res1 = malloc(size * sizeof(double));
    double start, end;
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
    start = clock();
    sp->mult(sp, vec, res, group, groupSize, verticeToGroup);
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
    vecDec(res, res1, group, groupSize);
    free(res1);
    return ((double)(end-start)/ CLOCKS_PER_SEC);
}

/*get a vector and initialize it values to 1*/
void initOneValVec(double *unitVec, int n, const int *group, int val) {
    int i, prev = *group, delta;
//    unitVec += prev;
    for (i = 0; i < n; ++i) {
//        *(unitVec + *(group + i)) = val;
//        *unitVec = val;
//        delta = *++group - prev;
//        prev = *group;
//        unitVec += delta;

        *unitVec = val;
        unitVec++;
    }
}


/* calculate the vector B^*v to res, by split the B^ into B and F vector as values of diag matrix*/
double multBRoof(spmat *sp, double *vec, const int *group, int groupSize, double *res, double *vecF, int *verticeToGroup,
               int debug) {
    int size = sp->n;
//    double *unitVec = malloc(size * sizeof(double));
    double *vecFCopy = malloc(size * sizeof(double));
    double total, start, end;
    if (vecFCopy == NULL) {
        printf("ERROR - memory allocation unsuccessful");
        exit(EXIT_FAILURE);
    }
    start = clock();
//    initOneValVec(unitVec, groupSize, group, 1);
//    multBv(sp, unitVec, group, vecF, groupSize, 0);
    copyDoubleVec(vecF, vecFCopy, group, groupSize);
    vecMult(vecFCopy, vec, group, groupSize);
    total = multBv(sp, vec, group, res, groupSize, debug, verticeToGroup);
//    if (debug == 1) {
//        printVector(vec, size, group);
//            printVector(result, size, group);
//        printf("PI multBv total time took %f\n", total);
//    }
    vecDec(res, vecFCopy, group, groupSize);
    free(vecFCopy);
    debug = 0;
//    free(unitVec);
    end = clock();
//    return ((double)(end-start)/ CLOCKS_PER_SEC);
    return total;
}


/*power iteration on B^ to calculate the leading eigenvalue, using matrix shifting*/
void powerIter(spmat *sp, double *b0, double shifting, int *group, int groupSize, double *result, double *vecF,
               int *verticeToGroup, int debug) {
    int flag = 1, i, idx;
    int size = sp->n;
    int counter = 0;
    double total = 0;
    while (flag == 1 && counter < 10000) {
        flag = 0;
        total +=multBRoof(sp, b0, group, groupSize, result, vecF, verticeToGroup, debug);
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
//    printf("took %d iterations\n", counter);
//    printf("mult A in PI took %f\n", total);
}

/*calculate the eigenvalue of the leading eigenVector found*/
double eigenValue(spmat *sp, double *vec, const int *group, int groupSize, double *vecF, int *verticeToGroup) {
    int size = sp->n;
    double *tmp = malloc(sizeof(double) * size);
    double res;
    multBRoof(sp, vec, group, groupSize, tmp, vecF, verticeToGroup, 0);
    res = dotDoubleProd(tmp, vec, group, groupSize);
    free(tmp);
    //   printf("eigen value is %f\n", res);
    return res;
}

/*modularity calculation by multiply +-1 vector with B^*/
double modularityCalc(spmat *sp, double *vec, int *group, int groupSize, double *vecF, int *verticeToGroup) {
    double res = 0;
    int size = sp->n;
    double *tmp = malloc(sizeof(double) * size);
    if (tmp == NULL) {
        printf("ERROR - memory allocation unsuccessful");
        exit(EXIT_FAILURE);
    }
    multBRoof(sp, vec, group, groupSize, tmp, vecF, verticeToGroup, 0);
    res = dotDoubleProd(tmp, vec, group, groupSize);
    free(tmp);
    return res / 2;
}
