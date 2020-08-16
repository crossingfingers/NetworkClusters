#include <stdio.h>
#include <stdlib.h>
#include "spmat.h"
#include "divider.h"
#include <assert.h>
#include <math.h>
#include <time.h>
#define IS_POSITIVE(X) ((X) > 0.00001)

void randomizeVec(int size, double *vec) {
    int i;
    srand(time(NULL));
    assert(vec != NULL);
    for (i = 0; i < size; i++) {
        vec[i] = rand();
    }
}


void normalize(int size, double *vec) {
    int i;
    double res = 0;
    for (i = 0; i < size; i++) {
        res += vec[i] * vec[i];
    }
    res = sqrt(res);
    for (i = 0; i < size; i++) {
        vec[i] /= res;
    }
}

void vecMult(const int *vec1, const double* vec2, double *res, int size){
    int i;
    for(i=0; i<size; ++i){
//        printf("%f\n", vec1[i]);
        res[i] = vec1[i] * vec2[i];
    }
}

void vecSum(double *vec, double* b0, double shifting, int n){
    int i;
    for(i=0; i<n; ++i){
        vec[i] += b0[i] * shifting;
    }
}

void vecDec(double *vec1, double *vec2, int n){
    int i;
    for(i=0; i<n; ++i){
        vec1[i] -= vec2[i];
    }
}

void printVector(double *vec, int n){
    int i;
    for(i = 0; i< n; ++i){
        printf("%d\t", vec[i]);
    }printf("\n");
}

void powerIter(spmat *sp, double *b0, double shifting, double *result) {
    int flag = 1, i;
    double * res1 = malloc(sp->n * sizeof(double));
    while (flag == 1) {
        flag = 0;
        vecMult(sp->k, b0, res1, sp->n);
        vecMult(sp->k, res1, res1, sp->n);
        sp->mult(sp, b0, result);
        vecDec(result, res1, sp->n);
        vecSum(result, b0, shifting, sp->n);
        normalize(sp->n, result);
        for(i=0; i<sp->n; i++){
//            if(IS_POSITIVE(fabs(result[i] - b0[i])))
            printf("%f\n",fabs(result[i] - b0[i]));
            if(fabs(result[i] - b0[i]) >= 0.00001)
                flag = 1;
            b0[i] = result[i];
        }
    }
    free(res1);
}


