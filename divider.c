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

void powerIter(spmat *sp, double *b0, double *result) {
    int flag = 1, i;
    while (flag == 1) {
        flag = 0;
        sp->mult(sp, b0, result);
        normalize(sp->n, result);
        for(i=0; i<sp->n; i++){
            if(IS_POSITIVE(fabs(result[i] - b0[i])))
                flag = 1;
            b0[i] = result[i];
        }
    }
}

void matrixShifting(spmat *A ,double *output){
    int i;
    int j;
    double sum;
    double max;
    for(i = 0; i< A->n)
}
