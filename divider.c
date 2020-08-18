#include <stdio.h>
#include <stdlib.h>
#include "spmat.h"
#include <assert.h>
#include <math.h>
#include "divider.h"
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

void vecMult(const int *vec1, const double *vec2, double *res, int size) {
    int i;
    for (i = 0; i < size; ++i) {
//        printf("%f\n", vec1[i]);
        res[i] = vec1[i] * vec2[i];
    }
}

void vecSum(double *vec, const double *b0, double shifting, int n) {
    int i;
    for (i = 0; i < n; ++i) {
        vec[i] += b0[i] * shifting;
    }
}

void vecDec(double *vec1, const double *vec2,int group, const int *groupid,int n) {
    int i;
    for (i = 0; i < n; ++i) {
        if(group == groupid[i])
            vec1[i] -= vec2[i];
        else
            vec1[i] = 0;
    }
}

void scalarMult(double *vec, double x, int group, const int *groupid, int n) {
    int i;
    for (i = 0; i < n; ++i) {
        if (group == groupid[i])
            vec[i] *= x;
        else
            vec[i] = 0;
    }
}

double dotProd(const int *vec1, const double *vec2, int group, const int *groupid, int n) {
    double res = 0;
    int i;
    for (i = 0; i < n; ++i) {
        if (groupid[i] == group)
            res += vec1[i] * vec2[i];
    }
    return res;
}

double dotDoubleProd(const double *vec1, const double *vec2, int group, const int *groupid, int n) {
    double res = 0;
    int i;
    for (i = 0; i < n; ++i) {
        if(group == groupid[i])
            res += vec1[i] * vec2[i];
    }
    return res;
}

void copyVec(const int *src, double *dst, int n) {
    int i;
    for (i = 0; i < n; ++i) {
        dst[i] = src[i];
    }
}

void multBv(spmat *sp, double *vec, int group, const int *groupid, double *res) {
    double dot;
    double *res1 = malloc(sp->n * sizeof(double));
    dot = dotProd(sp->k, vec, group, groupid, sp->n);
    if (sp->M == 0) {
        printf("ERROR - divide in zero");
        exit(EXIT_FAILURE);
    }
    copyVec(sp->k, res1, sp->n);
    scalarMult(res1, (double) dot / sp->M, group, groupid,sp->n);
    sp->mult(sp, vec, res);
    vecDec(res, res1, group, groupid, sp->n);
    free(res1);
}

void powerIter(spmat *sp, double *b0, double shifting, int group, const int *groupid, double *result) {
    int flag = 1, i;
    while (flag == 1) {
        flag = 0;
        multBv(sp, b0, group, groupid, result);
        vecSum(result, b0, shifting, sp->n);
        normalize(sp->n, result);
        for (i = 0; i < sp->n; i++) {
            if (IS_POSITIVE(fabs(result[i] - b0[i])))
                flag = 1;
            b0[i] = result[i];
        }
    }
}


void freeDivision(division *d) {
    free(d->groupid);
}

double modularityCalc(spmat *sp, double *vec, int group, const int *groupid) {
    double res = 0;
    double *tmp = malloc(sizeof(double) * sp->n);
    multBv(sp, vec, group, groupid,tmp);
    res = dotDoubleProd(tmp, vec, group, groupid, sp->n);
    return res / 2;
}

void split(struct _division *d, spmat *sp, double *vec, int group) {
    int flag;
    int newGroup = -1;
    int i;
    double delta = modularityCalc(sp, vec, group, d->groupid);
    d->Q += delta;
    flag = IS_POSITIVE(vec[0]) ? 1 : 0;
    for (i = 1; i < sp->n; ++i) {
        if(group != d->groupid[i])
            continue;
        if (IS_POSITIVE(vec[i]) != flag) {
            if (newGroup == -1) {
                newGroup = d->numOfGroups;
                d->numOfGroups+=1;
            }
            d->groupid[i] = newGroup;
        }
    }
}

void printGroups(division *d){
    int i;
    printf("number of groups = %d\n", d->numOfGroups);
    for(i=0; i< d->n; ++i){
        printf("(%d,%d)\t",i,d->groupid[i]);
    }printf("\nmodularity value: %f\n",d->Q);
}


division *allocateDivision(int n) {
    int i;
    division *d = malloc(sizeof(division));
    if(d == NULL){
        printf("ERROR - memory allocation unsuccessful");
        exit(EXIT_FAILURE);
    }
    d->n = n;
    d->printGroups = printGroups;
    d->split = split;
    d->groupid = malloc(sizeof(int) * n);
    d->Q = 0;
    if (d->groupid == NULL) {
        printf("ERROR - memory allocation unsuccessful");
        exit(EXIT_FAILURE);
    }
    d->numOfGroups = 1;
    for (i = 0; i < n; ++i) {
        d->groupid[i] = 0;
    }
    return d;
}

