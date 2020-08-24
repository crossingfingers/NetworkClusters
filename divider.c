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

void vecMult(double *vec1, const double *vec2, int group,const int *groupid, int size) {
    int i;
    for (i = 0; i < size; ++i) {
        if(group == groupid[i])
            vec1[i] = vec1[i] * vec2[i];
        else
            vec1[i] = 0;
    }
}

void vecSum(double *vec, const double *b0, double shifting, int group, const int *groupid, int n) {
    int i;
    for (i = 0; i < n; ++i) {
        if (group == groupid[i])
            vec[i] += b0[i] * shifting;
        else
            vec[i] = 0;
    }
}

void vecDec(double *vec1, const double *vec2, int group, const int *groupid, int n) {
    int i;
    for (i = 0; i < n; ++i) {
        if (group == groupid[i])
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
        if (group == groupid[i])
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

void normalize(int size, double *vec, int group, const int *groupid) {
    double res = dotDoubleProd(vec, vec, group, groupid, size);
    res = sqrt(res);
    scalarMult(vec, 1 / res, group, groupid, size);
}

void multBv(spmat *sp, double *vec, int group, const int *groupid, double *res) {
    double dot;
    int size = sp->n;
    double *res1 = malloc(size * sizeof(double));
    dot = dotProd(sp->k, vec, group, groupid, size);
    if (sp->M == 0) {
        printf("ERROR - divide in zero");
        exit(EXIT_FAILURE);
    }
    copyVec(sp->k, res1, size);
    scalarMult(res1, (double) dot / sp->M, group, groupid, size);
    sp->mult(sp, vec, res);
    vecDec(res, res1, group, groupid, size);
    free(res1);
}

void initUnitVec(double *unitVec, int n) {
    int i;
    for (i = 0; i < n; ++i) {
        unitVec[i] = 1;
    }
}

void printVector(double *vec, int n) {
    int i;
    for (i = 0; i < n; ++i) {
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

void multBRoof(spmat *sp, double *vec, int group, const int *groupid, double *res) {
    int size = sp->n;
    double *unitVec = malloc(size * sizeof(double));
    double *vecF = malloc(size * sizeof(double));
    if (unitVec == NULL || vecF == NULL) {
        printf("ERROR - memory allocation unsuccessful");
        exit(EXIT_FAILURE);
    }
    initUnitVec(unitVec, size);
    multBv(sp, unitVec, group, groupid, vecF);
    vecMult(vecF, vec, group, groupid, size);
//    printf("FFF\n");
//    printVector(vecF, size);
    multBv(sp, vec, group, groupid, res);
    vecDec(res, vecF, group, groupid, size);
    free(vecF);
    free(unitVec);
}

void powerIter(spmat *sp, double *b0, double shifting, int group, const int *groupid, double *result) {
    int flag = 1, i;
    int size = sp->n;
    while (flag == 1) {
        flag = 0;
        multBRoof(sp, b0, group, groupid, result);
        vecSum(result, b0, shifting, group, groupid, size);
        normalize(size, result, group, groupid);
        for (i = 0; i < size; i++) {
            if (IS_POSITIVE(fabs(result[i] - b0[i])))
                flag = 1;
            b0[i] = result[i];
        }
    }
}

double eigenValue(spmat *sp, double *vec, int group, const int *groupid) {
    int size = sp->n;
    double *tmp = malloc(sizeof(double) * size);
    double res;
    double div;
    multBRoof(sp, vec, group, groupid, tmp);
    res = dotDoubleProd(tmp, vec, group, groupid, size);
//    printVector(tmp, size);
    div = dotDoubleProd(vec, vec, group, groupid, size);
    if (div == 0) {
        printf("ERROR - divide in zero");
        exit(EXIT_FAILURE);
    }
    free(tmp);
    return res / div;
}

double modularityCalc(spmat *sp, double *vec, int group, const int *groupid) {
    double res = 0;
    int size = sp->n;
    double *tmp = malloc(sizeof(double) * size);
    multBRoof(sp, vec, group, groupid, tmp);
    res = dotDoubleProd(tmp, vec, group, groupid, size);
    free(tmp);
    return res / 2;
}

int split(struct _division *d, spmat *sp, double *vec, int group) {
    int flag;
    double delta;
    int newGroup = -1;
    int i;
    int size = sp->n;
    int *groupid = d->groupid;
    int *copyGroup = malloc(sizeof(int) * size);
    if (copyGroup == NULL) {
        printf("ERROR - memory allocation unsuccessful");
        exit(EXIT_FAILURE);
    }
    flag = IS_POSITIVE(vec[0]) ? 1 : 0;
    copyGroup[0] = groupid[0];
    vec[0] = 1;
    for (i = 1; i < size; ++i) {
        copyGroup[i] = groupid[i];
        if (group != groupid[i])
            continue;
        if (IS_POSITIVE(vec[i]) != flag) {
            if (newGroup == -1) {
                newGroup = d->numOfGroups;
                d->numOfGroups += 1;
            }
            groupid[i] = newGroup;
            d->nodesforGroup[newGroup]++;
            d->nodesforGroup[group]--;
            vec[i] = -1;
        } else {
            vec[i] = 1;
        }
    }
    delta = modularityCalc(sp, vec, group, copyGroup);
    if (!IS_POSITIVE(delta))
        return 0;
    d->Q += delta;
    free(copyGroup);
    return 1;
}

int divideToTwo(division *div, spmat *sp, int group) {
    printf("working on group: %d\n", group);
    int size = sp->n;
    double *b0 = malloc(sizeof(double) * size);
    double *res = malloc(sizeof(double) * size);
    if (b0 == NULL || res == NULL) {
        printf("ERROR - memory allocation unsuccessful");
        exit(EXIT_FAILURE);
    }
    randomizeVec(size, b0);
    powerIter(sp, b0, sp->matShifting(sp, group, div->groupid), group, div->groupid, res);
//    printVector(res, size);
    double eigen = eigenValue(sp, res, group, div->groupid);
    printf("eigen %f\n", eigen);
    if (!IS_POSITIVE(eigen))
        return 0;
    if (div->split(div, sp, res, group) == 0)
        return 0;
    div->printGroups(div);
    free(b0);
    free(res);
    return 1;
}

void findGroups(division *div, spmat *sp) {
    int flag;
    int last = div->numOfGroups - 1;
    while (last < div->numOfGroups) {
        flag = 1;
        while (flag == 1) {
            printf("last is %d\n", last);
            flag = divideToTwo(div, sp, last);
        }
        last++;
    }
    div->printGroups(div);
}

void printGroups(division *d) {
    int i;
    printf("number of groups = %d\n", d->numOfGroups);
    for (i = 0; i < d->n; ++i) {
        printf("(%d,%d)\t", i, d->groupid[i]);
    }
    printf("\nmodularity value: %f\n", d->Q);
}

void freeDivision(division *d) {
    free(d->groupid);
}

void writeDivision(struct _division *div, FILE *output) {
    int numOfGroups = div->numOfGroups;
    int *vertexForGroup = div->nodesforGroup;
    int *groupid = div->groupid;
    int i;
    int n = div->n;
    int j;
    fwrite(&numOfGroups, sizeof(int), 1, output);
    for (i = 0; i < numOfGroups; ++i) {
        fwrite(&vertexForGroup[i], sizeof(int), 1, output);
        for (j = 0; j < n; ++j) {
            if (groupid[j] == i)
                fwrite(&j, sizeof(int), 1, output);
        }
    }
}

division *allocateDivision(int n) {
    int i;
    division *d = malloc(sizeof(division));
    if (d == NULL) {
        printf("ERROR - memory allocation unsuccessful");
        exit(EXIT_FAILURE);
    }
    d->n = n;
    d->printGroups = printGroups;
    d->split = split;
    d->free = freeDivision;
    d->groupid = malloc(sizeof(int) * n);
    d->nodesforGroup = malloc(sizeof(int) * n);
    d->writeDivision = writeDivision;
    d->Q = 0;
    if (d->groupid == NULL || d->nodesforGroup == NULL) {
        printf("ERROR - memory allocation unsuccessful");
        exit(EXIT_FAILURE);
    }
    d->numOfGroups = 1;
    for (i = 0; i < n; ++i) {
        d->groupid[i] = 0;
        d->nodesforGroup[i] = 0;
    }
    d->nodesforGroup[0] = n;
    return d;
}

