#include <stdio.h>
#include <stdlib.h>
#include "spmat.h"
#include <assert.h>
#include <math.h>
#include "divider.h"
#include <time.h>

#define IS_POSITIVE(X) ((X) > 0.00001)


/*get a vector and the size of it and init any value to a random value*/
void randomizeVec(int size, double *vec) {
    int i;
    srand(time(NULL));
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
void vecSum(double *vec, const double *b0, double shifting,const int *group, int n) {
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
    int i,idx;
    for (i = 0; i < n; ++i) {
        idx = *group;
        vec[idx] *= x;
        group++;
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

/*this is a dot product of int vector with double vector*/
double dotProd(const int *vec1, const double *vec2, const int *group, int n) {
    double res = 0;
    int i, idx;
    for (i = 0; i < n; ++i) {
        idx = *(group+i);
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

void normalize(int size, double *vec,const int *group, int groupSize) {
    double res = dotDoubleProd(vec, vec, group, groupSize);
    res = sqrt(res);
    scalarMult(vec, 1 / res, group, groupSize);
}


/*the calculation of Bv by split B to A, and KiKj matrix*/
void multBv(spmat *sp, double *vec, const int *group, double *res, int groupSize, int debug) {
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
    sp->mult(sp, vec, res, group, groupSize);
//    if (debug == 1){
//        printf("res after mult A is : ");
//        printVector(res, size);
//    }
    vecDec(res, res1, group, groupSize);
    free(res1);
}

/*get a vector and initialize it values to 1*/
void initUnitVec(double *unitVec, int n, const int *group) {
    int i;
    for (i = 0; i < n; ++i) {
        *(unitVec + *(group+i)) = 1;
    }
}




/* calculate the vector B^*v to res, by split the B^ into B and F vector as values of diag matrix*/
void multBRoof(spmat *sp, double *vec, const int *group, int groupSize, double *res) {
    int size = sp->n;
    double *unitVec = malloc(size * sizeof(double));
    double *vecF = malloc(size * sizeof(double));
    if (unitVec == NULL || vecF == NULL) {
        printf("ERROR - memory allocation unsuccessful");
        exit(EXIT_FAILURE);
    }
    initUnitVec(unitVec, groupSize, group);
    multBv(sp, unitVec, group, vecF, groupSize, 0);
    vecMult(vecF, vec, group, groupSize);
    multBv(sp, vec, group, res, groupSize, 0);
    vecDec(res, vecF, group, groupSize);
    free(vecF);
    free(unitVec);
}


/*power iteration on B^ to calculate the leading eigenvalue, using matrix shifting*/
void powerIter(spmat *sp, double *b0, double shifting, int* group, int groupSize, double *result) {
    int flag = 1, i,idx;
    int size = sp->n;
    int counter = 0;
    while (flag == 1 && counter < 10000) {
        flag = 0;
        multBRoof(sp, b0, group, groupSize,result);
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
}

/*calculate the eigenvalue of the leading eigenVector found*/
double eigenValue(spmat *sp, double *vec, const int *group, int groupSize) {
    int size = sp->n;
    double *tmp = malloc(sizeof(double) * size);
    double res;
    multBRoof(sp, vec, group, groupSize , tmp);
    res = dotDoubleProd(tmp, vec, group, groupSize);
    free(tmp);
//    printf("eigen value is %f\n", res);
    return res;
}

/*modularity calculation by multiply +-1 vector with B^*/
double modularityCalc(spmat *sp, double *vec, int* group, int groupSize) {
    double res = 0;
    int size = sp->n;
    double *tmp = malloc(sizeof(double) * size);
    if (tmp == NULL) {
        printf("ERROR - memory allocation unsuccessful");
        exit(EXIT_FAILURE);
    }
    multBRoof(sp, vec, group, groupSize, tmp);
    res = dotDoubleProd(tmp, vec, group, groupSize);
    free(tmp);
    return res / 2;
}

double split(struct _division *d, spmat *sp, double *vec, int groupIdx) {
    int flag = -1;
    double delta;
    int newGroupIdx = -1;
    int i;
    int size = d->nodesforGroup[groupIdx];
    int **groups = d->groups;
    int *g = d->groups[groupIdx];
    int *tempGroup;
    int *g1Ptr, *g2Ptr;
    int counter = 0;
    /* make the leading eigen-vector, a +-1 vector*/
    for (i = 0; i < size; ++i) {
        if (flag == -1)
            flag = IS_POSITIVE(vec[g[i]]) ? 1 : 0;
        if (IS_POSITIVE(vec[g[i]]) != flag){
            vec[g[i]] = -1;
            counter ++;
        }
        else
            vec[g[i]] = 1;
    }
    delta = modularityCalc(sp, vec, g, size);
    if (!IS_POSITIVE(delta))
        return 0;
    d->Q += delta;
    /*create a new group for the -1 indexes in the +-1 vector*/
    newGroupIdx = d->numOfGroups;
    d->numOfGroups += 1;
    groups[newGroupIdx] = malloc(sizeof(int)*counter);
    tempGroup = malloc(sizeof(int)*(size-counter));
    if (groups[newGroupIdx] == NULL || tempGroup == NULL) {
        printf("ERROR - memory allocation unsuccessful");
        exit(EXIT_FAILURE);
    }
    g2Ptr = groups[newGroupIdx];
    g1Ptr = tempGroup;
    for (i = 0; i < size; ++i) {
        if (vec[g[i]] == -1) {
            *g2Ptr = g[i];
            d->vertexToGroup[g[i]] = newGroupIdx;
            g2Ptr++;

        }
        else{
            *g1Ptr = g[i];
            g1Ptr++;
        }
    }
    d->nodesforGroup[newGroupIdx] = counter;
    free(g);
    d->groups[groupIdx]=tempGroup;
    d->nodesforGroup[groupIdx] = size - counter;
    return delta;
}

int divideToTwo(division *div, spmat *sp, int groupIdx, double *res, double *b0) {
    int size = sp->n;
    double delta;
    randomizeVec(size, b0);
    int *group = *(div->groups+groupIdx);
    int groupSize = div->nodesforGroup[groupIdx];
    double *unitVec = malloc(size * sizeof(double));
    double *vecF = malloc(size * sizeof(double));
    if (unitVec == NULL || vecF == NULL) {
        printf("ERROR - memory allocation unsuccessful");
        exit(EXIT_FAILURE);
    }
    initUnitVec(unitVec, div->nodesforGroup[groupIdx], group);
    multBv(sp, unitVec, group, vecF, groupSize, 1);
//    printVector(vecF, size);
    //TODO fix mult Bv and shifting
//    printf("vec f is :");
//    printVector(vecF, size);
//    printf("shifting value is %f\n", sp->matShifting(sp, group, groupSize, div->vertexToGroup, groupIdx,vecF));
    powerIter(sp, b0, sp->matShifting(sp, group, groupSize, div->vertexToGroup, groupIdx, vecF), group, groupSize, res);
//    printf("HERE11\n");
    double eigen = eigenValue(sp, res, group, groupSize);
    if (!IS_POSITIVE(eigen))
        return 0;
    delta = split(div, sp, res, groupIdx);
    if (delta == 0)
        return 0;
//    div->printGroups(div);
    return 1;
}

void printGroups(division *d) {
    int i, j;
    printf("number of groups = %d\n", d->numOfGroups);
    for (i = 0; i < d->numOfGroups; ++i) {
        int *group = d->groups[i];
        for(j=0; j< d->nodesforGroup[i]; ++j){
            printf("(%d,%d)\t", group[j], i);
        }
    }
    printf("\nmodularity value: %f\n", d->Q);
}

void freeDivision(division *d) {
    int i;
    for (i = 0; i < d->numOfGroups; ++i) {
        free(d->groups[i]);
    }
    free(d->groups);
    free(d->nodesforGroup);
}

void writeDivision(struct _division *div, FILE *output) {
    int numOfGroups = div->numOfGroups;
    int *vertexForGroup = div->nodesforGroup;
    int **groups = div->groups;
    int i;
    fwrite(&numOfGroups, sizeof(int), 1, output);
    for (i = 0; i < numOfGroups; ++i) {
        fwrite(&vertexForGroup[i], sizeof(int), 1, output);
        fwrite(groups[i], sizeof(int), vertexForGroup[i], output);
    }
}


void findGroups(division *div, spmat *sp) {
    double delta;
    int size = sp->n;
    int last = div->numOfGroups - 1;
    double *b0 = malloc(sizeof(double) * size);
    double *res = malloc(sizeof(double) * size);
    if (b0 == NULL || res == NULL) {
        printf("ERROR - memory allocation unsuccessful");
        exit(EXIT_FAILURE);
    }
    while (last < div->numOfGroups) {
        delta = 1;
        while (delta == 1) {
            delta = divideToTwo(div, sp, last, res, b0);
//            divOptimization(div, last, delta, res, sp);
        }
        last++;
    }
    printf("modularity is %f\n",div->Q);
    free(b0);
    free(res);
}

division *allocateDivision(int n) {
    int i, *ptr;
    division *d = malloc(sizeof(division));
    if (d == NULL) {
        printf("ERROR - memory allocation unsuccessful");
        exit(EXIT_FAILURE);
    }
    d->n = n;
    d->printGroups = printGroups;
    d->free = freeDivision;
    d->groups = malloc(sizeof(int *) * n);
    d->nodesforGroup = malloc(sizeof(int) * n);
    d->vertexToGroup = malloc(sizeof(int) * n);
    d->writeDivision = writeDivision;
    d->findGroups = findGroups;
    d->Q = 0;
    d->groups[0] = malloc(sizeof(int) * n);
    if (d->groups == NULL || d->nodesforGroup == NULL || d->groups[0] == NULL || d->vertexToGroup == NULL) {
        printf("ERROR - memory allocation unsuccessful");
        exit(EXIT_FAILURE);
    }
    d->numOfGroups = 1;
    ptr = d->groups[0];
    for (i = 0; i < n; ++i) {
        d->vertexToGroup[i] = 0;
        *ptr = i;
        ptr++;
    }
    d->nodesforGroup[0] = n;
    return d;
}
