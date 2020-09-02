#include <stdio.h>
#include <stdlib.h>
#include "spmat.h"
#include <assert.h>
#include <math.h>
#include "divider.h"
#include <time.h>

#define IS_POSITIVE(X) ((X) > 0.00001)


void divOptimization(int *groupid, int size, int group, double *divVector, spmat *sp);

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
void vecMult(double *vec1, const double *vec2, int group, const int *groupid, int size) {
    int i;
    for (i = 0; i < size; ++i) {
        if (group == groupid[i])
            vec1[i] = vec1[i] * vec2[i];
        else
            vec1[i] = 0;
    }
}

/* gets to vectors and return the sum of vec with b0*shifting (shifting is the value from matrix shifting*/
void vecSum(double *vec, const double *b0, double shifting, int group, const int *groupid, int n) {
    int i;
    for (i = 0; i < n; ++i) {
        if (group == groupid[i])
            vec[i] += b0[i] * shifting;
        else
            vec[i] = 0;
    }
}

/*decrease vec1 by vec2*/
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

/*this is a dot product of int vector with double vector*/
double dotProd(const int *vec1, const double *vec2, int group, const int *groupid, int n) {
    double res = 0;
    int i;
    for (i = 0; i < n; ++i) {
        if (groupid[i] == group)
            res += vec1[i] * vec2[i];
    }
    return res;
}

/* dot product of two double vectors*/
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


/*the calculation of Bv by split B to A, and KiKj matrix*/
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
    sp->mult(sp, vec, res, group, groupid);
    vecDec(res, res1, group, groupid, size);
    free(res1);
}

/*get a vector and initialize it values to 1*/
void initUnitVec(double *unitVec, int n, int group, const int *groupid) {
    int i;
    for (i = 0; i < n; ++i) {
        if (group == groupid[i])
            unitVec[i] = 1;
        else
            unitVec[i] = 0;
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


/* calculate the vector B^*v to res, by split the B^ into B and F vector as values of diag matrix*/
void multBRoof(spmat *sp, double *vec, int group, const int *groupid, double *res) {
    int size = sp->n;
    double *unitVec = malloc(size * sizeof(double));
    double *vecF = malloc(size * sizeof(double));
    if (unitVec == NULL || vecF == NULL) {
        printf("ERROR - memory allocation unsuccessful");
        exit(EXIT_FAILURE);
    }
    initUnitVec(unitVec, size, group, groupid);
    multBv(sp, unitVec, group, groupid, vecF);
//    printf("vecF is \n");
//    printVector(vecF, size);
    vecMult(vecF, vec, group, groupid, size);
//    printf("vecF after mult with vec\n");
//    printVector(vecF,size);
//    printf("FFF\n");
//    printVector(vecF, size);
    multBv(sp, vec, group, groupid, res);
    vecDec(res, vecF, group, groupid, size);
    free(vecF);
    free(unitVec);
}


/*power iteration on B^ to calculate the leading eigenvalue, using matrix shifting*/
void powerIter(spmat *sp, double *b0, double shifting, int group, const int *groupid, double *result) {
    int flag = 1, i;
    int size = sp->n;
    int counter = 0;
    while (flag == 1 && counter < 10000) {
        flag = 0;
        multBRoof(sp, b0, group, groupid, result);
        vecSum(result, b0, shifting, group, groupid, size);
        normalize(size, result, group, groupid);
        for (i = 0; i < size; i++) {
            if (group != groupid[i])
                continue;
            if (IS_POSITIVE(fabs(result[i] - b0[i])))
                flag = 1;
            b0[i] = result[i];
        }
        counter++;
    }
//    printf("took %d iterations\n", counter);
}

/*calculate the eigenvalue of the leading eigenVector found*/
double eigenValue(spmat *sp, double *vec, int group, const int *groupid) {
    int size = sp->n;
    double *tmp = malloc(sizeof(double) * size);
    double res;
    multBRoof(sp, vec, group, groupid, tmp);
    res = dotDoubleProd(tmp, vec, group, groupid, size);
//    printVector(tmp, size);
//    div = dotDoubleProd(vec, vec, group, groupid, size);
//    printf("div in eigen calc is %f\n", div);
//    if (div == 0) {
//        printf("ERROR - divide in zero");
//        exit(EXIT_FAILURE);
//    }
    free(tmp);
    printf("eigen value is %f\n", res);
    return res;
}

/*modularity calculation by multiply +-1 vector with B^*/
double modularityCalc(spmat *sp, double *vec, int group, const int *groupid) {
    double res = 0;
    int size = sp->n;
    double *tmp = malloc(sizeof(double) * size);
    multBRoof(sp, vec, group, groupid, tmp);
    res = dotDoubleProd(tmp, vec, group, groupid, size);
    free(tmp);
    return res / 2;
}

double split(struct _division *d, spmat *sp, double *vec, int group) {
    int flag = -1;
    double delta;
    int newGroup = -1;
    int i;
    int size = sp->n;
    int *groupid = d->groupid;
    /* make the leading eigen-vector, a +-1 vector*/
    for (i = 0; i < size; ++i) {
        if (group != groupid[i])
            continue;
        if (flag == -1)
            flag = IS_POSITIVE(vec[i]) ? 1 : 0;
        if (IS_POSITIVE(vec[i]) != flag)
            vec[i] = -1;
        else
            vec[i] = 1;
    }
    divOptimization(groupid, d->nodesforGroup[group], group, vec, sp);
    delta = modularityCalc(sp, vec, group, groupid);
    if (!IS_POSITIVE(delta))
        return 0;

    delta = modularityCalc(sp, vec, group, groupid);
    d->Q += delta;
    /*create a new group for the -1 indexes in the +-1 vector*/
    for (i = 0; i < size; ++i) {
        if (group != groupid[i])
            continue;
        if (vec[i] == -1) {
            if (newGroup == -1) {
                newGroup = d->numOfGroups;
                d->numOfGroups += 1;
            }
            groupid[i] = newGroup;
            d->nodesforGroup[newGroup]++;
            d->nodesforGroup[group]--;
        }
    }
    return delta;
}

int divideToTwo(division *div, spmat *sp, int group, double *res, double *b0) {
//    printf("working on group: %d\n", group);
    int size = sp->n;
    double delta;
    randomizeVec(size, b0);
    int *groupid = div->groupid;
    double *unitVec = malloc(size * sizeof(double));
    double *vecF = malloc(size * sizeof(double));
    if (unitVec == NULL || vecF == NULL) {
        printf("ERROR - memory allocation unsuccessful");
        exit(EXIT_FAILURE);
    }
    initUnitVec(unitVec, size, group, groupid);
    multBv(sp, unitVec, group, groupid, vecF);
    powerIter(sp, b0, sp->matShifting(sp, group, groupid, vecF), group, groupid, res);
    double eigen = eigenValue(sp, res, group, div->groupid);
    if (!IS_POSITIVE(eigen))
        return 0;
    delta = split(div, sp, res, group);
    if (delta == 0)
        return 0;
    div->printGroups(div);
    return 1;
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

/*find vertice movement that maximizes Q*/
int
calcScore(double *Bvec, int *maxScoreIdx, int k, const double *divVec, spmat *sp, const int *unmoved, int group,
          int *groupID,
          double *score) {
    int i;
    int scoreIdx = 0;
    int maxIdx = 0;
    *maxScoreIdx = 0;
    for (i = 0; i < sp->n; i++) {
        if (group == groupID[i]) {
            if (unmoved[i]) {

                if (i == k) {
                    score[scoreIdx] = -1 * score[scoreIdx];  /*calcs new deltaQ */
                } else {

                    score[scoreIdx] = score[scoreIdx] -
                                      (4 * divVec[i] * divVec[k] * Bvec[i]);
                }
                if (score[scoreIdx] >= score[*maxScoreIdx]) {
                    *maxScoreIdx = scoreIdx;
                    maxIdx = i;
                }

            } /*keeps track of max Q */
            scoreIdx++;
        }
    }


    return maxIdx;

}

void initZeroVec(double *vec, int size, int group, const int *groupid) {
    int i;
    for (i = 0; i < size; ++i) {
        if (group == groupid[i])
            vec[i] = 0;
    }
}

/*optimizes division by moving one node from g1 to g2, saves division in res*/
void
optimize(double *divVec, int *unmoved, int *indices, double *improve, double *score, spmat *sp, int group,
         int *groupID, int size, double *vec) {
    if (sp->M == 0) {
        printf("ERROR - divide in zero");
        exit(EXIT_FAILURE);
    }
    int i;
    int scoreIdx = 0;
    int verticeIdx;
    int maxModularityIdx = size - 1;
    double deltaQ = 0;
    //TODO: FIX MEM RELEASE
    int *maxScoreIdx = malloc(sizeof(int));
    double *zeroVec = malloc(sizeof(double) * sp->n);
    double *bVec = malloc(sizeof(double) * sp->n);
    initZeroVec(zeroVec, sp->n, group, groupID);
    for (i = 0; i < sp->n; i++) { unmoved[i] = 1; }
    multBv(sp, divVec, group, groupID, vec);
    *maxScoreIdx = 0;

    for (i = 0; i < sp->n; i++) {
        if (groupID[i] == group) {
            score[scoreIdx] =
                    (-2) * ((divVec[i] * vec[i]) + (((double) sp->k[i] * (double) sp->k[i]) / (double) sp->M));
            if (score[scoreIdx] >= score[*maxScoreIdx]) {
                *maxScoreIdx = scoreIdx;
                verticeIdx = i;
            }
            scoreIdx++;
        }
    }

    scoreIdx = 1;
    unmoved[verticeIdx] = 0;
    divVec[verticeIdx] = (-1) * divVec[verticeIdx];
    indices[0] = verticeIdx;
    improve[0] = score[*maxScoreIdx];
    zeroVec[verticeIdx] = 1;
    multBv(sp, zeroVec, group, groupID, bVec);
    zeroVec[verticeIdx] = 0;


    while (scoreIdx < size) {

        verticeIdx = calcScore(bVec, maxScoreIdx, verticeIdx, divVec, sp, unmoved, group, groupID, score);
        divVec[verticeIdx] = (-1) * divVec[verticeIdx];
        indices[scoreIdx] = verticeIdx;
        improve[scoreIdx] = improve[scoreIdx - 1] + score[*maxScoreIdx];
        unmoved[verticeIdx] = 0;
        zeroVec[verticeIdx] = 1;
        multBv(sp, zeroVec, group, groupID, bVec);
        zeroVec[verticeIdx] = 0;
        scoreIdx++;
    }
    free(maxScoreIdx);

    for (i = 0; i < size; i++) {
        if (improve[i] >
            improve[maxModularityIdx]) { maxModularityIdx = i; }   /*saves max division state */}

    if ((maxModularityIdx == (size - 1)) || (improve[maxModularityIdx] < 0)) {
        deltaQ = 0;
        for (i = 0; i < size; i++) { divVec[indices[i]] = (-1) * divVec[indices[i]]; }
    } else {
        deltaQ = improve[maxModularityIdx];
        for (i = (size - 1); i > maxModularityIdx; i--) { divVec[indices[i]] = (-1) * divVec[indices[i]]; }
    }
    if (IS_POSITIVE(deltaQ)) {
        optimize(divVec, unmoved, indices, improve, score, sp, group, groupID, size, vec);
    }  /*finds max division state, if Q is positive, we try again */
}


void divOptimization(int *groupid, int size, int group, double *divVector, spmat *sp) {
    printf("optimization called\n");
    printf("mod calc before %f\n",modularityCalc(sp,divVector,group,groupid));
    int *unmoved = malloc(sizeof(int) * sp->n);
    int *indices = malloc(sizeof(int) * size);
    double *improve = malloc(sizeof(double) * size);
    double *score = malloc(size * sizeof(double));
    double *vec = malloc(sp->n * sizeof(double));
    optimize(divVector, unmoved, indices, improve, score, sp, group, groupid, size, vec);
    printf("mod calc after %f\n",modularityCalc(sp,divVector,group,groupid));
//TODO: fix memory release
// free(unmoved);
//    free(indices);
//    free(improve);
//    free(score);
//  free(vec);
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
        }
        last++;
    }
//    div->printGroups(div);
    free(b0);
    free(res);
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
    d->free = freeDivision;
    d->groupid = malloc(sizeof(int) * n);
    d->nodesforGroup = malloc(sizeof(int) * n);
    d->writeDivision = writeDivision;
    d->findGroups = findGroups;
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

//void resetUnmoved(int *groupid, int group, int *unmoved){
//
//}
//
//void divisionOptimizer(division *div, spmat *sp, double *vec, int group) {
//    double delta = 1;
//    int size = sp->n;
//    int *unmoved = malloc(sizeof(int) * size);
//    if (unmoved == NULL) {
//        printf("ERROR - memory allocation unsuccessful");
//        exit(EXIT_FAILURE);
//    }
//    while(IS_POSITIVE(delta)){
//        unmoved = resetUnmoved(div->groupid, group, unmoved);
//    }
//}

