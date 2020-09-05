#include <stdio.h>
#include <stdlib.h>
#include "spmat.h"
#include "divider.h"
#include <float.h>
#include "utils.h"


void resetUnmoved(int *unmoved, int groupSize) {
    int i;
    for (i = 0; i < groupSize; ++i) {
        unmoved[i] = 1;
    }
}


void updateScore(spmat *sp, const double *s, double *score, const int *group, int groupSize, const int *unmoved, int k,
                 double *zeroVec, double *res) {
    int i, idx;
    zeroVec[k] = 1;
    multBv(sp, zeroVec, group, res, groupSize, 0);
    for (i = 0; i < groupSize; ++i) {
        if (unmoved[i] == 0)
            continue;
        idx = group[i];
        if (k == idx)
            score[i] = -score[i];
        else {
            score[i] = score[i] - (4 * s[idx] * s[k] * res[idx]);
        }
    }
    zeroVec[k] = 0;
}

int findMaxIdx(const double *score, int groupSize, const int *unmoved) {
    double max = -DBL_MAX;
    int maxIdx=-1, i, flag = 0;
    for (i = 0; i < groupSize; ++i) {
        if (unmoved[i] == 0)
            continue;
        if (flag == 0 || max < score[i]) {
            flag = 1;
            max = score[i];
            maxIdx = i;
        }
    }
//    printf("maxidx is %d\n", maxIdx);
    return maxIdx;
}

/*int findMaxIdx(spmat *sp, double *s, double *score, int group, const int *groupid, const int *unmoved) {
    double q0 = modularityCalc(sp, s, group, groupid), max = -1;
    int i, size = sp->n, maxIdx = -1;
    for (i = 0; i < size; ++i) {
        if (group != groupid[i] || unmoved[i] == 0)
            continue;
        s[i] = -s[i];
        score[i] = modularityCalc(sp, s, group, groupid) - q0;
        s[i] = -s[i];
        if (max < score[i]) {
            max = score[i];
            maxIdx = i;
        }
    }
    return maxIdx;
}*/

double findMaxImprove(double *s, const double *improve, const int *indices, const int *group, int groupSize) {
    int i, maxIdx = -1, j;
    double delta, max;
    for (i = 0; i < groupSize; ++i) {
        if (maxIdx == -1) {
            max = improve[i];
            maxIdx = i;
        }
        if (improve[i] > max) {
            maxIdx = i;
            max = improve[i];
        }
    }
    for (i = groupSize - 1; i > maxIdx; --i) {
        j = indices[i];
        s[j] = -s[j];
    }
    if (maxIdx == groupSize - 1)
        delta = 0;
    else
        delta = improve[maxIdx];
    return delta;
}

void optimize(spmat *sp, double *s, int *group, int groupSize) {
    int size = sp->n;
    int i;
    int maxIdx, idx;
    double delta;
    int *unmoved = malloc(sizeof(int) * groupSize);
    int *indices = malloc(sizeof(int) * groupSize);
    double *score = malloc(sizeof(double) * groupSize);
    double *improve = malloc(sizeof(double) * groupSize);
    double *res = malloc(sizeof(double) * size);
    double *zeroVec = malloc(sizeof(double) * size);
    if (unmoved == NULL || indices == NULL || score == NULL || improve == NULL || res == NULL || zeroVec == NULL) {
        printf("ERROR - memory allocation unsuccessful");
        exit(EXIT_FAILURE);
    }

    initOneValVec(zeroVec, groupSize, group, 0);
    do {
        resetUnmoved(unmoved, groupSize);
        int square, M = sp->M;
        multBv(sp, s, group, res, groupSize, 0);
        for (i = 0; i < groupSize; ++i) {
            idx = group[i];
            square = sp->k[idx] * sp->k[idx];
            score[i] = -2 * ((s[idx] * res[idx]) + ((double) square / M));
        }
        for (i = 0; i < groupSize; ++i) {
            maxIdx = findMaxIdx(score, groupSize, unmoved);
            s[group[maxIdx]] = -s[group[maxIdx]];
            indices[i] = group[maxIdx];
            if (i == 0)
                improve[i] = score[maxIdx];
            else
                improve[i] = improve[i - 1] + score[maxIdx];
            updateScore(sp, s, score, group, groupSize, unmoved, group[maxIdx], zeroVec, res);
            unmoved[maxIdx] = 0;
        }
        delta = findMaxImprove(s, improve, indices, group, groupSize);
    } while (IS_POSITIVE(delta));
}

int getNewGroupSize(const double *s, const int *group, int groupSize) {
    int i, counter = 0;
    for (i = 0; i < groupSize; ++i) {
        if (s[group[i]] == -1)
            counter++;
    }
    return counter;
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
        if (IS_POSITIVE(vec[g[i]]) != flag) {
            vec[g[i]] = -1;
            counter++;
        } else
            vec[g[i]] = 1;
    }
//    optimize(sp, vec, g, size);
    counter = getNewGroupSize(vec, g, size);
    delta = modularityCalc(sp, vec, g, size);
    if (!IS_POSITIVE(delta))
        return 0;
    d->Q += delta;
    /*create a new group for the -1 indexes in the +-1 vector*/
    newGroupIdx = d->numOfGroups;
    d->numOfGroups += 1;
    groups[newGroupIdx] = malloc(sizeof(int) * counter);
    tempGroup = malloc(sizeof(int) * (size - counter));
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

        } else {
            *g1Ptr = g[i];
            g1Ptr++;
        }
    }
    d->nodesforGroup[newGroupIdx] = counter;
    free(g);
    d->groups[groupIdx] = tempGroup;
    d->nodesforGroup[groupIdx] = size - counter;
    printf("old group after split size is %d\n",d->nodesforGroup[groupIdx]);
    return delta;
}

int divideToTwo(division *div, spmat *sp, int groupIdx, double *res, double *b0) {
    printf("working on group %d\n", groupIdx);
//TODO make vecF only calculated here and send it to all functions !!!!
    int size = sp->n;
    double delta;
    randomizeVec(size, b0);
    int *group = *(div->groups + groupIdx);
    int groupSize = div->nodesforGroup[groupIdx];
    double *unitVec = malloc(size * sizeof(double));
    double *vecF = malloc(size * sizeof(double));
    if (unitVec == NULL || vecF == NULL) {
        printf("ERROR - memory allocation unsuccessful");
        exit(EXIT_FAILURE);
    }
    printf("group is: ");
    printIntVector(group, groupSize);
    initOneValVec(unitVec, div->nodesforGroup[groupIdx], group, 1);
    multBv(sp, unitVec, group, vecF, groupSize, 1);
//    printf("onevec is: ");
//    printVector(unitVec, size);

//    printf("vec f is: ");
//    printVector(vecF, size);
    printf("shifting value is %f\n", sp->matShifting(sp, group, groupSize, div->vertexToGroup, groupIdx,vecF));
    powerIter(sp, b0, sp->matShifting(sp, group, groupSize, div->vertexToGroup, groupIdx, vecF), group, groupSize,res);
//    printf("HERE11\n");
    double eigen = eigenValue(sp, res, group, groupSize);
    printf("eigen value is %f\n", eigen);
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
        for (j = 0; j < d->nodesforGroup[i]; ++j) {
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
    printf("modularity is %f\n", div->Q);
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
