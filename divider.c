#include <stdio.h>
#include <stdlib.h>
#include "spmat.h"
#include "divider.h"
#include <float.h>
#include "utils.h"


void resetUnmoved(int *unmoved, int groupSize) {
    int i;
    for (i = 0; i < groupSize; ++i) {
//        unmoved[i] = 1;
        *unmoved = 1;
        unmoved++;
    }
}


void updateScore(spmat *sp, double *s, double *score, const int *group, int groupSize, const int *unmoved, int k,
                 double *zeroVec, double *res, int *verticeToGroup) {
    int i;
    zeroVec[k] = 1;
    double *sCopy = s;
    multBv(sp, zeroVec, group, res, groupSize, 0, verticeToGroup);
    for (i = 0; i < groupSize; ++i) {
        if (*unmoved++ == 0) {
            score++;
            sCopy++;
            res++;
            continue;
        }
        if (k == i)
            *score = -(*score);
        else {
            *score -= (4 * (*sCopy) * s[k] * (*res));
        }
        score++;
        sCopy++;
        res++;
    }
    zeroVec[k] = 0;
}

int findMaxIdx(const double *score, int groupSize, const int *unmoved) {
    double max = -DBL_MAX;
    int maxIdx = -1, i, flag = 0;
    for (i = 0; i < groupSize; ++i) {
        if (*unmoved == 0) {
            unmoved++;
            score++;
            continue;
        }
        if (flag == 0 || max < *score) {
            flag = 1;
            max = *score;
            maxIdx = i;
        }
        score++;
        unmoved++;
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

double findMaxImprove(double *s, double *improve, int *indices, const int *group, int groupSize) {

    int i, maxIdx = -1, j;
    double delta, max;
    double *improveCopy=improve;
    int *indicesCopy=indices;
    int groupSizeMinusOne=groupSize-1;
    for (i = 0; i < groupSize; ++i) {
        if (maxIdx == -1) {
            max = *improveCopy;
            maxIdx = i;
        }
        if (*improveCopy > max) {
            maxIdx = i;
            max = *improveCopy;
        }
        improveCopy++;
    }
    indicesCopy+=groupSizeMinusOne;
    for (i = groupSizeMinusOne; i > maxIdx; --i) {
        j = *indicesCopy;
        s[j] = -s[j];
        indicesCopy--;
    }
    if (maxIdx == groupSizeMinusOne)
        delta = 0;
    else
        delta = improve[maxIdx];
    return delta;
}

void optimize(spmat *sp, double *s, int *group, int groupSize, int *verticeToGroup) {
    int size = sp->n;
    int i;
    int maxIdx;
    double delta;
    int *indicesCopy, *kCopy;
    double *scoreCopy, *improveCopy, *resCopy, *sCopy;
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
        indicesCopy = indices;
        kCopy = sp->k;
        scoreCopy = score;
        improveCopy = improve;
        indicesCopy = indices;
        resCopy = res;
        sCopy = s;
        int square, M = sp->M;

        multBv(sp, s, group, res, groupSize, 0, verticeToGroup);


        for (i = 0; i < groupSize; ++i) {
            square = (*kCopy) * (*kCopy);
            *scoreCopy = -2 * (((*sCopy) * (*resCopy)) + ((double) square / M));
            kCopy++;
            scoreCopy++;
            sCopy++;
            resCopy++;
        }
        for (i = 0; i < groupSize; ++i) {
            maxIdx = findMaxIdx(score, groupSize, unmoved);
            s[maxIdx] = -s[maxIdx];
            *indicesCopy = maxIdx;
            if (i == 0)
                *improveCopy = score[maxIdx];
            else
                *improveCopy = *(improveCopy - 1) + score[maxIdx];
            updateScore(sp, s, score, group, groupSize, unmoved, maxIdx, zeroVec, res, verticeToGroup);
            unmoved[maxIdx] = 0;
            indicesCopy++;
            improveCopy++;
        }
        delta = findMaxImprove(s, improve, indices, group, groupSize);
        printf("delta is: %f\n",delta);
    } while (IS_POSITIVE(delta));
}

int getNewGroupSize(const double *s, const int *group, int groupSize) {
    int i, counter = 0;
    for (i = 0; i < groupSize; ++i) {
        if (s[i] == -1)
            counter++;
    }
    return counter;
}

/* make the leading eigen-vector a +-1 vector*/
void createSVector(double *vec, int *g, int groupSize) {
    int i, flag = -1;
    for (i = 0; i < groupSize; ++i) {
        if (flag == -1)
            flag = IS_POSITIVE(vec[i]) ? 1 : 0;
        if (IS_POSITIVE(vec[i]) != flag) {
            vec[i] = -1;
//            counter++;
        } else
            vec[i] = 1;
    }
}

double split(struct _division *d, spmat *sp, networks *graphs, double *vec, int groupIdx, double *vecF) {
    double delta;
    int newGroupIdx = -1;
    int i;
    int size = d->nodesforGroup[groupIdx];
    int **groups = d->groups;
    int *g = d->groups[groupIdx];
    int *tempGroup;
    int *g1Ptr, *g2Ptr;
    int counter = 0;
    createSVector(vec, g, size);
    optimize(sp, vec, g, size, d->vertexToGroup);
    counter = getNewGroupSize(vec, g, size);
    delta = modularityCalc(sp, vec, g, size, vecF, d->vertexToGroup);
    if (!IS_POSITIVE(delta))
        return 0;
    d->Q += delta;
    /*create a new group for the -1 indexes in the +-1 vector*/
    if (counter != 0) {
        newGroupIdx = d->numOfGroups;
        d->numOfGroups += 1;
        sp->splitGraph(graphs, groupIdx, newGroupIdx, vec, g, size, size - counter, counter);
        groups[newGroupIdx] = malloc(sizeof(int) * counter);
        tempGroup = malloc(sizeof(int) * (size - counter));
        if (groups[newGroupIdx] == NULL || tempGroup == NULL) {
            printf("ERROR - memory allocation unsuccessful");
            exit(EXIT_FAILURE);
        }
        g2Ptr = groups[newGroupIdx];
        g1Ptr = tempGroup;
        for (i = 0; i < size; ++i) {
            if (vec[i] == -1) {
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
    }
    return delta;
}

int divideToTwo(division *div, spmat *sp, networks *graphs, int groupIdx, double *res, double *b0, double *vecF) {
//    printf("working on group %d\n", groupIdx);
//TODO make vecF only calculated here and send it to all functions !!!!
    int size = sp->n;
    double delta;
    int *group = *(div->groups + groupIdx);
    int groupSize = div->nodesforGroup[groupIdx];
    randomizeVec(size, b0, groupSize, group);
//    printf("group is: ");
//    printIntVector(group, groupSize);
//    printf("onevec is: ");
//    printVector(unitVec, size);

//    printf("vec f is: ");
//    printVector(vecF, size);
//    printf("shifting value is %f\n", sp->matShifting(sp, group, groupSize, div->vertexToGroup, groupIdx,vecF));
    powerIter(sp, b0, sp->matShifting(sp, group, groupSize, div->vertexToGroup, groupIdx, vecF), group, groupSize, res,
              vecF, div->vertexToGroup, 1);
//    printf("HERE11\n");
//    printVector(res, groupSize, group);
    double eigen = eigenValue(sp, res, group, groupSize, vecF, div->vertexToGroup);
//    printf("eigen value is %f\n", eigen);
    if (!IS_POSITIVE(eigen)) {
//        free(vecF);
//        free(unitVec);
        return 0;
    }
    delta = split(div, sp, graphs, res, groupIdx, vecF);
    if (delta == 0) {
//        free(vecF);
//        free(unitVec);
        return 0;
    }
//    div->printGroups(div);
//    free(vecF);
//    free(unitVec);
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


void findGroups(division *div, networks *graphs) {
    double delta;
    int size = graphs->n;
    spmat **mats = graphs->A;
    spmat *sp = *mats;
    int groupIdx = 0, *nodesForGroup = div->nodesforGroup, **groups = div->groups;
    double *b0 = malloc(sizeof(double) * size);
    double *res = malloc(sizeof(double) * size);
    double *unitVec = malloc(size * sizeof(double));
    double *vecF = malloc(size * sizeof(double));
    if (b0 == NULL || res == NULL || unitVec == NULL || vecF == NULL) {
        printf("ERROR - memory allocation unsuccessful");
        exit(EXIT_FAILURE);
    }
    initOneValVec(unitVec, size, groups[0], 1);
    while (groupIdx < div->numOfGroups) {
        delta = 1;
        while (delta == 1) {
            sp = *mats;
//            printVector(vecF, *nodesForGroup, *groups);
            multBv(sp, unitVec, *groups, vecF, *nodesForGroup, 0, div->vertexToGroup);
            delta = divideToTwo(div, sp, graphs, groupIdx, res, b0, vecF);
        }
        groupIdx++;
        groups++;
        nodesForGroup++;
        mats++;
    }
    printf("modularity is %f\n", div->Q);
    free(b0);
    free(res);
    free(vecF);
    free(unitVec);
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
