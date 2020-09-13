#include <stdio.h>
#include <stdlib.h>
#include "spmat.h"
#include "divider.h"
#include <float.h>
#include "utils.h"


void resetUnmoved(int *unmoved, int groupSize) {
    int i;
    for (i = 0; i < groupSize; ++i) {
        *unmoved = 1;
        unmoved++;
    }
}


void updateScore(spmat *sp, double *s, double *score, int groupSize, const int *unmoved, int k) {
    register int i, M = sp->M;
    register double sk = s[k];
//    multBv(sp, zeroVec, group, res, groupSize, 0, verticeToGroup);
    for (i = 0; i < groupSize; ++i) {
        if (*unmoved++ == 0) {
            score++;
            s++;
            continue;
        }
        if (k == i)
            *score = -*score;
        else {
//            *score -= (4 * *s * sk * *res);
            *score -= (4 * *s * sk * (sp->findAij(sp, i, k) - (double) (sp->k[i] * sp->k[k]) / M));
        }
        score++;
        s++;
    }

}

int findMaxIdx(const double *score, int groupSize, const int *unmoved) {
    register double max = -DBL_MAX;
    register int maxIdx = -1, i, flag = 0;
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

double findMaxImprove(double *s, double *improve, int *indices, int groupSize, int maxIdx) {
    register int i, j;
    double delta;
    for (i = groupSize - 1; i > maxIdx; --i) {
        j = indices[i];
        s[j] = -s[j];
    }
//    printf("max index is :%d, improve is: %f\n",maxIdx,improve[maxIdx]);
    if (maxIdx == groupSize - 1)
        delta = 0;
    else
        delta = improve[maxIdx];
    return delta;
}

void optimize(spmat *sp, double *s, int *group, int groupSize, int *verticeToGroup) {
    int size = sp->n;
    register int i, *k = sp->k, *kCopy;
    int maxIdx;
    register double delta, *sCopy, *resCopy, *scoreCopy, *improveCopy;
    register int *unmoved = malloc(sizeof(int) * groupSize);
    register int *indices = malloc(sizeof(int) * groupSize);
    register double *score = malloc(sizeof(double) * groupSize);
    register double *improve = malloc(sizeof(double) * groupSize);
    register double *res = malloc(sizeof(double) * size);
    register double maxImp, *prevImp = improve;
    register int maxImpIdx;
    if (unmoved == NULL || indices == NULL || score == NULL || improve == NULL || res == NULL) {
        printf("ERROR - memory allocation unsuccessful");
        exit(EXIT_FAILURE);
    }

    register int square, M = sp->M;
    do {
        resetUnmoved(unmoved, groupSize);
        maxImp = -DBL_MAX;
        maxImpIdx = -1;
        multBv(sp, s, group, res, groupSize, 0, verticeToGroup);
        sCopy = s;
        resCopy = res;
        scoreCopy = score;
        kCopy = k;
        for (i = 0; i < groupSize; ++i) {
            square = *kCopy * *kCopy;
            *scoreCopy = -2 * ((*sCopy * *resCopy) + ((double) square) / M);
            scoreCopy++;
            kCopy++;
            resCopy++;
            sCopy++;
        }
        for (i = 0; i < groupSize; ++i) {
            maxIdx = findMaxIdx(score, groupSize, unmoved);
//            printf("max index chosen is: %d, score : %f\n",maxIdx,score[maxIdx]);
            sCopy = s + maxIdx;
            *sCopy = -*sCopy;
            *indices = maxIdx;
            if (i == 0)
                *improve = score[maxIdx];
            else
                *improve = *(prevImp) + score[maxIdx];
//            printf("improve is :%f\n",*improve);
            if (*improve > maxImp) {
                maxImp = *improve;
                maxImpIdx = i;
            }
            updateScore(sp, s, score, groupSize, unmoved, maxIdx);
            unmoved[maxIdx] = 0;
            indices++;
            prevImp = improve;
            improve++;
        }

        indices -= i;
        improve -= i;
        delta = findMaxImprove(s, improve, indices, groupSize, maxImpIdx);
//        printf("delta is: %f\n",delta);
    } while (IS_POSITIVE(delta));
}

int getNewGroupSize(const double *s, int groupSize) {
    int i, counter = 0;
    for (i = 0; i < groupSize; ++i) {
        if (s[i] == -1)
            counter++;
    }
    return counter;
}

/* make the leading eigen-vector a +-1 vector*/
void createSVector(double *vec, int groupSize) {
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
    createSVector(vec, size);
    optimize(sp, vec, g, size, d->vertexToGroup);
    counter = getNewGroupSize(vec, size);
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

int divideToTwo(division *div, spmat *sp, networks *graphs, int groupIdx, double *res, double *b0, double *vecF,
                double shifting) {
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
    powerIter(sp, b0, sp->matShifting(sp, group, groupSize, vecF), group, groupSize, res,
              vecF, div->vertexToGroup, 1);
//    powerIter(sp, b0,shifting, group, groupSize, res,vecF, div->vertexToGroup, 1);
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
    double shifting = -1;
    int groupIdx = 0, *nodesForGroup = div->nodesforGroup, **groups = div->groups;
    double *b0 = malloc(sizeof(double) * size);
    double *res = malloc(sizeof(double) * size);
    double *unitVec = malloc(size * sizeof(double));
    double *vecF = malloc(size * sizeof(double));
    int counter = 0;
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
//            printf("counter %d\n", counter++);
            multBv(sp, unitVec, *groups, vecF, *nodesForGroup, 0, div->vertexToGroup);
            if (shifting == -1) {
                shifting = sp->matShifting(sp, div->groups[groupIdx], div->nodesforGroup[groupIdx], vecF);
            }
            delta = divideToTwo(div, sp, graphs, groupIdx, res, b0, vecF, shifting);
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
