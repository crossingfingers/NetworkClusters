#include <stdio.h>
#include <stdlib.h>
#include "spmat.h"
#include "divider.h"
#include <float.h>
#include "utils.h"

/**
@file divider.h
**Author:** Ofek Bransky & Gal Cohen
**Date:**  18.9.2020
## This is the Divider C file, maintains all main methods to find the graph subgroups
*/


/**resets Unmoved, an array that indicates if a vertice has been moved during optimization*/
void resetUnmoved(int *unmoved, int groupSize) {
    int i;
    for (i = 0; i < groupSize; ++i) {
        *unmoved = 1;
        unmoved++;
    }
}


/**calculates the difference in modularity if a certain vertice is moved
 * @param *score : an array keeping the score(modularity) of each vertice in the subgroup
 * @param M : the M value of the graph (sum of edges divided by two
 * @param *unmoved : an array that keeps track which vertice hasn't been moved
 * */
void
updateScore(spmat *sp, double *s, double *score, int groupSize, const int *unmoved, int k, int movedFlag, int *maxIdx) {
    register int i, M = sp->M;
    register double sk = s[k];
    int *Acol, Aval;
    double max = -DBL_MAX;
    int idx = -1;
    Acol = sp->getARowIterator(sp, k);
    for (i = 0; i < groupSize; ++i) {
        if (*unmoved++ == movedFlag) {
            score++;
            s++;
            if (Acol != NULL && *Acol <= i && sp->hasNextARow(sp, k, Acol))
                Acol++;
            continue;
        }
        Aval = (Acol != NULL && *Acol == i) ? 1 : 0;

        if (k == i)
            *score = -*score;
        else {
            *score -= (4 * *s * sk * (Aval - (double) (sp->k[i] * sp->k[k]) / M));
        }
        if (*score > max) {
            max = *score;
            idx = i;
        }
        if (Acol != NULL && *Acol <= i && sp->hasNextARow(sp, k, Acol)) {
            Acol++;
        }
        score++;
        s++;
    }
    *maxIdx = idx;
}

/**Finds the index of the vertice that if moved, will add maximum modularity
 * @param *score : an array keeping the score(modularity) of each vertice in the subgroup
 * @param *unmoved : an array that keeps track which vertice hasn't been moved
 * @return maxIdx : the index with the maximum modularity score
 * */
int findMaxIdx(const double *score, int groupSize, const int *unmoved, int movedFlag) {
    register double max = -DBL_MAX;
    register int maxIdx = -1, i, flag = 0;
    for (i = 0; i < groupSize; ++i) {
        if (*unmoved == movedFlag) {
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
    return maxIdx;
}

/**Reverts the group division to the optimal one, by moving back vertices that reduced the modularity
 * @param *improve : an array keeping the improvement in modularity after each vertice movement
 * @param *indices : keeps the order of vertices moved, during the optimization
 * @return delta : returns modularity improvement
 * */
double findMaxImprove(double *s, const double *improve, const int *indices, int groupSize, int maxIdx) {
    register int i, j;
    double delta;
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

/**Optimizes the group division to give max modularity
 * We move all the vertices from one group to another (ordered by maximum modularity)
 * We revert to the division that gives max modularity
 * if there was an improvement in modularity, we run the algorithm again
 * we stop when there is no improvement possible (max modularity is not positive)
 * @param *improve : an array keeping the improvement in modularity after each vertice movement
 * @param *indices : keeps the order of vertices moved, during the optimization
 * @param *unmoved : an array that keeps track which vertice hasn't been moved
 * @param *score : an array keeping the score(modularity) of each vertice in the subgroup
 * @param *res : a vector used for calculations
 * */
void optimize(struct _division *d, spmat *sp, double *s, int *group, int groupSize) {
    register int i, *k = sp->k;
    int maxIdx;
    int movedFlag = 1;
    register double delta;
    register int *unmoved = d->unmoved;
    register int *indices = d->indices;
    register double *score = d->score;
    register double *improve = d->improve;
    register double *res = d->res;

    if ((unmoved == NULL) || (indices == NULL) || (score == NULL) || (improve == NULL) || (res == NULL)) {
        printf("ERROR - memory allocation unsuccessful");
        exit(EXIT_FAILURE);
    }
    register double maxImp, *prevImp = improve, maxScore, *sMaxIdx;
    register int maxImpIdx;
    register int square, M = sp->M;
    resetUnmoved(unmoved, groupSize);

    /*runs until modularity improvement is not positive*/
    do {
        movedFlag = -movedFlag;
        maxImp = -DBL_MAX;
        maxImpIdx = -1;


        maxScore = -DBL_MAX;
        maxIdx = -1;
        multBv(sp, s, group, res, groupSize, 0);
        /*calculates initial score for all vertices*/
        for (i = 0; i < groupSize; ++i) {
            square = *k * *k;
            *score = -2 * ((*s * *res) + ((double) square) / M);
            if (*score >= maxScore) {
                maxScore = *score;
                maxIdx = i;
            }
            s++;
            k++;
            score++;
            res++;
        }
        s -= groupSize;
        k -= groupSize;
        res -= groupSize;
        score -= groupSize;

        /*runs until all vertices have been moved once (unmoved array is empty)*/
        for (i = 0; i < groupSize; ++i) {

            /*moves the vertice*/
            sMaxIdx = s + maxIdx;
            *sMaxIdx = -*sMaxIdx;

            /*updates movement order in *indices array*/
            *indices = maxIdx;
            /*updates modularity improvement in *improve array*/
            if (i == 0)
                *improve = score[maxIdx];
            else
                *improve = *(prevImp) + score[maxIdx];
            if (*improve > maxImp) {
                maxImp = *improve;
                maxImpIdx = i;
            }

            //TODO- check with gal about this part
            /*updates for all vertices modularity score after vertice movement*/
            updateScore(sp, s, score, groupSize, unmoved, maxIdx, movedFlag, &maxIdx);
            unmoved[maxIdx] = movedFlag;

            indices++;
            prevImp = improve;
            improve++;
        }

        indices -= i;
        improve -= i;
        /*reverts to division vector *s that gave max modularity improvement*/
        delta = findMaxImprove(s, improve, indices, groupSize, maxImpIdx);
//        printf("delta is: %f\n",delta);
    } while (IS_POSITIVE(delta));
}

/**Calculates the size of a subgroup divided by  vector s
 * @param s the division vector (splits the group into 2, one group with the value -1, one with 1
 * @param groupSize
 * @return counter the number of vertices in the group asigned with the val -1
 */
int getNewGroupSize(const double *s, int groupSize) {
    int i, counter = 0;
    for (i = 0; i < groupSize; ++i) {
        if (s[i] == -1)
            counter++;
    }
    return counter;
}

/** make the leading eigen-vector a +-1 vector*/
void createSVector(double *vec, int groupSize) {
    int i, flag = -1;
    for (i = 0; i < groupSize; ++i) {
        if (flag == -1)
            flag = IS_POSITIVE(vec[i]) ? 1 : 0;
        if (IS_POSITIVE(vec[i]) != flag) {
            vec[i] = -1;
        } else
            vec[i] = 1;
    }
}

/**
 * Method takes a division vector of a group, and splits:
 * Calls Optimize - to check if there is a better division and optimize modularity
 * the designated groups into 2 subgroups
 * Calls Sparse matrix split function to create a new sparse matrix for each subgroup
 * @param d : the division vector
 * @param sp : the Sparse matrix
 * @param graphs :an array containing all the sparse matrices for each subgroup
 * @param vec : a general vector created for calculations
 * @param groupIdx : the group index that will be split
 * @param vecF : vector containing column sums of the sparse matrix, useful for calculations
 * @return delta : modularity added to the graph after split is called
 */
double split(struct _division *d, spmat *sp, networks *graphs, double *vec, int groupIdx, double *vecF) {

    /*initiates 2 to subgroups*/
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
    /*calls modularity division optimization function*/
    optimize(d, sp, vec, g, size);
    counter = getNewGroupSize(vec, size);

    delta = modularityCalc(sp, vec, g, size, vecF);
    if (!IS_POSITIVE(delta))
        return 0;
    d->Q += delta;

    /*create a new group for the -1 indexes in the +-1 vector*/
    if (counter != 0) {
        newGroupIdx = d->numOfGroups;
        d->numOfGroups += 1;

        /*splits the sparse matrix into two new matrices*/
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

/**Finds an initial group division into two subgroups using power iteration
 * @param div : the division struct storing all the subdivisions
 * @param sp : the sparse matrix of the group to be split
 * @param graphs : an array containing the group sparse matrices
 * @param res : vector result used for calculations
 * @param b0 : initial random vector used for power iteration
 * @param vecF : sum of columns from B matrice
 * @param shifting : shifting value to get max positive eigenvalue
 * @return 1 if split was made, 0 if split with increase in modularity hasn't been found
 */
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
//    powerIter(sp, b0, sp->matShifting(sp, group, groupSize, vecF), group, groupSize, res,
//              vecF, 1);
    powerIter(sp, b0, shifting, group, groupSize, res, vecF, 1);
//    printf("HERE11\n");
//    printVector(res, groupSize, group);
    /*calculates eigen value of the division vector found*/
    double eigen = eigenValue(sp, res, group, groupSize, vecF);
//    printf("eigen value is %f\n", eigen);
    if (!IS_POSITIVE(eigen)) {
//        free(vecF);
//        free(unitVec);
        return 0;
    }
    /*calls split function*/
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

/** prints all subgroups of the graph
 * @param d : the division to be printed
 */
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

/**
 * frees all allocated memory in the division struct
 * @param d the division struct to be freed
 */
void freeDivision(division *d) {
    int i;
    for (i = 0; i < d->numOfGroups; ++i) {
        free(d->groups[i]);
    }
    free(d->groups);
    free(d->nodesforGroup);
    free(d->res);
    free(d->improve);
    free(d->score);
    free(d->indices);
    free(d->unmoved);
}

/**
 * Writes all subgroups found into binary file
 * @param div : struct containing all subgroups
 * @param output : the file to be written to
 */
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

/**
 * The main iterative method to find subgroups the the graph
 * called upon to split group into 2
 * runs until the modularity isn't increased
 * @param div : the division struct storing all groups
 * @param graphs : a struct containing all sparse matrices of the graphs
 * @return updates the struct division with the subgroups which give max modularity
 */
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
            multBv(sp, unitVec, *groups, vecF, *nodesForGroup, 0);
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

/** Allocates a division structure that contains all the subgroups
 * @param n : number of vertices in the graph
 * @return a pointer to the struct
 */
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
    d->unmoved = malloc(sizeof(int) * n);
    d->indices = malloc(sizeof(int) * n);
    d->score = malloc(sizeof(double) * n);
    d->improve = malloc(sizeof(double) * n);
    d->res = malloc(sizeof(double) * n);
    d->writeDivision = writeDivision;
    d->findGroups = findGroups;
    d->Q = 0;
    d->groups[0] = malloc(sizeof(int) * n);
    if (d->groups == NULL || d->nodesforGroup == NULL || d->groups[0] == NULL) {
        printf("ERROR - memory allocation unsuccessful");
        exit(EXIT_FAILURE);
    }
    d->numOfGroups = 1;
    ptr = d->groups[0];
    for (i = 0; i < n; ++i) {
        *ptr = i;
        ptr++;
    }
    d->nodesforGroup[0] = n;
    return d;
}
