#include <stdio.h>
#include <stdlib.h>
#include "divider.h"
#include <float.h>
#include "bmat.h"
#include <time.h>

/**
@file divider.c
**Author:** Ofek Bransky & Gal Cohen
**Date:**  18.9.2020
*Summary:
 * This is the Divider C file, mantains all main methods to find the graph subgroups, and contains two structs 'Networks' & 'Division'
 * 'Networks' - > a struct containing an array of B matrices
 * 'Division' -> a struct containing all data to find the communities in the input graph, including all subgroups found
 * Functions:
 * resetUnmoved - initializes an array to all '1' values (for optimization algorithm)
 * updateScore - updates the modularity score for each vertice after moving a vertice from one group to another (for optimization algorithm)
 * findMaxImprove - finds best group split to maximize modularity (for optimization algorithm)
 * calcInitialScore - updates initial modularity score for each vertice before movements (for optimization algorithm)
 * optimize - finds the best group split to maximize the modularity (for optimization algorithm)
 * getNewGroupSize - counts the size of a new group after receiving a split vector
 * createSVector - creates a +1 / -1 vector to split the group
 * splitGraph - splits the groups (found in the division struct)
 * split - splits all struct, based on a found division to maximize modularity
 * divideToTwo - finds a division into two subgroups
 * freeDivision - frees the division struct
 * writeDivision - writes the found groups into a file
 * findGroups - main function, finds the best division into subgroups (communities) to maximize modularity
 * allocateDivision - allocates the division struct
 * freeNetworks - frees the networks struct
 * freeDivision- frees the division struct
 * allocateNetworks - allocates the networks struct
 * readGraph - reads the initial graph from the input file
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
 * @param sp : the sparse matrix array
 * @param s : the division vector
 * @param score : an array keeping the score(modularity) of each vertice in the subgroup
 * @param *unmoved : an array that keeps track which vertice hasn't been moved
 * @param k : array containing vertice ranks
 * @param movedFlag : indicates which value is kept in unmoved for unmoved vertices
 * @param maxIdx : a pointer to update the max index found during the function
 * @return : updates all vertice scores after a vertice movement, and saves the max index in maxIdx
 * */
void updateScore(spmat *sp, double *s, double *score, const int *unmoved, int k, int movedFlag, int *maxIdx) {
    register int i, M = sp->M, groupSize = sp->n;
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
        if (IS_POSITIVE(*score - max)) {
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

/**Reverts the group division to the optimal one, by moving back vertices that reduced the modularity
 * @param s : the division vector
 * @param improve : an array keeping the improvement in modularity after each vertice movement
 * @param indices : an array that keeps the order of vertices moved, during the optimization
 * @param groupSize : group size
 * @param maxIdx : index of vertice with max modularity reached
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

/**
 * Updates the initial score before moving any vertices
 * @param k : an array containing all ranks for each vertice
 * @param M : sum of all vertice ranks
 * @param score : an array keeping the score(modularity) of each vertice in the subgroup
 * @param res : a vector used for calculations
 * @param s : the division vector
 * @param groupSize : size of the subgroup to be split
 * @param maxScore : a pointer to the max score found
 * @return updates all scores & saves the index with max modularity found in maxScore
 */
int calcInitialScore(int *k, int M, double *score, double *res, double *s, int groupSize, double *maxScore) {
    register int i, square;
    int maxIdx = -1;

    for (i = 0; i < groupSize; ++i) {
        square = *k * *k;
        *score = -2 * ((*s * *res) + ((double) square) / M);
        if (IS_POSITIVE(*score - *maxScore)) {
            *maxScore = *score;
            maxIdx = i;
        }
        s++;
        k++;
        score++;
        res++;
    }
    return maxIdx;

}

/**Optimizes the group division to give max modularity
 * We move all the vertices from one group to another (ordered by maximum modularity)
 * We revert to the division that gives max modularity
 * if there was an improvement in modularity, we run the algorithm again
 * we stop when there is no improvement possible (max modularity is not positive)
 * @param d : the division struct containing all elements needed for the function containing:
 *  * @param improve : an array keeping the improvement in modularity after each vertice movement
 *  * @param indices : an array that keeps the order of vertices moved, during the optimization
 *  * @param unmoved : an array that keeps track which vertice hasn't been moved
 *  * @param score : an array keeping the score(modularity) of each vertice in the subgroup
 *  * @param res : a vector used for calculations
 * @param s : the division vector

 * */
void optimize(division *d, BMat *B, double *s) {
    spmat *sp;
    register int i, *k, groupSize;
    int maxIdx;
    int movedFlag = 1;
    double delta, maxScore;
    register int *unmoved = d->unmoved;
    register int *indices = d->indices;
    register double *score = d->score;
    register double *improve = d->improve;
    register double *res = d->res;
    register double maxImp, *prevImp = improve, *sMaxIdx;
    register int maxImpIdx, M;
    sp = B->sp;
    M = sp->M;
    k = sp->k;
    groupSize = sp->n;
    resetUnmoved(unmoved, groupSize);

    /*runs until modularity improvement is not positive*/
    do {
        movedFlag = -movedFlag;
        maxImp = -DBL_MAX;
        maxImpIdx = -1;
        maxScore = -DBL_MAX;
        maxIdx = -1;
        multBv(sp, s, res);

        /*calculates initial score for all vertices*/
        maxIdx = calcInitialScore(k, M, score, res, s, groupSize, &maxScore);

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

            /*updates for all vertices modularity score after vertice movement*/
            unmoved[maxIdx] = movedFlag;
            updateScore(sp, s, score, unmoved, maxIdx, movedFlag, &maxIdx);
            indices++;
            prevImp = improve;
            improve++;
        }
        indices -= i;
        improve -= i;

        /*reverts to division vector *s that gave max modularity improvement*/
        delta = findMaxImprove(s, improve, indices, groupSize, maxImpIdx);
    } while (IS_POSITIVE(delta));
}

/**Calculates the size of a subgroup divided by vector s
 * @param s : the division vector (splits the group into 2, one group with the value -1, one with 1
 * @param groupSize : size of the subgroup
 * @return counter : the number of vertices in the group asigned with the val -1
 */
int getNewGroupSize(const double *s, int groupSize) {
    int i, counter = 0;
    for (i = 0; i < groupSize; ++i) {
        if (*s++ == -1)
            counter++;
    }
    return counter;
}

/** make the leading eigen-vector a +-1 vector*/
void createSVector(double *vec, int groupSize) {
    int i, flag;
    flag = IS_POSITIVE(*vec) ? 1 : 0;
    *vec = 1;
    vec++;
    for (i = 1; i < groupSize; ++i) {
        if (IS_POSITIVE(*vec) != flag)
            *vec = -1;
        else
            *vec = 1;
        vec++;
    }
}

/**
 * splits the B matrix designating a group into two
 * @param graphs : the networks struct
 * @param groupIdx : the group index to split
 * @param newGroupIdx : the new group index to insert the new group
 * @param s : the division vector
 * @param g : the original group
 * @param g1Size : new group 1 size
 * @param g2Size : new group 2 size
 * @return two new groups, saved in networks struct 'graphs'
 */
void splitGraph(networks *graphs, int groupIdx, int newGroupIdx, double *s, int *g, int g1Size, int g2Size) {
    BMat **BMats = graphs->B, *currB = BMats[groupIdx], *g2B;
    spmat **spMats, *currSp = currB->sp;
    spMats = currSp->splitGraph(currSp, s, g, g1Size, g2Size);
    g2B = allocateB();
    currB->sp = spMats[0];
    g2B->sp = spMats[1];
    free(spMats);
    BMats[newGroupIdx] = g2B;
}

/**
 * Method takes a division vector of a group, and splits:
 * Calls Optimize - to check if there is a better division and optimize modularity
 * the designated groups into 2 subgroups
 * Calls Sparse matrix split function to create a new sparse matrix for each subgroup
 * @param d : the division vector
 * @param B : the B matrix
 * @param graphs :an array containing all the sparse matrices for each subgroup
 * @param vec : a general vector created for calculations
 * @param groupIdx : the group index that will be split
 * @return delta : modularity added to the graph after split is called
 */
double split(struct _division *d, BMat *B, networks *graphs, double *vec, int groupIdx) {

    /*initiates 2 to subgroups*/
    double delta;
    spmat *sp = B->sp;
    double *vecF = B->vecF, eigen;
    int newGroupIdx = -1;
    int i;
    int size = d->nodesforGroup[groupIdx];
    int **groups = d->groups;
    int *g = d->groups[groupIdx];
    int *tempGroup;
    int *g1Ptr, *g2Ptr;
    int counter = 0;

    eigen = eigenValue(B, vec, graphs->tmp);
    if (!IS_POSITIVE(eigen))
        initOneValVec(vec, size, 1);
    else
        createSVector(vec, size);

    /*calls modularity division optimization function*/
    optimize(d, B, vec);
    counter = getNewGroupSize(vec, size);

    delta = modularityCalc(sp, vec, graphs->tmp, vecF);
    if (!IS_POSITIVE(delta))
        return 0;
    d->Q += delta;

    /*create a new group for the -1 indexes in the +-1 vector*/
    if (counter != 0) {
        newGroupIdx = d->numOfGroups;
        d->numOfGroups += 1;

        /*splits the sparse matrix into two new matrices*/
        splitGraph(graphs, groupIdx, newGroupIdx, vec, g, size - counter, counter);
        groups[newGroupIdx] = malloc(sizeof(int) * counter);
        tempGroup = malloc(sizeof(int) * (size - counter));
        if (groups[newGroupIdx] == NULL || tempGroup == NULL) {
            error(ALLOCERROR);
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
 * @param B : the B matrix of the group to be split
 * @param graphs : an array containing the group sparse matrices
 * @param groupIdx : the group index to be split (in the networks struct)
 * @param res : vector result used for calculations
 * @param b0 : initial random vector used for power iteration
 * @return 1 if split was made, 0 if split with increase in modularity hasn't been found
 */
int divideToTwo(division *div, BMat *B, networks *graphs, int groupIdx, double *res, double *b0) {
    double delta;
    int groupSize = div->nodesforGroup[groupIdx];
    randomizeVec(b0, groupSize);
    powerIter(B, b0, res);
    /*calculates eigen value of the division vector found*/

    /*calls split function*/
    delta = split(div, B, graphs, res, groupIdx);

    if (delta == 0) {
        return 0;
    }
    return 1;
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
    if (output == NULL) {
        error(FILEOUT);
        exit(EXIT_FAILURE);
    }
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
    double shifting = -1;
    BMat *B, **mats = graphs->B;
    int groupIdx = 0, *nodesForGroup = div->nodesforGroup, **groups = div->groups;
    double *b0;
    double *res;
    double *unitVec;
    double *vecF;
    spmat *sp;
    if (graphs->M == 0) {
        error(ZERODIV);
        exit(EXIT_FAILURE);
    }
    b0 = malloc(sizeof(double) * size);
    res = malloc(sizeof(double) * size);
    unitVec = malloc(size * sizeof(double));
    vecF = malloc(size * sizeof(double));
    if (b0 == NULL || res == NULL || unitVec == NULL || vecF == NULL) {
        error(ALLOCERROR);
        exit(EXIT_FAILURE);
    }

    initOneValVec(unitVec, size, 1);
    while (groupIdx < div->numOfGroups) {
        delta = 1;
        while (delta == 1) {
            B = *mats;
            sp = B->sp;
            multBv(sp, unitVec, vecF);
            if (shifting == -1)
                shifting = sp->matShifting(sp, div->nodesforGroup[groupIdx], vecF);
            B->shifting = shifting;
            B->vecF = vecF;
            delta = divideToTwo(div, B, graphs, groupIdx, res, b0);
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
        error(ALLOCERROR);
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

/**
 * Frees the networks- a struct containing all the sparse matrices of all subgroups
 * @param graphs : the array containing the sparse matrix pointers
 * @param numOfGroups : the number of allocated groups (each has a sparse matrix)
 */
void freeNetworks(networks *graphs, int numOfGroups) {
    int i;
    BMat **Bmats = graphs->B, *BMat;
    for (i = 0; i < numOfGroups; ++i) {
        BMat = *Bmats++;
        BMat->free(BMat);
        free(BMat);
    }
    free(graphs->B);
    free(graphs->tmp);
}

/**
 * Allocates the networks struct
 * @param n : size of the graph
 * @return a pointer to the struct
 */
networks *allocateNetworks(int n) {
    networks *graphs = malloc(sizeof(networks));
    if (graphs == NULL) {
        error(ALLOCERROR);
        exit(EXIT_FAILURE);
    }
    graphs->n = n;
    graphs->B = malloc(sizeof(BMat *) * n);
    graphs->tmp = malloc(sizeof(double) * n);
    if (graphs->B == NULL || graphs->tmp == NULL) {
        error(ALLOCERROR);
        exit(EXIT_FAILURE);
    }
    graphs->free = freeNetworks;
    return graphs;
}

/**
 * A cointainer to call the type of sparse matrix to be used (List/ Array)
 * @param input : the input file
 * @return : a pointer to the networks struct
 */
networks *readGraph(FILE *input) {
    BMat *B;
    networks *graphs;
    B = readGraphB(input);
    graphs = allocateNetworks(B->sp->n);
    graphs->B[0] = B;
    graphs->M = B->sp->M;
    return graphs;
}

