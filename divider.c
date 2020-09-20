#include <stdio.h>
#include <stdlib.h>
#include "divider.h"
#include <float.h>

/**
@file divider.c
**Author:** Ofek Bransky & Gal Cohen
**Date:**  18.9.2020
*Summary:
 * This is the Divider C file, maintains all main methods to find the graph subgroups, and contains two structs 'Networks' & 'Division'
 * 'Networks' - > a struct containing an array of B matrices
 * 'Division' -> a struct containing all data to find the communities in the input graph, including all subgroups found
 * Functions:
 * resetUnmoved - initializes an array to all '1' values (for optimization algorithm)
 * updateScore - updates the modularity score for each vertice after moving a vertice from one group to another (for optimization algorithm)
 * findMaxImprove - finds best group split to maximize modularity (for optimization algorithm)
 * initScore - updates initial modularity score for each vertice before movements (for optimization algorithm)
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
 *  updates all vertice scores after a vertice movement, and saves the max index in maxIdx
 * @param B : the B matrix
 * @param s : the division vector
 * @param score : an array keeping the score(modularity) of each vertice in the subgroup
 * @param unmoved : an array that keeps track which vertice hasn't been moved
 * @param k : array containing vertice ranks
 * @param movedFlag : indicates which value is kept in unmoved for unmoved vertices
 * @param M : sum of vertice ranks
 * @param kPtr : an array containing all ranks for each vertice
 * @return the new index which maximize the score
 * */
int updateScore(BMat *B, double *s, double *score, const int *unmoved, int k, int movedFlag, int M, int *kPtr) {
    register int i,Kk,Aval, groupSize = B->n;
    register double sk = s[k];
    int *Bcol;
    double max = -DBL_MAX;
    int idx = -1;
    Bcol = B->getBIterator(B, k);
    Kk=*(kPtr+k);
    for (i = 0; i < groupSize; ++i) {
        if (*unmoved++ == movedFlag) {
            score++;
            s++;
            kPtr++;
            if(B->iterHasNext(B, k, i, Bcol) > 0)
                Bcol++;
            continue;
        }
        if (k == i) {
            *score = -*score;
        }
        else {
            Aval = (Bcol != NULL && *Bcol == i) ? 1 : 0;
            *score -= (4 * *s * sk * (Aval - (double) (Kk * *kPtr) / M));
        }
        if (IS_POSITIVE(*score - max)) {
            max = *score;
            idx = i;
        }
        if(B->iterHasNext(B, k, i, Bcol) > 0)
            Bcol++;
        score++;
        s++;
        kPtr++;
    }
    return idx;
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
 * @return updates all scores & returns the index with max modularity found in maxScore
 */
int initScore(int *k, int M, double *score, double *res, double *s, int groupSize) {
    register int i, square;
    int maxIdx = -1;
    double maxScore = -DBL_MAX;
    for (i = 0; i < groupSize; ++i) {
        square = *k * *k;
        *score = -2 * ((*s * *res) + ((double) square) / M);
        if (IS_POSITIVE(*score - maxScore)) {
            maxScore = *score;
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
 * the vector s (group division) will be updated with the optimized split
 *  @param B : the B matrix
 *  @param d : the division struct containing all elements needed for the function containing:
 *  * @param improve : an array keeping the improvement in modularity after each vertice movement
 *  * @param indices : an array that keeps the order of vertices moved, during the optimization
 *  * @param unmoved : an array that keeps track which vertice hasn't been moved
 *  * @param score : an array keeping the score(modularity) of each vertice in the subgroup
 *  * @param res : a vector used for calculations
 * @param s : the division vector
 * @param M : sum of vertice ranks
 * */
void optimize(division *d, BMat *B, double *s, int M) {
    register int i, groupSize = B->n;
    int maxIdx;
    int movedFlag = 1;
    double delta;
    register int *unmoved = d->unmoved;
    register int *indices = d->indices;
    register double *score = d->score;
    register double *improve = d->improve;
    register double *res = d->res;
    register double maxImp, *prevImp = improve, *sMaxIdx;
    register int maxImpIdx;
    int *k = B->getKPtr(B);
    resetUnmoved(unmoved, groupSize);

    /*runs until modularity improvement is not positive*/
    do {
        movedFlag = -movedFlag;
        maxImp = -DBL_MAX;
        maxImpIdx = -1;
        B->Bv(B, s, res);

        /*calculates initial score for all vertices*/
        maxIdx = initScore(k, M, score, res, s, groupSize);

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
            maxIdx = updateScore(B, s, score, unmoved, maxIdx, movedFlag, M, k);
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
 * splits the B matrix designating a group into two new groups, saved in networks struct 'graphs'
 * @param graphs : the networks struct
 * @param groupIdx : the group index to split
 * @param newGroupIdx : the new group index to insert the new group
 * @param s : the division vector
 * @param g : the original group
 * @param g1Size : new group 1 size
 * @param g2Size : new group 2 size
 */
void splitGraph(networks *graphs, int groupIdx, int newGroupIdx, double *s, int *g, int g1Size, int g2Size) {
    BMat **BMats = graphs->B, *currB = BMats[groupIdx], *g2B;
    g2B = currB->splitGraphB(currB, s, g, g1Size, g2Size);
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
    double eigen;
    int newGroupIdx;
    int i;
    int size = d->nodesforGroup[groupIdx];
    int **groups = d->groups;
    int *g = d->groups[groupIdx];
    int *tempGroup;
    int *g1Ptr, *g2Ptr;
    int counter = 0;

    eigen = B->eigenValue(B, vec, graphs->tmp);
    if (!IS_POSITIVE(eigen))
        initOneValVec(vec, size, 1);
    else
        createSVector(vec, size);

    /*calls modularity division optimization function*/
    optimize(d, B, vec, graphs->M);
    counter = getNewGroupSize(vec, size);

    delta = B->modularityCalc(B, vec, graphs->tmp);
    if (!IS_POSITIVE(delta))
        return 0;

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
    B->powerIter(B, b0, res);
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
 * called upon to split group into 2, updates the struct division with the subgroups which give max modularity
 * runs until the modularity isn't increased
 * @param div : the division struct storing all groups
 * @param graphs : a struct containing all sparse matrices of the graphs
 */
void findGroups(division *div, networks *graphs) {
    double delta;
    int size = graphs->n;
    BMat *B, **mats = graphs->B;
    int groupIdx = 0, *nodesForGroup = div->nodesforGroup, **groups = div->groups;
    double *b0;
    double *res;
    double *unitVec;
    double *vecF;
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
            B->Bv(B, unitVec, vecF);
            B->updateFields(B, vecF);
            delta = divideToTwo(div, B, graphs, groupIdx, res, b0);
        }
        groupIdx++;
        groups++;
        nodesForGroup++;
        mats++;
    }
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
    graphs = allocateNetworks(B->n);
    graphs->B[0] = B;
    graphs->M = B->M;
    return graphs;
}

