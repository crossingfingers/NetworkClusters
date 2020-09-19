
#include "bmat.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/**
@file bmat.c
**Author:** Ofek Bransky & Gal Cohen
**Date:**  18.9.2020
* Summary:
* This is the B matrix C file, contains the B matrix struct, and various functions for different calculations with B
 * we maintain a B matrix for each subgroup (community) of vertices
 * Functions:
 * randomizeVec - returns a random vector
 * vecSum - sums 2 vectors
 * scalarMult - multiplies a vector with a scalar value
 * dotProd -calculates a dot product of a vector (1 int 1 double)
 * dotDoubleProd -calculates a dot product of a vector (2 double)
 * normalize- normalizes a vector
 * vecDecK - the calculation of multiplying B matrix with a vector v
 * multBv - multiplies a vector v with the B matrix
 * initOneValVec - initializes a vector, all values are equal to the value input of the function
 * vecDecFv - Decreases the vector F (sum of colums) from the vector result of Bv
 * multBRoof - multiplies a vector v with the B^ matrix (B roof matrix of a subgroup)
 * powerIter - power iteration algorithm
 * eigenValue - calculates the eigenvalue
 * modularityCalc - calculates the group's modularity
 * Bv - Outer function to multiply B matrix with vector v
 * getBIterator - outer function to get row iterator of B matrix values
 * getBValue- gets B matrix value in specific place
 * getNext - checks if theres a next value in the row
 * getKPtr - returns pointer to K (vertice rank) array of group
 * updateFields - updates B struct fields
 * splitGraphB - splits the B struct into two new B structs (one for each new subgroup)
 * freeB - frees a B matrix struct from memory
 * readGraphB - reads the graph into a initial B matrix
 * allocateB - allocates a B struct
*/

/**
 * inserts random variables into an initialized vector
 * @param vec : the initialized vector
 * @param groupSize : the number of values inserted, inserted into the first indexes
 */
void randomizeVec(double *vec, int groupSize) {
    int i;
    for (i = 0; i < groupSize; i++) {
        *vec = rand();
        vec++;
    }
}

/** gets two vectors and returns the sum of the vector with b0*shifting
 * @param vec : vec1 (and the result)
 * @param b0 : vec 2
 * @param shifting : the value from matrix shifting)
 * @param n : the group size
 * */
void vecSum(double *vec, const double *b0, double shifting, int n) {
    int i;
    for (i = 0; i < n; ++i) {
        *vec += *b0++ * shifting;
        vec++;
    }
}

/**multiplies a scalar value with a vector, updates the vector values received
 * @param vec : the vector
 * @param x : the scalar
 * @param n : the group size
 * */
void scalarMult(double *vec, double x, int n) {
    int i;
    for (i = 0; i < n; ++i) {
        *vec++ *= x;
    }
}

/**calculates the dot product of integer vector with double vector (size n), returns the dot product */
double dotProd(const int *vec1, const double *vec2, int n) {
    register double res;
    int i;
    res = 0;
    for (i = 0; i < n; ++i) {
        res += *vec1++ * *vec2++;
    }
    return res;
}

/**calculates the dot product of two double vectors (size n), returns the dot product*/
double dotDoubleProd(const double *vec1, const double *vec2, int n) {
    register double res = 0;
    int i;
    for (i = 0; i < n; ++i) {
        res += *vec1++ * *vec2++;
    }
    return res;
}

/**normalizes a vector (group size n)*/
void normalize(double *vec, int groupSize) {
    double res = dotDoubleProd(vec, vec, groupSize);
    res = sqrt(res);
    scalarMult(vec, 1 / res, groupSize);
}

/**function used for the calculation of multiplying B matrix with a vector v by spliting B :
 * to A (sparse matrix) and KiKj matrix (rank matrix), the result is saved in vec1
 * @param vec1 : the multiplication vector
 * @param sp : the sparse matrix
 * @param n : the size of the group to be multiplied
 * @param dotM : the dot product
 */
void vecDecK(double *vec1, spmat *sp, int n, double dotM) {
    int i, *k = sp->k;
    for (i = 0; i < n; ++i) {
        *vec1 -= *k++ * dotM;
        vec1++;
    }
}

/**
 * multiplies a vector and the modularity matrix B
 * @param sp : the sparse matrix
 * @param vec : the vector to be multiplied by
 * @param res : the vector result of the multiplication
 * @return : 1 on success
 */
double multBv(spmat *sp, double *vec, double *res) {
    double dot;
    int groupSize = sp->n;
    dot = dotProd(sp->k, vec, groupSize);
    sp->mult(sp, vec, res);
    vecDecK(res, sp, groupSize, (double) dot / sp->M);
    return 1;
}

/**
 * initializes a vector with all values identical
 * @param unitVec : the vector to be initialized
 * @param n : the size of the vector
 * @param val : the value to be inserted
 */
void initOneValVec(double *unitVec, int n, int val) {
    int i;
    for (i = 0; i < n; ++i) {
        *unitVec = val;
        unitVec++;
    }
}

/**
 * Decreases the vector F (sum of colums) from the vector result of Bv
 * @param res : the result of the calculation
 * @param vecF : the F vector
 * @param v : the result of BxV
 * @param n :size of the group
 */
void vecDecFv(double *res, double *vecF, double *v, int n) {
    int i;
    for (i = 0; i < n; ++i) {
        *res -= *vecF++ * *v++;
        res++;
    }
}

/**
 * multiplies a vector and the modularity matrix B-roof, of a specific subgroup
 * @param sp : the sparse matrix
 * @param vec : the vector to be multiplied by
 * @param res : the vector result of the multiplication
 * @param vecF : an array containing the sum of values for each column
 * @return returns 1 on success
 */
double
multBRoof(spmat *sp, double *vec, double *res, double *vecF) {
    int size = sp->n;
    double total;
    total = multBv(sp, vec, res);
    vecDecFv(res, vecF, vec, size);
    return total;
}

/**
 * Power iteration method, multiplying a random vector with the matrix B, to find max eigenvector,
 * saves the power iteration result into the result vector
 * @param B : the B matrix
 * @param b0  : initial random vector
 * @param result : the vector result of the power iteration
 */
void powerIter(BMat *B, double *b0, double *result) {
    register int flag = 1, i;
    spmat *sp = B->sp;
    int groupSize = sp->n;
    int counter = 0;
    double MAX_ITERS = 5000 * groupSize + 80000;
    double total = 0, shifting = B->shifting, *vecF = B->vecF;
    while (flag == 1 && counter < MAX_ITERS) {
        flag = 0;
        total += multBRoof(sp, b0, result, vecF);
        vecSum(result, b0, shifting, groupSize);
        normalize(result, groupSize);
        for (i = 0; i < groupSize; i++) {
            if (IS_POSITIVE(fabs(result[i] - b0[i])))
                flag = 1;
            b0[i] = result[i];
        }
        counter++;
    }
    if (counter >= MAX_ITERS) {
        error(PIERROR);
        exit(EXIT_FAILURE);
    }
}

/**
 * Calculates the eigenvalue of a vector
 * @param B : B matrix
 * @param vec : the eigenvector
 * @param tmp : a temporary vector to assist calculation
 * @return : the eigenvalue
 */
double eigenValue(BMat *B, double *vec, double *tmp) {
    double res;
    spmat *sp = B->sp;
    double *vecF = B->vecF;
    int groupSize = sp->n;
    multBRoof(sp, vec, tmp, vecF);
    res = dotDoubleProd(tmp, vec, groupSize);
    return res;
}

/**
 * calculates the modularity of a division vector for a subgroup
 * @param B : B matrix
 * @param vec : the eigenvector
 * @param tmp : a temporary vector to assist calculation
 * @param vecF : an array containing the sum of values for each column
 * @return the modularity of the division
 */
double modularityCalc(BMat *B, double *vec, double *tmp) {
    spmat *sp = B->sp;
    double res = 0, *vecF = B->vecF;
    int groupSize = sp->n;
    multBRoof(sp, vec, tmp, vecF);
    res = dotDoubleProd(tmp, vec, groupSize);
    return res / 2;
}

/**
 * Outer function to be used outside of Bmat module, calls Bv multiplication function
 * @param B : the B matrix
 * @param vec : the vector to be multiplied
 * @param res : the result, saved
 */
void Bv(BMat *B, double *vec, double *res){
    multBv(B->sp, vec, res);
}

/**
 * Outer function to be used outside of Bmat module,
 * returns an iterator to iterate upon row values of a B matrix
 * @param B : B matrix
 * @param i : row index
 * @return  : pointer to first row value
 */
int *getBIterator(BMat *B, int i){
    spmat *sp = B->sp;
    int *Acol = sp->getARowIterator(sp, i);
    return Acol;
}

/**
 * Outer function to be used outside, returns the value B from a B matrix, when iterated upon
 * used mostly for Optimization
 * @param B : B matrix
 * @param i : row index
 * @param j : the moved vertice
 * @param ptr :pointer to row values
 * @return the value of B in indexes i,j
 */
double getBValue(BMat *B, int i, int j, const int *ptr){
    spmat *sp = B->sp;
    int M = sp->M;
    int Aval = (ptr != NULL && *ptr == j) ? 1 : 0;
    return (Aval - (double) (sp->k[i] * sp->k[j]) / M);
}

/**
 * iterates upon a row of B matrix
 * @param B : B matrix
 * @param i : row index
 * @param j  : the last index in the row
 * @param ptr : pointer to first row value
 * @return the next value in the specific row
 */
int *getNext(BMat *B, int i, int j, int *ptr){
    spmat *sp = B->sp;
    if (ptr != NULL && *ptr <= j && sp->hasNextARow(sp, i, ptr)) {
        ptr++;
    }
    return ptr;
}

/**
 * Outer function, returns pointer to K array (vertice ranks)
 * @param B : the B matrix
 * @return pointer to K array
 */
int *getKPtr(BMat *B){
    return B->sp->k;
}

/**
 * updates B struct fields
 * @param B : B matrix
 * @param vecF : vector F (sum of columns)
 * @return updated B struct (with vector F and shifting value)
 */
void updateFields(BMat *B, double *vecF){
    if(B->shifting == -1)
        B->shifting = B->sp->matShifting(B->sp, vecF);
    B->vecF = vecF;
}

/**
 * Splits a B struct into two subgroups, including inner implementation of A sparse matrix
 * @param currB :current B struct
 * @param s : division vector
 * @param g : current group
 * @param g1Size : new group 1 size
 * @param g2Size : new group 2 size
 * @return : a pointer to the second group, first group is saved in original pointer
 */
BMat *splitGraphB(BMat *currB, double *s, int *g, int g1Size, int g2Size){
    BMat *g2B;
    spmat **spMats, *currSp = currB->sp;
    spMats = currSp->splitGraph(currSp, s, g, g1Size, g2Size);
    g2B = allocateB();
    currB->sp = spMats[0];
    currB->n = currB->sp->n;
    g2B->sp = spMats[1];
    g2B->n = g2B->sp->n;
    g2B->shifting = currB->shifting;
    free(spMats);
    return g2B;
}

/**
 * frees the B matrix allocated memory
 * @param B matrix
 */
void freeB(BMat *B) {
    B->sp->free(B->sp);
    free(B->sp);
}

/**
 * allocates a new B matrix
 * @return a pointer to the matrix
 */
BMat *allocateB() {
    BMat *B = malloc(sizeof(BMat));
    if (B == NULL) {
        error(ALLOCERROR);
        exit(EXIT_FAILURE);
    }
    B->free = freeB;
    B->Bv = Bv;
    B->modularityCalc = modularityCalc;
    B->eigenValue = eigenValue;
    B->powerIter = powerIter;
    B->getBIterator = getBIterator;
    B->getNext = getNext;
    B->getBValue = getBValue;
    B->getKPtr = getKPtr;
    B->updateFields = updateFields;
    B->splitGraphB = splitGraphB;
    return B;
}

/**
 * reads into a initial B struct the input graph from the binary file
 * @param input : the input file pointer
 * @return a pointer to the bmat struct
 */
BMat *readGraphB(FILE *input) {
    spmat *sp;
    BMat *B = allocateB();
    sp = readGraphA(input);
    B->sp = sp;
    B->n = sp->n;
    B->M = sp->M;
    B->shifting = -1;
    return B;
}