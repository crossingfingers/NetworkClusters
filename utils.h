//
// Created by gal21 on 05/09/2020.
//

#ifndef CPROJECT_UTILS_H
#define CPROJECT_UTILS_H

#include "spmat.h"

#define IS_POSITIVE(X) ((X) > 0.00001)

/**
@file utils.h
**Author:** Ofek Bransky & Gal Cohen
**Date:**  18.9.2020
## This is the utility header file, contains various functions for different calculations in the program
*/

/**
 * inserts random variables into an initialized vector
 * @param size : the size of the vector
 * @param vec : the initialized vector
 * @param groupSize : the number of values inserted, inserted into the first indexes
 * @param group : the subgroup array containing the vertices of the group
 */
void randomizeVec(int size, double *vec, int groupSize, int *group);

//TODO- remove debug
/**
 * multiplies a vector and the modularity matrix B
 * @param sp : the sparse matrix
 * @param vec : the vector to be multiplied by
 * @param res : the vector result of the multiplication
 * @param groupSize : the number of values inserted, inserted into the first indexes
 * @param group : the subgroup array containing the vertices of the group
 * @param debug
 * @return : the time taken to finish the method
 */
double multBv(spmat *sp, double *vec, const int *group, double *res, int groupSize, int debug);


/**
 * initializes a vector with all values identical
 * @param unitVec : the vector to be initialized
 * @param n : the size of the vector
 * @param group : the subgroup array containing the vertices of the group
 * @param val : the value to be inserted
 */
void initOneValVec(double *unitVec, int n, const int *group, int val);

//TODO- remove debug
/**
 * multiplies a vector and the modularity matrix B-roof, of a specific subgroup
 * @return
 * @param sp : the sparse matrix
 * @param vec : the vector to be multiplied by
 * @param res : the vector result of the multiplication
 * @param groupSize : the number of values inserted, inserted into the first indexes
 * @param group : the subgroup array containing the vertices of the group
 * @param vecF : an array containing the sum of values for each column
 * @param debug
 * @return : the time taken to finish the method
 */
double multBRoof(spmat *sp, double *vec, const int *group, int groupSize, double *res, double *vecF, int debug);

//TODO- remove debug
/**
 * Power iteration method, multiplying a random vector with the matrix B, to find max eigenvector
 * @param sp : sparse matrix
 * @param b0  : initial random vector
 * @param shifting  : shifting value to shift the matrix & find max positive eigenvalue
 * @param result : the vector result of the power iteration
 * @param groupSize : the number of values inserted, inserted into the first indexes
 * @param group : the subgroup array containing the vertices of the group
 * @param vecF : an array containing the sum of values for each column
 * @param debug
 */
void
powerIter(spmat *sp, double *b0, double shifting, int *group, int groupSize, double *result, double *vecF, int debug);

/**
 * Calculates the eigenvalue of a vector
 * @param sp : sparse matrix
 * @param vec : the eigenvector
 * @param group : the subgroup array containing the vertices of the group
 * @param groupSize : the number of values inserted, inserted into the first indexes
 * @param vecF : an array containing the sum of values for each column
 * @return : the eigenvalue
 */
double eigenValue(spmat *sp, double *vec, const int *group, int groupSize, double *vecF);

/**
 * calculates the modularity of a division vector for a subgroup
 * @param sp : sparse matrix
 * @param vec : the eigenvector
 * @param group : the subgroup array containing the vertices of the group
 * @param groupSize : the number of values inserted, inserted into the first indexes
 * @param vecF : an array containing the sum of values for each column
 * @return the modularity of the division
 */
double modularityCalc(spmat *sp, double *vec, int *group, int groupSize, double *vecF);

void printVector(double *vec, int n, const int *group);

void printIntVector(int *vec, int n);

#endif //CPROJECT_UTILS_H
