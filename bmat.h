
#ifndef CPROJECT_UTILS_H
#define CPROJECT_UTILS_H

#include "spmat.h"

/**
@file utils.h
**Author:** Ofek Bransky & Gal Cohen
**Date:**  18.9.2020
## This is the utility header file, contains various functions for different calculations in the program
*/




typedef struct _bmat{
    spmat *sp;
    double *vecF;
    double shifting;
    void (*free)(struct _bmat *B);
}BMat;

BMat *readGraphB(FILE *input);
BMat *allocateB();


/**
 * inserts random variables into an initialized vector
 * @param vec : the initialized vector
 * @param groupSize : the number of values inserted, inserted into the first indexes
 */
void randomizeVec( double *vec, int groupSize);

/**
 * multiplies a vector and the modularity matrix B
 * @param sp : the sparse matrix
 * @param vec : the vector to be multiplied by
 * @param res : the vector result of the multiplication
 * @param groupSize : the number of values inserted, inserted into the first indexes
 * @return : the time taken to finish the method
 */
double multBv(spmat *sp, double *vec, double *res);


/**
 * initializes a vector with all values identical
 * @param unitVec : the vector to be initialized
 * @param n : the size of the vector
 * @param val : the value to be inserted
 */
void initOneValVec(double *unitVec, int n, int val);

/**
 * multiplies a vector and the modularity matrix B-roof, of a specific subgroup
 * @return
 * @param sp : the sparse matrix
 * @param vec : the vector to be multiplied by
 * @param res : the vector result of the multiplication
 * @param groupSize : the number of values inserted, inserted into the first indexes
 * @param vecF : an array containing the sum of values for each column
 * @return the time taken to finish the method
 */
double multBRoof(spmat *sp, double *vec, double *res, double *vecF);

/**
 * Power iteration method, multiplying a random vector with the matrix B, to find max eigenvector
 * @param sp : sparse matrix
 * @param b0  : initial random vector
 * @param shifting  : shifting value to shift the matrix & find max positive eigenvalue
 * @param result : the vector result of the power iteration
 * @param groupSize : the number of values inserted, inserted into the first indexes
 * @param vecF : an array containing the sum of values for each column
 */
void
powerIter(BMat *B, double *b0, double *result);

/**
 * Calculates the eigenvalue of a vector
 * @param sp : sparse matrix
 * @param vec : the eigenvector
 * @param groupSize : the number of values inserted, inserted into the first indexes
 * @param vecF : an array containing the sum of values for each column
 * @return : the eigenvalue
 */
double eigenValue(BMat *B, double *vec,double *tmp);

/**
 * calculates the modularity of a division vector for a subgroup
 * @param sp : sparse matrix
 * @param vec : the eigenvector
 * @param group : the subgroup array containing the vertices of the group
 * @param groupSize : the number of values inserted, inserted into the first indexes
 * @param vecF : an array containing the sum of values for each column
 * @return the modularity of the division
 */
double modularityCalc(spmat *sp, double *vec,double *tmp, double *vecF);

/**
 *  A method to print the type of error event in runtime
 * @param errorCode : error to be printed
 */
void error(int errorCode);


#endif
