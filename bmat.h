
#ifndef CPROJECT_UTILS_H
#define CPROJECT_UTILS_H

#include "spmat.h"

/**
@file bmat.h
**Author:** Ofek Bransky & Gal Cohen
**Date:**  18.9.2020
* Summary:
* This is the B matrix header file, contains the B matrix struct, and various functions for different calculations with B
 * we maintain a B matrix for each subgroup (community) of vertices
 * the struct 'bmat' is a struct containing all components of the B matrix from the theoretical modularity calculations
 * Functions:
 * readGraphB - reads the graph into a initial B matrix
 * allocateB - allocates a B struct
 * randomizeVec - returns a random vector
 * multBv - multiplies a vector v with the B matrix
 * initOneValVec - initializes a vector, all values are equal to the value input of the function
 * multBRoof - multiplies a vector v with the B^ matrix (B roof matrix of a subgroup)
 * powerIter - power iteration algorithm
 * eigenValue - calculates the eigenvalue
 * modularityCalc - calculates the group's modularity
 * error - returns an error based on the input (prints the error)
*/



/**
 * The B struct
 * @param sp : a pointer to the Sparse matrix associated with the B matrix
 * @param vecF : a vector containing the sum of colums for each row in the matrix
 * @param shifting : the shifting value, we shift the matrix with this value to get positive eigenvalues
 * @function free : frees the struct
 */
typedef struct _bmat{
    spmat *sp;
    double *vecF;
    double shifting;
    void (*free)(struct _bmat *B);
}BMat;

/**
 * reads into a initial B struct the input graph from the binary file
 * @param input : the input file pointer
 * @return a pointer to the bmat struct
 */
BMat *readGraphB(FILE *input);

/**
 * allocates a new B matrix
 * @return a pointer to the matrix
 */
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
 * @return : returns 1 on success
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
 * @param sp : the sparse matrix
 * @param vec : the vector to be multiplied by
 * @param res : the vector result of the multiplication
 * @param vecF : an array containing the sum of values for each column
 * @return  returns 1 on success
 */
double multBRoof(spmat *sp, double *vec, double *res, double *vecF);

/**
 * Power iteration method, multiplying a random vector with the matrix B, to find max eigenvector
 * @param B : the B matrix
 * @param b0  : initial random vector
 * @param result : the vector result of the power iteration
 * @return the power iteration result into the result vector
 */
void
powerIter(BMat *B, double *b0, double *result);

/**
 * Calculates the eigenvalue of a vector
 * @param B : B matrix
 * @param vec : the eigenvector
 * @param tmp : a temporary vector to assist calculation
 * @return : the eigenvalue
 */
double eigenValue(BMat *B, double *vec,double *tmp);

/**
 * calculates the modularity of a division vector for a subgroup
 * @param sp : sparse matrix
 * @param vec : the eigenvector
 * @param tmp : a temporary vector to assist calculation
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
