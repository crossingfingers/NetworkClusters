//
// Created by gal21 on 05/09/2020.
//

#ifndef CPROJECT_UTILS_H
#define CPROJECT_UTILS_H

#include "spmat.h"

#define IS_POSITIVE(X) ((X) > 0.00001)

void randomizeVec(int size, double *vec, int groupSize, int *group);

double multBv(spmat *sp, double *vec, const int *group, double *res, int groupSize, int debug,int *verticeToGroup);

void initOneValVec(double *unitVec, int n, const int *group, int val);

double multBRoof(spmat *sp, double *vec, const int *group, int groupSize, double *res, double *vecF,int *verticeToGroup, int debug);

void powerIter(spmat *sp, double *b0, double shifting, int *group, int groupSize, double *result, double *vecF,int *verticeToGroup, int debug);

double eigenValue(spmat *sp, double *vec, const int *group, int groupSize, double *vecF,int *verticeToGroup);

double modularityCalc(spmat *sp, double *vec, int *group, int groupSize, double *vecF,int *verticeToGroup);

void printVector(double *vec, int n, const int *group);

void printIntVector(int *vec, int n);

#endif //CPROJECT_UTILS_H
