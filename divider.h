//
// Created by gal21 on 12/08/2020.
//

#ifndef CPROJECT_DIVIDER_H
#define CPROJECT_DIVIDER_H
#include "spmat.h"
typedef struct _division {
    int n;
    int *groupid;
    int numOfGroups;
    double Q;


    int (*split)(struct _division *d, spmat *A, double *vec, int group);

    void (*free)(struct _division *d);

    void (*printGroups)(struct _division *d);
    double (*modularityCalc)(spmat *A, double *vec, int group, const int *groupid);
    void (*divOptimization)(struct _division *div,int group,double q0,double *maxDiv, spmat *sp);
} division;

division *allocateDivision(int n);

void powerIter(spmat *A, double *b0, double shifting, int group, const int *groupid ,double *result);

void randomizeVec(int size, double *vec);

double modularityCalc(spmat *A, double *vec, int group, const int *groupid);

void divOptimization(division *div,int group,double q0,double *maxDiv, spmat *sp);

double eigenValue(spmat *A, double *vec, int group, const int* groupid);

int divideToTwo(division *div, spmat *sp, int group);

void findGroups(division *div, spmat *sp);

#endif //CPROJECT_DIVIDER_H
