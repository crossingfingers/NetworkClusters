//
// Created by gal21 on 12/08/2020.
//

#ifndef CPROJECT_DIVIDER_H
#define CPROJECT_DIVIDER_H

#include "spmat.h"

typedef struct _division {
    int n;
    int **groups;
    int numOfGroups;
    int *nodesforGroup;
    int *vertexToGroup;
    double Q;



    void (*free)(struct _division *d);

    void (*printGroups)(struct _division *d);

    void (*writeDivision)(struct _division *div, FILE *output);

    void (*findGroups)(struct _division *div, networks *graphs);

} division;

division *allocateDivision(int n);

//void powerIter(spmat *A, double *b0, double shifting, int group, const int *groupid ,double *result);

//void randomizeVec(int size, double *vec);

//double modularityCalc(spmat *A, double *vec, int group, const int *groupid);

//double eigenValue(spmat *A, double *vec, int group, const int* groupid);

//int divideToTwo(division *div, spmat *sp, int group);



#endif //CPROJECT_DIVIDER_H
