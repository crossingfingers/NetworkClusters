//
// Created by gal21 on 12/08/2020.
//

#ifndef CPROJECT_DIVIDER_H
#define CPROJECT_DIVIDER_H

typedef struct _division {
    int n;
    int *groupid;
    int numOfGroups;
    double Q;


    void (*split)(struct _division *d, spmat *A, double *vec, int group);

    void (*freeDivision)(struct _division *d);

    void (*printGroups)(struct _division *d);

} division;

division *allocateDivision(int n);

void powerIter(spmat *A, double *b0, double shifting, int group, const int *groupid ,double *result);

void randomizeVec(int size, double *vec);

double modularityCalc(spmat *A, double *vec, int group, const int *groupid);

#endif //CPROJECT_DIVIDER_H
