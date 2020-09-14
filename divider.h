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
    int *unmoved;
    int *indices;
    double *score;
    double *improve;
    double *res;
    double Q;



    void (*free)(struct _division *d);

    void (*printGroups)(struct _division *d);

    void (*writeDivision)(struct _division *div, FILE *output);

    void (*findGroups)(struct _division *div, networks *graphs);

} division;

division *allocateDivision(int n);


#endif //CPROJECT_DIVIDER_H
