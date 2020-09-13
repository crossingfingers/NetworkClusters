//
// Created by gal21 on 12/08/2020.
//

#ifndef CPROJECT_DIVIDER_H
#define CPROJECT_DIVIDER_H
#include "spmat.h"

/**
@file divider.h
**Author:** Ofek Bransky & Gal Cohen
**Date:**  18.9.2020
## This is the Divider header file, mantains all main methods to find the graph subgroups
*/



/**Division struct maintains the groups found during the algorithm run,
 *@param n : number of vertices in the graph
 *@param *groups : an array containing subgroups of the vertices
 *@param numOfGroups : the number of subgroups (communities) in the graph
 *@param *nodesofGroup : an array containing the number of nodes in each subgroup
 * */

typedef struct _division {
    int n;
    int **groups;
    int numOfGroups;
    int *nodesforGroup;
    double Q;

    /**frees the division struct*/
    void (*free)(struct _division *d);
    /**prints the division subgroups*/
    void (*printGroups)(struct _division *d);
    /**writes the division subgroups found*/
    void (*writeDivision)(struct _division *div, FILE *output);
    /**finds the subgroups, based on modularity*/
    void (*findGroups)(struct _division *div, networks *graphs);

} division;
/**function to allocate a division struct*/
division *allocateDivision(int n);

#endif //CPROJECT_DIVIDER_H
