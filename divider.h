
#ifndef CPROJECT_DIVIDER_H
#define CPROJECT_DIVIDER_H

#include "bmat.h"

/**
@file divider.h
**Author:** Ofek Bransky & Gal Cohen
**Date:**  18.9.2020
*Summary:
 * This is the Divider header file, mantains all main methods to find the graph subgroups
 * the struct 'networks' maintains an array to gather all data for the B matrices during the algorithm
 * the struct 'Division' maintains all data to find the communities in the input graph
 * the function 'allocateDivision' allocates the division struct
 * the function 'readGraph' reads the input graph into the program, gives an error if there are any problems during read
*/


/** a struct that contains all the B matrices in the program
 * @param B : a array containing pointers to a set of B matrices, each belonging to a subgroup
 * @paran n : number of vertices in the graph
 * @param M : sum of all vertice ranks in the graph
 * @param tmp : a temporary vector used during calculations (such as modularity calculation and eigenvalue calculation
 * @function free : frees the struct
 */
typedef struct _networks{
    BMat **B;
    int n;
    int M;
    double *tmp;
    void (*free)(struct _networks *graphs, int numOfGroups);
}networks;


/**Division struct maintains the groups found during the algorithm run,
 * @param groups : an array containing subgroups of the vertices
 * @param numOfGroups : the number of subgroups (communities) in the graph
 * @param nodesforGroup : an array containing the number of nodes in each subgroup
 * @param improve : an array keeping the improvement in modularity after each vertice movement
 * @param indices : an array that keeps the order of vertices moved, during the optimization
 * @param unmoved : an array that keeps track which vertice hasn't been moved
 * @param score : an array keeping the score(modularity) of each vertice in the subgroup
 * @param res : a vector used for calculations
 * */
typedef struct _division {
    int **groups;
    int numOfGroups;
    int *nodesforGroup;
    int *unmoved;
    int *indices;
    double *score;
    double *improve;
    double *res;

    /**frees the division struct*/
    void (*free)(struct _division *d);
    /**writes the division subgroups found*/
    void (*writeDivision)(struct _division *div, FILE *output);
    /**finds the subgroups, based on modularity*/
    void (*findGroups)(struct _division *div, networks *graphs);

} division;

/**reads the graph input into the program*/
networks *readGraph(FILE *input);

/**function to allocate a division struct*/
division *allocateDivision(int n);


#endif
