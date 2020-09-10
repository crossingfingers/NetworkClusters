#include <stdio.h>
#include <stdlib.h>
#include "spmat.h"
#include "divider.h"
#include <time.h>
//

int main(int argc, char **argv) {
    networks *graphs = NULL;
    division *div;
    FILE *input;
    FILE *output;
    double start;
    double end;

    srand(time(NULL));
    if (argc != 3) {
        printf("ERROR - there is not 2 arguments");
        exit(EXIT_FAILURE);
    }

    input = fopen(argv[1], "rb");
    output = fopen(argv[2], "wb");
    start = clock();
    graphs = readGraph(input,2);
/*
    int i;
    double res[]={0,0,0, 0,0,0, 0,0,0};
    double vec[]={1,2,3,4,5,6,7,8,9};
    int group[]={1,3,5};
    int verticeToGroup[]={0,1,0, 1,0,1, 0,0,0};

    graph->mult(graph,vec,res,group,3,verticeToGroup);

    for(i=0;i<9;i++){printf(" %f ",res[i]);}
    printf("\n");*/

//    graph->printSprase(graph);
    div = allocateDivision(graphs->n);
    spmat **A=graphs->A;
    A[0]->printSprase(A[0]);
    int g1[]={0,1,2,3,4};
    int g2[]={5,6,7,8};
    A[0]->splitGraph(graphs,0,1,g1,g2,5,4);
    A[0]->printSprase(A[0]);
    A[1]->printSprase(A[1]);


    //TODO change the graphs->A[0] to graphs. after implement the creation of A[group]
//    div->findGroups(div, graphs->A[0]);
//    div->writeDivision(div, output);
//    graphs->free(graphs, div->numOfGroups);
//    div->free(div);
//    end = clock();
    printf("took %f seconds\n", ((double)(end-start)/ CLOCKS_PER_SEC));
    free(graphs);
    free(div);
    fclose(input);
    fclose(output);
    return 0;
}
