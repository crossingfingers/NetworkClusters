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


    /*   int i;
      double res[]={0,0,0, 0,0,0, 0,0,0};
      double vec[]={2,4,6,4,5,6,7,8,9};
      int group[]={0,1,2,3,4,5,6,7,8};
      int verticeToGroup[]={0,1,0, 1,0,1, 0,0,0};

  //    graph->printSprase(graph);
     spmat **A=graphs->A;
       printf("\n==========Printing g=========\n");
       A[0]->printSprase(A[0]);
       int g1[]={0,2,4,6,7,8};
       int g2[]={1,3,5};
       double s1[]={1,-1,1,-1,1,-1,1,1,1};
       A[0]->splitGraph(graphs,0,1,s1,group,9,6,3);
       double s2[]={1,-1,1,-1,1,-1,0,0,0};
       int group2[]={0,2,4,6,7,8};
       A[0]->splitGraph(graphs,0,2,s2,group2,6,3,3);
       printf("\n==========Printing g0 (0,4,7)=========\n");
       A[0]->printSprase(A[0]);
       printf("\n==========Printing g1 (1,3,5)=========\n");
       A[1]->printSprase(A[1]);
        printf("\n==========Printing g2 (2,6,8)=========\n");
         A[2]->printSprase(A[2]);
         printf("\n flag\n");
         graphs->A[0]->mult(graphs->A[1],vec,res,g2,3,verticeToGroup);
        for(i=0;i<9;i++){printf(" %f ",res[i]);}
         printf("\n");
*/

//    double res2[]={0,0,0, 0,0,0, 0,0,0};
//    double vec2[]={1,2,3,4,5,6,7,8,9};
//    int verticeToGroup2[]={0,1,0, 1,0,1, 0,0,0};
//    graphs->A[1]->mult(graphs->A[1],vec2,res2,g2,3,verticeToGroup);
//    for(i=0;i<9;i++){printf(" %f ",res2[i]);}
//    printf("\n");
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
    div->findGroups(div, graphs);
    div->writeDivision(div, output);
    graphs->free(graphs, div->numOfGroups);
    div->free(div);
    end = clock();
    printf("took %f seconds\n", ((double)(end-start)/ CLOCKS_PER_SEC));
    free(graphs);
    free(div);
    fclose(input);
    fclose(output);
    return 0;
}
