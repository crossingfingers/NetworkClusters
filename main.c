#include <stdio.h>
#include <stdlib.h>
#include "spmat.h"
#include "divider.h"
//
//void printVector(double *vec, int n) {
//    int i;
//    for (i = 0; i < n; ++i) {
//        printf("%f\t", vec[i]);
//    }
//    printf("\n");
//}

spmat *readGraph(FILE *input) {
    spmat *graph;
    int i, size, elem, *row, j;
    unsigned int n;
    n = fread(&elem, sizeof(int), 1, input);
    if (n != 1) {
        printf("ERROR - mismatch reading value");
        exit(EXIT_FAILURE);
    }
    size = elem;
//  graph = spmat_allocate_list(size);
  graph = spmat_allocate_array(size,16);
//    printf("%d\n", *elem);
    row = malloc(sizeof(int) * size);
    if (row == NULL) {
        printf("ERROR - memory allocation unsuccessful");
        exit(EXIT_FAILURE);
    }
    for (i = 0; i < size; ++i) {
        n = fread(&elem, sizeof(int), 1, input);
        if (n != 1) {
            printf("ERROR - mismatch reading value");
            exit(EXIT_FAILURE);
        }
        n = fread(row, sizeof(int), elem, input);
        if (n != elem) {
            printf("ERROR - mismatch reading value");
            exit(EXIT_FAILURE);
        }

        graph->add_row(graph, row, i, elem);
       // graph->print_list(graph);
//        printf("%d - \t", i);
//        for (j = 0; j < elem; ++j) {
//            printf("%d\t", row[j]);
//        }
//        printf("\n");
    }
    free(row);
    graph->print_list(graph);
    return graph;
}

int main(int argc, char **argv) {
    spmat *graph = NULL;
    int i;
    division *div;
    FILE *input;
    FILE *output;
    if (argc != 3) {
        printf("ERROR - there is not 2 arguments");
        exit(EXIT_FAILURE);
    }
    input = fopen(argv[1], "rb");
  //  output = fopen(argv[2], "wb");
    graph = readGraph(input);
    printf("\nM is:%d\n",graph->M);
    double a[7]={1,2,3,4,5,6,7};
    for( i=0;i<7;i++)
    {printf("%f ",a[i]);}
    printf("\n");
    double re[7]={0,0,0,0,0,0,0};
    graph->mult(graph,a,re);
    for( i=0;i<7;i++)
    {printf("%f ",re[i]);}

  //  div = allocateDivision(graph->n);
 //   graph->print_list(graph);

//    double divVec[]={1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
//    divOptimization(div,0,0,divVec,graph);
     //findGroups(div, graph);
//
//    div->writeDivision(div, output);


 //   div->free(div);
    graph->free(graph);

   free(graph);
   // free(div);
    fclose(input);
   // fclose(output);
    return 0;
}
