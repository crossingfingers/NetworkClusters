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
    graph = spmat_allocate_list(size);
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
        printf("%d - \t", i);
        for (j = 0; j < elem; ++j) {
            printf("%d\t", row[j]);
        }
        printf("\n");
    }
    free(row);
    return graph;
}

int main(int argc, char **argv) {
    spmat *graph = NULL;
    division *div;
    FILE *input;
    if (argc != 3) {
        printf("ERROR - there is not 2 arguments");
        exit(EXIT_FAILURE);
    }
    input = fopen(argv[1], "rb");
    graph = readGraph(input);
    div = allocateDivision(graph->n);
    graph->print_list(graph);
    findGroups(div, graph);
    div->free(div);
    graph->free(graph);
    free(graph);
    free(div);
    return 0;
}
