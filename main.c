#include <stdio.h>
#include <stdlib.h>
#include "spmat.h"
#include "divider.h"

void printVector(double *vec, int n) {
    int i;
    for (i = 0; i < n; ++i) {
        printf("%f\t", vec[i]);
    }
    printf("\n");
}

int main(int argc, char **argv) {
    int elem[1];
    int size;
    int *row;
    double *res, *b0;
    spmat *graph;
    int i;
    int j;
    division *div;
    unsigned int n;
    FILE *input;
    if (argc != 3) {
        printf("ERROR - there is not 2 arguments");
        exit(EXIT_FAILURE);
    }
    input = fopen(argv[1], "rb");
    fread(elem, sizeof(int), 1, input);
    size = *elem;
    graph = spmat_allocate_list(size);
//    printf("%d\n", *elem);
    row = malloc(sizeof(int) * size);
    if (row == NULL) {
        printf("ERROR - memory allocation unsuccessful");
        exit(EXIT_FAILURE);
    }
    for (i = 0; i < size; ++i) {
        fread(elem, sizeof(int), 1, input);
        fread(row, sizeof(int), *elem, input);
        graph->add_row(graph, row, i, *elem);
        printf("%d - \t", i);
        for (j = 0; j < *elem; ++j) {
            printf("%d\t", row[j]);
        }
        printf("\n");
    }
    div = allocateDivision(graph->n);
    graph->print_list(graph);
    printf("%f\n", graph->matShifting(graph));
    b0 = malloc(sizeof(double) * graph->n);
    res = malloc(sizeof(double) * graph->n);
    randomizeVec(graph->n, b0);
    powerIter(graph, b0, graph->matShifting(graph), res);
    printVector(res, graph->n);
    div->split(div, graph, res);
    div->printGroups(div);
    return 0;
}
