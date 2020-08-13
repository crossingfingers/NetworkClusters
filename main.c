#include <stdio.h>
#include <stdlib.h>
#include "spmat.h"
#include "divider.c"
int main(int argc, char **argv) {
    int elem[1];
    int size;
    int *row;
    spmat *graph;
    int i;
    int j;
    unsigned int n;
    FILE * input;
    if (argc != 3){
        printf("ERROR - there is not 2 arguments");
        return 1;
    }
    input = fopen(argv[1], "rb");
    fread(elem, sizeof(int), 1, input);
    size = *elem;
    graph = spmat_allocate_list(size);
//    printf("%d\n", *elem);
    row = malloc(sizeof(int)*size);
    if (row == NULL){
        printf("ERROR - memory allocation unsuccessful");
        return 1;
    }
    for(i = 0; i<size; ++i){
        fread(elem, sizeof(int), 1, input);
        fread(row, sizeof(int), *elem, input);
        graph->add_row(graph, row, i, *elem);
        printf("%d - \t", i);
        for(j=0; j<*elem; ++j){
            printf("%d\t", row[j]);
        }
        printf("\n");
    }
    graph->print_list(graph);
    return 0;
}
