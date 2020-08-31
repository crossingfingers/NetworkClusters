#include <stdio.h>
#include <stdlib.h>
#include "spmat.h"
#include "divider.h"
#include <time.h>
//

int main(int argc, char **argv) {
    spmat *graph = NULL;
    division *div;
    FILE *input;
    FILE *output;
    double start;
    double end;
    if (argc != 3) {
        printf("ERROR - there is not 2 arguments");
        exit(EXIT_FAILURE);
    }
    start = clock();
    input = fopen(argv[1], "rb");
    output = fopen(argv[2], "wb");
    graph = readGraph(input);
//    graph->printSprase(graph);
    div = allocateDivision(graph->n);
    div->findGroups(div, graph);
    div->printGroups(div);
    div->writeDivision(div, output);
    div->free(div);
    graph->free(graph);
    end = clock();
    printf("took %f seconds\n", ((double)(end-start)/ CLOCKS_PER_SEC));
    free(graph);
    free(div);
    fclose(input);
    fclose(output);
    return 0;
}
