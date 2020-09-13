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

    /*Checks if program input arguments are valid*/
    srand(time(NULL));
    if (argc != 3) {
        printf("ERROR - there is not 2 arguments");
        exit(EXIT_FAILURE);
    }

    input = fopen(argv[1], "rb");
    output = fopen(argv[2], "wb");
    start = clock();
    graphs = readGraph(input);
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
