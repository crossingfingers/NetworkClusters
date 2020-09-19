#include <stdio.h>
#include <stdlib.h>
#include "divider.h"
#include <time.h>

/**
@file main.c
**Author:** Ofek Bransky & Gal Cohen
**Date:**  18.9.2020
 * main.c Summary:
 * This is the main file, reads the program input, runs the algorithms,
 * and saves the output to the specified directory.
 * if the specified path is incorrect or not found, an error will be given.

*/


/*Reads input, calculates groups, and saves to output */
int main(int argc, char **argv) {
    networks *graphs = NULL;
    division *div;
    FILE *input;
    FILE *output;

    /*Checks if program arguments are OK */
    srand(time(NULL));
    if (argc != 3) {
        error(ARGSERROR);
        exit(EXIT_FAILURE);
    }

    /*Reads input */
    input = fopen(argv[1], "rb");
    output = fopen(argv[2], "wb");
    graphs = readGraph(input);

    /* Allocates group division struct, runs findGroup Algorithm */
    div = allocateDivision(graphs->n);
    div->findGroups(div, graphs);

    /*Writes division to binary file */
    div->writeDivision(div, output);
    /*Frees all allocated data*/
    graphs->free(graphs, div->numOfGroups);
    div->free(div);
    free(graphs);
    free(div);
    fclose(input);
    fclose(output);
    return 0;
}
