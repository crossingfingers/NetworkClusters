#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void initMat(int **mat, int n) {
    int i;
    int j;
    for (i = 0; i < n; ++i) {
        mat[i] = malloc(sizeof(int) * n);
        for (j = 0; j < n; ++j) {
            mat[i][j] = 0;
        }
    }
}

int calcK(int **row, int n, int i) {
    int j;
    int sum = 0;
    for (j = 0; j < n; ++j) {
        sum += row[i][j];
    }
    return sum;
}

int main(int argc, char **argv) {
    int size = atoi(argv[2]);
    FILE *out = fopen(argv[1], "wb");
    int k;
    int i;
    int j;
    int prob = 3; /*the probability to get edge between two vertex*/
    int val;
    int **mat;
    srand(time(NULL));
    mat = malloc(sizeof(int *) * size);
    initMat(mat, size);
    for (i = 0; i < size; ++i) {
        for (j = 0; j < size; ++j) {
            if (j==i)
                continue;
            if (mat[i][j] != 1 && (rand() % (prob)) == 1) {
                mat[i][j] = 1;
                mat[j][i] = 1;
            }
        }
    }
    for (i = 0; i < size; ++i) {
        k = calcK(mat, size, i);
        printf("i=%d, k=%d  -  \t", i, k);
        fwrite(&k, sizeof(int), 1, out);
        for (j = 0; j < size; ++j) {
            if (mat[i][j] != 0) {
                val = j;
                printf("%d\t", val);
                fwrite(&val, sizeof(int), 1, out);
            }
        }
        free(mat[i]);
        printf("\n");
    }
    free(mat);
    fclose(out);
    return 0;

}
