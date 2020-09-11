#include <stdlib.h>
#include <stdio.h>
#include "spmat.h"
#include <math.h>

typedef struct _node {
    struct _node *next;
    int col_idx;
    int value;
} node;

void add_row_list(struct _spmat *A, int *row, int i, int k) {
    int j;
    node **rows = (node **) A->private;
    node *idx;
    node *temp = NULL;
    node *now = NULL;
    A->M += k;
    A->k[i] = k;
    if (k == 0) {
        rows[i] = malloc(sizeof(node));
        rows[i] = NULL;
    }
    for (j = 0; j < k; j++) {
        if (rows[i] == NULL) {
            rows[i] = malloc(sizeof(node));
            idx = rows[i];
            idx->value = 1;
            idx->next = NULL;
            idx->col_idx = row[j];
            now = idx;
        } else {
            temp = malloc(1 * sizeof(node));
            temp->next = NULL;
            temp->value = 1;
            temp->col_idx = row[j];
            now->next = temp;
            now = temp;
        }
    }
}


void print_list(struct _spmat *A) {
    int i;
    int j;
    node *curr;
    node **rows = (node **) A->private;
    for (i = 0; i < A->n; ++i) {
        printf("%d - \t", i);
        curr = rows[i];
        for (j = 0; j < A->k[i]; ++j) {
            printf("%d\t", curr->col_idx);
            curr = curr->next;
        }
        printf("\n");
    }
};

void free_list(struct _spmat *A) {
    node **rows = (node **) A->private;
    int i;
    node *head = NULL, *temp;
    for (i = 0; i < A->n; i++) {
        head = rows[i];
        while (head != NULL) {
            temp = head;
            head = head->next;
            free(temp);
        }
    }
    free(rows);
}


void mult_list(const struct _spmat *A, const double *v, double *result, const int *group,
               int groupSize, const int *groupToVertice) {
    double res;
    int i, j;
    node *now;
    node **rows = (node **) A->private;
    for (i = 0; i < groupSize; i++) {
        res = 0;
        now = rows[group[i]];
        j = 0;
        while (now != NULL) {
            if (now->col_idx == group[j]) {
                res += v[now->col_idx];
                j++;
                now = now->next;
            } else if (now->col_idx < group[j]) {
                now = now->next;
            } else {
                j++;
            }
        }
        result[group[i]] = res;
    }
}

double listShifting(spmat *A, const int *group, int groupSize, const int *vertexToGroup, int groupIdx, double *F) {
    double max = 0;
    double sum = 0;
    int val = 0;
    node *curr;
    int ki, kj;
    node **rows = (node **) A->private;
    int i, idxI;
    int idxJ;
    for (i = 0; i < groupSize; ++i) {
        idxI = group[i];
        curr = rows[idxI];
        sum = 0;
        ki = A->k[idxI];
        for (idxJ = 0; idxJ <= group[groupSize - 1]; ++idxJ) {
            kj = A->k[idxJ];
            if (curr != NULL && curr->col_idx == idxJ) {
                val = 1;
                curr = curr->next;
            } else {
                val = 0;
                if (curr != NULL && curr->col_idx < idxJ)
                    curr = curr->next;
            }
            if (vertexToGroup[idxJ] == groupIdx) {
                if (idxI != idxJ) {
                    sum += fabs((double) val - ((double) (ki * kj) / A->M));
                } else {
                    sum += fabs((double) val - ((double) (ki * kj) / A->M) - (double) F[i]);
                }
            }
        }
        max = (max >= sum) ? max : sum;
    }
    return max;
}

void initk(spmat *A) {
    int i;
    for (i = 0; i < A->n; ++i) {
        A->k[i] = 0;
    }
}


spmat *spmat_allocate_list(int n) {
    spmat *sp;
    node **rows = calloc(n, sizeof(node *));
    sp = malloc(sizeof(spmat));
    sp->n = n;
    sp->M = 0;
    sp->k = malloc(sizeof(int) * n);
    initk(sp);
    sp->add_row = add_row_list;
    sp->free = free_list;
    sp->mult = mult_list;
    sp->private = rows;
    sp->printSprase = print_list;
    sp->matShifting = listShifting;
    return sp;
}

typedef struct _array {
    int *colind;
    int *rowptr;
    int lastindex;
    int lastRowPtr;
    int nnz;
    int groupSize;
} array;

//TODO remove values array, no need....
void add_row_array(struct _spmat *A, int *row, int i, int k) {
    array *sparray = (array *) A->private;
    int *rowInput = row;
    int ci;
    A->M += k;
    A->k[i] = k;
    sparray->rowptr[sparray->lastRowPtr + 1] = sparray->rowptr[sparray->lastRowPtr] + k;
    sparray->lastRowPtr++;
/*updates values array*/
    for (ci = 0; ci < k; ci++) {
        sparray->colind[sparray->lastindex] = *rowInput;
        rowInput++;
        sparray->lastindex++;
    }
}

void mult_array(const struct _spmat *A, const double *vec, double *result, const int *group,
                int groupSize, const int *verticeToGroup) {

    array *sparray = (array *) A->private;
    int *rowPtr = sparray->rowptr;
    int *cols;
    int i;
    int colStart;
    int colEnd;
    double res;
    int size=0;
    for (i = 0; i < groupSize; ++i) {
        res = 0;
        colStart = *(rowPtr + i);

        cols = sparray->colind + colStart;
        colEnd = *(rowPtr + i + 1);
        size=colEnd-colStart;
        while (colStart < colEnd) {
            res += vec[*(cols)];
            colStart++;
            if (colStart < colEnd) { cols++; }
        }
        result[group[i]] = res;
    }
}

void free_array(struct _spmat *A) {
    array *sparray = (array *) A->private;
    free(sparray->rowptr);
    free(sparray->colind);
    free(A->private);
//    free(A->k);
}

void printMatrix(spmat *A) {

    int i;
    array *sparr = (array *) A->private;


    printf("\ncollarr:\n");
    for (i = 0; i < sparr->nnz; i++) {
        printf("%d ", sparr->colind[i]);
    }
    printf("\nrowptrarr:\n");
    for (i = 0; i <= sparr->groupSize; i++) {
        printf("%d ", sparr->rowptr[i]);
    }

}

void print_array(struct _spmat *A) {
    int i;
    int j;
    array *sparray = (array *) A->private;
    int *colindarr = (int *) sparray->colind;
    int counter = 0;
    for (i = 0; i < A->n; ++i) {
        printf("%d - \t", i);

        for (j = 0; j < A->k[i]; ++j) {
            printf("%d\t", colindarr[counter]);
            counter++;
        }
        printf("\n");
    }
    printMatrix(A);
    printf("\n");
}

spmat *spmat_allocate_array(int n, int groupSize, int nnz);

double arrayShifting(spmat *A, const int *group, int groupSize, const int *vertexToGroup, int groupIdx, double *F) {
    double max = 0;
    double sum;
    int val;
    array *sparray = A->private;
    int *rowPtr = sparray->rowptr;
    int *cols;
    int colStart;
    int colEnd;
    int i;
    int j;
    int ki;
    int kj;
    int M = A->M;
    int vertice1;
    double Fi;

    for (i = 0; i < groupSize; ++i) {
        sum = 0;
        vertice1 = group[i];
        ki = A->k[vertice1];
        Fi = (double) F[vertice1];
        colStart = *(rowPtr + i);
        cols = sparray->colind + colStart;
        colEnd = *(rowPtr + i + 1);
        for (j = 0; j < groupSize; ++j) {
            kj = A->k[group[j]];
            if ((*(cols) == group[j]) && (colStart < colEnd)) {
                val = 1;
                colStart++;
                if (colStart < colEnd) { cols++; }
            } else {
                val = 0;
            }
            if (vertice1 != group[j]) {
                sum += fabs((double) val - ((double) (ki * kj) / M));
            } else {
                sum += fabs(
                        (double) val - ((double) (ki * kj) / M) - Fi);
            }
        }
        max = (max >= sum) ? max : sum;
    }
    return max;
}

//g1(2,6,8)
//2: 0, 4:6, 6:4,7,8, 7:6,8, 8:6,7
int getNNZforGroup(array *arrAg, int *g1, int gSize) {
    int i;
    int j;
    int colIDX = 0;
    int currG = 0;
    int numOfVals;
    int NNZ = 0;
    for (i = 0; i < gSize; i++) {
        colIDX = arrAg->rowptr[i];
        numOfVals = (arrAg->rowptr[i + 1] - colIDX);
        j = 0;
        currG = 0;
        printf("vertex %d\n:", g1[i]);
        while ((j < numOfVals)) {
            printf(" %d ", arrAg->colind[colIDX + j]);
            if (g1[currG] == arrAg->colind[colIDX + j]) {
                NNZ++;
                j++;
                if ((currG < gSize))
                    currG++;
                printf("+1");
            } else {
                if ((currG < gSize) && (g1[currG] < arrAg->colind[colIDX + j])) { currG++; }
                else { j++; }
            }
        }
        printf("\n");
    }
    return NNZ;
}

void getNNZ(int *nnz1, int *nnz2, array *arrAg, int *g,int *g1, int *g2, int gSize, int g1Size, int g2Size) {
    int i, j;
    int *row = arrAg->rowptr;
    int *cols = arrAg->colind;
    int NNZ1 = 0;
    int NNZ2 = 0;
    int currG1 = 0;
    int currCol = 0;
    int currG2 = 0;
    int numOfVals = 0;
    int colIdx = 0;

    for (i = 0; i < gSize; i++) {
        j = 0;
        colIdx = row[i];
        numOfVals = row[i + 1] - colIdx;
        currCol = 0;
        if (g[i] == g1[currG1]) {
            while (j < numOfVals) {
                //if col[row[i]] is from group then add NNZ
                if (g1[currCol] == cols[colIdx + j]) {
                    NNZ1++;
                    j++;
                    if ((currCol < g1Size))
                        currCol++;
                } else {
                    if ((currCol < gSize) && (g1[currCol] < arrAg->colind[colIdx + j])) { currCol++; }
                    else { j++; }
                }
            }
            currG1++;
        } else {
            while (j < numOfVals) {
                //if col[row[i]] is from group then add NNZ
                if (g2[currCol] == cols[colIdx + j]) {
                    NNZ2++;
                    j++;
                    if ((currCol < g2Size))
                        currCol++;
                } else {
                    if ((currCol < gSize) && (g2[currCol] < arrAg->colind[colIdx + j])) { currCol++; }
                    else { j++; }
                }
            }
            currG2++;
        }

    }
    *nnz1 = NNZ1;
    *nnz2 = NNZ2;
}
void insertVals( array *arrAg,array *arrAg1, array *arrAg2, int *g,int *g1, int *g2, int gSize, int g1Size, int g2Size){
    int i, j;
    int *row = arrAg->rowptr;
    int *cols = arrAg->colind;
    int currG1row = 0;
    int currG1col=0;
    int currCol = 0;
    int currG2row = 0;
    int currG2col=0;

    int numOfVals = 0;
    int colIdx = 0;
    int numOfNewVal=0;
    for (i = 0; i < gSize; i++) {
        j = 0;
        numOfNewVal=0;
        colIdx = row[i];
        numOfVals = row[i + 1] - colIdx;
        currCol = 0;
        if (g[i] == g1[currG1row]) {
            while (j < numOfVals) {
                if (g1[currCol] == cols[colIdx + j]) {
                    arrAg1->colind[currG1col] = arrAg->colind[colIdx + j];
                    j++;
                    numOfNewVal++;
                    currG1col++;
                    if ((currCol < g1Size))
                        currCol++;
                } else {
                    if ((currCol < gSize) && (g1[currCol] < arrAg->colind[colIdx + j])) { currCol++; }
                    else { j++; }
                }
            }
            arrAg1->rowptr[currG1row + 1] = arrAg1->rowptr[currG1row] +numOfNewVal;
            currG1row++;
        }
        else {
            while (j < numOfVals) {
                if (g2[currCol] == cols[colIdx + j]) {
                    arrAg2->colind[currG2col] = arrAg->colind[colIdx + j];
                    j++;
                    numOfNewVal++;
                    currG2col++;
                    if ((currCol < g2Size))
                        currCol++;
                } else {
                    if ((currCol < g2Size) && (g2[currCol] < arrAg->colind[colIdx + j])) { currCol++; }
                    else { j++; }
                }
            }
            currG2row++;
            arrAg2->rowptr[currG2row + 1] = arrAg2->rowptr[currG2row] +numOfNewVal;

        }
    }
}
void insertValsToArr(array *arrAg, array *arrAg1, int *g1, int g1Size) {
    int i;
    int j;
    int numOfVals;
    int colIDX;
    int currG1;
    int currG = 0;
    for (i = 0; i < g1Size; i++) {
        numOfVals = arrAg->rowptr[g1[i] + 1] - arrAg->rowptr[g1[i]];
        arrAg1->rowptr[i + 1] = arrAg1->rowptr[i] +numOfVals;
        colIDX = arrAg->rowptr[g1[i]];
        j = 0;
        currG1 = 0;

        while (j < numOfVals) {
            if (arrAg->colind[colIDX + j] == g1[currG1]) {
                arrAg1->colind[currG] = arrAg->colind[colIDX + j];
                currG++;
                j++;
                arrAg1->rowptr[i + 1]++;
                if (currG1 < g1Size) { currG1++; }
            } else {
                if ((currG1 < g1Size) && (g1[currG1] < arrAg->colind[colIDX + j])) { currG1++; }
                else { j++; }
            }
        }
    }


}

void splitGraphArray(networks *graphs, int groupIdx, int newGroupIdx,int *g, int *g1, int *g2, int g1Size, int g2Size) {
    if (g2Size == 0) { return; }
    spmat **Amats = graphs->A;
    int n = Amats[0]->n;
    spmat *currSP = Amats[groupIdx];
    array *arrAg = (array *) currSP->private;
    spmat *Ag1;
    spmat *Ag2;
    int *g1NNZ=malloc(sizeof(int));
    int *g2NNZ=malloc(sizeof(int));
    getNNZ(g1NNZ,g2NNZ,arrAg,g,g1,g2,g1Size+g2Size,g1Size,g2Size);
    Ag1 = spmat_allocate_array(graphs->n, g1Size, *g1NNZ);
    Ag2 = spmat_allocate_array(graphs->n, g2Size, *g2NNZ);
    array *arrAg1 = Ag1->private;
    array *arrAg2 = Ag2->private;
    insertVals(arrAg,arrAg1,arrAg2,g,g1,g2,g1Size+g2Size,g1Size,g2Size);

    //TODO change the method of spmat_allocate to work without malloc of k. should be outside, they all use the same k, M too.
    free(Ag1->k);
    free(Ag2->k);
    Ag1->k = graphs->k;
    Ag2->k = graphs->k;
    Ag1->M = graphs->M;
    Ag2->M = graphs->M;
    currSP->free(currSP);
    free(g1NNZ);
    free(g2NNZ);
    free(currSP);
    Amats[groupIdx] = Ag1;
    Amats[newGroupIdx] = Ag2;
}

spmat *spmat_allocate_array(int n, int groupSize, int nnz) {
    spmat *sp;
    array *sparray = malloc(sizeof(array));
    sparray->colind = calloc(nnz, sizeof(int));
    sparray->rowptr = calloc(groupSize + 1, sizeof(int));
    sparray->lastindex = 0;
    sparray->lastRowPtr = 0;
    sparray->rowptr[0] = 0;
    sparray->nnz = nnz;
    sparray->groupSize = groupSize;
    sp = malloc(sizeof(spmat));
    sp->n = n;
    sp->add_row = add_row_array;
    sp->free = free_array;
    sp->mult = mult_array;
    sp->private = sparray;
    sp->printSprase = print_array;
    sp->splitGraph = splitGraphArray;
    sp->M = 0;
    sp->k = malloc(sizeof(int) * n);
    initk(sp);
    sp->matShifting = arrayShifting;
    return sp;
}

int find_nnz(FILE *input) {
    int res;
    int size;
    fread(&size, sizeof(int), 1, input);
    fseek(input, 0L, SEEK_END);
//    printf("size of file is %ld\n", ftell(input));
//    printf("size is %d\n", size);
    res = ftell(input) - (int) ((1 + size) * sizeof(int));
//    printf("res is %d\n", res);
    fseek(input, 0, SEEK_SET);
    return res / 4;
}

void freeNetworks(networks *graphs, int numOfGroups) {
    int i;
    spmat *sp, **mats = graphs->A;
    free(graphs->k);
    //TODO change 1 to numOfGroups!!
    for (i = 0; i < 1; ++i) {
        sp = *mats++;
        sp->free(sp);
        free(sp);
    }
}

networks *allocateNetworks(int n) {
    networks *graphs = malloc(sizeof(networks));
    if (graphs == NULL) {
        printf("ERROR - memory allocation unsuccessful");
        exit(EXIT_FAILURE);
    }
    graphs->n = n;
    graphs->A = malloc(sizeof(spmat *) * n);
    if (graphs->A == NULL) {
        printf("ERROR - memory allocation unsuccessful");
        exit(EXIT_FAILURE);
    }
    graphs->free = freeNetworks;
    return graphs;
}

networks *readArray(FILE *input) {
    spmat *graph;
    networks *graphs;
    int i, size, elem, *row, j;
    unsigned int n;
    int nnz;
    nnz = find_nnz(input);
//    printf("nnz is %d\n", nnz);
    n = fread(&elem, sizeof(int), 1, input);
    if (n != 1) {
        printf("ERROR - mismatch reading value");
        exit(EXIT_FAILURE);
    }
    size = elem;
//    printf("n: %d, nnz: %d\n", elem, nnz);
    graph = spmat_allocate_array(size, size, nnz);
    graphs = allocateNetworks(size);
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
    }
    free(row);
//    graph->printSprase(graph);
    graphs->M = graph->M;
    graphs->k = graph->k;
    graphs->A[0] = graph;
    return graphs;
}

networks *readList(FILE *input) {
    spmat *graph;
    networks *graphs;
    int i, size, elem, *row, j;
    unsigned int n;
    n = fread(&elem, sizeof(int), 1, input);
    if (n != 1) {
        printf("ERROR - mismatch reading value");
        exit(EXIT_FAILURE);
    }
    size = elem;
    graph = spmat_allocate_list(size);
    graphs = allocateNetworks(size);
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
    }
    free(row);
//    graph->printSprase(graph);
    graphs->M = graph->M;
    graphs->k = graph->k;
    graphs->A[0] = graph;
    return graphs;
}

networks *readGraph(FILE *input, int type) {
    if (type == 2) {
        return readArray(input);
    } else
        return readList(input);
}