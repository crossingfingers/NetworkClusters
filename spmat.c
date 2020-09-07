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
    int *values;
    int *colind;
    int *rowptr;
    int lastindex;
    int nnz;
} array;

void add_row_array(struct _spmat *A, int *row, int i, int k) {
    if (k) {
        array *sparray = (array *) A->private;
        int firstrowIdx;
        int ci;
        A->M += k;
        A->k[i] = k;
        if (!sparray->lastindex) { firstrowIdx = 0; }
        else { firstrowIdx = (sparray->rowptr[sparray->lastindex - 1]) + 1; }

/*updates values array*/
        for (ci = 0; ci < k; ci++) {
            sparray->values[sparray->lastindex] = 1;
            sparray->colind[sparray->lastindex] = *row;
            sparray->rowptr[sparray->lastindex] = firstrowIdx;
            sparray->lastindex++;
            row++;
        }

    }
}

void mult_array(const struct _spmat *A, const double *vec, double *result, const int *group,
                int groupSize, const int *groupToVertice) {

    array *sparray = (array *) A->private;
    int *rows = sparray->rowptr;
    int *cols = sparray->colind;
    double res;
    int rowIdx;
    int curr = 0;
    int nnz = sparray->nnz;
    int groupID = groupToVertice[group[0]];
    while (curr < nnz) {
        rowIdx = *(rows + curr);
        res = 0;
        while ((rowIdx == *(rows + curr)) && (curr < nnz)) {
            if ((*(groupToVertice + (*(cols + curr))) == groupID) &&
                (*(groupToVertice + (*(rows + curr))) == groupID)) {
                res += *(vec + (*(cols + curr)));
            }
            curr++;
        }
        *(result + rowIdx) = res;
    }
}

//       for (i = 0; i < A->n; i++) { result[i] = 0; }
//        for (curr = 0; curr < nnz; curr++) {
//            if ((groupToVertice[sparray->colind[curr]] == groupID) &&
//                (groupToVertice[sparray->rowptr[curr]] == groupID)) {
//                result[sparray->rowptr[curr]] += vec[sparray->colind[curr]];
//            }
//        }
//}
/* for (curr = 0; i < sparray->nnz; i++) {
     res = 0;
     rowIdx = sparray->rowptr[curr];
     while (rowIdx == sparray->rowptr[curr]) {
         if (groupToVertice[sparray->colind[curr]] == groupID)
             res += sparray->values[curr] * vec[sparray->colind[curr]];
         curr++;
     }
     result[rowIdx] = res;
 }*/

//    array *sparray = (array *) A->private;
//
//    int i;
//    int row = 0;
//    int col = 0;
//    int nnz = sparray->nnz;
//    int currRow;
//    printf("\nprinting vec:\n");
//    for(i=0;i<groupSize;i++)
//    {printf("%f ",vec[group[i]]);}
//    printf("\n");
//
//    for (i = 0; i < nnz; i++) {
//        if (sparray->rowptr[i] == group[row]) {
//            currRow = sparray->rowptr[i];
//            while (currRow == sparray->rowptr[i]) {
//                if (sparray->colind[i] == group[col]) {
//                    result[sparray->rowptr[i]] += vec[sparray->colind[i]];
//                    i++;
//                }
//                col++;
//            }
//            col=0;
//            row++;
//        }
//    }    printf("\nprinting vec:\n");
//    for(i=0;i<groupSize;i++)
//    {printf("%f ",result[group[i]]);}
//}

void free_array(struct _spmat *A) {
    array *sparray = (array *) A->private;
    free(A->k);
    free(sparray->values);
    free(sparray->rowptr);
    free(sparray->colind);
    free(A->private);

}

void printMatrix(spmat *A) {

    int i;
    array *sparr = (array *) A->private;
    printf("\nvalarr:\n");
    for (i = 0; i < sparr->nnz; i++) {
        printf("%d", sparr->values[i]);
    }
    printf("\ncollarr:\n");
    for (i = 0; i < sparr->nnz; i++) {
        printf("%d", sparr->colind[i]);
    }
    printf("\nrowptrarr:\n");
    for (i = 0; i <= sparr->nnz; i++) {
        printf("%d", sparr->rowptr[i]);
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
}

double arrayShifting(spmat *A, const int *group, int groupSize, const int *vertexToGroup, int groupIdx, double *F) {
    double max = 0;
    double sum;
    int val;
    array *sparray = A->private;
    int *rows = sparray->rowptr;
    int *cols = sparray->colind;
    int curr = 0;
    int i;
    int j;
    int ki;
    int kj;
    int M=A->M;
    double Fi;
    for (i = 0; i < groupSize; ++i) {
        sum = 0;
        ki = A->k[group[i]];
        Fi=(double) F[group[i]];
        for (j = 0; j < groupSize; ++j) {
            kj = A->k[group[j]];
            if ((*(rows + curr) == group[i]) && *(cols + curr) == group[j]) {
                val = 1;
                curr++;
            } else {
                val = 0;
            }
            if (group[i] != group[j]) {
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

spmat *spmat_allocate_array(int n, int nnz) {
    spmat *sp;
    array *sparray = malloc(sizeof(array));
    sparray->values = calloc(nnz, sizeof(int));
    sparray->colind = calloc(nnz, sizeof(int));
    sparray->rowptr = calloc(nnz + 1, sizeof(int));
    sparray->lastindex = 0;
    sparray->rowptr[nnz] = nnz;
    sparray->nnz = nnz;
    sp = malloc(sizeof(spmat));
    sp->n = n;
    sp->add_row = add_row_array;
    sp->free = free_array;
    sp->mult = mult_array;
    sp->private = sparray;
    sp->printSprase = print_array;
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

spmat *readGraph(FILE *input, int type) {
    spmat *graph;
    int i, size, elem, *row, j;
    unsigned int n;
    int nnz;
    if (type == 2) { nnz = find_nnz(input); }
//    printf("nnz is %d\n", nnz);
    n = fread(&elem, sizeof(int), 1, input);
    if (n != 1) {
        printf("ERROR - mismatch reading value");
        exit(EXIT_FAILURE);
    }
    size = elem;
    if (type == 1) {
        graph = spmat_allocate_list(size);
    }
    if (type == 2) {
        graph = spmat_allocate_array(size, nnz);
    }
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
    return graph;
}



