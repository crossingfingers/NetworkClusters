#include <stdlib.h>
#include <stdio.h>
#include "spmat.h"
#include <math.h>

typedef struct _node {
    struct _node *next;
    int col_idx;
    int value;
} node;

typedef struct _array {
    int *values;
    int *colind;
    int *rowptr;
    int lastindex;
    int nnz;
} array;

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

void add_row_array(struct _spmat *A, int *row, int i, int k) {
    if (k) {
        array *sparray = (array *) A->private;
        int firstrowidx;
        int ci;
        A->M += k;
        A->k[i] = k;
        if (!sparray->lastindex) { firstrowidx = 0; }
        else { firstrowidx = (sparray->rowptr[sparray->lastindex - 1]) + 1; }

/*updates values array*/
        for (ci = 0; ci < k; ci++) {
            sparray->values[sparray->lastindex] = 1;
            sparray->colind[sparray->lastindex] = *row;
            sparray->rowptr[sparray->lastindex] = firstrowidx;
            sparray->lastindex++;
            row++;
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

}

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

void free_array(struct _spmat *A) {
    array *sparray = (array *) A->private;
    free(A->k);
    free(sparray->values);
    free(sparray->rowptr);
    free(sparray->colind);
    free(A->private);

}

void mult_list(const struct _spmat *A, const double *v, double *result) {
    double res;
    int i;
    node *now;
    node **rows = (node **) A->private;
    for (i = 0; i < A->n; i++) {
        res = 0;
        now = rows[i];
        while (now != NULL) {
            res += now->value * v[now->col_idx];
            now = now->next;
        }
        result[i] = res;
    }
}

void mult_array(const struct _spmat *A, const double *vec, double *result) {
    array *sparray = (array *) A->private;
  int i=0;
  for(i=0;i<sparray->nnz;i++)
      result[sparray->rowptr[i]]+=sparray->values[i]*vec[sparray->colind[i]];
}

double listShifting(spmat *A, int group, const int *groupid) {
    double max = 0;
    double sum = 0;
    int val = 0;
    node *curr;
    node **rows = (node **) A->private;
    int i;
    int j;
    for (i = 0; i < A->n; ++i) {
        if (group != groupid[i])
            continue;
        curr = rows[i];
        sum = 0;
        for (j = 0; j < A->n; ++j) {
            if (curr != NULL && curr->col_idx == j) {
                val = 1;
                curr = curr->next;
            } else
                val = 0;
            if (groupid[j] == group)
                sum += fabs((double) val - ((double) (A->k[i] * A->k[j]) / A->M));
        }
        max = (max >= sum) ? max : sum;
    }
    return max;
}

double arrayShifting(spmat *A, int group, const int *groupid) {
    double max = 0;
    double sum;
    int counter = 0;
    int val = 0;
    int idx;
    array *sparray = (array *) A->private;
    int i;
    int j;
    int rowIDX=1;
    int vecIDX=0;
    for (i = 0; i < A->n; ++i) {
        sum = 0;
        rowIDX++;
        for (j = 0; j < A->n; ++j) {
            val = 0;
                if((sparray->colind[vecIDX]==j)&&(sparray->rowptr[vecIDX]<sparray->rowptr[rowIDX]))
                { val = 1;
                    vecIDX++;
                }
                sum += fabs((double) val - ((double) (A->k[i] * A->k[j]) / A->M));


        }
        max = (max >= sum) ? max : sum;
    }
    printf("\nmax shifting is %f\n", max);
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
    sp->print_list = print_list;
    sp->matShifting = listShifting;
    return sp;
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
    sp->M = 0;
    sp->k = malloc(sizeof(int) * n);
    initk(sp);
    sp->print_list = print_array;
    sp->matShifting = arrayShifting;
    return sp;
}

