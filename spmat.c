#include <stdlib.h>
#include <stdio.h>
#include "spmat.h"


typedef struct _list {
    struct _node *node;
    int k;
} list;

typedef struct _node {
    struct _node *next;
    int col_idx;
    int value;
} node;

//typedef struct _array{
//    double *values;
//    int *colind;
//    int *rowptr;
//    double *lastvalue;
//    int *lastcolind;
//    int *lastrowptr;
//    int nnz;
//}array;

void add_row_list(struct _spmat *A, const int *row, int i, int k) {
    int j;
    list **rows = (list **) A->private;
    node *idx;
    node *temp = NULL;
    node *now = NULL;
    A->M += k;
    if (k == 0) {
        rows[i] = malloc(sizeof(list));
        rows[i]->k = k;
        rows[i]->node = NULL;
    }
    for (j = 0; j < k; j++) {
        if (rows[i] == NULL) {
            rows[i] = malloc(sizeof(list));
            rows[i]->k = k;
            rows[i]->node = malloc(sizeof(node));
            idx = rows[i]->node;
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

//void add_row_array(struct _spmat *A, const double *row, int i){
//    double *valarr;
//    int *rowptrarr;
//    int *colindarr;
//    int ci;
//    array *sparray= (array *) A->private;
//    int temp=*sparray->lastrowptr;
//    sparray->lastrowptr++;
//    *sparray->lastrowptr=temp;
//
//     valarr=(double *)sparray->lastvalue;
//     rowptrarr=(int *)sparray->lastrowptr;
//     colindarr=(int *) sparray->lastcolind;
//
//    ci=i;
//
///*updates values array*/
//    for(ci=0;ci<A->n;ci++)
//    {   if((*row)!=0)
//        {
//        *valarr=*row;
//        *colindarr=ci;
//        (*rowptrarr)++;
//        valarr++;
//        colindarr++;
//        }
//        row++;
//
//    }
//
//   sparray->lastvalue= valarr;
//   sparray->lastrowptr=rowptrarr;
//   sparray->lastcolind=colindarr;
//}

void print_list(struct _spmat *A) {
    int i;
    int j;
    node *curr;
    list **rows = (list **) A->private;
    for (i = 0; i < A->n; ++i) {
        printf("%d - \t", i);
        curr = rows[i]->node;
        for (j = 0; j < rows[i]->k; ++j) {
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

//void free_array(struct _spmat *A) {
//    array *sparray= (array *) A->private;
//
//    free(sparray->colind);
//    free(sparray->rowptr);
//    free(sparray->values);
//    free(A->private);
//
//}

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

//void mult_array(const struct _spmat *A, const double *v, double *result) {
//    double sum = 0;
//    int vi=0;
//    int ri;
//    array *sparray = (array *) A->private;
//    double *val = (double *) sparray->values;
//    int *rowptr = (int *) sparray->rowptr;
//    int *colin = (int *) sparray->colind;
//
//
//    for(ri=1;ri<A->n+1;ri++) {
//        sum=0;
//       while(vi < rowptr[ri])  {
//            sum =sum+ val[vi] * v[colin[vi]];
//           vi++;
//        }
//       result[ri-1]=sum;
//    }
//
//}

int findNextRowList(spmat *A, int i) {
    list **rows = (list **) A->private;
    int k = rows[i]->k;
    while (k == 0 && i < A->n) {
        k = rows[++i]->k;
    }
    if (i < A->n) {
        return i;
    }
    return -1;
    spIterator *iter = A->iter;
    iter->rowIdx = i;
    iter->colIdx = rows[i]->node->col_idx;
    iter->kRow = k;
    iter->kCol = rows[iter->colIdx]->k;
    iter->private = (void *) rows[i]->node;
}

void setIterRowList(spmat *A, int i) {
    list **rows = (list **) A->private;
    spIterator *iter = A->iter;
    iter->rowIdx = i;
    iter->colIdx = rows[i]->node->col_idx;
    iter->kRow = rows[i]->k;
    iter->kCol = rows[iter->colIdx]->k;
    iter->private = (void *) rows[i]->node;
}

void initListIter(spmat *A) {
    setIterRowList(A, 0);
}

int listHasNext(spmat *A) {
    int i;
    node *curr = (node *) A->iter->private;
    if (curr->next != NULL) {
        return 1;
    }
    i = A->iter->rowIdx + 1;
    if (i < A->n) {
        if (findNextRowList(A, i) > 0)
            return 1;
    }
    return 0;
}

void iterListNext(spmat *A) {
    int i;
    list **rows = (list **) A->private;
    node *curr = (node *)A->iter->private;
    if (curr->next != NULL) {
        A->iter->private = (void *)curr->next;
        A->iter->colIdx = curr->next->col_idx;
        A->iter->kCol = rows[A->iter->colIdx]->k;
    }
    i = A->iter->rowIdx + 1;
    setIterRowList(A, findNextRowList(A, i));
}

spmat *spmat_allocate_list(int n) {
    spmat *sp;
    list **rows = calloc(n, sizeof(list *));
    sp = malloc(sizeof(spmat));
    sp->n = n;
    sp->M = 0;
    sp->add_row = add_row_list;
    sp->free = free_list;
    sp->mult = mult_list;
    sp->private = rows;
    sp->print_list = print_list;
    sp->iter = malloc(sizeof(spIterator));
    sp->initIterator = initListIter;
    sp->hasNext = listHasNext;
    sp->iterNext = iterListNext;
    return sp;
}



//spmat *spmat_allocate_array(int n, int nnz){
//    spmat *sp;
//    array *sparray=malloc(sizeof(array));
//    sparray->values=calloc(nnz, sizeof(double));
//    sparray->colind=calloc(nnz, sizeof(int));
//    sparray->rowptr=calloc(n+1, sizeof(int));
//    sparray->lastvalue= sparray->values;
//    sparray->lastcolind= sparray->colind;
//    sparray->lastrowptr=  sparray->rowptr;
//    sparray->nnz= nnz;
//    sp = malloc(sizeof(spmat));
//    sp->n = n;
//    sp->add_row = add_row_array;
//    sp->free = free_array;
//    sp->mult = mult_array;
//    sp->private = sparray;
//    return sp;
//}


