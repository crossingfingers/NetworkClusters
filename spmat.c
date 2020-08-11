#include <stdlib.h>
#include "spmat.h"


typedef struct _list {
    struct _node *node;
    int k;
}list;

typedef struct _node {
    struct _node *next;
    int col_idx;
    double value;
}node;

typedef struct _array{
    double *values;
    int *colind;
    int *rowptr;
    double *lastvalue;
    int *lastcolind;
    int *lastrowptr;
    int nnz;
}array;

void add_row_list(struct _spmat *A, const double *row, int i) {
    int j;
    list **rows = (list **) A->private;
    node* idx;
    node *temp = NULL;
    node *now = NULL;
    for (j = 0; j < A->n; j++) {
        if (row[j] != 0) {
            if (rows[i] == NULL) {
                rows[i] = malloc(sizeof(list));
//                rows[i]->node = malloc(sizeof(node));
                idx = rows[i]->node;
                rows[i]->node = row[j];
                rows[i]->next = NULL;
                rows[i]->col_idx = j;
                now = rows[i];
            } else {
                temp = malloc(1 * sizeof(node));
                temp->next = NULL;
                temp->value = row[j];
                temp->col_idx = j;
                now->next = temp;
                now = temp;
            }
        }
    }
}

void add_row_array(struct _spmat *A, const double *row, int i){
    double *valarr;
    int *rowptrarr;
    int *colindarr;
    int ci;
    array *sparray= (array *) A->private;
    int temp=*sparray->lastrowptr;
    sparray->lastrowptr++;
    *sparray->lastrowptr=temp;

     valarr=(double *)sparray->lastvalue;
     rowptrarr=(int *)sparray->lastrowptr;
     colindarr=(int *) sparray->lastcolind;

    ci=i;

/*updates values array*/
    for(ci=0;ci<A->n;ci++)
    {   if((*row)!=0)
        {
        *valarr=*row;
        *colindarr=ci;
        (*rowptrarr)++;
        valarr++;
        colindarr++;
        }
        row++;

    }

   sparray->lastvalue= valarr;
   sparray->lastrowptr=rowptrarr;
   sparray->lastcolind=colindarr;
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
    array *sparray= (array *) A->private;

    free(sparray->colind);
    free(sparray->rowptr);
    free(sparray->values);
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

void mult_array(const struct _spmat *A, const double *v, double *result) {
    double sum = 0;
    int vi=0;
    int ri;
    array *sparray = (array *) A->private;
    double *val = (double *) sparray->values;
    int *rowptr = (int *) sparray->rowptr;
    int *colin = (int *) sparray->colind;


    for(ri=1;ri<A->n+1;ri++) {
        sum=0;
       while(vi < rowptr[ri])  {
            sum =sum+ val[vi] * v[colin[vi]];
           vi++;
        }
       result[ri-1]=sum;
    }

}

/*
void print_Struct(struct _spmat *A, int size)
{
    array *sparray= (array *) A->private;
    double *valarr=(double *)sparray->values;
    int *rowptrarr=(int *)sparray->rowptr;
    int *colindarr=(int *) sparray->colind;
    double val;
    int vali;
    int i;
    for (i = 0; i < size; i++) {
        val=*valarr;
        printf("%f\t", val);
        valarr++;
    }
    printf("\n");

    for (i = 0; i < size; i++) {
        vali=*colindarr;
        printf("%d\t", vali);
        colindarr++;
    }
    printf("\n");
    for (i = 0; i < A->n+1; i++) {
        vali=*rowptrarr;
        printf("%d\t", vali);
    rowptrarr++;
    }
    printf("\n");

}
*/
spmat *spmat_allocate_list(int n) {
    spmat *sp;
    list **rows = calloc(n, sizeof(list *));
    sp = malloc(sizeof(spmat));
    sp->n = n;
    sp->add_row = add_row_list;
    sp->free = free_list;
    sp->mult = mult_list;
    sp->private = rows;
    return sp;
}



spmat *spmat_allocate_array(int n, int nnz){
    spmat *sp;
    array *sparray=malloc(sizeof(array));
    sparray->values=calloc(nnz, sizeof(double));
    sparray->colind=calloc(nnz, sizeof(int));
    sparray->rowptr=calloc(n+1, sizeof(int));
    sparray->lastvalue= sparray->values;
    sparray->lastcolind= sparray->colind;
    sparray->lastrowptr=  sparray->rowptr;
    sparray->nnz= nnz;
    sp = malloc(sizeof(spmat));
    sp->n = n;
    sp->add_row = add_row_array;
    sp->free = free_array;
    sp->mult = mult_array;
    sp->private = sparray;
    return sp;
}


