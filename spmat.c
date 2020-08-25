#include <stdlib.h>
#include <stdio.h>
#include "spmat.h"
#include <math.h>

typedef struct _node {
    struct _node *next;
    int col_idx;
    int value;
} node;

typedef struct _array{
    int *values;
    int *colind;
    int *rowptr;
    int *lastvalue;
    int *lastcolind;
    int *lastrowptr;
    int lastindex;
    int nnz;
}array;

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

/*works*/
void add_row_array(struct _spmat *A,  int *row, int i,int k){
    array *sparray= (array *) A->private;
    int *valarr=(int *)sparray->lastvalue;
    int *rowptrarr=(int *)sparray->lastrowptr;
    int *colindarr=(int *) sparray->lastcolind;
    int firstrowidx;
    int ci;
    A->M += k;
    A->k[i] = k;
    if(!sparray->lastindex){firstrowidx=0;}
    else{firstrowidx=(sparray->rowptr[sparray->lastindex-1])+1; }

/*updates values array*/
    for(ci=0;ci<k;ci++)
    {

        *valarr=1;
        *colindarr=*row;
        *rowptrarr=firstrowidx;


        valarr++;
        colindarr++;
        rowptrarr++;
        row++;
        sparray->lastindex++;

    }
   sparray->lastvalue= valarr;
   sparray->lastcolind=colindarr;
    sparray->lastrowptr=rowptrarr;
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

/*works*/
void print_array(struct _spmat *A) {
    int i;
    int j;
    array *sparray= (array *) A->private;
    double *valarr=(double *)sparray->values;
    int *rowptrarr=(int *)sparray->rowptr;
    int *colindarr=(int *) sparray->colind;
    int nnz=(int)sparray->nnz;

    for (i = 0; i < A->n; ++i) {
        printf("%d - \t", i);

        for (j = 0; j < A->k[i]; ++j) {
            printf("%d\t", *colindarr);
           colindarr++;

        }
        printf("\n");
    }
    printf("nnz val is %d",*(++rowptrarr));
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

/*works*/
void free_array(struct _spmat *A) {
    array *sparray= (array *) A->private;
   free(sparray->colind);
   free(sparray->rowptr);
   free(sparray->values);
   free(A->private);
   free(A->k);
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

double listShifting(spmat *A, int group, const int *groupid){
    double max = 0;
    double sum = 0;
    int val=0;
    node *curr;
    node **rows = (node **) A->private;
    int i;
    int j;
    for(i=0; i<A->n; ++i){
        if(group != groupid[i])
            continue;
        curr = rows[i];
        sum = 0;
        for(j = 0; j< A->n; ++j){
            if(curr != NULL && curr->col_idx == j){
                val = 1;
                curr = curr->next;
            }
            else
                val = 0;
            if(groupid[j] == group)
                sum += fabs((double)val - ((double)(A->k[i] * A->k[j])/A->M));
        }
        max = (max>=sum)? max:sum;
    }
    return max;
}



void mult_array(const struct _spmat *A, const double *v, double *result) {
    double sum = 0;
    int vi=0;
    int ri;
    array *sparray = (array *) A->private;
    double *val = (int *) sparray->values;
    int *rowptr = (int *) sparray->rowptr;
    int *colin = (int *) sparray->colind;

    for(ri=1;ri<A->n+1;ri++) {
        sum=0;
        while(vi < rowptr[ri])  {
            sum =+ (val[vi] * v[colin[vi]]);
            vi++;
        }
        result[ri-1]=sum;
    }

}
/*
need to fix
*/
double arrayShifting(spmat *A, int group, const int *groupid){
    double max = 0;
    double sum = 0;
    int val=0;
    node *curr;
    node **rows = (node **) A->private;
    int i;
    int j;
    for(i=0; i<A->n; ++i){
        if(group != groupid[i])
            continue;
        curr = rows[i];
        sum = 0;
        for(j = 0; j< A->n; ++j){
            if(curr != NULL && curr->col_idx == j){
                val = 1;
                curr = curr->next;
            }
            else
                val = 0;
            if(groupid[j] == group)
                sum += fabs((double)val - ((double)(A->k[i] * A->k[j])/A->M));
        }
        max = (max>=sum)? max:sum;
    }
    return max;
}


void initk(spmat* A){
    int i;
    for(i = 0; i<A->n; ++i){
        A->k[i] = 0;
    }
}

void printMatrix(spmat *A){

    int i;
    array *sparr=( array *) A->private;
    printf("\nvalarr:\n");
    for( i=0;i< sparr->nnz;i++)
    {  printf("%d",*sparr->values);
        sparr->values++;
    }
    printf("\ncollarr:\n");
    for( i=0;i< sparr->nnz;i++)
    {  printf("%d",*sparr->colind);
        sparr->colind++;
    }
    printf("\nrowptrarr:\n");
    for( i=0;i<= sparr->nnz;i++)
    {  printf("%d",*sparr->rowptr);
        sparr->rowptr++;
    }

}

spmat *spmat_allocate_list(int n) {
    spmat *sp;
    node **rows = calloc(n, sizeof(node *));
    sp = malloc(sizeof(spmat));
    sp->n = n;
    sp->M = 0;
    sp->k = malloc(sizeof(int)*n);
    initk(sp);
    sp->add_row = add_row_list;
    sp->free = free_list;
    sp->mult = mult_list;
    sp->private = rows;
    sp->print_list = print_list;
    sp->matShifting = listShifting;
    return sp;
}



spmat *spmat_allocate_array(int n, int nnz){
    spmat *sp;
    array *sparray=malloc(sizeof(array));
    sparray->values=calloc(nnz, sizeof(int));
    sparray->colind=calloc(nnz, sizeof(int));
    sparray->rowptr=calloc(nnz+1, sizeof(int));
    sparray->lastvalue= sparray->values;
    sparray->lastcolind= sparray->colind;
    sparray->lastrowptr=  sparray->rowptr;
    sparray->lastindex=0;
    sparray->rowptr[nnz]=nnz;
    sparray->nnz= nnz;
    sp = malloc(sizeof(spmat));
    sp->n = n;
    sp->add_row = add_row_array;
    sp->free = free_array;
    sp->mult = mult_array;
    sp->private = sparray;
    sp->M = 0;
    sp->k = malloc(sizeof(int)*n);
    initk(sp);
    sp->print_list = print_array;
    sp->matShifting = arrayShifting;
    return sp;
}


