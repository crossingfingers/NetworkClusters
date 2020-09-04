//
// Created by gal21 on 29/08/2020.
//

#include "stack.h"
#include <stdlib.h>
#include<stdio.h>

int isEmpty(stack *s){
    return s->size == 0;
}

stackNode* pop(stack *s){
    stackNode *tmp;
    if (s->size > 0){
        tmp = s->node;
        s->node = s->node->next;
        s->size--;
        return tmp;
    }
    return NULL;
}

int push(stack *s, stackNode *n){
    if(s->size > 0){
        n->next = s->node;
    }
    s->node = n;
    s->size++;
    return 1;
}

stack *initStack(){
    stack *s = malloc(sizeof(stack));
    if(s == NULL){
        printf("ERROR - memory allocation unsuccessful");
        exit(EXIT_FAILURE);
    }
    s->size=0;
    s->node=NULL;
    s->isEmpty = isEmpty;
    s->pop = pop;
    s->push = push;
    return s;
}
