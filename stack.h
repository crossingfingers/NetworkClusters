//
// Created by gal21 on 29/08/2020.
//

#ifndef CPROJECT_STACK_H
#define CPROJECT_STACK_H

typedef struct _stackNode{
    int start;
    int end;
    struct _stackNode *next;
}stackNode;

typedef struct _stack {
    stackNode *node;
    int size;
    int (*isEmpty)(struct _stack *s);
    stackNode* (*pop)(struct _stack *s);
    int (*push)(struct _stack *s, stackNode *n);
}stack;

stack* initStack();

#endif //CPROJECT_STACK_H
