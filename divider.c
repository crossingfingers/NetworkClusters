#include <stdio.h>
#include <stdlib.h>
#include "spmat.h"
#include <assert.h>
#include <math.h>
#include "divider.h"
#include <time.h>

#define IS_POSITIVE(X) ((X) > 0.00001)
#define IS_NEGATIVE(X) ((X) < -0.00001)

void randomizeVec(int size, double *vec) {
    int i;
    srand(time(NULL));
    assert(vec != NULL);
    for (i = 0; i < size; i++) {
        vec[i] = rand();
    }
}


void normalize(int size, double *vec) {
    int i;
    double res = 0;
    for (i = 0; i < size; i++) {
        res += vec[i] * vec[i];
    }
    res = sqrt(res);
    for (i = 0; i < size; i++) {
        vec[i] /= res;
    }
}

void vecMult(const int *vec1, const double *vec2, double *res, int size) {
    int i;
    for (i = 0; i < size; ++i) {
//        printf("%f\n", vec1[i]);
        res[i] = vec1[i] * vec2[i];
    }
}

void vecSum(double *vec, const double *b0, double shifting, int n) {
    int i;
    for (i = 0; i < n; ++i) {
        vec[i] += b0[i] * shifting;
    }
}

void vecDec(double *vec1, const double *vec2, int group, const int *groupid, int n) {
    int i;
    for (i = 0; i < n; ++i) {
        if (group == groupid[i])
            vec1[i] -= vec2[i];
        else
            vec1[i] = 0;
    }
}

void scalarMult(double *vec, double x, int group, const int *groupid, int n) {
    int i;
    for (i = 0; i < n; ++i) {
        if (group == groupid[i])
            vec[i] *= x;
        else
            vec[i] = 0;
    }
}

double dotProd(const int *vec1, const double *vec2, int group, const int *groupid, int n) {
    double res = 0;
    int i;
    for (i = 0; i < n; ++i) {
        if (groupid[i] == group)
            res += vec1[i] * vec2[i];
    }
    return res;
}

double dotDoubleProd(const double *vec1, const double *vec2, int group, const int *groupid, int n) {
    double res = 0;
    int i;
    for (i = 0; i < n; ++i) {
        if (group == groupid[i])
            res += vec1[i] * vec2[i];
    }
    return res;
}

void copyVec(const int *src, double *dst, int n) {
    int i;
    for (i = 0; i < n; ++i) {
        dst[i] = src[i];
    }
}

void multBv(spmat *sp, double *vec, int group, const int *groupid, double *res) {
    double dot;
    int size = sp->n;
    double *res1 = malloc(size * sizeof(double));
    dot = dotProd(sp->k, vec, group, groupid, size);
    if (sp->M == 0) {
        printf("ERROR - divide in zero");
        exit(EXIT_FAILURE);
    }
    copyVec(sp->k, res1, size);
    scalarMult(res1, (double) dot / sp->M, group, groupid, size);
    sp->mult(sp, vec, res);
    vecDec(res, res1, group, groupid, size);
    free(res1);
}

void initUnitVec(double *unitVec, int n) {
    int i;
    for (i = 0; i < n; ++i) {
        unitVec[i] = 1;
    }
}

void printVector(double *vec, int n) {
    int i;
    for (i = 0; i < n; ++i) {
        printf("%f\t", vec[i]);
    }
    printf("\n");
}

void powerIter(spmat *sp, double *b0, double shifting, int group, const int *groupid, double *result) {
    int flag = 1, i;
    int size = sp->n;
    double *unitVec = malloc(size * sizeof(double));
    double *tmp = malloc(size * sizeof(double));
    if (unitVec == NULL || tmp == NULL) {
        printf("ERROR - memory allocation unsuccessful");
        exit(EXIT_FAILURE);
    }
    initUnitVec(unitVec, size);
    multBv(sp, unitVec, group, groupid, tmp);
    while (flag == 1) {
        flag = 0;
        multBv(sp, b0, group, groupid, result);
        vecSum(result, b0, shifting, size);
        vecDec(result, tmp, group, groupid, size);
        normalize(size, result);
        for (i = 0; i < size; i++) {
            if (IS_POSITIVE(fabs(result[i] - b0[i])))
                flag = 1;
            b0[i] = result[i];
        }
    }
    free(tmp);
    free(unitVec);
}

double modularityCalc(spmat *sp, double *vec, int group, const int *groupid) {
    double res = 0;
    int size = sp->n;
    double *tmp = malloc(sizeof(double) * size);
    multBv(sp, vec, group, groupid, tmp);
    res = dotDoubleProd(tmp, vec, group, groupid, size);
    free(tmp);
    return res / 2;
}

void split(struct _division *d, spmat *sp, double *vec, int group) {
    int flag;
    double delta;
    int newGroup = -1;
    int i;
    flag = IS_POSITIVE(vec[0]) ? 1 : 0;
    vec[0] = 1;
    for (i = 1; i < sp->n; ++i) {
        if (group != d->groupid[i])
            continue;
        if (IS_POSITIVE(vec[i]) != flag) {
            if (newGroup == -1) {
                newGroup = d->numOfGroups;
                d->numOfGroups += 1;
            }
            d->groupid[i] = newGroup;
            vec[i] = -1;
        } else {
            vec[i] = 1;
        }
    }
    delta = modularityCalc(sp, vec, group, d->groupid);
    d->Q += delta;
}

void printGroups(division *d) {
    int i;
    printf("number of groups = %d\n", d->numOfGroups);
    for (i = 0; i < d->n; ++i) {
        printf("(%d,%d)\t", i, d->groupid[i]);
    }
    printf("\nmodularity value: %f\n", d->Q);
}

void freeDivision(division *d) {
    free(d->groupid);
}

division *allocateDivision(int n) {
    int i;
    division *d = malloc(sizeof(division));
    if (d == NULL) {
        printf("ERROR - memory allocation unsuccessful");
        exit(EXIT_FAILURE);
    }
    d->n = n;
    d->printGroups = printGroups;
    d->split = split;
    d->free = freeDivision;
    d->groupid = malloc(sizeof(int) * n);
    d->Q = 0;
    d->divOptimization=divOptimization;
    d->modularityCalc=modularityCalc;
    if (d->groupid == NULL) {
        printf("ERROR - memory allocation unsuccessful");
        exit(EXIT_FAILURE);
    }
    d->numOfGroups = 1;
    for (i = 0; i < n; ++i) {
        d->groupid[i] = 0;
    }
    return d;
}

void printdVector(double *vec, int n) {
    int i;
    for (i = 0; i < n; ++i) {
        printf("%f\t", vec[i]);
    }
    printf("\n");
}
int getDivSize(division *div, int group) {
    int i;
    int counter=0;

    for ( i = 0;i<div->n; i++)
    {
        if(div->groupid[i]==group) {counter++;}

    }  return counter;

}


/*find vertice movement that maximizes Q*/
int moveVertice(double q0,double *divVec,spmat *sp, const int *unmoved, double *maxDeltaQ,int group,const int *groupID,int size)
{

    int j;
    double deltaQ;
    int maxIndex;
    int flag=1;

    for(j=0;j<sp->n;j++)
    {   if(group==groupID[j])
        {
            if(unmoved[j])
            {

                divVec[j]= (-1*divVec[j]);  /*moves vertice j */
                deltaQ = modularityCalc(sp,divVec,group,groupID)-q0;  /*calcs new deltaQ */

                divVec[j]= (-1*divVec[j]);  /*moves back vertice j */

                if(flag){*maxDeltaQ=deltaQ;maxIndex=j;flag=0;  }

                if(deltaQ>=*maxDeltaQ){*maxDeltaQ=deltaQ;maxIndex=j;} /*keeps track of max Q */

            }
        }
    }

    return maxIndex;

}


void optimize(double q0, double *divVec,double *deltaQ,int *unmoved,int *indices, double *improve, spmat*sp,int group,int *groupID,int size)
{   int j;
    int i;
    int counter=0;
    int maxIndex;
    int maxScore;
    for(i=0;i<size;i++){unmoved[i]=1;}  /*keeps track who moved */


    for(j=0;j<sp->n;j++)
    {   if(groupID[j]==group){
        *deltaQ=0;
        maxIndex=moveVertice(q0,divVec,sp,unmoved,deltaQ,group,groupID,size);     /*finds vertice that maximizes deltaQ*/
        unmoved[maxIndex]=0;    /*moves vertice*/
        divVec[maxIndex]=(-1)*divVec[maxIndex];  /*moves vertice*/
        indices[counter]=maxIndex;


        if(counter==0){improve[counter]=(*deltaQ); maxScore=counter;}
        else { improve[counter] = (*deltaQ);}
        if (improve[counter] > improve[maxScore]) { maxScore = counter; }
            /*saves max division state */

            counter++;

    }

    }

    if(IS_POSITIVE(improve[maxScore]))
    {
        for(j=(size-1);j>maxScore;j--){ divVec[indices[j]]=(-1) * divVec[indices[j]];}
        optimize(improve[maxScore]+q0,divVec,deltaQ,unmoved,indices,improve,sp,group,groupID,size);
    }  /*finds max division state, if Q is positive, we try again */

    else {for(j=0;j<size;j++){divVec[indices[j]]=(-1)*divVec[indices[j]];}}



}
/*optimizes division by moving one node from g1 to g2, saves division in res*/

void divOptimization(division *div,int group,double q0,double *initialdiv, spmat *sp)
{
    int size=getDivSize(div,group);
    int *unmoved=malloc(sizeof(int)*div->n);
    double *deltaQ=malloc(sizeof(double)); /*DeltaQ result*/
    int *indices=malloc(sizeof(int)*size);
    double *improve=malloc(sizeof(double)*size);
    optimize(q0,initialdiv,deltaQ,unmoved,indices,improve,sp,group,div->groupid,size);
    div->split(div,sp,initialdiv,1);
    div->printGroups(div);
    free(deltaQ);
    free(unmoved);
    free(indices);
    free(improve);
 //   div->split(div,sp,initialdiv,group);
 //  div->printGroups(div);
}



