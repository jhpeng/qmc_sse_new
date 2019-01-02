#include <stdio.h>
#include <stdlib.h>

#include "DataTypes.h"

/*********Low level structure and function**********/
IntSequence *CreateIntSequence(int length)
{
    IntSequence* sequence = (IntSequence*)malloc(sizeof(IntSequence));
    sequence->length=length;
    sequence->data = (int*)malloc(length*sizeof(int));

    return sequence;
}

LongIntSequence *CreateLongIntSequence(int length)
{
    LongIntSequence* sequence = (LongIntSequence*)malloc(sizeof(LongIntSequence));
    sequence->length=length;
    sequence->data = (long int*)malloc(length*sizeof(long int));

    return sequence;
}

FloatSequence *CreateFloatSequence(int length)
{
    FloatSequence* sequence = (FloatSequence*)malloc(sizeof(FloatSequence));
    sequence->length=length;
    sequence->data = (float*)malloc(length*sizeof(float));

    return sequence;
}

DoubleSequence *CreateDoubleSequence(int length)
{
    DoubleSequence* sequence = (DoubleSequence*)malloc(sizeof(DoubleSequence));
    sequence->length=length;
    sequence->data = (double*)malloc(length*sizeof(double));

    return sequence;
}


void DestroyIntSequence(IntSequence* sequence)
{
    free(sequence->data);
    free(sequence);
}

void DestroyLongIntSequence(LongIntSequence* sequence)
{
    free(sequence->data);
    free(sequence);
}

void DestroyFloatSequence(FloatSequence* sequence)
{
    free(sequence->data);
    free(sequence);
}

void DestroyDoubleSequence(DoubleSequence* sequence)
{
    free(sequence->data);
    free(sequence);
}

/*********Medium level struct and function***********/
LatticeConf* CreateSquareLatticeConf(const int* shape, int dims)
{
    LatticeConf* lconf = (LatticeConf*)malloc(sizeof(LatticeConf));
    lconf->shape = (int*)malloc(dims*sizeof(int));
    lconf->dims  = dims;
    lconf->nsite = 1;
    for(int i=0;i<dims;++i){
        lconf->shape[i]=shape[i];
        lconf->nsite*=shape[i];
    }
    lconf->Nb=lconf->nsite*dims;
    lconf->sigma0=CreateIntSequence(lconf->nsite);
    lconf->sigmap=CreateIntSequence(lconf->nsite);

    return lconf;
}

OperatorSequence* CreateOperatorSequence(int length, int ndiff)
{
    OperatorSequence* ops = (OperatorSequence*)malloc(sizeof(OperatorSequence));
    ops->length=length;
    ops->noo=0;
    ops->ndiff=ndiff;
    ops->sequence = CreateIntSequence(length);
    ops->sort   = CreateIntSequence(length);
    for(int i=0;i<length;++i){
        ops->sequence->data[i]=-1;
        ops->sort->data[i]=-1;
    }

    return ops;
}

void DestroyLatticeConf(LatticeConf* lconf)
{
    DestroyIntSequence(lconf->sigma0);
    DestroyIntSequence(lconf->sigmap);
    free(lconf->shape);
    free(lconf);
}

void DestroyOperatorSequence(OperatorSequence* ops)
{
    DestroyIntSequence(ops->sequence);
    DestroyIntSequence(ops->sort);
    free(ops);
}

void LatticeConfSetMapping(LatticeConf* lconf, bond2sigma *mapping)
{
    lconf->mapping = mapping;
}

void LatticeConfApplyMapping(LatticeConf* lconf, int bond)
{
    lconf->mapping(&lconf->first,&lconf->second,lconf->shape,bond);
}

void OperatorSequenceSort(OperatorSequence* ops)
{
    int p;
    for(p=0;p<ops->length;++p) ops->sort->data[p]=-1;
    ops->noo=0;
    for(p=0;p<ops->length;++p){
        if(ops->sequence->data[p]!=-1){
            ops->sort->data[ops->noo]=p;
            ++ops->noo;
        }
    }
}

#if 1
int main()
{
    int i,ndiff=2,length=10000;
    int dims=2,shape[2]={32,64};

    for(i=0;i<10000;i++){
        IntSequence* ints = CreateIntSequence(length);
        LongIntSequence* lints = CreateLongIntSequence(length);
        FloatSequence* floats = CreateFloatSequence(length);
        DoubleSequence* doubles = CreateDoubleSequence(length);
        LatticeConf* lconf = CreateSquareLatticeConf(shape, dims);
        OperatorSequence* ops = CreateOperatorSequence(length,ndiff);
    
        OperatorSequenceSort(ops);

        DestroyIntSequence(ints);
        DestroyLongIntSequence(lints);
        DestroyFloatSequence(floats);
        DestroyDoubleSequence(doubles);
        DestroyLatticeConf(lconf);
        DestroyOperatorSequence(ops);
    }
}
#endif
