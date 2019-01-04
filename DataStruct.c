#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>

#include "DataStruct.h"

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
    lconf->first =CreateIntSequence(lconf->nsite);
    lconf->last  =CreateIntSequence(lconf->nsite);
    lconf->flip  =CreateIntSequence(lconf->nsite);

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

OperatorLoop* CreateOperatorLoop(int length)
{
    OperatorLoop* opl = (OperatorLoop*)malloc(sizeof(OperatorLoop));
    opl->length=length;
    opl->noo=0;
    opl->check=1;
    opl->loop = CreateIntSequence(2*length);
    opl->lvc  = CreateIntSequence(4*length);
    opl->cross= CreateIntSequence(4*length);
    for(int i=0;i<4*length;++i){
         opl->lvc->data[i]=-1;
         opl->cross->data[i]=-1;
    }
    for(int i=0;i<2*length;++i) opl->loop->data[i]=-1;

    return opl;
}

void DestroyLatticeConf(LatticeConf* lconf)
{
    DestroyIntSequence(lconf->sigma0);
    DestroyIntSequence(lconf->sigmap);
    DestroyIntSequence(lconf->first);
    DestroyIntSequence(lconf->last);
    DestroyIntSequence(lconf->flip);
    free(lconf->shape);
    free(lconf);
}

void DestroyOperatorSequence(OperatorSequence* ops)
{
    DestroyIntSequence(ops->sequence);
    DestroyIntSequence(ops->sort);
    free(ops);
}

void DestroyOperatorLoop(OperatorLoop* opl)
{
    DestroyIntSequence(opl->lvc);
    DestroyIntSequence(opl->loop);
    DestroyIntSequence(opl->cross);
    free(opl);
}


/*********High level struct and function***********/
SEPlaceHolder* CreateSEPlaceHolder()
{
    SEPlaceHolder* placeholder = (SEPlaceHolder*)malloc(sizeof(SEPlaceHolder));
    placeholder->initialize=0;
    placeholder->length=0;
    placeholder->isweep=0;
    placeholder->nsweep=0;
    placeholder->cutoff=0;
    placeholder->seed=0;
    placeholder->set_lattice=0;
    placeholder->set_length=0;
    placeholder->set_random=0;
    placeholder->beta=-1;
    placeholder->max_err=-1;
    return placeholder;
}

void DestroySEPlaceHolder(SEPlaceHolder* placeholder)
{
    if(placeholder->set_lattice) DestroyLatticeConf(placeholder->lconf);
    if(placeholder->set_length){
        DestroyOperatorSequence(placeholder->ops);
        DestroyOperatorLoop(placeholder->opl);
    }
    if(placeholder->set_random) gsl_rng_free(placeholder->rng);
    free(placeholder);
}

#if 0
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
        LinkedVertexConf* lvc = CreateLinkedVertexConf(length);
    

        DestroyIntSequence(ints);
        DestroyLongIntSequence(lints);
        DestroyFloatSequence(floats);
        DestroyDoubleSequence(doubles);
        DestroyLatticeConf(lconf);
        DestroyOperatorSequence(ops);
        DestroyLinkedVertexConf(lvc);
    }
}
#endif
