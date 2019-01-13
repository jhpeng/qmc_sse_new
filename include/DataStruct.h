#include <time.h>
#include <gsl/gsl_rng.h>

#ifndef DATASTRUCT_H
#define DATASTRUCT_H


/*********Low level structure and function**********/
typedef struct IntSequence{
    int length;
    int *data;
} IntSequence;

typedef struct LongIntSequence{
    int length;
    long int *data;
} LongIntSequence;

typedef struct FloatSequence{
    int length;
    float *data;
} FloatSequence;

typedef struct DoubleSequence{
    int length;
    double *data;
} DoubleSequence;

IntSequence *CreateIntSequence(int length);

LongIntSequence *CreateLongIntSequence(int length);

FloatSequence *CreateFloatSequence(int length);

DoubleSequence *CreateDoubleSequence(int length);

void DestroyIntSequence(IntSequence* sequence);

void DestroyLongIntSequence(LongIntSequence* sequence);

void DestroyFloatSequence(FloatSequence* sequence);

void DestroyDoubleSequence(DoubleSequence* sequence);


/*********Medium level struct and function***********/
//(&left,&right,shape,bond)
typedef void bond2sigma(int*, int*, const int*, int);

//the configuration of the lattice model
//dims : the dimension of the lattice model
//nsite : the whole system size 
//Nb : the total number of bond (nsite*dims for square lattice)
//left and right : the output of function mapping
//shape : the shape of lattice model nsite=shape[0]*shape[1]*...
//sigma0 : the spin configuration for | \alpha(0) > state
//sigmap : the spin configuration for | \alpha(p) > propagate state
//mapping : the mapping from bond to 2 leg
typedef struct LatticeConf{
    int dims;
    int nsite;
    int nflip;
    int Nb;
    int left;
    int right;
    int *shape;
    IntSequence *sigma0;
    IntSequence *sigmap;
    IntSequence *first;
    IntSequence *last;
    IntSequence *flip;
    DoubleSequence *J;
    bond2sigma *mapping;
} LatticeConf;

//the sequence of the operator H_{a(p),b(p)}
//SL={ [a(p),b(p)],[a(p+1),b(p+1)],... }
//p from 0 to length-1
//a(p)=0 is diagonal operator; a(p)=1 is off-diagonal operator; a(p)=2 is another
//a(p)=sequence->data[p]%ndiff
//b(p) is the operator on ID b(p) bond
//b(p)=sequence->data[p]/ndiff
typedef struct OperatorSequence{
    int length;
    int noo;
    int ndiff;
    IntSequence *sequence;
    IntSequence *sort;
} OperatorSequence;

typedef struct OperatorLoop{
    int length;
    int noo;
    int check;
    int v0;
    IntSequence *lvc;
    IntSequence *cross;
    IntSequence *loop;
} OperatorLoop;

LatticeConf* CreateSquareLatticeConf(const int* shape, int dims);

OperatorSequence* CreateOperatorSequence(int length, int ndiff);

OperatorLoop* CreateOperatorLoop(int length);

void DestroyLatticeConf(LatticeConf* lconf);

void DestroyOperatorSequence(OperatorSequence* ops);

void DestroyOperatorLoop(OperatorLoop* opl);


/*********High level struct and function***********/
typedef struct SEPlaceHolder{
    int initialize;
    int length;
    int isweep;
    int nsweep;
    int cutoff;
    int seed;
    int set_lattice;
    int set_length;
    int set_random;
    double beta;
    double max_err;
    gsl_rng* rng;
    LatticeConf* lconf;
    OperatorSequence* ops;
    OperatorLoop* opl;
} SEPlaceHolder;

typedef double measurement(SEPlaceHolder*, void* args);

typedef struct Observable{
    int nobs;
    int nave;
    int count;
    int set_obs;
    time_t start;
    time_t end;
    char** obs_name;
    measurement** measure;
    void** args;
    DoubleSequence* data;
    double* mean;
    double* var;
    double* err;
} Observable;

SEPlaceHolder* CreateSEPlaceHolder();

void DestroySEPlaceHolder(SEPlaceHolder* placeholder);
#endif
