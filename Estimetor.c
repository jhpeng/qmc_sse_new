#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "DataStruct.h"
#include "SEAlgorithm.h"

Observable* CreateObservable(int nobs, int nave)
{
    Observable* Obs = (Observable*)malloc(sizeof(Observable));
    Obs->nobs=nobs;
    Obs->nave=nave;
    Obs->count=0;
    Obs->set_obs=0;
    Obs->start=clock();
    Obs->obs_name = (char**)malloc(nobs*sizeof(char*));
    for(int i=0;i<nobs;i++) Obs->obs_name[i]=(char*)malloc(128*sizeof(char));

    Obs->measure = (measurement**)malloc(nobs*sizeof(measurement*));
    Obs->args = (void**)malloc(nobs*sizeof(void*));
    Obs->data = CreateDoubleSequence(nobs*nave);
    Obs->mean = (double*)malloc(nobs*sizeof(double));
    Obs->var  = (double*)malloc(nobs*sizeof(double));
    Obs->err  = (double*)malloc(nobs*sizeof(double));

    return Obs;
}

void DestroyObservable(Observable* Obs)
{
    free(Obs->mean);
    free(Obs->var);
    free(Obs->err);
    free(Obs->measure);
    free(Obs->args);
    for(int i=0;i<Obs->nobs;i++) free(Obs->obs_name[i]);
    free(Obs->obs_name);
    DestroyDoubleSequence(Obs->data);
    free(Obs);
}

void MeanAverage(Observable* Obs)
{
    int i_obs,i,nave;
    int nobs=Obs->nobs;
    double data;
    if(Obs->count<Obs->nave) nave=Obs->count;
    else nave=Obs->nave;

    for(i_obs=0;i_obs<nobs;i_obs++){
        Obs->mean[i_obs]=0;
        Obs->var[i_obs] =0;
        Obs->err[i_obs] =0;
        for(i=0;i<nave;++i){
            data = Obs->data->data[i*nobs+i_obs];
            Obs->mean[i_obs]+=data;
            Obs->var[i_obs]+=data*data;
        }
        Obs->mean[i_obs]=Obs->mean[i_obs]/nave;
        Obs->var[i_obs]=Obs->var[i_obs]/nave-Obs->mean[i_obs]*Obs->mean[i_obs];
        Obs->err[i_obs]=sqrt(Obs->var[i_obs]/nave);
    }
}

void ObservableSetMeasurement(
                    Observable* obs, 
                    measurement* measure, 
                    char* obs_name,
                    void* args)
{
    obs->measure[obs->set_obs] = measure;
    strcpy(obs->obs_name[obs->set_obs],obs_name);
    obs->args[obs->set_obs]=args;
    obs->set_obs++;
}

void ObservableDoMeasurement(
                    Observable* obs,
                    SEPlaceHolder* placeholder)
{
    int count;
    double sample;
    for(int i_obs=0;i_obs<obs->nobs;i_obs++){
        sample = obs->measure[i_obs](placeholder,obs->args[i_obs]);
        count = obs->count%obs->nave;
        obs->data->data[count*obs->nobs+i_obs] = sample;
    }
    obs->count++;
}

void ObservableShow(
                    Observable* obs,
                    SEPlaceHolder* placeholder,
                    int mode)
{
    MeanAverage(obs);
    //if model=0 : std output
    //if model=1 : output a text file
    //if model=2 : output a html file

    switch(mode){
        case 0:
            printf("Measurement ...\n");
            printf("length of sequence \t: %d \nnumber of operator \t: %d\n",placeholder->ops->length,placeholder->ops->noo);
            printf("obs name  \t| mean   \t| var    \t| err    \t|\n");
            char name[]={"         "};
            for(int i_obs=0;i_obs<obs->nobs;i_obs++){
                memcpy(name,obs->obs_name[i_obs],6);
                printf("%s  \t| %.4e \t| %.4e \t| %.4e \t|\n",name,obs->mean[i_obs],obs->var[i_obs],obs->err[i_obs]);
            }
            int isweep = placeholder->isweep+1;
            int nsweep = placeholder->nsweep;
            double ratio = (double)isweep/nsweep*100;
            printf("|<-------------------%.2lf------------------->| %d/%d \n",ratio,isweep,nsweep);
            printf("===========================================================================\n");
    }
}

double ObservableSpecificEnergy(
                    SEPlaceHolder* placeholder, 
                    void* args)
{
    int Nb = placeholder->lconf->Nb;
    int nsite = placeholder->lconf->nsite;
    int noo = placeholder->ops->noo;
    double beta = placeholder->beta;
    double senergy = -noo/beta/nsite + (double)Nb/nsite*0.25;

    return senergy;
}

double ObservableMagnetization(
                    SEPlaceHolder* placeholder,
                    void* args)
{
    int nsite = placeholder->lconf->nsite;
    double mz=0;
    for(int i=0;i<nsite;++i){
        mz+=(double)placeholder->lconf->sigma0->data[i];
    }
    mz = mz/nsite*0.5;

    return mz;
}

double ObservableSusceptibility(
                    SEPlaceHolder* placeholder,
                    void* args)
{
    int nsite = placeholder->lconf->nsite;
    double mz=0,beta=placeholder->beta;
    for(int i=0;i<nsite;++i){
        mz+=(double)placeholder->lconf->sigma0->data[i];
    }
    mz = mz*mz*beta/nsite*0.25;

    return mz;
}

double ObservableStiffnessX(
                    SEPlaceHolder* placeholder,
                    void* args)
{
    int i,type,bond,p,left,right,length=placeholder->ops->length;
    int ndiff=placeholder->ops->ndiff;
    int nsite=placeholder->lconf->nsite;
    double winding=0;
    double beta = placeholder->beta;
    for(i=0;i<nsite;++i){
        placeholder->lconf->sigmap->data[i] = placeholder->lconf->sigma0->data[i];
    }

    for(p=0;p<length;++p){
        type = placeholder->ops->sequence->data[p]%ndiff;
        if(type==1){
            bond = placeholder->ops->sequence->data[p]/ndiff;
            LatticeConfApplyMapping(placeholder->lconf,bond);
            right = placeholder->lconf->right;
            left  = placeholder->lconf->left;
            if(bond<nsite){
                winding += placeholder->lconf->sigmap->data[left];
            }
            placeholder->lconf->sigmap->data[left] *=-1;
            placeholder->lconf->sigmap->data[right] *=-1;
        }
    }

    return winding*winding/placeholder->lconf->shape[0]/placeholder->lconf->shape[0]/beta;
}
