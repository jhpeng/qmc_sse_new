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
                    char* obs_name)
{
    obs->measure[obs->set_obs] = measure;
    strcpy(obs->obs_name[obs->set_obs],obs_name);
    obs->set_obs++;
}
