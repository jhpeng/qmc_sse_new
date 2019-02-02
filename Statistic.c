#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "DataStruct.h"
#include "SEAlgorithm.h"
#include "Estimetor.h"

double StatisticAutocorrelation(Observable* obs, int i_obs, int length)
{
    if(length>(obs->nave/2)){
        printf("StatisticAutocorrelation : Do not have enought data to measure the correlation!\n");
        exit(-1);
    }
    int nobs=obs->nobs;
    int nave=obs->nave;

    MeanAverage(obs);
    
    double autoc=0;
    double data_t1;
    double data_t2;
    for(int i=0;i<(nave-length);++i){
        data_t1 = obs->data->data[i*nobs+i_obs];
        data_t2 = obs->data->data[(i+length)*nobs+i_obs];

        autoc += (data_t1-obs->mean[i_obs])*(data_t2-obs->mean[i_obs]);
    }

    autoc = autoc/nave/obs->var[i_obs];

    return autoc;
}
