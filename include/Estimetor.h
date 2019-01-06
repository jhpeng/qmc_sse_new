#include <time.h>
#include <DataStruct.h>

#ifndef ESTIMETOR_H
#define ESTIMETOR_H

Observable* CreateObservable(int nobs, int nave);

void DestroyObservable(Observable* Obs);

void MeanAverage(Observable* Obs);

void ObservableSetMeasurement(
                    Observable* obs, 
                    measurement* measure, 
                    char* obs_name,
                    void* args);

void ObservableDoMeasurement(
                    Observable* obs,
                    SEPlaceHolder* placeholder);

void ObservableShow(
                    Observable* obs,
                    SEPlaceHolder* placeholder,
                    int mode);

double ObservableSpecificEnergy(
                    SEPlaceHolder* placeholder, 
                    void* args);


double ObservableMagnetization(
                    SEPlaceHolder* placeholder,
                    void* args);

double ObservableSusceptibility(
                    SEPlaceHolder* placeholder,
                    void* args);

double ObservableStiffnessX(
                    SEPlaceHolder* placeholder,
                    void* args);



#endif
