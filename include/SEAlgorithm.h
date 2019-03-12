#include <stdlib.h>
#include <gsl/gsl_rng.h>

#include "DataStruct.h"

#ifndef SEALGORITHM_H
#define SEALGORITHM_H
void LatticeConfSetMapping(LatticeConf* lconf, bond2sigma *mapping);

void LatticeConfApplyMapping(LatticeConf* lconf, int bond);

void LatticeConfSynchronizeSigma(LatticeConf* lconf);

void OperatorSequenceSort(OperatorSequence* ops);

void SEOps2Lvc(
                    OperatorLoop* opl,
                    LatticeConf* lconf,
                    OperatorSequence* ops);

void SETraverseLoop(
                    OperatorLoop* opl,
                    LatticeConf* lconf,
                    OperatorSequence* ops);

void SELoopUpdate(
                    LatticeConf* lconf,
                    OperatorSequence* ops,
                    OperatorLoop* opl);

void CreateMappingList(bond2sigma* mapping, const int* shape, int Nb);

void DestroyMappingList();

void mapping_list(int* left, int* right, const int* shape, int bond);

void mapping_1d(int* left, int* right, const int* shape, int bond);

void mapping_2d(int *first, int *second, const int *size, int bond);

void SEPlaceHolderSetLattice(
                    SEPlaceHolder* placeholder, 
                    bond2sigma *mapping, 
                    const int* shape, 
                    int dims, 
                    int model);

void SEPlaceHolderSetLength(SEPlaceHolder* placeholder, int length, int ndiff);

void SEPlaceHolderSetRandomSeed(SEPlaceHolder* placeholder, int seed);

void SEPlaceHolderSetNsweep(SEPlaceHolder* placeholder, int nsweep, int cutoff);

void SEPlaceHolderSetBeta(SEPlaceHolder* placeholder, double beta);

void SEPlaceHolderSetError(SEPlaceHolder* placeholder, double max_err);

int SEPlaceHolderCheckSetting(SEPlaceHolder* placeholder);

void SEPlaceHolderLengthMonitor(SEPlaceHolder* placeholder, double buffer);

void SEPlaceHolderSetDisorder2D(
                    SEPlaceHolder* placeholder, 
                    double J);

void SEPlaceHolderSetHerringbond2D(
                    SEPlaceHolder* placeholder, 
                    double J);

void SEPlaceHolderSetHerringbondRandom2D(
                    SEPlaceHolder* placeholder, 
                    double Jc, 
                    double dJ, 
                    double p);

void SEPlaceHolderSetPlaquetteRandom2D(
                    SEPlaceHolder* placeholder, 
                    double Jc, 
                    double dJ, 
                    double p);

void SEPlaceHolderSetConfigurationalDisorder2D(
                    SEPlaceHolder* placeholder, 
                    double J);

void SEPlaceHolderBetaDoubling(SEPlaceHolder* placeholder);

#endif
