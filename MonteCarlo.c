#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>

#include "DataStruct.h"
#include "SEAlgorithm.h"
#include "Estimetor.h"
#include "Statistic.h"

void MCInitializeLatticeConf(SEPlaceHolder* placeholder)
{
    int nsite=placeholder->lconf->nsite;
    double dis;
    for(int i=0;i<nsite;++i){
        dis  = gsl_rng_uniform_pos(placeholder->rng);
        if(dis<0.5) placeholder->lconf->sigma0->data[i]=1;
        else placeholder->lconf->sigma0->data[i]=-1;
    }
}

void MCDiagonalOperatorUpdate(SEPlaceHolder* placeholder)
{
    int length=placeholder->length;
    int ndiff=placeholder->ops->ndiff;
    //OperatorSequenceSort(placeholder->ops);
    int noo=placeholder->ops->noo;
    int Nb=placeholder->lconf->Nb;
    int p,bond,left,right;
    double beta=placeholder->beta;
    double dis;

    LatticeConfSynchronizeSigma(placeholder->lconf);

    for(p=0;p<length;++p){
        if(placeholder->ops->sequence->data[p]==-1){
            bond = (int)(gsl_rng_uniform_pos(placeholder->rng)*Nb);
            dis  = gsl_rng_uniform_pos(placeholder->rng);
            LatticeConfApplyMapping(placeholder->lconf,bond);
            left=placeholder->lconf->left;
            right=placeholder->lconf->right;
            if(dis*2*(length-noo)<beta*placeholder->lconf->J->data[bond]*Nb && placeholder->lconf->sigmap->data[left]!=placeholder->lconf->sigmap->data[right]){
                placeholder->ops->sequence->data[p]=bond*ndiff;
                noo++;
            }
        }
        else if(placeholder->ops->sequence->data[p]%ndiff==0){
            bond = placeholder->ops->sequence->data[p]/ndiff;
            dis  = gsl_rng_uniform_pos(placeholder->rng);
            if(dis*beta*placeholder->lconf->J->data[bond]*Nb<2*(length-noo+1)){
                placeholder->ops->sequence->data[p]=-1;
                noo--;
            }
        }
        else if(placeholder->ops->sequence->data[p]%ndiff==1){
            bond = placeholder->ops->sequence->data[p]/ndiff;
            LatticeConfApplyMapping(placeholder->lconf,bond);
            left=placeholder->lconf->left;
            right=placeholder->lconf->right;
            placeholder->lconf->sigmap->data[left] *=-1;
            placeholder->lconf->sigmap->data[right]*=-1;
        }
    }

    OperatorSequenceSort(placeholder->ops);
}

void MCOffDiagOperatorUpdate(SEPlaceHolder* placeholder)
{
    double dis;
    SEOps2Lvc(placeholder->opl,placeholder->lconf,placeholder->ops);
    placeholder->opl->check=-1;
    while(placeholder->opl->check){
        SETraverseLoop(placeholder->opl,placeholder->lconf,placeholder->ops);
        dis = gsl_rng_uniform_pos(placeholder->rng);
        if(dis<0.5){
            SELoopUpdate(placeholder->lconf,placeholder->ops,placeholder->opl);
        }
    }
}

void MCFlipUpdate(SEPlaceHolder* placeholder)
{
    double dis;
    int nsite=placeholder->lconf->nsite;

    for(int i=0;i<nsite;++i){
        if(placeholder->lconf->last->data[i]==-1){
            dis = gsl_rng_uniform_pos(placeholder->rng);
            if(dis<0.5){
                placeholder->lconf->sigma0->data[i]*=-1;
            }
        }
    }
}

int MCCheckInnerProduct(SEPlaceHolder* placeholder)
{
    int left,right,i,p,bond,type;
    int noo=placeholder->ops->noo;
    int ndiff=placeholder->ops->ndiff;
    int nsite=placeholder->lconf->nsite;

    LatticeConfSynchronizeSigma(placeholder->lconf);

    for(i=0;i<noo;++i){
        p=placeholder->ops->sort->data[i];
        bond=placeholder->ops->sequence->data[p]/ndiff;
        type=placeholder->ops->sequence->data[p]%ndiff;
        if(type==1){
            LatticeConfApplyMapping(placeholder->lconf,bond);
            left =placeholder->lconf->left;
            right=placeholder->lconf->right;
            placeholder->lconf->sigmap->data[left] *=-1;
            placeholder->lconf->sigmap->data[right]*=-1;
        }
    }

    int check=0;
    for(i=0;i<nsite;++i) {
        if(placeholder->lconf->sigmap->data[i]!=placeholder->lconf->sigma0->data[i]){
            check=1;
            return check;
        }
    }

    return check;
}

void MCIsotropy2D(double beta, int* shape, int nsweep, int cutoff, int seed)
{
    int ndiff=2,length=50;
    int dims=2;
    double max_err=1.e-4;
    double buffer=1.3;
    char prefix[128];

    sprintf(prefix,"data/isotropy_shape_%d_%d_beta_%.1f",shape[0],shape[1],beta);

    int Nb=shape[0]*shape[1]*dims;
    CreateMappingList(mapping_2d,shape,Nb);

    SEPlaceHolder* placeholder = CreateSEPlaceHolder();
    SEPlaceHolderSetLattice(placeholder,mapping_list,shape,dims,0);
    SEPlaceHolderSetLength(placeholder,length,ndiff);
    SEPlaceHolderSetRandomSeed(placeholder, seed);
    SEPlaceHolderSetNsweep(placeholder, nsweep, cutoff);
    SEPlaceHolderSetBeta(placeholder, beta);
    SEPlaceHolderSetError(placeholder, max_err);
#if 0
    SEPlaceHolderCheckSetting(placeholder);
#endif
    MCInitializeLatticeConf(placeholder);

    int j=0;
    for(j=0;j<cutoff;j++){
        MCDiagonalOperatorUpdate(placeholder);
        MCOffDiagOperatorUpdate(placeholder);
        MCFlipUpdate(placeholder);
#if 0
#define CHECK_INNER_PRODUCT
#endif
#ifdef CHECK_INNER_PRODUCT
        int check = MCCheckInnerProduct(placeholder);
        if(check) {
            printf("faile passing the check inner product!\n");
            exit(-1);
        }
#endif
        SEPlaceHolderLengthMonitor(placeholder, buffer);
    }

    int nobs=7;
    int nave=nsweep;
    Observable *obs = CreateObservable(nobs,nave);
    ObservableSetMeasurement(obs,ObservableSpecificEnergy,"energy",NULL);
    ObservableSetMeasurement(obs,ObservableMagnetization,"magn_z",NULL);
    ObservableSetMeasurement(obs,ObservableSusceptibility,"susc_z",NULL);
#if 1
#define FAST_OBS
#endif
#ifndef FAST_OBS
    ObservableSetMeasurement(obs,ObservableStiffnessX,"stif_x",NULL);
    ObservableSetMeasurement(obs,ObservableAntiferroOrder1,"mz_1",NULL);
    ObservableSetMeasurement(obs,ObservableAntiferroOrder2,"mz_2",NULL);
    ObservableSetMeasurement(obs,ObservableAntiferroOrder4,"mz_4",NULL);
#endif
#ifdef FAST_OBS
    ObservableSetMeasurement(obs,ObservableFastStiffnessX,"stif_x",NULL);
    ObservableSetMeasurement(obs,ObservableFastAntiferroOrder1,"mz_1",NULL);
    ObservableSetMeasurement(obs,ObservableFastAntiferroOrder2,"mz_2",NULL);
    ObservableSetMeasurement(obs,ObservableFastAntiferroOrder4,"mz_4",NULL);
#endif

    j=0;
    placeholder->isweep=0;
    while(j<nsweep){
        MCDiagonalOperatorUpdate(placeholder);
        MCOffDiagOperatorUpdate(placeholder);
        MCFlipUpdate(placeholder);

#ifdef FAST_OBS
        ObservableFastPreCal(placeholder);
#endif
        ObservableDoMeasurement(obs,placeholder);
        if((j+1)%10000==0){
            //ObservableShow(obs,placeholder,NULL,0);
            //ObservableShow(obs,placeholder,prefix,1);
            ObservableShow(obs,placeholder,prefix,2);
            ObservableShow(obs,placeholder,prefix,3);
        }
#ifdef CHECK_INNER_PRODUCT
        int check = MCCheckInnerProduct(placeholder);
        if(check) {
            printf("faile passing the check inner product!\n");
            exit(-1);
        }
#endif

        SEPlaceHolderLengthMonitor(placeholder, buffer);
        placeholder->isweep++;
        j++;
    }

    DestroySEPlaceHolder(placeholder);
    DestroyMappingList();
}

void MCDisorder2D(double J, double beta, int* shape, int nsweep, int cutoff, int seed)
{
    int ndiff=2,length=50;
    int dims=2;
    double max_err=1.e-4;
    double buffer=1.3;
    char prefix[128];

    sprintf(prefix,"data/disorder_shape_%d_%d_J_%.4f_beta_%.1f",shape[0],shape[1],J,beta);

    int Nb=shape[0]*shape[1]*dims;
    CreateMappingList(mapping_2d,shape,Nb);

    SEPlaceHolder* placeholder = CreateSEPlaceHolder();
    SEPlaceHolderSetLattice(placeholder,mapping_list,shape,dims,0);
    SEPlaceHolderSetLength(placeholder,length,ndiff);
    SEPlaceHolderSetRandomSeed(placeholder, seed);
    SEPlaceHolderSetNsweep(placeholder, nsweep, cutoff);
    SEPlaceHolderSetBeta(placeholder, beta);
    SEPlaceHolderSetError(placeholder, max_err);
    SEPlaceHolderSetDisorder2D(placeholder,J);
    SEPlaceHolderCheckSetting(placeholder);

    MCInitializeLatticeConf(placeholder);

    int j=0;
    for(j=0;j<cutoff;j++){
        MCDiagonalOperatorUpdate(placeholder);
        MCOffDiagOperatorUpdate(placeholder);
        MCFlipUpdate(placeholder);
        SEPlaceHolderLengthMonitor(placeholder, buffer);
    }

    int nobs=8;
    int nave=nsweep;
    Observable *obs = CreateObservable(nobs,nave);
    ObservableSetMeasurement(obs,ObservableSpecificEnergy,"energy",NULL);
    ObservableSetMeasurement(obs,ObservableMagnetization,"magn_z",NULL);
    ObservableSetMeasurement(obs,ObservableSusceptibility,"susc_z",NULL);
    ObservableSetMeasurement(obs,ObservableFastStiffnessX,"stif_x",NULL);
    ObservableSetMeasurement(obs,ObservableFastStiffnessY,"stif_y",NULL);
    ObservableSetMeasurement(obs,ObservableFastAntiferroOrder1,"mz_1",NULL);
    ObservableSetMeasurement(obs,ObservableFastAntiferroOrder2,"mz_2",NULL);
    ObservableSetMeasurement(obs,ObservableFastAntiferroOrder4,"mz_4",NULL);

    j=0;
    placeholder->isweep=0;
    while(j<nsweep){
        MCDiagonalOperatorUpdate(placeholder);
        MCOffDiagOperatorUpdate(placeholder);
        MCFlipUpdate(placeholder);

        ObservableFastPreCal(placeholder);
        ObservableDoMeasurement(obs,placeholder);
        if((j+1)%10000==0){
            //ObservableShow(obs,placeholder,NULL,0);
            //ObservableShow(obs,placeholder,prefix,1);
            ObservableShow(obs,placeholder,prefix,2);
            ObservableShow(obs,placeholder,prefix,3);
        }

        SEPlaceHolderLengthMonitor(placeholder, buffer);
        placeholder->isweep++;
        j++;
    }

    DestroySEPlaceHolder(placeholder);
    DestroyMappingList();
}

void MCHerringbond2D(double J, double beta, int* shape, int nsweep, int nblock, int cutoff, int seed)
{
    int ndiff=2,length=50;
    int dims=2;
    double max_err=1.e-4;
    double buffer=1.3;
    char prefix[128];

    sprintf(prefix,"data/herringbond_shape_%d_%d_J_%.4f_beta_%.1f",shape[0],shape[1],J,beta);

    int Nb=shape[0]*shape[1]*dims;
    CreateMappingList(mapping_2d,shape,Nb);

    SEPlaceHolder* placeholder = CreateSEPlaceHolder();
    SEPlaceHolderSetLattice(placeholder,mapping_list,shape,dims,0);
    SEPlaceHolderSetLength(placeholder,length,ndiff);
    SEPlaceHolderSetRandomSeed(placeholder, seed);
    SEPlaceHolderSetNsweep(placeholder, nsweep, cutoff);
    SEPlaceHolderSetBeta(placeholder, beta);
    SEPlaceHolderSetError(placeholder, max_err);
    SEPlaceHolderSetHerringbond2D(placeholder,J);
    SEPlaceHolderCheckSetting(placeholder);

    MCInitializeLatticeConf(placeholder);

    int j=0;
    for(j=0;j<cutoff;j++){
        MCDiagonalOperatorUpdate(placeholder);
        MCOffDiagOperatorUpdate(placeholder);
        MCFlipUpdate(placeholder);
        SEPlaceHolderLengthMonitor(placeholder, buffer);
    }

    int nobs=8;
    int nave=nsweep;
    Observable *obs = CreateObservable(nobs,nave);
    ObservableSetMeasurement(obs,ObservableSpecificEnergy,"energy",NULL);
    ObservableSetMeasurement(obs,ObservableMagnetization,"magn_z",NULL);
    ObservableSetMeasurement(obs,ObservableSusceptibility,"susc_z",NULL);
    ObservableSetMeasurement(obs,ObservableFastStiffnessX,"stif_x",NULL);
    ObservableSetMeasurement(obs,ObservableFastStiffnessY,"stif_y",NULL);
    ObservableSetMeasurement(obs,ObservableFastAntiferroOrder1,"mz_1",NULL);
    ObservableSetMeasurement(obs,ObservableFastAntiferroOrder2,"mz_2",NULL);
    ObservableSetMeasurement(obs,ObservableFastAntiferroOrder4,"mz_4",NULL);

    for(int i_block=0;i_block<nblock;i_block++){
        placeholder->isweep=0;
        while(placeholder->isweep<nsweep){
            MCDiagonalOperatorUpdate(placeholder);
            MCOffDiagOperatorUpdate(placeholder);
            MCFlipUpdate(placeholder);

            ObservableFastPreCal(placeholder);
            ObservableDoMeasurement(obs,placeholder);

            SEPlaceHolderLengthMonitor(placeholder, buffer);
            placeholder->isweep++;
        }
        
        ObservableShow(obs,placeholder,prefix,4);
    }

    DestroySEPlaceHolder(placeholder);
    DestroyMappingList();
}

void MCHerringbond2DImproveSpeed(double J, double beta, int* shape, int nsweep, int nblock, int cutoff, int seed)
{
    int ndiff=2,length=50;
    int dims=2;
    double max_err=1.e-4;
    double buffer=1.3;
    char prefix[128];

    sprintf(prefix,"data/herringbond_improve_speed_shape_%d_%d_J_%.4f_beta_%.1f",shape[0],shape[1],J,beta);

    int Nb=shape[0]*shape[1]*dims;
    CreateMappingList(mapping_2d,shape,Nb);

    SEPlaceHolder* placeholder = CreateSEPlaceHolder();
    SEPlaceHolderSetLattice(placeholder,mapping_list,shape,dims,0);
    SEPlaceHolderSetLength(placeholder,length,ndiff);
    SEPlaceHolderSetRandomSeed(placeholder, seed);
    SEPlaceHolderSetNsweep(placeholder, nsweep, cutoff);
    SEPlaceHolderSetBeta(placeholder, beta);
    SEPlaceHolderSetError(placeholder, max_err);
    SEPlaceHolderSetHerringbond2D(placeholder,J);
    SEPlaceHolderCheckSetting(placeholder);

    MCInitializeLatticeConf(placeholder);

    int j=0;
    for(j=0;j<cutoff;j++){
        MCDiagonalOperatorUpdate(placeholder);
        MCOffDiagOperatorUpdate(placeholder);
        MCFlipUpdate(placeholder);
        SEPlaceHolderLengthMonitor(placeholder, buffer);
    }

    int nobs=6;
    int nave=nsweep;
    Observable *obs = CreateObservable(nobs,nave);
    ObservableSetMeasurement(obs,ObservableSpecificEnergy,"energy",NULL);
    ObservableSetMeasurement(obs,ObservableMagnetization,"magn_z",NULL);
    ObservableSetMeasurement(obs,ObservableSusceptibility,"susc_z",NULL);
    ObservableSetMeasurement(obs,ObservableFastAntiferroOrder1,"mz_1",NULL);
    ObservableSetMeasurement(obs,ObservableFastAntiferroOrder2,"mz_2",NULL);
    ObservableSetMeasurement(obs,ObservableFastAntiferroOrder4,"mz_4",NULL);

    for(int i_block=0;i_block<nblock;i_block++){
        placeholder->isweep=0;
        while(placeholder->isweep<nsweep){
            MCDiagonalOperatorUpdate(placeholder);
            MCOffDiagOperatorUpdate(placeholder);
            MCFlipUpdate(placeholder);

            ObservableImproveSpeedPreCal(placeholder);
            ObservableDoMeasurement(obs,placeholder);

            SEPlaceHolderLengthMonitor(placeholder, buffer);
            placeholder->isweep++;
        }
        
        ObservableShow(obs,placeholder,prefix,4);
    }

    DestroySEPlaceHolder(placeholder);
    DestroyMappingList();
}

void MCBetaIncreasePlaquetteDisorder2D(double J, double dJ, double p, double beta_i, double beta_f, double interval, int* shape, int nsweep, int cutoff, int seed)
{
    int ndiff=2,length=1000;
    int dims=2;
    double max_err=1.e-4;
    double buffer=1.5;
    char prefix[128];

    sprintf(prefix,"data/plaquette_disorder_shape_%d_%d_J_%.4f_bi_%.1f_p_%.4f_seed_%d",shape[0],shape[1],J,beta_i,p,seed);

    int Nb=shape[0]*shape[1]*dims;
    CreateMappingList(mapping_2d,shape,Nb);

    SEPlaceHolder* placeholder = CreateSEPlaceHolder();
    SEPlaceHolderSetLattice(placeholder,mapping_list,shape,dims,0);
    SEPlaceHolderSetLength(placeholder,length,ndiff);
    SEPlaceHolderSetRandomSeed(placeholder, seed);
    SEPlaceHolderSetNsweep(placeholder, nsweep, cutoff);
    SEPlaceHolderSetError(placeholder, max_err);
    SEPlaceHolderSetPlaquetteRandom2D(placeholder,J,dJ,p);

    MCInitializeLatticeConf(placeholder);

    double beta;
    for(beta=beta_i;beta<beta_f;beta+=interval){
        SEPlaceHolderSetBeta(placeholder, beta);
        SEPlaceHolderCheckSetting(placeholder);
        int j=0;
        for(j=0;j<cutoff;j++){
            MCDiagonalOperatorUpdate(placeholder);
            MCOffDiagOperatorUpdate(placeholder);
            MCFlipUpdate(placeholder);
            SEPlaceHolderLengthMonitor(placeholder, buffer);
        }

        int nobs=8;
        int nave=nsweep;
        Observable *obs = CreateObservable(nobs,nave);
        ObservableSetMeasurement(obs,ObservableSpecificEnergy,"energy",NULL);
        ObservableSetMeasurement(obs,ObservableMagnetization,"magn_z",NULL);
        ObservableSetMeasurement(obs,ObservableSusceptibility,"susc_z",NULL);
        ObservableSetMeasurement(obs,ObservableFastStiffnessX,"stif_x",NULL);
        ObservableSetMeasurement(obs,ObservableFastStiffnessY,"stif_y",NULL);
        ObservableSetMeasurement(obs,ObservableFastAntiferroOrder1,"mz_1",NULL);
        ObservableSetMeasurement(obs,ObservableFastAntiferroOrder2,"mz_2",NULL);
        ObservableSetMeasurement(obs,ObservableFastAntiferroOrder4,"mz_4",NULL);

        j=0;
        placeholder->isweep=0;
        while(j<nsweep){
            MCDiagonalOperatorUpdate(placeholder);
            MCOffDiagOperatorUpdate(placeholder);
            MCFlipUpdate(placeholder);

            ObservableFastPreCal(placeholder);
            ObservableDoMeasurement(obs,placeholder);

            //SEPlaceHolderLengthMonitor(placeholder, buffer);
            placeholder->isweep++;
            j++;
        }
        ObservableShow(obs,placeholder,prefix,4);
        DestroyObservable(obs);
    }

    DestroySEPlaceHolder(placeholder);
    DestroyMappingList();
}

void MCBetaIncreasePlaquetteDisorderImproveSpeed2D(double J, double dJ, double p, double beta_i, double beta_f, double interval, int* shape, int nsweep, int cutoff, int seed)
{
    int ndiff=2,length=1000;
    int dims=2;
    double max_err=1.e-4;
    double buffer=1.5;
    char prefix[128];

    sprintf(prefix,"data/plaquette_disorder_improve_speed_shape_%d_%d_J_%.4f_bi_%.1f_p_%.4f_seed_%d",shape[0],shape[1],J,beta_i,p,seed);

    int Nb=shape[0]*shape[1]*dims;
    CreateMappingList(mapping_2d,shape,Nb);

    SEPlaceHolder* placeholder = CreateSEPlaceHolder();
    SEPlaceHolderSetLattice(placeholder,mapping_list,shape,dims,0);
    SEPlaceHolderSetLength(placeholder,length,ndiff);
    SEPlaceHolderSetRandomSeed(placeholder, seed);
    SEPlaceHolderSetNsweep(placeholder, nsweep, cutoff);
    SEPlaceHolderSetError(placeholder, max_err);
    SEPlaceHolderSetPlaquetteRandom2D(placeholder,J,dJ,p);

    MCInitializeLatticeConf(placeholder);

    double beta;
    for(beta=beta_i;beta<beta_f;beta+=interval){
        SEPlaceHolderSetBeta(placeholder, beta);
        SEPlaceHolderCheckSetting(placeholder);
        int j=0;
        for(j=0;j<cutoff;j++){
            MCDiagonalOperatorUpdate(placeholder);
            MCOffDiagOperatorUpdate(placeholder);
            MCFlipUpdate(placeholder);
            SEPlaceHolderLengthMonitor(placeholder, buffer);
        }

        int nobs=6;
        int nave=nsweep;
        Observable *obs = CreateObservable(nobs,nave);
        ObservableSetMeasurement(obs,ObservableSpecificEnergy,"energy",NULL);
        ObservableSetMeasurement(obs,ObservableMagnetization,"magn_z",NULL);
        ObservableSetMeasurement(obs,ObservableSusceptibility,"susc_z",NULL);
        ObservableSetMeasurement(obs,ObservableFastAntiferroOrder1,"mz_1",NULL);
        ObservableSetMeasurement(obs,ObservableFastAntiferroOrder2,"mz_2",NULL);
        ObservableSetMeasurement(obs,ObservableFastAntiferroOrder4,"mz_4",NULL);

        j=0;
        placeholder->isweep=0;
        while(j<nsweep){
            MCDiagonalOperatorUpdate(placeholder);
            MCOffDiagOperatorUpdate(placeholder);
            MCFlipUpdate(placeholder);

            ObservableImproveSpeedPreCal(placeholder);
            ObservableDoMeasurement(obs,placeholder);

            //SEPlaceHolderLengthMonitor(placeholder, buffer);
            placeholder->isweep++;
            j++;
        }
        ObservableShow(obs,placeholder,prefix,4);
        DestroyObservable(obs);
    }

    DestroySEPlaceHolder(placeholder);
    DestroyMappingList();
}

void MCZeroTempPlaquetteDisorder2D(double J, double dJ, double p, double beta, int* shape, int nsweep, int cutoff, int ntime, int seed)
{
    int ndiff=2,length=1000;
    int dims=2;
    double max_err=1.e-4;
    double buffer=1.3;
    char prefix[128];

    sprintf(prefix,"data/zero_temp_plaquette_disorder_shape_%d_%d_J_%.4f_beta_%.2f_p_%.4f_seed_%d",shape[0],shape[1],J,beta,p,seed);

    int Nb=shape[0]*shape[1]*dims;
    CreateMappingList(mapping_2d,shape,Nb);

    SEPlaceHolder* placeholder = CreateSEPlaceHolder();
    SEPlaceHolderSetLattice(placeholder,mapping_list,shape,dims,0);
    SEPlaceHolderSetLength(placeholder,length,ndiff);
    SEPlaceHolderSetRandomSeed(placeholder, seed);
    SEPlaceHolderSetNsweep(placeholder, nsweep, cutoff);
    SEPlaceHolderSetError(placeholder, max_err);
    SEPlaceHolderSetPlaquetteRandom2D(placeholder,J,dJ,p);

    MCInitializeLatticeConf(placeholder);

    SEPlaceHolderSetBeta(placeholder, beta);
    SEPlaceHolderCheckSetting(placeholder);

    int nobs=8;
    int nave=nsweep;
    Observable *obs = CreateObservable(nobs,nave);
    ObservableSetMeasurement(obs,ObservableSpecificEnergy,"energy",NULL);
    ObservableSetMeasurement(obs,ObservableMagnetization,"magn_z",NULL);
    ObservableSetMeasurement(obs,ObservableSusceptibility,"susc_z",NULL);
    ObservableSetMeasurement(obs,ObservableFastStiffnessX,"stif_x",NULL);
    ObservableSetMeasurement(obs,ObservableFastStiffnessY,"stif_y",NULL);
    ObservableSetMeasurement(obs,ObservableFastAntiferroOrder1,"mz_1",NULL);
    ObservableSetMeasurement(obs,ObservableFastAntiferroOrder2,"mz_2",NULL);
    ObservableSetMeasurement(obs,ObservableFastAntiferroOrder4,"mz_4",NULL);

    for(int k=0;k<ntime*2;k++){
        int j=0;
        for(j=0;j<cutoff;j++){
            MCDiagonalOperatorUpdate(placeholder);
            MCOffDiagOperatorUpdate(placeholder);
            MCFlipUpdate(placeholder);
            SEPlaceHolderLengthMonitor(placeholder, buffer);
        }

        j=0;
        placeholder->isweep=0;
        while(j<nsweep){
            MCDiagonalOperatorUpdate(placeholder);
            MCOffDiagOperatorUpdate(placeholder);
            MCFlipUpdate(placeholder);

            ObservableFastPreCal(placeholder);
            ObservableDoMeasurement(obs,placeholder);

            SEPlaceHolderLengthMonitor(placeholder, buffer);
            placeholder->isweep++;
            j++;
        }
        ObservableShow(obs,placeholder,prefix,4);
        if(k%2==1) SEPlaceHolderBetaDoubling(placeholder);
    }

    if(0){
        int nc=20;
        for(int ic=1;ic<(nc+1);ic++){
            double ac = StatisticAutocorrelation(obs,6,ic);
            printf("%e ",ac);
        }
        printf("\n");
    }

    DestroyObservable(obs);

    DestroySEPlaceHolder(placeholder);
    DestroyMappingList();
}

void MCZeroTempHerringbondDisorder2D(double J, double dJ, double p, double beta, int* shape, int nsweep, int cutoff, int ntime, int seed)
{
    int ndiff=2,length=1000;
    int dims=2;
    double max_err=1.e-4;
    double buffer=1.3;
    char prefix[128];

    sprintf(prefix,"data/zero_temp_herringbond_disorder_shape_%d_%d_J_%.4f_beta_%.2f_p_%.4f_seed_%d",shape[0],shape[1],J,beta,p,seed);

    int Nb=shape[0]*shape[1]*dims;
    CreateMappingList(mapping_2d,shape,Nb);

    SEPlaceHolder* placeholder = CreateSEPlaceHolder();
    SEPlaceHolderSetLattice(placeholder,mapping_list,shape,dims,0);
    SEPlaceHolderSetLength(placeholder,length,ndiff);
    SEPlaceHolderSetRandomSeed(placeholder, seed);
    SEPlaceHolderSetNsweep(placeholder, nsweep, cutoff);
    SEPlaceHolderSetError(placeholder, max_err);
    SEPlaceHolderSetHerringbondRandom2D(placeholder,J,dJ,p);

    MCInitializeLatticeConf(placeholder);

    SEPlaceHolderSetBeta(placeholder, beta);
    SEPlaceHolderCheckSetting(placeholder);

    int nobs=8;
    int nave=nsweep;
    Observable *obs = CreateObservable(nobs,nave);
    ObservableSetMeasurement(obs,ObservableSpecificEnergy,"energy",NULL);
    ObservableSetMeasurement(obs,ObservableMagnetization,"magn_z",NULL);
    ObservableSetMeasurement(obs,ObservableSusceptibility,"susc_z",NULL);
    ObservableSetMeasurement(obs,ObservableFastStiffnessX,"stif_x",NULL);
    ObservableSetMeasurement(obs,ObservableFastStiffnessY,"stif_y",NULL);
    ObservableSetMeasurement(obs,ObservableFastAntiferroOrder1,"mz_1",NULL);
    ObservableSetMeasurement(obs,ObservableFastAntiferroOrder2,"mz_2",NULL);
    ObservableSetMeasurement(obs,ObservableFastAntiferroOrder4,"mz_4",NULL);

    for(int k=0;k<ntime*2;k++){
        int j=0;
        for(j=0;j<cutoff;j++){
            MCDiagonalOperatorUpdate(placeholder);
            MCOffDiagOperatorUpdate(placeholder);
            MCFlipUpdate(placeholder);
            SEPlaceHolderLengthMonitor(placeholder, buffer);
        }

        j=0;
        placeholder->isweep=0;
        while(j<nsweep){
            MCDiagonalOperatorUpdate(placeholder);
            MCOffDiagOperatorUpdate(placeholder);
            MCFlipUpdate(placeholder);

            ObservableFastPreCal(placeholder);
            ObservableDoMeasurement(obs,placeholder);

            SEPlaceHolderLengthMonitor(placeholder, buffer);
            placeholder->isweep++;
            j++;
        }
        ObservableShow(obs,placeholder,prefix,4);
        if(k%2==1) SEPlaceHolderBetaDoubling(placeholder);
    }

    if(0){
        int nc=20;
        for(int ic=1;ic<(nc+1);ic++){
            double ac = StatisticAutocorrelation(obs,6,ic);
            printf("%e ",ac);
        }
        printf("\n");
    }

    DestroyObservable(obs);

    DestroySEPlaceHolder(placeholder);
    DestroyMappingList();
}

void MCBetaIncreaseHerringboneDisorderImproveSpeed2D(double J, double dJ, double p, double beta_i, double beta_f, double interval, int* shape, int nsweep, int cutoff, int seed)
{
    int ndiff=2,length=1000;
    int dims=2;
    double max_err=1.e-4;
    double buffer=1.5;
    char prefix[128];

    sprintf(prefix,"data/herringbond_disorder_improve_speed_shape_%d_%d_J_%.4f_bi_%.1f_p_%.4f_seed_%d",shape[0],shape[1],J,beta_i,p,seed);

    int Nb=shape[0]*shape[1]*dims;
    CreateMappingList(mapping_2d,shape,Nb);

    SEPlaceHolder* placeholder = CreateSEPlaceHolder();
    SEPlaceHolderSetLattice(placeholder,mapping_list,shape,dims,0);
    SEPlaceHolderSetLength(placeholder,length,ndiff);
    SEPlaceHolderSetRandomSeed(placeholder, seed);
    SEPlaceHolderSetNsweep(placeholder, nsweep, cutoff);
    SEPlaceHolderSetError(placeholder, max_err);
    SEPlaceHolderSetHerringbondRandom2D(placeholder,J,dJ,p);

    MCInitializeLatticeConf(placeholder);

    double beta;
    for(beta=beta_i;beta<beta_f;beta+=interval){
        SEPlaceHolderSetBeta(placeholder, beta);
        SEPlaceHolderCheckSetting(placeholder);
        int j=0;
        for(j=0;j<cutoff;j++){
            MCDiagonalOperatorUpdate(placeholder);
            MCOffDiagOperatorUpdate(placeholder);
            MCFlipUpdate(placeholder);
            SEPlaceHolderLengthMonitor(placeholder, buffer);
        }

        int nobs=6;
        int nave=nsweep;
        Observable *obs = CreateObservable(nobs,nave);
        ObservableSetMeasurement(obs,ObservableSpecificEnergy,"energy",NULL);
        ObservableSetMeasurement(obs,ObservableMagnetization,"magn_z",NULL);
        ObservableSetMeasurement(obs,ObservableSusceptibility,"susc_z",NULL);
        ObservableSetMeasurement(obs,ObservableFastAntiferroOrder1,"mz_1",NULL);
        ObservableSetMeasurement(obs,ObservableFastAntiferroOrder2,"mz_2",NULL);
        ObservableSetMeasurement(obs,ObservableFastAntiferroOrder4,"mz_4",NULL);

        j=0;
        placeholder->isweep=0;
        while(j<nsweep){
            MCDiagonalOperatorUpdate(placeholder);
            MCOffDiagOperatorUpdate(placeholder);
            MCFlipUpdate(placeholder);

            ObservableImproveSpeedPreCal(placeholder);
            ObservableDoMeasurement(obs,placeholder);

            //SEPlaceHolderLengthMonitor(placeholder, buffer);
            placeholder->isweep++;
            j++;
        }
        ObservableShow(obs,placeholder,prefix,4);
        DestroyObservable(obs);
    }

    DestroySEPlaceHolder(placeholder);
    DestroyMappingList();
}

void MCBetaIncreaseConfigurationalDisorder2D(double J, double beta_i, double beta_f, double interval, int* shape, int nsweep, int cutoff, int seed)
{
    int ndiff=2,length=1000;
    int dims=2;
    double max_err=1.e-4;
    double buffer=1.5;
    char prefix[128];

    sprintf(prefix,"data/configurational_disorder_shape_%d_%d_J_%.4f_bi_%.1f_seed_%d",shape[0],shape[1],J,beta_i,seed);

    int Nb=shape[0]*shape[1]*dims;
    CreateMappingList(mapping_2d,shape,Nb);

    SEPlaceHolder* placeholder = CreateSEPlaceHolder();
    SEPlaceHolderSetLattice(placeholder,mapping_list,shape,dims,0);
    SEPlaceHolderSetLength(placeholder,length,ndiff);
    SEPlaceHolderSetRandomSeed(placeholder, seed);
    SEPlaceHolderSetNsweep(placeholder, nsweep, cutoff);
    SEPlaceHolderSetError(placeholder, max_err);
    SEPlaceHolderSetConfigurationalDisorder2D(placeholder,J);

    MCInitializeLatticeConf(placeholder);

    double beta;
    for(beta=beta_i;beta<beta_f;beta+=interval){
        SEPlaceHolderSetBeta(placeholder, beta);
        SEPlaceHolderCheckSetting(placeholder);
        int j=0;
        for(j=0;j<cutoff;j++){
            MCDiagonalOperatorUpdate(placeholder);
            MCOffDiagOperatorUpdate(placeholder);
            MCFlipUpdate(placeholder);
            SEPlaceHolderLengthMonitor(placeholder, buffer);
        }

        int nobs=8;
        int nave=nsweep;
        Observable *obs = CreateObservable(nobs,nave);
        ObservableSetMeasurement(obs,ObservableSpecificEnergy,"energy",NULL);
        ObservableSetMeasurement(obs,ObservableMagnetization,"magn_z",NULL);
        ObservableSetMeasurement(obs,ObservableSusceptibility,"susc_z",NULL);
        ObservableSetMeasurement(obs,ObservableFastStiffnessX,"stif_x",NULL);
        ObservableSetMeasurement(obs,ObservableFastStiffnessY,"stif_y",NULL);
        ObservableSetMeasurement(obs,ObservableFastAntiferroOrder1,"mz_1",NULL);
        ObservableSetMeasurement(obs,ObservableFastAntiferroOrder2,"mz_2",NULL);
        ObservableSetMeasurement(obs,ObservableFastAntiferroOrder4,"mz_4",NULL);

        j=0;
        placeholder->isweep=0;
        while(j<nsweep){
            MCDiagonalOperatorUpdate(placeholder);
            MCOffDiagOperatorUpdate(placeholder);
            MCFlipUpdate(placeholder);

            ObservableFastPreCal(placeholder);
            ObservableDoMeasurement(obs,placeholder);

            //SEPlaceHolderLengthMonitor(placeholder, buffer);
            placeholder->isweep++;
            j++;
        }
        ObservableShow(obs,placeholder,prefix,4);
        DestroyObservable(obs);
    }

    DestroySEPlaceHolder(placeholder);
    DestroyMappingList();
}

void MCGeneralSchemeAndLattice(int* shape, int mode, int lattice, double J, double dJ, double p, double beta, double beta_i, double beta_f, double interv, int thermal, int nsweep, int nblock, int ntime, int seed){
    int ndiff=2,length=1000;
    int dims=2;
    double max_err=1.e-4;
    double buffer=1.5;
    char prefix[128];

    sprintf(prefix,"data/lattice_%d_scheme_%d_shape_%d_%d_J_%.4f_dJ_%.4f_p_%.4f_seed_%d",lattice,mode,shape[0],shape[1],J,dJ,p,seed);

    int Nb=shape[0]*shape[1]*dims;
    CreateMappingList(mapping_2d,shape,Nb);

    SEPlaceHolder* placeholder = CreateSEPlaceHolder();
    SEPlaceHolderSetLattice(placeholder,mapping_list,shape,dims,0);
    SEPlaceHolderSetLength(placeholder,length,ndiff);
    SEPlaceHolderSetRandomSeed(placeholder, seed);
    SEPlaceHolderSetNsweep(placeholder, nsweep, thermal);
    SEPlaceHolderSetError(placeholder, max_err);

    if(lattice==0) SEPlaceHolderSetHerringbondRandom2D(placeholder,J,dJ,p);
    else if(lattice==1) SEPlaceHolderSetPlaquetteRandom2D(placeholder,J,dJ,p);
    else if(lattice==2) SEPlaceHolderSetConfigurationalDisorder2D(placeholder,J);

    MCInitializeLatticeConf(placeholder);

/* ------------------------------------------- **
** -------------- Normal Scheme -------------- **
** ------------------------------------------- */
    if(mode==0){
        SEPlaceHolderSetBeta(placeholder, beta);
        SEPlaceHolderCheckSetting(placeholder);

        int nobs=8;
        int nave=nsweep;
        Observable *obs = CreateObservable(nobs,nave);
        ObservableSetMeasurement(obs,ObservableSpecificEnergy,"energy",NULL);
        ObservableSetMeasurement(obs,ObservableMagnetization,"magn_z",NULL);
        //ObservableSetMeasurement(obs,ObservableFastStaggeredX,"stag_x",NULL);
        //ObservableSetMeasurement(obs,ObservableFastStaggeredY,"stag_y",NULL);
        ObservableSetMeasurement(obs,ObservableSusceptibility,"susc_z",NULL);
        ObservableSetMeasurement(obs,ObservableFastStiffnessX,"stif_x",NULL);
        ObservableSetMeasurement(obs,ObservableFastStiffnessY,"stif_y",NULL);
        ObservableSetMeasurement(obs,ObservableFastAntiferroOrder1,"mz_1",NULL);
        ObservableSetMeasurement(obs,ObservableFastAntiferroOrder2,"mz_2",NULL);
        ObservableSetMeasurement(obs,ObservableFastAntiferroOrder4,"mz_4",NULL);
        
        for(int j=0;j<thermal;j++){
            MCDiagonalOperatorUpdate(placeholder);
            MCOffDiagOperatorUpdate(placeholder);
            MCFlipUpdate(placeholder);
            SEPlaceHolderLengthMonitor(placeholder, buffer);
        }
        for(int k=0;k<nblock;++k){
            placeholder->isweep=0;
            for(int j=0;j<nsweep;++j){
                MCDiagonalOperatorUpdate(placeholder);
                MCOffDiagOperatorUpdate(placeholder);
                MCFlipUpdate(placeholder);

                ObservableFastPreCal(placeholder);
                ObservableDoMeasurement(obs,placeholder);

                placeholder->isweep++;
            }
            ObservableShow(obs,placeholder,prefix,4);
        }
        DestroyObservable(obs);
    }
/* ------------------------------------------- **
** -------------- Beta Doubling -------------- **
** ------------------------------------------- */
    else if(mode==1){
        SEPlaceHolderSetBeta(placeholder, beta);
        SEPlaceHolderCheckSetting(placeholder);

        int nobs=8;
        int nave=nsweep;
        Observable *obs = CreateObservable(nobs,nave);
        ObservableSetMeasurement(obs,ObservableSpecificEnergy,"energy",NULL);
        ObservableSetMeasurement(obs,ObservableMagnetization,"magn_z",NULL);
        //ObservableSetMeasurement(obs,ObservableFastStaggeredX,"stag_x",NULL);
        //ObservableSetMeasurement(obs,ObservableFastStaggeredY,"stag_y",NULL);
        ObservableSetMeasurement(obs,ObservableSusceptibility,"susc_z",NULL);
        ObservableSetMeasurement(obs,ObservableFastStiffnessX,"stif_x",NULL);
        ObservableSetMeasurement(obs,ObservableFastStiffnessY,"stif_y",NULL);
        ObservableSetMeasurement(obs,ObservableFastAntiferroOrder1,"mz_1",NULL);
        ObservableSetMeasurement(obs,ObservableFastAntiferroOrder2,"mz_2",NULL);
        ObservableSetMeasurement(obs,ObservableFastAntiferroOrder4,"mz_4",NULL);
        
        for(int k=0;k<ntime*2;++k){
            for(int j=0;j<thermal;j++){
                MCDiagonalOperatorUpdate(placeholder);
                MCOffDiagOperatorUpdate(placeholder);
                MCFlipUpdate(placeholder);
                SEPlaceHolderLengthMonitor(placeholder, buffer);
            }

            for(int i_b=0;i_b<nblock;++i_b){
                placeholder->isweep=0;
                for(int j=0;j<nsweep;++j){
                    MCDiagonalOperatorUpdate(placeholder);
                    MCOffDiagOperatorUpdate(placeholder);
                    MCFlipUpdate(placeholder);

                    ObservableFastPreCal(placeholder);
                    ObservableDoMeasurement(obs,placeholder);

                    placeholder->isweep++;
                }
                ObservableShow(obs,placeholder,prefix,4);
            }
            if(k%2==1) SEPlaceHolderBetaDoubling(placeholder);
        }
        DestroyObservable(obs);
    }
/* ------------------------------------------- **
** -------------- Beta Increase -------------- **
** ------------------------------------------- */
    else if(mode==2){
        int nobs=6;
        int nave=nsweep;
        Observable *obs = CreateObservable(nobs,nave);
        ObservableSetMeasurement(obs,ObservableSpecificEnergy,"energy",NULL);
        ObservableSetMeasurement(obs,ObservableMagnetization,"magn_z",NULL);
        ObservableSetMeasurement(obs,ObservableSusceptibility,"susc_z",NULL);
        ObservableSetMeasurement(obs,ObservableFastAntiferroOrder1,"mz_1",NULL);
        ObservableSetMeasurement(obs,ObservableFastAntiferroOrder2,"mz_2",NULL);
        ObservableSetMeasurement(obs,ObservableFastAntiferroOrder4,"mz_4",NULL);
        
        for(beta=beta_i;beta<beta_f;beta+=interv){
            SEPlaceHolderSetBeta(placeholder, beta);
            SEPlaceHolderCheckSetting(placeholder);

            for(int j=0;j<thermal;j++){
                MCDiagonalOperatorUpdate(placeholder);
                MCOffDiagOperatorUpdate(placeholder);
                MCFlipUpdate(placeholder);
                SEPlaceHolderLengthMonitor(placeholder, buffer);
            }

            placeholder->isweep=0;
            for(int j=0;j<nsweep;++j){
                MCDiagonalOperatorUpdate(placeholder);
                MCOffDiagOperatorUpdate(placeholder);
                MCFlipUpdate(placeholder);

                ObservableImproveSpeedPreCal(placeholder);
                ObservableDoMeasurement(obs,placeholder);

                placeholder->isweep++;
            }
            ObservableShow(obs,placeholder,prefix,4);
        }
        DestroyObservable(obs);
    }
/* ------------------------------------------- **
** ---- Beta Increase with specific heat ----- **
** ------------------------------------------- */
    else if(mode==3){
        int nobs=9;
        int nave=nsweep;
        Observable *obs = CreateObservable(nobs,nave);
        ObservableSetMeasurement(obs,ObservableMagnetization,"magn_z",NULL);
        ObservableSetMeasurement(obs,ObservableSusceptibility,"susc_z",NULL);
        ObservableSetMeasurement(obs,ObservableFastAntiferroOrder1,"mz_1",NULL);
        ObservableSetMeasurement(obs,ObservableFastAntiferroOrder2,"mz_2",NULL);
        ObservableSetMeasurement(obs,ObservableFastAntiferroOrder4,"mz_4",NULL);
        ObservableSetMeasurement(obs,ObservableNoo1,"noo1",NULL);
        ObservableSetMeasurement(obs,ObservableNoo2,"noo2",NULL);
        ObservableSetMeasurement(obs,ObservableFastStaggeredX,"stag_x",NULL);
        ObservableSetMeasurement(obs,ObservableFastStaggeredY,"stag_y",NULL);
        
        for(beta=beta_i;beta<beta_f;beta+=interv){
            SEPlaceHolderSetBeta(placeholder, beta);
            SEPlaceHolderCheckSetting(placeholder);

            for(int j=0;j<thermal;j++){
                MCDiagonalOperatorUpdate(placeholder);
                MCOffDiagOperatorUpdate(placeholder);
                MCFlipUpdate(placeholder);
                SEPlaceHolderLengthMonitor(placeholder, buffer);
            }

            placeholder->isweep=0;
            for(int j=0;j<nsweep;++j){
                MCDiagonalOperatorUpdate(placeholder);
                MCOffDiagOperatorUpdate(placeholder);
                MCFlipUpdate(placeholder);

                ObservableImproveSpeedPreCal(placeholder);
                ObservableDoMeasurement(obs,placeholder);

                placeholder->isweep++;
            }
            ObservableShow(obs,placeholder,prefix,4);
        }
        DestroyObservable(obs);
    }


    DestroySEPlaceHolder(placeholder);
    DestroyMappingList();
}
