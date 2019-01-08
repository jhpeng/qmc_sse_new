#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>

#include "DataStruct.h"
#include "SEAlgorithm.h"
#include "Estimetor.h"

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

void MCDiagonalOperatorUpdateIsotropy(SEPlaceHolder* placeholder)
{
    int length=placeholder->length;
    int ndiff=placeholder->ops->ndiff;
    OperatorSequenceSort(placeholder->ops);
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
            if(dis*2*(length-noo)<beta*Nb && placeholder->lconf->sigmap->data[left]!=placeholder->lconf->sigmap->data[right]){
                placeholder->ops->sequence->data[p]=bond*ndiff;
                noo++;
            }
        }
        else if(placeholder->ops->sequence->data[p]%ndiff==0){
            dis  = gsl_rng_uniform_pos(placeholder->rng);
            if(dis*beta*Nb<2*(length-noo+1)){
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

void MCIsotropy2D(double beta, int* shape, int nsweep, int cutoff, int seed)
{
    int ndiff=2,length=50;
    int dims=2;
    double max_err=1.e-4;
    double buffer=1.3;
    char prefix[128];

    sprintf(prefix,"data/isotropy_shape_%d_%d_beta_%.1f_0",shape[0],shape[1],beta);

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

    int nobs=7;
    int nave=nsweep;
    Observable *obs = CreateObservable(nobs,nave);
    ObservableSetMeasurement(obs,ObservableSpecificEnergy,"energy",NULL);
    ObservableSetMeasurement(obs,ObservableMagnetization,"magn_z",NULL);
    ObservableSetMeasurement(obs,ObservableSusceptibility,"susc_z",NULL);
    ObservableSetMeasurement(obs,ObservableStiffnessX,"stif_x",NULL);
    ObservableSetMeasurement(obs,ObservableAntiferroOrder1,"mz_1",NULL);
    ObservableSetMeasurement(obs,ObservableAntiferroOrder2,"mz_2",NULL);
    ObservableSetMeasurement(obs,ObservableAntiferroOrder4,"mz_4",NULL);

    int j=0;
    for(j=0;j<cutoff;j++){
        MCDiagonalOperatorUpdateIsotropy(placeholder);
        MCOffDiagOperatorUpdate(placeholder);
        MCFlipUpdate(placeholder);
        SEPlaceHolderLengthMonitor(placeholder, buffer);
    }
    j=0;
    placeholder->isweep=0;
    while(j<nsweep){
        MCDiagonalOperatorUpdateIsotropy(placeholder);
        MCOffDiagOperatorUpdate(placeholder);
        MCFlipUpdate(placeholder);
        ObservableDoMeasurement(obs,placeholder);
        if((j+1)%10000==0){
            //ObservableShow(obs,placeholder,NULL,0);
            //ObservableShow(obs,placeholder,prefix,1);
            ObservableShow(obs,placeholder,prefix,2);
        }

        SEPlaceHolderLengthMonitor(placeholder, buffer);
        placeholder->isweep++;
        j++;
    }

    DestroySEPlaceHolder(placeholder);
    DestroyMappingList();
}

#if 0
int main(int argn, char *argv[])
{
    int shape[2]={8,8};
    double beta=16;
    int nsweep=1000000,cutoff=20000;
    int seed=290318;

    int i;
    for(i=0;i<argn;i++){
        if(i==0) shape[0] = atoi(argv[1]);
        else if(i==1) shape[1] = atoi(argv[2]);
    }

    for(double i=0.5;i<2.5;i=i+0.5){
        beta = i*shape[0];
        seed = seed/shape[0]*i;
        MCIsotropy2D(beta,shape,nsweep,cutoff,seed);
    }
}
#endif

#if 1
int main(int argn, char *argv[])
{
    int shape[2]={32,32};
    double beta=16;
    int nsweep=1000000,cutoff=20000;
    int seed=2318;

    MCIsotropy2D(beta,shape,nsweep,cutoff,seed);
}
#endif
