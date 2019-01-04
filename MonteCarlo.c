#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>

#include "DataStruct.h"
#include "SEAlgorithm.h"

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

#if 1
int main()
{
    int ndiff=2,length=50;
    int dims=2,shape[2]={48,48};
    int nsweep=1000000,seed=2133124;
    int cutoff=2000;
    double beta=512;
    double max_err=1.e-4;
    double buffer=1.3;

    SEPlaceHolder* placeholder = CreateSEPlaceHolder();
    SEPlaceHolderSetLattice(placeholder,mapping_2d,shape,dims,0);
    SEPlaceHolderSetLength(placeholder,length,ndiff);
    SEPlaceHolderSetRandomSeed(placeholder, seed);
    SEPlaceHolderSetNsweep(placeholder, nsweep, cutoff);
    SEPlaceHolderSetBeta(placeholder, beta);
    SEPlaceHolderSetError(placeholder, max_err);
    SEPlaceHolderCheckSetting(placeholder);

    MCInitializeLatticeConf(placeholder);
    int j=0;
    while(j<1000){
        MCDiagonalOperatorUpdateIsotropy(placeholder);
        MCOffDiagOperatorUpdate(placeholder);
        MCFlipUpdate(placeholder);
        if(j%1000==0){
            printf("%d %d %d\n",j,placeholder->ops->noo,placeholder->length);
#if 0
            for(int i=0;i<length;++i) printf("%d ",placeholder->ops->sequence->data[i]);
            printf("\n");
#endif
        }

        SEPlaceHolderLengthMonitor(placeholder, buffer);
        j++;
    }

    DestroySEPlaceHolder(placeholder);
}
#endif
