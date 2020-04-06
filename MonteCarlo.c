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

static double MCNoo,MCNoo2;
static int MCNave;

double MCImproveNoo(
                SEPlaceHolder* placeholder,
                void* args)
{
    return MCNoo;
}

double MCImproveNoo2(
                SEPlaceHolder* placeholder,
                void* args)
{
    return MCNoo2;
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

    int index=0;
    int type;

    MCNoo = 0;
    MCNoo2 = 0;
    MCNave = 0;

    for(p=0;p<length;++p){
        type = placeholder->ops->sequence->data[p]%ndiff;
        if(placeholder->ops->sequence->data[p]==-1){
            bond = (int)(gsl_rng_uniform_pos(placeholder->rng)*Nb);
            LatticeConfApplyMapping(placeholder->lconf,bond);
            left=placeholder->lconf->left;
            right=placeholder->lconf->right;
            if(placeholder->lconf->sigmap->data[left]!=placeholder->lconf->sigmap->data[right]){
                dis  = gsl_rng_uniform_pos(placeholder->rng);
                if(dis*2*(length-noo)<beta*placeholder->lconf->J->data[bond]*Nb){
                    placeholder->ops->sequence->data[p]=bond*ndiff;
                    noo++;
                    placeholder->ops->sort->data[index]=p;
                    index++;
                }
            }

            MCNoo += (double)(noo);
            MCNoo2 += (double)(noo)*(double)(noo);
            MCNave += 1;
        }
        else if(type==0){
            bond = placeholder->ops->sequence->data[p]/ndiff;
            dis  = gsl_rng_uniform_pos(placeholder->rng);
            if(dis*beta*placeholder->lconf->J->data[bond]*Nb<2*(length-noo+1)){
                placeholder->ops->sequence->data[p]=-1;
                noo--;
            }
            else{
                placeholder->ops->sort->data[index]=p;
                index++;
            }

            MCNoo += (double)(noo);
            MCNoo2 += (double)(noo)*(double)(noo);
            MCNave += 1;
        }
        else if(type==1){
            bond = placeholder->ops->sequence->data[p]/ndiff;
            LatticeConfApplyMapping(placeholder->lconf,bond);
            left=placeholder->lconf->left;
            right=placeholder->lconf->right;
            placeholder->lconf->sigmap->data[left] *=-1;
            placeholder->lconf->sigmap->data[right]*=-1;
            placeholder->ops->sort->data[index]=p;
            index++;
        }
    }

    placeholder->ops->noo=noo;
    //OperatorSequenceSort(placeholder->ops);

    MCNoo = MCNoo/MCNave;
    MCNoo2 = MCNoo2/MCNave;
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

void MCSaveConf(SEPlaceHolder* placeholder, int lx, int ly, char* prefix){
    FILE* sfile = fopen(prefix,"w");
    for(int j=0;j<ly;++j){
        for(int i=0;i<lx;++i){
            fprintf(sfile,"%d ",placeholder->lconf->sigma0->data[j*lx+i]);
        }
        fprintf(sfile,"\n");
    }
    fclose(sfile);
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
        sprintf(prefix,"data/lattice_%d_scheme_%d_shape_%d_%d_J_%.4f_b_%.4f_dJ_%.4f_p_%.4f_seed_%d",lattice,mode,shape[0],shape[1],J,beta,dJ,p,seed);

        SEPlaceHolderSetBeta(placeholder, beta);
        SEPlaceHolderCheckSetting(placeholder);

        int nobs=14;
        int nave=nsweep;
        Observable *obs = CreateObservable(nobs,nave);
        ObservableSetMeasurement(obs,ObservableSpecificEnergy,"energy",NULL);
        ObservableSetMeasurement(obs,ObservableMagnetization,"magn_z",NULL);
        ObservableSetMeasurement(obs,ObservableFastStaggeredX,"stag_x",NULL);
        ObservableSetMeasurement(obs,ObservableFastStaggeredY,"stag_y",NULL);
        ObservableSetMeasurement(obs,ObservableSusceptibility,"susc_z",NULL);
        ObservableSetMeasurement(obs,ObservableFastStiffnessX,"stif_x",NULL);
        ObservableSetMeasurement(obs,ObservableFastStiffnessY,"stif_y",NULL);
        ObservableSetMeasurement(obs,ObservableFastAntiferroOrder1,"mz_1",NULL);
        ObservableSetMeasurement(obs,ObservableFastAntiferroOrder2,"mz_2",NULL);
        ObservableSetMeasurement(obs,ObservableFastAntiferroOrder4,"mz_4",NULL);
        ObservableSetMeasurement(obs,ObservableNoo1,"noo1",NULL);
        ObservableSetMeasurement(obs,ObservableNoo2,"noo2",NULL);
        //ObservableSetMeasurement(obs,MCImproveNoo,"noo1",NULL);
        //ObservableSetMeasurement(obs,MCImproveNoo2,"noo2",NULL);
        
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
                //SESaveConfiguration(placeholder,prefix);

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
/* ----------------------------------------------------- **
** ------------ Quantum Correlator Zero Temp ----------- **
** ------------------------------------------------------ */
    else if(mode==4){
        SEPlaceHolderSetBeta(placeholder, beta);
        SEPlaceHolderCheckSetting(placeholder);

        int nobs=13;
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
        ObservableSetMeasurement(obs,ObservableFastStaggeredX,"stag_x",NULL);
        ObservableSetMeasurement(obs,ObservableFastStaggeredY,"stag_y",NULL);
        ObservableSetMeasurement(obs,ObservableNoo1,"noo1",NULL);
        ObservableSetMeasurement(obs,ObservableNoo2,"noo2",NULL);
        ObservableSetMeasurement(obs,ObservableFastQuantumVariance,"var_q",NULL);
        
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

                    ObservableQuantumCorrelator(placeholder);
                    ObservableDoMeasurement(obs,placeholder);

                    if(k==(ntime*2-1) && (j%10)==9){
                        char conf_prefix[128];
                        sprintf(conf_prefix,"conf/conf_lattice_%d_scheme_%d_shape_%d_%d_J_%.4f_dJ_%.4f_p_%.4f_id_%d_seed_%d.txt",
                                                lattice,mode,shape[0],shape[1],J,dJ,p,(i_b*nsweep+j),seed);

                        MCSaveConf(placeholder,shape[0],shape[1],conf_prefix);
                    }

                    placeholder->isweep++;
                }
                ObservableShow(obs,placeholder,prefix,4);
            }
            if(k%2==1) SEPlaceHolderBetaDoubling(placeholder);
        }
        DestroyObservable(obs);
    }
/* ------------------------------------------- **
** ---- Beta Increase with specific heat ----- **
** ------------------------------------------- */
    else if(mode==5){
        
        for(beta=beta_i;beta<beta_f;beta+=interv){
            int nobs=11;
            int nave=(int)(nsweep*beta*beta);
            if(nave<thermal) nave=thermal;
            SEPlaceHolderSetNsweep(placeholder, nave, thermal);

            Observable *obs = CreateObservable(nobs,nave);
            ObservableSetMeasurement(obs,ObservableMagnetization,"magn_z",NULL);
            ObservableSetMeasurement(obs,ObservableSusceptibility,"susc_z",NULL);
            ObservableSetMeasurement(obs,ObservableFastAntiferroOrder1,"mz_1",NULL);
            ObservableSetMeasurement(obs,ObservableFastAntiferroOrder2,"mz_2",NULL);
            ObservableSetMeasurement(obs,ObservableFastAntiferroOrder4,"mz_4",NULL);
            ObservableSetMeasurement(obs,ObservableNoo1,"noo1",NULL);
            ObservableSetMeasurement(obs,ObservableNoo2,"noo2",NULL);
            //ObservableSetMeasurement(obs,MCImproveNoo,"noo1",NULL);
            //ObservableSetMeasurement(obs,MCImproveNoo2,"noo2",NULL);
            ObservableSetMeasurement(obs,ObservableFastStaggeredX,"stag_x",NULL);
            ObservableSetMeasurement(obs,ObservableFastStaggeredY,"stag_y",NULL);
            SEPlaceHolderSetBeta(placeholder, beta);
            SEPlaceHolderCheckSetting(placeholder);

            for(int j=0;j<thermal;j++){
                MCDiagonalOperatorUpdate(placeholder);
                MCOffDiagOperatorUpdate(placeholder);
                MCFlipUpdate(placeholder);
                SEPlaceHolderLengthMonitor(placeholder, buffer);
            }

            placeholder->isweep=0;
            for(int j=0;j<nave;++j){
                MCDiagonalOperatorUpdate(placeholder);
                MCOffDiagOperatorUpdate(placeholder);
                MCFlipUpdate(placeholder);

                ObservableImproveSpeedPreCal(placeholder);
                ObservableDoMeasurement(obs,placeholder);

                placeholder->isweep++;
            }
            ObservableShow(obs,placeholder,prefix,4);

            for(int j=0;j<thermal;j++){
                MCDiagonalOperatorUpdate(placeholder);
                MCOffDiagOperatorUpdate(placeholder);
                MCFlipUpdate(placeholder);
                SEPlaceHolderLengthMonitor(placeholder, buffer);
            }

            placeholder->isweep=0;
            for(int j=0;j<nave;++j){
                MCDiagonalOperatorUpdate(placeholder);
                MCOffDiagOperatorUpdate(placeholder);
                MCFlipUpdate(placeholder);

                ObservableImproveSpeedPreCal(placeholder);
                ObservableDoMeasurement(obs,placeholder);

                placeholder->isweep++;
            }
            ObservableShow(obs,placeholder,prefix,4);
            DestroyObservable(obs);
        }
    }


    DestroySEPlaceHolder(placeholder);
    DestroyMappingList();
}
