#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "DataStruct.h"
#include "SEAlgorithm.h"
#include "HTML.h"

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
                    char* prefix,
                    int mode)
{
    int isweep = placeholder->isweep+1;
    int nsweep = placeholder->nsweep;
    double ratio = (double)isweep/nsweep*100;
    char filename[128];
    char name[]={"         "};
    obs->end=clock();
    double dtime = difftime(obs->end,obs->start)/CLOCKS_PER_SEC;
    MeanAverage(obs);
    //if model=0 : std output
    //if model=1 : output a text file
    //if model=2 : output a html file

    if(mode==0){
            printf("Measurement ...\n");
            if(placeholder->lconf->dims==1){
                printf("the shape of lattice \t: [%d]\n",placeholder->lconf->shape[0]);
            }
            else if(placeholder->lconf->dims==2){
                printf("the shape of lattice \t: [%d,%d]\n",placeholder->lconf->shape[0],placeholder->lconf->shape[1]);
            }
            else if(placeholder->lconf->dims==3){
                printf("the shape of lattice \t: [%d,%d,%d]\n",placeholder->lconf->shape[0],placeholder->lconf->shape[1],placeholder->lconf->shape[2]);
            }
            printf("inverse temperature \t: %lf\n",placeholder->beta);
            printf("length of sequence \t: %d \nnumber of operator \t: %d\n",placeholder->ops->length,placeholder->ops->noo);
            printf("obs name  \t| mean   \t| var    \t| err    \t|\n");
            for(int i_obs=0;i_obs<obs->nobs;i_obs++){
                memcpy(name,obs->obs_name[i_obs],6);
                printf("%s  \t| %.4e \t| %.4e \t| %.4e \t|\n",name,obs->mean[i_obs],obs->var[i_obs],obs->err[i_obs]);
            }
            printf("|<-------------------%.2lf------------------->| %d/%d time : %.1lf(s) \n",ratio,isweep,nsweep,dtime);
            printf("===========================================================================\n");
     }
    else if(mode==1){
            sprintf(filename,"%s.out",prefix);
            FILE* outfile = fopen(filename,"w");
            
            fprintf(outfile,"Measurement ...\n");
            if(placeholder->lconf->dims==1){
                fprintf(outfile,"the shape of lattice \t: [%d]\n",placeholder->lconf->shape[0]);
            }
            else if(placeholder->lconf->dims==2){
                fprintf(outfile,"the shape of lattice \t: [%d,%d]\n",placeholder->lconf->shape[0],placeholder->lconf->shape[1]);
            }
            else if(placeholder->lconf->dims==3){
                fprintf(outfile,"the shape of lattice \t: [%d,%d,%d]\n",placeholder->lconf->shape[0],placeholder->lconf->shape[1],placeholder->lconf->shape[2]);
            }
            fprintf(outfile,"inverse temperature \t: %lf\n",placeholder->beta);
            fprintf(outfile,"length of sequence \t: %d \nnumber of operator \t: %d\n",placeholder->ops->length,placeholder->ops->noo);
            fprintf(outfile,"obs name  \t| mean   \t| var    \t| err    \t|\n");
            for(int i_obs=0;i_obs<obs->nobs;i_obs++){
                memcpy(name,obs->obs_name[i_obs],6);
                fprintf(outfile,"%s  \t| %.4e \t| %.4e \t| %.4e \t|\n",name,obs->mean[i_obs],obs->var[i_obs],obs->err[i_obs]);
            }
            fprintf(outfile,"|<-------------------%.2lf------------------->| %d/%d time : %.1lf(s) \n",ratio,isweep,nsweep,dtime);
            fprintf(outfile,"===========================================================================\n");
            fclose(outfile);
    }
    else if(mode==2){
        OutputHTML(obs,placeholder,prefix);
    }
    else if(mode==3){
        sprintf(filename,"%s.data",prefix);
        FILE* outfile = fopen(filename,"w");
        int i_obs;
        for(i_obs=0;i_obs<obs->nobs;++i_obs){
            fprintf(outfile,"%e ",obs->mean[i_obs]);
        }
        fprintf(outfile,"\n");
        for(i_obs=0;i_obs<obs->nobs;++i_obs){
            fprintf(outfile,"%e ",obs->err[i_obs]);
        }
        fprintf(outfile,"\n");
        fclose(outfile);
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
    int type,bond,p,left,right,length=placeholder->ops->length;
    int ndiff=placeholder->ops->ndiff;
    int nsite=placeholder->lconf->nsite;
    double winding=0;
    double beta = placeholder->beta;

    LatticeConfSynchronizeSigma(placeholder->lconf);

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

double ObservableAntiferroOrder1(
                    SEPlaceHolder* placeholder,
                    void* args)
{
    int i,j,id,p,bond;
    int left,right;
    int nsite  = placeholder->lconf->nsite;
    int length = placeholder->ops->length;
    double mz=0;
    double m1=0;

    LatticeConfSynchronizeSigma(placeholder->lconf);

    if(placeholder->lconf->dims==2){
        int xlen = placeholder->lconf->shape[0];
        mz=0;
        for(id=0;id<nsite;++id){
            i = id/xlen;
            j = id%xlen; 
            mz+= (((i+j)%2)*2-1)*placeholder->lconf->sigmap->data[id];
        }
        for(p=0;p<length;++p){
            m1+=fabs(mz);
            if(placeholder->ops->sequence->data[p]!=-1){
            if((placeholder->ops->sequence->data[p]%2)==1){
                bond = placeholder->ops->sequence->data[p]/2;
                LatticeConfApplyMapping(placeholder->lconf,bond);
                left  = placeholder->lconf->left;
                right = placeholder->lconf->right;
                i = left/xlen;
                j = left%xlen;
                mz-=4*(((i+j)%2)*2-1)*placeholder->lconf->sigmap->data[left];
                placeholder->lconf->sigmap->data[left]*=-1;
                placeholder->lconf->sigmap->data[right]*=-1;
            }
            }
        }
    }

    m1 = m1/length/nsite*0.5;

    return m1;
}

double ObservableAntiferroOrder2(
                    SEPlaceHolder* placeholder,
                    void* args)
{
    int i,j,id,p,bond;
    int left,right;
    int nsite  = placeholder->lconf->nsite;
    int length = placeholder->ops->length;
    double mz=0;
    double m2=0;

    LatticeConfSynchronizeSigma(placeholder->lconf);

    if(placeholder->lconf->dims==2){
        int xlen = placeholder->lconf->shape[0];
        mz=0;
        for(id=0;id<nsite;++id){
            i = id/xlen;
            j = id%xlen; 
            mz+= (((i+j)%2)*2-1)*placeholder->lconf->sigmap->data[id];
        }
        for(p=0;p<length;++p){
            m2+=mz*mz;
            if(placeholder->ops->sequence->data[p]!=-1){
            if((placeholder->ops->sequence->data[p]%2)==1){
                bond = placeholder->ops->sequence->data[p]/2;
                LatticeConfApplyMapping(placeholder->lconf,bond);
                left  = placeholder->lconf->left;
                right = placeholder->lconf->right;
                i = left/xlen;
                j = left%xlen;
                mz-=4*(((i+j)%2)*2-1)*placeholder->lconf->sigmap->data[left];
                placeholder->lconf->sigmap->data[left]*=-1;
                placeholder->lconf->sigmap->data[right]*=-1;
            }
            }
        }
    }

    m2 = m2/length/nsite/nsite*0.25;

    return m2;
}

double ObservableAntiferroOrder4(
                    SEPlaceHolder* placeholder,
                    void* args)
{
    int i,j,id,p,bond;
    int left,right;
    int nsite  = placeholder->lconf->nsite;
    int length = placeholder->ops->length;
    double mz=0;
    double m4=0;

    LatticeConfSynchronizeSigma(placeholder->lconf);

    if(placeholder->lconf->dims==2){
        int xlen = placeholder->lconf->shape[0];
        mz=0;
        for(id=0;id<nsite;++id){
            i = id/xlen;
            j = id%xlen; 
            mz+= (((i+j)%2)*2-1)*placeholder->lconf->sigmap->data[id];
        }
        for(p=0;p<length;++p){
            m4+=mz*mz*mz*mz;
            if(placeholder->ops->sequence->data[p]!=-1){
            if((placeholder->ops->sequence->data[p]%2)==1){
                bond = placeholder->ops->sequence->data[p]/2;
                LatticeConfApplyMapping(placeholder->lconf,bond);
                left  = placeholder->lconf->left;
                right = placeholder->lconf->right;
                i = left/xlen;
                j = left%xlen;
                mz-=4*(((i+j)%2)*2-1)*placeholder->lconf->sigmap->data[left];
                placeholder->lconf->sigmap->data[left]*=-1;
                placeholder->lconf->sigmap->data[right]*=-1;
            }
            }
        }
    }

    m4 = m4/length/nsite/nsite/nsite/nsite*0.0625;

    return m4;
}

static double obs_m1,obs_m2,obs_m4,obs_stifx;
void ObservableFastPreCal(
                    SEPlaceHolder* placeholder)
{
    int i,j,id,p,bond,type;
    int left,right;
    int nsite  = placeholder->lconf->nsite;
    int length = placeholder->ops->length;
    int ndiff  = placeholder->ops->ndiff;
    double beta = placeholder->beta;
    double mz=0;
    double m1=0;
    double m2=0;
    double m4=0;
    double winding=0;

    LatticeConfSynchronizeSigma(placeholder->lconf);

    if(placeholder->lconf->dims==2){
        int xlen = placeholder->lconf->shape[0];
        mz=0;
        for(id=0;id<nsite;++id){
            i = id/xlen;
            j = id%xlen; 
            mz+= (((i+j)%2)*2-1)*placeholder->lconf->sigmap->data[id];
        }
        for(p=0;p<length;++p){
            m1+=fabs(mz);
            m2+=mz*mz;
            m4+=mz*mz*mz*mz;
            if(placeholder->ops->sequence->data[p]!=-1){
                type = placeholder->ops->sequence->data[p]%ndiff;
                if(type==1){
                    bond = placeholder->ops->sequence->data[p]/ndiff;
                    LatticeConfApplyMapping(placeholder->lconf,bond);
                    left  = placeholder->lconf->left;
                    right = placeholder->lconf->right;
                    i = left/xlen;
                    j = left%xlen;
                    mz-=4*(((i+j)%2)*2-1)*placeholder->lconf->sigmap->data[left];
                    if(bond<nsite){
                        winding += placeholder->lconf->sigmap->data[left];
                    }
                    placeholder->lconf->sigmap->data[left]*=-1;
                    placeholder->lconf->sigmap->data[right]*=-1;
                }
            }
        }
    }

    obs_m1 = m1/length/nsite*0.5;
    obs_m2 = m2/length/nsite/nsite*0.25;
    obs_m4 = m4/length/nsite/nsite/nsite/nsite*0.0625;
    obs_stifx = winding*winding/placeholder->lconf->shape[0]/placeholder->lconf->shape[0]/beta;
}

double ObservableFastAntiferroOrder1(
                    SEPlaceHolder* placeholder,
                    void* args)
{
    return obs_m1;
}

double ObservableFastAntiferroOrder2(
                    SEPlaceHolder* placeholder,
                    void* args)
{
    return obs_m2;
}

double ObservableFastAntiferroOrder4(
                    SEPlaceHolder* placeholder,
                    void* args)
{
    return obs_m4;
}

double ObservableFastStiffnessX(
                    SEPlaceHolder* placeholder,
                    void* args)
{
    return obs_stifx;
}
