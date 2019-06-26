#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>

#include "MonteCarlo.h"

static int lx=8;
static int ly=8;
static int lz=0;
static int dims=2;
static int mode=0;
static int lattice=0;
static double J=1;
static double dJ=0.5;
static double p=0.5;
static double beta=4;
static double beta_i=1.0;
static double beta_f=3.0;
static double interv=0.5;
static int thermal=20000;
static int nsweep=2000;
static int nblock=50;
static int ntime=5;
static int seed=0;
static int help=0;


void SetupFromArgument(int argc, char** argv)
{
    int c;

    while((c=getopt(argc,argv,"hx:y:D:l:m:j:b:t:n:s:i:f:v:p:d:k:e:"))!=-1){
        switch(c){
            case 'h':
                help=1;
                printf("usage: \n");
                printf("\t-D <dimension> default 2\n");
                printf("\t-x <length of x> default 8\n");
                printf("\t-y <length of y> default 8\n");
                printf("\t-z <length of z> default 0\n");
                printf("\t-m <mode> default 0\n");
                printf("\t\tmode=0 : normal scheme\n");
                printf("\t\tmode=1 : zero temp scheme\n");
                printf("\t\tmode=2 : beta increase scheme\n");
                printf("\t-l <lattice> default 0\n");
                printf("\t\tlattice=0 : Herringbone\n");
                printf("\t\tlattice=1 : Plaquette\n");
                printf("\t\tlattice=2 : Configurational disorder\n");
                printf("\t-j <bond ratio> default 1\n");
                printf("\t-b <beta> default 4\n");
                printf("\t-i <beta_i>   default 1\n");
                printf("\t-f <beta_f>   default 3\n");
                printf("\t-v <interval> default 0.5\n");
                printf("\t-d <strong bond dJ> default 0.5\n");
                printf("\t-p <probability for strong bond> default 0.5\n");
                printf("\t-t <nsweep for thermal> default 20000\n");
                printf("\t-n <nsweep for estimetor> default 2000\n");
                printf("\t-k <nblock for estimetor> default 50\n");
                printf("\t-e <ntime for beta doubling> default 5\n");
                printf("\t-s <random seed> default 0\n");
                return;
            case 'x':
                lx=atoi(optarg);
                break;
            case 'y':
                ly=atoi(optarg);
                break;
            case 'z':
                lz=atoi(optarg);
                break;
            case 'D':
                dims=atoi(optarg);
                break;
            case 'm':
                mode=atoi(optarg);
                break;
            case 'l':
                lattice=atoi(optarg);
                break;
            case 'j':
                J=atof(optarg);
                break;
            case 'd':
                dJ=atof(optarg);
                break;
            case 'p':
                p=atof(optarg);
                break;
            case 'b':
                beta=atof(optarg);
                break;
            case 'i':
                beta_i=atof(optarg);
                break;
            case 'f':
                beta_f=atof(optarg);
                break;
            case 'v':
                interv=atof(optarg);
                break;
            case 't':
                thermal=atoi(optarg);
                break;
            case 'n':
                nsweep=atoi(optarg);
                break;
            case 'k':
                nblock=atoi(optarg);
                break;
            case 'e':
                ntime=atoi(optarg);
                break;
            case 's':
                seed=atoi(optarg);
                break;
            default:
            abort();
        }
    }
}

void Execute()
{
    if(help) return;
    if(dims==2){
        if(lz!=0){
            printf("Execute : For 2 dimension, lz must be zero\n");
            exit(-1);
        }
        else if(lx==0 || ly==0){
            printf("Execute : For 2 dimension, lx and ly can not be zero\n");
            exit(-1);
        }

        if(mode==0){
            int shape[2];
            shape[0]=lx;
            shape[1]=ly;
            MCDisorder2D(J,beta,shape,nsweep,thermal,seed);
        }
        else if(mode==1){
            int shape[2];
            shape[0]=lx;
            shape[1]=ly;
            MCHerringbond2D(J,beta,shape,nsweep,nblock,thermal,seed);
        }
        else if(mode==2){
            int shape[2];
            shape[0]=lx;
            shape[1]=ly;
            MCBetaIncreasePlaquetteDisorder2D(J,dJ,p,beta_i,beta_f,interv,shape,nsweep,thermal,seed);
        }
        else if(mode==3){
            int shape[2];
            shape[0]=lx;
            shape[1]=ly;
            MCBetaIncreaseConfigurationalDisorder2D(J,beta_i,beta_f,interv,shape,nsweep,thermal,seed);
        }
        else if(mode==4){
            int shape[2];
            shape[0]=lx;
            shape[1]=ly;
            MCZeroTempPlaquetteDisorder2D(J,dJ,p,beta,shape,nsweep,thermal,ntime,seed);
        }
        else if(mode==5){
            int shape[2];
            shape[0]=lx;
            shape[1]=ly;
            MCHerringbond2DImproveSpeed(J,beta,shape,nsweep,nblock,thermal,seed);
        }
        else if(mode==6){
            int shape[2];
            shape[0]=lx;
            shape[1]=ly;
            MCBetaIncreasePlaquetteDisorderImproveSpeed2D(J,dJ,p,beta_i,beta_f,interv,shape,nsweep,thermal,seed);
        }
        else if(mode==7){
            int shape[2];
            shape[0]=lx;
            shape[1]=ly;
            MCZeroTempHerringbondDisorder2D(J,dJ,p,beta,shape,nsweep,thermal,ntime,seed);
        }
        else if(mode==8){
            int shape[2];
            shape[0]=lx;
            shape[1]=ly;
            MCBetaIncreaseHerringboneDisorderImproveSpeed2D(J,dJ,p,beta_i,beta_f,interv,shape,nsweep,thermal,seed);
        }
        else{
            printf("Execute : Can not support mode=%d now\n",mode);
        }
    }
    else{
        printf("Execute : Can not support %d dimensions now\n",dims);
        exit(-1);
    }
}

#if 1
int main(int argc, char* argv[])
{
    SetupFromArgument(argc,argv);
    //Execute();
    int shape[2];
    shape[0]=lx;
    shape[1]=ly;
    
    if(help) return 0;
    else MCGeneralSchemeAndLattice(shape,mode,lattice,J,dJ,p,beta,beta_i,beta_f,interv,thermal,nsweep,nblock,ntime,seed);

    return 0;
}
#endif
