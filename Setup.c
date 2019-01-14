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
static double J=1;
static double beta=4;
static int thermal=20000;
static int nsweep=100000;
static int seed=0;
static int help=0;

void SetupFromArgument(int argc, char** argv)
{
    int c;

    while((c=getopt(argc,argv,"hx:y:d:m:j:b:t:n:s:"))!=-1){
        switch(c){
            case 'h':
                help=1;
                printf("usage: \n");
                printf("\t-d <dimension> default 2\n");
                printf("\t-x <length of x> default 8\n");
                printf("\t-y <length of y> default 8\n");
                printf("\t-z <length of z> default 0\n");
                printf("\t-m <mode> default 0\n");
                printf("\t\tmode=0 : disorder\n");
                printf("\t\tmode=1 : herringbond\n");
                printf("\t\tmode=2 : clean\n");
                printf("\t-j <bond ratio> default 1\n");
                printf("\t-b <beta> default 4\n");
                printf("\t-t <nsweep for thermal> default 20000\n");
                printf("\t-n <nsweep for estimetor> default 1000000\n");
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
            case 'd':
                dims=atoi(optarg);
                break;
            case 'm':
                mode=atoi(optarg);
                break;
            case 'j':
                J=atof(optarg);
                break;
            case 'b':
                beta=atof(optarg);
                break;
            case 't':
                thermal=atoi(optarg);
                break;
            case 'n':
                nsweep=atoi(optarg);
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
            MCHerringbond2D(J,beta,shape,nsweep,thermal,seed);
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

int main(int argc, char* argv[])
{
    SetupFromArgument(argc,argv);
    Execute();
    return 0;
}
