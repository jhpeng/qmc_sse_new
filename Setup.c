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
static int nblock=1;
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
                printf("\t\tmode=3 : beta increase scheme(specific heat)\n");
                printf("\t\tmode=4 : quantum correlator(zero temp scheme)\n");
                printf("\t\tmode=5 : beta increase scheme(specific heat, increase nsweep)\n");
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

#if 1
int main(int argc, char* argv[])
{
    SetupFromArgument(argc,argv);
    //Execute();
    int shape[2];

    //Test for the speed
    if(0){
        lx = 32;
        ly = 32;
        mode = 0;
        J = 2.4981;
        dJ = 0;
        beta = 10;
        thermal = 20000;
        nsweep = 0;
        nblock = 10;
        seed = 39829;
    }
    
    shape[0]=lx;
    shape[1]=ly;
    if(help) return 0;
    else MCGeneralSchemeAndLattice(shape,mode,lattice,J,dJ,p,beta,beta_i,beta_f,interv,thermal,nsweep,nblock,ntime,seed);

    return 0;
}
#endif
