#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>

#include "DataStruct.h"
#include "SEAlgorithm.h"


void LatticeConfSetMapping(LatticeConf* lconf, bond2sigma *mapping)
{
    lconf->mapping = mapping;
}

void LatticeConfApplyMapping(LatticeConf* lconf, int bond)
{
    lconf->mapping(&lconf->left,&lconf->right,lconf->shape,bond);
}

void LatticeConfSynchronizeSigma(LatticeConf* lconf)
{
    for(int i=0;i<lconf->nsite;++i) lconf->sigmap->data[i]=lconf->sigma0->data[i];
}

void OperatorSequenceSort(OperatorSequence* ops)
{
    int p;
    for(p=0;p<ops->length;++p) ops->sort->data[p]=-1;
    ops->noo=0;
    for(p=0;p<ops->length;++p){
        if(ops->sequence->data[p]!=-1){
            ops->sort->data[ops->noo]=p;
            ++ops->noo;
        }
    }
}

void SEOps2Lvc(
                    OperatorLoop* opl,
                    LatticeConf* lconf,
                    OperatorSequence* ops)
{
    if(ops->length!=opl->length){
        printf("ops2lvc : length error\n");
        exit(-1);
    }

    int i,p,bond,length=ops->length;
    int noo=ops->noo;
    int ndiff=ops->ndiff;
    int v0,v1,v2;
    
    for(i=0;i<lconf->nsite;++i){
        lconf->first->data[i]=-1;
        lconf->last->data[i] =-1;
    }
    for(i=0;i<4*length;++i){
        opl->lvc->data[i]=-1;
        opl->cross->data[i]=-1;
    }
    for(i=0;i<noo;++i){
        p  = ops->sort->data[i];
        v0 = 4*p;
        bond = ops->sequence->data[p]/ndiff;
        LatticeConfApplyMapping(lconf, bond);
        v1 = lconf->last->data[lconf->left];
        v2 = lconf->last->data[lconf->right];
        if(v1!=-1){
            opl->lvc->data[v0] = v1;
            opl->lvc->data[v1] = v0;
            lconf->last->data[lconf->left] = v0+2;
        }
        else{
            lconf->first->data[lconf->left] = v0;
            lconf->last->data[lconf->left]  = v0+2;
        }
    
        if(v2!=-1){
            opl->lvc->data[v0+1] = v2;
            opl->lvc->data[v2] = v0+1;
            lconf->last->data[lconf->right]  = v0+3;
        }
        else{
            lconf->first->data[lconf->right] = v0+1;
            lconf->last->data[lconf->right]  = v0+3;
        }
    }
    for(i=0;i<lconf->nsite;++i){
        if(lconf->last->data[i]!=-1){
            opl->lvc->data[lconf->last->data[i]]=lconf->first->data[i];
            opl->lvc->data[lconf->first->data[i]]=lconf->last->data[i];
            opl->cross->data[lconf->last->data[i]]=i;
            opl->cross->data[lconf->first->data[i]]=i;
        }
    }
}

void SETraverseLoop(
                    OperatorLoop* opl,
                    LatticeConf* lconf,
                    OperatorSequence* ops)
{
    if(opl->check==0) return;
    else if(opl->check==-1){
        opl->v0=0;
        opl->check=1;
    }
    if(ops->length!=opl->length){
        printf("traverse_loop : length error");
        exit(-1);
    }
    
    int v,length=ops->length;
    lconf->nflip=0;
    opl->noo=0;
    //for(int i=0;i<2*length;++i) opl->loop->data[i]=-1;
    //for(int i=0;i<lconf->nsite;++i) lconf->flip->data[i]=-1;
    for(v=opl->v0;v<4*length;++v){
        if(opl->lvc->data[v]!=-1){
            int check=1;
            int prev,next,flip;
            prev = v^1;
            while(check){
                opl->loop->data[opl->noo] = prev/4;
                opl->noo++;
                if(opl->cross->data[prev]!=-1){
                    flip=opl->cross->data[prev];
                    lconf->flip->data[lconf->nflip]=flip;
                    lconf->nflip++;
                }
                next = opl->lvc->data[prev];
                opl->lvc->data[prev]=-1;
                opl->lvc->data[prev^1]=-1;
                prev = next^1;
                if(opl->lvc->data[prev]==-1) check=0;
            }
            opl->v0=v+1;
            return;
        }
    }
    opl->check=0;
}

void SELoopUpdate(
                    LatticeConf* lconf,
                    OperatorSequence* ops,
                    OperatorLoop* opl)
{
    int i,p,type,bond;
    for(i=0;i<opl->noo;++i){
        p = opl->loop->data[i];
        type = ops->sequence->data[p]%ops->ndiff;
        bond = ops->sequence->data[p]/ops->ndiff;
        if(type==0 || type==1){
                type = type^1;
                ops->sequence->data[p] = bond*ops->ndiff+type;
        }
    }
    
    for(i=0;i<lconf->nflip;i++){
        lconf->sigma0->data[lconf->flip->data[i]] *=-1;
    }
}

void mapping_1d(int* left, int* right, const int* shape, int bond)
{
    if(bond==0){
        *left = shape[0]-1;
        *right = 0;
    }
    else if(bond>0 && bond<shape[0]){
        *left = bond-1;
        *right = bond;
    }
    else{
        printf("mapping_1d : input domain error\n");
        exit(-1);
    }
}

static int* BOND2LEFT;
static int* BOND2RIGHT;

void CreateMappingList(bond2sigma* mapping, const int* shape, int Nb)
{
    int left,right;
    BOND2LEFT =(int*)malloc(Nb*sizeof(int));
    BOND2RIGHT=(int*)malloc(Nb*sizeof(int));
    for(int i=0;i<Nb;++i){
        mapping(&left,&right,shape,i);
        BOND2LEFT[i] =left;
        BOND2RIGHT[i]=right;
    }
}

void DestroyMappingList()
{
    free(BOND2LEFT);
    free(BOND2RIGHT);
}

void mapping_list(int* left, int* right, const int* shape, int bond)
{
    *left  = BOND2LEFT[bond];
    *right = BOND2RIGHT[bond];
}

#if 1
#define FAST_MAPPING_2D
#endif
#ifdef FAST_MAPPING_2D
void mapping_2d(int* left, int* right, const int* shape, int bond)
{
    int i,j,k,t;
    int x=shape[0];
    int y=shape[1];
    int nsite=x*y;

    t=bond/nsite;
    k=bond%nsite;
    *left=k;

    i= k%x;
    j= k/x;


    i+=t^1;
    j+=t;
    i=i%x;
    j=j%y;
    *right=j*x+i;
}
#endif
#ifndef FAST_MAPPING_2D
void mapping_2d(int *first, int *second, const int *size, int bond)
{
    int M=size[0];
    int N=size[1];
    int nsite=M*N;
    int t=bond/nsite;

    if(t==0){
        if(bond%M==(M-1)){
            *first=bond;
            *second=M*(int)(bond/M);
        }
        else{
            *first=bond;
            *second=bond+1;
        }
    }
    else if(t==1){
        int k=bond-nsite;
        int p=k/M;
        if(p==(N-1)){
            *first=k;
            *second=k%M;
        }
        else{
            *first=k;
            *second=k+M;
        }
    }
    else{
        printf("bond is not satisfy the size! %d\n",bond);
        exit(-1);
    }
}
#endif

void SEPlaceHolderSetLattice(
                    SEPlaceHolder* placeholder, 
                    bond2sigma *mapping, 
                    const int* shape, 
                    int dims, 
                    int model)
{
    if(placeholder->set_lattice){
        printf("Create the new system and overwrite the original lconf\n");
        DestroyLatticeConf(placeholder->lconf);
    }

    if(model==0){
        placeholder->set_lattice=1;
        placeholder->lconf = CreateSquareLatticeConf(shape,dims);
        placeholder->lconf->mapping = mapping;
    }
    else{
        printf("No such model %d\n",model);
        exit(-1);
    }
}

void SEPlaceHolderSetLength(SEPlaceHolder* placeholder, int length, int ndiff)
{
    if(placeholder->set_length && placeholder->length<length){
        DestroyOperatorLoop(placeholder->opl);
        placeholder->opl = CreateOperatorLoop(length);
        
        OperatorSequence* ops = CreateOperatorSequence(length,ndiff);
        ops->noo = placeholder->ops->noo;
        for(int i=0;i<length;++i){
            if(i<placeholder->length){
                ops->sequence->data[i] = placeholder->ops->sequence->data[i];
                ops->sort->data[i] = placeholder->ops->sort->data[i];
            }
            else{
                ops->sequence->data[i] = -1;
                ops->sort->data[i] = -1;
            }
        }
        placeholder->length=length;
        DestroyOperatorSequence(placeholder->ops);
        placeholder->ops = ops;
        
    }
    else{
        placeholder->set_length=1;
        placeholder->length=length;
        placeholder->ops = CreateOperatorSequence(length,ndiff);
        placeholder->opl = CreateOperatorLoop(length);
    }
}

void SEPlaceHolderSetRandomSeed(SEPlaceHolder* placeholder, int seed)
{
    if(placeholder->set_random){
        gsl_rng_free(placeholder->rng);
        placeholder->rng = gsl_rng_alloc(gsl_rng_taus);
        gsl_rng_set(placeholder->rng,seed);
    }
    else{
        placeholder->set_random=1;
        placeholder->rng = gsl_rng_alloc(gsl_rng_taus);
        gsl_rng_set(placeholder->rng,seed);
    }
}

void SEPlaceHolderSetNsweep(SEPlaceHolder* placeholder, int nsweep, int cutoff)
{
    placeholder->nsweep=nsweep;
    placeholder->cutoff=cutoff;
}

void SEPlaceHolderSetBeta(SEPlaceHolder* placeholder, double beta)
{
    placeholder->beta=beta;
}

void SEPlaceHolderSetError(SEPlaceHolder* placeholder, double max_err)
{
    placeholder->max_err=max_err;
}

int SEPlaceHolderCheckSetting(SEPlaceHolder* placeholder)
{
    int check=1;
    if(placeholder->set_lattice^1){
        printf("SEPlaceHolderCheckSetting : please use SEPlaceHolderSetLattice\n");
        check=0;
    }
    
    if(placeholder->set_length^1){
        printf("SEPlaceHolderCheckSetting : please use SEPlaceHolderSetLength\n");
        check=0;
    }

    if(placeholder->set_random^1){
        printf("SEPlaceHolderCheckSetting : please use SEPlaceHolderSetRandomSeed\n");
        check=0;
    }

    if(placeholder->nsweep==0){
        printf("SEPlaceHolderCheckSetting : please use SEPlaceHolderSetNsweep\n");
        check=0;
    }

    if(placeholder->beta==-1){
        printf("SEPlaceHolderCheckSetting : please use SEPlaceHolderSetBeta\n");
        check=0;
    }
    
    if(placeholder->max_err==-1){
        printf("SEPlaceHolderCheckSetting : please use SEPlaceHolderSetError\n");
        check=0;
    }

    return check;
}

void SEPlaceHolderLengthMonitor(SEPlaceHolder* placeholder, double buffer)
{
    int noo=placeholder->ops->noo;
    int ndiff=placeholder->ops->ndiff;
    if(noo*buffer>placeholder->length){
        SEPlaceHolderSetLength(placeholder,(int)(noo*buffer),ndiff);
    }
}

void SEPlaceHolderSetDisorder2D(SEPlaceHolder* placeholder, double J)
{
    int Nb=placeholder->lconf->Nb;
    int nsite=placeholder->lconf->nsite;
    int bond,i,j;
    if(placeholder->lconf->dims!=2){
        printf("SEPlaceHolderSetDisorder2D : The dimension must be 2!\n");
        exit(-1);
    }
    else if(Nb!=(2*nsite)){
        printf("SEPlaceHolderSetDisorder2D : Nb must be 2*nsite!\n");
        exit(-1);
    }
    else if((placeholder->lconf->shape[0]%4)!=0 && (placeholder->lconf->shape[0]%4)!=0){
        printf("SEPlaceHolderSetDisorder2D : lx and ly must be multiple of 4!\n");
        exit(-1);
    }

    for(bond=0;bond<nsite;++bond){
        placeholder->lconf->J->data[bond]=1;
        i = (bond%placeholder->lconf->shape[0])%4;
        j = (bond/placeholder->lconf->shape[0])%4;
        if(i==0 && j==0) placeholder->lconf->J->data[bond]=J;
        else if(i==0 && j==1) placeholder->lconf->J->data[bond]=J;
        else if(i==2 && j==2) placeholder->lconf->J->data[bond]=J;
        else if(i==2 && j==3) placeholder->lconf->J->data[bond]=J;
    }
    for(bond=nsite;bond<(2*nsite);++bond){
        placeholder->lconf->J->data[bond]=1;
        i = ((bond-nsite)%placeholder->lconf->shape[0])%4;
        j = ((bond-nsite)/placeholder->lconf->shape[0])%4;
        if(i==0 && j==2) placeholder->lconf->J->data[bond]=J;
        else if(i==1 && j==2) placeholder->lconf->J->data[bond]=J;
        else if(i==2 && j==0) placeholder->lconf->J->data[bond]=J;
        else if(i==3 && j==0) placeholder->lconf->J->data[bond]=J;
    }
}

void SEPlaceHolderSetHerringbond2D(SEPlaceHolder* placeholder, double J)
{
    int Nb=placeholder->lconf->Nb;
    int nsite=placeholder->lconf->nsite;
    int bond,i,j;
    if(placeholder->lconf->dims!=2){
        printf("SEPlaceHolderSetDisorder2D : The dimension must be 2!\n");
        exit(-1);
    }
    else if(Nb!=(2*nsite)){
        printf("SEPlaceHolderSetDisorder2D : Nb must be 2*nsite!\n");
        exit(-1);
    }
    else if((placeholder->lconf->shape[0]%4)!=0 && (placeholder->lconf->shape[0]%4)!=0){
        printf("SEPlaceHolderSetDisorder2D : lx and ly must be multiple of 4!\n");
        exit(-1);
    }

    for(bond=0;bond<nsite;++bond){
        placeholder->lconf->J->data[bond]=1;
        i = (bond%placeholder->lconf->shape[0])%4;
        j = (bond/placeholder->lconf->shape[0])%4;
        if(i==1 && j==0) placeholder->lconf->J->data[bond]=J;
        else if(i==2 && j==1) placeholder->lconf->J->data[bond]=J;
        else if(i==3 && j==2) placeholder->lconf->J->data[bond]=J;
        else if(i==0 && j==3) placeholder->lconf->J->data[bond]=J;
    }
    for(bond=nsite;bond<(2*nsite);++bond){
        placeholder->lconf->J->data[bond]=1;
        i = ((bond-nsite)%placeholder->lconf->shape[0])%4;
        j = ((bond-nsite)/placeholder->lconf->shape[0])%4;
        if(i==0 && j==0) placeholder->lconf->J->data[bond]=J;
        else if(i==1 && j==1) placeholder->lconf->J->data[bond]=J;
        else if(i==2 && j==2) placeholder->lconf->J->data[bond]=J;
        else if(i==3 && j==3) placeholder->lconf->J->data[bond]=J;
    }
}

void SEPlaceHolderSetHerringbondRandom2D(SEPlaceHolder* placeholder, double Jc, double dJ, double p)
{
    int Nb=placeholder->lconf->Nb;
    int nsite=placeholder->lconf->nsite;
    int bond,i,j;
    double J;
    if(placeholder->lconf->dims!=2){
        printf("SEPlaceHolderSetDisorder2D : The dimension must be 2!\n");
        exit(-1);
    }
    else if(Nb!=(2*nsite)){
        printf("SEPlaceHolderSetDisorder2D : Nb must be 2*nsite!\n");
        exit(-1);
    }
    else if((placeholder->lconf->shape[0]%4)!=0 && (placeholder->lconf->shape[0]%4)!=0){
        printf("SEPlaceHolderSetDisorder2D : lx and ly must be multiple of 4!\n");
        exit(-1);
    }

    for(bond=0;bond<nsite;++bond){
        placeholder->lconf->J->data[bond]=1;
        i = (bond%placeholder->lconf->shape[0])%4;
        j = (bond/placeholder->lconf->shape[0])%4;
        if(gsl_rng_uniform_pos(placeholder->rng)<p) J=Jc*(1+dJ);
        else J=Jc*(1-dJ);
        if(i==1 && j==0) placeholder->lconf->J->data[bond]=J;
        else if(i==2 && j==1) placeholder->lconf->J->data[bond]=J;
        else if(i==3 && j==2) placeholder->lconf->J->data[bond]=J;
        else if(i==0 && j==3) placeholder->lconf->J->data[bond]=J;
    }
    for(bond=nsite;bond<(2*nsite);++bond){
        placeholder->lconf->J->data[bond]=1;
        i = ((bond-nsite)%placeholder->lconf->shape[0])%4;
        j = ((bond-nsite)/placeholder->lconf->shape[0])%4;
        if(gsl_rng_uniform_pos(placeholder->rng)<p) J=Jc*(1+dJ);
        else J=Jc*(1-dJ);
        if(i==0 && j==0) placeholder->lconf->J->data[bond]=J;
        else if(i==1 && j==1) placeholder->lconf->J->data[bond]=J;
        else if(i==2 && j==2) placeholder->lconf->J->data[bond]=J;
        else if(i==3 && j==3) placeholder->lconf->J->data[bond]=J;
    }
}

void SEPlaceHolderSetPlaquetteRandom2D(SEPlaceHolder* placeholder, double Jc, double dJ, double p)
{
    int Nb=placeholder->lconf->Nb;
    int nsite=placeholder->lconf->nsite;
    int bond,i,j;
    double J;
    if(placeholder->lconf->dims!=2){
        printf("SEPlaceHolderSetDisorder2D : The dimension must be 2!\n");
        exit(-1);
    }
    else if(Nb!=(2*nsite)){
        printf("SEPlaceHolderSetDisorder2D : Nb must be 2*nsite!\n");
        exit(-1);
    }
    else if((placeholder->lconf->shape[0]%2)!=0 && (placeholder->lconf->shape[0]%2)!=0){
        printf("SEPlaceHolderSetDisorder2D : lx and ly must be multiple of 2!\n");
        exit(-1);
    }

    for(bond=0;bond<nsite;++bond){
        placeholder->lconf->J->data[bond]=1;
        i = (bond%placeholder->lconf->shape[0])%4;
        j = (bond/placeholder->lconf->shape[0])%4;
        if(gsl_rng_uniform_pos(placeholder->rng)<p) J=Jc*(1+dJ);
        else J=Jc*(1-dJ);
        if(i==0) placeholder->lconf->J->data[bond]=J;
        else if(i==2) placeholder->lconf->J->data[bond]=J;
    }
    for(bond=nsite;bond<(2*nsite);++bond){
        placeholder->lconf->J->data[bond]=1;
        i = ((bond-nsite)%placeholder->lconf->shape[0])%4;
        j = ((bond-nsite)/placeholder->lconf->shape[0])%4;
        if(gsl_rng_uniform_pos(placeholder->rng)<p) J=Jc*(1+dJ);
        else J=Jc*(1-dJ);
        if(j==0) placeholder->lconf->J->data[bond]=J;
        else if(j==2) placeholder->lconf->J->data[bond]=J;
    }
}

void SEPlaceHolderSetConfigurationalDisorder2D(SEPlaceHolder* placeholder, double J)
{
    int Nb=placeholder->lconf->Nb;
    int nsite=placeholder->lconf->nsite;
    int bond,i,j;
    if(placeholder->lconf->dims!=2){
        printf("SEPlaceHolderSetDisorder2D : The dimension must be 2!\n");
        exit(-1);
    }
    else if(Nb!=(2*nsite)){
        printf("SEPlaceHolderSetDisorder2D : Nb must be 2*nsite!\n");
        exit(-1);
    }

    if(gsl_rng_uniform_pos(placeholder->rng)<0.5){
        for(bond=0;bond<nsite;++bond){
            placeholder->lconf->J->data[bond]=1;
            i = (bond%placeholder->lconf->shape[0])%2;
            j = (bond/placeholder->lconf->shape[0])%2;
            if(i==0) placeholder->lconf->J->data[bond]=J;
        }
    }
    else{
        for(bond=nsite;bond<(2*nsite);++bond){
            placeholder->lconf->J->data[bond]=1;
            i = (bond%placeholder->lconf->shape[0])%2;
            j = (bond/placeholder->lconf->shape[0])%2;
            if(j==0) placeholder->lconf->J->data[bond]=J;
        }
    }
}

#if 0
int main()
{
    int ndiff=2,length=12;
    int dims=1,shape[1]={8};
    int nsweep=1000000,seed=2133124;
    int cutoff=2000;
    double beta=1024;
    double max_err=1.e-4;

    SEPlaceHolder* placeholder = CreateSEPlaceHolder();
    SEPlaceHolderSetLattice(placeholder,mapping_1d,shape,dims,0);
    SEPlaceHolderSetLength(placeholder,length,ndiff);
    SEPlaceHolderSetRandomSeed(placeholder, seed);
    SEPlaceHolderSetNsweep(placeholder, nsweep, cutoff);
    SEPlaceHolderSetBeta(placeholder, beta);
    SEPlaceHolderSetError(placeholder, max_err);
    SEPlaceHolderCheckSetting(placeholder);

    placeholder->ops->sequence->data[ 0]=14;
    placeholder->ops->sequence->data[ 1]= 9;
    placeholder->ops->sequence->data[ 2]=-1;
    placeholder->ops->sequence->data[ 3]=13;
    placeholder->ops->sequence->data[ 4]= 4;
    placeholder->ops->sequence->data[ 5]=-1;
    placeholder->ops->sequence->data[ 6]=-1;
    placeholder->ops->sequence->data[ 7]= 6;
    placeholder->ops->sequence->data[ 8]=13;
    placeholder->ops->sequence->data[ 9]= 9;
    placeholder->ops->sequence->data[10]=-1;
    placeholder->ops->sequence->data[11]= 4;

    for(int i=0;i<length;++i) printf("%d ",placeholder->ops->sequence->data[i]);
    printf("\n");

    OperatorSequenceSort(placeholder->ops);
    SEOps2Lvc(placeholder->opl,placeholder->lconf,placeholder->ops);
    for(int i=0;i<length;++i){
        int l1 = placeholder->opl->lvc->data[4*i+0];
        int l2 = placeholder->opl->lvc->data[4*i+1];
        int l3 = placeholder->opl->lvc->data[4*i+2];
        int l4 = placeholder->opl->lvc->data[4*i+3];
        printf("%d %d %d %d\n",l1,l2,l3,l4);
    }

    int j=0;
    placeholder->opl->check=-1;
    while(placeholder->opl->check){
        SETraverseLoop(placeholder->opl,placeholder->lconf,placeholder->ops);
        for(int i=0;i<placeholder->opl->noo;++i){
            printf("%d ",placeholder->opl->loop->data[i]);
        }
        printf("\n");
        for(int i=0;i<placeholder->lconf->nflip;++i){
            printf("%d ",placeholder->lconf->flip->data[i]);
        }
        printf("\n");
        printf("----------------------\n");
    
        if(j==2) SELoopUpdate(placeholder->lconf,placeholder->ops,placeholder->opl);
        j++;
    }
    for(int i=0;i<length;++i) printf("%d ",placeholder->ops->sequence->data[i]);
    printf("\n");
}
#endif

#if 0
//Test memory leakage
int main()
{
    int i,ndiff=2,length=10000;
    int dims=2,shape[2]={32,64};

    for(i=0;i<1000000;i++){
        LatticeConf* lconf = CreateSquareLatticeConf(shape, dims);
        OperatorSequence* ops = CreateOperatorSequence(length,ndiff);
        OperatorLoop* opl = CreateOperatorLoop(length);

        DestroyLatticeConf(lconf);
        DestroyOperatorSequence(ops);
        DestroyOperatorLoop(opl);

#if 1
        SEPlaceHolder* placeholder = CreateSEPlaceHolder();
        SEPlaceHolderSetLattice(placeholder,shape,dims,0);
        SEPlaceHolderSetLength(placeholder,length,ndiff);
        DestroySEPlaceHolder(placeholder);
#endif
    }
}
#endif
