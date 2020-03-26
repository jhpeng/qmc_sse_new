#ifndef MONTECARLO_H
#define MONTECARLO_H

void MCGeneralSchemeAndLattice(
            int* shape, 
            int mode, 
            int lattice, 
            double J, 
            double dJ, 
            double p, 
            double beta, 
            double beta_i, 
            double beta_f, 
            double interv, 
            int thermal, 
            int nsweep, 
            int nblock, 
            int ntime, 
            int seed);
#endif
