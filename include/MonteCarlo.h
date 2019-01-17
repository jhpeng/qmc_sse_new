#ifndef MONTECARLO_H
#define MONTECARLO_H

void MCDisorder2D(
            double J, 
            double beta, 
            int* shape, 
            int nsweep, 
            int cutoff, 
            int seed);

void MCHerringbond2D(
            double J, 
            double beta, 
            int* shape, 
            int nsweep, 
            int cutoff, 
            int seed);

void MCBetaIncrease2D(
            double beta_i, 
            double beta_f, 
            double interval, 
            int* shape, 
            int nsweep, 
            int cutoff, 
            int seed);

#endif
