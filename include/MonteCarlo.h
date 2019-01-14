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


#endif
