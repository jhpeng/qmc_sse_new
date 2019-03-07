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
            int nblock,
            int cutoff, 
            int seed);

void MCHerringbond2DImproveSpeed(
            double J, 
            double beta, 
            int* shape, 
            int nsweep, 
            int nblock,
            int cutoff, 
            int seed);

void MCBetaIncreasePlaquetteDisorder2D(
            double J,
            double dJ,
            double p,
            double beta_i, 
            double beta_f, 
            double interval, 
            int* shape, 
            int nsweep, 
            int cutoff, 
            int seed);

void MCBetaIncreasePlaquetteDisorderImproveSpeed2D(
            double J,
            double dJ,
            double p,
            double beta_i, 
            double beta_f, 
            double interval, 
            int* shape, 
            int nsweep, 
            int cutoff, 
            int seed);

void MCBetaIncreaseConfigurationalDisorder2D(
            double J, 
            double beta_i, 
            double beta_f, 
            double interval, 
            int* shape, 
            int nsweep, 
            int cutoff, 
            int seed);

void MCZeroTempPlaquetteDisorder2D(
            double J, 
            double dJ, 
            double p, 
            double beta, 
            int* shape, 
            int nsweep, 
            int cutoff, 
            int seed);

#endif
