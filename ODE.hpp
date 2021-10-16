#ifndef _PSO_HPP_
#define _PSO_HPP_

/* Global PSO Parameters */
#define N_SPECIES 3
#define N_RATE_CONSTANTS 5
#define SAMPLE_SIZE_Y0 5000 // sample size of observed t = 0
#define SAMPLE_SIZE_X0 5000 // sample size of observed at time t
#define N_STEPS_1 50 // Blind PSO step size
#define N_STEPS_2 1000 // targeted PSO step size


/* Define which moments used in PSO estimation */
#define USE_VARIANCES 1  // True = 1, False = 0
#define USE_COVARIANCES 1

#endif