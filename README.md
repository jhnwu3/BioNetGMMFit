# **Sig_GMM**
Sig_GMM is a C++ software designed for rate constant estimation of CYTOF Snapshot Data. 
It takes into account both linear and nonlinear models of evolution for estimating rate constants.

## **Important Note: Operating System**
The program has only been compiled and tested on Linux based systems, using G++.

## **(Optional) Prerequisites for Compiling ** ##

### *Eigen*
Snapshot uses the Eigen 3.3.9 C++ linear algebra library for matrix computations. If on Ubuntu, you can do a quick install using:

    sudo apt install libeigen3-dev

otherwise you can see more detailed install instructions [here](https://eigen.tuxfamily.org/dox/GettingStarted.html)

### *Boost*
Snapshot uses the Boost 1.7.2 odeint C++ library for ODE estimations for nonlinear systems. To install the whole boost C++ library, you can try:

    sudo apt-get install libboost-all-dev

However, Snapshot only uses the C++ odeint library, so if storage space is an explicit concern, more
detailed install intructions can be found [here](https://www.boost.org/doc/libs/1_77_0/more/getting_started/unix-variants.html)


## **(Optional) Compilation** ##

If you wish to modify the code for your own use or if the binary is not sufficient, a Makefile has been provided in the /src directory. 
After entering the directory

    cd src

Run

    make

in order to recompile an executable.

## **Running the executable** ##

### *PSO Aside*

Although currently not available for estimating nonlinear systems, there is an optional two step procedure for the linear system, which
may improve estimates. If run time is a concern, one can simply turn off the second step "targeted PSO" by simply setting the number of steps
of the targeted to 0. 

### *Time Inputs*
Make sure to list your the times for time evolutions in time_steps.csv rowwise. For a single time evolution, only a single time point, the time that your second comparison sample was evolved to, is needed in the file.

However, especially in the nonlinear case where multiple time points and samples may be beneficial, simply list out the times that each sample was evolved to, as shown below.

    0.5
    2
    10
    20
    30

### *Data Inputs*
All data must be loaded in csv format. Make sure to load your "base" or initial sample values into
X_0.csv, which will be evolved using the differential equation model selected.

Make sure to specifiy all respective Yt and "evolved" sample file names in the Ytnames.csv file.

### *PSO Inputs*
To set the parameters you want for your estimation run, open up PSO.csv
and change the respective values. For instance, through excel 

| n_particles1 | 15 |
|--------------|----|

or using a default text editor

    n_particles1,15

sets the number of particles in blind pso to 15.

The default PSO parameters are listed below,

| Parameter                        | Value |
|----------------------------------|-------|
| Number of Particles Blind PSO    | 15    |
| Number of Steps Blind PSO        | 50    |
| Number of Particles Targeted PSO | 1     |
| Number of Steps Targeted PSO     | 5     |
| Using Only Second Moments?       | 0     |
| Using Only First Moments?        | 0     |
| Use Linear Model?                | 1     |
| Number of Runs                   | 1     |

By default, the PSO runs with all moments, with means, variances, and covariances. Currently, there are only two other options for specifying which estimators to use. For instance, set

    use_OnlySecMoments?,1

to use means + variances only.

## **Directory Structure** ##

## *sig_gmm*
The default directory contains all configuration csv files needed to run the program. 

### *src*
Contains all C++ files needed to recompile the program.

### *input*
This is where all the protein data is. The simulated data is provided as default for those wanting to use it based on the paper.




