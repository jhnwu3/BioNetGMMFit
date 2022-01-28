# **CyGMM**
CyGMM is a C++ software designed for rate constant estimation of CYTOF Snapshot Data. 
It takes into account both linear and nonlinear models of evolution for estimating rate constants.

## **Important Note: Operating System**
The program has only been compiled and tested on debian Linux based systems, specifically Ubuntu, using G++.
# *Getting Started*

## *Quickstart*
To quickly get started with one of the simulated examples, do:

1. In your terminal, pick a suitable directory for your liking and input
    
        git clone https://github.com/jhnwu3/CyGMM.git

    to quickly download the program to a specified directory. 


2. Install Eigen and Boost C++ libraries by running the install shell script by typing in the terminal

        source install.sh

    Make sure you run it with admin access.

3. Then, navigate to the previously created program folder

        cd Sig_GMM

4. By default, parameters for the 3 protein linear system are provided and simulated with a pre-defined evolution matrix defined in system.cpp in the main directory, hence to get started, simply run the run shell script to begin:

        source run.sh

5. All output (final estimate of rate constants) is recorded in out.txt

## *Prerequisites for Compiling* ##

### *Eigen*
Snapshot uses the Eigen 3.3.9 C++ linear algebra library for matrix computations. If on Ubuntu, you can do a quick install using:

    sudo apt install libeigen3-dev

otherwise you can see more detailed install instructions [here](https://eigen.tuxfamily.org/dox/GettingStarted.html)

### *Boost*
Snapshot uses the Boost 1.7.2 odeint C++ library for ODE estimations for nonlinear systems. To install the whole boost C++ library, you can try:

    sudo apt-get install libboost-all-dev

However, Snapshot only uses the C++ odeint library, so if storage space is an explicit concern, more
detailed install intructions can be found [here](https://www.boost.org/doc/libs/1_77_0/more/getting_started/unix-variants.html)

## **Compilation** ##

If you wish to modify the code for your own use or if the binary is not sufficient, a Makefile has been provided in the /src directory. 
After entering the directory

    cd src

Run

    make

in order to recompile an executable.

## **Running the executable** ##

To run the program, simply enter

    source run.sh

in your terminal. For more information about parameters and writing your own system, look below.

### *PSO Aside*
Although currently not available for estimating nonlinear systems, there is an optional two step procedure for the linear system, which may improve estimates. If run time is a concern, one can simply turn off the second step "targeted PSO" by simply setting the number of steps
of the targeted to 0. 

### *Data Inputs*
All data inputs are taken from the Data directory. By default, a set of randomly generated data points have been provided for the 3 species linear case for both X_0 and Y_0. For more run example data, look into the folder titled

    example

provided in this repo.

All data must be loaded in csv format. 

## *Loading in data*
If you want to load in X data file, make sure to delete or move any pre-existing csv file located in the 

    data/X 

directory and move or copy in your own X.csv file into the directory. Similarly, make sure to move
all Y_0 or Y_t data files into the directory listed as

    data/Y 

after moving or removing any previous Yt/Y0 files. Keep in mind, that the number of species of proteins you wish to simulate will correspond to the number of columns in each input X/Y file.

### *Rate Constant Inputs*
If you decide to simulate the rate constants and therefore simulate Y_t instead of manually inputting Yt files, make sure to define your set of rate constants in the "true_rates.csv" file. For instance, by default

    0.27678200
    0.83708059
    0.44321700
    0.04244124
    0.30464502

is defined in the file, which defines the true set of rate constants as "0.27678, 0.837, 0.44, 0.04, 0.30".

### *Time Inputs*
Make sure to list your the times for time evolutions in time_steps.csv rowwise. For a single time evolution, only a single time point, the time that your second comparison sample was evolved to, is needed in the file.

However, especially in the nonlinear case where multiple time points and samples may be beneficial, simply list out the times that each sample was evolved to, as shown below.

    0.5
    2
    10
    20
    30

## *Important Caveat*
One key thing to understand is every file in either the data/X or data/Y folders are read in alphabetical order. An error message and exit will output if the number of time steps do not match the number of Yt files. Make sure to label each file name in order of the time steps for proper loading.

Furthermore, if one chooses to simulate Y_t instead of inputing their own, keep in mind, it will specifically only choose the first Y file listed in the directory for use as Y_0.

### *Configuration Inputs*
To set the parameters you want for your estimation run, double click or open the

    Config.csv

file and change the respective values. For instance, through excel 

| n_particles1 | 15 |
|--------------|----|

or using a default text editor

    n_particles1,15

sets the number of particles in blind pso to 15.

The default PSO parameters are listed below,

| Parameter                        | Value |
|----------------------------------|-------|
| Number of Particles Blind PSO    | 1000  |
| Number of Steps Blind PSO        | 10    |
| Number of Particles Targeted PSO | 10    |
| Number of Steps Targeted PSO     | 1000  |
| Exclude Mixed Moments?           | 0     |
| Exclude Mixed and Second Moments?| 0     |
| Use Linear Model?                | 1     |
| Number of Runs                   | 1     |
| Simulate Y_t?                    | 1     |
| Use Matrix Inverse?              | 0     |
| Number of Rates                  | 5     |
| Sample Size                      | 5000  |
By default, the PSO runs with all moments, with means, variances, and covariances. Currently, there are only two other options for specifying which estimators to use. For instance, set

    Exclude Mixed Moments?,1

to use means and second moments only. All boolean options such as "Use Linear Model?" are set to on with 1, and set to off with 0.

### *System Parameters*

All ODE system parameters such as the number of protein species and rate constants are listed in system_parameters.csv. By default, the program runs with 3 protein species and 5 rate constants as well as respective X data sizes. *(Tentative feature to be added in)*, The program will autodetect the number of protein species based on the number of columns in the data csv files. 

(Tentative feature to be added in), will autodetect the total number of rows in the X/Y csv files. For now, a *sample size* must be defined in the Config.csv file.

## **Directory Structure** ##

### *sig_gmm*
The default directory contains all configuration csv files needed to run the program. 

### *src*
Contains all C++ files needed to recompile the program.

### *data*
This is where all the protein data is. The simulated data is provided as default for those wanting to use it based on the paper.




