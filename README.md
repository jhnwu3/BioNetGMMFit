# **CyGMM**
CyGMM is a C++ software designed for parameter estimation of CYTOF Snapshot Data. 
It takes into account both linear and nonlinear models of evolution for estimating parameters. 

## **Important Note: Operating System**
The program has only been compiled and tested on debian Linux based systems, specifically the latest version of Ubuntu.

Keep in mind, should you still want to pursue this on Windows 10, there is a Windows Linux Subsystem that runs at a similar level of performance that you can enable in the operating system. 
A Youtube video I've found helpful to get Linux up and running on Windows can be found [here](https://www.youtube.com/watch?v=A0eqZujVfYU). This code has been tested successfully on Windows Linux Subsystem. This is arguably a simpler solution than going through the trouble of CMake or Chocolatey to get this C++ code/program up and running on Windows.

Mac's Unix based system should feasibly still work well with this repo, but has been untested in this current iteration.

# *Getting Started*

## *Quickstart*
To quickly get started with one of the simulated examples, do:

1. In your terminal, pick a suitable directory for your liking and input
    
        git clone https://github.com/jhnwu3/CyGMM.git

    to quickly download the program to a specified directory or if you don't have git, simply click on the green code button and left click download. Images below show where to click to download the zipfile.
    ![Step 1](/img/REPO.png)
    ![Step 2](/img/DLST.png)

    make sure to unzip the directory before use.

2. Open the default (or WSL) terminal in the directory and install Eigen and Boost C++ libraries by running the install shell script by typing in the terminal

        source install.sh

    Make sure you run it with admin access.

3. Then, make sure you're still in the newly created repo directory

        cd CyGMM

4. By default, parameters for the 3 protein linear system are provided and simulated with a pre-defined evolution matrix defined in system.cpp in the main directory, hence to get started, simply run the run shell script to begin:

        source run.sh

5. All output (final estimate of rate constants) is recorded in **out.txt**

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
If you want to load in the X, control, or time 0 data file, make sure to delete or move any pre-existing csv file located in the 

    data/X 

directory and move or copy in your own X.csv file into the directory. Similarly, make sure to move
all Y_0 or Y_t data files into the directory listed as

    data/Y 

after moving or removing any previous Yt/Y0 files. Keep in mind, that the number of species of proteins you wish to simulate will correspond to the number of columns in each input X/Y file and the number of rows correspond to the cell count. 

For more examples please take a look at the /example directory. The real data X and Y are labeled in the 4_prot_real folder. All real data was taken from 
[here](https://dpeerlab.github.io/dpeerlab-website/dremi-data.html), specifically the first CD8, CD28, and CD3 naive time series (third column, 1st row).

### *Rate Constant Inputs*
If you decide to simulate the rate constants and therefore simulate Y_t instead of manually inputting Yt files, make sure to define your set of rate constants in the "true_rates.csv" file. For instance, by default

    0.27678200
    0.83708059
    0.44321700
    0.04244124
    0.30464502

is defined in the file, which defines the true set of rate constants as "0.27678, 0.837, 0.44, 0.04, 0.30".

### *Time Inputs*
Make sure to list your the times for time evolutions in time_steps.csv rowwise. For a single time evolution, only two time points, the end and start time of your evolution interval, is needed in the file.

However, especially in the nonlinear case where multiple time points and samples may be beneficial, simply list out each of the times evolved for rowwise, as shown below.
    0
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

| Parameter                        | Value | Explanation                                                              |
|----------------------------------|-------|--------------------------------------------------------------------------|
| Number of Particles Blind PSO    | 1000  | Sets number of particles in first stage PSO                              |
| Number of Steps Blind PSO        | 10    | Sets number of steps for first stage PSO                                 |
| Number of Particles Targeted PSO | 10    | Sets number of particles in second stage PSO                             |
| Number of Steps Targeted PSO     | 1000  | Sets number of steps for second stage PSO                                |
| Exclude Mixed Moments?           | 0     | 1 to use only means and variances, 0 otherwise                           |
| Exclude Mixed and Second Moments?| 0     | 1 to use only means, 0 otherwise                                         |
| Use Matrix Exp Model?            | 1     | 1 to use a matrix system, 0 to use nonlinear                             |
| Number of Runs                   | 1     | Sets total number of PSO runs for estimation                             |
| Simulate Y_t?                    | 1     | 1 to simulate Yt with a true rate vector, 0 to provide own Yt matrix     |
| Use Matrix Inverse?              | 0     | 1 to use C++'s Matrix Inverse, 0 otherwise                               |
| Number of Rates                  | 5     | Sets number of rates in rate vector                                      |
| Index of Held Rate Constant      | -1    | -1 to not hold a rate constant, else specified theta i is held constant  |
| Value of Held Rate Constant      | 0     | Value between 0 and 1 that a rate constant would be held at              |
| Hypercube Dimension              | 1.0   | Real Value Dimensions of Hypercube to be searched in PSO.                |
| Report Moments?                  | 1     | 1 to report predicted moments in out.txt                                 |
| Number of Nested Hypcubes        | 1     | Estimates rate constants in spaces of 2^n order, n = # of nested hypcubes|
| Bootstrap?                       | 1     | 1 to estimate 95% CI's, 0 otherwise                                      |   

By default, the PSO runs with all moments, with means, variances, and covariances. Currently, there are only two other options for specifying which estimators to use. For instance, set

    Exclude Mixed Moments?,1

to use means and second moments only while

    Exclude Mixed and Second Moments?, 1

will force the program to estimate rate constants using means only. All boolean options such as "Use Linear Model?" are set to on with 1, and set to off with 0. 

Finally, regarding holding parameter or rate constant values, these are currently only enabled for the nonlinear system where it's necessary for accurate estimation. 

### *Defining Your Own Linear or Nonlinear System*
Defining your own linear and nonlinear system will currently require writing a little bit of code. Unfortunately, there isn't a GUI to use, however, with enough of an understanding of interaction matrices (and how they relate to your system of equations), the syntax is fairly straightforward.

#### *Interaction Matrices*
Navigate to the 

    system.cpp

C++ file and you'll notice that there's a single function called interaction matrix and it should have a variable called intMatrix with numerous comments surrounding it. It should look something like this: 

![Step 1](/img/intMat.png)

In order to define your own linear matrix system, observe that we are assigning a matrix with specific values in code written as:

    intMatrix << 
        -k(2), k(2), 0,
		k(1), -k(1) - k(4), k(4),
		k(3), k(0), -k(0) - k(3);

now the only part you should be concerned about changing is 

    -k(2), k(2), 0,
		k(1), -k(1) - k(4), k(4),
		k(3), k(0), -k(0) - k(3);

on the next line after the << assignment operator. In this current iteration (possibly changed later on), you must pre-define the number of rates you're estimating or simulating in Config.csv.
In the example case of this linear system, we have 5 rate constants defined, hence

    Number of Rates, 5

Now, observe from above, we have a 3x3 interaction matrix for the 3 species linear system. Each element in the matrix is separated by a comma. In this case, we have 3 elements per row, hence we have 

    -k(2), k(2), 0,

3 commas, each comma right after every matrix element. Now, also notice that there is a k defined. This *k* is your rate constant vector that you are estimating. Understand that in C++, indexing starts at 0 instead of 1, hence the

    k(2)

inside the row

    -k(2), k(2), 0,

is the *second* matrix element of the first row referring to a positive rate constant k3 or put simply the third element in the k rate vector. Similarly,

    k(0)

and is equivalent to

    k1

Now let's take a look at forming compound expressions in an interaction matrix element, this row

    k(3), k(0), -k(0) - k(3);

with the element

    -k(0) - k(3)

is equivalent to 

    -k1 - k4

as an aside:
    - the + sign is not needed to denote positive values and is only needed for arithmetic operations i.e -k(0) + -k(3)
    - whitespace is not a concern
    - the compiler will take care of negative values such as -k(0) 


#### *Default System*
Currently, the nonlinear system is not very well documented. However, should one choose to go further into the code and attempt to define their own nonlinear system. It is important to understand how to interact with the boost odeint library. Generally speaking, outside of understanding the nuances between ( ) and [ ] syntaxes as we are feeding in an Eigen vector into a boost ode-solver, the syntax should still be fairly close to writing out a system of equations mathematically.

To get started, enter the src directory (you can also just double click the folder) 

    cd src

Open the the *system.hpp* file in any text or code editor, and you'll see code lined out as shown below:

![Nonlinear](/img/genSystem.png)

Observe that it resembles a differential system of equations. Defining a system is done within the circled region of text. You can simply ignore everything else surrounding the system of equations. 

Observe that the first element is listed as

    dcdt[0] = -(k(0) * c[0] * c[1]) + k(1) * c[2] + k(2) * c[2];

which, corresponds to the equation 

    dc1/dt = -k1 *c1*c2 + k2*c3 + k3*c3 

note the syntax where each equation ends with a semicolon ";" and observe that rate constant k1 is denoted by k(0) where the indexing starts at 0 again. 

Now, observe that all nonconstant terms, the derivative term dcdt, and its respective element c, are all indexed using [ ] instead of parenthesises. Again, indexed starting 0, we have dc1/dt and c1 respectively listed as

    dcdt[0], c[0]

In terms of mathematical operators, only basic, +,-,*,/ operators have been tested with the boost odeint library, but in theory math.exp(a) and other standard math libraries should work.

## **Directory Structure** ##

### *CyGMM*
Main Directory with general configuration files and system.cpp code files.

### *src*
Contains all C++ source and header files needed to recompile the program.

### *data*
This is where protein input data is read in from the program. Example data (used in the paper) is provided in the example folder.

### *example*
Contains various examples for use with the linear and nonlinear system provided by default in the program / code. There is a system.pdf that gives a brief summary explanation of the nonlinear ODE system used.





