# **CyGMM**
CyGMM is a C++ software designed for parameter estimation of CYTOF Snapshot Data. 
It takes into account both linear and nonlinear models of evolution for estimating parameters. 
# Table of Contents
1. [Quickstart Guide](#qstrt)
    1. [Dockers](#docker)
2. [Prerequisites](#paragraph1)
    1. [Eigen](#eig)
    2. [Boost](#bst)
    3. [libRoadRunner/BionetGen](#rr)
3. [Compilation](#compilation)
4. [Execution](#exe)
5. [Program Inputs](#pin)
    1. [Loading in Data](#indat)
    2. [Simulated Rate Constants](#rcons)
    3. [Time Inputs](#tim)
    4. [Program Configuration](#config)
6. [Defining Your Own System](#sys)
    1. [Interaction Matrices](#intmat)
    2. [Linear System](#linsys)
    3. [Nonlinear System](#nonlinsys)
## **Important Note: Operating System**
The program has only been compiled and tested on debian Linux based systems, specifically the latest version of Ubuntu.

Keep in mind, should you still want to pursue this on Windows 10, there is a Windows Linux Subsystem that runs at a similar level of performance that you can enable in the operating system. 
A Youtube video I've found helpful to get Linux up and running on Windows can be found [here](https://www.youtube.com/watch?v=A0eqZujVfYU). This code has been tested successfully on Windows Linux Subsystem. This is arguably a simpler solution than going through the trouble of CMake to get this C++ code/program up and running on Windows.

Mac's Unix based system should feasibly still work well with this repo, but has been untested in this current iteration.

## Quickstart  <a name="qstrt"></a>
To quickly get started with one of the simulated examples, do:

1. In your terminal, pick a suitable directory for your liking and input
    
        git clone https://github.com/jhnwu3/CyGMM.git

    to quickly download the program to a specified directory or if you don't have git, simply click on the green code button and left click download. Images below show where to click to download the zipfile.
    ![Step 1](/img/REPO.png)
    ![Step 2](/img/DLST.png)

    make sure to unzip the directory before use. Make sure you also have bionetgen installed, if not, look [here](https://bng-vscode-extension.readthedocs.io/en/latest/install.html).

2. Then, make sure you're still in the newly created repo directory

        cd CyGMM

3. By default, parameters for the 3 protein linear system are provided and simulated with a pre-defined evolution matrix defined in system.cpp in the main directory, hence to get started, simply run the run shell script to begin:

        source run.sh

4. All output (final estimate of rate constants) is recorded in **out.txt**

### Docker  <a name="docker"></a>
Should the quickstart statically compiled executable in the repository fail to run, there is a docker image that can be easily pulled and run on any operating system. 

1. Install docker or dockerhub [here](https://docs.docker.com/get-docker/) (this should be as easy as an executable)
2. Pull image from dockerhub, in your command line (bash, terminal, etc.), input 

    docker pull jhnwu3/cygmm:cygmm_run

3. Run docker image which currently should already have some pre-existing (birth-death) model inside.

    docker run -t cygmm_run

4.  Mounting Volumes to Write/Feed in Own Configuration Files (in the process of writing a python script to run this process with docker), do:

    docker run -v full_path_local_directory:cygmm/home -t cygmm_run

## Prequisites to Compiling <a name="prq"></a>

### *Eigen* <a name="eig"></a>
Snapshot uses the Eigen 3.3.9 C++ linear algebra library for matrix computations on Linux. If on Ubuntu, you can do a quick install using:

    sudo apt install libeigen3-dev

otherwise you can see more detailed install instructions [here](https://eigen.tuxfamily.org/dox/GettingStarted.html)

### *Boost* <a name="bst"></a>
Snapshot uses the Boost 1.7.2 odeint C++ library for ODE estimations for nonlinear systems. To install the whole boost C++ library, you can try:

    sudo apt-get install libboost-all-dev

However, Snapshot only uses the C++ odeint library, so if storage space is an explicit concern, more
detailed install intructions can be found [here](https://www.boost.org/doc/libs/1_77_0/more/getting_started/unix-variants.html)

### *libRoadRunner/bionetgen* <a name="rr"> </a>
The program uses a python library called bionetgen. Please make sure you have python installed from instructions [here](https://www.python.org/downloads/).
Please make sure to have bionetgen(https://bionetgen.org/) installed through 

    pip install bionetgen

and should you want to compile your own code from libroadrunner off of this, please look [here](https://libroadrunner.readthedocs.io/en/latest/Installation/installation.html)

## Compilation <a name="compilation"></a>

If you wish to modify the code for your own use or if the binary is not sufficient, a cmake has been provided in the /src directory. Fair warning this can be a tedious and bug-prone
process, that being said, assuming you have installed all of **Boost** and **Eigen** libraries through the above, then you can simply just download a fully built roadrunner + CyGMM library
[here](https://drive.google.com/drive/folders/1n6S7y2sf88mb_62evYO-RimbdT9FHMwh?usp=sharing). 

First unzip the folder, doing

    unzip 

After entering the directory /src

    cd src

Run
    cmake .
    make

in order to recompile an executable.

## Execution <a name="exe"></a>

To run the program, simply enter

    source run.sh

in your terminal. For more information about parameters and writing your own system, look below.

### *PSO Aside*
Currently, targeted PSO has only been provided should a user choose to use matrix exponentiation. 

## Program Inputs <a name="pin"></a>
All data inputs are taken from the Data directory. By default, a set of randomly generated data points have been provided for the 3 species linear case for both X_0 and Y_0. For more run example data, look into the folder titled

    example

provided in this repo.

All data must be loaded in csv format. 

### *Loading in data* <a name="indat"></a>
If you want to load in the X, control, or time 0 data file, make sure to delete or move any pre-existing csv file located in the 

    data/X 

directory and move or copy in your own X.csv file into the directory. Similarly, make sure to move
all Y_0 or Y_t data files into the directory listed as

    data/Y 

after moving or removing any previous Yt/Y0 files. Keep in mind, that the number of species of proteins you wish to simulate will correspond to the number of columns in each input X/Y file and the number of rows correspond to the cell count. 

For more examples please take a look at the /example directory. The real data X and Y are labeled in the 4_prot_real folder. All real data was taken from 
[here](https://dpeerlab.github.io/dpeerlab-website/dremi-data.html), specifically the first CD8, CD28, and CD3 naive time series (third column, 1st row).

### *Rate Constant Inputs* <a name="rcons"></a>
If you decide to simulate the rate constants and therefore simulate Y_t instead of manually inputting Yt files, make sure to define your set of rate constants in the "true_rates.csv" file. For instance, by default

    0.27678200
    0.83708059
    0.44321700
    0.04244124
    0.30464502

is defined in the file, which defines the true set of rate constants as "0.27678, 0.837, 0.44, 0.04, 0.30".

### *Time Inputs* <a name="tim"></a>
Make sure to list your the times for time evolutions in time_steps.csv rowwise. For a single time evolution, only two time points, the end and start time of your evolution interval, is needed in the file.

However, especially in the nonlinear case where multiple time points and samples may be beneficial, simply list out each of the times evolved for rowwise, as shown below.
    0
    0.5
    2
    10
    20
    30

### *Important Caveat* 
One key thing to understand is every file in either the data/X or data/Y folders are read in alphabetical order. An error message and exit will output if the number of time steps do not match the number of Yt files. Make sure to label each file name in order of the time steps for proper loading.

Furthermore, if one chooses to simulate Y_t instead of inputing their own, keep in mind, it will specifically only choose the first Y file listed in the directory for use as Y_0.

### *Configuration Inputs* <a name="config"></a>
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
| Use BNGL?                        | 1     | 1 to use bionetgen language for simulation, 0 otherwise                  |
| Use Deterministic?               | 1     | 1 to use CVode integrators, 0 to use roadrunner gillespie simulation     |
| Number of BNGL Steps             | 15    | Tuning Parameter for number of steps of integration                      |
| Seed                             | -1    | Used to seed the PSO, Off when seed < 0, On when seed > 0                |
| Parallel Number of Threads       | 8     | Number of threads to parallelize on.                                     |
By default, the PSO runs with all moments, with means, variances, and covariances. Currently, there are only two other options for specifying which estimators to use. For instance, set

    Exclude Mixed Moments?,1

to use means and second moments only while

    Exclude Mixed and Second Moments?, 1

will force the program to estimate rate constants using means only. All boolean options such as "Use Linear Model?" are set to on with 1, and set to off with 0. 

Finally, regarding holding parameter or rate constant values, these are currently only enabled for the nonlinear system where it's necessary for accurate estimation. 

## Defining Your Own Linear or Nonlinear System <a name="sys"></a>
There are two ways to define your own biological system, one is through bionetgen, by compiling bionetgen into sbml code shown in model.bngl, CyGMM is able to perform parameter estimation using bngl. For more information on bionetgen, see [here](https://webng.readthedocs.io/en/latest/). The most important caveat is to make sure the bngl file writes

    writeSBML()

at the end. 

However, should you want to write explicit odes, look below. Although this method provides substantial performance gains, it also requires you to recompile CyGMM.
 
### *Interaction Matrices* <a name="intmat"></a>
In order to minimize computation times, matrix exponentiation is commonly used to quickly solve a system of coupled linear odes, specifically in the form of 

![generalized linear system](/img/matExpSys.png)

where M = some Matrix where the elements of the M matrix are given by zero or linear sums of
the reaction rates k. P denotes the protein abundances, i their respective protein types, and observe that we can now solve it using matrix exponentiation in

![generalized linear solution](/img/matExpSoln.png)

To give an example, consider the protein system shown below,

![linear 3 example protein system for mat exp](/img/lin3ExpModel.png)

observe that we can model this using a system of linear differential equations,

![interaction diffeq system](/img/matExpSysEx.png)

and then observe that the coefficients in this differential equation system can be easily mapped to an interaction matrix M as shown below

![interaction Matrix](/img/matExpEx.png)

to take this into code, we can now bake this into our program. Navigate to the 

    system.cpp

C++ file and open it using your preferred text editor. You'll notice that there's a single function called interaction matrix and it should have a variable called intMatrix with numerous comments surrounding it. It should look something like this: 

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

    -k1 - k4   or -theta1 - theta4 (if using theta notation)

as an aside:
    - the + sign is not needed to denote positive values and is only needed for arithmetic operations i.e -k(0) + -k(3)
    - math operators are denoted as such, multiplication(*), division(/), addition(+), subtraction(-)
    - whitespace is not a concern
    - the compiler will take care of negative values such as -k(0) 

### *Linear System* <a name="linsys"></a>
Now suppose, your linear system cannot be modeled using an interaction matrix. A good example of a protein system that might have this characteristic is shown below,

![4 Protein Linear System](/img/linear4Sys.png)

with its respective linear system

![4 Protein Linear DiffEq](/img/linear4DiffEq.png)

Worry not, we can still solve this system using the conventional method of runge kutta solvers. To do this, simply navigate to the 

    /src/system.hpp 

header file (by first double clicking the src folder, and then double clicking again the system.hpp file where you can say open with notepad/etc.), and again open it with your preferred text editor (i.e notepad). Now observe it's direct mapping in code,

![linear4](/img/genSystem.png)

to get a further look, let us take a deeper look at one equation in code and how it maps directly to another equation in mathematic terms. Consider the equation,

    dcdt[0] = k(0) - k(4) * c[0];

and how it directly maps to the first equation in the linear system

    dP1/dt = k1 - k5P1

in this case, first observe that k(0) maps to k1, and P1 maps to c[0]. This is a key point, because there is two extremely important pieces of this syntax.

1. The code indexes from 0 not 1, hence k1 maps to k(0), k2 maps to k(1), etc.
2. The proteins indices are referenced through the use of [] brackets instead of () paranthesises while rate constants are referenced through () not [], hence why you see that it is c[0] not c(0) and k(0) not k[0]

Now also realize that the derivative dP1/dt is mapped to the variable dcdt, hence all specific dPi/dt are mapped by dcdt[i]. Once you have finished writing the equation at hand, make sure to end each equation with a semicolon ";" as shown at the end of the equation above.

###### *All arithmetic operators are the same as was defined in the Interaction Matrices Section.*

### *Nonlinear System* <a name="nonlinsys"></a>
Defining a nonlinear system can be done through the same process as the linear system, first consider some nonlinear protein system,

![Nonlinear Proteins](/img/nonlinearSysChemEqP.png)

This can be broken down into a system of differential equations as such

![Nonlinear Proteins](/img/nonlinearSysDiffEq.png)

Now, like the linear system, we can map this to code, to get started, enter the src directory (you can also just double click the folder) 

    cd src

Open the the *system.hpp* file in any text or code editor, and you'll see code lined out as shown below:

![Nonlinear](/img/nonlinearCode.png)

Observe that it resembles a differential system of equations. Defining a system is done within the circled region of text. You can simply ignore everything else surrounding the system of equations. 

Observe that the first element is listed as

    dcdt[0] = -(k(0) * c[0] * c[1]) + k(1) * c[2] + k(2) * c[2];

which, corresponds to the equation 

    dc1/dt = -k1 *c1*c2 + k2*c3 + k3*c3 

note the syntax where each equation ends with a semicolon ";" and observe that rate constant k1 is denoted by k(0) where the indexing starts at 0 again. 

Now, observe that all nonconstant terms, the derivative term dcdt, and its respective element c, are all indexed using [ ] instead of parenthesises. Again, indexed starting 0, we have dc1/dt and c1 respectively listed as

    dcdt[0], c[0]

In terms of mathematical operators, only basic, +,-,*,/ operators have been tested with the boost odeint library, but in theory math.exp(a) and other standard math libraries should work.

## **Directory Structure** ## <a name="dir"></a>

### *CyGMM*
Main Directory with general configuration files and system.cpp code files.

### *src*
Contains all C++ source and header files needed to recompile the program.

### *data*
This is where protein input data is read in from the program. Example data (used in the paper) is provided in the example folder.

### *example*
Contains various examples for use with the linear and nonlinear system provided by default in the program / code. There should be codes associated with them.





