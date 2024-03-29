# **BioNetGMMFit**
BioNetGMMFit (formerly CyGMM) is a C++ software designed for parameter estimation of BioNetGen models using Snapshot Data. 
It takes into account both deterministic as well as stochastic models of evolution for estimating parameters. 

[Methods Paper](https://www.biorxiv.org/content/10.1101/2022.03.17.484491v1)

[Software Paper](https://www.biorxiv.org/content/10.1101/2022.12.08.519526v1)

For a quick web demo to see what one can do, please see [here](https://bngmm.nchigm.org/)

For those who are unfamiliar with the intuition behind generalized method of moments (GMM), please look [here](https://github.com/jhnwu3/BioNetGMMFit/blob/main/example/gmm_tutorial.ipynb) for a quick python tutorial.

# Table of Contents
1. [Quickstart Guide](#qstrt)
    1. [Docker](#docker)
    2. [Static Binary (Recommended for Linux Users and Performance)](#statbin)
2. [Prerequisites](#paragraph1)
    1. [Eigen](#eig)
    2. [Boost](#bst)
    3. [libRoadRunner/BionetGen](#rr)
3. [Compilation](#compilation)
4. [Execution and Examples](#exe)
5. [Program Inputs](#pin)
    1. [Loading in Data](#indat)
    2. [Simulated Rate Constants](#rcons)
    3. [Time Inputs](#tim)
    4. [Program Configuration](#config)
6. [Defining Your Own System in BioNetGen](#bngl)

## **Important Note: Operating System**
The preferred beginner's method of using this software is through [Docker](https://docs.docker.com/get-started/) as the docker virtualization engine allows BioNetGMMFit to be compatible with all operating systems as well as reduces any operating system specific dependencies.

However should computational performance be a major concern, the non-docker program has only been compiled and tested on debian Linux based systems, specifically the latest version of Ubuntu.

Keep in mind, on Windows 10, there is a Windows Linux Subsystem that runs at a similar level of performance that you can enable in the operating system. A Youtube video I've found helpful to get Linux up and running on Windows can be found [here](https://www.youtube.com/watch?v=A0eqZujVfYU). This code has been tested successfully on Windows Linux Subsystem. This is arguably a simpler solution than going through the trouble of CMake to get this C++ code/program up and running on Windows. 

Mac's Unix based system should feasibly still work well with the BioNetGMMFit binary, but has been untested in this current iteration.

## Quickstart  <a name="qstrt"></a>

There are two ways to get BioNetGMMFit up and running, one if you're on Ubuntu and/or a Linux Based System and need maximum performance on an HPC, download the latest static binary from the static binary section, otherwise the preferred method is through docker.

To quickly get started with one of the simulated examples, do:

### Docker  <a name="docker"></a>
There is a docker image that can be easily pulled and run on any operating system from [here](https://hub.docker.com/r/jhnwu3/bngmm). 

1. Install docker [here](https://docs.docker.com/get-docker/) 

If you're on Windows, you will have to install WSL first, which will simply require you to open power shell as an administrator and inputting:

    wsl --install

For more information for Windows Linux Subsystem, see [here](https://learn.microsoft.com/en-us/windows/wsl/install). If you're not already familiar with a terminal/console, then please do see here [here](https://www.youtube.com/watch?v=A0eqZujVfYU) for how to setup a linux kernel on Windows. Setting up a linux based terminal will allow you to perform the following steps.

2. Pull image from dockerhub, in your command line (bash, linux terminal if on Windows, etc.), input 

    docker pull jhnwu3/bngmm:ui

3. Run docker image in terminal

    docker run -d -p 5000:5000  jhnwu3/bngmm:ui

    If you have docker desktop setup, you should be able to see the program's running logs as shown below.
    ![Step 3](/img/dockerdesktop.png)

4. Then, to see the web UI, simply go open your browser, and enter the URL as follows

    http://127.0.0.1:5000

    which should open up something like this below:
     ![Step 4](/img/WebUI.png)


5. Now, let's get started with a basic example, please download the folder from [here](https://github.com/jhnwu3/BioNetGMMFit-Example) or do

    git clone https://github.com/jhnwu3/BioNetGMMFit-Example.git

in one of your preferred directories, i.e a download folder, etc.

6. Now, to upload the necessary files, first click on the "Browse..." button, as shown below

     ![Step 6](/img/BrowseButton.png)

7. Then, a file explorer button should open up, this will change depending on which operating system you are on. If windows, it will look like below.
    
    ![Step 7](/img/FileUpload.png)

Simply highlight all of them by holding left click on all of the files and selecting them by clicking the open button.

![Step 7.1](/img/uploadButton.png)

Then just upload the files, by clicking the big blue "submit" button.

![Step 7.2](/img/Submit.png)

8. Once uploaded, you will see the uploaded files listed below.
![Step 8](/img/uploaded.png)

9. Now, match all uploaded files to each BNGMM output like below.
![Step 9](/img/UIConfig.png)

Please match:
    Config4pro.csv to the Configuration File.
    4proV2.bngl to BioNetGen
    time_steps4.csv to Time Steps
    t1m_processed.csv to X - Initial Abundances
    t2m_processed.csv to Y - Final Abundances

![Config UI](/img/FullyConfigured.png)

10. Click the big blue button Run CyGMM to run the program. 
![Run](/img/Run.png)

If you have docker desktop installed, its terminal view should show something like this.
![RunD](/img/dockerRun.png)

11. Wait approximately 20 - 30 seconds, and the browser page should update with something like this.

![Result](/img/uiResults.png)

### CMAKE Static Binary (Recommended for Linux/WSL and HPC Users) <a name="statbin"></a>

1. Get prerequisites for compiling.

        sudo apt-get install -y build-essential git cmake autoconf libtool pkg-config libncurses5-dev libeigen3-dev libboost-all-dev python3 python3-pip && \
        pip install bionetgen matplotlib flask

<!-- 0. Make sure to have git installed, see [here](https://github.com/git-guides/install-git) -->
2. In your bash terminal, pick a suitable directory for your liking and download [here](https://zenodo.org/record/7733865)

    <!-- to quickly download the program to a specified directory or if you don't have git, simply click on the green code button and left click download. Images below show where to click to download the zipfile. -->
    <!-- ![Step 1](/img/REPO.png)
    ![Step 2](/img/DLST.png) -->

    <!-- make sure to unzip the directory before use. Make sure you also have bionetgen installed, if not, look [here](https://bng-vscode-extension.readthedocs.io/en/latest/install.html). Or one can download a binary from Google Drive [here]() (might be missing). -->

    and make sure to unzip the zip file.

        unzip BNGMM_Build.zip

3. To run BNGMM (Ubuntu and WSL Ubuntu users only)

        cd BNGMM_Build/BNGMM_DockerBuild/BNGMM
        ./BNGMM -h 


4. If you have executable errors and most likely are on a separate linux distribution (i.e redhat), we will need to recompile the executable. So, first make sure you're in the build-release directory.

        cd /path/to/your/dir/BNGMM_build/BNGMM_DockerBuild/buildroadrunner/roadrunner/build-release

    where /path/to/your/dir is the directory that is housing your unzipped build folder. For instance, if you're in the downloads folder, it'd be something like 
        
        cd ~/Downloads/BNGMM_build/BNGMM_DockerBuild/buildroadrunner/roadrunner/build-release

5. cmake build to reset all build-targets to your local directory.

        cmake -DCMAKE_INSTALL_PREFIX="../install-Release" -DLLVM_INSTALL_PREFIX="full/path/to/dir/llvm13-ubuntu-gcc10-rel" -DRR_DEPENDENCIES_INSTALL_PREFIX="../../libroadrunner-deps/install-Release" -DCMAKE_BUILD_TYPE="Release" ..

        cmake --build . --target install --config Release

    For instance, if downloaded and unzipped in the Downloads directory, it might look like

        cmake -DCMAKE_INSTALL_PREFIX="../install-Release" -DLLVM_INSTALL_PREFIX="/home/user/Downloads/BNGMM_Build/BNGMM_DockerBuild/buildroadrunner/llvm13-ubuntu-gcc10-rel" -DRR_DEPENDENCIES_INSTALL_PREFIX="../../libroadrunner-deps/install-Release" -DCMAKE_BUILD_TYPE="Release" ..

        cmake --build . --target install --config Release 

    Now, one might run into errors in the process of compilation (i.e cannot find LLVM or something like that), which we will explain how to resolve below.

6. Unfortunately, the rebuilding of this distribution's libRoadRunner build package will most likely contain minor bugs with pathing, but fortunately, there's an easy fix. Please go to the new install-Release's cmake directory.

        cd ../install-Release/lib/cmake 

7. Now open the ImportRoadrunnerAndDependencies.cmake file with your preferred text editor.

        vi ImportRoadrunnerAndDependencies.cmake 

    where you will see something like the following below,

        find_package(Threads) # for libxml2, FindThreads.cmake is shipped with cmake
        find_package(LibLZMA) # for libxml2, LibLZMA.cmake is shipped with cmake
        find_package(zlib CONFIG REQUIRED)
        find_package(bzip2 CONFIG REQUIRED)
        find_package(iconv CONFIG REQUIRED)
        find_package(LibXml2 CONFIG REQUIRED)
        find_package(libsbml-static CONFIG REQUIRED)
        find_package(rr-libstruct CONFIG REQUIRED)
        find_package(clapack CONFIG REQUIRED)
        find_package(nleq1 CONFIG REQUIRED)
        find_package(nleq2 CONFIG REQUIRED)
        find_package(PocoFoundation CONFIG REQUIRED)
        find_package(PocoNet CONFIG REQUIRED)
        find_package(PocoXML CONFIG REQUIRED)
        find_package(Sundials CONFIG REQUIRED)
        find_package(LLVM REQUIRED)
        find_package(roadrunner-static CONFIG REQUIRED)
        find_package(roadrunner CONFIG REQUIRED)
        find_package(roadrunner_c_api CONFIG REQUIRED)

8. Edit this to edit the following line 

        find_package(libsbml-static CONFIG REQUIRED)
    
    To

        find_package(sbml-static CONFIG REQUIRED)


    and 

        find_package(sundials CONFIG REQUIRED)

    to 

        find_package(SUNDIALS CONFIG REQUIRED)

    or just copy and paste the following into the ImportRoadrunnerAndDependencies.cmake file

        find_package(Threads) # for libxml2, FindThreads.cmake is shipped with cmake
        find_package(LibLZMA) # for libxml2, LibLZMA.cmake is shipped with cmake
        find_package(zlib CONFIG REQUIRED)
        find_package(bzip2 CONFIG REQUIRED)
        find_package(iconv CONFIG REQUIRED)
        find_package(LibXml2 CONFIG REQUIRED)
        find_package(sbml-static CONFIG REQUIRED)
        find_package(rr-libstruct CONFIG REQUIRED)
        find_package(clapack CONFIG REQUIRED)
        find_package(nleq1 CONFIG REQUIRED)
        find_package(nleq2 CONFIG REQUIRED)
        find_package(PocoFoundation CONFIG REQUIRED)
        find_package(PocoNet CONFIG REQUIRED)
        find_package(PocoXML CONFIG REQUIRED)
        find_package(SUNDIALS CONFIG REQUIRED)
        find_package(LLVM REQUIRED)
        find_package(roadrunner-static CONFIG REQUIRED)
        find_package(roadrunner CONFIG REQUIRED)
        find_package(roadrunner_c_api CONFIG REQUIRED)

10. Go to the BNGMM src directory and compile the executable

        cd ../../../../../BNGMM/src
        cmake .
        make

    If it errors, and says something about the Eigen library (i.e cannot #include<Eigen/Dense>), taken from [stackoverflow](https://stackoverflow.com/questions/23284473/fatal-error-eigen-dense-no-such-file-or-directory), the solution would be to, depending on the Linux distribution and version,

        cd /usr/local/include
        sudo ln -sf eigen3/Eigen Eigen
        sudo ln -sf eigen3/unsupported unsupported

    or 

        cd /usr/include
        sudo ln -sf eigen3/Eigen Eigen
        sudo ln -sf eigen3/unsupported unsupported

11. Exit the src directory and run the options screen for all possible commands

        cd ..
        ./BNGMM -h 

12. As a final note, one can take out the BNGMM executable and throw it anywhere for use in any other directory. The rest of the directory can be discarded if space is a major concern.

 
## Prequisites to Compiling <a name="prq"></a>

### *Eigen* <a name="eig"></a>
BNGMMFit uses the Eigen 3.3.9 C++ linear algebra library for matrix computations on Linux. If on Ubuntu, you can do a quick install using:

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

and note that it is imperative that you have libRoadRunner installed as well, please look [here](https://libroadrunner.readthedocs.io/en/latest/Installation/installation.html)

## Compilation <a name="compilation"></a>

If you wish to modify the code for your own use or if the binary is not sufficient, a cmake has been provided in the /src directory. Fair warning this can be a tedious and bug-prone
process, that being said, assuming you have installed all of **Boost** and **Eigen** libraries through the above, then you can simply just download a fully built roadrunner + CyGMM library
[here (Defunct)](https://drive.google.com/drive/folders/1n6S7y2sf88mb_62evYO-RimbdT9FHMwh?usp=sharing). (outdated will remove section later)

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

    ./BNGMM

in your terminal. For more information about parameters and writing your own system, look below.

### Simple Example

Examples of parameter estimation tasks and models are in the example/ directory of the repo. In particular, the ones used in the [Software Paper](https://www.biorxiv.org/content/10.1101/2022.12.08.519526v1) (revision in progress as of 5/3/2023 so may be outdated currently) were of the linear 6 protein model under the "6_pro_lin_sim/" or "l6p_5p_sim/" directories, the yeast nonlinear model under "yeast/" directory, the 6 protein nonlinear model under "6_pro_nonlinear_slim/" directory, and the CD8+ T cell problem under "4_prot_CD3_CD8_CD28/". Please note that each directory should contain two directories, "X/" and "Y/" with their respective time snapshot data files, each should contain some .bngl file that would be ran for parameter estimation, and each should contain a time_steps.csv file that defines the time points of evolution that are being analyzed. Some directories also contain a true_rates.csv file that contain the rate constants (in order i.e topmost value corresponds to the first rate constant defined in the .bngl file) that were used for simulating the data (for examples containing simulated datasets). 

An example of one running BioNetGMMFit for a simple linear 3 protein problem (in CLI form) would be the following

    ./BNGMM -m example/3_prot_linear_sim/model.bngl -x example/3_prot_linear_sim/X -y example/3_prot_linear_sim/Y -t example/3_prot_linear_sim/time_steps.csv -r example/3_prot_linear_sim/true_rates.csv -c example/3_prot_linear_sim/Config.csv -o test/l3p/ --contour k1 k2

where it will run parameter estimation and save the results to the test/l3p/ directory. It will generate some contour plots as well and simulate the data at the given time steps in the time_steps.csv file. 

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

after moving or removing any previous Yt/Y0 files. Keep in mind, that the number of species of proteins you wish to simulate will correspond to the number of columns in each input X/Y file and the number of rows correspond to the cell count. Furthermore, note that if you have multiple time data, please make sure to label each time point with a "_tn" time tag such that the order of files are loaded in with the appropriate time point i.e _t2 for time point 2.

Furthermore, in docker to specify directories please do:

    docker run -t BNGMM -x <X-dir-of-initial-conditions> -y <Y-dir-of-SnapshotData>

or if statically compiled

    ./BNGMM -x <x-dir> -y <y-dir>

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
Make sure to list your times for time evolutions in time_steps.csv rowwise. For a single time evolution, only two time points, the end and start time of your evolution interval, is needed in the file.

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

### *Configuration i.e Hyperparameter Inputs* <a name="config"></a>
To set the parameters you want for your estimation run, double click or open the

    Config.csv

file and change the respective values. For instance, through excel 

| n_particles1 | 15 |
|--------------|----|

or using a default text editor

    n_particles1,15

sets the number of particles in blind pso to 15.

The default PSO parameters are listed below, Note that many of the hyperparameters were removed for simplicity, (in progress) removal of these from the configuration file will take.

| Parameter                        | Value | Explanation                                                              |
|----------------------------------|-------|--------------------------------------------------------------------------|
| Number of Particles PSO          | 1000  | Sets number of particles in PSO                                          |
| Number of Steps PSO              | 10    | Sets number of steps for  PSO                                            |
| Exclude Mixed Moments?           | 0     | 1 to use only means and variances, 0 otherwise                           |
| Exclude Mixed and Second Moments?| 0     | 1 to use only means, 0 otherwise                                         |
| Number of Runs                   | 1     | Sets total number of PSO runs for estimation                             |
| Simulate Y_t?                    | 1     | 1 to simulate Yt with a true rate vector, 0 to provide own Yt matrix     |
| Use Matrix Inverse?              | 0     | 1 to use C++'s Matrix Inverse, 0 otherwise                               |
| Number of Rates                  | 5     | Sets number of parameters to be estimated                                |
| Hypercube Dimension              | 1.0   | Real Value Bounds of Hypercube to be searched in PSO.                    |
| Report Moments?                  | 1     | 1 to report predicted moments in out.txt                                 |
| Bootstrap?                       | 1     | 1 to estimate 95% CI's, 0 otherwise                                      |   
| Use Deterministic?               | 1     | 1 to use CVode integrators, 0 to use roadrunner gillespie simulation     |
| Number of BNGL Steps             | 15    | Tuning Parameter for number of steps of integration                      |
| Seed                             | -1    | Used to seed the PSO, Off when seed < 0, On when seed > 0                |
| Parallel Number of Threads       | 8     | Number of threads to parallelize on.                                     |
| Initial Particle Best Weight     | 3.0   | How much historical weight (i.e last known particle position) to affect PSO step.|
| Global Best Weight               | 1.0   | How much weight best particle affects next PSO Step.                     |
| Particle Inertial Weight         | 6.0   | PSO Particle Inertia Component (to avoid local minima)                   |


By default, the PSO runs with all moments, with means, variances, and covariances. Currently, there are only two other options for specifying which estimators to use. For instance, set

    Exclude Mixed Moments?,1

to use means and second moments only while

    Exclude Mixed and Second Moments?, 1

will force the program to estimate rate constants using means only. All boolean options such as "Use Linear Model?" are set to on with 1, and set to off with 0. 

Finally, regarding holding parameter or rate constant values, these are currently only enabled for the nonlinear system where it's necessary for accurate estimation. 

## ** Defining a Model in BioNetGen ** <a name="bngl"></a>
One key note of importance is that this piece of code is not a full merger of BioNetGen and CyGMM (PSO + GMM). The reality is that all BNGL written is directly converted into sbml such that there are two key notes of importance: 

1. BioNetGMMFit can only do mechanistic (ODE) and stochastic (gillespie) models.
2. All models defined must have the following lines at the end of the BNGL file.

    generate_network()
    writeSBML()

With that being said, let's delve into a quick example. Start by creating a model.bngl file. If you're using VScode, this is relatively painless as shown below.

    ![create](/img/bnglStart.png)

Once created, open the file and start your model up with:

    begin model

    end model

Then, first define parameters to be estimated:

    begin model
        begin parameters
            k1 0.1
            k2 0.1
            k3 0.95
            k4 0.17
            k5 0.05
        end parameters
    end model

Then, define the species of observed data. In this case, we have 4. Note that the values defined here are placeholder values. Furthermore, note that if the species are used in reactions, their respective observable definitions are required . (If there are still species in your system that are not measured, but you wish to give them some value see [here]())

    begin model
        begin parameters
            k1 0.1
            k2 0.1
            k3 0.95
            k4 0.17
            k5 0.05
        end parameters
        begin species
            pCD3z() 192.7959
            pSLP76() 1463.265
            pErk() 5.251
            pS6() 435.2968
        end species
    end model

Now let's define the reaction rules.

    begin model
            begin parameters
                k1 0.1
                k2 0.1
                k3 0.95
                k4 0.17
                k5 0.05
            end parameters
            begin species
                pCD3z() 192.7959
                pSLP76() 1463.265
                pErk() 5.251
                pS6() 435.2968
            end species
            begin reaction rules
                0 -> pCD3z() k1
                pCD3z() -> pSLP76() + pCD3z() k2
                pSLP76() -> pErk() + pSLP76() k3
                pErk() -> pS6() + pErk() k4   
                pCD3z() -> 0 k5
                pSLP76() -> 0 k5
                pErk() -> 0 k5
                pS6() -> 0 k5
            end reaction rules
    end model

Finally, make sure to write the generate_network(), writeSBML() at the end to convert to SBML for use.

    begin model
            begin parameters
                k1 0.1
                k2 0.1
                k3 0.95
                k4 0.17
                k5 0.05
            end parameters
            begin species
                pCD3z() 192.7959
                pSLP76() 1463.265
                pErk() 5.251
                pS6() 435.2968
            end species
            begin reaction rules
                0 -> pCD3z() k1
                pCD3z() -> pSLP76() + pCD3z() k2
                pSLP76() -> pErk() + pSLP76() k3
                pErk() -> pS6() + pErk() k4   
                pCD3z() -> 0 k5
                pSLP76() -> 0 k5
                pErk() -> 0 k5
                pS6() -> 0 k5
            end reaction rules
    end model
    generate_network()
    writeSBML()

### Defining Proteins of Interest <a name="poi"></a>

Please note that should one need to define proteins that have not been observed with some initial values. They can simply just define their presence with their values normally as shown above. However, the user will need to specify a text .txt file to BioNetGMMFit with the proteins that have been observed in order. In this case, make sure to include in the .txt file their observables name. For instance, if I only observed pCD3z() and pS6() in the data such that pCD3z() was column 1 and pS6() was column 2 in the data files, we would specify in the poi.txt file:

    pCD3z()
    pS6()

and then call BioNetGMMFit with the -p tag as such.

    ./BNGMM -p poi.txt

Example BioNetGen Models and their .bngl files can be found in example/ such as the one seen with /6_pro_nonlin_sim_slim/6pro.bngl.
More documentation of BioNetGen can be found [here](http://bionetgen.org/)

Other command line arguments can be shown in the table below.

![cliargs](/img/CommandLineArguments.png)

## **Directory Structure (In Progress)** ## <a name="dir"></a>

### *BNGMM*
Main Directory with general configuration files and system.cpp code files.

### *src*
Contains all C++ source and header files needed to recompile the program.

### *example*
Contains various examples for use with the linear and nonlinear system provided by default in the program / code. There should be codes associated with them. The 4 protein CD8 T Cell data was sourced from [here](https://dpeerlab.github.io/dpeerlab-website/dremi-data.html).



### *(not updated to match most recent C++ version)* Mounting Volumes and Running Command Line Version of BNGMM Docker
  
Mounting Volumes to Write/Feed in Own Configuration Files (in the process of writing a python script to run this process with docker), do:

    docker run --rm -v $PWD:/data jhnwu3/bngmm -c /data/Config4pro.csv -m /data/4proV2.bngl -t /data/time_steps4.csv -x /data/X -y /data/Y -o /data

Note each directory and observe which each parameter corresponds to.

    For more information see [here](https://docs.docker.com/storage/volumes/) on how to mount a drive to be able to input configuration and data files.

