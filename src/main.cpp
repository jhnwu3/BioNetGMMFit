/* Have one boolean variable (number) to decide whether or not to use linear or nonlinear */
/* 
List of different functions between linear and nonlinear.
#1: Different recomputation of weight matrices depending on if you're using inversion or not.
#2: Different updating mechanism for rpoints (new position vector)
#3: Different Moment Vector mechanic.
#4: Different number of weight matrices

Parameters of PSO

particles, steps, using cross moments or no, X_0 and Y_0 -> nSpecies, number of rate constants in model used.

Later Todo: -> modularize code such that we can just swap in models.

Note: Nonlinear.hpp should contain all the boost ODE stuff and the nonlinear weight function and PSO.
Note: linear.hpp/cpp will store literally the linear PSO function and their respective different functions.
Note: main.hpp will contain the inclusion of all the needed libraries such that all other .hpp files can link directly to it (single point)
*/
#include "main.hpp" // all necessary libraries
#include "linear.hpp" // only linear model
#include "nonlinear.hpp" // only nonlinear model
#include "fileIO.hpp" // reading input and output
#include "calc.hpp" // common clac functions
int main(){
    auto t1 = std::chrono::high_resolution_clock::now();

    /* some default values for PSO just in case. */
    int nParts1 = 25; // first part PSO
    int nSteps1 = 50;
    int nParts2 = 10; // second part PSO
    int nSteps2 = 1000;
    int useOnlySecMom = 1;
    int useOnlyFirstMom = 1;
    int useLinear = 0;
    int xDataSize = 0;
    int yDataSize = 0;
    int nSpecies = 0;
    int nRates = 0;
    int nRuns = 0;
    
    cout << "Reading in data!" << endl;
    if(readCsvPSO(nParts1, nSteps1, nParts2, nSteps2, useOnlySecMom, useOnlyFirstMom, useLinear, nRuns) != 0 || 
        readCsvDataParam(nSpecies, nRates, xDataSize, yDataSize) != 0){
        cout << "failed to effectively read in parameters!" << endl;
        return EXIT_FAILURE;
    }

    cout << "nParts1:" << nParts1 << endl;
    cout << "nSteps1:" << nSteps1 << endl;
    cout << "nParts2:" << nParts2 << endl;
    cout << "nSteps2:" << nSteps2 << endl;

    MatrixXd X_0 = csvToMatrix("../data/testXr.csv", xDataSize);
    MatrixXd Y_0 = csvToMatrix("../data/testYr.csv", yDataSize);
    cout << "X_0:" << "(" << X_0.rows() << "," << X_0.cols() << ")" << endl;
    cout << "Y_0:" << "(" << Y_0.rows() << "," << Y_0.cols() << ")" << endl;
    int nMoments = (X_0.cols() * (X_0.cols() + 3)) / 2;
    cout << "moments:" << nMoments << endl;
    /* Temp Initial Conditions */
    // MatrixXd X_0(xDataSize, nSpecies);
    // MatrixXd Y_0(yDataSize, nSpecies);
    // X_0 = txtToMatrix("input/knewX.0.txt", xDataSize, nSpecies);
    // Y_0 = txtToMatrix("input/knewY.0.txt", yDataSize, nSpecies);
    
    /* Fixed Initial Conditions */

    // cout << "Using starting row of data:" << startRow << " and " << N << " data pts!" << endl;
    // cout << "first row X0:" << X_0.row(0) << endl;
    // cout << "final row X0:" << X_0.row(N - 1) << endl << endl << endl << endl;

    if(useOnlySecMom){  // these will be added to the options sheet later.
        cout << "USING NONMIXED MOMENTS!!" << endl;
        nMoments = 2 * X_0.cols();
    }
    if(useOnlyFirstMom){
        cout << "USING ONLY MEANS!" << endl;
        nMoments = X_0.cols();
    }
    MatrixXd GBMAT;
    if(useLinear == 1){
        GBMAT = linearModel(nParts1, nSteps1, nParts2, nSteps2, X_0, Y_0, nRates, nMoments);
    }else{
        GBMAT = nonlinearModel(nParts1, nSteps1, nParts2, nSteps2, X_0, Y_0, nRates, nRuns, nMoments);
    }

    auto t2 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count();
    cout << "CODE FINISHED RUNNING IN " << duration << " s TIME!" << endl;
    return 0;
}