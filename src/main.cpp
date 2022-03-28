/* 
Author: John W., Dr. Stewart
Summary of Src File: 
CyGMM -> <Insert Statement About Rate Constant Estimation>

Coder's Note:
If you're wondering why the rungekutta solution method is all in main and not segmented in nonlinear.cpp is because currently
the way in how Boost's ODE struct/class system for solving ODEs results in major performance decreases when inside a function, very possibly
due to something with memory or how functions work from the compiler.

*/

#include "main.hpp" // all necessary libraries
#include "linear.hpp" // only linear model
#include "nonlinear.hpp" // only nonlinear model
#include "fileIO.hpp" // reading input and output
#include "calc.hpp" // common calc functions
#include "system.hpp" // user defined ode systems
int main(){
    auto t1 = std::chrono::high_resolution_clock::now();
    cout << "Program Begin:" << endl;
    /* Input Parameters for Program */
    int nParts = 25; // first part PSO
    int nSteps = 50;
    int nParts2 = 10; // second part PSO
    int nSteps2 = 100;
    int useOnlySecMom = 1;
    int useOnlyFirstMom = 1;
    int useLinear = 0;
    int nRates = 0;
    int nRuns = 0;
    int simulateYt = 1;
    int useInverse = 0; // currently just inverse only occurs in linear model.
    int heldTheta = -1;
    int reportMoments = 1;
    int nest = 0;
    int bootstrap = 0;
    double heldThetaVal = 0;
    double hyperCubeScale = 1.0;
    VectorXd times = readCsvTimeParam();
    if(times.size() < 1){
        cout << "Error! Unable to read in timesteps properly or number of time steps inputted is equal to 0" << endl;
        exit(1);
    }
    
    if(readCsvPSO(nParts, nSteps, nParts2, nSteps2, useOnlySecMom, useOnlyFirstMom, useLinear, nRuns, simulateYt, useInverse, nRates, heldTheta, heldThetaVal, reportMoments, hyperCubeScale, nest, bootstrap) !=0 ){
        cout << "failed to effectively read in parameters!" << endl;
        return EXIT_FAILURE;
    }

    MatrixXd X_0;
    MatrixXd ogX_0;
    X_0 = readX("../data/X");
    X_0 = filterZeros(X_0);
    ogX_0 = X_0;
    cout << "After removing all negative rows, X has " << X_0.rows() << " rows." << endl;
    // cout << "---------" << endl << X_0 << endl << "--------" << endl;
    int nMoments = (X_0.cols() * (X_0.cols() + 3)) / 2;
    if(useOnlySecMom){  // these will be added to the options sheet later.
        nMoments = 2 * X_0.cols();
    }
    if(useOnlyFirstMom){
        nMoments = X_0.cols();
    }
    printParameters(nParts, nSteps, nParts2, nSteps2, useOnlySecMom, useOnlyFirstMom, useLinear, nRuns, simulateYt, useInverse, nRates, heldTheta, heldThetaVal, reportMoments, hyperCubeScale, times, nMoments);
    MatrixXd GBMAT;
    MatrixXd GBVECS = MatrixXd::Zero(nRuns, nRates + 1);
    if(useLinear == 1){
        for(int r = 0; r < nRuns; ++r){
            GBMAT = linearModel(nParts, nSteps, nParts2, nSteps2, X_0, nRates, nMoments, times, simulateYt, useInverse);
            GBVECS.row(r) = GBMAT.row(GBMAT.rows() - 1);
        }
    }else{
        /*---------------------- Nonlinear Setup ------------------------ */

        double dt = 1.0; // nonlinear time evolution variables
        /* Explicit Boundary Parameters */
        double squeeze = 0.500, sdbeta = 0.10; // how much to shrink PSO search over time (how much variability each position is iterated upon)
        double boundary = 0.001;
        /* SETUP */
        double sf2 = 1; // factor that can be used to regularize particle weights (global, social, inertial)
        double epsi = 0.02;
        double nan = 0.005;
        /* PSO params */
        double sfp = 3.0, sfg = 1.0, sfe = 6.0; // initial particle historical weight, global weight social, inertial
        double sfi = sfe, sfc = sfp, sfs = sfg; // below are the variables being used to reiterate weights
        double alpha = 0.2;
        int hone = 28; 
        int startRow = 0;
        double low = 0.0, high = 1.0; // boundaries for PSO rate estimation, 0 to 1.0
        random_device RanDev;
        mt19937 gen(RanDev());
        vector<MatrixXd> weights;
        MatrixXd PBMAT(nParts, nRates + 1); // particle best matrix + 1 for cost component
        MatrixXd POSMAT(nParts, nRates); // Position matrix as it goees through it in parallel

        /* Solve for Y_t (mu). */
        VectorXd tru;
        vector<MatrixXd> Yt3Mats;
        vector<MatrixXd> ogYt3Mats;
        vector<VectorXd> Yt3Vecs;
        Controlled_RK_Stepper_N controlledStepper;
        double trukCost = 0;
        if(simulateYt == 1){
            cout << "------ SIMULATING YT! ------" << endl;
            tru = readRates(nRates);
            MatrixXd Y_0 = readY("../data/Y")[0];
            Y_0 = filterZeros(Y_0);
            cout << "After removing all negative rows, Y has " << Y_0.rows() << " rows." << endl;
            for(int t = 1; t < times.size(); t++){ // start at t1, because t0 is now in the vector
                Nonlinear_ODE trueSys(tru);
                Protein_Components Yt(times(t), nMoments, Y_0.rows(), X_0.cols());
                Moments_Mat_Obs YtObs(Yt);
                for (int i = 0; i < Y_0.rows(); ++i) {
                    State_N y0 = convertInit(Y_0.row(i));
                    Yt.index = i;
                    integrate_adaptive(controlledStepper, trueSys, y0, times(0), times(t), dt, YtObs); 
                }
                Yt.mVec /= Y_0.rows();
                Yt3Mats.push_back(Yt.mat);
                Yt3Vecs.push_back(Yt.mVec);
            }
            cout << "---------------------------" << endl;
        }else{
            Yt3Mats = readY("../data/Y");
            if(Yt3Mats.size() + 1 != times.size()){
                cout << "Error, number of Y_t files read in do not match the number of timesteps!" << endl;
                exit(1);
            }
            // filter all zeroes and compute moments vectors for cost calcs
            for(int i = 0; i < Yt3Mats.size(); i++){
                Yt3Mats[i] = filterZeros(Yt3Mats[i]);
                cout << "Yt Means:" << Yt3Mats[i].colwise().mean() << endl;
                cout << "After removing all negative rows, Y"<< i << " has " << Yt3Mats[i].rows() << " rows." << endl;
                Yt3Vecs.push_back(momentVector(Yt3Mats[i], nMoments));
            }
            ogYt3Mats = Yt3Mats;
            // cout << "Y_0" << endl << Yt3Mats[0] << endl;
        }

        cout << "Computing Weight Matrices!" << endl;
        /* Compute initial wolfe weights */
        for(int y = 0; y < Yt3Mats.size(); ++y){ 
            weights.push_back(wolfWtMat(Yt3Mats[y], nMoments, false));
        }
        cout << weights[0] << endl;
       
        for(int run = 0; run < nRuns; ++run){ // for multiple runs aka bootstrapping (for now)
            VectorXd nestedHolds = VectorXd::Zero(nRates);
            if (run > 0 && bootstrap == 1){
                for(int y = 0; y < Yt3Mats.size(); ++y){ 
                    weights[y] = wolfWtMat(Yt3Mats[y], nMoments, false);
                }
            }
            for(int ne = 0; ne < nest; ne++){ // now we will include the ability to nest hypercubes. Every nested hypercube doubles the hypercube width.
                // make sure to reset GBMAT, POSMAT, AND PBMAT every run
                double sfi = sfe, sfc = sfp, sfs = sfg; // below are the variables being used to reiterate weights
                GBMAT = MatrixXd::Zero(0,0); // iterations of global best vectors
                MatrixXd PBMAT = MatrixXd::Zero(nParts, nRates + 1); // particle best matrix + 1 for cost component
                MatrixXd POSMAT = MatrixXd::Zero(nParts, nRates); // Position matrix as it goees through it in parallel
                
                /* Initialize Global Best  */
                VectorXd seed = VectorXd::Zero(nRates);
                for (int i = 0; i < nRates; i++) {seed(i) = rndNum(low,high);}
                if(heldTheta > -1){seed(heldTheta) = heldThetaVal;
                seed << 1.54103,	0.0417282,	0.311026,	1.1217,	0.384331,	0.425928;
                }
                
                /* Evolve initial Global Best and Calculate a Cost*/
                double costSeedK = 0;
                
                if(ne > 0){
                    for(int i = 0; i < nRates; i++){
                        if(nestedHolds(i)  != 0){
                            seed(i) = nestedHolds(i) / hyperCubeScale;
                        }
                    }
                }
                for(int t = 1; t < times.size(); t++){
                    Protein_Components Xt(times(t), nMoments, X_0.rows(), X_0.cols());
                    Moments_Mat_Obs XtObs(Xt);
                    Nonlinear_ODE sys(hyperCubeScale * seed);
                    for (int i = 0; i < X_0.rows(); ++i) {
                        State_N c0 = convertInit(X_0.row(i));
                        Xt.index = i;
                        integrate_adaptive(controlledStepper, sys, c0, times(0), times(t), dt, XtObs);
                    }
                    Xt.mVec /= X_0.rows();  
                    costSeedK += costFunction(Yt3Vecs[t - 1], Xt.mVec, weights[t - 1]); // indexed one off for each weight matrix.
                }
                cout << "PSO Seeded At:"<< seed.transpose() << "| cost:" << costSeedK << endl;
                
                double gCost = costSeedK; //initialize costs and GBMAT
                VectorXd GBVEC = seed;
                
                GBMAT.conservativeResize(GBMAT.rows() + 1, nRates + 1);
                for (int i = 0; i < nRates; i++) {
                    GBMAT(GBMAT.rows() - 1, i) = seed(i);
                }
                GBMAT(GBMAT.rows() - 1, nRates) = gCost;
                double probabilityToTeleport = 3.0/4.0; 
                /* Blind PSO begins */
                cout << "PSO Estimation Has Begun, This may take some time..." << endl;
                for(int step = 0; step < nSteps; step++){
                #pragma omp parallel for 
                    for(int particle = 0; particle < nParts; particle++){
                        /* initialize all particle rate constants with unifDist */
                        if(step == 0){
                            /* initialize all particles with random rate constant positions */
                            for(int i = 0; i < nRates; i++){
                                POSMAT(particle, i) = rndNum(low, high);
                            }
                            if(heldTheta > -1){POSMAT.row(particle)(heldTheta) = heldThetaVal;
                            POSMAT.row(particle)  << 1.54103,	0.0417282,	0.311026,	1.1217,	0.384331,	0.425928;
                            }
                            
                            double cost = 0;    
                            if(step == 0 && ne > 0){
                                for(int i = 0 ; i < nRates; i++){
                                    if(nestedHolds(i) != 0){
                                        POSMAT.row(particle)(i) = nestedHolds(i) / hyperCubeScale; 
                                    }
                                }
                            }
                            
                            for(int t = 1; t < times.size(); ++t){
                                Nonlinear_ODE initSys(POSMAT.row(particle) * hyperCubeScale);
                                Protein_Components XtPSO(times(t), nMoments, X_0.rows(), X_0.cols());
                                Moments_Mat_Obs XtObsPSO(XtPSO);
                                for(int i = 0; i < X_0.rows(); ++i){
                                    State_N c0 = convertInit(X_0.row(i));
                                    XtPSO.index = i;
                                    integrate_adaptive(controlledStepper, initSys, c0, times(0), times(t), dt, XtObsPSO);
                                }
                                XtPSO.mVec /= X_0.rows();
                                cost += costFunction(Yt3Vecs[t - 1], XtPSO.mVec, weights[t - 1]);
                            }
                            
                            /* instantiate PBMAT */
                            for(int i = 0; i < nRates; i++){
                                PBMAT(particle, i) = POSMAT(particle, i);
                            }
                            PBMAT(particle, nRates) = cost; // add cost to final column
                        }else{ 
                            /* step into PSO */
                            double w1 = sfi * rndNum(low,high) / sf2, w2 = sfc * rndNum(low,high) / sf2, w3 = sfs * rndNum(low,high) / sf2;
                            double sumw = w1 + w2 + w3; 
                            w1 = w1 / sumw; w2 = w2 / sumw; w3 = w3 / sumw;
                    
                            VectorXd rpoint = adaptVelocity(POSMAT.row(particle), particle, epsi, nan, hone);
                            VectorXd PBVEC(nRates);
                            for(int i = 0; i < nRates; ++i){PBVEC(i) = PBMAT(particle, i);}
                            
                            POSMAT.row(particle) = (w1 * rpoint + w2 * PBVEC + w3 * GBVEC); // update position of particle
                        
                            if(heldTheta > -1){POSMAT.row(particle)(heldTheta) = heldThetaVal;}
                            double cost = 0;

                            for(int i = 0 ; i < nRates; i++){
                                if(nestedHolds(i) != 0){
                                    POSMAT(particle, i) = nestedHolds(i) / hyperCubeScale; 
                                }
                            }
                            for(int t = 1; t < times.size(); ++t){
                                /*solve ODEs and recompute cost */
                                Protein_Components XtPSO(times(t), nMoments, X_0.rows(), X_0.cols());
                                Moments_Mat_Obs XtObsPSO1(XtPSO);
                                Nonlinear_ODE stepSys(POSMAT.row(particle) * hyperCubeScale);
                                for(int i = 0; i < X_0.rows(); ++i){
                                    State_N c0 = convertInit(X_0.row(i));
                                    XtPSO.index = i;
                                    integrate_adaptive(controlledStepper, stepSys, c0, times(0), times(t), dt, XtObsPSO1);
                                }
                                XtPSO.mVec /= X_0.rows();
                                cost += costFunction(Yt3Vecs[t - 1], XtPSO.mVec, weights[t - 1]);
                            }
                        
                            /* update gBest and pBest */
                        #pragma omp critical
                        {
                            if(cost < PBMAT(particle, nRates)){ // particle best cost
                                for(int i = 0; i < nRates; i++){
                                    PBMAT(particle, i) = POSMAT.row(particle)(i);
                                }
                                PBMAT(particle, nRates) = cost;
                                if(cost < gCost){
                                    gCost = cost;
                                    GBVEC = POSMAT.row(particle);
                                }   
                            }
                        }
                        }
                    }
                    GBMAT.conservativeResize(GBMAT.rows() + 1, nRates + 1); // Add to GBMAT after resizing
                    for (int i = 0; i < nRates; i++) {GBMAT(GBMAT.rows() - 1, i) = GBVEC(i);}
                    GBMAT(GBMAT.rows() - 1, nRates) = gCost;
                    sfi = sfi - (sfe - sfg) / nSteps;   // reduce the inertial weight after each step 
                    sfs = sfs + (sfe - sfg) / nSteps;
                }

                

                // double hypercube size and rescan global best vector to see what needs to be held constant.
                if(nest > 1){
                    for(int i = 0; i < GBVEC.size(); i++){
                        if(hyperCubeScale * GBVEC(i) < hyperCubeScale - epsi){
                            nestedHolds(i) = hyperCubeScale * GBVEC(i);
                        }
                    }
                }
                for(int i = 0; i < nRates; i++){
                    if(nestedHolds(i) != 0){
                        GBVECS(run, i) = nestedHolds(i); // GBVECS is either something that was held constant.
                    }else{ 
                        GBVECS(run, i) = hyperCubeScale * GBVEC(i); // or is now something scaled with GBVEC.
                    } 
                }
                GBVECS(run, nRates) = gCost;
                if(nest > 1){
                    hyperCubeScale *= 2.0;
                }
            } // nested loop

            /* bootstrap X0 and Y matrices if more than 1 run is specified */
            if(nRuns > 1 && bootstrap == 1){
                X_0 = bootStrap(ogX_0);
                for(int y = 0; y < Yt3Mats.size(); ++y){
                    Yt3Mats[y] = bootStrap(ogYt3Mats[y]);
                    Yt3Vecs[y] = momentVector(Yt3Mats[y], nMoments);
                }
                cout << "bootstrap means" << endl << "X_0:" << X_0.colwise().mean() << endl << "Yt:" << Yt3Mats[0].colwise().mean() << endl;
            }
            cout << "hypercubescale before reset:" << hyperCubeScale << endl;
            if(nest > 1){
                for(int ne = 0; ne < nest; ne++){ // reset cube for each run
                    hyperCubeScale /= 2.0;
                }
            }
            cout << "hypercubescale after reset:" << hyperCubeScale << endl;
            cout << "nested hypercubes:" << nestedHolds.transpose() << endl; 
           
        } // run loop
        // when done, find what the original max hypercube size was from nesting
        for(int ne = 1; ne < nest; ne++){ 
            hyperCubeScale *= 2.0;
        }
        cout << "hypercubescale after final nesting:" << hyperCubeScale << endl;
        if(simulateYt == 1){cout << "Simulated Truth:" << tru.transpose() << endl;}
        if(reportMoments == 1){
            for(int t = 1; t < times.size(); ++t){
                /*solve ODEs and recompute cost */
                Protein_Components XtPSO(times(t), nMoments, X_0.rows(), X_0.cols());
                Moments_Mat_Obs XtObsPSO1(XtPSO);
                Nonlinear_ODE stepSys(GBVECS.colwise().mean());
                for(int i = 0; i < X_0.rows(); ++i){
                    State_N c0 = convertInit(X_0.row(i));
                    XtPSO.index = i;
                    integrate_adaptive(controlledStepper, stepSys, c0, times(0), times(t), dt, XtObsPSO1);
                }
                XtPSO.mVec/=X_0.rows();
                cout << "GBVEC:" << GBVECS.colwise().mean() << endl;
                cout << "Simulated Xt Moments for time " << times(t) << ":" << XtPSO.mVec.transpose() << endl;
                // cout << "Final Evolved Matrix" << endl << XtPSO.mat << endl;
            }
        }
    }
    cout << endl << "-------------- All Run Estimates: -------------------" << endl;
    cout << GBVECS << endl;
    /* Compute 95% CI's with basic z=1.96 normal distribution assumption for now if n>1 */
    if(nRuns > 1){ computeConfidenceIntervals(GBVECS, 1.96, nRates);}
    
    auto tB = std::chrono::high_resolution_clock::now();
    auto bDuration = std::chrono::duration_cast<std::chrono::seconds>(tB - t1).count();
    cout << "CODE FINISHED RUNNING IN " << bDuration << " s TIME!" << endl;
    return 0;
}