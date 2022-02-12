/* 
Author: John W., Dr. Stewart
Summary of Src File: 
- main, main function that calls all fileIO requirements

*/

#include "main.hpp" // all necessary libraries
#include "linear.hpp" // only linear model
#include "nonlinear.hpp" // only nonlinear model
#include "fileIO.hpp" // reading input and output
#include "calc.hpp" // common calc functions
#include "system.hpp" // user defined ode systems
int main(){
    auto t1 = std::chrono::high_resolution_clock::now();
    double scaleFactor = 20;
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
    double heldThetaVal = 0;
    VectorXd times = readCsvTimeParam();
    if(times.size() < 1){
        cout << "Error! Unable to read in timesteps properly or number of time steps inputted is equal to 0" << endl;
        exit(1);
    }
    
    cout << "Reading in data!" << endl;
    if(readCsvPSO(nParts, nSteps, nParts2, nSteps2, useOnlySecMom, useOnlyFirstMom, useLinear, nRuns, simulateYt, useInverse, nRates, heldTheta, heldThetaVal) !=0 ){
        cout << "failed to effectively read in parameters!" << endl;
        return EXIT_FAILURE;
    }

    MatrixXd X_0;
    X_0 = readX("../data/X");
    // X_0 = X_0 / scaleFactor;
    cout << "Reading in " << X_0.rows() << " elements from X data directory" << endl;
    int nMoments = (X_0.cols() * (X_0.cols() + 3)) / 2;
    if(useOnlySecMom){  // these will be added to the options sheet later.
        cout << "USING NONMIXED MOMENTS!!" << endl;
        nMoments = 2 * X_0.cols();
    }
    if(useOnlyFirstMom){
        cout << "USING ONLY MEANS!" << endl;
        nMoments = X_0.cols();
    }
    cout << "Moments:" << nMoments << endl;

    MatrixXd GBMAT;
    MatrixXd GBVECS = MatrixXd::Zero(nRuns, nRates + 1);
    if(useLinear == 1){
        for(int r = 0; r < nRuns; ++r){
            GBMAT = linearModel(nParts, nSteps, nParts2, nSteps2, X_0, nRates, nMoments, times, simulateYt, useInverse);
            GBVECS.row(r) = GBMAT.row(GBMAT.rows() - 1);
        }
    }else{
        /*---------------------- Nonlinear Setup ------------------------ */

        double t0 = 0, dt = 1.0; // nonlinear time evolution variables
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
        cout << "--------- Parameters ---------" << endl;
        cout << "X Size:" << X_0.rows() << " with:" << nMoments << " moments." << endl;
        cout << "Using Times:" << times.transpose() << endl;
        cout << "Bounds for Uniform Distribution (" << low << "," << high << ")"<< endl;
        cout << "Blind PSO --> nParts:" << nParts << " Nsteps:" << nSteps << endl;
        cout << "Targeted PSO --> nParts:" <<  nParts2 << " Nsteps:" << nSteps2 << endl;
        cout << "sdbeta:" << sdbeta << endl;
        cout << "------------------------------" << endl;

        MatrixXd PBMAT(nParts, nRates + 1); // particle best matrix + 1 for cost component
        MatrixXd POSMAT(nParts, nRates); // Position matrix as it goees through it in parallel

        /* Solve for Y_t (mu). */
        struct K tru;
        tru.k = readRates(nRates); // read in rates.

        vector<MatrixXd> Yt3Mats;
        vector<VectorXd> Yt3Vecs;
        Controlled_RK_Stepper_N controlledStepper;
        double trukCost = 0;
        if(simulateYt == 1){
            cout << "------ SIMULATING YT! ------" << endl;
            MatrixXd Y_0 = readY("../data/Y")[0];
            for(int t = 0; t < times.size(); t++){
                Nonlinear_ODE trueSys(tru);
                Protein_Components Yt(times(t), nMoments, Y_0.rows(), X_0.cols());
                Moments_Mat_Obs YtObs(Yt);
                for (int i = 0; i < Y_0.rows(); ++i) {
                    State_N y0 = convertInit(Y_0.row(i));
                    Yt.index = i;
                    integrate_adaptive(controlledStepper, trueSys, y0, t0, times(t), dt, YtObs);
                }
                Yt.mVec /= Y_0.rows();
                Yt3Mats.push_back(Yt.mat);
                Yt3Vecs.push_back(Yt.mVec);
            }
            cout << "---------------------------" << endl;
        }else{
            Yt3Mats = readY("../data/Y");
            // Yt3Mats[0] = Yt3Mats[0] / scaleFactor;// temporarily normalize it
            if(Yt3Mats.size() != times.size()){
                cout << "Error, number of Y_t files read in do not match the number of timesteps!" << endl;
                exit(1);
            }
            for(int i = 0; i < Yt3Mats.size(); i++){
                Yt3Vecs.push_back(momentVector(Yt3Mats[i], nMoments));
            }
        }

        cout << "Computing Weight Matrices!" << endl;
        /* Compute initial wolfe weights */
        for(int t = 0; t < times.size(); ++t){
            weights.push_back(wolfWtMat(Yt3Mats[t], nMoments, false));
        }
        cout << weights[0] << endl;
        
        for(int run = 0; run < nRuns; ++run){
            // make sure to reset GBMAT, POSMAT, AND PBMAT every run
            double sfi = sfe, sfc = sfp, sfs = sfg; // below are the variables being used to reiterate weights
            GBMAT = MatrixXd::Zero(0,0); // iterations of global best vectors
            MatrixXd PBMAT = MatrixXd::Zero(nParts, nRates + 1); // particle best matrix + 1 for cost component
            MatrixXd POSMAT = MatrixXd::Zero(nParts, nRates); // Position matrix as it goees through it in parallel
            
            /* Initialize Global Best  */
            double holdTheta2 = 0.1;
            struct K seed;
            seed.k = VectorXd::Zero(nRates); 
            for (int i = 0; i < nRates; i++) { seed.k(i) = rndNum(low,high);}
            if(heldTheta > -1){seed.k(heldTheta) = heldThetaVal;}
            
            /* Evolve initial Global Best and Calculate a Cost*/
            double costSeedK = 0;
            for(int t = 0; t < times.size(); t++){
                Protein_Components Xt(times(t), nMoments, X_0.rows(), X_0.cols());
                Moments_Mat_Obs XtObs(Xt);
                Nonlinear_ODE sys(seed);
                for (int i = 0; i < X_0.rows(); ++i) {
                    State_N c0 = convertInit(X_0.row(i));
                    Xt.index = i;
                    integrate_adaptive(controlledStepper, sys, c0, t0, times(t), dt, XtObs);
                }
                Xt.mVec /= X_0.rows();  
                costSeedK += costFunction(Yt3Vecs[t], Xt.mVec, weights[t]);
            }
            cout << "seedk:"<< seed.k.transpose() << "| cost:" << costSeedK << endl;
            
            double gCost = costSeedK; //initialize costs and GBMAT
            VectorXd GBVEC = seed.k;
            
            GBMAT.conservativeResize(GBMAT.rows() + 1, nRates + 1);
            for (int i = 0; i < nRates; i++) {
                GBMAT(GBMAT.rows() - 1, i) = seed.k(i);
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
                        POSMAT(particle, 1) = holdTheta2;
                        struct K pos;
                        pos.k = VectorXd::Zero(nRates);
                        for(int i = 0; i < nRates; i++){
                            pos.k(i) = POSMAT(particle, i);
                        }
                        
                        double cost = 0;
                        for(int t = 0; t < times.size(); ++t){
                            Nonlinear_ODE initSys(pos);
                            Protein_Components XtPSO(times(t), nMoments, X_0.rows(), X_0.cols());
                            Moments_Mat_Obs XtObsPSO(XtPSO);
                            for(int i = 0; i < X_0.rows(); ++i){
                                //State_N c0 = gen_multi_norm_iSub();
                                State_N c0 = convertInit(X_0.row(i));
                                XtPSO.index = i;
                                integrate_adaptive(controlledStepper, initSys, c0, t0, times(t), dt, XtObsPSO);
                            }
                            XtPSO.mVec/= X_0.rows();
                            cost += costFunction(Yt3Vecs[t], XtPSO.mVec, weights[t]);
                        }
                        
                        /* instantiate PBMAT */
                        for(int i = 0; i < nRates; i++){
                            PBMAT(particle, i) = POSMAT(particle, i);
                        }
                        PBMAT(particle, nRates) = cost; // add cost to final column
                    }else{ 
                        /* step into PSO */
                        double w1 = sfi * rndNum(low,high) / sf2, w2 = sfc * rndNum(low,high) / sf2, w3 = sfs * rndNum(low,high) / sf2;
                        double sumw = w1 + w2 + w3; //w1 = inertial, w2 = pbest, w3 = gbest
                        w1 = w1 / sumw; w2 = w2 / sumw; w3 = w3 / sumw;
                        struct K pos;
                        pos.k = VectorXd::Zero(nRates);
                        pos.k = POSMAT.row(particle);
                        VectorXd rpoint = adaptVelocity(pos.k, particle, epsi, nan, hone);
                        VectorXd PBVEC(nRates);
                        for(int i = 0; i < nRates; ++i){PBVEC(i) = PBMAT(particle, i);}
                        
                        pos.k = w1 * rpoint + w2 * PBVEC + w3 * GBVEC; // update position of particle
                        
                        // if(rndNum(low,high) < probabilityToTeleport){ // hard coded grid re-search for an adaptive component
                        //     pos.k(0) = rndNum(low,high);
                        //     pos.k(1) = rndNum(low,high);
                        //     pos.k(4) = rndNum(low,high);
                        // }
                        if(heldTheta > -1){
                            pos.k(heldTheta) = heldThetaVal;
                        }
                        POSMAT.row(particle) = pos.k;
                        double cost = 0;
                        VectorXd temp;
                        for(int t = 0; t < times.size(); ++t){
                            /*solve ODEs and recompute cost */
                            Protein_Components XtPSO(times(t), nMoments, X_0.rows(), X_0.cols());
                            Moments_Mat_Obs XtObsPSO1(XtPSO);
                            Nonlinear_ODE stepSys(pos);
                            for(int i = 0; i < X_0.rows(); i++){
                                State_N c0 = convertInit(X_0.row(i));
                                XtPSO.index = i;
                                integrate_adaptive(controlledStepper, stepSys, c0, t0, times(t), dt, XtObsPSO1);
                            }
                            XtPSO.mVec/=X_0.rows();
                            temp = XtPSO.mVec;
                            cost += costFunction(Yt3Vecs[t], XtPSO.mVec, weights[t]);
                        }
                    
                        /* update gBest and pBest */
                    #pragma omp critical
                    {
                        if(cost < PBMAT(particle, nRates)){ // particle best cost
                            for(int i = 0; i < nRates; i++){
                                PBMAT(particle, i) = pos.k(i);
                            }
                            PBMAT(particle, nRates) = cost;
                            if(cost < gCost){
                                gCost = cost;
                                GBVEC = pos.k;
                                
                                cout << "Moment Estimates:" << temp.transpose() << endl;
                                
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

            for(int i = 0; i < nRates; i++){GBVECS(run, i) = GBVEC(i);}
            GBVECS(run, nRates) = gCost;


            /* bootstrap X0 and Y matrices if more than 1 run is specified */
            if(nRuns > 1){
                X_0 = bootStrap(X_0);
                for(int t = 0; t < times.size(); ++t){
                    Yt3Mats[t] = bootStrap(Yt3Mats[t]);
                    Yt3Vecs[t] = momentVector(Yt3Mats[t], nMoments);
                }
            }
        }
        if(simulateYt == 1){cout << "Simulated Truth:" << tru.k.transpose() << endl;}
    }
    cout << "All Run Results:" << endl;
    cout << GBVECS << endl;
    /* Compute 95% CI's with basic z=1.96 normal distribution assumption for now if n>1 */
    if(nRuns > 1){       
        computeConfidenceIntervals(GBVECS, 1.96, nRates);
    }
    auto tB = std::chrono::high_resolution_clock::now();
    auto bDuration = std::chrono::duration_cast<std::chrono::seconds>(tB - t1).count();
    cout << "CODE FINISHED RUNNING IN " << bDuration << " s TIME!" << endl;
    return 0;
}