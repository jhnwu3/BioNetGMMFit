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
    cout << "Program Begin:" << endl;
    /* some default values for PSO just in case. */
    int nParts = 25; // first part PSO
    int nSteps = 50;
    int nParts2 = 10; // second part PSO
    int nSteps2 = 100;
    int useOnlySecMom = 1;
    int useOnlyFirstMom = 1;
    int useLinear = 0;
    int sampleSize = 0;
    int nRates = 0;
    int nRuns = 0;
    int simulateYt = 1;
    int useInverse = 0; // currently just inverse only occurs in linear model.
    VectorXd times = readCsvTimeParam();
    if(times.size() < 1){
        cout << "Error! Unable to read in timesteps properly or number of time steps inputted is equal to 0" << endl;
        exit(1);
    }
    
    cout << "Reading in data!" << endl;
    if(readCsvPSO(nParts, nSteps, nParts2, nSteps2, useOnlySecMom, useOnlyFirstMom, useLinear, nRuns, simulateYt, useInverse, nRates, sampleSize) !=0 ){
        cout << "failed to effectively read in parameters!" << endl;
        return EXIT_FAILURE;
    }

    MatrixXd X_0;
    X_0 = readX("../data/X", sampleSize);
    cout << "X:" << X_0.rows() << endl;
    int nMoments = (X_0.cols() * (X_0.cols() + 3)) / 2;
    if(useOnlySecMom){  // these will be added to the options sheet later.
        cout << "USING NONMIXED MOMENTS!!" << endl;
        nMoments = 2 * X_0.cols();
    }
    if(useOnlyFirstMom){
        cout << "USING ONLY MEANS!" << endl;
        nMoments = X_0.cols();
    }
    cout << "moments:" << nMoments << endl;

    MatrixXd GBMAT;
    if(useLinear == 1){
        GBMAT = linearModel(nParts, nSteps, nParts2, nSteps2, X_0, nRates, nMoments, times, simulateYt);
    }else{
        /*---------------------- Nonlinear Setup ------------------------ */
        /* Variables (global) */
        double t0 = 0, dt = 1.0; // time variables
        double squeeze = 0.500, sdbeta = 0.10; 
        double boundary = 0.001;
        /* SETUP */
        int useDiag = 0;
        int sf1 = 1;
        int sf2 = 1;
        double epsi = 0.02;
        double nan = 0.005;
        /* PSO params */
        double sfp = 3.0, sfg = 1.0, sfe = 6.0; // initial particle historical weight, global weight social, inertial
        double sfi = sfe, sfc = sfp, sfs = sfg; // below are the variables being used to reiterate weights
        double alpha = 0.2;
        int hone = 28; 
        int startRow = 0;
        //nMoments = 2*N_SPECIES; // mean + var only!
        VectorXd wmatup(4);
        wmatup << 0.15, 0.35, 0.60, 0.9;
        double uniLowBound = 0.0, uniHiBound = 1.0;
        random_device RanDev;
        mt19937 gen(RanDev());
        uniform_real_distribution<double> unifDist(uniLowBound, uniHiBound);
        
        vector<MatrixXd> weights;

        for(int i = 0; i < times.size(); i++){
            weights.push_back(MatrixXd::Identity(nMoments, nMoments));
        }
        
        cout << "Using two part PSO " << "Sample Size:" << X_0.rows() << " with:" << nMoments << " moments." << endl;
        cout << "Using Times:" << times.transpose() << endl;
        cout << "Bounds for Uniform Distribution (" << uniLowBound << "," << uniHiBound << ")"<< endl;
        cout << "Blind PSO --> nParts:" << nParts << " Nsteps:" << nSteps << endl;
        cout << "Targeted PSO --> nParts:" <<  nParts2 << " Nsteps:" << nSteps2 << endl;
        cout << "sdbeta:" << sdbeta << endl;
        // cout << "wt:" << endl << wt << endl;

        MatrixXd GBMAT(0, 0); // iterations of global best vectors
        MatrixXd PBMAT(nParts, nRates + 1); // particle best matrix + 1 for cost component
        MatrixXd POSMAT(nParts, nRates); // Position matrix as it goees through it in parallel

        /* Solve for Y_t (mu). */
        cout << "Loading in Truk!" << endl;
        struct K tru;
        tru.k = readRates(nRates);
        // tru.k << 0.1, 0.1, 0.95, 0.17, 0.05, 0.18;

        cout << "Calculating Yt!" << endl;
        vector<MatrixXd> Yt3Mats;
        vector<VectorXd> Yt3Vecs;
        Controlled_RK_Stepper_N controlledStepper;
        double trukCost = 0;
        if(simulateYt == 1){
            cout << "SIMULATING YT!" << endl;
            MatrixXd Y_0 = readY("../data/Y", sampleSize)[0];
            for(int t = 0; t < times.size(); t++){
                Nonlinear_ODE6 trueSys(tru);
                Protein_Components Yt(times(t), nMoments, Y_0.rows(), X_0.cols());
                Protein_Components Xt(times(t), nMoments, X_0.rows(), X_0.cols());
                Moments_Mat_Obs YtObs(Yt);
                Moments_Mat_Obs XtObs(Xt);
                for (int i = 0; i < sampleSize; ++i) {
                    State_N y0 = convertInit(Y_0.row(i));
                    State_N x0 = convertInit(X_0.row(i));
                    Yt.index = i;
                    Xt.index = i;
                    integrate_adaptive(controlledStepper, trueSys, y0, t0, times(t), dt, YtObs);
                    integrate_adaptive(controlledStepper, trueSys, x0, t0, times(t), dt, XtObs);
                }
                Yt.mVec /= Y_0.rows();
                Xt.mVec /= X_0.rows();
                trukCost += calculate_cf2(Yt.mVec,Xt.mVec, weights[t]);
                Yt3Mats.push_back(Yt.mat);
                Yt3Vecs.push_back(Yt.mVec);
            }
        }else{
            Yt3Mats = readY("../data/Y", sampleSize);
            if(Yt3Mats.size() != times.size()){
                cout << "Error, number of Y_t files read in do not match the number of timesteps!" << endl;
                exit(1);
            }
            for(int i = 0; i < Yt3Mats.size(); i++){
                Yt3Vecs.push_back(moment_vector(Yt3Mats[i], nMoments));
            }
        }

        /* Compute initial wolfe weights */
        for(int t = 0; t < times.size(); ++t){
            weights[t] = ytWtMat(Yt3Mats[t], nMoments, false);
        }

        MatrixXd GBVECS = MatrixXd::Zero(nRuns, nRates + 1);
        for(int run = 0; run < nRuns; ++run){
            // make sure to reset GBMAT, POSMAT, AND PBMAT every run
            double sfi = sfe, sfc = sfp, sfs = sfg; // below are the variables being used to reiterate weights
            GBMAT = MatrixXd::Zero(0,0); // iterations of global best vectors
            MatrixXd PBMAT = MatrixXd::Zero(nParts, nRates + 1); // particle best matrix + 1 for cost component
            MatrixXd POSMAT = MatrixXd::Zero(nParts, nRates); // Position matrix as it goees through it in parallel
            
            /* Instantiate seedk aka global costs */
            double holdTheta2 = 0.1;
            struct K seed;
            seed.k = VectorXd::Zero(nRates); 
            //seed.k = testVec;
            for (int i = 0; i < nRates; i++) { 
                seed.k(i) = unifDist(gen);
            }
            seed.k(1) = holdTheta2;
            double costSeedK = 0;
            for(int t = 0; t < times.size(); t++){
                Protein_Components Xt(times(t), nMoments, X_0.rows(), X_0.cols());
                Moments_Mat_Obs XtObs(Xt);
                Nonlinear_ODE6 sys(seed);
                for (int i = 0; i < X_0.rows(); ++i) {
                    State_N c0 = convertInit(X_0.row(i));
                    Xt.index = i;
                    integrate_adaptive(controlledStepper, sys, c0, t0, times(t), dt, XtObs);
                }
                Xt.mVec /= X_0.rows();  
                cout << "XtmVec:" << Xt.mVec.transpose() << endl;
                costSeedK += calculate_cf2(Yt3Vecs[t], Xt.mVec, weights[t]);
            }

            cout << "seedk:"<< seed.k.transpose() << "| cost:" << costSeedK << endl;
            
            double gCost = costSeedK; //initialize costs and GBMAT
            // global values
            VectorXd GBVEC = seed.k;
            
            GBMAT.conservativeResize(GBMAT.rows() + 1, nRates + 1);
            for (int i = 0; i < nRates; i++) {
                GBMAT(GBMAT.rows() - 1, i) = seed.k(i);
            }
            GBMAT(GBMAT.rows() - 1, nRates) = gCost;
            double probabilityToTeleport = 3.0/4.0; 
            /* Blind PSO begins */
            cout << "PSO begins!" << endl;
            for(int step = 0; step < nSteps; ++step){
            #pragma omp parallel for 
                for(int particle = 0; particle < nParts; particle++){
                    random_device pRanDev;
                    mt19937 pGenerator(pRanDev());
                    uniform_real_distribution<double> pUnifDist(uniLowBound, uniHiBound);
                    /* instantiate all particle rate constants with unifDist */
                    if(step == 0){
                        /* temporarily assign specified k constants */
                        for(int i = 0; i < nRates; i++){
                            POSMAT(particle, i) = pUnifDist(pGenerator);
                        }
                        POSMAT(particle, 1) = holdTheta2;
                        struct K pos;
                        pos.k = VectorXd::Zero(nRates);
                        for(int i = 0; i < nRates; i++){
                            pos.k(i) = POSMAT(particle, i);
                        }
                        
                        double cost = 0;
                        for(int t = 0; t < times.size(); ++t){
                            Nonlinear_ODE6 initSys(pos);
                            Protein_Components XtPSO(times(t), nMoments, X_0.rows(), X_0.cols());
                            Moments_Mat_Obs XtObsPSO(XtPSO);
                            for(int i = 0; i < X_0.rows(); ++i){
                                //State_N c0 = gen_multi_norm_iSub();
                                State_N c0 = convertInit(X_0.row(i));
                                XtPSO.index = i;
                                integrate_adaptive(controlledStepper, initSys, c0, t0, times(t), dt, XtObsPSO);
                            }
                            XtPSO.mVec/= X_0.rows();
                            cost += calculate_cf2(Yt3Vecs[t], XtPSO.mVec, weights[t]);
                        }
                        
                        
                        /* instantiate PBMAT */
                        for(int i = 0; i < nRates; i++){
                            PBMAT(particle, i) = POSMAT(particle, i);
                        }
                        PBMAT(particle, nRates) = cost; // add cost to final column
                    }else{ 
                        /* using new rate constants, instantiate particle best values */
                        /* step into PSO */
                        double w1 = sfi * pUnifDist(pGenerator)/ sf2, w2 = sfc * pUnifDist(pGenerator) / sf2, w3 = sfs * pUnifDist(pGenerator)/ sf2;
                        double sumw = w1 + w2 + w3; //w1 = inertial, w2 = pbest, w3 = gbest
                        w1 = w1 / sumw; w2 = w2 / sumw; w3 = w3 / sumw;
                        struct K pos;
                        pos.k = VectorXd::Zero(nRates);
                        pos.k = POSMAT.row(particle);
                        VectorXd rpoint = adaptVelocity(pos.k, particle, epsi, nan, hone);
                        VectorXd PBVEC(nRates);
                        for(int i = 0; i < nRates; i++){
                            PBVEC(i) = PBMAT(particle, i);
                        }
                        
                        pos.k = w1 * rpoint + w2 * PBVEC + w3 * GBVEC; // update position of particle
                        
                        if(pUnifDist(pGenerator) < probabilityToTeleport){ // hard coded grid re-search for an adaptive component
                            pos.k(0) = pUnifDist(pGenerator);
                            pos.k(1) = pUnifDist(pGenerator);
                            pos.k(4) = pUnifDist(pGenerator);
                        }
                        // pos.k(4) = 0.05;
                        pos.k(1) = holdTheta2;
                        POSMAT.row(particle) = pos.k;
                        double cost = 0;
                        for(int t = 0; t < times.size(); t++){
                            /*solve ODEs and recompute cost */
                            Protein_Components XtPSO(times(t), nMoments, X_0.rows(), X_0.cols());
                            Moments_Mat_Obs XtObsPSO1(XtPSO);
                            Nonlinear_ODE6 stepSys(pos);
                            for(int i = 0; i < X_0.rows(); i++){
                                State_N c0 = convertInit(X_0.row(i));
                                XtPSO.index = i;
                                integrate_adaptive(controlledStepper, stepSys, c0, t0, times(t), dt, XtObsPSO1);
                            }
                            XtPSO.mVec/=X_0.rows();
                            cost += calculate_cf2(Yt3Vecs[t], XtPSO.mVec, weights[t]);
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

            for(int i = 0; i < nRates; i++){
                GBVECS(run, i) = GBVEC(i);
            }
            GBVECS(run, nRates) = gCost;
        }
        trukCost = 0;
    }
    cout << "Final Estimate:" << GBMAT.row(GBMAT.rows() - 1) << endl;
    return 0;
}