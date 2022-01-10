#include "nonlinear.hpp"
State_N convertInit(const VectorXd &v1){
    vector<double> v2;
    v2.resize(v1.size());
    VectorXd::Map(&v2[0], v1.size()) = v1;
    return v2;
}
VectorXd adaptVelocity(const VectorXd& posK, int seed, double epsi, double nan, int hone) {
    
    VectorXd rPoint;
    rPoint = posK;
    std::random_device rand_dev;
    std::mt19937 generator(rand_dev());
    vector<int> rand;
    uniform_real_distribution<double> unifDist(0.0, 1.0);
    for (int i = 0; i < posK.size(); i++) {
        rand.push_back(i);
    }
    shuffle(rand.begin(), rand.end(), generator); // shuffle indices as well as possible. 
    int ncomp = rand.at(0);
    VectorXd wcomp(ncomp);
    shuffle(rand.begin(), rand.end(), generator);
    for (int i = 0; i < ncomp; i++) {
        wcomp(i) = rand.at(i);
    }
    
    for (int smart = 0; smart < ncomp; smart++) {
        int px = wcomp(smart);
        double pos = rPoint(px);
        if (pos > 1.0 - nan) {
            cout << "overflow!" << endl;
            pos -= epsi;
        }else if (pos < nan) {
            cout << "underflow!"<< pos << endl;
            pos += epsi;
            cout << "pos" << posK.transpose() << endl; 
        }
        double alpha = hone * pos; // Component specific
        double beta = hone - alpha; // pos specific
    // cout << "alpha:" << alpha << "beta:" << beta << endl;
        std::gamma_distribution<double> aDist(alpha, 1); // beta distribution consisting of gamma distributions
        std::gamma_distribution<double> bDist(beta, 1);

        double x = aDist(generator);
        double y = bDist(generator);

        rPoint(px) = (x / (x + y)); 
    }
    
    return rPoint;
}

MatrixXd nonlinearModel(int nParts, int nSteps, int nParts2, int nSteps2, const MatrixXd& X_0, const MatrixXd& Y_0, int nRates, int nRuns, int nMoments){
    auto t1 = std::chrono::high_resolution_clock::now();
    /*---------------------- Setup ------------------------ */
    VectorXd times = readCsvTimeParam();
    /* Variables (global) */
    double t0 = 0, dt = 1.0; // time variables
    int nTimeSteps = times.size();
    int Npars = nRates;
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
    int N = X_0.rows();
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

    for(int i = 0; i < nTimeSteps; i++){
        weights.push_back(MatrixXd::Identity(nMoments, nMoments));
    }
    
    cout << "Using two part PSO " << "Sample Size:" << N << " with:" << nMoments << " moments." << endl;
    cout << "Using Times:" << times.transpose() << endl;
    cout << "Bounds for Uniform Distribution (" << uniLowBound << "," << uniHiBound << ")"<< endl;
    cout << "Blind PSO --> nParts:" << nParts << " Nsteps:" << nSteps << endl;
    cout << "Targeted PSO --> nParts:" <<  nParts2 << " Nsteps:" << nSteps2 << endl;
    cout << "sdbeta:" << sdbeta << endl;
    // cout << "wt:" << endl << wt << endl;

    MatrixXd GBMAT(0, 0); // iterations of global best vectors
    MatrixXd PBMAT(nParts, Npars + 1); // particle best matrix + 1 for cost component
    MatrixXd POSMAT(nParts, Npars); // Position matrix as it goees through it in parallel

    /* Solve for Y_t (mu). */
    cout << "Loading in Truk!" << endl;
    struct K tru;
    tru.k = VectorXd::Zero(Npars);
    tru.k << 0.1, 0.1, 0.95, 0.17, 0.05, 0.18;

    cout << "Calculating Yt!" << endl;
    vector<MatrixXd> Yt3Mats;
    vector<VectorXd> Yt3Vecs;
    vector<VectorXd> Xt3Vecs;
    Controlled_RK_Stepper_N controlledStepper;
    double trukCost = 0;
    for(int t = 0; t < nTimeSteps; t++){
        Nonlinear_ODE6 trueSys(tru);
        Protein_Components Yt(times(t), nMoments, N, X_0.cols());
        Protein_Components Xt(times(t), nMoments, N, X_0.cols());
        Moments_Mat_Obs YtObs(Yt);
        Moments_Mat_Obs XtObs(Xt);
        for (int i = 0; i < N; ++i) {
            //State_N c0 = gen_multi_norm_iSub(); // Y_0 is simulated using norm dist.
            State_N y0 = convertInit(Y_0.row(i));
            State_N x0 = convertInit(X_0.row(i));
            Yt.index = i;
            Xt.index = i;
            integrate_adaptive(controlledStepper, trueSys, y0, t0, times(t), dt, YtObs);
            integrate_adaptive(controlledStepper, trueSys, x0, t0, times(t), dt, XtObs);
        }
        Yt.mVec /= N;
        Xt.mVec /= N;
        trukCost += calculate_cf2(Yt.mVec,Xt.mVec, weights[t]);
        Xt3Vecs.push_back(Xt.mVec);
        Yt3Mats.push_back(Yt.mat);
        Yt3Vecs.push_back(Yt.mVec);
    }
    cout << "truk cost:"<< trukCost << endl;


    /* Compute initial wolfe weights */
    for(int t = 0; t < nTimeSteps; ++t){
        weights[t] = ytWtMat(Yt3Mats[t], nMoments, false);
    }

    MatrixXd GBVECS = MatrixXd::Zero(nRuns, Npars + 1);
    for(int run = 0; run < nRuns; ++run){
        // make sure to reset GBMAT, POSMAT, AND PBMAT EVERY RUN!
        double sfi = sfe, sfc = sfp, sfs = sfg; // below are the variables being used to reiterate weights
        MatrixXd GBMAT = MatrixXd::Zero(0,0); // iterations of global best vectors
        MatrixXd PBMAT = MatrixXd::Zero(nParts, Npars + 1); // particle best matrix + 1 for cost component
        MatrixXd POSMAT = MatrixXd::Zero(nParts, Npars); // Position matrix as it goees through it in parallel
        
        /* Instantiate seedk aka global costs */
        double holdTheta2 = 0.1;
        struct K seed;
        seed.k = VectorXd::Zero(Npars); 
        //seed.k = testVec;
        for (int i = 0; i < Npars; i++) { 
            seed.k(i) = unifDist(gen);
        }
        // seed.k(4) = tru.k(4);
        seed.k(1) = holdTheta2;
        // seed.k <<    0.094531 , 0.99 , 0.938388 , 0.170400 , 0.0517104 , 0.180564;
        // holdTheta2 = seed.k(1);
        // seed.k = tru.k;
        double costSeedK = 0;
        for(int t = 0; t < nTimeSteps; t++){
            Protein_Components Xt(times(t), nMoments, N, X_0.cols());
            Moments_Mat_Obs XtObs(Xt);
            Nonlinear_ODE6 sys(seed);
            for (int i = 0; i < N; ++i) {
                //State_N c0 = gen_multi_norm_iSub();
                State_N c0 = convertInit(X_0.row(i));
                Xt.index = i;
                integrate_adaptive(controlledStepper, sys, c0, t0, times(t), dt, XtObs);
            }
            Xt.mVec /= N;  
            cout << "XtmVec:" << Xt.mVec.transpose() << endl;
            costSeedK += calculate_cf2(Yt3Vecs[t], Xt.mVec, weights[t]);
        }

        cout << "seedk:"<< seed.k.transpose() << "| cost:" << costSeedK << endl;
        
        double gCost = costSeedK; //initialize costs and GBMAT
        // global values
        VectorXd GBVEC = seed.k;
        
        GBMAT.conservativeResize(GBMAT.rows() + 1, Npars + 1);
        for (int i = 0; i < Npars; i++) {
            GBMAT(GBMAT.rows() - 1, i) = seed.k(i);
        }
        GBMAT(GBMAT.rows() - 1, Npars) = gCost;
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
                    for(int i = 0; i < Npars; i++){
                        POSMAT(particle, i) = pUnifDist(pGenerator);
                    }
                    // POSMAT(particle, 4) = 0.05;
                    POSMAT(particle, 1) = holdTheta2;
                    // POSMAT.row(particle) = seed.k;
                    struct K pos;
                    pos.k = VectorXd::Zero(Npars);
                    for(int i = 0; i < Npars; i++){
                        pos.k(i) = POSMAT(particle, i);
                    }
                    
                    double cost = 0;
                    for(int t = 0; t < nTimeSteps; t++){
                        Nonlinear_ODE6 initSys(pos);
                        Protein_Components XtPSO(times(t), nMoments, N, X_0.cols());
                        Moments_Mat_Obs XtObsPSO(XtPSO);
                        for(int i = 0; i < N; i++){
                            //State_N c0 = gen_multi_norm_iSub();
                            State_N c0 = convertInit(X_0.row(i));
                            XtPSO.index = i;
                            integrate_adaptive(controlledStepper, initSys, c0, t0, times(t), dt, XtObsPSO);
                        }
                        XtPSO.mVec/=N;
                        cost += calculate_cf2(Yt3Vecs[t], XtPSO.mVec, weights[t]);
                    }
                    
                    
                    /* instantiate PBMAT */
                    for(int i = 0; i < Npars; i++){
                        PBMAT(particle, i) = POSMAT(particle, i);
                    }
                    PBMAT(particle, Npars) = cost; // add cost to final column
                }else{ 
                    /* using new rate constants, instantiate particle best values */
                    /* step into PSO */
                    double w1 = sfi * pUnifDist(pGenerator)/ sf2, w2 = sfc * pUnifDist(pGenerator) / sf2, w3 = sfs * pUnifDist(pGenerator)/ sf2;
                    double sumw = w1 + w2 + w3; //w1 = inertial, w2 = pbest, w3 = gbest
                    w1 = w1 / sumw; w2 = w2 / sumw; w3 = w3 / sumw;
                    //w1 = 0.05; w2 = 0.90; w3 = 0.05;
                    struct K pos;
                    pos.k = VectorXd::Zero(Npars);
                    pos.k = POSMAT.row(particle);
                    VectorXd rpoint = adaptVelocity(pos.k, particle, epsi, nan, hone);
                    VectorXd PBVEC(Npars);
                    for(int i = 0; i < Npars; i++){
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
                    for(int t = 0; t < nTimeSteps; t++){
                        /*solve ODEs and recompute cost */
                        Protein_Components XtPSO(times(t), nMoments, N, X_0.cols());
                        Moments_Mat_Obs XtObsPSO1(XtPSO);
                        Nonlinear_ODE6 stepSys(pos);
                        for(int i = 0; i < N; i++){
                            State_N c0 = convertInit(X_0.row(i));
                            XtPSO.index = i;
                            integrate_adaptive(controlledStepper, stepSys, c0, t0, times(t), dt, XtObsPSO1);
                        }
                        XtPSO.mVec/=N;
                        cost += calculate_cf2(Yt3Vecs[t], XtPSO.mVec, weights[t]);
                    }
                
                    /* update gBest and pBest */
                #pragma omp critical
                {
                    if(cost < PBMAT(particle, Npars)){ // particle best cost
                        for(int i = 0; i < Npars; i++){
                            PBMAT(particle, i) = pos.k(i);
                        }
                        PBMAT(particle, Npars) = cost;
                        if(cost < gCost){
                            gCost = cost;
                            GBVEC = pos.k;
                        }   
                    }
                }
                }
            }
            GBMAT.conservativeResize(GBMAT.rows() + 1, Npars + 1); // Add to GBMAT after resizing
            for (int i = 0; i < Npars; i++) {GBMAT(GBMAT.rows() - 1, i) = GBVEC(i);}
            GBMAT(GBMAT.rows() - 1, Npars) = gCost;
            sfi = sfi - (sfe - sfg) / nSteps;   // reduce the inertial weight after each step 
            sfs = sfs + (sfe - sfg) / nSteps;
            cout << "current:" << GBVEC.transpose()<<" "<< gCost << endl;
        }

        cout << "GBMAT from blind PSO:" << endl << endl;
        cout << GBMAT << endl << endl;
        cout << "truk: " << tru.k.transpose() << endl;
        double dist = calculate_cf1(tru.k, GBVEC);
        cout << "total difference b/w truk and final GBVEC" << dist << endl << endl; // compute difference
        auto tB = std::chrono::high_resolution_clock::now();
        auto bDuration = std::chrono::duration_cast<std::chrono::seconds>(tB - t1).count();
        cout << "blind PSO FINISHED RUNNING IN " << bDuration << " s TIME!" << endl;
    
        for(int i = 0; i < Npars; i++){
            GBVECS(run, i) = GBVEC(i);
        }
        GBVECS(run, Npars) = gCost;
    }
    trukCost = 0;
    for(int t = 0; t < nTimeSteps; t++){
        trukCost += calculate_cf2(Yt3Vecs[t], Xt3Vecs[t], weights[t]);
    }

    cout << "truk: " << tru.k.transpose() << " with trukCost with new weights:" << trukCost << endl;
    
    auto t2 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count();
    cout << "CODE FINISHED RUNNING IN " << duration << " s TIME!" << endl;

    return GBMAT; // just to close the program at the end.
}

