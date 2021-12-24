#include "linear.hpp"

VectorXd moment_vector(const MatrixXd &sample, int nMoments){
    VectorXd moments(nMoments);
    VectorXd mu = sample.colwise().mean();
    VectorXd variances(sample.cols());

    int nVar = sample.cols();// check to make sure if we are using variances to compute means, variances, etc. 
    if(nMoments < sample.cols()){
        nVar = 0;
    }
    // dont forget to subtract mu.
    for(int c = 0; c < nVar; c++){
        variances(c) = (sample.col(c).array() - sample.col(c).array().mean()).square().sum() / ((double) sample.col(c).array().size() - 1);
    }

    int nCross = nMoments - 2*sample.cols();
    VectorXd covariances(0);
    if(nCross > 0){
        int n = 0;
        covariances.conservativeResize(nCross);
        for (int i = 0; i < sample.cols(); i++) {
            for (int j = i + 1; j < sample.cols(); j++) {
                covariances(n) = ((sample.col(i).array() - sample.col(i).array().mean()) * (sample.col(j).array() - sample.col(j).array().mean())).sum() / ( sample.rows() - 1);
                n++;
            }
        }
    }

    // concatenate all moment vectors needed.
    for(int i = 0; i < nMoments; i++){
        if(i < sample.cols()){
            moments(i) = mu(i);
        }else if (i >= sample.cols() && i < 2 * sample.cols()){
            moments(i) = variances(i - sample.cols());
        }else if (i >= 2 * sample.cols()){
            moments(i) = covariances(i - (2 * sample.cols()));
        }
    }
    return moments;

}


/* placeholder function -> will generalize for greater proteins if we ever get a better matrix. */
MatrixXd evolutionMatrix(VectorXd &k, double tf, int nSpecies){
    MatrixXd M(nSpecies, nSpecies);
    cout << "(" << M.rows() << "," << M.cols() << ")" << endl;
	M << -k(2), k(2), 0,
		k(1), -k(1) - k(4), k(4),
		k(3), k(0), -k(0) - k(3);

	MatrixXd MT(nSpecies, nSpecies);
	MT = tf * M.transpose();
   
	MatrixXd EMT(nSpecies, nSpecies);
	EMT = MT.exp();
    return EMT;
}


VectorXd linearVelVec(const VectorXd& posK, int seed, double epsi, double nan, int hone) {
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


/*
IMPORTANT NOTE FOR FUTURE CHANGES! -> Y0 will turn into Yt as Y0 technically doesn't exist, but regardless, we'll leave it here
for pseudo-tests as this project is still in its infancy.
*/

MatrixXd linearModel(int nParts, int nSteps, int nParts2, int nSteps2, MatrixXd& X_0, MatrixXd &Y_0, int nRates, int nMoments) {
  
    /* Variables (global) */
    VectorXd times = readCsvTimeParam();
    cout << "times:" << times.transpose() << endl << endl;
    int midPt = times.size() / 2;
    double tf = times(midPt);
    double t0 = 0, dt = 1.0; // time variables
    int Npars = nRates, nSpecies = X_0.cols();
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
    // int nParts = 5000; // first part PSO
    // int nSteps = 5;
    // int nParts2 = 5; // second part PSO
    // int nSteps2 = 5000;
    int hone = 4;
    //nMoments = 2*N_SPECIES; // mean + var only!
    VectorXd wmatup(4);
    wmatup << 0.15, 0.30, 0.45, 0.60;
    double uniLowBound = 0.0, uniHiBound = 1.0;
    random_device RanDev;
    mt19937 gen(RanDev());
    uniform_real_distribution<double> unifDist(uniLowBound, uniHiBound);
    
    MatrixXd weight = MatrixXd::Identity(nMoments, nMoments);

    MatrixXd GBMAT(0, 0); // iterations of global best vectors
    MatrixXd PBMAT(nParts, Npars + 1); // particle best matrix + 1 for cost component
    MatrixXd POSMAT(nParts, Npars); // Position matrix as it goees through it in parallel
    MatrixXd X_t = MatrixXd::Zero(X_0.rows(), X_0.cols());
    MatrixXd Y_t = MatrixXd::Zero(Y_t.rows(), Y_t.cols());
    VectorXd YtmVec(nMoments);
    VectorXd XtmVec(nMoments);


    /* Solve for Y_t (mu). */
    cout << "Loading in Truk!" << endl;
    VectorXd trueK = VectorXd::Zero(Npars);
    trueK <<  0.27678200, 0.83708059, 0.44321700, 0.04244124, 0.30464502; // Bill k
    cout << "Calculating Yt!" << endl;

    Y_t = (evolutionMatrix(trueK, tf, nSpecies) * Y_0.transpose()).transpose();
    cout << "moment vec" << endl;
    YtmVec = moment_vector(Y_t, nMoments);
    /* Instantiate seedk aka global costs */
    VectorXd seed;
    seed = VectorXd::Zero(Npars); 
    for (int i = 0; i < Npars; i++) { 
        seed(i) = unifDist(gen);
    }

    double costSeedK = 0;
    X_t = (evolutionMatrix(seed, tf, nSpecies) * X_0.transpose()).transpose();
    XtmVec = moment_vector(X_t, nMoments);
    costSeedK = calculate_cf2(YtmVec, XtmVec, weight);
    cout << "seedk:"<< seed.transpose()<< "| cost:" << costSeedK << endl;
    cout << "XtmVec:" << XtmVec.transpose() << endl;
    cout << "YtmVec:" << YtmVec.transpose() << endl;
    double gCost = costSeedK; //initialize costs and GBMAT
    // global values
    VectorXd GBVEC = seed;
    
    GBMAT.conservativeResize(GBMAT.rows() + 1, Npars + 1);
    for (int i = 0; i < Npars; i++) {
        GBMAT(GBMAT.rows() - 1, i) = seed(i);
    }
    GBMAT(GBMAT.rows() - 1, Npars) = gCost;
    
    /* Blind PSO begins */
    cout << "PSO begins!" << endl;
    for(int step = 0; step < nSteps; step++){
        
    // #pragma omp parallel for 
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
                
                VectorXd pos;
                pos = VectorXd::Zero(Npars);
                for(int i = 0; i < Npars; i++){
                    pos(i) = POSMAT(particle, i);
                }
                double cost = 0;
                X_t = (evolutionMatrix(pos, tf, nSpecies) * X_0.transpose()).transpose();
                XtmVec = moment_vector(X_t, nMoments);
                cost = calculate_cf2(YtmVec, XtmVec, weight);
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
                VectorXd pos;
                pos = VectorXd::Zero(Npars);
                pos = POSMAT.row(particle);
                VectorXd rpoint = linearVelVec(pos, particle, epsi, nan, hone);
                VectorXd PBVEC(Npars);
                for(int i = 0; i < Npars; i++){
                    PBVEC(i) = PBMAT(particle, i);
                }
                pos = w1 * rpoint + w2 * PBVEC + w3 * GBVEC; // update position of particle
                POSMAT.row(particle) = pos;

                double cost = 0;
                X_t = (evolutionMatrix(pos, tf, nSpecies) * X_0.transpose()).transpose();
                XtmVec = moment_vector(X_t, nMoments);
                cost = calculate_cf2(YtmVec, XtmVec, weight);
               
                /* update gBest and pBest */
            //     #pragma omp critical
            //    {
                if(cost < PBMAT(particle, Npars)){ // particle best cost
                    for(int i = 0; i < Npars; i++){
                        PBMAT(particle, i) = pos(i);
                    }
                    PBMAT(particle, Npars) = cost;
                    if(cost < gCost){
                        gCost = cost;
                        GBVEC = pos;
                    }   
                }
            //   }
            }
        }
        GBMAT.conservativeResize(GBMAT.rows() + 1, Npars + 1); // Add to GBMAT after resizing
        for (int i = 0; i < Npars; i++) {GBMAT(GBMAT.rows() - 1, i) = GBVEC(i);}
        GBMAT(GBMAT.rows() - 1, Npars) = gCost;
        sfi = sfi - (sfe - sfg) / nSteps;   // reduce the inertial weight after each step 
        sfs = sfs + (sfe - sfg) / nSteps;
    }

    cout << "GBMAT from blind PSO:" << endl << endl;
    cout << GBMAT << endl << endl;
    cout << "truk: " << trueK.transpose() << endl;
    double dist = calculate_cf1(trueK, GBVEC);
    cout << "total difference b/w truk and final GBVEC" << dist << endl << endl; // compute difference

    /*** targeted PSO ***/
    POSMAT.conservativeResize(nParts2, Npars); // resize matrices to fit targetted PSO
    PBMAT.conservativeResize(nParts2, Npars + 1);

    cout << "targeted PSO has started!" << endl; 
    sfp = 3.0, sfg = 1.0, sfe = 6.0; // initial particle historical weight, global weight social, inertial
    sfi = sfe, sfc = sfp, sfs = sfg; // below are the variables being used to reiterate weights
    double nearby = sdbeta;
    VectorXd chkpts = wmatup * nSteps2;
    for(int step = 0; step < nSteps2; step++){
        if(step == 0 || step == chkpts(0) || step == chkpts(1) || step == chkpts(2) || step == chkpts(3)){ /* update wt   matrix || step == chkpts(0) || step == chkpts(1) || step == chkpts(2) || step == chkpts(3) */
           
            cout << "GBVEC AND COST:" << GBMAT.row(GBMAT.rows() - 1) << endl;
            nearby = squeeze * nearby;
            /* reinstantiate gCost */
            VectorXd gPos = GBVEC;
            
            double cost = 0;
            X_t = (evolutionMatrix(gPos, tf, nSpecies) * X_0.transpose()).transpose();
            XtmVec = moment_vector(X_t, nMoments);
            weight = customWtMat(Y_t, X_t, nMoments, N, false, true);
            cout << "Updated Weight Matrix!" << endl;
            cost = calculate_cf2(YtmVec, XtmVec, weight);
            
            gCost = cost;
            hone += 4;
            GBMAT.conservativeResize(GBMAT.rows() + 1, Npars + 1);
            for (int i = 0; i < Npars; i++) {GBMAT(GBMAT.rows() - 1, i) = gPos(i);}
            GBMAT(GBMAT.rows() - 1, Npars) = gCost;
        }
    // #pragma omp parallel for 
        for(int particle = 0; particle < nParts2; particle++){
            random_device pRanDev;
            mt19937 pGenerator(pRanDev());
            uniform_real_distribution<double> pUnifDist(uniLowBound, uniHiBound);
        
            if(step == 0 || step == chkpts(0) || step == chkpts(1) || step == chkpts(2) || step == chkpts(3)){
                /* reinitialize particles around global best */
                for(int edim = 0; edim < Npars; edim++){
                    int wasflipped = 0;
                    double tmean = GBVEC(edim);
                    if (GBVEC(edim) > 0.5) {
                        tmean = 1 - GBVEC(edim);
                        wasflipped = 1;
                    }
                    double myc = (1 - tmean) / tmean;
                    double alpha = myc / ((1 + myc) * (1 + myc) * (1 + myc)*nearby*nearby);
                    double beta = myc * alpha;

                    std::gamma_distribution<double> aDist(alpha, 1);
                    std::gamma_distribution<double> bDist(beta, 1);

                    double x = aDist(pGenerator);
                    double y = bDist(pGenerator);
                    double myg = x / (x + y);

                    if (wasflipped == 1) {
                        wasflipped = 0;
                        myg = 1 - myg;
                    }
                    POSMAT(particle, edim) = myg;
                }

                /* Write new POSMAT into Ks to be passed into system */
                VectorXd pos = VectorXd::Zero(Npars);
                for(int i = 0; i < Npars; i++){
                    pos(i) = POSMAT(particle, i);
                }
                //VectorXd XtPSO3 = VectorXd::Zero(nMoments);
                double cost = 0;
                X_t = (evolutionMatrix(pos, tf, nSpecies) * X_0.transpose()).transpose();
                XtmVec = moment_vector(X_t, nMoments);
                cost = calculate_cf2(YtmVec, XtmVec, weight);
                
                /* initialize PBMAT */
                for(int i = 0; i < Npars; i++){
                    PBMAT(particle, i) = POSMAT(particle, i);
                }
                PBMAT(particle, Npars) = cost; // add cost to final column
            }else{ 
                /* using new rate constants, initialize particle best values */
                /* step into PSO */
                double w1 = sfi * pUnifDist(pGenerator)/ sf2, w2 = sfc * pUnifDist(pGenerator) / sf2, w3 = sfs * pUnifDist(pGenerator)/ sf2;
                double sumw = w1 + w2 + w3; //w1 = inertial, w2 = pbest, w3 = gbest
                w1 = w1 / sumw; w2 = w2 / sumw; w3 = w3 / sumw;
                //w1 = 0.05; w2 = 0.90; w3 = 0.05;
                VectorXd pos = POSMAT.row(particle);
                VectorXd rpoint = linearVelVec(pos, particle, epsi, nan, hone);
                VectorXd PBVEC(Npars);
                for(int i = 0; i < Npars; i++){
                    PBVEC(i) = PBMAT(particle, i);
                }
                pos = w1 * rpoint + w2 * PBVEC + w3 * GBVEC; // update position of particle
                POSMAT.row(particle) = pos; // back into POSMAT
                
                double cost = 0;
                /* solve ODEs with new system and recompute cost */
                X_t = (evolutionMatrix(pos, tf, nSpecies) * X_0.transpose()).transpose();
                XtmVec = moment_vector(X_t, nMoments);
                cost = calculate_cf2(YtmVec, XtmVec, weight);
                
                /* update pBest and gBest */
                // #pragma omp critical
                // {
                if(cost < PBMAT(particle, Npars)){ // update particle best 
                    for(int i = 0; i < Npars; i++){
                        PBMAT(particle, i) = pos(i);
                    }
                    PBMAT(particle, Npars) = cost;
                    if(cost < gCost){ // update global 
                        gCost = cost;
                        GBVEC = pos;
                    }   
                }
                // }
            }
        }
        GBMAT.conservativeResize(GBMAT.rows() + 1, Npars + 1); // Add to GBMAT after each step.
        for (int i = 0; i < Npars; i++) {GBMAT(GBMAT.rows() - 1, i) = GBVEC(i);}
        GBMAT(GBMAT.rows() - 1, Npars) = gCost;

        sfi = sfi - (sfe - sfg) / nSteps2;   // reduce the inertial weight after each step 
        sfs = sfs + (sfe - sfg) / nSteps2;

    }
    cout << "GBMAT after targeted PSO:" << endl << GBMAT << endl;
    dist = calculate_cf1(trueK, GBVEC);
    cout << "total difference b/w truk and final GBVEC:" << dist << endl; // compute difference
    
	MatrixXd GBMATWithSteps(GBMAT.rows(), GBMAT.cols() + 1);
	VectorXd globalIterations(GBMAT.rows());
	for(int i = 0; i < GBMAT.rows(); i++){
		globalIterations(i) = i;
	}
	matrixToCsv(GBMATWithSteps, "GBMAT");


    return GBMAT; // just to close the program at the end.
}

