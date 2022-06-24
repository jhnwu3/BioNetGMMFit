#include "linear.hpp"

/*
    Summary:
        Takes a sample in Matrix form and converts it into a nMoment x 1 moment vector.
    Input: 
        sample - matrix of cell abundances
        nMoments - number of moments (i.e 9 moments for 3 proteins, means, var, covars)
    Output:
        moment vector - nMoments x 1 size.



 */
VectorXd momentVector(const MatrixXd &sample, int nMoments){
    VectorXd moments(nMoments);
    VectorXd mu = sample.colwise().mean();
    VectorXd variances(sample.cols());

    int nVar = sample.cols();// check to make sure if we are using variances to compute means, variances, etc. 
    if(nMoments < sample.cols()){
        nVar = 0;
    }

    // Compute sample variances
    for(int c = 0; c < nVar; c++){
        variances(c) = (sample.col(c).array() - sample.col(c).array().mean()).square().sum() / ((double) sample.col(c).array().size() - 1);
    }
    if(sample.rows() < 2){
        variances = VectorXd::Zero(sample.cols());
    }

    // again only compute covariances, if number of moments allow for it
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

    // Now after all computations, add to moment vector
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


/* 
    Summary:
        Takes the interaction matrix and "evolves" it using the specified interaction matrix in system.cpp
    Input:
        k - parameter vector that is being "estimated" (i.e rate constants)
        tf - time of evolution
        nSpecies - number of species in the model


 */
MatrixXd evolutionMatrix(VectorXd &k, double tf, int nSpecies){
    MatrixXd M(nSpecies, nSpecies);
    M = interactionMatrix(nSpecies, k);
	
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
    if(seed > 0){generator.seed(seed);}
    vector<int> rand;
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
            pos -= epsi;
        }else if (pos < nan) {
            pos += epsi;
        }
        double alpha = hone * pos; // Component specific
        double beta = hone - alpha; // pos specific
        std::gamma_distribution<double> aDist(alpha, 1); // beta distribution consisting of gamma distributions
        std::gamma_distribution<double> bDist(beta, 1);

        double x = aDist(generator);
        double y = bDist(generator);

        rPoint(px) = (x / (x + y)); 
    }
    
    return rPoint;
}


/*
    Summary:
        Takes program parameters and computes using the described linear model in system.cpp and computes a parameter estimate.   
    Input:
        nParts - Number of Blind PSO particles
        nSteps - Number of Blind PSO steps
        nParts2 - Number of Targeted PSO particles
        nSteps2 - Number of Targeted PSO step counts
        X_0 - Control Matrix loaded in
        nRates - number of rates that are being simulated
        nMoments - number of moments, i.e 9 moments for both means + variances + covariances for 3 proteins
        times - Vector of time steps that you'll be plugging in
        simulateYt - integer/boolean for simulating Yt so you only need to plug in Y_0 (used to test ground truths)
    Output:
        GBMAT - Global Best Matrix with each row being a new "global" best parameter estimate in the PSO, the final row is the final estimate.
*/

MatrixXd linearModel(int nParts, int nSteps, int nParts2, int nSteps2, MatrixXd& X_0, int nRates, int nMoments, const VectorXd &times, int simulateYt, int useInverse, int argc, char ** argv, int rngSeed) {
  
    int midPt = times.size() / 2; // take only the midpoint of all listed time points for now for evolution
    double tf = times(midPt);
    double t0 = 0, dt = 1.0; // time variables
    int Npars = nRates, nSpecies = X_0.cols();
    double squeeze = 0.96, sdbeta = 0.05; 
    double boundary = 0.001;
    bool willInvert = useInverse == 1;
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
    int hone = 4;
    VectorXd wmatup(4);
    wmatup << 0.15, 0.30, 0.45, 0.60;
    double low = 0.0, high = 1.0;
    random_device RanDev;
    mt19937 gen(RanDev());
    if(rngSeed > 0){
        gen.seed(rngSeed);
    }
    uniform_real_distribution<double> unifDist(low, high);
    MatrixXd weight = MatrixXd::Identity(nMoments, nMoments);
    MatrixXd GBMAT(0, 0); // iterations of global best vectors
    MatrixXd PBMAT(nParts, Npars + 1); // particle best matrix + 1 for cost component
    MatrixXd POSMAT(nParts, Npars); // Position matrix as it goees through it in parallel
    MatrixXd Y_t = MatrixXd::Zero(Y_t.rows(), Y_t.cols());
    VectorXd YtmVec(nMoments);
    /* Solve or load Y_t  */
    VectorXd trueK = readRates(nRates, getTrueRatesPath(argc, argv)); 
    if(simulateYt == 1){
        MatrixXd Y_0 = readY("data/Y")[0];
        Y_0 = filterZeros(Y_0);
        cout << "Simulating Yt!" << endl;
        cout << "with evolution matrix:" << endl << evolutionMatrix(trueK, tf, nSpecies) << endl;
        Y_t = (evolutionMatrix(trueK, tf, nSpecies) * Y_0.transpose()).transpose();
        YtmVec = momentVector(Y_t, nMoments);
    }else{
        Y_t = readY("data/Y")[0];
        Y_t = filterZeros(Y_t);
        YtmVec = momentVector(Y_t, nMoments);
    }
    weight = wolfWtMat(Y_t, nMoments, willInvert); // wolf weights

    /* Initialize seedk aka global costs */
    VectorXd seed;
    seed = VectorXd::Zero(Npars); 
    for (int i = 0; i < Npars; i++) { 
        seed(i) = unifDist(gen);
    }
    cout << "Seeded at " << seed.transpose() << endl;
    double costSeedK = 0;
    MatrixXd seedXt = (evolutionMatrix(seed, tf, nSpecies) * X_0.transpose()).transpose();
    VectorXd seedMoms = momentVector(seedXt, nMoments);
    costSeedK = costFunction(YtmVec, seedMoms, weight);

    /* Initialize the start of the global best matrix */
    double gCost = costSeedK; 
    VectorXd GBVEC = seed;
    GBMAT.conservativeResize(GBMAT.rows() + 1, Npars + 1);
    for (int i = 0; i < Npars; i++) {GBMAT(GBMAT.rows() - 1, i) = seed(i);}
    GBMAT(GBMAT.rows() - 1, Npars) = gCost;
    
    /* Blind PSO begins */
    cout << "PSO has begun!" << endl;
    for(int step = 0; step < nSteps; step++){
    #pragma omp parallel for 
        for(int particle = 0; particle < nParts; particle++){
            random_device pRanDev;
            mt19937 pGenerator(pRanDev());
            uniform_real_distribution<double> pUnifDist(low, high);
            int pSeed = -1;
            if(rngSeed > 0){
                pSeed = particle + step + rngSeed;
                pGenerator.seed(pSeed);
            }
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
                MatrixXd X_t = (evolutionMatrix(pos, tf, nSpecies) * X_0.transpose()).transpose();
                cost = costFunction(YtmVec, momentVector(X_t, nMoments), weight);
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
                VectorXd rpoint = linearVelVec(pos, pSeed, epsi, nan, hone);
                VectorXd PBVEC(Npars);
                for(int i = 0; i < Npars; i++){PBVEC(i) = PBMAT(particle, i);}
                pos = w1 * rpoint + w2 * PBVEC + w3 * GBVEC; // update position of particle
                POSMAT.row(particle) = pos;

                double cost = 0;
                MatrixXd X_t = (evolutionMatrix(pos, tf, nSpecies) * X_0.transpose()).transpose();
                cost = costFunction(YtmVec, momentVector(X_t, nMoments), weight);
               
                /* update gBest and pBest */
                #pragma omp critical
               {
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
              }
            }
        }
        GBMAT.conservativeResize(GBMAT.rows() + 1, Npars + 1); // Add to GBMAT after resizing
        for (int i = 0; i < Npars; i++) {GBMAT(GBMAT.rows() - 1, i) = GBVEC(i);}
        GBMAT(GBMAT.rows() - 1, Npars) = gCost;
        sfi = sfi - (sfe - sfg) / nSteps;   // reduce the inertial weight after each step 
        sfs = sfs + (sfe - sfg) / nSteps;
    }

    /*** targeted PSO ***/
    POSMAT.conservativeResize(nParts2, Npars); // resize matrices to fit targetted PSO
    PBMAT.conservativeResize(nParts2, Npars + 1);

    cout << "Targeted PSO has started!" << endl; 
    sfp = 3.0, sfg = 1.0, sfe = 6.0; // initial particle historical weight, global weight social, inertial
    sfi = sfe, sfc = sfp, sfs = sfg; // below are the variables being used to reiterate weights
    double nearby = sdbeta;
    VectorXd chkpts = wmatup * nSteps2;
    for(int step = 0; step < nSteps2; step++){
        if(step == 0 || step == chkpts(0) || step == chkpts(1) || step == chkpts(2) || step == chkpts(3)){ /* update wt   matrix || step == chkpts(0) || step == chkpts(1) || step == chkpts(2) || step == chkpts(3) */
            nearby = squeeze * nearby;
            /* reinstantiate gCost */
            VectorXd gPos = GBVEC;
            double cost = 0;
            MatrixXd X_t = (evolutionMatrix(gPos, tf, nSpecies) * X_0.transpose()).transpose();
            weight = dasWtMat(Y_t, X_t, nMoments, N, willInvert);
            cost = costFunction(YtmVec, momentVector(X_t, nMoments), weight);
            gCost = cost;
            hone += 4;
            GBMAT.conservativeResize(GBMAT.rows() + 1, Npars + 1);
            for (int i = 0; i < Npars; i++) {GBMAT(GBMAT.rows() - 1, i) = gPos(i);}
            GBMAT(GBMAT.rows() - 1, Npars) = gCost;
        }
    #pragma omp parallel for 
        for(int particle = 0; particle < nParts2; particle++){
            random_device pRanDev;
            mt19937 pGenerator(pRanDev());
            uniform_real_distribution<double> pUnifDist(low, high);
            int pSeed = -1;
            if(rngSeed > 0){
                pSeed = particle + step + rngSeed;
                pGenerator.seed(pSeed);
            }
            if(step == 0 || step == chkpts(0) || step == chkpts(1) || step == chkpts(2) || step == chkpts(3)){
                /* reinitialize particles around global best */
                for(int edim = 0; edim < Npars; edim++){
                    int wasflipped = 0;
                    double tmean = GBVEC(edim);
                    /* If we hit a boundary (close to 1), flip the mean back towards the center of the beta distribution */
                    if (GBVEC(edim) > 0.5) {
                        tmean = 1 - GBVEC(edim);
                        wasflipped = 1;
                    }
                    /* Compute Specific parameters for beta dist */
                    double myc = (1 - tmean) / tmean;
                    double alpha = myc / ((1 + myc) * (1 + myc) * (1 + myc)*nearby*nearby);
                    double beta = myc * alpha;

                    std::gamma_distribution<double> aDist(alpha, 1);
                    std::gamma_distribution<double> bDist(beta, 1);
                    /* analytic solution for beta distribution */
                    double x = aDist(pGenerator);
                    double y = bDist(pGenerator);
                    double myg = x / (x + y);

                    if (wasflipped == 1) {
                        wasflipped = 0;
                        myg = 1 - myg;
                    }
                    /* set new particle starting position for targeted PSO */
                    POSMAT(particle, edim) = myg;
                }

                /* Write new POSMAT into Ks to be passed into system */
                VectorXd pos = VectorXd::Zero(Npars);
                for(int i = 0; i < Npars; i++){
                    pos(i) = POSMAT(particle, i);
                }
                //VectorXd XtPSO3 = VectorXd::Zero(nMoments);
                double cost = 0;
                MatrixXd X_t = (evolutionMatrix(pos, tf, nSpecies) * X_0.transpose()).transpose();
                cost = costFunction(YtmVec, momentVector(X_t, nMoments), weight);
                
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
                VectorXd rpoint = linearVelVec(pos, pSeed, epsi, nan, hone);
                VectorXd PBVEC(Npars);
                for(int i = 0; i < Npars; i++){
                    PBVEC(i) = PBMAT(particle, i);
                }
                pos = w1 * rpoint + w2 * PBVEC + w3 * GBVEC; // update position of particle
                POSMAT.row(particle) = pos; // back into POSMAT
                
                double cost = 0;
                /* solve ODEs with new system and recompute cost */
                MatrixXd X_t = (evolutionMatrix(pos, tf, nSpecies) * X_0.transpose()).transpose();
                cost = costFunction(YtmVec, momentVector(X_t, nMoments), weight);
                
                /* update pBest and gBest */
                #pragma omp critical
                {
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
                }
            }
        }
        GBMAT.conservativeResize(GBMAT.rows() + 1, Npars + 1); // Add to GBMAT after each step.
        for (int i = 0; i < Npars; i++) {GBMAT(GBMAT.rows() - 1, i) = GBVEC(i);}
        GBMAT(GBMAT.rows() - 1, Npars) = gCost;

        sfi = sfi - (sfe - sfg) / nSteps2;   // reduce the inertial weight after each step 
        sfs = sfs + (sfe - sfg) / nSteps2;
    }
    
    if(simulateYt == 1){
        cout << "Simulation Ground Truth:" << trueK.transpose() << endl;
    }
    return GBMAT; // just to close the program at the end.
}

