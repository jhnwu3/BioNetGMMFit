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
#include "sbml.hpp"
#include "param.hpp"
#include "cli.hpp"
int main(int argc, char** argv){
    if(helpCall(argc, argv)){return EXIT_SUCCESS;}

    auto t1 = std::chrono::high_resolution_clock::now();
    cout << "Program Begin:" << endl;
    cout << "** Please Make Sure That All Inputted Files are in the UNIX Line Formatting to Prevent Bugs! **" << endl;
    const string bngl = ".bngl"; // suffixes for sbml/bngl file types
    const string sbml = "_sbml.xml";
    
    /* Input Parameters for Program */
    Parameters parameters = Parameters(getConfigPath(argc, argv));

    /* Important! Max Thread Count Test */
    omp_set_num_threads(parameters.nThreads);

    /* Run an update for bngl */
    string modelPath = getModelPath(argc, argv);
    string baseModelFile = modelPath.substr(modelPath.find_last_of("/\\") + 1);
    std::string::size_type const p(baseModelFile.find_last_of('.'));
    std::string file_without_extension = baseModelFile.substr(0, p);
    string sbmlModel = "sbml/"+ file_without_extension + sbml;
    const string bnglCall = "bionetgen run -i" + modelPath + " -o sbml";
    if(parameters.useSBML > 0 && system(bnglCall.c_str()) < 0){
        cout << "Error Running BioNetGen -> Make sure you have installed it through pip install bionetgen or check program permissions!" << endl;
    }

    /* Read time steps*/
    VectorXd times = readCsvTimeParam(getTimeStepsPath(argc, argv));
    if(times.size() < 2){
        cout << "Error! Unable to read in timesteps properly, need two or more time steps!" << endl;
        exit(1);
    }

    /* Read First Matrix */
    MatrixXd X_0;
    MatrixXd ogX_0; // copy of initial X values if bootstrap is used.
    X_0 = readX(getXPath(argc, argv));
    X_0 = filterZeros(X_0);
    ogX_0 = X_0;
    cout << "After removing all negative rows, X has " << X_0.rows() << " rows." << endl;
    cout << "If dimensions are unexpected of input data, please make sure to check csv files contains all possible values in each row/column." << endl;
    /* Decide Number of Moments to Use in Estimation */
    int nMoments = (X_0.cols() * (X_0.cols() + 3)) / 2;
    if(parameters.useOnlySecMom){  // these will be added to the options sheet later.
        nMoments = 2 * X_0.cols();
    }
    if(parameters.useOnlyFirstMom){
        nMoments = X_0.cols();
    }
    parameters.printParameters(nMoments, times);

    /* RoadRunner Configuration and Simulation Variables */
    RoadRunner r = RoadRunner(sbmlModel);
    vector<int> specifiedProteins;
    vector<string> parameterNames = r.getGlobalParameterIds(); // parameter names check
    vector<string> speciesNames =  getSpeciesNames(sbmlModel);
    cout << "--------------------------------------------------------" << endl;
    if(r.getNumberOfIndependentSpecies() +  r.getNumberOfDependentSpecies() != X_0.cols() && parameters.useSBML > 0){
        if(X_0.cols() > r.getNumberOfIndependentSpecies() +  r.getNumberOfDependentSpecies() ){
            cout << "Error Too Many Species/Columns in X.csv file! Please remove some before continuing!" << endl;
            cout << "Expected:" <<  r.getNumberOfIndependentSpecies() +  r.getNumberOfDependentSpecies() << " Got:" << X_0.cols() << endl;
            return EXIT_FAILURE;
        }
        cout << "Number of Species Defined in BNGL does not match number of columns in data files! Now listing all species in system and respective indices in order!" << endl;
        cout << "Note: User can supply a \"./CyGMM -p protein_observed.txt \" to specify explicitly which proteins are observed in data. Please make sure names are in order from top to bottom matching left to right in data csv file." << endl;
        for(int i = 0; i < speciesNames.size(); i++){
            cout << "(" << i << ") " << speciesNames[i] << endl; 
        }
        if(proPathExists(argc, argv)){
            specifiedProteins = specifySpeciesFromProteinsList(getProPath(argc,argv), speciesNames, X_0.cols());
            if(specifiedProteins.size() < 1){
                cout << "Error, no proteins specified!" << endl;
                exit(EXIT_FAILURE);
            }
            cout << "From Above List of Indexed Species, We are using..." << endl;
            for(int i = 0; i < specifiedProteins.size(); i++){
                cout << "(" <<specifiedProteins[i] <<") "<< specifiedProteins[i] <<endl;
            }
        }else{
            cout << "No Proteins Specified Using Argument \"-p Proteins.txt\", Thus Using First " << X_0.cols() << " Species Listed" << endl;
        }
        cout << "Proteins Not Observed Will Default to Initial Values Defined in .BNGL File" << endl;
        // check if some txt -p file has been used and use that otherwise have user manually select.
    }else if(parameters.useSBML > 0){
        cout << "------- Matching Columns of X Data files to Ids -------" << endl;
        for(int i = 0; i < speciesNames.size(); i++){
            cout << speciesNames[i] << " to column:"<< i << " with first value:" << X_0(0,i) << endl;
        }
    }
    cout << "--------------------------------------------------------" << endl;
    if(parameters.useDet > 0){
        r.setIntegrator("cvode");
    }else{
        r.setIntegrator("gillespie");
    }
    SimulateOptions opt;
    opt.steps = parameters.odeSteps;
    double theta[parameters.nRates];// static array to be constantly used with road runner model parameters.

    MatrixXd GBMAT;
    MatrixXd GBVECS = MatrixXd::Zero(parameters.nRuns, parameters.nRates + 1);
    if(parameters.useLinear > 0){
        for(int r = 0; r < parameters.nRuns; ++r){
            GBMAT = linearModel(parameters.nParts, parameters.nSteps, parameters.nParts2, parameters.nSteps2, X_0, parameters.nRates, nMoments, times, parameters.simulateYt, parameters.useInverse, argc, argv, parameters.seed);
            GBVECS.row(r) = GBMAT.row(GBMAT.rows() - 1);
        }
    }else{
        /*---------------------- Nonlinear Setup PSO ------------------------ */
        double dt = 1.0; // nonlinear time evolution variables
        /* Explicit Boundary Parameters */
        double squeeze = 0.500, sdbeta = 0.10; // how much to shrink PSO search over time (how much variability each position is iterated upon)
        double boundary = 0.001;
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
        vector<MatrixXd> weights;
        MatrixXd PBMAT(parameters.nParts, parameters.nRates + 1); // particle best matrix + 1 for cost component
        MatrixXd POSMAT(parameters.nParts, parameters.nRates); // Position matrix as it goees through it in parallel

        /* RNG seeding just in case */
        random_device RanDev;
        mt19937 gen(RanDev());
        uniform_real_distribution<double> unifDist(low, high);
        if(parameters.seed > 0){
            gen.seed(parameters.seed);
        }
        /* Solve for Y_t (mu). */
        VectorXd tru;
        vector<MatrixXd> Yt3Mats;
        vector<MatrixXd> ogYt3Mats;
        vector<VectorXd> Yt3Vecs;
        Controlled_RK_Stepper_N controlledStepper;
        double trukCost = 0;
        if(parameters.simulateYt > 0){
            cout << "------ SIMULATING YT! ------" << endl;
            tru = readRates(parameters.nRates, getTrueRatesPath(argc, argv));
            cout << "Read in Rates:" << tru.transpose() << endl;
            MatrixXd Y_0 = readY(getYPath(argc, argv))[0];
            Y_0 = filterZeros(Y_0);
            cout << "Note: We will only be using the first Yt file read in for this simulation!" << endl;
            cout << "After removing all negative rows, Y has " << Y_0.rows() << " rows." << endl;
            cout << "Time Point \t\t Moments" << endl;
            for(int t = 1; t < times.size(); t++){ // start at t1, because t0 is now in the vector
                if(parameters.useSBML > 0){
                    MatrixXd YtMat = MatrixXd::Zero(Y_0.rows(), Y_0.cols());
                    opt.start = times(0);
                    opt.duration = times(t);
                    for(int i = 0; i < tru.size(); ++i){
                        theta[i] = tru(i);  
                    }
                    r.getModel()->setGlobalParameterValues(tru.size(), 0, theta); // set new global parameter values here.
                    for(int i = 0; i < Y_0.rows(); ++i){
                        if(specifiedProteins.size() > 0){
                            vector<double> init = r.getFloatingSpeciesInitialConcentrations();
                            for(int p = 0; p < specifiedProteins.size(); p++){
                                init[specifiedProteins[p]] = Y_0(i,p);
                            }
                            r.changeInitialConditions(init);
                        }else{
                            r.changeInitialConditions(convertInit(Y_0.row(i)));
                        }
                        const DoubleMatrix res = *r.simulate(&opt);
                        for(int j = 0; j < Y_0.cols(); ++j){
                            YtMat(i,j) = res[res.numRows() - 1][j + 1];
                        }
                    }
                    Yt3Vecs.push_back(momentVector(YtMat, nMoments));
                    Yt3Mats.push_back(YtMat);
                }else{
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
                cout << t << " "<< Yt3Vecs[t-1].transpose() << endl;
            }
            cout << "--------------------------------------------------------" << endl;
        }else{
            Yt3Mats = readY(getYPath(argc, argv));
            if(Yt3Mats.size() + 1 != times.size()){
                cout << "Error, number of Y_t files read in do not match the number of timesteps!" << endl;
                exit(1);
            }
            // filter all zeroes and compute moments vectors for cost calcs
            for(int i = 0; i < Yt3Mats.size(); i++){
                Yt3Mats[i] = filterZeros(Yt3Mats[i]);
                cout << "After removing all negative rows, Y"<< i << " has " << Yt3Mats[i].rows() << " rows." << endl;
                Yt3Vecs.push_back(momentVector(Yt3Mats[i], nMoments));
            }
            ogYt3Mats = Yt3Mats;
        }
        cout << "Yt Means For First Time Step:" << Yt3Mats[0].colwise().mean() << endl;
        cout << "Computing Weight Matrices!" << endl;
        /* Compute initial wolfe weights */
        if (Yt3Mats[0].cols() != X_0.cols()){
            cout << "Error, mismatch in number of species/columns between X and Y!" << endl;
            cout << "X:" << X_0.cols() << " Y:" << Yt3Mats[0].cols() << endl;
            exit(1);
        }
        for(int y = 0; y < Yt3Mats.size(); ++y){ 
            weights.push_back(wolfWtMat(Yt3Mats[y], nMoments, false));
        }
        cout << weights[0] << endl;
        for(int run = 0; run < parameters.nRuns; ++run){ // for multiple runs aka bootstrapping (for now)
            VectorXd nestedHolds = VectorXd::Zero(parameters.nRates);
            if (run > 0 && parameters.bootstrap > 0){
                for(int y = 0; y < Yt3Mats.size(); ++y){ 
                    weights[y] = wolfWtMat(Yt3Mats[y], nMoments, false);
                }
            }
            for(int ne = 0; ne < parameters.nest; ne++){ // now we will include the ability to nest hypercubes. Every nested hypercube doubles the hypercube width.
                // make sure to reset GBMAT, POSMAT, AND PBMAT every run
                double sfi = sfe, sfc = sfp, sfs = sfg; // below are the variables being used to reiterate weights
                GBMAT = MatrixXd::Zero(0,0); // iterations of global best vectors
                MatrixXd PBMAT = MatrixXd::Zero(parameters.nParts, parameters.nRates + 1); // particle best matrix + 1 for cost component
                MatrixXd POSMAT = MatrixXd::Zero(parameters.nParts, parameters.nRates); // Position matrix as it goees through it in parallel
                
                /* Initialize Global Best  */
                VectorXd seed = VectorXd::Zero(parameters.nRates);
                for (int i = 0; i < parameters.nRates; i++) {seed(i) = unifDist(gen);}
                // seed << 0.1,  0.1,  0.95,  0.17, 0.05,  0.18;
                if(parameters.heldTheta > -1){seed(parameters.heldTheta) = parameters.heldThetaVal;}
                
                /* Evolve initial Global Best and Calculate a Cost*/
                double costSeedK = 0;
                if(ne > 0){
                    for(int i = 0; i < parameters.nRates; i++){
                        if(nestedHolds(i)  != 0){
                            seed(i) = nestedHolds(i) / parameters.hyperCubeScale;
                        }
                    }
                }
                for(int t = 1; t < times.size(); t++){
                    if(parameters.useSBML > 0){
                        MatrixXd XtMat = MatrixXd::Zero(X_0.rows(), X_0.cols());
                        for(int i = 0; i < seed.size(); ++i){
                            theta[i] = seed(i) * parameters.hyperCubeScale;
                        }
                        r.getModel()->setGlobalParameterValues(seed.size(),0,theta); // set new global parameter values here.
                        opt.start = times(0);
                        opt.duration = times(t);
                        for(int i = 0; i < X_0.rows(); ++i){
                            if(specifiedProteins.size() > 0){
                                vector<double> init = r.getFloatingSpeciesInitialConcentrations();
                                for(int p = 0; p < specifiedProteins.size(); p++){
                                    init[specifiedProteins[p]] = X_0(i,p);
                                }
                                r.changeInitialConditions(init);
                            }else{ 
                                r.changeInitialConditions(convertInit(X_0.row(i)));
                            }
                            const DoubleMatrix res = *r.simulate(&opt);
                            for(int j = 0; j < X_0.cols(); ++j){
                                XtMat(i,j) = res[res.numRows() - 1][j + 1];
                            }
                        }
                        VectorXd XtmVec = momentVector(XtMat, nMoments);
                        costSeedK += costFunction(Yt3Vecs[t - 1], XtmVec, weights[t - 1]); 
                    }else{
                        Protein_Components Xt(times(t), nMoments, X_0.rows(), X_0.cols());
                        // Xt = evolveSystem(seed, X_0, nMoments, times(t), dt, times(0));
                        Moments_Mat_Obs XtObs(Xt);
                        Nonlinear_ODE sys(parameters.hyperCubeScale * seed);
                        for (int i = 0; i < X_0.rows(); ++i) {
                            State_N c0 = convertInit(X_0.row(i));
                            Xt.index = i;
                            integrate_adaptive(controlledStepper, sys, c0, times(0), times(t), dt, XtObs);
                        }
                        Xt.mVec /= X_0.rows();  
                        cout << "mVec using boost:" << momentVector(Xt.mat, nMoments).transpose() << endl;
                        costSeedK += costFunction(Yt3Vecs[t - 1], Xt.mVec, weights[t - 1]); // indexed one off for each weight matrix.
                    }
                }
                cout << "PSO Seeded At:"<< seed.transpose() << "| cost:" << costSeedK << endl;
                
                double gCost = costSeedK; //initialize costs and GBMAT
                VectorXd GBVEC = seed;
                
                GBMAT.conservativeResize(GBMAT.rows() + 1, parameters.nRates + 1);
                for (int i = 0; i < parameters.nRates; i++) {
                    GBMAT(GBMAT.rows() - 1, i) = seed(i);
                }
                GBMAT(GBMAT.rows() - 1, parameters.nRates) = gCost;
                double probabilityToTeleport = 3.0/4.0; 
                /* Blind PSO begins */
                cout << "PSO Estimation Has Begun, This may take some time..." << endl;
                for(int step = 0; step < parameters.nSteps; step++){
                #pragma omp parallel for 
                    for(int particle = 0; particle < parameters.nParts; particle++){
                        /* initialize all particle rate constants with unifDist */
                        random_device pRanDev;
                        mt19937 pGen(pRanDev());
                        uniform_real_distribution<double> pUnifDist(low, high);
                        int pSeed = -1;
                        if(parameters.seed > 0){
                            pSeed = particle + step + parameters.seed;
                            pGen.seed(pSeed);
                        }
                        if(step == 0){

                            /* initialize all particles with random rate constant positions */
                            for(int i = 0; i < parameters.nRates; i++){POSMAT(particle, i) = pUnifDist(pGen);}
                            if(parameters.heldTheta > -1){POSMAT.row(particle)(parameters.heldTheta) = parameters.heldThetaVal;}
                            
                            double cost = 0;    
                            if(step == 0 && ne > 0){
                                for(int i = 0 ; i < parameters.nRates; i++){
                                    if(nestedHolds(i) != 0){
                                        POSMAT.row(particle)(i) = nestedHolds(i) / parameters.hyperCubeScale; 
                                    }
                                }
                            }
                            
                            for(int t = 1; t < times.size(); ++t){
                                if(parameters.useSBML > 0){
                                    MatrixXd XtMat = MatrixXd::Zero(X_0.rows(), X_0.cols());
                                    RoadRunner paraModel = r;
                                    double parallelTheta[parameters.nRates];
                                    for(int i = 0; i < parameters.nRates; ++i){
                                        parallelTheta[i] = POSMAT(particle, i) * parameters.hyperCubeScale;
                                    }
                                    paraModel.getModel()->setGlobalParameterValues(parameters.nRates,0, parallelTheta); // set new global parameter values here.
                                    SimulateOptions pOpt = opt;
                                    pOpt.start = times(0);
                                    pOpt.duration = times(t);                   
                                    for(int i = 0; i < X_0.rows(); ++i){
                                        if(specifiedProteins.size() > 0){
                                            vector<double> init = paraModel.getFloatingSpeciesInitialConcentrations();
                                            for(int p = 0; p < specifiedProteins.size(); p++){
                                                init[specifiedProteins[p]] = X_0(i,p);
                                            }
                                            paraModel.changeInitialConditions(init);
                                        }else{
                                            paraModel.changeInitialConditions(convertInit(X_0.row(i)));
                                        }
                                        const DoubleMatrix res = *paraModel.simulate(&pOpt);
                                        for(int j = 0; j < X_0.cols(); ++j){
                                            XtMat(i,j) = res[res.numRows() - 1][j + 1];
                                        }
                                    }
                                    VectorXd XtmVec = momentVector(XtMat, nMoments);
                                    cost += costFunction(Yt3Vecs[t - 1], XtmVec, weights[t - 1]); 
                                }else{
                                    Nonlinear_ODE initSys(POSMAT.row(particle) * parameters.hyperCubeScale);
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
                                
                            }
                            
                            /* instantiate PBMAT */
                            for(int i = 0; i < parameters.nRates; i++){
                                PBMAT(particle, i) = POSMAT(particle, i);
                            }
                            PBMAT(particle, parameters.nRates) = cost; // add cost to final column
                        }else{ 
                            /* step into PSO */
                            double w1 = sfi * pUnifDist(pGen) / sf2, w2 = sfc * pUnifDist(pGen) / sf2, w3 = sfs * pUnifDist(pGen) / sf2;
                            double sumw = w1 + w2 + w3; 
                            w1 = w1 / sumw; w2 = w2 / sumw; w3 = w3 / sumw;
                    
                            VectorXd rpoint = adaptVelocity(POSMAT.row(particle), pSeed, epsi, nan, hone);
                            VectorXd PBVEC(parameters.nRates);
                            for(int i = 0; i < parameters.nRates; ++i){PBVEC(i) = PBMAT(particle, i);}
                            POSMAT.row(particle) = (w1 * rpoint + w2 * PBVEC + w3 * GBVEC); // update position of particle
                        
                            if(parameters.heldTheta > -1){POSMAT.row(particle)(parameters.heldTheta) = parameters.heldThetaVal;}
                            double cost = 0;

                            for(int i = 0 ; i < parameters.nRates; i++){
                                if(nestedHolds(i) != 0){
                                    POSMAT(particle, i) = nestedHolds(i) / parameters.hyperCubeScale; 
                                }
                            }
                            for(int t = 1; t < times.size(); ++t){
                                if(parameters.useSBML > 0){
                                    MatrixXd XtMat = MatrixXd::Zero(X_0.rows(), X_0.cols());
                                    RoadRunner paraModel = r;
                                    double parallelTheta[parameters.nRates];
                                    for(int i = 0; i < parameters.nRates; ++i){
                                        parallelTheta[i] = POSMAT(particle, i) * parameters.hyperCubeScale;
                                    }
                                    paraModel.getModel()->setGlobalParameterValues(parameters.nRates,0, parallelTheta); // set new global parameter values here.
                                    SimulateOptions pOpt = opt;
                                    pOpt.start = times(0);
                                    pOpt.duration = times(t);                   
                                    for(int i = 0; i < X_0.rows(); ++i){
                                        if(specifiedProteins.size() > 0){
                                            vector<double> init = paraModel.getFloatingSpeciesInitialConcentrations();
                                            for(int p = 0; p < specifiedProteins.size(); p++){
                                                init[specifiedProteins[p]] = X_0(i,p);
                                            }
                                            paraModel.changeInitialConditions(init);
                                        }else{
                                            paraModel.changeInitialConditions(convertInit(X_0.row(i)));
                                        }
                                        const DoubleMatrix res = *paraModel.simulate(&pOpt);
                                        for(int j = 0; j < X_0.cols(); ++j){
                                            XtMat(i,j) = res[res.numRows() - 1][j + 1];
                                        }
                                    }
                                    VectorXd XtmVec = momentVector(XtMat, nMoments);
                                    cost += costFunction(Yt3Vecs[t - 1], XtmVec, weights[t - 1]); 
                                }else{
                                    /*solve ODEs and recompute cost */
                                    Protein_Components XtPSO(times(t), nMoments, X_0.rows(), X_0.cols());
                                    Moments_Mat_Obs XtObsPSO1(XtPSO);
                                    Nonlinear_ODE stepSys(POSMAT.row(particle) * parameters.hyperCubeScale);
                                    for(int i = 0; i < X_0.rows(); ++i){
                                        State_N c0 = convertInit(X_0.row(i));
                                        XtPSO.index = i;
                                        integrate_adaptive(controlledStepper, stepSys, c0, times(0), times(t), dt, XtObsPSO1);
                                    }
                                    XtPSO.mVec /= X_0.rows();
                                    cost += costFunction(Yt3Vecs[t - 1], XtPSO.mVec, weights[t - 1]);
                                }
                            }
                        
                            /* update gBest and pBest */
                        #pragma omp critical
                        {
                            if(cost < PBMAT(particle, parameters.nRates)){ // particle best cost
                                for(int i = 0; i < parameters.nRates; i++){
                                    PBMAT(particle, i) = POSMAT.row(particle)(i);
                                }
                                PBMAT(particle, parameters.nRates) = cost;
                                if(cost < gCost){
                                    gCost = cost;
                                    GBVEC = POSMAT.row(particle);
                                }   
                            }
                        }
                        }
                    }
                    GBMAT.conservativeResize(GBMAT.rows() + 1, parameters.nRates + 1); // Add to GBMAT after resizing
                    for (int i = 0; i < parameters.nRates; i++) {GBMAT(GBMAT.rows() - 1, i) = GBVEC(i);}
                    GBMAT(GBMAT.rows() - 1, parameters.nRates) = gCost;
                    sfi = sfi - (sfe - sfg) / parameters.nSteps;   // reduce the inertial weight after each step 
                    sfs = sfs + (sfe - sfg) / parameters.nSteps;
                }

                

                // double hypercube size and rescan global best vector to see what needs to be held constant.
                if(parameters.nest > 1){
                    for(int i = 0; i < GBVEC.size(); i++){
                        if(parameters.hyperCubeScale * GBVEC(i) < parameters.hyperCubeScale - epsi){
                            nestedHolds(i) = parameters.hyperCubeScale * GBVEC(i);
                        }
                    }
                }
                for(int i = 0; i < parameters.nRates; i++){
                    if(nestedHolds(i) != 0){
                        GBVECS(run, i) = nestedHolds(i); // GBVECS is either something that was held constant.
                    }else{ 
                        GBVECS(run, i) = parameters.hyperCubeScale * GBVEC(i); // or is now something scaled to GBVEC.
                    } 
                }
                GBVECS(run, parameters.nRates) = gCost;
                if(parameters.nest > 1){
                    parameters.hyperCubeScale *= 2.0;
                }
            } // nested loop

            /* bootstrap X0 and Y matrices if more than 1 run is specified */
            if(parameters.nRuns > 1 && parameters.bootstrap > 0){
                X_0 = bootStrap(ogX_0);
                for(int y = 0; y < Yt3Mats.size(); ++y){
                    Yt3Mats[y] = bootStrap(ogYt3Mats[y]);
                    Yt3Vecs[y] = momentVector(Yt3Mats[y], nMoments);
                }
                cout << "bootstrap means" << endl << "X_0:" << X_0.colwise().mean() << endl << "Yt:" << Yt3Mats[0].colwise().mean() << endl;
            }
            if(parameters.nest > 1){
                for(int ne = 0; ne < parameters.nest; ne++){ // reset cube for each run
                    parameters.hyperCubeScale /= 2.0;
                }
            }
            cout << "Held Parameter Estimates:" << nestedHolds.transpose() << endl; 
           
        } // run loop
        // when done, find what the original max hypercube size was from nesting
        for(int ne = 1; ne < parameters.nest; ne++){ 
            parameters.hyperCubeScale *= 2.0;
        }
        cout << "Hypercubescale Max:" << parameters.hyperCubeScale << endl;
        if(parameters.simulateYt > 0){cout << "Simulated Truth:" << tru.transpose() << endl;}
        if(parameters.reportMoments > 0){
            cout << "Average Estimate:" << GBVECS.colwise().mean() << endl;
            cout << "Simulated Xt Moments For Various Times:" << endl;  
            for(int t = 1; t < times.size(); ++t){
                if(parameters.useSBML > 0){
                    VectorXd avgRunPos = VectorXd::Zero(parameters.nRates);
                    for(int i = 0; i < avgRunPos.size(); ++i){
                        avgRunPos(i) = GBVECS.colwise().mean()(i);
                    }
                    MatrixXd XtMat = MatrixXd::Zero(X_0.rows(), X_0.cols());
                    for(int i = 0; i < avgRunPos.size(); ++i){
                        theta[i] = avgRunPos(i);
                    }
                    r.getModel()->setGlobalParameterValues(avgRunPos.size(),0,theta); // set new global parameter values here.
                    opt.start = times(0);
                    opt.duration = times(t);
                    for(int i = 0; i < X_0.rows(); ++i){
                        if(specifiedProteins.size() > 0){
                            vector<double> init = r.getFloatingSpeciesInitialConcentrations();
                            for(int p = 0; p < specifiedProteins.size(); p++){
                                init[specifiedProteins[p]] = X_0(i,p);
                            }
                            r.changeInitialConditions(init);
                        }else{    
                            r.changeInitialConditions(convertInit(X_0.row(i)));
                        }
                        const DoubleMatrix res = *r.simulate(&opt);
                        for(int j = 0; j < X_0.cols(); ++j){
                            XtMat(i,j) = res[res.numRows() - 1][j + 1];
                        }
                    }
                    VectorXd XtmVec = momentVector(XtMat, nMoments);
                    cout << times(t) << " " << XtmVec.transpose() << endl;
                }else{
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
                    cout << times(t) << " " << XtPSO.mVec.transpose() << endl;
                }
            }
        }
    }
    cout << endl << "-------------- All Run Estimates: -------------------" << endl;
    for (int i = 0; i < GBVECS.cols() - 1; i++){cout << parameterNames[i] << "\t\t";}
    cout << "cost" << endl;
    cout << GBVECS << endl;
    /* Compute 95% CI's with basic z=1.96 normal distribution assumption for now if n>1 */
    if(parameters.nRuns > 1){computeConfidenceIntervals(GBVECS, 1.96, parameters.nRates);}
    
    auto tB = std::chrono::high_resolution_clock::now();
    auto bDuration = std::chrono::duration_cast<std::chrono::seconds>(tB - t1).count();
    cout << "CODE FINISHED RUNNING IN " << bDuration << " s TIME!" << endl;
    return EXIT_SUCCESS;
}