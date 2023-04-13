/* 
Author: John W., Dr. Stewart
Summary of Src File: 
CyGMM -> <Insert Statement About Rate Constant Estimation>

Coder's Note:
If you're wondering why the rungekutta solution method is all in main and not segmented in nonlinear.cpp is because currently
the way in how Boost's ODE struct/class system for solving ODEs results in major performance decreases when inside a function, very possibly
due to something with memory or how functions work from the compiler. Also libRR is substantially slower than boost ode solvers and Eigen.expm.

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
#include "graph.hpp"
int main(int argc, char** argv){
    auto t1 = std::chrono::high_resolution_clock::now();
    /* Input Parameters for Program */
    if(helpCall(argc, argv)){return EXIT_SUCCESS;}
    cout << "Program Begin:" << endl;
    cout << "** Please Make Sure That All Inputted Files are in the UNIX Line Formatting to Prevent Bugs! To see the full list of commands with BNGMM, please do ./BNGMM -h **" << endl;
    Parameters parameters = Parameters(getConfigPath(argc, argv));
    parameters.useSBML = int(modelPathExists(argc,argv));
    if(outPathExists(argc, argv)){
        parameters.outPath = getOutputPath(argc,argv);       
        fs::create_directory(parameters.outPath);
    }
    /* Read time steps*/
    VectorXd times = readCsvTimeParam(getTimeStepsPath(argc, argv));
    /* Important! Max Thread Count Test */
    omp_set_num_threads(parameters.nThreads);

    /* Read First Matrix */
    MatrixXd x0;
    MatrixXd ogx0; // copy of initial X values if bootstrap is used.
    vector<MatrixXd> xt3Mats;
    x0 = readX(getXPath(argc, argv));
    ogx0 = x0;
    /* Decide Number of Moments to Use in Estimation */
    int nMoments = (x0.cols() * (x0.cols() + 3)) / 2;
    if(parameters.useOnlySecMom > 0){  // these will be added to the options sheet later.
        nMoments = 2 * x0.cols();
    }
    if(parameters.useOnlyFirstMom > 0){
        nMoments = x0.cols();
    }
    parameters.printParameters(nMoments, times);
    parameters.nMoments = nMoments;
    /* Parameter Estimate Matrices */
    MatrixXd GBMAT; // per PSO
    MatrixXd GBVECS = MatrixXd::Zero(parameters.nRuns, parameters.nRates + 1); // for each run
    
    /*---------------------- Nonlinear Setup PSO ------------------------ */
    double dt = 1.0; // nonlinear time evolution variables
    /* Explicit Boundary Parameters */
    double squeeze = 0.500, sdbeta = 0.10; // how much to shrink PSO search over time (how much variability each position is iterated upon)
    double sf2 = 1; // factor that can be used to regularize particle weights (global, social, inertial)
    double epsi = 0.02;
    double nan = 0.005;
    /* PSO params */
    double sfp = parameters.pBestWeight, sfg = parameters.globalBestWeight, sfe = parameters.pInertia; // initial particle historical weight, global weight social, in ertial
    int hone = 28; 
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
    vector<MatrixXd> yt3Mats;
    vector<MatrixXd> ogYt3Mats;
    vector<VectorXd> yt3Vecs; // vector of observed moments for each time point (each element in in this vector is a vector)
    vector<MatrixXd> allMomentsAcrossTime; // we can store all moment vectors that are generated through each run in matrices for each time point.
    for(int t = 1; t < times.size();  ++t){
        allMomentsAcrossTime.push_back(MatrixXd::Zero(parameters.nRuns, parameters.nMoments));
    }
    if(parameters.simulateYt < 1){        
        yt3Mats = readY(getYPath(argc, argv));
        if(yt3Mats.size() + 1 != times.size()){
            cout << "Error, number of Y_t files read in do not match the number of timesteps!" << endl;
            exit(1);
        }
        // filter all zeroes and compute moments vectors for cost calcs
        for(int i = 0; i < yt3Mats.size(); i++){
            yt3Mats[i] = filterZeros(yt3Mats[i]);
            cout << "After removing all negative rows, Y"<< i << " has " << yt3Mats[i].rows() << " rows." << endl;
            yt3Vecs.push_back(momentVector(yt3Mats[i], nMoments));
            cout << "t" << times(i+1) << " moments:"<< yt3Vecs[i].transpose() << endl;
        }
        ogYt3Mats = yt3Mats;
    }
    /* Default Bionetgen Mode */
    vector<string> parameterNames;
    if(parameters.useSBML > 0){
        const string bngl = ".bngl"; // suffixes for sbml/bngl file types
        const string sbml = "_sbml.xml";

        /* Run an update for bngl */
        string modelPath = getModelPath(argc, argv);
        string file_without_extension = getFileNameWithoutExtensions(modelPath);
        string sbmlModel = "sbml/"+ file_without_extension + sbml;
        const string bnglCall = "bionetgen run -i" + modelPath + " -o sbml";
        Grapher graph = Grapher(parameters.outPath, file_without_extension, getTrueRatesPath(argc, argv), times, parameters.nRates);
        if(system(bnglCall.c_str()) < 0){
            cout << "Error Running BioNetGen -> Make sure you have installed it through pip install bionetgen or check program permissions!" << endl;
            return EXIT_FAILURE;
        }
        /* RoadRunner Configuration and Simulation Variables */
        RoadRunner r = RoadRunner(sbmlModel);
        vector<int> specifiedProteins;
        parameterNames = r.getGlobalParameterIds(); // parameter names check
        parameterNames = vector<string>(parameterNames.begin(),parameterNames.begin() + parameters.nRates);
        parameterNames.push_back("cost");
        vector<string> speciesNames =  getSpeciesNames(sbmlModel);
        cout << "--------------------------------------------------------" << endl;
        if(r.getNumberOfIndependentSpecies() +  r.getNumberOfDependentSpecies() != x0.cols()){
            if(x0.cols() > r.getNumberOfIndependentSpecies() +  r.getNumberOfDependentSpecies() ){
                cout << "Error Too Many Species/Columns in X.csv file! Please remove some before continuing!" << endl;
                cout << "Expected:" <<  r.getNumberOfIndependentSpecies() +  r.getNumberOfDependentSpecies() << " Got:" << x0.cols() << endl;
                return EXIT_FAILURE;
            }
            cout << "Number of Species Defined in BNGL does not match number of columns in data files! Now listing all species in system and respective indices in order!" << endl;
            cout << "Note: User can supply a \"./BNGMM -p protein_observed.txt \" to specify explicitly which proteins are observed in data. Please make sure names are in order from top to bottom matching left to right in data csv file." << endl;
            for(int i = 0; i < speciesNames.size(); i++){
                cout << "(" << i << ") " << speciesNames[i] << endl; 
            }
            if(proPathExists(argc, argv)){
                specifiedProteins = specifySpeciesFromProteinsList(getProPath(argc,argv), speciesNames, x0.cols());
                if(specifiedProteins.size() < 1){
                    cout << "Error, no proteins specified!" << endl;
                    exit(EXIT_FAILURE);
                }
                cout << "From Above List of Indexed Species, We are using..." << endl;
                for(int i = 0; i < specifiedProteins.size(); i++){
                    cout << "(" <<specifiedProteins[i] <<") "<< specifiedProteins[i] <<endl;
                }
            }else{
                cout << "No Proteins Specified Using Argument \"-p Proteins.txt\", Thus Using First " << x0.cols() << " Species Listed" << endl;
            }
            cout << "Proteins Not Observed Will Default to Initial Values Defined in .BNGL File" << endl;
            // check if some txt -p file has been used and use that otherwise have user manually select.
        }else{
            cout << "------- Matching Columns of X Data files to Ids -------" << endl;
            for(int i = 0; i < speciesNames.size(); i++){
                cout << speciesNames[i] << " to column:"<< i << " with first value:" << x0(0,i) << endl;
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
                yt3Vecs.push_back(momentVector(YtMat, nMoments));
                yt3Mats.push_back(YtMat);
                cout << times(t) << " "<< yt3Vecs[t-1].transpose() << endl;
            }
            cout << "--------------------------------------------------------" << endl;
        }

        MatrixXd heldTheta;
        /* HeldRates if it happens*/ 
        if (holdRates(argc, argv)){
            cout << "Held Rate Constants With Respect to Their Indices:" << endl;
            const string hr = getHeldRatesDir(argc, argv);
            cout << hr << endl;
            heldTheta = heldThetas(parameters.nRates, hr);
            cout << heldTheta << endl;
        }

        /* Compute initial wolfe weights */
        if (yt3Mats[0].cols() != x0.cols()){
            cout << "Error, mismatch in number of species/columns between X and Y!" << endl;
            cout << "X:" << x0.cols() << " Y:" << yt3Mats[0].cols() << endl;
            return EXIT_FAILURE;
        }
        for(int y = 0; y < yt3Mats.size(); ++y){ 
            weights.push_back(wolfWtMat(yt3Mats[y], nMoments, parameters.useInverse > 0));
            cout << "--------------------------------------------------------" << endl;
            cout << "Computed GMM Weight Matrix:" << endl;
            cout << weights[y] << endl;
            // matrixToCsv(weights[y], parameters.outPath + file_without_extension + "_weight_t" + to_string_with_precision(times(y+1),2));
            cout << "--------------------------------------------------------" << endl << endl;
        }

        /* Contour Function - ONLY RUNS IF SIMULATED OR IF SEEDED */
        if(contour(argc, argv) && (seedRates(argc, argv) || parameters.simulateYt > 0 )){
            cout << "--------------------------------------------------------" << endl;
            cout << "Generating Contour Files With" << endl;
            int stepSize = 25;
            vector<int> pairwise1;
            vector<int> pairwise2;
            for(int i = 0; i < parameters.nRates - 1; ++i){
                pairwise1.push_back(i);
                pairwise2.push_back(i+1);
            }

            VectorXd contourTheta;
            if(parameters.simulateYt > 0){
                contourTheta = tru;
            }else if(seedRates(argc, argv)){
                contourTheta = readSeed(parameters.nRates, getSeededRates(argc,argv));
            }
            cout << "Theta:" << contourTheta.transpose() << endl;
            #pragma omp parallel for schedule(dynamic)
                for(int thta = 0; thta < parameters.nRates-1; ++thta){
                    vector<string> contourLabels;
                    MatrixXd contour = MatrixXd::Zero(stepSize*stepSize, 3); // # of pairwise params x Coordinate + Cost
                    int fIdx = pairwise1[thta];
                    int sIdx = pairwise2[thta];
                    contourLabels.push_back(parameterNames[fIdx]);
                    contourLabels.push_back(parameterNames[sIdx]);
                    contourLabels.push_back("Cost");
                    int cdex = 0;
                    RoadRunner paraMod = r;
                    double firstTheta = 0;
                    for(int idxs = 0; idxs < stepSize; ++idxs){
                        double secondTheta = 0;
                        for(int jdx = 0; jdx < stepSize; ++jdx){
                            double gmm = 0;
                            for(int t = 1; t < times.size(); t++){
                                MatrixXd XtMat = MatrixXd::Zero(x0.rows(), x0.cols());
                                for(int idx = 0; idx < contourTheta.size(); ++idx){
                                    theta[idx] = contourTheta(idx);
                                }
                                theta[fIdx] = firstTheta;
                                theta[sIdx] = secondTheta;
                                paraMod.getModel()->setGlobalParameterValues(contourTheta.size(),0,theta); // set new global parameter values here.
                                opt.start = times(0);
                                opt.duration = times(t);
                                for(int i = 0; i < x0.rows(); ++i){
                                    if(specifiedProteins.size() > 0){
                                        vector<double> init = paraMod.getFloatingSpeciesInitialConcentrations();
                                        for(int p = 0; p < specifiedProteins.size(); p++){
                                            init[specifiedProteins[p]] = x0(i,p);
                                        }
                                        r.changeInitialConditions(init);
                                    }else{ 
                                        r.changeInitialConditions(convertInit(x0.row(i)));
                                    }
                                    const DoubleMatrix res = *paraMod.simulate(&opt);
                                    for(int j = 0; j < x0.cols(); ++j){
                                        XtMat(i,j) = res[res.numRows() - 1][j + 1];
                                    }
                                }
                                VectorXd XtmVec = momentVector(XtMat, nMoments);
                                gmm += costFunction(yt3Vecs[t - 1], XtmVec, weights[t - 1]); 
                            }
                            contour(cdex,0) = firstTheta;
                            contour(cdex,1) = secondTheta;
                            contour(cdex,2) = gmm;
                            secondTheta += 1.0 / stepSize;
                            // contour()
                            cdex++; 
                        }
                        firstTheta += 1.0 / stepSize;
                    }
                    #pragma omp critical
                    {
                        matrixToCsvWithLabels(contour, contourLabels, parameters.outPath + file_without_extension + "_contour" + to_string(fIdx) + "_" + to_string(sIdx));
                    }
                }
            graph.graphContours(parameters.nRates);
            cout << "--------------------------------------------------------" << endl;
        }
        
        /*------------ PSO SECTION ------------*/
        for(int run = 0; run < parameters.nRuns; ++run){ // for multiple runs aka bootstrapping (for now)
            if (run > 0 && parameters.bootstrap > 0){
                for(int y = 0; y < yt3Mats.size(); ++y){ 
                    weights[y] = wolfWtMat(yt3Mats[y], nMoments, parameters.useInverse > 0);
                }
            }
            // make sure to reset GBMAT, POSMAT, AND PBMAT every run
            double sfi = sfe, sfc = sfp, sfs = sfg; // below are the variables being used to reiterate weights
            GBMAT = MatrixXd::Zero(0,0); // iterations of global best vectors
            MatrixXd PBMAT = MatrixXd::Zero(parameters.nParts, parameters.nRates + 1); // particle best matrix + 1 for cost component
            MatrixXd POSMAT = MatrixXd::Zero(parameters.nParts, parameters.nRates); // Position matrix as it goees through it in parallel
            
            /* Initialize Global Best  */
            VectorXd seed = VectorXd::Zero(parameters.nRates);
            for (int i = 0; i < parameters.nRates; i++) {seed(i) = unifDist(gen);}
            // seed << 0.1,  0.1,  0.95,  0.17, 0.05,  0.18;
            if(holdRates(argc,argv)){
                for(int i = 0; i < heldTheta.rows(); ++i){
                    if (heldTheta(i,0) != 0){
                        seed(i) = heldTheta(i,1);
                    }
                }
            }
            if(seedRates(argc, argv)){
                seed = readSeed(parameters.nRates, getSeededRates(argc,argv));
            }
            
            /* Evolve initial Global Best and Calculate a Cost*/
            double costSeedK = 0;
            for(int t = 1; t < times.size(); t++){
                MatrixXd XtMat = MatrixXd::Zero(x0.rows(), x0.cols());
                for(int i = 0; i < seed.size(); ++i){
                    theta[i] = parameters.hyperCubeScale * seed(i);
                }
                r.getModel()->setGlobalParameterValues(seed.size(),0,theta); // set new global parameter values here.
                opt.start = times(0);
                opt.duration = times(t);
                for(int i = 0; i < x0.rows(); ++i){
                    if(specifiedProteins.size() > 0){
                        vector<double> init = r.getFloatingSpeciesInitialConcentrations();
                        for(int p = 0; p < specifiedProteins.size(); p++){
                            init[specifiedProteins[p]] = x0(i,p);
                        }
                        r.changeInitialConditions(init);
                    }else{ 
                        r.changeInitialConditions(convertInit(x0.row(i)));
                    }
                    const DoubleMatrix res = *r.simulate(&opt);
                    for(int j = 0; j < x0.cols(); ++j){
                        XtMat(i,j) = res[res.numRows() - 1][j + 1];
                    }
                }
                VectorXd XtmVec = momentVector(XtMat, nMoments);
                costSeedK += costFunction(yt3Vecs[t - 1], XtmVec, weights[t - 1]); 
            }
            cout << "PSO Seeded At:"<< seed.transpose() << "| cost:" << costSeedK << endl;
            
            double gCost = costSeedK; //initialize costs and GBMAT
            VectorXd GBVEC = seed;
            VectorXd scaledGBVEC = parameters.hyperCubeScale * GBVEC;
            
            GBMAT.conservativeResize(GBMAT.rows() + 1, parameters.nRates + 1);
            for (int i = 0; i < parameters.nRates; i++) {
                GBMAT(GBMAT.rows() - 1, i) = seed(i);
            }
            GBMAT(GBMAT.rows() - 1, parameters.nRates) = gCost;
            double probabilityToTeleport = 3.0/4.0; 
            /* Blind PSO begins */
            cout << "PSO Estimation Has Begun, This may take some time..." << endl;
            for(int step = 0; step < parameters.nSteps; step++){
            #pragma omp parallel for schedule(dynamic)
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
                    VectorXd scaledPos = VectorXd::Zero(parameters.nRates);
                    if(step == 0){
                        /* initialize all particles with random rate constant positions */
                        for(int i = 0; i < parameters.nRates; i++){
                            POSMAT(particle, i) = pUnifDist(pGen);
                            scaledPos(i)= parameters.hyperCubeScale * POSMAT(particle,i);
                        }
                        if(holdRates(argc,argv)){
                            for(int i = 0; i < heldTheta.rows(); ++i){
                                if (heldTheta(i,0) != 0){
                                    scaledPos(i) = heldTheta(i,1);
                                }
                            }
                        }
                        
                        double cost = 0;    
                        
                        for(int t = 1; t < times.size(); ++t){
                            MatrixXd XtMat = MatrixXd::Zero(x0.rows(), x0.cols());
                            RoadRunner paraModel = r;
                            double parallelTheta[parameters.nRates];

                            for(int i = 0; i < parameters.nRates; ++i){
                                parallelTheta[i] = scaledPos(i);
                            }

                            paraModel.getModel()->setGlobalParameterValues(parameters.nRates,0, parallelTheta); // set new global parameter values here.
                            SimulateOptions pOpt = opt;
                            pOpt.start = times(0);
                            pOpt.duration = times(t);                   
                            for(int i = 0; i < x0.rows(); ++i){
                                if(specifiedProteins.size() > 0){
                                    vector<double> init = paraModel.getFloatingSpeciesInitialConcentrations();
                                    for(int p = 0; p < specifiedProteins.size(); p++){
                                        init[specifiedProteins[p]] = x0(i,p);
                                    }
                                    paraModel.changeInitialConditions(init);
                                }else{
                                    paraModel.changeInitialConditions(convertInit(x0.row(i)));
                                }
                                const DoubleMatrix res = *paraModel.simulate(&pOpt);
                                for(int j = 0; j < x0.cols(); ++j){
                                    XtMat(i,j) = res[res.numRows() - 1][j + 1];
                                }
                            }
                            VectorXd XtmVec = momentVector(XtMat, nMoments);
                            cost += costFunction(yt3Vecs[t - 1], XtmVec, weights[t - 1]); 
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
                        scaledPos = parameters.hyperCubeScale*POSMAT.row(particle);
                        if(holdRates(argc,argv)){
                            for(int i = 0; i < heldTheta.rows(); ++i){
                                if (heldTheta(i,0) != 0){
                                    scaledPos(i) = heldTheta(i,1);
                                    POSMAT.row(particle)(i) = scaledPos(i);
                                }
                            }
                        }
                        // if(holdRates(argc, argv)){POSMAT.row(particle)(parameters.heldTheta) = parameters.heldThetaVal;}
                        double cost = 0;
                    
                        for(int t = 1; t < times.size(); ++t){ 
                            MatrixXd XtMat = MatrixXd::Zero(x0.rows(), x0.cols());
                            RoadRunner paraModel = r;
                            double parallelTheta[parameters.nRates];
                            for(int i = 0; i < parameters.nRates; ++i){
                                parallelTheta[i] = scaledPos(i);
                            }
                            paraModel.getModel()->setGlobalParameterValues(parameters.nRates,0, parallelTheta); // set new global parameter values here.
                            SimulateOptions pOpt = opt;
                            pOpt.start = times(0);
                            pOpt.duration = times(t);                   
                            for(int i = 0; i < x0.rows(); ++i){
                                if(specifiedProteins.size() > 0){
                                    vector<double> init = paraModel.getFloatingSpeciesInitialConcentrations();
                                    for(int p = 0; p < specifiedProteins.size(); p++){
                                        init[specifiedProteins[p]] = x0(i,p);
                                    }
                                    paraModel.changeInitialConditions(init);
                                }else{
                                    paraModel.changeInitialConditions(convertInit(x0.row(i)));
                                }
                                const DoubleMatrix res = *paraModel.simulate(&pOpt);
                                for(int j = 0; j < x0.cols(); ++j){
                                    XtMat(i,j) = res[res.numRows() - 1][j + 1];
                                }
                            }
                            VectorXd XtmVec = momentVector(XtMat, nMoments);
                            cost += costFunction(yt3Vecs[t - 1], XtmVec, weights[t - 1]);   
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
                                scaledGBVEC = parameters.hyperCubeScale * GBVEC;
                            }   
                        }
                    }
                    }
                }
                GBMAT.conservativeResize(GBMAT.rows() + 1, parameters.nRates + 1); // Add to GBMAT after resizing
                for (int i = 0; i < parameters.nRates; i++) {GBMAT(GBMAT.rows() - 1, i) = scaledGBVEC(i);} // ideally want to save
                GBMAT(GBMAT.rows() - 1, parameters.nRates) = gCost;
                sfi = sfi - (sfe - sfg) / parameters.nSteps;   // reduce the inertial weight after each step 
                sfs = sfs + (sfe - sfg) / parameters.nSteps;
            }
            cout << "----------------PSO Best Each Iterations----------------" << endl;
            cout << GBMAT << endl;
            cout << "--------------------------------------------------------" << endl;


            for(int i = 0; i < parameters.nRates; i++){
                GBVECS(run, i) = scaledGBVEC(i); // or is now something scaled to GBVEC. 
            }
            GBVECS(run, parameters.nRates) = gCost;
            cout << GBVECS.row(run) << endl;
            /* bootstrap X0 and Y matrices if more than 1 run is specified */
            if(parameters.nRuns > 1 && parameters.bootstrap > 0){
                x0 = bootStrap(ogx0);
                for(int y = 0; y < yt3Mats.size(); ++y){
                    yt3Mats[y] = bootStrap(ogYt3Mats[y]);
                    yt3Vecs[y] = momentVector(yt3Mats[y], nMoments);
                }
                cout << "bootstrap means" << endl << "x0:" << x0.colwise().mean() << endl << "Yt:" << yt3Mats[0].colwise().mean() << endl;
            }

        } // run loop

        VectorXd leastCostRunPos = VectorXd::Zero(parameters.nRates);
        int indexOfLeastCost = 0;
        double currentCost = GBVECS(indexOfLeastCost, GBVECS.cols() - 1);
        for(int i = 0; i < GBVECS.rows(); ++i){
            if(currentCost > GBVECS(i, GBVECS.cols() - 1)){
                indexOfLeastCost = i;
                currentCost = GBVECS(i,GBVECS.cols() - 1);
            }
        }
        for(int i = 0; i < leastCostRunPos.size(); ++i){
            leastCostRunPos(i) = GBVECS(indexOfLeastCost, i);
        }
        
        vectorToCsv(leastCostRunPos, parameters.outPath + file_without_extension + "_leastCostEstimate");
        for(int t = 1; t < times.size(); ++t){
            MatrixXd XtMat = MatrixXd::Zero(x0.rows(), x0.cols());
            for(int i = 0; i < leastCostRunPos.size(); ++i){
                theta[i] = leastCostRunPos(i);
            }
            r.getModel()->setGlobalParameterValues(leastCostRunPos.size(),0,theta); // set new global parameter values here.
            opt.start = times(0);
            opt.duration = times(t);
            for(int i = 0; i < x0.rows(); ++i){
                if(specifiedProteins.size() > 0){
                    vector<double> init = r.getFloatingSpeciesInitialConcentrations();
                    for(int p = 0; p < specifiedProteins.size(); p++){
                        init[specifiedProteins[p]] = x0(i,p);
                    }
                    r.changeInitialConditions(init);
                }else{    
                    r.changeInitialConditions(convertInit(x0.row(i)));
                }
                const DoubleMatrix res = *r.simulate(&opt);
                for(int j = 0; j < x0.cols(); ++j){
                    XtMat(i,j) = res[res.numRows() - 1][j + 1];
                }
            }
            VectorXd XtmVec = momentVector(XtMat, nMoments);
            xt3Mats.push_back(XtMat);    
            reportLeastCostMoments(XtmVec,yt3Vecs[t-1],times(t), parameters.outPath + file_without_extension); // FIND BEST FIT.
            if(parameters.reportMoments > 0){
                cout << "--------------------------------------------------------" << endl;
                cout << "For Least Cost Estimate:" << leastCostRunPos.transpose() << endl;
                cout << "RSS (NOT GMM) COST FROM DATASET:" << costFunction(XtmVec, yt3Vecs[t-1], MatrixXd::Identity(nMoments, nMoments)) << endl;
                cout << "t                  moments" << endl;
                cout << times(t) << " " << XtmVec.transpose() << endl;
                cout << "--------------------------------------------------------" << endl;
            }
        }

        for(int n = 0; n < GBVECS.rows(); ++n ){
            for(int t = 1; t < times.size(); ++t){
                MatrixXd XtMat = MatrixXd::Zero(x0.rows(), x0.cols());
                for(int j = 0; j < GBVECS.cols() - 1; ++j){
                    theta[j] = GBVECS(n,j);
                }
                r.getModel()->setGlobalParameterValues(leastCostRunPos.size(),0,theta); // set new global parameter values here.
                opt.start = times(0);
                opt.duration = times(t);
                for(int i = 0; i < x0.rows(); ++i){
                    if(specifiedProteins.size() > 0){
                        vector<double> init = r.getFloatingSpeciesInitialConcentrations();
                        for(int p = 0; p < specifiedProteins.size(); p++){
                            init[specifiedProteins[p]] = x0(i,p);
                        }
                        r.changeInitialConditions(init);
                    }else{    
                        r.changeInitialConditions(convertInit(x0.row(i)));
                    }
                    const DoubleMatrix res = *r.simulate(&opt);
                    for(int j = 0; j < x0.cols(); ++j){
                        XtMat(i,j) = res[res.numRows() - 1][j + 1];
                    }
                }
                VectorXd XtmVec = momentVector(XtMat, nMoments);
                allMomentsAcrossTime[t-1].row(n) = XtmVec; 
            }
        }

        /* Save Data for Plotting */
        reportAllMoments(allMomentsAcrossTime,yt3Vecs,times, parameters.outPath + file_without_extension);
        for(int y = 0; y < yt3Mats.size(); ++y){ 
            matrixToCsvWithLabels(yt3Mats[y], speciesNames, parameters.outPath + file_without_extension + "Yt" + to_string_with_precision(times(y + 1), 2));
        }
        for(int t = 1; t < times.size(); t++){
            matrixToCsvWithLabels(xt3Mats[t-1], speciesNames, parameters.outPath + file_without_extension + "Xt" + to_string_with_precision(times(t),2));
        }
        matrixToCsvWithLabels(GBVECS, parameterNames, parameters.outPath + file_without_extension + "_estimates");

        // Graphing Time
        /* Necessary Graphing Initialization */
        cout << "Plotting R^2 Plot and Confidence Intervals!" << endl;
        graph.graphMoments(xt3Mats[0].cols());
        graph.graphConfidenceIntervals(parameters.simulateYt > 0 );

        if(forecast(argc, argv)){
            VectorXd futureT = readCsvTimeParam(getForecastedTimes(argc, argv));
            VectorXd avgMu = GBVECS.colwise().mean();
            MatrixXd futurecast = MatrixXd::Zero(futureT.size(), nMoments + 1);
            MatrixXd observedData = MatrixXd::Zero(times.size(), nMoments + 1);
            for(int t = 1; t < times.size(); ++t){
                observedData(t - 1, 0) = times(t);
                for(int mom = 1; mom < nMoments + 1; ++mom){
                    observedData(t - 1, mom) = yt3Vecs[t - 1](mom - 1);
                }
            }
            matrixToCsv(observedData, parameters.outPath + file_without_extension + "_observed");
            /* Calculate New Moments */
            cout << "--------------- Forecasted Moments in Time: ----------" << endl;
            for(int t = 0; t < futureT.size(); ++t){
                MatrixXd XtMat = MatrixXd::Zero(x0.rows(), x0.cols());
                for(int j = 0; j < GBVECS.cols() - 1; ++j){
                    theta[j] = avgMu(j);
                }
                r.getModel()->setGlobalParameterValues(avgMu.size() - 1, 0, theta); // set new global parameter values here.
                opt.start = futureT(0);
                opt.duration = futureT(t);
                for(int i = 0; i < x0.rows(); ++i){
                    if(specifiedProteins.size() > 0){
                        vector<double> init = r.getFloatingSpeciesInitialConcentrations();
                        for(int p = 0; p < specifiedProteins.size(); p++){
                            init[specifiedProteins[p]] = x0(i,p);
                        }
                        r.changeInitialConditions(init);
                    }else{    
                        r.changeInitialConditions(convertInit(x0.row(i)));
                    }
                    const DoubleMatrix res = *r.simulate(&opt);
                    for(int j = 0; j < x0.cols(); ++j){
                        XtMat(i,j) = res[res.numRows() - 1][j + 1];
                    }
                }
                VectorXd XtmVec = momentVector(XtMat, nMoments);
                cout << futureT(t) <<" " << XtmVec.transpose() << endl;
                futurecast(t,0) = futureT(t);
                for(int mom = 1; mom < nMoments + 1; ++mom){
                    futurecast(t,mom) = XtmVec(mom - 1);
                }
            }
            cout << "------------------------------------------------------" << endl;
            matrixToCsv(futurecast,  parameters.outPath + file_without_extension + "_forecast");
            graph.graphForecasts(x0.cols());
        }

    /* 
    ******************************************************************************************************************************
    ********************* Should user decide to compile their own C++ ode system.************************************************* 
    ******************************************************************************************************************************
    ******************************************************************************************************************************
    */
    }else{
        cout << "Error No .BNGL Model Specified! Please specify a model by \"./BNGMM -m model.bngl\". To get more possible BNGMM parameters, please do \"./BNGMM -h\". Exiting!" << endl;
        return EXIT_FAILURE;
    }
    /* 
    ******************************************************************************************************************************
    ************************************************************** PSO END ************************************************* 
    ******************************************************************************************************************************
    ******************************************************************************************************************************
    */
    cout << endl << "--------------- All Run Estimates: -------------------" << endl;
    if(parameters.useSBML > 0){
        for (int i = 0; i < parameterNames.size(); i++){cout << parameterNames[i] << " ";}
        cout << endl;
    }
    cout << GBVECS << endl;
    /* Compute 95% CI's with basic z=1.96 normal distribution assumption for now if n>1 */
    if(parameters.nRuns > 1){computeConfidenceIntervals(GBVECS, 1.96, parameters.nRates);}

    auto tB = std::chrono::high_resolution_clock::now();
    auto bDuration = std::chrono::duration_cast<std::chrono::seconds>(tB - t1).count();
    cout << "CODE FINISHED RUNNING IN " << bDuration << " s TIME!" << endl;
    return EXIT_SUCCESS;
}