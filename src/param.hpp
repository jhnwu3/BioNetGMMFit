#ifndef _PARAM_HPP_
#define _PARAM_HPP_
#include "main.hpp"
#include "fileIO.hpp"
//
class Parameters{
    public:
        int nParts; // first part PSO
        int nSteps;
        int nParts2; // second part PSO
        int nSteps2;
        int useOnlySecMom;
        int useOnlyFirstMom;
        int useLinear;
        int nRates;
        int nRuns;
        int simulateYt;
        int useInverse; // currently just inverse only occurs in linear model.
        int heldTheta;
        int reportMoments;
        int nest;
        int bootstrap;
        int useSBML;
        int useDet;
        int odeSteps;
        int seed;
        int nThreads;
        int nMoments;
        double heldThetaVal;
        double hyperCubeScale;
        string outPath;
        double pBestWeight;
        double globalBestWeight;
        double pInertia;
        Parameters(const string &path){
            cout << "Reading in Parameters from Configuration File!" << endl;
            std::ifstream input(path);
            if(!input.is_open()){
                throw std::runtime_error("Could Not Open Configuration file! Please specify a configuration file by \"./BNGMM -c /path/to/config.csv\"! Call \"./BNGMM -h\" for more information!");
                exit(EXIT_FAILURE);
            }
            vector<double> params;
            string line;
            while(std::getline(input, line)){
                std::stringstream ss(line); // make a string stream from the line such that you can isolate each word even further.
                string col;
                while(std::getline(ss, col, ',')){
                    if(isNumber(col) || isDouble(col)){ // only add into parameter vector if actually an int.
                        params.push_back(std::stod(col)); 
                    }
                }
            }
            nParts = params.at(0);
            nSteps = params.at(1);
            useOnlySecMom = params.at(2);
            useOnlyFirstMom = params.at(3);
            nRuns = params.at(4);
            simulateYt = params.at(5);
            useInverse = params.at(6);
            nRates = params.at(7);
            heldTheta = params.at(8);
            heldThetaVal = params.at(9);
            hyperCubeScale = params.at(10);
            reportMoments = params.at(11);
            bootstrap = params.at(12);
            useDet = params.at(13);
            odeSteps = params.at(14);
            seed = params.at(15);
            nThreads = params.at(16); 
            pBestWeight = params.at(17);
            globalBestWeight = params.at(18);
            pInertia = params.at(19);
            
            useSBML = 0;
            outPath = "";
            input.close();
        }
        void printParameters(int nMoments, const VectorXd& times){
            cout << "---------------------  Parameters  --------------------" << endl;
            if(useLinear){
                cout << "Using Matrix Interaction Model Instead of Runge Kutta Solvers!" << endl;
            }
            cout << "Total Number of Runs:" << nRuns << endl;
            cout << "Number of Moments:" << nMoments << endl;
            if(useOnlyFirstMom){
                cout << "Using Only Means!" << endl;
            }else if(useOnlySecMom){
                cout << "Using Only Means and Second Moments!" << endl;
            }
            if(heldTheta > - 1){
                cout << "Theta Held Index:" << heldTheta << " held value:" << heldThetaVal << endl;
            }
            cout << "Hyper Cube Width:" << hyperCubeScale << endl;
            cout << "Using Times:" << times.transpose() << endl;
            cout << "Blind PSO --> nParts:" << nParts << " Nsteps:" << nSteps << endl;
            cout << "Targeted PSO --> nParts:" <<  nParts2 << " Nsteps:" << nSteps2 << endl;
            cout << "Number of Rates:" << nRates << endl;
            cout << "Particle Best Weight:" << pBestWeight << " Global Best Weight:"<< globalBestWeight << " Particle Inertia:" << pInertia << endl;
            if(useSBML){
                cout << "Redirecting Model to SBML/BNGL" << endl;
                if(useDet > 0){
                    cout << "Modeling With Deterministic ODEs" << endl;
                }else{
                    cout << "Modeling with Gillespie" << endl;
                }
                cout << "Number of Steps of Integration Determined:" << odeSteps << endl;
            }
            if(useInverse){
                cout << "Using Matrix Inverse!" << endl;
            }
            if(seed > 0){
                cout << "Now seeding with value " << seed << endl;
            }
            if(bootstrap > 0){
                cout << "Enabled Bootstrapping!" << endl;
            }
            cout << "Outputting data to:" << outPath << endl;
            cout << "--------------------------------------------------------" << endl;
        }
};
#endif