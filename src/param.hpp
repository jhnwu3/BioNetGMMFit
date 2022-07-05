#ifndef _PARAM_HPP_
#define _PARAM_HPP_
#include "main.hpp"
#include "fileIO.hpp"
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
        double heldThetaVal;
        double hyperCubeScale;
        Parameters(const string &path){
            cout << "Reading in Parameters from Configuration File!" << endl;
            std::ifstream input(path);
            if(!input.is_open()){
                throw std::runtime_error("Could Not Open Configuration file");
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
            nParts2 = params.at(2);
            nSteps2 = params.at(3);
            useOnlySecMom = params.at(4);
            useOnlyFirstMom = params.at(5);
            useLinear = params.at(6);
            nRuns = params.at(7);
            simulateYt = params.at(8);
            useInverse = params.at(9);
            nRates = params.at(10);
            heldTheta = params.at(11);
            heldThetaVal = params.at(12);
            hyperCubeScale = params.at(13);
            reportMoments = params.at(14);
            nest = params.at(15);
            bootstrap = params.at(16);
            useSBML = params.at(17);
            useDet = params.at(18);
            odeSteps = params.at(19);
            seed = params.at(20);
            nThreads = params.at(21); 
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
            if(useSBML){
                cout << "Redirecting Model to SBML" << endl;
                if(useDet > 0){
                    cout << "and Modeling With Deterministic ODEs" << endl;
                }else{
                    cout << "using gillespie" << endl;
                }
            }
            if(seed > 0){
                cout << "Now seeding with value " << seed << endl;
            }
            cout << "--------------------------------------------------------" << endl;
        }
};
#endif