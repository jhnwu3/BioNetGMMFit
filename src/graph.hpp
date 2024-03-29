#ifndef _GRAPH_HPP_
#define _GRAPH_HPP_
/*
Author: John Wu
Summary: Functions Used to Generate Graphs. Note that every file path is essentially a vector of strings that gets read in and outputted by a python graphing script.
Python makes life easier.
 */

#include "main.hpp"
#include "fileIO.hpp"
class Grapher{
    public:
        string generalPath;
        string estFile;
        string modelName;
        string leastCostEstFile;
        string trueRatesFile;
        string forecastFile;
        string observedFile;
        // string contourFile; // these may not exist, but it just assumes a set of paths automatically
        vector<string> leastCostMoments;
        vector<string> allEstimatedMoments;
        vector<string> observedDataFiles;
        vector<string> estimatedDataFiles;
        vector<string> observedMoments; 
        Grapher(const string & path, const string & modelFileName, const string & trueRatesPath, const VectorXd &times, int nRates){
            generalPath = path + modelFileName;
            modelName = modelFileName;
            estFile = generalPath + "_estimates.csv";
            leastCostEstFile = generalPath + "_leastCostEstimate.csv";
            for(int t=1; t < times.size(); ++t){
                leastCostMoments.push_back(generalPath + "t" + to_string_with_precision(times(t), 2) + "_leastCostMoments.csv");
                observedDataFiles.push_back(generalPath +"Yt" + to_string_with_precision(times(t),2));
                estimatedDataFiles.push_back(generalPath + "Xt" + to_string_with_precision(times(t),2));
                allEstimatedMoments.push_back(generalPath + "XtMoments" + to_string_with_precision(times(t),2));
                observedMoments.push_back(generalPath + "YtMoments" + to_string_with_precision(times(t),2));
            }
            
            // contourFile = generalPath + "_contour";
            trueRatesFile = trueRatesPath;
            forecastFile = generalPath + "_forecast.csv";
            observedFile = generalPath + "_observed.csv";
        }

    void graphMoments(int nSpecies){
        int status = 0;
        string whichMoments = " -m " + to_string(nSpecies);

        for(int i = 0; i < leastCostMoments.size(); ++i){
            string cmd = "python3 graph.py -f " + leastCostMoments[i] + " -g Moments -n \"Fit of Predicted Moments\"" + whichMoments;
            status = system(cmd.c_str());
        }

        if (status < 0)
            std::cout << "Error: " << strerror(errno) << '\n';
        else
        {
            if (WIFEXITED(status))
                std::cout << "Program returned normally, exit code " << WEXITSTATUS(status) << '\n';
            else
                std::cout << "Program exited abnormaly\n";
        }
    }
    void graphConfidenceIntervals(bool simulated){
        string cmd = "python3 graph.py -f " + estFile + " -g CI";
        if(simulated){
            cmd = "python3 graph.py -f " + estFile + " -g CI_truth -r " + trueRatesFile + " -n \"Parameter Estimates\"";
        }
        int status = system(cmd.c_str());
        if (status < 0)
            std::cout << "Error: " << strerror(errno) << '\n';
        else
        {
            if (WIFEXITED(status))
                std::cout << "Program returned normally, exit code " << WEXITSTATUS(status) << '\n';
            else
                std::cout << "Program exited abnormaly\n";
        }
    }

    void graphForecasts(int nSpecies){
        string whichMoments = " -m " + to_string(nSpecies);
        string cmd = "python3 graph.py -f " + forecastFile + " " + observedFile + " -g forecast" + whichMoments;

        int status = system(cmd.c_str());
        if (status < 0)
            std::cout << "Error: " << strerror(errno) << '\n';
        else
        {
            if (WIFEXITED(status))
                std::cout << "Program returned normally, exit code " << WEXITSTATUS(status) << '\n';
            else
                std::cout << "Program exited abnormaly\n";
        }
    }

    void graphContours(int nRates, const string & contourFile){
        string whichMoments = " -m " + to_string(nRates);
        string titleFig = " -n " + modelName;
        string cmd = "python3 graph.py -f "+ contourFile + " -g contour" + whichMoments + titleFig;

        int status = system(cmd.c_str());
        if (status < 0)
            std::cout << "Error: " << strerror(errno) << '\n';
        else
        {
            if (WIFEXITED(status))
                std::cout << "Program returned normally, exit code " << WEXITSTATUS(status) << '\n';
            else
                std::cout << "Program exited abnormaly\n";
        }
    }
    
};


#endif