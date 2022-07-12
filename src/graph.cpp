#include "graph.hpp"
#include "fileIO.hpp"



void graphMoments(const MatrixXd &X, const MatrixXd &Y, double t, vector<string> &labels, const string &name){
    matrixToCsvWithLabels(X, labels, "frontend/graph/" + name + "_X_t" + to_string(t) + ".csv"); // graphing Xt vs. Yt for a specific time point
    matrixToCsvWithLabels(Y, labels, "frontend/graph/" + name + "_Y_t" + to_string(t) + ".csv"); 
    string pythonCall = "python3 graph.py"
}
void graphConfidenceIntervals(const MatrixXd &estimate, vector<string> &labels, const string &name){
    matrixToCsvWithLabels(estimate, labels, "frontend/graph/" + name + "_estimates.csv"); 
}
