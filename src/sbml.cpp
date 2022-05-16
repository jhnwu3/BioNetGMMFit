#include "sbml.hpp"

VectorXd simulateSBML(int useDet, double ti, double tf, const VectorXd &c0, const VectorXd &k){
    RoadRunner r = RoadRunner("sbml/model_sbml.xml");
    vector<double> init = convertInit(c0);
    r.changeInitialConditions(init);
    SimulateOptions opt;
    opt.start = ti;
    opt.duration = tf;
    opt.steps = 2;

    // to apparently change the values of parameters in the model, we must first feed the vector into a double array.
    double kVals[k.size()];
    for(int i = 0; i < k.size(); ++i){
        kVals[i] = k(i);
    }
    r.getModel()->setGlobalParameterValues(k.size(),0,kVals); // set new global parameter values here.
    if(useDet > 0){
        r.setIntegrator("cvode");
    }else{
        r.setIntegrator("gillespie");
    }
    const DoubleMatrix res = *r.simulate(&opt);
    VectorXd evolved = VectorXd::Zero(res.numCols() - 1);
    for(int i = 1; i < res.numCols(); i++){
        evolved(i - 1) = res[res.numRows() - 1][i];
    }
    return evolved;
}
