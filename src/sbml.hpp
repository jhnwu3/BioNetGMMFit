
#ifndef _SBML_HPP_
#define _SBML_HPP_
#include "main.hpp"
#include "nonlinear.hpp"
#include "tinyxml2.h"

VectorXd simulateSBML(int useDet, double ti, double tf, const VectorXd &c0, const VectorXd &k);
vector<string> getSpeciesNames(const string& path);
vector<int> specifySpeciesFromProteinsList(const string& path, vector<string> &species, int nObs);
#endif
