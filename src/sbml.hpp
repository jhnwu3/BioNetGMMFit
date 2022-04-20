
#ifndef _SBML_HPP_
#define _SBML_HPP_
#include "main.hpp"
#include "nonlinear.hpp"
VectorXd simulateSBML(int useDet, double ti, double tf, const VectorXd &c0, const VectorXd &k);
#endif
