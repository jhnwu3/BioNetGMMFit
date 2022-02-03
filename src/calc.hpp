#ifndef _CALC_HPP_
#define _CALC_HPP_
#include "main.hpp"
bool isInvertible(const MatrixXd& m);
double rndNum(double low, double high);
double costFunction(const VectorXd& trueVec, const  VectorXd& estVec, const MatrixXd& w);
MatrixXd ytWtMat(const MatrixXd& Yt, int nMoments, bool useBanks);
MatrixXd customWtMat(const MatrixXd& Yt, const MatrixXd& Xt, int nMoments, int N, bool useBanks);
#endif 