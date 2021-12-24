#ifndef _CALC_HPP_
#define _CALC_HPP_
#include "main.hpp"
double calculate_cf1(const VectorXd& trueVec, const VectorXd& estVec);
double calculate_cf2(const VectorXd& trueVec, const  VectorXd& estVec, const MatrixXd& w);
MatrixXd ytWtMat(const MatrixXd& Yt, int nMoments, bool useBanks);
MatrixXd customWtMat(const MatrixXd& Yt, const MatrixXd& Xt, int nMoments, int N, bool useBanks, bool useInverse);
#endif 