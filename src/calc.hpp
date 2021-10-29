#ifndef _CALC_HPP_
#define _CALC_HPP_
#include "main.hpp"
double calculate_cf1(const VectorXd& trueVec, const VectorXd& estVec);
double calculate_cf2(const VectorXd& trueVec, const  VectorXd& estVec, const MatrixXd& w);
MatrixXd customWtMat(const MatrixXd& Yt, const MatrixXd& Xt, int nMoments, int N);
#endif _CALC_HPP_