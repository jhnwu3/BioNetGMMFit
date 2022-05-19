#ifndef _LINEAR_HPP_
#define _LINEAR_HPP_
#include "main.hpp"
#include "calc.hpp"
#include "fileIO.hpp"
#include "system.hpp"
#include "cli.hpp"

VectorXd momentVector(const MatrixXd &sample, int nMoments);
MatrixXd evolutionMatrix(VectorXd &k, double tf, int nSpecies);
VectorXd linearVelVec(const VectorXd& posK, int seed, double epsi, double nan, int hone);
MatrixXd linearModel(int nParts, int nSteps, int nParticles2, int nSteps2, MatrixXd& X_0, int nRates, int nMoments, const VectorXd &times, int simulateYt, int useInverse, int argc, char** argv, int rngSeed);

#endif