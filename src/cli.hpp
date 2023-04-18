#ifndef _CLI_HPP_
#define _CLI_HPP_
#include "main.hpp"
string getConfigPath(int argc, char **argv);
string getTimeStepsPath(int argc, char **argv);
string getTrueRatesPath(int argc, char **argv);
string getXPath(int argc, char **argv);
string getYPath(int argc, char **argv);
string getModelPath(int argc, char **argv);
string getProPath(int argc, char **argv);
string getOutputPath(int argc, char **argv);
bool modelPathExists(int argc, char **argv);
bool outPathExists(int argc, char **argv);
bool graphingEnabled(int argc, char **argv);
bool proPathExists(int argc, char **argv);
bool helpCall(int argc, char **argv);
bool holdRates(int argc, char **argv);
bool seedRates(int argc, char **argv);
bool forecast(int argc, char **argv);
bool contour(int argc, char **argv);

string getContourTheta1(int argc, char**argv);
string getContourTheta2(int argc, char**argv);

string getHeldRatesDir(int argc, char **argv);
string getSeededRates(int argc, char **argv);
string getForecastedTimes(int argc, char **argv);
#endif

