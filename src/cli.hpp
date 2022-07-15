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
#endif

