#ifndef _GRAPH_HPP_
#define _GRAPH_HPP_
/*
Author: John Wu
Summary: Header file for all file input and output functions used to read in inputs from Config.csv and data/X or data/Y files.
 */

#include "main.hpp"

void graphMoments(const MatrixXd &X, const MatrixXd &Y, double t, const string &name);
void graphConfidenceIntervals(const MatrixXd &estimates, vector<string> &labels, const string &name);


#endif