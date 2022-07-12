#ifndef _FILEIO_HPP_
#define _FILEIO_HPP_
/*
Author: John Wu
Summary: Header file for all file input and output functions used to read in inputs from Config.csv and data/X or data/Y files.
 */

#include "main.hpp"

/* String Processing Functions */
bool isNumber(const std::string& str); 
bool isDouble(const std::string& s);
string removeWhiteSpace(string current);
string findDouble(string line, int startPos);

/* General Input File Functions*/
MatrixXd txtToMatrix(const string& fileName, int rows, int cols);
MatrixXd csvToMatrix(const std::string & path);

/* Specific Input File Functions for data/X and data/Y */
MatrixXd readX(const std::string &path);
vector<MatrixXd> readY(const std::string & path);

/* General Output Functions*/
void matrixToCsv(const MatrixXd& mat, const string& fileName);
void matrixToCsvWithLabels(const MatrixXd& mat,  vector<string> &labels, const string& fileName);
/* Reading in time and rate parameters */
VectorXd readCsvTimeParam(const string &path);
VectorXd readRates(int nRates, const string &path);

#endif