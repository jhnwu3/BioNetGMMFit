
#ifndef _FILEIO_HPP_
#define _FILEIO_HPP_
#include "main.hpp"
bool isNumber(const std::string& str);
string removeWhiteSpace(string current);
string findDouble(string line, int startPos);
MatrixXd txtToMatrix(const string& fileName, int rows, int cols);

// /* 
//     path - returns whatever.
//     returns a row x column matrix/vector from data of csv file.
// */
// template<typename M>
MatrixXd csvToMatrix(const std::string & path, int fileSize);
MatrixXd readX(const std::string &path, int xSize);
vector<MatrixXd> readY(const std::string & path, int ySize);
void matrixToCsv(const MatrixXd& mat, const string& fileName);
int readCsvPSO(int &nPart1, int &nSteps1, int &nPart2, int &nSteps2, int &useOnlySecMom, int &useOnlyFirstMom, int &useLinear, int &nRuns, int &simulateYt, int &useInverse, int &nRates, int &sampleSize);
int readCsvDataParam(int &nSpecies, int &nRates, int &xSize, int &ySize);
VectorXd readCsvTimeParam();
VectorXd readRates(int nRates);
#endif