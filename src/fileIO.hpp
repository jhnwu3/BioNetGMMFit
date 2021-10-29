
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
// M load_csv (const std::string & path) {
//     std::ifstream indata;
//     indata.open(path);
//     std::string line;
//     std::vector<double> values;
//     uint rows = 0;
//     while (std::getline(indata, line)) {
//         std::stringstream lineStream(line);
//         std::string cell;
//         while (std::getline(lineStream, cell, ',')) {
//             values.push_back(std::stod(cell));
//         }
//         ++rows;
//     }
//     return Map<const Matrix<typename M::Scalar, M::RowsAtCompileTime, M::ColsAtCompileTime, RowMajor>>(values.data(), rows, values.size()/rows);
// }

void matrixToCsv(const MatrixXd& mat, const string& fileName);
int readCsvPSO(int &nPart1, int &nSteps1, int &nPart2, int &nSteps2, int &useMixMom, int &useLinear);
int readCsvDataParam(int &xDataSize, int &yDataSize, int &nSpecies, int &nRates);
VectorXd readCsvTimeParam();
#endif