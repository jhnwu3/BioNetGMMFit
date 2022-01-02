#include <iostream>
#include <fstream>
#include <boost/math/distributions.hpp>
#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>
#include <random>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>
#include <cmath>
#include <chrono>
#include <omp.h>
#include <Eigen/StdVector>
#include <boost/numeric/odeint/external/openmp/openmp.hpp>
#include "fileIO.hpp"

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;


bool isNumber(const std::string& s)
{
    int it = 0;
    while (it < s.length() && std::isdigit(s.at(it))) ++it; 
    return (it == (s.length() - 1));
}

bool isDouble(const std::string& s)
{
    int decimals = 0;
    int it = 0;
    for(int i = 0; i < s.length(); i++){
        if(s.at(i) == '.'){
            decimals++;
        }
    }
    if(decimals != 1){
        return false;
    }
    while (it < s.length() && std::isdigit(s.at(it))) ++it; 
    return (it == (s.length() - 2)); // -2, 1 for null, 1 for the .
}

string removeWhiteSpace(string current)
{
  string myNewString = "";
  string temp = "x";
  for (char c : current)
  {
    if (temp.back() != ' ' && c == ' ')
    {
      myNewString.push_back(' ');
    }
    temp.push_back(c);
    if (c != ' ')
    {
      myNewString.push_back(c);
    }
  }
  return myNewString;
}

string findDouble(string line, int startPos) {
    string doble;
    int i = startPos;
    int wDist = 0;
    while (i < line.length() && !isspace(line.at(i)) && line.at(i) != '\t') {
        i++;
        wDist++;
    }
    doble = line.substr(startPos, wDist);

    return doble;
}

MatrixXd txtToMatrix(const string& fileName, int rows, int cols) {
    MatrixXd mat(rows, cols);
    ifstream in(fileName);
    // use first row to determine how many columns to read.
    for (int i = 0; i < rows; i++) {
        string line;
        if (in.is_open()) {
            getline(in, line);
            line = removeWhiteSpace(line);
            int wordPos = 0;
            for (int j = 0; j < cols; j++) {
                string subs = findDouble(line, wordPos);
                mat(i, j) = stod(subs);
                wordPos += subs.length() + 1;
            }
        }
        else {
            cout << "Error: File Closed!" << endl;
        }

    }
    in.close();

    return mat;
}

MatrixXd csvToMatrix (const std::string & path) {
    std::ifstream indata;
    indata.open(path);
    if(!indata.is_open()){
        throw std::runtime_error("Invalid Sample File Name!");
        exit(EXIT_FAILURE);
    }
    std::string line;
    std::vector<double> values;
    unsigned int rows = 0;
    while (std::getline(indata, line)) {
        std::stringstream lineStream(line);
        std::string cell;
        while (std::getline(lineStream, cell, ',')) {
            values.push_back(std::stod(cell));
        }
        ++rows;
    }
    MatrixXd mat = MatrixXd::Zero(rows, values.size()/rows);
    int i = 0;
    for(int r = 0; r < rows; r++){
        for(int c = 0; c < mat.cols(); c++){
            mat(r,c) = values[i];
            i++;
        }
    }
    return mat;
}

void matrixToCsv(const MatrixXd& mat, const string& fileName){ // prints matrix to csv
    ofstream plot;
    string csvFile = fileName + ".csv";
	plot.open(csvFile);

    for(int i = 0; i < mat.rows(); i++){
        for(int j = 0; j < mat.cols(); j++){
            if(j == 0){
                plot << mat(i,j);
            }else{
                plot << "," << mat(i,j);
            }
        }
        plot << endl;
    }
    plot.close();
}

// Reads PSO Parameters File
int readCsvPSO(int &nPart1, int &nSteps1, int &nPart2, int &nSteps2, int &useOnlySecMom, int &useOnlyFirstMom, int &useLinear, int &nRuns){
    ifstream input("../PSO.csv");
    if(!input.is_open()){
        throw std::runtime_error("Could not open PSO file");
        return EXIT_FAILURE;
    }
    cout << "csvPSO" << endl;
    vector<int> params;
    string line;
    while(std::getline(input, line)){
        std::stringstream ss(line); // make a string stream from the line such that you can isolate each word even further.
        string col;
        while(std::getline(ss, col, ',')){
            cout << "col:" << col << "isNumber:" << isNumber(col) << endl;
            if(isNumber(col)){ // only add into parameter vector if actually an int.
                params.push_back(std::stoi(col)); 
            }
        }
    }
    cout << "params:" << params.size() << endl;
    nPart1 = params.at(0);
    nSteps1 = params.at(1);
    nPart2 = params.at(2);
    nSteps2 = params.at(3);
    useOnlySecMom = params.at(4);
    useOnlyFirstMom = params.at(5);
    useLinear = params.at(6);
    nRuns = params.at(7);
    input.close();
    return 0;
}

// Reads Input Data Parameters.
int readCsvDataParam(int &nSpecies, int &nRates){
    ifstream input("../data_parameters.csv");
    cout << "csvData" << endl;
    if(!input.is_open()){
        throw std::runtime_error("Could not open data parameters file");
        return EXIT_FAILURE;
    }
    vector<int> params;
    string line;
    while(std::getline(input, line)){
        std::stringstream ss(line); // make a string stream from the line such that you can isolate each word even further.
        string col;
        while(std::getline(ss, col, ',')){
            if(isNumber(col)){ // only add into parameter vector if actually an int.
                params.push_back(std::stoi(col)); 
            }
        }
    }
    nSpecies = params.at(0);
    nRates = params.at(1);
    input.close();
    return 0;
}

// Reads Time Step Parameters.
VectorXd readCsvTimeParam(){

    ifstream input("../time_steps.csv");
    if(!input.is_open()){
        throw std::runtime_error("Could not open data parameters file");
        exit;
    }
    vector<double> params;
    string line;
    while(std::getline(input, line)){
        std::stringstream ss(line); // make a string stream from the line such that you can isolate each word even further.
        string col;
        while(std::getline(ss, col, ',')){
            if(isNumber(col) || isDouble(col)){ // only add into parameter vector if actually an int.
                params.push_back(std::stod(col)); 
            }
        }
    }
    VectorXd times(params.size());
    for(int i = 0; i < params.size(); i++){
        times(i) = params.at(i);
    }
    input.close();
    return times;
}
