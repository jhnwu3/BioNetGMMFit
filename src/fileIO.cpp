#include "fileIO.hpp"

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

/* Returns true if string is an integer*/
bool isNumber(const std::string& s)
{
    int it = 0;
    bool allDigits = true;
    while (it < s.length()){
        if (!(std::isdigit(s.at(it))) && !(s.at(it) == '\r') && !(s.at(it) == '\0') && !(s.at(it)==' ') && !(s.at(it)=='-')){ // check for control characters to check for numbers correctly.
            allDigits = false;
        }
        ++it; 
    } 
    return allDigits;
}

/* Returns true if string is a Double */
bool isDouble(const std::string& s)
{
    bool allDigits = true;
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
     while (it < s.length()){
        if (!(std::isdigit(s.at(it))) && !(s.at(it) == '\r') && !(s.at(it) == '\0') && !(s.at(it) == '.')  && !(s.at(it)==' ') && !(s.at(it)=='-')){ // check for control characters to check for numbers correctly.
            allDigits = false;
        }
        ++it; 
    } 
    return allDigits; // -2, 1 for null, 1 for the .
}

/* Remove whitespace from string */
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

/* Finds next double in string position */
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
string vectorToString(const VectorXd & v){
    string vec = "";
    for(int t = 0; t < v.size(); t++){
        vec += "_" + to_string(((int) v(t)));
    }
    return vec;
}
string getFileNameWithoutExtensions(const string & path){
    string baseModelFile = path.substr(path.find_last_of("/\\") + 1);
    string::size_type const p(baseModelFile.find_last_of('.'));
    string file_without_extension = baseModelFile.substr(0, p);
    return file_without_extension;
}
/* 
    Summary:
        Takes a txt file and converts it into a matrix 
    Input:
        fileName - name of file (with file extension) that's being converted
        rows - number of rows 
        cols - No. of Cols
    Output:
        A matrix of rows x cols defined in inputs.
    
*/
MatrixXd txtToMatrix(const string& fileName, int rows, int cols) {
    MatrixXd mat(rows, cols);
    std::ifstream in(fileName);
    // use first row to determine how many columns to read.
    for (int i = 0; i < rows; i++) {
        string line;
        if (in.is_open()) {
            std::getline(in, line);
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

/* 
    Summary:
        Same idea as above, but now can read an entire file and no longer need to specifiy number of columns.
    Input:
        path - file name with directory including the .csv extension
        fileSize - number of rows of matrix being outputted
    Output:
        Matrix of file with only the number of rows of fileSize.

*/
MatrixXd csvToMatrix (const std::string & path){
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
    indata.close();
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

/* 
    Summary:
        Reads a single file from the data/X directory into a matrix.
    Input:
        path - path of X directory
        xSize - number of rows of matrix to be returned
    Output:
        matrix of xSize rows.
*/
MatrixXd readX(const std::string &path){
    int nFile = 0;
    MatrixXd X_0;
    cout << "------ Reading in X_0! ------" << endl;
    try{
        for(const auto & entry : fs::directory_iterator(path)){
            cout << entry << endl;
            X_0 = csvToMatrix(entry.path().string());
            ++nFile;
        }
    }catch(...){
        cout << "Error with X directory path, please make sure to input a valid path that is a DIRECTORY not a file, received path:" << path << endl;
        exit(-1);
    }
    if(nFile < 1){
        cout << "Error! No X csv file detected in " << path << "!" << endl;
        exit(-1);
    }
    if(nFile > 1){
        cout << "Multiple X files detected, reading only the last X file read in." << endl;
    }
    cout << "Reading in (rows,columns): (" << X_0.rows() <<"," << X_0.cols() << ") from X data directory" << endl;
    X_0 = filterZeros(X_0);
    cout << "After removing all negative rows, X has " << X_0.rows() << " rows." << endl;
    cout << "If dimensions are unexpected of input data, please make sure to check csv files contains all possible values in each row/column." << endl;
    cout << "---------------------------" << endl;
    return X_0;
}
/*
    Summary:
        Reads all files in data/Y directory into a list of matrices
    Input:
        path - path of Y directory
        ySize - number of rows of matrices to be returned
    Output:
        vector of matrices of ySize. note: vector == list.
*/
vector<MatrixXd> readY(const std::string & path){
    vector<MatrixXd> Y;
    cout << "------ Reading in Yt! ------" << endl;
    try{
        for(const auto & entry : fs::directory_iterator(path)){
            cout << entry << endl;
            Y.push_back(csvToMatrix(entry.path().string()));
            cout << "Read in " <<Y[Y.size() - 1].rows() << " rows! " << Y[Y.size() - 1].cols() << " columns!"<< endl;
        }
    }catch(...){
        cout << "Error! Unable to Read in values from the true Y directory. Make sure to specifiy a directory path not an explicit file name. Error with path:"<< path << endl;
        exit(-1);
    }
    if(Y.size() < 1){
        cout << "Error! 0 Y Files read in!" << endl;
        exit(-1);
    }
    cout << "---------------------------" << endl;
    return Y;
}
/* 
    Summary:
        Converts a matrix mat into a csv file, defined by the string variable fileName.
 */
void matrixToCsv(const MatrixXd& mat, const string& fileName){ // prints matrix to csv
    std::ofstream plot;
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

void matrixToCsvWithLabels(const MatrixXd& mat,  vector<string> &labels, const string& fileName){ // prints matrix to csv
    std::ofstream plot;
    string csvFile = fileName + ".csv";
	plot.open(csvFile);
    for(int j = 0; j < mat.cols(); j++){
        if(j == 0){
            plot << labels[j];
        }else{
            plot << "," << labels[j];
        }
    }
    plot << endl;
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
void vectorToCsv(const VectorXd& v, const string& fileName){
    std::ofstream plot;
    string csvFile = fileName + ".csv";
	plot.open(csvFile);
    plot << v << endl;
    plot.close();
}
void reportLeastCostMoments(const VectorXd & est, const VectorXd & obs, double t, const string& fileName){
    std::ofstream plot;
    string csvFile = fileName + "t" + to_string_with_precision(t, 2)+ "_leastCostMoments.csv";
	plot.open(csvFile);
    plot << "Estimated,Observed" << endl;
    for(int i = 0; i < est.size(); ++i){
        plot << est(i) << "," << obs(i) << endl;
    } 
    plot.close();
}
void reportAllMoments(vector<MatrixXd> & x, vector<VectorXd> & y, const VectorXd& times, const string& fileName){

    for(int t = 0; t < times.size() - 1; ++t){
        matrixToCsv(x[t], fileName + "XtMoments" + to_string_with_precision(times(t + 1),2));
        vectorToCsv(y[t], fileName + "YtMoments" + to_string_with_precision(times(t + 1),2));
    } 

}
// Reads Time Step Parameters.
VectorXd readCsvTimeParam(const string &path){

    std::ifstream input(path);
    if(!input.is_open()){
        throw std::runtime_error("Could not open time steps csv file!");
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
    if(times.size() < 2){
        cout << "Error! Unable to read in timesteps properly, need two or more time steps!" << endl;
        exit(EXIT_FAILURE);
    }
    return times;
}

VectorXd readRates(int nRates, const string &path){
    std::ifstream input(path);
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
    VectorXd rates(nRates);
    if(params.size() < nRates){
        cout << "Error, number of rates in parameters do not match number of true rates simulated!" << endl;
        exit(1);
    }
    for(int i = 0; i < nRates; i++){
        rates(i) = params.at(i);
    }
    input.close();
    return rates;
}


MatrixXd heldThetas(int nRates, const string& path){
    cout << path << endl;
    std::ifstream input(path);
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
                // cout << std::stod(col) << endl;
            }
        }
    }
    cout << "HI" << endl;
    MatrixXd rates = MatrixXd::Zero(nRates, 2);
    if(params.size() % 2 != 0){
        cout << "Error, missing some columns in held rates file!" << endl;
        exit(1);
    }
    // go every other element 
    for(int i = 0; i < params.size(); i++){
        if (i % 2 == 0){
            rates(int(params.at(i)), 0) = 1;
        }else{
            rates(int(params.at(i - 1)), 1) = params.at(i);
        }
    }

    input.close();
    return rates;
}

VectorXd readSeed(int nRates, const string& path){
    std::ifstream input(path);
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
    VectorXd rates(nRates);
    if(params.size() < nRates){
        cout << "Error, number of rates in parameters do not match number of true rates simulated!" << endl;
        exit(1);
    }
    for(int i = 0; i < nRates; i++){
        rates(i) = params.at(i);
    }
    input.close();
    return rates;
}