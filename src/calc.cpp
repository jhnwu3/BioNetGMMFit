#include "calc.hpp"
#include "fileIO.hpp"

/* Cost Function, by default with an identity weight matrix = square of differences, however, formally it is in the form of:
    (true vector - estimated vector)' * weight * (true vector - estimated vector)
 */

bool isInvertible(const MatrixXd& m){
    return (m.determinant() != 0);
}

double rndNum(double low, double high){
    random_device RanDev;
    mt19937 gen(RanDev());
    uniform_real_distribution<double> unifDist(low, high);
    return unifDist(gen);
}


double costFunction(const VectorXd& trueVec, const  VectorXd& estVec, const MatrixXd& w) {
    double cost = 0;
    VectorXd diff(trueVec.size());
    diff = trueVec - estVec;
    cost = diff.transpose() * w * (diff.transpose()).transpose();
    return cost;
}
/*TODO: Rename to wolfe weights */
MatrixXd wolfWtMat(const MatrixXd& Yt, int nMoments, bool useInverse){
    /* first moment differences */
    MatrixXd fmdiffs = MatrixXd::Zero(Yt.rows(), Yt.cols());
    for(int i = 0; i < Yt.cols(); i++){
        fmdiffs.col(i) = Yt.col(i).array() - Yt.col(i).array().mean();
    }
    /* second moment difference computations - @todo make it variable later */
    MatrixXd smdiffs(Yt.rows(), Yt.cols());
    for(int i = 0; i < Yt.cols(); i++){
        smdiffs.col(i) = (Yt.col(i).array() * Yt.col(i).array()) - (Yt.col(i).array().mean() * Yt.col(i).array().mean());
    }
    /* If no cross moments, then have a check for it */
    int nCross = nMoments - (2 * Yt.cols());
    if (nCross < 0){
        nCross = 0;
    }
    MatrixXd cpDiff(Yt.rows(), nCross);

    /* cross differences */
    if(nCross > 0){
        int upperDiag = 0;
        for(int i = 0; i < Yt.cols(); i++){
            for(int j = i + 1; j < Yt.cols(); j++){
                cpDiff.col(upperDiag) = (Yt.col(i).array() * Yt.col(j).array()) - (Yt.col(i).array().mean() * Yt.col(j).array().mean());
                upperDiag++;
            }
        }
    }

    MatrixXd aDiff(Yt.rows(), nMoments);
    for(int i = 0; i < Yt.rows(); i++){
        for(int moment = 0; moment < nMoments; moment++){
            if(moment < Yt.cols()){
                aDiff(i, moment) = fmdiffs(i, moment);
            }else if (moment >= Yt.cols() && moment < 2 * Yt.cols()){
                aDiff(i, moment) = smdiffs(i, moment - Yt.cols());
            }else if (moment >= 2 * Yt.cols()){
                aDiff(i, moment) = cpDiff(i, moment - (2 * Yt.cols()));
            }
        }
    }
    double cost = 0;
    VectorXd variances(nMoments);
    if(aDiff.rows() > 1){ // are there 2 or more cells? Ok, we can compute a variance, otherwise, let's default to an identity matrix.
        for(int i = 0; i < nMoments; i++){
            variances(i) = (aDiff.col(i).array() - aDiff.col(i).array().mean()).square().sum() / ((double) aDiff.col(i).array().size() - 1);
        }
    }else{
        for(int i = 0; i < nMoments; i++){
            variances(i) = 0;
        }
    }
  
    MatrixXd wt = MatrixXd::Zero(nMoments, nMoments);
    
    if(useInverse){
         // compute covariances for differences.
        for(int i = 0; i < nMoments; i++){
            wt(i,i) = variances(i); // cleanup code and make it more vectorized later.
        }
        for(int i = 0; i < nMoments; i++){
            for(int j = i + 1; j < nMoments; j++){
                wt(i,j) = ((aDiff.col(i).array() - aDiff.col(i).array().mean()).array() * (aDiff.col(j).array() - aDiff.col(j).array().mean()).array() ).sum() / ((double) aDiff.col(i).array().size() - 1); 
                wt(j,i) = wt(i,j); // across diagonal
            }
        }
        // [2 0 <- variances on diagonal
        //  0 4]
        // [1/2 0 <- inverse of that
        //   0  1/4 ]
        // [2 1
        //  -1 4 ]
        //  [0.asdasd  0.asddas ]
        //    0.asdas
        // 
        wt = wt.colPivHouseholderQr().solve(MatrixXd::Identity(wt.cols(), wt.cols()));

    }else{
        for(int i = 0; i < nMoments; i++){
            if(variances(i) == 0){ variances(i) = 1; } // error check for invalid variancees
            wt(i,i) = 1.0 / variances(i); // cleanup code and make it more vectorized later.
        }
    }
    
    return wt;
}

/* TODO: Rename to Das Weights */
MatrixXd dasWtMat(const MatrixXd& Yt, const MatrixXd& Xt, int nMoments, int N, bool useInverse){
    if(Yt.rows() != Xt.rows() || Yt.cols() != Xt.cols()){
        cout << "Error! Dimension mismatch between X and Y! Calculation of Das Weights cancelled!" << endl;
        return MatrixXd::Identity(nMoments, nMoments);
    }
    /* first moment differences */
    MatrixXd fmdiffs = Yt - Xt; 
    /* second moment difference computations - @todo make it variable later */
    MatrixXd smdiffs(N, Yt.cols());
    for(int i = 0; i < Yt.cols(); i++){
        smdiffs.col(i) = (Yt.col(i).array() * Yt.col(i).array()) - (Xt.col(i).array() * Xt.col(i).array());
    }
    /* If no cross moments, then have a check for it */
    int nCross = nMoments - 2 * Yt.cols();
    if (nCross < 0){
        nCross = 0;
    }
    MatrixXd cpDiff(N, nCross);
    
    /* cross differences */
    if(nCross > 0){
        int upperDiag = 0;
        for(int i = 0; i < Yt.cols(); i++){
            for(int j = i + 1; j < Yt.cols(); j++){
                cpDiff.col(upperDiag) = (Yt.col(i).array() * Yt.col(j).array()) - (Xt.col(i).array() * Xt.col(j).array());
                upperDiag++;
            }
        }
    }
    MatrixXd aDiff(N, nMoments);
    for(int i = 0; i < N; i++){
        for(int moment = 0; moment < nMoments; moment++){
            if(moment < Yt.cols()){
                aDiff(i, moment) = fmdiffs(i, moment);
            }else if (moment >= Yt.cols() && moment < 2 * Yt.cols()){
                aDiff(i, moment) = smdiffs(i, moment - Yt.cols());
            }else if (moment >= 2 * Yt.cols()){
                aDiff(i, moment) = cpDiff(i, moment - (2 * Yt.cols()));
            }
        }
    }
    
    MatrixXd wt = MatrixXd::Zero(nMoments, nMoments);
    if(useInverse){
        // compute covariances for differences.
        VectorXd variances(nMoments);
        for(int i = 0; i < nMoments; i++){
            variances(i) = (aDiff.col(i).array() - aDiff.col(i).array().mean()).square().sum() / ((double) aDiff.col(i).array().size() - 1);
        }
        for(int i = 0; i < nMoments; i++){
            wt(i,i) = variances(i); // cleanup code and make it more vectorized later.
        }
        for(int i = 0; i < nMoments; i++){
            for(int j = i + 1; j < nMoments; j++){
                wt(i,j) = ((aDiff.col(i).array() - aDiff.col(i).array().mean()).array() * (aDiff.col(j).array() - aDiff.col(j).array().mean()).array() ).sum() / ((double) aDiff.col(i).array().size() - 1); 
                wt(j,i) = wt(i,j); // across diagonal
            }
        }
        wt = wt.completeOrthogonalDecomposition().solve(MatrixXd::Identity(nMoments, nMoments));
    }else{
        
        VectorXd variances(nMoments);
        for(int i = 0; i < nMoments; i++){
            variances(i) = (aDiff.col(i).array() - aDiff.col(i).array().mean()).square().sum() / ((double) aDiff.col(i).array().size() - 1);
        }
        for(int i = 0; i < nMoments; i++){
            wt(i,i) = 1 / variances(i); // cleanup code and make it more vectorized later.
        }
        cout << "Weights:"<< endl;
        cout << wt << endl;
        
    }
    
    return wt;
}

MatrixXd bootStrap(const MatrixXd& sample){

    MatrixXd bSample = MatrixXd::Zero(sample.rows(), sample.cols());
    std::random_device rand_dev;
    std::mt19937 generator(rand_dev());
    std::uniform_int_distribution<> unif(0, sample.rows() - 1);

    /* Now resample to create boostrapped data sample */
    for(int i = 0; i < sample.rows(); ++i){
        bSample.row(i) = sample.row(unif(generator));
    }
    return bSample;
}

VectorXd cwiseVar(const MatrixXd& sample){
    VectorXd variances(sample.cols());
    for(int c = 0; c < sample.cols(); ++c){
        variances(c) = (sample.col(c).array() - sample.col(c).array().mean()).square().sum() / ((double) sample.col(c).array().size() - 1);
    }
    return variances;
}

void computeConfidenceIntervals(const MatrixXd& sample, double z, int nRates){
    cout << "------- 95 Percent Confidence Intervals -------" << endl;
    /* Cheap Way to Compute Means and Standard Deviation */
    VectorXd estMu = sample.colwise().mean();
    VectorXd estSigma = cwiseVar(sample).array().sqrt();
    cout << "Rates | Standard Deviation" << endl;
    for(int r = 0; r < nRates; ++r){
        cout << estMu(r) << "   |   " << estSigma(r) << endl;
    }
    VectorXd delta = estSigma / sqrt(sample.rows()); 
    cout << "Confidence Intervals for Each Rate:" << endl;
    for(int r = 0; r < nRates; ++r){
        cout << "Theta" << r << ": [" << estMu(r) - z*delta(r) << "," <<estMu(r) + z * delta(r) << "]" << endl;
    }
    cout << "-----------------------------------------------" << endl;
}

bool rowIsAllPositive(const VectorXd &x){
    bool allPos = true;
    for(int i = 0; i < x.size(); i++){
        if (x(i) <  0){ //|| x(0) > 60 || x(1) > 70 || x(2) > 25 || x(3) > 70
            return false;
        }
    }
    return allPos;
}

MatrixXd filterZeros(const MatrixXd &X){
    MatrixXd x_filtered = MatrixXd::Zero(0,0);
    for(int i = 0; i < X.rows(); i++){
        if(rowIsAllPositive(X.row(i))){
            x_filtered.conservativeResize(x_filtered.rows() + 1, X.cols());
            x_filtered.row(x_filtered.rows() - 1) = X.row(i);
        }
    }
    return x_filtered;
}
// pos should be the ground truth parameter set that remains.

// MatrixXd generatePairwiseContour(const RoadRunner &model, SimulateOptions &opt, const VectorXd &pos, const MatrixXd &x0, const vector<VectorXd> &yt, const vector<MatrixXd> &wts, int theta1, int theta2, int stepSize){


//     MatrixXd
// }