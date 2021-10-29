#include "calc.hpp"

double calculate_cf1(const VectorXd& trueVec, const VectorXd& estVec) {
    double cost = 0;
    VectorXd diff(trueVec.size());
    diff = trueVec - estVec;
    cost = diff.transpose() * diff.transpose().transpose();
    return cost;
}
double calculate_cf2(const VectorXd& trueVec, const  VectorXd& estVec, const MatrixXd& w) {
    double cost = 0;
    VectorXd diff(trueVec.size());
    diff = trueVec - estVec;
    cost = diff.transpose() * w * (diff.transpose()).transpose();
    return cost;
}


MatrixXd customWtMat(const MatrixXd& Yt, const MatrixXd& Xt, int nMoments, int N){
    
    bool useInverse = true;
    if(Yt.cols() > 4){
        useInverse = false;
    }
    /* first moment differences */
    MatrixXd fmdiffs = Yt - Xt; 
    /* second moment difference computations - @todo make it variable later */
    MatrixXd smdiffs(N, Yt.cols());
    for(int i = 0; i < smdiffs.cols(); i++){
        smdiffs.col(i) = (Yt.col(i).array() - Yt.col(i).array().mean()) * (Yt.col(i).array() - Yt.col(i).array().mean()) - ((Xt.col(i).array() -Xt.col(i).array().mean()) *(Xt.col(i).array() - Xt.col(i).array().mean()));
    }

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
                cpDiff.col(upperDiag) = ((Yt.col(i).array() - Yt.col(i).array().mean()) * (Yt.col(j).array() -Yt.col(j).array().mean())) - ((Xt.col(i).array() - Xt.col(i).array().mean()) * (Xt.col(j).array() - Xt.col(j).array().mean()));
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
    double cost = 0;
    VectorXd means = aDiff.colwise().mean();
    VectorXd variances(nMoments);
    for(int i = 0; i < nMoments; i++){
        variances(i) = (aDiff.col(i).array() - aDiff.col(i).array().mean()).square().sum() / ((double) aDiff.col(i).array().size() - 1);
    }
    int rank = nMoments;
    MatrixXd wt = MatrixXd::Identity(rank, rank);
    if(!useInverse){
        for(int i = 0; i < rank; i++){
            wt(i,i) = 1 / variances(i); // cleanup code and make it more vectorized later.
        }
    }else{
        for(int i = 0; i < aDiff.rows(); i++){
            wt += aDiff.row(i).transpose() * aDiff.row(i);
        }
        wt = (wt / aDiff.rows()).inverse();
    }

    cout << "new wt mat:" << endl << wt << endl;

    return wt;
}