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

MatrixXd ytWtMat(const MatrixXd& Yt, int nMoments, bool useBanks){
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
    VectorXd means = aDiff.colwise().mean();
    VectorXd variances(nMoments);
    for(int i = 0; i < nMoments; i++){
        variances(i) = (aDiff.col(i).array() - aDiff.col(i).array().mean()).square().sum() / ((double) aDiff.col(i).array().size() - 1);
    }
    VectorXd covariances(nMoments - 1);
    
    for(int i = 0; i < nMoments - 1; i++){
        int j = i + 1;
        covariances(i) = ( (aDiff.col(i).array() - aDiff.col(i).array().mean()).array() * (aDiff.col(j).array() - aDiff.col(j).array().mean()).array() ).sum() / ((double) aDiff.col(i).array().size() - 1);
    }

    MatrixXd wt = MatrixXd::Zero(nMoments, nMoments);
   
    
    if(useBanks){
        for(int i = 0; i < nMoments; i++){
        wt(i,i) = variances(i); // cleanup code and make it more vectorized later.
        }
        for(int i = 0; i < nMoments - 1; i++){
            int j = i + 1;
            wt(i,j) = covariances(i);
            wt(j,i) = covariances(i);
        }
        cout << "Weights Before Inversion:" << endl << wt << endl;
        wt = wt.llt().solve(MatrixXd::Identity(nMoments, nMoments));
        cout << "Weights:" << endl;
        cout << wt << endl;
    }else{
        for(int i = 0; i < nMoments; i++){
            wt(i,i) = 1 / variances(i); // cleanup code and make it more vectorized later.
        }
        cout << "Weights:"<< endl;
        cout << wt << endl;
    }
    return wt;
}

MatrixXd customWtMat(const MatrixXd& Yt, const MatrixXd& Xt, int nMoments, int N, bool useBanks, bool useInverse){
    
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
    double cost = 0;
    VectorXd means = aDiff.colwise().mean();
    VectorXd variances(nMoments);
    for(int i = 0; i < nMoments; i++){
        variances(i) = (aDiff.col(i).array() - aDiff.col(i).array().mean()).square().sum() / ((double) aDiff.col(i).array().size() - 1);
    }
    
    MatrixXd wt = MatrixXd::Zero(nMoments, nMoments);
    if(useInverse){
        // compute covariances for differences.
        for(int i = 0; i < nMoments; i++){
            for(int j = i + 1; j < nMoments; j++){
                wt(i,j) = ((aDiff.col(i).array() - aDiff.col(i).array().mean()).array() * (aDiff.col(j).array() - aDiff.col(j).array().mean()).array() ).sum() / ((double) aDiff.col(i).array().size() - 1); 
                wt(j,i) = wt(i,j); // across diagonal
            }
        }

        wt = (wt / aDiff.rows()).inverse();
        cout << "wt:" << endl << wt << endl << endl;
    }else{
        if(useBanks){
            VectorXd covariances(nMoments - 1);
    
            for(int i = 0; i < nMoments - 1; i++){
                int j = i + 1;
                covariances(i) = ((aDiff.col(i).array() - aDiff.col(i).array().mean()).array() * (aDiff.col(j).array() - aDiff.col(j).array().mean()).array() ).sum() / ((double) aDiff.col(i).array().size() - 1);
            }
            for(int i = 0; i < nMoments; i++){
                wt(i,i) = variances(i); // cleanup code and make it more vectorized later.
            }
            for(int i = 0; i < nMoments - 1; i++){
                int j = i + 1;
                wt(i,j) = covariances(i);
                wt(j,i) = covariances(i);
            }
            cout << "Weights Before Inversion:" << endl << wt << endl;
            wt = wt.llt().solve(MatrixXd::Identity(nMoments, nMoments));
            cout << "Weights:" << endl;
            cout << wt << endl;
        }else{
            for(int i = 0; i < nMoments; i++){
                wt(i,i) = 1 / variances(i); // cleanup code and make it more vectorized later.
            }
            cout << "Weights:"<< endl;
            cout << wt << endl;
        }
    }
    
    return wt;
}