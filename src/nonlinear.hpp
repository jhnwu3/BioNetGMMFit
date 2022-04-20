#ifndef _NONLINEAR_HPP_
#define _NONLINEAR_HPP_
#include "main.hpp"
#include "calc.hpp"
#include "fileIO.hpp"
#include "system.hpp"
/* MVN Generator Struct, currently unused in either model, but is useful for generating values from multivariate normal distributions */
struct Multi_Normal_Random_Variable
{
    Multi_Normal_Random_Variable(Eigen::MatrixXd const& covar)
        : Multi_Normal_Random_Variable(Eigen::VectorXd::Zero(covar.rows()), covar)
    {}

    Multi_Normal_Random_Variable(Eigen::VectorXd const& mean, Eigen::MatrixXd const& covar)
        : mean(mean)
    {
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenSolver(covar);
        transform = eigenSolver.eigenvectors() * eigenSolver.eigenvalues().cwiseSqrt().asDiagonal();
    }

    Eigen::VectorXd mean;
    Eigen::MatrixXd transform;

    Eigen::VectorXd operator()() const
    {
        static std::mt19937 gen{ std::random_device{}() }; //std::random_device {} ()
        static std::normal_distribution<> dist;

        return mean + transform * Eigen::VectorXd{ mean.size() }.unaryExpr([&](auto x) { return dist(gen); });
    }
};

/* Stepper Typedefs for ode solvers. */
typedef std::vector<double> State_N;
typedef runge_kutta_cash_karp54< State_N > Error_RK_Stepper_N;
typedef controlled_runge_kutta< Error_RK_Stepper_N > Controlled_RK_Stepper_N;

/* Nonlinear rate constant struct for use in Nonlinear class */
struct K
{
    VectorXd k;
};

/* Observer structs used to observe ode solver over time. Keep in mind, there should be no reason to change this for ode-solving  */
struct Protein_Components {
    int index;
    MatrixXd mat;
    VectorXd mVec;
    double timeToRecord;
    Protein_Components(double tf, int mom, int n, int nSpecies) {
        mVec = VectorXd::Zero(mom);
        mat = MatrixXd::Zero(n, nSpecies);
        timeToRecord = tf;
    }
};

struct Moments_Mat_Obs
{
    struct Protein_Components& dComp;
    Moments_Mat_Obs(struct Protein_Components& dCom) : dComp(dCom) {}
    void operator()(State_N const& c, const double t) const
    {
        if (t == dComp.timeToRecord) {
            int upperDiag = 2 * dComp.mat.cols();
            for (int i = 0; i < dComp.mat.cols(); i++) {
                dComp.mVec(i) += c[i];
                dComp.mat(dComp.index, i) = c[i];
                for (int j = i; j < dComp.mat.cols(); j++) {
                    if (i == j && (dComp.mat.cols() + i) < dComp.mVec.size()) { // diagonal elements
                        dComp.mVec(dComp.mat.cols() + i) += c[i] * c[j]; // 2nd moments
                    }
                    else if (upperDiag < dComp.mVec.size()){
                        dComp.mVec(upperDiag) += c[i] * c[j]; // cross moments
                        upperDiag++;
                    }
                }
            }
        }
    }
};

State_N convertInit(const VectorXd &v1);
VectorXd adaptVelocity(const VectorXd& posK, int seed, double epsi, double nan, int hone);
MatrixXd nonlinearModel(int nParts, int nSteps, int nParts2, int nSteps2, const MatrixXd& X_0, const MatrixXd &Y_0, int nRates, int nRuns, int nMoments);


#endif