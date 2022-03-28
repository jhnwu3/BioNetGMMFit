#ifndef _SYSTEM_HPP_
#define _SYSTEM_HPP_
#include "main.hpp"
#include "linear.hpp"
#include "nonlinear.hpp"
MatrixXd interactionMatrix(int nSpecies, const VectorXd &k);

/* Nonlinear ODE System to be defined for nonlinear system. Simply refer to the comments for each term. */
class Nonlinear_ODE
{
    // estimate vector
    VectorXd k;

public:
    Nonlinear_ODE(VectorXd G) : k(G) {}

    void operator() (const State_N& c, State_N& dcdt, double t)
    {
        /* WRITE YOUR SYSTEM OF DIFFERENTIAL EQUATIONS BELOW, BEWARE OF THE "()"" INDEXING FOR YOUR PARAMETERS AND "[]"" FOR YOUR PROTEINS */

        dcdt[0] = k(0) - k(4) * c[0]; // pcd3z
        dcdt[1] = k(1) * c[0] - k(4) * c[1]; // pslp
        dcdt[2] = k(2) * c[1] - k(4) * c[2]; // perk
        dcdt[3] = k(3) * c[2] - k(4) * c[3]; //ps6
    }
};
#endif