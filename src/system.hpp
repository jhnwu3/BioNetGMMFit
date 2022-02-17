#ifndef _SYSTEM_HPP_
#define _SYSTEM_HPP_
#include "main.hpp"
#include "linear.hpp"
#include "nonlinear.hpp"
MatrixXd interactionMatrix(int nSpecies, const VectorXd &k);

/* Nonlinear ODE System to be defined for nonlinear system. Simply refer to the comments for each term. */
class Nonlinear_ODE
{
    // struct K rate;
    VectorXd k;

public:
    Nonlinear_ODE(VectorXd G) : k(G) {}

    void operator() (const State_N& c, State_N& dcdt, double t)
    {
        // dcdt[0] = -(rate.k(0) * c[0] * c[1])  // Syk = dc1/dt = k1 *c1*c2 + k2*c3 + k3*c3 
        //     + rate.k(1) * c[2]
        //     + rate.k(2) * c[2];

        // dcdt[1] = -(rate.k(0) * c[0] * c[1]) // Vav
        //     + rate.k(1) * c[2]
        //     + rate.k(5) * c[5];

        // dcdt[2] = rate.k(0) * c[0] * c[1] // Syk-Vav
        //     - rate.k(1) * c[2]
        //     - rate.k(2) * c[2];

        // dcdt[3] = rate.k(2) * c[2] //pVav
        //     - rate.k(3) * c[3] * c[4]
        //     + rate.k(4) * c[5];

        // dcdt[4] = -(rate.k(3) * c[3] * c[4]) // SHP1 
        //     + rate.k(4) * c[5]
        //     + rate.k(5) * c[5];

        // dcdt[5] = rate.k(3) * c[3] * c[4]  // SHP1-pVav
        //     - rate.k(4) * c[5]
        //     - rate.k(5) * c[5];

        // dcdt[0] = rate.k(0) - rate.k(5) * c[0];
        // dcdt[1] = rate.k(1) * c[0] - rate.k(4) * c[1];
        // dcdt[2] = rate.k(2) * c[1] - rate.k(4) * c[2];
        // dcdt[3] = rate.k(3) * c[2] - rate.k(4) * c[3];
        dcdt[0] = k(0) - k(5) * c[0];
        dcdt[1] = k(1) * c[0] - k(4) * c[1];
        dcdt[2] = k(2) * c[1] - k(4) * c[2];
        dcdt[3] = k(3) * c[2] - k(4) * c[3];
    }
};
#endif