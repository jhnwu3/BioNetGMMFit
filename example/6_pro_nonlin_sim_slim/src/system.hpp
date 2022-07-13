#ifndef _SYSTEM_HPP_
#define _SYSTEM_HPP_
#include "main.hpp"
#include "linear.hpp"
#include "nonlinear.hpp"
MatrixXd interactionMatrix(int nSpecies, const VectorXd &k);

/* Nonlinear ODE System to be defined for nonlinear system. Simply refer to the comments for each term. */
class Nonlinear_ODE
{
    struct K rate;

public:
    Nonlinear_ODE(struct K G) : rate(G) {}

    void operator() (const State_N& c, State_N& dcdt, double t)
    {
       dcdt[0] = -(k(0) * c[0] * c[1])  // Syk = dc1/dt = k1 c1c2 + k2c3 + k3c3 
            + k(1) * c[2]
            + k(2) * c[2];

        dcdt[1] = -(k(0) * c[0] * c[1]) // Vav
            + k(1) * c[2]
            + k(5) * c[5];

        dcdt[2] = k(0) * c[0] * c[1] // Syk-Vav
            - k(1) * c[2]
            - k(2) * c[2];

        dcdt[3] = k(2) * c[2] //pVav
            - k(3) * c[3] * c[4]
            + k(4) * c[5];

        dcdt[4] = -(k(3) * c[3] * c[4]) // SHP1 
            + k(4) * c[5]
            + k(5) * c[5];

        dcdt[5] = k(3) * c[3] * c[4]  // SHP1-pVav
            - k(4) * c[5]
            - k(5) * c[5];
    }
};
#endif