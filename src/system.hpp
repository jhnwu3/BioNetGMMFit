#ifndef _SYSTEM_HPP_
#define _SYSTEM_HPP_
#include "main.hpp"
#include "linear.hpp"
#include "nonlinear.hpp"
MatrixXd interactionMatrix(int nSpecies, const VectorXd &k);
typedef std::vector<double> State_N;
typedef runge_kutta_cash_karp54< State_N > Error_RK_Stepper_N;
typedef controlled_runge_kutta< Error_RK_Stepper_N > Controlled_RK_Stepper_N;
/* Nonlinear ODE System to be defined for nonlinear system. Simply refer to the comments for each term. */
class Nonlinear_ODE
{
    // estimate vector
    VectorXd k;

public:
    Nonlinear_ODE(VectorXd G) : k(G) {}

    void operator() (const State_N& c, State_N& dcdt, double t)
    {
    /* Linear 3 */
        dcdt[0] = -k(2)*c[0] + k(1) * c[1] + k(3) * c[2]; // p1
        dcdt[1] = k(2) * c[0] + (-k(1) - k(4))*c[1] + k(0)*c[2] ; // p2
        dcdt[2] = k(4) * c[1] + (-k(0) - k(3))*c[2]; // p3
    }
};
#endif

/* WRITE YOUR SYSTEM OF DIFFERENTIAL EQUATIONS BELOW, BEWARE OF THE "()"" INDEXING FOR YOUR PARAMETERS AND "[]"" FOR YOUR PROTEINS */
        // dcdt[0] = -(k(0) * c[0] * c[1])  // Syk = dc1/dt = k1 c1c2 + k2c3 + k3c3 
        //     + k(1) * c[2]
        //     + k(2) * c[2];

        // dcdt[1] = -(k(0) * c[0] * c[1]) // Vav
        //     + k(1) * c[2]
        //     + k(5) * c[5];

        // dcdt[2] = k(0) * c[0] * c[1] // Syk-Vav
        //     - k(1) * c[2]
        //     - k(2) * c[2];

        // dcdt[3] = k(2) * c[2] //pVav
        //     - k(3) * c[3] * c[4]
        //     + k(4) * c[5];

        // dcdt[4] = -(k(3) * c[3] * c[4]) // SHP1 
        //     + k(4) * c[5]
        //     + k(5) * c[5];

        // dcdt[5] = k(3) * c[3] * c[4]  // SHP1-pVav
        //     - k(4) * c[5]
        //     - k(5) * c[5];

    /* Linear 4 */
    //     dcdt[0] = k(0) - k(4) * c[0]; // pcd3z
    //     dcdt[1] = k(1) * c[0] - k(4) * c[1]; // pslp
    //     dcdt[2] = k(2) * c[1] - k(4) * c[2]; // perk
    //     dcdt[3] = k(3) * c[2] - k(4) * c[3]; //ps6