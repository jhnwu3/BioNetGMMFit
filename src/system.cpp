#include "system.hpp"
/*
Author: John Wu
Summary: Interaction Matrix used in the linear system.
 */
MatrixXd interactionMatrix(int nSpecies, const VectorXd &k){
    MatrixXd intMatrix = MatrixXd::Zero(nSpecies, nSpecies);
    /*--------- Define Interaction Matrix Here! (Below) ---------*/
    intMatrix << 
        -k(2), k(2), 0,
		k(1), -k(1) - k(4), k(4),
		k(3), k(0), -k(0) - k(3);

    /* 
        Make sure to define the matrix in the form like this with the semicolon ; at the end and make sure to have a comma between each matrix element.
       
       -k(2), k(2), 0,
		k(1), -k(1) - k(4), k(4),
		k(3), k(0), -k(0) - k(3);
        
        observe that between each comma (",") is a component in the system of linear differential equations. 

        ** IMPORTANT!!!!! **
        C++ indexing starts from 0, so the first rate constant is represented by k(0) not k(1)!


        <insert mathematical representation of an ode here, will ask Dr. Stewart about this, but I can almost guess>


    */
    
    /*--------- Define Interaction Matrix here! (Above) ---------*/
    return intMatrix;
}