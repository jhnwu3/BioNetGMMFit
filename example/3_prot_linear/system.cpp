#include "src/system.hpp"
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

    */
    
    /*--------- Define Interaction Matrix here! (Above) ---------*/
    return intMatrix;
}