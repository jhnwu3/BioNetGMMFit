// PSO_S.cpp : Dr. Stewart's 3 variable linear system R script in C++ format.
//

#include <iostream>
#include <fstream>
#include <boost/math/distributions.hpp>
#include <random>
#include <vector>
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>
#include <cmath>
#include <chrono>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;
using namespace boost;
using namespace boost::math;

int main() {
	
	auto t1 = std::chrono::high_resolution_clock::now();
	/*---------------------- Setup ------------------------ */
	int bsi = 1, Nterms = 9, useEqual = 0, Niter = 1, Biter = 1; 

	/* Variables (global) */
	
	int wasflipped = 0, Nprots = 3, Npars = 5;
	double squeeze = 0.96, sdbeta = 0.05;

	/* SETUP */
	int useDiag = 0;
	int sf1 = 1;
	int sf2 = 1;
	
	int N = 10000;

	int Nparts_1 = 5000;
	int Nsteps_1 = 5;

	int Nparts_2 = 5;
	int Nsteps_2 = 5000;

	// note for coder: wmatup is a list 
	vector<double> wmatup; 
	wmatup.push_back(0.15);
	wmatup.push_back(0.30);
	wmatup.push_back(0.45);
	wmatup.push_back(0.60);

	double dp = 1, sfp = 3, sfg = 1, sfe = 6;

	vector<double> k;
	k.push_back(0.27678200 / sf1);
	k.push_back(0.83708059 / sf1);
	k.push_back(0.44321700 / sf1);
	k.push_back(0.04244124 / sf1);
	k.push_back(0.30464502 / sf1);

	vector<double> truk; // make a copy of a vector/ array/list 

	for (unsigned int i = 0; i < k.size(); i++) {
		truk.push_back(k.at(i));
	}

	// print truk values to a .par file w/ 5 columns? 
	ofstream truk_file("truk.par");
	for (unsigned int i = 0; i < truk.size(); i++) {
		truk_file << " " << truk.at(i);
	}
	truk_file.close();

	
	double mu_x = 1.47, mu_y = 1.74, mu_z = 1.99; // true means for MVN(theta)

	double var_x = 0.77, var_y = 0.99, var_z = 1.11; // true variances for MVN(theta);

	double rho_xy = 0.10, rho_xz = 0.05, rho_yz = 0.10; // true correlations for MVN

	double sigma_x = sqrt(var_x), sigma_y = sqrt(var_y), sigma_z = sqrt(var_z);

	double cov_xy = rho_xy * sigma_x * sigma_y;
	double cov_xz = rho_xz * sigma_x * sigma_z;
	double cov_yz = rho_yz * sigma_y * sigma_z;

	/* sigma matrices */
	MatrixXd sigma_12(1, 2);
	sigma_12 << cov_xz, cov_yz;
	
	MatrixXd sigma_22(2, 2);
	sigma_22 << var_x, cov_xy,
		cov_xy, var_y;
	
	MatrixXd sigma_21(2, 1);
	sigma_21 = sigma_12.transpose();

	/* conditional variances of proteins*/
	double cvar_ygx = (1 - (rho_xy * rho_xy)) * (sigma_y * sigma_y);

	double cvar_zgxy = var_z - (sigma_12 * sigma_22.inverse() * sigma_21)(0,0); // note: since matrix types are incompatible witgh doubles, we must do matrix math first, then convert to double.

	int t = 3; // num steps away from initial state

	// first instantiate a matrix, then performance matrix math to form MT
	MatrixXd M(3,3);
	M << -k.at(2), k.at(2), 0,
		k.at(1), -k.at(1) - k.at(4), k.at(4),
		k.at(3), k.at(0), -k.at(0) - k.at(3);

	MatrixXd MT(3, 3);
	MT = t * M.transpose();
   
	MatrixXd EMT(3, 3);
	EMT = MT.exp();
	
	/* actually global variables that are being recalculated in PSO */
	double omp_1, omp_2, omp_3, ovp_1 = 0, ovp_2 = 0, ovp_3 = 0, ocov_12, ocov_13, ocov_23;
	double pmp_1, pmp_2, pmp_3, pvp_1 = 0, pvp_2 = 0, pvp_3 = 0, pcov_12, pcov_13, pcov_23;
	double cost_seedk, cost_gbest, cost_sofar;
	MatrixXd X_0(N, 3);
	MatrixXd X_0_obs(N, 3);
	MatrixXd Y_t_obs(N, 3);
	MatrixXd Y_t(N, 3);
	VectorXd pmpV(3);

	MatrixXd GBMAT;
	MatrixXd w_mat(9, 9);

	VectorXd gbest(Npars), best_sofar(Npars);

	VectorXd x(N); //note the data sample x is a list of 10000 RV from normal dist
	VectorXd pa_x(N);
	VectorXd y(N);
	VectorXd pa_y(N);
	VectorXd z(N); // z big questions about how to get the data values for it. It takes in a far bigger value???
	VectorXd pa_z(N);

	VectorXd all_terms(9);
	VectorXd term_vec(9);

	/* IMPORTANT THAT YOU INSTANTIATE THE RANDOM GENERATOR LIKE THIS!*/
	std::random_device rand_dev;
	std::mt19937 generator(rand_dev());
	uniform_real_distribution<double> unifDist(0.0, 1.0);

	for (int q = 1; q <= Niter; q++) {

		int dpFlag = 1;   // Unsure what this print statement does, will ask later.
		if (q % 10 == 0) {
			cout << "Working on replicate " << q << "\n";
		}

		if (bsi == 0 || q == 1) {
			/* Simulate Y(t) and X(0) */
			
			
			std::normal_distribution<double> xNorm(mu_x, sigma_x);

			for (int i = 0; i < N; i++) {
				x(i) = (xNorm(generator));
				pa_x(i) = (exp(x(i)));
			}
			
			
			for (int i = 0; i < x.size(); i++) {
				std::normal_distribution<double> yNorm(mu_y + sigma_y * rho_xy * (x(i) - mu_x) / sigma_x, sqrt(cvar_ygx));
				y(i) = (yNorm(generator));
				pa_y(i) = (exp(y(i)));
			}
			
			/* matrix math for the z random vals. */
			MatrixXd rbind(2, N); // first calculate a 2xN rbind matrix
			for (int i = 0; i < x.size(); i++) {
				rbind(0, i) = x(i) - mu_x;
				rbind(1, i) = y(i) - mu_y;
			}
			MatrixXd zMean(1, N); // calculate the vector of means
			zMean = sigma_12 * sigma_22.inverse() * rbind;
			for (int i = 0; i < zMean.size(); i++) {
				zMean(0, i) = zMean(0, i) + mu_z;
			}
			// finally actually calculate z and pa_z vectors
			for (int i = 0; i < N; i++) {
				std::normal_distribution<double> zNorm(zMean(0,i), sqrt(cvar_zgxy));
				z(i) = (zNorm(generator));
				pa_z(i) = (exp(z(i)));
			}

			/* Create Y.0 */
			MatrixXd Y_0(N, 3);
			/*for (int i = 0; i < N; i++) {
				// fill it up from vectors
				Y_0(i, 0) = pa_x(i);
				Y_0(i, 1) = pa_y(i);
				Y_0(i, 2) = pa_z(i);
			}*/
			Y_0.col(0) = pa_x;
			Y_0.col(1) = pa_y;
			Y_0.col(2) = pa_z;
			
			Y_t = (EMT * Y_0.transpose()).transpose();
			
			if (bsi == 1 && q == 1) {
				Y_t_obs = Y_t;
			}

			/*  # Compute the observed means, variances, and covariances
				# Add random noise to Y.t
				# trusd < -apply(Y.t, 2, sd)
				# error     <- t(matrix(rnorm(N*Nprots,rep(0,Nprots),trusd*0.01),nrow=3))
				# Y.t < -Y.t + error */

			/* means */
			
			VectorXd ompV = Y_t.colwise().mean();
		
			omp_1 = ompV(0);
			omp_2 = ompV(1);
			omp_3 = ompV(2);

			/* variances - actually have to manually calculate it, no easy library  */
			ovp_1 = (Y_t.col(0).array() - Y_t.col(0).array().mean()).square().sum() / ((double)Y_t.col(0).array().size() - 1);
			ovp_2 = (Y_t.col(1).array() - Y_t.col(1).array().mean()).square().sum() / ((double)Y_t.col(1).array().size() - 1);
			ovp_3 = (Y_t.col(2).array() - Y_t.col(2).array().mean()).square().sum() / ((double)Y_t.col(2).array().size() - 1);

			/* covariances - also requires manual calculation*/
			double sum12 = 0, sum13 = 0, sum23 = 0;
			for (int n = 0; n < N; n++)
			{
				sum12 += (Y_t(n, 0) - omp_1) * (Y_t(n, 1) - omp_2);
				sum13 += (Y_t(n, 0) - omp_1) * (Y_t(n, 2) - omp_3);
				sum23 += (Y_t(n, 1) - omp_2) * (Y_t(n, 2) - omp_3);

			}
			int N_SUBTRACT_ONE = N - 1;
			ocov_12 = sum12 / N_SUBTRACT_ONE;
			ocov_13 = sum13 / N_SUBTRACT_ONE;
			ocov_23 = sum23 / N_SUBTRACT_ONE;
			
			for (int i = 0; i < N; i++) {
				x(i) = (xNorm(generator));
				pa_x(i) = (exp(x(i)));
			}

		    
			for (int i = 0; i < N; i++) {
				std::normal_distribution<double> yNorm(mu_y + sigma_y * rho_xy * (x(i) - mu_x) / sigma_x, sqrt(cvar_ygx));
				y(i) = (yNorm(generator));
				pa_y(i) = (exp(y(i)));
			}

			/* matrix math for the z random vals. */
			MatrixXd r1bind(2, N); // first calculate a 2xN rbind matrix
			for (int i = 0; i < N; i++) {
				r1bind(0, i) = x(i) - mu_x;
				r1bind(1, i) = y(i) - mu_y;
			}
			MatrixXd z1Mean(1, N); // calculate the vector of means
			z1Mean = sigma_12 * sigma_22.inverse() * r1bind;
			for (int i = 0; i < z1Mean.size(); i++) {
				z1Mean(0, i) = z1Mean(0, i) + mu_z;
			}
			// finally actually calculate z and pa_z vectors
			for (int i = 0; i < N; i++) {
				std::normal_distribution<double> zNorm(z1Mean(0, i), sqrt(cvar_zgxy));
				z(i) = (zNorm(generator));
				pa_z(i) = (exp(z(i)));
			}

			
			/*for (int i = 0; i < N; i++) {
				// fill it up from vectors
				X_0(i, 0) = pa_x(i);
				X_0(i, 1) = pa_y(i);
				X_0(i, 2) = pa_z(i);
			}*/
			X_0.col(0) = pa_x;
			X_0.col(1) = pa_y;
			X_0.col(2) = pa_z;
			
			if (bsi == 1 && q == 1) {// save the simulated CYTOF data time 0
				X_0_obs = X_0;
			}

		}
		
		if (bsi == 1 && q > 1) {

			/* create shuffled indices based on uniform rand dist */
			vector<int> bindices;
			for (int i = 0; i < N; i++) {
				bindices.push_back(i);
			}
			shuffle(bindices.begin(), bindices.end(), generator); // shuffle indices as well as possible. 

			MatrixXd X_0(N, 3);
			MatrixXd Y_t(N, 3);
			/* shuffle all of the values in the observed matrices*/
			for (int i = 0; i < N; i++) {
				X_0(i , 0) = X_0_obs(bindices.at(i), 0);
				X_0(i, 1) = X_0_obs(bindices.at(i), 1);
				X_0(i, 2) = X_0_obs(bindices.at(i), 2);
				Y_t(i, 0) = Y_t_obs(bindices.at(i), 0);
				Y_t(i, 1) = Y_t_obs(bindices.at(i), 1);
				Y_t(i, 2) = Y_t_obs(bindices.at(i), 2);
			}
			
			/* re-calc new omp, ovp, and ocovs, which should be the same???*/
			VectorXd ompV = Y_t.colwise().mean();
			
			omp_1 = ompV(0);
			omp_2 = ompV(1);
			omp_3 = ompV(2);

			
			/* variances - actually have to manually calculate it, no easy library  */
		
			ovp_1 = (Y_t.col(0).array() - Y_t.col(0).array().mean()).square().sum() / ((double)Y_t.col(0).array().size() - 1);
			ovp_2 = (Y_t.col(1).array() - Y_t.col(1).array().mean()).square().sum() / ((double)Y_t.col(1).array().size() - 1);
			ovp_3 = (Y_t.col(2).array() - Y_t.col(2).array().mean()).square().sum() / ((double)Y_t.col(2).array().size() - 1);

	

			/* covariances - also requires manual calculation*/
			double sum12 = 0, sum13 = 0, sum23 = 0;
			for (int n = 0; n < N; n++)
			{
				sum12 += (Y_t(n, 0) - omp_1) * (Y_t(n, 1) - omp_2);
				sum13 += (Y_t(n, 0) - omp_1) * (Y_t(n, 2) - omp_3);
				sum23 += (Y_t(n, 1) - omp_2) * (Y_t(n, 2) - omp_3);

			}
			int N_SUBTRACT_ONE = N - 1;
			ocov_12 = sum12 / N_SUBTRACT_ONE;
			ocov_13 = sum13 / N_SUBTRACT_ONE;
			ocov_23 = sum23 / N_SUBTRACT_ONE;

		}
		
		// Initialize variables to start the layered particle swarms
		int Nparts = Nparts_1;
		int Nsteps = Nsteps_1;

		VectorXd vectorOfOnes(9);		
		for (int i = 0; i < 9; i++) {
			vectorOfOnes(i) = 1;
		}
		w_mat = vectorOfOnes.asDiagonal(); //initialize weight matrix
		
		VectorXd seedk(Npars); //initialize global best
		for (int i = 0; i < Npars; i++) { seedk(i) = unifDist(generator) /sf2; }

		/*Compute cost of seedk */
		for (int i = 0; i < Npars; i++) { k.at(i) = seedk(i); }

		MatrixXd HM(3, 3);
		HM << -k.at(2), k.at(2), 0,
			k.at(1), -k.at(1) - k.at(4), k.at(4),
			k.at(3), k.at(0), -k.at(0) - k.at(3);

		MatrixXd HMT(3, 3);
		HMT = t * HM.transpose();
		
		MatrixXd EHMT;
		EHMT = HMT.exp();
		
		MatrixXd Q(N,3);
		Q = (EHMT * X_0.transpose()).transpose();
		
		//re-calc new omp, ovp, and ocovs, which should be the same???
	    pmpV = Q.colwise().mean();

		pmp_1 = pmpV(0);
		pmp_2 = pmpV(1);
		pmp_3 = pmpV(2);
	
		// variances - actually have to manually calculate it, no easy library
		pvp_1 = (Q.col(0).array() - Q.col(0).array().mean()).square().sum() / ((double) Q.col(0).array().size() - 1);
		pvp_2 = (Q.col(1).array() - Q.col(1).array().mean()).square().sum() / ((double) Q.col(1).array().size() - 1);
		pvp_3 = (Q.col(2).array() - Q.col(2).array().mean()).square().sum() / ((double) Q.col(2).array().size() - 1);

		// covariances - also requires manual calculation 
		double sum12 = 0, sum13 = 0, sum23 = 0;
		
		for (int n = 0; n < Q.rows(); n++)
		{
			sum12 += (Q(n, 0) - pmp_1) * (Q(n, 1) - pmp_2);	
			sum13 += (Q(n, 0) - pmp_1) * (Q(n, 2) - pmp_3);
			sum23 += (Q(n, 1) - pmp_2) * (Q(n, 2) - pmp_3);

		}
		double N_SUBTRACT_ONE = Q.rows() - 1.0;
		
		pcov_12 = sum12 / N_SUBTRACT_ONE;
		pcov_13 = sum13 / N_SUBTRACT_ONE;
		pcov_23 = sum23 / N_SUBTRACT_ONE;

		double term_1 = pmp_1 - omp_1, 
			term_2 = pmp_2 - omp_2, 
			term_3 = pmp_3 - omp_3, 
			term_4 = pvp_1 - ovp_1, 
			term_5 = pvp_2 - ovp_2,
			term_6 = pvp_3 - ovp_3,
			term_7 = pcov_12 - ocov_12, 
			term_8 = pcov_13 - ocov_13,
			term_9 = pcov_23 - ocov_23;
		// note to self: I'm using vectorXd from now on b/c it plays way better than the vectors built into C++ unless ofc there are strings that we need to input.
	
		all_terms << term_1, term_2, term_3, term_4, term_5, term_6, term_7, term_8, term_9;
	
		
		term_vec = all_terms;
		
		
		cost_seedk = term_vec.transpose() * w_mat * (term_vec.transpose()).transpose();

		// instantiate values 
		gbest = seedk;
		best_sofar = seedk;
		cost_gbest = cost_seedk;
		cost_sofar = cost_seedk;

		GBMAT.conservativeResize(1, 6);
		
		// will probably find a better method later, but will for now just create temp vec to assign values.
		
		VectorXd cbind(gbest.size() + 1);
		cbind << gbest, cost_gbest;
		GBMAT.row(GBMAT.rows() - 1) = cbind;
	
		
		double nearby = sdbeta;
		MatrixXd POSMAT(Nparts, Npars);
		
		for (int pso = 1; pso <= Biter + 1 ; pso++) {
			cout << "PSO:" << pso << endl;

			if (pso < Biter + 1) {
				for (int i = 0; i < Nparts; i++) {
					// row by row in matrix using uniform dist.
					for (int n = 0; n < Npars; n++) {
						POSMAT(i, n) = unifDist(generator) / sf2;
					}
				}
				
			}
			
			if (pso == Biter + 1) {
				Nparts = Nparts_2;
				Nsteps = Nsteps_2;

				GBMAT.conservativeResize(GBMAT.rows() + 1, 6);
				cbind << best_sofar, cost_sofar;
				GBMAT.row(GBMAT.rows() - 1) = cbind;

				gbest = best_sofar;
				cost_gbest = cost_sofar;

				// reset POSMAT? 
				POSMAT.resize(Nparts, Npars);
				POSMAT.setZero();
	
				for (int init = 0; init < Nparts; init++) {
					for (int edim = 0; edim < Npars; edim++) {
						double tmean = gbest(edim);
						if (gbest(edim) > 0.5) {
							tmean = 1 - gbest(edim);
							wasflipped = 1;
						}
						double myc = (1 - tmean) / tmean;
						double alpha = myc / ((1 + myc) * (1 + myc) * (1 + myc)*nearby*nearby);
						double beta = myc * alpha;

						std::gamma_distribution<double> aDist(alpha, 1);
						std::gamma_distribution<double> bDist(beta, 1);

						double x = aDist(generator);
						double y = bDist(generator);
						double myg = x / (x + y);
						// sample from beta dist - this can be quite inefficient and taxing, there is another way with a gamma dist (THAT NEEDS TO BE REINVESTIGATED), but works so far. 
						//beta_distribution<double> betaDist(alpha, beta);
						//double randFromUnif = unifDist(generator);
						//double myg = quantile(betaDist, randFromUnif);

						if (wasflipped == 1) {
							wasflipped = 0;
							myg = 1 - myg;
						}
						POSMAT(init, edim) = myg;
					}
				}
				
			} 
			
			// initialize PBMAT 
			MatrixXd PBMAT = POSMAT; // keep track of ea.particle's best, and it's corresponding cost
			
			PBMAT.conservativeResize(POSMAT.rows(), POSMAT.cols() + 1);
			for (int i = 0; i < PBMAT.rows(); i++) { PBMAT(i, PBMAT.cols() - 1) = 0; } // add the 0's on far right column
			
			for (int h = 0; h < Nparts; h++) {
				for (int init = 0; init < Npars; init++) { k.at(init) = PBMAT(h, init); }

				HM << -k.at(2), k.at(2), 0,
					k.at(1), -k.at(1) - k.at(4), k.at(4),
					k.at(3), k.at(0), -k.at(0) - k.at(3);

				HMT = t * HM.transpose();
				EHMT = HMT.exp(); // NOTE MATRIX EXPONENTIATION ACTUALLY TAKES A REALLY LONG TIME TO RUN
				Q = (EHMT * X_0.transpose()).transpose();

				pmpV = Q.colwise().mean();
				
				pmp_1 = pmpV(0);
				pmp_2 = pmpV(1);
				pmp_3 = pmpV(2);

				// variances - below is manual calculation  
				pvp_1 = (Q.col(0).array() - Q.col(0).array().mean()).square().sum() / ((double)Q.col(0).array().size() - 1);
				pvp_2 = (Q.col(1).array() - Q.col(1).array().mean()).square().sum() / ((double)Q.col(1).array().size() - 1);
				pvp_3 = (Q.col(2).array() - Q.col(2).array().mean()).square().sum() / ((double)Q.col(2).array().size() - 1);
				
				// covariances - also requires manual calculation 
				double sum12 = 0, sum13 = 0, sum23 = 0;

				for (int n = 0; n < Q.rows(); n++)
				{
					sum12 += (Q(n, 0) - pmp_1) * (Q(n, 1) - pmp_2);
					sum13 += (Q(n, 0) - pmp_1) * (Q(n, 2) - pmp_3);
					sum23 += (Q(n, 1) - pmp_2) * (Q(n, 2) - pmp_3);
				}
			    N_SUBTRACT_ONE = Q.rows() - 1.0;

				pcov_12 = sum12 / N_SUBTRACT_ONE;
				pcov_13 = sum13 / N_SUBTRACT_ONE;
				pcov_23 = sum23 / N_SUBTRACT_ONE;
				
				 term_1 = pmp_1 - omp_1,
					term_2 = pmp_2 - omp_2,
					term_3 = pmp_3 - omp_3,
					term_4 = pvp_1 - ovp_1,
					term_5 = pvp_2 - ovp_2,
					term_6 = pvp_3 - ovp_3,
					term_7 = pcov_12 - ocov_12,
					term_8 = pcov_13 - ocov_13,
					term_9 = pcov_23 - ocov_23;
				// note to self: I'm using vectorXd from now on b/c it plays way better than the vectors built into C++ unless ofc there are strings that we need to input.
				
				all_terms << term_1, term_2, term_3, term_4, term_5, term_6, term_7, term_8, term_9;
				term_vec = all_terms; 

				PBMAT(h, Npars) = term_vec.transpose() * w_mat * (term_vec.transpose()).transpose();
			}

			
			// ALL SWARMS BEGIN TO MOVE HERE 
			double sfi = sfe;
			double sfc = sfp;
			double sfs = sfg;

			for (int iii = 0; iii < Nsteps; iii++) { //REMEMBER IF THERE IS ITERATION WITH iii MAKE SURE TO SUBTRACT ONE

				
				if (pso == (Biter + 1)) {
					vector<int> chkpts;
					
					for (unsigned int i = 0; i < wmatup.size(); i++) {	
						chkpts.push_back(wmatup.at(i)* Nsteps);
					}
					
					if (iii == chkpts.at(0) || iii == chkpts.at(1) || iii == chkpts.at(2) || iii == chkpts.at(3)) {
						nearby = squeeze * nearby;

						for (int i = 0; i < Npars; i++) { // best estimate of k to compute w.mat
							k.at(i) = gbest(i);
						}

						HM << -k.at(2), k.at(2), 0,
							k.at(1), -k.at(1) - k.at(4), k.at(4),
							k.at(3), k.at(0), -k.at(0) - k.at(3);

						HMT = t * HM.transpose();
						EHMT = HMT.exp();
						Q = (EHMT * X_0.transpose()).transpose();

						MatrixXd fmdiffs(Q.rows(), 3);
						fmdiffs = Y_t - Q;
						
						VectorXd mxt(3);
						mxt = Q.colwise().mean();
				
						VectorXd myt(3); 
						myt = Y_t.colwise().mean();
					
						MatrixXd residxt(Q.rows(), Q.cols());
						residxt.col(0) = mxt.row(0).replicate(N, 1);
						residxt.col(1) = mxt.row(1).replicate(N, 1);
						residxt.col(2) = mxt.row(2).replicate(N, 1);
						residxt = Q - residxt;

						MatrixXd residyt(Y_t.rows(), Y_t.cols());
						residyt.col(0) = myt.row(0).replicate(N, 1);
						residyt.col(1) = myt.row(1).replicate(N, 1);
						residyt.col(2) = myt.row(2).replicate(N, 1);
						residyt = Y_t - residyt;

						MatrixXd smdiffs(N, 3);
						smdiffs = (residyt.array() * residyt.array()) - (residxt.array()* residxt.array());

						MatrixXd cprxt(N, 3);
						cprxt.col(0) = residxt.col(0).array() * residxt.col(1).array();
						cprxt.col(1) = residxt.col(0).array() * residxt.col(2).array();
						cprxt.col(2) = residxt.col(1).array() * residxt.col(2).array();

						MatrixXd cpryt(N, 3);
						cpryt.col(0) = residyt.col(0).array() * residyt.col(1).array();
						cpryt.col(1) = residyt.col(0).array() * residyt.col(2).array();
						cpryt.col(2) = residyt.col(1).array() * residyt.col(2).array();

						MatrixXd cpdiffs(N, 3);
						cpdiffs = cpryt - cprxt;

						MatrixXd Adiffs(N, 9);
						Adiffs << fmdiffs, smdiffs, cpdiffs; // concatenate

						MatrixXd g_mat(N, Nterms);
						g_mat = Adiffs;
						for (int m = 0; m < N; m++) { w_mat = w_mat + g_mat.row(m).transpose() * g_mat.row(m); }
						w_mat = w_mat / N;
						w_mat = w_mat.inverse();

						if (useDiag == 1) { w_mat = w_mat.diagonal().diagonal(); }

						// CALCULATE MEANS, VARIANCES, AND COVARIANCES
						pmpV = Q.colwise().mean();

						pmp_1 = pmpV(0);
						pmp_2 = pmpV(1);
						pmp_3 = pmpV(2);

						// variances - manually calculate it, no easy library 
						pvp_1 = (Q.col(0).array() - Q.col(0).array().mean()).square().sum() / ((double)Q.col(0).array().size() - 1);
						pvp_2 = (Q.col(1).array() - Q.col(1).array().mean()).square().sum() / ((double)Q.col(1).array().size() - 1);
						pvp_3 = (Q.col(2).array() - Q.col(2).array().mean()).square().sum() / ((double)Q.col(2).array().size() - 1);

						// covariances - manual calculation 
						double sum12 = 0, sum13 = 0, sum23 = 0;

						for (int n = 0; n < Q.rows(); n++)
						{
							sum12 += (Q(n, 0) - pmp_1) * (Q(n, 1) - pmp_2);
							sum13 += (Q(n, 0) - pmp_1) * (Q(n, 2) - pmp_3);
							sum23 += (Q(n, 1) - pmp_2) * (Q(n, 2) - pmp_3);

						}
						N_SUBTRACT_ONE = Q.rows() - 1.0;

						pcov_12 = sum12 / N_SUBTRACT_ONE;
						pcov_13 = sum13 / N_SUBTRACT_ONE;
						pcov_23 = sum23 / N_SUBTRACT_ONE;

						term_1 = pmp_1 - omp_1,
							term_2 = pmp_2 - omp_2,
							term_3 = pmp_3 - omp_3,
							term_4 = pvp_1 - ovp_1,
							term_5 = pvp_2 - ovp_2,
							term_6 = pvp_3 - ovp_3,
							term_7 = pcov_12 - ocov_12,
							term_8 = pcov_13 - ocov_13,
							term_9 = pcov_23 - ocov_23;
						// note to self: I'm using vectorXd from now on b/c it plays way better than the vectors built into C++ unless ofc there are strings that we need to input.
						
						all_terms << term_1, term_2, term_3, term_4, term_5, term_6, term_7, term_8, term_9;
						term_vec = all_terms;
						
						cost_gbest = term_vec.transpose() * w_mat * term_vec.transpose().transpose();
						
						GBMAT.conservativeResize(GBMAT.rows() + 1, GBMAT.cols());
						VectorXd cbind1(GBMAT.cols());
						cbind1 << gbest, cost_gbest;
						GBMAT.row(GBMAT.rows() - 1) = cbind1;
						
						POSMAT.resize(Nparts,Npars); //reset to 0???
						POSMAT.setZero();
						
						for (int init = 0; init < Nparts; init++) {
							for (int edim = 0; edim < Npars; edim++) {
								double tmean = gbest(edim);
								if (gbest(edim) > 0.5) {
									tmean = 1 - gbest(edim);
									wasflipped = 1;
									
								}
								double myc = (1 - tmean) / tmean;
								double alpha = myc / ((1 + myc) * (1 + myc) * (1 + myc) * nearby * nearby);
								double beta = myc * alpha;

								// sample from beta dist - this can be quite inefficient and taxing, there is another way with a gamma dist (THAT NEEDS TO BE REINVESTIGATED), but works so far.
								std::gamma_distribution<double> aDist(alpha, 1);
								std::gamma_distribution<double> bDist(beta, 1);

								double x = aDist(generator);
								double y = bDist(generator);
								double myg = x/(x+y);
								if (wasflipped == 1) {
									wasflipped = 0;
									myg = 1 - myg;
								}
								POSMAT(init, edim) = myg;
							}
						}
						
						MatrixXd cbindMat(POSMAT.rows(), POSMAT.cols() + 1); // keep track of each particle's best and it's corresponding cost
						cbindMat << POSMAT, VectorXd::Zero(POSMAT.rows());

						for (int h = 0; h < Nparts; h++) {
							for (int init = 0; init < Npars; init++) { k.at(init) = PBMAT(h, init); } 

							HM << -k.at(2), k.at(2), 0,
								k.at(1), -k.at(1) - k.at(4), k.at(4),
								k.at(3), k.at(0), -k.at(0) - k.at(3);

							HMT = t * HM.transpose();
							EHMT = HMT.exp();
							Q = (EHMT * X_0.transpose()).transpose();


							// CALCULATE MEANS, VARIANCES, AND COVARIANCES
							pmpV = Q.colwise().mean();

							pmp_1 = pmpV(0);
							pmp_2 = pmpV(1);
							pmp_3 = pmpV(2);

							// variances - manually calculate it, no easy library 
							pvp_1 = (Q.col(0).array() - Q.col(0).array().mean()).square().sum() / ((double)Q.col(0).array().size() - 1);
							pvp_2 = (Q.col(1).array() - Q.col(1).array().mean()).square().sum() / ((double)Q.col(1).array().size() - 1);
							pvp_3 = (Q.col(2).array() - Q.col(2).array().mean()).square().sum() / ((double)Q.col(2).array().size() - 1);
							// covariances - manual calculation 
							double sum12 = 0, sum13 = 0, sum23 = 0;

							for (int n = 0; n < Q.rows(); n++)
							{
								sum12 += (Q(n, 0) - pmp_1) * (Q(n, 1) - pmp_2);
								sum13 += (Q(n, 0) - pmp_1) * (Q(n, 2) - pmp_3);
								sum23 += (Q(n, 1) - pmp_2) * (Q(n, 2) - pmp_3);

							}
							N_SUBTRACT_ONE = Q.rows() - 1.0;

							pcov_12 = sum12 / N_SUBTRACT_ONE;
							pcov_13 = sum13 / N_SUBTRACT_ONE;
							pcov_23 = sum23 / N_SUBTRACT_ONE;

							term_1 = pmp_1 - omp_1,
								term_2 = pmp_2 - omp_2,
								term_3 = pmp_3 - omp_3,
								term_4 = pvp_1 - ovp_1,
								term_5 = pvp_2 - ovp_2,
								term_6 = pvp_3 - ovp_3,
								term_7 = pcov_12 - ocov_12,
								term_8 = pcov_13 - ocov_13,
								term_9 = pcov_23 - ocov_23;
							// note to self: I'm using vectorXd from now on b/c it plays way better than the vectors built into C++ unless ofc there are strings that we need to input.

							all_terms << term_1, term_2, term_3, term_4, term_5, term_6, term_7, term_8, term_9;
							term_vec = all_terms;
							PBMAT(h, Npars) = term_vec.transpose() * w_mat * term_vec.transpose().transpose();
						}
					}
				}

				//cout << "line 802" << endl;
				for (int jjj = 0; jjj < Nparts; jjj++) {

					double w1 = sfi * unifDist(generator) /sf2, w2 = sfc*  unifDist(generator) / sf2, w3 = sfs * unifDist(generator) / sf2;
					double sumw = w1 + w2 + w3;

					w1 = w1 / sumw;
					w2 = w2 / sumw;
					w3 = w3 / sumw;

					// R -sample ~ shuffle
					
					vector<int> seqOneToFive;
					seqOneToFive.clear();
					for (int i = 0; i < Npars; i++) {
						seqOneToFive.push_back(i);
					}
					shuffle(seqOneToFive.begin(), seqOneToFive.end(), generator); // shuffle indices as well as possible. 
					int ncomp = seqOneToFive.at(0);
					VectorXd wcomp(ncomp);
					shuffle(seqOneToFive.begin(), seqOneToFive.end(), generator);
					for (int i = 0; i < ncomp; i++) {
						wcomp(i) = seqOneToFive.at(i);
					}
				
					VectorXd rpoint = POSMAT.row(jjj);

					for (int smart = 0; smart < ncomp; smart++) {
						int px = wcomp(smart);
						double pos = rpoint(px);
						double alpha = 4 * pos;
						double beta = 4 - alpha;

						std::gamma_distribution<double> aDist(alpha, 1);
						std::gamma_distribution<double> bDist(beta, 1);
						
						double x = aDist(generator);
						double y = bDist(generator);
						
						rpoint(px) = (x/(x+y)) / sf2;
					}
					
					VectorXd PBMATV(5);
					PBMATV << PBMAT(jjj, 0), PBMAT(jjj, 1), PBMAT(jjj, 2), PBMAT(jjj, 3), PBMAT(jjj, 4);
					POSMAT.row(jjj) = w1 * rpoint + w2 * PBMATV + w3 * gbest;


					
					for (int i = 0; i < Npars; i++) { k.at(i) = POSMAT(jjj, i); }

					HM << -k.at(2), k.at(2), 0,
						k.at(1), -k.at(1) - k.at(4), k.at(4),
						k.at(3), k.at(0), -k.at(0) - k.at(3);

					HMT = t * HM.transpose();
					EHMT = HMT.exp();
					Q = (EHMT * X_0.transpose()).transpose();


					// CALCULATE MEANS, VARIANCES, AND COVARIANCES
					pmpV = Q.colwise().mean();
				
					pmp_1 = pmpV(0);
					pmp_2 = pmpV(1);
					pmp_3 = pmpV(2);
 
					// variances - below is best way to calculate column wise
					pvp_1 = (Q.col(0).array() - Q.col(0).array().mean()).square().sum() / ((double)Q.col(0).array().size() - 1);
					pvp_2 = (Q.col(1).array() - Q.col(1).array().mean()).square().sum() / ((double)Q.col(1).array().size() - 1);
					pvp_3 = (Q.col(2).array() - Q.col(2).array().mean()).square().sum() / ((double)Q.col(2).array().size() - 1);

					// covariances - manual calculation 
					double sum12 = 0, sum13 = 0, sum23 = 0;

					for (int n = 0; n < Q.rows(); n++)
					{
						sum12 += (Q(n, 0) - pmp_1) * (Q(n, 1) - pmp_2);
						sum13 += (Q(n, 0) - pmp_1) * (Q(n, 2) - pmp_3);
						sum23 += (Q(n, 1) - pmp_2) * (Q(n, 2) - pmp_3);

					}
					N_SUBTRACT_ONE = Q.rows() - 1.0;

					pcov_12 = sum12 / N_SUBTRACT_ONE;
					pcov_13 = sum13 / N_SUBTRACT_ONE;
					pcov_23 = sum23 / N_SUBTRACT_ONE;

					term_1 = pmp_1 - omp_1,
						term_2 = pmp_2 - omp_2,
						term_3 = pmp_3 - omp_3,
						term_4 = pvp_1 - ovp_1,
						term_5 = pvp_2 - ovp_2,
						term_6 = pvp_3 - ovp_3,
						term_7 = pcov_12 - ocov_12,
						term_8 = pcov_13 - ocov_13,
						term_9 = pcov_23 - ocov_23;
					// note to self: I'm using vectorXd from now on b/c it plays way better than the vectors built into C++ unless ofc there are strings that we need to input.

					all_terms << term_1, term_2, term_3, term_4, term_5, term_6, term_7, term_8, term_9;
					term_vec = all_terms;

					// USE THE MOST RECENT ESTIMATE OF WMAT UNLESS USEEQUAL == 1
					double cost_newpos;
					if (useEqual == 0) {
						
						cost_newpos = term_vec.transpose() * w_mat * term_vec.transpose().transpose();
						
					}
					if (useEqual == 1) {
						cost_newpos = term_vec.transpose() * vectorOfOnes.asDiagonal() * term_vec.transpose().transpose();
						
					}
					if (cost_newpos < PBMAT(jjj, Npars)) {
						//cout << "line 913" << endl;
						VectorXd POSMAT_cost_newpos(POSMAT.cols());
						POSMAT_cost_newpos = POSMAT.row(jjj);
						POSMAT_cost_newpos.conservativeResize(PBMAT.cols());
						POSMAT_cost_newpos(PBMAT.cols() - 1) = cost_newpos;
						PBMAT.row(jjj) = POSMAT_cost_newpos;

						if (cost_newpos < cost_gbest) {
							gbest = POSMAT.row(jjj);
							cost_gbest = cost_newpos;
						}
					}
					
				}


				sfi = sfi - (sfe - sfg) / Nsteps; // reduce the inertial weight after each step
				sfs = sfs + (sfe - sfg) / Nsteps; // increase social weight after each step


				// CHECK IF NEW GLOBAL BEST HAS BEEN FOUND
				int neflag = 0;
				int lastrow = GBMAT.rows() - 1;
				
				for (int ne = 0; ne < Npars; ne++) {
					if (GBMAT(lastrow, ne) != gbest(ne)) {
						neflag = 1;
					}
				}

				if (GBMAT(lastrow, Npars) != cost_gbest) {
					neflag = 1;
				}

				//IF NEW GLOBAL BEST HAS BEEN FOUND, THEN UPDATE GBMAT
				if (neflag == 1) {
					//cout << "line 943" << endl;
					GBMAT.conservativeResize(GBMAT.rows() + 1, GBMAT.cols()); //rbind method currently.... where cbind is the "column bind vector"
					cbind << gbest, cost_gbest;
					GBMAT.row(GBMAT.rows() - 1) = cbind;
				}
				
			}

			if (pso < (Biter + 1)) {
				for (int init = 0; init < Npars; init++) { // best estimate of k to compute w.mat
					k.at(init) = gbest(init);
				}

				HM << -k.at(2), k.at(2), 0,
					k.at(1), -k.at(1) - k.at(4), k.at(4),
					k.at(3), k.at(0), -k.at(0) - k.at(3);

				HMT = t * HM.transpose();
				EHMT = HMT.exp();
				Q = (EHMT * X_0.transpose()).transpose();
				
				MatrixXd fmdiffs(N, 3);
				fmdiffs = Y_t - Q;

				VectorXd mxt(3);
				mxt = Q.colwise().mean();

				VectorXd myt(3);
				myt = Y_t.colwise().mean();

				MatrixXd residxt(Q.rows(), Q.cols());
				residxt.col(0) = mxt.row(0).replicate(N, 1);
				residxt.col(1) = mxt.row(1).replicate(N, 1);
				residxt.col(2) = mxt.row(2).replicate(N, 1);
				residxt = Q - residxt;

				MatrixXd residyt(Y_t.rows(), Y_t.cols());
				residyt.col(0) = myt.row(0).replicate(N, 1);
				residyt.col(1) = myt.row(1).replicate(N, 1);
				residyt.col(2) = myt.row(2).replicate(N, 1);
				residyt = Y_t - residyt;

				MatrixXd smdiffs(N, 3);
				smdiffs = (residyt.array() * residyt.array()) - (residxt.array() * residxt.array());
				
				MatrixXd cprxt(N, 3);
				cprxt.col(0) = residxt.col(0).array() * residxt.col(1).array();
				cprxt.col(1) = residxt.col(0).array() * residxt.col(2).array();
				cprxt.col(2) = residxt.col(1).array() * residxt.col(2).array();

				MatrixXd cpryt(N, 3);
				cpryt.col(0) = residyt.col(0).array() * residyt.col(1).array();
				cpryt.col(1) = residyt.col(0).array() * residyt.col(2).array();
				cpryt.col(2) = residyt.col(1).array() * residyt.col(2).array();

				MatrixXd cpdiffs(N, 3);
				cpdiffs = cpryt - cprxt;

				MatrixXd Adiffs(N, 9);

				Adiffs << fmdiffs, smdiffs, cpdiffs; // concatenate

				MatrixXd g_mat(N, Nterms);
				g_mat = Adiffs;

				w_mat.setZero();
				for (int m = 0; m < N; m++) { w_mat = w_mat + (g_mat.row(m).transpose()) * g_mat.row(m); }
				w_mat = (w_mat/N).inverse();
			
				if (useDiag == 1) { w_mat = w_mat.diagonal().diagonal(); }

				// Update cost_gbest with w_mat

				for (int init = 0; init < Npars; init++) { // recompute the cost for seedk using this w.mat
					k.at(init) = gbest(init);
				}

				HM << -k.at(2), k.at(2), 0,
					k.at(1), -k.at(1) - k.at(4), k.at(4),
					k.at(3), k.at(0), -k.at(0) - k.at(3);

				HMT = t * HM.transpose();
				EHMT = HMT.exp();
				Q = (EHMT * X_0.transpose()).transpose();


				// CALCULATE MEANS, VARIANCES, AND COVARIANCES
				pmpV = Q.colwise().mean();

				pmp_1 = pmpV(0);
				pmp_2 = pmpV(1);
				pmp_3 = pmpV(2);

				pvp_1 = 0;
				pvp_2 = 0;
				pvp_3 = 0;
				// variances - manually calculate it, no easy library 
				pvp_1 = (Q.col(0).array() - Q.col(0).array().mean()).square().sum() / ((double)Q.col(0).array().size() - 1);
				pvp_2 = (Q.col(1).array() - Q.col(1).array().mean()).square().sum() / ((double)Q.col(1).array().size() - 1);
				pvp_3 = (Q.col(2).array() - Q.col(2).array().mean()).square().sum() / ((double)Q.col(2).array().size() - 1);
				// covariances - manual calculation 
				double sum12 = 0, sum13 = 0, sum23 = 0;

				for (int n = 0; n < Q.rows(); n++)
				{
					sum12 += (Q(n, 0) - pmp_1) * (Q(n, 1) - pmp_2);
					sum13 += (Q(n, 0) - pmp_1) * (Q(n, 2) - pmp_3);
					sum23 += (Q(n, 1) - pmp_2) * (Q(n, 2) - pmp_3);

				}
				double N_SUBTRACT_ONE = Q.rows() - 1.0;

				pcov_12 = sum12 / N_SUBTRACT_ONE;
				pcov_13 = sum13 / N_SUBTRACT_ONE;
				pcov_23 = sum23 / N_SUBTRACT_ONE;

				term_1 = pmp_1 - omp_1,
					term_2 = pmp_2 - omp_2,
					term_3 = pmp_3 - omp_3,
					term_4 = pvp_1 - ovp_1,
					term_5 = pvp_2 - ovp_2,
					term_6 = pvp_3 - ovp_3,
					term_7 = pcov_12 - ocov_12,
					term_8 = pcov_13 - ocov_13,
					term_9 = pcov_23 - ocov_23;
				// note to self: I'm using vectorXd from now on b/c it plays way better than the vectors built into C++ unless ofc there are strings that we need to input.

				all_terms << term_1, term_2, term_3, term_4, term_5, term_6, term_7, term_8, term_9;
				term_vec = all_terms;

				cost_gbest = term_vec.transpose() * w_mat * (term_vec.transpose()).transpose();

				GBMAT.conservativeResize(GBMAT.rows() + 1, GBMAT.cols());
				VectorXd cbind(gbest.size() + 1);
				cbind << gbest, cost_gbest;
				GBMAT.row(GBMAT.rows() - 1) = cbind;

				if (cost_gbest < cost_sofar) {
					best_sofar = gbest;
					cost_sofar = cost_gbest;
				}
			}

			if (bsi == 0 || q == 1) {
				if (pso < (Biter + 1)) {
					cout << "blindpso+cost.est:" << endl << best_sofar << endl << cost_sofar << endl << endl;
				}
				if (pso == (Biter + 1)) {
					cout << "igmme + cost.est:"<< endl << gbest << endl << cost_gbest  << endl << endl;
					cout << "w.mat" << w_mat.transpose()  << endl << endl;
				}
			}
			if (bsi == 1 && q > 1) {
				if (pso == 1) {
					cout << gbest << cost_gbest << "ubsreps+cost.mat" << endl << endl;
				}
				if (pso > 1) {
					cout << gbest << cost_gbest << "wbsreps+cost.mat" << endl << endl;
				} 
			} 
		}  // end loop over PSO layers */
		

	} // end loop over NIter simulations
	cout << "GBMAT: " << endl;
	cout << GBMAT << endl;
	auto t2 = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count();
	cout << "CODE FINISHED RUNNING IN "<< duration<< " s TIME!" << endl;


	return 0; // just to close the program at the end.
}







/* generate random individual values */

/*double uniformRandDouble(double range_from, double range_to) {
	std::random_device rand_dev;
	std::mt19937 generator(rand_dev());
	std::uniform_int_distribution<double>    distr(range_from, range_to);
	return distr(generator);
}*/














// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
