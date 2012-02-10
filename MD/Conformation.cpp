/*
 * Conformation.cpp
 *
 *  Created on: Feb 7, 2012
 *      Author: alfaceor
 */

#include "Conformation.h"

Conformation::Conformation() {
	int d=3;
	for(int i=0;i<N;i++){
		for(int j=0;j<d;j++){
			if(j==0) chain[i].vec_r[j] = i*chain[i].zigma;
			else	 chain[i].vec_r[j] = 0;
		}
	}
}

Conformation::~Conformation() {
	// TODO Auto-generated destructor stub
}

void Conformation::GenDeltas(){
	// Generate a top diagonal matrix x_ij
	for (int i=0;i<N-1;i++){
		for (int j=i+1; j<N;j++){
			deltaX[i][j]	= 	chain[i].vec_r[0] - chain[j].vec_r[0];
			deltaY[i][j]	= 	chain[i].vec_r[1] - chain[j].vec_r[1];
			deltaZ[i][j]	= 	chain[i].vec_r[2] - chain[j].vec_r[2];
			deltaR2[i][j]	=	deltaX[i][j]*deltaX[i][j]+
								deltaY[i][j]*deltaY[i][j]+
								deltaZ[i][j]*deltaZ[i][j];
			deltaR[i][j]	=	sqrt(deltaR2[i][j]);
		}
	}

}

double Conformation::GaussianForce(double mean, double var){

		const gsl_rng_type * T;
		gsl_rng * r;
		gsl_rng_env_setup();
		T = gsl_rng_default;
		r = gsl_rng_alloc(T);
		// FIXME: complete this function to generate random forces.
		double randomNumber;
		for (int i=0;i<10;i++){
			//randomNumber = (rand() % 10) / 10.0;

			randomNumber = gsl_ran_gaussian(r,1);
			printf("%f\n",randomNumber);
		}
		return randomNumber;

}

double Conformation::Force_cc(int i,int j, double r_ij, double epsilon){
	double force;
	double r_inv6,r_inv8;
	// calculate r_ij e seus exponenetes
	r_inv6=1.0/pow(r_ij,6);
	r_inv8=r_inv6/pow(r_ij,2);
	double x_i=0.3, x_j=0.5;

	if(r_ij<=1){
		printf("r_ij <= 1");
		force = 12.0*epsilon*(r_inv6-1)*r_inv8*(x_i - x_j);
		return force;
	}else{
		printf("r_ij > 1");
		force = 0;
		return force;
	}
}

void Conformation::nextStep(){
	// TODO make the evolution of the systems
	// calculate the force
	// update the positions and velocities
	//
}

void Conformation::print_r(){
	for(int i=0;i<N;i++){
		chain[i].print_r(); printf("\n");
	}
}

