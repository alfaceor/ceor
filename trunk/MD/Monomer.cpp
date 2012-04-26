/*
 * Monomer.cpp
 *
 *  Created on: Feb 7, 2012
 *      Author: alfaceor
 */

#include "Monomer.h"

Monomer::Monomer() {
	zigma=1;
	mass=1;
	hydro=1;
	gamma=0.001;
	for (int d=0; d<DIM; d++){
		total_force[d]=0.0;
		total_force_old[d]=0.0;
		vec_r[d]=0.0;
		vec_v[d]=0.0;
	}
}

Monomer::~Monomer() {
	// TODO Auto-generated destructor stub
}

void Monomer::actualizeVec_v(double dt){
	for (int d=0; d<DIM; d++){
		//vec_v[d] += 0.5*dt*(total_force[d]+total_force_old[d])/mass;
		vec_v[d] = (vec_r[d] - vec_r_old[d])/dt;
	}
}

void Monomer::actualizeVec_r(double dt){
	for (int d=0; d<DIM; d++){
		//vec_r[d] += dt*vec_v[d]+0.5*(dt*dt)*total_force[d]/mass;
		vec_r_old[d] = vec_r[d];
		vec_r[d] += dt*total_force[d];
	}
}

void Monomer::addPositionNoise(double dt, double KT, gsl_rng *r){
	double etha;
	for (int d=0; d<DIM; d++){
		// Generate the gaussian random number
		etha=gsl_ran_gaussian(r,1);
		vec_r[d] += sqrt(2*KT*dt)*etha;
	}

}

void Monomer::print_r(){
	for (int i=0;i<3;i++){
		printf("%f\t",vec_r[i]);
	}
}

void Monomer::print_v(){
	for (int i=0;i<3;i++){
		printf("%f\t",vec_v[i]);
	}
}
