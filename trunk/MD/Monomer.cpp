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
	total_force[0]=0.0;
	total_force[1]=0.0;
	total_force[2]=0.0;
	vec_r[0]=0; vec_r[1]=0; vec_r[2]=0;
	vec_v[0]=0.1; vec_v[1]=0; vec_v[2]=0;
}

Monomer::~Monomer() {
	// TODO Auto-generated destructor stub
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
