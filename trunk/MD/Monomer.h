/*
 * Monomer.h
 *
 *  Created on: Feb 7, 2012
 *      Author: alfaceor
 */

#ifndef MONOMER_H_
#define MONOMER_H_

#include "util_functions.h"

class Monomer {
public:
	Monomer();
	virtual ~Monomer();
	double vec_r[DIM];	// position vector
	double vec_v[DIM];	// velocity vector
	double zigma;		// monomer diameter
	double hydro;	// Hydropathy Index (negative polar && positive nonpolar)
	double mass;
	double gamma;	// Viscosius constant
	double total_force[DIM];
	double total_force_old[DIM];
	void actualizeVec_r(double dt);
	void actualizeVec_v(double dt);
	void print_r();
	void print_v();
};

#endif /* MONOMER_H_ */
