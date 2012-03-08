/*
 * Monomer.h
 *
 *  Created on: Feb 7, 2012
 *      Author: alfaceor
 */

#ifndef MONOMER_H_
#define MONOMER_H_

#include "stdio.h"
#define DIM 3
class Monomer {
public:
	Monomer();
	virtual ~Monomer();
	double vec_r[DIM];	// position vector
	double vec_v[DIM];	// velocity vector
	double zigma;		// monomer diameter
	double hydro;	// Hydropathy Index (negative polar && positive nonpolar)
	double mass;
	double total_force[DIM];
	void print_r();
	void print_v();
};

#endif /* MONOMER_H_ */
