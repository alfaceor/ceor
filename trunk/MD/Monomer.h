/*
 * Monomer.h
 *
 *  Created on: Feb 7, 2012
 *      Author: alfaceor
 */

#ifndef MONOMER_H_
#define MONOMER_H_

#include "stdio.h"

class Monomer {
public:
	Monomer();
	virtual ~Monomer();
	double vec_r[3];	// position vector
	double vec_v[3];	// velocity vector
	double zigma;		// monomer diameter
	double hydroIndex;	// Hydropathy Index (negative polar && positive nonpolar)
	double mass;
	void print_r();
	void print_v();
};

#endif /* MONOMER_H_ */
