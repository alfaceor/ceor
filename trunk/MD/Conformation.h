/*
 * Conformation.h
 *
 *  Created on: Feb 7, 2012
 *      Author: alfaceor
 */

#ifndef CONFORMATION_H_
#define CONFORMATION_H_

#include "Monomer.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>

#define DIM 3

class Conformation {

public:
	Conformation(int N);
	Conformation(int N, char *basename);
	virtual ~Conformation();
	Monomer *chain;	// Polymer chain
	double Energy;		// System Energy
	double Rg;			// Radius of Gyration
	double D;			// End to End
	double *deltaR2;
	void calculateDeltaR2();
	void calculateBondForces(double epsi, double q);
	void calculateHydroForces(double epsi, double Ec);
	void calculateTotalForces(double epsi, double q, double Ec);
	void integratorVerlet(double dt);
	void cleanForces();

	void nextStep();	// time evolution, TODO: CLEAN force value
	void print_r();		// Print data
//	void print_pdbfile(char *filename);
	void print_pdb_line(int serial);
private:
	int N;
	char *basename;


};

#endif /* CONFORMATION_H_ */
