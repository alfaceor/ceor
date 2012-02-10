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
#define N 3

class Conformation {
public:
	Conformation();
	virtual ~Conformation();
	Monomer chain[N];	// Polymer chain
	double Energy;		// System Energy
	double Rg;			// Radius of Gyration
	double D;			// End to End

	double deltaX[N][N];
	double deltaY[N][N];
	double deltaZ[N][N];
	double deltaR2[N][N];	// r_ij^2
	double deltaR[N][N];	// r_ij

	void GenDeltas();	// Generate a diagonal matrix x_ij

	double GaussianForce(double mean,double var);
	double Force_cc(int i,int j, double r_ij, double epsilon);
	void nextStep();	// time evolution
	void print_r();		// Print data
	void print_pdbfile(char filename); //
};

#endif /* CONFORMATION_H_ */
