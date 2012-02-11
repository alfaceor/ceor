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
#define N 3	// FIXME: make a general chain!!!!!!

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
	double force[N][3];

	void GenDeltas();	// Generate a diagonal matrix x_ij

	double GaussianForce(double mean,double var);
	void Force_cc(double epsilon, double q);
	double Force_cc2(int i,int j, double r_ij, double epsilon);
	void nextStep();	// time evolution, TODO: CLEAN force value
	void print_r();		// Print data
//	void print_pdbfile(char *filename);
	void print_pdb_line(int serial);

};

#endif /* CONFORMATION_H_ */
