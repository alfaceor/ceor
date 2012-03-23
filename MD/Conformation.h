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
	double KinecticEnergy;	// Kinetic energy
	double PotentialEnergy;	// Potentatial energy
	double Energy;			// System Energy
	double Rg;				// Radius of Gyration
	double D;				// End to End
	double *deltaR2;
	void calculateDeltaR2();
	void calculateBondForces(double epsi, double q);
	void calculateBondPotential(double epsi, double q);
	void calculateHydroForces(double epsi, double Ec);
	void calculateHydroPotential(double epsi, double Ec);
	void calculateTotalForces(double epsi, double q, double Ec);
	void calculateTotalEnergy(double epsi, double q, double Ec);
//	void integratorVerlet(double dt);
	void actualizePositions(double dt);
	void actualizeVelocities(double dt);
	void calculateKineticEnergy();
	void cleanForces();
	void cleanEnergyValues();
	void randomPositions();
	void gaussianRandomVelocities();

	void nextStep();	// time evolution, TODO: CLEAN force value
	void print_r();		// Print data
	void print_pdb_conformation(FILE *fp, int time_model);
//	void print_pdbfile(char *filename);
	void print_pdb_line(FILE *fp,int serial, double x, double y, double z,char *resName,double tempFactor);

private:
	int N;
	char *basename;


};

#endif /* CONFORMATION_H_ */
