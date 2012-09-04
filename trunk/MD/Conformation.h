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
	Conformation(int N, char *hydroChain, double temp, char *basename);
	virtual ~Conformation();
	Monomer *chain;	// Polymer chain
	double KinecticEnergy;	// Kinetic energy
	double PotentialEnergy;	// Potentatial energy
	double Energy;			// System Energy
	double Rg;				// Radius of Gyration: the rms distance of each atom to the centroid.
	double HRg;
	double PRg;
	double CenterMass[DIM];
	double D;				// End to End
	double *deltaR2;
	double Temperature;		// to get the temperature in Kelvin divide by Boltzmann constant
	void calculateDeltaR2();
	void calculateCenterMass();
	void calculateRg();
	void calculateD();
	void set_D_to(double Dnew);
	void calculateBondForces(double epsi, double q);
	void calculateBondPotential(double epsi, double q);
	void calculateHydroForces(double epsi, double Ec);
	void calculateHydroPotential(double epsi, double Ec);
	void calculateDampingForces();
	void calculateRandomForces();
	void calculateTotalForces(double epsi, double q, double Ec);
	void calculateTotalEnergy(double epsi, double q, double Ec);
//	void integratorVerlet(double dt);

	void setTemperature(double temp);
	void calculateTemperature();

	void actualizePositions(double dt);
	void actualizePositionsFixedEnds(double dt);
	void addPosition2DNoise(double dt, double KT, gsl_rng *r);
	void addPosition3DNoise(double dt, double KT, gsl_rng *r);
	void addPosition3DNoiseFixedEnds(double dt, double KT, gsl_rng *r);
	void actualizeVelocities(double dt);
	void calculateKineticEnergy();
	void cleanForces();
	void cleanEnergyValues();
	void randomPositions();
	void gaussianRandomVelocities(double temp);
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
