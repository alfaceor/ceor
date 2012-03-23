//============================================================================
// Name        : MD.cpp
// Author      : Carlos Olivares
// Version     : 0.5
// Copyright   : GPL
// Description : Molecular dynamics simulation
//============================================================================

#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include "Conformation.h"

int main(int argc, char* argv[]) {
	// input parameters

	const int M	=	13;
	double epsi	=	 1.0;
	double q	=	 0.1;
	double Ec	=	-1.0;
	double temp	=	0.04;
	double dt	=	0.001;
	int total_time = 200000;

	char *filename;	filename="data_inicial.dat";
	char *pdbfile;	pdbfile	="protein.pdb";

	// Simulation
	Conformation protein(M,filename);
	protein.calculateTotalForces(epsi,q,Ec);
	FILE *fp;
	fp = fopen(pdbfile,"w");

	int time	= 0;
	while(time<total_time){
		protein.actualizePositions(dt);
		protein.calculateTotalForces(epsi,q,Ec);
		protein.actualizeVelocities(dt);
		protein.calculateTotalEnergy(epsi,q,Ec);

		if (time % 10000 == 0){
			protein.print_pdb_conformation(fp,time);
			printf("%d\t%f\t%f\t%f\n",time,protein.Energy, protein.KinecticEnergy, protein.PotentialEnergy);
		}
		time++;
	}
	fclose(fp);

	printf("END SIMULATION\n");

	return 1;
}
