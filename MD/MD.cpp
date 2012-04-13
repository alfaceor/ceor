//============================================================================
// Name        : MD.cpp
// Author      : Carlos Olivares
// Version     : 0.5
// Copyright   : GPL
// Description : Molecular dynamics simulation
//============================================================================

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include "Conformation.h"


int main(int argc, char* argv[]) {
	// input parameters
	if(argc != 9){
		printf("Usage: MD [chain] [temperature] [total_time] [dt] [epsi] [q] [Ec] [print_each]\n");
		return 0;
	}

	const int M		=	strlen(argv[1]);//13;
	char *hydroChain;
	hydroChain =	argv[1];
	double temp		=	atof(argv[2]);	//0.04;
	int total_time	=	atoi(argv[3]);//200000;
	double dt	=	atof(argv[4]);//0.001;
	double epsi	=	atof(argv[5]);	//1.0;
	double q	=	atof(argv[6]); //0.1;
	double Ec	=	atof(argv[7]); //-1.0;
	int print_each	=	atoi(argv[8]);

	char filename_pattern[90];
	char filename_pdb[100];
	char filename_dat[100];
	char ext_pdb[]=".pdb";
	char ext_dat[]=".dat";

	strcpy(filename_pattern,"md-N");
	for (int i=1;i<argc;i++){
		strcat(filename_pattern,"_");
		strcat(filename_pattern,argv[i]);
	}

	strcat(filename_pdb,filename_pattern);
	strcat(filename_dat,filename_pattern);
	strcat(filename_pdb,ext_pdb);
	strcat(filename_dat,ext_dat);

	//----------------- Simulation
	Conformation protein(M, hydroChain, temp, filename_dat);
	protein.calculateTotalForces(epsi,q,Ec);
	FILE *fp;
	fp = fopen(filename_pdb,"w");
	printf("#%s\t%s\t%s\t%s\t%s\t%s\n","time","Energy", "KinecticEnergy", "PotentialEnergy", "Rg", "D");
	int ttime	= 0;
	while(ttime<total_time){
		protein.actualizePositions(dt);
		protein.calculateTotalForces(epsi,q,Ec);
		protein.actualizeVelocities(dt);
		protein.calculateTotalEnergy(epsi,q,Ec);

		if (ttime % print_each == 0){
			protein.print_pdb_conformation(fp,ttime);
			protein.calculateRg();
			protein.calculateD();
			printf("%d\t%f\t%f\t%f\t%f\t%f\n",ttime,protein.Energy, protein.KinecticEnergy, protein.PotentialEnergy, protein.Rg, protein.D);
		}
		ttime++;
	}
	fclose(fp);

	printf("# END SIMULATION\n");

	return EXIT_SUCCESS;
}
