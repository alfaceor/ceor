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
#include <time.h>

#include "Conformation.h"


int main(int argc, char* argv[]) {
	// input parameters
	if(argc != 10){
		printf("Usage: MD [rev_prefix] [chain] [temperature] [total_time] [dt] [epsi] [q] [Ec] [print_each]\n");
		return EXIT_FAILURE;
	}

	char *prefix_file, *hydroChain;
	prefix_file	=	argv[1];
	hydroChain	=	argv[2];
	const int M		=	strlen(argv[2]);//13;
	double temp		=	atof(argv[3]);	//0.04;
	int total_time	=	atoi(argv[4]);//200000;
	double dt	=	atof(argv[5]);//0.001;
	double epsi	=	atof(argv[6]);	//1.0;
	double q	=	atof(argv[7]); //0.1;
	double Ec	=	atof(argv[8]); //-1.0;
	int print_each	=	atoi(argv[9]);

	char filename_pattern[90];
	char filename_pdb[100];
	char filename_dat[100];
	char filename_ini[100];
	char ext_pdb[]=".pdb";
	char ext_dat[]=".dat";
	char ext_ini[]=".ini";

	strcpy(filename_pattern,prefix_file);
	for (int i=2;i<argc;i++){
		strcat(filename_pattern,"_");
		strcat(filename_pattern,argv[i]);
	}

	strcpy(filename_pdb,filename_pattern);
	strcpy(filename_dat,filename_pattern);
	strcpy(filename_ini,filename_pattern);
	strcat(filename_pdb,ext_pdb);
	strcat(filename_dat,ext_dat);
	strcat(filename_ini,ext_ini);

	//----------------- Simulation
	Conformation protein(M, hydroChain, temp, filename_ini);
	//protein.calculateTotalForces(epsi,q,Ec);
	FILE *fp_pdb, *fp_dat;
	fp_pdb = fopen(filename_pdb,"w");
	fp_dat = fopen(filename_dat,"w");
	fprintf(fp_dat,"#%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n","time","Energy", "KinecticEnergy", "PotentialEnergy", "Rg", "D","HRg","PRg");
	int ttime	= 0;

	const gsl_rng_type *T; T = gsl_rng_default;
	gsl_rng *r; r = gsl_rng_alloc(T);
	long unsigned int seed;
	seed = time(NULL); // Get the time of the system as seed
	gsl_rng_set(r,seed);

	while(ttime<total_time){
		protein.calculateTotalForces(epsi,q,Ec);
		protein.actualizePositions(dt);
		protein.addPosition3DNoise(dt,temp,r);
		protein.actualizeVelocities(dt);
		protein.calculateTotalEnergy(epsi,q,Ec);

		if (ttime % print_each == 0){
			//protein.print_pdb_conformation(fp_pdb,ttime);
			protein.calculateRg();
			protein.calculateD();
			fprintf(fp_dat,"%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",ttime,protein.Energy, protein.KinecticEnergy, protein.PotentialEnergy, protein.Rg, protein.D, protein.HRg, protein.PRg);
		}
		ttime++;
	}
	fclose(fp_pdb);
	fclose(fp_dat);

	printf("# END SIMULATION\n");

	return EXIT_SUCCESS;
}
