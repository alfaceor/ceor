//============================================================================
// Name        : MD.cpp
// Author      : Carlos Olivares
// Version     : 0.5
// Copyright   : GPL
// Description : Molecular dynamics simulation
//============================================================================

// TODO: Add function to calculate the matrix contacts
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
		printf("Usage: MD2d [rev_prefix] [chain] [temperature] [total_time] [dt] [epsi] [q] [Ec] [print_each]\n");
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

	// FIXME: find a better way to put the cutoff
	char strdcutoff[]="1.2";
	double dcutoff = atof(strdcutoff); // cutoff for the Contact matrix

	// TODO: Add fin file
	char filename_pattern[90];
	char filename_pdb[100];
	char filename_dat[100];
	char filename_ini[100];
	char filename_con[100];
	char filename_seed[100];
	char ext_pdb[]=".pdb";
	char ext_dat[]=".dat";
	char ext_ini[]=".ini";
	char ext_con[]=".con";
	char ext_seed[]=".seed";

	strcpy(filename_pattern,prefix_file);
	for (int i=2;i<argc;i++){
		strcat(filename_pattern,"_");
		strcat(filename_pattern,argv[i]);
	}

	strcpy(filename_pdb,filename_pattern);
	strcpy(filename_dat,filename_pattern);
	strcpy(filename_ini,filename_pattern);
	strcpy(filename_con,filename_pattern);
	strcpy(filename_seed,filename_pattern);

	strcat(filename_pdb,ext_pdb);
	strcat(filename_dat,ext_dat);
	strcat(filename_ini,ext_ini);
	strcat(filename_seed,ext_seed);

	strcat(filename_con,"_d-");
	strcat(filename_con,strdcutoff);	//FIXME: dcutoff as string
	strcat(filename_con,ext_con);

	// Seed for Random numbers
	const gsl_rng_type *T; T = gsl_rng_default;
	gsl_rng *r; r = gsl_rng_alloc(T);
	long unsigned int seed;
	seed = (long unsigned int) time(NULL); // Get the time of the system as seed
	gsl_rng_set(r,seed);


	//protein.calculateTotalForces(epsi,q,Ec);
	FILE *fp_pdb, *fp_dat, *fp_con, *fp_seed;
	fp_pdb = fopen(filename_pdb,"w");
	fp_dat = fopen(filename_dat,"w");
	fp_con = fopen(filename_con,"w");
	fp_seed = fopen(filename_seed,"w");

	// FIXME: to save space in the hard disk i'm not printing Total Energy, KinecticEnergy, HRg, PRg
	fprintf(fp_dat,"#%s\t%s\t%s\t%s\t%s\t%s\n",
		"time", "PotentialEnergy", "Rg",
		"Z","HRg","PRg");
	fprintf(fp_con,"#%s\t%s\t%s\t%s\t%s\n",
		"time","HHcontacts", "HPcontacts",
		"PPcontacts", "ALLcontacts");
	fprintf(fp_seed,"seed %lu\n",seed);

	// ========> BEGIN Simulation -----------------
	Conformation protein(M, hydroChain, temp, filename_ini);
	int ttime	= 0;
	int ttaux = total_time/100;
	protein.print_pdb_unfoldedchain(fp_pdb); // First Positions just estetic function

	while(ttime<total_time){
		protein.calculateTotalForces(epsi,q,Ec);
		protein.actualizePositions(dt);
		protein.addPosition2DNoise(dt,temp,r);
		protein.actualizeVelocities(dt);
		protein.calculateTotalEnergy(epsi,q,Ec);

		if ( ttime % print_each == 0 ){
			protein.calculateRg();
			protein.calculateD();
			fprintf(fp_dat,"%d\t%f\t%f\t%f\t%f\t%f\n",
				ttime, protein.PotentialEnergy, protein.Rg,
				protein.D, protein.HRg, protein.PRg);
			// FIXME: Review the number of contacts, something is WRONG!
			protein.binarizeDeltaR2(dcutoff);
			protein.calculateContacts();
			fprintf(fp_con,"%d\t%f\t%f\t%f\t%f\n",
				ttime, protein.HHcontacts, protein.HPcontacts,
				protein.PPcontacts, protein.ALLcontacts);

		}
		if (ttime%ttaux==0) protein.print_pdb_conformation(fp_pdb,ttime);
		ttime++;
	}
	protein.print_pdb_conformation(fp_pdb,ttime); // last positions

	//========>  END Simulation -----------------

	long unsigned int finaltime;
	finaltime = (long unsigned int) time(NULL);
	fprintf(fp_seed,"End %lu",finaltime);

	fclose(fp_pdb);
	fclose(fp_dat);
	fclose(fp_con);
	fclose(fp_seed);
	// first is printed the existence of an ini file
	printf("%s # END_SIMULATION %lu %lu\n",filename_pattern,seed,finaltime);

	return EXIT_SUCCESS;
}
