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

/*
// Python External scripts
void pyplot_time_evolution(char *filename_dat){
	char *strcmd = (char *) malloc((strlen(filename_dat)+50)*sizeof(char));
	strcpy(strcmd,"grafica_evolucion_temporal_Dfixed.py ");
	strcat(strcmd,filename_dat);
//	free(strcmd);
	system(strcmd);
}

void pyplot_transition_matrix(char *filename_dat){
	char *strcmd = (char *) malloc((strlen(filename_dat)+50)*sizeof(char));
	strcpy(strcmd,"TransitionMatrix.py ");
	strcat(strcmd,filename_dat);
//	free(strcmd);
	system(strcmd);
}


// Shell tools
void shell_compress_pdbfiles(char *filename_pattern){
	char *strcmd = (char *) malloc((strlen(filename_pattern)+50)*sizeof(char));
	strcpy(strcmd,"tar cvzf pdbfiles");
	strcat(strcmd,filename_pattern);
	strcat(strcmd,".tar.gz *.pdb");
	system(strcmd);
//	free(strcmd);
	system("rm *.pdb");
}
*/


/*************************** MAIN ***************************/

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

	char filename_pattern[100];
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

	strcpy(filename_ini,filename_pattern);
	strcat(filename_ini,ext_ini);

	//----------------- Simulation
	Conformation protein(M, hydroChain, temp, filename_ini);
	protein.calculateD();	// calculate initial D
	double Dmax = protein.D;
	double Dmin = 1.0;
	int ttime	= 0;

	const gsl_rng_type *T; T = gsl_rng_default;
	gsl_rng *r; r = gsl_rng_alloc(T);
	long unsigned int seed;
	seed = time(NULL); // Get the time of the system as seed
	gsl_rng_set(r,seed);

		// open files to write data
		FILE *fp_pdb, *fp_dat;

		// copy the filename pattern
		strcpy(filename_pdb,filename_pattern);
		strcpy(filename_dat,filename_pattern);

		// add number sufix, transition sufix and file extension
		strcat(filename_pdb,ext_pdb);
		strcat(filename_dat,ext_dat);

		// open files to write data
		fp_pdb = fopen(filename_pdb,"w");
		fp_dat = fopen(filename_dat,"w");

		fprintf(fp_dat,"#%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n","time","Energy", "KinecticEnergy", "PotentialEnergy", "Rg", "D","HRg","PRg");
		// FIXME: TENGO que hacer una translacion y luego una rotacion

		protein.print_pdb_conformation(fp_pdb,-2);
		ttime=0;
		//protein.displace(-protein.chain[0].vec_r[0],-protein.chain[0].vec_r[1],-protein.chain[0].vec_r[2]);
		protein.print_pdb_conformation(fp_pdb,-1);
		// calcular el angulo a rotar en polar
		protein.alingWithDaxis();

		protein.print_pdb_conformation(fp_pdb,0);

		while(ttime<total_time){
			protein.calculateTotalForces(epsi,q,Ec);
			protein.actualizePositionsFixedEnds(dt);		// FIXME: fix positions to the ends
			protein.addPosition3DNoiseFixedEnds(dt,temp,r);	// FIXME: remove ends from perturbation
			protein.actualizeVelocitiesFixedEnds(dt);
			protein.calculateTotalEnergy(epsi,q,Ec);

			if (ttime % print_each == 0){
				protein.print_pdb_conformation(fp_pdb,ttime);
				protein.calculateRg();
				fprintf(fp_dat,"%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",dt*ttime,protein.Energy, protein.KinecticEnergy, protein.PotentialEnergy, protein.Rg, protein.D, protein.HRg, protein.PRg);
			}
			ttime++;
		}


		// close D fixed files
		fclose(fp_pdb);
		fclose(fp_dat);

	printf("# END SIMULATION\n");
//	shell_compress_pdbfiles(filename_pattern);

	return EXIT_SUCCESS;

}
