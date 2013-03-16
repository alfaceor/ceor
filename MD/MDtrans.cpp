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

/**********
// Python External scripts
void pyplot_time_evolution(char *filename_dat){
	char *strcmd = (char *) malloc((strlen(filename_dat)+50)*sizeof(char));
	strcpy(strcmd,"grafica_evolucion_temporal_Dtrans.py ");
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

**********/

/*************************** MAIN ***************************/

int main(int argc, char* argv[]) {
	// input parameters
	if(argc != 12){
		printf("Usage: MD [rev_prefix] [chain] [temperature] [total_time] [dt] [epsi] [q] [Ec] [print_each] [Dfin] [Drate]\n");
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
	double Dfin	= atof(argv[10]);//1.0;
	double Drate	=	atof(argv[11]); // 1

	// FIXME: find a better way to put the cutoff
	char strdcutoff[]="1.2";
	double dcutoff = atof(strdcutoff); // cutoff for the Contact matrix

	char filename_pattern[100];
	char filename_pdb[100];
	char filename_dat[100];
	char filename_ini[100];
	char filename_con[100];
	char ext_pdb[]=".pdb";
	char ext_dat[]=".dat";
	char ext_ini[]=".ini";
	char ext_con[]=".con";

	strcpy(filename_pattern,prefix_file);
	for (int i=2;i<argc;i++){ // Don't print the Drate, Dfin
		strcat(filename_pattern,"_");
		strcat(filename_pattern,argv[i]);
	}

	strcpy(filename_ini,filename_pattern);
	strcat(filename_ini,ext_ini);

	//----------------- Simulation
	Conformation protein(M, hydroChain, temp, filename_ini);
	protein.calculateD();	// calculate initial D
	double Dini = protein.D;
	printf("M= %d,protein.N,Dini = %f\n",M,protein.D);

	// TODO: Validate the values of Dfin and Drate

	int ttime	= 0;

	const gsl_rng_type *T; T = gsl_rng_default;
	gsl_rng *r; r = gsl_rng_alloc(T);
	long unsigned int seed;
	seed = time(NULL); // Get the time of the system as seed
	gsl_rng_set(r,seed);

	// copy the filename pattern
	strcpy(filename_pdb,filename_pattern);
	strcpy(filename_dat,filename_pattern);
	strcpy(filename_con,filename_pattern);

	// add number sufix, transition sufix and file extension
	strcat(filename_pdb,ext_pdb);
	strcat(filename_dat,ext_dat);

	strcat(filename_con,"_d-");
	strcat(filename_con,strdcutoff);	//FIXME: dcutoff as string
	strcat(filename_con,ext_con);

	// open files to write data
	FILE *fp_pdb, *fp_dat,*fp_con;
	fp_pdb = fopen(filename_pdb,"w");
	fp_dat = fopen(filename_dat,"w");
	fp_con = fopen(filename_con,"w");

	protein.alingWithDaxis();

	fprintf(fp_dat,"#%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n","time","Energy", "KinecticEnergy", "PotentialEnergy", "Rg", "Z","HRg","PRg");
	fprintf(fp_con,"#%s\t%s\t%s\t%s\t%s\n","time","HHcontacts", "HPcontacts", "PPcontacts", "ALLcontacts");
	ttime	= 0;
	// Disminuye el valor de D
	if ( Dini > Dfin ){
		if (Drate > 0){ // If the Drate is positive then the D value will be decrease
			// Run the program
			while(protein.D > Dfin){ // D must decrease
				protein.calculateTotalForces(epsi,q,Ec);
				protein.actualizePositionsFixedEnds(dt);
				protein.addPosition3DNoiseFixedEnds(dt,temp,r);
				protein.actualizeVelocitiesFixedEnds(dt);
				protein.calculateTotalEnergy(epsi,q,Ec);
				protein.set_D_to(protein.D - Drate*dt);

				if (ttime % print_each == 0){
					protein.print_pdb_conformation(fp_pdb,ttime);
					protein.calculateRg();
					fprintf(fp_dat,"%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",ttime,protein.Energy, protein.KinecticEnergy, protein.PotentialEnergy, protein.Rg, protein.D, protein.HRg, protein.PRg);

					protein.binarizeDeltaR2(dcutoff);
					protein.calculateContacts();
					fprintf(fp_con,"%d\t%f\t%f\t%f\t%f\n",ttime, protein.HHcontacts, protein.HPcontacts, protein.PPcontacts, protein.ALLcontacts);

				}
				ttime++;
			}

		} else{
			printf("Bad paramaters: Drate must be POSITIVE\n");
			return EXIT_FAILURE;
		}
	} else if ( Dini < Dfin ){	// Aumenta el valor de D
		if (Drate < 0){
			// Run the program
			while(protein.D < Dfin){ // D must increase until reach Dfin
				protein.calculateTotalForces(epsi,q,Ec);
				protein.actualizePositionsFixedEnds(dt);
				protein.addPosition3DNoiseFixedEnds(dt,temp,r);
				protein.actualizeVelocitiesFixedEnds(dt);
				protein.calculateTotalEnergy(epsi,q,Ec);
				protein.set_D_to(protein.D - Drate*dt);

				if (ttime % print_each == 0){
					protein.print_pdb_conformation(fp_pdb,ttime);
					protein.calculateRg();
					fprintf(fp_dat,"%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",ttime,protein.Energy, protein.KinecticEnergy, protein.PotentialEnergy, protein.Rg, protein.D, protein.HRg, protein.PRg);
					protein.binarizeDeltaR2(dcutoff);
					protein.calculateContacts();
					fprintf(fp_con,"%d\t%f\t%f\t%f\t%f\n",ttime, protein.HHcontacts, protein.HPcontacts, protein.PPcontacts, protein.ALLcontacts);

				}
				ttime++;
			}

		} else{
			printf("Bad paramaters: Drate must be NEGATIVE\n");
			return EXIT_FAILURE;
		}
	}
	else{
		printf("fuck u!!");
	}

	// TODO: there must be a more efficient way to do that.
	// With this part of code the protein will end with D = Dfin
	protein.set_D_to(Dfin);

	protein.calculateTotalForces(epsi,q,Ec);
	protein.actualizePositionsFixedEnds(dt);
	protein.addPosition3DNoiseFixedEnds(dt,temp,r);
	protein.actualizeVelocitiesFixedEnds(dt);
	protein.calculateTotalEnergy(epsi,q,Ec);
	ttime++;

	protein.print_pdb_conformation(fp_pdb,ttime);
	protein.calculateRg();

	fprintf(fp_dat,"%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",ttime,protein.Energy, protein.KinecticEnergy, protein.PotentialEnergy, protein.Rg, protein.D, protein.HRg, protein.PRg);

	protein.binarizeDeltaR2(dcutoff);
	protein.calculateContacts();
	fprintf(fp_con,"%d\t%f\t%f\t%f\t%f\n",ttime, protein.HHcontacts, protein.HPcontacts, protein.PPcontacts, protein.ALLcontacts);
	// close transition files
	fclose(fp_pdb);
	fclose(fp_dat);
	fclose(fp_con);

	// time evolution plot
	// pyplot_time_evolution(filename_dat);
	//pyplot_transition_matrix(filename_dat);
	printf("# END SIMULATION\n");
	//shell_compress_pdbfiles(filename_pattern);

	return EXIT_SUCCESS;
}
