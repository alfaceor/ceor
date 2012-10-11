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

#define MAXCHARLEN 200

class MDfilenames{
public:
  MDfilenames(const int argc,char **argv);
  void setDfixed_filename(double Dvalue);
  void onejobfile();
  char filename_pattern[MAXCHARLEN];
  char filename_pdb[MAXCHARLEN];
  char filename_dat[MAXCHARLEN];
  char filename_ini[MAXCHARLEN];
  char ext_pdb[5];
  char ext_dat[5];
  char ext_ini[5];

  char letter_prefix[MAXCHARLEN];
  char token[5];
};

// TODO: Aqui debo el numero de argumentos que entran
MDfilenames::MDfilenames(const int argc, char **argv){
  // Init file extensions
  strcpy(ext_pdb,".pdb");
  strcpy(ext_dat,".dat");
  strcpy(ext_ini,".ini");

  strcpy(letter_prefix,"_ict____");
   /*
  i: idjob
  c: chain
  t: temperature
  */
  strcpy(token,"__");
  strcpy(filename_pattern,argv[0]);

  for (int i=1; i<argc; i++){
    token[1]=letter_prefix[i];
    strcat(filename_pattern,token);
    strcat(filename_pattern,argv[i]);
  }
  strcpy(filename_ini,filename_pattern);
  strcat(filename_ini,ext_ini);
}

void MDfilenames::onejobfile(){
	strcpy(filename_pdb,filename_pattern);strcat(filename_pdb,ext_pdb);
	strcpy(filename_dat,filename_pattern);strcat(filename_dat,ext_dat);
}

void MDfilenames::setDfixed_filename(double Dvalue){
	strcpy(filename_pdb,filename_pattern);
	strcpy(filename_dat,filename_pattern);
	char strDvalue[MAXCHARLEN]; sprintf(strDvalue,"_D%3.3f",Dvalue);
	strcat(filename_pdb,strDvalue); strcat(filename_pdb,ext_pdb);
	strcat(filename_dat,strDvalue); strcat(filename_dat,ext_dat);
}


void write_header_dat(FILE *fp_dat, char *exec_line, int argc, char **argv){
	// header INPUT data
	fprintf(fp_dat,"# %s",exec_line);
	fprintf(fp_dat,"# ");
	for (int i=0; i<argc; i++){
		fprintf(fp_dat,"%s\t",argv[i]);
	}
	fprintf(fp_dat,"\n");
	fprintf(fp_dat,"#%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n","time","Energy", "KinecticEnergy", "PotentialEnergy", "Rg", "D","HRg","PRg");
}



/*************************************/
/***************  MAIN ***************/
/*************************************/
int main(int argc, char* argv[]) {
	// INPUT PARAMS
	char exec_line[]="MD\t[prefix_file]\t[chain]\t[temp]\t[total_time]\t[dt]\t[epsi]\t[q]\t[Ec]\t[print_each]\n";
	int input_args=10;
	if(argc != input_args){
		printf("Usage: %s",exec_line);
		printf("min arguments = %d\n",input_args);
		return EXIT_FAILURE;
	}

	char *prefix_file, *hydroChain;
	prefix_file	=	argv[1];
	hydroChain	=	argv[2];
	const int M		=	strlen(argv[2]);//13;
	double temp		=	atof(argv[3]);	//0.04;
	int total_time	=	atoi(argv[4]);	//200000;
	double dt	=	atof(argv[5]);	//0.001;
	double epsi	=	atof(argv[6]);	//1.0;
	double q	=	atof(argv[7]);	//0.1;
	double Ec	=	atof(argv[8]);	//-1.0;
	int print_each	=	atoi(argv[9]);

	// FIXME: Times as arguments
	int ttime	= 0;
	double Drate=1;
	double ttrans=total_time;

	// D values to analize
	double Dwell[3]={3.0, 2.0, 1.0};

	// Initiate randomness
	const gsl_rng_type *T; T = gsl_rng_default;
	gsl_rng *r; r = gsl_rng_alloc(T);
	long unsigned int seed;
	seed = time(NULL); // Get the time of the system as seed
	gsl_rng_set(r,seed);

	MDfilenames nmfiles(4,argv);
	//----------------- Simulation
	Conformation protein(M, hydroChain, temp, nmfiles.filename_ini);
	FILE *fp_pdb, *fp_dat;
	protein.calculateD();
	for (int j=0;j<3;j++){
		// Hacer una transicion hasta el Dwell[0] a un Drate dado
		while(protein.D > Dwell[j]){
			protein.calculateTotalForces(epsi,q,Ec);
			protein.actualizePositionsFixedEnds(dt);		// FIXME: fix positions to the ends
			protein.addPosition3DNoiseFixedEnds(dt,temp,r);	// FIXME: remove ends from perturbation
			protein.actualizeVelocitiesFixedEnds(dt);
			protein.calculateTotalEnergy(epsi,q,Ec);
			protein.set_D_to(protein.D - Drate*dt);
			printf("%f\n",protein.D);
		}

		nmfiles.setDfixed_filename(Dwell[j]);
		fp_pdb = fopen(nmfiles.filename_pdb,"w");
		fp_dat = fopen(nmfiles.filename_dat,"w");

		// header INPUT data
		write_header_dat(fp_dat,exec_line,argc,argv);
		// Dejar la simulacion ahi correr por un tiempo largo
		ttime = 0;
		int taux = ttrans/100;
		time_t tstart = time(0);
		time_t tend;
		while(ttime < ttrans){ // FIXME: LONG TIME
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
//			if (ttime % taux == 0){
//				tend = time(0);
//				printf("%lf\n",tend);
//			}
			ttime++;
		};
		fclose(fp_pdb);
		fclose(fp_dat);
	}

	printf("# END SIMULATION\n");

	return EXIT_SUCCESS;
}
