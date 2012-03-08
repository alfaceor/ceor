//============================================================================
// Name        : MD.cpp
// Author      : Carlos Olivares
// Version     : 0.5
// Copyright   : GPL
// Description : Molecular dynamics simulation, Ansi-style
//============================================================================

#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>


#include "Conformation.h"

/**
 * @brief Function to generate initial random conditions
 * @param chain_r 2d array for save positions
 * @param chain_v 2d array for save velocities
 * @param hydro 1d array for save hydrophobicity values
 * @param filename filename of the file to save this initial values
 * @param M Number of particles (first dimension of the 2d array for positions and velocities)
 * */
int generateInitialConditions(double *chain_r,double *chain_v, double *hydro, char *filename, const int M){
	FILE *fop;
	fop = fopen(filename,"w");
	double sigma=1.0;

	// File header
	fprintf(fop,"N = %d\n",M);
	fprintf(fop,"x\ty\tz\tvx\tvy\tvz\thydro\n");

	// random variables
	gsl_rng *r1,*r2;
	const gsl_rng_type * T;
	T = gsl_rng_default;
	r1 = gsl_rng_alloc(T);
	r2 = gsl_rng_alloc(T);
	double randomaux;
	for (int i=0;i<M;i++){
		for (int d=0;d<3;d++){
			// vector force
//			total_force[i*DIM+d]=0;
			// vector positions and velocities
			chain_r[i*DIM+d] = 0.0;
			chain_v[i*DIM+d] = 0.0;
			if (d==0){
				randomaux = 0.001*gsl_rng_uniform_pos(r1);
				chain_r[i*DIM+d] = randomaux;
			}
			// vector velocities in any direction.
			randomaux = 0.5*(2*gsl_ran_gaussian(r2,sigma) - 1);
			chain_v[i*DIM+d] = randomaux;
			// la misma cadena que usa Lois 2008
		}
		//-------------------------
		if(i==0 || i==4 || i==8 || i==12) hydro[i]=-1.0;
		else			hydro[i]=1.0;
		chain_r[i*DIM+0] += sigma*i;
		fprintf(fop,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", chain_r[i*DIM+0], chain_r[i*DIM+1], chain_r[i*DIM+2], chain_v[i*DIM+0], chain_v[i*DIM+1], chain_v[i*DIM+2], hydro);
//		printf("%f\t%f\n",chain_r[i*N+0],chain_v[i*N+0]);
	}

	fclose(fop);
	return true;
}

int main(int argc, char* argv[]) {
	char *filename;
	filename="data_inicial.dat";
	const int M=13;

	Conformation protein(M,filename);
	printf("--------------");
	protein.calculateDeltaR2();
	for (int i=0;i<M;i++){
		for (int j=0; j<M; j++){
			printf("%lf",protein.deltaR2[i*DIM+j]);
		}
		printf("\n");
	}


//	protein.chain[1].print_r();

	/**
	 * TODO:
	 * -Calculate the distance matrix. deltaR a triangular matrix.
	 * -Calculate the force (aceleration) to change the velocity
	 * -Then change the position.
	 * */
	// Calculate

	return EXIT_SUCCESS;
}
