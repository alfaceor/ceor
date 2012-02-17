//============================================================================
// Name        : Tests.cpp
// Author      : carlos
// Version     :
// Copyright   : GPL
// Description : Hello World in C, Ansi-style
//============================================================================

#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include "forces.h"
#include "pdb_utils.h"


gsl_rng *r1,*r2;

int main(int argc, char* argv[]) {
	const int N = 13; // = atoi(argv[1]);
	const gsl_rng_type * T;
	T = gsl_rng_default;
	r1 = gsl_rng_alloc(T);
	r2 = gsl_rng_alloc(T);

	// phi_cc function
//	double dx=0.025;
//	double x;
	int time=0;
	int total_time=50;
	double dt=0.01;
	double total_force[N][3];
	double total_force_old[N][3];
	double chain_r[N][3];		// monomer positions
	double chain_v[N][3];		// monomer velocities

	int    hydro[N];	// = {1,1,-1,1,-1,-1,1,1,1,-1,1,1,1};
	double deltaR2[N][N];
	double epsi=1.0, Ec=-1.0, q=0.1;
	double sigma = 1.0;	// diametro


	//----- Initial positions and velocities
	double random;
	for (int i=0;i<N;i++){
		for (int d=0;d<N;d++){
			// vector force
			total_force[i][d]=0;
			// vector positions
			random = 0.001*gsl_rng_uniform_pos(r1);
			chain_r[i][d] = random;
			// vector velocities
			random = 0.5*(2*gsl_ran_gaussian(r2,sigma) - 1);
			chain_v[i][d] = random;
		}
		chain_r[i][0] += sigma*i;
	}

	// FIXME: Nose pq en el bucle no puede inicializar la variable hydro
	hydro[0]  =  1;
	hydro[1]  =  1;
	hydro[2]  = -1;
	hydro[3]  =  1;
	hydro[4]  = -1;
	hydro[5]  = -1;
	hydro[6]  =  1;
	hydro[7]  =  1;
	hydro[8]  =  1;
	hydro[9]  = -1;
	hydro[10] =  1;
	hydro[11] =  1;
	hydro[12] =  1;

	/****************************************/

FILE *fp;
fp = fopen("test.pdb", "w");
  while(time<total_time){

		fprintf(fp,"MODEL\t%d\n",time);
		char resname[]="   ";
		for (int i=0;i<N;i++){
			if(hydro[i] == 1)	resname[2]='C';
			else				resname[2]='N';
	//		printf("\n--- %d --- %s \n",hydro[i],resname);
			print_pdb_line(fp,i+1,chain_r[i][0],chain_r[i][1],chain_r[i][2],resname);
		}

	//----- Actualizacion de las posiciones
	for (int i=0;i<N;i++){
		for (int d=0;d<3;d++){
			// r(t+dt) = r(t) + dt*v(t)+0.5*(dt*dt)*a(t)
			chain_r[i][d] += dt*chain_v[i][d] + 0.5*(dt*dt)*total_force[i][d];
			total_force_old[i][d]= total_force[i][d];
		}
	}

	//----- Calculo de los r_ij^2
	for (int k=0; k<N-1; k++){
		for (int l=1; l<N; l++){
			deltaR2[k][l] = 0;
			for (int d=0;d<3;d++){
				deltaR2[k][l] += (chain_r[k][d]-chain_r[l][d])*(chain_r[k][d]-chain_r[l][d]);
			}
			deltaR2[l][k] = deltaR2[k][l];
		}
		deltaR2[k][k] = 0;
	}
//	for (int w=0;w<N;w++){
//		for (int n=0;n<N;n++){
//			printf("%f\t",deltaR2[w][n]);
//		}
//		printf("\n");
//	}


	double auxvar;
	double auxforce[3];

	//----- Calculo de las fuerzas sobre las particulas.
	for (int m=0;m<N-1;m++){
		auxvar = force_cc_r(epsi,q,deltaR2[m][m+1]);
		for (int d=0;d<3;d++){
			auxforce[d] = auxvar*(chain_r[m][d]-chain_r[m+1][d]);
			total_force[m][d]	+=  auxforce[d];
			total_force[m+1][d]	+= -auxforce[d];
		}

		for(int i=m+2;i<N;i++){
			auxvar = force_hydro(epsi,Ec,deltaR2[m][i],hydro[m]+hydro[i]);
			for (int d=0;d<3;d++){
				auxforce[d] = auxvar*(chain_r[m][d]-chain_r[i][d]);
				total_force[m][d] +=  auxforce[d];
				total_force[i][d] += -auxforce[d];
			}
		}
	}
	//----- Calcular la fuerza nuevamente y velocidades.
	for (int i=0;i<N;i++){
		for(int d=0;d<3;d++){
			// v(t+dt) = v(t) + 0.5*dt*(a(t)+a(t+dt))
			chain_v[i][d] += 0.5*dt*(total_force[i][d]+total_force_old[i][d]);
		}
	}

	//----- Print data.


	fprintf(fp,"ENDMDL\n");

	time++;
  }
	fclose(fp);
	return EXIT_SUCCESS;
}
