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


void InitialValues(double *total_force,double *chain_r, double *chain_v, double sigma,const int N, const int DIM){
	gsl_rng *r1,*r2;
	const gsl_rng_type * T;
	T = gsl_rng_default;
	r1 = gsl_rng_alloc(T);
	r2 = gsl_rng_alloc(T);
	double randomaux;
	for (int i=0;i<N;i++){
		for (int d=0;d<3;d++){
			// vector force
			total_force[i*DIM+d]=0;
			// vector positions and velocities
			chain_r[i*DIM+d] = 0.0;
			chain_v[i*DIM+d] = 0.0;
			if (d==0){
				randomaux = 0.001*gsl_rng_uniform_pos(r1);
				chain_r[i*DIM+d] = randomaux;
				// vector velocities
				randomaux = 0.5*(2*gsl_ran_gaussian(r2,sigma) - 1);
				chain_v[i*DIM+d] = randomaux;
			}
		}
		chain_r[i*DIM+0] += sigma*i;
//		printf("%f\t%f\n",chain_r[i*N+0],chain_v[i*N+0]);
	}

}

void CalculateR2(double *deltaR2,double *chain_r, const int N, const int DIM){
//		printf("----------- r**2 -----------\n");
	double auxdelta;
	for (int k=0; k<N; k++){
		for (int l=0; l<N; l++){
			deltaR2[k*N+l] = 0.0;
			for (int d=0;d<DIM;d++){
				auxdelta = (chain_r[k*DIM+d]-chain_r[l*DIM+d]);
				deltaR2[k*N+l] += auxdelta*auxdelta;
//				printf("%d=%f\t",d,auxdelta);
			}
//			printf("-> %d%d\t",k,l);
			deltaR2[l*N+k] = deltaR2[k*N+l];
		}
//		printf("%f\n",deltaR2[k][0]);
//		printf("\n");
//		deltaR2[k][k] = 0.0;
	}
//	printf("*******************\n");

}

int main(int argc, char* argv[]) {
	const int N = 9; // = atoi(argv[1]);
	const int DIM = 3;

	// phi_cc function
//	double dx=0.025;
//	double x;
	int time=0;
	int total_time= 112;//atoi(argv[1]);
	double dt=0.001;
	double total_force[N][DIM];
	double total_force_old[N][DIM];
	double chain_r[N][DIM];		// monomer positions
	double chain_v[N][DIM];		// monomer velocities

	int    hydro[N];	// = {1,1,-1,1,-1,-1,1,1,1,-1,1,1,1};
	double deltaR2[N][N];
	double epsi=1.0, Ec=-0.7, q=0.2;
	double sigma = 1.0;	// diametro


	//----- Initial positions and velocities
	InitialValues(total_force[0], chain_r[0], chain_v[0],sigma, N, DIM);

	// FIXME: Nose pq en el bucle no puede inicializar la variable hydro
	hydro[0]  =  1;
	hydro[1]  =  1;
	hydro[2]  = -1;
	hydro[3]  = -1;
	hydro[4]  =  1;
	hydro[5]  = -1;
	hydro[6]  =  1;
	hydro[7]  = -1;
	hydro[8]  =  1;
//	hydro[9]  = -1;
//	hydro[10] =  1;
//	hydro[11] =  1;
//	hydro[12] =  1;

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

	//----- Calculo de los r_ij^2
	CalculateR2(deltaR2[0],chain_r[0], N, DIM);

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

//		FIXME: las fuerzas hidrophobicas o hidrofilicas
//		for(int i=m+2;i<N;i++){
//			auxvar = force_hydro(epsi,Ec,deltaR2[m][i],hydro[m]+hydro[i]);
//			for (int d=0;d<3;d++){
//				auxforce[d] = auxvar*(chain_r[m][d]-chain_r[i][d]);
//				total_force[m][d] +=  auxforce[d];
//				total_force[i][d] += -auxforce[d];
//			}
//		}
	}

	//----- Actualizacion de las posiciones
	for (int i=0;i<N;i++){
		for (int d=0;d<3;d++){
			// r(t+dt) = r(t) + dt*v(t)+0.5*(dt*dt)*a(t)
			chain_r[i][d] += dt*chain_v[i][d] + 0.5*(dt*dt)*total_force[i][d];
			chain_v[i][d] += 0.5*dt*(total_force[i][d]+total_force_old[i][d]);
			total_force_old[i][d]= total_force[i][d];
		}
	}

	//----- Print data.
	fprintf(fp,"ENDMDL\n");

	time++;
  }
	fclose(fp);
	return EXIT_SUCCESS;
}
