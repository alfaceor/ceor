/*
 * Conformation.cpp
 *
 *  Created on: Feb 7, 2012
 *      Author: alfaceor
 */

#include "Conformation.h"

Conformation::Conformation(int N, char *basename)
{
	// Initialize chain
	this->N = N;
	this->chain = new Monomer[N]();
	this->deltaR2 = new double[N*N]();
	this->basename = basename;

	FILE *fp;
	fp=fopen(basename,"r");
	if(fp==NULL){
		printf("filename= %s doesnt exist!\n",basename);
		printf("creating new file %s with random data...\n",basename);
		fp=fopen(basename,"w");

		//---------- random variables
		gsl_rng *r1,*r2;
		const gsl_rng_type * T;
		T = gsl_rng_default;
		r1 = gsl_rng_alloc(T);
		r2 = gsl_rng_alloc(T);
		double randomaux;
		for (int i=0;i<N;i++){
			for (int d=0;d<DIM;d++){
				// vector force
				this->chain[i].total_force[d]=0.0;
				// vector positions and velocities
				this->chain[i].vec_r[d] = 0.0;
				this->chain[i].vec_v[d] = 0.0;
				if (d==0){
					randomaux = 0.001*gsl_rng_uniform_pos(r1);
					this->chain[i].vec_r[d] = randomaux;
				}
				// vector velocities in any direction.
				randomaux = 0.5*(2*gsl_ran_gaussian(r2,this->chain[i].zigma) - 1);
				this->chain[i].vec_v[d] = randomaux;
			}
			//-------------------------
			// la misma cadena que usa Lois 2008
			if(i==0 || i==4 || i==8 || i==12) this->chain[i].hydro = -1.0;
			else			this->chain[i].hydro =  1.0;
			this->chain[i].vec_r[0] += this->chain[i].zigma*i;

			// save data in a file.
			fprintf(fp,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t\n",
					this->chain[i].vec_r[0], this->chain[i].vec_r[1], this->chain[i].vec_r[2],
					this->chain[i].vec_v[0], this->chain[i].vec_v[1], this->chain[i].vec_v[2],
					this->chain[i].hydro);
		}

	}else{
		// read a file input data
		for (int i=0; i<N; i++){
			for (int d=0; d<DIM; d++){
				fscanf(fp,"%lf",&this->chain[i].vec_r[d]);
			}
			for (int d=0; d<DIM; d++){
				fscanf(fp,"%lf",&this->chain[i].vec_v[d]);
			}
			fscanf(fp,"%lf",&this->chain[i].hydro);
		}
	}

	fclose(fp);
}

Conformation::~Conformation() {
	// TODO Auto-generated destructor stub
}

void Conformation::calculateDeltaR2(){
	double auxdelta;
	for (int k=0; k<N; k++){
		for (int l=0; l<N; l++){
			deltaR2[k*N+l] = 0.0;
			for (int d=0;d<DIM;d++){
				auxdelta = (chain[k].vec_r[d] - chain[l].vec_r[d]);
				deltaR2[k*N+l] += auxdelta*auxdelta;
			}
			deltaR2[l*N+k] = deltaR2[k*N+l];
		}
	}
}

void Conformation::calculateBondForces(double epsi, double q){
	double auxvar;
	double auxforce[DIM];
	for (int m=0;m<N-1;m++){
		auxvar = force_cc_r(epsi,q,deltaR2[m*N+(m+1)]);
		for (int d=0;d<DIM;d++){
			auxforce[d] = auxvar*(chain[m].vec_r[d] - chain[m+1].vec_r[d]);
			chain[m].total_force[d]		+=  auxforce[d];
			chain[m+1].total_force[d]	+= -auxforce[d];
		}
	}
}

void Conformation::calculateHydroForces(double epsi,double Ec){
	double auxvar;
	double auxforce[DIM];
	for (int m=0;m<N-1;m++){
//		FIXME: las fuerzas hidrophobicas o hidrofilicas
		for (int i=m+2;i<N;i++){
			auxvar = force_hydro(epsi,Ec,deltaR2[m*N+i],chain[m].hydro+chain[i].hydro);
			for (int d=0;d<3;d++){
				auxforce[d] = auxvar*(chain[m].vec_r[d]-chain[i].vec_r[d]);
				chain[m].total_force[d] +=  auxforce[d];
				chain[i].total_force[d] += -auxforce[d];
			}
		}
	}
}

void Conformation::calculateTotalForces(double epsi, double q, double Ec){
	cleanForces();
	calculateDeltaR2();
	calculateBondForces(epsi,q);
	calculateHydroForces(epsi,Ec);
}

void Conformation::integratorVerlet(double dt){
	for (int i=0;i<N;i++){
		chain[i].actualizeVec_v(dt);
	}
}

void Conformation::cleanForces(){
	for(int i=0; i<N; i++){
		for (int d=0; d<DIM; d++){
			chain[i].total_force[d]=0.0;
		}
	}

}

//double Conformation::GaussianForce(double mean, double var){
//
//		const gsl_rng_type * T;
//		gsl_rng * r;
//		gsl_rng_env_setup();
//		T = gsl_rng_default;
//		r = gsl_rng_alloc(T);
//		// FIXME: complete this function to generate random forces.
//		double randomNumber;
//		for (int i=0;i<10;i++){
//			//randomNumber = (rand() % 10) / 10.0;
//
//			randomNumber = gsl_ran_gaussian(r,1);
//			printf("%f\n",randomNumber);
//		}
//		return randomNumber;
//
//}

void Conformation::print_r(){
	for(int i=0;i<N;i++){
		this->chain[i].print_r(); printf("\n");
	}
}


void Conformation::print_pdb_line(int serial){
	char recordname[]="HETATM";	// 1 - 6        Record name    "HETATM"
//  int  serial    =1;			// 7 - 11       Integer   Atom serial number.
	char name    []="    ";		//13 - 16       Atom      Atom name.
	char altLoc  []=" ";		//17            character Alternate location indicator.
	char resName []=" MG";		//18 - 20       Residue name  Residue name.
	char chainID []="A";		//22            character     Chain identifier.
	char resSeq  []="1";
	char iCode   []=" ";
	double x       = 11.3;	//chain[serial].vec_r[0];
	double y       = 10.1;	//chain[serial].vec_r[1];
	double z       = 0.0;	//chain[serial].vec_r[2];
	double occupancy =1.00;		// FIXME: WHY this value?
	double tempFactor =27.36;	// FIXME:
	char element []="  ";
	char charge  []="  ";

  char pdb_line[80];
  const char atom_line_iformat[]=
    "%6s%5d %4s%1s%3s %1s%4s%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s";
  sprintf(pdb_line,atom_line_iformat, recordname, serial,name, altLoc, resName, chainID, resSeq, iCode,x, y, z, occupancy, tempFactor, element,charge);
  printf("%s\n",pdb_line);
}




