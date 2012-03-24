/*
 * Conformation.cpp
 *
 *  Created on: Feb 7, 2012
 *      Author: alfaceor
 */

#include "Conformation.h"

Conformation::Conformation(int N, char *hydroChain, double temp, char *basename)
{
	// Initialize chain
	this->N = N;
	this->chain		= new Monomer[N]();
	this->deltaR2	= new double[N*N]();
	this->basename	= basename;
	this->Energy	= 0.0;
	this->KinecticEnergy =0.0;
	this->PotentialEnergy=0.0;

	FILE *fp;
	fp=fopen(basename,"r");
	//FIXME: temporal flag remove later
//	if(fp==NULL){
	if(true){
		printf("filename= %s doesnt exist!\n",basename);
		printf("creating new file %s with random data...\n",basename);
		fp=fopen(basename,"w");
		//---------- random variables
		randomPositions();
		gaussianRandomVelocities(temp); // temp

		for (int i=0;i<N;i++){
			//-------------------------
			// la misma cadena que usa Lois 2008
			if(hydroChain[i]=='1'){
				this->chain[i].hydro = -1.0;
			}else{
				this->chain[i].hydro =  1.0;
			}
//			if(i==0 || i==4 || i==8 || i==12) this->chain[i].hydro = -1.0;
//			else			this->chain[i].hydro =  1.0;

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

void Conformation::calculateBondPotential(double epsi,double q){
	double auxvar=0.0;
	for (int i=0;i<N-1;i++){
		auxvar	+=	phi_cc(epsi,q,deltaR2[i*N+(i+1)]);
	}
	PotentialEnergy+=auxvar;
}

void Conformation::calculateHydroForces(double epsi,double Ec){
	double auxvar;
	double auxforce[DIM];
	for (int m=0;m<N-2;m++){
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

void Conformation::calculateHydroPotential(double epsi,double Ec){
	double auxvar=0.0;
	for(int m=0; m<N-2; m++){
		for (int i=m+2;i<N;i++){
			auxvar += potential_hydro(epsi,Ec,deltaR2[m*N+i],chain[m].hydro+chain[i].hydro);
		}
	}
	PotentialEnergy+=auxvar;
}

void Conformation::calculateTotalForces(double epsi, double q, double Ec){
	cleanForces();
	calculateDeltaR2();
	calculateBondForces(epsi,q);
	calculateHydroForces(epsi,Ec);
}

void Conformation::calculateTotalEnergy(double epsi,double q, double Ec){
	cleanEnergyValues();
	calculateBondPotential(epsi,q);
	calculateHydroPotential(epsi,Ec);
	calculateKineticEnergy();
	Energy = KinecticEnergy + PotentialEnergy;
}



void Conformation::actualizePositions(double dt){
	for (int i=0;i<N;i++){
		chain[i].actualizeVec_r(dt);
	}
}

void Conformation::actualizeVelocities(double dt){
	for (int i=0;i<N;i++){
		chain[i].actualizeVec_v(dt);
	}
}

void Conformation::calculateKineticEnergy(){
	double kinectic_energy=0.0;
	for (int i=0; i<N; i++){
		for (int d=0; d<DIM; d++){
			kinectic_energy+=chain[i].mass*(chain[i].vec_v[d]*chain[i].vec_v[d]);
		}
	}
	KinecticEnergy=0.5*kinectic_energy;
}

void Conformation::cleanForces(){
	for(int i=0; i<N; i++){
		for (int d=0; d<DIM; d++){
			chain[i].total_force[d]=0.0;
		}
	}
}

void Conformation::cleanEnergyValues(){
	this->Energy = 0.0;
	this->KinecticEnergy = 0.0;
	this->PotentialEnergy= 0.0;
}

void Conformation::randomPositions(){
	gsl_rng *r1;
	const gsl_rng_type * T;
	T = gsl_rng_default;
	r1 = gsl_rng_alloc(T);
	double randomaux;

	for (int i=0; i<N; i++){
		for (int d=0;d<DIM;d++){
			// vector positions and velocities
			chain[i].vec_r[d] = 0.0;
			if (d==0){
				randomaux = 0.001*gsl_rng_uniform_pos(r1);
				chain[i].vec_r[d] = randomaux;
			}
		}
		chain[i].vec_r[0] += chain[i].zigma*i;
	}

}

void Conformation::gaussianRandomVelocities(double temp){
	gsl_rng *r2;
	gsl_rng *r1;
	const gsl_rng_type * T;
	T = gsl_rng_default;
	r2 = gsl_rng_alloc(T);
	r1 = gsl_rng_alloc(T);
	double randomaux;
	double randomang;
//	double temp = 0.04;
	double avg_vel[DIM];

	for (int i=0;i<N;i++){
		randomaux = gsl_ran_gaussian(r2,1);
		randomang = M_PI*gsl_rng_uniform(r1);
		for (int d=0;d<DIM;d++){
			// vector positions and velocities
			chain[i].vec_v[d] = 0.0;
			// vector velocities in any direction.

			if (d==0){
				chain[i].vec_v[d] = randomaux*cos(randomang);
			}else if (d==1){
				chain[i].vec_v[d] = randomaux*sin(randomang);
			}else{
				chain[i].vec_v[d] = 0.0;
			}
			avg_vel[d] += chain[i].vec_v[d];
		}
	}
	// TOTAL MOMENTUM = 0 (ZERO)
	for (int i=0;i<N;i++){
		for (int d=0; d<DIM; d++){
			if(i == 0){ avg_vel[d] = avg_vel[d]/N; }
			chain[i].vec_v[d] -= avg_vel[d];
		}
	}

	calculateKineticEnergy();

	double scale_temp=sqrt(3.0*temp/KinecticEnergy);

	for (int i=0;i<N;i++){
		for (int d=0; d<DIM; d++){
			chain[i].vec_v[d] = scale_temp*chain[i].vec_v[d];
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

void Conformation::print_pdb_conformation(FILE *fp,int time_model){
	char name[]="    ";
	fprintf(fp,"MODEL\t%d\n",time_model);
	for (int i=0;i<N;i++){
		if(chain[i].hydro == 1){	name[0]='Z';name[1]='n';}
		else{				name[0]='N';name[1]=' ';}
		print_pdb_line(fp,i+1,chain[i].vec_r[0],chain[i].vec_r[1],chain[i].vec_r[2],name,chain[i].vec_v[0]);
	}
	fprintf(fp,"ENDMDL\n");
}

void Conformation::print_pdb_line(FILE *fp,int serial, double x, double y, double z,char *name,double tempFactor){
	char recordname[]="HETATM";	// 1 - 6        Record name    "HETATM"
//	int  serial    =1;						// 7 - 11       Integer   Atom serial number.
//	char name    []="    ";		//13 - 16       Atom      Atom name.
	char altLoc  []=" ";				//17            character Alternate location indicator.
	char resName []="   ";			//18 - 20       Residue name  Residue name.
	char chainID []="A";				//22            character     Chain identifier.
	char resSeq  []="1";
	char iCode   []=" ";
//	double x       = 13.872;
//	double y       = -2.555;
//	double z       = -29.045;
	double occupancy	=  1.00;
//	double tempFactor	= 27.36;
	char element []="  ";
	char charge  []="  ";

//	char pdb_line[80];
//  const char atom_line_iformat[]=
//		"%6s%5d %4s%1s%3s %1s%4s%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s";
//
//	sprintf(pdb_line,atom_line_iformat, recordname, serial, name, altLoc, resName, chainID, resSeq, iCode, x, y, z, occupancy, tempFactor, element,charge );
	fprintf(fp,"%6s%5d %4s%1s%3s %1s%4s%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s", recordname, serial, name, altLoc, resName, chainID, resSeq, iCode, x, y, z, occupancy, tempFactor, element,charge );
	fprintf(fp,"\n");
}




//void Conformation::print_pdb_line(int serial){
//	char recordname[]="HETATM";	// 1 - 6        Record name    "HETATM"
////  int  serial    =1;			// 7 - 11       Integer   Atom serial number.
//	char name    []="    ";		//13 - 16       Atom      Atom name.
//	char altLoc  []=" ";		//17            character Alternate location indicator.
//	char resName []=" MG";		//18 - 20       Residue name  Residue name.
//	char chainID []="A";		//22            character     Chain identifier.
//	char resSeq  []="1";
//	char iCode   []=" ";
//	double x       = 11.3;	//chain[serial].vec_r[0];
//	double y       = 10.1;	//chain[serial].vec_r[1];
//	double z       = 0.0;	//chain[serial].vec_r[2];
//	double occupancy =1.00;		// FIXME: WHY this value?
//	double tempFactor =27.36;	// FIXME:
//	char element []="  ";
//	char charge  []="  ";
//
//  char pdb_line[80];
//  const char atom_line_iformat[]=
//    "%6s%5d %4s%1s%3s %1s%4s%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s";
//  sprintf(pdb_line,atom_line_iformat, recordname, serial,name, altLoc, resName, chainID, resSeq, iCode,x, y, z, occupancy, tempFactor, element,charge);
//  printf("%s\n",pdb_line);
//}




