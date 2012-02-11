/*
 * Conformation.cpp
 *
 *  Created on: Feb 7, 2012
 *      Author: alfaceor
 */

#include "Conformation.h"

Conformation::Conformation() {
	int d=3;
	for(int i=0;i<N;i++){
		for(int j=0;j<d;j++){
			if(j==0) chain[i].vec_r[j] = i*chain[i].zigma;
			else	 chain[i].vec_r[j] = 0;
		}
	}
	// TODO: initialize force matrix
}

Conformation::~Conformation() {
	// TODO Auto-generated destructor stub
}

void Conformation::GenDeltas(){
	// Generate a top diagonal matrix x_ij
	for (int i=0;i<N-1;i++){
		for (int j=i+1; j<N;j++){
			deltaX[i][j]	= 	chain[i].vec_r[0] - chain[j].vec_r[0];
			deltaY[i][j]	= 	chain[i].vec_r[1] - chain[j].vec_r[1];
			deltaZ[i][j]	= 	chain[i].vec_r[2] - chain[j].vec_r[2];
			// element by element matrix product.
			deltaR2[i][j]	=	deltaX[i][j]*deltaX[i][j]+
								deltaY[i][j]*deltaY[i][j]+
								deltaZ[i][j]*deltaZ[i][j];
			deltaR[i][j]	=	sqrt(deltaR2[i][j]);
		}
	}

}

double Conformation::GaussianForce(double mean, double var){

		const gsl_rng_type * T;
		gsl_rng * r;
		gsl_rng_env_setup();
		T = gsl_rng_default;
		r = gsl_rng_alloc(T);
		// FIXME: complete this function to generate random forces.
		double randomNumber;
		for (int i=0;i<10;i++){
			//randomNumber = (rand() % 10) / 10.0;

			randomNumber = gsl_ran_gaussian(r,1);
			printf("%f\n",randomNumber);
		}
		return randomNumber;

}
// TODO: Calculate the force over the particle k.
// So the force with the interaction beteween k-1 and k plus k and k+1
void Conformation::Force_cc(double epsilon, double q){

	for (int k=1;k<N-1;k++){
		if(deltaR[k-1][k]<=1){
		// XXX: k-1 and k
			double deltaR6 = pow(deltaR2[k-1][k],3);
			double deltaR8 = deltaR6*deltaR2[k-1][k];
			// Acumulate axis force projection.
			for (int x=0;x<3;x++)	// FIXME: IS NOT
				force[k][x] += -12*epsilon*(1/deltaR6 - 1)*deltaX[k-1][k]/deltaR8;
		}else{
			// TODO: Logarithmic force

			force[k][x] += 2*epsilon*(deltaX[k-1][k]/deltaR[k-1][k])*(1/((q/deltaR[k-1][k])-(deltaR[k-1][k]/q)));
		}

		// XXX: k and k+1
		if(deltaR[k-1][k]<=1){
		// k-1 and k
			double deltaR6 = pow(deltaR2[k-1][k],3);
			double deltaR8 = deltaR6*deltaR2[k-1][k];
			// Acumulate axis force projection.
			for (int x=0;x<3;x++)	// FIXME: IS NOT
				force[k][x] += -12*epsilon*(1/deltaR6 - 1)*deltaX[k-1][k]/deltaR8;
		}else{
			// TODO: Logarithmic force

			force[k][x] += 2*epsilon*(deltaX[k-1][k]/deltaR[k-1][k])*(1/((q/deltaR[k-1][k])-(deltaR[k-1][k]/q)));
		}


		// epsilon*(1/pow(deltaR2[k][k+1],3) - 1)*deltaX[k][k+1]/pow(deltaR2[k][k+1],4);
	}

}


double Conformation::Force_cc2(int i,int j, double r_ij, double epsilon){
	double force;
	double r_inv6,r_inv8;
	// calculate r_ij e seus exponenetes
	r_inv6=1.0/pow(r_ij,6);
	r_inv8=r_inv6/pow(r_ij,2);
	double x_i=0.3, x_j=0.5;

	if(r_ij<=1){
		printf("r_ij <= 1");
		force = 12.0*epsilon*(r_inv6-1)*r_inv8*(x_i - x_j);
		return force;
	}else{
		printf("r_ij > 1");
		force = 0;
		return force;
	}
}

void Conformation::nextStep(){
	// TODO make the evolution of the systems
	// calculate the force
	// update the positions and velocities
	//
}

void Conformation::print_r(){
	for(int i=0;i<N;i++){
		chain[i].print_r(); printf("\n");
	}
}

//void Conformation::print_pdbfile(char *filename)
//{
//	FILE *fp;
//	int rating = 9;
//	if (fp = fopen(filename, "w"))
//	{
//		fprintf(fp, "WWW\n");
//		fprintf(fp, "Topic: computer programming\n");
//		fprintf(fp, "Rating out of 10 : %d \n",rating );
//		fclose(fp);
//	}
//	else
//	printf("Error opening d:/website.txt\n");
//
//}

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




