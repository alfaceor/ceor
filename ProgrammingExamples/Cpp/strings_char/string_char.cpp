#include <stdlib.h>
#include <stdio.h>
int main(int argc, char* argv[]){

	const int M	=	atoi(argv[1]);
	double epsi	=	 1.0;
	double q	=	 0.1;
	double Ec	=	-1.0;
	double temp	=	0.04;
	double dt	=	0.001;
	int total_time = 200000;

	char *filename;	filename="data_inicial.dat";
	char *pdbfile;	pdbfile	="protein.pdb";

	FILE *fp;
	fp = fopen(pdbfile,"w");
	fprintf(fp,"joder");
	fclose(fp);

}
