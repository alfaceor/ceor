#include <stdlib.h>
#include <stdio.h>
#include <string.h>

void filenames(){

}

int main(int argc, char* argv[]){

	const int M	=	strlen(argv[1]);
	double epsi	=	 atof(argv[2]);
	int total_time = atoi(argv[3]);
	double q	=	 0.1;
	double Ec	=	-1.0;
	double temp	=	0.04;
	double dt	=	0.001;
	char ext_pdb[]=".pdb";
	char ext_dat[]=".dat";
	char filename_pattern[40];
	char filename_pdb[50];
	char filename_dat[50];

	strcpy(filename_pattern,"md-N");
	for (int i=1;i<argc;i++){
		strcat(filename_pattern,"_");
		strcat(filename_pattern,argv[i]);
	}

	strcat(filename_pdb,filename_pattern);
	strcat(filename_dat,filename_pattern);
	strcat(filename_pdb,ext_pdb);
	strcat(filename_dat,ext_dat);

	// char *pdbfile; strcpy(pdbfile,"protein.pdb");
	printf("fuck you!\n");	
	printf("%s\t%d\n",filename_pattern,M);

	FILE *fp;
	fp = fopen(filename_pattern,"w");
	fprintf(fp,"joder N=%d\n",M);
	fprintf(fp,"joder\n");
	fprintf(fp,"joder\n");
	fclose(fp);

}
