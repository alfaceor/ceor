
#include <stdio.h>
#include <stdlib.h>

void print_array2d(double *array,const int N,const int M){
	for (int i=0;i<N;i++){
		for (int j=0;j<M;j++){
			printf("%f\t",*(array+i*M+j));
		}
		printf("\n");
	}
}


void print_array(double *array,const int N){
	for (int i=0;i<N;i++){
		printf("%f\n",*(array+i));
	}
}


int main(int argc, char* argv[]) {
	FILE *fp;
	char *filename;
	filename = "data.dat";
	fp = fopen(filename,"r");

	double number;
	double arraynum[6];
	double arraynum2d[6][2];
	for (int i=0;i<6;i++){
		for(int j=0;j<2;j++){
			arraynum2d[i][j]=10.0*i+j;
		}
	}

	printf("read file\n");
	
	for (int i=0;i<6;i++){
		for (int j=0;j<2;j++)
		// fscanf(fp,"%lf",&number);
			fscanf(fp,"%lf",&arraynum2d[i][j]);
//		printf("number = %f\n",arraynum[i]);
	}
	
	print_array2d(arraynum2d[0],6,2);
	fclose(fp);
	return EXIT_SUCCESS;
}

