
#include <stdio.h>
#include <stdlib.h>

void print_array(double *array,const int N){
	for (int i=0;i<N;i++){
		printf("%f\n",*(array+i));
	}
}


int main(int argc, char* argv[]) {
	FILE *fp;
	fp = fopen("data.dat","r");

	double number;
	double arraynum[6];
	printf("read file\n");
	
	for (int i=0;i<6;i++){
		// fscanf(fp,"%lf",&number);
		fscanf(fp,"%lf",&arraynum[i]);
//		printf("number = %f\n",arraynum[i]);
	}
	
	print_array(arraynum,6);
	fclose(fp);
	return EXIT_SUCCESS;
}

