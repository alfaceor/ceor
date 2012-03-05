
#include <stdio.h>
#include <stdlib.h>


int main(int argc, char* argv[]) {
	FILE *fp;
	fp = fopen("data.dat","r");
	double number;
	printf("read file\n");
	for (int i=0;i<6;i++){
		fscanf(fp,"%lf",&number);
		printf("number = %f\n",number);
	}
	fclose(fp);
	return EXIT_SUCCESS;
}

