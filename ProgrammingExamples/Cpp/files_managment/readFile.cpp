
#include <stdio.h>
#include <stdlib.h>


int main(int argc, char* argv[]) {
	FILE *fp;
	fp = fopen("data.dat","r");
	int number;
	printf("read file\n");
	for (int i=0;i<6;i++){
		fscanf(fp,"%d",&number);
		printf("number = %d\n",number);
	}
	fclose(fp);
	return EXIT_SUCCESS;
}

