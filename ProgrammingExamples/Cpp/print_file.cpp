#include <stdio.h> 


void print_pdb(char *filename){
FILE *fp;
int rating = 9;
//char *filename;
//filename = "carlos.txt";
if (fp = fopen(filename, "w"))
{
fprintf(fp, "WWW\n");
fprintf(fp, "Topic: computer programming\n");
fprintf(fp, "Rating out of 10 : %d \n",rating );
fclose(fp);
}
else
printf("Error opening d:/website.txt\n");


}

int main()
{
	print_pdb("joder");
	return 0;
}

