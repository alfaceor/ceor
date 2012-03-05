#include <stdio.h>

/*************************************/

/***
// DON'T WORK
void func(int x, int y, int a[x][y]){
	for (int i=0;i<x;i++){
		for (int j=0;j<y;j++){
			a[i][j] = i+j;
		}
	}
}
***/

void func02(int *ptr_A,int N){
	for (int i=0;i<N;i++){
		for (int j=0;j<N;j++){
			printf("%d\t",ptr_A[i*N+j]);
		}
		printf("\n");
	}	
}

void func03(int *ptr_A,int N){
	for (int i=0;i<N;i++){
		for (int j=0;j<N;j++){
			ptr_A[i*N+j] = i-j;
		}
	}	
}

int main(){
	printf("hello world\n");
	const int N=5;
	int array_var[N][N];
	int *ptr_var;
	
	for (int i=0;i<N;i++){
		for (int j=0;j<N;j++){
			array_var[i][j] = i*j;
		}
	}
	
	// ptr_var=array_var[0];
	func02(array_var[0],N);
	func03(array_var[0],N);
	printf("\n*****************\n");
	func02(array_var[0],N);
	printf("\n++++++++++++\n");
	for (int i=0;i<N;i++){
		for (int j=0;j<N;j++){
			printf("%d\t",array_var[i][j]);
		}
		printf("\n");
	}
	
}
