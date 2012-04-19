#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>


int main(int argc, char* argv[]) {
	gsl_rng *r1, *r2;
	
	const gsl_rng_type * T1, * T2;
	
	T1	= gsl_rng_default;
	T2	= gsl_rng_default;
	r1	= gsl_rng_alloc(T1);
	r2	= gsl_rng_alloc(T2);
	
	double randomaux1,randomaux2;
	int N=atoi(argv[1]);
	double sigma=atof(argv[2]);
	for (int i=0; i<N; i++){
		randomaux1	= gsl_ran_gaussian(r1,sigma);
		randomaux2	= gsl_ran_gaussian(r2,sigma);
		printf("%f\t%f\n",randomaux1,randomaux2);
	}
	
}
