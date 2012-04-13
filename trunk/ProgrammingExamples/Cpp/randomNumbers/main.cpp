#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>


int main(int argc, char* argv[]) {
	gsl_rng *r1;
	
	const gsl_rng_type * T;
	
	T		= gsl_rng_default;
	r1	= gsl_rng_alloc(T);
	
	double randomaux;
	int N=atoi(argv[1]);
	double sigma=atof(argv[2]);
	for (int i=0; i<N; i++){
		randomaux	= gsl_ran_gaussian(r1,sigma);
		printf("%f\n",randomaux);
	}
	
}
