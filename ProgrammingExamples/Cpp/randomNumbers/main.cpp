#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

void generate2randomNumbers(int N, double sigma){
	gsl_rng *r1, *r2;
	const gsl_rng_type * T1, * T2;
	T1	= gsl_rng_default;
	T2	= gsl_rng_default;
	r1	= gsl_rng_alloc(T1);
	r2	= gsl_rng_alloc(T2);
	gsl_rng_set(r1,10);
	gsl_rng_set(r2,4);
	double randomaux1,randomaux2;
	for (int i=0; i<N; i++){
		randomaux1	= gsl_ran_gaussian(r1,sigma);
		randomaux2	= gsl_ran_gaussian(r2,sigma);
		printf("%f\t%f\n",randomaux1,randomaux2);
	}
	gsl_rng_free(r1);
	gsl_rng_free(r2);
}

void generateRandomNumber(gsl_rng *r, double sigma){
	double randomaux;
	randomaux	= gsl_ran_gaussian(r,sigma);
	printf("%f\n",randomaux);
}

int main(int argc, char* argv[]) {
	if (argv[1] == NULL || argv[2] == NULL){
		printf("Usage: ./main <N> <sigma>\n");
		return EXIT_FAILURE;
	}
	int N=atoi(argv[1]);
	double sigma=atof(argv[2]);
	generate2randomNumbers(N,sigma);

	const gsl_rng_type * T1;
	gsl_rng *r1;
//	gsl_rng_env_setup();
	unsigned long int seed;
//	seed = time (NULL) * getpid();   // FIXME: SEED 
	seed=10;
	
	T1	= gsl_rng_default;
	r1	= gsl_rng_alloc(T1);
	gsl_rng_set(r1, seed);                  // set seed
	printf("generator type: %s\n", gsl_rng_name(r1));
	printf("seed = %lu\n",gsl_rng_default_seed);
	printf ("first value = %lu\n", gsl_rng_get (r1));

	for (int i=0; i<N; i++){
		generateRandomNumber(r1,sigma);
	}

	printf("------------\n");

	gsl_rng_free(r1);
}
