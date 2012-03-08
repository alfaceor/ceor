//============================================================================
// Name        : MD.cpp
// Author      : Carlos Olivares
// Version     : 0.5
// Copyright   : GPL
// Description : Molecular dynamics simulation
//============================================================================

#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include "Conformation.h"

int main(int argc, char* argv[]) {
	// input parameters
	char *filename;
	filename="data_inicial.dat";
	const int M	=	13;
	double epsi	=	 1.0;
	double q	=	 0.1;
	double Ec	=	-0.8;

	// Simulation
	Conformation protein(M,filename);
	protein.calculateTotalForces(epsi,q,Ec);

	return EXIT_SUCCESS;
}
