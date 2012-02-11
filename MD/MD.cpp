//============================================================================
// Name        : MD.cpp
// Author      : Carlos Olivares
// Version     :
// Copyright   : GPL
// Description : Hello World in C, Ansi-style
//============================================================================

#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>

#include "Conformation.h"

int main(void) {
	Conformation protein;
	/**
	 * TODO:
	 * -Calculate the distance matrix. deltaR a triangular matrix.
	 * -Calculate the force (aceleration) to change the velocity
	 * -Then change the position.
	 * */
	// Calculate
	protein.nextStep();
	printf("protein coordenates\n");
	protein.print_r();

	/********************/

	return EXIT_SUCCESS;
}
