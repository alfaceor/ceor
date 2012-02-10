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
	Conformation polimero;
	polimero.nextStep();
	printf("positions polimero\n");
	polimero.print_r();
	return EXIT_SUCCESS;
}
