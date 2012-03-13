/*
 * util_functions.h
 *
 *  Created on: Mar 8, 2012
 *      Author: alfaceor
 */

#ifndef UTIL_FUNCTIONS_H_
#define UTIL_FUNCTIONS_H_

#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#define DIM 3

double phi_cc(double epsilon, double q,double  r2);
double phi_att(double epsilon,double E_c,double r2);
double phi_rep(double epsilon,double r2);
double norm2(double *vector1, double *vector2, int d);
double force_cc_r(double epsilon, double q, double r2);
double force_att_r(double epsilon, double Ec, double r2);
double force_rep_r(double epsilon, double r2);
double force_hydro(double epsilon, double Ec, double r2, int hydroValue);



#endif /* UTIL_FUNCTIONS_H_ */
