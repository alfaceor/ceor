/*
 * forces.h
 *
 *  Created on: Feb 16, 2012
 *      Author: alfaceor
 */

#ifndef FORCES_H_
#define FORCES_H_

#include <gsl/gsl_math.h>

double phi_cc(double epsilon, double q,double  r2);
double phi_att(double epsilon,double E_c,double r2);
double phi_rep(double epsilon,double r2);
double norm2(double *vector1, double *vector2, int d);
double force_cc_r(double epsilon, double q, double r2);
double force_att_r(double epsilon, double Ec, double r2);
double force_rep_r(double epsilon, double r2);
double force_hydro(double epsilon, double Ec, double r2, int hydroValue);


#endif /* FORCES_H_ */
