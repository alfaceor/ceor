/*
 * util_functions.cpp
 *
 *  Created on: Mar 8, 2012
 *      Author: alfaceor
 */
#include "util_functions.h"

// r2 = r*r
double phi_cc(double epsilon, double q,double  r2){
	if (r2<=1){
		double r6  = 1/(r2*r2*r2);
		double r12 = r6*r6;
		return (epsilon*(r12 -2*r6 +1));
	}else{
		double r = sqrt(r2);
		return (-epsilon*log(1-((r-1)*(r-1))/(q*q)));
	}
}

double phi_att(double epsilon,double E_c,double r2){
	double r6  = 1/(r2*r2*r2);
	double r12 = r6*r6;
	return (-epsilon*E_c*(r12-2*r6));
}

double phi_rep(double epsilon,double r2){
	double r6  = 1/(r2*r2*r2);
	double r12 = r6*r6;
	if (r2<=1){
		return ( epsilon*(r12 - 2*r6 + 1) );
	}else{
		return 0;
	}
}

double potential_hydro(double epsilon,double E_c, double r2, int hydroValue){
	if(hydroValue==-2){
			// fuerza atractiva
			return(phi_att(epsilon,E_c,r2));
		}else{
			// fuerza repulsiva
			return(phi_rep(epsilon,r2));
		}
}

// Calculate quadratic norm
double norm2(double *vector1, double *vector2, int d){
	double norm2=0;
	for (int i=0;i<d;i++){
		norm2 += (vector1[i]-vector2[i])*(vector1[i]-vector2[i]);
	}
	return norm2;
}

// FORCES
// FIXME: Complete function
double force_cc_r(double epsilon, double q, double r2){
	double force;
	double r6 = 1.0/(r2*r2*r2);
	double r8 = r6*(1.0/(r2));
	if(r2<=1){
		force = ( 12*epsilon*r8*(r6-1) );
	}else{
		double r   = sqrt(r2);
		double aux = (r-1)/q;
		double denominador = r*q*(1.0/aux - aux);
		force =  -2*epsilon/denominador;
	}
	return force;
}

double force_att_r(double epsilon, double Ec, double r2){
	double r6 = 1.0/(r2*r2*r2);
	double r8 = r6*(1.0/(r2));

	return (-12*epsilon*Ec*r8*(r6-1));
}

double force_rep_r(double epsilon, double r2){
	double force;
	double r6 = 1.0/(r2*r2*r2);
	double r8 = r6*(1.0/(r2));
	if (r2<1){
		force = 12*epsilon*r8*(r6-1);
	}else{
		force = 0;
	}
	return force;
}

double force_hydro(double epsilon, double Ec, double r2, int hydroValue){
	if(hydroValue==-2){
		// fuerza atractiva
		return(force_att_r(epsilon,Ec,r2));
	}else{
		// fuerza repulsiva
		return(force_rep_r(epsilon,r2));
	}
}



