/*
 * pdb_utils.h
 *
 *  Created on: Feb 16, 2012
 *      Author: alfaceor
 */

#ifndef PDB_UTILS_H_
#define PDB_UTILS_H_

#include <stdio.h>


void print_pdb_line(FILE *fp,int serial, double x, double y, double z,char *resName);

#endif /* PDB_UTILS_H_ */
