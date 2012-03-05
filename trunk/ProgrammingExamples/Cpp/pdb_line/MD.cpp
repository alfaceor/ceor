#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>

#include "Conformation.h"
/*
void print_pdb_line(int serial){
	char recordname[]="HETATM";	// 1 - 6        Record name    "HETATM"
//  int  serial    =1;			// 7 - 11       Integer   Atom serial number.
	char name    []="    ";		//13 - 16       Atom      Atom name.
	char altLoc  []=" ";		//17            character Alternate location indicator.
	char resName []=" MG";		//18 - 20       Residue name  Residue name.
	char chainID []="A";		//22            character     Chain identifier.
	char resSeq  []="1";
	char iCode   []=" ";
	double x       = 11.3;	//chain[serial].vec_r[0];
	double y       = 10.1;	//chain[serial].vec_r[1];
	double z       = 0.0;	//chain[serial].vec_r[2];
	double occupancy =1.00;		// FIXME: WHY this value?
	double tempFactor =27.36;	// FIXME:
	char element []="  ";
	char charge  []="  ";

  char pdb_line[80];
  const char atom_line_iformat[]=
    "%6s%5d %4s%1s%3s %1s%4s%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s";
  sprintf(pdb_line,atom_line_iformat, recordname, serial,name, altLoc, resName, chainID, resSeq, iCode,x, y, z, occupancy, tempFactor, element,charge);
  printf("%s\n",pdb_line);
}
*/

int main(void) {
//	char resname[]="  N";
	Conformation polimero;
//	polimero.nextStep();
//	printf("positions polimero\n");
//	//polimero.print_r();
//	resname[2]='O';
//	printf("resnmae =%s",resname);
//	polimero.print_pdb_line(1);
	/**i******************/
 polimero.print_pdb_line(1);
	/********************/
	return EXIT_SUCCESS;
}
