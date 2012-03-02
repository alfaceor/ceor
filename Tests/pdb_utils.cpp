/*
 * pdb_utils.cpp
 *
 *  Created on: Feb 16, 2012
 *      Author: alfaceor
 */

/**
Description: write a pdb file for a HETATM
author: Carlos Olivares
e-mail: alfaceor@gmail.com
license: creative commons
Reference: http://www.wwpdb.org/documentation/format33/sect9.html#HETATM

COLUMNS       DATA  TYPE     FIELD         DEFINITION
-----------------------------------------------------------------------
 1 - 6        Record name    "HETATM"
 7 - 11       Integer        serial        Atom serial number.
12       --- SPACE
13 - 16       Atom           name          Atom name.
17            Character      altLoc        Alternate location indicator.
18 - 20       Residue name   resName       Residue name.
21       --- SPACE
22            Character      chainID       Chain identifier.
23 - 26       Integer        resSeq        Residue sequence number.
27            AChar          iCode         Code for insertion of residues.
28 - 30  --- SPACE
31 - 38       Real(8.3)      x             Orthogonal coordinates for X.
39 - 46       Real(8.3)      y             Orthogonal coordinates for Y.
47 - 54       Real(8.3)      z             Orthogonal coordinates for Z.
55 - 60       Real(6.2)      occupancy     Occupancy.
61 - 66       Real(6.2)      tempFactor    Temperature factor.
67 - 76  --- SPACE
77 - 78       LString(2)     element       Element symbol; right-justified.
79 - 80       LString(2)     charge        Charge on the atom.


EXAMPLE

         1         2         3         4         5         6         7         8
12345678901234567890123456789012345678901234567890123456789012345678901234567890
HETATM 8237 MG    MG A1001      13.872  -2.555 -29.045  1.00 27.36          MG

HETATM 3835 FE   HEM A   1      17.140   3.115  15.066  1.00 14.14          FE
HETATM 8238  S   SO4 A2001      10.885 -15.746 -14.404  1.00 47.84           S
HETATM 8239  O1  SO4 A2001      11.191 -14.833 -15.531  1.00 50.12           O
HETATM 8240  O2  SO4 A2001       9.576 -16.338 -14.706  1.00 48.55           O
HETATM 8241  O3  SO4 A2001      11.995 -16.703 -14.431  1.00 49.88           O
HETATM 8242  O4  SO4 A2001      10.932 -15.073 -13.100  1.00 49.91           O

*/

#include "pdb_utils.h"

void print_pdb_line(FILE *fp,int serial, double x, double y, double z,char *name,double tempFactor){
	char recordname[]="HETATM";	// 1 - 6        Record name    "HETATM"
//	int  serial    =1;						// 7 - 11       Integer   Atom serial number.
//	char name    []="    ";		//13 - 16       Atom      Atom name.
	char altLoc  []=" ";				//17            character Alternate location indicator.
	char resName []="   ";			//18 - 20       Residue name  Residue name.
	char chainID []="A";				//22            character     Chain identifier.
	char resSeq  []="1";
	char iCode   []=" ";
//	double x       = 13.872;
//	double y       = -2.555;
//	double z       = -29.045;
	double occupancy	=  1.00;
//	double tempFactor	= 27.36;
	char element []="  ";
	char charge  []="  ";

//	char pdb_line[80];
//  const char atom_line_iformat[]=
//		"%6s%5d %4s%1s%3s %1s%4s%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s";
//
//	sprintf(pdb_line,atom_line_iformat, recordname, serial, name, altLoc, resName, chainID, resSeq, iCode, x, y, z, occupancy, tempFactor, element,charge );
	fprintf(fp,"%6s%5d %4s%1s%3s %1s%4s%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s", recordname, serial, name, altLoc, resName, chainID, resSeq, iCode, x, y, z, occupancy, tempFactor, element,charge );
	fprintf(fp,"\n");
}

