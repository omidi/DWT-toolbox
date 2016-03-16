/*
 * Alignment.h
 *
 *  Created on: Aug 18, 2010
 *      Author: omidi
 */

#ifndef ALIGNMENT_H_
#define ALIGNMENT_H_

#include "Motif.h"
#include "constants.h"
#include <iostream>
#include <stdlib.h>
#include <boost/algorithm/string.hpp>
#include <math.h>
#include <vector>
#include <fstream>
#include <exception>

namespace std {

class Alignment {
	std::vector <Motif> alignment;
//	Motif * alignment;
	std::vector <std::vector<float> > n_i;	// to hold the frequencies for each column
//	float ** n_i;
	std::vector <std::vector<std::vector<std::vector<float> > > > n_ij;		// to hold the mutual columns frequency
	float n_rows;
	unsigned short int n_cols;
	string TF;
	string filename;
public:
	Alignment(string);
	Alignment();
	virtual ~Alignment();
	vector<Motif> read_alignment(string);
//	Motif* read_alignment(string);
	Motif operator[](unsigned int);
	float nrows();
	unsigned short int ncols();
	float get_n_i(char, unsigned short int);
	float get_n_i(unsigned short int, unsigned short int);	// same version but without convert( function
	float get_n_ij(char, char, unsigned short int, unsigned short int);
	float get_n_ij(unsigned short int, unsigned short int, unsigned short int, unsigned short int);	// same version but without convert( function
	string get_TF();
	void initialize_n_i(unsigned short int);
	void initialize_n_ij(unsigned short int);
};

}

#endif /* ALIGNMENT_H_ */
