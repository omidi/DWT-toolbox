/*
 * ParameterFile.h
 *
 *  Created on: Feb 28, 2011
 *      Author: omidi
 */

#ifndef PARAMETERFILE_H_
#define PARAMETERFILE_H_


#include <iostream>
#include <vector>
#include <iomanip>
#include <string>
#include <fstream>
#include <stdlib.h>
#include "constants.h"

using namespace std;


class ParameterFile {
	string file_name;
	string resulting_file_name;	// the name of the file in which the results of sliding window will be written in.
	std::vector<float> background;
	bool background_is_set;
	bool rho_is_set;
	double precision;	// by default is 5.0
	double lowest_score;	// the minimum score for printing out or saving into the file, in LOG space
	double rho; 	// rho is a control parameter that lead to flexibility of dependency structures, with rho=1 the algorithm tries to fit spanning trees, and
					// with rho~0, algorithm work in no-dependency-at-all fashion. One should find the rho in which maximizes the likelihood function.
	double lambda;
public:
	ParameterFile();
	ParameterFile(string);
	virtual ~ParameterFile();
	std::vector<float> give_background();
	double give_precision();
	bool is_background();
	bool is_rho();
	bool is_resulting_file_name();
	double give_lowest_score();
	double give_rho();
	string give_resulting_file_name();
	double give_lambda();
	void initialize_background (ifstream&);
	void initialize_minimum_score (ifstream&);
	void initialize_resulting_file_name (ifstream&);
	void initialize_rho(ifstream&);
	void initialize_precision(ifstream&);
	void initialize_lambda(ifstream&);
	void change_minimum_score (double);
	void change_precision(double);
	void change_background(std::vector <float>);
	void change_rho(double);
	void change_resulting_file_name(string);
};

#endif /* PARAMETERFILE_H_ */
