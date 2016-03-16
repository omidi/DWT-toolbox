/*
 * Score.h
 *
 *  Created on: Aug 20, 2010
 *      Author: omidi
 */

#ifndef SCORE_H_
#define SCORE_H_

#include "Alignment.h"
#include <boost/math/special_functions/gamma.hpp>
#include <vector>
#include "Decomposition.h"
#include "BackgroundModel.h"
#include "ParameterFile.h"

using namespace std;

class Score {
	Alignment alignment;
	std::vector <std::vector<double> > logR;
	std::vector <std::vector<double> > rescaled_R;
	double laplacian_determinant;
	string tf;
	double PRECISION;	// for re-scaling the numbers in LogR (default 5.0)
	double LAMBDA;
	ParameterFile parameters;
	BackgroundModel background;
	double alpha_exponent;	// in order to re-scale all logRs by the same alpha_exponent for the original logR matrix
public:
	Score();
	Score(Alignment&);
	Score(Alignment&, string);
	Score(Alignment&, string, ParameterFile);
	virtual ~Score();
	std::vector <std::vector<double> > calculate_logR();
	std::vector <std::vector<double> > shift(std::vector <std::vector<double> >);
	std::vector <std::vector<double> > shift(std::vector <std::vector<double> >, double);	// this method is written in order to use it for fitting rho
	std::vector <std::vector<double> > re_scale(std::vector <std::vector<double> >, double);
	double find_alpha_exponent(std::vector <std::vector<double> >); // in order to find the alpha_exponent of a given matrix
	double determinant();
	inline double calculate_determinant(std::vector <std::vector<double> >, int); // int parameter helps to remove as many columns and rows as we like
	std::vector <std::vector<double> > laplacian_of(std::vector <std::vector<double> >);	// finding laplacian matrix of logR which has been initialized before in calculate_logR
	double score_of_sequence(string);	// calculate the log-odd of the sequence according to the alignment
	std::vector <std::vector<double> > logR_plus_seq(string);	// here, logR of the alignment in addition to the new sequence will be calculated
	string get_TF();
	double transf (double, double, double, double);
	double wm_score(string);
	std::vector<std::vector<double> > contract_edge(std::vector <std::vector<double> >, unsigned int, unsigned int);
	double posterior(unsigned int, unsigned int);
	double fit_rho();
	std::vector <std::vector<double> > add_one(std::vector <std::vector<double> >);
};


#endif /* SCORE_H_ */
