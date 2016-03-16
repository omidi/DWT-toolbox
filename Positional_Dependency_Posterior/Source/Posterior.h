/*
 * Posterior.h
 *
 *  Created on: Feb 8, 2012
 *      Author: omidi
 */

#ifndef POSTERIOR_H_
#define POSTERIOR_H_

#include "Score.h"
#include "Decomposition.h"


using namespace std;

class Posterior {
	std::vector <std::vector<double> > logR;
public:
//	Posterior();
	Posterior(Score);
	virtual ~Posterior();
	std::vector <std::vector<double> > contract_edge(std::vector <std::vector<double> >,unsigned int,unsigned int);
	void calculate_posterior(Score);
	double calculate_determinant(std::vector <std::vector<double> >, int);
};

#endif /* POSTERIOR_H_ */
