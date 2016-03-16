/*
 * BackgroundModel.h
 *
 *  Created on: Aug 26, 2010
 *      Author: omidi
 */

#ifndef BACKGROUNDMODEL_H_
#define BACKGROUNDMODEL_H_

#include <string>
#include <iostream>
#include <vector>
#include "constants.h"


using namespace std;

class BackgroundModel {
private:
	float A,C,G,T;
public:
	BackgroundModel();
	BackgroundModel(string);	// TODO: takes a filename that contains the background probabilities
	BackgroundModel(std::vector <float>);
	virtual ~BackgroundModel();
	double uniform_background(int);	// The simplest form of background 1/(4^n) in log scale, means: -n*log(4)
	double non_uniform_background(string);	// takes the sequence and gives back the background probability of the sequence.
	void change_background (std::vector <float>);
};

#endif /* BACKGROUNDMODEL_H_ */
