/*
 * BackgroundModel.cpp
 *
 *  Created on: Aug 26, 2010
 *      Author: omidi
 */

#include "BackgroundModel.h"
#include <math.h>
#include <iostream>


BackgroundModel::BackgroundModel() {
/*
 * If the background probabilities are not provided,
 * all of them will be set to 0.25, or uniform probability which shows
 * the state of total ignorance.
 */
	A = log(0.25);
	C = log(0.25);
	G = log(0.25);
	T = log(0.25);
}


BackgroundModel::BackgroundModel (std::vector <float> probs){
	A = log(probs[convert('A')]);
	C = log(probs[convert('C')]);
	G = log(probs[convert('G')]);
	T = log(probs[convert('T')]);
}

BackgroundModel::~BackgroundModel() {
	// TODO Auto-generated destructor stub
}

double BackgroundModel::uniform_background(int n){
	return  -n*log(4);
}

double BackgroundModel::non_uniform_background(string seq){
	double prob = 0.0;

	for (string::iterator s = seq.begin();
			s != seq.end(); ++s){
		switch (*s){
		case 'A':
			prob += A;
			break;
		case 'C':
			prob += C;
			break;
		case 'G':
			prob += G;
			break;
		case 'T':
			prob += T;
			break;
		}
	}
	return prob;
}

void BackgroundModel::change_background(std::vector <float> probs){
	A = log(probs[convert('A')]);
	C = log(probs[convert('C')]);
	G = log(probs[convert('G')]);
	T = log(probs[convert('T')]);

}
