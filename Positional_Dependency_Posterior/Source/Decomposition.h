/*
 * Decomposition.h
 *
 *  Created on: Aug 23, 2010
 *      Author: omidi
 */

#ifndef DECOMPOSITION_H_
#define DECOMPOSITION_H_

#include <math.h>
#include "constants.h"
using namespace std;


class Decomposition {
public:
	Decomposition();
	virtual ~Decomposition();
	std::vector<double>  QR_decomposition(std::vector<std::vector<double> >, int);	// returns the main diagonal elements of R
	double determinant (std::vector<std::vector<double> >, int);
};

#endif /* DECOMPOSITION_H_ */
