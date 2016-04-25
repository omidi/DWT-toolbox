/*
 * Motif.h
 *
 *  Created on: Aug 18, 2010
 *      Author: omidi
 */

#ifndef MOTIF_H_
#define MOTIF_H_

#include <string>
using namespace std;

class Motif {
	string motif;
	float probability;
public:
	Motif();
	virtual ~Motif();
	Motif(string,float);
	char operator[](int);
	float score();
	int len();
	string sequence();
};

#endif /* MOTIF_H_ */
