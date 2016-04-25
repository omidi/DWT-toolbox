/*
 * Motif.cpp
 *
 *  Created on: Aug 18, 2010
 *      Author: omidi
 */

#include "Motif.h"

Motif::Motif(string str, float prob) {
	motif = str;
	probability = prob;
}

char Motif::operator [](int i){
	return toupper(motif[i]);
}


float Motif::score(){
	return probability;
}

int Motif::len(){
	return motif.length();
}

string Motif::sequence(){
	return motif;
}

Motif::Motif() {
	// TODO Auto-generated constructor stub

}

Motif::~Motif() {
	// TODO Auto-generated destructor stub
}
