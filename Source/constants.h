/*
 * constants.h
 *
 *  Created on: Aug 18, 2010
 *      Author: omidi
 */

#ifndef CONSTANTS_H_
#define CONSTANTS_H_

#include <boost/math/special_functions/gamma.hpp>



#define MAX(a,b) \
       ({ typeof (a) _a = (a); \
           typeof (b) _b = (b); \
         _a > _b ? _a : _b; })

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))


enum Nucleotides {A=65, C=67, G=71, T=84};	// The numbers that assigned are ASCII codes of ACGT
// this can be used in statements like "if (vec[1] == A)"
#define ALPH_NUM	4
// #define LAMBDA		0.125
#define LGAMMA_4LAMBDA	(2*(lgamma(4*LAMBDA)))
#define LGAMMA_LAMBDA	(lgamma(LAMBDA))
#define LOG_10 			(log(10.))
#define TINY 			(1e-20)

double FMAX(double,double);	// find the biggest parameter of types double

double ABS(double);	// find the Absolute value of a parameter of type double

short int convert(char);	// converting nucleotides to appropriate numbers for referring in vectors

//void sequence_file_score(string, Score);	// Open a file and read each sequence in each line and find the score of each sequence according to the alignment that assigned to the parameter of Score


#endif /* CONSTANTS_H_ */
