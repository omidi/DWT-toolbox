/*
 * Window.h
 *
 *  Created on: Oct 28, 2010
 *      Author: omidi
 */

#ifndef WINDOW_H_
#define WINDOW_H_

#include <string>
#include <iomanip>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/replace.hpp>
#include "Score.h"
#include "ParameterFile.h"

using namespace std;

class Window {
private:
	string filename;
	Score score;
	unsigned int l;
	double lowest_score;
	ParameterFile parameters;
public:
	Window();
	Window(string, Score , int);
	Window(string, Score , int, double);
	Window(string, ParameterFile, Score , int);
	virtual ~Window();
	inline void slide(string);	// here, the object slides a window along the genome, as an argument it takes the name of the output file
	inline void regional_slide(string, string, ofstream&);	//  takes a region of genome and slides the window on it
	inline char convert(char);	//  gives the complement nucleotide of the parameter
	inline string complement(string);	// gives the complement of the given site
	inline char char_is_valid(char);
	inline bool is_valid(string);
};

#endif /* WINDOW_H_ */
