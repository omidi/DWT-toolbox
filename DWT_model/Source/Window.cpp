/*
 * Window.cpp
 *
 *  Created on: Oct 28, 2010
 *      Author: omidi
 */

#include "Window.h"

Window::Window() {
	// TODO Auto-generated constructor stub
}

Window::Window(string str, Score sc , int len){
	filename = str;
	score = sc;
	l = len;
	slide(sc.get_TF() + "_result");
	lowest_score = -50.;
}

Window::Window(string str, Score sc , int len, double minimum_score){
	filename = str;
	score = sc;
	l = len;
	lowest_score = minimum_score;
	slide(sc.get_TF() + "_result");
}

Window::Window(string str, ParameterFile params, Score sc ,int len){
	filename = str;
	score = sc;
	l = len;
	lowest_score = params.give_lowest_score();

	if (!params.is_resulting_file_name())
		slide(sc.get_TF() + "_result");
	else
		slide(params.give_resulting_file_name().c_str());
}


Window::~Window() {

}

inline string Window::complement(string str){
	string c_str = str;
	unsigned short int len = str.length();

	for (unsigned short int i = 0; i < len / 2 ; i++){
		c_str [i] = convert(str[(len-i-1)]);
		c_str [(len-i-1)] = convert(str[i]);
	}
	if (len%2)	// if the length of sequence was odd. otherwise, it doesn't do anything.
		c_str [len/2] = convert(str[len/2]);

	return c_str;
}

inline char Window::convert(char c){
	switch (c){
	case 'A':
		return 'T';
	case 'T':
		return 'A';
	case 'C':
		return 'G';
	case 'G':
		return 'C';
	default:
		return '*';
	}
}

inline char Window::char_is_valid(char c){
	switch (c){
	case 'A':
		return true;
	case 'T':
		return true;
	case 'C':
		return true;
	case 'G':
		return true;
	default:
		return false;
	}
}

inline void Window::slide(string result_file_name){
	ifstream FILE (filename.c_str(), ios::in);
	if(!FILE.is_open()){
		cerr << "The file is not open!" << endl;
	}
	else{
		ofstream RESULT (result_file_name.c_str(), ios::out);
		string sequence, specie;
		while (true){
			FILE >> specie;

			if (FILE.eof())	// if it's not here, we read a sequence twice! for example if the last sequence is AAAAAA, it considers the sequence twice!
				break;
//			FILE >> chr;
//			FILE >> position;
			FILE >> sequence;
            specie = specie.substr(1);
			regional_slide (specie, sequence, RESULT);

		}
		RESULT.close();
	}
}

inline bool Window::is_valid(string str){
	for (unsigned short int i = 0; i < str.length(); i++){
		if (!char_is_valid(str[i]))
			return false;
	}

	return true;
}


inline void Window::regional_slide(string specie, string sequence, ofstream& RESULT){
	double temp = .0;	// this variable holds the score of the sequence to prevent re-calculation
	double wm_score = .0;
	string position="0";
	string c_sequence;	// this variable holds the complement of the sequence
	if (sequence.length() >= l)  // this is to check whether the sequence is at least as big as the motif length
		for (unsigned short int i = 0 ; i <= (sequence.length() - l) ; i++){
			if (is_valid(sequence.substr(i,l)))
			{
				try{
					temp = score.score_of_sequence(sequence.substr(i,l));
					wm_score = score.wm_score(sequence.substr(i,l));
				}
				catch (exception& e){
					cout << e.what() << endl;
					cout << sequence << endl;
				}
				if (temp > lowest_score && wm_score > lowest_score){
					RESULT << specie << '\t';
					RESULT << atoi (position.substr(0 , position.find_first_of('-')).c_str()) + i << '\t' << atoi (position.substr(0 , position.find_first_of('-')).c_str()) + l + i << '\t';
					RESULT << '+' << '\t';
					RESULT << score.get_TF() << '\t';
					RESULT << sequence.substr(i,l) << '\t';
//					RESULT << fixed << showpoint << setprecision(4) << wm_score << '\t';		// here it writes the 'Weight Matrix' score of the sequence
					RESULT << fixed << showpoint << setprecision(4) << temp << '\t';		// here it writes the score of the sequence
					RESULT << '\n';
				}
				// sequence on the opposite strand, reverse compliment of the last sequence.
				c_sequence = complement(sequence.substr(i,l));
				temp = score.score_of_sequence(c_sequence);
				wm_score = score.wm_score(c_sequence);
				if (temp > lowest_score && wm_score > lowest_score){
					RESULT << specie << '\t';
					RESULT << atoi (position.substr(0 , position.find_first_of('-')).c_str()) + i << '\t'  << atoi (position.substr(0 , position.find_first_of('-')).c_str()) + l + i << '\t';
					RESULT << '-' << '\t';
					RESULT << score.get_TF() << '\t';
					RESULT << c_sequence << '\t';	// the complement of the sequence
//					RESULT << fixed << showpoint << setprecision(4) << wm_score << '\t';		// here it writes the 'Weight Matrix' score of the sequence
					RESULT << fixed << showpoint << setprecision(4) << temp << '\t';		// here it writes the score of the sequence
					RESULT << '\n';
				}
			}
		}
}
