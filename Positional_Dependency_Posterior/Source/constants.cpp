/*
 * functions.cpp
 *
 *  Created on: Aug 19, 2010
 *      Author: omidi
 */

#include "constants.h"

using namespace std;


double FMAX(double a, double b){
	return (a > b)? a : b;
}

double ABS(double a){
	return (a<0)?-a:a;
}

short int convert(char ch){
	switch(ch){
	case A:
		return 0;
	case C:
		return 1;
	case G:
		return 2;
	case T:
		return 3;
	}
	return -1;
}

//void sequence_file_score(string file, Score score){
//	time_t start, end;
//	time(&start);
//	ifstream FILE(file.c_str(), ios::in);
//	string result_filename = file + "_result";
//	ofstream RESULT(result_filename.c_str());
//	if (!FILE) {
//		cerr << "There is no such a file or directory: " << file << endl;
//		exit(1);
//	} else {
//		while(! FILE.eof()){
//			string seq;
//			FILE >> seq;
//			RESULT << seq << '\t' << score.score_of_sequence(seq) << endl;
//		}
//		FILE.close();
//		RESULT.close();
//	}
//	time(&end);
//	cout << "Time for computation is: "<< difftime(end, start);
//}


