//============================================================================
// Name        : DWT_model.cpp
// Author      : Saeed Omidi
// Version     :
// Copyright   : Your copyright notice
// Description : Positional Dependency based on DWT model in C++, Ansi-style
//============================================================================

#include <iostream>
#include <boost/math/special_functions/gamma.hpp>
#include <time.h>
#include <vector>
#include <iomanip>
#include "Motif.h"
#include "constants.h"
#include "Alignment.h"
#include "Decomposition.h"
#include "Score.h"
#include "Window.h"
#include "ParameterFile.h"


using namespace std;

int main(int argc, char* argv[]){
	if (argc != 2){
		cerr << "Calculates the pairwise dependency posterior between all pairs of positions." << endl;
		cerr << "e.g. CTCF.dwt" << endl;
		exit(1);
	}

	string dwt_model_file = argv[1];
	Alignment alignment(dwt_model_file);
	ParameterFile parameters;
	Score S(alignment, alignment.get_TF(), parameters);
	for(unsigned int i=0; i<alignment.ncols(); i++){
		for(unsigned int j=0;j<i;j++){
			printf("%d\t%d\t%.4lf\n", j, i, S.posterior(i,j));
		}
	}
	return 0;
}
