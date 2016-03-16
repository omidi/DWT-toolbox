//============================================================================
// Name        : DWT_model.cpp
// Author      : Saeed Omidi
// Version     :
// Copyright   : Your copyright notice
// Description : Dinucleotide Weighted Tensor model in C++, Ansi-style
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
	if (argc != 4){
		cerr << "The DWT model, input sequences in FASTA format, and a parameter file." << endl;
		cerr << "e.g. CTCF.dwt sequences.fasta param_file.txt" << endl;
		exit(1);
	}

	string dwt_model_file = argv[1];
	Alignment alignment(dwt_model_file);
	ParameterFile parameters(argv[3]);
	Score S(alignment, alignment.get_TF(), parameters);
	Window win(argv[2] , parameters, S, alignment.ncols());
	return 0;
}
