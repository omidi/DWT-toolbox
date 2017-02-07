/*
 * ParameterFile.cpp
 *
 *  Created on: Feb 28, 2011
 *      Author: omidi
 */

#include "ParameterFile.h"

ParameterFile::ParameterFile() {
	background.resize(4, 0.25);
	lowest_score = -50;
	file_name = "NULL";
	resulting_file_name = "NULL";
	background_is_set = false;
	precision = 15.;
	lambda = 0.125; // the default value is set to 1/8 to make have WM score
	rho = 1.;
	rho_is_set = false;
}

ParameterFile::ParameterFile (string file){
	file_name = file;
	ifstream  param_file (file_name.c_str(), ios::in);
	background_is_set = false;	// if background is not mentioned in param_file
	background.resize(4, 0.25);	// if background is not mentioned in param_file
	lowest_score = -80.0;	// if lowest_score is not mentioned in param_file
	resulting_file_name = "NULL";
	precision = 15.;
	lambda = 0.125; // the default value is set to 1/8 to make have WM score

	rho = 1.;	// the default value of rho is 1.
	rho_is_set = false;

	if (!param_file.is_open()){
		cerr << "The parameter file could not be opened." << endl;
		exit(1);
	} else{
		do{
			string temp;
			param_file >> temp;
			if (temp.find('#') < temp.length()){
				if (temp.find("background_model") < temp.length()) {
					initialize_background(param_file);
					background_is_set = true;
				}
				if (temp.find("output_file") < temp.length()) {
					initialize_resulting_file_name(param_file);
				}
			}
			if (temp.find("//") < temp.length()){
				while(temp.find('\n') < temp.length()) {
					param_file >> temp;
				}
			}

		}
		while (! param_file.eof());
	}

}

ParameterFile::~ParameterFile() {
	// TODO Auto-generated destructor stub
}

void ParameterFile::initialize_lambda(ifstream & param_file){
	string temp;
	param_file >> temp;
	lambda = atof(temp.c_str());
	if (lambda == 0) {
		cout << "The psudo-count is set to zero." << endl;
	}
	if (lambda < 0) {
		cerr << "Negative value for Lambda (psudo-count) is not valid!" << endl;
		exit(0);
	}
}


void ParameterFile::initialize_precision(ifstream & param_file){
	string temp;
	param_file >> temp;
	precision = atof(temp.c_str());
	if (precision <= 0) {
		cerr << "The precision value in parameter file is not valid and will be set to 5.0" << endl;
		precision = 5.;
	}
}

void ParameterFile::initialize_rho(ifstream & param_file){
		string temp;
		param_file >> temp;
		rho = atof(temp.c_str());
		if (rho < 0 || rho > 1) {
			cerr << "The rho value in parameter file is not valid and will be set to 1.0" << endl;
			rho = 1.;
		}
}

void ParameterFile::initialize_minimum_score(ifstream & param_file){
	string temp;
	param_file >> temp;
	lowest_score = atof(temp.c_str());
}

void ParameterFile::initialize_resulting_file_name(ifstream & param_file){
	string temp;
	param_file >> temp;
	if (temp.length() > 0)
		resulting_file_name = temp.c_str();
	else
		cerr << "The string provided for the resulting file name is not valid" << endl;
}

void ParameterFile::initialize_background(ifstream & param_file){
//	background.resize(4);
	for (unsigned short int i = 0 ; i < 4 ; i++)
		background[0] = .0;

	string temp;
	short int index = 0;
	float probability = .0;

	for (unsigned short int i = 0; i < 4; i++){
		param_file >> temp;
		index = convert(temp.c_str()[0]);
		if (index == -1){
			cerr << "Watch out your background probabilities are wrong!" << endl;
			exit(1);
		}
		param_file >> temp;
		probability = atof(temp.c_str());
		background[index] = probability;
	}

	// check if the background is normalized
	float sum = 0.;
	for (std::vector<float>::iterator i = background.begin();
			i < background.end(); i++)
		sum += *i;
	// if (sum != 1.0){
	// 	cerr << "Watch out your background is not normalized correctly." << endl;
	// 	exit(1);
	// }
}

bool ParameterFile::is_background() {
	return background_is_set;
}

bool ParameterFile::is_rho() {
	return rho_is_set;
}

bool ParameterFile::is_resulting_file_name() {
	if (resulting_file_name == "NULL")
		return false;
	else
		return true;
}

double ParameterFile::give_lambda() {
	return lambda;
}

double ParameterFile::give_rho(){
	return rho;
}

double ParameterFile::give_precision(){
	return precision;
}

string ParameterFile::give_resulting_file_name(){
	return resulting_file_name;
}

std::vector <float> ParameterFile::give_background(){
	return background;
}

double ParameterFile::give_lowest_score(){
	return lowest_score;
}

void ParameterFile::change_minimum_score(double min){
	lowest_score = min;
}

void ParameterFile::change_resulting_file_name(string new_name){
	if (new_name.length() > 0) {
		resulting_file_name = new_name.c_str();
	} else {
		cerr << "The string provided for the resulting file name is not valid" << endl;
	}
}

void ParameterFile::change_rho(double new_rho){
	if (!(new_rho < 0 && new_rho > 1)){
		rho = new_rho;
	} else {
		cout << "The new rho value is not valid. It might be between 0 and 1." << endl;
	}
	rho_is_set = true;
}

void ParameterFile::change_precision(double new_val) {
	if (new_val <= 0) {
		cerr << "The precision value provided is not valid." << endl;
	} else {
		precision = new_val;
	}
}

void ParameterFile::change_background(std::vector <float> probs){
	if (probs.size() == 4){
		background [convert('A')] = probs[convert('A')];
		background [convert('C')] = probs[convert('C')];
		background [convert('G')] = probs[convert('G')];
		background [convert('T')] = probs[convert('T')];
	} else{
		cerr << "The vector which was provided for changing background probabilities is inappropriate." << endl;
	}

	// check if the background is normalized
	float sum = 0.;
	for (std::vector<float>::iterator i = background.begin();
			i < background.end(); i++)
		sum += *i;
//	if (sum != 1.0){
//		cerr << "Watch out your background is not normalized correctly." << endl;
//		exit(1);
//	}

}
