/*
 * Alignment.cpp
 *
 *  Created on: Aug 18, 2010
 *      Author: omidi
 */

#include "Alignment.h"
using namespace std;

Alignment::Alignment() {
	//cerr << "A filename should be provided." << endl;
}

Alignment::~Alignment() {
	// TODO Auto-generated destructor stub
}

Alignment::Alignment(string file) {
	alignment = read_alignment(file);
	filename = file;
}

vector<Motif> Alignment::read_alignment(string file){
//Motif* Alignment::read_alignment(string file){
	vector<Motif> M;
//	Motif * M;
	ifstream File(file.c_str(), ios::in);
	if (!File) {
		cerr << "There is no such a file or directory: " << file << endl;
		exit(1);
	} else {
		string line;
		vector<string> elems;
		int line_number = 0;
		int num_positions = 0;
		while(!File.eof()){
			getline(File, line);
			if(line.find("NA") == 0){
				boost::split(elems, line, boost::is_any_of("\t "));
				TF = elems[1];
			}
			if((line.find("\\") == 0) || (line.find("PO") == 0)){
				continue;
			}
			boost::split(elems, line, boost::is_any_of("\t "));
			if(elems.size() == 18){
				line_number++;
			}
		}
		num_positions = int((1 + sqrt(8*line_number+1))/2);
		initialize_n_i(num_positions);
		initialize_n_ij(num_positions);
		File.close();
		ifstream File(file.c_str(), ios::in);
		int index = 0;
		int line_count = 0;
		while(!File.eof()){
			getline(File, line);
			if(line.find("NA") == 0){
				continue;
			}
			if((line.find("\\") == 0) || (line.find("PO") == 0)){
				continue;
			}
			boost::split(elems, line, boost::is_any_of("\t "));
			if(elems.size() == 18){
				index = 2;
				for(int a=0; a<4; a++){
					for(int b=0; b<4; b++){
						n_ij[a][b][atoi(elems[1].c_str()) - 1][atoi(elems[0].c_str()) - 1] =
								n_ij[a][b][atoi(elems[0].c_str()) - 1][atoi(elems[1].c_str()) - 1] = atof(elems[index].c_str());
						index++;
					}
				}
				index = 2;
				for(int a=0; a<4; a++){
					n_i[a][atoi(elems[0].c_str()) - 1] = atof(elems[index].c_str()) +
							atof(elems[index+1].c_str()) +
							atof(elems[index+2].c_str()) +
							atof(elems[index+3].c_str());
					index += 4;
					}
				line_count++;
				if(line_count==line_number){
					index = 2;
					for(int a=0; a<4; a++){
						n_i[a][atoi(elems[1].c_str()) - 1] = atof(elems[index].c_str()) +
								atof(elems[index+4].c_str()) +
								atof(elems[index+8].c_str()) +
								atof(elems[index+12].c_str());
						index++;
						}
				}
				n_rows = float(n_i[0][0] + n_i[1][0] + n_i[2][0] + n_i[3][0]);
				n_cols = num_positions;
			}
		}
		File.close();
	}
	return M;
}


string Alignment::get_TF(){
	return TF;
}

Motif Alignment::operator[](unsigned int index){
	try {
		return alignment[index];
	} catch(exception& e){
		cerr << "Exception occurred: "<< e.what() <<endl;
	}
	return Motif("ERROR OCCURRED!", 1);
}

float Alignment::nrows(){
	return n_rows;
}

unsigned short int Alignment::ncols(){
	return n_cols;
}

float Alignment::get_n_i(char alpha, unsigned short int i){
	return n_i[convert(toupper(alpha))][i];	// I should add an exception here in order to prevent it from index out of range
}

float Alignment::get_n_i(unsigned short int alpha, unsigned short int i){
	return n_i[alpha][i];
}

float Alignment::get_n_ij(char alpha, char beta, unsigned short int i, unsigned short int j){
	return n_ij[convert(toupper(alpha))][convert(toupper(beta))][i][j];	// same as get_n_i();
}

float Alignment::get_n_ij(unsigned short int alpha, unsigned short int beta, unsigned short int i, unsigned short int j){
	return n_ij[alpha][beta][i][j];
}


void Alignment::initialize_n_i(unsigned short int cols){
	n_i.resize(ALPH_NUM);
	for(unsigned int i = 0; i < ALPH_NUM ; i++){
		n_i[i].resize(cols, 0);
	}
//	n_i = (float **) malloc(ALPH_NUM*sizeof(float));
//	for(unsigned int i=0; i < ALPH_NUM ; i++){
//		n_i[i] = (float*) malloc(cols*sizeof(float));
//		for(unsigned int j = 0; j < cols ; j++){
//			n_i[i][j] = .0;
//		}
//	}
}

void Alignment::initialize_n_ij(unsigned short int cols){
	n_ij.resize(ALPH_NUM);
	for(unsigned int i=0; i< ALPH_NUM; i++){
		n_ij[i].resize(ALPH_NUM);
		for(unsigned int j=0; j<ALPH_NUM ; j++){
			n_ij[i][j].resize(cols);
			for(unsigned int k = 0; k<cols ; k++){
				n_ij[i][j][k].resize(cols, 0);
			}
		}
	}
}






