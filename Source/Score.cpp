/*
 * Score.cpp
 *
 *  Created on: Aug 20, 2010
 *      Author: omidi
 */

#include "Score.h"

Score::Score() {
	// TODO Auto-generated constructor stub
}

Score::~Score() {
	// TODO Auto-generated destructor stub
}

Score::Score(Alignment& alg){
	PRECISION = 15.;	// the default value of PRECISION
	LAMBDA = .125;
	alignment = alg;
	alpha_exponent = 1.;
	logR = calculate_logR();
	std::vector <std::vector<double> > temp = shift(logR);
	alpha_exponent = find_alpha_exponent(temp);
	rescaled_R = re_scale(temp, alpha_exponent);

	laplacian_determinant = calculate_determinant(laplacian_of(rescaled_R), rescaled_R.size() - 1);

	tf = "TF-Unknown";

	ParameterFile p;
	parameters = p;

	BackgroundModel bg (parameters.give_background());
	background = bg;
}

Score::Score(Alignment& alg, string tf_name){
	PRECISION = 15.;	// the default value of PRECISION
	LAMBDA = .125;
	alignment = alg;
	alpha_exponent = 1.;
	logR = calculate_logR();

	std::vector <std::vector<double> > temp = shift(logR);
	alpha_exponent = find_alpha_exponent(temp);
	laplacian_determinant = calculate_determinant(laplacian_of(rescaled_R), rescaled_R.size() - 1);

	tf = tf_name;

	ParameterFile p;
	parameters = p;

	BackgroundModel bg (parameters.give_background());
	background = bg;
}

Score::Score(Alignment& alg, string tf_name, ParameterFile params){
	alpha_exponent = 1.;
	alignment = alg;
	parameters = params;
	PRECISION = parameters.give_precision();
	LAMBDA = parameters.give_lambda();
	tf = tf_name;
	BackgroundModel bg (parameters.give_background());
	background = bg;

	logR = calculate_logR();
	if (!parameters.is_rho()) {
		double rho = 0.;
		laplacian_determinant = -10000000000000000;	// to make sure it will be changed at the first run of algorithm
		alpha_exponent = find_alpha_exponent(logR);
		rescaled_R = re_scale(logR, alpha_exponent);
		double fitted_rho = 1.;
		double det_temp = 0.;
        for (short int i = 0; i < 101; i++){
			det_temp = calculate_determinant(laplacian_of(shift(rescaled_R, rho)), logR.size() - 1);
			if (det_temp > laplacian_determinant){
				laplacian_determinant = det_temp;
				fitted_rho = rho;
			}
			rho += 0.01;
        }
        parameters.change_rho(fitted_rho);
	}
	else
	{
		alpha_exponent = find_alpha_exponent(logR);
		rescaled_R = re_scale(logR, alpha_exponent);
		laplacian_determinant = calculate_determinant(laplacian_of(shift(rescaled_R,parameters.give_rho())), rescaled_R.size()-1);
	}
}


double Score::determinant(){
  return laplacian_determinant;
}


double Score::fit_rho(){
    double rho = 0.;
    double alpha = find_alpha_exponent(logR);
    double determinant = -10000000000000000;        // to make sure it will be changed at the first run of algorithm
    double fitted_rho = 1.;
    double det_temp = 0.;
    for (short int i = 0; i < 101; i++){
            det_temp = calculate_determinant(laplacian_of(shift(re_scale(logR, alpha), rho)), logR.size() - 1);
            cout << rho << '\t' << det_temp << endl;
            if (det_temp > determinant){
                    determinant = det_temp;
                    fitted_rho = rho;
            }
            rho += 0.01;
}
return fitted_rho;
}


std::vector <std::vector<double> > Score::calculate_logR(){
	std::vector <std::vector<double> > M;
	M.resize(alignment.ncols());
	for(unsigned int i = 0; i < alignment.ncols(); i++){
		M[i].resize(alignment.ncols(), .0);
	}

	double logP_indep = 0.;
	double logP_dep = 0;

	for (unsigned short int i=0; i < alignment.ncols() ; i++){
	  for (unsigned short int j=0; j < i ; j++){
		logP_indep = 0.0;
	    logP_dep = lgamma(alignment.nrows()+16*LAMBDA) - lgamma(16*LAMBDA);
		  for(unsigned short int alpha = 0; alpha < ALPH_NUM; alpha++){
		    logP_indep += lgamma(alignment.get_n_i(alpha, i) + 4*LAMBDA);	// Gamma(n_alpha^i + Lambda)
		    logP_indep += lgamma(alignment.get_n_i(alpha, j) + 4*LAMBDA);	// Gamma(n_beta^j + Lambda)
		    logP_indep -= LGAMMA_4LAMBDA; // Gamma(4Lambda)^2 _ defined in constants.h

		    for(unsigned short int beta = 0; beta < ALPH_NUM; beta++){
		      logP_dep += lgamma(alignment.get_n_ij(alpha,beta,i,j) + LAMBDA);	// Gamma(n_ij^alpha,beta + Lambda)
		      logP_dep -= LGAMMA_LAMBDA;	// Gamma(Lambda) _ defined in constants.h
		    }
		  }
		  M[i][j] = M[j][i] = (logP_dep-logP_indep); // last version
		}

	}
	return M;
}


std::vector <std::vector<double> > Score::shift(std::vector <std::vector<double> > logR, double rho){
	/*
	 * This method is a polymorphism of the shift() method and it will be used for fitting rho.
	 * This method is used in case the rho value is not set inside parameter file.
	 */

	std::vector <std::vector<double> > M;
	M.resize(alignment.ncols());
	for(unsigned int i = 0; i < alignment.ncols(); i++){
		M[i].resize(alignment.ncols(), .0);
	}
	  for (unsigned int l=0 ; l < alignment.ncols() ;l++){
		for (unsigned int k=0; k < l ; k++){
			M[l][k] = rho*exp(logR[l][k])  - rho + 1.;
			M[k][l] = M[l][k];
		 }
	  }
	  return M;
}

std::vector <std::vector<double> > Score::add_one(std::vector <std::vector<double> > logR){
	/*
	 * This method is a polymorphism of the shift() method and it will be used for fitting rho.
	 * This method is used in case the rho value is not set inside parameter file.
	 */

	std::vector <std::vector<double> > M;
	M.resize(alignment.ncols());
	for(unsigned int i = 0; i < alignment.ncols(); i++){
		M[i].resize(alignment.ncols(), .0);
	}
	  for (unsigned int l=0 ; l < alignment.ncols() ;l++){
		for (unsigned int k=0; k < l ; k++){
			M[l][k] = (logR[l][k]) + 1.;
			M[k][l] = M[l][k];
		 }
	  }
	  return M;
}


std::vector <std::vector<double> > Score::shift(std::vector <std::vector<double> > logR){
	/*
	 * It's a kind of re-scaling. Here the minimum entry in the matrix will be subtracted from all the entries of the matrix
	 * In addition the average of rows and column will be reduced.
	 */
	double rho = parameters.give_rho();
	std::vector <std::vector<double> > M;
	M.resize(alignment.ncols());
	for(unsigned int i = 0; i < alignment.ncols(); i++){
		M[i].resize(alignment.ncols(), .0);
	}

	  for (unsigned int l=0 ; l < alignment.ncols() ;l++){
		for (unsigned int k=0; k < l ; k++){

		  // rho part
		  if (M[k][l] > 0) {
			  double temp = rho*exp(logR[k][l]) + (1.0 - rho);
			  M[k][l] = (temp);
		  } else {
			  double temp = log(rho) + logR[k][l] + log (1.0 + ((1.0 - rho)/rho)*exp(-logR[k][l]));
			  M[k][l] = exp(temp);
		  }
		  M[l][k] = M[k][l];
		 }
	  }

	  return M;
}


double Score::find_alpha_exponent(std::vector <std::vector<double> > logR) {
	double max = logR[0][1];
	double min = logR[0][1];
	double alpha = 1.;

	for (unsigned int i = 0; i < logR.size(); i++){
		for(unsigned int j = 0; j < i; j++){
			if (logR[i][j] > max)
				max = logR[i][j];
			if (logR[i][j] < min)
				min = logR[i][j];
		}
	}
	if (max  > PRECISION){
		alpha = PRECISION / max;
	}
	return alpha;
}

std::vector <std::vector<double> > Score::re_scale(std::vector <std::vector<double> > logR, double alpha_exponent){
/*
 * The re-scaling formula is as following on:
 * 		R -> R^alpha
 * 		where alpha = (K*log10)/(log_10(max) - log_10(min))
 * 		K determines the precision, and should be set based on the numerical capability of the machine, in original paper this value
 * 		set as K=5. In this code K refers as PRECISION (constant) class member.
 * 		the function of this re-scaling method is to shrink the data values in a way the their relative difference will be conserved, however in
 * 		smaller scale.
 * 		For more info refers to "Tractable Bayesian Learning of Tree Augmented Naive Bayes Classifiers", Jesus Cerquides and Ramon de Mantaras, 2003
 * */
	std::vector <std::vector<double> > M;
	M.resize(alignment.ncols());
	for(unsigned int i = 0; i < alignment.ncols(); i++){
		M[i].resize(alignment.ncols(), .0);
	}
	for (unsigned int i = 0; i < logR.size(); i++){
		for(unsigned int j = 0; j < i; j++){
				M[j][i] = M[i][j] = alpha_exponent*logR[i][j];
		}
	}
	return M;
}

double Score::transf(double log_x, double min_logW, double max_logW, double alpha){
	  double x_transf=pow(10., alpha*(log_x-min_logW) - PRECISION + 2);
	  if (isinf(x_transf)==1){
		  return 1.0;
	  }
	  else
		  if (isinf(x_transf)==-1){
		  return .0;
	 	}
	  return x_transf;//not in log_scale!
}

inline double Score::calculate_determinant(std::vector <std::vector<double> > M, int n){
	Decomposition D;
	double determinant = 1.0;

	determinant = D.determinant(M, n);
	return determinant;
}

std::vector<std::vector<double> > Score::contract_edge(std::vector <std::vector<double> > logR, unsigned int node1, unsigned int node2){
	std::vector <std::vector<double> > M;
	M.resize(logR.size()-1);
	for(unsigned int i = 0; i < logR.size()-1; i++){
		M[i].resize(logR.size(), .0);
	}
  for(unsigned int i=0;i<logR.size();i++){
    if((i!=node1) && (i!=node2)){
      for(unsigned int j=0;j<i;j++){
	if((j!=node1) && (j!=node2)){
	  unsigned int index1=i;
	  unsigned int index2=j;
	  //take out rows/columns given by node1 and node2
	  if(index1>=node1){index1--;}
	  if(index1>=node2){index1--;}
	  if(index2>=node1){index2--;}
	  if(index2>=node2){index2--;}
	  if((index1>=0) && (index2>=0)){
	    M[index1][index2]=logR[i][j];
	    M[index2][index1]=logR[i][j];
	  }
	}
      }
    }
  }
  //add contracted edge in the last column
  unsigned int last_index=logR.size()-2;
  for(unsigned int i=0;i<logR.size();i++){
    if((i!=node1) && (i!=node2)){
      unsigned int current_index=i;
      if(current_index>=node1){current_index--;}
      if(current_index>=node2){current_index--;}
      if(current_index>=0){
	M[last_index][current_index]=logR[node1][i]+logR[node2][i];
	M[current_index][last_index]=M[last_index][current_index];
      }
    }
  }
  return M;
}


double Score::posterior(unsigned int i, unsigned int j){
    if(i>logR.size() || j>logR.size()){
            return -1.;
    }
    double posterior = 0.;
    std::vector <std::vector<double> > M = contract_edge(shift(rescaled_R, parameters.give_rho()), i, j);
    double numerator_determinant = calculate_determinant(laplacian_of(M), M.size()-1);
    posterior = log(shift(rescaled_R, parameters.give_rho())[i][j] - 1 + parameters.give_rho()) + numerator_determinant;
    posterior -= laplacian_determinant;    
    return exp(posterior);
}

std::vector <std::vector<double> > Score::laplacian_of(std::vector <std::vector<double> > M){
	std::vector <std::vector<double> > L;
	L.resize(M.size());
	for(unsigned short int i = 0; i < M.size(); i++)		L[i].resize(M.size(), .0);
	for(unsigned short int i = 0; i < M.size(); i++){
		double sum = .0;
		for(unsigned short int j = 0; j< M.size(); j++){
			sum += M[i][j];
			L[i][j] = -M[i][j];
		}
		L[i][i] = sum ; //	I didn't use "- M[i][i]" because logR[i][i] is 0
	}
	return L;
}

double Score::score_of_sequence(string seq){
	/*
	 * in this method score of a given sequence will be calculated.
	 * the formula for this calculation is:
	 * P(s|S) = [det(L(new_logR)) / det(L(old_logR))] * \prod_i^|seq| (n^i_si + 4*LAMBDA)
	 */
	double score = 0.0;
	std::vector <std::vector<double> > temp = logR_plus_seq(seq);
	std::vector <std::vector<double> > R_new = shift(temp, parameters.give_rho());

	double determinant_new_R = calculate_determinant(laplacian_of(R_new), R_new.size() - 1);

	for(unsigned short int i = 0; i < seq.size(); i++){
		score += log(alignment.get_n_i(seq[i], i) + 4*LAMBDA);	// multiply by (n_{s_i}^i + 4*LAMBDA)	(This part changed by Erik to change the pseudo-count part)
//		score -= log(alignment.nrows() + 16*LAMBDA);	// multiply by 1/(n + 16*LAMBDA)
	}
	score -= alignment.ncols()*log(alignment.nrows() + 16*LAMBDA); // multiply by 1/(n + 16*LAMBDA)^l which l is the length of TF binding sequence

	score += determinant_new_R;	    // det(L(min(R_new))) / det(L(min(R_old)))
	score -= laplacian_determinant;	// determinant is always positive
									// because determinant method gives its absolute value
	score -= background.non_uniform_background(seq);
	return score;	// log-odd of the given sequence
}

double Score::wm_score(string seq){
	double score = 0.0;
	for(unsigned short int i = 0; i < seq.size(); i++){
		score += log(alignment.get_n_i(seq[i], i) + 4*LAMBDA);	// multiply by (n_{s_i}^i + 4*LAMBDA)
	}
	score -= alignment.ncols()*log(alignment.nrows() + 16*LAMBDA); // multiply by 1/(n + 16*LAMBDA)^l which l is the length of TF binding sequence
	score -= background.non_uniform_background(seq);
	return score;
}

std::vector <std::vector<double> > Score::logR_plus_seq(string seq){
	/* instead of counting again number of appearances for each letter and generate new n_i and n_ij matrices,
	 * here, logR of the alignment plus new sequence will be calculated according to the following formula:
	 * logR_ij^s = log(R_ij) + log(n^ij_{si,sj} + LAMBDA) - log(n^i_si + 4*LAMBDA) - log(n^j_sj + 4*LAMBDA) - log(n + 16*LAMBDA) + 2log(n+4*LAMBDA)
	 * where, logR_ij^s is an entry of the new dependency matrix and R_ij is an entry of previously computed dependency matrix R
	 * and n^ij, n^i are frequency matrix that we have them from an object of Alignment class.
	*/

	// the original logR replaced with rescaled_R in following
	std::vector <std::vector<double> > logR_new;
	logR_new.resize(rescaled_R.size());	// initialization of logR_new matrix
	for(unsigned short int i = 0; i < rescaled_R.size(); i++)		logR_new[i].resize(rescaled_R.size(), .0);

	for(unsigned short int i = 0; i < seq.size(); i++){
		for(unsigned short int j = 0; j < i; j++){
		  logR_new[i][j] = rescaled_R[i][j] + log(alignment.get_n_ij(seq[j], seq[i], j, i) + LAMBDA);
		  logR_new[i][j] -= log(alignment.nrows()+16*LAMBDA);
		  logR_new[i][j] -= log(alignment.get_n_i(seq[j], j) + 4*LAMBDA);
		  logR_new[i][j] -= log(alignment.get_n_i(seq[i], i) + 4*LAMBDA);
		  logR_new[i][j] += 2*log(alignment.nrows() + 16*LAMBDA);

		  logR_new[j][i] = logR_new[i][j];
		}
	}

	return logR_new;
}

string Score::get_TF(){
	return tf;
}


