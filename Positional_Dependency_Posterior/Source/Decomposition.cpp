/*
 * Decomposition.cpp
 *
 *  Created on: Aug 23, 2010
 *      Author: omidi
 */

#include "Decomposition.h"
#include "constants.h"

Decomposition::Decomposition() {
	// TODO Auto-generated constructor stub

}

Decomposition::~Decomposition() {
	// TODO Auto-generated destructor stub
}


std::vector<double>  Decomposition::QR_decomposition(std::vector<std::vector<double> > A, int n){
	std::vector<double> R(n, .0);
	std::vector<double> c(n, .0);
	double scale, sigma, tau;
	for(unsigned short int k = 0; k < n ; k++){
		scale = 0.0;
		for(unsigned short int i=k; i<n; i++) scale = FMAX(scale, ABS(A[i][k]));
		if(scale == 0.0){
			c[k] = R[k] = .0;
		} else {	// form !_k and Q_k.A
			for (int i = k; i < n; i++)	A[i][k] /= scale;
			double sum = .0;
			for(unsigned short int i = k; i < n; i++){ // calculating the norm of column 'k'
				sum += pow(A[i][k], 2.0);
			}
			sigma = SIGN(sqrt(sum), A[k][k]);
			A[k][k] += sigma;
			c[k] = sigma * A[k][k];
			R[k] = scale*sigma;
			for(unsigned short int j = k+1; j<n;j++){	// calculate Q_k.A
				double sum = .0;
				for(unsigned short int i = k; i<n; i++){
					sum += A[i][k]*A[i][j];
				}
				tau = sum/c[k];
				for(unsigned short int i = k; i<n; i++){
					A[i][j] -= tau*A[i][k];
				}
			}
		}
	}
	return R;
}


double Decomposition::determinant(std::vector<std::vector<double> > A, int n){ 	// From Lukas Burger's code

	//calculates the determinant in logscale, without the sign
	  //sign will be printed to STDERR
	  //by setting n, one can specify whether one wants to calculate
	  //the determinant of the whole matrix or a submatrix a[0....n][0....n]
	  vector<int> indx(n);
	  int i,j,k;
	  int imax = 0;
	  double big,dum,sum,temp;
	  vector<double> vv(n);
	  double d;
	  double det=0;
	  d=1.0;
	  for(i=0;i<n;i++){
	    big=0.0;
	    for(j=0;j<n;j++)
	      if((temp=fabs(A[i][j]))>big) big=temp;
	    if(big==0.0){cerr<<"Singular matrix in routine ludcmp"<<endl; exit(1);}
	    vv[i]=1.0/big;
	  }
	  for(j=0;j<n;j++){
	    for(i=0;i<j;i++){
	      sum=A[i][j];
	      for(k=0;k<i;k++) sum-=A[i][k]*A[k][j];
	      A[i][j]=sum;
	    }
	    big=0.0;
	    for(i=j;i<n;i++){
	      sum=A[i][j];
	      for(k=0;k<j;k++)
		sum-=A[i][k]*A[k][j];
	      A[i][j]=sum;
	      if((dum=vv[i]*fabs(sum))>=big){
		big=dum;
		imax=i;
	      }
	    }
	    if(j!=imax){
	      for(k=0;k<n;k++){
		dum=A[imax][k];
		A[imax][k]=A[j][k];
		A[j][k]=dum;
	      }
	      d=-d;
	      vv[imax]=vv[j];
	    }
	    indx[j]=imax;
	    if(A[j][j]==0.0){
	      A[j][j]=TINY;
	      cerr<<"watch out, matrix is singular\n";
	    }
	    if(j!=(n-1)){
	      dum=1.0/(A[j][j]);
	      for(i=j+1;i<n;i++) A[i][j]*=dum;
	    }
	  }
//		for(unsigned int i=0;i<=A.size()-1;i++){
//			for(unsigned int j=0;j<=i;j++){
//				cout << j << '\t' << i << '\t' << A[i][j] << endl;
//			}
//		}

	  double signum=d;//the sign of the determinant
	  //double min_diag=200;
	  for(int j=0;j<n;j++){
	    if(A[j][j]<0){
	      signum*=-1;
	    }
	    // if(fabs(a[j][j])<min_diag){
	    //  min_diag=fabs(a[j][j]);
	    //}
	    det+=log(fabs(A[j][j]));
	  }
	  //cerr<<"signum is "<<signum<<endl;
	  // cerr<<"det is "<<det<<endl;
	  //cerr<<signum;
	  // cout<<det<<" "<<signum<<" "<<min_diag<<endl;
	  return det;
}
