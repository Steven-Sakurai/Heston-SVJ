/*
	need to get Yobs, i.e. simulation
*/
//#include <armadillo>
#include <iostream>
#include <numeric>

#include "InitialMoments.hpp"
#include "Coef.hpp"

// N stands for the number of obs of x = log(S)
double loglike(const stdVec& par, const stdVec& y) {
	// obs
	int N  = y.size();
	// M_{k, i}^{(l)}
	stdMat filterMoment(N+1, stdVec(4*nbasis)); 
	
	stdVec initm(4*nbasis);
	initialMoments(par, initm, false);
	filterMoment[0] = initm;
	
	stdVec Li(N);
	
	stdVec alpha(4*nbasis);
	stdMat beta(4*nbasis, stdVec(4*nbasis));

	int i, j;

	stdVec nodes(23);
	findNodes(par, nodes);
	for(i = 0; i < N; ++i) {
		BasisCoef(alpha, beta, y[i+1], y[i], par, nodes);
		Li[i] = inner_product(filterMoment[i].begin(), filterMoment[i].end(), alpha.begin(), 0.0);
		/*
		if(i == 1) {
			for(int l = 0; l < 4*nbasis; ++l)
				cout << beta[21][l] << endl;
		}
		*/
		//cout << Li[i] << endl;
		for(j = 0; j < 4*nbasis; ++j) { 
			filterMoment[i+1][j] = inner_product(beta[j].begin(), beta[j].end(), filterMoment[i].begin(), 0.0) / Li[i];
		}
	}
	//cout << "filterMoment at stage 1: " << endl;
	//for(int i = 0; i < 4*nbasis; ++i)
	//	cout << filterMoment[1][i] << endl;
	//cout << "1st and 2nd updated likelihood: " << endl;
	//cout << Li[0] << " " << Li[1] << endl;
	double ret = 0;
	for(int i = 0; i < N; ++i) {
		if(Li[i] <= 0)
			continue;
		ret -= log(Li[i]);
	}
	return ret;
}
