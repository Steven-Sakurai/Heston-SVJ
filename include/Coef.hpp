#pragma once
#include <cmath>
#include <algorithm>
#include <functional>

using namespace std;
#include "Functor.hpp"

// Originally struct intercoefficient InterCoef_cpp()
// this part is just copy and paste since there's no elegant way to do this...
// test this function in `test_coef.cpp`
stdVec interCoef(function<double(double)> f, double x, double interval) {
	// use stdVec to store coef, I think it's better than using a struct...
	stdVec ret(4);
	double f1 = f(x);
	double f2 = f(x+ interval);
	double f3 = f(x+ 2*interval);
	double f4 = f(x+ 3*interval);

	ret[0] = f1 + ((11*f1 - 18*f2 + 9*f3 - 
		2*f4)*x)/(6.*interval) + ((f1 - (5*f2)/2. + 2*f3 - 
		f4/2.)*pow(x,2))/pow(interval,2) + ((f1 - 3*f2 + 
		3*f3 - f4)*pow(x,3))/(6.*pow(interval,3));

	ret[1] = (-11*f1 + 18*f2 - 9*f3 +
	 2*f4)/(6.*interval) + ((-2*f1 + 5*f2 - 
	 	4*f3 + f4)*x)/pow(interval,2) + ((-f1 + 3*f2 - 
	 	3*f3 + f4)*pow(x,2))/(2.*pow(interval,3));

	ret[2] = (f1 - (5*f2)/2. + 2*f3 - 
		f4/2.)/pow(interval,2) + ((f1 - 3*f2 + 3*f3 - 
			f4)*x)/(2.*pow(interval,3));
	
	ret[3] = (-f1 + 3*f2 - 3*f3 + f4)/(6.*pow(interval,3));
	return ret;
}

/*
	Main Dish!
		Arguments:
		1. marginal transition expansion coef \alpha_{k, i}^{(l)}, vector of length 4*nbasis
		2. marginal moment expansion coef \beta_{k, i}^{(l)}, square matrix of size (4*nbasis)^2
		calculate 1 & 2 using
		3. (i-1) th y0
		4. i th y
		5. parameters
*/

void findNodes(const stdVec& par, stdVec& interNodes) {
	// create the (nbasis+1) nodes we want to interpolate
	interNodes = {-0.8, -0.75, -0.7, -0.65, -0.6, -0.55, -0.5, -0.45, -0.4, -0.35,
		 -0.3, -0.25, -0.2, -0.15, -0.1, -0.05, 0, 0.05, 0.1, 0.2, 0.3, 0.4};
	int i = 0;
	for(; i < 22; ++i) 
		interNodes[i] = exp(log(par[2]) - interNodes[i]*log(sqrt( par[2]*par[3]*par[3] / (2*par[1]) ))); 	
}


void BasisCoef(stdVec& alpha, stdMat& beta, double y, double y0, const stdVec& par, const stdVec& interNodes) {
	// determine the interpolation nodes, from the stationary gamma distribution
	// same as method in InitialMoments.hpp
	// first calculate alpha easily 
	static auto p1_y = bind(prePy, par, y, y0, placeholders::_1); // the marginal density functional object, i.e. prestfun
	static double stepSize;
	static stdVec tmp(4);
	int i = 0;
	for(; i < nbasis; ++i) {
		stepSize = (interNodes[i+1] - interNodes[i])/3.0;
		tmp = interCoef(p1_y, interNodes[i], stepSize);
		alpha[i] = tmp[0];
		alpha[i+nbasis] = tmp[1];
		alpha[i+2*nbasis] = tmp[2];
		alpha[i+3*nbasis] = tmp[3];
	}

	int k = 0, j;
	// traverse v first then v0, row first
	for(; k < nbasis; ++k) {
		// initialize the four functor
		auto I0 = bind(Integ0, par, y, y0, interNodes[k], interNodes[k+1], placeholders::_1);
		auto I1 = bind(Integ1, par, y, y0, interNodes[k], interNodes[k+1], placeholders::_1);
		auto I2 = bind(Integ2, par, y, y0, interNodes[k], interNodes[k+1], placeholders::_1);
		auto I3 = bind(Integ3, par, y, y0, interNodes[k], interNodes[k+1], placeholders::_1);

		for(j = max(k-3, 0); j < min(k+3, nbasis); ++j) {
			//if(abs(k-j) > 3)
			//	continue;
			stepSize = (interNodes[j+1] - interNodes[j])/3.0;
			
			tmp = interCoef(I0, interNodes[j], stepSize);
			beta[k][j] = tmp[0];          beta[k][j+nbasis] = tmp[1];
			beta[k][j+2*nbasis] = tmp[2]; beta[k][j+3*nbasis] = tmp[3];

			tmp = interCoef(I1, interNodes[j], stepSize);
			beta[k+nbasis][j] = tmp[0];          beta[k+nbasis][j+nbasis] = tmp[1];
			beta[k+nbasis][j+2*nbasis] = tmp[2]; beta[k+nbasis][j+3*nbasis] = tmp[3];

			tmp = interCoef(I2, interNodes[j], stepSize);
			beta[k+2*nbasis][j] = tmp[0];          beta[k+2*nbasis][j+nbasis] = tmp[1];
			beta[k+2*nbasis][j+2*nbasis] = tmp[2]; beta[k+2*nbasis][j+3*nbasis] = tmp[3];

			tmp = interCoef(I3, interNodes[j], stepSize);
			beta[k+3*nbasis][j] = tmp[0];          beta[k+3*nbasis][j+nbasis] = tmp[1];
			beta[k+3*nbasis][j+2*nbasis] = tmp[2]; beta[k+3*nbasis][j+3*nbasis] = tmp[3];
		}
	}
}





