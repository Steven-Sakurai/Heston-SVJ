// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH)]]
#include <RcppArmadillo.h>
#include <boost/math/special_functions/gamma.hpp>
using namespace arma;
using boost::math::gamma_p;

#include <cmath>
#include <vector>
#include <algorithm>
#include <functional>
#include <numeric>



using namespace std;
typedef vector< double > stdVec;
typedef vector< vector< double > > stdMat;

#define Pi 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679

#define Delta 1.0/252.0
#define nbasis 21

inline double tran_y_mean(const stdVec& par, double& y0, double& v0) {
	return (y0 + (par[0]-0.5*v0)*Delta + par[5]*Delta*par[6]);
}

inline double tran_y_var(const stdVec& par, double& v0) {
	return (v0*Delta + par[5]*Delta*(par[7] + pow(par[6], 2)));
}

inline double tran_v_mean(const stdVec& par, double& v0) {
	return (v0 + par[1]*(par[2] - v0)*Delta);
}

inline double tran_v_var(const stdVec& par, double& v0) {
	return (pow(par[3], 2)*v0*Delta);
}

/*
	The one-step transition density decomposition: p(y, v) = p(y) * p(v|y) 
*/
inline double cond_v_var(const stdVec& par, double& v0)
{
	return (pow(par[3], 2)*v0*Delta * (1 - pow(par[4], 2)));
}

inline double cond_v_mean(const stdVec& par, double& y, double& y0, double& v0)
{
	return ((v0 + par[1]*(par[2] - v0)*Delta) + par[3]*par[4]*(y - (y0 + (par[0]-0.5*v0)*Delta + par[5]*Delta*par[6])));
}


inline double prePy(const stdVec& par, double y, double y0, double v0) {
	double a = pow(y - tran_y_mean(par, y0, v0), 2);
	double b = tran_y_var(par, v0);
	return exp(-a/(2*b)) / sqrt(2*Pi*b);
}

inline double pdf_c(const stdVec& par, double y, double y0, double v, double v0) {
	double mu = cond_v_mean(par, y, y0, v0);
	double var = cond_v_var(par, v0);
	return exp(-pow(v-mu, 2)/(2*var)) / sqrt(2*Pi*var);
}

inline double cdf_c(const stdVec& par, double y, double y0, double v, double v0) {
	double mu = cond_v_mean(par, y, y0, v0);
	double var = cond_v_var(par, v0);
	return 1 - 0.5*erfc((v-mu)/sqrt(var)/sqrt(2));
}

inline double Integ0(const stdVec& par, double y, double y0, double v1, double v2, double v0) {
	return prePy(par, y, y0, v0) * (cdf_c(par, y, y0, v2, v0) - cdf_c(par, y, y0, v1, v0));
}

inline double Integ1(const stdVec& par, double y, double y0, double v1, double v2, double v0) {
	double mu_c = cond_v_mean(par, y, y0, v0);
	double var_c = cond_v_var(par, v0);
	return Integ0(par, y, y0, v1, v2, v0)*mu_c  - 
		prePy(par, y, y0, v0)*var_c*(pdf_c(par, y, y0, v2, v0) - 
			pdf_c(par, y, y0, v1, v0));
}

inline double Integ2(const stdVec& par, double y, double y0, double v1, double v2, double v0) {
	double mu_c = cond_v_mean(par, y, y0, v0);
	double var_c = cond_v_var(par, v0);
	return Integ0(par, y, y0, v1, v2, v0)*var_c + Integ1(par, y, y0, v1, v2, v0)*mu_c - 
		prePy(par, y, y0, v0)*var_c*(pdf_c(par, y, y0, v2, v0)*v2 - pdf_c(par, y, y0, v1, v0)*v1);
}

inline double Integ3(const stdVec& par, double y, double y0, double v1, double v2, double v0) {
	double mu_c = cond_v_mean(par, y, y0, v0);
	double var_c = cond_v_var(par, v0);
	return 2*Integ1(par, y, y0, v1, v2, v0)*var_c + 
		Integ2(par, y, y0, v1, v2, v0)*mu_c - 
		prePy(par, y, y0, v0)*var_c*(pdf_c(par, y, y0, v2, v0)*v2*v2 - 
			pdf_c(par, y, y0, v1, v0)*v1*v1);
}

stdVec interCoef(function<double(double)> f, double x, double interval) {
	// use stdVec to store coef, I think it's better than using a struct...
	double f1 = f(x);
	double f2 = f(x+ interval);
	double f3 = f(x+ 2*interval);
	double f4 = f(x+ 3*interval);
	double c1 = ((11*f1 - 18*f2 + 9*f3 - 
		2*f4)*x)/(6.*interval);
	double c2 = ((f1 - (5*f2)/2. + 2*f3 - 
		f4/2.)*pow(x,2))/pow(interval,2);
	double c3 = ((f1 - 3*f2 + 3*f3 - 
			f4))/(6.*pow(interval,3));

	stdVec ret = {f1 + c1 + c2 + c3*pow(x, 3), 
				- c1/x - 2*c2/x - 3*c3*pow(x, 2),
				 c2/pow(x, 2) + 3*c3*x, 
				 -c3};
	return ret;
}

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

void initialMoments(const stdVec& par, stdVec& initialMoments) {
    double mu = par[0];
    double kappa = par[1];
    double theta = par[2];
    double xi = par[3];
    double rho = par[4];
	// scale parameter of gamma dis
	double beta = xi*xi/(2*kappa);
	// taking log of the mean and the standard deviation of gamma dis
	double stLogMean = log(theta);
	double stLogSD = log(sqrt(theta*xi*xi/(2*kappa)));
	// create the nodes we want to interpolate
	stdVec interNodes = {-0.8, -0.75, -0.7, -0.65, -0.6, -0.55, -0.5, -0.45, -0.4, -0.35,
		 -0.3, -0.25, -0.2, -0.15, -0.1, -0.05, 0, 0.05, 0.1, 0.2, 0.3, 0.4};
	for(int i = 0; i < 22; ++i) 
		interNodes[i] = exp(stLogMean - interNodes[i]*stLogSD) / (xi*xi/(2*kappa));
	// shape parameter of gamma dis
	double shape = 2*kappa*theta/(xi*xi);
	// calculate using the beautiful property of gamma distribution	
	for(int i = 0; i < 21; ++i)
		initialMoments[i] = gamma_p(shape, interNodes[i+1]) - gamma_p(shape, interNodes[i]);
	for(int i = 0; i < 21; ++i)
		initialMoments[21+i] = shape*beta*(gamma_p(shape+1, interNodes[i+1]) - gamma_p(shape+1, interNodes[i]));
	for(int i = 0; i < 21; ++i)
		initialMoments[42+i] = shape*(shape+1)*beta*beta * (gamma_p(shape+2, interNodes[i+1]) - gamma_p(shape+2, interNodes[i]));
	for(int i = 0; i < 21; ++i)
		initialMoments[63+i] = shape*(shape+1)*(shape+2)*beta*beta*beta * (gamma_p(shape+3, interNodes[i+1]) - gamma_p(shape+3, interNodes[i]));
}

// [[Rcpp::export]]
double loglike(const rowvec& par_, const rowvec& y_) {
	// obs
	stdVec par = conv_to< stdVec >::from(par_);
	stdVec y = conv_to< stdVec >::from(y_);
	int N  = y.size();
	// M_{k, i}^{(l)}
	stdMat filterMoment(N+1, stdVec(4*nbasis)); 
	
	stdVec initm(4*nbasis);
	initialMoments(par, initm);
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
		for(j = 0; j < 4*nbasis; ++j) { 
			filterMoment[i+1][j] = inner_product(beta[j].begin(), beta[j].end(), filterMoment[i].begin(), 0.0) / Li[i];
		}
	}
	double ret = 0;
	for(int i = 0; i < N; ++i) {
		if(Li[i] <= 0)
			continue;
		ret -= log(Li[i]);
	}
	return ret;
}