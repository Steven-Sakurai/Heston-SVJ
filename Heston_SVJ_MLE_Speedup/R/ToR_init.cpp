#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;

// [[Rcpp::depends(BH)]]
#include <boost/math/special_functions/gamma.hpp>
using boost::math::gamma_p;

// [[Rcpp::export]]
rowvec initialMoments(rowvec par) {
	rowvec initialMoments(4*21);
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
	rowvec interNodes = {-0.80, -0.75, -0.70, -0.65, -0.60, -0.55, -0.50, -0.45, -0.40, -0.35, -0.30, -0.25, -0.20, -0.15, -0.10, -0.05, 0, 0.05, 0.10, 0.20, 0.30, 0.40};
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
	return initialMoments;
}
