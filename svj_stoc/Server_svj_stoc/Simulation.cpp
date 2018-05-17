// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
using namespace Rcpp;

#include <vector>
#include <random>
#include <cmath>

using namespace std;
typedef vector< double > stdVec;

// [[Rcpp::export]]
NumericVector SimSVJ(NumericVector par, int len, double y0, double v0, double Delta, unsigned int seed1, unsigned int seed2, unsigned int seed3, unsigned int seed4) {
	NumericVector y(len);
	stdVec v(len);

	default_random_engine gen1{seed1};
	default_random_engine gen2{seed2};
	default_random_engine gen3{seed3};
	default_random_engine gen4{seed4};

	double mu = par[0], kappa = par[1], theta = par[2], xi = par[3], rho = par[4], 
		l0 = par[5], l1 = par[6], mu_j = par[7], sd_j = sqrt(par[8]);

	normal_distribution<double> norm{0,1};
        
    y[0] = y0; v[0] = v0;
	for(int i = 1; i < len; ++i) {
		double rnorm1 = norm(gen1);
		double rnorm2 = rho*rnorm1 + sqrt(1-rho*rho)*norm(gen2);

		v[i] = v[i-1] + kappa*(theta - v[i-1])*Delta + xi*sqrt(v[i-1]*Delta)*rnorm1 + 0.25*xi*xi*Delta*(rnorm1*rnorm1 - 1);
		if(v[i] < 0) v[i] = -v[i];
			
		y[i] = y[i-1] + (mu - (0.5)*v[i-1])*Delta + sqrt(v[i-1]*Delta)*rnorm2;    
        // add jumps
        poisson_distribution<int> pois((l0 + l1*v[i-1])*Delta);
        int njumps = pois(gen3);
        while(njumps != 0) {
            y[i] += norm(gen4)*sd_j + mu_j;
            --njumps;
        }
	}
	return y;
}