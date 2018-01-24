#ifndef _CORSTDNORM
#define _CORSTDNORM

#include <armadillo>
#include <iostream>
#include <vector>

using namespace std;
using namespace arma;

typedef vector<double> stdvec;
typedef std::vector< std::vector<double> > stdvecvec;

stdvecvec mvrnormArma(int nrows, stdvec std_mu, stdvecvec std_sigma) {
	arma_rng::set_seed_random();
	int ncols = std_mu.size(); 
	vec mu = conv_to< vec >::from(std_mu);
	mat sigma = zeros<mat>(nrows, ncols);
	for(int i = 0; i < nrows; ++i)
		sigma.row(i) = conv_to< rowvec >::from(std_sigma[i]);
	stdvecvec ret(nrows);
    mat Y = randn(nrows, ncols); 
    Y = repmat(mu, 1, nrows).t() + Y * arma::chol(sigma);
	for(int i = 0; i < nrows; ++i) {
		ret[i] = conv_to< stdvec >::from(Y.row(i));
	}
	return ret;
}

mat mvrnormArma(int nrows, vec mu, mat sigma) {
	arma_rng::set_seed_random();
	int ncols = sigma.n_cols;
	mat Y = randn(nrows, ncols);
    return repmat(mu, 1, nrows).t() + Y * arma::chol(sigma);
}
#endif