// [[Rcpp::depends(Rcpp)]]
// [[Rcpp::depends(BH)]]
#include <Rcpp.h>
using namespace Rcpp;
#include <boost/math/special_functions/gamma.hpp>
using boost::math::gamma_p;

#include <cmath>
#include <vector>
#include <algorithm>
#include <functional>
#include <numeric>

using namespace std;
typedef vector< double > stdVec;
typedef vector< vector< double > > stdMat;

static const double Pi = M_PI;

const double Delta = 1.0/252.0; // 1/252

const int nbeta = 3;

const int nbasis = 21;
const stdVec myNodes = {-0.8, -0.75, -0.7, -0.65, -0.6, -0.55, -0.5, -0.45, -0.4, -0.35, -0.3, -0.25, -0.2, -0.15, -0.1, -0.05, 0, 0.05, 0.1, 0.2, 0.3, 0.4};

void findNodes(const stdVec& par, stdVec& interNodes) {
	for(int i = 0; i < nbasis+1; ++i) 
		interNodes[i] = exp(log(par[2]) - myNodes[i]*log(sqrt(par[2]*par[3]*par[3]/(2*par[1])))); 
}


const int njumps = 2;
/*
	Transition density parameter
	0 stands for no jump
	1 stands for 1 jump
*/

/*
	mu, kappa, theta, xi, rho, lambda, l0,  l1,  mu_s, sigma_s^2
	0    1       2     3   4      5     6    7    8       9
*/

inline double tran_y_mean_0(const stdVec& par, double y0, double v0) {
	return y0 + (par[0] - (0.5)*v0) * Delta;
}

inline double tran_y_var_0(const stdVec& par, double v0) {
	return v0*Delta;
}

inline double tran_y_mean_n(const stdVec& par, double y0, double v0, int n) {
	return y0 + (par[0] - (0.5)*v0)*Delta + par[7]*n;
}

inline double tran_y_var_n(const stdVec& par, double v0, int n) {
	return v0*Delta + par[8]*n;
}


inline double tran_v_mean(const stdVec& par, double v0) {
	return v0 + par[1]*(par[2] - v0)*Delta;
}

inline double tran_v_var(const stdVec& par, double v0) {
	return par[3]*par[3]*v0*Delta;
}

/*
	transition density decomposition: p(y, v) = p(y) * p(v|y) 
*/

inline double cond_v_var(const stdVec& par, double v0)
{
	return par[3]*par[3]*v0*Delta * (1.0 - pow(par[4], 2));
}

inline double cond_v_mean_0(const stdVec& par, double y, double y0, double v0)
{
	return tran_v_mean(par, v0) + 
		sqrt((par[3]*par[3]*v0*Delta)/(v0*Delta))*par[4]*(y - tran_y_mean_0(par, y0, v0));
}

inline double cond_v_mean_n(const stdVec& par, double y, double y0, double v0, int n)
{
	return tran_v_mean(par, v0) +
		sqrt((par[3]*par[3]*v0*Delta)/(v0*Delta + par[8]*n))*par[4]*(y - tran_y_mean_n(par, y0, v0, n));
}


/*

*/


inline double prePy_0(const stdVec& par, double y, double y0, double v0) {
	double mu0 = tran_y_mean_0(par, y0, v0);
	double var0 = tran_y_var_0(par, v0);
	return exp(-(y-mu0)*(y-mu0) / (2.0*var0)) / sqrt(2*Pi*var0);
}

inline double prePy_n(const stdVec& par, double y, double y0, double v0, int n) {
	double mu1 = tran_y_mean_n(par, y0, v0, n);
	double var1 = tran_y_var_n(par, v0, n);
	return exp(-(y- mu1)*(y-mu1)/(2.0*var1)) / sqrt(2*Pi*var1);
}

// p(v|y) conditional distribution
inline double pdf_c_0(const stdVec& par, double y, double y0, double v, double v0) {
	double mu = cond_v_mean_0(par, y, y0, v0);
	double var = cond_v_var(par, v0);
	return ( exp(-(v-mu)*(v-mu) / (2.0*var) ) / sqrt(2*Pi*var) );
}

inline double pdf_c_n(const stdVec& par, double y, double y0, double v, double v0, int n) {
	double mu = cond_v_mean_n(par, y, y0, v0, n);
	double var = cond_v_var(par, v0);
	return ( exp(-(v-mu)*(v-mu) / (2.0*var) ) / sqrt(2*Pi*var) );
}

inline double cdf_c_0(const stdVec& par, double y, double y0, double v, double v0) {
	double mu = cond_v_mean_0(par, y, y0, v0);
	double var = cond_v_var(par, v0);
	return 1.0 - 0.5 * erfc((v-mu)/sqrt(2.0*var));
}

inline double cdf_c_n(const stdVec& par, double y, double y0, double v, double v0, int n) {
	double mu = cond_v_mean_n(par, y, y0, v0, n);
	double var = cond_v_var(par, v0);
	return 1.0 - 0.5 * erfc((v-mu)/sqrt(2.0*var));
}

inline double prePy(const stdVec& par, double y, double y0, double v0) {
	double ret = prePy_0(par, y, y0, v0);
	double fact = (par[5] + par[6]*v0)*Delta;
	for(int i = 1; i <= njumps; ++i) {
		ret += prePy_n(par, y, y0, v0, i) * fact;
		fact *= (par[5] + par[6]*v0)*Delta / (double) i;
	}
	return exp(-(par[5] + par[6]*v0)*Delta) * ret;
}

// lagrange cubic interpolation, coefficients stored in coef
void interCoef(double f[4], double x, double interval, stdVec& coef) {
	double f1 = f[0];
	double f2 = f[1];
	double f3 = f[2];
	double f4 = f[3];

	coef[0] = f1 + ((11*f1 - 18*f2 + 9*f3 - 
		2*f4)*x)/(6.*interval) + ((f1 - (5*f2)/2. + 2*f3 - 
		f4/2.)*x*x)/(interval*interval) + ((f1 - 3*f2 + 
		3*f3 - f4)*(x*x*x))/(6.*(interval*interval*interval));

	coef[1] = (-11*f1 + 18*f2 - 9*f3 +
	 2*f4)/(6.*interval) + ((-2*f1 + 5*f2 - 
	 	4*f3 + f4)*x)/(interval*interval) + ((-f1 + 3*f2 - 
	 	3*f3 + f4)*x*x)/(2.*(interval*interval*interval));

	coef[2] = (f1 - (5*f2)/2. + 2*f3 - 
		f4/2.)/(interval*interval) + ((f1 - 3*f2 + 3*f3 - 
			f4)*x)/(2.*(interval*interval*interval));
	
	coef[3] = (-f1 + 3*f2 - 3*f3 + f4)/(6.*(interval*interval*interval));
}

void get4Integral(const stdVec& par, double y, double y0, double v1, double v2, double v0, double *res0, double *res1, double *res2, double *res3) {
	double t0 = 0, t1 = 0, t2 = 0, t3 = 0;

	double mu_c = cond_v_mean_0(par, y, y0, v0);
	double var_c = cond_v_var(par, v0);
	double theprePy = prePy_0(par, y, y0, v0);
	double thepdf_c2 = pdf_c_0(par, y, y0, v2, v0);
	double thepdf_c1 = pdf_c_0(par, y, y0, v1, v0);
	double thecdf_c2 = cdf_c_0(par, y, y0, v2, v0);
	double thecdf_c1 = cdf_c_0(par, y, y0, v1, v0);
	t0 = theprePy * (thecdf_c2 - thecdf_c1);
	t1 = t0*mu_c - theprePy*var_c*(thepdf_c2 - thepdf_c1);
	t2 = t0*var_c + t1*mu_c - theprePy*var_c*(thepdf_c2*v2 - thepdf_c1*v1);
	t3 = 2*t1*var_c + t2*mu_c - theprePy*var_c*(thepdf_c2*v2*v2 - thepdf_c1*v1*v1);

	double fact = (par[5] + par[6]*v0)*Delta;
	double norm = 1;
	for(int i = 1; i <= njumps; ++i) {
		double mu_c_n = cond_v_mean_n(par, y, y0, v0, i);
		double var_c_n = cond_v_var(par, v0);
		double theprePy_n = prePy_n(par, y, y0, v0, i);
		double thepdf_c2_n = pdf_c_n(par, y, y0, v2, v0, i);
		double thepdf_c1_n = pdf_c_n(par, y, y0, v1, v0, i);
		double thecdf_c2_n = cdf_c_n(par, y, y0, v2, v0, i);
		double thecdf_c1_n = cdf_c_n(par, y, y0, v1, v0, i);

		double plus0 = theprePy_n * (thecdf_c2_n - thecdf_c1_n);
		double plus1 = plus0*mu_c_n - theprePy_n*var_c_n*(thepdf_c2_n - thepdf_c1_n);
		double plus2 = plus0*var_c_n + plus1*mu_c_n - theprePy_n*var_c_n*(thepdf_c2_n*v2 - thepdf_c1_n*v1);
		double plus3 = 2*plus1*var_c_n + plus2*mu_c_n - theprePy_n*var_c_n*(thepdf_c2_n*v2*v2 - thepdf_c1_n*v1*v1);

		norm += fact;
		t0 += plus0 * fact;
		t1 += plus1 * fact;
		t2 += plus2 * fact;
		t3 += plus3 * fact;
		fact *= (par[5] + par[6]*v0)*Delta / (double) i;
	}
	t0 *= 1.0/norm;
	t1 *= 1.0/norm;
	t2 *= 1.0/norm;
	t3 *= 1.0/norm;
	*res0 = t0, *res1 = t1, *res2 = t2, *res3 = t3;
}


// find the coef vector alpha, expansion of marginal transition density
// and the coef matrix beta, expansion of B_k
void BasisCoef(stdVec& alpha, stdMat& beta, double y, double y0, const stdVec& par, const stdVec& interNodes) {
	// calculate alpha first
	for(int k = 0; k < nbasis; ++k) {
		double v1 = interNodes[k];
		double v2 = interNodes[k+1];

		stdVec tmp(4);
		double stepSize = (v2 - v1)/3.0;
			
		double p1_y[4] = {0};
		p1_y[0] = prePy(par, y, y0, v1);
		p1_y[1] = prePy(par, y, y0, v1+stepSize);
		p1_y[2] = prePy(par, y, y0, v1+stepSize*2);
		p1_y[3] = prePy(par, y, y0, v2);

		interCoef(p1_y, interNodes[k], stepSize, tmp);
		alpha[k] = tmp[0];
		alpha[k+nbasis] = tmp[1];
		alpha[k+2*nbasis] = tmp[2];
		alpha[k+3*nbasis] = tmp[3];

		double I0[4] = {0};
		double I1[4] = {0};
		double I2[4] = {0};
		double I3[4] = {0};

		for(int j = 0; j < nbasis; ++j) {
			if(abs(k-j) > nbeta)
				continue;
			stepSize = (interNodes[j+1] - interNodes[j])/3.0;
			
			for(int index = 0; index < 4; ++index) {
				double v0 = interNodes[j]+stepSize*index;
				get4Integral(par, y, y0, v1, v2, v0, &I0[index], &I1[index], &I2[index], &I3[index]);
			}

			interCoef(I0, interNodes[j], stepSize, tmp);
			beta[k][j] = tmp[0];          beta[k][j+nbasis] = tmp[1];
			beta[k][j+2*nbasis] = tmp[2]; beta[k][j+3*nbasis] = tmp[3];

			interCoef(I1, interNodes[j], stepSize, tmp);
			beta[k+nbasis][j] = tmp[0];          beta[k+nbasis][j+nbasis] = tmp[1];
			beta[k+nbasis][j+2*nbasis] = tmp[2]; beta[k+nbasis][j+3*nbasis] = tmp[3];

			interCoef(I2, interNodes[j], stepSize, tmp);
			beta[k+2*nbasis][j] = tmp[0];          beta[k+2*nbasis][j+nbasis] = tmp[1];
			beta[k+2*nbasis][j+2*nbasis] = tmp[2]; beta[k+2*nbasis][j+3*nbasis] = tmp[3];

			interCoef(I3, interNodes[j], stepSize, tmp);
			beta[k+3*nbasis][j] = tmp[0];          beta[k+3*nbasis][j+nbasis] = tmp[1];
			beta[k+3*nbasis][j+2*nbasis] = tmp[2]; beta[k+3*nbasis][j+3*nbasis] = tmp[3];
		}
	}			
}

// calculate the initial moments

void initialMoments(const stdVec& par, stdVec& initialMoments, const stdVec& nodes) {
	// create vector to hold the scaled interval nodes
	stdVec interNodes(nbasis+1);
	// scale and shape parameter of gamma distribution
	double beta = 0.5*par[3]*par[3]/par[1];
	double shape = 2*par[1]*par[2]/pow(par[3], 2);
	// scale the nodes by 1/scale
	transform(nodes.begin(), nodes.end(), interNodes.begin(), [=](double x) { return x/beta; });
	
	//#pragma omp parallel for num_threads(1)
		for(int i = 0; i < nbasis; ++i) {
			initialMoments[i] = gamma_p(shape, interNodes[i+1]) - gamma_p(shape, interNodes[i]);
			initialMoments[nbasis+i] = shape*beta*(gamma_p(shape+1, interNodes[i+1]) - gamma_p(shape+1, interNodes[i]));
			initialMoments[2*nbasis+i] = shape*(shape+1)*beta*beta * (gamma_p(shape+2, interNodes[i+1]) - gamma_p(shape+2, interNodes[i]));
			initialMoments[3*nbasis+i] = shape*(shape+1)*(shape+2)*beta*beta*beta * (gamma_p(shape+3, interNodes[i+1]) - gamma_p(shape+3, interNodes[i]));
		}
}

/*
	negative log-likelihood function
	objective
*/

// [[Rcpp::export]]
double nll(NumericVector par_, NumericVector y_) {
	// convert armadillo vector to std vector
	stdVec par = as< stdVec >(par_);
	stdVec y = as< stdVec >(y_);
	int N = y.size() - 1;

	// calculate the interval
	stdVec nodes(nbasis+1);
	findNodes(par, nodes);

	// calculate the initial moments
	stdVec initm(4*nbasis);
	initialMoments(par, initm, nodes);
	stdMat filterMoment(N+1, stdVec(4*nbasis)); 
	filterMoment[0] = initm;
	
	stdVec Li(N);
	stdVec alpha(4*nbasis, 0);
	stdMat beta(4*nbasis, stdVec(4*nbasis, 0));

	int i, j;
	double ret = 0;
	// online induction
	for(i = 0; i < N; ++i) {
		BasisCoef(alpha, beta, y[i+1], y[i], par, nodes);
		Li[i] = inner_product(filterMoment[i].begin(), filterMoment[i].end(), alpha.begin(), 0.0);
		ret -= log(Li[i]);
		
		for(j = 0; j < 4*nbasis; ++j) { 
			filterMoment[i+1][j] = inner_product(beta[j].begin(), beta[j].end(), filterMoment[i].begin(), 0.0) / Li[i];
		}
	}
	return ret;
}
