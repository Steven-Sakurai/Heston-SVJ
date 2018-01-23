#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>

#include <boost/math/special_functions/gamma.hpp>

using boost::math::gamma_p;
using namespace std;

typedef vector< double > stdVec;
typedef vector< vector<double> > stdMat;

int main() {
    double mu = 0.05;
    double kappa = 2;
    double theta = 0.2;
    double xi = 0.25;
    double rho = -0.8;
	
	// scale parameter of gamma dis
	double beta = xi*xi/(2*kappa);
	
	// taking log of the mean and the standard deviation of gamma dis
	double stLogMean = log(theta);
	double stLogSD = log(sqrt(theta*xi*xi/(2*kappa)));
	
	stdVec interNodes = {-0.8, -0.75, -0.7, -0.65, -0.6, -0.55, -0.5, -0.45, -0.4, -0.35,
		 -0.3, -0.25, -0.2, -0.15, -0.1, -0.05, 0, 0.05, 0.1, 0.2, 0.3, 0.4};
	
	for(int i = 0; i < 22; ++i) 
		interNodes[i] = exp(stLogMean - interNodes[i]*stLogSD) / (xi*xi/(2*kappa));
	
	double shape = 2*kappa*theta/(xi*xi);
	
	stdMat initialMoments(4, stdVec(21, 0));
	
	for(int i = 0; i < 21; ++i)
		initialMoments[0][i] = gamma_p(shape, interNodes[i+1]) - gamma_p(shape, interNodes[i]);
	for(int i = 0; i < 21; ++i)
		initialMoments[1][i] = shape*beta*(gamma_p(shape+1, interNodes[i+1]) - gamma_p(shape+1, interNodes[i]));
	for(int i = 0; i < 21; ++i)
		initialMoments[2][i] = shape*(shape+1)*beta*beta * (gamma_p(shape+2, interNodes[i+1]) - gamma_p(shape+2, interNodes[i]));
	for(int i = 0; i < 21; ++i)
		initialMoments[3][i] = shape*(shape+1)*(shape+2)*beta*beta*beta * (gamma_p(shape+3, interNodes[i+1]) - gamma_p(shape+3, interNodes[i]));
	
	cout << "0th order basis:" << setprecision(16) << endl;
	for(int i = 0; i < 21; ++i) {
		cout << initialMoments[0][i] << endl;
	}
	
	cout << "1th order basis:" << endl;
	for(int i = 0; i < 21; ++i) {
		cout << initialMoments[1][i] << endl;
	}
		
	cout << "2th order basis:" << endl;
	for(int i = 0; i < 21; ++i) {
		cout << initialMoments[2][i] << endl;
	}
	
	cout << "3th order basis:" << endl;
	for(int i = 0; i < 21; ++i) {
		cout << initialMoments[3][i] << endl;
	}
	
	
	/*
    cout << "test incomplete gamma lower function in boost: " << endl;
	double z1 = 2*kappa*theta/(xi*xi);
	double a1 = 0.154221*10*kappa*theta/(xi*xi)*pow(1/kappa*theta*xi*xi, 0.375);
	double z2 = 2*kappa*theta/(xi*xi);
	double a2 = 0.151572*10*kappa*theta/(xi*xi)*pow(1/kappa*theta*xi*xi, 0.4);
    cout << std::setprecision(16) << boost::math::gamma_p(z1, a1) - boost::math::gamma_p(z2, a2) << endl;
	*/
	return 0;
}