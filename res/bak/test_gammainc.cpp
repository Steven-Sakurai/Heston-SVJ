#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>

#include <boost/math/special_functions/gamma.hpp>

using boost::math;

#include "Parameter.hpp"

using namespace std;

typedef vector< double > stdVec;

int main() {
    double mu = 0.05;
    double kappa = 2;
    double theta = 0.2;
    double xi = 0.25;
    double rho = -0.8;

    cout << "test incomplete gamma lower function in boost: " << endl;
	double z1 = 2*kappa*theta/(xi*xi);
	double a1 = 0.154221*10*kappa*theta/(xi*xi)*pow(1/kappa*theta*xi*xi, 0.375);
	double z2 = 2*kappa*theta/(xi*xi);
	double a2 = 0.151572*10*kappa*theta/(xi*xi)*pow(1/kappa*theta*xi*xi, 0.4);
    cout << std::setprecision(16) << boost::math::gamma_p(z1, a1) - boost::math::gamma_p(z2, a2) << endl;
	return 0;
}

