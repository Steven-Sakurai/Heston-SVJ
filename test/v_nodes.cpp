#include "Coef.hpp"
#include <iostream>


int main() {
	stdVec par = {0.05, 2, 0.2, 0.25, -0.8};

    double mu = par[0];
    double kappa = par[1];
    double theta = par[2];
    double xi = par[3];
    double rho = par[4];

	double stLogMean = log(theta);
	double stLogSD = log(sqrt(theta*xi*xi/(2*kappa)));
	// create the nodes we want to interpolate
	stdVec interNodes = {-0.8, -0.75, -0.7, -0.65, -0.6, -0.55, -0.5, -0.45, -0.4, -0.35,
		 -0.3, -0.25, -0.2, -0.15, -0.1, -0.05, 0, 0.05, 0.1, 0.2, 0.3, 0.4};
	for(int i = 0; i < 22; ++i) 
		interNodes[i] = exp(stLogMean - interNodes[i]*stLogSD);
	for(int i = 0; i < interNodes.size(); ++i)
		cout << interNodes[i] << " " << endl;

	double v0 = interNodes[0], v1 = interNodes[1];
	double y0 = 4.6052, y1 = 4.6195;


}
