#include "Coef.hpp"

#include <iostream>
#include <iomanip>

int main() {
    
	stdVec alpha(84);
	stdMat beta(84, stdVec(84));
    double y0 = 4.6052, y1 = 4.6195;

    int N = 1000;
    for(int j = 0; j < N; ++j)
	    BasisCoef(alpha, beta, y0, y1);
	
    cout << "alpha:" << setprecision(16) << setw(16) << endl;
    for(int i = 0; i < 84; ++i) {
        cout << alpha[i] << '\t';
        if((i+1) % 4 == 0)
            cout << endl;
    }
    cout << "length: " << alpha.size() << endl;
	/*
	int i, j;
	double stLogMean = log(par[2]);
	double stLogSD = log(sqrt( par[2]*par[3]*par[3] / (2*par[1]) ));
	stdVec interNodes = {-0.8, -0.75, -0.7, -0.65, -0.6, -0.55, -0.5, -0.45, -0.4, -0.35,
		 -0.3, -0.25, -0.2, -0.15, -0.1, -0.05, 0, 0.05, 0.1, 0.2, 0.3, 0.4};
	for(i = 0; i < 22; ++i) 
		interNodes[i] = exp(stLogMean - interNodes[i]*stLogSD) / (par[3]*par[3]/(2*par[1]));
	double y0 = 4.6052, y = 4.6195;
	i = 0, j = 0;
	Integral0 I0(y, y0, interNodes[i], interNodes[i+1]);
	double stepSize = (interNodes[j+1] - interNodes[j])/3.0;
	stdVec tmp = interCoef(I0, interNodes[j], stepSize);
	for(i = 0; i < 4; ++i)
		cout << setprecision(16) << setw(16) << I0(interNodes[i]) << endl;
	*/
	
	cout << "beta:" << endl;
	for(int i = 0; i < 84; ++i) {
		for(int j = 0; j < 84; ++j) {
				cout << beta[i][j] << '\t';
		}
		cout << endl;
    }
	
}
