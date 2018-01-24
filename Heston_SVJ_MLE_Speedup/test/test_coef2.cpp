#include "Coef.hpp"

#include <iostream>
#include <iomanip>

int main() {
	stdVec par = {0.05, 2.0, 0.2, 0.25, -0.8};
    stdVec alpha(84);
	//stdMat beta(84, stdVec(84));
    double y0 = 4.6052, y1 = 4.6195;
    BasisCoef(alpha, y1, y0, par);
    cout << setprecision(16) << setw(16) << endl;
    for(int i = 0; i < 84; ++i) {
        cout << alpha[i] << '\t';
        if((i+1) % 4 == 0)
            cout << endl;
    }
    cout << alpha.size() << endl;
}
