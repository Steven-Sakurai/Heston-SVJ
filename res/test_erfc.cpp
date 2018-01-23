#include <cmath>
#include <iostream>

using namespace std;

inline double normalCdf(double x) {
    return (1 - 0.5*erfc(x/sqrt(2)));
}

int main() {
    cout << normalCdf(-1) << '\t' << normalCdf(-0.5) << '\t' << normalCdf(1) << '\t' << normalCdf(2) << endl;
}
