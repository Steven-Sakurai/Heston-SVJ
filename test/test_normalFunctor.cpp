#include "Functor.hpp"

#include <iostream>

using namespace std;

int main() {
    double mu = 0.0;
    double var = 1.0;
    NormalPdf my_fun(mu, var);
    cout << my_fun(0) << endl;
}
