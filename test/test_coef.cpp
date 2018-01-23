#include "Coef.hpp"
#include <iostream>
#include <cmath>

class MyFun: public Functor {
public:
    MyFun() {}
    ~MyFun() {}
    double operator() (double x) const { 
        return 7 - x + 4*pow(x, 2) + 3.1*pow(x, 3);
    }
};

int main() {
    NormalPdf my_fun(5.0, 1.0);
    cout << my_fun(5.0) << endl;
    MyFun myfun;
    stdVec coef = interCoef(myfun, -2, 1);
    for(int i = 0; i < 4; ++i)
        cout << coef[i] << endl;
    return 0;
}
