#ifndef __OPTION_H
#define __OPTION_H

#include "payoff.hpp"

class Option {
public:
    PayOff* pay_off;
    double K;
    double r;
    double T;

    Option(double _K, double _r, double _T, PayOff* _pay_off):
        K(_K), r(_r), T(_T), pay_off(_pay_off) {}
    virtual ~Option() {}
};

#endif
