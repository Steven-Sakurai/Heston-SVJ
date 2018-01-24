#ifndef __PAY_OFF_HPP
#define __PAY_OFF_HPP

class PayOff {
public:
    PayOff() {}// Default (no parameter) constructor
    virtual ~PayOff() {}; // Virtual destructor
    // Overloaded () operator, turns the PayOff into an abstract function object
    virtual double operator() (const double& S) const = 0;
};

class PayOffCall : public PayOff {
private:
	double K;
public:
	PayOffCall(const double& K_) {
		K = K_;
	}
    virtual ~PayOffCall() {}
    virtual double operator() (const double& S) const {
        return S > K ? (S-K) : 0;
    }
};

class PayOffPut : public PayOff {
private:
	double K;
public:
	PayOffPut(const double& K_) {
		K = K_;
	}
    virtual ~PayOffPut() {}
    virtual double operator() (const double& S) const {
        return S < K ? (K-S) : 0;
    }
};

#endif
