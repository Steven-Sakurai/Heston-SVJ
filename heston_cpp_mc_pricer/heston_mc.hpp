#ifndef __HESTON_MC_H
#define __HESTON_MC_H

#include <cmath>
#include <vector>
#include "option.hpp"

typedef std::vector<double> stdvec;

// The HestonEuler class stores the necessary information
class HestonEuler {
private:
 	Option* pOption; // K, r, T
 	double kappa;
 	double theta;
 	double xi;
 	double rho;

public:
 	HestonEuler(Option* _pOption, double _kappa, double _theta, double _xi, double _rho):
		pOption(_pOption), kappa(_kappa), theta(_theta), xi(_xi), rho(_rho) {}
	virtual ~HestonEuler() {}

	// Calculate the volatility path
	// vol_draw is one part of the correlated B.M.
	void calc_vol_path(const stdvec& vol_draws, stdvec& vol_path) {
	    size_t vec_size = vol_draws.size();
	    double dt = (double) pOption->T/vec_size;
	    // Iterate through the correlated random draws vector and
	    // 'Full Truncation'
	    for (int i=1; i < vec_size; i++) {
	      double v_max = vol_path[i-1] >= 0.0 ? vol_path[i-1] : 0;
	      vol_path[i] = vol_path[i-1] + kappa * dt * (theta - v_max) + xi * sqrt(v_max * dt) * vol_draws[i-1];
	    }
	}

	// Calculate the asset price path
	// spot_draws is the other part of the correlated B.M.
 	void calc_spot_path(const stdvec& spot_draws, const stdvec& vol_path, stdvec& spot_path) {
		size_t vec_size = spot_draws.size();
		double dt = (double) pOption->T / vec_size;
		
		// 'Full Truncation'
		for (int i=1; i<vec_size; i++) {
			double v_max = vol_path[i-1] >= 0.0 ? vol_path[i-1] : 0;
			spot_path[i] = spot_path[i-1] * exp( (pOption->r - 0.5*v_max)*dt + sqrt(v_max*dt)*spot_draws[i-1]);
		}
 	}
};

#endif
