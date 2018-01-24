#include <iostream>

#include "payoff.hpp"
#include "option.hpp"
#include "CorStdNorm.hpp"
#include "heston_mc.hpp"


void generate_normal_correlation_paths(double rho, stdvec& spot_normals, stdvec& cor_normals, int N) {
	vec stdmu = {0, 0};
	mat stdsigma = {{1, rho}, {rho, 1}};
	mat tmp_res = mvrnormArma(N, stdmu, stdsigma);
	spot_normals = conv_to< stdvec >::from(tmp_res.col(0));
	cor_normals = conv_to< stdvec >::from(tmp_res.col(1));
}

int main(int argc, char **argv) {

    unsigned num_sims = 100000;   // Number of simulated asset paths
    unsigned num_intervals = 1000;  // Number of intervals for the asset path to be sampled 

    double S_0 = 100.0;    // Initial spot price
    double K = 100.0;      // Strike price
    double r = 0.0319;     // Risk-free rate
    double v_0 = 0.010201; // Initial volatility 
    double T = 1.00;       // One year until expiry

    double rho = -0.7;     // Correlation of asset and volatility
    double kappa = 6.21;   // Mean-reversion rate
    double theta = 0.019;  // Long run average volatility
    double xi = 0.61;      // "Vol of vol"

    // Create the PayOff, Option and Heston objects
    PayOff* pPayOffCall = new PayOffCall(K);
    Option* pOption = new Option(K, r, T, pPayOffCall);
    HestonEuler hest_euler(pOption, kappa, theta, xi, rho);

    // Store the correlated normal random numbers
    stdvec spot_draws(num_intervals, 0.0);  
    stdvec vol_draws(num_intervals, 0.0);   
    // Paths
    stdvec spot_prices(num_intervals, S_0);  // Vector of initial spot prices
    stdvec vol_prices(num_intervals, v_0);   // Vector of initial vol prices

    // Monte Carlo options pricing
    double payoff_sum = 0.0;
    for (unsigned i=0; i < num_sims; i++) {
        if((i+1)%10000 == 0) {
            cout << "Calculating path " << i+1 << " of " << num_sims << endl; 
        }
        generate_normal_correlation_paths(rho, spot_draws, vol_draws, num_intervals);
        hest_euler.calc_vol_path(vol_draws, vol_prices);
        hest_euler.calc_spot_path(spot_draws, vol_prices, spot_prices);
        payoff_sum += (*(pOption->pay_off))(spot_prices[num_intervals-1]);
    }

    double option_price = (double) payoff_sum / num_sims * exp(-r*T);
    cout << "Option Price: " << option_price << endl;

    // Free memory
    delete pOption;
    delete pPayOffCall;

    return 0;
}
