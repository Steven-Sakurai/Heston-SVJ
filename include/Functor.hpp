#pragma once

#include <cmath>
#include "Jump.hpp"

double prePy(const stdVec& par, double y, double y0, double v0) {
	double a = pow(y - tran_y_mean(par, y0, v0), 2);
	double b = tran_y_var(par, v0);
	return exp(-a/(2*b)) / sqrt(2*Pi*b);
}

double pdf_c(const stdVec& par, double y, double y0, double v, double v0) {
	double mu = cond_v_mean(par, y, y0, v0);
	double var = cond_v_var(par, v0);
	return exp(-pow(v-mu, 2)/(2*var)) / sqrt(2*Pi*var);
}

double cdf_c(const stdVec& par, double y, double y0, double v, double v0) {
	double mu = cond_v_mean(par, y, y0, v0);
	double var = cond_v_var(par, v0);
	return 1 - 0.5*erfc((v-mu)/sqrt(var)/sqrt(2));
}

double Integ0(const stdVec& par, double y, double y0, double v1, double v2, double v0) {
	return prePy(par, y, y0, v0) * (cdf_c(par, y, y0, v2, v0) - cdf_c(par, y, y0, v1, v0));
}

double Integ1(const stdVec& par, double y, double y0, double v1, double v2, double v0) {
	double mu_c = cond_v_mean(par, y, y0, v0);
	double var_c = cond_v_var(par, v0);
	return Integ0(par, y, y0, v1, v2, v0)*mu_c  - 
		prePy(par, y, y0, v0)*var_c*(pdf_c(par, y, y0, v2, v0) - 
			pdf_c(par, y, y0, v1, v0));
}

double Integ2(const stdVec& par, double y, double y0, double v1, double v2, double v0) {
	double mu_c = cond_v_mean(par, y, y0, v0);
	double var_c = cond_v_var(par, v0);
	return Integ0(par, y, y0, v1, v2, v0)*var_c + Integ1(par, y, y0, v1, v2, v0)*mu_c - 
		prePy(par, y, y0, v0)*var_c*(pdf_c(par, y, y0, v2, v0)*v2 - pdf_c(par, y, y0, v1, v0)*v1);
}

double Integ3(const stdVec& par, double y, double y0, double v1, double v2, double v0) {
	double mu_c = cond_v_mean(par, y, y0, v0);
	double var_c = cond_v_var(par, v0);
	return 2*Integ1(par, y, y0, v1, v2, v0)*var_c + 
		Integ2(par, y, y0, v1, v2, v0)*mu_c - 
		prePy(par, y, y0, v0)*var_c*(pdf_c(par, y, y0, v2, v0)*v2*v2 - 
			pdf_c(par, y, y0, v1, v0)*v1*v1);
}

