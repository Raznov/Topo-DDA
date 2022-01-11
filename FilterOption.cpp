#include "definition.h"

FilterOption::FilterOption(double beta_min_, double beta_max_, double ita_, string beta_type_, double rfilter_, bool fixit_, int MAX_ITERATION_FIXED_) {
	beta_min = beta_min_;
	beta_max = beta_max_;
	beta = beta_min;
	ita = ita_;
	beta_type = beta_type_;
	rfilter = rfilter_;
	fixit = fixit_;
	MAX_ITERATION_FIXED = MAX_ITERATION_FIXED_;
	if (fixit) {
		cout << "MAX_ITERATION_FIXED: " << MAX_ITERATION_FIXED << endl;
	}
}

void FilterOption::update_beta(const int iteration, const int Max_iteration) {
	int MAX_iteration_real = Max_iteration;
	if (fixit) {
		MAX_iteration_real = MAX_ITERATION_FIXED;
	}
	
	if (beta_type == "exp") {
		beta=exp_update(iteration, MAX_iteration_real, beta_min, beta_max);
		return;
	}
	else if(beta_type == "piecewise") {
		beta = piecewise_update(iteration, MAX_iteration_real, beta_min, beta_max);
		return;
	}
	else if (beta_type == "linear") {
		beta = linear_update(iteration, MAX_iteration_real, beta_min, beta_max);
		return;
	}
	else {
		cout << "ERROR: FilterOption::update_beta: beta_type not defined" << endl;
		throw 1;
		return;
	}
}

double FilterOption::SmoothDensity(double input) {
	double result = 0.0;
	if (input <= ita && input >= 0.0) {
		return ita * (exp(-beta * (1 - input / ita)) - (1 - input / ita) * exp(-beta));
	}
	else if (input > ita && input <= 1.0) {
		return (1 - ita) * (1 - exp(-beta * (input - ita) / (1 - ita)) + (input - ita) / (1 - ita) * exp(-beta)) + ita;
	}
	else {
		cout << "ERROR: FilterOption::SmoothDensity(double input)--input out of range" << endl;
		throw 1;
	}
}

double FilterOption::get_beta() {
	return beta;
}

double FilterOption::get_ita() {
	return ita;
}

double FilterOption::get_rfilter() {
	return rfilter;
}