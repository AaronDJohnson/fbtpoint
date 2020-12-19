#ifndef FMINSOL_H
#define FMINSOL_H

#include <complex>
#include <vector>

using namespace std;

complex<double> find_fn(int n, vector<complex<double>> &f_vec,
	vector<int> &n_vec, complex<double> nu, double lambda, double omega,
	double aa, int em, int ess=-2);

#endif