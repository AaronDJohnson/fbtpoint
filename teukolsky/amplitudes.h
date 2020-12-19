#ifndef AMPS_H
#define AMPS_H

#include <complex>
#include <vector>

using namespace std;

complex<double> Binc(vector<complex<double>> &f_vec,
	vector<int> &n_vec, vector<complex<double>> &f_neg_vec,
	vector<int> &n_neg_vec, complex<double> nu, double lambda, double omega,
	double aa, int em, double ess=-2, double tol=1e-15);

#endif