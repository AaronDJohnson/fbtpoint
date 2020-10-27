#ifndef HYPGEO_H
#define HYPGEO_H

#include <complex>
#include <vector>

using namespace std;

complex<double> find_hypgeo(int n, vector<complex<double>> &hyp_vec, 
	vector<int> &nf_vec, complex<double> a, complex<double> b,
	complex<double> c, complex<double> d);

complex<double> find_hypgeo_deriv(int n, vector<complex<double>> &hyp_vec, 
	vector<int> &nf_vec, complex<double> a, complex<double> b,
	complex<double> c, complex<double> d);

#endif