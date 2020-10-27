#ifndef SPFUNC_H
#define SPFUNC_H

#include <complex>

using namespace std;

complex<double> gamma(complex<double> x);

complex<double> rf(complex<double> x, complex<double> n);

complex<double> hypgeo_2f1(complex<double> a, complex<double> b,
	complex<double> c, complex<double> z);

#endif