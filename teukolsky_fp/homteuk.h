#ifndef HOMTEUK_H
#define HOMTEUK_H

#include <complex>
#include <vector>

using namespace std;

void find_R(double &re_R_in, double &im_R_in, double &re_Rd_in,
			double &im_Rd_in, double &re_Rdd_in, double &im_Rdd_in,
			double r, double re_nu, double im_nu, double eigen, double aa,
			double omega, int em, int ess=-2, double tol=1e-15, double M=1.);

void find_Bin(double &re_B_in, double &im_B_in, double re_nu, double im_nu,
			  double eigen, double aa, double omega,
			  int em, int ess=-2, double tol=1e-15);

#endif