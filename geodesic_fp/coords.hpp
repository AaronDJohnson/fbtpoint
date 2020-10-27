#ifndef COORDS_H
#define COORDS_H

void calc_equatorial_coords(double &t, double &r, double &theta, double &phi,
	double psi, double ups_r, double ups_theta,
	double ups_phi, double gamma, double r1, double r2, double r3,
	double r4, double zp, double zm, double En, double Lz, double aa,
	double slr, double ecc, double qt0, double qr0, double qz0,
	double qphi0);

void calc_circular_eq_coords(double &t, double &r, double &theta, double &phi,
				 			 double psi, double En, double Lz, double aa,
                                 double slr, double M);

#endif