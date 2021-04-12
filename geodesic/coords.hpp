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


void calc_gen_coords_mino(double &t, double &r, double &theta, double &phi,
	double mino_t, double ups_r, double ups_theta, double ups_phi, double gamma,
    double r1, double r2, double r3, double r4, double zp, double zm, double En,
    double Lz, double aa, double qphi0, double qr0, double qz0, double qt0);


// void calc_gen_coords_psi(double &t, double &r, double &theta, double &phi,
// 	double psi, double ups_r, double ups_theta, double ups_phi, double gamma,
//     double r1, double r2, double r3, double r4, double zp, double zm, double En,
// 	double Lz, double Q, double aa, double slr, double ecc, double x,
// 	double qphi0, double qr0, double qz0, double qt0);


// double calc_wr(double psi, double ups_r, double En, double Lz,
// 			   double Q, double aa, double slr, double ecc, double x);


double calc_dwr_dpsi(double psi, double ups_r, double En, double Lz,
					 double Q, double aa, double slr, double ecc);


// double calc_wtheta(double chi, double ups_theta, double zp, double zm,
// 				   double En, double Lz, double aa, double slr, double x);


double calc_dwtheta_dchi(double chi, double zp, double zm);



#endif