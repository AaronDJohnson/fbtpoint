#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include <math.h>
#include <boost/math/special_functions/jacobi_elliptic.hpp>

// calculates the coordinates as fast as possible to machine precision
// this can be extended to multiple precision if needed using arblib

// declare const parameters -- this could lead to speed increases

// TODO: there is a problem with ecc = 0 in this method


double am(double x, double m) {
	if (m == 0) {
		return x;
	} else if (x == 0) {
		return 0;
	} else {
		// function to calculate jacobi amplitude for real input
		// "borrowed" from sagemath
		double sn, cn, dn;
		if (fabs(m) == 1) {
			double tanhx = tanh(x);
			return asin(tanhx);
		} else if (fabs(m) > 1){
			gsl_sf_elljac_e(x, m, &sn, &cn, &dn);
			return atan2(sn, cn);
		}
		else{
			double K = gsl_sf_ellint_Kcomp(sqrt(m), GSL_PREC_DOUBLE);
			if (fabs(x) <= K) {
				gsl_sf_elljac_e(x, m, &sn, &cn, &dn);
				return atan2(sn, cn);
			} else {
				// Do argument reduction on x to end up with z = x - 2nK, with
				// abs(z) <= K
				double tK = 2 * K;
				int n = floor(x / tK);
				double tnK = n * tK;
				double npi = n * M_PI;
				double z = x - tnK;

				gsl_sf_elljac_e(z, m, &sn, &cn, &dn);
				return atan2(sn, cn) + npi;
			}
		}
	}
}


double calc_radius(double psi, double slr, double ecc) {
	double r = slr / (1 + ecc * cos(psi));
	// printf("slr = %.17g \n", slr);
	// printf("ecc = %.17g \n", ecc);
	// printf("psi = %.17g \n", psi);
	return r;
}


double calc_lambda_r(double r, double r1, double r2, double r3,
	double r4, double En) {

	// compute function to change from t[lambda] -> t[psi] with psi an angle
	// associated with radial motion. 

	double yr, F_asin, kr;
	kr = ((r1 - r2)*(r3 - r4))/((r1 - r3)*(r2 - r4));
	// there is some problem here because of the high precision part of the
	// code. r can be less than r2 (peribothron).
	// TODO: fix the code, but for now a temporary fix:
	if (r < r2) {
		r = r2;
		yr = sqrt(((r - r2)*(r1 - r3))/((r1 - r2)*(r - r3)));
	}
	else {
		yr = sqrt(((r - r2)*(r1 - r3))/((r1 - r2)*(r - r3)));
	}
	// printf("r = %.17g \n", r);
	// printf("r2 = %.17g \n", r2);
	F_asin = gsl_sf_ellint_F(asin(yr), sqrt(kr), GSL_PREC_DOUBLE);
	return (2*F_asin)/(sqrt(1 - En * En)*sqrt((r1 - r3)*(r2 - r4)));
}


double calc_lambda_psi(double psi, double ups_r, double r1,
	double r2, double r3, double r4, double En, double slr, double ecc) {

	double lam_r, lam_r1, r, turns, res, kr;
	kr = ((r1 - r2)*(r3 - r4))/((r1 - r3)*(r2 - r4));
	r = calc_radius(psi, slr, ecc);
	lam_r = 2 * M_PI / ups_r;  // radial period
	lam_r1 = calc_lambda_r(r2, r1, r2, r3, r4, En);
	turns = floor(psi / (2 * M_PI));
	if (fmod(psi, 2 * M_PI) <= M_PI) {
		// printf("psi = %.17g\n", psi);
		res = calc_lambda_r(r, r1, r2, r3, r4, En) - lam_r1;
		// printf("res = %.17g\n", res);
	} else {
		res = lam_r1 - calc_lambda_r(r, r1, r2, r3, r4, En);
	}
	// printf("two\n");

	return lam_r * turns + res;
}


double calc_lambda_chi(double chi, double zp, double zm, double En, double aa) {
	// compute function to change from t[lambda] -> t[chi] with chi an angle
	// associated with polar motion. 
	double beta, k, k2, prefactor, ellipticK_k, ellipticF;
	beta = aa * aa * (1 - En * En);
	k = sqrt(zm / zp);
	k2 = k * k;
	prefactor = 1 / sqrt(beta * zp);
	ellipticK_k = gsl_sf_ellint_Kcomp(sqrt(k2), GSL_PREC_DOUBLE);
	ellipticF = gsl_sf_ellint_F(M_PI / 2 - chi, sqrt(k2), GSL_PREC_DOUBLE);

	return prefactor * (ellipticK_k - ellipticF);
}

// routines used to find t, r, theta, phi


double calc_theta_chi(double chi, double zm) {
	return acos(zm * cos(chi));
}


double calc_rq(double qr, double r1, double r2,
	double r3, double r4) {

	double arg, sn2, sn, cn, dn, kr;

	kr = ((r1 - r2)*(r3 - r4))/((r1 - r3)*(r2 - r4));

	arg = (qr*gsl_sf_ellint_Kcomp(sqrt(kr), GSL_PREC_DOUBLE))/M_PI;
	gsl_sf_elljac_e(arg, kr, &sn, &cn, &dn);
	sn2 = sn * sn;

	return ((-(r2*(r1 - r3)) + (r1 - r2)*r3*sn2)/ (-r1 + r3 + (r1 - r2)*sn2));
}


double calc_zq(double qz, double zp, double zm, double En, double aa) {
	double ktheta, sn, dn, cn, ellipticK_k, arg1, arg2;
	ktheta = (aa*aa *(1 - En*En)*zm*zm)/zp*zp;
	ellipticK_k = gsl_sf_ellint_Kcomp(sqrt(ktheta), GSL_PREC_DOUBLE);
	arg1 = (2*(M_PI/2. + qz)*ellipticK_k)/M_PI;
	arg2 = ktheta;
	gsl_sf_elljac_e(arg1, arg2, &sn, &cn, &dn);

	return zm * sn;
}


double calc_psi_r(double qr, double r1, double r2, double r3, double r4) {
	double kr, amplitude;
	kr = ((r1 - r2)*(r3 - r4))/((r1 - r3)*(r2 - r4));
	amplitude = am((qr * gsl_sf_ellint_Kcomp(sqrt(kr), GSL_PREC_DOUBLE))/M_PI,
				   kr);
	return amplitude;
}


double calc_t_r(double qr, double r1, double r2, double r3, double r4,
	double En, double Lz, double aa, double M=1) {
	double psi_r, kr, rp, rm, hr, hp, hm, aa2, En2, sin_psi_r, ellipe_k, ince_k;
	double ellippi_hm, ellippi_hr, ellippi_hp, inc_pi_hm, inc_pi_hr, inc_pi_hp;
	aa2 = aa * aa;
	En2 = En * En;
	psi_r = calc_psi_r(qr, r1, r2, r3, r4);

	kr = ((r1 - r2)*(r3 - r4))/((r1 - r3)*(r2 - r4));

	rp = M + sqrt(-aa*aa + M*M);
	rm = M - sqrt(-aa*aa + M*M);

	hr = (r1 - r2)/(r1 - r3);
	hp = ((r1 - r2)*(r3 - rp))/((r1 - r3)*(r2 - rp));
	hm = ((r1 - r2)*(r3 - rm))/((r1 - r3)*(r2 - rm));

	sin_psi_r = pow(sin(psi_r), 2);
	ince_k = gsl_sf_ellint_E(psi_r, sqrt(kr), GSL_PREC_DOUBLE);
	ellipe_k = gsl_sf_ellint_Ecomp(sqrt(kr), GSL_PREC_DOUBLE);
	ellippi_hm = gsl_sf_ellint_Pcomp(sqrt(kr), -hm, GSL_PREC_DOUBLE);
	ellippi_hp = gsl_sf_ellint_Pcomp(sqrt(kr), -hp, GSL_PREC_DOUBLE);
	ellippi_hr = gsl_sf_ellint_Pcomp(sqrt(kr), -hr, GSL_PREC_DOUBLE);

	// printf("ellippi_hm = %.17g \n", ellippi_hm);

	inc_pi_hm = gsl_sf_ellint_P(psi_r, sqrt(kr), -hm, GSL_PREC_DOUBLE);
	inc_pi_hp = gsl_sf_ellint_P(psi_r, sqrt(kr), -hp, GSL_PREC_DOUBLE);
	inc_pi_hr = gsl_sf_ellint_P(psi_r, sqrt(kr), -hr, GSL_PREC_DOUBLE);


	return  (-((En * ((-4*(r2 - r3)*(-(((-2*aa2 + (4 - (aa*Lz)/En)*rm)*
			((qr*ellippi_hm)/M_PI - inc_pi_hm))/
			((r2 - rm)*(r3 - rm))) + 
			((-2*aa2 + (4 - (aa*Lz)/En)*rp)*
			((qr*ellippi_hp)/M_PI - inc_pi_hp))/
			((r2 - rp)*(r3 - rp))))/(-rm + rp) + 
			4*(r2 - r3)*((qr*ellippi_hr)/M_PI - inc_pi_hr) + 
			(r2 - r3)*(r1 + r2 + r3 + r4)*
			((qr*ellippi_hr)/M_PI - inc_pi_hr) + 
			(r1 - r3)*(r2 - r4)*((qr*ellipe_k)/M_PI - ince_k + 
			(hr*cos(psi_r)*sin(psi_r)*sqrt(1 - kr*sin_psi_r))/
			(1 - hr*sin_psi_r))))/sqrt((1 - En2)*(r1 - r3)*(r2 - r4))));

}


double calc_phi_r(double qr, double r1, double r2, double r3, double r4,
	double En, double Lz, double aa, double M=1){

	double psi_r, kr, rp, rm, hp, hm, aa2, En2;
	double ellpi_hp, ellpi_hm, incpi_hp, incpi_hm;
	aa2 = aa * aa;
	En2 = En * En;
	psi_r = calc_psi_r(qr, r1, r2, r3, r4);
	kr = ((r1 - r2)*(r3 - r4))/((r1 - r3)*(r2 - r4));
	rp = M + sqrt(-aa2 + M*M);
	rm = M - sqrt(-aa2 + M*M);

	hp = ((r1 - r2)*(r3 - rp))/((r1 - r3)*(r2 - rp));
	hm = ((r1 - r2)*(r3 - rm))/((r1 - r3)*(r2 - rm));

	ellpi_hp = gsl_sf_ellint_Pcomp(sqrt(kr), -hp, GSL_PREC_DOUBLE);
	ellpi_hm = gsl_sf_ellint_Pcomp(sqrt(kr), -hm, GSL_PREC_DOUBLE);
	incpi_hp = gsl_sf_ellint_P(psi_r, sqrt(kr), -hp, GSL_PREC_DOUBLE);
	incpi_hm = gsl_sf_ellint_P(psi_r, sqrt(kr), -hm, GSL_PREC_DOUBLE);

	return ((2*aa*En*(-(((r2 - r3)*(-((aa*Lz)/En) + 2*rm)*
			((qr*ellpi_hm)/M_PI - incpi_hm))/
			((r2 - rm)*(r3 - rm))) + 
			((r2 - r3)*(-((aa*Lz)/En) + 2*rp)*
			((qr*ellpi_hp)/M_PI - incpi_hp))/
			((r2 - rp)*(r3 - rp))))/
			(sqrt((1 - En2)*(r1 - r3)*(r2 - r4))*(-rm + rp)));

}


double calc_psi_z(double qz, double zp, double zm, double En, double aa){
	double ktheta, ellK_k;
	ktheta = (aa*aa *(1 - En*En)*zm*zm)/zp*zp;
	ellK_k = gsl_sf_ellint_Kcomp(sqrt(ktheta), GSL_PREC_DOUBLE);
	return am((2*(M_PI/2. + qz) * ellK_k)/M_PI, ktheta);
}


double calc_t_z(double qz, double zp, double zm, double En, double aa){
	double psi_z, ktheta, elle_k, ince_k;
	psi_z = calc_psi_z(qz, zp, zm, En, aa);
	ktheta = (aa*aa *(1 - En*En)*zm*zm)/zp*zp;
	elle_k = gsl_sf_ellint_Ecomp(sqrt(ktheta), GSL_PREC_DOUBLE);
	ince_k = gsl_sf_ellint_E(psi_z, sqrt(ktheta), GSL_PREC_DOUBLE);
	return ((En*zp*((2*(M_PI/2. + qz)*elle_k)/M_PI - ince_k))/
			(1 - En*En));
}


double calc_phi_z(double qz, double zp, double zm, double En, double Lz,
	double aa){

	double psi_z, ktheta, ellpi_k, incpi_k;
	psi_z = calc_psi_z(qz, zp, zm, En, aa);
	ktheta = (aa*aa * (1 - En*En)*zm*zm)/zp*zp;
	ellpi_k = gsl_sf_ellint_Pcomp(sqrt(ktheta), -zm*zm, GSL_PREC_DOUBLE);
	incpi_k = gsl_sf_ellint_P(psi_z, sqrt(ktheta), -zm*zm, GSL_PREC_DOUBLE);
	return (-((Lz*((2*(M_PI/2. + qz)*ellpi_k)/M_PI - 
			incpi_k))/zp));
}


double calc_Ct(double qr0, double qz0, double r1, double r2, double r3,
	double r4, double zp, double zm, double En, double Lz, double aa) {

	double t_r, t_z;
	t_r = calc_t_r(qr0, r1, r2, r3, r4, En, Lz, aa);
	t_z = calc_t_z(qz0, zp, zm, En, aa);
	return t_r + t_z;
}


double calc_Cz(double qr0, double qz0, double r1, double r2, double r3, double r4,
	double zp, double zm, double En, double Lz, double aa) {

	double phi_r, phi_z;
	phi_r = calc_phi_r(qr0, r1, r2, r3, r4, En, Lz, aa);
	phi_z = calc_phi_z(qz0, zp, zm, En, Lz, aa);
	return phi_r + phi_z;
}


double calc_t(double mino_t, double ups_r, double ups_theta, double gamma,
	double qt0, double qr0, double qz0, double r1, double r2, double r3,
	double r4, double zp, double zm, double En, double Lz, double aa){

	double eta_t, eta_r, eta_z, t_r, t_z, Ct;
	eta_t = qt0 + gamma * mino_t;
	eta_r = qr0 + ups_r * mino_t;
	eta_z = qz0 + ups_theta * mino_t;
	if (r1 == r2) {
		t_r = 0;
	} else {
		t_r = calc_t_r(eta_r, r1, r2, r3, r4, En, Lz, aa);
		// printf("t_r = %.17g \n", t_r);
	}
	if (zm == 0){
		t_z = 0;
	} else {
		t_z = calc_t_z(eta_z, zp, zm, En, aa);
	}
	if (qr0 == 0 && qz0 == 0) {
		Ct = 0;
	} else {
		Ct = calc_Ct(qr0, qz0, r1, r2, r3, r4, zp, zm, En, Lz, aa);
	}
	return eta_t + t_r + t_z - Ct;
}


double calc_r(double mino_t, double ups_r, double qr0, double r1, double r2,
	double r3, double r4){

	double eta;
	eta = ups_r * mino_t + qr0;
	return calc_rq(eta, r1, r2, r3, r4);
}


double calc_theta(double mino_t, double ups_theta, double qz0, double zp,
	double zm, double En, double aa){

	double eta;
	eta = ups_theta * mino_t + qz0;
	return acos(calc_zq(eta, zp, zm, En, aa));
}


double calc_phi(double mino_t, double ups_r, double ups_theta, double ups_phi,
	double qphi0, double qr0, double qz0, double r1, double r2, double r3,
	double r4, double zp, double zm, double En, double Lz, double aa){

	double eta_phi, eta_r, eta_theta, phi_r, phi_z, Cz;
	eta_phi = ups_phi * mino_t + qphi0;
	eta_r = ups_r * mino_t + qr0;
	eta_theta = ups_theta * mino_t + qz0;
	if (r1 == r2) {
		phi_r = 0;
	} else {
		phi_r = calc_phi_r(eta_r, r1, r2, r3, r4, En, Lz, aa);
	}

	if (zm == 0) {
		phi_z = 0;
	}
	else {
		phi_z = calc_phi_z(eta_theta, zp, zm, En, Lz, aa);
	}

	if (qr0 == 0 && qz0 == 0){
		Cz = 0;
	} else {
		Cz = calc_Cz(qr0, qz0, r1, r2, r3, r4, zp, zm, En, Lz, aa);
	}
	return eta_phi + phi_r + phi_z - Cz;
}


void calc_circular_eq_coords(double &t, double &r, double &theta, double &phi,
							 double psi, double En, double Lz, double aa, double slr, double M) {
    // lam_psi is used to cross check with circular orbits in BHPTK
    // lam_psi = (psi / (sqrt(1 - En**2 + 3*(1 - En**2) + 2*(1 - En**2 - 1/slr) + 
    //            (aa**2*(1 - En**2) + Lz**2)/slr**2 - 4/slr)*slr))
    // print(lam_psi)

	double x, Vr, Vphi, Vt, J, denom;
	double aa2 = aa * aa;
	double slr2 = slr * slr;
	x = Lz - aa * En;
	double x2 = x * x;

    Vr = aa2 + x2 + 2 * aa * x * En - 6 * M * x2 / slr;
    Vphi = x + aa * En - 2 * M * x / slr;
    Vt = aa2 * En - 2 * aa * M * x / slr + En * slr2;
    J = 1 - 2 * M / slr + aa2 / slr2;
    denom = J * sqrt(Vr);
	t = Vt / denom * psi;
    phi = Vphi / denom * psi;
    r = slr;
	theta = M_PI / 2;
}


void calc_equatorial_coords(double &t, double &r, double &theta, double &phi,
	double psi, double ups_r, double ups_theta,
	double ups_phi, double gamma, double r1, double r2, double r3,
	double r4, double zp, double zm, double En, double Lz, double aa,
	double slr, double ecc, double qt0, double qr0, double qz0,
	double qphi0) {

	double lam_psi;
	if (zm != 0){
		printf("The orbit specified is not equatorial.\n");
	}
	r = calc_radius(psi, slr, ecc);
	lam_psi = calc_lambda_psi(psi, ups_r, r1, r2, r3, r4, En, slr, ecc);
	t = calc_t(lam_psi, ups_r, ups_theta, gamma, qt0, qr0, qz0, r1, r2, r3, r4,
			   zp, zm, En, Lz, aa);
	theta = M_PI / 2;  // TODO(aaron): check generic theta
	phi = calc_phi(lam_psi, ups_r, ups_theta, ups_phi, qphi0, qr0, qz0, r1, r2,
				   r3, r4, zp, zm, En, Lz, aa);

}


// int main() {
// 	double t, r, theta, phi;
// 	double psi = M_PI / 2;
// 	double En = 0.95656754033756958293025762715202237237836692590871, Lz = 3.7823473723611689443109490343829688125249006073227, Q = 0;
// 	double r1 = 11.111111111111111111111111111111111111111111111111, r2 = 9.0909090909090909090909090909090909090909090909091, r3 = 3.3333333333333333333333333333333333333333333333337, r4 = 0;
// 	double zp = 2.77210885820935, zm = 0;
// 	double ups_r = 2.08698701855696, ups_theta = 2.77210885820935;
// 	double ups_phi = 3.080504407723, gamma = 50.492536214173;

// 	double aa = 0, slr = 10.0, ecc = 0.1;

// 	double qt0 = 0;
// 	double qr0 = 0;
// 	double qz0 = 0;
// 	double qphi0 = 0;

// 	//double ellipticK_k = gsl_sf_ellint_Kcomp(sqrt(ktheta), GSL_PREC_DOUBLE);

// 	calc_equatorial_coords(t, r, theta, phi,
// 	psi, ups_r, ups_theta, ups_phi, gamma, r1, r2, r3, r4, zp, zm, En, Lz, aa,
// 	slr, ecc, qt0, qr0, qz0, qphi0);

// 	printf("t = %.17g \n", t);
// 	printf("r = %.17g \n", r);
// 	printf("theta = %.17g \n", theta);
// 	printf("phi = %.17g \n", phi);



// 	return 0;
// }

