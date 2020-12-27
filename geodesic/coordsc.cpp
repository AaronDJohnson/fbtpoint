#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include <math.h>
#include <complex>
#include <boost/math/special_functions/jacobi_elliptic.hpp>
#include <boost/math/special_functions/ellint_1.hpp>


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


double calc_J(double chi, double En, double Lz, double Q, double aa, double slr, double ecc) {
    // """
    // Schmidt's J function.

    // Parameters:
    //     chi (float): radial angle
    //     En (float): energy
    //     Lz (float): angular momentum
    //     Q (float): Carter constant
    //     aa (float): spin
    //     slr (float): semi-latus rectum
    //     ecc (float): eccentricity

    // Returns:
    //     J (float)
    // """
    double En2 = En * En;
    double ecc2 = ecc * ecc;
    double aa2 = aa * aa;
    double Lz2 = Lz * Lz;
    double slr2 = slr * slr;

    double eta = 1 + ecc * cos(chi);
    double eta2 = eta * eta;

    double J = (
        (1 - ecc2) * (1 - En2)
        + 2 * (1 - En2 - (1 - ecc2) / slr) * eta
        + (
            ((3 + ecc2) * (1 - En2)) / (1 - ecc2)
            + ((1 - ecc2) * (aa2 * (1 - En2) + Lz2 + Q)) / slr2
            - 4 / slr
        )
        * eta2
    );

    return J;
}


double calc_wr(double psi, double ups_r, double En, double Lz,
			   double Q, double aa, double slr, double ecc, double x) {
	// note that this has only been checked for psi in [0, pi]
	double aa2 = aa * aa;
	double slr2 = slr * slr;
	double ecc2 = ecc * ecc;
	double En2 = En * En;
	double Lz2 = Lz * Lz;
	double a1 = (1 - ecc2) * (1 - En2);
    double b1 = 2 * (1 - En2 - (1 - ecc2) / slr);
    double c1 = (((3 + ecc2)*(1 - En2))/(1 - ecc2) - 4/slr + 
          		 ((1 - ecc2) * (aa2*(1 - En2) + Lz2 + Q)) / slr2);

	double b12 = b1 * b1
	if (psi == M_PI) {
		return M_PI;
	} else {
		complex<double> phi {0, asinh(sqrt((a1 - (-1 + ecc)*(b1 + c1 - c1*ecc))/
       						 (a1 + b1 + c1 - c1*ecc2 + sqrt(b12 - 4*a1*c1)*ecc2)))* tan(psi/2.))};

	 	double m = ((a1 + b1 + c1 - c1*ecc2 + sqrt((b12 - 4*a1*c1)*ecc2))/
  				    (a1 + b1 + c1 - c1*ecc2 - sqrt((b12 - 4*a1*c1)*ecc2)));
		double k = sqrt(m)
		complex<double> ellint_f = ellint_1(k, phi);
		complex<double> res {0, ((-2*(1 - ecc2) * ellint_f * ups_r * power(cos(psi/2.),2)*
				sqrt(2 + (2*(a1 - (-1 + ecc)*(b1 + c1 - c1*ecc))*power(tan(psi/2.),2))/
					(a1 + b1 + c1 - c1*ecc2 - sqrt(b12 - 4*a1*c1)*ecc2)))*
				sqrt(1 + ((a1 - (-1 + ecc)*(b1 + c1 - c1*ecc))*power(tan(psi/2.),2))/
					(a1 + b1 + c1 - c1*ecc2 + sqrt((b12 - 4*a1*c1)*ecc2))))/
			(sqrt((a1 - (-1 + ecc)*(b1 + c1 - c1*ecc))/
				(a1 + b1 + c1 - c1*ecc2 + sqrt((b12 - 4*a1*c1)*ecc2)))*
				slr*sqrt(2*a1 + 2*b1 + 2*c1 + c1*ecc2 + 2*(b1 + 2*c1)*ecc*cos(psi) + 
				c1*ecc2*cos(2*psi))))};
		double re_res = real(res);
		
		return re_res;
	}
}


double calc_dwr_dpsi(double psi, double ups_r, double En, double Lz,
					 double Q, double aa, double slr, double ecc) {
    double J = calc_J(psi, En, Lz, Q, aa, slr, ecc);
	double ecc2 = ecc * ecc;
    return (1 - ecc2) / slr * ups_r / sqrt(J);
}


double calc_wtheta(double chi, double ups_theta, double zp, double zm,
				   double En, double Lz, double aa, double slr, double x) {
    // """
    // w_theta = ups_theta * lambda as a function of polar angle chi

    // Parameters:
    //     chi (float): polar angle
    //     ups_theta (float): theta mino frequency
    //     En (float): energy
    //     Lz (float): angular momentum
    //     aa (float): spin
    //     slr (float): semi-latus rectum
    //     x (float): inclination

    // Returns:
    //     w_theta (float)
    // """
    double pi = M_PI;
    if (chi >= 0 && chi <= pi / 2) {
        return ups_theta * calc_lambda_chi(chi, zp, zm, En, Lz, aa, slr, x);
	} else if (chi > pi / 2 && chi <= pi) {
        return pi - ups_theta * calc_lambda_chi(pi - chi, zp, zm, En, Lz, aa, slr, x);
	} else if (chi > pi && chi <= 3 * pi / 2) {
        return pi + ups_theta * calc_lambda_chi(chi - pi, zp, zm, En, Lz, aa, slr, x);
	} else if (chi > 3 * pi / 2 && chi <= 2 * pi) {
        return 2 * pi - ups_theta * calc_lambda_0(2 * pi - chi, zp, zm, En, Lz, aa, slr, x);
	} else {
        // print("Something went wrong in calc_wtheta!")
        return 0.0;  // this case should not occur, but is required by C++
	}
}


double calc_dwtheta_dchi(double chi, double zp, double zm) {
    // """
    // derivative of w_theta

    // Parameters:
    //     chi (float): polar angle
    //     zp (float): polar root
    //     zm (float): polar root

    // Returns:
    //     dw_dtheta (float)
    // """
    double k = sqrt(zm / zp);
    double ellipticK_k = gsl_sf_ellint_Kcomp(k, GSL_PREC_DOUBLE);
    return M_PI / (2 * ellipticK_k) * (1 / (1 - k * k * power(cos(chi),2)));
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


void calc_gen_coords_mino(double &t, double &r, double &theta, double &phi,
	double mino_t, double ups_r, double ups_theta, double ups_phi, double gamma,
    double r1, double r2, double r3, double r4, double zp, double zm, double En,
    double Lz, double aa, double qphi0, double qr0, double qz0, double qt0) {

	t = calc_t(mino_t, ups_r, ups_theta, gamma, qt0, qr0, qz0, r1, r2, r3, r4, zp, zm, En, Lz, aa);
	r = calc_r(mino_t, ups_r, qr0, r1, r2, r3, r4);
	theta = calc_theta(mino_t, ups_theta, qz0, zp, zm, En, aa);
	phi = calc_phi(mino_t, ups_r, ups_theta, ups_phi, qphi0, qr0, qz0, r1, r2, r3, r4, zp, zm, En, Lz, aa);
}


void calc_gen_coords_psi(double &t, double &r, double &theta, double &phi,
	double psi, double ups_r, double ups_theta, double ups_phi, double gamma,
    double r1, double r2, double r3, double r4, double zp, double zm, double En,
	double Lz, double Q, double aa, double slr, double ecc, double x,
	double qphi0, double qr0, double qz0, double qt0) {

	double wr = calc_wr(psi, ups_r, En, Lz, Q, aa, slr, ecc, x);
	double mino_t = wr / ups_r;
	t = calc_t(mino_t, ups_r, ups_theta, gamma, qt0, qr0, qz0, r1, r2, r3, r4, zp, zm, En, Lz, aa);
	r = calc_r(mino_t, ups_r, qr0, r1, r2, r3, r4);
	theta = calc_theta(mino_t, ups_theta, qz0, zp, zm, En, aa);
	phi = calc_phi(mino_t, ups_r, ups_theta, ups_phi, qphi0, qr0, qz0, r1, r2, r3, r4, zp, zm, En, Lz, aa);
}


