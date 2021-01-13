#include "fminsol.h"
#include "sp_func.h"
#include <gsl/gsl_sf.h>
// #include <gsl/gsl_complex.h>
#include <complex>
#include <vector>
#include <cmath>

using namespace std;


complex<double> f_sum_up(vector<complex<double>> &f_vec,
	vector<int> &n_vec, complex<double> nu, double lambda, double omega,
	double aa, int em, double tol=1e-15) {
	//
	// old
	complex<double> f_sum, fn;
	double re_err, im_err;
	int n = 0; 
	fn = find_fn(n, f_vec, n_vec, nu, lambda, omega, aa, em);
	f_sum += fn;

	//new
	n = 1;
	fn = find_fn(n, f_vec, n_vec, nu, lambda, omega, aa, em);
	f_sum += fn;

	re_err = fabs(real(fn));
    im_err = fabs(imag(fn));

    while (re_err > tol or im_err > tol) {
    	n += 1;
    	fn = find_fn(n, f_vec, n_vec, nu, lambda, omega, aa, em);
    	f_sum += fn;

    	re_err = fabs(real(fn));
    	im_err = fabs(imag(fn));
    	if (n == 1000) {
            printf("f_sum_up did not converge.\n");
            return 0.;
    	}
    }
    return f_sum;
}


complex<double> f_sum_dn(vector<complex<double>> &f_vec,
	vector<int> &n_vec, complex<double> nu, double lambda, double omega,
	double aa, int em, double tol=1e-15) {
	//
	// old
	complex<double> f_sum, fn;
	double re_err, im_err;

	//new
	int n = -1;
	fn = find_fn(n, f_vec, n_vec, nu, lambda, omega, aa, em);
	// printf("fn = %.17g + %.17gi \n", real(fn), imag(fn));
	f_sum += fn;

	re_err = fabs(real(fn));
    im_err = fabs(imag(fn));

    while (re_err > tol or im_err > tol) {
    	n -= 1;
    	fn = find_fn(n, f_vec, n_vec, nu, lambda, omega, aa, em);
    	// printf("fn = %.17g + %.17gi \n", real(fn), imag(fn));
    	f_sum += fn;

    	re_err = fabs(real(fn));
    	im_err = fabs(imag(fn));
    	if (n == 1000) {
            printf("f_sum_up did not converge.\n");
            return 0.;
    	}
    }
    return f_sum;
}


complex<double> f_sum_Knu_up(vector<complex<double>> &f_vec,
	vector<int> &n_vec, complex<double> nu, double lambda, double omega,
	double aa, int em, double ess=-2, double tol=1e-15) {

	int n = 0;
	double epsilon = 2 * omega;
	complex<double> prefactor;

	double kappa = sqrt(1 - aa * aa);
	double tau = (epsilon - em * aa) / kappa;

	complex<double> i_eps {0, epsilon};
	complex<double> i_tau {0, tau};
	complex<double> f_new, f_sum;

	prefactor = ((pow(-1, n) * gamma(1. + n + i_eps + ess + nu) *
                gamma(1. + n + nu + nu) * gamma(1. + n + nu + i_tau)) /
                (gsl_sf_fact(n) * gamma(1. + n - i_eps - ess + nu) *
                gamma(1. + n + nu - i_tau)));
	// old
	f_sum = prefactor;
	n = 1;
	// new
	f_new = find_fn(n, f_vec, n_vec, nu, lambda, omega, aa, em);
	prefactor = ((pow(-1, n) * gamma(1. + n + i_eps + ess + nu) *
                gamma(1. + n + nu + nu) * gamma(1. + n + nu + i_tau)) /
                (gsl_sf_fact(n) * gamma(1. + n - i_eps - ess + nu) *
                gamma(1. + n + nu - i_tau)));
	f_sum += f_new * prefactor;

	double re_err = fabs(real(f_new));
	double im_err = fabs(imag(f_new));

	while (re_err > tol or im_err > tol) {
		n++;
		f_new = find_fn(n, f_vec, n_vec, nu, lambda, omega, aa, em);
		prefactor = ((pow(-1, n) * gamma(1. + n + i_eps + ess + nu) *
                gamma(1. + n + nu + nu) * gamma(1. + n + nu + i_tau)) /
                (gsl_sf_fact(n) * gamma(1. + n - i_eps - ess + nu) *
                gamma(1. + n + nu - i_tau)));
		f_sum += f_new * prefactor;

		re_err = fabs(real(f_new));
		im_err = fabs(imag(f_new));
		if (n == 1000) {
			printf("f_sum_up did not converge.");
			return 0.;
		}
	}
	return f_sum;
}


complex<double> f_sum_Knu_dn(vector<complex<double>> &f_vec,
	vector<int> &n_vec, complex<double> nu, double lambda, double omega,
	double aa, int em, double ess=-2, double tol=1e-15) {

	int n = 0;
	double epsilon = 2 * omega;
	complex<double> prefactor;

	double kappa = sqrt(1 - aa * aa);
	double tau = (epsilon - em * aa) / kappa;

	complex<double> i_eps {0, epsilon};
	complex<double> i_tau {0, tau};
	complex<double> f_new, f_sum;

	prefactor = ((pow(-1, n) * rf(1. - i_eps + ess + nu, n)) /
                (gsl_sf_fact(-n) * rf(1. + i_eps - ess + nu, n) *
                 rf(2. + nu + nu, n)));
	// old
	f_sum = prefactor;

	// new
	n = -1;
	f_new = find_fn(n, f_vec, n_vec, nu, lambda, omega, aa, em);
	prefactor = ((pow(-1, n) * rf(1. - i_eps + ess + nu, n)) /
                (gsl_sf_fact(-n) * rf(1. + i_eps - ess + nu, n) *
                 rf(2. + nu + nu, n)));
	f_sum += f_new * prefactor;

	double re_err = fabs(real(f_new));
	double im_err = fabs(imag(f_new));

	while (re_err > tol or im_err > tol) {
		n--;
		f_new = find_fn(n, f_vec, n_vec, nu, lambda, omega, aa, em);
		prefactor = ((pow(-1, n) * rf(1. - i_eps + ess + nu, n)) /
                (gsl_sf_fact(-n) * rf(1. + i_eps - ess + nu, n) *
                 rf(2. + nu + nu, n)));
		f_sum += f_new * prefactor;

		re_err = fabs(real(f_new));
		im_err = fabs(imag(f_new));
		if (n == 1000) {
			printf("f_sum_dn did not converge.");
			return 0.;
		}
	}
	return f_sum;
}


complex<double> Knu(vector<complex<double>> &f_vec,
	vector<int> &n_vec, complex<double> nu, double lambda, double omega,
	double aa, int em, double ess=-2, double tol=1e-15) {

	double epsilon = 2 * omega;
	double kappa = sqrt(1 - aa * aa);
	double tau = (epsilon - em * aa) / kappa;
	double epsilonp = 0.5 * (tau + epsilon);

	complex<double> fSumKnu1Up, fSumKnu1Dn;
	complex<double> i {0, 1};

	fSumKnu1Up = f_sum_Knu_up(f_vec, n_vec, nu, lambda, omega, aa, em);
	// printf("fSumKnu1Up = %.17g + %.17gi \n", real(fSumKnu1Up), imag(fSumKnu1Up));
	fSumKnu1Dn = f_sum_Knu_dn(f_vec, n_vec, nu, lambda, omega, aa, em);
	// printf("fSumKnu1Dn = %.17g + %.17gi \n", real(fSumKnu1Dn), imag(fSumKnu1Dn));

	complex<double> res;
	// printf("check = %.17g + %.17gi \n", real(exp(i * epsilon * kappa)), imag(exp(i * epsilon * kappa)));
	res = ((exp(i * epsilon * kappa) * fSumKnu1Up *
            pow(epsilon * kappa, ess - nu) * gamma(2. + 2. * nu) *
            gamma(1. - ess - 2. * i * epsilonp)) /
            (pow(2., nu) * fSumKnu1Dn * gamma(1. + i * epsilon - ess + nu) *
            gamma(1. + i * epsilon + ess + nu) *
            gamma(1. + nu + i * tau)));

	return res;

}


complex<double> Aplus(vector<complex<double>> &f_vec,
	vector<int> &n_vec, complex<double> nu, double lambda, double omega,
	double aa, int em, double ess=-2, double tol=1e-15) {

	complex<double> f_up, f_dn, f_sum, res;
	complex<double> i {0, 1};

	double epsilon = 2 * omega;

	f_up = f_sum_up(f_vec, n_vec, nu, lambda, omega, aa, em); // good
	f_dn = f_sum_dn(f_vec, n_vec, nu, lambda, omega, aa, em);
	f_sum = f_up + f_dn;

	// printf("f_up = %.17g + %.17gi \n", real(f_up), imag(f_up));
	// printf("f_dn = %.17g + %.17gi \n", real(f_dn), imag(f_dn));

	res = ((pow(2, (-1. - i * epsilon + ess)) *
            exp(-(epsilon * M_PI) / 2. + 0.5 * i *
                   (1. - ess + nu) * M_PI) *
            gamma(1. + i * epsilon - ess + nu)) /
            gamma(1. - i * epsilon + ess + nu)) * f_sum;
	return res;
}


complex<double> csc(complex<double> x) {
	return 1. / sin(x);
}

// TODO: link and use this function in python.
// currently this is only programmed, but not used in the
// main python program
// NOTE: this is not going to be a massive speed up though
complex<double> Binc(vector<complex<double>> &f_vec,
	vector<int> &n_vec, vector<complex<double>> &f_neg_vec,
	vector<int> &n_neg_vec, complex<double> nu, double lambda, double omega,
	double aa, int em, double ess=-2, double tol=1e-15) {

	double epsilon = 2 * omega;
	double kappa = sqrt(1 - aa * aa);

	complex<double> A, Knu1, Knu2, res;
	complex<double> i {0, 1};

	A = Aplus(f_vec, n_vec, nu, lambda, omega, aa, em);
	// printf("A = %.17g + %.17gi \n", real(A), imag(A)); // good
	Knu1 = Knu(f_vec, n_vec, nu, lambda, omega, aa, em);
	// printf("Knu1 = %.17g + %.17gi \n", real(Knu1), imag(Knu1));

	// FIXME: this following function uses the old nu vector for f.
	// we need to make a second set of vectors or invert them or something.
	// there is a symmetry here which could be exploited which is not currently.
	Knu2 = Knu(f_neg_vec, n_neg_vec, -(1. + nu), lambda, omega, aa, em);
	// printf("Knu2 = %.17g + %.17gi \n", real(Knu2), imag(Knu2));
	// printf("A = %.17g + %.17gi \n", real(A), imag(A));

	// complex<double> check;
	complex<double> eps_c {epsilon, 0};
	// check = exp(i * epsilon * ((-1. + kappa) * 0.5 + log(epsilon))) * omega;
	// check = log(eps_c);

	res = A * ((Knu1 - (i * Knu2 * csc((-i * eps_c + ess + nu) *
                M_PI) *
                sin((i * eps_c - ess + nu) * M_PI)) /
        exp(i * nu * M_PI)) /
        (exp(i * eps_c * ((-1. + kappa) * 0.5 + log(eps_c))) * omega));
	// printf("check = %.17g + %.17gi \n", real(check), imag(check));
	// printf("Binc = %.17g + %.17gi \n", real(res), imag(res));
	return res;
}


// int main() {
// 	int em = 2;
// 	complex<double> nu {1.77980305023701826085742951665338224628827792617851607563424529604734, 0};
//     double aa = 0.998;
//     double epsilon = 0.6457703535574504596266418475864428068406741427634647973592090;
//     double omega = epsilon / 2;
//     double lambda = 1.87728029208475374768756357658320995667770697758914259921550427;
//     vector<complex<double>> f_vec, f_neg_vec;
//     vector<int> n_vec, n_neg_vec;
//     complex<double> B_inc, fn;

//     int n = -10;

//     // fn = find_fn(n, f_vec, n_vec, nu, lambda, omega, aa, em);
//     // printf("fn = %.17g + %.17gi \n", real(fn), imag(fn));
//     B_inc = Binc(f_vec, n_vec, f_neg_vec, n_neg_vec, nu, lambda, omega, aa, em);
//     printf("B_inc = %.17g + %.17gi \n", real(B_inc), imag(B_inc));
// }


















