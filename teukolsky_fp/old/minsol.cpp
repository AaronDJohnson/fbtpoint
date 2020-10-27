#include <complex>
#include <vector>

using namespace std;

complex<double> f_alpha(int n, double kappa, double tau, double epsilon,
	complex<double> nu, int ess=-2) {

	complex<double> i {0, 1};

	return ((i*epsilon*kappa*(1.0 - i*epsilon + nu + (n + ess) * 1.0)*
     		(1.0 + i*epsilon + nu + (n + ess) * 1.0)*(1.0 + n + nu + i*tau))/
   			((1.0 + n + nu)*(3.0 + 2.0*n + 2.0*nu)));
}


complex<double> f_beta(int n, double epsilon, complex<double> nu,
	double lambda, double aa, int em, int ess=-2) {

	double epsilon2 = epsilon * epsilon;
	double ess2 = ess * ess;

	return (epsilon2 - ess*(1 + ess) + (1.*n + nu)*(1. + n + nu) +
	epsilon*(epsilon - em*aa) + (epsilon*(epsilon2 + ess2)*(epsilon - em*aa))/
    ((1.*n + nu)*(1. + n + nu)) - lambda);
}


complex<double> f_gamma(int n, double kappa, double tau, double epsilon,
	complex<double> nu, int ess=-2) {

	complex<double> i {0, 1};

	return ((i*epsilon*kappa*(-i*epsilon - 1.*ess + 1.*n + nu)*
     		(i*epsilon - 1.*ess + 1.*n + nu)*(1.*n + nu - i*tau))/
   			((1.*n + nu)*(-1. + 2.*n + 2.*nu)));
}


complex<double> num_up(int i, double kappa, double tau, double epsilon,
	complex<double> nu) {

	complex<double> alpha, gamma;
	alpha = f_alpha(i - 1, kappa, tau, epsilon, nu);
	gamma = f_gamma(i, kappa, tau, epsilon, nu);
	return -alpha * gamma;
}


complex<double> den_up(int i, double epsilon, complex<double> nu,
	double lambda, double aa, double em) {

	complex<double> beta;

	beta = f_beta(i, epsilon, nu, lambda, aa, em);
	return beta;
}


complex<double> f_cont_frac_up(int n0, double kappa, double tau,
	double lambda, double epsilon, complex<double> nu, double aa, int em,
	int ess=-2, double tol=1e-15) {
	// printf("n0 = %d \n", n0);
	int j;
	complex<double> ak, bk, Am2, Bm2, Am1, Bm1, A, B, x, xm1;
	double re_err, im_err;
	xm1 = 1e30;  // large value to keep the function from converging first iter

	ak = num_up(n0, kappa, tau, epsilon, nu);
	bk = den_up(n0, epsilon, nu, lambda, aa, em);

	Am2 = 1;
	Bm2 = 0;
	Am1 = 0;
	Bm1 = 1;
	A = bk * Am1 + ak * Am2;
	B = bk * Bm1 + ak * Bm2;
	x = A / B;
	// printf("A = %.17g \n", A);
	// printf("B = %.17g \n", B);

	j = n0;
	re_err = fabs(real(x) - real(xm1)) / fabs(real(xm1));
	im_err = fabs(imag(x) - imag(xm1)) / fabs(imag(xm1));
	// printf("rel_err = %.17g \n", rel_err);
	while (re_err > tol or im_err > tol) {
		// update values
		xm1 = x;
		Am2 = Am1;
		Am1 = A;
		Bm2 = Bm1;
		Bm1 = B;

		// find new values
		ak = num_up(j + 1, kappa, tau, epsilon, nu);
		bk = den_up(j + 1, epsilon, nu, lambda, aa, em);

		A = bk * Am1 + ak * Am2;
		B = bk * Bm1 + ak * Bm2;
		x = A / B;
		// printf("x = %.17g \n", x);

		j++;
		// printf("j = %d \n", j);
		// printf("rel_err = %.17g \n", rel_err);
		re_err = fabs(real(x) - real(xm1)) / fabs(real(xm1));
		im_err = fabs(imag(x) - imag(xm1)) / fabs(imag(xm1));
		if (j == 10000) {
			// infinite loop prevention
			printf("Continued fraction did not converge.\n");
			return 0.;
		}
	}

	return x;
}


vector<complex<double>> find_minsol_up(double kappa, double tau, double lambda,
	double epsilon, complex<double> nu, double aa, int em, double tol=1e-15) {
	vector<complex<double>> f_up;

	complex<double> one {1, 0}, res, R, alpha;
	double re_err, im_err;

	f_up.push_back(one);
	res = 1e30;  // large number so that loop doesn't converge first iteration
	int n0 = 1;
	re_err = fabs(real(f_up[n0-1]) - real(res)) / real(f_up[n0-1]);
	im_err = fabs(imag(f_up[n0-1]) - imag(res)) / imag(f_up[n0-1]);
	while (re_err > tol or im_err > tol) {
		res = f_up[n0 - 1];
		alpha = f_alpha(n0, kappa, tau, epsilon, nu);
		R = f_cont_frac_up(n0, kappa, tau, lambda, epsilon, nu, aa, em);
		f_up.push_back(R * f_up[n0 - 1] / alpha);

		re_err = fabs(real(f_up.back()) - real(res)) / real(f_up.back());
		im_err = fabs(imag(f_up.back()) - imag(res)) / imag(f_up.back());
		n0++;
		if (n0 == 10000) {
			printf("minsol_up did not converge.");
			return f_up;
		}
	}
	return f_up;
}



int main() {
	double kappa, tau, lambda, epsilon, aa;
	int em = 2;
	complex<double> nu {-2.779803050237018260857429516653382246288277, 0};
	aa = 0.998;
	epsilon = 0.6457703535574504596266418475864428068406741427634647973592090;
	lambda = 1.87728029208475374768756357658320995667770697758914259921550427;
	kappa = sqrt(1 - aa*aa);
	tau = (epsilon - em * aa) / kappa;
	vector<complex<double>> f_up;
	f_up = find_minsol_up(kappa, tau, lambda, epsilon, nu, aa, em);
	for (int i=0; i<=f_up.size(); i++) {
		printf("f_up[%d] = (%.17g, %.17gi) \n", i, real(f_up[i]), imag(f_up[i]));
	}











}









