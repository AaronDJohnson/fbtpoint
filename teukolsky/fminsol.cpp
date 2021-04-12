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

	return (-(i*epsilon*kappa*(-i*epsilon - 1.*ess + 1.*n + nu)*
     		(i*epsilon - 1.*ess + 1.*n + nu)*(1.*n + nu - i*tau))/
   			((1.*n + nu)*(-1. + 2.*n + 2.*nu)));
}


complex<double> numer_up(int n, double kappa, double tau, double epsilon,
	complex<double> nu, int ess=-2) {

	complex<double> alpha, gamma;
	alpha = f_alpha(n - 1, kappa, tau, epsilon, nu, ess);
	gamma = f_gamma(n, kappa, tau, epsilon, nu, ess);
	return -alpha * gamma;
}


complex<double> numer_dn(int n, double kappa, double tau, double epsilon,
	complex<double> nu, int ess=-2) {

	complex<double> alpha, gamma;
	alpha = f_alpha(n, kappa, tau, epsilon, nu, ess);
	gamma = f_gamma(n + 1, kappa, tau, epsilon, nu, ess);
	return -alpha * gamma;
}


complex<double> denom(int n, double epsilon, complex<double> nu,
	double lambda, double aa, double em) {

	complex<double> beta;

	beta = f_beta(n, epsilon, nu, lambda, aa, em);
	return beta;
}


complex<double> cont_frac_up(int n0, double kappa, double tau,
	double lambda, double epsilon, complex<double> nu, double aa, int em,
	int ess=-2, double tol=1e-15) {
	// printf("n0 = %d \n", n0);
	int j;
	complex<double> ak, bk, Am2, Bm2, Am1, Bm1, A, B, res;
	double re_err, im_err;
	// large value to keep the function from converging first iter
	complex<double> res_old {1e30, 1e30};

	ak = numer_up(n0, kappa, tau, epsilon, nu);
	// printf("ak = %.17g + %.17gi \n", real(ak), imag(ak));
	bk = denom(n0, epsilon, nu, lambda, aa, em);
	// printf("bk = %.17g + %.17gi \n", real(bk), imag(bk));

	Am2 = 1;
	Bm2 = 0;
	Am1 = 0;
	Bm1 = 1;
	A = bk * Am1 + ak * Am2;
	B = bk * Bm1 + ak * Bm2;
	res = A / B;
	// printf("A = %.17g \n", A);
	// printf("B = %.17g \n", B);

	j = n0 + 1;
	re_err = fabs(real(res) - real(res_old));
	im_err = fabs(imag(res) - imag(res_old));
	// printf("rel_err = %.17g \n", rel_err);
	while (re_err > tol or im_err > tol) {
		// update values
		res_old = res;
		Am2 = Am1;
		Am1 = A;
		Bm2 = Bm1;
		Bm1 = B;

		// find new values
		ak = numer_up(j, kappa, tau, epsilon, nu);
		// printf("ak = %.17g + %.17gi \n", real(ak), imag(ak));
		bk = denom(j, epsilon, nu, lambda, aa, em);
		// printf("bk = %.17g + %.17gi \n", real(bk), imag(bk));
		A = bk * Am1 + ak * Am2;
		B = bk * Bm1 + ak * Bm2;
		res = A / B;
		// printf("x = %.17g \n", x);
		// printf("j = %d \n", j);
		// printf("rel_err = %.17g \n", rel_err);
		re_err = fabs(real(res) - real(res_old));
		im_err = fabs(imag(res) - imag(res_old));
		j++;
		if (j == 10000) {
			// infinite loop prevention
			printf("Continued fraction did not converge.\n");
			return 0.;
		}
	}

	return res;
}


complex<double> cont_frac_dn(int n0, double kappa, double tau,
	double lambda, double epsilon, complex<double> nu, double aa, int em,
	int ess=-2, double tol=1e-15) {
	// printf("n0 = %d \n", n0);
	int j;
	complex<double> ak, bk, Am2, Bm2, Am1, Bm1, A, B, res;
	double re_err, im_err;
	// large value to keep the function from converging first iter
	complex<double> res_old {1e30, 1e30};

	ak = numer_dn(n0, kappa, tau, epsilon, nu);
	bk = denom(n0, epsilon, nu, lambda, aa, em);

	Am2 = 1;
	Bm2 = 0;
	Am1 = 0;
	Bm1 = 1;
	A = bk * Am1 + ak * Am2;
	B = bk * Bm1 + ak * Bm2;
	res = A / B;
	// printf("A = %.17g \n", A);
	// printf("B = %.17g \n", B);

	j = n0 - 1;
	re_err = fabs(real(res) - real(res_old));
	im_err = fabs(imag(res) - imag(res_old));
	// printf("rel_err = %.17g \n", rel_err);
	while (re_err > tol or im_err > tol) {
		// update values
		res_old = res;
		Am2 = Am1;
		Am1 = A;
		Bm2 = Bm1;
		Bm1 = B;

		// find new values
		ak = numer_dn(j, kappa, tau, epsilon, nu);
		bk = denom(j, epsilon, nu, lambda, aa, em);

		A = bk * Am1 + ak * Am2;
		// printf("A = %.17g + %.17gi \n", real(A), imag(A));
		B = bk * Bm1 + ak * Bm2;
		// printf("B = %.17g + %.17gi \n", real(B), imag(B));
		res = A / B;
		// printf("res = %.17g + %.17gi \n", real(res), imag(res));
		// printf("x = %.17g \n", x);
		// printf("j = %d \n", j);
		// printf("rel_err = %.17g \n", rel_err);
		re_err = fabs(real(res) - real(res_old));
		// printf("re_err = %.17g + %.17gi \n", real(re_err), imag(re_err));
		im_err = fabs(imag(res) - imag(res_old));
		// printf("im_err = %.17g + %.17gi \n", real(im_err), imag(im_err));
		j--;
		if (j == 10000) {
			// infinite loop prevention
			printf("Continued fraction did not converge.\n");
			return 0.;
		}
	}

	return res;
}


complex<double> minsol(int nf, vector<complex<double> > &f_vec,
	vector<int> &n_vec, double kappa, double tau, double lambda, double epsilon,
	complex<double> nu, double aa, int em, int ess=-2, double tol=1e-15) {
	// the input vector must have nf = 0 when minsol = 1
	vector<int>::iterator it = find(n_vec.begin(), n_vec.end(), nf);
	if (it != n_vec.end()) {
		// this will check to see if the vector already has the result
		// would it be faster just to recompute everything?
		// Or perhaps we only need to check the last value
		// because everything in between will have been computed.
		// There are definitely potential speed increases here.
		int index = distance(n_vec.begin(), it);
		return f_vec[index];
	} else {
		if (nf == 0) {
			n_vec.push_back(0);
			f_vec.push_back(1); // check that this ends up in the right location
			return 1;
		} else if (nf > 0) {
			complex<double> frac_up, f_last, alpha, f;
			n_vec.push_back(nf);
			f_last = f_vec.back();
			frac_up = cont_frac_up(nf, kappa, tau, lambda, epsilon, nu, aa, em);
			// printf("frac_up = %.17g + %.17gi \n", real(frac_up), imag(frac_up));
			alpha = f_alpha(nf - 1, kappa, tau, epsilon, nu);
			f = f_last * frac_up / alpha;
			f_vec.push_back(f);
			return f;
		} else {
			complex<double> frac_dn, f_last, gamma, f;
			n_vec.insert(n_vec.begin(), nf);
			f_last = f_vec.front(); // 0th term
			// printf("f_last = %.17g + %.17gi \n", real(f_last), imag(f_last));
			frac_dn = cont_frac_dn(nf, kappa, tau, lambda, epsilon, nu, aa, em);
			gamma = f_gamma(nf + 1, kappa, tau, epsilon, nu);
			// printf("frac_dn = %.17g + %.17gi \n", real(frac_dn), imag(frac_dn));
			f = f_last * frac_dn / gamma;
			f_vec.insert(f_vec.begin(), f);
			return f;
		}
	}
}


complex<double> find_fn(int n, vector<complex<double> > &f_vec,
	vector<int> &n_vec, complex<double> nu, double lambda, double omega,
	double aa, int em, int ess=-2) {
	//
	double epsilon = 2 * omega;
	double kappa = sqrt(1 - aa * aa);
	double tau = (epsilon - em * aa) / kappa;
	complex<double> f;
	// printf("n = %i \n", n);
	if (n == 0) {
		f = minsol(n, f_vec, n_vec, kappa, tau, lambda,
				   epsilon, nu, aa, em);
		return f;
	} else if (n > 0) {
		for (int i = 0; i <= n; i++) {
			// printf("%i \n", i);
			f = minsol(i, f_vec, n_vec, kappa, tau, lambda,
				   	   epsilon, nu, aa, em);
			// printf("f = %.17g + %.17gi \n", real(f), imag(f));
		}
		return f;
	} else {
		for (int i = 0; i <= abs(n); i++) {
			// printf("%i \n", i);
			f = minsol(-i, f_vec, n_vec, kappa, tau, lambda,
				   	   epsilon, nu, aa, em);
			// printf("f = %.17g + %.17gi \n", real(f), imag(f));
		}
		return f;
	}

}


// int main() {
// 	double kappa, tau, lambda, epsilon, aa;
// 	int em = 2;
// 	int n = -3;
// 	complex<double> nu {-2.779803050237018260857429516653382246288277, 0};
// 	aa = 0.998;
// 	epsilon = 0.6457703535574504596266418475864428068406741427634647973592090;
// 	lambda = 1.87728029208475374768756357658320995667770697758914259921550427;
// 	kappa = sqrt(1 - aa*aa);
// 	tau = (epsilon - em * aa) / kappa;
// 	vector<complex<double> > f_up;
// 	// f_up = find_minsol_up(kappa, tau, lambda, epsilon, nu, aa, em);
// 	complex<double> num;
// 	vector<complex<double>> f_vec;
// 	vector<int> n_vec;
// 	double omega = 0.5 * epsilon;
// 	num = find_fn(n, f_vec, n_vec, nu, lambda, omega, aa, em);
// 	//num = find_fn(n + 1, f_vec, n_vec, nu, lambda, omega, aa, em);
// 	printf("num = %.17g + %.17gi \n", real(num), imag(num));


// }



