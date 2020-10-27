#include "hypgeo.h"
#include "fminsol.h"
#include "amplitudes.h"

#include <complex>
#include <vector>

using namespace std;

//TODO: ess and tol are not listed in some of these parts
// go back through and make sure even KWargs are listed in every
// function call

complex<double> R_2F1(vector<complex<double>> hyp_vec, vector<int> nf_vec,
	vector<complex<double>> f_vec, vector<int> n_vec, complex<double> &y,
	double r, complex<double> nu, double lambda, double aa,
	double omega, int em, int ess=-2, double tol=1e-15) {
    /*
    B.2 in Throwe's thesis
    */
    double kappa = sqrt(1 - aa*aa);
    double rp = 1 + kappa;
    double epsilon = 2 * omega;
    double tau = (epsilon - em * aa) / kappa;
    double x = omega * (rp - r) / (epsilon * kappa);

    complex<double> i {0, 1};

    // exponents:
    complex<double> a1 = i * epsilon * kappa * x;
    complex<double> a2 = 0. -ess - i * (epsilon + tau) / 2.;
    complex<double> a3 = - 1. + i * (epsilon + tau) / 2. - nu;

    complex<double> c = 1. - ess - i * epsilon - i * tau;
    complex<double> d = -x / (1 - x);
    complex<double> prefactor = exp(a1) * pow(-x, a2) * pow(1 - x, a3);

    complex<double> a = 1. + nu - i * tau;
    complex<double> b = 1. - ess + nu - i * epsilon;

    complex<double> hypgeo, f, y0, y_up;

    double re_err, im_err;

    int n = 0;
    hypgeo = find_hypgeo(n, hyp_vec, nf_vec, a, b, c, d);
    f = find_fn(n, f_vec, n_vec, nu, lambda, omega, aa, em);
    y0 = hypgeo * pow(1 - x, -n) * f;

    n++;
    hypgeo = find_hypgeo(n, hyp_vec, nf_vec, a, b, c, d);
    f = find_fn(n, f_vec, n_vec, nu, lambda, omega, aa, em);
    y_up = hypgeo * pow(1 - x, -n) * f + y0;

    re_err = fabs(real(y_up) - real(y0)) / fabs(real(y0));
    im_err = fabs(imag(y_up) - imag(y0)) / fabs(imag(y0));

    // up
    while (re_err > tol or im_err > tol) {
        n++;
        hypgeo = find_hypgeo(n, hyp_vec, nf_vec, a, b, c, d);
    	f = find_fn(n, f_vec, n_vec, nu, lambda, omega, aa, em);
        y0 = y_up;

        y_up += hypgeo * pow(1 - x, -n) * f;

        re_err = fabs(real(y_up) - real(y0)) / fabs(real(y0));
        im_err = fabs(imag(y_up) - imag(y0)) / fabs(imag(y0));
        if (n == 1000) {
            printf("R_2F1_up did not converge \n");
            return 0;
        }
    }
    // down

    complex<double> y_dn;

    n = -1;
    hypgeo = find_hypgeo(n, hyp_vec, nf_vec, a, b, c, d);
    f = find_fn(n, f_vec, n_vec, nu, lambda, omega, aa, em);
    y0 = hypgeo * pow(1 - x, -n) * f;

    n--;
    hypgeo = find_hypgeo(n, hyp_vec, nf_vec, a, b, c, d);
    f = find_fn(n, f_vec, n_vec, nu, lambda, omega, aa, em);
    y_dn = hypgeo * pow(1 - x, -n) * f + y0;

    re_err = fabs(real(y_dn) - real(y0)) / fabs(real(y0));
    im_err = fabs(imag(y_dn) - imag(y0)) / fabs(imag(y0));
    while (re_err > tol or im_err > tol) {
        n--;
        hypgeo = find_hypgeo(n, hyp_vec, nf_vec, a, b, c, d);
    	f = find_fn(n, f_vec, n_vec, nu, lambda, omega, aa, em);
        y0 = y_dn;

        y_dn += hypgeo * pow(1 - x, -n) * f;

        re_err = fabs(real(y_dn) - real(y0)) / fabs(real(y0));
        im_err = fabs(imag(y_dn) - imag(y0)) / fabs(imag(y0));
        if (n == 1000) {
            printf("R_2F1_up did not converge \n");
            return 0;
        }
    }
    y = y_up + y_dn;
    complex<double> Rin = prefactor * y;

    // printf("Rin = %.17g + i %.17g \n", real(Rin), imag(Rin));
    return Rin;
}


complex<double> dR_2F1dr(vector<complex<double>> hyp_vec,
	vector<int> nf_vec,
	vector<complex<double>> hyp_der_vec,
	vector<int> nf_der_vec,
	vector<complex<double>> f_vec, vector<int> n_vec,
	complex<double> &Rin,
	double r, complex<double> nu, double lambda, double aa,
	double omega, int em, int ess=-2, double tol=1e-15) {

    double kappa = sqrt(1 - aa*aa);
    double rp = 1 + kappa;
    double epsilon = 2 * omega;
    double tau = (epsilon - em * aa) / kappa;
    double x = omega * (rp - r) / (epsilon * kappa);

    complex<double> i {0, 1};

    // exponents:
    complex<double> a1 = i * epsilon * kappa;
    complex<double> a2 = 0. -ess - i * (epsilon + tau) / 2.;
    complex<double> a3 = -1. -nu + i * (epsilon + tau) / 2.;

    complex<double> c = 1. - ess - i * epsilon - i * tau;
    complex<double> d = -x / (1 - x);
    complex<double> a = 1. - i * tau + nu;
    complex<double> b = 1. - ess - i * epsilon + nu;

    complex<double> yr;

    Rin = R_2F1(hyp_vec, nf_vec, f_vec, n_vec, yr, r, nu, lambda, aa,
    			omega, em);
    // printf("Rin = %.17g + i %.17g \n", real(Rin), imag(Rin));

    complex<double> prefactor, hypgeo, hypgeo_deriv, f, y0, y_up;
    double re_err, im_err;
    
    int n = 0;
    hypgeo = find_hypgeo(n, hyp_vec, nf_vec, a, b, c, d);
    hypgeo_deriv = find_hypgeo_deriv(n, hyp_der_vec, nf_der_vec, a, b, c, d);
    prefactor = (-1. / c * pow(1 - x, -n - 2) *
            	(c * (x - 1.) * (0. + n) * hypgeo +
            	(0. + n + a) * (0. + n + b) * hypgeo_deriv));

    f = find_fn(n, f_vec, n_vec, nu, lambda, omega, aa, em);
    y0 = prefactor * f;

    // printf("y0 = %.17g + i %.17g \n", real(y0), imag(y0));

    n++;
    hypgeo = find_hypgeo(n, hyp_vec, nf_vec, a, b, c, d);
    hypgeo_deriv = find_hypgeo_deriv(n, hyp_der_vec, nf_der_vec, a, b, c, d);
    prefactor = (-1. / c * pow(1 - x, -n - 2) *
            	(c * (x - 1.) * (0. + n) * hypgeo +
            	(0. + n + a) * (0. + n + b) * hypgeo_deriv));
    f = find_fn(n, f_vec, n_vec, nu, lambda, omega, aa, em);
    y_up = prefactor * f + y0;

    re_err = fabs(real(y_up) - real(y0)) / fabs(real(y0));
    im_err = fabs(imag(y_up) - imag(y0)) / fabs(imag(y0));

    // up
    while (re_err > tol or im_err > tol) {
        n++;
        hypgeo = find_hypgeo(n, hyp_vec, nf_vec, a, b, c, d);
	    hypgeo_deriv = find_hypgeo_deriv(n, hyp_der_vec, nf_der_vec, a, b, c, d);
	    prefactor = (-1. / c * pow(1 - x, -n - 2) *
	            	(c * (x - 1.) * (0. + n) * hypgeo +
	            	(0. + n + a) * (0. + n + b) * hypgeo_deriv));
        f = find_fn(n, f_vec, n_vec, nu, lambda, omega, aa, em);
        y0 = y_up;

        y_up += prefactor * f;

        re_err = fabs(real(y_up) - real(y0)) / fabs(real(y0));
   		im_err = fabs(imag(y_up) - imag(y0)) / fabs(imag(y0));
        if (n == 1000) {
            printf("R_2F1_up did not converge \n");
            return 0;
        }
    }
    // printf("y_up = %.17g + i %.17g \n", real(y_up), imag(y_up));
    // down
    complex<double> y_dn;

    n = -1;
    hypgeo = find_hypgeo(n, hyp_vec, nf_vec, a, b, c, d);
    hypgeo_deriv = find_hypgeo_deriv(n, hyp_der_vec, nf_der_vec, a, b, c, d);
    prefactor = (-1. / c * pow(1 - x, -n - 2) *
            	(c * (x - 1.) * (0. + n) * hypgeo +
            	(0. + n + a) * (0. + n + b) * hypgeo_deriv));
    f = find_fn(n, f_vec, n_vec, nu, lambda, omega, aa, em);
    y0 = prefactor * f;
    // printf("y0 = %.17g + i %.17g \n", real(y0), imag(y0));

    n--;
    hypgeo = find_hypgeo(n, hyp_vec, nf_vec, a, b, c, d);
    hypgeo_deriv = find_hypgeo_deriv(n, hyp_der_vec, nf_der_vec, a, b, c, d);
    prefactor = (-1. / c * pow(1 - x, -n - 2) *
            	(c * (x - 1.) * (0. + n) * hypgeo +
            	(0. + n + a) * (0. + n + b) * hypgeo_deriv));
    f = find_fn(n, f_vec, n_vec, nu, lambda, omega, aa, em);
    y_dn = prefactor * f + y0;
    // printf("y_dn = %.17g + i %.17g \n", real(y_dn), imag(y_dn));

    re_err = fabs(real(y_dn) - real(y0)) / fabs(real(y0));
   	im_err = fabs(imag(y_dn) - imag(y0)) / fabs(imag(y0));
    while (re_err > tol or im_err > tol) {
        n--;
        // printf("%i \n", n);
        hypgeo = find_hypgeo(n, hyp_vec, nf_vec, a, b, c, d);
	    hypgeo_deriv = find_hypgeo_deriv(n, hyp_der_vec, nf_der_vec, a, b, c, d);
	    prefactor = (-1. / c * pow(1 - x, -n - 2) *
	            	(c * (x - 1.) * (0. + n) * hypgeo +
	            	(0. + n + a) * (0. + n + b) * hypgeo_deriv));
	    f = find_fn(n, f_vec, n_vec, nu, lambda, omega, aa, em);
	    y0 = y_dn;

        y_dn += prefactor * f;
        // printf("y_dn = %.17g + i %.17g \n", real(y_dn), imag(y_dn));

        re_err = fabs(real(y_dn) - real(y0)) / fabs(real(y0));
   		im_err = fabs(imag(y_dn) - imag(y0)) / fabs(imag(y0));
        if (n == 1000) {
            printf("R_2F1_up did not converge \n");
            return 0;
        }
    }
    // printf("y_dn = %.17g + i %.17g \n", real(y_dn), imag(y_dn));

    complex<double> y, logR_prime, Rdin;
    y = y_up + y_dn;

    logR_prime = a1 + (a2 / x) + (a3 / (x - 1)) + (y / yr);

    Rdin = -omega / (epsilon * kappa) * logR_prime * Rin;
    // printf("Rdin = %.17g + i %.17g \n", real(Rdin), imag(Rdin));
    return Rdin;
}


complex<double> potential(double r, double aa, int em, double omega,
	double lambda, int ess=-2, double M=1.) {
    double delta = r * r - 2 * M * r + aa * aa;
    double K = -aa * em + (aa * aa + r * r) * omega;

    complex<double> i {0, 1};

    double ess_d = ess;
    complex<double> V;

    V = (lambda -
         4. * i * r * ess_d * omega -
         (-2. * i * (-M + r) * ess_d * K + K * K) / delta
         );
    return V;
}


complex<double> dR_2F1dr2(vector<complex<double>> hyp_vec,
	vector<int> nf_vec,
	vector<complex<double>> hyp_der_vec,
	vector<int> nf_der_vec,
	vector<complex<double>> f_vec, vector<int> n_vec,
	complex<double> &Rin, complex<double> &Rdin,
	double r, complex<double> nu, double lambda, double aa,
	double omega, int em, int ess=-2, double tol=1e-15, double M=1.) {
    /*
    Uses the Teukolsky equation to find the second derivative
    */
    complex<double> Rdd;
    double delta = r * r - 2 * M * r + aa * aa;
    double delta_prime = 2 * r - 2 * M;
    complex<double> V = potential(r, aa, em, omega, lambda);
    Rdin = dR_2F1dr(hyp_vec, nf_vec,
				  hyp_der_vec, nf_der_vec,
				  f_vec, n_vec, Rin, r, nu, lambda, aa, omega, em, ess, tol);
    Rdd = (Rin * V + Rdin * delta_prime) / delta;

    // printf("Rdd = %.17g + i %.17g \n", real(Rdd), imag(Rdd));
    return Rdd;
}


void find_R(double &re_R_in, double &im_R_in, double &re_Rd_in,
	double &im_Rd_in, double &re_Rdd_in, double &im_Rdd_in,
	double r, double re_nu, double im_nu, double eigen, double aa,
	double omega, int em, int ess=-2, double tol=1e-15, double M=1.) {

	double lambda = eigen; // python doesn't like lambda

	complex<double> nu {re_nu, im_nu};
	complex<double> R_in, Rd_in, Rdd_in;
	vector<complex<double>> hyp_vec, hyp_der_vec, f_vec;
	vector<int> nf_vec, nf_der_vec, n_vec;
    Rdd_in = dR_2F1dr2(hyp_vec, nf_vec, hyp_der_vec, nf_der_vec, f_vec, n_vec,
    				   R_in, Rd_in, r, nu, lambda, aa, omega, em);
    // printf("R_in = %.17g + i%.17g \n", real(R_in), imag(R_in));
    re_R_in = real(R_in);
    // printf("re_R_in = %.17g \n", re_R_in);
    im_R_in = imag(R_in);
    // printf("im_R_in = %.17g \n", im_R_in);
    re_Rd_in = real(Rd_in);
    im_Rd_in = imag(Rd_in);
    re_Rdd_in = real(Rdd_in);
    im_Rdd_in = imag(Rdd_in);
    // printf("Rdd_in = %.17g + i%.17g \n", real(Rdd_in), imag(Rdd_in));
}


void find_Bin(double &re_B_in, double &im_B_in, double re_nu, double im_nu,
	double eigen, double aa, double omega,
	int em, int ess=-2, double tol=1e-15) {

	double lambda = eigen;
	complex<double> nu {re_nu, im_nu};
	complex<double> B_in;
	vector<complex<double>> f_vec, f_neg_vec;
	vector<int> n_vec, n_neg_vec;

	B_in = Binc(f_vec, n_vec, f_neg_vec, n_neg_vec, nu, lambda, omega, aa, em);

	re_B_in = real(B_in);
	im_B_in = imag(B_in);

}


// int main() {
// 	int em = 2;
//     complex<double> nu = 1.783400403135214919625587973236499545319854788778001927189583309161564163708964633825665147355791331;
//     double re_nu = real(nu);
//     double im_nu = imag(nu);
//     double aa = 0.998;
//     double epsilon = 2 * 0.3209735919183500952842127182168827323371877619227890953225746693737469758370741001675147452797702759;
//     double eigen = 1.889708166284304609494898318227660972426238914985137355645957343847354541537472115043317557726883335;
//     double kappa = sqrt(1 - aa * aa);
//     double r = 2.727272727272727272727272727272727272727272727272727272727272727272727272727272727272727272727272727;
//     double omega = epsilon / 2.;

//     // double re_R_in, im_R_in, re_Rd_in, im_Rd_in, re_Rdd_in, im_Rdd_in;
    
//     // find_R(re_R_in, im_R_in, re_Rd_in, im_Rd_in, re_Rdd_in, im_Rdd_in, r, re_nu, im_nu, eigen, aa, omega, em);

//     double re_B_in, im_B_in;
//     find_Bin(re_B_in, im_B_in, re_nu, im_nu, eigen, aa, omega, em);
//     printf("B_in = %.17g + i%.17g \n", re_B_in, im_B_in);
// }








