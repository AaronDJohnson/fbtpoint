#include "sp_func.h"
#include <complex>
#include <vector>

using namespace std;

complex<double> F000(complex<double> a, complex<double> b, complex<double> c,
	complex<double> z, complex<double> Fm1m10, complex<double> Fm2m20) {

	complex<double> alpha, beta, gamma, ratio, prefactor, res;

	alpha = (-1. + a - c) * (-1. + b - c) * (-1. + a + b - c);
    beta = ((-2. + a + b - c) * (2. * (-1. + a) * (-1. + b) +
            (-1. + a) * (-1. + a) * z - (-1. + a) * c * (1. + z) +
            (-b + c) * (c + z - (-1. + b) * z)));
    gamma = ((-1. + a) * (-1. + b) * (-3. + a + b - c) * (-1. + z) * (-1. + z));

    ratio = (Fm1m10 / Fm2m20);

    prefactor = - 1. / gamma * (alpha - beta * ratio);

    res = prefactor * Fm2m20;
    return res;
}


complex<double> Fm2m2(complex<double> a, complex<double> b, complex<double> c,
	complex<double> z, complex<double> F, complex<double> Fm1m10) {

	complex<double> alpha, beta, gamma, ratio, prefactor, res;
    alpha = (-1. + a - c) * (-1. + b - c) * (-1. + a + b - c);
    beta = ((-2. + a + b - c) * (2. * (-1. + a) * (-1. + b) +
            (-1. + a) * (-1. + a) * z - (-1. + a) * c * (1. + z) +
            (-b + c) * (c + z - (-1. + b) * z)));
    gamma = ((-1. + a) * (-1. + b) * (-3. + a + b - c) * (-1. + z) * (-1. + z));

    ratio = (F / Fm1m10);

    prefactor = 1. / alpha * (beta - gamma * ratio);

    res = prefactor * Fm1m10;
    return res;
}


complex<double> F001(complex<double> a, complex<double> b, complex<double> c,
	complex<double> z, complex<double> Fm1m11, complex<double> Fm2m21) {
    /*
    Uses Gauss' contiguous relations to solve for 2F1(0, 0, 0)
    in terms of 2F1(-1, -1, 0) and 2F1(-2, -2, 0)
    */
    return ((-2. + a + b - c) * (((-2. + a - c) * (2. - b + c) *
            (Fm2m21)) / (-4. + a + b - c) -
           (3. + c + a * (-2. + z) - ((-1. + a) * (-1. + a - c) * (-1. + z)) /
            (-2. + a + b - c) +
            ((-2. + a) * (2. - a + c) * (-1. + z)) / (-4. + a + b - c) - b * z) *
        	Fm1m11)) / ((-1. + a) * (-1. + b) * pow(-1. + z, 2));
}


complex<double> Fm2m21(complex<double> a, complex<double> b, complex<double> c,
	complex<double> z, complex<double> F00p1, complex<double> Fm1m11) {
	complex<double> res;
    res = ((-4. + a + b - c) * ((-3. - c - a * (-2. + z) +
           ((-1. + a) * (-1. + a - c) * (-1. + z)) / (-2. + a + b - c) -
           ((-2. + a) * (2. - a + c) * (-1. + z)) / (-4. + a + b - c) + b * z) *
        	Fm1m11 -
          ((-1. + a) * (-1. + b) * pow(-1. + z, 2) * F00p1) /
        	(-2. + a + b - c))) / ((2. - a + c) * (2. - b + c));
    return res;
}


complex<double> hypgeo_up(int n, complex<double> a, complex<double> b,
	complex<double> c, complex<double> d, complex<double> hypgeo_m1,
	complex<double> hypgeo_m2) {
	double nf = n;
    return F000(a + nf, b + nf, c, d, hypgeo_m1, hypgeo_m2);
}


complex<double> hypgeo_dn(int n, complex<double> a, complex<double> b,
	complex<double> c, complex<double> d, complex<double> hypgeo_m1,
	complex<double> hypgeo_m2) {
	double nf = n;
    return Fm2m2(a + nf, b + nf, c, d, hypgeo_m1, hypgeo_m2);
}


complex<double> hypgeo(int nf, vector<complex<double>> &hyp_vec, 
	vector<int> &nf_vec, complex<double> a, complex<double> b,
	complex<double> c, complex<double> d) {

	vector<int>::iterator it = find(nf_vec.begin(), nf_vec.end(), nf);
	if (it != nf_vec.end()) {
		int index = distance(nf_vec.begin(), it);
		return hyp_vec[index];
	} else {
		if (nf == -1) {
			complex<double> res = hypgeo_2f1(a - 1., b - 1., c, d);
			nf_vec.insert(nf_vec.begin(), nf);
			hyp_vec.insert(hyp_vec.begin(), res);
			return res;
	    } else if (nf == -2) {
	    	complex<double> res = hypgeo_2f1(a - 2., b - 2., c, d);
			nf_vec.insert(nf_vec.begin(), nf);
			hyp_vec.insert(hyp_vec.begin(), res);
	    	return res;
	    } else if (nf >= 0) {
	    	complex<double> hyp_m1 = hyp_vec[hyp_vec.size() - 1];
	    	complex<double> hyp_m2 = hyp_vec[hyp_vec.size() - 2];
	    	complex<double> res = hypgeo_up(nf, a, b, c, d, hyp_m1, hyp_m2);
	    	nf_vec.push_back(nf);
	    	hyp_vec.push_back(res);
	        return res;
	    } else {
	    	complex<double> hyp_m1 = hyp_vec[1];
	    	complex<double> hyp_m2 = hyp_vec[0];
	    	complex<double> res = hypgeo_dn(nf + 2, a, b, c, d, hyp_m1, hyp_m2);
	    	nf_vec.insert(nf_vec.begin(), nf);
			hyp_vec.insert(hyp_vec.begin(), res);
	        return res;
	    }
	}
}


complex<double> find_hypgeo(int n, vector<complex<double>> &hyp_vec, 
	vector<int> &nf_vec, complex<double> a, complex<double> b,
	complex<double> c, complex<double> d) {

	// make sure the vectors contain -1 and -2:
	complex<double> hyp_m1, hyp_m2, hypgeo_new;
	hyp_m1 = hypgeo(-1, hyp_vec, nf_vec, a, b, c, d);
	hyp_m2 = hypgeo(-2, hyp_vec, nf_vec, a, b, c, d);
    if (n == -1) {
        return hyp_m1;
    } else if (n == -2) {
        return hyp_m2;
    } else if (n >= 0) {
        hypgeo_new = hypgeo(0, hyp_vec, nf_vec, a, b, c, d);
        for (int i=1; i<=n; i++) {
            hypgeo_new = hypgeo(i, hyp_vec, nf_vec, a, b, c, d);
        }
        return hypgeo_new;
    } else {
        for (int i=3; i<=-n; i++) {
            hypgeo_new = hypgeo(-i, hyp_vec, nf_vec, a, b, c, d);
        }
        return hypgeo_new;
    }
}


complex<double> hypgeo_deriv(int nf, vector<complex<double>> &hyp_vec, 
	vector<int> &nf_vec, complex<double> a, complex<double> b,
	complex<double> c, complex<double> d) {

	vector<int>::iterator it = find(nf_vec.begin(), nf_vec.end(), nf);
	if (it != nf_vec.end()) {
		int index = distance(nf_vec.begin(), it);
		return hyp_vec[index];
	} else {
		if (nf == -1) {
			complex<double> res = hypgeo_2f1(a, b, c + 1., d);
			nf_vec.insert(nf_vec.begin(), nf);
			hyp_vec.insert(hyp_vec.begin(), res);
			return res;
	    } else if (nf == -2) {
	    	complex<double> res = hypgeo_2f1(a - 1., b - 1., c + 1., d);
			nf_vec.insert(nf_vec.begin(), nf);
			hyp_vec.insert(hyp_vec.begin(), res);
	    	return res;
	    } else if (nf >= 0) {
	    	complex<double> hyp_m1 = hyp_vec[hyp_vec.size() - 1];
	    	complex<double> hyp_m2 = hyp_vec[hyp_vec.size() - 2];
	    	complex<double> res = hypgeo_up(nf, a + 1., b + 1., c + 1., d,
	    									hyp_m1, hyp_m2);
	    	nf_vec.push_back(nf);
	    	hyp_vec.push_back(res);
	        return res;
	    } else {
	    	complex<double> hyp_m1 = hyp_vec[1];
	    	complex<double> hyp_m2 = hyp_vec[0];
	    	complex<double> res = hypgeo_dn(nf + 2, a + 1., b + 1., c + 1., d,
	    									hyp_m1, hyp_m2);
	    	nf_vec.insert(nf_vec.begin(), nf);
			hyp_vec.insert(hyp_vec.begin(), res);
	        return res;
	    }
	}
}


complex<double> find_hypgeo_deriv(int n, vector<complex<double>> &hyp_vec, 
	vector<int> &nf_vec, complex<double> a, complex<double> b,
	complex<double> c, complex<double> d) {

	// make sure the vectors contain -1 and -2:
	complex<double> hyp_m1, hyp_m2, hypgeo_new;
	hyp_m1 = hypgeo_deriv(-1, hyp_vec, nf_vec, a, b, c, d);
	hyp_m2 = hypgeo_deriv(-2, hyp_vec, nf_vec, a, b, c, d);
    if (n == -1) {
        return hyp_m1;
    } else if (n == -2) {
        return hyp_m2;
    } else if (n >= 0) {
        hypgeo_new = hypgeo_deriv(0, hyp_vec, nf_vec, a, b, c, d);
        for (int i=1; i<=n; i++) {
            hypgeo_new = hypgeo_deriv(i, hyp_vec, nf_vec, a, b, c, d);
        }
        return hypgeo_new;
    } else {
        for (int i=3; i<=-n; i++) {
            hypgeo_new = hypgeo_deriv(-i, hyp_vec, nf_vec, a, b, c, d);
        }
        return hypgeo_new;
    }
}


// int main() {
//     int ess = -2;
//     int em = 2;
//     complex<double> nu = 1.779803050237018260857429516653382246288277926178516075634245296047;
//     double aa = 0.998;
//     double epsilon = 0.6457703535574504596266418475864428068406741427634647973592090;
//     double eigen = 1.87728029208475374768756357658320995667770697758914259921550427;
//     double kappa = sqrt(1 - aa * aa);
//     double r = 3;
//     double omega = epsilon / 2.;

//     double rp = 1 + kappa;
//     double tau = (epsilon - em * aa) / kappa;
//     double x = omega * (rp - r) / (epsilon * kappa);

//     complex<double> i {0, 1};
//     complex<double> a1 = i * epsilon * kappa;
//     complex<double> a2 = 0. -ess - 0.5 * i * (epsilon + tau);
//     complex<double> a3 = -nu - 1. + 0.5 * i * (epsilon + tau);
//     complex<double> c = 1. - ess - i * epsilon - i * tau;
//     complex<double> d = -x / (1. - x);

//     complex<double> a = nu + 1. - i * tau;
//     complex<double> b = 1. - ess + nu - i * epsilon;

//     vector<complex<double>> hyp_vec;
//     vector<int> nf_vec;

//     complex<double> h2f1;
//     h2f1 = find_hypgeo_deriv(10, hyp_vec, nf_vec, a, b, c, d);

//     printf("h2f1 = %.17g + %.17g i \n", real(h2f1), imag(h2f1));

// }










