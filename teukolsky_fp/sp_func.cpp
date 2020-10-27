#include "arb.h"
#include "acb.h"
#include "acb_hypgeom.h"
#include <complex>

using namespace std;

complex<double> gamma(complex<double> x) {
	double res_real, res_imag;
	// break x into real and imag parts:
	double re_x = real(x);
	double im_x = imag(x);
	// printf("re_x = %.17g \n", re_x);
	// initialize special types:
	acb_t arb_input, arb_result;
	slong prec = 55; // floating point precision (1e-17)
	// set values to zero and allocate memory
	acb_init(arb_input);
	acb_init(arb_result);
	// set values to double values
	acb_set_d_d(arb_input, re_x, im_x);
	// compute gamma function
	acb_gamma(arb_result, arb_input, prec);
	// acb_printn(ay, 16, 0);
	// printf("\n");

	// allocate result values
	arb_t re_y, im_y;
	arb_init(re_y);
	arb_init(im_y);

	acb_get_real(re_y, arb_result);
	acb_get_imag(im_y, arb_result);

	// arb_printn(re_y, 16, 0);
	// printf("\n");

	res_real = arf_get_d(arb_midref(re_y), ARF_RND_NEAR);
	res_imag = arf_get_d(arb_midref(im_y), ARF_RND_NEAR);

	acb_clear(arb_input);
	acb_clear(arb_result);
	arb_clear(re_y);
	arb_clear(im_y);

	complex<double> res {res_real, res_imag};
	return res;
}


complex<double> rf(complex<double> x, complex<double> n) {
	double res_real, res_imag;

	double re_x = real(x);
	double im_x = imag(x);

	double re_n = real(n);
	double im_n = imag(n);

	acb_t arb_input_x, arb_input_n, arb_result;
	slong prec = 55;

	acb_init(arb_input_x);
	acb_init(arb_input_n);
	acb_init(arb_result);

	acb_set_d_d(arb_input_x, re_x, im_x);
	acb_set_d_d(arb_input_n, re_n, im_n);

	acb_rising(arb_result, arb_input_x, arb_input_n, prec);

	arb_t re_res, im_res;
	arb_init(re_res);
	arb_init(im_res);

	acb_get_real(re_res, arb_result);
	acb_get_imag(im_res, arb_result);

	res_real = arf_get_d(arb_midref(re_res), ARF_RND_NEAR);
	res_imag = arf_get_d(arb_midref(im_res), ARF_RND_NEAR);

	acb_clear(arb_input_x);
	acb_clear(arb_input_n);
	acb_clear(arb_result);
	arb_clear(re_res);
	arb_clear(im_res);

	complex<double> res {res_real, res_imag};
	return res;
}


complex<double> hypgeo_2f1(complex<double> a, complex<double> b,
	complex<double> c, complex<double> z) {

	double res_real, res_imag;

	double re_a = real(a);
	double im_a = imag(a);

	double re_b = real(b);
	double im_b = imag(b);

	double re_c = real(c);
	double im_c = imag(c);

	double re_z = real(z);
	double im_z = imag(z);


	acb_t arb_input_a, arb_input_b, arb_input_c, arb_input_z, arb_result;
	slong prec = 55;

	acb_init(arb_input_a);
	acb_init(arb_input_b);
	acb_init(arb_input_c);
	acb_init(arb_input_z);
	acb_init(arb_result);

	acb_set_d_d(arb_input_a, re_a, im_a);
	acb_set_d_d(arb_input_b, re_b, im_b);
	acb_set_d_d(arb_input_c, re_c, im_c);
	acb_set_d_d(arb_input_z, re_z, im_z);

	acb_hypgeom_2f1(arb_result, arb_input_a, arb_input_b, arb_input_c,
					arb_input_z, 0, prec);

	arb_t re_res, im_res;
	arb_init(re_res);
	arb_init(im_res);

	acb_get_real(re_res, arb_result);
	acb_get_imag(im_res, arb_result);

	res_real = arf_get_d(arb_midref(re_res), ARF_RND_NEAR);
	res_imag = arf_get_d(arb_midref(im_res), ARF_RND_NEAR);

	acb_clear(arb_input_a);
	acb_clear(arb_input_b);
	acb_clear(arb_input_c);
	acb_clear(arb_input_z);
	acb_clear(arb_result);
	arb_clear(re_res);
	arb_clear(im_res);

	complex<double> res {res_real, res_imag};
	return res;
}


// int main() {
// 	complex<double> a {3, 5};
// 	complex<double> b {-1./3., 1./2.};
// 	complex<double> c {1, 1};
// 	complex<double> z {0.25, 0};
// 	complex<double> res;
// 	res = hypgeo_2f1(a, b, c, z);
// 	double re_res, im_res;
// 	re_res = real(res);
// 	im_res = imag(res);
// 	printf("res = %.17g + %.17g i \n", re_res, im_res);
// }
























