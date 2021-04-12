import numpy as np


def trapezoidal_rule(f, a, b, tol=1e-8):
    """
    The trapezoidal rule is known to be very accurate for
    oscillatory integrals integrated over their period.

    See papers on spectral integration (it's just the composite trapezoidal rule....)

    TODO (aaron): f is memoized to get the already computed points quickly.
    Ideally, we should put this into a C++ function and call it with Cython. (Maybe someday)
    """
    # endpoints first:
    num = 2
    dx = b - a
    res0 = 1e30
    res1 = 0.5 * dx * (f(b) + f(a))
    delta_res = res0 - res1
    re_err = np.abs(np.real(delta_res))
    im_err = np.abs(np.imag(delta_res))
    while re_err > tol or im_err > tol:
        res0 = res1
        num = 2 * num - 1
        # print(num)
        x = np.linspace(a, b, num=num)
        res = 0
        dx = (x[1] - x[0])
        res += f(x[0])
        for i in range(1, len(x) - 1):
            res += 2 * f(x[i])
        res += f(x[-1]) 
        res1 = 0.5 * dx * res
        delta_res = res1 - res0
        re_err = np.abs(np.real(delta_res))
        im_err = np.abs(np.imag(delta_res))
        if num > 100000:
            print('Integral failed to converge with', num, 'points.')
            return np.nan, np.nan, np.nan
    return res1, re_err, im_err