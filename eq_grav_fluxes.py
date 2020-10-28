# import sys
# sys.path.append('/Users/aaronjohnson/Desktop/pykerr/')

from mpmath import mpf, mp

import numpy as np

from geodesic_fp.constants import fp_calc_constants
from geodesic_fp.frequencies import fp_radial_roots, fp_polar_roots, fp_boyer_freqs
from geodesic_fp.frequencies import fp_mino_freqs, fp_find_omega

from swsh_fp.swsh import calc_swsh_eq
# from swsh.accurate.leaver import swsh_eigen, swsh_constants

from renormnu_hp.functions import calc_nu
from teukolsky_fp.homteuk import py_find_Bin

from flux_fp.flux import flux_inf
from data import make_folder, save_data, save_wave_data


def equatorial_fluxes(slr, ecc, aa, x, ell, em, kay, ess=-2, verbose=True):
    """
    """
    # TODO (aaron): clean this up and hide some of these functions elsewhere
    # TODO (aaron): ecc = 0 sometimes doesn't work (for the FT data in Throwe's
    #               thesis for example)
    #               The easiest way is to change to the circ calc from BHPTK
    foldername = make_folder(slr, ecc, aa)
    if verbose:
        print("----------------")
        print("mode parameters:")
        print("l =", ell)
        print("m =", em)
        print("k =", kay)
        print("----------------")

    En, Lz, Q = fp_calc_constants(aa, slr, ecc, 1)
    r1, r2, r3, r4 = fp_radial_roots(En, Q, aa, slr, ecc)
    zp, zm = fp_polar_roots(En, Lz, aa, slr, x)
    ups_r, ups_theta, ups_phi, gamma = fp_mino_freqs(r1, r2, r3, r4, En, Lz, Q,
                                                     aa, slr, ecc, x)
    omega_r, omega_theta, omega_phi = fp_boyer_freqs(ups_r, ups_theta, ups_phi,
                                                     gamma, aa, slr, ecc, x)

    omega = fp_find_omega(omega_r, omega_theta, omega_phi, em, 0, kay)
    nu = calc_nu(aa, slr, ecc, x, ell, kay, em, 0)

    eigen, Slm, Slmd, Slmdd = calc_swsh_eq(aa, omega, ell, em, ess)

    re_nu = np.real(nu)
    im_nu = np.imag(nu)
    Bin = py_find_Bin(re_nu, im_nu, eigen, aa, omega, em)

    energy_inf, Z = flux_inf(nu, Bin, eigen, slr, ecc, aa, ups_r, ups_theta,
                             ups_phi, gamma, omega, em, Lz, En, Slm, Slmd,
                             Slmdd, omega_r, r1, r2, r3, r4, zp, zm)
    if np.isnan(energy_inf):  # TODO (aaron): investigate why this can be nan
        energy_inf = 0
    print("E_inf = ", 2 * energy_inf)
    save_data(foldername, ell, em, kay, 2 * energy_inf)
    save_wave_data(foldername, ell, em, kay, omega, Z, Slm)


def mode_search(slr, ecc, aa, ellmax, kaymax, ess=-2, verbose=True):
    for ell in range(2, ellmax + 1):
        for em in range(1, ell + 1):
            for kay in range(-kaymax, kaymax + 1):
                equatorial_fluxes(slr, ecc, aa, ell, em, kay, ess, verbose)


if __name__ == "__main__":
    slr = 7
    ecc = 0.1
    aa = 0.998
    ell = 2
    em = 2
    kay = 0
    x = 1

    # mode_search(slr, ecc, aa, 2, 5)
    equatorial_fluxes(slr, ecc, aa, x, ell, em, kay, verbose=True)








