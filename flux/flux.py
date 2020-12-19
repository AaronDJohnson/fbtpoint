import os, sys, inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)

from geodesic.coords import py_calc_equatorial_coords
from source.tsource import eq_source
from teukolsky.homteuk import py_find_R

from scipy.integrate import romberg, quad

import numpy as np

from mpmath import mp

# TODO (aaron): set up special case of circular orbits (no integral in this case)


def eq_find_z(nu, Bin, eigen, slr, ecc, aa, ups_r, ups_theta, ups_phi, gamma,
              omega, em, Lz, En, Slm, Slmd, Slmdd, omega_r, r1, r2, r3, r4, zp,
              zm, ess=-2):

    # print(type(nu))
    # print(type(Bin))
    # print(type(eigen))
    # print(type(slr))
    # print(type(ecc))
    # print(type(aa))
    # print(type(omega))
    # print(type(Lz))
    # print(type(En))
    # print(type(omega_r))
    # print(type(r1))
    # print(type(zp))

    def find_psi_integrand(psi):
        t, r, __, phi = py_calc_equatorial_coords(psi, ups_r, ups_theta,
                                                  ups_phi, gamma, r1,
                                                  r2, r3, r4, zp, zm, En, Lz,
                                                  aa, slr, ecc)
        re_nu = np.real(nu)
        im_nu = np.imag(nu)
        Rin, dRdr, dRdr2 = py_find_R(r, re_nu, im_nu, aa, omega, em, eigen)
        J, V_t, V_r, I_plus = eq_source(psi, 1, slr, ecc, aa, omega, em, Lz,
                                        En, Slm, Slmd, Slmdd, Rin, dRdr, dRdr2)
        _, _, _, I_minus = eq_source(psi, -1, slr, ecc, aa, omega, em, Lz,
                                     En, Slm, Slmd, Slmdd, Rin, dRdr, dRdr2)

        result = (
            V_t / (J * np.sqrt(V_r)) *
            (I_plus * np.exp(1j * omega * t - 1j * em * phi) +
             I_minus * np.exp(-1j * omega * t + 1j * em * phi))
        )
        return result

    # coeff = abs(find_psi_integrand(0)) / abs(find_psi_integrand(np.pi))
    # print(coeff)

    def integrand(zeta):
        psi = coeff * np.log(1 + zeta / coeff)
        dchi_dzeta = np.exp(-psi / coeff)

        t, r, __, phi = py_calc_equatorial_coords(psi, ups_r, ups_theta,
                                                  ups_phi, gamma, r1,
                                                  r2, r3, r4, zp, zm, En, Lz,
                                                  aa, slr, ecc)

        re_nu = np.real(nu)
        im_nu = np.imag(nu)
        Rin, dRdr, dRdr2 = py_find_R(r, re_nu, im_nu, aa, omega, em, eigen)
        J, V_t, V_r, I_plus = eq_source(psi, 1, slr, ecc, aa, omega, em, Lz,
                                        En, Slm, Slmd, Slmdd, Rin, dRdr, dRdr2)
        _, _, _, I_minus = eq_source(psi, -1, slr, ecc, aa, omega, em, Lz,
                                     En, Slm, Slmd, Slmdd, Rin, dRdr, dRdr2)

        result = (
            V_t / (J * np.sqrt(V_r)) *
            (I_plus * np.exp(1j * omega * t - 1j * em * phi) +
             I_minus * np.exp(-1j * omega * t + 1j * em * phi))
        ) * dchi_dzeta
        # print(result)
        return result

    def integrand_re(zeta):
        return np.real(integrand(zeta))

    def integrand_im(zeta):
        return np.imag(integrand(zeta))

    def cheap_re(zeta):
        return np.real(find_psi_integrand(zeta))

    def cheap_im(zeta):
        return np.imag(find_psi_integrand(zeta))

    # a = 0
    # mp.dps += 50
    # TODO: fix overflow in the following line
    # b = float((np.exp(np.pi / coeff) - 1) * coeff)
    # print(b)
    # print("b = ", b)
    # re_res, re_err = quad(integrand_re, a, b)
    # im_res, im_err = quad(integrand_im, a, b)

    re_res, re_err = quad(cheap_re, 0, np.pi)
    im_res, im_err = quad(cheap_im, 0, np.pi)

    print(re_err)
    print(im_err)

    # re_res = romberg(integrand_re, a, b, divmax=20)
    # im_res = romberg(integrand_im, a, b, divmax=20)

    Z = omega_r / (2 * 1j * omega * Bin) * (re_res + 1j * im_res)
    return Z


def flux_inf(nu, Bin, eigen, slr, ecc, aa, ups_r, ups_theta, ups_phi,
             gamma, omega, em, Lz, En, Slm, Slmd, Slmdd, omega_r, r1, r2,
             r3, r4, zp, zm, ess=-2):
    Z = eq_find_z(nu, Bin, eigen, slr, ecc, aa, ups_r, ups_theta, ups_phi,
                  gamma, omega, em, Lz, En, Slm, Slmd, Slmdd, omega_r, r1, r2,
                  r3, r4, zp, zm, ess)
    energy = abs(Z)**2 / (4 * np.pi * omega**2)
    return energy, Z
