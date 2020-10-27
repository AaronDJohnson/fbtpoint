"""
This file covers relevant cases for non-extremal spin SC/Kerr holes (a < 1).

This file is translated from ConstantsOfMotion.m from the Black Hole
Perturbation Toolkit by Niels Warburton et al.
"""

# !/usr/bin/env python
# -*- coding: utf-8 -*- #

from mpmath import sqrt, mp
from sys import exit

# ------------------------------------------------------------------------------
#  SC orbits (a = 0)
# ------------------------------------------------------------------------------


def calc_sc_constants(slr, ecc, x):
    """
    Energy, angular momentum, and carter constant calculation.

    Schwarzschild case (spin parameter a = 0)

    Parameters:
        slr (mpf): semi-latus rectum [6, inf)
        ecc (mpf): eccentricity [0, 1)
        x (mpf): inclination value given by cos(theta_inc) (0, 1]
                   negative x -> retrograde
                   positive x -> prograde

    Returns:
        En (mpf): energy
        Lz (mpf): angular momentum
        Q (mpf): Carter constant
    """
    try:
        prec = mp.prec
        mp.prec += 10

        ecc2 = ecc * ecc
        slr2 = slr * slr
        x2 = x * x

        mp.prec += 10
        En = sqrt(
            (-4 * ecc2 + (-2 + slr)**2) /
            (slr * (-3 - ecc2 + slr))
        )
        Lz = (slr * x) / sqrt(-3 - ecc2 + slr)
        Q = (slr2 * (-1 + x2)) / (3 + ecc2 - slr)

        return En, Lz, Q
    finally:
        mp.prec = prec


# ------------------------------------------------------------------------------
#  Kerr equatorial orbits (x = +/- 1) -> sign is pro/retrograde motion
# ------------------------------------------------------------------------------


def eq_energy(aa, slr, ecc, x):
    """
    Energy for the Kerr equatorial case (inclination value x = 1).

    Parameters:
        aa (mpf): spin parameter (0, 1)
        slr (mpf): semi-latus rectum [~2, inf)
        ecc (mpf): eccentricity [0, 1)
        x (mpf): inclination value given by cos(theta_inc)
                   x = -1 -> retrograde
                   x = +1 -> prograde

    Returns:
        En (mpf): energy
    """
    try:
        prec = mp.prec
        mp.prec += 10

        ecc2 = ecc * ecc
        eta = ecc2 - 1
        eta2 = eta * eta
        aa2 = aa * aa
        aa4 = aa2 * aa2
        aa6 = aa4 * aa2
        mp.prec += 10
        return (
            sqrt(
                1 - (-eta * (1 + (eta * (aa2 * (1 + 3 * ecc2 + slr) +
                     slr * (-3 - ecc2 + slr -
                     2 * sqrt((aa6 * eta2 +
                              aa2 * (-4 * ecc2 + (-2 + slr)**2) * slr**2 +
                              2 * aa4 * slr * (-2 + slr + ecc2 * (2 + slr))) /
                             (slr**3 * x**2)) * x))) /
                    (-4 * aa2 * eta2 + (3 + ecc2 - slr)**2 * slr))) / slr)
        )
    finally:
        # print('here')
        mp.prec = prec


def eq_ang_momentum(aa, slr, ecc, x):
    """
    Angular momentum for the Kerr equatorial case (inclination value x = 1).

    Parameters:
        aa (mpf): spin parameter (0, 1)
        slr (mpf): semi-latus rectum [6, inf)
        ecc (mpf): eccentricity [0, 1)
        x (mpf): inclination value given by cos(theta_inc)
                   x = -1 -> retrograde
                   x = +1 -> prograde

    Returns:
        Lz (mpf): angular momentum
    """
    try:
        prec = mp.prec
        mp.prec += 10

        ecc2 = ecc * ecc
        eta = ecc2 - 1
        eta2 = eta * eta
        aa2 = aa * aa
        aa4 = aa2 * aa2
        aa6 = aa4 * aa2
        slr2 = slr * slr
        slr3 = slr2 * slr
        x2 = x * x
        mp.prec += 10
        num_root = slr * (-3 - ecc2 + slr - 2 *
                          sqrt((aa6 * eta2 + aa2 *
                               (-4 * ecc2 + (-2 + slr)**2) *
                               slr2 + 2 * aa4 * slr *
                               (-2 + slr + ecc2 * (2 + slr))) /
                               (slr3 * x2)) * x)
        denom = (-4 * aa2 * eta2 + (3 + ecc2 - slr)**2 * slr)
        return (
            slr * x *
            sqrt((aa2 * (1 + 3 * ecc2 + slr) + num_root) / (denom * x2)) +
            aa * sqrt(1 - (-eta * (1 + (eta * (aa2 * (1 + 3 * ecc2 + slr) +
                      num_root)) / denom)) / slr)
        )
    finally:
        mp.prec = prec


def calc_eq_constants(aa, slr, ecc, x):
    """
    Calculate equatorial constants in one function.

    Note that for the equatorial case, the Carter constant, Q = 0.

    Parameters:
        aa (mpf): spin parameter (0, 1)
        slr (mpf): semi-latus rectum [6, inf)
        ecc (mpf): eccentricity [0, 1)
        x (mpf): inclination value given by cos(theta_inc)
                   x = -1 -> retrograde
                   x = +1 -> prograde

    Returns:
        En (mpf): energy
        Lz (mpf): angular momentum
        Q (mpf): Carter constant
    """
    En = eq_energy(aa, slr, ecc, x)
    Lz = eq_ang_momentum(aa, slr, ecc, x)
    Q = 0  # equatorial case
    return En, Lz, Q

# ------------------------------------------------------------------------------
#  Kerr spherical orbits (ecc = 0)
# ------------------------------------------------------------------------------


def spherical_energy(aa, slr, x):
    """
    Compute energy for the Kerr spherical case (ecc = 0).

    Parameters:
        aa (mpf): spin parameter (0, 1)
        slr (mpf): semi-latus rectum [6, inf)
        x (mpf): inclination value given by cos(theta_inc) (0, 1]
                   negative x -> retrograde
                   positive x -> prograde

    Returns:
        En (mpf): energy
    """
    try:
        prec = mp.prec
        mp.prec += 10

        slrm2 = (-2 + slr)**2
        slr2 = slr * slr
        slr3 = slr * slr2
        slr4 = slr * slr3
        slr5 = slr3 * slr2
        slr7 = slr5 * slr2
        x2 = x * x
        x3 = x2 * x
        aa2 = aa * aa
        aa3 = aa2 * aa
        aa4 = aa2 * aa2
        aa5 = aa4 * aa
        aa6 = aa4 * aa2
        x2m1 = x2 - 1
        x2m12 = x2m1 * x2m1

        mp.prec += 10
        return (
            sqrt(((-3 + slr) * slrm2 * slr5 -
                 2 * aa5 * x * x2m1 * sqrt(slr3 + aa2 * slr * x2m1) +
                 aa4 * slr2 * x2m1 * (4 - 5 * slr * (-1 + x2) +
                 3 * slr2 * x2m1) - aa6 * x2m12 *
                 (x2 + slr2 * x2m1 - slr * (1 + 2 * x2)) +
                 aa2 * slr3 * (4 - 4 * x2 + slr * (12 - 7 * x2) -
                               3 * slr3 * (-1 + x2) +
                               slr2 * (-13 + 10 * x2)) +
                 aa * (-2 * slr**4.5 * x * sqrt(slr2 + aa2 * x2m1) +
                       4 * slr3 * x * sqrt(slr3 + aa2 * slr * x2m1)) +
                 2 * aa3 * (2 * slr * x * x2m1 *
                            sqrt(slr3 + aa2 * slr * x2m1) -
                            x3 * sqrt(slr7 + aa2 * slr5 * x2m1))) /
                 ((slr2 - aa2 * x2m1) *
                 ((-3 + slr)**2 * slr4 - 2 * aa2 * slr2 *
                  (3 + 2 * slr - 3 * x2 + slr2 * x2m1) +
                  aa4 * x2m1 * (-1 + x2 + slr2 * x2m1 -
                  2 * slr * (1 + x2)))))
        )
    finally:
        mp.prec = prec


def spherical_ang_momentum(En, aa, slr, x):
    """
    Compute angular momentum for the Kerr spherical case (ecc = 0).

    Parameters:
        En (mpf): energy
        aa (mpf): spin parameter (0, 1)
        slr (mpf): semi-latus rectum [6, inf)
        x (mpf): inclination value given by cos(theta_inc) (0, 1]
                   negative x -> retrograde
                   positive x -> prograde

    Returns:
        Lz (mpf): angular momentum
    """
    try:
        prec = mp.prec
        mp.prec += 10

        aa2 = aa * aa
        slr2 = slr * slr
        x2 = x * x
        mp.prec += 10
        g = 2 * aa * slr
        d = (aa2 + (-2 + slr) * slr) * (slr2 - aa2 * (-1 + x2))
        h = ((-2 + slr) * slr - aa2 * (-1 + x2)) / x2
        f = (slr**4 + aa2 * (slr * (2 + slr) - (aa2 + (-2 + slr) * slr) *
             (-1 + x2)))

        return (-(En * g) +
                sqrt((-(d * h) + En**2 * (g**2 + f * h)) / x2) * x) / h
    finally:
        mp.prec = prec

# ------------------------------------------------------------------------------
#  Generic orbit constant calculation
# ------------------------------------------------------------------------------


def calc_delta(r, aa):
    """
    Calculate ubiquitous function on Kerr spacetimes.

    Parameters:
        r (mpf): radius
        aa (mpf): spin parameter (0, 1)

    Returns:
        delta (mpf)
    """
    try:
        prec = mp.prec
        mp.prec += 20
        # return r * r - 2 * r + aa * aa
        return r * (r - 2) + aa * aa
    finally:
        mp.prec = prec


def calc_f(r, zm, aa):
    """
    DESCRIPTION TBD (FUNCTION USED IN GENERAL CONSTANTS).

    Parameters:
        r (mpf): radius
        zm (mpf): 1 - x * x
        aa (mpf): spin parameter (0, 1)

    Returns:
        f (mpf)
    """
    try:
        prec = mp.prec
        mp.prec += 10
        r2 = r * r
        mp.prec += 10
        r4 = r2 * r2
        aa2 = aa * aa
        zm2 = zm * zm

        delta = calc_delta(r, aa)
        mp.prec += 10
        return 2 * aa2 * r + aa2 * r2 + r4 + aa2 * zm2 * delta
    finally:
        mp.prec = prec


def calc_g(r, aa):
    """
    DESCRIPTION TBD (FUNCTION USED IN GENERAL CONSTANTS).

    Parameters:
        r (mpf): radius
        aa (mpf): spin parameter (0, 1)

    Returns:
        g (mpf)
    """
    try:
        prec = mp.prec
        mp.prec += 10
        return 2 * aa * r
    finally:
        mp.prec = prec


def calc_h(r, zm, aa):
    """
    DESCRIPTION TBD (FUNCTION USED IN GENERAL CONSTANTS).

    Parameters:
        r (mpf): radius
        zm (mpf): 1 - x * x
        aa (mpf): spin parameter (0, 1)

    Returns:
        h (mpf)
    """
    try:
        prec = mp.prec
        mp.prec += 10
        zm2 = zm * zm
        delta = calc_delta(r, aa)
        mp.prec += 10
        return (-2 + r) * r + (zm2 * delta) / (1 - zm2)
    finally:
        mp.prec = prec


def calc_d(r, zm, aa):
    """
    DESCRIPTION TBD (FUNCTION USED IN GENERAL CONSTANTS).

    Parameters:
        r (mpf): radius
        zm (mpf): 1 - x * x
        aa (mpf): spin parameter (0, 1)

    Returns:
        d (mpf)
    """
    try:
        prec = mp.prec
        mp.prec += 10
        r2 = r * r
        aa2 = aa * aa
        zm2 = zm * zm

        delta = calc_delta(r, aa)
        mp.prec += 10
        return (r2 + aa2 * zm2) * delta
    finally:
        mp.prec = prec


def gen_energy(zm, aa, slr, ecc, x):
    """
    Compute energy for generic orbit case.

    Parameters:
        zm (mpf): 1 - x * x
        aa (mpf): spin parameter (0, 1)
        slr (mpf): semi-latus rectum [~2, inf)
        ecc (mpf): eccentricity [0, 1)
        x (mpf): inclination value given by cos(theta_inc)
                   x < 0 -> retrograde
                   x > 0 -> prograde

    Returns:
        En (mpf): energy
    """
    try:
        prec = mp.prec
        mp.prec += 20
        r1 = slr / (1 - ecc)
        r2 = slr / (1 + ecc)

        dr1 = calc_d(r1, zm, aa)
        dr2 = calc_d(r2, zm, aa)
        gr1 = calc_g(r1, aa)
        gr2 = calc_g(r2, aa)
        hr1 = calc_h(r1, zm, aa)
        hr2 = calc_h(r2, zm, aa)
        fr1 = calc_f(r1, zm, aa)
        fr2 = calc_f(r2, zm, aa)

        mp.prec += 10
        kappa = dr1 * hr2 - hr1 * dr2
        epsilon = dr1 * gr2 - gr1 * dr2
        rho = fr1 * hr2 - hr1 * fr2
        eta = fr1 * gr2 - gr1 * fr2
        sigma = gr1 * hr2 - hr1 * gr2

        mp.prec += 10
        kappa2 = kappa * kappa
        epsilon2 = epsilon * epsilon
        rho2 = rho * rho
        x2 = x * x

        mp.prec += 10
        En = sqrt(
            (kappa * rho + 2 * epsilon * sigma -
             2 * sqrt((sigma * (-(eta * kappa2) + epsilon * kappa * rho +
                      epsilon2 * sigma)) / x2) * x) / (rho2 + 4 * eta * sigma))

        return En
    finally:
        mp.prec = prec


def gen_ang_momentum(En, aa, slr, ecc, x):
    """
    Compute angular momentum for generic orbit case.

    Parameters:
        En (mpf): energy
        aa (mpf): spin parameter (0, 1)
        slr (mpf): semi-latus rectum [6, inf)
        ecc (mpf): eccentricity [0, 1)
        x (mpf): inclination value given by cos(theta_inc)
                   x < 0 -> retrograde
                   x > 0 -> prograde

    Returns:
        z_ang_mom (mpf): angular momentum
    """
    try:
        prec = mp.prec
        mp.prec += 10
        r1 = slr / (1 - ecc)
        zm = sqrt(1 - x * x)

        # // delta_r1 = calc_delta(r1, aa)
        mp.prec += 10
        fr1 = calc_f(r1, zm, aa)
        gr1 = calc_g(r1, aa)
        hr1 = calc_h(r1, zm, aa)
        dr1 = calc_d(r1, zm, aa)

        mp.prec += 10
        En2 = En * En
        x2 = x * x
        mp.prec += 10
        Lz = ((-(En * gr1) + x * sqrt((-(dr1 * hr1) + En2 *
                     (gr1 * gr1 + fr1 * hr1)) / x2)) / hr1)
        return Lz
    finally:
        mp.prec = prec


def gen_carter_const(En, Lz, aa, slr, ecc, x):
    """
    Compute Carter constant for generic orbit case.

    Parameters:
        En (mpf): energy
        Lz (mpf): angular momentum
        aa (mpf): spin parameter (0, 1)
        slr (mpf): semi-latus rectum [6, inf)
        ecc (mpf): eccentricity [0, 1)
        x (mpf): inclination value given by cos(theta_inc)
                   x < 0 -> retrograde
                   x > 0 -> prograde

    Returns:
        Q (mpf): Carter constant
    """
    try:
        prec = mp.prec
        mp.prec += 10
        zm = sqrt(1 - x * x)
        mp.prec += 10
        zm2 = zm * zm
        mp.prec += 10
        return (zm2 * (aa * aa * (1 - En * En) + Lz * Lz /
                (1 - zm2)))
    finally:
        mp.prec = prec


def calc_gen_constants(aa, slr, ecc, x):
    """
    Call generic orbit constant calculating functions in one function.

    Parameters:
        aa (mpf): spin parameter (0, 1)
        slr (mpf): semi-latus rectum [6, inf)
        ecc (mpf): eccentricity [0, 1)
        x (mpf): inclination value given by cos(theta_inc)
                   x < 0 -> retrograde
                   x > 0 -> prograde

    Returns:
        En (mpf): energy
        Lz (mpf): angular momentum
        Q (mpf): Carter constant
    """
    try:
        prec = mp.prec
        mp.prec += 30
        zm = sqrt(1 - x * x)

        En = gen_energy(zm, aa, slr, ecc, x)
        Lz = gen_ang_momentum(En, aa, slr, ecc, x)
        Q = gen_carter_const(En, Lz, aa, slr, ecc, x)
        return En, Lz, Q
    finally:
        mp.prec = prec


def pol_energy(aa, slr, ecc):
    """
    Calculate energy for polar Kerr case.

    Parameters:
        aa (mpf): spin parameter (0, 1)
        slr (mpf): semi-latus rectum [6, inf)
        ecc (mpf): eccentricity

    Returns:
        En (mpf): energy
    """
    try:
        prec = mp.prec
        mp.prec += 10
        aa2 = aa * aa
        aa4 = aa2 * aa2
        ecc2 = ecc * ecc
        ecc4 = ecc2 * ecc2
        slr2 = slr * slr
        slr4 = slr2 * slr2
        mp.prec += 10
        return sqrt(-((slr * (aa4 * (-1 + ecc2)**2 +
                    (-4 * ecc2 + (-2 + slr)**2) * slr2 + 
                    2 * aa2 * slr * (-2 + slr + ecc2 * (2 + slr)))) /
                    (aa4 * (-1 + ecc2)**2 * (-1 + ecc2 - slr) +
                     (3 + ecc2 - slr) * slr4 - 2 * aa2 * slr2 *
                     (-1 - ecc4 + slr + ecc2*(2 + slr)))))
    finally:
        mp.prec = prec


def pol_carter_const(aa, slr, ecc):
    """
    Calculate energy for polar Kerr case.

    Parameters:
        aa (mpf): spin parameter (0, 1)
        slr (mpf): semi-latus rectum [6, inf)
        ecc (mpf): eccentricity

    Returns:
        Q (mpf): Carter constant
    """
    try:
        prec = mp.prec
        mp.prec += 10
        aa2 = aa * aa
        aa4 = aa2 * aa2
        ecc2 = ecc * ecc
        ecc4 = ecc2 * ecc2
        slr2 = slr * slr
        slr4 = slr2 * slr2
        mp.prec += 10
        return -((slr2 * (aa4 * (-1 + ecc2)**2 + slr4 + 
                 2 * aa2 * slr * (-2 + slr + ecc2 * (2 + slr)))) /
                 (aa4 * (-1 + ecc2)**2 *
                  (-1 + ecc2 - slr) + (3 + ecc2 - slr) * slr4 - 
                 2 * aa2 * slr2 * (-1 - ecc4 + slr + ecc2 * (2 + slr))))
    finally:
        mp.prec = prec


def calc_sph_constants(aa, slr, x):
    """
    Call spherical orbit (ecc = 0) functions.

    Parameters:
        aa (mpf): spin parameter (0, 1)
        slr (mpf): semi-latus rectum [6, inf)
        x (mpf): inclination value given by cos(theta_inc)
                   x < 0 -> retrograde
                   x > 0 -> prograde

    Returns:
        En (mpf): energy
        Lz (mpf): angular momentum
        Q (mpf): Carter constant
    """
    ecc = 0

    En = spherical_energy(aa, slr, x)
    Lz = spherical_ang_momentum(En, aa, slr, x)
    Q = gen_carter_const(En, Lz, aa, slr, ecc, x)

    return En, Lz, Q


def calc_pol_constants(aa, slr, ecc):
    """
    Call Kerr polar orbit (x = 0) functions.

    Parameters:
        aa (mpf): spin parameter [0, 1)
        slr (mpf): semi-latus rectum
        ecc (mpf): eccentricity [0, 1)

    Returns:
        En (mpf): energy
        Lz (mpf): angular momentum
        Q (mpf): Carter constant
    """
    Lz = 0

    En = pol_energy(aa, slr, ecc)
    Q = pol_carter_const(aa, slr, ecc)

    return En, Lz, Q


def mp_calc_constants(aa, slr, ecc, x):
    """
    Choose which function to call based on input parameters.

    This version uses mpmath for extended precision.

    Parameters:
        aa (mpf): spin parameter (0, 1)
        slr (mpf): semi-latus rectum [6, inf)
        ecc (mpf): eccentricity [0, 1)
        x (mpf): inclination value given by cos(theta_inc)
                   x < 0 -> retrograde
                   x > 0 -> prograde

    Returns:
        En (mpf): energy
        Lz (mpf): angular momentum
        Q (mpf): Carter constant
    """
    if aa == 0:
        if slr < 6:
            print("slr value is too small for the SC case.")
            exit("Error in input values.")
        return calc_sc_constants(slr, ecc, x)
    elif x == 0:
        return calc_pol_constants(aa, slr, ecc)
    elif x**2 == 1:
        return calc_eq_constants(aa, slr, ecc, x)
    elif ecc == 0:
        return calc_sph_constants(aa, slr, x)
    else:
        return calc_gen_constants(aa, slr, ecc, x)
