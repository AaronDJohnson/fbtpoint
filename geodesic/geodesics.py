"""
This file takes the functions from each file and puts them in terms of
input variables aa, slr, ecc, x
"""

# !/usr/bin/env python
# -*- coding: utf-8 -*- #

# import relative files unless running in the folder
try:
    from .constants import fp_calc_constants
    from .frequencies import fp_radial_roots, fp_polar_roots, fp_boyer_freqs
    from .frequencies import fp_mino_freqs, fp_find_omega, fp_mino_freqs
    from .coords import py_calc_equatorial_coords
    from .py_coords import mp_calc_equatorial_coords
    from .py_coords import mp_calc_circular_eq_coords
except:
    from constants import fp_calc_constants
    from frequencies import fp_radial_roots, fp_polar_roots, fp_boyer_freqs
    from frequencies import fp_mino_freqs, fp_find_omega, fp_mino_freqs
    from coords import py_calc_equatorial_coords
    from py_coords import mp_calc_equatorial_coords
    from py_coords import mp_calc_circular_eq_coords


def calc_constants(aa, slr, ecc, x):
    """
    Compute adiabatic constants.

    Parameters:
        aa (float): SMBH spin
        slr (float): semi-latus rectum
        ecc (float): eccentricity
        x (float): cos of the inclination

    Returns:
        En (float): energy
        Lz (float): angular momentum
        Q (float): Carter constant
    """
    En, Lz, Q = fp_calc_constants(aa, slr, ecc, x)
    return En, Lz, Q


def radial_roots(aa, slr, ecc, x):
    """
    Compute radial roots.

    Parameters:
        aa (float): SMBH spin
        slr (float): semi-latus rectum
        ecc (float): eccentricity
        x (float): cos of the inclination

    Returns:
        r1 (float): apastron
        r2 (float): periastron
        r3 (float): radial root 3
        r4 (float): radial root 4
    """
    En, __, Q = fp_calc_constants(aa, slr, ecc, x)
    r1, r2, r3, r4 = fp_radial_roots(En, Q, aa, slr, ecc)
    return r1, r2, r3, r4


def polar_roots(aa, slr, ecc, x):
    """
    Compute polar roots.

    Parameters:
        aa (float): SMBH spin
        slr (float): semi-latus rectum
        ecc (float): eccentricity
        x (float): cos of the inclination

    Returns:
        zp (float): polar root
        zm (float): polar root
    """
    En, Lz, __ = fp_calc_constants(aa, slr, ecc, x)
    zp, zm = fp_polar_roots(En, Lz, aa, slr, x)
    return zp, zm


def mino_freqs(aa, slr, ecc, x, M=1):
    """
    Compute Mino frequencies.

    Parameters:
        aa (float): SMBH spin
        slr (float): semi-latus rectum
        ecc (float): eccentricity
        x (float): cos of the inclination

    Returns:
        ups_r (float): radial Mino frequency
        ups_theta (float): polar Mino frequency
        ups_phi (float): azimuthal Mino frequency
        gamma (float): temporal Mino frequency
    """
    En, Lz, Q = fp_calc_constants(aa, slr, ecc, x)
    r1, r2, r3, r4 = fp_radial_roots(En, Q, aa, slr, ecc, M)
    ups_r, ups_theta, ups_phi, gamma = fp_mino_freqs(r1, r2, r3, r4, En, Lz,
                                                        Q, aa, slr, ecc, x)
    return ups_r, ups_theta, ups_phi, gamma


def boyer_freqs(aa, slr, ecc, x, M=1):
    """
    Compute Boyer-Lindquist frequencies.

    Parameters:
        aa (float): SMBH spin
        slr (float): semi-latus rectum
        ecc (float): eccentricity
        x (float): cos of the inclination

    Returns:
        Omega_r (float): radial Boyer-Lindquist frequency
        Omega_theta (float): polar Boyer-Lindquist frequency
        Omega_phi (float): azimuthal Boyer-Lindquist frequency
    """
    En, Lz, Q = fp_calc_constants(aa, slr, ecc, x)
    r1, r2, r3, r4 = fp_radial_roots(En, Q, aa, slr, ecc, M)
    ups_r, ups_theta, ups_phi, gamma = fp_mino_freqs(r1, r2, r3, r4, En, Lz,
                                                     Q, aa, slr, ecc, x)
    omega_r, omega_theta, omega_phi = fp_boyer_freqs(ups_r, ups_theta, ups_phi,
                                                     gamma, aa, slr, ecc, x, M)
    return omega_r, omega_theta, omega_phi


def find_omega(en, em, kay, aa, slr, ecc, x, M=1):
    """
    Compute gravitational wave frequency omega.

    Parameters:
        en (int): radial mode
        em (int): azimuthal mode
        kay (int): polar mode
        aa (float): SMBH spin
        slr (float): semi-latus rectum
        ecc (float): eccentricity
        x (float): cos of the inclination

    Returns:
        omega (float): gravitational wave frequency
    """
    En, Lz, Q = fp_calc_constants(aa, slr, ecc, x)
    r1, r2, r3, r4 = fp_radial_roots(En, Q, aa, slr, ecc, M)
    ups_r, ups_theta, ups_phi, gamma = fp_mino_freqs(r1, r2, r3, r4, En, Lz,
                                                     Q, aa, slr, ecc, x)
    omega_r, omega_theta, omega_phi = fp_boyer_freqs(ups_r, ups_theta, ups_phi,
                                                     gamma, aa, slr, ecc, x, M)
    omega = en * omega_r + em * omega_phi + kay * omega_theta
    return omega


def find_equatorial_coords(psi, aa, slr, ecc, x, method='fast', M=1):
    """
    Finds the coordinates for a geodesic on the equator

    Parameters:
        psi (float): radial orbital parameter
        aa (float): SMBH spin
        slr (float): semi-latus rectum
        ecc (float): eccentricity
        x (float): cos of the inclination
    
    Returns:
        t (float): time geodesic coordinate
        r (float): radial geodesic coordinate
        theta (float): polar geodesic coordinate
        phi (float): azimuthal geodesic coordinate
    """

    En, Lz, Q = fp_calc_constants(aa, slr, ecc, x)
    r1, r2, r3, r4 = fp_radial_roots(En, Q, aa, slr, ecc, M)
    zp, zm = fp_polar_roots(En, Lz, aa, slr, x)
    ups_r, ups_theta, ups_phi, gamma = fp_mino_freqs(r1, r2, r3, r4, En, Lz,
                                                     Q, aa, slr, ecc, x)

    if method == 'fast':
        t, r, theta, phi = py_calc_equatorial_coords(psi, ups_r, ups_theta,
                                                     ups_phi, gamma, r1, r2, r3,
                                                     r4, zp, zm, En, Lz, aa,
                                                     slr, ecc)
    elif method == 'slow':
        if ecc == 0:
            t, r, theta, phi = mp_calc_circular_eq_coords(psi, En, Lz, aa, slr)
        else:
            t, r, theta, phi = mp_calc_equatorial_coords(psi, ups_r, ups_theta, ups_phi, gamma, r1, r2, r3,
                                                         r4, zp, zm, En, Lz, aa, slr, ecc)

    return t, r, theta, phi
