import numpy as np
from .spspec import fast_spheroidal_harmonics

# TODO: add test cases

def calc_swsh_eq(phi, aa, omega, ell, em, ess=-2):
    """
    Finds Slm(pi/2) based on a spectral decomposition

    Normalization is that from Glampedakis and Kennefick (2002)

    Inputs:
        aa (float): spin parameter (0, 1)
        omega (float): gravitational wave frequency
        ell (int): swsh mode
        em (int): azimuthal mode
        ess (int): spin number

    Returns:
        eigen (float): eigenvalue
        Slm (float): spin-weighted spheroidal harmonic (SWSH)
        dSlm (float): first derivative of SWSH
        d2Slm (float): second derivative of SWSH
    """
    theta = np.pi / 2
    eigen, Slm, dSlm, d2Slm = fast_spheroidal_harmonics(theta, phi, aa, omega, ell, em, ess)
    return eigen, Slm, dSlm, d2Slm


def calc_swsh(theta, phi, aa, omega, ell, em, ess=-2):
    """
    Finds Slm(theta) based on a spectral decomposition

    Normalization is that from Glampedakis and Kennefick (2002)

    Inputs:
        aa (float): spin parameter (0, 1)
        omega (float): gravitational wave frequency
        ell (int): swsh mode
        em (int): azimuthal mode
        ess (int): spin number

    Returns:
        eigen (float): eigenvalue
        Slm (complex float): spin-weighted spheroidal harmonic (SWSH)
        dSlm (complex float): first derivative of SWSH
        d2Slm (complex float): second derivative of SWSH
    """
    eigen, Slm, dSlm, d2Slm = fast_spheroidal_harmonics(theta, phi, aa, omega, ell, em, ess)
    return eigen, Slm, dSlm, d2Slm
