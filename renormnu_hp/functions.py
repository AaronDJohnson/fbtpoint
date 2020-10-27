from mpmath import mpf, mp
from .geo_const import mp_calc_constants
from .geo_freq import mp_boyer_freqs, mp_radial_roots, mp_mino_freqs
from .swsh_leaver import swsh_eigen, swsh_constants
from .renormnu import find_nu

def calc_nu(aa, slr, ecc, x, ell, en, em, kay, digits=100, ess=-2, M=1):
    """
    Find renormalized angular momentum using monodromy method. All functions here are
    computed to high precision for nu to be computed properly.

    Parameters:
        aa (float): SMBH spin
        slr (float): semi-latus rectum
        ecc (float): eccentricity
        x (float): cos of the inclination
        ell (int): SWSH mode
        en (int): radial mode
        em (int): azimuthal mode
        kay (int): polar mode
        digits (int): number of digits of accuracy requested

    Returns:
        nu (mpf): renormalized angular momentum
    """
    mp.dps = digits
    aa = mpf(str(aa))
    slr = mpf(str(slr))
    ecc = mpf(str(ecc))
    x = mpf(str(x))
    En, Lz, Q = mp_calc_constants(aa, slr, ecc, x)
    r1, r2, r3, r4 = mp_radial_roots(En, Q, aa, slr, ecc, M)
    ups_r, ups_theta, ups_phi, gamma = mp_mino_freqs(r1, r2, r3, r4, En, Lz,
                                                     Q, aa, slr, ecc, x)
    omega_r, omega_theta, omega_phi = mp_boyer_freqs(ups_r, ups_theta, ups_phi,
                                                     gamma, aa, slr, ecc, x, M)
    omega = en * omega_r + em * omega_phi + kay * omega_theta
    c, km, kp, nInv = swsh_constants(aa, omega, ell, em, ess)
    __, eigen = swsh_eigen(c, km, kp, ell, em, nInv, ess)
    return float(find_nu(aa, omega, eigen, ell, em))


    
    