try:
    from .spspec import fast_spheroidal_harmonics_eq
except:
    from spspec import fast_spheroidal_harmonics_eq

# TODO: test the case aa = 0
# TODO: include a check that stops people from using this function for
#       non-equatorial orbits

def calc_swsh_eq(aa, omega, ell, em, ess=-2):
    """
    Finds Slm(theta) based on a spectral decomposition

    Normalization is that from Glampedakis and Kennefick (2002)

    Parameters:
        aa (float): spin parameter (0, 1)
        omega (float): gravitational wave frequency
        ell (int): swsh mode
        em (int): azimuthal mode
        ess (int): 

    Returns:
        eigen (float): eigenvalue
        Slm (float): spin-weighted spheroidal harmonic (SWSH)
        dSlm (float): first derivative of SWSH
        d2Slm (float): second derivative of SWSH
    """
    eigen, Slm, dSlm, d2Slm = fast_spheroidal_harmonics_eq(aa, omega, ell, em, ess)
    return eigen, Slm, dSlm, d2Slm