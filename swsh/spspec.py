import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigs
from scipy.special import comb, factorial

#-------------------------------------------------------------------------------
# create and solve the spectral method matrix to get eigenvalues and coeffs
# for Ylms
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# TODO List:
#
# * make sure that the spectral matrix is symmetric for every situation
# * solve only for the necessary eigenvalue (closest to A0) and its vector
# * port everything to C++ to make this fast and arbitrary precision
#    - this should use Eigen and GMP
# * check that the solution is correct for different cases
# * figure out what to do in the SC case (this should be easy)
#    - the relevant case is a = 0 => a * omega = 0
#    - this will reduce down to sYlm and A0
# * think about whether nmax is large enough (incorrect if not!)
#-------------------------------------------------------------------------------


def swsh_constants(aa, omega, ell, em, ess=-2):
    km = abs(em - ess)
    kp = abs(em + ess)
    c = aa * omega
    if em == 1:
        ell = ell - 1
    return c, km, kp


def kHat(c, ell, em, ess=-2):
    # print(ell)
    # print(np.where(ell == 0)[0])
    # if ell == 0 and em == 0:
    #     return pow(c, 2 / 3)
    # else:
    ell2 = ell * ell
    em2 = em * em
    ell3 = ell2 * ell
    ess2 = ess * ess
    c2 = c * c

    return (-(ell*(1 + ell)) + (2*em*ess2*c)/(ell + ell2) + 
             ((1 + (2*(ell + ell2 - 3*em2)*(ell + ell2 - 3*ess2))/
             (ell*(-3 + ell + 8*ell2 + 4*ell3)))*c2)/3.)


def k2(c, ell, em, ess=-2):
    c2 = c * c
    ellmem = ell - em
    ellpem = ell + em
    ellmess = ell - ess
    ellpess = ell + ess

    ellmemp1 = ellmem + 1
    ellmemp2 = ellmem + 2
    ellpemp1 = ellpem + 1
    ellpemp2 = ellpem + 2
    ellmessp1 = ellmess + 1
    ellmessp2 = ellmess + 2
    ellpessp1 = ellpess + 1
    ellpessp2 = ellpess + 2
    ell2p1 = 2 * ell + 1
    ell2p5 = 2 * ell + 5

    fraction = ((ellmemp1*ellmemp2*ellpemp1*ellpemp2*
            ellmessp1*ellmessp2*ellpessp1*ellpessp2)/
            (ell2p1*ell2p5))

    return ((c2*np.sqrt(fraction))/((1 + ell)*(2 + ell)*(3 + 2*ell)))


def kTilde2(c, ell, em, ess=-2):
    ess2 = ess * ess
    # if (ell == 0 and em == 0):
    #     return (-2*c*ess*np.sqrt(1 - ess2))/np.sqrt(3)
    # else:
    ell2 = ell * ell
    em2 = em * em
    return ((-2*c*(2*ell + ell2 + c*em)*ess*
            np.sqrt(((1 + 2*ell + ell2 - em2)*(1 + 2*ell + ell2 - ess2))/
            (3 + 8*ell + 4*ell2)))/(ell*(2 + 3*ell + ell2)))


def sparse_spectral_matrix(c, ell, em, ess=-2):
    """
    Is this symmetric?
    """
    lmin = max(abs(ess), abs(em))

    nmax = 50 + np.ceil(abs((3*c)/2. - c*c/250.))
    if nmax % 2 == 0:
        nmax += 1
    nmax = int(nmax)

    nmin = min(ell - lmin, nmax)
    i = np.arange(0, nmax + nmin + 1)  # iterator
    spec_ii = -kHat(c, ell - nmin - 1 + i + 1, em, ess)  # good

    i = np.delete(i, 0)
    spec_ijm1 = -kTilde2(c, ell - nmin + i - 2 + 1, em, ess)  # good
    spec_ijp1 = -kTilde2(c, ell - nmin + i - 2 + 1, em, ess)  # good

    i = np.delete(i, 0)
    spec_ijp2 = -k2(c, ell - nmin + i - 3 + 1, em, ess)  # good
    spec_ijm2 = -k2(c, ell - nmin - 3 + i + 1, em, ess)  # good

    diagonals = [spec_ijm2, spec_ijm1, spec_ii, spec_ijp1, spec_ijp2]
    # print(diagonals)
    band_matrix = sparse.diags(diagonals, [-2, -1, 0, 1, 2])
    return band_matrix


def solve_sparse_matrix(c, ell, em, ess=-2):
    matrix = sparse_spectral_matrix(c, ell, em, ess)
    # print(matrix.todense())
    A0 = ell*(ell + 1) - ess*(ess + 1)
    eigenvalues, eigenvectors = eigs(matrix, 1, sigma=A0)
    # print(eigenvalues)
    # print(eigenvectors.T)
    xmin = min(np.real(eigenvalues)) - ess * (ess + 1)
    eigen = xmin + c * c - 2 * c * em
    # ind_max = np.argmax(np.abs(np.real(eigenvectors)))
    # print(ind_max)
    if eigenvectors[0] < 0:  # largest value is always [0]
        # print('mult by -1')
        eigenvectors = eigenvectors * -1
    coeffs = eigenvectors.T[0]
    # print(coeffs)
    # coeffs = np.where(abs(eigenvectors.T[0])>1e-16, eigenvectors.T[0], 0)
    # print(eigenvalues)
    # print(eigenvectors.T)
    # print(coeffs)

    return eigen, coeffs

#-------------------------------------------------------------------------------
# We'll need Ylm and it's first derivative for the Slms
#-------------------------------------------------------------------------------


def sphericalY(theta, phi, ell, em, ess=-2):
    prefactor = (((-1)**em*np.sqrt(((1 + 2*ell)*factorial(ell - em)*
        factorial(ell + em))/ (factorial(ell - ess)*factorial(ell + ess))))/
        (2.*np.sqrt(np.pi)) * np.exp(1j * em * phi))
    start = max(em - ess, 0)
    end = min(ell - ess, ell + em)
    j = np.arange(start, end + 1)
    binoms = comb(ell - ess, j) * comb(ell + ess, j + ess - em)
    phase = (-1)**(ell - j - ess)
    sin = np.sin(theta / 2)**(2 * ell - 2 * j - ess + em)
    cos = np.cos(theta / 2)**(2 * j + ess - em)
    total = np.sum(binoms * phase * sin * cos)
    return prefactor * total


def sphericalY_dtheta(theta, phi, ell, em, ess=-2):
    prefactor = (((-1)**em*np.sqrt(((1 + 2*ell)*factorial(ell - em)*
        factorial(ell + em))/ (factorial(ell - ess)*factorial(ell + ess))))/
        (2.*np.sqrt(np.pi)) * np.exp(1j * em * phi))
    start = max(em - ess, 0)
    end = min(ell - ess, ell + em)
    j = np.arange(start, end + 1)
    binoms = comb(ell - ess, j) * comb(ell + ess, j + ess - em)
    phase = (-1)**(ell - j - ess)
    deriv = (np.cos(theta/2.)**(-em + ess + 2*j)*(ell + em - ess - 2*j +
            ell*np.cos(theta))*
            (1/np.sin(theta))*np.sin(theta/2.)**(2*ell + em - ess - 2*j))
    total = np.sum(binoms * phase * deriv)
    return prefactor * total


# def sphericalY_d2theta(theta, phi, ell, em, ess=-2):
#     """
#     This function does not currently work. Check deriv.
#     """
#     prefactor = (((-1)**em*np.sqrt(((1 + 2*ell)*factorial(ell - em)*
#         factorial(ell + em))/ (factorial(ell - ess)*factorial(ell + ess))))/
#         (2.*np.sqrt(np.pi)) * np.exp(1j * em * phi))
#     start = max(em - ess, 0)
#     end = min(ell - ess, ell + em)
#     j = np.arange(start, end + 1)
#     binoms = comb(ell - ess, j) * comb(ell + ess, j + ess - em)
#     phase = (-1)**(ell - j - ess)
#     deriv = ((np.cos(theta/2.)**(-em + ess + 2*j)*
#             (-2*ell + 2*(ell + em - ess - 2*j + (-1 + ell)*np.cos(theta))*
#             (ell + em - ess - 2*j + ell*np.cos(theta))*(1/np.sin(theta))**2)*
#             np.sin(theta/2.)**(2*ell + em - ess - 2*j))/2.)
#     total = np.sum(binoms * phase * deriv)
#     return prefactor * total

#-------------------------------------------------------------------------------
# Define the Slm functions and use the DE to get the second derivative
#-------------------------------------------------------------------------------

def spectral_Slm(theta, phi, coeffs, c, ell, em, ess=-2):
    lmin = max(abs(ess), abs(em))
    nmax = int(50 + np.ceil(abs((3*c)/2. - c*c/250.)))
    if nmax % 2 == 1:
        nmax += 1
    nmin = min(ell - lmin, nmax)
    total = 0.
    for j in range(-nmin, nmax):
        Ylm = sphericalY(theta, 0, ell + j, em, ess)
        total += coeffs[j + nmin] * Ylm
    return np.sqrt(2*np.pi) * total * np.exp(1j * em * phi)


def spectral_dSlm(theta, phi, coeffs, c, ell, em, ess=-2):
    lmin = max(abs(ess), abs(em))
    nmax = int(50 + np.ceil(abs((3*c)/2. - c*c/250.)))
    if nmax % 2 == 1:
        nmax += 1
    nmin = min(ell - lmin, nmax)
    total = 0.
    for j in range(-nmin, nmax):
        dYlm = sphericalY_dtheta(theta, 0, ell + j, em, ess)
        total += coeffs[j + nmin] * dYlm
    return np.sqrt(2*np.pi) * total * np.exp(1j * em * phi)


# def spectral_d2Slm(theta, phi, coeffs, c, ell, em, ess=-2):
#     lmin = max(abs(ess), abs(em))
#     nmax = int(50 + np.ceil(abs((3*c)/2. - c*c/250.)))
#     if nmax % 2 == 1:
#         nmax += 1
#     nmin = min(ell - lmin, nmax)
#     total = 0.
#     for j in range(-nmin, nmax):
#         d2Ylm = sphericalY_d2theta(theta, 0, ell + j, em, ess)
#         total += coeffs[j + nmin] * d2Ylm
#     return np.sqrt(2*np.pi) * total * np.exp(1j * em * phi)


def find_d2Slm(x, Slm, dSlm, eigen, c, em, ess=-2):
    Alm = -c**2 + eigen + 2 * c * em
    V = ((c * x)*(c * x) - 2 * c * ess * x + ess + Alm -
         (em + ess * x)*(em + ess * x) / (1 - x*x))
    d2Slm = ((2 * x) * dSlm - (V * Slm)) / (1 - x*x)
    return d2Slm


def fast_spheroidal_harmonics_eq(aa, omega, ell, em, ess=-2):
    theta = np.pi / 2
    phi = 0
    x = 0  # this is the argument of the SWSH (not the inclination)
    c, __, __ = swsh_constants(aa, omega, ell, em, ess)
    # print(aa, omega, ell, em, ess)
    eigen, coeffs = solve_sparse_matrix(c, ell, em)
    Slm = spectral_Slm(theta, phi, coeffs, c, ell, em, ess)
    dSlm = spectral_dSlm(theta, phi, coeffs, c, ell, em, ess)
    d2Slm = find_d2Slm(x, Slm, dSlm, eigen, c, em, ess)
    # d2Slm = spectral_d2Slm(theta, phi, coeffs, c, ell, em, ess)
    return eigen, np.real(Slm), np.real(dSlm), np.real(d2Slm)

