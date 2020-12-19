from mpmath import mp, mpf, exp, gamma, sign
from .swsh_contfrac import contfracK
import numpy as np


def alpha_gamma(aa, omega, ell, em, nmax, ess=-2):
    """
    Computes alpha and gamma for use in the SWSHs

    Parameters:
        aa (float): spin
        omega (float): omega_mkn
        ell (int)
        em (int)
        nmax (int)

    Keyword Args:
        ess (int)

    Returns:
        alpha (float array): nmax x 1
        gamma (float array): nmax x 1
    """
    km = abs(em - ess) / 2
    kp = abs(em + ess) / 2
    c = aa * omega
    # eta = km + kp

    alpha = np.zeros(nmax + 1)
    alpha = mpf("1") * alpha
    gamma = np.zeros(nmax + 1)
    gamma = mpf("1") * gamma

    for j in range(nmax + 1):
        alpha[j] = (-2 * (j + 1) * (j + 2 * km + 1))
        gamma[j] = (2 * c * (j + km + kp + ess))
    return alpha, gamma


def find_beta(Alm, aa, omega, ell, em, nmax, ess=-2):
    """
    Computes beta for use in the SWSHs

    Parameters:
        Alm (float): separation value used to find eigenvalue
        aa (float): spin
        omega (float): omega_mkn
        ell (int): mode
        em (int): mode
        nmax (int): number of array values

    Keyword Args:
        ess (int): mode

    Returns:
        beta (float array)

    """
    km = abs(em - ess) / 2
    kp = abs(em + ess) / 2
    c = aa * omega
    # eta = km + kp

    beta = np.zeros(nmax + 1)
    beta = mpf("1") * beta
    for j in range(nmax + 1):
        beta[j] = (
            j * (j - 1) + 2 * j * (km + kp + 1 - 2 * c) -
            (2 * c * (2 * km + ess + 1) - (km + kp) * (km + kp + 1)) -
            (c**2 + ess * (ess + 1) + Alm)
        )

    return beta


def alpha_gamma2(aa, omega, ell, em, nmax, ess=-2):
    """
    Computes alpha2 and gamma2 for use in the SWSHs

    Parameters:
        aa (float): spin
        omega (float): omega_mkn
        ell (int)
        em (int)
        nmax (int)

    Keyword Args:
        ess (int)

    Returns:
        alpha (float array): nmax x 1
        gamma (float array): nmax x 1
    """
    km = abs(em - ess)
    kp = abs(em + ess)
    c = aa * omega
    # eta = km + kp

    alpha = np.zeros(nmax + 1)
    gamma = np.zeros(nmax + 1)
    alpha = mpf("1") * alpha
    gamma = mpf("1") * gamma

    for j in range(nmax + 1):
        alpha[j] = ((-4 * c * (j + kp + 1) * (j + km + 1) *
                     (j + (kp + km) / 2 + 1 + ess)) /
                    ((2 * j + kp + km + 2) * (2 * j + kp + km + 3))
                    )
        gamma[j] = ((4 * c * j * (j + kp + km) *
                     (j + (kp + km) / 2 - ess)) /
                    ((2 * j + kp + km - 1) * (2 * j + kp + km)))
    # print("alpha2 = ", alpha)
    # print("gamma2 = ", gamma)
    return alpha, gamma


def beta2(Alm, aa, omega, ell, em, nmax, ess=-2):
    """
    Computes beta2 for use in the SWSHs

    Parameters:
        Alm (float): separation value used to find eigenvalue
        aa (float): spin
        omega (float): omega_mkn
        ell (int): mode
        em (int): mode
        nmax (int): number of array values

    Keyword Args:
        ess (int): mode

    Returns:
        beta (float array)

    """
    km = abs(em - ess)
    kp = abs(em + ess)
    # print("km = ", km)
    # print("kp = ", kp)

    c = aa * omega
    # print("gamma = ", c)
    # eta = km + kp

    beta = np.zeros(nmax + 1)
    beta = mpf("1") * beta
    if ess != 0:
        for j in range(nmax + 1):  # reindexed from j to j-1 to work
            beta[j] = (
                Alm + ess * (ess + 1) + c**2 - (j + (kp + km) / 2) *
                (j + (kp + km) / 2 + 1) + (8 * em * ess**2 * c) /
                ((2 * j + kp + km) * (2 * j + kp + km + 2))
            )
    else:
        for j in range(nmax + 1):  # reindexed from j to j-1 to work
            beta[j] = (
                Alm + ess * (ess + 1) + c**2 - (j + (kp + km) / 2) *
                (j + (kp + km) / 2 + 1)
            )
    # print("beta2 = ", beta)
    return beta


def find_eigenvalue(aa, omega, ell, em, nmax, ess=-2):
    """
    Uses Secant method to calculate the zeros of the continued fraction for
    eigenvalues of the SWSH equations (Leaver's method).

    Part of this is taken from the BPTK

    Parameters:
        aa (float): spin
        omega (float): omega_mkn
        ell (int)
        em (int)
        nmax (int)

    Keyword Args:
        ess (int)

    Returns:
        xmin (float)
        eigen (float)
        alpha (float)
        gamma (float)
    """
    # TODO: Fix the error here. This doesn't show up in C++ somehow
    # This is a hotfix for some issue with computing inversion n in Python
    if em == 1:
        ell = ell - 1
        nInv = int(ell - abs(em))  # inversion n
    else:
        nInv = int(ell - abs(em))  # inversion n
    # print("nInv = ", nInv)

    def LHS(Alm):
        b_LHS = beta2(Alm, aa, omega, ell, em, nmax, ess)
        alpha2, gamma2 = alpha_gamma2(aa, omega, ell, em, nmax, ess)
        # print("alpha = ", alpha2)
        # print("gamma = ", gamma2)
        if nInv == 0:
            LHS_result = b_LHS[0]
            # print(Alm)
            # print("LHS = ", LHS_result)
        elif nInv == 1:
            LHS_result = (
                b_LHS[1] - alpha2[0] * gamma2[1] / b_LHS[0]
            )
            # print(Alm)
            # print(b_LHS[1])
            # print(alpha2[0])
            # print(b_LHS[0])
            # print("LHS = ", LHS_result)
        else:
            b_LHS_0 = b_LHS[nInv]
            a_LHS = np.zeros(nInv)
            a_LHS = mpf("1") * a_LHS
            bot = np.zeros(nInv)
            bot = mpf("1") * bot
            for j in range(0, nInv):
                a_LHS[j] = (-alpha2[(nInv) - (j + 1)] *
                            gamma2[(nInv) - (j + 1) + 1])
                # print((nInv) - (j + 1), a_LHS[j])
                bot[j] = b_LHS[nInv - (j + 1)]
                # print((nInv) - (j + 1), bot[j])
            LHS_result = b_LHS_0 + contfracK(a_LHS, bot)
        # print("Alm1 = ", Alm)
        # print("LHS_result1 = ", LHS_result)
        return LHS_result

    def RHS(Alm):
        b_RHS = beta2(Alm, aa, omega, ell, em, nmax)
        b_RHS = b_RHS[nInv + 1:-1]
        # print("b_RHS = ", b_RHS)
        alpha2, gamma2 = alpha_gamma2(aa, omega, ell, em, nmax, ess)
        a_RHS = np.zeros(len(b_RHS))
        a_RHS = mpf("1") * a_RHS
        for j in range(len(b_RHS)):
            a_RHS[j] = -alpha2[j + nInv] * gamma2[j + (nInv + 1)]
        # print("a_RHS = ", a_RHS)
        RHS_result = -contfracK(a_RHS, b_RHS)
        # print("Alm = ", Alm)
        # print("RHS_result1 = ", RHS_result)
        return RHS_result

    def fraction(Alm):
        # print("Alm = ", Alm)
        LHS_result = LHS(Alm)
        RHS_result = RHS(Alm)
        # print("LHS = ", LHS_result)
        # print("RHS = ", RHS_result)
        result = LHS_result - RHS_result
        # print("res = ", result)
        return result

    A = ell * (ell + 1) - ess * (ess + 1)  # initial guess
    c = aa * omega
    # println("c = ", c)
    alpha, gamma = alpha_gamma(aa, omega, ell, em, nmax, ess)
    # print('alpha = ', alpha)
    # print('gamma = ', gamma)
    if c == 0.0:
        return A, A, alpha, gamma

    x0 = A
    xmin = mp.findroot(fraction, x0, tol=mpf(str(10**-(mp.dps))))  # secant method
    # print("xmin2 = ", xmin)
    eigen = xmin + c**2 - 2 * c * em  # eigenvalues
    # print(eigen)
    return xmin, eigen, alpha, gamma


##-----------------------------------------------------------------------------
##
## functions used to determine the eigenvalue using Leaver's method
##
##-----------------------------------------------------------------------------

def swsh_constants(aa, omega, ell, em, ess=-2):
    km = abs(em - ess)
    kp = abs(em + ess)
    c = aa * omega
    if em == 1:
        ell = ell - 1
    nInv = ell - abs(em)
    return c, km, kp, nInv


def alpha_eigen(n, c, km, kp, ess=-2):
    return ((-4*c*(1 + km + n)*(1 + kp + n)*(1 + (km + kp)/2. + n + ess))/
             ((2 + km + kp + 2*n)*(3 + km + kp + 2*n)))


def beta_eigen(n, Alm, c, km, kp, em, ess=-2):
    return (Alm + c * c + ess*(1 + ess) -
            ((km + kp)/2. + n)*(1 + (km + kp)/2. + n) + 
            (8*c*ess*ess*em)/((km + kp + 2*n)*(2 + km + kp + 2*n)))


def gamma_eigen(n, c, km, kp, ess=-2):
    return ((4*c*n*(km + kp + n)*(-ess + (km + kp)/2. + n))/
             ((-1 + km + kp + 2*n)*(km + kp + 2*n)))


def numer_RHS(n, c, km, kp, ess=-2):
    alpha = alpha_eigen(n - 1, c, km, kp, ess)
    gamma = gamma_eigen(n, c, km, kp, ess)
    return -alpha * gamma


def denom_RHS(n, Alm, c, km, kp, em, ess=-2):
    return beta_eigen(n, Alm, c, km, kp, em, ess)


def cont_frac_eigen_RHS(n0, Alm, c, km, kp, em, ess=-2,
                        tol=mpf(str(10**(-mp.dps)))):
    xm1 = 1e30  # large value to keep the function from converging first iter
    ak = numer_RHS(n0, c, km, kp, ess)
    bk = denom_RHS(n0, Alm, c, km, kp, em, ess)

    Am2 = 1
    Bm2 = 0
    Am1 = 0
    Bm1 = 1
    A = bk * Am1 + ak * Am2
    B = bk * Bm1 + ak * Bm2
    x = A / B
    # // printf("A = %.17g \n", A);
    # // printf("B = %.17g \n", B);

    j = 1
    rel_err = abs(x - xm1) / abs(xm1)
    # // printf("rel_err = %.17g \n", rel_err);
    while (rel_err > tol):
        # // update values
        xm1 = x
        Am2 = Am1
        Am1 = A
        Bm2 = Bm1
        Bm1 = B

        # // find new values
        ak = numer_RHS(n0 + j, c, km, kp, ess)
        bk = denom_RHS(n0 + j, Alm, c, km, kp, em, ess)

        A = bk * Am1 + ak * Am2
        B = bk * Bm1 + ak * Bm2
        x = A / B
        # // printf("x = %.17g \n", x);

        j = j + 1
        # // printf("j = %d \n", j);
        # // printf("rel_err = %.17g \n", rel_err);
        rel_err = abs(x - xm1) / abs(xm1)
        if (j == 10000):
            # // infinite loop prevention
            print("Continued fraction did not converge.\n")
            return 0.

    return x


def numer_LHS(n, c, km, kp, nInv, ess=-2):
    # // printf("n = %d \n", n);
    alpha = alpha_eigen(nInv - n, c, km, kp, ess)
    # // printf("alpha = %.17g \n", alpha);
    gamma = gamma_eigen(nInv - n + 1, c, km, kp, ess)
    # // printf("gamma = %.17g \n", gamma);
    return -alpha * gamma
    # // printf("numer_LHS = %.17g \n", result);


def denom_LHS(n, Alm, c, km, kp, em, nInv, ess=-2):
    return beta_eigen(nInv - n, Alm, c, km, kp, em, ess)


def cont_frac_K_LHS(Alm, c, km, kp, em, nInv, ess=-2):
    n = 1
    numer = numer_LHS(nInv - n + 1, c, km, kp, nInv, ess)
    denom = denom_LHS(nInv - n + 1, Alm, c, km, kp, em, nInv, ess)
    result = numer / denom
    # // printf("numer = %.17g \n", numer);
    # // printf("denom = %.17g \n", denom);
    # // printf("result = %.17g \n", result);
    for n in range(2, nInv + 1):  # check this!
        old_result = result
        # // printf("nInv - n = %d \n", nInv - n);
        numer = numer_LHS(nInv - n + 1, c, km, kp, nInv, ess)
        denom = denom_LHS(nInv - n + 1, Alm, c, km, kp, em, nInv, ess)
        result = 1 / (old_result + denom) * numer
    return result


def fraction_eqn(Alm, c, km, kp, em, nInv, ess=-2, tol=mpf(str(10**(-mp.dps)))):
    b_LHS_0 = beta_eigen(0, Alm, c, km, kp, em, ess)

    # // LHS
    if (nInv == 0):
        LHS_result = b_LHS_0
    elif (nInv == 1):
        LHS_result = (beta_eigen(1, Alm, c, km, kp, em, ess) -
            alpha_eigen(0, c, km, kp, ess) * gamma_eigen(1, c, km, kp, ess) /
            b_LHS_0)
    else:
        LHS_result = (beta_eigen(nInv, Alm, c, km, kp, em, ess) +
                      cont_frac_K_LHS(Alm, c, km, kp, em, nInv, ess))
    # // printf("b_LHS_0 = %.17g \n", beta_eigen(nInv, Alm, c, km, kp, em, ess));
    # // printf("LHS_result = %.17g \n", LHS_result);
    # // RHS
    RHS_result = -cont_frac_eigen_RHS(nInv + 1, Alm, c, km, kp, em, ess, tol)

    # // printf("RHS_result = %.17g \n", RHS_result);
    # print("Alm =", Alm)
    # print("LHS_result =", LHS_result)
    # print("RHS_result =", RHS_result)
    return LHS_result - RHS_result


def swsh_eigen(c, km, kp, ell, em, nInv, ess=-2):
    tol = mpf(str(10**(-mp.dps)))  # set tolerance internally
    Alm0 = ell * (ell + 1) - ess * (ess + 1)
    if (c == 0):
        return Alm0, Alm0

    def fraction(Alm):
        return fraction_eqn(Alm, c, km, kp, em, nInv, ess, tol)

    xmin = mp.findroot(fraction, Alm0, tol=tol)
    eigen = xmin + c * c - 2 * c * em

    # print("xmin = ", xmin)
    return xmin, eigen


##
## TODO(aaron): Fix the following to work with the above
## everything below isn't currently used, but could be
## repurposed for a high precision code
##

# def alpha_gamma(aa, omega, ell, em, nmax, ess=-2):
#     km = np.abs(em - ess) / 2
#     kp = np.abs(em + ess) / 2
#     c = aa * omega
#     # eta = km + kp

#     alpha = np.zeros(nmax + 1)
#     alpha = mpf("1") * alpha
#     gamma = np.zeros(nmax + 1)
#     gamma = mpf("1") * gamma

#     for j in range(nmax + 1):
#         alpha[j] = (-2 * (j + 1) * (j + 2 * km + 1))
#         gamma[j] = (2 * c * (j + km + kp + ess))

#     return alpha, gamma


# def find_beta(Alm, aa, omega, ell, em, nmax, ess=-2):
#     km = np.abs(em - ess) / 2
#     kp = np.abs(em + ess) / 2
#     c = aa * omega

#     beta = np.zeros(nmax + 1)
#     beta = mpf("1") * beta
#     for j in range(nmax + 1):
#         beta[j] = (
#             j * (j - 1) + 2 * j * (km + kp + 1 - 2 * c) -
#             (2 * c * (2 * km + ess + 1) - (km + kp) * (km + kp + 1)) -
#             (c**2 + ess * (ess + 1) + Alm)
#         )

#     return beta


# def f_coeffs(aa, omega, ell, em, nmax, ess=-2):
#     c = aa * omega
#     # print(omega)
#     km = abs(em - ess) / 2
#     kp = abs(em + ess) / 2
#     Alm, eigen = swsh_eigen(c, km, kp, ell, em, ess)
#     alpha, gamma1 = alpha_gamma(aa, omega, ell, em, nmax, ess)
#     # print(gamma1)
#     # print(Alm, aa, omega, ell, em)
#     b = find_beta(Alm, aa, omega, ell, em, nmax, ess)
#     an = np.zeros(nmax) * mpf(1)
#     an = mpf("1") * an
#     an[0] = mpf("1")
#     # print(an[0])
#     # print(b[0])
#     # print(alpha[0])
#     an[1] = -b[0] * an[0] / alpha[0]
#     # print(an[1])
#     for j in range(2, nmax):
#         # print(b[j - 1])
#         # print(an[j - 1])
#         # print(gamma1[j - 1])
#         # print(an[j - 2])
#         # print(alpha[j - 1])
#         an[j] = ((-b[j - 1] * an[j - 1] - gamma1[j - 1] * an[j - 2]) /
#                  alpha[j - 1])
#         # print((-b[j - 1] * an[j - 1] - gamma1[j - 1] * an[j - 2]))
#         # print(alpha[j - 1])
#         # print(str(j), an[j])

#     # print(an)

#     an_dotprods = np.zeros(nmax)
#     an_dotprods = mpf("1") * an_dotprods
#     for j in range(nmax):
#         a1 = 1 + j + 2 * km
#         a1 = mpf(str(a1))
#         b1 = j + 2 * (1 + kp + km)
#         b1 = mpf(str(b1))
#         c1 = 4 * c
#         c1 = mpf(str(c1))

#         prefactor = (2**j * rf(j + 2 * (1 + km + kp), -2 * kp - 1) *
#                      hyp1f1(a1, b1, c1))
#         # print(prefactor)
#         # prefactor is good!
#         # print(str(j), prefactor)
#         an_forward = an[0:j + 1]
#         # print(str(j), an_forward)
#         an_dotprods[j] = prefactor * np.dot(an_forward, an_forward[::-1])
#         # print(np.dot(an_forward, an_forward[::-1]))
#         # print(str(j), np.dot(an_forward, an_forward[::-1]))
#         # print(hyp1f1(a1, b1, c1))
#     # print(an_dotprods)
#     norm = (mpf("2")**(1 + 2 * km + 2 * kp) * exp(-mpf("2") * c) *
#             gamma(1 + 2 * kp) * np.sum(an_dotprods))**(1 / 2)
#     # print(an_dotprods)
#     # print('thing = ', np.sum(an_dotprods))
#     # print("sum = ", np.sum(an_dotprods))
#     # print("norm = ", norm)
#     f = an / norm
#     # print(an)
#     # print(norm)
#     # for j in range(len(f)):
#         # print(f[j])
#     # print(f)
#     return f, Alm, eigen


# def harmonics(x, aa, omega, ell, em, coeffs, nmax, ess=-2):
#     c = aa * omega
#     km = abs(em - ess) / 2
#     kp = abs(em + ess) / 2
#     vectors = np.zeros(nmax)
#     vectors = mpf("1") * vectors
#     for j in range(nmax):
#         vectors[j] = (1 + x)**j

#     harmonics = (
#         exp(c * x) * (1 + x)**km * (1 - x)**kp *
#         np.dot(vectors, coeffs)
#     )
#     integrand = np.abs(harmonics)**2
#     return integrand, harmonics


# def Slm_leaver(x, aa, omega, ell, em, coeffs, nmax, ess=-2):
#     """
#     Return the Spheroidal Eigenfunctions. This agrees with BPTK, but the
#     precision check in their SWSH package doesn't work.
#     """
#     c = aa * omega
#     km = abs(em - ess) / 2
#     kp = abs(em + ess) / 2
#     vectors = np.zeros(nmax)
#     vectors = mpf("1") * vectors
#     for j in range(nmax):
#         vectors[j] = (1 + x)**j
#         # print(vectors[j])
#     Slm = exp(c * x) * (1 + x)**km * (1 - x)**kp * np.dot(vectors, coeffs)
#     Slm = Slm * sign(omega)  # remove sign ambiguity
#     return Slm


# def dSlm_leaver(x, aa, omega, ell, em, coeffs, Slm, nmax, ess=-2):
#     """
#     Calculates the derivative of the SWSH functions using an analytical
#     derivative of the Leaver's method Slm.
#     """
#     c = aa * omega
#     km = abs(em - ess) / 2
#     kp = abs(em + ess) / 2
#     vectors = np.zeros(nmax)
#     vectors = mpf("1") * vectors
#     for j in range(nmax):
#         vectors[j] = (1 + x)**(j)
#     y = np.dot(coeffs, vectors)
#     vectors = np.zeros(nmax)
#     vectors = mpf("1") * vectors
#     for j in range(nmax):
#         vectors[j] = j * (1 + x)**(j - 1)
#     y_prime = np.dot(coeffs, vectors)

#     lnS_prime = c + km * (1 / (1 + x)) - kp * (1 / (1 - x)) + y_prime / y
#     S_prime = Slm * lnS_prime
#     return S_prime


# def d2Slm_leaver(x, aa, omega, ell, em, Alm, Slm, Slmd, ess=-2):
#     """
#     Uses the second order ODE to calculate the second derivative of Slm.
#     """
#     c = aa * omega
#     # km = abs(em - ess) / 2
#     # kp = abs(em + ess) / 2

#     V = ((c * x)**2 - 2 * c * ess * x + ess + Alm - (em + ess * x)**2 /
#          (1 - x**2))
#     S_doubleprime = -((2 * x) * Slmd - (V * Slm)) / (1 - x**2)
#     return S_doubleprime


# def find_swsh_eq_leaver(aa, omega, ell, em, ess=-2):
#     """
#     Culmination of all routines in this file to make these easily callable.
#     This will return eigenvalues, S, dSdth, dSdth2.
#     """
#     nmax = 100  # TODO: find a way to remove this

#     x = mpf(0)
#     # x = cos(pi / 2)  # specialize to equatorial orbits
#     coeffs, Alm, eigen = f_coeffs(aa, omega, ell, em, nmax, ess)
#     # print(coeffs)
#     Slm = Slm_leaver(x, aa, omega, ell, em, coeffs, nmax, ess)
#     Slmd = dSlm_leaver(x, aa, omega, ell, em, coeffs, Slm, nmax, ess)
#     Slmdd = d2Slm_leaver(x, aa, omega, ell, em, Alm, Slm, Slmd, ess)
#     return eigen, -Slm, Slmd, Slmdd


# def swsh_eigenvalue(aa, omega, ell, em, ess=-2):

#     c, km, kp, nInv = swsh_constants(aa, omega, ell, em, ess)
#     __, eigen = swsh_eigen(c, km, kp, ell, em, nInv, ess)

#     return eigen


# if __name__ == "__main__":
#     mp.dps = 50
#     aa = mpf("0.998")
#     omega = mpf("0.1220182085786500841797685807337234729735480")
#     ell = 2
#     em = 2
#     tol=mpf(str(10**(-mpf(mp.dps))))
#     eigen = swsh_eigenvalue(aa, omega, ell, em)

