import numpy as np
from mpmath import mpc, sqrt, gamma, cos, acos, re, im, mp, pi

#-------------------------------------------------------------------------------
# TODO List:
# * Make this function depend on some tolerance rather than a static number of
#   points (currently set at 200 -- this is too many)
#
# * Something is making this cutoff at 55 digits of accuracy. While this is
#   enough for our purposes, it would be nice to do better.
#
# * Find out where this is from so that we can use numpy/scipy
#   functions to make this a little faster and maybe drop high precision
#
#-------------------------------------------------------------------------------

def monodromy_nu(aa, omega, eigen, ell, em, nmax, ess=-2):
    """
    This is the monodromy method used in the BPTK. It seems to work well.

    Small speed gains to be had by porting this to C++ plus GMP OR using the
    result from the arb package that has already been programmed
    """

    nhalf = nmax // 2 + 1
    epsilon = 2 * omega
    kappa = sqrt(1 - aa**2)
    tau = (epsilon - em * aa) / kappa

    # parameters for confluent Heun eqn.
    gamma_ch = 1 - ess - 1j * epsilon - 1j * tau
    delta_ch = 1 + ess + 1j * epsilon - 1j * tau
    epsilon_ch = 2 * 1j * epsilon * kappa
    alpha_ch = 2 * epsilon * kappa * (1j - 1j * ess - epsilon + tau)
    q_ch = -(-ess * (1 + ess) + epsilon**2 + 1j * (-1 + 2 * ess) * epsilon *
             kappa - eigen - tau * (1j + tau))

    mu_1c = alpha_ch / epsilon_ch - (gamma_ch + delta_ch)
    mu_2c = -(alpha_ch / epsilon_ch)

    # recurrence relations here
    a1_vec = np.zeros(nmax + 3, dtype=complex) * mpc("1")
    a1_vec[0] = mpc("0")
    a1_vec[1] = mpc("1")

    a2_vec = np.zeros(nmax + 3, dtype=complex) * mpc("1")
    a2_vec[0] = mpc("0")
    a2_vec[1] = mpc("1")
        
    for i in range(2, nmax + 3):
        n = i - 1
        a1_vec[i] = (((alpha_ch - (-1 + n + delta_ch) * epsilon_ch) *
                      (alpha_ch - (-2 + n + gamma_ch + delta_ch) *
                       epsilon_ch)) /
                     (n * epsilon_ch) * a1_vec[i - 2] - 1 /
                     (n * epsilon_ch**2) *
                     (alpha_ch**2 + alpha_ch *
                      epsilon_ch *
                      (1 - 2 * n - gamma_ch - delta_ch + epsilon_ch) +
                      epsilon_ch**2 *
                      (n**2 - q_ch + n *
                       (-1 + gamma_ch + delta_ch -
                        epsilon_ch) + epsilon_ch - delta_ch * epsilon_ch)) *
                     a1_vec[i - 1])
        a2_vec[i] = (-(((alpha_ch + (-2 + n) * epsilon_ch) *
                        (alpha_ch +
                         (-1 + n - gamma_ch) * epsilon_ch)) /
                       (n * epsilon_ch)) * a2_vec[i - 2] +
                     1 / (n * epsilon_ch**2) *
                     (alpha_ch**2 + (n**2 - q_ch + gamma_ch +
                                     delta_ch - n *
                                     (1 + gamma_ch + delta_ch - epsilon_ch) -
                                     epsilon_ch) *
                      epsilon_ch**2 + alpha_ch * epsilon_ch *
                      (-1 + 2 * n - gamma_ch -
                       delta_ch + epsilon_ch)) * a2_vec[i - 1])
        # ratio_a1 = abs(a1_vec[i]) / abs(a1_vec[i - 1])
        # ratio_a2 = abs(a2_vec[i]) / abs(a2_vec[i - 1])
        # if ratio_a1 < tol and ratio_a2 < tol:
        #     break
        # elif i == nmax + 2:
        #     print('a1 or a2 arrays failed to converge in renorm_nu. ' +
        #           'Increase nmax')

    a1_vec = a1_vec[1:-1]
    a2_vec = a2_vec[1:-1]

    pochhammerp1m2 = np.zeros(nmax + 3, dtype=complex) * mpc("1")
    pochhammerp1m2[0] = mpc("1")

    pochhammerm1p2 = np.zeros(nmax + 3, dtype=complex) * mpc("1")
    pochhammerm1p2[0] = mpc("1")

    for i in range(1, nmax + 3):
        n = i
        pochhammerp1m2[i] = (-mu_2c + mu_1c + n - 1) * pochhammerp1m2[i - 1]
        pochhammerm1p2[i] = (mu_2c - mu_1c + n - 1) * pochhammerm1p2[i - 1]

    pochhammerp1m2 = pochhammerp1m2[0:-2]
    pochhammerm1p2 = pochhammerm1p2[0:-2]

    # summations:
    a1sum = np.zeros(nhalf, dtype=complex) * mpc("1")
    a2sum = np.zeros(nhalf, dtype=complex) * mpc("1")

    for j in range(nhalf):
        a1sum[j] = a1_vec[j] * pochhammerp1m2[nmax - j]
        a2sum[j] = (-1)**j * a2_vec[j] * pochhammerm1p2[nmax - j]

    a1 = gamma(-mu_2c + mu_1c) * np.sum(a1sum)
    a2 = gamma(mu_2c - mu_1c) * np.sum(a2sum)

    # print("a1 = ", a1)
    # print("a2 = ", a2)

    c2pn = (cos(pi * (mu_1c - mu_2c)) + (2 * pi**2) / (a1 * a2) *
            (-1)**(nmax - 1) * a1_vec[-1] * a2_vec[-1])

    # print('c2pn = ', c2pn)

    # print(np.imag(acos(np.real(c2pn)) / (2 * pi) * 1j))

    # nu_0 = ell - acos(c2pn) / (2 * pi)
    # print(acos(np.real(c2pn)) / (2 * pi))
    # print("nu_0 = ", nu_0)

    if re(c2pn) >= -1 and re(c2pn) <= 1:
        # region = "re"
        nu_0 = ell - acos(re(c2pn)) / (2 * pi)
        # print("nu_0 = ", nu_0)
    elif re(c2pn) < -1:
        # region = "half_int"
        nu_0 = 1/2 - im(acos(re(c2pn)) / (2 * pi)) * 1j
        # print("nu_0 = ", nu_0)
    elif re(c2pn) > 1:
        # region = "integer"
        nu_0 = -1j * im(acos(re(c2pn)) / (2 * pi))
    # print("nu_0 = ", nu_0)
    # print(region)
    return nu_0


def find_nu(aa, omega, eigen, ell, em, ess=-2, tol=10**(-mp.dps)):
    nmax = 25
    nu0 = 1e30
    nu = monodromy_nu(aa, omega, eigen, ell, em, nmax, ess)
    while abs(nu - nu0) / abs(nu0) > tol:
        nu0 = nu
        nmax *= 2
        nu = monodromy_nu(aa, omega, eigen, ell, em, nmax, ess)
        if nmax >= 600:
            print("nu failed to converge.")
            return 0
    return nu











