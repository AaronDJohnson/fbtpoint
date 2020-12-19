import numpy as np


def eq_source(chi, sign, slr, ecc, aa, omega, em, Lz, En,
              Slm, Slmd, Slmdd, Rin, dRdr, dRdr2):
    """
    This finds the I that are used in the Z function. 'sign' specifies
    I+ or I-.
    """
    Xn = -aa * En + Lz
    eta = aa * omega - em
    irootpi = 1 / np.sqrt(np.pi)
    sigma = 1 + ecc * np.cos(chi)
    aa2 = aa * aa
    sigma2 = sigma * sigma
    slr2 = slr * slr
    Xn2 = Xn * Xn
    u = sigma / slr
    denom = (1 - 2 * u + aa2 * u**2)

    J = (1 - 2 / slr * sigma + aa2 / slr2 * sigma2)
    # print('chi = ' + str(chi))
    # print('J = ' + str(J))
    V_t = (
        aa2 * En - 2 * aa * Xn / slr * sigma + En * slr2 /
        sigma2)
    # print('Vt = ' + str(V_t))
    V_r = Xn2 + aa2 + 2 * aa * Xn * En - 2 * Xn2 / slr * (2 + sigma)

    # print('Vr = ' + str(V_r))
    # print(chi, sign)
    # print('J, V_t, V_r = ' + str([J, V_t, (V_r)**(1 / 2)]))

    Cnn = (
        J / (4 * slr**4 * V_t) *
        (slr2 * En - aa * Xn * sigma2 + ecc * slr * sign * np.sin(chi) *
         np.sqrt(V_r))**2
    )
    # print('Cnn = ' + str(Cnn))
    # print(chi, sign)
    Cmn = (
        1j * Xn * J / (2**(3 / 2) * slr**3 * V_t) * sigma *
        (slr2 * En - aa * Xn * sigma2 + ecc * slr *
         sign * np.sin(chi) * np.sqrt(V_r))
    )
    # print('Cmn = ' + str(Cmn))
    Cmm = (-Xn2 * J / (2 * slr2 * V_t) * sigma2)
    # print('Cmm = ' + str(Cmm))
    Amn0 = (
        2 * irootpi * Cmn / (u * denom**2) *
        (2 * aa2 * u**3 + (1j * aa * eta - 4) * u**2 + 2 * u +
         1j * omega) * (Slmd + eta * Slm)
    )
    # print('Amn0 = ' + str(Amn0))
    Amm0 = (
        1 / np.sqrt(2) * irootpi * Cmm * Slm / (u**2 * denom**2) *
        (-2 * 1j * aa**3 * eta * u**5 + aa * eta *
         (6 * 1j + aa * eta) * u**4 - 4 * 1j * aa * eta *
         u**3 + 2 * omega * (1j + aa * eta) * u**2 - 2 * 1j * omega * u +
         omega**2)
    )
    # print('Amm0 = ' + str(Amm0))
    Amn1 = (
        2 * irootpi * Cmn / (u * denom) * (Slmd + eta * Slm)
    )
    # print('Amn1= ' + str(Amn1))
    Amm1 = (
        -np.sqrt(2) * irootpi * Cmm * Slm / (u**2 * denom) *
        (aa**2 * u**3 + (1j * aa * eta - 2) * u**2 + u + 1j * omega)
    )
    # print('Amm1= ' + str(Amm1))
    Amm2 = (
        -1 / np.sqrt(2) * irootpi * Cmm * Slm / u**2
    )
    # print('Amm2= ' + str(Amm2))
    Ann0 = (
        -np.sqrt(2) * irootpi * Cnn / denom**2 *
        (-2 * 1j * aa * (Slmd + eta * Slm) * u + Slmdd + 2 * eta *
         Slmd + (eta**2 - 2) * Slm)
    )
    # print('Ann0=' + str(Ann0))
    # print('Slmdd = ' + str(Slmdd))

    I_inf = (
        Rin * (Ann0 + Amn0 + Amm0) - dRdr * (Amn1 + Amm1) + dRdr2 * Amm2
    )
    # print('chi = ', chi)
    # print('I_inf = ' + str(I_inf))
    return J, V_t, V_r, I_inf
