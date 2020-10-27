from mpmath import sqrt, ellipk, ellippi, ellipe, mp


def mp_radial_roots(En, Q, aa, slr, ecc, M=1):
    """
    Roots of the radial geodesic eq, these correspond to turning points of the
    orbits.

    Parameters:
        En (float): energy
        Q (float): Carter constant
        aa (float): spin parameter (0, 1)
        slr (float): semi-latus rectum [6, inf)
        ecc (float): eccentricity [0, 1)
    
    Keyword Args:
        M (float): mass of the large body

    Returns:
        r1 (float): apastron
        r2 (float): periastron
        r3 (float)
        r4 (float)
    """
    try:
        prec = mp.prec
        mp.prec += 10
        En2 = En * En

        mp.prec += 10
        r1 = slr / (1 - ecc)
        r2 = slr / (1 + ecc)

        mp.prec += 10
        AplusB = (2 * M) /  (1 - En2) - (r1 + r2)
        AB = (aa * aa * Q) / ((1 - En2) * r1 * r2)
        mp.prec += 10
        r3 = (AplusB + sqrt((AplusB * AplusB - 4 * AB))) / 2
        r4 = 0
        if r3 != 0:
            r4 = AB / r3
        return r1, r2, r3, r4
    finally:
        mp.prec = prec


def mp_polar_roots(En, Lz, aa, slr, x):
    """
    Roots of the polar geodesic eq, these correspond to turning points of the
    orbits.

    Parameters:
        En (float): energy
        Lz (float): angular momentum
        aa (float): spin parameter (0, 1)
        slr (float): semi-latus rectum [6, inf)
        x (float): inclination value given by cos(theta_inc) (0, 1]
                   negative x -> retrograde
                   positive x -> prograde

    Returns:
        zp (float)
        zm (float)
    """
    L2 = Lz * Lz

    zm = sqrt(1 - x * x)
    if abs(zm) == 1:
        zp = 0
    else:
        zp = sqrt(aa * aa * (1 - En * En) + L2 / (1 - zm * zm))
    return zp, zm


def mino_freqs_sc(slr, ecc, x):
    """
    Mino frequencies for the SC case (aa = 0)

    Parameters:
        slr (float): semi-latus rectum
        ecc (float): eccentricity
        x (float): inclincation

    Returns:
        ups_r (float): radial Mino frequency
        ups_theta (float): polar Mino frequency
        ups_phi (float): azimuthal Mino frequency
        gamma (float): time Mino frequency
    """
    pi = mp.pi
    ups_r = ((pi*sqrt(-((slr*(-6 + 2*ecc + slr))/(3 + ecc**2 - slr))))/
             (2*ellipk((4*ecc)/(-6 + 2*ecc + slr))))
    ups_theta = slr/sqrt(-3 - ecc**2 + slr)
    ups_phi = (slr*x)/(sqrt(-3 - ecc**2 + slr)*abs(x))
    gamma = ((sqrt((-4*ecc**2 + (-2 + slr)**2)/(slr*(-3 - ecc**2 + slr)))*
            (8 + (-(((-4 + slr)*slr**2*(-6 + 2*ecc + slr)*
            ellipe((4*ecc)/(-6 + 2*ecc + slr)))/(-1 + ecc**2)) + 
            (slr**2*(28 + 4*ecc**2 - 12*slr + slr**2)*
            ellipk((4*ecc)/(-6 + 2*ecc + slr)))/(-1 + ecc**2) - 
            (2*(6 + 2*ecc - slr)*(3 + ecc**2 - slr)*slr**2*
            ellippi((2*ecc*(-4 + slr))/((1 + ecc)*(-6 + 2*ecc + slr)),
            (4*ecc)/(-6 + 2*ecc + slr)))/((-1 + ecc)*(1 + ecc)**2) + 
            (4*(-4 + slr)*slr*(2*(1 + ecc)*ellipk((4*ecc)/(-6 + 2*ecc + slr)) + 
            (-6 - 2*ecc + slr)*
            ellippi((2*ecc*(-4 + slr))/((1 + ecc)*(-6 + 2*ecc + slr)),
            (4*ecc)/(-6 + 2*ecc + slr))))/(1 + ecc) + 
            2*(-4 + slr)**2*((-4 + slr)*ellipk((4*ecc)/(-6 + 2*ecc + slr)) - 
            ((6 + 2*ecc - slr)*slr*
            ellippi((16*ecc)/(12 + 8*ecc - 4*ecc**2 - 8*slr + slr**2),
            (4*ecc)/(-6 + 2*ecc + slr)))/(2 + 2*ecc - slr)))/
            ((-4 + slr)**2*ellipk((4*ecc)/(-6 + 2*ecc + slr)))))/2.)

    return ups_r, ups_theta, ups_phi, gamma


def mino_freqs_kerr(r1, r2, r3, r4, En, Lz, Q, aa, slr, ecc, x, M=1):
    """
    Mino frequencies for the Kerr case (aa != 0)

    Parameters:
        r1 (float): radial root
        r2 (float): radial root
        r3 (float): radial root
        r4 (float): radial root
        En (float): energy
        Lz (float): angular momentum
        Q (float): Carter constant
        aa (float): spin
        slr (float): semi-latus rectum
        ecc (float): eccentricity
        x (float): inclincation

    Keyword Args:
        M (float): mass

    Returns:
        ups_r (float): radial Mino frequency
        ups_theta (float): polar Mino frequency
        ups_phi (float): azimuthal Mino frequency
        gamma (float): time Mino frequency
    """
    # En, Lz, Q = calc_eq_constants(aa, slr, ecc)
    # r1, r2, r3, r4 = radial_roots(En, Q, aa, slr, ecc, M)
    pi = mp.pi
    L2 = Lz * Lz
    aa2 = aa * aa
    En2 = En * En
    M2 = M * M

    # polar pieces
    zm = 1 - x * x
    # a2zp = (L2 + aa2*(-1 + En2)*(-1 + zm))/((-1 + En2)*(-1 + zm))
    eps0zp = -((L2 + aa2*(-1 + En2)*(-1 + zm))/(L2*(-1 + zm)))
    zmOverzp = (aa2*(-1 + En2)*(-1 + zm)*zm)/(L2 + aa2*(-1 + En2)*(-1 + zm))

    kr = sqrt(((r1 - r2)*(r3 - r4))/((r1 - r3)*(r2 - r4)))
    ktheta = sqrt(zmOverzp)

    kr2 = kr * kr
    ktheta2 = ktheta * ktheta

    ellipticK_r = ellipk(kr2)
    ellipticK_theta = ellipk(ktheta2)

    rp = M + sqrt(M2 - aa2)
    rm = M - sqrt(M2 - aa2)
    hr = (r1 - r2)/(r1 - r3)
    hp = ((r1 - r2)*(r3 - rp))/((r1 - r3)*(r2 - rp))
    hm = ((r1 - r2)*(r3 - rm))/((r1 - r3)*(r2 - rm))

    ellipticPi_hmkr = ellippi(hm, kr2)
    ellipticPi_hpkr = ellippi(hp, kr2)
    ellipticPi_hrkr = ellippi(hr, kr2)
    ellipticPi_zmktheta = ellippi(zm, ktheta2)
    ellipticE_kr = ellipe(kr2)
    ellipticE_ktheta = ellipe(ktheta2)


    ups_r = (pi * sqrt((1 - En2)*(r1 - r3)*(r2 - r4)))/(2 * ellipticK_r)
    ups_theta = (sqrt(eps0zp) * Lz * pi)/(2. * ellipticK_theta)
    ups_phi = ((2 * aa * ups_r * (-(((-(aa * Lz) + 2 * En * M * rm)*
            (ellipticK_r - ((r2 - r3) * ellipticPi_hmkr)/(r2 - rm)))/
            (r3 - rm)) + ((-(aa * Lz) + 2*En * M * rp)*
            (ellipticK_r - ((r2 - r3) * ellipticPi_hpkr)/(r2 - rp)))/
            (r3 - rp)))/(pi * sqrt((1 - En2)*(r1 - r3)*(r2 - r4))*(-rm + rp)) + 
            (2 * ups_theta * ellipticPi_zmktheta) / (sqrt(eps0zp)*pi))
    gamma = (4*En*M2 + (2*En*ups_theta*(L2 + aa2*(-1 + En2)*(-1 + zm))*
            (-ellipticE_ktheta + ellipticK_theta))/
            ((-1 + En2)*sqrt(eps0zp)*Lz*pi*(-1 + zm)) + 
            (2*ups_r*((2*M*(-(((-2*aa2*En*M + (-(aa*Lz) + 4*En*M2)*rm)*
            (ellipticK_r - ((r2 - r3) * ellipticPi_hmkr)/(r2 - rm)))/
            (r3 - rm)) + ((-2*aa2*En*M + (-(aa*Lz) + 4*En*M2)*rp)*
            (ellipticK_r - ((r2 - r3)*ellipticPi_hpkr)/(r2 - rp)))/
            (r3 - rp)))/(-rm + rp) + 
            2*En*M*(r3*ellipticK_r + (r2 - r3)*ellipticPi_hrkr) + 
            (En*((r1 - r3)*(r2 - r4)*ellipticE_kr + 
            (-(r1*r2) + r3*(r1 + r2 + r3))*ellipticK_r + 
            (r2 - r3)*(r1 + r2 + r3 + r4)*ellipticPi_hrkr))/2.))/
            (pi*sqrt((1 - En2)*(r1 - r3)*(r2 - r4))))

    return ups_r, ups_theta, ups_phi, gamma


def mp_mino_freqs(r1, r2, r3, r4, En, Lz, Q, aa, slr, ecc, x):
    """
    Mino frequency calculation using mpmath with arbitrary precision

    Parameters:
        r1 (float): radial root
        r2 (float): radial root
        r3 (float): radial root
        r4 (float): radial root
        En (float): energy
        Lz (float): angular momentum
        Q (float): Carter constant
        aa (float): spin parameter (0, 1)
        slr (float): semi-latus rectum [6, inf)
        ecc (float): eccentricity [0, 1)
        x (float): inclination value given by cos(theta_inc) (0, 1]
                   negative x -> retrograde
                   positive x -> prograde

    Returns:
        ups_r (float): radial frequency
        ups_theta (float): theta frequency
        ups_phi (float): phi frequency
        gamma (float): time frequency
    """
    if aa == 0:
        ups_r, ups_theta, ups_phi, gamma = mino_freqs_sc(slr, ecc, x)
        return ups_r, ups_theta, ups_phi, gamma
    else:
        ups_r, ups_theta, ups_phi, gamma = mino_freqs_kerr(r1, r2, r3, r4, En,
                                                        Lz, Q, aa, slr, ecc, x)
        return ups_r, ups_theta, ups_phi, gamma


def mp_boyer_freqs(ups_r, ups_theta, ups_phi, gamma, aa, slr, ecc, x, M=1):
    """
    Boyer frequency calculation using mpmath with arbitrary precision

    Parameters:
        ups_r (float): radial frequency
        ups_theta (float): theta frequency
        ups_phi (float): phi frequency
        gamma (float): time frequency
        aa (float): spin parameter (0, 1)
        slr (float): semi-latus rectum [6, inf)
        ecc (float): eccentricity [0, 1)
        x (float): inclination value given by cos(theta_inc) (0, 1]
                   negative x -> retrograde
                   positive x -> prograde

    Keyword Args:
        M (float): mass of large body

    Returns:
        omega_r (float): radial boyer lindquist frequency
        omega_theta (float): theta boyer lindquist frequency
        omega_phi (float): phi boyer lindquist frequency
    """

    omega_r = ups_r / gamma
    omega_theta = ups_theta / gamma
    omega_phi = ups_phi / gamma

    return omega_r, omega_theta, omega_phi


def mp_find_omega(omega_r, omega_theta, omega_phi, em, kay, en, M=1):
    """
    Uses Boyer frequencies to find omega_mkn

    Parameters:
        omega_r (float): radial boyer lindquist frequency
        omega_theta (float): theta boyer lindquist frequency
        omega_phi (float): phi boyer lindquist frequency
        em (int): phi mode
        kay (int): theta mode
        en (int): radial mode

    Keyword Args:
        M (float): mass of large body

    Returns:
        omega (float): frequency of motion
    """

    return en * omega_r + em * omega_phi + kay * omega_theta








