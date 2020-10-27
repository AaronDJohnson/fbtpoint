# from constants import fp_calc_constants
# from frequencies import fp_mino_freqs, fp_boyer_freqs, fp_polar_roots
# from frequencies import fp_find_omega, fp_radial_roots

from mpmath import cos, sin, acos, asin, mp, pi
from mpmath import sqrt, floor, quad, ellipfun
from mpmath import ellipk, ellipf, ellipe, ellippi

# ------------------------------------------------------------------------------
# TODO List:
#
# * test all known cases and see if we get the correct answer
# * create automatic testing profile and input solutions
# * port this to C++ so that we get answers fast!
#   - this may not lead to arbitrary precision (do we need high precision?)
#   - if we do need arbitrary precision, implement spectral integration
# * make math functions use numpy/scipy at low precision
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# These may or may not be needed
# ------------------------------------------------------------------------------


def jacobi_am_f(x, m):
    """Borrowed from sagemath."""
    ctx = mp
    prec = ctx.prec
    try:
        # print(x)
        # print(m)
        x = ctx.convert(x)
        m = ctx.convert(m)
        if not isinstance(x, ctx.mpf) or not isinstance(m, ctx.mpf):
            raise ValueError('arguments must be real')
        if abs(m) == 1:
            # gd(x)
            ctx.prec += 10
            tanhx = ctx.tanh(x)
            ctx.prec += 10
            return ctx.asin(tanhx)
        elif abs(m) > 1:
            ctx.prec += 10
            # Real values needed for atan2; as per "Handbook of Elliptic
            # Integrals for Engineers and Scientists" 121.02, sn is real for
            # real x. The imaginary components can thus be safely discarded.
            snx = ctx.ellipfun('sn', x, m).real
            cnx = ctx.ellipfun('cn', x, m).real
            ctx.prec += 10
            return ctx.atan2(snx, cnx)
        else:
            ctx.prec += 10
            K = ctx.ellipk(m)
            if abs(x) <= K:
                snx = ctx.ellipfun('sn', x, m).real
                cnx = ctx.ellipfun('cn', x, m).real
                ctx.prec += 10
                return ctx.atan2(snx, cnx)
            else:
                # Do argument reduction on x to end up with z = x - 2nK, with
                # abs(z) <= K
                ctx.prec += 10
                tK = 2 * K
                ctx.prec += 10
                n = ctx.floor(x / tK)
                ctx.prec += 10
                tnK = n * tK
                npi = n * ctx.pi()
                ctx.prec += 10
                z = x - tnK
                ctx.prec += 10
                # z (and therefore sn(z, m) and cn(z, m)) is real because K(m)
                # is real for abs(m) <= 1.
                snz = ctx.ellipfun('sn', z, m).real
                cnz = ctx.ellipfun('cn', z, m).real
                ctx.prec += 10
                return ctx.atan2(snz, cnz) + npi
    finally:
        ctx.prec = prec


def calc_J(chi, En, Lz, Q, aa, slr, ecc):
    """
    Schmidt's J function.

    Parameters:
        chi (float): radial angle
        En (float): energy
        Lz (float): angular momentum
        Q (float): Carter constant
        aa (float): spin
        slr (float): semi-latus rectum
        ecc (float): eccentricity

    Returns:
        J (float)
    """
    En2 = En * En
    ecc2 = ecc * ecc
    aa2 = aa * aa
    Lz2 = Lz * Lz
    slr2 = slr * slr

    eta = 1 + ecc * cos(chi)
    eta2 = eta * eta

    J = ((1 - ecc2)*(1 - En2) + 2*(1 - En2 - (1 - ecc2)/slr)*eta + 
        (((3 + ecc2)*(1 - En2))/(1 - ecc2) + 
        ((1 - ecc2)*(aa2*(1 - En2) + Lz2 + Q))/slr2 - 4/slr)*
        eta2)

    return J


def calc_wr(psi, ups_r, En, Lz, Q, aa, slr, ecc):
    """
    w_r = ups_r * lambda as a function of the radial angle

    Parameters:
        psi (float): radial angle
        ups_r (float): radial mino frequency
        En (float): energy
        Lz (float): angular momentum
        Q (float): Carter constant
        aa (float): spin
        slr (float): semi-latus rectum
        ecc (float): eccentricity

    Returns:
        w_r (float)
    """

    # TODO: this function has an analytical integral solution
    # but the integral doesn't work at pi where the integral equals pi
    def wr_integrand(chi):
        sqrt_J = sqrt(calc_J(chi, En, Lz, Q, aa, slr, ecc))
        result = ups_r / sqrt_J
        return result
    answer = quad(wr_integrand, [0, psi])
    return (1 - ecc**2) / slr * answer


def calc_dwr_dpsi(psi, ups_r, En, Lz, Q, aa, slr, ecc):
    """
    Derivative of w_r with respect to psi

    Parameters:
        psi (float): radial angle
        ups_r (float): radial mino frequency
        En (float): energy
        Lz (float): angular momentum
        Q (float): Carter constant
        aa (float): spin
        slr (float): semi-latus rectum
        ecc (float): eccentricity

    Returns:
        dwr_dpsi (float)
    """
    sqrt_J = sqrt(calc_J(psi, En, Lz, Q, aa, slr, ecc))
    return (1 - ecc * ecc) / slr * ups_r / sqrt_J


def calc_lambda_0(chi, zp, zm, En, Lz, aa, slr, x):
    """
    Mino time as a function of polar angle, chi

    Parameters:
        chi (float): polar angle
        zp (float): polar root
        zm (float): polar root
        En (float): energy
        Lz (float): angular momentum
        aa (float): spin
        slr (float): semi-latus rectum
        x (float): inclination

    Returns:
        lambda_0 (float)

    """
    beta = aa * aa * (1 - En * En)
    k = sqrt(zm / zp)
    k2 = k * k
    prefactor = 1 / sqrt(beta * zp)
    ellipticK_k = ellipk(k2)
    ellipticF = ellipf(pi / 2 - chi, k2)

    return prefactor * (ellipticK_k - ellipticF)


def calc_wtheta(chi, ups_theta, En, Lz, aa, slr, zp, zm, x):
    """
    w_theta = ups_theta * lambda as a function of polar angle chi

    Parameters:
        chi (float): polar angle
        ups_theta (float): theta mino frequency
        En (float): energy
        Lz (float): angular momentum
        aa (float): spin
        slr (float): semi-latus rectum
        x (float): inclination

    Returns:
        w_theta (float)
    """
    if (chi >= 0 and chi <= pi / 2):
        return ups_theta * calc_lambda_0(chi, zp, zm, En, Lz, aa, slr, x)
    elif (chi > pi / 2 and chi <= pi):
        return pi - ups_theta * calc_lambda_0(pi - chi, zp, zm, En, Lz, aa,
                                              slr, x)
    elif (chi > pi and chi <= 3 * pi / 2):
        return pi + ups_theta * calc_lambda_0(chi - pi, zp, zm, En, Lz, aa,
                                              slr, x)
    elif (chi > 3 * pi / 2 and chi <= 2 * pi):
        return (2 * pi - ups_theta *
                calc_lambda_0(2 * pi - chi, zp, zm, En, Lz, aa, slr, x))
    else:
        print("Something went wrong in calc_wtheta!")
        return 0.  # this case should not occur, but is required by C++


def calc_dwtheta_dchi(chi, zp, zm):
    """
    derivative of w_theta

    Parameters:
        chi (float): polar angle
        zp (float): polar root
        zm (float): polar root

    Returns:
        dw_dtheta (float)
    """
    k = sqrt(zm / zp)
    ellipticK_k = ellipk(k**2)
    return pi / (2 * ellipticK_k) * (1 / (1 - k * k * cos(chi)**2))

#------------------------------------------------------------------------------
# Needed!
#------------------------------------------------------------------------------

def calc_radius(psi, slr, ecc):
    """
    r coordinate in terms of radial angle, psi

    Parameters:
        psi (float): radial angle
        slr (float): semi-latus rectum
        ecc (float): eccentricity

    Returns:
        radius (float)
    """
    return slr / (1 + ecc * cos(psi))


def calc_theta(chi, zm):
    """
    theta coordinate in terms of polar angle chi

    Parameters:
        chi (float): polar angle
        zm (float): polar root

    Returns:
        theta (float)
    """
    return acos(zm * cos(chi))


# def calc_Vt_theta(theta, En, aa):
#     return -En * aa**2 * sin(theta)**2


# def calc_Vt_r(r, En, Lz, aa):
#     return ((En*(aa**2 + r**2)**2)/(aa**2 - 2*r + r**2) + 
#             aa*Lz*(1 - (aa**2 + r**2)/(aa**2 - 2*r + r**2)))


# def calc_Vphi_theta(theta, Lz):
#     return Lz * csc(theta)**2


# def calc_Vphi_r():
#     return (-((aa**2*Lz)/(aa**2 - 2*r + r**2)) + aa*En*(-1 + (aa**2 + r**2)/
#              (aa**2 - 2*r + r**2)))


# def calc_delta_t_theta(ups_theta, zp, zm, En, Lz, aa, slr, x):
#     k = 2
#     def integrand(chi):
#         print('chi = ', chi)
#         theta = calc_theta(chi, zm)
#         # print('theta = ', theta)
#         dwtheta_dchi = calc_dwtheta_dchi(chi, zp, zm)
#         print('dwtheta_dchi = ', dwtheta_dchi)
#         wtheta = calc_wtheta(chi, ups_theta, En, Lz, aa, slr, x)
#         print('wtheta = ', wtheta)
#         Vt_theta = calc_Vt_theta(theta, En, aa)
#         print('Vt_theta = ', Vt_theta)
#         return dwtheta_dchi * Vt_theta * cos(k * wtheta)
#     answer = 4 / (pi * k * ups_theta) * quadgl(integrand, [pi/2, pi])

#     print('answer = ', answer)


# def am(x, k):
#     # TODO (aaron): change this to a different algorithm for speed
#     def integrand(x):
#         dn = ellipfun('dn')
#         return dn(x, k)
#     return quadgl(integrand, [0, x])


def am(x, m):
    """
    Jacobi amplitude function for real input parameters

    (JacobiAmplitude[] in Mathematica)

    Parameters:
        x (float)
        m (float)

    Returns:
        am (float)
    """
    if m == 0:
        return x
    elif x == 0:
        return 0
    else:
        return jacobi_am_f(x, m)


def calc_rq(qr, r1, r2, r3, r4):
    """
    function used in computing radial geodesic coordinates

    Parameters:
        qr (float)
        r1 (float): radial root
        r2 (float): radial root
        r3 (float): radial root
        r4 (float): radial root

    Returns:
        rq (float)
    """
    kr = ((r1 - r2) * (r3 - r4)) / ((r1 - r3) * (r2 - r4))

    sn = ellipfun('sn')
    return ((-(r2*(r1 - r3)) +
            (r1 - r2)*r3*sn((qr*ellipk(kr))/pi,kr)**2)/
            (-r1 + r3 + (r1 - r2)*sn((qr*ellipk(kr))/pi,kr)**2))


def calc_zq(qz, zp, zm, En, aa):
    """
    function used in computing polar geodesic coordinates

    Parameters:
        qz (float)
        zp (float): polar root
        zm (float): polar root
        En (float): energy
        aa (float): spin

    Returns:
        zq (float)
    """
    ktheta = (aa**2*(1 - En**2)*zm**2)/zp**2
    sn = ellipfun('sn')

    return zm*sn((2*(pi/2. + qz)*ellipk(ktheta))/pi, ktheta)


def calc_psi_r(qr, r1, r2, r3, r4):
    """
    radial geodesic angle

    Parameters:
        qr (float)
        r1 (float): radial root
        r2 (float): radial root
        r3 (float): radial root
        r4 (float): radial root

    Returns:
        psi_r (float)
    """
    # print(qr)
    kr = ((r1 - r2)*(r3 - r4))/((r1 - r3)*(r2 - r4))
    # print('ellipk(kr) = ', ellipk(kr))
    return am((qr*ellipk(kr))/pi, kr)


# for circular orbits, the following are zero:

def calc_t_r(qr, r1, r2, r3, r4, En, Lz, aa, M=1):
    """
    delta_t_r in Drasco and Hughes (2003?)

    Parameters:
        qr (float)
        r1 (float): radial root
        r2 (float): radial root
        r3 (float): radial root
        r4 (float): radial root
        En (float): energy
        Lz (float): angular momentum
        aa (float): spin

    Keyword Args:
        M (float): mass

    Returns:
        t_r (float)
    """
    psi_r = calc_psi_r(qr, r1, r2, r3, r4)

    kr = ((r1 - r2)*(r3 - r4))/((r1 - r3)*(r2 - r4))

    rp = M + sqrt(-aa**2 + M**2)
    rm = M - sqrt(-aa**2 + M**2)

    hr = (r1 - r2)/(r1 - r3)
    hp = ((r1 - r2)*(r3 - rp))/((r1 - r3)*(r2 - rp))
    hm = ((r1 - r2)*(r3 - rm))/((r1 - r3)*(r2 - rm))

    return  (-((En * ((-4*(r2 - r3)*(-(((-2*aa**2 + (4 - (aa*Lz)/En)*rm)*
            ((qr*ellippi(hm,kr))/pi - ellippi(hm,psi_r,kr)))/
            ((r2 - rm)*(r3 - rm))) + 
            ((-2*aa**2 + (4 - (aa*Lz)/En)*rp)*
            ((qr*ellippi(hp,kr))/pi - ellippi(hp,psi_r,kr)))/
            ((r2 - rp)*(r3 - rp))))/(-rm + rp) + 
            4*(r2 - r3)*((qr*ellippi(hr,kr))/pi - ellippi(hr,psi_r,kr)) + 
            (r2 - r3)*(r1 + r2 + r3 + r4)*
            ((qr*ellippi(hr,kr))/pi - ellippi(hr,psi_r,kr)) + 
            (r1 - r3)*(r2 - r4)*((qr*ellipe(kr))/pi - ellipe(psi_r,kr) + 
            (hr*cos(psi_r)*sin(psi_r)*sqrt(1 - kr*sin(psi_r)**2))/
            (1 - hr*sin(psi_r)**2))))/sqrt((1 - En**2)*(r1 - r3)*(r2 - r4))))


def calc_phi_r(qr, r1, r2, r3, r4, En, Lz, aa, M=1):
    """
    delta_phi_r in Drasco and Hughes (2003?)

    Parameters:
        qr (float)
        r1 (float): radial root
        r2 (float): radial root
        r3 (float): radial root
        r4 (float): radial root
        En (float): energy
        Lz (float): angular momentum
        aa (float): spin

    Keyword Args:
        M (float): mass

    Returns:
        phi_r (float)
    """
    psi_r = calc_psi_r(qr, r1, r2, r3, r4)
    kr = ((r1 - r2)*(r3 - r4))/((r1 - r3)*(r2 - r4))
    rp = M + sqrt(-aa**2 + M**2)
    rm = M - sqrt(-aa**2 + M**2)

    hp = ((r1 - r2)*(r3 - rp))/((r1 - r3)*(r2 - rp))
    hm = ((r1 - r2)*(r3 - rm))/((r1 - r3)*(r2 - rm))

    return ((2*aa*En*(-(((r2 - r3)*(-((aa*Lz)/En) + 2*rm)*
            ((qr*ellippi(hm,kr))/pi - ellippi(hm,psi_r,kr)))/
            ((r2 - rm)*(r3 - rm))) + 
            ((r2 - r3)*(-((aa*Lz)/En) + 2*rp)*
            ((qr*ellippi(hp,kr))/pi - ellippi(hp,psi_r,kr)))/
            ((r2 - rp)*(r3 - rp))))/
            (sqrt((1 - En**2)*(r1 - r3)*(r2 - r4))*(-rm + rp)))


def calc_psi_z(qz, zp, zm, En, aa):
    """
    angle used in polar geodesic calculations

    Parameters:
        qz (float)
        zp (float): polar root
        zm (float): polar root
        En (float): energy
        aa (float): spin

    Returns:
        psi_z (float)
    """
    ktheta = (aa**2*(1 - En**2)*zm**2)/zp**2
    return am((2*(pi/2. + qz)*ellipk(ktheta))/pi, ktheta)

# for equatorial orbits, the following are zero:

def calc_t_z(qz, zp, zm, En, aa):
    """
    delta_t_theta in Drasco and Hughes (2003?)

    Parameters:
        qz (float)
        zp (float): polar root
        zm (float): polar root
        En (float): energy
        aa (float): spin

    Returns:
        t_z (float)
    """
    psi_z = calc_psi_z(qz, zp, zm, En, aa)
    ktheta = (aa**2*(1 - En**2)*zm**2)/zp**2
    return ((En*zp*((2*(pi/2. + qz)*ellipe(ktheta))/pi - ellipe(psi_z,ktheta)))/
            (1 - En**2))


def calc_phi_z(qz, zp, zm, En, Lz, aa):
    """
    delta_phi_theta in Drasco and Hughes (2003?)

    Parameters:
        qz (float)
        zp (float): polar root
        zm (float): polar root
        En (float): energy
        Lz (float): angular momentum
        aa (float): spin

    Returns:
        phi_z (float)
    """
    psi_z = calc_psi_z(qz, zp, zm, En, aa)
    ktheta = (aa**2*(1 - En**2)*zm**2)/zp**2
    return (-((Lz*((2*(pi/2. + qz)*ellippi(zm**2,ktheta))/pi - 
            ellippi(zm**2,psi_z,ktheta)))/zp))


def calc_Ct(qr0, qz0, r1, r2, r3, r4, zp, zm, En, Lz, aa):
    """
    phase constant so that Mino time starts at 0

    Parameters:
        qr0 (float): initial radial phase
        qz0 (float): initial polar phase
        r1 (float): radial root
        r2 (float): radial root
        r3 (float): radial root
        r4 (float): radial root
        zp (float): polar root
        zm (float): polar root
        En (float): energy
        Lz (float): angular momentum
        aa (float): spin

    Returns:
        Ct (float)
    """
    t_r = calc_t_r(qr0, r1, r2, r3, r4, En, Lz, aa)
    t_z = calc_t_z(qz0, zp, zm, En, aa)
    return t_r + t_z


def calc_Cz(qr0, qz0, r1, r2, r3, r4, zp, zm, En, Lz, aa):
    """
    phase constant so that Mino time starts at 0

    Parameters:
        qr0 (float): initial radial phase
        qz0 (float): initial polar phase
        r1 (float): radial root
        r2 (float): radial root
        r3 (float): radial root
        r4 (float): radial root
        zp (float): polar root
        zm (float): polar root
        En (float): energy
        Lz (float): angular momentum
        aa (float): spin

    Returns:
        Cz (float)
    """
    phi_r = calc_phi_r(qr0, r1, r2, r3, r4, En, Lz, aa)
    phi_z = calc_phi_z(qz0, zp, zm, En, Lz, aa)
    return phi_r + phi_z


def calc_t(mino_t, ups_r, ups_theta, gamma, qt0, qr0, qz0, r1, r2, r3, r4, zp,
           zm, En, Lz, aa):
    """
    time geodesic coordinate

    Parameters:
        mino_t (float): Mino time
        ups_r (float): radial Mino frequency
        ups_theta (float): theta Mino frequency
        gamma (float): time Mino frequency
        qt0 (float): initial time phase
        qr0 (float): initial radial phase
        qz0 (float): initial theta phase
        r1 (float): radial root
        r2 (float): radial root
        r3 (float): radial root
        r4 (float): radial root
        zp (float): polar root
        zm (float): polar root
        En (float): energy
        Lz (float): angular momentum
        aa (float): spin

    Returns:
        t (float)
    """
    eta_t = qt0 + gamma * mino_t
    eta_r = qr0 + ups_r * mino_t
    eta_z = qz0 + ups_theta * mino_t
    if r1 == r2:
        t_r = 0
    else:
        t_r = calc_t_r(eta_r, r1, r2, r3, r4, En, Lz, aa)
    if zm == 0:
        t_z = 0
    else:
        t_z = calc_t_z(eta_z, zp, zm, En, aa)
    # print('t_z = ', t_z)
    if qr0 == 0 and qz0 == 0:
        Ct = 0
    else:
        Ct = calc_Ct(qr0, qz0, r1, r2, r3, r4, zp, zm, En, Lz, aa)
    return eta_t + t_r + t_z - Ct


def calc_r(mino_t, ups_r, qr0, r1, r2, r3, r4):
    """
    radius in terms of Mino time

    Parameters:
        mino_t (float): Mino time
        ups_r (float): Mino radial frequency
        qr0 (float): inital radial phase
        r1 (float): radial root
        r2 (float): radial root
        r3 (float): radial root
        r4 (float): radial root

    Returns:
        r (float)
    """
    eta = ups_r * mino_t + qr0
    return calc_rq(eta, r1, r2, r3, r4)


# def calc_theta(mino_t, ups_theta, qz0, zp, zm, En, aa):
#     """
#     theta in terms of Mino time

#     Parameters:
#         mino_t (float): Mino time
#         ups_theta (float): Mino theta frequency
#         qz0 (float): inital polar phase
#         zp (float): polar root
#         zm (float): polar root
#         En (float): energy
#         aa (float): spin

#     Returns:
#         theta (float)
#     """
#     eta = ups_theta * mino_t + qz0
#     return acos(calc_zq(eta, zp, zm, En, aa))


def calc_phi(mino_t, ups_r, ups_theta, ups_phi, qphi0, qr0, qz0, r1, r2, r3,
             r4, zp, zm, En, Lz, aa):
    """
    phi in terms of Mino time

    Parameters:
        mino_t (float): Mino time
        ups_r (float): Mino radial frequency
        ups_theta (float): Mino theta frequency
        ups_phi (float): Mino phi frequency
        qphi0 (float): initial phi phase
        qr0 (float): initial radial phase
        qz0 (float): initial polar phase
        r1 (float): radial root
        r2 (float): radial root
        r3 (float): radial root
        r4 (float): radial root
        zp (float): polar root
        zm (float): polar root
        En (float): energy
        Lz (float): angular momentum
        aa (float): spin

    Returns:
        phi (float)
    """
    eta_phi = ups_phi * mino_t + qphi0
    eta_r = ups_r * mino_t + qr0
    eta_theta = ups_theta * mino_t + qz0
    if r1 == r2:
        phi_r = 0
    else:
        phi_r = calc_phi_r(eta_r, r1, r2, r3, r4, En, Lz, aa)
    # print('phi_r = ', phi_r)
    if zm == 0:
        phi_z = 0
    else:
        phi_z = calc_phi_z(eta_theta, zp, zm, En, Lz, aa)
    if qr0 == 0 and qz0 == 0:
        Cz = 0
    else:
        Cz = calc_Cz(qr0, qz0, r1, r2, r3, r4, zp, zm, En, Lz, aa)
    return eta_phi + phi_r + phi_z - Cz


def calc_lambda_r(r, r1, r2, r3, r4, En):
    """
    Mino time as a function of r (which in turn is a function of psi)

    Parameters:
        r (float): radius
        r1 (float): radial root
        r2 (float): radial root
        r3 (float): radial root
        r4 (float): radial root
        En (float): energy

    Returns:
        lambda (float)
    """
    kr = ((r1 - r2)*(r3 - r4))/((r1 - r3)*(r2 - r4))
    if r1 == r2:
        # circular orbit
        print('Circular orbits currently do not work.')
        return 0
    else:
        if r1 < r2:
            r1 = r2
        yr = sqrt(((r - r2) * (r1 - r3)) / ((r1 - r2) * (r - r3)))
    F_asin = ellipf(asin(yr), kr)
    return (2*F_asin)/(sqrt(1 - En * En)*sqrt((r1 - r3)*(r2 - r4)))


def calc_lambda_psi(psi, ups_r, r1, r2, r3, r4, En, slr, ecc):
    """
    changes lambda(r) -> lambda(psi) by computing lambda(r(psi))

    Parameters:
        psi (float): radial angle
        ups_r (float): radial Mino frequency
        r1 (float): radial root
        r2 (float): radial root
        r3 (float): radial root
        r4 (float): radial root
        En (float): energy
        slr (float): semi-latus rectum
        ecc (float): eccentricity

    Returns:
        lambda_psi (float)
    """
    # print('psi = ', psi)
    r = calc_radius(psi, slr, ecc)
    lam_r = 2 * pi / ups_r  # radial period
    # print('lam_r = ', lam_r)
    lam_r1 = calc_lambda_r(r2, r1, r2, r3, r4, En)
    # print('lam_r1 = ', lam_r1)
    turns = floor(psi / (2 * pi))
    # print('turns = ', turns)
    if (psi % (2 * pi)) <= pi:
        # print('here')
        # print('r = ', r)
        res = calc_lambda_r(r, r1, r2, r3, r4, En) - lam_r1
        # print(res)
    else:
        res = lam_r1 - calc_lambda_r(r, r1, r2, r3, r4, En)
    # print('res = ', res)

    return r, lam_r * turns + res


# def mp_calc_circular_eq_coords_old(psi, En, Lz, aa, slr):
#     """
#     THIS FUNCTION DOESN'T WORK PROPERLY
#     Computes circular equatorial coords in a convenient function

#     Parameters:
#         psi (float): radial angle
#         En (float): energy
#         Lz (float): z direction of angular momentum
#         aa (float): spin
#         slr (float): semi-latus rectum

#     Returns:
#         t (float): time coordinate
#         r (float): radial coordinate
#         theta (float): polar coordinate
#         phi (float): azimuthal coordinate
#     """
#     lam_psi = (psi / (sqrt(1 - En**2 + 3*(1 - En**2) + 2*(1 - En**2 - 1/slr) + 
#                (aa**2*(1 - En**2) + Lz**2)/slr**2 - 4/slr)*slr))
#     # print(lam_psi)
#     t = ((lam_psi*(aa**3*sqrt(2*aa + (-3 + slr)*sqrt(slr))*slr**2 + 
#            aa*sqrt(2*aa + (-3 + slr)*sqrt(slr))*(-2 + slr)*slr**3 + 
#            aa**2*sqrt((2*aa + (-3 + slr)*sqrt(slr))*slr**7) - 
#            2*sqrt((2*aa + (-3 + slr)*sqrt(slr))*slr**9) + 
#            sqrt((2*aa + (-3 + slr)*sqrt(slr))*slr**11))) /
#          ((2*aa + (-3 + slr)*sqrt(slr))*slr**0.75*(aa**2 + (-2 + slr)*slr)))
#     theta = pi / 2
#     phi = (lam_psi*slr**1.25)/sqrt(2*aa + (-3 + slr)*sqrt(slr))
#     return float(t), float(slr), float(theta), float(phi)


def mp_calc_circular_eq_coords(psi, En, Lz, aa, slr, M=1):
    """
    Computes circular equatorial coords in a convenient function

    This function is based on the phi and t functions in Glampedakis
    and Kennfick (2002).

    Parameters:
        psi (float): radial angle
        En (float): energy
        Lz (float): z direction of angular momentum
        aa (float): spin
        slr (float): semi-latus rectum

    Returns:
        t (float): time coordinate
        r (float): radial coordinate
        theta (float): polar coordinate
        phi (float): azimuthal coordinate
    """
    # lam_psi is used to cross check with circular orbits in BHPTK
    # lam_psi = (psi / (sqrt(1 - En**2 + 3*(1 - En**2) + 2*(1 - En**2 - 1/slr) + 
    #            (aa**2*(1 - En**2) + Lz**2)/slr**2 - 4/slr)*slr))
    # print(lam_psi)

    x = Lz - aa * En
    Vr = aa**2 + x**2 + 2 * aa * x * En - 6 * M * x**2 / slr
    Vphi = x + aa * En - 2 * M * x / slr
    Vt = aa**2 * En - 2 * aa * M * x / slr + En * slr**2
    J = 1 - 2 * M / slr + aa**2 / slr**2
    denom = J * Vr**(1/2)
    t = Vt / denom * psi
    phi = Vphi / denom * psi
    return float(t), float(slr), float(pi / 2), float(phi)



def mp_calc_equatorial_coords(psi, ups_r, ups_theta, ups_phi, gamma, r1, r2, r3,
                              r4, zp, zm, En, Lz, aa, slr, ecc, qt0=0,
                              qr0=0, qz0=0, qphi0=0):
    """
    Computes all equatorial coordinates in a convenient function

    Parameters:
        psi (float): radial angle
        ups_r (float): radial Mino frequency
        ups_theta (float): theta Mino frequency
        ups_phi (float): phi Mino frequency
        gamma (float): time Mino frequency
        r1 (float): radial root
        r2 (float): radial root
        r3 (float): radial root
        r4 (float): radial root
        zp (float): polar root
        zm (float): polar root
        En (float): energy
        Lz (float): angular momentum
        aa (float): spin
        slr (float): semi-latus rectum
        ecc (float): eccentricity

    Keyword Args:
        qt0 (float): initial time phase
        qr0 (float): initial radial phase
        qz0 (float): initial theta phase
        qphi0 (float): initial phi phase

    Returns:
        t (float): time coordinate
        r (float): radial coordinate
        theta (float): polar coordinate
        phi (float): azimuthal coordinate
    """
    if zm != 0:
        print('The orbit specified is not equatorial.')
    r, lam_psi = calc_lambda_psi(psi, ups_r, r1, r2, r3, r4, En, slr, ecc)
    t = calc_t(lam_psi, ups_r, ups_theta, gamma, qt0, qr0, qz0, r1, r2, r3, r4,
               zp, zm, En, Lz, aa)
    theta = pi / 2
    phi = calc_phi(lam_psi, ups_r, ups_theta, ups_phi, qphi0, qr0, qz0, r1, r2,
                   r3, r4, zp, zm, En, Lz, aa)
    return float(t), float(r), float(theta), float(phi)


# if __name__ == '__main__':
#     aa = 0.998
#     slr = 3
#     ecc = 0.1
#     x = 1

#     En, Lz, Q = mp_calc_constants(aa, slr, ecc, x)
#     # print(En, Lz, Q)
#     r1, r2, r3, r4 = mp_radial_roots(En, Q, aa, slr, ecc)
#     # print(r1, r2, r3, r4)
#     zp, zm = mp_polar_roots(En, Lz, aa, slr, x)
#     # print(zp, zm)
#     ups_r, ups_theta, ups_phi, gamma = mp_mino_freqs(r1, r2, r3, r4, En, Lz,
#                                                      Q, aa, slr, ecc, x)
#     # print(ups_r, ups_theta, ups_phi, gamma)
#     psi = 3.082771039571993143494951148712902333298005500564555505306750660854449711396384940176778119898643597005105
#     t, r, theta, phi = mp_calc_equatorial_coords(psi, ups_r, ups_theta,
#                    ups_phi, gamma, r1, r2, r3, r4, zp, zm, En, Lz, aa, slr, ecc)
#     print(t)
#     print(r)
#     print(theta)
#     print(phi)





