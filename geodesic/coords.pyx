cdef extern from "coords.hpp":
    void calc_equatorial_coords(double &t, double &r, double &theta, double &phi,
            double psi, double ups_r, double ups_theta,
            double ups_phi, double gamma, double r1, double r2, double r3,
            double r4, double zp, double zm, double En, double Lz, double aa,
            double slr, double ecc, double qt0, double qr0, double qz0,
            double qphi0)
    void calc_circular_eq_coords(double &t, double &r, double &theta, double &phi,
				                 double psi, double En, double Lz, double aa,
                                 double slr, double M)
    void calc_gen_coords_mino(double &t, double &r, double &theta, double &phi,
        double mino_t, double ups_r, double ups_theta, double ups_phi, double gamma,
        double r1, double r2, double r3, double r4, double zp, double zm, double En,
        double Lz, double aa, double qphi0, double qr0, double qz0, double qt0)
    void calc_gen_coords_psi(double &t, double &r, double &theta, double &phi,
        double psi, double ups_r, double ups_theta, double ups_phi, double gamma,
        double r1, double r2, double r3, double r4, double zp, double zm, double En,
        double Lz, double Q, double aa, double slr, double ecc, double x,
        double qphi0, double qr0, double qz0, double qt0)
    double calc_wr(double psi, double ups_r, double En, double Lz,
                double Q, double aa, double slr, double ecc, double x)
    double calc_dwr_dpsi(double psi, double ups_r, double En, double Lz,
                        double Q, double aa, double slr, double ecc)
    double calc_wtheta(double chi, double ups_theta, double zp, double zm,
                    double En, double Lz, double aa, double slr, double x)
    double calc_dwtheta_dchi(double chi, double zp, double zm)


def py_calc_equatorial_coords(double psi, double ups_r, double ups_theta,
            double ups_phi, double gamma, double r1, double r2, double r3,
            double r4, double zp, double zm, double En, double Lz, double aa,
            double slr, double ecc, double M=1):
    """
    This assumes that all the initial phases are zero
    """

    cdef double t = 0.0
    cdef double r = 0.0
    cdef double theta = 0.0
    cdef double phi = 0.0
    if ecc == 0:
        calc_circular_eq_coords(t, r, theta, phi, psi, En, Lz, aa, slr, M)
    else:
        calc_equatorial_coords(t, r, theta, phi, psi, ups_r, ups_theta,
                                ups_phi, gamma, r1, r2, r3, r4, zp, zm, En, Lz, aa,
                                slr, ecc, 0, 0, 0, 0)
    return t, r, theta, phi


def py_calc_gen_coords_mino(double mino_t, double ups_r, double ups_theta, double ups_phi, double gamma,
        double r1, double r2, double r3, double r4, double zp, double zm, double En,
        double Lz, double aa, double qphi0=0, double qr0=0, double qz0=0, double qt0=0):
    cdef double t = 0.0
    cdef double r = 0.0
    cdef double theta = 0.0
    cdef double phi = 0.0
    calc_gen_coords_mino(t, r, theta, phi,
        mino_t, ups_r, ups_theta, ups_phi, gamma,
        r1, r2, r3, r4, zp, zm, En,
        Lz, aa, qphi0, qr0, qz0, qt0)
    return t, r, theta, phi


def py_calc_gen_coords_psi(double psi, double ups_r, double ups_theta, double ups_phi, double gamma,
        double r1, double r2, double r3, double r4, double zp, double zm, double En,
        double Lz, double Q, double aa, double slr, double ecc, double x,
        double qphi0=0, double qr0=0, double qz0=0, double qt0=0):

    cdef double t = 0.0
    cdef double r = 0.0
    cdef double theta = 0.0
    cdef double phi = 0.0

    calc_gen_coords_psi(t, r, theta, phi, psi, ups_r, ups_theta, ups_phi, gamma,
                        r1, r2, r3, r4, zp, zm, En, Lz, Q, aa, slr, ecc, x, qphi0, qr0, qz0, qt0)
    return t, r, theta, phi


def py_calc_wr(double psi, double ups_r, double En, double Lz,
               double Q, double aa, double slr, double ecc, double x):
    return calc_wr(psi, ups_r, En, Lz,
                   Q, aa, slr, ecc, x)


def py_calc_dwr_dpsi(double psi, double ups_r, double En, double Lz,
                        double Q, double aa, double slr, double ecc):

    return calc_dwr_dpsi(psi, ups_r, En, Lz, Q, aa, slr, ecc)


def py_calc_wtheta(double chi, double ups_theta, double zp, double zm,
                    double En, double Lz, double aa, double slr, double x):

    return calc_wtheta(chi, ups_theta, zp, zm, En, Lz, aa, slr, x)


def py_calc_dwtheta_dchi(double chi, double zp, double zm):

    return calc_dwtheta_dchi(chi, zp, zm)
