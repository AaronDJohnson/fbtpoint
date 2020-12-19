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



def py_calc_equatorial_coords(double psi, double ups_r, double ups_theta,
            double ups_phi, double gamma, double r1, double r2, double r3,
            double r4, double zp, double zm, double En, double Lz, double aa,
            double slr, double ecc, double M=1):

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
    
    
