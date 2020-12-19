cdef extern from "homteuk.h":
    void find_R(double &re_R_in, double &im_R_in, double &re_Rd_in,
                double &im_Rd_in, double &re_Rdd_in, double &im_Rdd_in,
                double r, double re_nu, double im_nu, double eigen, double aa,
                double omega, int em)

    void find_Bin(double &re_B_in, double &im_B_in, double re_nu, double im_nu,
                  double eigen, double aa, double omega,
                  int em)


def py_find_R(double r, double re_nu, double im_nu, double aa, double omega,
    int em, double eigen, int ess=-2, double M=1.):

    cdef double re_R_in = 0.0
    cdef double im_R_in = 0.0
    cdef double re_Rd_in = 0.0
    cdef double im_Rd_in = 0.0
    cdef double re_Rdd_in = 0.0
    cdef double im_Rdd_in = 0.0

    find_R(re_R_in, im_R_in, re_Rd_in, im_Rd_in, re_Rdd_in, im_Rdd_in,
           r, re_nu, im_nu, eigen, aa, omega, em)
    # print(re_R_in)
    # print(im_R_in)
    R_in = re_R_in + 1j * im_R_in
    Rd_in = re_Rd_in + 1j * im_Rd_in
    Rdd_in = re_Rdd_in + 1j * im_Rdd_in
    return R_in, Rd_in, Rdd_in


def py_find_Bin(double re_nu, double im_nu, double eigen, double aa,
    double omega, int em):

    cdef double re_B_in = 0.0
    cdef double im_B_in = 0.0

    find_Bin(re_B_in, im_B_in, re_nu, im_nu, eigen, aa, omega, em)

    B_in = re_B_in + 1j * im_B_in
    return B_in
