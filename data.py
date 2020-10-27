import os
import numpy as np


def make_folder(slr, ecc, aa, verbose=False):
    f_name = ("data/slr=" + str(float(slr)) + "_ecc=" + str(ecc) + "_aa=" +
              str(aa))
    try:
        os.makedirs(f_name)
        print("Making folder for parameters.")
    except FileExistsError:
        if verbose:
            print("Folder for these parameters exists!")
    return f_name


def save_data(folder_path, ell, em, kay, energy):
    filename = "/l=" + str(ell) + "_m=" + str(em)
    abs_path = folder_path + filename
    if os.path.isfile(abs_path):
        with open(abs_path, "ab+") as f:
            np.savetxt(f, np.c_[int(ell), int(em), int(kay), energy])
    else:
        with open(abs_path, "ab+") as f:
            np.savetxt(f, np.c_[int(ell), int(em), int(kay), energy],
                       header="l m k E_inf")


def save_kmax_data(folder_path, ell, em, kay_guess, kay_max):
    filename = "/kmax_data"
    abs_path = folder_path + filename
    with open(abs_path, "ab+") as f:
        np.savetxt(f, np.c_[int(ell), int(em), int(kay_guess), int(kay_max)],
                   fmt='%5d')


def save_wave_data(folder_path, ell, em, kay, omega, Z, S):
    re_Z = np.real(Z)
    im_Z = np.imag(Z)
    re_S = np.real(S)
    im_S = np.imag(S)
    filename = "/l=" + str(ell) + "_m=" + str(em) + "_wave"
    abs_path = folder_path + filename
    if os.path.isfile(abs_path):
        with open(abs_path, "ab+") as f:
            np.savetxt(f, np.c_[int(ell), int(em), int(kay),
                                re_Z, im_Z, re_S, im_S, omega])
    else:
        with open(abs_path, "ab+") as f:
            np.savetxt(f, np.c_[int(ell), int(em), int(kay),
                                re_Z, im_Z, re_S, im_S, omega],
                       header="l m k real(Z) imag(Z) real(S) imag(S)")
