try:
    from .coords import py_calc_equatorial_coords, py_calc_gen_coords_mino
    from .coords import py_calc_gen_coords_psi, py_calc_wr, py_calc_dwr_dpsi
    from .coords import py_calc_wtheta, py_calc_dwtheta_dchi
except:
    from coords import py_calc_equatorial_coords, py_calc_gen_coords_mino
    from coords import py_calc_gen_coords_psi, py_calc_wr, py_calc_dwr_dpsi
    from coords import py_calc_wtheta, py_calc_dwtheta_dchi


def calc_coordinates():
    