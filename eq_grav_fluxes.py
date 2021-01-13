import numpy as np
from orbit.orbits import Orbit


low_ecc = Orbit(aa=0.9, slr=6, ecc=0.7, x=1)
low_ecc.set_modes(ell=2, em=2, en=50)
low_ecc.energy_inf()








