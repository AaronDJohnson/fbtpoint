# SWSH

* This package computes the spin-weighted spheroidal harmonics (SWSHs) following the method of the [Black Hole Perturbation Toolkit][https://bhptoolkit.org/].
* The SWSHs are computed using a spectral decomposition in terms of the spin-weighted spherical harmonics.
* This package also computes the first and second theta derivatives as needed for frequency based teukolsky codes.
* At a later date, Leaver's method may be implemented at high precision. This code exists currently, but is not as robust, and needs to be tested for higher (l, m) values.