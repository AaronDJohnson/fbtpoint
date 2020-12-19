# fbtpoint
 Frequency based, hypergeometric Teukolsky solving code for black hole perturbation theory

 ## Currently this does equatorial orbits, but it will be generalized

 ## This work makes use of the black hole perturbation toolkit

 * bhptoolkit.org

 This code is a work in progress. Currently, the geodesic and renormnu functions run well and have been tested on some data. I plan to add many more test functions with data from the literature.

 Geodesic and Teukolsky folders have Cython pieces in them. These pieces use functions from `gsl`, `cblas`, `boost` and `arb` which can be found and installed from 

 * http://arblib.org/setup.html

 via the `conda-forge` script which is highly recommended.

 Once the dependencies are installed properly, you can run `python setup.py build_ext --inplace` in a terminal to build the Cython pieces. Again, in the coming weeks, I hope to have the package installing itself easily.

 # Goals
 * Setup scripts to install package with `pip`
 * Generalize to non-equatorial orbits
 * Create a specific "highly-eccentric" mode that integrates more carefully