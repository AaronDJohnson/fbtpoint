from distutils.core import setup, Extension
from Cython.Build import cythonize

ext = Extension("homteuk",
    sources=["homteuk.pyx", "homteukc.cpp", "fminsol.cpp", "hypgeo.cpp", "amplitudes.cpp", "sp_func.cpp"],
    library_dirs=['/usr/local/include/gsl', '/usr/local/opt/openblas'],
    libraries=["gsl", "arb", "cblas"],
    include_dirs=['/usr/local/include', '/usr/local/opt'],
    language="c++", extra_compile_args=["-std=c++17"])
setup(name="find_R", ext_modules=cythonize([ext]))