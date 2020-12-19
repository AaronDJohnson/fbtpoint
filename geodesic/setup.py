from distutils.core import setup, Extension
from Cython.Build import cythonize

ext = Extension("coords", sources=["coords.pyx", "coordsc.cpp"],
    library_dirs=['/usr/local/include/gsl', '/usr/local/opt/openblas'], libraries=["gsl","cblas"],
    include_dirs=['/usr/local/include', '/usr/local/opt'],
    language="c++", extra_compile_args=["-std=c++17"])
setup(name="eq_coords", ext_modules=cythonize([ext]))