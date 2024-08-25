from setuptools import setup
from Cython.Build import cythonize
import numpy

setup(
    ext_modules=cythonize([
         './tdavec/TDAvec.pyx'
         ],  # Python code file with primes() function
        annotate=True),                 # enables generation of the html annotation file
        include_dirs=[numpy.get_include()]
)
