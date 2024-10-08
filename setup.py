from setuptools import setup, Extension
from Cython.Build import cythonize  

import numpy

extensions = [
    Extension("tdavec.TDAvec", ["tdavec/TDAvec.pyx"],
              include_dirs=[numpy.get_include()],  # Include NumPy's header files
              )
]

setup(
    name = "tdavec",
    version = "0.1.0",
    author = "Umar Islambekov, Aleksei Luchinsky",
    packages = ["tdavec"],
    ext_modules = cythonize(extensions),
    setup_requires=["Cython", "numpy", "ripser"],
    zip_safe = False
)