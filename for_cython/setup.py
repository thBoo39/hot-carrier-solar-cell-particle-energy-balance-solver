from distutils.core import setup
from Cython.Build import cythonize

setup(name='myProject',
      ext_modules=cythonize("*.pyx"))
