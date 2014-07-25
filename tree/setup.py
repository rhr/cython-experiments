from distutils.core import setup
from Cython.Build import cythonize

setup(
  name = 'tree',
  ext_modules = cythonize("tree.pyx"),
)
