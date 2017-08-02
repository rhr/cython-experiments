from setuptools import setup
from setuptools.extension import Extension
from Cython.Build import cythonize

extensions = [
    Extension(
        'cyexpokit',
        ['cyexpokit.c'],
        libraries = ['lapack', 'blas', 'm'],
        extra_objects = ['dexpm_c.o', 'expokit.o', 'mataid.o', 'clock.o'],
        define_macros=[('CYTHON_TRACE', '1')]
        ),
    Extension(
        'odeiv',
        ['odeiv.c'],
        libraries = ['gsl', 'blas']
        ),
    Extension('tree', ['tree.c']),
    Extension('mcmc', ['mcmc.c'])
]

setup(
  name = 'cython-experiments',
  ext_modules = extensions
)
