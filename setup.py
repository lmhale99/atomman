from setuptools import setup
from Cython.Build import cythonize

# Compile cython code files
extensions = cythonize([
    "atomman/core/dmag.pyx",
    "atomman/core/dvect.pyx",
    "atomman/core/nlist.pyx",
    "atomman/defect/slip_vector.pyx",
    "atomman/defect/Strain.pyx",
])

setup(
    name = 'atomman',
    ext_modules = extensions
)
