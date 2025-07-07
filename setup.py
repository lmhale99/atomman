import os
from setuptools import setup, find_packages
from Cython.Build import cythonize

# Compile cython code files
extensions = cythonize([
    "atomman/core/dmag.pyx",
    "atomman/core/dvect.pyx",
    "atomman/core/nlist.pyx",
    "atomman/defect/slip_vector.pyx",
    "atomman/defect/Strain.pyx",
])

def getversion():
    """Fetches version information from VERSION file"""
    with open(os.path.join('atomman', 'VERSION'), encoding='UTF-8') as version_file:
        version = version_file.read().strip()
    return version

def getreadme():
    """Fetches description from README.rst file"""
    with open('README.rst', encoding='UTF-8') as readme_file:
        return readme_file.read()

setup(
    name = 'atomman',
    version = getversion(),
    description = 'Atomistic Manipulation Toolkit',
    long_description = getreadme(),
    ext_modules = extensions,
    classifiers = [
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'Natural Language :: English',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Topic :: Scientific/Engineering :: Physics'
    ],
    keywords = [
        'atom',
        'atomic',
        'atomistic',
        'molecular dynamics'
    ],
    url = 'https://github.com/usnistgov/atomman/',
    author = 'Lucas Hale',
    author_email = 'lucas.hale@nist.gov',
    packages = find_packages(),
    install_requires = [
        'xmltodict',
        'DataModelDict',
        'numericalunits',
        'numpy>=1.15',
        'matplotlib',
        'scipy',
        'pandas',
        'cython',
        'requests',
        'toolz',
        'potentials>=0.4.1',
        'yabadaba>=0.3.2'
    ],
    include_package_data = True,
    zip_safe = False
)
