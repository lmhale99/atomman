import os
from setuptools import setup, find_packages

def getversion():
    with open(os.path.join('atomman', 'VERSION')) as version_file:
        version = version_file.read().strip()
    return version

def getreadme():
    with open('README.rst') as f:
        return f.read()
    
setup(name = 'atomman',
      version = getversion(),
      description = 'Atomistic Manipulation Toolkit',
      long_description = getreadme(),
      classifiers = [
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Natural Language :: English',
        'Programming Language :: Python :: 2.7',
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
        'numpy', 
        'matplotlib',
        'scipy',
        'pandas',
        'numericalunits'
      ],
      include_package_data = True,
      zip_safe = False)