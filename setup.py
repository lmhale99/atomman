from setuptools import setup, find_packages

def readme():
    with open('README.rst') as f:
        return f.read()
    

setup(name='atomman',
      version='0.6',
      description='Atomistic Manipulation Toolkit',
      long_description=readme(),
      classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Natural Language :: English',
        'Programming Language :: Python :: 2.7',
        'Topic :: Scientific/Engineering :: Physics'
      ],
      keywords=['atom', 'atomic', 'atomistic', 'molecular dynamics'], 
      url='https://github.com/usnistgov/atomman',
      author='Lucas Hale',
      author_email='lucas.hale@nist.gov',
      packages=find_packages(),
      install_requires=[
        'xmltodict',
        'DataModelDict',
        'numpy', 
        'matplotlib',
        'scipy',
        'numericalunits'
      ],
      include_package_data=True,
      zip_safe=False)