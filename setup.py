from setuptools import setup

def readme():
    with open('README.rst') as f:
        return f.read()
    

setup(name='atomman',
      version='0.5',
      description='Atomistic Manipulation Toolkit',
      long_description=readme(),
      classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Natural Language :: English',
        'Programming Language :: Python :: 2.7',
        'Topic :: Scientific/Engineering :: Physics'
      ],
      keywords='atom atomic atomistic molecular dynamics', 
      url='git@github.com:lmhale99/atomman.git',
      author='Lucas Hale',
      author_email='lucas.hale@nist.gov',
      license='NIST',
      packages=['atomman'],
      install_requires=[
        'mendeleev==0.1.0',
        'numpy', 
        'matplotlib',
        'scipy',
        'numericalunits',
        'xmltodict'
      ],
      include_package_data=True,
      zip_safe=False)