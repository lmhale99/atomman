# coding: utf-8

# Standard Python libraries
from importlib import resources

import atomman as am

def test_version():
    if hasattr(resources, 'files'):
        assert resources.files('atomman').joinpath('VERSION').is_file()
    else:
        assert resources.is_resource('atomman', 'VERSION')
    print('Testing atomman version', am.__version__)
