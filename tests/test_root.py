# coding: utf-8

# Standard Python libraries
from importlib import resources

import atomman as am

def test_version():
    assert resources.is_resource('atomman', 'VERSION')
    print('Testing atomman version', am.__version__)
