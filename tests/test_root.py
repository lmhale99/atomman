# coding: utf-8

# Standard Python libraries
import os

# https://docs.pytest.org/en/latest/
import pytest

# http://www.numpy.org/
import numpy as np

import atomman as am

def test_rootdir():
    assert os.path.isdir(am.rootdir)

def test_version():
    assert os.path.isfile(os.path.join(os.path.join(am.rootdir, 'VERSION')))
    try:
        print('Testing atomman version', am.__version__)
    except:
        assert False