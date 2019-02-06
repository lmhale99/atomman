import os
import pytest

import atomman as am

def test_rootdir():
    assert os.path.isdir(am.rootdir)
    testrootdir = os.path.abspath(os.path.join(__file__, '..', '..', 'atomman'))
    assert os.path.samefile(am.rootdir, testrootdir)

def test_version():
    assert os.path.isfile(os.path.join(os.path.join(am.rootdir, 'VERSION')))
    try:
        print('Testing atomman version', am.__version__)
    except:
        assert False