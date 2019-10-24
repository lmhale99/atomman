# coding: utf-8

# https://docs.pytest.org/en/latest/
import pytest

from atomman.tools import filltemplate

def test_filltemplate():
    template1 = 'this is my {val} to test'
    template2 = 'another <val> to test'
    template3 = 'yet another !@#$val?@ to test'

    var = {'val': 'template', 'notused':'STEVE!'}

    assert (filltemplate(template1, var, '{', '}')
            == 'this is my template to test')

    assert (filltemplate(template2, var, '<', '>')
            == 'another template to test')

    assert (filltemplate(template3, var, '!@#$', '?@')
            == 'yet another template to test')