import pytest
import atomman as am

def test_filltemplate():
    template1 = 'this is my {val} to test'
    template2 = 'another <val> to test'
    template3 = 'yet another !@#$val?@ to test'

    var = {'val': 'template', 'notused':'STEVE!'}

    assert (am.tools.filltemplate(template1, var, '{', '}')
            == 'this is my template to test')

    assert (am.tools.filltemplate(template2, var, '<', '>')
            == 'another template to test')

    assert (am.tools.filltemplate(template3, var, '!@#$', '?@')
            == 'yet another template to test')