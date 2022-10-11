# coding: utf-8

# https://docs.pytest.org/en/latest/
import pytest

import atomman as am

class Test_crystalsystem(object):

    @property
    def boxcubic(self):
        return am.Box(a=3, b=3, c=3, alpha=90, beta=90, gamma=90)
    
    @property
    def boxhexagonal(self):
        return am.Box(a=3, b=3, c=8, alpha=90, beta=90, gamma=120)

    @property
    def boxtetragonal(self):
        return am.Box(a=3, b=3, c=8, alpha=90, beta=90, gamma=90)
    
    @property
    def boxrhombohedral(self):
        return am.Box(a=3, b=3, c=3, alpha=84, beta=84, gamma=84)
    
    @property
    def boxorthorhombic(self):
        return am.Box(a=3, b=4, c=5, alpha=90, beta=90, gamma=90)
    
    @property
    def boxmonoclinic(self):
        return am.Box(a=3, b=4, c=5, alpha=90, beta=105, gamma=90)

    @property
    def boxtriclinic(self):
        return am.Box(a=3, b=4, c=5, alpha=74, beta=87, gamma=105)
    
    def test_iscubic(self):
        assert self.boxcubic.iscubic()
        assert not self.boxhexagonal.iscubic()
        assert not self.boxtetragonal.iscubic()
        assert not self.boxrhombohedral.iscubic()
        assert not self.boxorthorhombic.iscubic()
        assert not self.boxmonoclinic.iscubic()
        assert not self.boxtriclinic.iscubic()

    def test_ishexagonal(self):
        assert not self.boxcubic.ishexagonal()
        assert self.boxhexagonal.ishexagonal()
        assert not self.boxtetragonal.ishexagonal()
        assert not self.boxrhombohedral.ishexagonal()
        assert not self.boxorthorhombic.ishexagonal()
        assert not self.boxmonoclinic.ishexagonal()
        assert not self.boxtriclinic.ishexagonal()

    def test_istetragonal(self):
        assert not self.boxcubic.istetragonal()
        assert not self.boxhexagonal.istetragonal()
        assert self.boxtetragonal.istetragonal()
        assert not self.boxrhombohedral.istetragonal()
        assert not self.boxorthorhombic.istetragonal()
        assert not self.boxmonoclinic.istetragonal()
        assert not self.boxtriclinic.istetragonal()

    def test_isrhombohedral(self):
        assert not self.boxcubic.isrhombohedral()
        assert not self.boxhexagonal.isrhombohedral()
        assert not self.boxtetragonal.isrhombohedral()
        assert self.boxrhombohedral.isrhombohedral()
        assert not self.boxorthorhombic.isrhombohedral()
        assert not self.boxmonoclinic.isrhombohedral()
        assert not self.boxtriclinic.isrhombohedral()

    def test_isorthorhombic(self):
        assert not self.boxcubic.isorthorhombic()
        assert not self.boxhexagonal.isorthorhombic()
        assert not self.boxtetragonal.isorthorhombic()
        assert not self.boxrhombohedral.isorthorhombic()
        assert self.boxorthorhombic.isorthorhombic()
        assert not self.boxmonoclinic.isorthorhombic()
        assert not self.boxtriclinic.isorthorhombic()

    def test_ismonoclinic(self):
        assert not self.boxcubic.ismonoclinic()
        assert not self.boxhexagonal.ismonoclinic()
        assert not self.boxtetragonal.ismonoclinic()
        assert not self.boxrhombohedral.ismonoclinic()
        assert not self.boxorthorhombic.ismonoclinic()
        assert self.boxmonoclinic.ismonoclinic()
        assert not self.boxtriclinic.ismonoclinic()

    def test_istriclinic(self):
        assert not self.boxcubic.istriclinic()
        assert not self.boxhexagonal.istriclinic()
        assert not self.boxtetragonal.istriclinic()
        assert not self.boxrhombohedral.istriclinic()
        assert not self.boxorthorhombic.istriclinic()
        assert not self.boxmonoclinic.istriclinic()
        assert self.boxtriclinic.istriclinic()

    def test_identifyfamily(self):
        assert self.boxcubic.identifyfamily() == 'cubic'
        assert self.boxhexagonal.identifyfamily() == 'hexagonal'
        assert self.boxtetragonal.identifyfamily() == 'tetragonal'
        assert self.boxrhombohedral.identifyfamily() == 'rhombohedral'
        assert self.boxorthorhombic.identifyfamily() == 'orthorhombic'
        assert self.boxmonoclinic.identifyfamily() == 'monoclinic'
        assert self.boxtriclinic.identifyfamily() == 'triclinic'