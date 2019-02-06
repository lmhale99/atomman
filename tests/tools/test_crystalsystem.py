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
        assert am.tools.iscubic(self.boxcubic)
        assert not am.tools.iscubic(self.boxhexagonal)
        assert not am.tools.iscubic(self.boxtetragonal)
        assert not am.tools.iscubic(self.boxrhombohedral)
        assert not am.tools.iscubic(self.boxorthorhombic)
        assert not am.tools.iscubic(self.boxmonoclinic)
        assert not am.tools.iscubic(self.boxtriclinic)

    def test_ishexagonal(self):
        assert not am.tools.ishexagonal(self.boxcubic)
        assert am.tools.ishexagonal(self.boxhexagonal)
        assert not am.tools.ishexagonal(self.boxtetragonal)
        assert not am.tools.ishexagonal(self.boxrhombohedral)
        assert not am.tools.ishexagonal(self.boxorthorhombic)
        assert not am.tools.ishexagonal(self.boxmonoclinic)
        assert not am.tools.ishexagonal(self.boxtriclinic)

    def test_istetragonal(self):
        assert not am.tools.istetragonal(self.boxcubic)
        assert not am.tools.istetragonal(self.boxhexagonal)
        assert am.tools.istetragonal(self.boxtetragonal)
        assert not am.tools.istetragonal(self.boxrhombohedral)
        assert not am.tools.istetragonal(self.boxorthorhombic)
        assert not am.tools.istetragonal(self.boxmonoclinic)
        assert not am.tools.istetragonal(self.boxtriclinic)

    def test_isrhombohedral(self):
        assert not am.tools.isrhombohedral(self.boxcubic)
        assert not am.tools.isrhombohedral(self.boxhexagonal)
        assert not am.tools.isrhombohedral(self.boxtetragonal)
        assert am.tools.isrhombohedral(self.boxrhombohedral)
        assert not am.tools.isrhombohedral(self.boxorthorhombic)
        assert not am.tools.isrhombohedral(self.boxmonoclinic)
        assert not am.tools.isrhombohedral(self.boxtriclinic)

    def test_isorthorhombic(self):
        assert not am.tools.isorthorhombic(self.boxcubic)
        assert not am.tools.isorthorhombic(self.boxhexagonal)
        assert not am.tools.isorthorhombic(self.boxtetragonal)
        assert not am.tools.isorthorhombic(self.boxrhombohedral)
        assert am.tools.isorthorhombic(self.boxorthorhombic)
        assert not am.tools.isorthorhombic(self.boxmonoclinic)
        assert not am.tools.isorthorhombic(self.boxtriclinic)

    def test_ismonoclinic(self):
        assert not am.tools.ismonoclinic(self.boxcubic)
        assert not am.tools.ismonoclinic(self.boxhexagonal)
        assert not am.tools.ismonoclinic(self.boxtetragonal)
        assert not am.tools.ismonoclinic(self.boxrhombohedral)
        assert not am.tools.ismonoclinic(self.boxorthorhombic)
        assert am.tools.ismonoclinic(self.boxmonoclinic)
        assert not am.tools.ismonoclinic(self.boxtriclinic)

    def test_istriclinic(self):
        assert not am.tools.istriclinic(self.boxcubic)
        assert not am.tools.istriclinic(self.boxhexagonal)
        assert not am.tools.istriclinic(self.boxtetragonal)
        assert not am.tools.istriclinic(self.boxrhombohedral)
        assert not am.tools.istriclinic(self.boxorthorhombic)
        assert not am.tools.istriclinic(self.boxmonoclinic)
        assert am.tools.istriclinic(self.boxtriclinic)

    def test_identifyfamily(self):
        assert am.tools.identifyfamily(self.boxcubic) == 'cubic'
        assert am.tools.identifyfamily(self.boxhexagonal) == 'hexagonal'
        assert am.tools.identifyfamily(self.boxtetragonal) == 'tetragonal'
        assert am.tools.identifyfamily(self.boxrhombohedral) == 'rhombohedral'
        assert am.tools.identifyfamily(self.boxorthorhombic) == 'orthorhombic'
        assert am.tools.identifyfamily(self.boxmonoclinic) == 'monoclinic'
        assert am.tools.identifyfamily(self.boxtriclinic) == 'triclinic'