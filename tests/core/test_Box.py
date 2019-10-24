# coding: utf-8

# https://docs.pytest.org/en/latest/
import pytest

# http://www.numpy.org/
import numpy as np

import atomman as am
import atomman.unitconvert as uc

class Test_Box:
    def test_default(self):
        box = am.Box()
        assert np.allclose(box.avect, [1., 0., 0.])
        assert np.allclose(box.bvect, [0., 1., 0.])
        assert np.allclose(box.cvect, [0., 0., 1.])
        assert np.allclose(box.origin, [0., 0., 0.])
        assert np.allclose(box.vects, [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]])
        assert np.isclose(box.a, 1.0)
        assert np.isclose(box.b, 1.0)
        assert np.isclose(box.c, 1.0)
        assert np.isclose(box.alpha, 90.0)
        assert np.isclose(box.beta , 90.0)
        assert np.isclose(box.gamma, 90.0)
        assert np.isclose(box.xlo, 0.0)
        assert np.isclose(box.xhi, 1.0)
        assert np.isclose(box.ylo, 0.0)
        assert np.isclose(box.yhi, 1.0)
        assert np.isclose(box.zlo, 0.0)
        assert np.isclose(box.zhi, 1.0)
        assert np.isclose(box.lx, 1.0)
        assert np.isclose(box.ly, 1.0)
        assert np.isclose(box.lz, 1.0)
        assert np.isclose(box.xy, 0.0)
        assert np.isclose(box.xz, 0.0)
        assert np.isclose(box.yz, 0.0)
        assert np.isclose(box.volume, 1.0)

    def test_cubic(self):
        a = 3.25
        box = am.Box.cubic(a)
        assert np.isclose(box.a, a)
        assert np.isclose(box.b, a)
        assert np.isclose(box.c, a)
        assert np.isclose(box.alpha, 90.0)
        assert np.isclose(box.beta , 90.0)
        assert np.isclose(box.gamma, 90.0) 

    def test_hexagonal(self):
        a = 3.25
        c = 11.53
        box = am.Box.hexagonal(a, c)
        assert np.isclose(box.a, a)
        assert np.isclose(box.b, a)
        assert np.isclose(box.c, c)
        assert np.isclose(box.alpha, 90.0)
        assert np.isclose(box.beta , 90.0)
        assert np.isclose(box.gamma, 120.0) 

    def test_tetragonal(self):
        a = 3.25
        c = 11.53
        box = am.Box.tetragonal(a, c)
        assert np.isclose(box.a, a)
        assert np.isclose(box.b, a)
        assert np.isclose(box.c, c)
        assert np.isclose(box.alpha, 90.0)
        assert np.isclose(box.beta , 90.0)
        assert np.isclose(box.gamma, 90.0) 

    def test_trigonal(self):
        a = 3.25
        alpha = 96.3
        box = am.Box.trigonal(a, alpha)
        assert np.isclose(box.a, a)
        assert np.isclose(box.b, a)
        assert np.isclose(box.c, a)
        assert np.isclose(box.alpha, alpha)
        assert np.isclose(box.beta , alpha)
        assert np.isclose(box.gamma, alpha) 

    def test_orthorhombic(self):
        a = 3.25
        b = 5.26
        c = 11.53
        box = am.Box.orthorhombic(a, b, c)
        assert np.isclose(box.a, a)
        assert np.isclose(box.b, b)
        assert np.isclose(box.c, c)
        assert np.isclose(box.alpha, 90.0)
        assert np.isclose(box.beta , 90.0)
        assert np.isclose(box.gamma, 90.0) 

    def test_monoclinic(self):
        a = 3.25
        b = 5.26
        c = 11.53
        beta = 101.4
        box = am.Box.monoclinic(a, b, c, beta)
        assert np.isclose(box.a, a)
        assert np.isclose(box.b, b)
        assert np.isclose(box.c, c)
        assert np.isclose(box.alpha, 90.0)
        assert np.isclose(box.beta , beta)
        assert np.isclose(box.gamma, 90.0) 

    def test_triclinic(self):
        a = 3.25
        b = 5.26
        c = 11.53
        alpha = 96.3
        beta = 101.4
        gamma = 78.4
        box = am.Box.triclinic(a, b, c, alpha, beta, gamma)
        assert np.isclose(box.a, a)
        assert np.isclose(box.b, b)
        assert np.isclose(box.c, c)
        assert np.isclose(box.alpha, alpha)
        assert np.isclose(box.beta , beta)
        assert np.isclose(box.gamma, gamma) 
        
    def test_set_vectors(self):
        avect = [3.2, 0.0, 0.0]
        bvect = [0.1, 3.2, 0.0]
        cvect = [-0.2, 0.15, 3.2]
        origin = [3.1, 1.4, 4.1]
        box = am.Box(avect=avect, bvect=bvect, cvect=cvect, origin=origin)
        assert np.allclose(box.avect, avect)
        assert np.allclose(box.bvect, bvect)
        assert np.allclose(box.cvect, cvect)
        assert np.allclose(box.origin, origin)
        
        avect = [6.2, 0.0, 0.0]
        bvect = [0.0, 4.2, 0.0]
        cvect = [0.0, 0.1, 5.2]
        box.set_vectors(avect=avect, bvect=bvect, cvect=cvect)
        assert np.allclose(box.avect, avect)
        assert np.allclose(box.bvect, bvect)
        assert np.allclose(box.cvect, cvect)
        assert np.allclose(box.origin, np.zeros(3))

    def test_set_abc(self):
        a = 4.3
        b = 3.2
        c = 8.1
        alpha = 110
        origin = [3.1, 1.4, 4.1]
        box = am.Box(a=a, b=b, c=c, alpha=alpha, origin=origin)
        assert np.isclose(box.a, a)
        assert np.isclose(box.b, b)
        assert np.isclose(box.c, c)
        assert np.isclose(box.alpha, alpha)
        assert np.isclose(box.beta , 90.0)
        assert np.isclose(box.gamma, 90.0)
        assert np.allclose(box.origin, origin)
        
        a = 24.5
        b = 236.2
        c = 42.5
        beta = 101
        box = am.Box(a=a, b=b, c=c, beta=beta)
        assert np.isclose(box.a, a)
        assert np.isclose(box.b, b)
        assert np.isclose(box.c, c)
        assert np.isclose(box.alpha, 90.0)
        assert np.isclose(box.beta , beta)
        assert np.isclose(box.gamma, 90.0)
        assert np.allclose(box.origin, np.zeros(3))

    def test_set_lengths(self):
        lx=42
        ly=57
        lz=112
        xz=15
        origin=[1,2,3]
        box = am.Box(lx=lx, ly=ly, lz=lz, xz=xz, origin=origin)
        assert np.isclose(box.lx, lx)
        assert np.isclose(box.ly, ly)
        assert np.isclose(box.lz, lz)
        assert np.isclose(box.xz, xz)
        assert np.isclose(box.xy, 0.0)
        assert np.isclose(box.yz, 0.0)
        assert np.allclose(box.origin, origin)

        lx=163
        ly=347
        lz=235
        xy=19
        box.set_lengths(lx=lx, ly=ly, lz=lz, xy=xy)
        assert np.isclose(box.lx, lx)
        assert np.isclose(box.ly, ly)
        assert np.isclose(box.lz, lz)
        assert np.isclose(box.xz, 0.0)
        assert np.isclose(box.xy, xy)
        assert np.isclose(box.yz, 0.0)
        assert np.allclose(box.origin, np.zeros(3))

    def test_set_hi_los(self):
        xlo=-1
        xhi=5
        ylo=-2.1
        yhi=5
        zlo=0.1
        zhi=3.1
        xy=0.5
        box = am.Box(xlo=xlo, xhi=xhi, ylo=ylo, yhi=yhi, zlo=zlo, zhi=zhi, xy=xy)
        assert np.isclose(box.xlo, xlo)
        assert np.isclose(box.xhi, xhi)
        assert np.isclose(box.ylo, ylo)
        assert np.isclose(box.yhi, yhi)
        assert np.isclose(box.zlo, zlo)
        assert np.isclose(box.zhi, zhi)
        assert np.isclose(box.xz, 0.0)
        assert np.isclose(box.xy, xy)
        assert np.isclose(box.yz, 0.0)

        xlo=-143
        xhi=125
        ylo=-124
        yhi=364
        zlo=-172
        zhi=235
        yz=10
        box = am.Box(xlo=xlo, xhi=xhi, ylo=ylo, yhi=yhi, zlo=zlo, zhi=zhi, yz=yz)
        assert np.isclose(box.xlo, xlo)
        assert np.isclose(box.xhi, xhi)
        assert np.isclose(box.ylo, ylo)
        assert np.isclose(box.yhi, yhi)
        assert np.isclose(box.zlo, zlo)
        assert np.isclose(box.zhi, zhi)
        assert np.isclose(box.xz, 0.0)
        assert np.isclose(box.xy, 0.0)
        assert np.isclose(box.yz, yz)
    
    def test_is_lammps_norm(self):
        vects = np.array([[123, 0, 0], [4, 142, 0], [7, -9, 145]])
        box = am.Box(vects=vects)
        assert box.is_lammps_norm()

        vects = np.array([[123, 1, 0], [4, 142, 0], [7, -9, 145]])
        box = am.Box(vects=vects)
        assert not box.is_lammps_norm()

        vects = np.array([[123, 0, 1], [4, 142, 0], [7, -9, 145]])
        box = am.Box(vects=vects)
        assert not box.is_lammps_norm()

        vects = np.array([[123, 0, 0], [4, 142, 1], [7, -9, 145]])
        box = am.Box(vects=vects)
        assert not box.is_lammps_norm()

        vects = np.array([[-123, 0, 0], [4, 142, 0], [7, -9, 145]])
        box = am.Box(vects=vects)
        assert not box.is_lammps_norm()

        vects = np.array([[123, 0, 0], [4, -142, 0], [7, -9, 145]])
        box = am.Box(vects=vects)
        assert not box.is_lammps_norm()

        vects = np.array([[123, 0, 0], [4, 142, 0], [7, -9, -145]])
        box = am.Box(vects=vects)
        assert not box.is_lammps_norm()

    def test_model(self):
        box = am.Box.cubic(5.42)
        model = box.model()
        assert np.allclose(uc.value_unit(model['box']['avect']), np.array([5.42, 0.0, 0.0]))
        assert np.allclose(uc.value_unit(model['box']['bvect']), np.array([0.0, 5.42, 0.0]))
        assert np.allclose(uc.value_unit(model['box']['cvect']), np.array([0.0, 0.0, 5.42]))
        assert np.allclose(uc.value_unit(model['box']['origin']), np.array([0.0, 0.0, 0.0]))
        box = am.Box(model=model)
        assert np.allclose(box.vects, np.array([[5.42, 0.0, 0.0],
                                                [0.0, 5.42, 0.0], 
                                                [0.0, 0.0, 5.42]]))