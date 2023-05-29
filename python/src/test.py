#
# Copyright (c) 2014 Mitsuaki Kawamura
#
# Permission is hereby granted, free of charge, to any person obtaining a
# copy of this software and associated documentation files (the
# "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject to
# the following conditions:
#
# The above copyright notice and this permission notice shall be included
# in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
# OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
# CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
# TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
# SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#
import math
import numpy
import libtetrabz


def test_occ(nb, bvec, vbz, eig, mat):
    """

    :return:
    """
    #
    wght = libtetrabz.occ(bvec, eig - 0.5)
    #
    val = numpy.empty(2, dtype=numpy.float_)
    for ib in range(nb):
        val[ib] = numpy.sum(wght[:, :, :, ib] * mat[:, :, :])
    print("# libtetrabz_occ")
    print('     %15.5e %15.5e' % (4.0 * numpy.pi / 5.0, val[0] * vbz))
    print('     %15.5e %15.5e' % (numpy.pi / (5.0 * math.sqrt(2.0)), val[1] * vbz))
    print("")


def test_fermieng(nb, bvec, vbz, eig, mat):
    """

    :return:
    """
    #
    nelec = (4.0 * numpy.pi / 3.0 + math.sqrt(2.0) * numpy.pi / 3.0) / vbz

    ef, wght, iterations = libtetrabz.fermieng(bvec, eig, nelec)
    #
    val = numpy.empty(2, dtype=numpy.float_)
    for ib in range(nb):
        val[ib] = numpy.sum(wght[:, :, :, ib] * mat[:, :, :])
    #
    print("# libtetrabz_fermieng")
    print("     %15.5e %15.5e" % (0.5, ef))
    print("     %15.5e %15.5e" % (4.0 * numpy.pi / 5.0, val[0] * vbz))
    print("     %15.5e %15.5e" % (numpy.pi / (5.0 * math.sqrt(2.0)), val[1] * vbz))
    print("")


def test_dos(nb, bvec, vbz, eig, mat):
    """

    :return:
    """
    ne = 5
    #
    e0 = numpy.empty(ne, dtype=numpy.float_)
    val0 = numpy.empty([nb, ne], dtype=numpy.float_)
    val = numpy.empty([nb, ne], dtype=numpy.float_)
    for ie in range(ne):
        e0[ie] = 0.2 * (ie + 1)
        val0[0, ie] = 4.0 * numpy.pi * e0[ie]**3
        if e0[ie] > 1.0 / math.sqrt(2.0):
            val0[1, ie] = math.sqrt(2.0) * numpy.pi * math.sqrt(-1.0 + 2.0 * e0[ie]**2)**3
        else:
            val0[1, ie] = 0.0
    e0[0:ne] = e0[0:ne]**2 * 0.5
    #
    wght = libtetrabz.dos(bvec, eig, e0)
    #
    for ib in range(nb):
        for ie in range(ne):
            val[ib, ie] = numpy.sum(wght[:, :, :, ib, ie] * mat[:, :, :])
    #
    print("# libtetrabz_dos")
    for ib in range(nb):
        for ie in range(ne):
            print("     %15.5e %15.5e" % (val0[ib, ie], val[ib, ie] * vbz))
    print("")


def test_intdos(nb, bvec, vbz, eig, mat):
    """

    :return:
    """
    ne = 5
    e0 = numpy.empty(ne, dtype=numpy.float_)
    val0 = numpy.empty([nb, ne], dtype=numpy.float_)
    val = numpy.empty([nb, ne], dtype=numpy.float_)
    #
    for ie in range(ne):
        e0[ie] = 0.2 * (ie + 1)
        val0[0, ie] = 4.0 * numpy.pi * e0[ie]**5 / 5.0
        if e0[ie] > 1.0 / math.sqrt(2.0):
            val0[1, ie] = numpy.pi * math.sqrt(-1.0 + 2.0 * e0[ie]**2)**5 / (5.0 * math.sqrt(2.0))
        else:
            val0[1, ie] = 0.0
    e0[0:ne] = e0[0:ne]**2 * 0.5
    #
    wght = libtetrabz.intdos(bvec, eig, e0)
    #
    for ib in range(nb):
        for ie in range(ne):
            val[ib, ie] = numpy.sum(wght[:, :, :, ib, ie] * mat[:, :, :])
    #
    print("# libtetrabz_intdos")
    for ib in range(nb):
        for ie in range(ne):
            print("     %15.5e %15.5e" % (val0[ib, ie], val[ib, ie] * vbz))
    print("")


def test_dblstep(nb, bvec, vbz, eig1, eig2, mat):
    """

    :param bvec:
    :param vbz:
    :param eig1:
    :param eig2:
    :param mat:
    :return:
    """
    #
    wght = libtetrabz.dblstep(bvec, eig1 - 0.5, eig2 - 0.5)
    #
    val = numpy.empty([nb, nb], dtype=numpy.float_)
    for ib in range(nb):
        for jb in range(nb):
            val[ib, jb] = numpy.sum(wght[:, :, :, ib, jb] * mat[:, :, :])
    #
    print("# libtetrabz_dblstep")
    print("     %15.5e %15.5e" % (49.0 * numpy.pi / 320.0, val[0, 0] * vbz))
    print("     %15.5e %15.5e" % (0.0, val[0, 1] * vbz))
    print("     %15.5e %15.5e" % (numpy.pi * (512.0 * math.sqrt(2.0) - 319.0) / 10240.0, val[1, 0] * vbz))
    print("     %15.5e %15.5e" % (0.0, val[1, 1] * vbz))
    print("")


def test_dbldelta(nb, bvec, vbz, eig1, eig2, mat):
    """

    :param bvec:
    :param vbz:
    :param eig1:
    :param eig2:
    :param mat:
    :return:
    """
    #
    wght = libtetrabz.dbldelta(bvec, eig1 - 0.5, eig2 - 0.5)
    #
    val = numpy.empty([nb, nb], dtype=numpy.float_)
    for ib in range(nb):
        for jb in range(nb):
            val[ib, jb] = numpy.sum(wght[:, :, :, ib, jb] * mat[:, :, :])
    print("# libtetrabz_dbldelta")
    print("     %15.5e %15.5e" % (2.0 * numpy.pi, val[0, 0] * vbz))
    print("     %15.5e %15.5e" % (0.0, val[0, 1] * vbz))
    print("     %15.5e %15.5e" % (numpy.pi, val[1, 0] * vbz))
    print("     %15.5e %15.5e" % (0.0, val[1, 1] * vbz))
    print("")


def test_polstat(nb, bvec, vbz, eig1, eig2, mat):
    """

    :param bvec:
    :param vbz:
    :param eig1:
    :param eig2:
    :param mat:
    :return:
    """
    wght = libtetrabz.polstat(bvec, eig1 - 0.5, eig2 - 0.5)
    #
    val = numpy.empty([nb, nb], dtype=numpy.float_)
    for ib in range(nb):
        for jb in range(nb):
            val[ib, jb] = numpy.sum(wght[:, :, :, ib, jb] * mat[:, :, :])
    #
    print("# libtetrabz_polstat")
    print("     %15.5e %15.5e" % (numpy.pi * (68.0 + 45.0 * math.log(3.0)) / 96.0, val[0, 0] * vbz))
    print("     %15.5e %15.5e" % (numpy.pi * 8.0 / 5.0, val[0, 1] * vbz))
    print("     %15.5e %15.5e" % (numpy.pi * (228.0 + 22.0 * math.sqrt(2.0) - 96.0 * math.log(2.0)
                                              + 192.0 * math.log(4.0 + math.sqrt(2.0))
                                              - 3.0 * math.log(1.0 + 2.0 * math.sqrt(2.0))) / 1536.0,
                                  val[1, 0] * vbz))
    print("     %15.5e %15.5e" % (numpy.pi * math.sqrt(8.0) / 5.0, val[1, 1] * vbz))
    print("")


def test_fermigr(nb, bvec, vbz, eig1, eig2, mat):
    """

    :param bvec:
    :param vbz:
    :param eig1:
    :param eig2:
    :param mat:
    :return:
    """
    ne = 3
    #
    e0 = numpy.empty(ne, dtype=numpy.float_)
    val0 = numpy.empty([nb, nb, ne], dtype=numpy.float_)
    val = numpy.empty([nb, nb, ne], dtype=numpy.float_)
    e0[0] = 1.0 / 3.0
    e0[1] = 2.0 / 3.0
    e0[2] = 1.0
    val0[0, 0, 0] = 4.0 * numpy.pi / 9.0
    val0[0, 0, 1] = 1295.0 * numpy.pi / 2592.0
    val0[0, 0, 2] = 15.0 * numpy.pi / 32.0
    val0[1, 0, 0] = 5183.0 * numpy.pi / 41472.0
    val0[1, 0, 1] = 4559.0 * numpy.pi / 41472.0
    val0[1, 0, 2] = 0.0
    val0[0:nb, 1, 0:3] = 0.0
    #
    wght = libtetrabz.fermigr(bvec, eig1 - 0.5, eig2 - 0.5, e0)
    #
    for ib in range(nb):
        for jb in range(nb):
            for ie in range(ne):
                val[ib, jb, ie] = numpy.sum(wght[:, :, :, ib, jb, ie] * mat[:, :, :])
    #
    print("# libtetrabz_fermigr")
    for ib in range(nb):
        for jb in range(nb):
            for ie in range(ne):
                print("     %15.5e %15.5e" % (val0[ib, jb, ie], val[ib, jb, ie] * vbz))
    print("")


def test_polcmplx(nb, bvec, vbz, eig1, eig2, mat):
    """

    :param bvec:
    :param vbz:
    :param eig1:
    :param eig2:
    :param mat:
    :return:
    """
    ne = 3
    e0 = numpy.empty(ne, dtype=numpy.complex_)
    val0 = numpy.empty([nb, nb, ne], dtype=numpy.complex_)
    val = numpy.empty([nb, nb, ne], dtype=numpy.complex_)
    e0[0] = -2.0 + 1.0j
    e0[1] = 0.0 + 2.0j
    e0[2] = 1.0 - 0.5j
    val0[0, 0, 0] = -0.838243341280338 - 0.734201894333234j
    val0[0, 0, 1] = 0.270393588876530 - 0.771908416949610j
    val0[0, 0, 2] = 0.970996830573510 + 0.302792326476720j
    val0[1, 0, 0] = -0.130765724778920 - 0.087431218706638j
    val0[1, 0, 1] = 0.030121954547245 - 0.135354254293510j
    val0[1, 0, 2] = 0.178882244951203 + 0.064232167683425j
    val0[0, 1, 0:3] = (8.0 * numpy.pi) / (5.0 * (1.0 + 2.0 * e0[0:3]))
    val0[1, 1, 0:3] = (math.sqrt(8.0) * numpy.pi) / (5.0 * (1.0 + 4.0 * e0[0:3]))
    #
    wght = libtetrabz.polcmplx(bvec, eig1 - 0.5, eig2 - 0.5, e0)
    #
    for ib in range(nb):
        for jb in range(nb):
            for ie in range(ne):
                val[ib, jb, ie] = numpy.sum(wght[:, :, :, ib, jb, ie] * mat[:, :, :])
    #
    print("# libtetrabz_polcmplx")
    for ib in range(nb):
        for jb in range(nb):
            for ie in range(ne):
                print("     %15.5e %15.5e" % (val0[ib, jb, ie].real, val[ib, jb, ie].real * vbz))
                print("     %15.5e %15.5e" % (val0[ib, jb, ie].imag, val[ib, jb, ie].imag * vbz))
    print("")


def test():
    ng0 = 8
    ng = numpy.array([ng0, ng0, ng0])
    nb = 2
    bvec = numpy.array([[3.0, 0.0, 0.0],
                        [0.0, 3.0, 0.0],
                        [0.0, 0.0, 3.0]])
    vbz = abs(numpy.linalg.det(bvec))

    #
    eig1 = numpy.empty([ng[0], ng[1], ng[2], nb], dtype=numpy.float_)
    eig2 = numpy.empty([ng[0], ng[1], ng[2], nb], dtype=numpy.float_)
    mat = numpy.empty([ng[0], ng[1], ng[2]], dtype=numpy.float_)
    #
    kvec = numpy.empty(3, dtype=numpy.float_)
    for i0 in range(ng[1]):
        for i1 in range(ng[1]):
            for i2 in range(ng[2]):
                #
                kvec[0:3] = numpy.array([i0, i1, i2]) / ng[0:3]
                for ii in range(3):
                    if kvec[ii] > 0.5:
                        kvec[ii] += -1
                kvec[0:3] = kvec.dot(bvec)
                #
                eig1[i0, i1, i2, 0] = 0.5 * kvec.dot(kvec)
                eig1[i0, i1, i2, 1] = eig1[i0, i1, i2, 0] + 0.25
                #
                kvec[0] += 1.0
                eig2[i0, i1, i2, 0] = 0.5 * kvec.dot(kvec)
                eig2[i0, i1, i2, 1] = eig1[i0, i1, i2, 0] + 0.5
                #
    #
    for i0 in range(ng[1]):
        for i1 in range(ng[1]):
            for i2 in range(ng[2]):
                #
                kvec[0:3] = numpy.array([i0, i1, i2]) / ng[0:3]
                for ii in range(3):
                    if kvec[ii] > 0.5:
                        kvec[ii] += -1
                kvec[0:3] = kvec.dot(bvec)
                #
                mat[i0, i1, i2] = kvec.dot(kvec)
    #
    print("#          Ideal          Result")
    #
    test_occ(nb, bvec, vbz, eig1, mat)
    #
    test_fermieng(nb, bvec, vbz, eig1, mat)
    #
    test_dos(nb, bvec, vbz, eig1, mat)
    #
    test_intdos(nb, bvec, vbz, eig1, mat)
    #
    test_dblstep(nb, bvec, vbz, eig1, eig2, mat)
    #
    test_dbldelta(nb, bvec, vbz, eig1, eig2, mat)
    #
    test_polstat(nb, bvec, vbz, eig1, eig2, mat)
    #
    test_fermigr(nb, bvec, vbz, eig1, eig2, mat)
    #
    test_polcmplx(nb, bvec, vbz, eig1, eig2, mat)
    #


test()
