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
from .libtetrabzc import dos_c
from .libtetrabzc import intdos_c
from .libtetrabzc import occ_c
from .libtetrabzc import fermieng_c
from .libtetrabzc import polstat_c
from .libtetrabzc import fermigr_c
from .libtetrabzc import dbldelta_c
from .libtetrabzc import dblstep_c
from .libtetrabzc import polcmplx_c
import numpy


def occ(bvec=numpy.array([1.0, 0.0, 0.0]), eig=numpy.array([0.0])):
    """

    :return:
    """
    #
    ng = numpy.array(eig.shape[0:3])
    nk = ng.prod(0)
    nb = eig.shape[3]
    eig_c = eig.reshape(nk*nb).tolist()
    #
    wght_c = libtetrabzc.occ_c(ng[0], ng[1], ng[2], nk, nb,
                               bvec[0, 0], bvec[0, 1], bvec[0, 2],
                               bvec[1, 0], bvec[1, 1], bvec[1, 2],
                               bvec[2, 0], bvec[2, 1], bvec[2, 2], eig_c)
    #
    wght = numpy.array(wght_c).reshape([ng[0], ng[1], ng[2], nb])
    return wght


def fermieng(bvec=numpy.array([1.0, 0.0, 0.0]), eig=numpy.array([0.0]), nelec=0.0):
    """

    :return:
    """
    #
    ng = numpy.array(eig.shape[0:3])
    nk = ng.prod(0)
    nb = eig.shape[3]
    eig_c = eig.reshape(nk*nb).tolist()

    wght_c = libtetrabzc.fermieng_c(ng[0], ng[1], ng[2], nk, nb,
                                    bvec[0, 0], bvec[0, 1], bvec[0, 2],
                                    bvec[1, 0], bvec[1, 1], bvec[1, 2],
                                    bvec[2, 0], bvec[2, 1], bvec[2, 2], eig_c, nelec)
    #
    ef = wght_c[nk * nb]
    iteration = wght_c[nk*nb + 1]
    wght = numpy.array(wght_c)[0:nk*nb].reshape([ng[0], ng[1], ng[2], nb])
    return ef, wght, iteration


def dos(bvec=numpy.array([1.0, 0.0, 0.0]), eig=numpy.array([0.0]), e0=numpy.array([0.0])):
    """

    :return:
    """
    ng = numpy.array(eig.shape[0:3])
    nk = ng.prod(0)
    nb = eig.shape[3]
    ne = e0.shape[0]
    eig_c = eig.reshape(nk*nb).tolist()
    e0_c = e0.tolist()
    #
    wght_c = libtetrabzc.dos_c(ng[0], ng[1], ng[2], nk, nb, ne,
                               bvec[0, 0], bvec[0, 1], bvec[0, 2],
                               bvec[1, 0], bvec[1, 1], bvec[1, 2],
                               bvec[2, 0], bvec[2, 1], bvec[2, 2], eig_c, e0_c)
    #
    wght = numpy.array(wght_c).reshape([ng[0], ng[1], ng[2], nb, ne])
    return wght


def intdos(bvec=numpy.array([1.0, 0.0, 0.0]), eig=numpy.array([0.0]), e0=numpy.array([0.0])):
    """

    :return:
    """
    ng = numpy.array(eig.shape[0:3])
    nk = ng.prod(0)
    nb = eig.shape[3]
    ne = e0.shape[0]
    eig_c = eig.reshape(nk*nb).tolist()
    e0_c = e0.tolist()
    #
    wght_c = libtetrabzc.intdos_c(ng[0], ng[1], ng[2], nk, nb, ne,
                                  bvec[0, 0], bvec[0, 1], bvec[0, 2],
                                  bvec[1, 0], bvec[1, 1], bvec[1, 2],
                                  bvec[2, 0], bvec[2, 1], bvec[2, 2], eig_c, e0_c)
    #
    wght = numpy.array(wght_c).reshape([ng[0], ng[1], ng[2], nb, ne])
    return wght


def dblstep(bvec=numpy.array([1.0, 0.0, 0.0]), eig1=numpy.array([0.0]), eig2=numpy.array([0.0])):
    """

    :param bvec:
    :param eig1:
    :param eig2:
    :return:
    """
    #
    ng = numpy.array(eig1.shape[0:3])
    nk = ng.prod(0)
    nb = eig1.shape[3]
    eig1_c = eig1.reshape(nk*nb).tolist()
    eig2_c = eig2.reshape(nk*nb).tolist()

    wght_c = libtetrabzc.dblstep_c(ng[0], ng[1], ng[2], nk, nb,
                                   bvec[0, 0], bvec[0, 1], bvec[0, 2],
                                   bvec[1, 0], bvec[1, 1], bvec[1, 2],
                                   bvec[2, 0], bvec[2, 1], bvec[2, 2], eig1_c, eig2_c)
    #
    wght = numpy.array(wght_c).reshape([ng[0], ng[1], ng[2], nb, nb])
    return wght


def dbldelta(bvec=numpy.array([1.0, 0.0, 0.0]), eig1=numpy.array([0.0]), eig2=numpy.array([0.0])):
    """

    :param bvec:
    :param eig1:
    :param eig2:
    :return:
    """
    #
    ng = numpy.array(eig1.shape[0:3])
    nk = ng.prod(0)
    nb = eig1.shape[3]
    eig1_c = eig1.reshape(nk * nb).tolist()
    eig2_c = eig2.reshape(nk * nb).tolist()

    wght_c = libtetrabzc.dbldelta_c(ng[0], ng[1], ng[2], nk, nb,
                                    bvec[0, 0], bvec[0, 1], bvec[0, 2],
                                    bvec[1, 0], bvec[1, 1], bvec[1, 2],
                                    bvec[2, 0], bvec[2, 1], bvec[2, 2], eig1_c, eig2_c)
    #
    wght = numpy.array(wght_c).reshape([ng[0], ng[1], ng[2], nb, nb])
    return wght


def polstat(bvec=numpy.array([1.0, 0.0, 0.0]), eig1=numpy.array([0.0]), eig2=numpy.array([0.0])):
    """

    :param bvec:
    :param eig1:
    :param eig2:
    :return:
    """
    ng = numpy.array(eig1.shape[0:3])
    nk = ng.prod(0)
    nb = eig1.shape[3]
    eig1_c = eig1.reshape(nk * nb).tolist()
    eig2_c = eig2.reshape(nk * nb).tolist()

    wght_c = libtetrabzc.polstat_c(ng[0], ng[1], ng[2], nk, nb,
                                   bvec[0, 0], bvec[0, 1], bvec[0, 2],
                                   bvec[1, 0], bvec[1, 1], bvec[1, 2],
                                   bvec[2, 0], bvec[2, 1], bvec[2, 2], eig1_c, eig2_c)
    #
    wght = numpy.array(wght_c).reshape([ng[0], ng[1], ng[2], nb, nb])
    return wght


def fermigr(bvec=numpy.array([1.0, 0.0, 0.0]),
            eig1=numpy.array([0.0]), eig2=numpy.array([0.0]), e0=numpy.array([0.0])):
    """

    :param bvec:
    :param eig1:
    :param eig2:
    :param e0:
    :return:
    """
    ng = numpy.array(eig1.shape[0:3])
    nk = ng.prod(0)
    nb = eig1.shape[3]
    ne = e0.shape[0]
    eig1_c = eig1.reshape(nk * nb).tolist()
    eig2_c = eig2.reshape(nk * nb).tolist()
    e0_c = e0.tolist()

    wght_c = libtetrabzc.fermigr_c(ng[0], ng[1], ng[2], nk, nb, ne,
                                   bvec[0, 0], bvec[0, 1], bvec[0, 2],
                                   bvec[1, 0], bvec[1, 1], bvec[1, 2],
                                   bvec[2, 0], bvec[2, 1], bvec[2, 2], eig1_c, eig2_c, e0_c)
    #
    wght = numpy.array(wght_c).reshape([ng[0], ng[1], ng[2], nb, nb, ne])
    return wght


def polcmplx(bvec=numpy.array([1.0, 0.0, 0.0]),
             eig1=numpy.array([0.0]), eig2=numpy.array([0.0]), e0=numpy.array([0.0+0.0j])):
    """

    :param bvec:
    :param eig1:
    :param eig2:
    :param e0:
    :return:
    """
    ng = numpy.array(eig1.shape[0:3])
    nk = ng.prod(0)
    nb = eig1.shape[3]
    ne = e0.shape[0]
    eig1_c = eig1.reshape(nk * nb).tolist()
    eig2_c = eig2.reshape(nk * nb).tolist()

    e0_1 = numpy.empty([ne, 2], dtype=numpy.float_)
    e0_1[:, 0] = e0[:].real
    e0_1[:, 1] = e0[:].imag
    e0_c = e0_1.reshape(ne * 2).tolist()

    wght_c = libtetrabzc.polcmplx_c(ng[0], ng[1], ng[2], nk, nb, ne,
                                    bvec[0, 0], bvec[0, 1], bvec[0, 2],
                                    bvec[1, 0], bvec[1, 1], bvec[1, 2],
                                    bvec[2, 0], bvec[2, 1], bvec[2, 2], eig1_c, eig2_c, e0_c)
    #
    wght_1 = numpy.array(wght_c).reshape([ng[0], ng[1], ng[2], nb, nb, ne, 2])
    #
    wght = wght_1[:, :, :, :, :, :, 0] + 1.0j * wght_1[:, :, :, :, :, :, 1]
    return wght
