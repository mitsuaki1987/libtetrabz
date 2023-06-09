API usage
=========

You can call a function in this library as follows:

.. code-block:: python

    import libtetrabz
..

Total energy, charge density, occupations
-----------------------------------------

.. math::
    \begin{align}
    \sum_{n}
    \int_{\rm BZ} \frac{d^3 k}{V_{\rm BZ}}
    \theta(\varepsilon_{\rm F} -
    \varepsilon_{n k}) X_{n k}
    \end{align}
..

.. code-block:: python

    wght = libtetrabz.occ(bvec, eig)
..

Parameters

    .. code-block:: python

        bvec = numpy.array([[b0x, b0y, b0z], [b1x, b1y, b1z], [b2x, b2y, b2z]])
    ..

        Reciprocal lattice vectors (arbitrary unit).
        Because they are used to choose the direction of tetrahedra,
        only their ratio is used.

    .. code-block:: python
                   
        eig = numpy.empty([ng0, ng1, ng2, nb], dtype=numpy.float_)
    ..

        The energy eigenvalue :math:`\varepsilon_{n k}` in arbitrary unit.
        Where `ng0, ng1, ng2` are the number of grid for each direction of
        the reciprocal lattice vectors and `nb` is the number of bands.

Return

    .. code-block:: python

        wght = numpy.empty([ng0, ng1, ng2, nb], dtype=numpy.float_)
    ..
   
      The integration weights.

Fermi energy and occupations
----------------------------

.. math::

   \begin{align}
   \sum_{n}
   \int_{\rm BZ} \frac{d^3 k}{V_{\rm BZ}}
   \theta(\varepsilon_{\rm F} -
   \varepsilon_{n k}) X_{n k} 
   \end{align}

.. code-block:: python

    ef, wght, iteration = libtetrabz.fermieng(bvec, eig, nelec)
..

Parameters

    .. code-block:: python

        bvec = numpy.array([[b0x, b0y, b0z], [b1x, b1y, b1z], [b2x, b2y, b2z]])
    ..

        Reciprocal lattice vectors (arbitrary unit).
        Because they are used to choose the direction of tetrahedra,
        only their ratio is used.

    .. code-block:: python

        eig = numpy.empty([ng0, ng1, ng2, nb], dtype=numpy.float_)
    ..

        The energy eigenvalue :math:`\varepsilon_{n k}` in arbitrary unit.
        Where `ng0, ng1, ng2` are the number of grid for each direction of
        the reciprocal lattice vectors and `nb` is the number of bands.

    .. code-block:: python

        nelec = # float
    ..

        The number of electrons.

Return

    .. code-block:: python

        ef = # float
    ..

        The Fermi energy. Unit is the same as `eig`.

    .. code-block:: python

        wght = numpy.empty([ng0, ng1, ng2, nb], dtype=numpy.float_)
    ..

      The integration weights.

Partial density of states
-------------------------

.. math::

   \begin{align}
   \sum_{n}
   \int_{\rm BZ} \frac{d^3 k}{V_{\rm BZ}}
   \delta(\omega - \varepsilon_{n k})
   X_{n k}(\omega) 
   \end{align}

.. code-block:: python

    wght = libtetrabz.dos(bvec, eig, e0)
..

Parameters

    .. code-block:: python

        bvec = numpy.array([[b0x, b0y, b0z], [b1x, b1y, b1z], [b2x, b2y, b2z]])
    ..

        Reciprocal lattice vectors (arbitrary unit).
        Because they are used to choose the direction of tetrahedra,
        only their ratio is used.

    .. code-block:: python

        eig = numpy.empty([ng0, ng1, ng2, nb], dtype=numpy.float_)
    ..

        The energy eigenvalue :math:`\varepsilon_{n k}` in arbitrary unit.
        Where `ng0, ng1, ng2` are the number of grid for each direction of
        the reciprocal lattice vectors and `nb` is the number of bands.

    .. code-block:: python

        e0 = numpy.empty(ne, dtype=numpy.float_)
    ..

        The energy point :math:`\omega` in the same unit of `eig`.
        Where `ne` is the number of energy points.

Return

    .. code-block:: python

        wght = numpy.empty([ng0, ng1, ng2, nb, ne], dtype=numpy.float_)
    ..

      The integration weights.

Integrated density of states
----------------------------

.. math::

   \begin{align}
   \sum_{n}
   \int_{\rm BZ} \frac{d^3 k}{V_{\rm BZ}}
   \theta(\omega - \varepsilon_{n k})
   X_{n k}(\omega) 
   \end{align}

.. code-block:: python

    wght = libtetrabz.intdos(bvec, eig, e0)
..

Parameters

    .. code-block:: python

        bvec = numpy.array([[b0x, b0y, b0z], [b1x, b1y, b1z], [b2x, b2y, b2z]])
    ..

        Reciprocal lattice vectors (arbitrary unit).
        Because they are used to choose the direction of tetrahedra,
        only their ratio is used.

    .. code-block:: python

        eig = numpy.empty([ng0, ng1, ng2, nb], dtype=numpy.float_)
    ..

        The energy eigenvalue :math:`\varepsilon_{n k}` in arbitrary unit.
        Where `ng0, ng1, ng2` are the number of grid for each direction of
        the reciprocal lattice vectors and `nb` is the number of bands.

    .. code-block:: python

        e0 = numpy.empty(ne, dtype=numpy.float_)
    ..

        The energy point :math:`\omega` in the same unit of `eig`.
        Where `ne` is the number of energy points.

Return

    .. code-block:: python

        wght = numpy.empty([ng0, ng1, ng2, nb, ne], dtype=numpy.float_)
    ..

      The integration weights.

Nesting function and Fr&oumlhlich parameter
-------------------------------------------

.. math::

   \begin{align}
   \sum_{n n'}
   \int_{\rm BZ} \frac{d^3 k}{V_{\rm BZ}}
   \delta(\varepsilon_{\rm F} -
   \varepsilon_{n k}) \delta(\varepsilon_{\rm F} - \varepsilon'_{n' k})
   X_{n n' k}
   \end{align}

.. code-block:: python

    wght = libtetrabz.dbldelta(bvec, eig1, eig2)
..

Parameters

    .. code-block:: python

        bvec = numpy.array([[b0x, b0y, b0z], [b1x, b1y, b1z], [b2x, b2y, b2z]])
    ..

        Reciprocal lattice vectors (arbitrary unit).
        Because they are used to choose the direction of tetrahedra,
        only their ratio is used.

    .. code-block:: python

        eig1 = numpy.empty([ng0, ng1, ng2, nb], dtype=numpy.float_)
    ..

        The energy eigenvalue :math:`\varepsilon_{n k}` in arbitrary unit.
        Where `ng0, ng1, ng2` are the number of grid for each direction of
        the reciprocal lattice vectors and `nb` is the number of bands.

    .. code-block:: python

        eig2 = numpy.empty([ng0, ng1, ng2, nb], dtype=numpy.float_)
    ..

        Another energy eigenvalue :math:`\varepsilon'_{n k}` in the same unit of `eig1`.
        Typically it is assumed to be :math:`\varepsilon_{n k+q}`.

Return

    .. code-block:: python

        wght = numpy.empty([ng0, ng1, ng2, nb, nb], dtype=numpy.float_)
    ..

      The integration weights.

A part of DFPT calculation
--------------------------

.. math::

   \begin{align}
   \sum_{n n'}
   \int_{\rm BZ} \frac{d^3 k}{V_{\rm BZ}}
   \theta(\varepsilon_{\rm F} -
   \varepsilon_{n k}) \theta(\varepsilon_{n k} - \varepsilon'_{n' k}) 
   X_{n n' k}
   \end{align}

.. code-block:: python

    wght = libtetrabz.dblstep(bvec, eig1, eig2)
..

Parameters

    .. code-block:: python

        bvec = numpy.array([[b0x, b0y, b0z], [b1x, b1y, b1z], [b2x, b2y, b2z]])
    ..

        Reciprocal lattice vectors (arbitrary unit).
        Because they are used to choose the direction of tetrahedra,
        only their ratio is used.

    .. code-block:: python

        eig1 = numpy.empty([ng0, ng1, ng2, nb], dtype=numpy.float_)
    ..

        The energy eigenvalue :math:`\varepsilon_{n k}` in arbitrary unit.
        Where `ng0, ng1, ng2` are the number of grid for each direction of
        the reciprocal lattice vectors and `nb` is the number of bands.

    .. code-block:: python

        eig2 = numpy.empty([ng0, ng1, ng2, nb], dtype=numpy.float_)
    ..

        Another energy eigenvalue :math:`\varepsilon'_{n k}` in the same unit of `eig1`.
        Typically it is assumed to be :math:`\varepsilon_{n k+q}`.

Return

    .. code-block:: python

        wght = numpy.empty([ng0, ng1, ng2, nb, nb], dtype=numpy.float_)
    ..

      The integration weights.

Static polarization function
----------------------------

.. math::

   \begin{align}
   \sum_{n n'}
   \int_{\rm BZ} \frac{d^3 k}{V_{\rm BZ}}
   \frac{\theta(\varepsilon_{\rm F} - \varepsilon_{n k})
   \theta(\varepsilon'_{n' k} - \varepsilon_{\rm F})}
   {\varepsilon'_{n' k} - \varepsilon_{n k}}
   X_{n n' k} 
   \end{align}

.. code-block:: python

    wght = libtetrabz.polstat(bvec, eig1, eig2)
..

Parameters

    .. code-block:: python

        bvec = numpy.array([[b0x, b0y, b0z], [b1x, b1y, b1z], [b2x, b2y, b2z]])
    ..

        Reciprocal lattice vectors (arbitrary unit).
        Because they are used to choose the direction of tetrahedra,
        only their ratio is used.

    .. code-block:: python

        eig1 = numpy.empty([ng0, ng1, ng2, nb], dtype=numpy.float_)
    ..

        The energy eigenvalue :math:`\varepsilon_{n k}` in arbitrary unit.
        Where `ng0, ng1, ng2` are the number of grid for each direction of
        the reciprocal lattice vectors and `nb` is the number of bands.

    .. code-block:: python

        eig2 = numpy.empty([ng0, ng1, ng2, nb], dtype=numpy.float_)
    ..

        Another energy eigenvalue :math:`\varepsilon'_{n k}` in the same unit of `eig1`.
        Typically it is assumed to be :math:`\varepsilon_{n k+q}`.

Return

    .. code-block:: python

        wght = numpy.empty([ng0, ng1, ng2, nb, nb], dtype=numpy.float_)
    ..

      The integration weights.

Phonon linewidth
----------------

.. math::

   \begin{align}
   \sum_{n n'}
   \int_{\rm BZ} \frac{d^3 k}{V_{\rm BZ}}
   \theta(\varepsilon_{\rm F} -
   \varepsilon_{n k}) \theta(\varepsilon'_{n' k} - \varepsilon_{\rm F})
   \delta(\varepsilon'_{n' k} - \varepsilon_{n k} - \omega)
   X_{n n' k}(\omega) 
   \end{align}

.. code-block:: python

    wght = libtetrabz.fermigr(bvec, eig1, eig2, e0)
..

Parameters

    .. code-block:: python

        bvec = numpy.array([[b0x, b0y, b0z], [b1x, b1y, b1z], [b2x, b2y, b2z]])
    ..

        Reciprocal lattice vectors (arbitrary unit).
        Because they are used to choose the direction of tetrahedra,
        only their ratio is used.

    .. code-block:: python

        eig1 = numpy.empty([ng0, ng1, ng2, nb], dtype=numpy.float_)
    ..

        The energy eigenvalue :math:`\varepsilon_{n k}` in arbitrary unit.
        Where `ng0, ng1, ng2` are the number of grid for each direction of
        the reciprocal lattice vectors and `nb` is the number of bands.

    .. code-block:: python

        eig2 = numpy.empty([ng0, ng1, ng2, nb], dtype=numpy.float_)
    ..

        Another energy eigenvalue :math:`\varepsilon'_{n k}` in the same unit of `eig1`.
        Typically it is assumed to be :math:`\varepsilon_{n k+q}`.

    .. code-block:: python

        e0 = numpy.empty(ne, dtype=numpy.float_)
    ..

        The energy point :math:`\omega` in the same unit of `eig`.
        Where `ne` is the number of energy points.

Return

    .. code-block:: python

        wght = numpy.empty([ng0, ng1, ng2, nb, nb, ne], dtype=numpy.float_)
    ..

      The integration weights.

Polarization function (complex frequency)
-----------------------------------------

.. math::

   \begin{align}
   \sum_{n n'}
   \int_{\rm BZ} \frac{d^3 k}{V_{\rm BZ}}
   \frac{\theta(\varepsilon_{\rm F} - \varepsilon_{n k})
   \theta(\varepsilon'_{n' k} - \varepsilon_{\rm F})}
   {\varepsilon'_{n' k} - \varepsilon_{n k} + i \omega}
   X_{n n' k}(\omega) 
   \end{align}

.. code-block:: python

    wght = libtetrabz.polcmplex(bvec, eig1, eig2, e0)
..

Parameters

    .. code-block:: python

        bvec = numpy.array([[b0x, b0y, b0z], [b1x, b1y, b1z], [b2x, b2y, b2z]])
    ..

        Reciprocal lattice vectors (arbitrary unit).
        Because they are used to choose the direction of tetrahedra,
        only their ratio is used.

    .. code-block:: python

        eig1 = numpy.empty([ng0, ng1, ng2, nb], dtype=numpy.float_)
    ..

        The energy eigenvalue :math:`\varepsilon_{n k}` in arbitrary unit.
        Where `ng0, ng1, ng2` are the number of grid for each direction of
        the reciprocal lattice vectors and `nb` is the number of bands.

    .. code-block:: python

        eig2 = numpy.empty([ng0, ng1, ng2, nb], dtype=numpy.float_)
    ..

        Another energy eigenvalue :math:`\varepsilon'_{n k}` in the same unit of `eig1`.
        Typically it is assumed to be :math:`\varepsilon_{n k+q}`.

    .. code-block:: python

        e0 = numpy.empty(ne, dtype=numpy.complex_)
    ..

        The energy point :math:`\omega` in the same unit of `eig`.
        Where `ne` is the number of energy points.

Return

    .. code-block:: python

        wght = numpy.empty([ng0, ng1, ng2, nb, nb, ne], dtype=numpy.complex)
    ..

      The integration weights.
