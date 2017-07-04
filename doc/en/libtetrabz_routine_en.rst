Subroutines
===========

You can call a subroutine in this library as follows:

.. code-block:: fortran

   USE libtetrabz, ONLY : libtetrabz_occ
   
   CALL libtetrabz_occ(ltetra,bvec,nb,nge,eig,ngw,wght)
        
Every subroutine has a name starts from ``libtetrabz_``.

For the C program, it can be used as follows:

.. code-block:: c

   #include "libtetrabz.h"

   libtetrabz_occ(&ltetra,bvec,&nb,nge,eig,ngw,wght)
        
Variables should be passed as pointers.
Arrays should be declared as one dimensional arrays.
Also, the communicator argument for the routine should be
transformed from the C/C++'s one to the fortran's one as follows.

.. code-block:: c

      comm_f = MPI_Comm_c2f(comm_c);

Total energy, charge density, occupations
-----------------------------------------

.. math::

   \begin{align}
   \sum_{n k} \theta(\varepsilon_{\rm F} -
   \varepsilon_{n k}) X_{n k}
   \end{align}

.. code-block:: fortran

    CALL libtetrabz_occ(ltetra,bvec,nb,nge,eig,ngw,wght,comm)

Parameters

   .. code-block:: fortran
                   
      INTEGER,INTENT(IN) :: ltetra
   ..
   
      Specify the type of the tetrahedron method.
      1 :math:`\cdots` the linear tetrahedron method.
      2 :math:`\cdots` the optimized tetrahedron method :ref:`[1] <ref>`.

   .. code-block:: fortran
                   
      REAL(8),INTENT(IN) :: bvec(3,3)
   ..
   
      Reciprocal lattice vectors (arbitrary unit). 
      Because they are used to choose the direction of tetrahedra,
      only their ratio is used.

   .. code-block:: fortran
                   
      INTEGER,INTENT(IN) :: nb
   ..
   
      The number of bands.

   .. code-block:: fortran
                   
      INTEGER,INTENT(IN) :: nge(3)
   ..
   
      Specify the :math:`k`\ -grid
      for input orbital energy.

   .. code-block:: fortran
                   
      REAL(8),INTENT(IN) :: eig(nb,nge(1),nge(2),nge(3))
   ..
   
      The orbital energy measured from the Fermi energy
      ( :math:`\varepsilon_{\rm F} = 0` ).

   .. code-block:: fortran
                   
      INTEGER,INTENT(IN) :: ngw(3)
   ..
   
      Specify the :math:`k`\ -grid for output integration weights.
      You can make ``ngw`` :math:`\neq` ``nge`` (See :ref:`app`).

   .. code-block:: fortran
                   
      REAL(8),INTENT(OUT) :: wght(nb,ngw(1),ngw(2),ngw(3))
   ..
   
      The integration weights.
      
   .. code-block:: fortran
                   
      INTEGER,INTENT(IN),OPTIONAL :: comm
   ..

      Optional argument. Communicators for MPI such as ``MPI_COMM_WORLD``.
      Only for MPI / Hybrid parallelization.
      For C compiler without MPI, just pass ``NULL`` to omit this argment.

Fermi energy and occupations
----------------------------

.. math::

   \begin{align}
   \sum_{n k} \theta(\varepsilon_{\rm F} -
   \varepsilon_{n k}) X_{n k} 
   \end{align}

.. code-block:: fortran

    CALL libtetrabz_fermieng(ltetra,bvec,nb,nge,eig,ngw,wght,ef,nelec,comm)
        
Parameters

   .. code-block:: fortran
                   
      INTEGER,INTENT(IN) :: ltetra
   ..
   
      Specify the type of the tetrahedron method.
      1 :math:`\cdots` the linear tetrahedron method.
      2 :math:`\cdots` the optimized tetrahedron method :ref:`[1] <ref>`.

   .. code-block:: fortran
                   
      REAL(8),INTENT(IN) :: bvec(3,3)
   ..
   
      Reciprocal lattice vectors (arbitrary unit). 
      Because they are used to choose the direction of tetrahedra,
      only their ratio is used.

   .. code-block:: fortran
                   
      INTEGER,INTENT(IN) :: nb
   ..
   
      The number of bands.

   .. code-block:: fortran
                   
      INTEGER,INTENT(IN) :: nge(3)
   ..
   
      Specify the :math:`k`\ -grid
      for input orbital energy.

   .. code-block:: fortran
                   
      REAL(8),INTENT(IN) :: eig(nb,nge(1),nge(2),nge(3))
   ..
   
      The orbital energy measured from the Fermi energy
      ( :math:`\varepsilon_{\rm F} = 0` ).

   .. code-block:: fortran
                   
      INTEGER,INTENT(IN) :: ngw(3)
   ..
   
      Specify the :math:`k`\ -grid for output integration weights.
      You can make ``ngw`` :math:`\neq` ``nge`` (See :ref:`app`).

   .. code-block:: fortran
                   
      REAL(8),INTENT(OUT) :: wght(nb,ngw(1),ngw(2),ngw(3))
   ..
   
      The integration weights.

   .. code-block:: fortran
                         
      REAL(8),INTENT(OUT) :: ef
   ..

      The Fermi energy.

   .. code-block:: fortran
                         
      REAL(8),INTENT(IN) :: nelec
   ..

      The number of (valence) electrons per spin.

   .. code-block:: fortran
                   
      INTEGER,INTENT(IN),OPTIONAL :: comm
   ..

      Optional argument. Communicators for MPI such as ``MPI_COMM_WORLD``.
      Only for MPI / Hybrid parallelization.
      For C compiler without MPI, just pass ``NULL`` to omit this argment.

Partial density of states
-------------------------

.. math::

   \begin{align}
   \sum_{n k} \delta(\omega - \varepsilon_{n k})
   X_{n k}(\omega) 
   \end{align}

.. code-block:: fortran

   CALL libtetrabz_dos(ltetra,bvec,nb,nge,eig,ngw,wght,ne,e0,comm)
        
Parameters

   .. code-block:: fortran
                   
      INTEGER,INTENT(IN) :: ltetra
   ..
   
      Specify the type of the tetrahedron method.
      1 :math:`\cdots` the linear tetrahedron method.
      2 :math:`\cdots` the optimized tetrahedron method :ref:`[1] <ref>`.

   .. code-block:: fortran
                   
      REAL(8),INTENT(IN) :: bvec(3,3)
   ..
   
      Reciprocal lattice vectors (arbitrary unit). 
      Because they are used to choose the direction of tetrahedra,
      only their ratio is used.

   .. code-block:: fortran
                   
      INTEGER,INTENT(IN) :: nb
   ..
   
      The number of bands.

   .. code-block:: fortran
                   
      INTEGER,INTENT(IN) :: nge(3)
   ..
   
      Specify the :math:`k`\ -grid
      for input orbital energy.

   .. code-block:: fortran
                   
      REAL(8),INTENT(IN) :: eig(nb,nge(1),nge(2),nge(3))
   ..
   
      The orbital energy measured from the Fermi energy
      ( :math:`\varepsilon_{\rm F} = 0` ).

   .. code-block:: fortran
                   
      INTEGER,INTENT(IN) :: ngw(3)
   ..
   
      Specify the :math:`k`\ -grid for output integration weights.
      You can make ``ngw`` :math:`\neq` ``nge`` (See :ref:`app`).

   .. code-block:: fortran
                   
      REAL(8),INTENT(OUT) :: wght(ne,nb,ngw(1),ngw(2),ngw(3))
   ..
   
      The integration weights.

   .. code-block:: fortran
                         
      INTEGER,INTENT(IN) :: ne
   ..
   
      The number of energy where DOS is calculated.

   .. code-block:: fortran
                         
      REAL(8),INTENT(IN) :: e0(ne)
   ..

      Energies where DOS is calculated.

   .. code-block:: fortran
                   
      INTEGER,INTENT(IN),OPTIONAL :: comm
   ..

      Optional argument. Communicators for MPI such as ``MPI_COMM_WORLD``.
      Only for MPI / Hybrid parallelization.
      For C compiler without MPI, just pass ``NULL`` to omit this argment.

Integrated density of states
----------------------------

.. math::

   \begin{align}
   \sum_{n k} \theta(\omega - \varepsilon_{n k})
   X_{n k}(\omega) 
   \end{align}

.. code-block:: fortran

   CALL libtetrabz_intdos(ltetra,bvec,nb,nge,eig,ngw,wght,ne,e0,comm)
        
Parameters

   .. code-block:: fortran
                   
      INTEGER,INTENT(IN) :: ltetra
   ..
   
      Specify the type of the tetrahedron method.
      1 :math:`\cdots` the linear tetrahedron method.
      2 :math:`\cdots` the optimized tetrahedron method :ref:`[1] <ref>`.

   .. code-block:: fortran
                   
      REAL(8),INTENT(IN) :: bvec(3,3)
   ..
   
      Reciprocal lattice vectors (arbitrary unit). 
      Because they are used to choose the direction of tetrahedra,
      only their ratio is used.

   .. code-block:: fortran
                   
      INTEGER,INTENT(IN) :: nb
   ..
   
      The number of bands.

   .. code-block:: fortran
                   
      INTEGER,INTENT(IN) :: nge(3)
   ..
   
      Specify the :math:`k`\ -grid
      for input orbital energy.

   .. code-block:: fortran
                   
      REAL(8),INTENT(IN) :: eig(nb,nge(1),nge(2),nge(3))
   ..
   
      The orbital energy measured from the Fermi energy
      ( :math:`\varepsilon_{\rm F} = 0` ).

   .. code-block:: fortran
                   
      INTEGER,INTENT(IN) :: ngw(3)
   ..
   
      Specify the :math:`k`\ -grid for output integration weights.
      You can make ``ngw`` :math:`\neq` ``nge`` (See :ref:`app`).

   .. code-block:: fortran
                   
      REAL(8),INTENT(OUT) :: wght(ne,nb,ngw(1),ngw(2),ngw(3))
   ..
   
      The integration weights.

   .. code-block:: fortran
                         
      INTEGER,INTENT(IN) :: ne
   ..
   
      The number of energy where DOS is calculated.

   .. code-block:: fortran
                         
      REAL(8),INTENT(IN) :: e0(ne)
   ..

      Energies where DOS is calculated.

   .. code-block:: fortran
                   
      INTEGER,INTENT(IN),OPTIONAL :: comm
   ..

      Optional argument. Communicators for MPI such as ``MPI_COMM_WORLD``.
      Only for MPI / Hybrid parallelization.
      For C compiler without MPI, just pass ``NULL`` to omit this argment.

Nesting function and Fr&oumlhlich parameter
-------------------------------------------

.. math::

   \begin{align}
   \sum_{n n' k} \delta(\varepsilon_{\rm F} -
   \varepsilon_{n k}) \delta(\varepsilon_{\rm F} - \varepsilon'_{n' k})
   X_{n n' k}
   \end{align}

.. code-block:: fortran

    CALL libtetrabz_dbldelta(ltetra,bvec,nb,nge,eig1,eig2,ngw,wght,comm)
        
Parameters

   .. code-block:: fortran
                   
      INTEGER,INTENT(IN) :: ltetra
   ..
   
      Specify the type of the tetrahedron method.
      1 :math:`\cdots` the linear tetrahedron method.
      2 :math:`\cdots` the optimized tetrahedron method :ref:`[1] <ref>`.

   .. code-block:: fortran
                   
      REAL(8),INTENT(IN) :: bvec(3,3)
   ..
   
      Reciprocal lattice vectors (arbitrary unit). 
      Because they are used to choose the direction of tetrahedra,
      only their ratio is used.

   .. code-block:: fortran
                   
      INTEGER,INTENT(IN) :: nb
   ..
   
      The number of bands.

   .. code-block:: fortran
                   
      INTEGER,INTENT(IN) :: nge(3)
   ..
   
      Specify the :math:`k`\ -grid
      for input orbital energy.

   .. code-block:: fortran
                   
      REAL(8),INTENT(IN) :: eig1(nb,nge(1),nge(2),nge(3))
   ..
   
      The orbital energy measured from the Fermi energy
      ( :math:`\varepsilon_{\rm F} = 0` ).
      Do the same with ``eig2``.

   .. code-block:: fortran
                   
      REAL(8),INTENT(IN) :: eig2(nb,nge(1),nge(2),nge(3))
   ..

      Another orbital energy.
      E.g. :math:`\varepsilon_{k + q}` on a shifted grid.

   .. code-block:: fortran
                   
      INTEGER,INTENT(IN) :: ngw(3)
   ..
   
      Specify the :math:`k`\ -grid for output integration weights.
      You can make ``ngw`` :math:`\neq` ``nge`` (See :ref:`app`).

   .. code-block:: fortran
                   
      REAL(8),INTENT(OUT) :: wght(nb,nb,ngw(1),ngw(2),ngw(3))
   ..
   
      The integration weights.

   .. code-block:: fortran
                   
      INTEGER,INTENT(IN),OPTIONAL :: comm
   ..

      Optional argument. Communicators for MPI such as ``MPI_COMM_WORLD``.
      Only for MPI / Hybrid parallelization.
      For C compiler without MPI, just pass ``NULL`` to omit this argment.

A part of DFPT calculation
--------------------------

.. math::

   \begin{align}
   \sum_{n n' k} \theta(\varepsilon_{\rm F} -
   \varepsilon_{n k}) \theta(\varepsilon_{n k} - \varepsilon'_{n' k}) 
   X_{n n' k}
   \end{align}

.. code-block:: fortran

    CALL libtetrabz_dblstep(ltetra,bvec,nb,nge,eig1,eig2,ngw,wght,comm)
        
Parameters

   .. code-block:: fortran
                   
      INTEGER,INTENT(IN) :: ltetra
   ..
   
      Specify the type of the tetrahedron method.
      1 :math:`\cdots` the linear tetrahedron method.
      2 :math:`\cdots` the optimized tetrahedron method :ref:`[1] <ref>`.

   .. code-block:: fortran
                   
      REAL(8),INTENT(IN) :: bvec(3,3)
   ..
   
      Reciprocal lattice vectors (arbitrary unit). 
      Because they are used to choose the direction of tetrahedra,
      only their ratio is used.

   .. code-block:: fortran
                   
      INTEGER,INTENT(IN) :: nb
   ..
   
      The number of bands.

   .. code-block:: fortran
                   
      INTEGER,INTENT(IN) :: nge(3)
   ..
   
      Specify the :math:`k`\ -grid
      for input orbital energy.

   .. code-block:: fortran
                   
      REAL(8),INTENT(IN) :: eig1(nb,nge(1),nge(2),nge(3))
   ..
   
      The orbital energy measured from the Fermi energy
      ( :math:`\varepsilon_{\rm F} = 0` ).
      Do the same with ``eig2``.

   .. code-block:: fortran
                   
      REAL(8),INTENT(IN) :: eig2(nb,nge(1),nge(2),nge(3))
   ..

      Another orbital energy.
      E.g. :math:`\varepsilon_{k + q}` on a shifted grid.

   .. code-block:: fortran
                   
      INTEGER,INTENT(IN) :: ngw(3)
   ..
   
      Specify the :math:`k`\ -grid for output integration weights.
      You can make ``ngw`` :math:`\neq` ``nge`` (See :ref:`app`).

   .. code-block:: fortran
                   
      REAL(8),INTENT(OUT) :: wght(nb,nb,ngw(1),ngw(2),ngw(3))
   ..
   
      The integration weights.

   .. code-block:: fortran
                   
      INTEGER,INTENT(IN),OPTIONAL :: comm
   ..

      Optional argument. Communicators for MPI such as ``MPI_COMM_WORLD``.
      Only for MPI / Hybrid parallelization.
      For C compiler without MPI, just pass ``NULL`` to omit this argment.

Static polarization function
----------------------------

.. math::

   \begin{align}
   \sum_{n n' k} \frac{\theta(\varepsilon_{\rm F} - \varepsilon_{n k})
   \theta(\varepsilon'_{n' k} - \varepsilon_{\rm F})}
   {\varepsilon'_{n' k} - \varepsilon_{n k}}
   X_{n n' k} 
   \end{align}

.. code-block:: fortran

    CALL libtetrabz_polstat(ltetra,bvec,nb,nge,eig1,eig2,ngw,wght,comm)
        
Parameters

   .. code-block:: fortran
                   
      INTEGER,INTENT(IN) :: ltetra
   ..
   
      Specify the type of the tetrahedron method.
      1 :math:`\cdots` the linear tetrahedron method.
      2 :math:`\cdots` the optimized tetrahedron method :ref:`[1] <ref>`.

   .. code-block:: fortran
                   
      REAL(8),INTENT(IN) :: bvec(3,3)
   ..
   
      Reciprocal lattice vectors (arbitrary unit). 
      Because they are used to choose the direction of tetrahedra,
      only their ratio is used.

   .. code-block:: fortran
                   
      INTEGER,INTENT(IN) :: nb
   ..
   
      The number of bands.

   .. code-block:: fortran
                   
      INTEGER,INTENT(IN) :: nge(3)
   ..
   
      Specify the :math:`k`\ -grid
      for input orbital energy.

   .. code-block:: fortran
                   
      REAL(8),INTENT(IN) :: eig1(nb,nge(1),nge(2),nge(3))
   ..
   
      The orbital energy measured from the Fermi energy
      ( :math:`\varepsilon_{\rm F} = 0` ).
      Do the same with ``eig2``.

   .. code-block:: fortran
                   
      REAL(8),INTENT(IN) :: eig2(nb,nge(1),nge(2),nge(3))
   ..

      Another orbital energy.
      E.g. :math:`\varepsilon_{k + q}` on a shifted grid.

   .. code-block:: fortran
                   
      INTEGER,INTENT(IN) :: ngw(3)
   ..
   
      Specify the :math:`k`\ -grid for output integration weights.
      You can make ``ngw`` :math:`\neq` ``nge`` (See :ref:`app`).

   .. code-block:: fortran
                   
      REAL(8),INTENT(OUT) :: wght(nb,nb,ngw(1),ngw(2),ngw(3))
   ..
   
      The integration weights.

   .. code-block:: fortran
                   
      INTEGER,INTENT(IN),OPTIONAL :: comm
   ..

      Optional argument. Communicators for MPI such as ``MPI_COMM_WORLD``.
      Only for MPI / Hybrid parallelization.
      For C compiler without MPI, just pass ``NULL`` to omit this argment.

Phonon linewidth
----------------

.. math::

   \begin{align}
   \sum_{n n' k} \theta(\varepsilon_{\rm F} -
   \varepsilon_{n k}) \theta(\varepsilon'_{n' k} - \varepsilon_{\rm F})
   \delta(\varepsilon'_{n' k} - \varepsilon_{n k} - \omega)
   X_{n n' k}(\omega) 
   \end{align}

.. code-block:: fortran

    CALL libtetrabz_fermigr(ltetra,bvec,nb,nge,eig1,eig2,ngw,wght,ne,e0,comm)
        
Parameters

   .. code-block:: fortran
                   
      INTEGER,INTENT(IN) :: ltetra
   ..
   
      Specify the type of the tetrahedron method.
      1 :math:`\cdots` the linear tetrahedron method.
      2 :math:`\cdots` the optimized tetrahedron method :ref:`[1] <ref>`.

   .. code-block:: fortran
                   
      REAL(8),INTENT(IN) :: bvec(3,3)
   ..
   
      Reciprocal lattice vectors (arbitrary unit). 
      Because they are used to choose the direction of tetrahedra,
      only their ratio is used.

   .. code-block:: fortran
                   
      INTEGER,INTENT(IN) :: nb
   ..
   
      The number of bands.

   .. code-block:: fortran
                   
      INTEGER,INTENT(IN) :: nge(3)
   ..
   
      Specify the :math:`k`\ -grid
      for input orbital energy.

   .. code-block:: fortran
                   
      REAL(8),INTENT(IN) :: eig1(nb,nge(1),nge(2),nge(3))
   ..
   
      The orbital energy measured from the Fermi energy
      ( :math:`\varepsilon_{\rm F} = 0` ).
      Do the same with ``eig2``.

   .. code-block:: fortran
                   
      REAL(8),INTENT(IN) :: eig2(nb,nge(1),nge(2),nge(3))
   ..

      Another orbital energy.
      E.g. :math:`\varepsilon_{k + q}` on a shifted grid.

   .. code-block:: fortran
                   
      INTEGER,INTENT(IN) :: ngw(3)
   ..
   
      Specify the :math:`k`\ -grid for output integration weights.
      You can make ``ngw`` :math:`\neq` ``nge`` (See :ref:`app`).

   .. code-block:: fortran
                   
      REAL(8),INTENT(OUT) :: wght(ne,nb,nb,ngw(1),ngw(2),ngw(3))
   ..
   
      The integration weights.

   .. code-block:: fortran
                         
      INTEGER,INTENT(IN) :: ne
   ..
   
      The number of branches of the phonon.

   .. code-block:: fortran
                         
      REAL(8),INTENT(IN) :: e0(ne)
   ..
   
      Phonon frequencies.

   .. code-block:: fortran
                   
      INTEGER,INTENT(IN),OPTIONAL :: comm
   ..

      Optional argument. Communicators for MPI such as ``MPI_COMM_WORLD``.
      Only for MPI / Hybrid parallelization.
      For C compiler without MPI, just pass ``NULL`` to omit this argment.

Polarization function (complex frequency)
-----------------------------------------

.. math::

   \begin{align}
   \sum_{n n' k} \frac{\theta(\varepsilon_{\rm F} - \varepsilon_{n k})
   \theta(\varepsilon'_{n' k} - \varepsilon_{\rm F})}
   {\varepsilon'_{n' k} - \varepsilon_{n k} + i \omega}
   X_{n n' k}(\omega) 
   \end{align}

.. code-block:: fortran

    CALL libtetrabz_polcmplx(ltetra,bvec,nb,nge,eig1,eig2,ngw,wght,ne,e0,comm)
        
Parameters

   .. code-block:: fortran
                   
      INTEGER,INTENT(IN) :: ltetra
   ..
   
      Specify the type of the tetrahedron method.
      1 :math:`\cdots` the linear tetrahedron method.
      2 :math:`\cdots` the optimized tetrahedron method :ref:`[1] <ref>`.

   .. code-block:: fortran
                   
      REAL(8),INTENT(IN) :: bvec(3,3)
   ..
   
      Reciprocal lattice vectors (arbitrary unit). 
      Because they are used to choose the direction of tetrahedra,
      only their ratio is used.

   .. code-block:: fortran
                   
      INTEGER,INTENT(IN) :: nb
   ..
   
      The number of bands.

   .. code-block:: fortran
                   
      INTEGER,INTENT(IN) :: nge(3)
   ..
   
      Specify the :math:`k`\ -grid
      for input orbital energy.

   .. code-block:: fortran
                   
      REAL(8),INTENT(IN) :: eig1(nb,nge(1),nge(2),nge(3))
   ..
   
      The orbital energy measured from the Fermi energy
      ( :math:`\varepsilon_{\rm F} = 0` ).
      Do the same with ``eig2``.

   .. code-block:: fortran
                   
      REAL(8),INTENT(IN) :: eig2(nb,nge(1),nge(2),nge(3))
   ..

      Another orbital energy.
      E.g. :math:`\varepsilon_{k + q}` on a shifted grid.

   .. code-block:: fortran
                   
      INTEGER,INTENT(IN) :: ngw(3)
   ..
   
      Specify the :math:`k`\ -grid for output integration weights.
      You can make ``ngw`` :math:`\neq` ``nge`` (See :ref:`app`).

   .. code-block:: fortran
                   
      COMPLEX(8),INTENT(OUT) :: wght(ne,nb,nb,ngw(1),ngw(2),ngw(3))
   ..
   
      The integration weights.

   .. code-block:: fortran
                         
      INTEGER,INTENT(IN) :: ne
   ..
   
      The number of imaginary frequencies where
      polarization functions are calculated.

   .. code-block:: fortran
                         
      COMPLEX(8),INTENT(IN) :: e0(ne)
   ..
   
      Complex frequencies where
      polarization functions are calculated.

   .. code-block:: fortran
                   
      INTEGER,INTENT(IN),OPTIONAL :: comm
   ..

      Optional argument. Communicators for MPI such as ``MPI_COMM_WORLD``.
      Only for MPI / Hybrid parallelization.
      For C compiler without MPI, just pass ``NULL`` to omit this argment.

