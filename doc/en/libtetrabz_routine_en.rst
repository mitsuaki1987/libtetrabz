Subroutines
===========

You can call a subroutine in this library as follows:

.. code-block:: fortran

   use libtetrabz, only : libtetrabz_occ
   
   call libtetrabz_occ(ltetra,bvec,nb,nge,eig,ngw,wght)
        
Every subroutine has a name starts from ``libtetrabz_``.
MPI version has a name starts from ``libtetrabz_mpi_``
and it requires ``use libtetrabz_mpi``.
The difference of arguments between the serial
version and the MPI version is an integer input argument ``comm``
which specifies the communicator.

For the C program, it can be used as follows:

.. code-block:: fortran

   #include "libtetrabz.h"

   libtetrabz_mp_libtetrabz_occ_(&ltetra,bvec,&nb,nge,eig,ngw,wght)
        
The name of a function in C becomes
``libtetrabz_mp_`` + fortran subroutine name + ``_``.
For the MPI version, ``libtetrabz_mpi.h``
should be included, and the name of a function becomes
``libtetrabz_mpi_mp_`` + fortran subroutine name + ``_``.
Variables should be passed as pointers.
Arrays should be declared as one dimensional arrays.

Total energy, charge density, occupations
-----------------------------------------

.. math::

   \begin{align}
   \sum_{n k} \theta(\varepsilon_{\rm F} -
   \varepsilon_{n k}) X_{n k}
   \end{align}

.. code-block:: fortran

    call libtetrabz_occ(ltetra,bvec,nb,nge,eig,ngw,wght)
    call libtetrabz_mpi_occ(ltetra,comm,bvec,nb,nge,eig,ngw,wght)

Parameters

   .. code-block:: fortran
                   
      integer,intent(in) :: ltetra
   ..
   
      Specify the type of the tetrahedron method.
      1 :math:`\cdots` the linear tetrahedron method.
      2 :math:`\cdots` the optimized tetrahedron method :ref:`[1] <ref>`.

   .. code-block:: fortran
                   
      integer,intent(in) :: comm
   ..
   
      Only for MPI version.
      Specify the communicator.

   .. code-block:: fortran
                   
      real(8),intent(in) :: bvec(3,3)
   ..
   
      Reciprocal lattice vectors (arbitrary unit). 
      Because they are used to choose the direction of tetrahedra,
      only their ratio is used.

   .. code-block:: fortran
                   
      integer,intent(in) :: nb
   ..
   
      The number of bands.

   .. code-block:: fortran
                   
      integer,intent(in) :: nge(3)
   ..
   
      Specify the :math:`k`\ -grid
      for input orbital energy.

   .. code-block:: fortran
                   
      real(8),intent(in) :: eig(nb,nge(1),nge(2),nge(3))
   ..
   
      The orbital energy measured from the Fermi energy
      ( :math:`\varepsilon_{\rm F} = 0` ).

   .. code-block:: fortran
                   
      integer,intent(in) :: ngw(3)
   ..
   
      Specify the :math:`k`\ -grid for output integration weights.
      You can make ``ngw`` :math:`\neq` ``nge`` (See :ref:`app`).

   .. code-block:: fortran
                   
      real(8),intent(out) :: wght(nb,ngw(1),ngw(2),ngw(3))
   ..
   
      The integration weights.
      
Fermi energy and occupations
----------------------------

.. math::

   \begin{align}
   \sum_{n k} \theta(\varepsilon_{\rm F} -
   \varepsilon_{n k}) X_{n k} 
   \end{align}

.. code-block:: fortran

    call libtetrabz_fermieng(ltetra,bvec,nb,nge,eig,ngw,wght,ef,nelec)
    call libtetrabz_mpi_fermieng(ltetra,comm,bvec,nb,nge,eig,ngw,wght,ef,nelec)
        
Parameters

   .. code-block:: fortran
                   
      integer,intent(in) :: ltetra
   ..
   
      Specify the type of the tetrahedron method.
      1 :math:`\cdots` the linear tetrahedron method.
      2 :math:`\cdots` the optimized tetrahedron method :ref:`[1] <ref>`.

   .. code-block:: fortran
                   
      integer,intent(in) :: comm
   ..
   
      Only for MPI version.
      Specify the communicator.

   .. code-block:: fortran
                   
      real(8),intent(in) :: bvec(3,3)
   ..
   
      Reciprocal lattice vectors (arbitrary unit). 
      Because they are used to choose the direction of tetrahedra,
      only their ratio is used.

   .. code-block:: fortran
                   
      integer,intent(in) :: nb
   ..
   
      The number of bands.

   .. code-block:: fortran
                   
      integer,intent(in) :: nge(3)
   ..
   
      Specify the :math:`k`\ -grid
      for input orbital energy.

   .. code-block:: fortran
                   
      real(8),intent(in) :: eig(nb,nge(1),nge(2),nge(3))
   ..
   
      The orbital energy measured from the Fermi energy
      ( :math:`\varepsilon_{\rm F} = 0` ).

   .. code-block:: fortran
                   
      integer,intent(in) :: ngw(3)
   ..
   
      Specify the :math:`k`\ -grid for output integration weights.
      You can make ``ngw`` :math:`\neq` ``nge`` (See :ref:`app`).

   .. code-block:: fortran
                   
      real(8),intent(out) :: wght(nb,ngw(1),ngw(2),ngw(3))
   ..
   
      The integration weights.

   .. code-block:: fortran
                         
      real(8),intent(out) :: ef
   ..

      The Fermi energy.

   .. code-block:: fortran
                         
      real(8),intent(in) :: nelec
   ..

      The number of (valence) electrons per spin.

Partial density of states
-------------------------

.. math::

   \begin{align}
   \sum_{n k} \delta(\omega - \varepsilon_{n k})
   X_{n k}(\omega) 
   \end{align}

.. code-block:: fortran

   call libtetrabz_dos(ltetra,bvec,nb,nge,eig,ngw,wght,ne,e0)
   call libtetrabz_mpi_dos(ltetra,comm,bvec,nb,nge,eig,ngw,wght,ne,e0)
        
Parameters

   .. code-block:: fortran
                   
      integer,intent(in) :: ltetra
   ..
   
      Specify the type of the tetrahedron method.
      1 :math:`\cdots` the linear tetrahedron method.
      2 :math:`\cdots` the optimized tetrahedron method :ref:`[1] <ref>`.

   .. code-block:: fortran
                   
      integer,intent(in) :: comm
   ..
   
      Only for MPI version.
      Specify the communicator.

   .. code-block:: fortran
                   
      real(8),intent(in) :: bvec(3,3)
   ..
   
      Reciprocal lattice vectors (arbitrary unit). 
      Because they are used to choose the direction of tetrahedra,
      only their ratio is used.

   .. code-block:: fortran
                   
      integer,intent(in) :: nb
   ..
   
      The number of bands.

   .. code-block:: fortran
                   
      integer,intent(in) :: nge(3)
   ..
   
      Specify the :math:`k`\ -grid
      for input orbital energy.

   .. code-block:: fortran
                   
      real(8),intent(in) :: eig(nb,nge(1),nge(2),nge(3))
   ..
   
      The orbital energy measured from the Fermi energy
      ( :math:`\varepsilon_{\rm F} = 0` ).

   .. code-block:: fortran
                   
      integer,intent(in) :: ngw(3)
   ..
   
      Specify the :math:`k`\ -grid for output integration weights.
      You can make ``ngw`` :math:`\neq` ``nge`` (See :ref:`app`).

   .. code-block:: fortran
                   
      real(8),intent(out) :: wght(ne,nb,ngw(1),ngw(2),ngw(3))
   ..
   
      The integration weights.

   .. code-block:: fortran
                         
      integer,intent(in) :: ne
   ..
   
      The number of energy where DOS is calculated.

   .. code-block:: fortran
                         
      real(8),intent(in) :: e0(ne)
   ..

      Energies where DOS is calculated.

Nesting function and Fr&oumlhlich parameter
-------------------------------------------

.. math::

   \begin{align}
   \sum_{n n' k} \delta(\varepsilon_{\rm F} -
   \varepsilon_{n k}) \delta(\varepsilon_{\rm F} - \varepsilon'_{n' k})
   X_{n n' k}
   \end{align}

.. code-block:: fortran

    call libtetrabz_doubledelta(ltetra,bvec,nb,nge,eig1,eig2,ngw,wght)
    call libtetrabz_mpi_doubledelta(ltetra,comm,bvec,nb,nge,eig1,eig2,ngw,wght)
        
Parameters

   .. code-block:: fortran
                   
      integer,intent(in) :: ltetra
   ..
   
      Specify the type of the tetrahedron method.
      1 :math:`\cdots` the linear tetrahedron method.
      2 :math:`\cdots` the optimized tetrahedron method :ref:`[1] <ref>`.

   .. code-block:: fortran
                   
      integer,intent(in) :: comm
   ..
   
      Only for MPI version.
      Specify the communicator.

   .. code-block:: fortran
                   
      real(8),intent(in) :: bvec(3,3)
   ..
   
      Reciprocal lattice vectors (arbitrary unit). 
      Because they are used to choose the direction of tetrahedra,
      only their ratio is used.

   .. code-block:: fortran
                   
      integer,intent(in) :: nb
   ..
   
      The number of bands.

   .. code-block:: fortran
                   
      integer,intent(in) :: nge(3)
   ..
   
      Specify the :math:`k`\ -grid
      for input orbital energy.

   .. code-block:: fortran
                   
      real(8),intent(in) :: eig1(nb,nge(1),nge(2),nge(3))
   ..
   
      The orbital energy measured from the Fermi energy
      ( :math:`\varepsilon_{\rm F} = 0` ).
      Do the same with ``eig2``.

   .. code-block:: fortran
                   
      real(8),intent(in) :: eig2(nb,nge(1),nge(2),nge(3))
   ..

      Another orbital energy.
      E.g. :math:`\varepsilon_{k + q}` on a shifted grid.

   .. code-block:: fortran
                   
      integer,intent(in) :: ngw(3)
   ..
   
      Specify the :math:`k`\ -grid for output integration weights.
      You can make ``ngw`` :math:`\neq` ``nge`` (See :ref:`app`).

   .. code-block:: fortran
                   
      real(8),intent(out) :: wght(nb,nb,ngw(1),ngw(2),ngw(3))
   ..
   
      The integration weights.

A part of DFPT calculation
--------------------------

.. math::

   \begin{align}
   \sum_{n n' k} \theta(\varepsilon_{\rm F} -
   \varepsilon_{n k}) \theta(\varepsilon_{n k} - \varepsilon'_{n' k}) 
   X_{n n' k}
   \end{align}

.. code-block:: fortran

    call libtetrabz_occstep(ltetra,bvec,nb,nge,eig1,eig2,ngw,wght)
    call libtetrabz_mpi_occstep(ltetra,comm,bvec,nb,nge,eig1,eig2,ngw,wght)
        
Parameters

   .. code-block:: fortran
                   
      integer,intent(in) :: ltetra
   ..
   
      Specify the type of the tetrahedron method.
      1 :math:`\cdots` the linear tetrahedron method.
      2 :math:`\cdots` the optimized tetrahedron method :ref:`[1] <ref>`.

   .. code-block:: fortran
                   
      integer,intent(in) :: comm
   ..
   
      Only for MPI version.
      Specify the communicator.

   .. code-block:: fortran
                   
      real(8),intent(in) :: bvec(3,3)
   ..
   
      Reciprocal lattice vectors (arbitrary unit). 
      Because they are used to choose the direction of tetrahedra,
      only their ratio is used.

   .. code-block:: fortran
                   
      integer,intent(in) :: nb
   ..
   
      The number of bands.

   .. code-block:: fortran
                   
      integer,intent(in) :: nge(3)
   ..
   
      Specify the :math:`k`\ -grid
      for input orbital energy.

   .. code-block:: fortran
                   
      real(8),intent(in) :: eig1(nb,nge(1),nge(2),nge(3))
   ..
   
      The orbital energy measured from the Fermi energy
      ( :math:`\varepsilon_{\rm F} = 0` ).
      Do the same with ``eig2``.

   .. code-block:: fortran
                   
      real(8),intent(in) :: eig2(nb,nge(1),nge(2),nge(3))
   ..

      Another orbital energy.
      E.g. :math:`\varepsilon_{k + q}` on a shifted grid.

   .. code-block:: fortran
                   
      integer,intent(in) :: ngw(3)
   ..
   
      Specify the :math:`k`\ -grid for output integration weights.
      You can make ``ngw`` :math:`\neq` ``nge`` (See :ref:`app`).

   .. code-block:: fortran
                   
      real(8),intent(out) :: wght(nb,nb,ngw(1),ngw(2),ngw(3))
   ..
   
      The integration weights.

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

    call libtetrabz_polstat(ltetra,bvec,nb,nge,eig1,eig2,ngw,wght)
    call libtetrabz_mpi_occstep(ltetra,comm,bvec,nb,nge,eig1,eig2,ngw,wght)
        
Parameters

   .. code-block:: fortran
                   
      integer,intent(in) :: ltetra
   ..
   
      Specify the type of the tetrahedron method.
      1 :math:`\cdots` the linear tetrahedron method.
      2 :math:`\cdots` the optimized tetrahedron method :ref:`[1] <ref>`.

   .. code-block:: fortran
                   
      integer,intent(in) :: comm
   ..
   
      Only for MPI version.
      Specify the communicator.

   .. code-block:: fortran
                   
      real(8),intent(in) :: bvec(3,3)
   ..
   
      Reciprocal lattice vectors (arbitrary unit). 
      Because they are used to choose the direction of tetrahedra,
      only their ratio is used.

   .. code-block:: fortran
                   
      integer,intent(in) :: nb
   ..
   
      The number of bands.

   .. code-block:: fortran
                   
      integer,intent(in) :: nge(3)
   ..
   
      Specify the :math:`k`\ -grid
      for input orbital energy.

   .. code-block:: fortran
                   
      real(8),intent(in) :: eig1(nb,nge(1),nge(2),nge(3))
   ..
   
      The orbital energy measured from the Fermi energy
      ( :math:`\varepsilon_{\rm F} = 0` ).
      Do the same with ``eig2``.

   .. code-block:: fortran
                   
      real(8),intent(in) :: eig2(nb,nge(1),nge(2),nge(3))
   ..

      Another orbital energy.
      E.g. :math:`\varepsilon_{k + q}` on a shifted grid.

   .. code-block:: fortran
                   
      integer,intent(in) :: ngw(3)
   ..
   
      Specify the :math:`k`\ -grid for output integration weights.
      You can make ``ngw`` :math:`\neq` ``nge`` (See :ref:`app`).

   .. code-block:: fortran
                   
      real(8),intent(out) :: wght(nb,nb,ngw(1),ngw(2),ngw(3))
   ..
   
      The integration weights.

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

    call libtetrabz_fermigr(ltetra,bvec,nb,nge,eig1,eig2,ngw,wght,ne,e0)
    call libtetrabz_mpi_fermigr(ltetra,comm,bvec,nb,nge,eig1,eig2,ngw,wght,ne,e0)
        
Parameters

   .. code-block:: fortran
                   
      integer,intent(in) :: ltetra
   ..
   
      Specify the type of the tetrahedron method.
      1 :math:`\cdots` the linear tetrahedron method.
      2 :math:`\cdots` the optimized tetrahedron method :ref:`[1] <ref>`.

   .. code-block:: fortran
                   
      integer,intent(in) :: comm
   ..
   
      Only for MPI version.
      Specify the communicator.

   .. code-block:: fortran
                   
      real(8),intent(in) :: bvec(3,3)
   ..
   
      Reciprocal lattice vectors (arbitrary unit). 
      Because they are used to choose the direction of tetrahedra,
      only their ratio is used.

   .. code-block:: fortran
                   
      integer,intent(in) :: nb
   ..
   
      The number of bands.

   .. code-block:: fortran
                   
      integer,intent(in) :: nge(3)
   ..
   
      Specify the :math:`k`\ -grid
      for input orbital energy.

   .. code-block:: fortran
                   
      real(8),intent(in) :: eig1(nb,nge(1),nge(2),nge(3))
   ..
   
      The orbital energy measured from the Fermi energy
      ( :math:`\varepsilon_{\rm F} = 0` ).
      Do the same with ``eig2``.

   .. code-block:: fortran
                   
      real(8),intent(in) :: eig2(nb,nge(1),nge(2),nge(3))
   ..

      Another orbital energy.
      E.g. :math:`\varepsilon_{k + q}` on a shifted grid.

   .. code-block:: fortran
                   
      integer,intent(in) :: ngw(3)
   ..
   
      Specify the :math:`k`\ -grid for output integration weights.
      You can make ``ngw`` :math:`\neq` ``nge`` (See :ref:`app`).

   .. code-block:: fortran
                   
      real(8),intent(out) :: wght(ne,nb,nb,ngw(1),ngw(2),ngw(3))
   ..
   
      The integration weights.

   .. code-block:: fortran
                         
      integer,intent(in) :: ne
   ..
   
      The number of branches of the phonon.

   .. code-block:: fortran
                         
      real(8),intent(in) :: e0(ne)
   ..
   
      Phonon frequencies.

Polarization function (imaginary frequency)
-------------------------------------------

.. math::

   \begin{align}
   \sum_{n n' k} \frac{\theta(\varepsilon_{\rm F} - \varepsilon_{n k})
   \theta(\varepsilon'_{n' k} - \varepsilon_{\rm F})}
   {\varepsilon'_{n' k} - \varepsilon_{n k} + i \omega}
   X_{n n' k}(\omega) 
   \end{align}

.. code-block:: fortran

    call libtetrabz_polimg(ltetra,bvec,nb,nge,eig1,eig2,ngw,wght,ne,e0)
    call libtetrabz_mpi_polimg(ltetra,comm,bvec,nb,nge,eig1,eig2,ngw,wght,ne,e0)
        
Parameters

   .. code-block:: fortran
                   
      integer,intent(in) :: ltetra
   ..
   
      Specify the type of the tetrahedron method.
      1 :math:`\cdots` the linear tetrahedron method.
      2 :math:`\cdots` the optimized tetrahedron method :ref:`[1] <ref>`.

   .. code-block:: fortran
                   
      integer,intent(in) :: comm
   ..
   
      Only for MPI version.
      Specify the communicator.

   .. code-block:: fortran
                   
      real(8),intent(in) :: bvec(3,3)
   ..
   
      Reciprocal lattice vectors (arbitrary unit). 
      Because they are used to choose the direction of tetrahedra,
      only their ratio is used.

   .. code-block:: fortran
                   
      integer,intent(in) :: nb
   ..
   
      The number of bands.

   .. code-block:: fortran
                   
      integer,intent(in) :: nge(3)
   ..
   
      Specify the :math:`k`\ -grid
      for input orbital energy.

   .. code-block:: fortran
                   
      real(8),intent(in) :: eig1(nb,nge(1),nge(2),nge(3))
   ..
   
      The orbital energy measured from the Fermi energy
      ( :math:`\varepsilon_{\rm F} = 0` ).
      Do the same with ``eig2``.

   .. code-block:: fortran
                   
      real(8),intent(in) :: eig2(nb,nge(1),nge(2),nge(3))
   ..

      Another orbital energy.
      E.g. :math:`\varepsilon_{k + q}` on a shifted grid.

   .. code-block:: fortran
                   
      integer,intent(in) :: ngw(3)
   ..
   
      Specify the :math:`k`\ -grid for output integration weights.
      You can make ``ngw`` :math:`\neq` ``nge`` (See :ref:`app`).

   .. code-block:: fortran
                   
      real(8),intent(out) :: wght(2,ne,nb,nb,ngw(1),ngw(2),ngw(3))
   ..
   
      The integration weights.

   .. code-block:: fortran
                         
      integer,intent(in) :: ne
   ..
   
      The number of imaginary frequencies where
      polarization functions are calculated.

   .. code-block:: fortran
                         
      real(8),intent(in) :: e0(ne)
   ..
   
      Imaginary frequencies where
      polarization functions are calculated.
