Installation
============

Files in this package
---------------------

``doc/manual_en.html`` : Manual in English (This file)

``doc/manual_jp.html`` : Manual in Japanese

``examples/`` : Sample programs using ``libtetrabz``

``src/`` : source files

``Makefile``, ``make.sys`` : Make file and the configuration file

Prerequisite
------------

``libtetrabz`` requires the following

-  fortran compiler

-  MPI library (If you use MPI/Hybrid version)

Installation guide
------------------

#. Download ``.tar.gz`` file from following web page.

   http://osdn.jp/projects/libtetrabz/releases/
               
#. Uncompress ``.tar.gz`` file and enter the generated directory.

   .. code-block:: bash

      $ tar xzvf xzvf libtetrabz_1.0.1.tar.gz
      $ cd libtetrabz
               

#. Modify ``make.sys`` for your platform. Adjust the following
   variables:

   ``TOPDIR`` : Path to the ``libtetrabz`` top directory

   ``F90`` : Fortran compile command (for serial version) (gfortran,
   ifort, frt, etc.)

   ``MPIF90`` : MPI-fortran compile command (mpif90, mpiifort, mpifrt,
   etc.)

   ``FFLAGS`` : fortran compile option

#. Compile.

   .. code-block:: bash

      $ make 

   If it succeed, you obtain the following files in ``src/``

   ::

       libtetrabz.a
       libtetrabz.mod
       libtetrabz_mpi.a
       libtetrabz_mpi.mod
               

   \* If you do not need the MPI version, ``make`` serial version only.

   .. code-block:: bash

      $ cd src
      $ make libtetrabz.a
               
Sample programs using ``libtetrabz`` are also compiled in ``example/`` .

``example/dos.x`` : Compute DOS of a tight-binding model in the cubic
lattice. The source code is ``dos.f90``

   .. figure:: ../dos.png
               :scale: 50

               Density of states of the tight-binding model on the
               cubic lattice calculated by using ``dos.x``.
               The solid line indicates the
               result converged about the number of :math:`k`.
               " :math:`+` " and " :math:`\times` " indicate
               results by using the linear tetrahedron method and the optimized
               tetrahedron method,
               respectively with :math:`8\times8\times8 k` grid.

``example/lindhard.x`` : Compute the Lindhard function. The source code
is ``lindhard.f90``

   .. figure:: ../lindhard.png
               :scale: 50

               (solid line) The analytical result of the Lindhard
               function. " :math:`+` " and " :math:`\times` " indicate results by using the linear
               tetrahedron method and the optimized tetrahedron method, respectively
               with :math:`8\times8\times8 k` grid.
