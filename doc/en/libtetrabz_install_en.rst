Installation
============

Important files and directories
-------------------------------

- ``doc/`` : Directory for manuals
   - ``doc/index.html`` : Index page 
- ``src/`` : Directory for the sources of the library
- ``example/`` : Directory for the sample program
- ``test/`` : Directory for tests
- ``configure`` : Configuration script for the build

Prerequisite
------------

``libtetrabz`` requires the following

-  fortran and C compiler

-  MPI library (If you use MPI/Hybrid version)

Installation guide
------------------

#. Download ``.tar.gz`` file from following web page.

   http://osdn.jp/projects/libtetrabz/releases/
               
#. Uncompress ``.tar.gz`` file and enter the generated directory.

   .. code-block:: bash

      $ tar xzvf xzvf libtetrabz_1.0.1.tar.gz
      $ cd libtetrabz

#. Configure the build environment.
   
   .. code-block:: bash

      $ ./configure --prefix=install_dir

   Then, this script checks the compiler and the libraries required for the installation,
   and creates Makefiles.
   ``install_dir`` indicates the full path of the directory where the library is installed
   (you should replace it according to your case).
   If none is specified, ``/use/local/`` is chosen for storing libraries
   by ``make install``  (Therefore, if one is not the admin, ``install_dir`` must be specified to
   the different directory).
   ``configure`` has many options, and they are used according to the environment etc.
   For more details, please see :ref:`configoption`.

#. After ``configure`` finishes successfully and Makefiles are generated,
   please type

   .. code-block:: bash

      $ make

   to build libraries. Then please type

   .. code-block:: bash

      $ make install

   to store libraries and the sample program to ``install_dir/lib`` and ``install_dir/bin``, respectively.
   Although one can use libraries and the sample program without ``make install``,
   they are a little different to the installed one.

#. Add the libtetrabz's library directory (``install_dir/lib``) to the
   search path of the dynamically linked program (environment variable ``LD_LIBRARY_PATH``).

   .. code-block:: bash

      $ export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:install_dir/lib
               
#. Sample programs using ``libtetrabz`` are also compiled in ``example/`` .

   ``example/dos.x`` : Compute DOS of a tight-binding model in the cubic
   lattice. The source code is ``dos.f90``

      .. figure:: ../figs/dos.png
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

      .. figure:: ../figs/lindhard.png
         :scale: 50

         (solid line) The analytical result of the Lindhard
         function. " :math:`+` " and " :math:`\times` " indicate results by using the linear
         tetrahedron method and the optimized tetrahedron method, respectively
         with :math:`8\times8\times8 k` grid.

.. _configoption:

Options for configure
---------------------

``configure`` has many options and environment variables.
They can be specified at once. E.g.

.. code-block:: bash

  $ ./configure --prefix=/home/komega/ --with-mpi=yes FC=mpif90

All options and variables have default values.
We show a part of them as follows:

``---prefix``

   Default: ``---prefix=/usr/local/``.
   Specify the directory where the library etc. are installed.

``--with-mpi``

   Default: ``--with-mpi=no`` (without MPI).
   Whether use MPI (``--with-mpi=yes``), or not.

``--with-openmp``

   Default: ``--with-openmp=yes`` (with OpenMP).
   Whether use OpenMP or not (``--with-openmp=no``).

``--enable-shared``

   Default: ``--enable-shared``.
   Whether generate shared library.

``--enable-static``

   Default: ``--enable-static``.
   Whether generate static library.

``FC``, ``C``

   Default: The fortran/C compiler chosen automatically from those in the system.
   When ``--with-mpi`` is specified, the corresponding MPI compiler
   (such as ``mpif90``) is searched.
   If ``FC`` printed the end of the standard-output of ``configure`` is not
   what you want, please set it manually as ``./configure FC=gfortran``.

``--help``

   Display all options including above, and stop without configuration.
