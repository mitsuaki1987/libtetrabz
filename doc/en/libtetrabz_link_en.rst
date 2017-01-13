Linking libtetrabz
==================

**e. g. / For intel fortran**

   .. code-block:: bash

      $ ifort program.f90 -L {libtetrabz path}/src/ -I {libtetrabz path}/src -ltetrabz -fopenmp
      $ mpiifort program.f90 -L {libtetrabz path}/src/  -I {libtetrabz path}/src -ltetrabz_mpi -fopenmp
          

**e. g. / For intel C**

   .. code-block:: bash

      $ icc -lifcore program.f90 -L libtetrabz path}/src/ -I libtetrabz path}/src -ltetrabz -fopenmp
      $ mpiicc -lifcore program.f90 -L libtetrabz path}/src/  -I libtetrabz path}/src -ltetrabz_mpi -fopenmp
          

where ``{libtetrabz path}`` means the path to the top directory of ``libtetrabz``.
If you perform following setting, you do not need ``-L`` and ``-I`` option.

-  Copy

   ``.a`` files to a directory in ``$LIBRALLY_PATH``.

-  Copy

   ``.mod`` files to a directory in ``$INCLUDE``.
