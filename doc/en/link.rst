Linking libtetrabz
==================

**e. g. / For intel fortran**

   .. code-block:: bash

      $ ifort program.f90 -L install_dir/lib -I install_dir/include -ltetrabz -fopenmp
      $ mpiifort program.f90 -L install_dir/lib -I install_dir/include -ltetrabz -fopenmp

**e. g. / For intel C**

   .. code-block:: bash

      $ icc -lifcore program.f90 -L install_dir/lib -I install_dir/include -ltetrabz -fopenmp
      $ mpiicc -lifcore program.f90 -L install_dir/lib -I install_dir/include -ltetrabz_mpi -fopenmp
