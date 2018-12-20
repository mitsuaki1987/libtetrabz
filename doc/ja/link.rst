ライブラリのリンク方法
======================

**例/ intel fortran の場合**

   .. code-block:: bash

      $ ifort program.f90 -L install_dir/lib -I install_dir/include -ltetrabz -fopenmp
      $ mpiifort program.f90 -L install_dir/lib -I install_dir/include -ltetrabz -fopenmp
          
**例/ intel C の場合**

   .. code-block:: bash

      $ icc -lifcore program.f90 -L install_dir/lib -I install_dir/include -ltetrabz -fopenmp
      $ mpiicc -lifcore program.f90 -L install_dir/lib -I install_dir/include -ltetrabz_mpi -fopenmp
