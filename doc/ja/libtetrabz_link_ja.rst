ライブラリのリンク方法
======================

**例/ intel fortran の場合**

   .. code-block:: bash

      $ ifort program.f90 -L {libtetrabzパス}/src/ -I {libtetrabzパス}/src -ltetrabz -fopenmp
      $ mpiifort program.f90 -L {libtetrabzパス}/src/ -I {libtetrabzパス}/src -ltetrabz_mpi -fopenmp
          
**例/ intel C の場合**

   .. code-block:: bash

      $ icc -lifcore program.f90 -L {libtetrabzパス}/src/ -I {libtetrabzパス}/src -ltetrabz -fopenmp
      $ mpiicc -lifcore program.f90 -L {libtetrabzパス}/src/ -I {libtetrabzパス}/src -ltetrabz_mpi -fopenmp
          
``-L``, ``-I`` オプションについては,
以下のように設定すればつける必要は無い.

-  ``.a`` ファイルを環境変数
   ``LIBRALLY_PATH`` に含まれているディレクトリへコピーする.

-  ``.mod`` ファイルを環境変数
   ``INCLUDE`` に含まれているディレクトリへコピーする.
