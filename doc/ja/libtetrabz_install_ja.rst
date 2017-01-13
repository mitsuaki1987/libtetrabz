インストール方法
================

このパッケージに含まれているファイル
------------------------------------

``doc/manual_en.html`` : マニュアル(英語)

``doc/manual_jp.html`` : マニュアル(日本語)(今開いて見ているファイル)

``examples/`` : ライブラリ使用例

``src/`` : ソースファイルディレクトリ

``Makefile``, ``make.sys`` :
Makeファイルおよびコンパイル環境設定ファイル

要件
----

以下のものが必要となる.

-  fortran コンパイラ

-  MPI ライブラリ (MPI/ハイブリッド並列版を利用する場合)

インストール手順
----------------

#. 以下の場所から ``.tar.gz`` ファイルをダウンロードする.

   http://osdn.jp/projects/libtetrabz/releases/
               
#. ダウンロードした ``.tar.gz`` ファイルを展開し,
   出来たディレクトリに入る.

   .. code-block:: bash

      $ tar xzvf libtetrabz_1.0.1.tar.gz
      $ cd libtetrabz
               

#. 自分の環境に合わせて ``make.sys`` の以下の変数を書き換える.

   ``TOPDIR`` : 展開してできたディレクトリの絶対パス

   ``F90`` : シリアル用fortran コンパイルコマンド (gfortran, ifort, frt等)

   ``MPIF90`` : MPI用fortran コンパイルコマンド (mpif90, mpiifort, mpifrt 等)

   ``FFLAGS`` : fortranコンパイルオプション

#. 次のコマンドを実行しコンパイルする.

   .. code-block:: bash

      $ make 

   コンパイルが成功すると ``src/`` に以下のファイルが生成される.

   ::

       src/libtetrabz.a
       src/libtetrabz.mod
       src/libtetrabz_mpi.a
       src/libtetrabz_mpi.mod
               
   .. note::
       
      シリアル版のみ ``make`` したい場合は次の様にする.

      .. code-block:: bash

         $ cd src 
         $ make libtetrabz.a
               
また,
``example/`` 以下のライブラリ使用例のプログラムもコンパイルされる.

``example/dos.x`` :
立法格子シングルバンドタイトバインディングモデルのDOSを 計算する.
ソースコードは ``dos.f90``

   .. figure:: ../figs/dos.png
               :scale: 50

               dos.xを使って計算した,
               立法格子タイトバインディング模型の状態密度.
               実線は十分多くの :math:`k` 点を利用した時の結果.
               " :math:`+` ",
               " :math:`\times` "はそれぞれ
               :math:`8\times8\times8 k` グリッドでの
               線形テトラへドロン法および最適化テトラへドロン法の結果.

``example/lindhard.x`` : リントハルト関数を計算する.
ソースコードは ``lindhard.f90``

   .. figure:: ../figs/lindhard.png
               :scale: 50

               lindhard.xを使って計算したLindhard関数.
               実線は解析的な結果.
               " :math:`+` ", " :math:`\times` "はそれぞれ
               :math:`8\times8\times8 k`
               グリッドでの線形テトラへドロン法および最適化テトラへドロン法の結果.

