インストール方法
================

主なファイルとディレクトリ
--------------------------

- ``doc/`` : マニュアルのディレクトリ
   - ``doc/index.html`` : 目録ページ 
- ``src/`` : ライブラリのソースコードのディレクトリ
- ``example/`` : ライブラリ使用例のディレクトリ
- ``test/`` : ライブラリのビルドのテスト用ディレクトリ
- ``configure`` : ビルド環境設定用スクリプト

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

      $ tar xzvf libtetrabz-version.tar.gz
      $ cd libtetrabz-version
               
最もシンプルには次のとおりである.

.. code-block:: bash

   $ ./configure --prefix=install_dir

これにより, ビルドに必要なコンパイラやライブラリ等の環境のチェックが行われ,
Makefile等が作成される.
ただし ``install_dir`` はインストール先のディレクトリの絶対パスとする (以後各自のディレクトリ名で読み替えること).
なにも指定しないと ``/use/local/`` が設定され, 後述の ``make install`` で
``/usr/local/lib`` 内にライブラリが置かれる (したがって, 管理者権限がない場合には ``install_dir`` を
別の場所に指定しなければならない).
``configure`` にはこの他にも様々なオプションがあり, 必要に応じて用途や環境に合わせてそれらを使用する.
詳しくは :ref:`configoption` を参照.

``configure`` の実行が正常に行われ, ``Makefile`` が生成された後は

.. code-block:: bash

   $ make

とタイプしてライブラリ等のビルドを行う.これが成功したのちに

.. code-block:: bash

   $ make install

とすると, ライブラリが ``install_dir/lib`` に置かれる.
``make install`` をしなくても, ビルドをしたディレクトリ内にあるライブラリやミニアプリを使うことは可能であるが, 使い勝手がやや異なる.

共有リンクを行ったプログラムの実行時にライブラリを探しにいけるよう,
環境変数 ``LD_LIBRARY_PATH`` にlibtetrabzをインストールしたディレクトリを追加する.

.. code-block:: bash

   $ export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:install_dir/lib

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

.. _configoption:

configureのオプション
---------------------

configureには多数のオプションと変数があり, それらを組み合わせて指定する.
指定しない場合にはデフォルト値が使われる.

.. code-block:: bash

  $ ./configure --prefix=/home/libtetrabz/ --with-mpi=yes FC=mpif90

おもなものを次に挙げる.

``---prefix``

   デフォルト: ``---prefix=/usr/local/``.
   ライブラリ等のインストールを行うディレクトリツリーを指定する.

``--with-mpi``

   デフォルト: ``--with-mpi=no`` (MPIを用いない).
   MPIを用いるか (``--with-mpi=yes``), 否かを指定する.

``--with-openmp``

   デフォルト: ``--with-openmp=yes`` (OpenMPを用いる).
   OpenMPを用いるか否か (``--with-openmp=no``) を指定する.

``--enable-shared``

   デフォルト: ``--enable-shared``.
   共有ライブラリを作成するか否か

``--enable-static``

   デフォルト: ``--enable-static``.
   静的ライブラリを作成するか否か.

``FC``

   デフォルト: システムにインストールされているfortranコンパイラをスキャンして,
   自動的に設定する. ``--with-mpi`` を指定した時にはそれに応じたコマンド
   (``mpif90`` 等)を自動で探し出して設定する. 
   ``configure`` の最後に出力される ``FC`` が望んだものでは無かった場合には
   ``./configure FC=gfortran`` のように手で指定する.

``--help``

   このオプションを指定した時には, ビルドの環境設定は行われず,
   上記を含めたすべてのオプションを表示する.
