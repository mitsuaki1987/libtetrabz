各サブルーチンの説明
====================

以下のサブルーチンを任意のプログラム内で

.. code-block:: fortran

   use libtetrabz, only : libtetrabz_occ
    
   call libtetrabz_occ(ltetra,bvec,nb,nge,eig,ngw,wght)
        
のように呼び出して使用できる.
サブルーチン名はすべて ``libtetrabz_`` からはじまる.
MPI版については ``libtetrabz_mpi_`` からはじまる.
またMPI版では ``libtetrabz_mpi`` モジュールをつかう.
シリアル版と MPI版の引数の違いはコミニュケータを指定する整数 ``comm`` だけで,
他は同じである.

C言語で書かれたプログラムから呼び出す場合には次のようにする.

.. code-block:: fortran

   #include "libtetrabz.h"
   
   libtetrabz_mp_libtetrabz_occ_(&ltetra,bvec,&nb,nge,eig,ngw,wght)
        
fortranサブルーチン名の前に ``libtetrabz_mp_`` を,
うしろに ``_`` をつけたものが Cでの関数名となる.
MPI版では ``libtetrabz_mpi.h`` をインクルードし,
fortranサブルーチン名の前には ``libtetrabz_mpi_mp_`` をつける.
変数はすべてポインタとして渡す.
配列はすべて1次元配列として定義し一番左の添字が内側のループとなるようにする.

全エネルギー, 電荷密度等(占有率の計算)
--------------------------------------

.. math::

   \begin{align}
   \sum_{n k} \theta(\varepsilon_{\rm F} -
   \varepsilon_{n k}) X_{n k}
   \end{align}

.. code-block:: fortran

    call libtetrabz_occ(ltetra,bvec,nb,nge,eig,ngw,wght)
    call libtetrabz_mpi_occ(ltetra,comm,bvec,nb,nge,eig,ngw,wght)

パラメーター
    
   .. code-block:: fortran
                   
      integer,intent(in) :: ltetra
   ..
   
      テトラへドロン法の種類を決める.
      1 :math:`\cdots` 線形テトラへドロン法,
      2 :math:`\cdots` 最適化線形テトラへドロン法

   .. code-block:: fortran
                   
      integer,intent(in) :: comm
   ..
   
      MPI 版のみ. コミニュケータ.

   .. code-block:: fortran
                   
      real(8),intent(in) :: bvec(3,3)
   ..
   
      逆格子ベクトル. 単位は任意で良い.
      逆格子の形によって四面体の切り方を決めるため,
      それらの長さの比のみが必要であるため.

   .. code-block:: fortran
                   
      integer,intent(in) :: nb
   ..
   
      バンド本数

   .. code-block:: fortran
                   
      integer,intent(in) :: nge(3)
   ..
   
      軌道エネルギーのメッシュ数.

   .. code-block:: fortran
                   
      real(8),intent(in) :: eig(nb,nge(1),nge(2),nge(3))
   ..
   
      軌道エネルギー.
      Fermiエネルギーを基準とすること( :math:`\varepsilon_{\rm F} = 0` ).

   .. code-block:: fortran
                   
      integer,intent(in) :: ngw(3)
   ..
   
      ``ngw(3)`` : (入力, 整数配列) 積分重みの :math:`k` メッシュ.
      ``nge`` と違っていても構わない(:ref:`app` 参照).

   .. code-block:: fortran
                   
      real(8),intent(out) :: wght(nb,ngw(1),ngw(2),ngw(3))
   ..
   
      ``wght(nb,ngw(1),ngw(2),ngw(3))`` : (出力, 実数配列) 積分重み

Fermi エネルギー(占有率も同時に計算する)
----------------------------------------

.. math::

   \begin{align}
   \sum_{n k} \theta(\varepsilon_{\rm F} -
   \varepsilon_{n k}) X_{n k} 
   \end{align}

.. code-block:: fortran

    call libtetrabz_fermieng(ltetra,bvec,nb,nge,eig,ngw,wght,ef,nelec)
    call libtetrabz_mpi_fermieng(ltetra,comm,bvec,nb,nge,eig,ngw,wght,ef,nelec)
        
パラメーター
    
   .. code-block:: fortran
                   
      integer,intent(in) :: ltetra
   ..
   
      テトラへドロン法の種類を決める.
      1 :math:`\cdots` 線形テトラへドロン法,
      2 :math:`\cdots` 最適化線形テトラへドロン法

   .. code-block:: fortran
                         
      integer,intent(in) :: comm
   ..
   
      ``comm`` : (入力, 整数) MPI 版のみ. コミニュケータ.

   .. code-block:: fortran
                         
      real(8),intent(in) :: bvec(3,3)
   ..
   
      逆格子ベクトル. 単位は任意で良い.
      逆格子の形によって四面体の切り方を決めるため,
      それらの長さの比のみが必要であるため.

   .. code-block:: fortran
                         
      integer,intent(in) :: nb
   ..
   
      バンド本数

   .. code-block:: fortran
                         
      integer,intent(in) :: nge(3)
   ..
   
      軌道エネルギーのメッシュ数.

   .. code-block:: fortran
                         
      real(8),intent(in) :: eig(nb,nge(1),nge(2),nge(3))
   ..
   
      軌道エネルギー.

   .. code-block:: fortran
                         
      integer,intent(in) :: nge(3)
   ..
   
      軌道エネルギーのメッシュ数.
      
   .. code-block:: fortran
                         
      integer,intent(in) :: ngw(3)
   ..
   
      積分重みの :math:`k` メッシュ.
      ``nge`` と違っていても構わない(:ref:`app` 参照).

   .. code-block:: fortran
                         
      real(8),intent(out) :: wght(nb,ngw(1),ngw(2),ngw(3))
   ..
   
      積分重み

   .. code-block:: fortran
                         
      real(8),intent(out) :: ef
   ..
   
      Fermi エネルギー

   .. code-block:: fortran
                         
      real(8),intent(in) :: nelec
   ..
   
      スピンあたりの(荷)電子数

(部分)状態密度
--------------

.. math::

   \begin{align}
   \sum_{n k} \delta(\omega - \varepsilon_{n k})
   X_{n k}(\omega) 
   \end{align}

.. code-block:: fortran

   call libtetrabz_dos(ltetra,bvec,nb,nge,eig,ngw,wght,ne,e0)
   call libtetrabz_mpi_dos(ltetra,comm,bvec,nb,nge,eig,ngw,wght,ne,e0)
        
パラメーター

   .. code-block:: fortran
                         
      integer,intent(in) :: ltetra
   ..
   
      テトラへドロン法の種類を決める.
      1 :math:`\cdots` 線形テトラへドロン法,
      2 :math:`\cdots` 最適化線形テトラへドロン法

   .. code-block:: fortran
                         
      integer,intent(in) :: comm
   ..
   
      MPI 版のみ. コミニュケーター.

   .. code-block:: fortran
                         
      real(8),intent(in) :: bvec(3,3)
   ..
   
      逆格子ベクトル. 単位は任意で良い.
      逆格子の形によって四面体の切り方を決めるため,
      それらの長さの比のみが必要であるため.

   .. code-block:: fortran
                         
      integer,intent(in) :: nb
   ..
   
      バンド本数

   .. code-block:: fortran
                         
      integer,intent(in) :: nge(3)
   ..
   
      軌道エネルギーの :math:`k` メッシュ数.

   .. code-block:: fortran
                         
      real(8),intent(in) :: eig(nb,nge(1),nge(2),nge(3))
   ..
   
      軌道エネルギー.

   .. code-block:: fortran
                         
      integer,intent(in) :: ngw(3)
   ..
   
      積分重みの :math:`k` メッシュ.
      ``nge`` と違っていても構わない(:ref:`app` 参照).

   .. code-block:: fortran
                         
      real(8),intent(out) :: wght(ne,nb,ngw(1),ngw(2),ngw(3))
   ..
   
      積分重み

   .. code-block:: fortran
                         
      integer,intent(in) :: ne
   ..
   
      状態密度を計算するエネルギー点数

   .. code-block:: fortran
                         
      real(8),intent(in) :: e0(ne)
   ..
   
      状態密度を計算するエネルギー

ネスティング関数, Fröhlich パラメーター
---------------------------------------

.. math::

   \begin{align}
   \sum_{n n' k} \delta(\varepsilon_{\rm F} -
   \varepsilon_{n k}) \delta(\varepsilon_{\rm F} - \varepsilon'_{n' k})
   X_{n n' k}
   \end{align}

.. code-block:: fortran

    call libtetrabz_doubledelta(ltetra,bvec,nb,nge,eig1,eig2,ngw,wght)
    call libtetrabz_mpi_doubledelta(ltetra,comm,bvec,nb,nge,eig1,eig2,ngw,wght)
        
パラメーター

   .. code-block:: fortran
                         
      integer,intent(in) :: ltetra
   ..
   
      テトラへドロン法の種類を決める.
      1 :math:`\cdots` 線形テトラへドロン法,
      2 :math:`\cdots` 最適化線形テトラへドロン法

   .. code-block:: fortran
                         
      integer,intent(in) :: comm
   ..
   
      ``comm`` : (入力, 整数) MPI 版のみ. コミニュケータ.

   .. code-block:: fortran
                         
      real(8),intent(in) :: bvec(3,3)
   ..
   
      逆格子ベクトル. 単位は任意で良い.
      逆格子の形によって四面体の切り方を決めるため,
      それらの長さの比のみが必要であるため.

   .. code-block:: fortran
                         
      integer,intent(in) :: nb
   ..
   
      バンド本数

   .. code-block:: fortran
                         
      integer,intent(in) :: nge(3)
   ..
   
      軌道エネルギーの :math:`k` メッシュ数.

   .. code-block:: fortran
                         
      real(8),intent(in) :: eig1(nb,nge(1),nge(2),nge(3))
   ..
   
      軌道エネルギー.
      Fermi エネルギーを基準とすること( :math:`\varepsilon_{\rm F}=0` ).
      ``eig2`` も同様.

   .. code-block:: fortran
                         
      real(8),intent(in) :: eig2(nb,nge(1),nge(2),nge(3))
   ..
   
      軌道エネルギー.
      移行運動量の分だけグリッドをずらしたものなど.

   .. code-block:: fortran
                         
      integer,intent(in) :: ngw(3)
   ..
   
      積分重みの :math:`k` メッシュ.
      ``nge`` と違っていても構わない(:ref:`app` 参照).

   .. code-block:: fortran
                         
      real(8),intent(out) :: wght(nb,nb,ngw(1),ngw(2),ngw(3))
   ..
   
      積分重み

DFPT 計算の一部
---------------

.. math::

   \begin{align}
   \sum_{n n' k} \theta(\varepsilon_{\rm F} -
   \varepsilon_{n k}) \theta(\varepsilon_{n k} - \varepsilon'_{n' k}) 
   X_{n n' k}
   \end{align}

.. code-block:: fortran

    call libtetrabz_occstep(ltetra,bvec,nb,nge,eig1,eig2,ngw,wght)
    call libtetrabz_mpi_occstep(ltetra,comm,bvec,nb,nge,eig1,eig2,ngw,wght)
        
パラメーター

   .. code-block:: fortran
                         
      integer,intent(in) :: ltetra
   ..
   
      テトラへドロン法の種類を決める.
      1 :math:`\cdots` 線形テトラへドロン法,
      2 :math:`\cdots` 最適化線形テトラへドロン法

   .. code-block:: fortran
                         
      integer,intent(in) :: comm
   ..
   
      MPI 版のみ. コミニュケータ.

   .. code-block:: fortran
                         
      real(8),intent(in) :: bvec(3,3)
   ..
   
      逆格子ベクトル. 単位は任意で良い.
      逆格子の形によって四面体の切り方を決めるため,
      それらの長さの比のみが必要であるため.

   .. code-block:: fortran
                         
      integer,intent(in) :: nb
   ..
   
      バンド本数

   .. code-block:: fortran
                         
      integer,intent(in) :: nge(3)
   ..
   
      軌道エネルギーのメッシュ数.

   .. code-block:: fortran
                         
      real(8),intent(in) :: eig1(nb,nge(1),nge(2),nge(3))
   ..
   
      軌道エネルギー.
      Fermi エネルギーを基準とすること
      ( :math:`\varepsilon_{\rm F}=0` ). ``eig2`` も同様.

   .. code-block:: fortran
                         
      real(8),intent(in) :: eig2(nb,nge(1),nge(2),nge(3))
   ..
   
      軌道エネルギー.
      移行運動量の分だけグリッドをずらしたものなど.

   .. code-block:: fortran
                         
      integer,intent(in) :: ngw(3)
   ..
   
      積分重みの :math:`k` メッシュ. ``nge``
      と違っていても構わない(:ref:`app` 参照).

   .. code-block:: fortran
                         
      real(8),intent(out) :: wght(nb,nb,ngw(1),ngw(2),ngw(3))
   ..
   
      積分重み

独立分極関数(静的,  :math:`\omega=0` )
--------------------------------------

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
        
パラメーター

   .. code-block:: fortran
                         
      integer,intent(in) :: ltetra
   ..
   
      テトラへドロン法の種類を決める.
      1 :math:`\cdots` 線形テトラへドロン法,
      2 :math:`\cdots` 最適化線形テトラへドロン法

   .. code-block:: fortran
                         
      integer,intent(in) :: comm
   ..
   
      MPI 版のみ. コミニュケータ.

   .. code-block:: fortran
                   
      real(8),intent(in) :: bvec(3,3)
   ..
   
      逆格子ベクトル. 単位は任意で良い.
      逆格子の形によって四面体の切り方を決めるため,
      それらの長さの比のみが必要であるため.

   .. code-block:: fortran
                         
      integer,intent(in) :: nb
   ..
   
      バンド本数

   .. code-block:: fortran
                         
      integer,intent(in) :: nge(3)
   ..
   
      軌道エネルギーのメッシュ数.

   .. code-block:: fortran
                         
      real(8),intent(in) :: eig1(nb,nge(1),nge(2),nge(3))
   ..
   
      軌道エネルギー.
      Fermi エネルギーを基準とすること
      ( :math:`\varepsilon_{\rm F}=0` ). ``eig2`` も同様.

   .. code-block:: fortran
                         
      real(8),intent(in) :: eig2(nb,nge(1),nge(2),nge(3))
   ..
   
      軌道エネルギー.
      移行運動量の分だけグリッドをずらしたものなど.

   .. code-block:: fortran
                         
      integer,intent(in) :: ngw(3)
   ..
   
      積分重みの :math:`k` メッシュ.
      ``nge`` と違っていても構わない(:ref:`app` 参照).

   .. code-block:: fortran
                         
      real(8),intent(out) :: wght(nb,nb,ngw(1),ngw(2),ngw(3))
   ..
   
      積分重み

フォノン線幅等
--------------

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
        
パラメーター

   .. code-block:: fortran
                         
      integer,intent(in) :: ltetra
   ..
   
      テトラへドロン法の種類を決める.
      1 :math:`\cdots` 線形テトラへドロン法,
      2 :math:`\cdots` 最適化線形テトラへドロン法

   .. code-block:: fortran
                         
      integer,intent(in) :: comm
   ..
   
      MPI 版のみ. コミニュケータ.

   .. code-block:: fortran
                         
      real(8),intent(in) :: bvec(3,3)
   ..
   
      逆格子ベクトル. 単位は任意で良い.
      逆格子の形によって四面体の切り方を決めるため,
      それらの長さの比のみが必要であるため.

   .. code-block:: fortran
                   
      integer,intent(in) :: nb
   ..
   
      バンド本数

   .. code-block:: fortran
                         
      integer,intent(in) :: nge(3)
   ..
   
      軌道エネルギーのメッシュ数.

   .. code-block:: fortran
                         
      real(8),intent(in) :: eig1(nb,nge(1),nge(2),nge(3))
   ..
   
      軌道エネルギー.
      Fermi エネルギーを基準とすること
      ( :math:`\varepsilon_{\rm F}=0` ). ``eig2`` も同様.

   .. code-block:: fortran
                         
      real(8),intent(in) :: eig2(nb,nge(1),nge(2),nge(3))
   ..
   
      軌道エネルギー.
      移行運動量の分だけグリッドをずらしたものなど.

   .. code-block:: fortran
                         
      integer,intent(in) :: ngw(3)
   ..
   
      積分重みの :math:`k` メッシュ.
      ``nge`` と違っていても構わない(:ref:`app` 参照).

   .. code-block:: fortran
                         
      real(8),intent(out) :: wght(ne,nb,nb,ngw(1),ngw(2),ngw(3))
   ..
   
      積分重み

   .. code-block:: fortran
                         
      integer,intent(in) :: ne
   ..
   
      フォノンモード数

   .. code-block:: fortran
                         
      real(8),intent(in) :: e0(ne)
   ..
   
      フォノン振動数

分極関数(虚振動数)
------------------

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
        
パラメーター

   .. code-block:: fortran
                         
      integer,intent(in) :: ltetra
   ..
   
      テトラへドロン法の種類を決める.
      1 :math:`\cdots` 線形テトラへドロン法,
      2 :math:`\cdots` 最適化線形テトラへドロン法

   .. code-block:: fortran
                         
      integer,intent(in) :: comm
   ..
   
      MPI 版のみ. コミニュケータ.

   .. code-block:: fortran
                         
      real(8),intent(in) :: bvec(3,3)
   ..
   
      逆格子ベクトル. 単位は任意で良い.
      逆格子の形によって四面体の切り方を決めるため,
      それらの長さの比のみが必要であるため.

   .. code-block:: fortran
                   
      integer,intent(in) :: nb
   ..
   
      バンド本数

   .. code-block:: fortran
                         
      integer,intent(in) :: nge(3)
   ..
   
      軌道エネルギーのメッシュ数.

   .. code-block:: fortran
                         
      real(8),intent(in) :: eig1(nb,nge(1),nge(2),nge(3))
   ..
   
      軌道エネルギー.
      Fermi エネルギーを基準とすること
      ( :math:`\varepsilon_{\rm F}=0` ). ``eig2`` も同様.

   .. code-block:: fortran
                         
      real(8),intent(in) :: eig2(nb,nge(1),nge(2),nge(3))
   ..
   
      軌道エネルギー.
      移行運動量の分だけグリッドをずらしたものなど.

   .. code-block:: fortran
                         
      integer,intent(in) :: ngw(3)
   ..
   
      積分重みの :math:`k` メッシュ.
      ``nge`` と違っていても構わない(:ref:`app` 参照).

   .. code-block:: fortran
                         
      real(8),intent(out) :: wght(2,ne,nb,nb,ngw(1),ngw(2),ngw(3))
   ..
   
      積分重み .
      1番目の次元は実部と虚部を格納する.

   .. code-block:: fortran
                         
      integer,intent(in) :: ne
   ..
   
      計算を行う虚振動数の点数

   .. code-block:: fortran
                         
      real(8),intent(in) :: e0(ne)
   ..
   
      計算を行う虚振動数

