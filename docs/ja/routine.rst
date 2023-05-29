各サブルーチンの説明
====================

以下のサブルーチンを任意のプログラム内で

.. code-block:: fortran

   USE libtetrabz, ONLY : libtetrabz_occ
    
   CALL libtetrabz_occ(ltetra,bvec,nb,nge,eig,ngw,wght)
        
のように呼び出して使用できる.
サブルーチン名はすべて ``libtetrabz_`` からはじまる.

C言語で書かれたプログラムから呼び出す場合には次のようにする.

.. code-block:: c

   #include "libtetrabz.h"
   
   libtetrabz_occ(&ltetra,bvec,&nb,nge,eig,ngw,wght)
        
変数はすべてポインタとして渡す.
配列はすべて1次元配列として定義し一番左の添字が内側のループとなるようにする.
またMPI/ハイブリッド並列のときにライブラリに渡すコミュニケーター変数を, 
次のようにC/C++のものからfortranのものに変換する。

.. code-block:: c

   comm_f = MPI_Comm_c2f(comm_c);

全エネルギー, 電荷密度等(占有率の計算)
--------------------------------------

.. math::

   \begin{align}
   \sum_{n}
   \int_{\rm BZ} \frac{d^3 k}{V_{\rm BZ}}
   \theta(\varepsilon_{\rm F} -
   \varepsilon_{n k}) X_{n k}
   \end{align}

.. code-block:: fortran

    CALL libtetrabz_occ(ltetra,bvec,nb,nge,eig,ngw,wght,comm)

パラメーター
    
   .. code-block:: fortran
                   
      INTEGER,INTENT(IN) :: ltetra
   ..
   
      テトラへドロン法の種類を決める.
      1 :math:`\cdots` 線形テトラへドロン法,
      2 :math:`\cdots` 最適化線形テトラへドロン法 :ref:`[1] <ref>`

   .. code-block:: fortran
                   
      REAL(8),INTENT(IN) :: bvec(3,3)
   ..
   
      逆格子ベクトル. 単位は任意で良い.
      逆格子の形によって四面体の切り方を決めるため,
      それらの長さの比のみが必要であるため.

   .. code-block:: fortran
                   
      INTEGER,INTENT(IN) :: nb
   ..
   
      バンド本数

   .. code-block:: fortran
                   
      INTEGER,INTENT(IN) :: nge(3)
   ..
   
      軌道エネルギーのメッシュ数.

   .. code-block:: fortran
                   
      REAL(8),INTENT(IN) :: eig(nb,nge(1),nge(2),nge(3))
   ..
   
      軌道エネルギー.
      Fermiエネルギーを基準とすること( :math:`\varepsilon_{\rm F} = 0` ).

   .. code-block:: fortran
                   
      INTEGER,INTENT(IN) :: ngw(3)
   ..
   
      ``ngw(3)`` : (入力, 整数配列) 積分重みの :math:`k` メッシュ.
      ``nge`` と違っていても構わない(:ref:`app` 参照).

   .. code-block:: fortran
                   
      REAL(8),INTENT(OUT) :: wght(nb,ngw(1),ngw(2),ngw(3))
   ..
   
      ``wght(nb,ngw(1),ngw(2),ngw(3))`` : (出力, 実数配列) 積分重み

   .. code-block:: fortran
                   
      INTEGER,INTENT(IN),OPTIONAL :: comm
   ..
   
      オプショナル引数. 
      MPIのコミニュケーター( ``MPI_COMM_WORLD`` など)を入れる.
      libtetrabz を内部でMPI/Hybrid並列するときのみ入力する.
      C言語では使用しないときには ``NULL`` を入れる.

Fermi エネルギー(占有率も同時に計算する)
----------------------------------------

.. math::

   \begin{align}
   \sum_{n}
   \int_{\rm BZ} \frac{d^3 k}{V_{\rm BZ}}
   \theta(\varepsilon_{\rm F} -
   \varepsilon_{n k}) X_{n k} 
   \end{align}

.. code-block:: fortran

    CALL libtetrabz_fermieng(ltetra,bvec,nb,nge,eig,ngw,wght,ef,nelec,comm)
        
パラメーター
    
   .. code-block:: fortran
                   
      INTEGER,INTENT(IN) :: ltetra
   ..
   
      テトラへドロン法の種類を決める.
      1 :math:`\cdots` 線形テトラへドロン法,
      2 :math:`\cdots` 最適化線形テトラへドロン法 :ref:`[1] <ref>`

   .. code-block:: fortran
                         
      REAL(8),INTENT(IN) :: bvec(3,3)
   ..
   
      逆格子ベクトル. 単位は任意で良い.
      逆格子の形によって四面体の切り方を決めるため,
      それらの長さの比のみが必要であるため.

   .. code-block:: fortran
                         
      INTEGER,INTENT(IN) :: nb
   ..
   
      バンド本数

   .. code-block:: fortran
                         
      INTEGER,INTENT(IN) :: nge(3)
   ..
   
      軌道エネルギーのメッシュ数.

   .. code-block:: fortran
                         
      REAL(8),INTENT(IN) :: eig(nb,nge(1),nge(2),nge(3))
   ..
   
      軌道エネルギー.

   .. code-block:: fortran
                         
      INTEGER,INTENT(IN) :: nge(3)
   ..
   
      軌道エネルギーのメッシュ数.
      
   .. code-block:: fortran
                         
      INTEGER,INTENT(IN) :: ngw(3)
   ..
   
      積分重みの :math:`k` メッシュ.
      ``nge`` と違っていても構わない(:ref:`app` 参照).

   .. code-block:: fortran
                         
      REAL(8),INTENT(OUT) :: wght(nb,ngw(1),ngw(2),ngw(3))
   ..
   
      積分重み

   .. code-block:: fortran
                         
      REAL(8),INTENT(OUT) :: ef
   ..
   
      Fermi エネルギー

   .. code-block:: fortran
                         
      REAL(8),INTENT(IN) :: nelec
   ..
   
      スピンあたりの(荷)電子数

   .. code-block:: fortran
                   
      INTEGER,INTENT(IN),OPTIONAL :: comm
   ..
   
      オプショナル引数. 
      MPIのコミニュケーター( ``MPI_COMM_WORLD`` など)を入れる.
      libtetrabz を内部でMPI/Hybrid並列するときのみ入力する.
      C言語では使用しないときには ``NULL`` を入れる.

(部分)状態密度
--------------

.. math::

   \begin{align}
   \sum_{n}
   \int_{\rm BZ} \frac{d^3 k}{V_{\rm BZ}}
   \delta(\omega - \varepsilon_{n k})
   X_{n k}(\omega) 
   \end{align}

.. code-block:: fortran

   CALL libtetrabz_dos(ltetra,bvec,nb,nge,eig,ngw,wght,ne,e0,comm)
        
パラメーター

   .. code-block:: fortran
                         
      INTEGER,INTENT(IN) :: ltetra
   ..
   
      テトラへドロン法の種類を決める.
      1 :math:`\cdots` 線形テトラへドロン法,
      2 :math:`\cdots` 最適化線形テトラへドロン法 :ref:`[1] <ref>`

   .. code-block:: fortran
                         
      REAL(8),INTENT(IN) :: bvec(3,3)
   ..
   
      逆格子ベクトル. 単位は任意で良い.
      逆格子の形によって四面体の切り方を決めるため,
      それらの長さの比のみが必要であるため.

   .. code-block:: fortran
                         
      INTEGER,INTENT(IN) :: nb
   ..
   
      バンド本数

   .. code-block:: fortran
                         
      INTEGER,INTENT(IN) :: nge(3)
   ..
   
      軌道エネルギーの :math:`k` メッシュ数.

   .. code-block:: fortran
                         
      REAL(8),INTENT(IN) :: eig(nb,nge(1),nge(2),nge(3))
   ..
   
      軌道エネルギー.

   .. code-block:: fortran
                         
      INTEGER,INTENT(IN) :: ngw(3)
   ..
   
      積分重みの :math:`k` メッシュ.
      ``nge`` と違っていても構わない(:ref:`app` 参照).

   .. code-block:: fortran
                         
      REAL(8),INTENT(OUT) :: wght(ne,nb,ngw(1),ngw(2),ngw(3))
   ..
   
      積分重み

   .. code-block:: fortran
                         
      INTEGER,INTENT(IN) :: ne
   ..
   
      状態密度を計算するエネルギー点数

   .. code-block:: fortran
                         
      REAL(8),INTENT(IN) :: e0(ne)
   ..
   
      状態密度を計算するエネルギー

   .. code-block:: fortran
                   
      INTEGER,INTENT(IN),OPTIONAL :: comm
   ..
   
      オプショナル引数. 
      MPIのコミニュケーター( ``MPI_COMM_WORLD`` など)を入れる.
      libtetrabz を内部でMPI/Hybrid並列するときのみ入力する.
      C言語では使用しないときには ``NULL`` を入れる.

積分状態密度
------------

.. math::

   \begin{align}
   \sum_{n}
   \int_{\rm BZ} \frac{d^3 k}{V_{\rm BZ}}
   \theta(\omega - \varepsilon_{n k})
   X_{n k}(\omega) 
   \end{align}

.. code-block:: fortran

   CALL libtetrabz_intdos(ltetra,bvec,nb,nge,eig,ngw,wght,ne,e0,comm)
        
パラメーター

   .. code-block:: fortran
                         
      INTEGER,INTENT(IN) :: ltetra
   ..
   
      テトラへドロン法の種類を決める.
      1 :math:`\cdots` 線形テトラへドロン法,
      2 :math:`\cdots` 最適化線形テトラへドロン法 :ref:`[1] <ref>`

   .. code-block:: fortran
                         
      REAL(8),INTENT(IN) :: bvec(3,3)
   ..
   
      逆格子ベクトル. 単位は任意で良い.
      逆格子の形によって四面体の切り方を決めるため,
      それらの長さの比のみが必要であるため.

   .. code-block:: fortran
                         
      INTEGER,INTENT(IN) :: nb
   ..
   
      バンド本数

   .. code-block:: fortran
                         
      INTEGER,INTENT(IN) :: nge(3)
   ..
   
      軌道エネルギーの :math:`k` メッシュ数.

   .. code-block:: fortran
                         
      REAL(8),INTENT(IN) :: eig(nb,nge(1),nge(2),nge(3))
   ..
   
      軌道エネルギー.

   .. code-block:: fortran
                         
      INTEGER,INTENT(IN) :: ngw(3)
   ..
   
      積分重みの :math:`k` メッシュ.
      ``nge`` と違っていても構わない(:ref:`app` 参照).

   .. code-block:: fortran
                         
      REAL(8),INTENT(OUT) :: wght(ne,nb,ngw(1),ngw(2),ngw(3))
   ..
   
      積分重み

   .. code-block:: fortran
                         
      INTEGER,INTENT(IN) :: ne
   ..
   
      状態密度を計算するエネルギー点数

   .. code-block:: fortran
                         
      REAL(8),INTENT(IN) :: e0(ne)
   ..
   
      状態密度を計算するエネルギー

   .. code-block:: fortran
                   
      INTEGER,INTENT(IN),OPTIONAL :: comm
   ..
   
      オプショナル引数. 
      MPIのコミニュケーター( ``MPI_COMM_WORLD`` など)を入れる.
      libtetrabz を内部でMPI/Hybrid並列するときのみ入力する.
      C言語では使用しないときには ``NULL`` を入れる.

ネスティング関数, Fröhlich パラメーター
---------------------------------------

.. math::

   \begin{align}
   \sum_{n n'}
   \int_{\rm BZ} \frac{d^3 k}{V_{\rm BZ}}
   \delta(\varepsilon_{\rm F} -
   \varepsilon_{n k}) \delta(\varepsilon_{\rm F} - \varepsilon'_{n' k})
   X_{n n' k}
   \end{align}

.. code-block:: fortran

    CALL libtetrabz_dbldelta(ltetra,bvec,nb,nge,eig1,eig2,ngw,wght,comm)
        
パラメーター

   .. code-block:: fortran
                         
      INTEGER,INTENT(IN) :: ltetra
   ..
   
      テトラへドロン法の種類を決める.
      1 :math:`\cdots` 線形テトラへドロン法,
      2 :math:`\cdots` 最適化線形テトラへドロン法 :ref:`[1] <ref>`

   .. code-block:: fortran
                         
      REAL(8),INTENT(IN) :: bvec(3,3)
   ..
   
      逆格子ベクトル. 単位は任意で良い.
      逆格子の形によって四面体の切り方を決めるため,
      それらの長さの比のみが必要であるため.

   .. code-block:: fortran
                         
      INTEGER,INTENT(IN) :: nb
   ..
   
      バンド本数

   .. code-block:: fortran
                         
      INTEGER,INTENT(IN) :: nge(3)
   ..
   
      軌道エネルギーの :math:`k` メッシュ数.

   .. code-block:: fortran
                         
      REAL(8),INTENT(IN) :: eig1(nb,nge(1),nge(2),nge(3))
   ..
   
      軌道エネルギー.
      Fermi エネルギーを基準とすること( :math:`\varepsilon_{\rm F}=0` ).
      ``eig2`` も同様.

   .. code-block:: fortran
                         
      REAL(8),INTENT(IN) :: eig2(nb,nge(1),nge(2),nge(3))
   ..
   
      軌道エネルギー.
      移行運動量の分だけグリッドをずらしたものなど.

   .. code-block:: fortran
                         
      INTEGER,INTENT(IN) :: ngw(3)
   ..
   
      積分重みの :math:`k` メッシュ.
      ``nge`` と違っていても構わない(:ref:`app` 参照).

   .. code-block:: fortran
                         
      REAL(8),INTENT(OUT) :: wght(nb,nb,ngw(1),ngw(2),ngw(3))
   ..
   
      積分重み

   .. code-block:: fortran
                   
      INTEGER,INTENT(IN),OPTIONAL :: comm
   ..
   
      オプショナル引数. 
      MPIのコミニュケーター( ``MPI_COMM_WORLD`` など)を入れる.
      libtetrabz を内部でMPI/Hybrid並列するときのみ入力する.
      C言語では使用しないときには ``NULL`` を入れる.

DFPT 計算の一部
---------------

.. math::

   \begin{align}
   \sum_{n n'}
   \int_{\rm BZ} \frac{d^3 k}{V_{\rm BZ}}
   \theta(\varepsilon_{\rm F} -
   \varepsilon_{n k}) \theta(\varepsilon_{n k} - \varepsilon'_{n' k}) 
   X_{n n' k}
   \end{align}

.. code-block:: fortran

    CALL libtetrabz_dblstep(ltetra,bvec,nb,nge,eig1,eig2,ngw,wght,comm)
        
パラメーター

   .. code-block:: fortran
                         
      INTEGER,INTENT(IN) :: ltetra
   ..
   
      テトラへドロン法の種類を決める.
      1 :math:`\cdots` 線形テトラへドロン法,
      2 :math:`\cdots` 最適化線形テトラへドロン法 :ref:`[1] <ref>`

   .. code-block:: fortran
                         
      REAL(8),INTENT(IN) :: bvec(3,3)
   ..
   
      逆格子ベクトル. 単位は任意で良い.
      逆格子の形によって四面体の切り方を決めるため,
      それらの長さの比のみが必要であるため.

   .. code-block:: fortran
                         
      INTEGER,INTENT(IN) :: nb
   ..
   
      バンド本数

   .. code-block:: fortran
                         
      INTEGER,INTENT(IN) :: nge(3)
   ..
   
      軌道エネルギーのメッシュ数.

   .. code-block:: fortran
                         
      REAL(8),INTENT(IN) :: eig1(nb,nge(1),nge(2),nge(3))
   ..
   
      軌道エネルギー.
      Fermi エネルギーを基準とすること
      ( :math:`\varepsilon_{\rm F}=0` ). ``eig2`` も同様.

   .. code-block:: fortran
                         
      REAL(8),INTENT(IN) :: eig2(nb,nge(1),nge(2),nge(3))
   ..
   
      軌道エネルギー.
      移行運動量の分だけグリッドをずらしたものなど.

   .. code-block:: fortran
                         
      INTEGER,INTENT(IN) :: ngw(3)
   ..
   
      積分重みの :math:`k` メッシュ. ``nge``
      と違っていても構わない(:ref:`app` 参照).

   .. code-block:: fortran
                         
      REAL(8),INTENT(OUT) :: wght(nb,nb,ngw(1),ngw(2),ngw(3))
   ..
   
      積分重み

   .. code-block:: fortran
                   
      INTEGER,INTENT(IN),OPTIONAL :: comm
   ..
   
      オプショナル引数. 
      MPIのコミニュケーター( ``MPI_COMM_WORLD`` など)を入れる.
      libtetrabz を内部でMPI/Hybrid並列するときのみ入力する.
      C言語では使用しないときには ``NULL`` を入れる.

独立分極関数(静的)
------------------

.. math::

   \begin{align}
   \sum_{n n'}
   \int_{\rm BZ} \frac{d^3 k}{V_{\rm BZ}}
   \frac{\theta(\varepsilon_{\rm F} - \varepsilon_{n k})
   \theta(\varepsilon'_{n' k} - \varepsilon_{\rm F})}
   {\varepsilon'_{n' k} - \varepsilon_{n k}}
   X_{n n' k} 
   \end{align}

.. code-block:: fortran

    CALL libtetrabz_polstat(ltetra,bvec,nb,nge,eig1,eig2,ngw,wght,comm)
        
パラメーター

   .. code-block:: fortran
                         
      INTEGER,INTENT(IN) :: ltetra
   ..
   
      テトラへドロン法の種類を決める.
      1 :math:`\cdots` 線形テトラへドロン法,
      2 :math:`\cdots` 最適化線形テトラへドロン法 :ref:`[1] <ref>`

   .. code-block:: fortran
                   
      REAL(8),INTENT(IN) :: bvec(3,3)
   ..
   
      逆格子ベクトル. 単位は任意で良い.
      逆格子の形によって四面体の切り方を決めるため,
      それらの長さの比のみが必要であるため.

   .. code-block:: fortran
                         
      INTEGER,INTENT(IN) :: nb
   ..
   
      バンド本数

   .. code-block:: fortran
                         
      INTEGER,INTENT(IN) :: nge(3)
   ..
   
      軌道エネルギーのメッシュ数.

   .. code-block:: fortran
                         
      REAL(8),INTENT(IN) :: eig1(nb,nge(1),nge(2),nge(3))
   ..
   
      軌道エネルギー.
      Fermi エネルギーを基準とすること
      ( :math:`\varepsilon_{\rm F}=0` ). ``eig2`` も同様.

   .. code-block:: fortran
                         
      REAL(8),INTENT(IN) :: eig2(nb,nge(1),nge(2),nge(3))
   ..
   
      軌道エネルギー.
      移行運動量の分だけグリッドをずらしたものなど.

   .. code-block:: fortran
                         
      INTEGER,INTENT(IN) :: ngw(3)
   ..
   
      積分重みの :math:`k` メッシュ.
      ``nge`` と違っていても構わない(:ref:`app` 参照).

   .. code-block:: fortran
                         
      REAL(8),INTENT(OUT) :: wght(nb,nb,ngw(1),ngw(2),ngw(3))
   ..
   
      積分重み

   .. code-block:: fortran
                   
      INTEGER,INTENT(IN),OPTIONAL :: comm
   ..
   
      オプショナル引数. 
      MPIのコミニュケーター( ``MPI_COMM_WORLD`` など)を入れる.
      libtetrabz を内部でMPI/Hybrid並列するときのみ入力する.
      C言語では使用しないときには ``NULL`` を入れる.

フォノン線幅等
--------------

.. math::

   \begin{align}
   \sum_{n n'}
   \int_{\rm BZ} \frac{d^3 k}{V_{\rm BZ}}
   \theta(\varepsilon_{\rm F} -
   \varepsilon_{n k}) \theta(\varepsilon'_{n' k} - \varepsilon_{\rm F})
   \delta(\varepsilon'_{n' k} - \varepsilon_{n k} - \omega)
   X_{n n' k}(\omega) 
   \end{align}

.. code-block:: fortran

    CALL libtetrabz_fermigr(ltetra,bvec,nb,nge,eig1,eig2,ngw,wght,ne,e0,comm)
        
パラメーター

   .. code-block:: fortran
                         
      INTEGER,INTENT(IN) :: ltetra
   ..
   
      テトラへドロン法の種類を決める.
      1 :math:`\cdots` 線形テトラへドロン法,
      2 :math:`\cdots` 最適化線形テトラへドロン法 :ref:`[1] <ref>`

   .. code-block:: fortran
                         
      REAL(8),INTENT(IN) :: bvec(3,3)
   ..
   
      逆格子ベクトル. 単位は任意で良い.
      逆格子の形によって四面体の切り方を決めるため,
      それらの長さの比のみが必要であるため.

   .. code-block:: fortran
                   
      INTEGER,INTENT(IN) :: nb
   ..
   
      バンド本数

   .. code-block:: fortran
                         
      INTEGER,INTENT(IN) :: nge(3)
   ..
   
      軌道エネルギーのメッシュ数.

   .. code-block:: fortran
                         
      REAL(8),INTENT(IN) :: eig1(nb,nge(1),nge(2),nge(3))
   ..
   
      軌道エネルギー.
      Fermi エネルギーを基準とすること
      ( :math:`\varepsilon_{\rm F}=0` ). ``eig2`` も同様.

   .. code-block:: fortran
                         
      REAL(8),INTENT(IN) :: eig2(nb,nge(1),nge(2),nge(3))
   ..
   
      軌道エネルギー.
      移行運動量の分だけグリッドをずらしたものなど.

   .. code-block:: fortran
                         
      INTEGER,INTENT(IN) :: ngw(3)
   ..
   
      積分重みの :math:`k` メッシュ.
      ``nge`` と違っていても構わない(:ref:`app` 参照).

   .. code-block:: fortran
                         
      REAL(8),INTENT(OUT) :: wght(ne,nb,nb,ngw(1),ngw(2),ngw(3))
   ..
   
      積分重み

   .. code-block:: fortran
                         
      INTEGER,INTENT(IN) :: ne
   ..
   
      フォノンモード数

   .. code-block:: fortran
                         
      REAL(8),INTENT(IN) :: e0(ne)
   ..
   
      フォノン振動数

   .. code-block:: fortran
                   
      INTEGER,INTENT(IN),OPTIONAL :: comm
   ..
   
      オプショナル引数. 
      MPIのコミニュケーター( ``MPI_COMM_WORLD`` など)を入れる.
      libtetrabz を内部でMPI/Hybrid並列するときのみ入力する.
      C言語では使用しないときには ``NULL`` を入れる.

分極関数(複素振動数)
--------------------

.. math::

   \begin{align}
   \sum_{n n'}
   \int_{\rm BZ} \frac{d^3 k}{V_{\rm BZ}}
   \frac{\theta(\varepsilon_{\rm F} - \varepsilon_{n k})
   \theta(\varepsilon'_{n' k} - \varepsilon_{\rm F})}
   {\varepsilon'_{n' k} - \varepsilon_{n k} + i \omega}
   X_{n n' k}(\omega) 
   \end{align}

.. code-block:: fortran

    CALL libtetrabz_polcmplx(ltetra,bvec,nb,nge,eig1,eig2,ngw,wght,ne,e0,comm)
        
パラメーター

   .. code-block:: fortran
                         
      INTEGER,INTENT(IN) :: ltetra
   ..
   
      テトラへドロン法の種類を決める.
      1 :math:`\cdots` 線形テトラへドロン法,
      2 :math:`\cdots` 最適化線形テトラへドロン法 :ref:`[1] <ref>`

   .. code-block:: fortran
                         
      REAL(8),INTENT(IN) :: bvec(3,3)
   ..
   
      逆格子ベクトル. 単位は任意で良い.
      逆格子の形によって四面体の切り方を決めるため,
      それらの長さの比のみが必要であるため.

   .. code-block:: fortran
                   
      INTEGER,INTENT(IN) :: nb
   ..
   
      バンド本数

   .. code-block:: fortran
                         
      INTEGER,INTENT(IN) :: nge(3)
   ..
   
      軌道エネルギーのメッシュ数.

   .. code-block:: fortran
                         
      REAL(8),INTENT(IN) :: eig1(nb,nge(1),nge(2),nge(3))
   ..
   
      軌道エネルギー.
      Fermi エネルギーを基準とすること
      ( :math:`\varepsilon_{\rm F}=0` ). ``eig2`` も同様.

   .. code-block:: fortran
                         
      REAL(8),INTENT(IN) :: eig2(nb,nge(1),nge(2),nge(3))
   ..
   
      軌道エネルギー.
      移行運動量の分だけグリッドをずらしたものなど.

   .. code-block:: fortran
                         
      INTEGER,INTENT(IN) :: ngw(3)
   ..
   
      積分重みの :math:`k` メッシュ.
      ``nge`` と違っていても構わない(:ref:`app` 参照).

   .. code-block:: fortran
                         
      COMPLEX(8),INTENT(OUT) :: wght(ne,nb,nb,ngw(1),ngw(2),ngw(3))
   ..
   
      積分重み .

   .. code-block:: fortran
                         
      INTEGER,INTENT(IN) :: ne
   ..
   
      計算を行う虚振動数の点数

   .. code-block:: fortran
                         
      COMPLEX(8),INTENT(IN) :: e0(ne)
   ..
   
      計算を行う複素振動数

   .. code-block:: fortran
                   
      INTEGER,INTENT(IN),OPTIONAL :: comm
   ..
   
      オプショナル引数. 
      MPIのコミニュケーター( ``MPI_COMM_WORLD`` など)を入れる.
      libtetrabz を内部でMPI/Hybrid並列するときのみ入力する.
      C言語では使用しないときには ``NULL`` を入れる.

