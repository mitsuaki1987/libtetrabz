.. _app:

補遺
====

逆補間
------

積分

.. math::

   \begin{align}
   \langle X \rangle = \sum_{k} X_k w(\varepsilon_k)
   \end{align}

を計算するとする. このとき,

-   :math:`w` は :math:`\varepsilon_k` に敏感な関数(階段関数 :math:`\cdot` デルタ関数等)であり,
    なるべく細かいグリッド上の :math:`\varepsilon_k` が必要である.

-   :math:`X_k` を求めるための計算コストが :math:`\varepsilon_k` の計算コストよりかなり大きい.

という場合には :math:`X_k` のグリッドを補間により増やす方法が有効である.
それは,

#.  :math:`\varepsilon_k` を細かい :math:`k` グリッド上で計算する.

#.  :math:`X_k` を粗いグリッド上で計算し, それを補間(線形補間, 多項式補間,
    スプライン補間など)して細かいグリッド上での値を得る.

.. math::
   
   \begin{align}
   X_k^{\rm dense} = \sum_{k'}^{\rm coarse}
   F_{k k'} X_{k'}^{\rm coarse}
   \end{align}

#. 細かい :math:`k` グリッドで上記の積分を行う.

.. math::
   
   \begin{align}
   \langle X \rangle = \sum_{k}^{\rm dense}
   X_k^{\rm dense} w_k^{\rm dense}
   \end{align}

という流れで行われる.

さらに,
この計算と同じ結果を得るように粗いグリッド上での積分重み
:math:`w_k^{\rm coarse}` を  :math:`w_k{\rm dense}` から求める(逆補間)ことも可能である
( :ref:`[2] <ref>` のAppendix).
すなわち,

.. math::
   
   \begin{align}
   \sum_k^{\rm dense} X_k^{\rm dense} w_k^{\rm dense}
   = \sum_k^{\rm coarse} X_k^{\rm coarse} w_k^{\rm coarse}
   \end{align}

が満たされる事を要請すると

.. math::

   \begin{align}
   w_k^{\rm coarse} = \sum_k^{\rm dense} F_{k' k}
   w_{k'}^{\rm dense}
   \end{align}

となる. この場合の計算手順は,

#. 細かい :math:`k` グリッド上の  :math:`\varepsilon_k` から
   :math:`w_k^{\rm dense}` を計算する.

#. 逆補間により :math:`w_k^{\rm coarse}` を求める.

#. 粗いグリッド上での :math:`X_k` との積和を行う.

となる. このライブラリ内の全ルーチンはこの逆補間の機能を備えており,
軌道エネルギーの :math:`k` グリッドと重み関数の :math:`k` グリッドを
異なる値にすると逆補間された :math:`w_k^{\rm coarse}` が出力される.

2重デルタ関数の積分
-------------------

