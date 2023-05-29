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
:math:`w_k^{\rm coarse}` を  :math:`w_k{\rm dense}` から求める
**逆補間** ( :ref:`[2] <ref>` のAppendix)も可能である.
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

For the integration

.. math::

   \begin{align}
   \sum_{n n' k} \delta(\varepsilon_{\rm F} -
   \varepsilon_{n k}) \delta(\varepsilon_{\rm F} - \varepsilon'_{n' k})
   X_{n n' k}
   \end{align}

first, we cut out one or two triangles where
:math:`\varepsilon_{n k} = \varepsilon_{\rm F}` from a tetrahedron
and evaluate :math:`\varepsilon_{n' k+q}` at the corners of each triangles as

.. math::

   \begin{align}
   \varepsilon'^{k+q}_{i} = \sum_{j=1}^4 F_{i j}(
   \varepsilon_1^{k}, \cdots, \varepsilon_{4}^{k}, \varepsilon_{\rm F}) 
   \epsilon_{j}^{k+q}.
   \end{align}
   
Then we calculate :math:`\delta(\varepsilon_{n' k+q} - \varepsilon{\rm F})`
in each triangles and obtain weights of corners.
This weights of corners are mapped into those of corners of the original tetrahedron as

.. math::
   
   \begin{align}
   W_{i} = \sum_{j=1}^3 \frac{S}{\nabla_k \varepsilon_k}F_{j i}(
   \varepsilon_{1}^k, \cdots, \varepsilon_{4}^k, \varepsilon_{\rm F}) 
   W'_{j}.
   \end{align}

:math:`F_{i j}` and :math:`\frac{S}{\nabla_k \varepsilon_k}` are calculated as follows 
(:math:`a_{i j} \equiv (\varepsilon_i - \varepsilon_j)/(\varepsilon_{\rm F} - \varepsilon_j)`):

.. _dbldeltapng:

.. figure:: ../figs/dbldelta.png
   :scale: 100

   How to divide a tetrahedron 
   in the case of :math:`\epsilon_1 \leq \varepsilon_{\rm F} \leq \varepsilon_2` (a), 
   :math:`\varepsilon_2 \leq \varepsilon_{\rm F} \leq \varepsilon_3` (b), and
   :math:`\varepsilon_3 \leq \varepsilon_{\rm F} \leq \varepsilon_4` (c).

- When :math:`\varepsilon_1 \leq \varepsilon_{\rm F} \leq \varepsilon_2 \leq \varepsilon_3 \leq\varepsilon_4`
  [Fig. :num:`dbldeltapng` (a)], 

   .. math::
   
      \begin{align}
      F &= 
      \begin{pmatrix}
      a_{1 2} & a_{2 1} &       0 & 0 \\
      a_{1 3} &       0 & a_{3 1} & 0 \\
      a_{1 4} &       0 &       0 & a_{4 1}
      \end{pmatrix}, 
      \qquad
      \frac{S}{\nabla_k \varepsilon_k} = \frac{3 a_{2 1} a_{3 1} a_{4 1}}{\varepsilon_{\rm F} - \varepsilon_1}
      \end{align}
  
- When :math:`\varepsilon_1 \leq \varepsilon_2 \leq \varepsilon_{\rm F} \leq \varepsilon_3 \leq\varepsilon_4`
  [Fig. :num:`dbldeltapng` (b)], 

   .. math::
   
      \begin{align}
      F &= 
      \begin{pmatrix}
      a_{1 3} &       0 & a_{3 1} & 0 \\
      a_{1 4} &       0 &       0 & a_{4 1} \\
      0 & a_{2 4} &       0 & a_{4 2}
      \end{pmatrix}, 
      \qquad
      \frac{S}{\nabla_k \varepsilon_k} = \frac{3 a_{3 1} a_{4 1} a_{2 4}}{\varepsilon_{\rm F} - \varepsilon_1}
      \end{align}
  
   .. math::
   
      \begin{align}
      F &= 
      \begin{pmatrix}
      a_{1 3} &       0 & a_{3 1} & 0 \\
      0 & a_{2 3} & a_{3 2} & 0 \\
      0 & a_{2 4} &       0 & a_{4 2}
      \end{pmatrix}, 
      \qquad
      \frac{S}{\nabla_k \varepsilon_k} = \frac{3 a_{2 3} a_{3 1} a_{4 2}}{\varepsilon_{\rm F} - \varepsilon_1}
      \end{align}

- When :math:`\varepsilon_1 \leq \varepsilon_2 \leq \varepsilon_3 \leq \varepsilon_{\rm F} \leq \varepsilon_4`
  [Fig. :num:`dbldeltapng` (c)], 

   .. math::
   
      \begin{align}
      F &= 
      \begin{pmatrix}
      a_{1 4} &       0 &       0 & a_{4 1} \\
      a_{1 3} & a_{2 4} &       0 & a_{4 2} \\
      a_{1 2} &       0 & a_{3 4} & a_{4 3}
      \end{pmatrix}, 
      \qquad
      \frac{S}{\nabla_k \varepsilon_k} = \frac{3 a_{1 4} a_{2 4} a_{3 4}}{\varepsilon_1 - \varepsilon_{\rm F}}
      \end{align}

Weights on each corners of the triangle are computed as follows
[(:math:`a'_{i j} \equiv (\varepsilon'_i - \varepsilon'_j)/(\varepsilon_{\rm F} - \varepsilon'_j)`)]:

- When :math:`\varepsilon'_1 \leq \varepsilon_{\rm F} \leq \varepsilon'_2 \leq \varepsilon'_3` [Fig. :num:`dbldeltapng` (d)], 

   .. math::
   
      \begin{align}
      W'_1 = L (a'_{1 2} + a'_{1 3}), \qquad
      W'_2 = L a'_{2 1}, \qquad
      W'_3 = L a'_{3 1}, \qquad
      L \equiv \frac{a'_{2 1} a'_{3 1}}{\varepsilon_{\rm F} - \varepsilon'_{1}}
      \end{align}

- When :math:`\varepsilon'_1 \leq \varepsilon'_2 \leq \varepsilon_{\rm F} \leq \varepsilon'_3` [Fig. :num:`dbldeltapng` (e)], 

   .. math::
   
      \begin{align}
      W'_1 = L a'_{1 3}, \qquad
      W'_2 = L a'_{2 3}, \qquad
      W'_3 = L (a'_{3 1} + a'_{3 2}), \qquad
      L \equiv \frac{a'_{1 3} a'_{2 3}}{\varepsilon'_{3} - \varepsilon_{\rm F}} 
      \end{align}
