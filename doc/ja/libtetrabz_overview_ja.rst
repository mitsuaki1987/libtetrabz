はじめに
========

この文書ではテトラへドロン法ライブラリ ``libtetrabz`` についての解説を行っている.
``libtetrabz`` は線形テトラへドロン法もしくは最適化線形テトラへドロン法 :ref:`[1] <ref>`
を用いて全エネルギーや電荷密度, 部分状態密度,
分極関数等を計算するためのライブラリ群である.
このライブラリには, 軌道エネルギーをインプットとして,

.. math::

   \begin{align}
   \sum_{n n' k} F(\varepsilon_{n k}, \varepsilon_{n' k+q})X_{n n' k}
   = \sum_{n n' k} w_{n n' k} X_{n n' k}
   \end{align}

のような積分における, 積分重み :math:`w_{n n' k}` を出力するサブルーチンを,
各種計算について取り揃えている. 具体的には以下の計算に対応している.

.. math::

   \begin{align}
   \sum_{n k}
   \theta(\varepsilon_{\rm F} - \varepsilon_{n k})
   X_{n k}
   \end{align}

.. math::
 
   \begin{align}
   \sum_{n k}
   \delta(\omega - \varepsilon_{n k})
   X_{n k}(\omega)
   \end{align}

.. math::

   \begin{align}
   \sum_{n n' k}
   \delta(\varepsilon_{\rm F} - \varepsilon_{n k})
   \delta(\varepsilon_{\rm F} - \varepsilon'_{n' k})
   X_{n n' k}
   \end{align}

.. math::

   \begin{align}
   \sum_{n n' k}
   \theta(\varepsilon_{\rm F} - \varepsilon_{n k})
   \theta(\varepsilon_{n k} - \varepsilon'_{n' k})
   X_{n n' k}
   \end{align}

.. math::

   \begin{align}
   \sum_{n n' k}
   \frac{
   \theta(\varepsilon_{\rm F} - \varepsilon_{n k})
   \theta(\varepsilon'_{n' k} - \varepsilon_{\rm F})}
   {\varepsilon'_{n' k} - \varepsilon_{n k}}
   X_{n n' k}
   \end{align}

.. math::

   \begin{align}
   \sum_{n n' k}
   \theta(\varepsilon_{\rm F} - \varepsilon_{n k})
   \theta(\varepsilon'_{n' k} - \varepsilon_{\rm F})
   \delta(\varepsilon'_{n' k} - \varepsilon_{n k} - \omega)
   X_{n n' k}(\omega)
   \end{align}

.. math::

   \begin{align}
   \sum_{n n' k}
   \frac{
   \theta(\varepsilon_{\rm F} - \varepsilon_{n k})
   \theta(\varepsilon'_{n' k} - \varepsilon_{\rm F})}
   {\varepsilon'_{n' k} - \varepsilon_{n k} + i \omega}
   X_{n n' k}(\omega) 
   \end{align}

