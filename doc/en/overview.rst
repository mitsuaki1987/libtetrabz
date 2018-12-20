Introduction
============

This document explains a tetrahedron method library ``libtetrabz``.
``libtetrabz`` is a library to calculate the total energy, the charge
density, partial density of states, response functions, etc. in a solid
by using the optimized tetrahedron method :ref:`[1] <ref>`.
Subroutines in this library receive the orbital (Kohn-Sham) energies as an input and
calculate weights :math:`w_{n n' k}` for integration such as

.. math::

   \begin{align}
   \sum_{n n'}
   \int_{\rm BZ} \frac{d^3 k}{V_{\rm BZ}}
   F(\varepsilon_{n k}, \varepsilon_{n' k+q})X_{n n' k}
   = \sum_{n n'} \sum_{k}^{N_k} w_{n n' k} X_{n n' k}
   \end{align}

``libtetrabz`` supports following Brillouin-zone integrations

.. math::

   \begin{align}
   \sum_{n}
   \int_{\rm BZ} \frac{d^3 k}{V_{\rm BZ}}
   \theta(\varepsilon_{\rm F} - \varepsilon_{n k})
   X_{n k}
   \end{align}

.. math::
 
   \begin{align}
   \sum_{n}
   \int_{\rm BZ} \frac{d^3 k}{V_{\rm BZ}}
   \delta(\omega - \varepsilon_{n k})
   X_{n k}(\omega)
   \end{align}

.. math::

   \begin{align}
   \sum_{n n'}
   \int_{\rm BZ} \frac{d^3 k}{V_{\rm BZ}}
   \delta(\varepsilon_{\rm F} - \varepsilon_{n k})
   \delta(\varepsilon_{\rm F} - \varepsilon'_{n' k})
   X_{n n' k}
   \end{align}

.. math::

   \begin{align}
   \sum_{n n'}
   \int_{\rm BZ} \frac{d^3 k}{V_{\rm BZ}}
   \theta(\varepsilon_{\rm F} - \varepsilon_{n k})
   \theta(\varepsilon_{n k} - \varepsilon'_{n' k})
   X_{n n' k}
   \end{align}

.. math::

   \begin{align}
   \sum_{n n'}
   \int_{\rm BZ} \frac{d^3 k}{V_{\rm BZ}}
   \frac{
   \theta(\varepsilon_{\rm F} - \varepsilon_{n k})
   \theta(\varepsilon'_{n' k} - \varepsilon_{\rm F})}
   {\varepsilon'_{n' k} - \varepsilon_{n k}}
   X_{n n' k}
   \end{align}

.. math::

   \begin{align}
   \sum_{n n'}
   \int_{\rm BZ} \frac{d^3 k}{V_{\rm BZ}}
   \theta(\varepsilon_{\rm F} - \varepsilon_{n k})
   \theta(\varepsilon'_{n' k} - \varepsilon_{\rm F})
   \delta(\varepsilon'_{n' k} - \varepsilon_{n k} - \omega)
   X_{n n' k}(\omega)
   \end{align}

.. math::

   \begin{align}
   \sum_{n n'}
   \int_{\rm BZ} \frac{d^3 k}{V_{\rm BZ}}
   \frac{
   \theta(\varepsilon_{\rm F} - \varepsilon_{n k})
   \theta(\varepsilon'_{n' k} - \varepsilon_{\rm F})}
   {\varepsilon'_{n' k} - \varepsilon_{n k} + i \omega}
   X_{n n' k}(\omega) 
   \end{align}

