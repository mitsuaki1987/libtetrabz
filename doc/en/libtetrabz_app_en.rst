.. _app:

Appendix: Inverse interpolation
===============================

We consider an integration as follows:

.. math::

   \begin{align}
   \langle X \rangle = \sum_{k} X_k w(\varepsilon_k)
   \end{align}

If this integration has conditions that

-  :math:`w(\varepsilon_k)` is sensitive to :math:`\varepsilon_k` (e. g. the
   stepfunction, the delta function, etc.) and requires
   :math:`\varepsilon_k` on a dense :math:`k` grid, and

-  the numerical cost to obtain :math:`X_k` is much larger than the cost for
   :math:`\varepsilon_k` (e. g. the polarization function),

it is efficient to interpolate :math:`X_k` into a denser :math:`k` grid and
evaluate that integration in a dense :math:`k` grid. This method is performed
as follows:

#. Calculate :math:`\varepsilon_k` on a dense :math:`k` grid.

#. Calculate :math:`X_k` on a coarse :math:`k` grid and obtain that on a dense :math:`k`
   grid by using the linear interpolation, the polynomial interpolation,
   the spline interpolation, etc.

.. math::
   
   \begin{align}
   X_k^{\rm dense} = \sum_{k'}^{\rm coarse}
   F_{k k'} X_{k'}^{\rm coarse}
   \end{align}

#. Evaluate that integration in the dense :math:`k` grid.

.. math::
   
   \begin{align}
   \langle X \rangle = \sum_{k}^{\rm dense}
   X_k^{\rm dense} w_k^{\rm dense}
   \end{align}

**The inverse interpolation method** arrows as to obtain the same result
to above without interpolating :math:`X_k` into a dense :math:`k` grid. In this
method, we map the integration weight on a dense :math:`k` grid into that on a
coarse :math:`k` grid (inverse interpolation). Therefore, if we require

.. math::
   
   \begin{align}
   \sum_k^{\rm dense} X_k^{\rm dense} w_k^{\rm dense}
   = \sum_k^{\rm coarse} X_k^{\rm coarse} w_k^{\rm coarse}
   \end{align}

we obtain

.. math::

   \begin{align}
   w_k^{\rm coarse} = \sum_k^{\rm dense} F_{k' k}
   w_{k'}^{\rm dense}
   \end{align}

The numerical procedure for this method is as follows:

#. Calculate the integration weight on a dense :math:`k` grid
   :math:`w_k^{\rm dense}` from :math:`\varepsilon_k` on a dense :math:`k` grid.

#. Obtain the integration weight on a coarse :math:`k` grid :math:`w_k^{\rm
   coarse}` by using the inverse interpolation method.

#. Evaluate that integration in a coarse :math:`k` grid where :math:`X_k` was
   calculated.

All routines in ``libtetrabz`` can perform the inverse interpolation
method; if we make :math:`k` grids for the orbital energy (``nge``) and the
integration weight (``ngw``) different, we obtain :math:`w_k^{\rm coarse}`
calculated by using the inverse interpolation method.

