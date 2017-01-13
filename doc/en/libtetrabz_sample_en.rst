Piece of sample code
====================

This sample shows the calculation of the charge density.

.. math::

   \begin{align}
   \rho(r) = 2 \sum_{n k} \theta(\varepsilon_{\rm F} - \varepsilon_{n k})
   |\varphi_{n k}(r)|^2
   \end{align}

.. code-block:: fortran

    subroutin calc_rho(nr,nb,ng,nelec,bvec,eig,ef,phi,rho)
      !
      use libtetrabz, only : libtetrabz_fermieng
      implicit none
      !
      integer,intent(in) :: nr ! number of r
      integer,intent(in) :: nb ! number of bands
      integer,intent(in) :: ng(3)
      ! k-point mesh
      real(8),intent(in) :: nelec ! number of electrons per spin
      real(8),intent(in) :: bvec(3,3) ! reciplocal lattice vector
      real(8),intent(in) :: eig(nb,ng(1),ng(2),ng(3)) ! Kohn-Sham eigenvalues
      real(8),intent(out) :: ef ! Fermi energy
      complex(8),intent(in) :: phi(nr,nb,ng(1),ng(2),ng(3)) ! Kohn-Sham orbitals
      real(8),intent(out) :: rho(nr) ! Charge density
      !
      integer :: ib, i1, i2, i3, ltetra
      real(8) :: wght(nb,ng(1),ng(2),ng(3))
      !
      ltetra = 2
      !
      call libtetrabz_fermieng(ltetra,bvec,nb,ng,eig,ng,wght,ef,nelec)
      !
      rho(1:nr) = 0d0
      do i1 = 1, ng(3)
         do i2 = 1, ng(2)
            do i1 = 1, ng(1)
               do ib = 1, nb
                  rho(1:nr) = rho(1:nr) + 2d0 * wght(ib,i1,i2,i3) &
                  &     * dble(conjg(phi(1:nr,ib,i1,i2,i3)) * phi(1:nr,ib,i1,i2,i3))
               end do
            end do
         end do
      end do
      !
    end subroutin calc_rho
        

