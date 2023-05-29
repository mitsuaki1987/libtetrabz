サンプルコード(抜粋)
====================

以下では電荷密度

.. math::

   \begin{align}
   \rho(r) = 2 \sum_{n k} \theta(\varepsilon_{\rm F} - \varepsilon_{n k})
   |\varphi_{n k}(r)|^2
   \end{align}

を計算するサブルーチンを示す.
   
.. code-block:: fortran

    SUBROUTINE calc_rho(nr,nb,ng,nelec,bvec,eig,ef,phi,rho)
      !
      USE libtetrabz, ONLY : libtetrabz_fermieng
      IMPLICIT NONE
      !
      INTEGER,INTENT(IN) :: nr ! number of r
      INTEGER,INTENT(IN) :: nb ! number of bands
      INTEGER,INTENT(IN) :: ng(3)
      ! k-point mesh
      REAL(8),INTENT(IN) :: nelec ! number of electrons per spin
      REAL(8),INTENT(IN) :: bvec(3,3) ! reciplocal lattice vector
      REAL(8),INTENT(IN) :: eig(nb,ng(1),ng(2),ng(3)) ! Kohn-Sham eigenvalues
      REAL(8),INTENT(OUT) :: ef ! Fermi energy
      COMPLEX(8),INTENT(IN) :: phi(nr,nb,ng(1),ng(2),ng(3)) ! Kohn-Sham orbitals
      REAL(8),INTENT(OUT) :: rho(nr) ! Charge density
      !
      INTEGER :: ib, i1, i2, i3, ltetra
      REAL(8) :: wght(nb,ng(1),ng(2),ng(3))
      !
      ltetra = 2
      !
      CALL libtetrabz_fermieng(ltetra,bvec,nb,ng,eig,ng,wght,ef,nelec)
      !
      rho(1:nr) = 0d0
      DO i1 = 1, ng(3)
         DO i2 = 1, ng(2)
            DO i1 = 1, ng(1)
               DO ib = 1, nb
                  rho(1:nr) = rho(1:nr) + 2d0 * wght(ib,i1,i2,i3) &
                  &     * DBLE(CONJG(phi(1:nr,ib,i1,i2,i3)) * phi(1:nr,ib,i1,i2,i3))
               END DO
            END DO
         END DO
      END DO
      !
    END SUBROUTINE calc_rho
        

