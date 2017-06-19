MODULE libtetrabz_doubledelta_mod
  !
  IMPLICIT NONE
  !
CONTAINS
!
! Compute doubledelta
!
SUBROUTINE libtetrabz_mpi_doubledelta(ltetra0,comm0,bvec,nb0,nge,eig1,eig2,ngw,wght0)
  !
  USE mpi, ONLY : MPI_IN_PLACE, MPI_DOUBLE_PRECISION, MPI_SUM
  USE libtetrabz_vals, ONLY : ltetra, ng, nb, nk0, indx1, indx2, indx3
  USE libtetrabz_common, ONLY : libtetrabz_initialize, libtetrabz_interpol_weight
  USE libtetrabz_doubledelta_mod, ONLY : libtetrabz_doubledelta_main
  !
  USE libtetrabz_mpi_routines, ONLY : comm, libtetrabz_mpi_kgrid
  !
  INTEGER,INTENT(IN) :: ltetra0, comm0, nge(3), ngw(3), nb0
  REAL(8),INTENT(IN) :: bvec(3,3)
  REAL(8),INTENT(IN) :: eig1(nb0,PRODUCT(nge(1:3))), eig2(nb0,PRODUCT(nge(1:3)))
  REAL(8),INTENT(OUT) :: wght0(nb0,nb0,PRODUCT(ngw(1:3)))
  !
  INTEGER :: ierr, nn
  REAL(8),ALLOCATABLE :: wght1(:,:)
  !
  ltetra = ltetra0
  comm = comm0
  nb = nb0
  ng(1:3) = nge(1:3)
  nn = nb * nb
  !
  CALL libtetrabz_initialize(bvec)
  CALL libtetrabz_mpi_kgrid()
  !
  ALLOCATE(wght1(nn, nk0))
  CALL libtetrabz_doubledelta_main(eig1,eig2,wght1)
  !
  CALL libtetrabz_interpol_weight(nn,ngw,nge,wght0,wght1)
  !
  DEALLOCATE(wght1, indx1, indx2, indx3)
  !
  CALL MPI_allREDUCE(MPI_IN_PLACE, wght0, nn * PRODUCT(ngw(1:3)), &
  &                  MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)
  !
END SUBROUTINE libtetrabz_mpi_doubledelta
!
! Main SUBROUTINE for Delta(E1) * Delta(E2)
!
SUBROUTINE libtetrabz_doubledelta_main(eig1,eig2,ddel)
  !
  USE libterabz_val, ONLY : nb, nk, nk0, fst, lst, indx1, indx2, wlsm
  IMPLICIT NONE
  !
  REAL(8),INTENT(IN) :: eig1(nb,nk), eig2(nb,nk)
  REAL(8),INTENT(OUT) :: dbldelta(nb,nb,nk0)
  !
  INTEGER :: it, ib, ii, nn
  REAL(8) :: e(4), V, &
  &          ei1(4,nb), ei2(3), ej1(4,nb), ej2(3,nb), &
  &          w1(nb,4), w2(nb,3)
  !
  dbldelta(1:nb,1:nb,1:nk0) = 0d0
  !
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP & SHARED(fst,lst,nb,nn,indx1,indx2,wlsm,eig1,eig2,w0,dbldelta) &
  !$OMP & PRIVATE(ib,it,ii,e,a,tmp,tmp2,w1,w2,ei,ej,ej2,V)
  !
  DO it = fst, lst
     !
     DO ib = 1, nb
        ei1(1:4,ib) = MATMUL(wlsm(1:4,1:20), eig1(ib, indx1(1:20,it)))
        ej1(1:4,ib) = MATMUL(wlsm(1:4,1:20), eig2(ib, indx1(1:20,it)))
     END DO
     !
     !$OMP DO
     DO ib = 1, nb
        !
        w1(1:nb,1:4) = 0d0
        e(1:4) = ei1(1:4, ib)
        CALL libtetrabz_sort(4,e,sind)
        !
        IF(e(1) < 0d0 .AND. 0d0 <= e(2)) THEN
           !
           CALL libtetrabz_tsmall_a1(e,V,tsmall)
           !
           IF(V > thr) THEN
              !
              ej2(1:3,1:nb) = MATMUL(tsmall(1:3,1:4), ej1(sind(1:4),1:nb))
              CALL libtetrabz_doubledelta2(ej2,w2)
              w1(1:nb,sind(1:4)) = w1(1:nb,            sind(1:4)) &
              &       + V * MATMUL(w2(1:nb,1:3), tsmall(1:3,1:4))
              !
           END IF
           !
        ELSE IF( e(2) < 0d0 .AND. 0d0 <= e(3)) THEN
           !
           CALL libtetrabz_tsmall_b1(e,V,tsmall)
           !
           IF(V > thr) THEN
              !
              ej2(1:3,1:nb) = MATMUL(tsmall(1:3,1:4), ej1(sind(1:4),1:nb))
              CALL libtetrabz_doubledelta2(ej2,w2)
              w1(1:nb,sind(1:4)) = w1(1:nb,            sind(1:4)) &
              &       + V * MATMUL(w2(1:nb,1:3), tsmall(1:3,1:4))
              !
           END IF
           !
           CALL libtetrabz_tsmall_b2(e,V,tsmall)
           !
           IF(V > thr) THEN
              !
              ej2(1:3,1:nb) = MATMUL(tsmall(1:3,1:4), ej1(sind(1:4),1:nb))
              CALL libtetrabz_doubledelta2(ej2,w2)
              w1(1:nb,sind(1:4)) = w1(1:nb,            sind(1:4)) &
              &       + V * MATMUL(w2(1:nb,1:3), tsmall(1:3,1:4))
              !
           END IF
           !
        ELSE IF(e(3) < 0d0 .AND. 0d0 < e(4)) THEN
           !
           CALL libtetrabz_tsmall_c1(e,V,tsmall)
           !
           IF(V > thr) THEN
              !
              ej2(1:3,1:nb) = MATMUL(tsmall(1:3,1:4), ej1(sind(1:4),1:nb))
              CALL libtetrabz_doubledelta2(ej2,w2)
              w1(1:nb,sind(1:4)) = w1(1:nb,            sind(1:4)) &
              &       + V * MATMUL(w2(1:nb,1:3), tsmall(1:3,1:4))
              !
           END IF
           !
        END IF
        !
        dbldelta(1:nb,ib,indx2(1:20,it)) = dbldelta(1:nb,ib,      indx2(1:20,it)) &
        &                               + MATMUL(w1(1:nb,1:4), wlsm(1:4,1:20))
        !
     END DO ! ib
     !$OMP END DO NOWAIT
     !
  END DO ! it
  !
  !$OMP END PARALLEL
  !
  dbldelta(1:nb,1:nb,1:nk0) = dbldelta(1:nb,1:nb,1:nk0) / DBLE(6 * nk)
  !
END SUBROUTINE libtetrabz_doubledelta_main
!
! 2nd step of tetrahedra method.
!
SUBROUTINE libtetrabz_doubledelta2(ej1,w1)
  !
  USE libterabz_val, ONLY : nb
  IMPLICIT NONE
  !
  REAL(8),INTENT(IN) :: ej(3,nb)
  REAL(8),INTENT(INOUT) :: w(nb,3)
  !
  INTEGER :: ib, ii, sind(3)
  REAL(8) :: e(3), a(3,3), V
  !
  DO ib = 1, nb
     !
     IF(maxval(ABS(ej(ib,1:3))) < 1d-10) STOP "Nesting !!"
     !
     w1(ib,1:3) = 0d0
     e(1:3) = ej(1:3, ib)
     CALL libtetrabz_sort(3,e,sind)
     !
     DO ii = 1, 3
        a(1:3,ii) = (0d0 - e(ii)) / (e(1:3) - e(ii))
     END DO
     !
     IF((e(1) < 0d0 .AND. 0d0 <= e(2)) .OR. (e(1) <= 0d0 .AND. 0d0 < e(2))) THEN
        !
        !V = a(2,1) * a(3,1) / (0d0 - e(1)) 
        V = a(2,1)           / (e(3) - e(1)) 
        !
        w(ib,sind(1)) = V * (a(1,2) + a(1,3))
        w(ib,sind(2)) = V * a(2,1)
        w(ib,sind(3)) = V * a(3,1)
        !
     ELSE IF((e(2) <= 0d0 .AND. 0d0 < e(3)) .OR. (e(2) < 0d0 .AND. 0d0 <= e(3))) THEN
        !
        !V = a(1,3) * a(2,3) / (e(3) - 0d0) 
        V = a(1,3)           / (e(3) - e(2)) 
        !
        w(ib,sind(1)) = V * a(1,3)
        w(ib,sind(2)) = V * a(2,3)
        w(ib,sind(3)) = V * (a(3,1) + a(3,2))
        !
     END IF
     !
  END DO ! ib
  !
END SUBROUTINE libtetrabz_doubledelta2
!
END MODULE libtetrabz_doubledelta_mod
