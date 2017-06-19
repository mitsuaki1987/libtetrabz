MODULE libetrabz_dblstep_mod
  !
  IMPLICIT NONE
  !
CONTAINS
!
! Compute Occ * Step
!
SUBROUTINE libtetrabz_mpi_doublestep(ltetra0,comm0,bvec,nb0,nge,eig1,eig2,ngw,wght0)
  !
  USE mpi, ONLY : MPI_IN_PLACE, MPI_DOUBLE_PRECISION, MPI_SUM
  USE libtetrabz_vals, ONLY : ltetra, ng, nb, nk0, indx1, indx2, indx3
  USE libtetrabz_common, ONLY : libtetrabz_initialize, libtetrabz_interpol_weight
  USE libtetrabz_doublestep_mod, ONLY : libtetrabz_doublestep_main
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
  !
  CALL libtetrabz_initialize(bvec)
  CALL libtetrabz_mpi_kgrid()
  nn = nb * nb
  !
  ALLOCATE(wght1(nn, nk0))
  CALL libtetrabz_doublestep_main(eig1,eig2,wght1)
  !
  CALL libtetrabz_interpol_weight(nn,ngw,nge,wght0,wght1)
  !
  DEALLOCATE(wght1, indx1, indx2, indx3)
  !
  CALL MPI_allREDUCE(MPI_IN_PLACE, wght0, nn * PRODUCT(ngw(1:3)), &
  &                  MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)
  !
END SUBROUTINE libtetrabz_mpi_doublestep
!
! Main SUBROUTINE for Theta(- E1) * Theta(E1 - E2)
!
SUBROUTINE libtetrabz_dblstep_main(eig1,eig2,dblstep)
  !
  USE libterabz_val, ONLY : nb, nk, nk0, fst, lst, indx1, indx2, wlsm
  IMPLICIT NONE
  !
  REAL(8),INTENT(IN) :: eig1(nb,nk), eig2(nb,nk)
  REAL(8),INTENT(OUT) :: dblstep(nb,nb,nk)
  !
  INTEGER :: it, ib
  REAL(8) :: e(4), V, thr = 1d-10, &
  &          ei1(4,nb), ej1(4,nb), ei2(4), ej2(4,nb), w1(nb,4), w2(nb,4)
  !
  dblstep(1:nb,1:nb,1:nk0) = 0d0
  !
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP & SHARED(fst,lst,nb,indx1,indx2,wlsm,eig1,eig2,thr,dblstep) &
  !$OMP & PRIVATE(ib,it,e,a,tmp,tmp2,w1,w2,ei1,ej1,de,V)
  !
  DO it = fst, lst
     !
     DO ib = 1, nb
        ei1(1:4, ib) = MATMUL(wlsm(1:4,1:20) * eig1(ib, indx1(1:20,it)))
        ej1(1:4, ib) = MATMUL(wlsm(1:4,1:20) * eig2(ib, indx1(1:20,it)))
     END DO
     !
     !$OMP DO
     DO ib = 1, nb
        !
        w1(1:nb,1:4) = 0d0
        e(1:4) = ei1(1:4, ib)
        CALL libtetrabz_sort(4,e,sind)
        !
        IF(e(1) <= 0d0 .AND. 0d0 < e(2) ) THEN
           !
           CALL libtetrabz_tsmall_a1(e,V,tsmall)
           !
           IF(V > thr) THEN
              !
              ei2(1:4     ) = MATMUL(tsmall(1:4,1:4), ei1(sind(1:4),  ib))
              ej2(1:4,1:nb) = MATMUL(tsmall(1:4,1:4), ej1(sind(1:4),1:nb))
              CALL libtetrabz_dblstep2(ei2,ej2,w2)
              w1(1:nb,sind(1:4)) = w1(1:nb,              sind(1:4)) &
              &       + V * MATMUL(w2(1:nb, 1:4 ), tsmall(1:4,1:4))
              !
           END IF
           !
        ELSE IF(e(2) <= 0d0 .AND. 0d0 < e(3)) THEN
           !
           CALL libtetrabz_tsmall_b1(e,V,tsmall)
           !
           IF(V > thr) THEN
              !
              ei2(1:4     ) = MATMUL(tsmall(1:4,1:4), ei1(sind(1:4),  ib))
              ej2(1:4,1:nb) = MATMUL(tsmall(1:4,1:4), ej1(sind(1:4),1:nb))
              CALL libtetrabz_dblstep2(ei2,ej2,w2)
              w1(1:nb,sind(1:4)) = w1(1:nb,            sind(1:4)) &
              &       + V * MATMUL(w2(1:nb,1:4), tsmall(1:4,1:4))
              !
           END IF
           !
           CALL libtetrabz_tsmall_b2(e,V,tsmall)
           !
           IF(V > thr) THEN
              !
              ei2(1:4     ) = MATMUL(tsmall(1:4,1:4), ei1(sind(1:4),  ib))
              ej2(1:4,1:nb) = MATMUL(tsmall(1:4,1:4), ej1(sind(1:4),1:nb))
              CALL libtetrabz_dblstep2(ei2,ej2,w2)
              w1(1:nb,sind(1:4)) = w1(1:nb,            sind(1:4)) &
              &       + V * MATMUL(w2(1:nb,1:4), tsmall(1:4,1:4))
              !
           END IF
           !
           CALL libtetrabz_tsmall_b3(e,V,tsmall)
           !
           IF(V > thr) THEN
              !
              ei2(1:4     ) = MATMUL(tsmall(1:4,1:4), ei1(sind(1:4),  ib))
              ej2(1:4,1:nb) = MATMUL(tsmall(1:4,1:4), ej1(sind(1:4),1:nb))
              CALL libtetrabz_dblstep2(ei2,ej2,w2)
              w1(1:nb,sind(1:4)) = w1(1:nb,            sind(1:4)) &
              &       + V * MATMUL(w2(1:nb,1:4), tsmall(1:4,1:4))
              !
           END IF
           !
        ELSE IF( e(3) <= 0d0 .AND. 0d0 < e(4)) THEN
           !
           CALL libtetrabz_tsmall_c1(e,V,tsmall)
           !
           IF(V > thr) THEN
              !
              ei2(1:4     ) = MATMUL(tsmall(1:4,1:4), ei1(sind(1:4),  ib))
              ej2(1:4,1:nb) = MATMUL(tsmall(1:4,1:4), ej1(sind(1:4),1:nb))
              CALL libtetrabz_dblstep2(ei2,ej2,w2)
              w1(1:nb,sind(1:4)) = w1(1:nb,            sind(1:4)) &
              &       + V * MATMUL(w2(1:nb,1:4), tsmall(1:4,1:4))
              !
           END IF
           !
           CALL libtetrabz_tsmall_c2(e,V,tsmall)
           !
           IF(V > thr) THEN
              !
              ei2(1:4     ) = MATMUL(tsmall(1:4,1:4), ei1(sind(1:4),  ib))
              ej2(1:4,1:nb) = MATMUL(tsmall(1:4,1:4), ej1(sind(1:4),1:nb))
              CALL libtetrabz_dblstep2(ei2,ej2,w2)
              w1(1:nb,sind(1:4)) = w1(1:nb,            sind(1:4)) &
              &       + V * MATMUL(w2(1:nb,1:4), tsmall(1:4,1:4))
              !
           END IF
           !
           CALL libtetrabz_tsmall_c3(e,V,tsmall)
           !
           IF(V > thr) THEN
              !
              ei2(1:4     ) = MATMUL(tsmall(1:4,1:4), ei1(sind(1:4),  ib))
              ej2(1:4,1:nb) = MATMUL(tsmall(1:4,1:4), ej1(sind(1:4),1:nb))
              CALL libtetrabz_dblstep2(ei2,ej2,w2)
              w1(1:nb,sind(1:4)) = w1(1:nb,            sind(1:4)) &
              &       + V * MATMUL(w2(1:nb,1:4), tsmall(1:4,1:4))
              !
           END IF
           !
        ELSE IF(e(4) <= 0d0) THEN
           !
           ei2(1:4     ) = ei(1:4,  ib)
           ej2(1:4,1:nb) = ej(1:4,1:nb)
           CALL libtetrabz_dblstep2(ei2,ej2,w2)
           w1(1:nb,1:4) = w1(1:nb,1:4) + w2(1:nb,1:4)
           !
        ELSE
           !
           CYCLE
           !
        END IF
        !
        dblstep(1:nb,ib,indx2(1:20,it)) = dblstep(1:nb,ib,       indx2(1:20,it)) &
        &                             + MATMUL(w1(1:nb, 1:4), wlsm(1:4,1:20))
        !
     END DO ! ib
     !$OMP END DO NOWAIT
     !
  END DO ! it
  !
  !$OMP END PARALLEL
  !
  dblstep(1:nb,1:nb,1:nk0) = dblstep(1:nb,1:nb,1:nk0) / DBLE(6 * nk)
  !
END SUBROUTINE libtetrabz_dblstep_main
!
! Tetrahedra method for theta( - de)
!
SUBROUTINE libtetrabz_dblstep2(ei1,ej1,w1)
  !
  USE libterabz_val, ONLY : nb
  IMPLICIT NONE
  !
  REAL(8),INTENT(IN) :: ei1(4), ej1(4,nb)
  REAL(8),INTENT(OUT) :: w1(nb,4)
  !
  INTEGER :: ib
  REAL(8) :: V, w2(4), thr = 1d-8, e(4)
  !
  w2(1:4) = 0.25d0
  !
  DO ib = 1, nb
     !
     w1(ib,1:4) = 0d0
     e(1:4) = ei(1:4) - ej(1:4,ib)
     CALL libtetrabz_sort(4,e,sind)
     !
     IF(ABS(e(1)) < thr .AND. ABS(e(4)) < thr) THEN
        !
        ! Theta(0) = 0.5
        !
        V = 0.5d0
        w1(ib,1:4) = w1(ib,1:4) + V * w2(1:4)
        !
     ELSE IF((e(1) <= 0d0 .AND. 0d0 < e(2)) .OR. (e(1) < 0d0 .AND. 0d0 <= e(2))) THEN
        !
        CALL libtetrabz_tsmall_a1(e,V,tsmall)
        w1(ib,sind(1:4)) = w1(ib,sind(1:4)) + V * MATMUL(w2(1:4), tsmall(1:4,1:4))
        !
     ELSE IF((e(2) <= 0d0 .AND. 0d0 < e(3)) .OR. (e(2) < 0d0 .AND. 0d0 <= e(3))) THEN
        !
        CALL libtetrabz_tsmall_b1(e,V,tsmall)
        w1(ib,sind(1:4)) = w1(ib,sind(1:4)) + V * MATMUL(w2(1:4), tsmall(1:4,1:4))
        !
        CALL libtetrabz_tsmall_b2(e,V,tsmall)
        w1(ib,sind(1:4)) = w1(ib,sind(1:4)) + V * MATMUL(w2(1:4), tsmall(1:4,1:4))
        !
        CALL libtetrabz_tsmall_b3(e,V,tsmall)
        w1(ib,sind(1:4)) = w1(ib,sind(1:4)) + V * MATMUL(w2(1:4), tsmall(1:4,1:4))
        !
     ELSE IF((e(3) <= 0d0 .AND. 0d0 < e(4)) .OR. (e(3) < 0d0 .AND. 0d0 <= e(4))) THEN
        !
        CALL libtetrabz_tsmall_c1(e,V,tsmall)
        w1(ib,sind(1:4)) = w1(ib,sind(1:4)) + V * MATMUL(w2(1:4), tsmall(1:4,1:4))
        !
        CALL libtetrabz_tsmall_c2(e,V,tsmall)
        w1(ib,sind(1:4)) = w1(ib,sind(1:4)) + V * MATMUL(w2(1:4), tsmall(1:4,1:4))
        !
        CALL libtetrabz_tsmall_c3(e,V,tsmall)
        w1(ib,sind(1:4)) = w1(ib,sind(1:4)) + V * MATMUL(w2(1:4), tsmall(1:4,1:4))
        !
     ELSE IF(e(4) <= 0d0) THEN
        !
        w1(ib,1:4) = w1(ib,1:4) + w2(1:4)
        !
     END IF
     !
  END DO ! ib = 1, nb
  !
END SUBROUTINE libtetrabz_doublestep2
