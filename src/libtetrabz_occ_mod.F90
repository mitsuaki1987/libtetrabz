MODULE libterabz_occ_mod
  !
  IMPLICIT NONE
  !
CONTAINS
!
! Compute occupation
!
SUBROUTINE libtetrabz_occ(ltetra0,comm0,bvec,nb0,nge,eig,ngw,wght0)
  !
  USE mpi, ONLY : MPI_IN_PLACE, MPI_DOUBLE_PRECISION, MPI_SUM
  USE libtetrabz_vals,    ONLY : ltetra, ng, nb, nk0, indx1, indx2, indx3
  USE libtetrabz_common,  ONLY : libtetrabz_initialize, libtetrabz_interpol_weight
  USE libtetrabz_occ_mod, ONLY : libtetrabz_occ_main
  !
  USE libtetrabz_routines, ONLY : comm, libtetrabz_kgrid
  !
  INTEGER,INTENT(IN) :: ltetra0, comm0, nge(3), ngw(3), nb0
  REAL(8),INTENT(IN) :: bvec(3,3)
  REAL(8),INTENT(IN) :: eig(nb0,PRODUCT(nge(1:3)))
  REAL(8),INTENT(OUT) :: wght0(nb0,PRODUCT(ngw(1:3)))
  !
  INTEGER :: ierr, nn
  REAL(8),ALLOCATABLE :: wght1(:,:)
  !
  ltetra = ltetra0
  comm = comm0
  nb = nb0
  ng(1:3) = nge(1:3)
  nn = nb
  !
  CALL libtetrabz_initialize(bvec)
  CALL libtetrabz_kgrid()
  !
  ALLOCATE(wght1(nn, nk0))
  CALL libtetrabz_occ_main(0d0,eig,wght1)
  !
  CALL libtetrabz_interpol_weight(nn,ngw,nge,wght0,wght1)
  !
  DEALLOCATE(wght1, indx1, indx2, indx3)
  !
  CALL MPI_allREDUCE(MPI_IN_PLACE, wght0, nn * PRODUCT(ngw(1:3)), &
  &                  MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)
  !
END SUBROUTINE libtetrabz_occ
!
! Calculate Fermi energy
!
SUBROUTINE libtetrabz_fermieng(ltetra0,comm0,bvec,nb0,nge,eig,ngw,wght0,ef,nelec)
  !
  USE mpi, ONLY : MPI_IN_PLACE, MPI_DOUBLE_PRECISION, MPI_SUM
  USE libtetrabz_vals,    ONLY : ltetra, ng, nb, nk, nk0, indx1, indx2, indx3
  USE libtetrabz_common,  ONLY : libtetrabz_initialize, libtetrabz_interpol_weight
  USE libtetrabz_occ_mod, ONLY : libtetrabz_occ_main
  !
  USE libtetrabz_routines, ONLY : comm, libtetrabz_kgrid
  !
  INTEGER,INTENT(IN) :: ltetra0, comm0, nge(3), ngw(3), nb0
  REAL(8),INTENT(IN) :: bvec(3,3), nelec
  REAL(8),INTENT(IN) :: eig(nb0,PRODUCT(nge(1:3)))
  REAL(8),INTENT(OUT) :: ef
  REAL(8),INTENT(OUT) :: wght0(nb0,PRODUCT(ngw(1:3)))
  !
  INTEGER :: ierr, nn
  REAL(8),ALLOCATABLE :: wght1(:,:)
  !
  INTEGER :: iter, maxiter = 300
  REAL(8) :: elw, eup, sumkmid, eps= 1.0d-10
  !
  ltetra = ltetra0
  comm = comm0
  nb = nb0
  ng(1:3) = nge(1:3)
  nn = nb
  !
  CALL libtetrabz_initialize(bvec)
  CALL libtetrabz_kgrid()
  !
  ALLOCATE(wght1(nn, nk0))
  !
  elw = MINVAL(eig(1:nb,1:nk))
  eup = MAXVAL(eig(1:nb,1:nk))
  !
  ! Bisection method
  !
  DO iter = 1, maxiter
     !
     ef = (eup + elw) / 2.d0
     !
     ! Calc. # of electrons 
     !
     CALL libtetrabz_occ_main(ef, eig,wght1)
     !
     sumkmid = SUM(wght1(1:nb,1:nk0))
     CALL MPI_allREDUCE(MPI_IN_PLACE, sumkmid, 1, &
     &                  MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)
     !
     ! convergence check
     !
     IF(ABS(sumkmid - nelec) < eps) THEN
        EXIT
     ELSE IF(sumkmid < nelec) THEN
        elw = ef
     ELSE
        eup = ef
     ENDIF
     !
  END DO  ! iter
  !
  IF(iter >= maxiter) STOP "libtetrabz_omp_fermieng"
  !
  CALL libtetrabz_interpol_weight(nn,ngw,nge,wght0,wght1)
  !
  DEALLOCATE(wght1, indx1, indx2, indx3)
  !
  CALL MPI_allREDUCE(MPI_IN_PLACE, wght0, nn * PRODUCT(ngw(1:3)), &
  &                  MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)
  !
END SUBROUTINE libtetrabz_fermieng
!
! Main SUBROUTINE for occupation : Theta(EF - E1)
!
SUBROUTINE libtetrabz_occ_main(ef,eig,occ)
  !
  USE libterabz_val, ONLY : nb, nk, nk0, fst, lst, wlsm, indx1, indx2
  IMPLICIT NONE
  !
  REAL(8),INTENT(IN) :: ef, eig(nb,nk)
  REAL(8),INTENT(OUT) :: occ(nb,nk0)
  !
  INTEGER :: it, ib, ii, sind(4)
  REAL(8) :: e(4), tsmall(4,4), V, w1(4), w2(4), ei1(4,nb)
  !
  occ(1:nb,1:nk0) = 0d0
  w2(1:4) = 0.25d0
  !
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP & SHARED(fst,lst,nb,wlsm,indx1,indx2,eig,ef,w2,occ) &
  !$OMP & PRIVATE(it,ii,ib,ei1,tsmall,e,V,w1)
  !
  DO it = fst, lst
     !
     DO ib = 1, nb
        ei1(1:4,ib) = MATMUL(wlsm(1:4,1:20) * eig(ib,indx1(1:20,it)))
     END DO
     !
     !$OMP DO
     DO ib = 1, nb
        !
        e(1:4) = ei1(1:4, ib) - ef
        CALL libtetrabz_sort(4,e,sind)
        w1(1:4) = 0d0
        !
        IF(e(1) <= 0d0 .AND. 0d0 < e(2)) THEN
           !
           CALL libtetrabz_tsmall_a1(e,V,tsmall)
           w1(sind(1:4)) = w1(sind(1:4)) &
           &  + V * MATMUL(w2(     1:4 ), tsmall(1:4,1:4))
           !
        ELSE IF(e(2) <= 0d0 .AND. 0d0 < e(3)) THEN
           !
           CALL libtetrabz_tsmall_a1(e,V,tsmall)
           w1(sind(1:4)) = w1(sind(1:4)) &
           &  + V * MATMUL(w2(     1:4 ), tsmall(1:4,1:4))
           !
           CALL libtetrabz_tsmall_b2(e,V,tsmall)
           w1(sind(1:4)) = w1(sind(1:4)) &
           &  + V * MATMUL(w2(     1:4 ), tsmall(1:4,1:4))
           !
           CALL libtetrabz_tsmall_b3(e,V,tsmall)
           w1(sind(1:4)) = w1(sind(1:4)) &
           &  + V * MATMUL(w2(     1:4 ), tsmall(1:4,1:4))
           !
        ELSE IF(e(3) <= 0d0 .AND. 0d0 < e(4)) THEN
           !
           CALL libtetrabz_tsmall_c1(e,V,tsmall)
           w1(sind(1:4)) = w1(sind(1:4)) &
           &  + V * MATMUL(w2(     1:4 ), tsmall(1:4,1:4))
           !
           CALL libtetrabz_tsmall_c2(e,V,tsmall)
           w1(sind(1:4)) = w1(sind(1:4)) &
           &  + V * MATMUL(w2(     1:4 ), tsmall(1:4,1:4))
           !
           CALL libtetrabz_tsmall_c3(e,V,tsmall)
           w1(sind(1:4)) = w1(sind(1:4)) &
           &  + V * MATMUL(w2(     1:4 ), tsmall(1:4,1:4))
           !
        ELSE IF(e(4) <= 0d0) THEN
           !
           w1(1:4) = w2(1:4)
           !
        ELSE
           !
           CYCLE
           !
        END IF
        !
        occ(ib,indx2(1:20,it)) = occ(ib,indx2(1:20,it)) &
        &                      + MATMUL(w1(1:4), wlsm(1:4,1:20))
        !
     END DO ! ib
     !$OMP END DO NOWAIT
     !
  END DO ! it
  !
  !$OMP END PARALLEL
  !
  occ(1:nb,1:nk0) = occ(1:nb,1:nk0) / DBLE(6 * nk)
  !
END SUBROUTINE libtetrabz_occ1
  !
END MODULE libterabz_occ_mod
