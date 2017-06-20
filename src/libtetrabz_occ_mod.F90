MODULE libtetrabz_occ_mod
  !
  IMPLICIT NONE
  !
CONTAINS
!
! Compute occupation
!
SUBROUTINE libtetrabz_occ(ltetra,bvec,nb,nge,eig,ngw,wght0,comm0) BIND(C)
  !
#if defined(__MPI)
  USE mpi, ONLY : MPI_IN_PLACE, MPI_DOUBLE_PRECISION, MPI_SUM
#endif
  USE ISO_C_BINDING
  USE libtetrabz_val,   ONLY : nk_local, ik_global, ik_local, kvec, lmpi, linterpol, comm
  USE libtetrabz_common, ONLY : libtetrabz_initialize, libtetrabz_interpol_indx
  !
  INTEGER(C_INT),INTENT(IN) :: ltetra, nge(3), ngw(3), nb
  REAL(C_DOUBLE),INTENT(IN) :: bvec(3,3)
  REAL(C_DOUBLE),INTENT(IN) :: eig(nb,PRODUCT(nge(1:3)))
  REAL(C_DOUBLE),INTENT(OUT) :: wght0(nb,PRODUCT(ngw(1:3)))
  INTEGER(C_INT),INTENT(IN),OPTIONAL :: comm0
  !
  INTEGER :: ii, ik, kintp(4)
  REAL(8) :: wintp(4)
  REAL(8),ALLOCATABLE :: wght1(:,:)
#if defined(__MPI)
  INTEGER :: ierr
#endif
  !
  lmpi = .FALSE.
  IF(PRESENT(comm0)) THEN
     comm = comm0
#if defined(__MPI)
     lmpi = .TRUE.
#endif
  END IF
  !
  CALL libtetrabz_initialize(ltetra,bvec,nge,ngw,nb)
  !
  IF(lmpi .OR. linterpol) THEN
     !
     ALLOCATE(wght1(nb,nk_local))
     CALL libtetrabz_occ_main(0d0,eig,wght1)
     !
     ! Interpolation
     !
     wght0(1:nb,1:PRODUCT(ngw(1:3))) = 0d0
     DO ik = 1, nk_local
        CALL libtetrabz_interpol_indx(ngw,kvec(1:3,ik),kintp,wintp)
        DO ii = 1, 4
           wght0(1:nb,kintp(ii)) = wght0(1:nb,       kintp(ii)) &
           &                     + wght1(1:nb, ik) * wintp(ii)
        END DO
     END DO ! ik = 1, nk_local
     DEALLOCATE(wght1, kvec)
     !
#if defined(__MPI)
     IF(lmpi) &
     &  CALL MPI_allREDUCE(MPI_IN_PLACE, wght0, nb * PRODUCT(ngw(1:3)), &
     &                     MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)
#endif
     !
  ELSE
     CALL libtetrabz_occ_main(0d0,eig,wght0)
  END IF
  !
  DEALLOCATE(ik_global, ik_local)
  !
END SUBROUTINE libtetrabz_occ
!
! Calculate Fermi energy
!
SUBROUTINE libtetrabz_fermieng(ltetra,bvec,nb,nge,eig,ngw,wght0,ef,nelec,comm0) BIND(C)
  !
#if defined(__MPI)
  USE mpi, ONLY : MPI_IN_PLACE, MPI_DOUBLE_PRECISION, MPI_SUM
#endif
  USE ISO_C_BINDING
  USE libtetrabz_val,    ONLY : nk_local, ik_global, ik_local, kvec, comm, lmpi, linterpol
  USE libtetrabz_common,  ONLY : libtetrabz_initialize, libtetrabz_interpol_indx
  !
  INTEGER(C_INT),INTENT(IN) :: ltetra, nge(3), ngw(3), nb
  REAL(8),INTENT(IN) :: bvec(3,3), nelec
  REAL(8),INTENT(IN) :: eig(nb,PRODUCT(nge(1:3)))
  REAL(8),INTENT(OUT) :: ef
  REAL(8),INTENT(OUT) :: wght0(nb,PRODUCT(ngw(1:3)))
  INTEGER(C_INT),INTENT(IN),OPTIONAL :: comm0
  !
  INTEGER :: ii, ik, kintp(4)
  REAL(8) :: wintp(4)
  REAL(8),ALLOCATABLE :: wght1(:,:)
#if defined(__MPI)
  INTEGER :: ierr
#endif
  !
  INTEGER :: iter, maxiter = 300
  REAL(8) :: elw, eup, sumkmid, eps= 1.0d-10
  !
  lmpi = .FALSE.
  IF(PRESENT(comm0)) THEN
     comm = comm0
#if defined(__MPI)
     lmpi = .TRUE.
#endif
  END IF
  !
  CALL libtetrabz_initialize(ltetra,bvec,nge,ngw,nb)
  !
  IF(lmpi .OR. linterpol) ALLOCATE(wght1(nb, nk_local))
  !
  elw = MINVAL(eig(1:nb,1:PRODUCT(nge(1:3))))
  eup = MAXVAL(eig(1:nb,1:PRODUCT(nge(1:3))))
  !
  ! Bisection method
  !
  DO iter = 1, maxiter
     !
     ef = (eup + elw) / 2.d0
     !
     ! Calc. # of electrons 
     !
     IF(lmpi .AND. linterpol) THEN
        CALL libtetrabz_occ_main(ef, eig,wght1)
        sumkmid = SUM(wght1(1:nb,1:nk_local))
        !
#if defined(__MPI)
        IF(lmpi) &
        &  CALL MPI_allREDUCE(MPI_IN_PLACE, sumkmid, 1, &
        &                     MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)
#endif
        !
     ELSE
        CALL libtetrabz_occ_main(ef, eig,wght0)
        sumkmid = SUM(wght0(1:nb,1:nk_local))
     END IF
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
  IF(lmpi .OR. linterpol) THEN
     !
     ! Interpolation
     !
     wght0(1:nb,1:PRODUCT(ngw(1:3))) = 0d0
     DO ik = 1, nk_local
        CALL libtetrabz_interpol_indx(ngw,kvec(1:3,ik),kintp,wintp)
        DO ii = 1, 4
           wght0(1:nb,kintp(ii)) = wght0(1:nb,       kintp(ii)) &
           &                     + wght1(1:nb, ik) * wintp(ii)
        END DO
     END DO ! ik = 1, nk_local
     DEALLOCATE(wght1, kvec)
     !
#if defined(__MPI)
     IF(lmpi) &
     &  CALL MPI_allREDUCE(MPI_IN_PLACE, wght0, nb * PRODUCT(ngw(1:3)), &
     &                     MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)
#endif
     !
  END IF ! (lmpi .OR. linterpol)
  !
  DEALLOCATE(ik_global, ik_local)
  !
END SUBROUTINE libtetrabz_fermieng
!
! Main SUBROUTINE for occupation : Theta(EF - E1)
!
SUBROUTINE libtetrabz_occ_main(ef,eig,occ)
  !
  USE libtetrabz_val, ONLY : nb, nk_local, nt_local, wlsm, ik_global, ik_local, nkBZ
  IMPLICIT NONE
  !
  REAL(8),INTENT(IN) :: ef, eig(nb,nkBZ)
  REAL(8),INTENT(OUT) :: occ(nb,nk_local)
  !
  INTEGER :: it, ib, indx(4)
  REAL(8) :: e(4), tsmall(4,4), V, w1(4), w2(4), ei1(4,nb)
  !
  occ(1:nb,1:nk_local) = 0d0
  w2(1:4) = 0.25d0
  !
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP & SHARED(nt_local,nb,wlsm,ik_local,ik_global,eig,ef,w2,occ) &
  !$OMP & PRIVATE(it,ib,ei1,tsmall,e,V,w1,indx)
  !
  DO it = 1, nt_local
     !
     DO ib = 1, nb
        ei1(1:4,ib) = MATMUL(wlsm(1:4,1:20), eig(ib,ik_global(1:20,it)))
     END DO
     !
     !$OMP DO
     DO ib = 1, nb
        !
        e(1:4) = ei1(1:4, ib) - ef
        CALL libtetrabz_sort(4,e,indx)
        w1(1:4) = 0d0
        !
        IF(e(1) <= 0d0 .AND. 0d0 < e(2)) THEN
           !
           CALL libtetrabz_tsmall_a1(e,V,tsmall)
           w1(indx(1:4)) = w1(indx(1:4)) &
           &  + V * MATMUL(w2(     1:4 ), tsmall(1:4,1:4))
           !
        ELSE IF(e(2) <= 0d0 .AND. 0d0 < e(3)) THEN
           !
           CALL libtetrabz_tsmall_a1(e,V,tsmall)
           w1(indx(1:4)) = w1(indx(1:4)) &
           &  + V * MATMUL(w2(     1:4 ), tsmall(1:4,1:4))
           !
           CALL libtetrabz_tsmall_b2(e,V,tsmall)
           w1(indx(1:4)) = w1(indx(1:4)) &
           &  + V * MATMUL(w2(     1:4 ), tsmall(1:4,1:4))
           !
           CALL libtetrabz_tsmall_b3(e,V,tsmall)
           w1(indx(1:4)) = w1(indx(1:4)) &
           &  + V * MATMUL(w2(     1:4 ), tsmall(1:4,1:4))
           !
        ELSE IF(e(3) <= 0d0 .AND. 0d0 < e(4)) THEN
           !
           CALL libtetrabz_tsmall_c1(e,V,tsmall)
           w1(indx(1:4)) = w1(indx(1:4)) &
           &  + V * MATMUL(w2(     1:4 ), tsmall(1:4,1:4))
           !
           CALL libtetrabz_tsmall_c2(e,V,tsmall)
           w1(indx(1:4)) = w1(indx(1:4)) &
           &  + V * MATMUL(w2(     1:4 ), tsmall(1:4,1:4))
           !
           CALL libtetrabz_tsmall_c3(e,V,tsmall)
           w1(indx(1:4)) = w1(indx(1:4)) &
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
        occ(ib,ik_local(1:20,it)) = occ(ib,ik_local(1:20,it)) &
        &                + MATMUL(w1(1:4), wlsm(1:4,1:20))
        !
     END DO ! ib
     !$OMP END DO NOWAIT
     !
  END DO ! it
  !
  !$OMP END PARALLEL
  !
  occ(1:nb,1:nk_local) = occ(1:nb,1:nk_local) / DBLE(6 * nkBZ)
  !
END SUBROUTINE libtetrabz_occ_main
  !
END MODULE libtetrabz_occ_mod
