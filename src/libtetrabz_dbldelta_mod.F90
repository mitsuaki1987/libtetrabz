MODULE libtetrabz_doubledelta_mod
  !
  IMPLICIT NONE
  !
CONTAINS
!
! Compute doubledelta
!
SUBROUTINE libtetrabz_doubledelta(ltetra,comm0,bvec,nb,nge,eig1,eig2,ngw,wght0) BIND(C)
  !
#if defined(__MPI)
  USE mpi, ONLY : MPI_IN_PLACE, MPI_DOUBLE_PRECISION, MPI_SUM
#endif
  USE ISO_C_BINDING
  USE libtetrabz_val,   ONLY : nk_local, ik_global, ik_local, kvec, lmpi, linterpol, comm
  USE libtetrabz_common, ONLY : libtetrabz_initialize, libtetrabz_interpol_indx
  IMPLICIT NONE
  !
  INTEGER(C_INT),INTENT(IN) :: ltetra, nge(3), ngw(3), nb
  REAL(C_DOUBLE),INTENT(IN) :: bvec(3,3)
  REAL(C_DOUBLE),INTENT(IN) :: eig1(nb,PRODUCT(nge(1:3))), eig2(nb,PRODUCT(nge(1:3)))
  REAL(C_DOUBLE),INTENT(OUT) :: wght0(nb,nb,PRODUCT(ngw(1:3)))
  INTEGER(C_INT),INTENT(IN),OPTIONAL :: comm0
  !
  INTEGER :: ii, ik, kintp(4)
  REAL(8) :: wintp(4)
  REAL(8),ALLOCATABLE :: wght1(:,:,:)
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
     ALLOCATE(wght1(nb,nb,nk_local))
     CALL libtetrabz_doubledelta_main(eig1,eig2,wght1)
     !
     ! Interpolation
     !
     wght0(1:nb,1:nb,1:PRODUCT(ngw(1:3))) = 0d0
     DO ik = 1, nk_local
        CALL libtetrabz_interpol_indx(ngw,kvec(1:3,ik),kintp,wintp)
        DO ii = 1, 4
           wght0(1:nb,1:nb,kintp(ii)) = wght0(1:nb,1:nb,       kintp(ii)) &
           &                          + wght1(1:nb,1:nb, ik) * wintp(ii)
        END DO
     END DO ! ik = 1, nk_local
     DEALLOCATE(wght1, kvec)
     !
#if defined(__MPI)
     IF(lmpi) &
     &  CALL MPI_allREDUCE(MPI_IN_PLACE, wght0, nb * nb * PRODUCT(ngw(1:3)), &
     &                     MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)
#endif
     !
  ELSE
     CALL libtetrabz_doubledelta_main(eig1,eig2,wght0)
  END IF
  !
  DEALLOCATE(ik_global, ik_local)
  !
END SUBROUTINE libtetrabz_doubledelta
!
! Main SUBROUTINE for Delta(E1) * Delta(E2)
!
SUBROUTINE libtetrabz_doubledelta_main(eig1,eig2,dbldelta)
  !
  USE libtetrabz_val, ONLY : nb, nk_local, nt_local, wlsm, ik_global, ik_local, nkBZ
  IMPLICIT NONE
  !
  REAL(8),INTENT(IN) :: eig1(nb,nkBZ), eig2(nb,nkBZ)
  REAL(8),INTENT(OUT) :: dbldelta(nb,nb,nk_local)
  !
  INTEGER :: it, ib, indx(4)
  REAL(8) :: e(4), V, thr = 1d-10, ei1(4,nb), ei2(3), ej1(4,nb), ej2(3,nb), &
  &          w1(nb,4), w2(nb,3), tsmall(3,4)
  !
  dbldelta(1:nb,1:nb,1:nk_local) = 0d0
  !
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP & SHARED(nt_local,nb,ik_global,ik_local,wlsm,eig1,eig2,dbldelta) &
  !$OMP & PRIVATE(ib,it,e,w1,w2,ei1,ej1,ei2,ej2,V,indx)
  !
  DO it = 1, nt_local
     !
     DO ib = 1, nb
        ei1(1:4,ib) = MATMUL(wlsm(1:4,1:20), eig1(ib, ik_global(1:20,it)))
        ej1(1:4,ib) = MATMUL(wlsm(1:4,1:20), eig2(ib, ik_global(1:20,it)))
     END DO
     !
     !$OMP DO
     DO ib = 1, nb
        !
        w1(1:nb,1:4) = 0d0
        e(1:4) = ei1(1:4, ib)
        CALL libtetrabz_sort(4,e,indx)
        !
        IF(e(1) < 0d0 .AND. 0d0 <= e(2)) THEN
           !
           CALL libtetrabz_tsmall_a1(e,V,tsmall)
           !
           IF(V > thr) THEN
              !
              ej2(1:3,1:nb) = MATMUL(tsmall(1:3,1:4), ej1(indx(1:4),1:nb))
              CALL libtetrabz_doubledelta2(ej2,w2)
              w1(1:nb,indx(1:4)) = w1(1:nb,            indx(1:4)) &
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
              ej2(1:3,1:nb) = MATMUL(tsmall(1:3,1:4), ej1(indx(1:4),1:nb))
              CALL libtetrabz_doubledelta2(ej2,w2)
              w1(1:nb,indx(1:4)) = w1(1:nb,            indx(1:4)) &
              &       + V * MATMUL(w2(1:nb,1:3), tsmall(1:3,1:4))
              !
           END IF
           !
           CALL libtetrabz_tsmall_b2(e,V,tsmall)
           !
           IF(V > thr) THEN
              !
              ej2(1:3,1:nb) = MATMUL(tsmall(1:3,1:4), ej1(indx(1:4),1:nb))
              CALL libtetrabz_doubledelta2(ej2,w2)
              w1(1:nb,indx(1:4)) = w1(1:nb,            indx(1:4)) &
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
              ej2(1:3,1:nb) = MATMUL(tsmall(1:3,1:4), ej1(indx(1:4),1:nb))
              CALL libtetrabz_doubledelta2(ej2,w2)
              w1(1:nb,indx(1:4)) = w1(1:nb,            indx(1:4)) &
              &       + V * MATMUL(w2(1:nb,1:3), tsmall(1:3,1:4))
              !
           END IF
           !
        END IF
        !
        dbldelta(1:nb,ib,ik_local(1:20,it)) = dbldelta(1:nb,ib,   ik_local(1:20,it)) &
        &                                  + MATMUL(w1(1:nb,1:4), wlsm(1:4,1:20))
        !
     END DO ! ib
     !$OMP END DO NOWAIT
     !
  END DO ! it
  !
  !$OMP END PARALLEL
  !
  dbldelta(1:nb,1:nb,1:nk_local) = dbldelta(1:nb,1:nb,1:nk_local) / DBLE(6 * nkBZ)
  !
END SUBROUTINE libtetrabz_doubledelta_main
!
! 2nd step of tetrahedra method.
!
SUBROUTINE libtetrabz_doubledelta2(ej,w)
  !
  USE libtetrabz_val, ONLY : nb
  IMPLICIT NONE
  !
  REAL(8),INTENT(IN) :: ej(3,nb)
  REAL(8),INTENT(INOUT) :: w(nb,3)
  !
  INTEGER :: ib, ii, indx(3)
  REAL(8) :: e(3), V, a(3,3)
  !
  DO ib = 1, nb
     !
     IF(maxval(ABS(ej(ib,1:3))) < 1d-10) STOP "Nesting !!"
     !
     w(ib,1:3) = 0d0
     e(1:3) = ej(1:3, ib)
     CALL libtetrabz_sort(3,e,indx)
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
        w(ib,indx(1)) = V * (a(1,2) + a(1,3))
        w(ib,indx(2)) = V * a(2,1)
        w(ib,indx(3)) = V * a(3,1)
        !
     ELSE IF((e(2) <= 0d0 .AND. 0d0 < e(3)) .OR. (e(2) < 0d0 .AND. 0d0 <= e(3))) THEN
        !
        !V = a(1,3) * a(2,3) / (e(3) - 0d0) 
        V = a(1,3)           / (e(3) - e(2)) 
        !
        w(ib,indx(1)) = V * a(1,3)
        w(ib,indx(2)) = V * a(2,3)
        w(ib,indx(3)) = V * (a(3,1) + a(3,2))
        !
     END IF
     !
  END DO ! ib
  !
END SUBROUTINE libtetrabz_doubledelta2
!
END MODULE libtetrabz_doubledelta_mod
