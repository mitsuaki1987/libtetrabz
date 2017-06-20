MODULE libtetrabz_polstat_mod
  !
  IMPLICIT NONE
  !
CONTAINS
!
! Compute Static polalization function
!
SUBROUTINE libtetrabz_polstat(ltetra,comm0,bvec,nb,nge,eig1,eig2,ngw,wght0) BIND(C)
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
     CALL libtetrabz_polstat_main(eig1,eig2,wght1)
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
     CALL libtetrabz_polstat_main(eig1,eig2,wght0)
  END IF
  !
  DEALLOCATE(ik_global, ik_local)
  !
END SUBROUTINE libtetrabz_polstat
!
! Main SUBROUTINE for polalization function : Theta(- E1) * Theta(E2) / (E2 - E1)
!
SUBROUTINE libtetrabz_polstat_main(eig1,eig2,pols)
  !
  USE libtetrabz_val, ONLY : nb, nk_local, nt_local, wlsm, ik_global, ik_local, nkBZ
  IMPLICIT NONE
  !
  REAL(8),INTENT(IN) :: eig1(nb,nkBZ), eig2(nb,nkBZ)
  REAL(8),INTENT(OUT) :: pols(nb,nb,nk_local)
  !
  INTEGER :: it, ib, indx(4)
  REAL(8) :: e(4), a(4,4), V, thr = 1d-10, tsmall(4,4), &
  &          ei1(4,nb), ei2(4), ej1(4,nb), ej2(4,nb), &
  &          w1(nb,4), w2(nb,4)
  !
  pols(1:nb,1:nb,1:nk_local) = 0d0
  !
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP & SHARED(nt_local,nb,ik_global,ik_local,wlsm,eig1,eig2,pols,thr) &
  !$OMP & PRIVATE(ib,it,e,w1,w2,ei1,ei2,ej1,ej2,V,indx)
  !
  DO it = 1, nt_local
     !
     DO ib = 1, nb
        ei1(1:4, ib) = MATMUL(wlsm(1:4,1:20), eig1(ib, ik_global(1:20,it)))
        ej1(1:4, ib) = MATMUL(wlsm(1:4,1:20), eig2(ib, ik_global(1:20,it)))
     END DO
     !
     !$OMP DO
     DO ib = 1, nb
        !
        w1(1:nb,1:4) = 0d0
        e(1:4) = ei1(1:4, ib)
        CALL libtetrabz_sort(4,e,indx)
        !
        IF(e(1) <= 0d0 .AND. 0d0 < e(2)) THEN
           !
           CALL libtetrabz_tsmall_a1(e,V,tsmall)
           !
           IF(V > thr) THEN
              !
              ei2(1:4     ) = MATMUL(tsmall(1:4,1:4), ei1(indx(1:4),  ib))
              ej2(1:4,1:nb) = MATMUL(tsmall(1:4,1:4), ej1(indx(1:4),1:nb))
              CALL libtetrabz_polstat2(ei2,ej2,w2)
              w1(1:nb,indx(1:4)) = w1(1:nb,            indx(1:4)) &
              &       + V * MATMUL(w2(1:nb,1:4), tsmall(1:4,1:4))
              !
           END IF
           !
        ELSE IF(e(2) <= 0d0 .AND. 0d0 < e(3)) THEN
           !
           CALL libtetrabz_tsmall_b1(e,V,tsmall)
           !
           IF(V > thr) THEN
              !
              ei2(1:4     ) = MATMUL(tsmall(1:4,1:4), ei1(indx(1:4),  ib))
              ej2(1:4,1:nb) = MATMUL(tsmall(1:4,1:4), ej1(indx(1:4),1:nb))
              CALL libtetrabz_polstat2(ei2,ej2,w2)
              w1(1:nb,indx(1:4)) = w1(1:nb,            indx(1:4)) &
              &       + V * MATMUL(w2(1:nb,1:4), tsmall(1:4,1:4))
              !
           END IF
           !
           CALL libtetrabz_tsmall_b2(e,V,tsmall)
           !
           IF(V > thr) THEN
              !
              ei2(1:4     ) = MATMUL(tsmall(1:4,1:4), ei1(indx(1:4),  ib))
              ej2(1:4,1:nb) = MATMUL(tsmall(1:4,1:4), ej1(indx(1:4),1:nb))
              CALL libtetrabz_polstat2(ei2,ej2,w2)
              w1(1:nb,indx(1:4)) = w1(1:nb,            indx(1:4)) &
              &       + V * MATMUL(w2(1:nb,1:4), tsmall(1:4,1:4))
              !
           END IF
           !
           CALL libtetrabz_tsmall_b3(e,V,tsmall)
           !
           IF(V > thr) THEN
              !
              ei2(1:4     ) = MATMUL(tsmall(1:4,1:4), ei1(indx(1:4),  ib))
              ej2(1:4,1:nb) = MATMUL(tsmall(1:4,1:4), ej1(indx(1:4),1:nb))
              CALL libtetrabz_polstat2(ei2,ej2,w2)
              w1(1:nb,indx(1:4)) = w1(1:nb,            indx(1:4)) &
              &       + V * MATMUL(w2(1:nb,1:4), tsmall(1:4,1:4))
              !
           END IF
           !
        ELSE IF(e(3) <= 0d0 .AND. 0d0 < e(4)) THEN
           !
           CALL libtetrabz_tsmall_c1(e,V,tsmall)
           !
           IF(V > thr) THEN
              !
              ei2(1:4     ) = MATMUL(tsmall(1:4,1:4), ei1(indx(1:4),  ib))
              ej2(1:4,1:nb) = MATMUL(tsmall(1:4,1:4), ej1(indx(1:4),1:nb))
              CALL libtetrabz_polstat2(ei2,ej2,w2)
              w1(1:nb,indx(1:4)) = w1(1:nb,            indx(1:4)) &
              &       + V * MATMUL(w2(1:nb,1:4), tsmall(1:4,1:4))
              !
           END IF
           !
           CALL libtetrabz_tsmall_c2(e,V,tsmall)
           !
           IF(V > thr) THEN
              !
              ei2(1:4     ) = MATMUL(tsmall(1:4,1:4), ei1(indx(1:4),  ib))
              ej2(1:4,1:nb) = MATMUL(tsmall(1:4,1:4), ej1(indx(1:4),1:nb))
              CALL libtetrabz_polstat2(ei2,ej2,w2)
              w1(1:nb,indx(1:4)) = w1(1:nb,            indx(1:4)) &
              &       + V * MATMUL(w2(1:nb,1:4), tsmall(1:4,1:4))
              !
           END IF
           !
           CALL libtetrabz_tsmall_c3(e,V,tsmall)
           !
           IF(V > thr) THEN
              !
              ei2(1:4     ) = MATMUL(tsmall(1:4,1:4), ei1(indx(1:4),  ib))
              ej2(1:4,1:nb) = MATMUL(tsmall(1:4,1:4), ej1(indx(1:4),1:nb))
              CALL libtetrabz_polstat2(ei2,ej2,w2)
              w1(1:nb,indx(1:4)) = w1(1:nb,            indx(1:4)) &
              &       + V * MATMUL(w2(1:nb,1:4), tsmall(1:4,1:4))
              !
           END IF
           !
        ELSE IF(e(4) <= 0d0) THEN
           !
           ei2(1:4     ) = ei1(1:4,  ib)
           ej2(1:4,1:nb) = ej1(1:4,1:nb)
           CALL libtetrabz_polstat2(ei2,ej2,w2)
           w1(1:nb,1:4) = w1(1:nb,1:4) + w2(1:nb,1:4)
           !
        ELSE
           !
           CYCLE
           !
        END IF
        !
        pols(1:nb,ib,ik_local(1:20,it)) = pols(1:nb,ib,   ik_local(1:20,it)) &
        &                          + MATMUL(w1(1:nb,1:4), wlsm(1:4,1:20))
        !
     END DO ! ib
     !$OMP END DO NOWAIT
     !
  END DO ! it
  !
  !$OMP END PARALLEL
  !
  pols(1:nb,1:nb,1:nk_local) = pols(1:nb,1:nb,1:nk_local) / DBLE(6 * nkBZ)
  !
END SUBROUTINE libtetrabz_polstat_main
!
! Tetrahedra method for theta( - E2)
!
SUBROUTINE libtetrabz_polstat2(ei1,ej1,w1)
  !
  USE libtetrabz_val, ONLY : nb
  IMPLICIT NONE
  !
  REAL(8),INTENT(IN) :: ei1(4), ej1(4,nb)
  REAL(8),INTENT(INOUT) :: w1(nb,4)
  !
  INTEGER :: ib, indx(4)
  REAL(8) :: V, w2(4), thr = 1d-8, &
  &          e(4), ei2(4), ej2(4), tsmall(4,4)
  !
  DO ib = 1, nb
     !
     w1(ib,1:4) = 0d0
     e(1:4) = - ej1(1:4, ib)
     CALL libtetrabz_sort(4,e,indx)
     !
     IF((e(1) <= 0d0 .AND. 0d0 < e(2)) .OR. (e(1) < 0d0 .AND. 0d0 <= e(2))) THEN
        !
        CALL libtetrabz_tsmall_a1(e,V,tsmall)
        !
        IF(V > thr) THEN
           !
           ei2(1:4) = MATMUL(tsmall(1:4,1:4), ei1(indx(1:4)   ))
           ej2(1:4) = MATMUL(tsmall(1:4,1:4), ej1(indx(1:4),ib))
           CALL libtetrabz_polstat3(ei2,ej2,w2)
           w1(ib,indx(1:4)) = w1(ib,                    indx(1:4)) &
           &                + V * MATMUL(w2(1:4), tsmall(1:4,1:4))
           !
        END IF
        !
     ELSE IF((e(2) <= 0d0 .AND. 0d0 < e(3)) .OR. (e(2) < 0d0 .AND. 0d0 <= e(3))) THEN
        !
        CALL libtetrabz_tsmall_b1(e,V,tsmall)
        !
        IF(V > thr) THEN
           !
           ei2(1:4) = MATMUL(tsmall(1:4,1:4), ei1(indx(1:4)   ))
           ej2(1:4) = MATMUL(tsmall(1:4,1:4), ej1(indx(1:4),ib))
           CALL libtetrabz_polstat3(ei2,ej2,w2)
           w1(ib,indx(1:4)) = w1(ib,indx(1:4)) &
           &          + V * MATMUL(w2(        1:4 ), tsmall(1:4,1:4))
           !
        END IF
        !
        CALL libtetrabz_tsmall_b2(e,V,tsmall)
        !
        IF(V > thr) THEN
           !
           ei2(1:4) = MATMUL(tsmall(1:4,1:4), ei1(indx(1:4)   ))
           ej2(1:4) = MATMUL(tsmall(1:4,1:4), ej1(indx(1:4),ib))
           CALL libtetrabz_polstat3(ei2,ej2,w2)
           w1(ib,indx(1:4)) = w1(ib,                    indx(1:4)) &
           &                + V * MATMUL(w2(1:4), tsmall(1:4,1:4))
           !
        END IF
        !
        CALL libtetrabz_tsmall_b3(e,V,tsmall)
        !
        IF(V > thr) THEN
           !
           ei2(1:4) = MATMUL(tsmall(1:4,1:4), ei1(indx(1:4)   ))
           ej2(1:4) = MATMUL(tsmall(1:4,1:4), ej1(indx(1:4),ib))
           CALL libtetrabz_polstat3(ei2,ej2,w2)
           w1(ib,indx(1:4)) = w1(ib,                    indx(1:4)) &
           &                + V * MATMUL(w2(1:4), tsmall(1:4,1:4))
           !
        END IF
        !
     ELSE IF((e(3) <= 0d0 .AND. 0d0 < e(4)) .OR. (e(3) < 0d0 .AND. 0d0 <= e(4))) THEN
        !
        CALL libtetrabz_tsmall_c1(e,V,tsmall)
        !
        IF(V > thr) THEN
           !
           ei2(1:4) = MATMUL(tsmall(1:4,1:4), ei1(indx(1:4)   ))
           ej2(1:4) = MATMUL(tsmall(1:4,1:4), ej1(indx(1:4),ib))
           CALL libtetrabz_polstat3(ei2,ej2,w2)
           w1(ib,indx(1:4)) = w1(ib,                    indx(1:4)) &
           &                + V * MATMUL(w2(1:4), tsmall(1:4,1:4))
           !
        END IF
        !
        CALL libtetrabz_tsmall_c2(e,V,tsmall)
        !
        IF(V > thr) THEN
           !
           ei2(1:4) = MATMUL(tsmall(1:4,1:4), ei1(indx(1:4)   ))
           ej2(1:4) = MATMUL(tsmall(1:4,1:4), ej1(indx(1:4),ib))
           CALL libtetrabz_polstat3(ei2,ej2,w2)
           w1(ib,indx(1:4)) = w1(ib,                    indx(1:4)) &
           &                + V * MATMUL(w2(1:4), tsmall(1:4,1:4))
           !
        END IF
        !
        CALL libtetrabz_tsmall_c3(e,V,tsmall)
        !
        IF(V > thr) THEN
           !
           ei2(1:4) = MATMUL(tsmall(1:4,1:4), ei1(indx(1:4)   ))
           ej2(1:4) = MATMUL(tsmall(1:4,1:4), ej1(indx(1:4),ib))
           CALL libtetrabz_polstat3(ei2,ej2,w2)
           w1(ib,indx(1:4)) = w1(ib,                    indx(1:4)) &
           &                + V * MATMUL(w2(1:4), tsmall(1:4,1:4))
           !
        END IF
        !
     ELSE IF( e(4) <= 0d0 ) THEN
        !
        ei2(1:4) = ei1(1:4)
        ej2(1:4) = ej1(1:4,ib)
        CALL libtetrabz_polstat3(ei2,ej2,w2)
        w1(ib,1:4) = w1(ib,1:4) + V * w2(1:4)
        !
     END IF
     !
  END DO ! ib = 1, nb
  !
END SUBROUTINE libtetrabz_polstat2
!
! Tetarahedra method for delta(om - ep + e)
!
SUBROUTINE libtetrabz_polstat3(ei1,ej1,w1)
  !
  IMPLICIT NONE
  !
  REAL(8),INTENT(IN) :: ei1(4), ej1(4)
  REAL(8),INTENT(INOUT) :: w1(4)
  !
  INTEGER :: ii, indx(4)
  REAL(8) :: e(4), ln(4), thr, thr2
  !
  e(1:4) = ej1(1:4) - ei1(1:4)
  CALL libtetrabz_sort(4,e,indx)
  !
  thr = MAXVAL(e(1:4)) * 1d-3
  thr2 = 1d-8
  !
  DO ii = 1, 4
     IF(e(ii) < thr2) THEN
        IF(ii == 3) THEN
           STOP "  Nesting ! "
        END IF
        ln(ii) = 0d0
        e(ii) = 0d0
     ELSE
        ln(ii) = LOG(e(ii))
     END IF
  END DO
  !
  IF(ABS(e(4) - e(3)) < thr ) THEN
     IF(ABS(e(4) - e(2)) < thr ) THEN
        IF(ABS(e(4) - e(1)) < thr ) THEN
           !
           ! e(4) = e(3) = e(2) = e(1)
           !
           w1(indx(4)) = 0.25d0 / e(4)
           w1(indx(3)) = w1(indx(4))
           w1(indx(2)) = w1(indx(4))
           w1(indx(1)) = w1(indx(4))
           !
        ELSE
           !
           ! e(4) = e(3) = e(2)
           !
           w1(indx(4)) = libtetrabz_polstat_1211(e(4),e(1),ln(4),ln(1))
           w1(indx(3)) = w1(indx(4))
           w1(indx(2)) = w1(indx(4))
           w1(indx(1)) = libtetrabz_polstat_1222(e(1),e(4),ln(1),ln(4))
           !
           IF(ANY(w1(1:4) < 0d0)) THEN
              WRITE(*,'(100e15.5)') e(1:4)
              WRITE(*,'(100e15.5)') w1(indx(1:4))
              STOP "weighting 4=3=2"
           END IF
           !
        END IF
     ELSE IF(ABS(e(2) - e(1)) < thr) THEN
        !
        ! e(4) = e(3), e(2) = e(1)
        !
        w1(indx(4)) = libtetrabz_polstat_1221(e(4),e(2), ln(4),ln(2))
        w1(indx(3)) = w1(indx(4))
        w1(indx(2)) = libtetrabz_polstat_1221(e(2),e(4), ln(2),ln(4))
        w1(indx(1)) = w1(indx(2))
        !
        IF(ANY(w1(1:4) < 0d0)) THEN
           WRITE(*,'(100e15.5)') e(1:4)
           WRITE(*,'(100e15.5)') w1(indx(1:4))
           STOP "weighting 4=3 2=1"
        END IF
        !
     ELSE
        !
        ! e(4) = e(3)
        !
        w1(indx(4)) = libtetrabz_polstat_1231(e(4),e(1),e(2),ln(4),ln(1),ln(2))
        w1(indx(3)) = w1(indx(4))
        w1(indx(2)) = libtetrabz_polstat_1233(e(2),e(1),e(4),ln(2),ln(1),ln(4))
        w1(indx(1)) = libtetrabz_polstat_1233(e(1),e(2),e(4),ln(1),ln(2),ln(4))
        !
        IF(ANY(w1(1:4) < 0d0)) THEN
           WRITE(*,'(100e15.5)') e(1:4)
           WRITE(*,'(100e15.5)') w1(indx(1:4))
           STOP "weighting 4=3"
        END IF
        !
     END IF
  ELSE IF(ABS(e(3) - e(2)) < thr) THEN
     IF(ABS(e(3) - e(1)) < thr) THEN
        !
        ! e(3) = e(2) = e(1)
        !
        w1(indx(4)) = libtetrabz_polstat_1222(e(4),e(3), ln(4),ln(3))
        w1(indx(3)) = libtetrabz_polstat_1211(e(3),e(4), ln(3),ln(4))
        w1(indx(2)) = w1(indx(3))
        w1(indx(1)) = w1(indx(3))
        !
        IF(ANY(w1(1:4) < 0d0)) THEN
           WRITE(*,'(100e15.5)') e(1:4)
           WRITE(*,'(100e15.5)') w1(indx(1:4))
           STOP "weighting 3=2=1"
        END IF
        !
     ELSE
        !
        ! e(3) = e(2)
        !
        w1(indx(4)) = libtetrabz_polstat_1233(e(4),e(1),e(3),ln(4),ln(1),ln(3))
        w1(indx(3)) = libtetrabz_polstat_1231(e(3),e(1),e(4),ln(3),ln(1),ln(4))
        w1(indx(2)) = w1(indx(3))
        w1(indx(1)) = libtetrabz_polstat_1233(e(1),e(4),e(3),ln(1),ln(4),ln(3))
        !
        IF(ANY(w1(1:4) < 0d0)) THEN
           WRITE(*,'(100e15.5)') e(1:4)
           WRITE(*,'(100e15.5)') w1(indx(1:4))
           STOP "weighting 3=2"
        END IF
        !
     END IF
  ELSE IF(ABS(e(2) - e(1)) < thr) THEN
     !
     ! e(2) = e(1)
     !
     w1(indx(4)) = libtetrabz_polstat_1233(e(4),e(3),e(2),ln(4),ln(3),ln(2))
     w1(indx(3)) = libtetrabz_polstat_1233(e(3),e(4),e(2),ln(3),ln(4),ln(2))
     w1(indx(2)) = libtetrabz_polstat_1231(e(2),e(3),e(4),ln(2),ln(3),ln(4))
     w1(indx(1)) = w1(indx(2))
     !
     IF(ANY(w1(1:4) < 0d0)) THEN
        WRITE(*,'(100e15.5)') e(1:4)
        WRITE(*,'(100e15.5)') w1(indx(1:4))
        STOP "weighting 2=1"
     END IF
     !
  ELSE
     !
     ! Different each other.
     !
     w1(indx(4)) = libtetrabz_polstat_1234(e(4),e(1),e(2),e(3),ln(4),ln(1),ln(2),ln(3))
     w1(indx(3)) = libtetrabz_polstat_1234(e(3),e(1),e(2),e(4),ln(3),ln(1),ln(2),ln(4))
     w1(indx(2)) = libtetrabz_polstat_1234(e(2),e(1),e(3),e(4),ln(2),ln(1),ln(3),ln(4))
     w1(indx(1)) = libtetrabz_polstat_1234(e(1),e(2),e(3),e(4),ln(1),ln(2),ln(3),ln(4))
     !
     IF(ANY(w1(1:4) < 0d0)) THEN
        WRITE(*,'(100e15.5)') e(1:4)
        WRITE(*,'(100e15.5)') w1(indx(1:4))
        STOP "weighting"
     END IF
     !
  END IF
  !
END SUBROUTINE libtetrabz_polstat3
!
! Results of Integration (1-x-y-z)/(g0+(g1-g0)x+(g2-g0)y+(g3-g0))
!  for 0<x<1, 0<y<1-x, 0<z<1-x-y
!
! 1, Different each other
!
FUNCTION libtetrabz_polstat_1234(g1,g2,g3,g4,lng1,lng2,lng3,lng4) RESULT(w)
  !
  IMPLICIT NONE
  !
  REAL(8),INTENT(IN) :: g1,g2,g3,g4,lng1,lng2,lng3,lng4
  REAL(8) :: w
  !
  REAL(8) :: w2, w3, w4
  !
  w2 = ((lng2 - lng1)/(g2 - g1)*g2 - 1d0)*g2/(g2 - g1)
  w3 = ((lng3 - lng1)/(g3 - g1)*g3 - 1d0)*g3/(g3 - g1)
  w4 = ((lng4 - lng1)/(g4 - g1)*g4 - 1d0)*g4/(g4 - g1)
  w2 = ((w2 - w3)*g2)/(g2 - g3)
  w4 = ((w4 - w3)*g4)/(g4 - g3)
  w = (w4 - w2)/(g4 - g2)
  !
END FUNCTION libtetrabz_polstat_1234
!
! 2, g4 = g1
!
FUNCTION libtetrabz_polstat_1231(g1,g2,g3,lng1,lng2,lng3) RESULT(w)
  !
  IMPLICIT NONE
  !
  REAL(8),INTENT(IN) :: g1,g2,g3,lng1,lng2,lng3
  REAL(8) :: w
  !
  REAL(8) :: w2, w3
  !
  w2 = ((lng2 - lng1)/(g2 - g1)*g2 - 1d0)*g2**2/(g2 - g1) - g1/( &
  &   2d0)
  w2 = w2/(g2 - g1)
  w3 = ((lng3 - lng1)/(g3 - g1)*g3 - 1d0)*g3**2/(g3 - g1) - g1/( &
  &   2d0)
  w3 = w3/(g3 - g1)
  w = (w3 - w2)/(g3 - g2)
  !
END FUNCTION libtetrabz_polstat_1231
!
! 3, g4 = g3
!
FUNCTION libtetrabz_polstat_1233(g1,g2,g3,lng1,lng2,lng3) RESULT(w)
  !
  IMPLICIT NONE
  !
  REAL(8),INTENT(IN) :: g1,g2,g3,lng1,lng2,lng3
  REAL(8) :: w
  !
  REAL(8) :: w2, w3
  !
  w2 = (lng2 - lng1)/(g2 - g1)*g2 - 1d0
  w2 = (g2*w2)/(g2 - g1)
  w3 = (lng3 - lng1)/(g3 - g1)*g3 - 1d0
  w3 = (g3*w3)/(g3 - g1)
  w2 = (w3 - w2)/(g3 - g2)
  w3 = (lng3 - lng1)/(g3 - g1)*g3 - 1d0
  w3 = 1d0 - (2d0*w3*g1)/(g3 - g1)
  w3 = w3/(g3 - g1)
  w = (g3*w3 - g2*w2)/(g3 - g2)
  !
END FUNCTION libtetrabz_polstat_1233
!
! 4, g4 = g1 and g3 = g2
!
FUNCTION libtetrabz_polstat_1221(g1,g2,lng1,lng2) RESULT(w)
  !
  IMPLICIT NONE
  !
  REAL(8),INTENT(IN) :: g1, g2, lng1, lng2
  REAL(8) :: w
  !
  w = 1d0 - (lng2 - lng1)/(g2 - g1)*g1
  w = -1d0 + (2d0*g2*w)/(g2 - g1)
  w = -1d0 + (3d0*g2*w)/(g2 - g1)
  w = w/(2d0*(g2 - g1))
  !
END FUNCTION libtetrabz_polstat_1221
!
! 5, g4 = g3 = g2
!
FUNCTION libtetrabz_polstat_1222(g1,g2,lng1,lng2) RESULT(w)
  !
  IMPLICIT NONE
  !
  REAL(8),INTENT(IN) :: g1, g2, lng1, lng2
  REAL(8) :: w
  !
  w = (lng2 - lng1)/(g2 - g1)*g2 - 1d0
  w = (2d0*g1*w)/(g2 - g1) - 1d0
  w = (3d0*g1*w)/(g2 - g1) + 1d0
  w = w/(2d0*(g2 - g1))
  !
END FUNCTION libtetrabz_polstat_1222
!
! 6, g4 = g3 = g1
!
FUNCTION libtetrabz_polstat_1211(g1,g2,lng1,lng2) RESULT(w)
  !
  IMPLICIT NONE
  !
  REAL(8),INTENT(IN) :: g1,g2,lng1,lng2
  REAL(8) :: w
  !
  w = -1d0 + (lng2 - lng1)/(g2 - g1)*g2
  w = -1d0 + (2d0*g2*w)/(g2 - g1)
  w = -1d0 + (3d0*g2*w)/(2d0*(g2 - g1))
  w = w/(3d0*(g2 - g1))
  !
END FUNCTION libtetrabz_polstat_1211
!
END MODULE libtetrabz_polstat_mod
