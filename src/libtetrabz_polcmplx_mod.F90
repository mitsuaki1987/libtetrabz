MODULE libtetrabz_polcmplx_mod
  !
  IMPLICIT NONE
  !
CONTAINS
!
! Compute Polarization of imaginary frequency
!
SUBROUTINE libtetrabz_polcmplx(ltetra,comm0,bvec,nb,nge,eig1,eig2,ngw,wght0,ne,e0) BIND(C)
  !
#if defined(__MPI)
  USE mpi, ONLY : MPI_IN_PLACE, MPI_DOUBLE_PRECISION, MPI_SUM
#endif
  USE ISO_C_BINDING
  USE libtetrabz_val,   ONLY : nk_local, ik_global, ik_local, kvec, lmpi, linterpol, comm
  USE libtetrabz_common, ONLY : libtetrabz_initialize, libtetrabz_interpol_indx
  IMPLICIT NONE
  !
  INTEGER(C_INT),INTENT(IN) :: ltetra, nge(3), ngw(3), nb, ne
  REAL(C_DOUBLE),INTENT(IN) :: bvec(3,3)
  COMPLEX(C_DOUBLE_COMPLEX),INTENT(IN) :: e0(ne)
  REAL(C_DOUBLE),INTENT(IN) :: eig1(nb,PRODUCT(nge(1:3))), eig2(nb,PRODUCT(nge(1:3)))
  COMPLEX(C_DOUBLE_COMPLEX),INTENT(OUT) :: wght0(ne,nb,nb,PRODUCT(ngw(1:3)))
  INTEGER(C_INT),INTENT(IN),OPTIONAL :: comm0
  !
  INTEGER :: ii, ik, kintp(4)
  REAL(8) :: wintp(4)
  COMPLEX(8),ALLOCATABLE :: wght1(:,:,:,:)
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
  CALL libtetrabz_initialize(ltetra,bvec,nge,ngw,nb,ne)
  !
  IF(lmpi .OR. linterpol) THEN
     !
     ALLOCATE(wght1(ne,nb,nb,nk_local))
     CALL libtetrabz_polcmplx_main(eig1,eig2,e0,wght1)
     !
     ! Interpolation
     !
     wght0(1:ne,1:nb,1:nb,1:PRODUCT(ngw(1:3))) = 0d0
     DO ik = 1, nk_local
        CALL libtetrabz_interpol_indx(ngw,kvec(1:3,ik),kintp,wintp)
        DO ii = 1, 4
           wght0(1:ne,1:nb,1:nb,kintp(ii)) = wght0(1:ne,1:nb,1:nb,       kintp(ii)) &
           &                               + wght1(1:ne,1:nb,1:nb, ik) * wintp(ii)
        END DO
     END DO ! ik = 1, nk_local
     DEALLOCATE(wght1, kvec)
     !
#if defined(__MPI)
     IF(lmpi) &
     &  CALL MPI_allREDUCE(MPI_IN_PLACE, wght0, ne * nb * nb * PRODUCT(ngw(1:3)), &
     &                     MPI_DOUBLE_COMPLEX, MPI_SUM, comm, ierr)
#endif
     !
  ELSE
     CALL libtetrabz_polcmplx_main(eig1,eig2,e0,wght0)
  END IF
  !
END SUBROUTINE libtetrabz_polcmplx
!
! Main SUBROUTINE for Polaization (Imaginaly axis) : Theta(- E1) * Theta(E2) / (E2 - E1 - iw)
!
SUBROUTINE libtetrabz_polcmplx_main(eig1,eig2,e0,polcmplx)
  !
  USE libtetrabz_val, ONLY : nb, nk_local, nt_local, wlsm, ik_global, ik_local, nkBZ, ne
  IMPLICIT NONE
  !
  COMPLEX(8),INTENT(IN) :: e0(ne)
  REAL(8),INTENT(IN) :: eig1(nb,nkBZ), eig2(nb,nkBZ)
  COMPLEX(8),INTENT(OUT) :: polcmplx(ne*nb,nb,nk_local)
  !
  INTEGER :: it, ib, indx(4)
  REAL(8) :: V, thr = 1d-8, e(4), tsmall(4,4), &
  &          ei1(4,nb), ei2(4), ej1(4,nb), ej2(4,nb)
  COMPLEX(8) :: w1(ne*nb,4), w2(ne*nb,4)
  !
  polcmplx(1:ne*nb,1:nb,1:nk_local) = 0d0
  !
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP & SHARED(nt_local,nb,ne,ik_global,ik_local,wlsm,eig1,eig2,polcmplx,thr,e0) &
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
              CALL libtetrabz_polcmplx2(e0,ei2,ej2,w2)
              w1(1:ne*nb,indx(1:4)) = w1(1:ne*nb,            indx(1:4)) &
              &          + V * MATMUL(w2(1:ne*nb,1:4), tsmall(1:4,1:4))
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
              CALL libtetrabz_polcmplx2(e0,ei2,ej2,w2)
              w1(1:ne*nb,indx(1:4)) = w1(1:ne*nb,            indx(1:4)) &
              &          + V * MATMUL(w2(1:ne*nb,1:4), tsmall(1:4,1:4))
              !
           END IF
           !
           CALL libtetrabz_tsmall_b2(e,V,tsmall)
           !
           IF(V > thr) THEN
              !
              ei2(1:4     ) = MATMUL(tsmall(1:4,1:4), ei1(indx(1:4),  ib))
              ej2(1:4,1:nb) = MATMUL(tsmall(1:4,1:4), ej1(indx(1:4),1:nb))
              CALL libtetrabz_polcmplx2(e0,ei2,ej2,w2)
              w1(1:ne*nb,indx(1:4)) = w1(1:ne*nb,            indx(1:4)) &
              &          + V * MATMUL(w2(1:ne*nb,1:4), tsmall(1:4,1:4))
              !
           END IF
           !
           CALL libtetrabz_tsmall_b3(e,V,tsmall)
           !
           IF(V > thr) THEN
              !
              ei2(1:4     ) = MATMUL(tsmall(1:4,1:4), ei1(indx(1:4),  ib))
              ej2(1:4,1:nb) = MATMUL(tsmall(1:4,1:4), ej1(indx(1:4),1:nb))
              CALL libtetrabz_polcmplx2(e0,ei2,ej2,w2)
              w1(1:ne*nb,indx(1:4)) = w1(1:ne*nb,            indx(1:4)) &
              &          + V * MATMUL(w2(1:ne*nb,1:4), tsmall(1:4,1:4))
              !
           END IF
           !
        ELSE IF(e(3) <= 0d0 .AND. 0d0 < e(4)) THEN
           !
           CALL libtetrabz_tsmall_b1(e,V,tsmall)
           !
           IF(V > thr) THEN
              !
              ei2(1:4     ) = MATMUL(tsmall(1:4,1:4), ei1(indx(1:4),  ib))
              ej2(1:4,1:nb) = MATMUL(tsmall(1:4,1:4), ej1(indx(1:4),1:nb))
              CALL libtetrabz_polcmplx2(e0,ei2,ej2,w2)
              w1(1:ne*nb,indx(1:4)) = w1(1:ne*nb,            indx(1:4)) &
              &          + V * MATMUL(w2(1:ne*nb,1:4), tsmall(1:4,1:4))
              !
           END IF
           !
           CALL libtetrabz_tsmall_b2(e,V,tsmall)
           !
           IF(V > thr) THEN
              !
              ei2(1:4     ) = MATMUL(tsmall(1:4,1:4), ei1(indx(1:4),  ib))
              ej2(1:4,1:nb) = MATMUL(tsmall(1:4,1:4), ej1(indx(1:4),1:nb))
              CALL libtetrabz_polcmplx2(e0,ei2,ej2,w2)
              w1(1:ne*nb,indx(1:4)) = w1(1:ne*nb,            indx(1:4)) &
              &          + V * MATMUL(w2(1:ne*nb,1:4), tsmall(1:4,1:4))
              !
           END IF
           !
           CALL libtetrabz_tsmall_b3(e,V,tsmall)
           !
           IF(V > thr) THEN
              !
              ei2(1:4     ) = MATMUL(tsmall(1:4,1:4), ei1(indx(1:4),  ib))
              ej2(1:4,1:nb) = MATMUL(tsmall(1:4,1:4), ej1(indx(1:4),1:nb))
              CALL libtetrabz_polcmplx2(e0,ei2,ej2,w2)
              w1(1:ne*nb,indx(1:4)) = w1(1:ne*nb,            indx(1:4)) &
              &          + V * MATMUL(w2(1:ne*nb,1:4), tsmall(1:4,1:4))
              !
           END IF
           !
        ELSE IF( e(4) <= 0d0 ) THEN
           !
           ei2(1:4     ) = ei1(1:4,  ib)
           ej2(1:4,1:nb) = ej1(1:4,1:nb)
           CALL libtetrabz_polcmplx2(e0,ei2,ej2,w2)
           w1(1:ne*nb,1:4) = w1(1:ne*nb,1:4) + w2(1:ne*nb,1:4)
           !
        ELSE
           !
           CYCLE
           !
        END IF
        !
        polcmplx(1:ne*nb,ib,ik_local(1:20,it)) = polcmplx(1:ne*nb,ib,   ik_local(1:20,it)) &
        &                                     + MATMUL(w1(1:ne*nb,1:4), wlsm(1:4,1:20))
        !
     END DO ! ib = 1, nb
     !$OMP END DO NOWAIT
     !
  END DO ! it
  !
  !$OMP END PARALLEL
  !
  polcmplx(1:ne*nb,1:nb,1:nk_local) = polcmplx(1:ne*nb,1:nb,1:nk_local) / DBLE(6 * nkBZ)
  !
END SUBROUTINE libtetrabz_polcmplx_main
!
! Tetrahedra method for theta( - E2)
!
SUBROUTINE libtetrabz_polcmplx2(e0,ei1,ej1,w1)
  !
  USE libtetrabz_val, ONLY : nb, ne
  IMPLICIT NONE
  !
  COMPLEX(8),INTENT(IN) :: e0(ne)
  REAL(8),INTENT(IN) :: ei1(4), ej1(4,nb)
  COMPLEX(8),INTENT(OUT) :: w1(ne,nb,4)
  !
  INTEGER :: ib, indx(4)
  REAL(8) :: V, thr = 1d-8, e(4), ei2(4), ej2(4), tsmall(4,4)
  COMPLEX(8) :: w2(ne,4)
  !
  DO ib = 1, nb
     !
     w1(1:ne,ib,1:4) = 0d0
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
           CALL libtetrabz_polcmplx3(e0,ei2,ej2,w2)
           w1(1:ne,ib,indx(1:4)) = w1(1:ne,ib,         indx(1:4)) &
           &          + V * MATMUL(w2(1:ne,1:4), tsmall(1:4,1:4))
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
           CALL libtetrabz_polcmplx3(e0,ei2,ej2,w2)
           w1(1:ne,ib,indx(1:4)) = w1(1:ne,ib,         indx(1:4)) &
           &            + V * MATMUL(w2(1:ne,1:4), tsmall(1:4,1:4))
           !
        END IF
        !
        CALL libtetrabz_tsmall_b2(e,V,tsmall)
        !
        IF(V > thr) THEN
           !
           ei2(1:4) = MATMUL(tsmall(1:4,1:4), ei1(indx(1:4)   ))
           ej2(1:4) = MATMUL(tsmall(1:4,1:4), ej1(indx(1:4),ib))
           CALL libtetrabz_polcmplx3(e0,ei2,ej2,w2)
           w1(1:ne,ib,indx(1:4)) = w1(1:ne,ib,         indx(1:4)) &
           &          + V * MATMUL(w2(1:ne,1:4), tsmall(1:4,1:4))
           !
        END IF
        !
        CALL libtetrabz_tsmall_b3(e,V,tsmall)
        !
        IF(V > thr) THEN
           !
           ei2(1:4) = MATMUL(tsmall(1:4,1:4), ei1(indx(1:4)   ))
           ej2(1:4) = MATMUL(tsmall(1:4,1:4), ej1(indx(1:4),ib))
           CALL libtetrabz_polcmplx3(e0,ei2,ej2,w2)
           w1(1:ne,ib,indx(1:4)) = w1(1:ne,ib,         indx(1:4)) &
           &          + V * MATMUL(w2(1:ne,1:4), tsmall(1:4,1:4))
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
           CALL libtetrabz_polcmplx3(e0,ei2,ej2,w2)
           w1(1:ne,ib,indx(1:4)) = w1(1:ne,ib,         indx(1:4)) &
           &          + V * MATMUL(w2(1:ne,1:4), tsmall(1:4,1:4))
           !
        END IF
        !
        CALL libtetrabz_tsmall_c2(e,V,tsmall)
        !
        IF(V > thr) THEN
           !
           ei2(1:4) = MATMUL(tsmall(1:4,1:4), ei1(indx(1:4)   ))
           ej2(1:4) = MATMUL(tsmall(1:4,1:4), ej1(indx(1:4),ib))
           CALL libtetrabz_polcmplx3(e0,ei2,ej2,w2)
           w1(1:ne,ib,indx(1:4)) = w1(1:ne,ib,         indx(1:4)) &
           &          + V * MATMUL(w2(1:ne,1:4), tsmall(1:4,1:4))
           !
        END IF
        !
        CALL libtetrabz_tsmall_c3(e,V,tsmall)
        !
        IF(V > thr) THEN
           !
           ei2(1:4) = MATMUL(tsmall(1:4,1:4), ei1(indx(1:4)   ))
           ej2(1:4) = MATMUL(tsmall(1:4,1:4), ej1(indx(1:4),ib))
           CALL libtetrabz_polcmplx3(e0,ei2,ej2,w2)
           w1(1:ne,ib,indx(1:4)) = w1(1:ne,ib,         indx(1:4)) &
           &          + V * MATMUL(w2(1:ne,1:4), tsmall(1:4,1:4))
           !
        END IF
        !
     ELSE IF(e(4) <= 0d0) THEN
        !
        ei2(1:4) = ei1(1:4)
        ej2(1:4) = ej1(1:4,ib)
        CALL libtetrabz_polcmplx3(e0,ei2,ej2,w2)
        w1(1:ne,ib,1:4) = w1(1:ne,ib,1:4) + w2(1:ne,1:4)
        !
     END IF
     !
  END DO
  !
END SUBROUTINE libtetrabz_polcmplx2
!
! Tetarahedra method for delta(om - ep + e)
!
SUBROUTINE libtetrabz_polcmplx3(e0,ei1,ej1,w1)
  !
  USE libtetrabz_val, ONLY : ne
  IMPLICIT NONE
  !
  COMPLEX(8),INTENT(IN) :: e0(ne)
  REAL(8),INTENT(IN) :: ei1(4), ej1(4)
  COMPLEX(8),INTENT(OUT) :: w1(ne,4)
  !
  INTEGER :: ie, ierr, indx(4)
  REAL(8) :: w2(2,4), de(4), e(4), thr
  !
  de(1:4) = ej1(1:4) - ei1(1:4)
  CALL libtetrabz_sort(4,de,indx)
  !
  DO ie = 1, ne
     !
     e(1:4) = (de(1:4) - DBLE(e0(ie))) / AIMAG(e0(ie))
     !thr = maxval(de(1:4)) * 1d-3
     thr = max(1d-3,  maxval(e(1:4)) * 1d-2)
     !
     IF(ABS(e(4) - e(3)) < thr ) THEN
        IF(ABS(e(4) - e(2)) < thr ) THEN
           IF(ABS(e(4) - e(1)) < thr ) THEN
              !
              ! e(4) = e(3) = e(2) = e(1)
              !
              w2(1,4) = 0.25d0 * e(4) / ((1d0 + e(4)**2))
              w2(2,4) = 0.25d0        / ((1d0 + e(4)**2))
              w2(1:2,3) = w2(1:2,4)
              w2(1:2,2) = w2(1:2,4)
              w2(1:2,1) = w2(1:2,4)
              !
           ELSE
              !
              ! e(4) = e(3) = e(2)
              !
              w2(1:2,4) = libtetrabz_polcmplx_1211(e(4),e(1))
              w2(1:2,3) = w2(1:2,4)
              w2(1:2,2) = w2(1:2,4)
              w2(1:2,1) = libtetrabz_polcmplx_1222(e(1),e(4))
              !
              IF(ANY(w2(1:2,1:4) < 0d0)) THEN
                 WRITE(*,*) ie
                 WRITE(*,'(100e15.5)') e(1:4)
                 WRITE(*,'(2e15.5)') w2(1:2,1:4)
                 STOP "weighting 4=3=2"
              END IF
              !
           END IF
        ELSE IF(ABS(e(2) - e(1)) < thr ) THEN
           !
           ! e(4) = e(3), e(2) = e(1)
           !
           w2(1:2,4) = libtetrabz_polcmplx_1221(e(4),e(2))
           w2(1:2,3) = w2(1:2,4)
           w2(1:2,2) = libtetrabz_polcmplx_1221(e(2),e(4))
           w2(1:2,1) = w2(1:2,2)
           !
           IF(ANY(w2(1:2,1:4) < 0d0)) THEN
              WRITE(*,*) ie
              WRITE(*,'(100e15.5)') e(1:4)
              WRITE(*,'(2e15.5)') w2(1:2,1:4)
              STOP "weighting 4=3 2=1"
           END IF
           !
        ELSE
           !
           ! e(4) = e(3)
           !
           w2(1:2,4) = libtetrabz_polcmplx_1231(e(4),e(1),e(2))
           w2(1:2,3) = w2(1:2,4)
           w2(1:2,2) = libtetrabz_polcmplx_1233(e(2),e(1),e(4))
           w2(1:2,1) = libtetrabz_polcmplx_1233(e(1),e(2),e(4))
           !
           IF(ANY(w2(1:2,1:4) < 0d0)) THEN
              WRITE(*,*) ie
              WRITE(*,'(100e15.5)') e(1:4)
              WRITE(*,'(2e15.5)') w2(1:2,1:4)
              STOP "weighting 4=3"
           END IF
           !
        END IF
     ELSE IF(ABS(e(3) - e(2)) < thr) THEN
        IF(ABS(e(3) - e(1)) < thr) THEN
           !
           ! e(3) = e(2) = e(1)
           !
           w2(1:2,4) = libtetrabz_polcmplx_1222(e(4),e(3))
           w2(1:2,3) = libtetrabz_polcmplx_1211(e(3),e(4))
           w2(1:2,2) = w2(1:2,3)
           w2(1:2,1) = w2(1:2,3)
           !
           IF(ANY(w2(1:2,1:4) < 0d0)) THEN
              WRITE(*,*) ie
              WRITE(*,'(100e15.5)') e(1:4)
              WRITE(*,'(2e15.5)') w2(1:2,1:4)
              STOP "weighting 3=2=1"
           END IF
           !
        ELSE
           !
           ! e(3) = e(2)
           !
           w2(1:2,4) = libtetrabz_polcmplx_1233(e(4),e(1),e(3))
           w2(1:2,3) = libtetrabz_polcmplx_1231(e(3),e(1),e(4))
           w2(1:2,2) = w2(1:2,3)
           w2(1:2,1) = libtetrabz_polcmplx_1233(e(1),e(4),e(3))
           !
           IF(ANY(w2(1:2,1:4) < 0d0)) THEN
              WRITE(*,*) ie
              WRITE(*,'(100e15.5)') e(1:4)
              WRITE(*,'(2e15.5)') w2(1:2,1:4)
              STOP "weighting 3=2"
           END IF
           !
        END IF
     ELSE IF(ABS(e(2) - e(1)) < thr) THEN
        !
        ! e(2) = e(1)
        !
        w2(1:2,4) = libtetrabz_polcmplx_1233(e(4),e(3),e(2))
        w2(1:2,3) = libtetrabz_polcmplx_1233(e(3),e(4),e(2))
        w2(1:2,2) = libtetrabz_polcmplx_1231(e(2),e(3),e(4))
        w2(1:2,1) = w2(1:2,2)
        !
        IF(ANY(w2(1:2,1:4) < 0d0)) THEN
           WRITE(*,*) ie
           WRITE(*,'(100e15.5)') e(1:4)
           WRITE(*,'(2e15.5)') w2(1:2,1:4)
           STOP "weighting 2=1"
        END IF
        !
     ELSE
        !
        ! Different each other.
        !
        w2(1:2,4) = libtetrabz_polcmplx_1234(e(4),e(1),e(2),e(3))
        w2(1:2,3) = libtetrabz_polcmplx_1234(e(3),e(1),e(2),e(4))
        w2(1:2,2) = libtetrabz_polcmplx_1234(e(2),e(1),e(3),e(4))
        w2(1:2,1) = libtetrabz_polcmplx_1234(e(1),e(2),e(3),e(4))
        !
        IF(ANY(w2(1:2,1:4) < 0d0)) THEN
           WRITE(*,*) ie
           WRITE(*,'(100e15.5)') e(1:4)
           WRITE(*,'(2e15.5)') w2(1:2,1:4)
           STOP "weighting"
        END IF
        !
     END IF
     !
     w1(ie,indx(1:4)) = CMPLX(w2(1,1:4) / AIMAG(e0(ie)), w2(2,1:4) / (- AIMAG(e0(ie))))
     !
  END DO ! ie
  !
END SUBROUTINE libtetrabz_polcmplx3
!
! Results of Integration (1-x-y-z)/(g0+(g1-g0)x+(g2-g0)y+(g3-g0))
!  for 0<x<1, 0<y<1-x, 0<z<1-x-y
!
! 1, Different each other
!
FUNCTION libtetrabz_polcmplx_1234(g1,g2,g3,g4) RESULT(w)
  !
  IMPLICIT NONE
  !
  REAL(8),INTENT(IN) :: g1, g2, g3, g4
  REAL(8) :: w(2)
  !
  REAL(8) :: w2, w3, w4
  !
  ! Real
  !
  w2 = 2d0*(3d0*g2**2 - 1d0)*(ATAN(g2) - ATAN(g1)) + (g2**2 - &
  &      3d0)*g2*LOG((1d0 + g2**2)/( 1d0 + g1**2))
  w2 = -2d0*(g2**2 - 1d0) + w2/(g2 - g1 )
  w2 = w2/(g2 - g1 )
  w3 = 2d0*(3d0*g3**2 - 1d0)*(ATAN(g3) - ATAN(g1)) + (g3**2 -  &
  &      3d0)*g3*LOG((1d0 + g3**2)/( 1d0 + g1**2))
  w3 = -2d0*(g3**2 - 1d0) + w3/(g3 - g1 )
  w3 = w3/(g3 - g1 )
  w4 = 2d0*(3d0*g4**2 - 1d0)*(ATAN(g4) - ATAN(g1)) + (g4**2 -  &
  &      3d0)*g4*LOG((1d0 + g4**2)/( 1d0 + g1**2))
  w4 = -2d0*(g4**2 - 1d0) + w4/(g4 - g1 )
  w4 = w4/(g4 - g1 )
  w2 = (w2 - w3)/(g2 - g3)
  w4 = (w4 - w3)/(g4 - g3)
  w(1) = (w4 - w2)/(2d0*(g4 - g2))
  !
  ! Imaginal
  !
  w2 = 2d0*(3d0 - g2**2)* &
  &    g2*(ATAN(g2) - ATAN(g1)) + (3d0*g2**2 - 1d0)* &
  &    LOG((1d0 + g2**2)/(1d0 + g1**2))
  w2 = 4d0*g2 - w2/(g2 - g1)
  w2 = w2/(g2 - g1)
  w3 = 2d0*(3d0 - g3**2)* &
  &    g3*(ATAN(g3) - ATAN(g1)) + (3d0*g3**2 - 1d0)* &
  &    LOG((1d0 + g3**2)/(1d0 + g1**2))
  w3 = 4d0*g3 - w3/(g3 - g1)
  w3 = w3/(g3 - g1)
  w4 = 2d0*(3d0 - g4**2)* &
  &    g4*(ATAN(g4) - ATAN(g1)) + (3d0*g4**2 - 1d0)* &
  &    LOG((1d0 + g4**2)/(1d0 + g1**2))
  w4 = 4d0*g4 - w4/(g4 - g1)
  w4 = w4/(g4 - g1)
  w2 = (w2 - w3)/(g2 - g3)
  w4 = (w4 - w3)/(g4 - g3)
  w(2) = (w4 - w2)/(2d0*(g4 - g2))
  !
END FUNCTION libtetrabz_polcmplx_1234
!
! 2, g4 = g1
!
FUNCTION libtetrabz_polcmplx_1231(g1,g2,g3) RESULT(w)
  !
  IMPLICIT NONE
  !
  REAL(8),INTENT(IN) :: g1, g2, g3
  REAL(8) :: w(2)
  !
  REAL(8) :: w2, w3
  !
  ! Real
  !
  w2 = 2d0*(-1d0 + 3d0*g2**2)*(ATAN(g2) - ATAN(g1)) +  &
  &   g2*(-3d0 + g2**2)*LOG((1d0 + g2**2)/(1d0 + g1**2))
  w2 = 2d0*(1d0 - g2**2) + w2/(g2 - g1)
  w2 = -g1 + w2/(g2 - g1)
  w2 = w2/(g2 - g1)
  w3 = 2d0*(-1d0 + 3d0*g3**2)*(ATAN(g3) - ATAN(g1)) +  &
  &   g3*(-3d0 + g3**2)*LOG((1d0 + g3**2)/(1d0 + g1**2))
  w3 = 2d0*(1 - g3**2) + w3/(g3 - g1)
  w3 = -g1 + w3/(g3 - g1)
  w3 = w3/(g3 - g1)
  w(1) = (w3 - w2)/(2d0*(g3 - g2))
  !
  ! Imaginal
  !
  w2 = 2d0* &
  &    g2*(3d0 - g2**2)*(ATAN(g2) - ATAN(g1)) + (-1d0 + 3d0*g2**2)* &
  &    LOG((1d0 + g2**2)/(1d0 + g1**2))
  w2 = 4d0*g2 - w2/(g2 - g1)
  w2 = 1 + w2/(g2 - g1)
  w2 = w2/(g2 - g1)
  w3 = 2d0* &
  &    g3*(3d0 - g3**2)*(ATAN(g3) - ATAN(g1)) + (-1d0 + 3d0*g3**2)* &
  &    LOG((1d0 + g3**2)/(1d0 + g1**2))
  w3 = 4d0*g3 - w3/(g3 - g1)
  w3 = 1 + w3/(g3 - g1)
  w3 = w3/(g3 - g1)
  w(2) = (w3 - w2)/(2d0*(g3 - g2))
  !
END FUNCTION libtetrabz_polcmplx_1231
!
! 3, g4 = g3
!
FUNCTION libtetrabz_polcmplx_1233(g1, g2, g3) RESULT(w)
  !
  IMPLICIT NONE
  !
  REAL(8),INTENT(IN) :: g1, g2, g3
  REAL(8) :: w(2)
  !
  REAL(8) :: w2, w3
  !
  ! Real
  !
  w2 = 2d0*(1d0 - 3d0*g2**2)*(ATAN(g2) - ATAN(g1)) +  &
  &   g2*(3d0 - g2**2)*LOG((1d0 + g2**2)/(1d0 + g1**2))
  w2 = 2d0*(1 - g2**2) - w2/(g2 - g1)
  w2 = w2/(g2 - g1)
  w3 = 2d0*(1d0 - 3d0*g3**2)*(ATAN(g3) - ATAN(g1)) +  &
  &   g3*(3d0 - g3**2)*LOG((1d0 + g3**2)/(1d0 + g1**2))
  w3 = 2d0*(1 - g3**2) - w3/(g3 - g1)
  w3 = w3/(g3 - g1)
  w2 = (w3 - w2)/(g3 - g2)
  w3 = 4d0*(1d0 - 3d0*g1*g3)*(ATAN(g3) - ATAN(g1)) + (3d0*g1 +  &
  &      3d0*g3 - 3d0*g1*g3**2 + g3**3) * LOG((1d0 + g3**2)/( &
  &     1d0 + g1**2))
  w3 = -4d0*(1d0 - g1**2) + w3/(g3 - g1)
  w3 = 4d0*g1 + w3/(g3 - g1)
  w3 = w3/(g3 - g1)
  w(1) = (w3 - w2)/(2d0*(g3 - g2))
  !
  ! Imaginal
  !
  w2 = 2d0* &
  &    g2*(3d0 - g2**2)*(ATAN(g2) - ATAN(g1)) + (-1d0 + 3d0*g2**2)* &
  &    LOG((1d0 + g2**2)/(1d0 + g1**2))
  w2 = 4d0*g2 - w2/(g2 - g1)
  w2 = w2/(g2 - g1)
  w3 = 2d0* &
  &    g3*(3d0 - g3**2)*(ATAN(g3) - ATAN(g1)) + (-1d0 + 3d0*g3**2)* &
  &    LOG((1d0 + g3**2)/(1d0 + g1**2))
  w3 = 4d0*g3 - w3/(g3 - g1)
  w3 = w3/(g3 - g1)
  w2 = (w3 - w2)/(g3 - g2)
  w3 = (3d0*g1 - 3d0*g1*g3**2 + 3d0*g3 + g3**3)*(ATAN(g3) -  &
  &      ATAN(g1)) + (3d0*g1*g3 - 1d0)* &
  &    LOG((1d0 + g3**2)/(1d0 + g1**2))
  w3 = w3/(g3 - g1) - 4d0*g1
  w3 = w3/(g3 - g1) - 2d0
  w3 = (2d0*w3)/(g3 - g1)
  w(2) = (w3 - w2)/(2d0*(g3 - g2))
  !
END FUNCTION libtetrabz_polcmplx_1233
!
! 4, g4 = g1 and g3 = g2
!
FUNCTION libtetrabz_polcmplx_1221(g1,g2) RESULT(w)
  !
  IMPLICIT NONE
  !
  REAL(8),INTENT(IN) :: g1, g2
  REAL(8) :: w(2)
  !
  ! Real
  !
  w(1) = -2d0*(-1d0 + 2d0*g1*g2 + g2**2)*(ATAN(g2) -  &
  &      ATAN(g1)) + (g1 + 2d0*g2 - g1*g2**2)* &
  &    LOG((1d0 + g2**2)/(1d0 + g1**2))
  w(1) = 2d0*(-1d0 + g1**2) + w(1)/(g2 - g1)
  w(1) = 3d0*g1 + w(1)/(g2 - g1)
  w(1) = 2d0 + (3d0*w(1))/(g2 - g1)
  w(1) = w(1)/(2d0*(g2 - g1))
  !
  ! Imaginal
  !
  w(2) = 2d0*(g1 + 2d0*g2 - g1*g2**2)*(ATAN(g2) -  &
  &      ATAN(g1)) + (-1d0 + 2d0*g1*g2 + g2**2)* &
  &    LOG((1 + g2**2)/(1 + g1**2))
  w(2) = -4d0*g1 + w(2)/(g2 - g1)
  w(2) = -3d0 + w(2)/(g2 - g1)
  w(2) = (3d0*w(2))/(2d0*(g2 - g1)**2)
  !
END FUNCTION libtetrabz_polcmplx_1221
!
! 5, g4 = g3 = g2
!
FUNCTION libtetrabz_polcmplx_1222(g1,g2) RESULT(w)
  !
  IMPLICIT NONE
  !
  REAL(8),INTENT(IN) :: g1, g2
  REAL(8) :: w(2)
  !
  ! Real
  !
  w(1) = 2d0*(-1d0 + g1**2 + 2d0*g1*g2)*(ATAN(g2) -  &
  &      ATAN(g1)) + (-2d0*g1 - g2 + g1**2*g2) * LOG((1d0 + g2**2)/( &
  &     1d0 + g1**2))
  w(1) = 2d0*(1d0 - g1**2) + w(1)/(g2 - g1)
  w(1) = g1 - w(1)/(g2 - g1)
  w(1) = 1d0 - (3d0*w(1))/(g2 - g1)
  w(1) = w(1)/(2d0*(g2 - g1))
  !
  ! Imaginal
  !
  w(2) = 2d0*(-2d0*g1 - g2 + g1**2*g2)*(ATAN(g2) - ATAN(g1)) + (1d0 - &
  &       g1**2 - 2d0*g1*g2) * LOG((1d0 + g2**2)/(1d0 + g1**2))
  w(2) = 4d0*g1 + w(2)/(g2 - g1)
  w(2) = 1d0 + w(2)/(g2 - g1)
  w(2) = (3d0*w(2))/(2d0*(g2 - g1)**2)
  !
END FUNCTION libtetrabz_polcmplx_1222
!
! 6, g4 = g3 = g1
!
FUNCTION libtetrabz_polcmplx_1211(g1,g2) RESULT(w)
  !
  IMPLICIT NONE
  !
  REAL(8),INTENT(IN) :: g1, g2
  REAL(8) :: w(2)
  !
  ! Real
  !
  w(1) = 2d0*(3d0*g2**2 - 1d0)*(ATAN(g2) - ATAN(g1)) +  &
  &   g2*(g2**2 - 3d0)*LOG((1d0 + g2**2)/(1d0 + g1**2))
  w(1) = 2d0*(1d0 - g1**2) + w(1)/(g2 - g1)
  w(1) = -5d0*g1 + w(1)/(g2 - g1)
  w(1) = -11d0 + (3d0*w(1))/(g2 - g1)
  w(1) = w(1)/(6d0*(g2 - g1))
  !
  ! Imaginal
  !
  w(2) = 2d0*g2*(-3d0 + g2**2)*(ATAN(g2) - ATAN(g1)) + (1d0 -  &
  &      3d0*g2**2)*LOG((1d0 + g2**2)/(1d0 + g1**2))
  w(2) = 4d0*g2 + w(2)/(g2 - g1)
  w(2) = 1d0 + w(2)/(g2 - g1)
  w(2) = w(2)/(2d0*(g2 - g1)**2)
  !
END FUNCTION libtetrabz_polcmplx_1211
!
END MODULE libtetrabz_polcmplx_mod
