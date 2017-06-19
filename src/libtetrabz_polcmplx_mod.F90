MODULE libtetrabz_polcmplx_mod
  !
  IMPLICIT NONE
  !
CONTAINS
!
! Compute Polarization of imaginary frequency
!
SUBROUTINE libtetrabz_mpi_polcmplx(ltetra0,comm0,bvec,nb0,nge,eig1,eig2,ngw,wght0,ne0,e0)
  !
  USE mpi, ONLY : MPI_IN_PLACE, MPI_DOUBLE_PRECISION, MPI_SUM
  USE libtetrabz_vals, ONLY : ltetra, ng, nb, nk0, indx1, indx2, indx3, ne
  USE libtetrabz_common, ONLY : libtetrabz_initialize, libtetrabz_interpol_weight
  USE libtetrabz_polcmplx_mod, ONLY : libtetrabz_polcmplx_main
  !
  USE libtetrabz_mpi_routines, ONLY : comm, libtetrabz_mpi_kgrid
  !
  INTEGER,INTENT(IN) :: ltetra0, comm0, nge(3), ngw(3), nb0, ne0
  REAL(8),INTENT(IN) :: bvec(3,3), e0(ne0)
  REAL(8),INTENT(IN) :: eig1(nb0,PRODUCT(nge(1:3))), eig2(nb0,PRODUCT(nge(1:3)))
  REAL(8),INTENT(OUT) :: wght0(2,ne0,nb0,nb0,PRODUCT(ngw(1:3)))
  !
  INTEGER :: ierr, nn
  REAL(8),ALLOCATABLE :: wght1(:,:)
  !
  ltetra = ltetra0
  comm = comm0
  nb = nb0
  ng(1:3) = nge(1:3)
  ne = ne0
  nn = 2 * ne * nb * nb
  !
  CALL libtetrabz_initialize(bvec)
  CALL libtetrabz_mpi_kgrid()
  !
  ALLOCATE(wght1(nn, nk0))
  CALL libtetrabz_polcmplx_main(eig1,eig2,e0,wght1)
  !
  CALL libtetrabz_interpol_weight(nn,ngw,nge,wght0,wght1)
  !
  DEALLOCATE(wght1, indx1, indx2, indx3)
  !
  CALL MPI_allREDUCE(MPI_IN_PLACE, wght0, nn * PRODUCT(ngw(1:3)), &
  &                  MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)
  !
END SUBROUTINE libtetrabz_mpi_polcmplx
!
! Main SUBROUTINE for Polaization (Imaginaly axis) : Theta(- E1) * Theta(E2) / (E2 - E1 - iw)
!
SUBROUTINE libtetrabz_polcmplx_main(eig1,eig2,e0,poli)
  !
  USE libterabz_val, ONLY : nb, nk, nk0, ne, fst, lst, indx1, indx2, wlsm
  IMPLICIT NONE
  !
  REAL(8),INTENT(IN) :: eig1(nb,nk), eig2(nb,nk), e0(ne)
  COMPLEX(8),INTENT(OUT) :: polcmplx(ne*nb,nb,nk0)
  !
  INTEGER :: it, ib
  REAL(8) :: V, thr = 1d-8, &
  &          ei1(4,nb), ei2(4), ej1(4,nb), ej2(4,nb)
  COMPLEX(8) :: w1(ne*nb,4), w2(ne*nb,4)
  !
  polcmplx(1:ne*nb,1:nb,1:nk0) = 0d0
  !
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP & SHARED(fst,lst,nb,ne,nn,indx1,indx2,wlsm,eig1,eig2,w0,polcmplx,thr,e0) &
  !$OMP & PRIVATE(ib,it,ii,e,a,tmp,tmp2,w1,w2,ei,ei2,ej,ej2,V)
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
        IF(e(1) <= 0d0 .AND. 0d0 < e(2)) THEN
           !
           CALL libtetrabz_tsmall_a1(e,V,tsmall)
           !
           IF(V > thr) THEN
              !
              ei2(1:4     ) = MATMUL(tsmall(1:4,1:4), ei1(sind(1:4),  ib))
              ej2(1:4,1:nb) = MATMUL(tsmall(1:4,1:4), ej1(sind(1:4),1:nb))
              CALL libtetrabz_polcmplx2(e0,ei2,ej2,w2)
              w1(1:ne*nb,sind(1:4)) = w1(1:ne*nb,            sind(1:4)) &
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
              ei2(1:4     ) = MATMUL(tsmall(1:4,1:4), ei1(sind(1:4),  ib))
              ej2(1:4,1:nb) = MATMUL(tsmall(1:4,1:4), ej1(sind(1:4),1:nb))
              CALL libtetrabz_polcmplx2(e0,ei2,ej2,w2)
              w1(1:ne*nb,sind(1:4)) = w1(1:ne*nb,            sind(1:4)) &
              &          + V * MATMUL(w2(1:ne*nb,1:4), tsmall(1:4,1:4))
              !
           END IF
           !
           CALL libtetrabz_tsmall_b2(e,V,tsmall)
           !
           IF(V > thr) THEN
              !
              ei2(1:4     ) = MATMUL(tsmall(1:4,1:4), ei1(sind(1:4),  ib))
              ej2(1:4,1:nb) = MATMUL(tsmall(1:4,1:4), ej1(sind(1:4),1:nb))
              CALL libtetrabz_polcmplx2(e0,ei2,ej2,w2)
              w1(1:ne*nb,sind(1:4)) = w1(1:ne*nb,            sind(1:4)) &
              &          + V * MATMUL(w2(1:ne*nb,1:4), tsmall(1:4,1:4))
              !
           END IF
           !
           CALL libtetrabz_tsmall_b3(e,V,tsmall)
           !
           IF(V > thr) THEN
              !
              ei2(1:4     ) = MATMUL(tsmall(1:4,1:4), ei1(sind(1:4),  ib))
              ej2(1:4,1:nb) = MATMUL(tsmall(1:4,1:4), ej1(sind(1:4),1:nb))
              CALL libtetrabz_polcmplx2(e0,ei2,ej2,w2)
              w1(1:ne*nb,sind(1:4)) = w1(1:ne*nb,            sind(1:4)) &
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
              ei2(1:4     ) = MATMUL(tsmall(1:4,1:4), ei1(sind(1:4),  ib))
              ej2(1:4,1:nb) = MATMUL(tsmall(1:4,1:4), ej1(sind(1:4),1:nb))
              CALL libtetrabz_polcmplx2(e0,ei2,ej2,w2)
              w1(1:ne*nb,sind(1:4)) = w1(1:ne*nb,            sind(1:4)) &
              &          + V * MATMUL(w2(1:ne*nb,1:4), tsmall(1:4,1:4))
              !
           END IF
           !
           CALL libtetrabz_tsmall_b2(e,V,tsmall)
           !
           IF(V > thr) THEN
              !
              ei2(1:4     ) = MATMUL(tsmall(1:4,1:4), ei1(sind(1:4),  ib))
              ej2(1:4,1:nb) = MATMUL(tsmall(1:4,1:4), ej1(sind(1:4),1:nb))
              CALL libtetrabz_polcmplx2(e0,ei2,ej2,w2)
              w1(1:ne*nb,sind(1:4)) = w1(1:ne*nb,            sind(1:4)) &
              &          + V * MATMUL(w2(1:ne*nb,1:4), tsmall(1:4,1:4))
              !
           END IF
           !
           CALL libtetrabz_tsmall_b3(e,V,tsmall)
           !
           IF(V > thr) THEN
              !
              ei2(1:4     ) = MATMUL(tsmall(1:4,1:4), ei1(sind(1:4),  ib))
              ej2(1:4,1:nb) = MATMUL(tsmall(1:4,1:4), ej1(sind(1:4),1:nb))
              CALL libtetrabz_polcmplx2(e0,ei2,ej2,w2)
              w1(1:ne*nb,sind(1:4)) = w1(1:ne*nb,            sind(1:4)) &
              &          + V * MATMUL(w2(1:ne*nb,1:4), tsmall(1:4,1:4))
              !
           END IF
           !
        ELSE IF( e(4) <= 0d0 ) THEN
           !
           ei2(1:4     ) = ei1(1:4,  ib))
           ej2(1:4,1:nb) = ej1(1:4,1:nb))
           CALL libtetrabz_polcmplx2(e0,ei2,ej2,w2)
           w1(1:ne*nb,1:4) = w1(1:ne*nb,1:4) + w2(1:ne*nb,1:4)
           !
        ELSE
           !
           CYCLE
           !
        END IF
        !
        polcmplx(1:ne*nb,ib,indx2(1:20,it)) = polcmplx(1:ne*nb,ib,      indx2(1:20,it)) &
        &                                  + MATMUL(w1(1:ne*nb,1:4), wlsm(1:4,1:20))
        !
     END DO ! ib = 1, nb
     !$OMP END DO NOWAIT
     !
  END DO ! it
  !
  !$OMP END PARALLEL
  !
  polcmplx(1:ne*nb,1:nb,1:nk0) = polcmplx(1:ne*nb,1:nb,1:nk0) / DBLE(6 * nk)
  !
END SUBROUTINE libtetrabz_polcmplx_main
!
! Tetrahedra method for theta( - E2)
!
SUBROUTINE libtetrabz_polcmplx2(e0,ei1,ej1,w1)
  !
  USE libterabz_val, ONLY : nb, ne
  IMPLICIT NONE
  !
  REAL(8),INTENT(IN) :: e0(ne), ei1(4), ej1(4,nb)
  COMPLEX(8),INTENT(OUT) :: w1(ne,nb,4)
  !
  INTEGER :: ib
  REAL(8) :: V, thr = 1d-8, e(4)
  COMPLEX(8) :: w2(ne,4)
  !
  DO ib = 1, nb
     !
     w1(1:ne,ib,1:4) = 0d0
     e(1:4) = - ej1(1:4, ib)
     CALL libtetrabz_sort(4,e,sind)
     !
     IF((e(1) <= 0d0 .AND. 0d0 < e(2)) .OR. (e(1) < 0d0 .AND. 0d0 <= e(2))) THEN
        !
        CALL libtetrabz_tsmall_a1(e,V,tsmall)
        !
        IF(V > thr) THEN
           !
           ei2(1:4) = MATMUL(tsmall(1:4,1:4), ei1(sind(1:4)   ))
           ej2(1:4) = MATMUL(tsmall(1:4,1:4), ej1(sind(1:4),ib))
           CALL libtetrabz_polcmplx3(e0,ei2,ej2,w2)
           w1(1:ne,ib,sind(1:4)) = w1(1:ne,ib,         sind(1:4)) &
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
           ei2(1:4) = MATMUL(tsmall(1:4,1:4), ei1(sind(1:4)   ))
           ej2(1:4) = MATMUL(tsmall(1:4,1:4), ej1(sind(1:4),ib))
           CALL libtetrabz_polcmplx3(e0,ei2,ej2,w2)
           w1(1:ne,ib,sind(1:4)) = w1(1:ne,ib,         sind(1:4)) &
           &            + V * MATMUL(w2(1:ne,1:4), tsmall(1:4,1:4))
           !
        END IF
        !
        CALL libtetrabz_tsmall_b2(e,V,tsmall)
        !
        IF(V > thr) THEN
           !
           ei2(1:4) = MATMUL(tsmall(1:4,1:4), ei1(sind(1:4)   ))
           ej2(1:4) = MATMUL(tsmall(1:4,1:4), ej1(sind(1:4),ib))
           CALL libtetrabz_polcmplx3(e0,ei2,ej2,w2)
           w1(1:ne,ib,sind(1:4)) = w1(1:ne,ib,         sind(1:4)) &
           &          + V * MATMUL(w2(1:ne,1:4), tsmall(1:4,1:4))
           !
        END IF
        !
        CALL libtetrabz_tsmall_b3(e,V,tsmall)
        !
        IF(V > thr) THEN
           !
           ei2(1:4) = MATMUL(tsmall(1:4,1:4), ei1(sind(1:4)   ))
           ej2(1:4) = MATMUL(tsmall(1:4,1:4), ej1(sind(1:4),ib))
           CALL libtetrabz_polcmplx3(e0,ei2,ej2,w2)
           w1(1:ne,ib,sind(1:4)) = w1(1:ne,ib,         sind(1:4)) &
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
           ei2(1:4) = MATMUL(tsmall(1:4,1:4), ei1(sind(1:4)   ))
           ej2(1:4) = MATMUL(tsmall(1:4,1:4), ej1(sind(1:4),ib))
           CALL libtetrabz_polcmplx3(e0,ei2,ej2,w2)
           w1(1:ne,ib,sind(1:4)) = w1(1:ne,ib,         sind(1:4)) &
           &          + V * MATMUL(w2(1:ne,1:4), tsmall(1:4,1:4))
           !
        END IF
        !
        CALL libtetrabz_tsmall_c2(e,V,tsmall)
        !
        IF(V > thr) THEN
           !
           ei2(1:4) = MATMUL(tsmall(1:4,1:4), ei1(sind(1:4)   ))
           ej2(1:4) = MATMUL(tsmall(1:4,1:4), ej1(sind(1:4),ib))
           CALL libtetrabz_polcmplx3(e0,ei2,ej2,w2)
           w1(1:ne,ib,sind(1:4)) = w1(1:ne,ib,         sind(1:4)) &
           &          + V * MATMUL(w2(1:ne,1:4), tsmall(1:4,1:4))
           !
        END IF
        !
        CALL libtetrabz_tsmall_c3(e,V,tsmall)
        !
        IF(V > thr) THEN
           !
           ei2(1:4) = MATMUL(tsmall(1:4,1:4), ei1(sind(1:4)   ))
           ej2(1:4) = MATMUL(tsmall(1:4,1:4), ej1(sind(1:4),ib))
           CALL libtetrabz_polcmplx3(e0,ei2,ej2,w2)
           w1(1:ne,ib,sind(1:4)) = w1(1:ne,ib,         sind(1:4)) &
           &          + V * MATMUL(w2(1:ne,1:4), tsmall(1:4,1:4))
           !
        END IF
        !
     ELSE IF( e(4) <= 0d0 ) THEN
        !
        ei2(1:4) = MATMUL(tsmall(1:4,1:4), ei1(sind(1:4)   ))
        ej2(1:4) = MATMUL(tsmall(1:4,1:4), ej1(sind(1:4),ib))
        CALL libtetrabz_polcmplx3(e0,ei2,ej2,w2)
        w1(1:ne,ib,sind(1:4)) = w1(1:ne,ib,         sind(1:4)) &
        &          + V * MATMUL(w2(1:ne,1:4), tsmall(1:4,1:4))
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
  USE libterabz_val, ONLY : ne
  IMPLICIT NONE
  !
  COMPLEX(8),INTENT(IN) :: e0(ne)
  REAL(8),INTENT(IN) :: ei1(4), ej1(4)
  COMPLEX(8),INTENT(OUT) :: w1(ne,4)
  !
  INTEGER :: ii, ie, ierr
  REAL(8) :: w2(2,4), de(4), e(4), thr
  !
  de(1:4) = ej1(1:4) - ei1(1:4)
  CALL libtetrabz_sort(4,de,sind)
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
     w1(ie,sind(1:4)) = CMPLX(w2(1,1:4) / AIMAG(e0(ie)), w2(2,1:4) / (- AIMAG(e0(ie))))
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
