MODULE libtetrabz_dos_mod
  !
  IMPLICIT NONE
  !
CONTAINS
!
! Compute DOS
!
SUBROUTINE libtetrabz_mpi_dos(ltetra0,comm0,bvec,nb0,nge,eig,ngw,wght0,ne0,e0)
  !
  USE mpi, ONLY : MPI_IN_PLACE, MPI_DOUBLE_PRECISION, MPI_SUM
  USE libtetrabz_vals, ONLY : ltetra, ng, nb, nk0, indx1, indx2, indx3, ne
  USE libtetrabz_common, ONLY : libtetrabz_initialize, libtetrabz_interpol_weight
  USE libtetrabz_dos_mod, ONLY : libtetrabz_dos_main
  !
  USE libtetrabz_mpi_routines, ONLY : comm, libtetrabz_mpi_kgrid
  !
  INTEGER,INTENT(IN) :: ltetra0, comm0, nge(3), ngw(3), nb0, ne0
  REAL(8),INTENT(IN) :: bvec(3,3), e0(ne0)
  REAL(8),INTENT(IN) :: eig(nb0,PRODUCT(nge(1:3)))
  REAL(8),INTENT(OUT) :: wght0(ne0,nb0,PRODUCT(ngw(1:3)))
  !
  INTEGER :: ierr, nn
  REAL(8),ALLOCATABLE :: wght1(:,:)
  !
  ltetra = ltetra0
  comm = comm0
  nb = nb0
  ng(1:3) = nge(1:3)
  ne = ne0
  nn = ne * nb
  !
  CALL libtetrabz_initialize(bvec)
  CALL libtetrabz_mpi_kgrid()
  !
  ALLOCATE(wght1(nn, nk0))
  CALL libtetrabz_dos_main(eig,e0,wght1)
  !
  CALL libtetrabz_interpol_weight(nn,ngw,nge,wght0,wght1)
  !
  DEALLOCATE(wght1, indx1, indx2, indx3)
  !
  CALL MPI_allREDUCE(MPI_IN_PLACE, wght0, nn * PRODUCT(ngw(1:3)), &
  &                  MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)
  !
END SUBROUTINE libtetrabz_mpi_dos
!
! Main SUBROUTINE for Dos : Delta(E - E1)
!
SUBROUTINE libtetrabz_dos_main(eig,e0,dos)
  !
  USE libterabz_val, ONLY : nb, nk, nk0, ne, fst, lst, wlsm, indx1, indx2
  IMPLICIT NONE
  !
  REAL(8),INTENT(IN) :: eig(nb,nk), e0(ne)
  REAL(8),INTENT(OUT) :: dos(ne,nb,nk0)
  !
  INTEGER :: ib, it, ie, sind(4)
  REAL(8) :: e(4), tmp(5,4), V, w1(ne,4), w2(3), ei1(4,nb)
  !
  dos(1:ne, 1:nb, 1:nk0) = 0d0
  w2(1:3) = 1d0 / 3d0
  !
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP & SHARED(fst,lst,nb,ne,indx1,indx2,wlsm,eig,w2,e0,dos) &
  !$OMP & PRIVATE(ib,it,ie,e,w1,ei1,V,sind)
  !
  DO it = fst, lst
     !
     DO ib = 1, nb
        ei1(1:4,ib) = MATMUL(wlsm(1:4,1:20), eig(ib,indx1(1:20,it)))
     END DO
     !
     !$OMP DO
     DO ib = 1, nb
        !
        w1(1:ne,1:4) = 0d0
        e(1:4) = ei1(1:4, ib)
        CALL libtetrabz_sort(4,e,sind)
        !
        DO ie = 1, ne
           !
           IF(e(1) < e0(ie) .AND. e0(ie) <= e(2)) THEN
              !
              CALL libtetrabz_triangle_a1(e(1:4) - e0(ie),V,tsmall)
              w1(ie,sind(1:4)) = w1(ie,sind(1:4)) &
              &  + V * MATMUL(w2(1:3), tsmall(1:3,1:4))
              !
           ELSE IF(e(2) < e0(ie) .AND. e0(ie) <= e(3)) THEN
              !
              CALL libtetrabz_triangle_b1(e(1:4) - e0(ie),V,tsmall)
              w1(ie,sind(1:4)) = w1(ie,sind(1:4)) &
              &  + V * MATMUL(w2(1:3), tsmall(1:3,1:4))
              !
              CALL libtetrabz_triangle_b1(e(1:4) - e0(ie),V,tsmall)
              w1(ie,sind(1:4)) = w1(ie,sind(1:4)) &
              &  + V * MATMUL(w2(1:3), tsmall(1:3,1:4))
              !
           ELSE IF(e(3) < e0(ie) .AND. e0(ie) < e(4)) THEN
              !
              CALL libtetrabz_triangle_c1(e(1:4) - e0(ie),V,tsmall)
              w1(ie,sind(1:4)) = w1(ie,sind(1:4)) &
              &  + V * MATMUL(w2(1:3), tsmall(1:3,1:4))
              !
           END IF
           !
        END DO ! ie
        !
        dos(1:ne,ib,indx2(1:20,it)) = dos(1:ne,ib,       indx2(1:20,it)) &
        &                     + MATMUL(w1(1:ne, 1:4), wlsm(1:4,1:20))
        !
     END DO ! ib
     !$OMP END DO NOWAIT
     !
  END DO ! it
  !
  !$OMP END PARALLEL
  !
  dos(1:ne,1:nb,1:nk0) = dos(1:ne,1:nb,1:nk0) / DBLE(6 * nk)
  !
END SUBROUTINE libtetrabz_dos_main
!
END MODULE libtetrabz_dos_mod
