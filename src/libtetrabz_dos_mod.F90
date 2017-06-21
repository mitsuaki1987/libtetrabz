MODULE libtetrabz_dos_mod
  !
  IMPLICIT NONE
  !
CONTAINS
!
! Compute DOS
!
SUBROUTINE libtetrabz_dos(ltetra,bvec,nb,nge,eig,ngw,wght0,ne,e0,comm0) BIND(C)
  !
#if defined(__MPI)
  USE mpi, ONLY : MPI_DOUBLE_PRECISION, MPI_IN_PLACE, MPI_SUM
#endif
  USE ISO_C_BINDING
  USE libtetrabz_val,    ONLY : comm, ik_global, ik_local, kvec, linterpol, lmpi, nk_local
  USE libtetrabz_common, ONLY : libtetrabz_initialize, libtetrabz_interpol_indx
  IMPLICIT NONE
  !
  INTEGER(C_INT),INTENT(IN) :: ltetra, nb, nge(3), ngw(3), ne
  REAL(C_DOUBLE),INTENT(IN) :: bvec(3,3), eig(nb,PRODUCT(nge(1:3))), e0(ne)
  REAL(C_DOUBLE),INTENT(OUT) :: wght0(ne,nb,PRODUCT(ngw(1:3)))
  INTEGER(C_INT),INTENT(IN),OPTIONAL :: comm0
  !
  INTEGER :: ik, ii, kintp(4)
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
  CALL libtetrabz_initialize(ltetra,bvec,nge,ngw,nb,ne)
  !
  IF(linterpol .OR. lmpi) THEN
     !
     ALLOCATE(wght1(ne,nb,nk_local))
     CALL libtetrabz_dos_main(eig,e0,wght1)
     !
     ! Interpolation
     !
     wght0(1:ne,1:nb,1:PRODUCT(ngw(1:3))) = 0d0
     DO ik = 1, nk_local
        CALL libtetrabz_interpol_indx(ngw,kvec(1:3,ik),kintp,wintp)
        DO ii = 1, 4
           wght0(1:ne,1:nb,kintp(ii)) = wght0(1:ne,1:nb,       kintp(ii)) &
           &                          + wght1(1:ne,1:nb, ik) * wintp(ii)
        END DO
     END DO ! ik = 1, nk_local
     DEALLOCATE(wght1, kvec)
     !
#if defined(__MPI)
     IF(lmpi) &
     &  CALL MPI_allREDUCE(MPI_IN_PLACE, wght0, ne * nb * PRODUCT(ngw(1:3)), &
     &                     MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)
#endif
     !
  ELSE
     CALL libtetrabz_dos_main(eig,e0,wght0)
  END IF
  !
  DEALLOCATE(ik_global, ik_local)
  !
END SUBROUTINE libtetrabz_dos
!
! Main SUBROUTINE for Dos : Delta(E - E1)
!
SUBROUTINE libtetrabz_dos_main(eig,e0,dos)
  !
  USE libtetrabz_val, ONLY : ik_global, ik_local, nb, ne, nkBZ, nk_local, nt_local, wlsm
  IMPLICIT NONE
  !
  REAL(8),INTENT(IN) :: eig(nb,nkBZ), e0(ne)
  REAL(8),INTENT(OUT) :: dos(ne,nb,nk_local)
  !
  INTEGER :: ib, it, ie, indx(4)
  REAL(8) :: e(4), ei1(4,nb), tsmall(3,4), V, w1(ne,4)
  !
  dos(1:ne, 1:nb, 1:nk_local) = 0d0
  !
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP & SHARED(dos,eig,e0,ik_global,ik_local,nb,ne,nt_local,wlsm) &
  !$OMP & PRIVATE(e,ei1,ib,ie,indx,it,tsmall,V,w1)
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
        w1(1:ne,1:4) = 0d0
        e(1:4) = ei1(1:4, ib)
        CALL libtetrabz_sort(4,e,indx)
        !
        DO ie = 1, ne
           !
           IF(e(1) < e0(ie) .AND. e0(ie) <= e(2)) THEN
              !
              CALL libtetrabz_triangle_a1(e(1:4) - e0(ie),V,tsmall)
              w1(ie,indx(1:4)) = w1(ie,indx(1:4)) + V * SUM(tsmall(1:3,1:4), 1) / 3d0
              !
           ELSE IF(e(2) < e0(ie) .AND. e0(ie) <= e(3)) THEN
              !
              CALL libtetrabz_triangle_b1(e(1:4) - e0(ie),V,tsmall)
              w1(ie,indx(1:4)) = w1(ie,indx(1:4)) + V * SUM(tsmall(1:3,1:4), 1) / 3d0
              !
              CALL libtetrabz_triangle_b1(e(1:4) - e0(ie),V,tsmall)
              w1(ie,indx(1:4)) = w1(ie,indx(1:4)) + V * SUM(tsmall(1:3,1:4), 1) / 3d0
              !
           ELSE IF(e(3) < e0(ie) .AND. e0(ie) < e(4)) THEN
              !
              CALL libtetrabz_triangle_c1(e(1:4) - e0(ie),V,tsmall)
              w1(ie,indx(1:4)) = w1(ie,indx(1:4)) + V * SUM(tsmall(1:3,1:4), 1) / 3d0
              !
           END IF
           !
        END DO ! ie
        !
        dos(1:ne,ib,ik_local(1:20,it)) = dos(1:ne,ib,    ik_local(1:20,it)) &
        &                        + MATMUL(w1(1:ne, 1:4), wlsm(1:4,1:20))
        !
     END DO ! ib
     !$OMP END DO NOWAIT
     !
  END DO ! it
  !
  !$OMP END PARALLEL
  !
  dos(1:ne,1:nb,1:nk_local) = dos(1:ne,1:nb,1:nk_local) / DBLE(6 * nkBZ)
  !
END SUBROUTINE libtetrabz_dos_main
!
END MODULE libtetrabz_dos_mod
