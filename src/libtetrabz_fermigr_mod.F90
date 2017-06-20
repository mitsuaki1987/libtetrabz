MODULE libtetrabz_fermigr_mod
  !
  IMPLICIT NONE
  !
CONTAINS
!
! Compute Fermi's golden rule
!
SUBROUTINE libtetrabz_fermigr(ltetra,comm0,bvec,nb,nge,eig1,eig2,ngw,wght0,ne,e0) BIND(C)
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
  REAL(C_DOUBLE),INTENT(IN) :: bvec(3,3), e0(ne)
  REAL(C_DOUBLE),INTENT(IN) :: eig1(nb,PRODUCT(nge(1:3))), eig2(nb,PRODUCT(nge(1:3)))
  REAL(C_DOUBLE),INTENT(OUT) :: wght0(ne,nb,nb,PRODUCT(ngw(1:3)))
  INTEGER(C_INT),INTENT(IN),OPTIONAL :: comm0
  !
  INTEGER :: ii, ik, kintp(4)
  REAL(8) :: wintp(4)
  REAL(8),ALLOCATABLE :: wght1(:,:,:,:)
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
     CALL libtetrabz_fermigr_main(eig1,eig2,e0,wght1)
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
     &                     MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)
#endif
     !
  ELSE
     CALL libtetrabz_fermigr_main(eig1,eig2,e0,wght0)
  END IF
  !
  DEALLOCATE(ik_global, ik_local)
  !
END SUBROUTINE libtetrabz_fermigr
!
! Main SUBROUTINE for Fermi's Gorlden rule : Theta(- E1) * Theta(E2) * Delta(E2 - E1 - w)
!
SUBROUTINE libtetrabz_fermigr_main(eig1,eig2,e0,fgr)
  !
  USE libtetrabz_val, ONLY : nb, nk_local, nt_local, wlsm, ik_global, ik_local, nkBZ, ne
  IMPLICIT NONE
  !
  REAL(8),INTENT(IN) :: eig1(nb,nkBZ), eig2(nb,nkBZ), e0(ne)
  REAL(8),INTENT(OUT) :: fgr(ne*nb,nb,nk_local)
  !
  INTEGER :: it, ib, indx(4)
  REAL(8) :: e(4), V, thr = 1d-10, tsmall(4,4), &
  &          ei1(4,nb), ei2(4), ej1(4,nb), ej2(4,nb), &
  &          w1(ne*nb,4), w2(ne*nb,4)
  !
  fgr(1:ne*nb,1:nb,1:nk_local) = 0d0
  !
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP & SHARED(nt_local,nb,ne,ik_global,ik_local,wlsm,eig1,eig2,fgr,thr,e0) &
  !$OMP & PRIVATE(ib,it,e,w1,w2,ei1,ej1,ei2,ej2,V,indx)
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
        w1(1:ne*nb,1:4) = 0d0
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
              CALL libtetrabz_fermigr2(e0,ei2,ej2,w2)
              w1(1:ne*nb,indx(1:4)) = w1(1:ne*nb,indx(1:4)) &
              &          + V * MATMUL(w2(1:ne*nb,     1:4 ), tsmall(1:4,1:4))
              !
           END IF
           !
        ELSE IF( e(2) <= 0d0 .AND. 0d0 < e(3)) THEN
           !
           CALL libtetrabz_tsmall_b1(e,V,tsmall)
           !
           IF(V > thr) THEN
              !
              ei2(1:4     ) = MATMUL(tsmall(1:4,1:4), ei1(indx(1:4),  ib))
              ej2(1:4,1:nb) = MATMUL(tsmall(1:4,1:4), ej1(indx(1:4),1:nb))
              CALL libtetrabz_fermigr2(e0,ei2,ej2,w2)
              w1(1:ne*nb,indx(1:4)) = w1(1:ne*nb,indx(1:4)) &
              &          + V * MATMUL(w2(1:ne*nb,     1:4 ), tsmall(1:4,1:4))
              !
           END IF
           !
           CALL libtetrabz_tsmall_b2(e,V,tsmall)
           !
           IF(V > thr) THEN
              !
              ei2(1:4     ) = MATMUL(tsmall(1:4,1:4), ei1(indx(1:4),  ib))
              ej2(1:4,1:nb) = MATMUL(tsmall(1:4,1:4), ej1(indx(1:4),1:nb))
              CALL libtetrabz_fermigr2(e0,ei2,ej2,w2)
              w1(1:ne*nb,indx(1:4)) = w1(1:ne*nb,indx(1:4)) &
              &          + V * MATMUL(w2(1:ne*nb,     1:4 ), tsmall(1:4,1:4))
              !
           END IF
           !
           CALL libtetrabz_tsmall_b3(e,V,tsmall)
           !
           IF(V > thr) THEN
              !
              ei2(1:4     ) = MATMUL(tsmall(1:4,1:4), ei1(indx(1:4),  ib))
              ej2(1:4,1:nb) = MATMUL(tsmall(1:4,1:4), ej1(indx(1:4),1:nb))
              CALL libtetrabz_fermigr2(e0,ei2,ej2,w2)
              w1(1:ne*nb,indx(1:4)) = w1(1:ne*nb,indx(1:4)) &
              &          + V * MATMUL(w2(1:ne*nb,     1:4 ), tsmall(1:4,1:4))
              !
           END IF
           !
        ELSE IF( e(3) <= 0d0 .AND. 0d0 < e(4)) THEN
           !
           CALL libtetrabz_tsmall_c1(e,V,tsmall)
           !
           IF(V > thr) THEN
              !
              ei2(1:4     ) = MATMUL(tsmall(1:4,1:4), ei1(indx(1:4),  ib))
              ej2(1:4,1:nb) = MATMUL(tsmall(1:4,1:4), ej1(indx(1:4),1:nb))
              CALL libtetrabz_fermigr2(e0,ei2,ej2,w2)
              w1(1:ne*nb,indx(1:4)) = w1(1:ne*nb,indx(1:4)) &
              &          + V * MATMUL(w2(1:ne*nb,     1:4 ), tsmall(1:4,1:4))
              !
           END IF
           !
           CALL libtetrabz_tsmall_c2(e,V,tsmall)
           !
           IF(V > thr) THEN
              !
              ei2(1:4     ) = MATMUL(tsmall(1:4,1:4), ei1(indx(1:4),  ib))
              ej2(1:4,1:nb) = MATMUL(tsmall(1:4,1:4), ej1(indx(1:4),1:nb))
              CALL libtetrabz_fermigr2(e0,ei2,ej2,w2)
              w1(1:ne*nb,indx(1:4)) = w1(1:ne*nb,indx(1:4)) &
              &          + V * MATMUL(w2(1:ne*nb,     1:4 ), tsmall(1:4,1:4))
              !
           END IF
           !
           CALL libtetrabz_tsmall_c3(e,V,tsmall)
           !
           IF(V > thr) THEN
              !
              ei2(1:4     ) = MATMUL(tsmall(1:4,1:4), ei1(indx(1:4),  ib))
              ej2(1:4,1:nb) = MATMUL(tsmall(1:4,1:4), ej1(indx(1:4),1:nb))
              CALL libtetrabz_fermigr2(e0,ei2,ej2,w2)
              w1(1:ne*nb,indx(1:4)) = w1(1:ne*nb,indx(1:4)) &
              &          + V * MATMUL(w2(1:ne*nb,     1:4 ), tsmall(1:4,1:4))
              !
           END IF
           !
        ELSE IF(e(4) <= 0d0) THEN
           !
           ei2(1:4     ) = MATMUL(tsmall(1:4,1:4), ei1(1:4,  ib))
           ej2(1:4,1:nb) = MATMUL(tsmall(1:4,1:4), ej1(1:4,1:nb))
           CALL libtetrabz_fermigr2(e0,ei2,ej2,w2)
           w1(1:ne*nb,1:4) = w1(1:ne*nb,1:4) &
           &    + V * MATMUL(w2(1:ne*nb,1:4), tsmall(1:4,1:4))
           !
        END IF
        !
        fgr(1:ne*nb,ib,ik_local(1:20,it)) = fgr(1:ne*nb,ib,   ik_local(1:20,it)) &
        &                           + MATMUL(w1(1:ne*nb,1:4), wlsm(1:4,1:20))
        !
     END DO ! ib = 1, nb
     !$OMP END DO NOWAIT
     !
  END DO ! it
  !
  !$OMP END PARALLEL
  !
  fgr(1:ne*nb,1:nb,1:nk_local) = fgr(1:ne*nb,1:nb,1:nk_local) / DBLE(6 * nkBZ)
  !
END SUBROUTINE libtetrabz_fermigr_main
!
! Tetrahedra method for theta( - E2)
!
SUBROUTINE libtetrabz_fermigr2(e0,ei1,ej1,w1)
  !
  USE libtetrabz_val, ONLY : nb, ne
  IMPLICIT NONE
  !
  REAL(8),INTENT(IN) :: e0(ne), ei1(4), ej1(4,nb)
  REAL(8),INTENT(OUT) :: w1(ne,nb,4)
  !
  INTEGER :: ib, indx(4)
  REAL(8) :: V, ei2(4), ej2(4), w2(ne,4), thr = 1d-8, e(4), tsmall(4,4)
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
           CALL libtetrabz_fermigr3(e0,ei2,ej2,w2)
           w1(1:ne,ib,indx(1:4)) = w1(1:ne,ib,indx(1:4)) &
           &          + V * MATMUL(w2(1:ne,        1:4 ), tsmall(1:4,1:4))
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
           CALL libtetrabz_fermigr3(e0,ei2,ej2,w2)
           w1(1:ne,ib,indx(1:4)) = w1(1:ne,ib,indx(1:4)) &
           &          + V * MATMUL(w2(1:ne,        1:4 ), tsmall(1:4,1:4))
           !
        END IF
        !
        CALL libtetrabz_tsmall_b2(e,V,tsmall)
        !
        IF(V > thr) THEN
           !
           ei2(1:4) = MATMUL(tsmall(1:4,1:4), ei1(indx(1:4)   ))
           ej2(1:4) = MATMUL(tsmall(1:4,1:4), ej1(indx(1:4),ib))
           CALL libtetrabz_fermigr3(e0,ei2,ej2,w2)
           w1(1:ne,ib,indx(1:4)) = w1(1:ne,ib,indx(1:4)) &
           &          + V * MATMUL(w2(1:ne,        1:4 ), tsmall(1:4,1:4))
           !
        END IF
        !
        CALL libtetrabz_tsmall_b3(e,V,tsmall)
        !
        IF(V > thr) THEN
           !
           ei2(1:4) = MATMUL(tsmall(1:4,1:4), ei1(indx(1:4)   ))
           ej2(1:4) = MATMUL(tsmall(1:4,1:4), ej1(indx(1:4),ib))
           CALL libtetrabz_fermigr3(e0,ei2,ej2,w2)
           w1(1:ne,ib,indx(1:4)) = w1(1:ne,ib,indx(1:4)) &
           &          + V * MATMUL(w2(1:ne,        1:4 ), tsmall(1:4,1:4))
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
           CALL libtetrabz_fermigr3(e0,ei2,ej2,w2)
           w1(1:ne,ib,indx(1:4)) = w1(1:ne,ib,indx(1:4)) &
           &          + V * MATMUL(w2(1:ne,        1:4 ), tsmall(1:4,1:4))
           !
        END IF
        !
        CALL libtetrabz_tsmall_c2(e,V,tsmall)
        !
        IF(V > thr) THEN
           !
           ei2(1:4) = MATMUL(tsmall(1:4,1:4), ei1(indx(1:4)   ))
           ej2(1:4) = MATMUL(tsmall(1:4,1:4), ej1(indx(1:4),ib))
           CALL libtetrabz_fermigr3(e0,ei2,ej2,w2)
           w1(1:ne,ib,indx(1:4)) = w1(1:ne,ib,indx(1:4)) &
           &          + V * MATMUL(w2(1:ne,        1:4 ), tsmall(1:4,1:4))
           !
        END IF
        !
        CALL libtetrabz_tsmall_c3(e,V,tsmall)
        !
        IF(V > thr) THEN
           !
           ei2(1:4) = MATMUL(tsmall(1:4,1:4), ei1(indx(1:4)   ))
           ej2(1:4) = MATMUL(tsmall(1:4,1:4), ej1(indx(1:4),ib))
           CALL libtetrabz_fermigr3(e0,ei2,ej2,w2)
           w1(1:ne,ib,indx(1:4)) = w1(1:ne,ib,indx(1:4)) &
           &          + V * MATMUL(w2(1:ne,        1:4 ), tsmall(1:4,1:4))
           !
        END IF
        !
     ELSE IF(e(4) <= 0d0) THEN
        !
        ei2(1:4) = ei1(indx(1:4)   )
        ej2(1:4) = ej1(indx(1:4),ib)
        CALL libtetrabz_fermigr3(e0,ei2,ej2,w2)
        w1(1:ne,ib,1:4) = w1(1:ne,ib,1:4) + w2(1:ne,1:4)
        !
     END IF
     !
  END DO ! ib = 1, nb
  !
END SUBROUTINE libtetrabz_fermigr2
!
!
!
SUBROUTINE libtetrabz_fermigr3(e0,ei1,ej1,w1)
  !
  USE libtetrabz_val, ONLY : ne
  IMPLICIT NONE
  !
  REAL(8),INTENT(IN) :: e0(ne), ei1(4), ej1(4)
  REAL(8),INTENT(OUT) :: w1(ne,4)
  !
  INTEGER :: ie, indx(4)
  REAL(8) :: e(4), V, tsmall(3,4), w2(3)
  !
  w2(1:3) = 1d0 / 3d0
  !
  w1(1:ne,1:4) = 0d0
  e(1:4) = ej1(1:4) - ei1(1:4)
  CALL libtetrabz_sort(4,e,indx)
  !
  DO ie = 1, ne
     !
     IF(e(1) < e0(ie) .AND. e0(ie) <= e(2)) THEN
        !
        CALL libtetrabz_triangle_a1(e(1:4) - e0(ie),V,tsmall)
        w1(ie,indx(1:4)) = w1(ie,indx(1:4)) &
        &  + V * MATMUL(w2(1:3), tsmall(1:3,1:4))
        !
     ELSE IF(e(2) < e0(ie) .AND. e0(ie) <= e(3)) THEN
        !
        CALL libtetrabz_triangle_b1(e(1:4) - e0(ie),V,tsmall)
        w1(ie,indx(1:4)) = w1(ie,indx(1:4)) &
        &  + V * MATMUL(w2(1:3), tsmall(1:3,1:4))
        !
        CALL libtetrabz_triangle_b2(e(1:4) - e0(ie),V,tsmall)
        w1(ie,indx(1:4)) = w1(ie,indx(1:4)) &
        &  + V * MATMUL(w2(1:3), tsmall(1:3,1:4))
        !
     ELSE IF(e(3) < e0(ie) .AND. e0(ie) < e(4)) THEN
        !
        !
        CALL libtetrabz_triangle_c1(e(1:4) - e0(ie),V,tsmall)
        w1(ie,indx(1:4)) = w1(ie,indx(1:4)) &
        &  + V * MATMUL(w2(1:3), tsmall(1:3,1:4))
        !
     END IF
     !
  END DO ! ie
  !
END SUBROUTINE libtetrabz_fermigr3
!
END MODULE libtetrabz_fermigr_mod
