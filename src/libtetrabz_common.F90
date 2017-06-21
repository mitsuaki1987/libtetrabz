MODULE libtetrabz_common
  !
  IMPLICIT NONE
  !
CONTAINS
!
! define shortest diagonal line & define type of tetragonal
!
SUBROUTINE libtetrabz_initialize(ltetra,bvec,nge,ngw,nb0,ne0)
  !
  USE libtetrabz_val, ONLY : wlsm, nb, ne, ng, linterpol, nkBZ
  IMPLICIT NONE
  !
  INTEGER,INTENT(IN) :: ltetra, nb0, nge(3), ngw(3)
  REAL(8),INTENT(IN) :: bvec(3,3)
  INTEGER,INTENT(IN),OPTIONAL :: ne0
  !
  INTEGER :: itype, i1, i2, i3, it, divvec(4,4), ivvec0(4), ivvec(3,20,6)
  REAL(8) :: l(4), bvec2(3,3), bvec3(3,4)
  !
  nb = nb0
  ng(1:3) = nge(1:3)
  nkBZ = PRODUCT(ng(1:3))
  linterpol = .NOT. ALL(nge(1:3) == ngw(1:3))
  IF(PRESENT(ne0)) ne = ne0
  !
  DO i1 = 1, 3
     bvec2(1:3,i1) = bvec(1:3,i1) / DBLE(ng(i1))
  END DO
  !
  bvec3(1:3,1) = -bvec2(1:3,1) + bvec2(1:3,2) + bvec2(1:3,3)
  bvec3(1:3,2) =  bvec2(1:3,1) - bvec2(1:3,2) + bvec2(1:3,3)
  bvec3(1:3,3) =  bvec2(1:3,1) + bvec2(1:3,2) - bvec2(1:3,3)
  bvec3(1:3,4) =  bvec2(1:3,1) + bvec2(1:3,2) + bvec2(1:3,3)
  !
  ! length of delta bvec
  !
  DO i1 = 1, 4
     l(i1) = DOT_PRODUCT(bvec3(1:3,i1),bvec3(1:3,i1))
  END DO
  !
  itype = MINLOC(l(1:4),1)
  !
  ! start & last
  !
  ivvec0(1:4) = (/ 0, 0, 0, 0 /)
  !
  divvec(1:4,1) = (/ 1, 0, 0, 0 /)
  divvec(1:4,2) = (/ 0, 1, 0, 0 /)
  divvec(1:4,3) = (/ 0, 0, 1, 0 /)
  divvec(1:4,4) = (/ 0, 0, 0, 1 /)
  !
  ivvec0(itype) = 1
  divvec(itype, itype) = - 1
  !
  ! Corners of tetrahedra
  !
  it = 0
  DO i1 = 1, 3
     DO i2 = 1, 3
        IF(i2 == i1) CYCLE
        DO i3 = 1, 3
           IF(i3 == i1 .OR. i3 == i2) CYCLE
           !
           it = it + 1
           !
           ivvec(1:3,1,it) = ivvec0(1:3)
           ivvec(1:3,2,it) = ivvec(1:3,1,it) + divvec(1:3,i1)
           ivvec(1:3,3,it) = ivvec(1:3,2,it) + divvec(1:3,i2)
           ivvec(1:3,4,it) = ivvec(1:3,3,it) + divvec(1:3,i3)
           !
        END DO
     END DO
  END DO
  !
  ! Additional points
  !
  ivvec(1:3, 5,1:6) = 2 * ivvec(1:3,1,1:6) - ivvec(1:3,2,1:6)
  ivvec(1:3, 6,1:6) = 2 * ivvec(1:3,2,1:6) - ivvec(1:3,3,1:6)
  ivvec(1:3, 7,1:6) = 2 * ivvec(1:3,3,1:6) - ivvec(1:3,4,1:6)
  ivvec(1:3, 8,1:6) = 2 * ivvec(1:3,4,1:6) - ivvec(1:3,1,1:6)
  !
  ivvec(1:3, 9,1:6) = 2 * ivvec(1:3,1,1:6) - ivvec(1:3,3,1:6)
  ivvec(1:3,10,1:6) = 2 * ivvec(1:3,2,1:6) - ivvec(1:3,4,1:6)
  ivvec(1:3,11,1:6) = 2 * ivvec(1:3,3,1:6) - ivvec(1:3,1,1:6)
  ivvec(1:3,12,1:6) = 2 * ivvec(1:3,4,1:6) - ivvec(1:3,2,1:6)
  !
  ivvec(1:3,13,1:6) = 2 * ivvec(1:3,1,1:6) - ivvec(1:3,4,1:6)
  ivvec(1:3,14,1:6) = 2 * ivvec(1:3,2,1:6) - ivvec(1:3,1,1:6)
  ivvec(1:3,15,1:6) = 2 * ivvec(1:3,3,1:6) - ivvec(1:3,2,1:6)
  ivvec(1:3,16,1:6) = 2 * ivvec(1:3,4,1:6) - ivvec(1:3,3,1:6)
  !
  ivvec(1:3,17,1:6) =  ivvec(1:3,4,1:6) - ivvec(1:3,1,1:6) + ivvec(1:3,2,1:6)
  ivvec(1:3,18,1:6) =  ivvec(1:3,1,1:6) - ivvec(1:3,2,1:6) + ivvec(1:3,3,1:6)
  ivvec(1:3,19,1:6) =  ivvec(1:3,2,1:6) - ivvec(1:3,3,1:6) + ivvec(1:3,4,1:6)
  ivvec(1:3,20,1:6) =  ivvec(1:3,3,1:6) - ivvec(1:3,4,1:6) + ivvec(1:3,1,1:6)
  !
  IF(ltetra == 1) THEN
     !
     !WRITE(*,*) "[libtetrabz] Linear tetrahedron method is used."
     !
     wlsm(1:4,1:20) = 0.0d0
     wlsm(1,1) = 1.0d0
     wlsm(2,2) = 1.0d0
     wlsm(3,3) = 1.0d0
     wlsm(4,4) = 1.0d0
     !
  ELSE IF(ltetra == 2) THEN
     !
     !WRITE(*,*) "[libtetrabz] Improved tetrahedron method is used."
     !
     wlsm(1, 1: 4) = DBLE((/1440,    0,   30,    0/))
     wlsm(2, 1: 4) = DBLE((/   0, 1440,    0,   30/))
     wlsm(3, 1: 4) = DBLE((/  30,    0, 1440,    0/))
     wlsm(4, 1: 4) = DBLE((/   0,   30,    0, 1440/))
     !
     wlsm(1, 5: 8) = DBLE((/ -38,    7,   17,  -28/))
     wlsm(2, 5: 8) = DBLE((/ -28,  -38,    7,   17/))
     wlsm(3, 5: 8) = DBLE((/  17,  -28,  -38,    7/))
     wlsm(4, 5: 8) = DBLE((/   7,   17,  -28,  -38/))
     !
     wlsm(1, 9:12) = DBLE((/ -56,    9,  -46,    9/))
     wlsm(2, 9:12) = DBLE((/   9,  -56,    9,  -46/))
     wlsm(3, 9:12) = DBLE((/ -46,    9,  -56,    9/))
     wlsm(4, 9:12) = DBLE((/   9,  -46,    9,  -56/))
     !
     wlsm(1,13:16) = DBLE((/ -38,  -28,   17,    7/))
     wlsm(2,13:16) = DBLE((/   7,  -38,  -28,   17/))
     wlsm(3,13:16) = DBLE((/  17,    7,  -38,  -28/))
     wlsm(4,13:16) = DBLE((/ -28,   17,    7,  -38/))
     !
     wlsm(1,17:20) = DBLE((/ -18,  -18,   12,  -18/))
     wlsm(2,17:20) = DBLE((/ -18,  -18,  -18,   12/))
     wlsm(3,17:20) = DBLE((/  12,  -18,  -18,  -18/))
     wlsm(4,17:20) = DBLE((/ -18,   12,  -18,  -18/))
     !
     wlsm(1:4,1:20) = wlsm(1:4,1:20) / 1260d0
     !
  ELSE
     !
     WRITE(*,*) "[libtetrabz] STOP! ltetrta is invalid."
     STOP
     !
  END IF
  !
  CALL libtetrabz_kgrid(ivvec)
  !
END SUBROUTINE libtetrabz_initialize
!
! Initialize grid
!
SUBROUTINE libtetrabz_kgrid(ivvec)
  !
  USE libtetrabz_val, ONLY : nk_local, ik_global, ik_local, kvec, ng, nt_local, lmpi, linterpol, nkBZ
  IMPLICIT NONE
  !
  INTEGER,INTENT(IN) :: ivvec(3,20,6)
  !
  INTEGER :: it, i1, i2, i3, ii, ikv(3), nt, ik, nt_front, loc2glob(nkBZ)
  !
  CALL libtetrabz_fst_and_lst(6 * nkBZ, nt_front, nt_local)
  ALLOCATE(ik_global(20, nt_local), ik_local(20, nt_local))
  !
  ! k-index for energy (Global index)
  !
  nt = 0
  DO i3 = 1, ng(3)
     DO i2 = 1, ng(2)
        DO i1 = 1, ng(1)
           !
           DO it = 1, 6
              !
              nt = nt + 1
              IF(nt <= nt_front .OR. nt_front + nt_local < nt) CYCLE
              !
              DO ii = 1, 20
                 !
                 ikv(1:3) = (/i1, i2, i3/) + ivvec(1:3,ii,it) - 1
                 ikv(1:3) = MODULO(ikv(1:3), ng(1:3))
                 !
                 ik_global(ii,nt - nt_front) = 1 + ikv(1) + ng(1) * ikv(2) + ng(1) * ng(2) * ikv(3)
                 !
              END DO
              !
           END DO
           !
        END DO
     END DO
  END DO
  !
  ! k-index for weight (Local index)
  !
  ik_local(1:20,1:nt_local) = ik_global(1:20,1:nt_local)
  IF((.NOT. lmpi) .AND. (.NOT. linterpol)) THEN
     nk_local = PRODUCT(ng(1:3))
     RETURN
  END IF
  !
  nk_local = 0
  DO nt = 1, nt_local
     DO ii = 1, 20
        !
        IF(ik_local(ii,nt) <= nk_local) CYCLE
        !
        nk_local = nk_local + 1
        loc2glob(nk_local) = ik_local(ii,nt)
        WHERE(ik_local(1:20,1:nt_local) == loc2glob(nk_local)) &
        &  ik_local(1:20,1:nt_local) = nk_local
        !
     END DO
  END DO
  !
  ! k-vector in the fractional coordinate
  !
  ALLOCATE(kvec(3,nk_local))
  DO ik = 1, nk_local
     ! loc2glob(ik) - 1 = i1 + ng(1) * i2 + ng(1) * ng(2) * i3
     i1 = MOD(loc2glob(ik) - 1, ng(1))
     i2 = MOD((loc2glob(ik) - 1) / ng(1), ng(2))
     i3 = (loc2glob(ik) - 1) / (ng(1) * ng(2))
     kvec(1:3,ik) = DBLE((/i1, i2, i3/)) / DBLE(ng(1:3))
  END DO
  !
END SUBROUTINE libtetrabz_kgrid
!
! Compute cnt and dsp
!
SUBROUTINE libtetrabz_fst_and_lst(nt,nt_front,nt_local)
  !
#if defined(__MPI)
  USE mpi, ONLY : mpi_comm_size, mpi_comm_rank
  USE libtetrabz_val, ONLY : comm, lmpi
#endif
  IMPLICIT NONE
  !
  INTEGER,INTENT(IN) :: nt
  INTEGER,INTENT(OUT) :: nt_front, nt_local
  !
  INTEGER :: petot, my_rank
#if defined(__MPI)
  INTEGER :: ierr
#endif
  !
  petot = 1
  my_rank = 0
#if defined(__MPI)
  IF(lmpi) THEN
     CALL MPI_COMM_SIZE(comm, petot, ierr)
     CALL MPI_COMM_RANK(comm, my_rank, ierr)
  END if
#endif
  !
  IF(my_rank < MOD(nt, petot)) THEN
     nt_local = nt / petot + 1
     nt_front = my_rank * nt_local
  ELSE
     nt_local = nt / petot
     nt_front = my_rank * nt_local + MOD(nt, petot)
  END IF
  !
END SUBROUTINE libtetrabz_fst_and_lst
!
! Simple sort
!
SUBROUTINE libtetrabz_sort(n,key,indx)
  !
  IMPLICIT NONE
  !
  integer,INTENT(IN) :: n
  REAL(8),INTENT(inout) :: key(n)
  INTEGER,INTENT(OUT) :: indx(n)
  !
  INTEGER :: i, i0, indx0
  REAL(8) :: key0
  !
  DO i = 1, n
     indx(i) = i
  END DO
  !
  DO i = 1, n - 1
     key0 = MINVAL(key(i+1:n))
     i0   = MINLOC(key(i+1:n),1) + i
     IF(key(i) > key0) THEN
        key(i0) = key(i)
        key(i) = key0
        !
        indx0 = indx(i0)
        indx(i0) = indx(i)
        indx(i) = indx0
     END IF
  END DO
  !
END SUBROUTINE libtetrabz_sort
!
! Linear interpolation
!
SUBROUTINE libtetrabz_interpol_indx(ng,kvec,kintp,wintp)
  !
  IMPLICIT NONE
  !
  INTEGER,INTENT(in) :: ng(3)
  REAL(8),INTENT(in) :: kvec(3)
  INTEGER,INTENT(out) :: kintp(4)
  REAL(8),INTENT(out) :: wintp(4)
  !
  INTEGER :: ikv0(3), ikv1(3), dikv(3), ii
  REAL(8) :: x(3)
  !
  ! Search nearest neighbor grid points.
  !
  x(1:3) = kvec(1:3) * DBLE(ng(1:3))
  ikv0(1:3) = NINT(x(1:3))
  dikv(1:3) = ikv0(1:3) - FLOOR(x(1:3))
  dikv(1:3) = 1 - 2 * dikv(1:3)
  x(1:3) = ABS(x(1:3) - DBLE(ikv0(1:3)))
  !
  ! Interpolation k & weights
  !
  DO ii = 1, 3
     wintp(ii) = x(ii)
     ikv1(1:3) = ikv0(1:3)
     ikv1(ii) = ikv1(ii) + dikv(ii)
     kintp(ii) = 1 + ikv1(1) + ng(1) * ikv1(2) + ng(1) * ng(2) * ikv1(3)
  END DO
  wintp(4) = 1d0 - SUM(x(1:3))
  kintp(4) = 1 + ikv0(1) + ng(1) * ikv0(2) + ng(1) * ng(2) * ikv0(3)
  !
end subroutine libtetrabz_interpol_indx
!
! Cut small tetrahedron A1
!
SUBROUTINE libtetrabz_tsmall_a1(e,V,tsmall)
  !
  IMPLICIT NONE
  !
  REAL(8),INTENT(IN) :: e(4)
  REAL(8),INTENT(OUT) :: V
  REAL(8),INTENT(OUT) :: tsmall(4,4)
  !
  INTEGER :: ii
  REAL(8) :: a(4,4)
  !
  DO ii = 1, 4
     a(1:4,ii) = (0d0 - e(ii)) / (e(1:4) - e(ii))
  END DO
  !
  V = 0.25d0 * a(2,1) * a(3,1) * a(4,1)
  !
  tsmall(1, 1:4) = (/   1d0,    0d0,    0d0,    0d0/)
  tsmall(2, 1:4) = (/a(1,2), a(2,1),    0d0,    0d0/)
  tsmall(3, 1:4) = (/a(1,3),    0d0, a(3,1),    0d0/)
  tsmall(4, 1:4) = (/a(1,4),    0d0,    0d0, a(4,1)/)
  !
END SUBROUTINE libtetrabz_tsmall_a1
!
! Cut small tetrahedron B1
!
SUBROUTINE libtetrabz_tsmall_b1(e,V,tsmall)
  !
  IMPLICIT NONE
  !
  REAL(8),INTENT(IN) :: e(4)
  REAL(8),INTENT(OUT) :: V
  REAL(8),INTENT(OUT) :: tsmall(4,4)
  !
  INTEGER :: ii
  REAL(8) :: a(4,4)
  !
  DO ii = 1, 4
     a(1:4,ii) = (0d0 - e(ii)) / (e(1:4) - e(ii))
  END DO
  !
  V = a(3,1) * a(4,1) * a(2,4)
  !
  tsmall(1, 1:4) = (/   1d0,    0d0,    0d0,    0d0/)
  tsmall(2, 1:4) = (/a(1,3),    0d0, a(3,1),    0d0/)
  tsmall(3, 1:4) = (/a(1,4),    0d0,    0d0, a(4,1)/)
  tsmall(4, 1:4) = (/   0d0, a(2,4),    0d0, a(4,2)/)
  !
END SUBROUTINE libtetrabz_tsmall_b1
!
! Cut small tetrahedron B2
!
SUBROUTINE libtetrabz_tsmall_b2(e,V,tsmall)
  !
  IMPLICIT NONE
  !
  REAL(8),INTENT(IN) :: e(4)
  REAL(8),INTENT(OUT) :: V
  REAL(8),INTENT(OUT) :: tsmall(4,4)
  !
  INTEGER :: ii
  REAL(8) :: a(4,4)
  !
  DO ii = 1, 4
     a(1:4,ii) = (0d0 - e(ii)) / (e(1:4) - e(ii))
  END DO
  !
  V = a(3,2) * a(4,2)
  !
  tsmall(1, 1:4) = (/1d0,    0d0,    0d0,    0d0/)
  tsmall(2, 1:4) = (/0d0,    1d0,    0d0,    0d0/)
  tsmall(3, 1:4) = (/0d0, a(2,3), a(3,2),    0d0/)
  tsmall(4, 1:4) = (/0d0, a(2,4),    0d0, a(4,2)/)
  !
END SUBROUTINE libtetrabz_tsmall_b2
!
! Cut small tetrahedron B3
!
SUBROUTINE libtetrabz_tsmall_b3(e,V,tsmall)
  !
  IMPLICIT NONE
  !
  REAL(8),INTENT(IN) :: e(4)
  REAL(8),INTENT(OUT) :: V
  REAL(8),INTENT(OUT) :: tsmall(4,4)
  !
  INTEGER :: ii
  REAL(8) :: a(4,4)
  !
  DO ii = 1, 4
     a(1:4,ii) = (0d0 - e(ii)) / (e(1:4) - e(ii))
  END DO
  !
  V = a(2,3) * a(3,1) * a(4,2)
  !
  tsmall(1, 1:4) = (/   1d0,    0d0,    0d0,    0d0/)
  tsmall(2, 1:4) = (/a(1,3),    0d0, a(3,1),    0d0/)
  tsmall(3, 1:4) = (/   0d0, a(2,3), a(3,2),    0d0/)
  tsmall(4, 1:4) = (/   0d0, a(2,4),    0d0, a(4,2)/)
  !
END SUBROUTINE libtetrabz_tsmall_b3
!
! Cut small tetrahedron C1
!
SUBROUTINE libtetrabz_tsmall_c1(e,V,tsmall)
  !
  IMPLICIT NONE
  !
  REAL(8),INTENT(IN) :: e(4)
  REAL(8),INTENT(OUT) :: V
  REAL(8),INTENT(OUT) :: tsmall(4,4)
  !
  INTEGER :: ii
  REAL(8) :: a(4,4)
  !
  DO ii = 1, 4
     a(1:4,ii) = (0d0 - e(ii)) / (e(1:4) - e(ii))
  END DO
  !
  V = a(4,3)
  !
  tsmall(1, 1:4) = (/1d0, 0d0,    0d0,    0d0/)
  tsmall(2, 1:4) = (/0d0, 1d0,    0d0,    0d0/)
  tsmall(3, 1:4) = (/0d0, 0d0,    1d0,    0d0/)
  tsmall(4, 1:4) = (/0d0, 0d0, a(3,4), a(4,3)/)
  !
END SUBROUTINE libtetrabz_tsmall_c1
!
! Cut small tetrahedron C2
!
SUBROUTINE libtetrabz_tsmall_c2(e,V,tsmall)
  !
  IMPLICIT NONE
  !
  REAL(8),INTENT(IN) :: e(4)
  REAL(8),INTENT(OUT) :: V
  REAL(8),INTENT(OUT) :: tsmall(4,4)
  !
  INTEGER :: ii
  REAL(8) :: a(4,4)
  !
  DO ii = 1, 4
     a(1:4,ii) = (0d0 - e(ii)) / (e(1:4) - e(ii))
  END DO
  !
  V = a(3,4) * a(4,2)
  !
  tsmall(1, 1:4) = (/1d0,    0d0,    0d0,    0d0/)
  tsmall(2, 1:4) = (/0d0,    1d0,    0d0,    0d0/)
  tsmall(3, 1:4) = (/0d0, a(2,4),    0d0, a(4,2)/)
  tsmall(4, 1:4) = (/0d0,    0d0, a(3,4), a(4,3)/)
  !
END SUBROUTINE libtetrabz_tsmall_c2
!
! Cut small tetrahedron C3
!
SUBROUTINE libtetrabz_tsmall_c3(e,V,tsmall)
  !
  IMPLICIT NONE
  !
  REAL(8),INTENT(IN) :: e(4)
  REAL(8),INTENT(OUT) :: V
  REAL(8),INTENT(OUT) :: tsmall(4,4)
  !
  INTEGER :: ii
  REAL(8) :: a(4,4)
  !
  DO ii = 1, 4
     a(1:4,ii) = (0d0 - e(ii)) / (e(1:4) - e(ii))
  END DO
  !
  V = a(3,4) * a(2,4) * a(4,1)
  !
  tsmall(1, 1:4) = (/   1d0,    0d0,    0d0,    0d0/)
  tsmall(2, 1:4) = (/a(1,4),    0d0,    0d0, a(4,1)/)
  tsmall(3, 1:4) = (/   0d0, a(2,4),    0d0, a(4,2)/)
  tsmall(4, 1:4) = (/   0d0,    0d0, a(3,4), a(4,3)/)
  !
END SUBROUTINE libtetrabz_tsmall_c3
!
! Cut triangle A1
!
SUBROUTINE libtetrabz_triangle_a1(e,V,tsmall)
  !
  IMPLICIT NONE
  !
  REAL(8),INTENT(IN) :: e(4)
  REAL(8),INTENT(OUT) :: V
  REAL(8),INTENT(OUT) :: tsmall(3,4)
  !
  INTEGER :: ii
  REAL(8) :: a(4,4)
  !
  DO ii = 1, 4
     a(1:4,ii) = (0d0 - e(ii)) / (e(1:4) - e(ii))
  END DO
  !
  !V = 3d0 * a(2,1) * a(3,1) * a(4,1) / (0d0 - e(1))
  V = 3d0 * a(2,1) * a(3,1)           / (e(4) - e(1))
  !
  tsmall(1,1:4) = (/a(1,2), a(2,1),    0d0,    0d0/)
  tsmall(2,1:4) = (/a(1,3),    0d0, a(3,1),    0d0/) 
  tsmall(3,1:4) = (/a(1,4),    0d0,    0d0, a(4,1)/)
  !
END SUBROUTINE libtetrabz_triangle_a1
!
! Cut triangle B1
!
SUBROUTINE libtetrabz_triangle_b1(e,V,tsmall)
  !
  IMPLICIT NONE
  !
  REAL(8),INTENT(IN) :: e(4)
  REAL(8),INTENT(OUT) :: V
  REAL(8),INTENT(OUT) :: tsmall(3,4)
  !
  INTEGER :: ii
  REAL(8) :: a(4,4)
  !
  DO ii = 1, 4
     a(1:4,ii) = (0d0 - e(ii)) / (e(1:4) - e(ii))
  END DO
  !
  !V = 3d0 * a(3,1) * a(4,1) * a(2,4) / (0d0 - e(1))
  V = 3d0           * a(4,1) * a(2,4) / (e(3) - e(1))
  !
  tsmall(1,1:4) = (/a(1,3),    0d0, a(3,1),    0d0/)
  tsmall(2,1:4) = (/a(1,4),    0d0,    0d0, a(4,1)/)
  tsmall(3,1:4) = (/   0d0, a(2,4),    0d0, a(4,2)/)
  !
END SUBROUTINE libtetrabz_triangle_b1
!
! Cut triangle B2
!
SUBROUTINE libtetrabz_triangle_b2(e,V,tsmall)
  !
  IMPLICIT NONE
  !
  REAL(8),INTENT(IN) :: e(4)
  REAL(8),INTENT(OUT) :: V
  REAL(8),INTENT(OUT) :: tsmall(3,4)
  !
  INTEGER :: ii
  REAL(8) :: a(4,4)
  !
  DO ii = 1, 4
     a(1:4,ii) = (0d0 - e(ii)) / (e(1:4) - e(ii))
  END DO
  !
  !V = 3d0 * a(2,3) * a(3,1) * a(4,2) / (0d0 - e(1))
  V = 3d0 * a(2,3)           * a(4,2) / (e(3) - e(1))
  !
  tsmall(1,1:4) = (/a(1,3),    0d0, a(3,1),    0d0/)
  tsmall(2,1:4) = (/   0d0, a(2,3), a(3,2),    0d0/)
  tsmall(3,1:4) = (/   0d0, a(2,4),    0d0, a(4,2)/)
  !
END SUBROUTINE libtetrabz_triangle_b2
!
! Cut triangle C1
!
SUBROUTINE libtetrabz_triangle_c1(e,V,tsmall)
  !
  IMPLICIT NONE
  !
  REAL(8),INTENT(IN) :: e(4)
  REAL(8),INTENT(OUT) :: V
  REAL(8),INTENT(OUT) :: tsmall(3,4)
  !
  INTEGER :: ii
  REAL(8) :: a(4,4)
  !
  DO ii = 1, 4
     a(1:4,ii) = (0d0 - e(ii)) / (e(1:4) - e(ii))
  END DO
  !
  !V = 3d0 * a(1,4) * a(2,4) * a(3,4) / (e(4) - 0d0)
  V = 3d0 * a(1,4) * a(2,4)           / (e(4) - e(3))
  !
  tsmall(1,1:4) = (/a(1,4),    0d0,    0d0, a(4,1)/)
  tsmall(2,1:4) = (/   0d0, a(2,4),    0d0, a(4,2)/)
  tsmall(3,1:4) = (/   0d0,    0d0, a(3,4), a(4,3)/)
  !
END SUBROUTINE libtetrabz_triangle_c1
!
END MODULE libtetrabz_common
