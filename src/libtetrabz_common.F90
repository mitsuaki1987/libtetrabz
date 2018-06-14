!
! Copyright (c) 2014 Mitsuaki Kawamura
!
! Permission is hereby granted, free of charge, to any person obtaining a
! copy of this software and associated documentation files (the
! "Software"), to deal in the Software without restriction, including
! without limitation the rights to use, copy, modify, merge, publish,
! distribute, sublicense, and/or sell copies of the Software, and to
! permit persons to whom the Software is furnished to do so, subject to
! the following conditions:
! 
! The above copyright notice and this permission notice shall be included
! in all copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
! OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
! MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
! IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
! CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
! TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
! SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
!
MODULE libtetrabz_common
  !
  IMPLICIT NONE
  !
  PRIVATE
  PUBLIC libtetrabz_initialize, libtetrabz_sort, libtetrabz_interpol_indx, &
  &      libtetrabz_tsmall_a1, libtetrabz_tsmall_b1, libtetrabz_tsmall_b2, libtetrabz_tsmall_b3, &
  &      libtetrabz_tsmall_c1, libtetrabz_tsmall_c2, libtetrabz_tsmall_c3, &
  &      libtetrabz_triangle_a1, libtetrabz_triangle_b1, &
  &      libtetrabz_triangle_b2, libtetrabz_triangle_c1, &
  &      libtetrabz_mpisum_d, libtetrabz_mpisum_dv, libtetrabz_mpisum_zv
  !
CONTAINS
!
! define shortest diagonal line & define type of tetragonal
!
SUBROUTINE libtetrabz_initialize(ltetra,nge,ngw,bvec,linterpol,wlsm,nk_local,nt_local,nkBZ,ik_global,ik_local,kvec,comm)
  !
  IMPLICIT NONE
  !
  INTEGER,INTENT(IN) :: ltetra, nge(3), ngw(3)
  REAL(8),INTENT(IN) :: bvec(3,3)
  LOGICAL,INTENT(OUT) :: linterpol
  REAL(8),INTENT(OUT) :: wlsm(4,20)
  INTEGER,INTENT(OUT) :: nk_local, nt_local, nkBZ
  INTEGER,INTENT(OUT),ALLOCATABLE :: ik_global(:,:), ik_local(:,:)
  REAL(8),INTENT(OUT),ALLOCATABLE :: kvec(:,:)
  INTEGER,INTENT(IN),OPTIONAL :: comm
  !
  INTEGER :: itype, i1, i2, i3, it, divvec(4,4), ivvec0(4), ivvec(3,20,6)
  REAL(8) :: l(4), bvec2(3,3), bvec3(3,4)
  !
  nkBZ = PRODUCT(nge(1:3))
  linterpol = .NOT. ALL(nge(1:3) == ngw(1:3))
  !
  DO i1 = 1, 3
     bvec2(1:3,i1) = bvec(1:3,i1) / DBLE(nge(i1))
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
  IF (PRESENT(comm)) THEN
     CALL libtetrabz_kgrid(linterpol,ivvec,nge,nkBZ,nk_local,nt_local,ik_global,ik_local,kvec,comm)
  ELSE
     CALL libtetrabz_kgrid(linterpol,ivvec,nge,nkBZ,nk_local,nt_local,ik_global,ik_local,kvec)
  END IF
  !
END SUBROUTINE libtetrabz_initialize
!
! Initialize grid
!
SUBROUTINE libtetrabz_kgrid(linterpol,ivvec,ng,nkBZ,nk_local,nt_local,ik_global,ik_local,kvec,comm)
  !
  IMPLICIT NONE
  !
  LOGICAL,INTENT(INOUT) :: linterpol
  INTEGER,INTENT(IN) :: ivvec(3,20,6), ng(3), nkBZ
  INTEGER,INTENT(OUT) :: nk_local, nt_local
  INTEGER,INTENT(OUT),ALLOCATABLE :: ik_global(:,:), ik_local(:,:)
  REAL(8),INTENT(OUT),ALLOCATABLE :: kvec(:,:)
  INTEGER,INTENT(IN),OPTIONAL :: comm
  !
  INTEGER :: it, i1, i2, i3, ii, ikv(3), nt, ik, nt_front, loc2glob(nkBZ), glob2loc(nkBZ)
  !
  IF(PRESENT(comm)) THEN
     CALL libtetrabz_divideMPI(comm,6 * nkBZ,nt_front,nt_local)
     linterpol = linterpol .OR. (6*nkBZ /= nt_local)
  ELSE
     nt_front = 0
     nt_local = 6 * nkBZ
  END IF
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
  IF(.NOT. linterpol) THEN
     nk_local = nkBZ
     ik_local(1:20,1:nt_local) = ik_global(1:20,1:nt_local)     
     RETURN
  END IF
  !
  glob2loc(1:nkBZ) = 0
  nk_local = 0
  DO nt = 1, nt_local
     DO ii = 1, 20
        !
        IF(glob2loc(ik_global(ii,nt)) /= 0) THEN
           ik_local(ii,nt) = glob2loc(ik_global(ii,nt))
        ELSE
           !
           nk_local = nk_local + 1
           loc2glob(nk_local) = ik_global(ii,nt)
           glob2loc(ik_global(ii,nt)) = nk_local
           ik_local(ii,nt) = nk_local
           !
        END IF
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
SUBROUTINE libtetrabz_divideMPI(comm,nt,nt_front,nt_local)
  !
#if defined(__MPI)
  USE mpi, ONLY : MPI_COMM_SIZE, MPI_COMM_RANK
#endif
  IMPLICIT NONE
  !
  INTEGER,INTENT(IN) :: comm, nt
  INTEGER,INTENT(OUT) :: nt_front, nt_local
  !
  INTEGER :: petot = 1, my_rank = 0
#if defined(__MPI)
  INTEGER :: ierr
  CALL MPI_COMM_SIZE(comm, petot, ierr)
  CALL MPI_COMM_RANK(comm, my_rank, ierr)
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
END SUBROUTINE libtetrabz_divideMPI
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
  INTEGER,INTENT(out) :: kintp(8)
  REAL(8),INTENT(out) :: wintp(8)
  !
  INTEGER :: ii
  REAL(8) :: xv(3)
  integer :: ikv0(3), ikv1(3), i1, i2, i3
  !
  xv(1:3) = kvec(1:3) * DBLE(ng(1:3))
  ikv0(1:3) = FLOOR(xv(1:3))
  xv(1:3) = xv(1:3) - DBLE(ikv0(1:3))
  !
  ii = 0
  DO i1 = 0, 1
     DO i2 = 0, 1
        DO i3 = 0, 1
           !
           ii = ii + 1
           !
           ikv1(1:3) = ikv0(1:3) + (/i1, i2, i3/)
           ikv1(1:3) = MODULO(ikv1(1:3), ng(1:3))
           kintp(ii) = 1 + ikv1(1) + ng(1) * ikv1(2) + ng(1) * ng(2) * ikv1(3)
           !
           wintp(ii) = xv(1)**i1 * (1.0d0 - xv(1))**(1 - i1) &
           &         * xv(2)**i2 * (1.0d0 - xv(2))**(1 - i2) &
           &         * xv(3)**i3 * (1.0d0 - xv(3))**(1 - i3) 
           !
        END DO
     END DO
  END DO
  !
END SUBROUTINE libtetrabz_interpol_indx
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
  V = a(2,1) * a(3,1) * a(4,1)
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
! MPI_Allreduce for double scaler
!
SUBROUTINE libtetrabz_mpisum_d(comm,scaler)
  !
#if defined(__MPI)
  USE mpi, ONLY : MPI_DOUBLE_PRECISION, MPI_IN_PLACE, MPI_SUM
#endif
  IMPLICIT NONE
  !
  INTEGER :: comm
  REAL(8) :: scaler
  !
#if defined(__MPI)
  INTEGER :: ierr
  !
  CALL MPI_allREDUCE(MPI_IN_PLACE, scaler, 1, &
  &                  MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)
#endif
  !
END SUBROUTINE libtetrabz_mpisum_d
!
! MPI_Allreduce for double vector
!
SUBROUTINE libtetrabz_mpisum_dv(comm,ndim,vector)
  !
#if defined(__MPI)
  USE mpi, ONLY : MPI_DOUBLE_PRECISION, MPI_IN_PLACE, MPI_SUM
#endif
  IMPLICIT NONE
  !
  INTEGER :: comm, ndim
  REAL(8) :: vector(ndim)
  !
#if defined(__MPI)
  INTEGER :: ierr
  !
  CALL MPI_allREDUCE(MPI_IN_PLACE, vector, ndim, &
  &                  MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)
#endif
  !
END SUBROUTINE libtetrabz_mpisum_dv
!
! MPI_Allreduce for double complex vector
!
SUBROUTINE libtetrabz_mpisum_zv(comm,ndim,vector)
  !
#if defined(__MPI)
  USE mpi, ONLY : MPI_DOUBLE_COMPLEX, MPI_IN_PLACE, MPI_SUM
#endif
  IMPLICIT NONE
  !
  INTEGER :: comm, ndim
  COMPLEX(8) :: vector(ndim)
  !
#if defined(__MPI)
  INTEGER :: ierr
  !
  CALL MPI_allREDUCE(MPI_IN_PLACE, vector, ndim, &
  &                  MPI_DOUBLE_COMPLEX, MPI_SUM, comm, ierr)
#endif
  !
END SUBROUTINE libtetrabz_mpisum_zv
!
END MODULE libtetrabz_common
