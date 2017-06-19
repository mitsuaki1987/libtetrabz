MODULE libtetrabz_routines
  !
  IMPLICIT NONE
  !
CONTAINS
!
! define shortest diagonal line & define type of tetragonal
!
SUBROUTINE libtetrabz_initialize(bvec)
  !
  USE libterabz_val, ONLY : ltetra, ng, wlsm, ivvec, nt, nk, ng
  IMPLICIT NONE
  !
  REAL(8),INTENT(IN) :: bvec(3,3)
  !
  INTEGER :: itype, i1, i2, i3, it, ii, divvec(4,4), ivvec0(4)
  REAL(8) :: l(4), bvec2(3,3), bvec3(3,4)
  !
  nk = PRODUCT(ng(1:3))
  nt = nk * 6
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
END SUBROUTINE libtetrabz_initialize
  for (i = 0; i < n; ++i) swap[i] = i;

  for (i = 0; i < n - 1; ++i) {
    for (j = i + 1; j < n; ++j) {
      if (key[swap[j]] < key[swap[i]]) {
        /*
         Swap
        */
        k = swap[j];
        swap[j] = swap[i];
        swap[i] = k;
      }/*if (sortee[j][0] < sortee[i][0])*/
    }/*for (j = i + 1; j < n; ++j)*/
  }/*for (i = 0; i < n - 1; ++i)*/
!
! Simple sort
!
SUBROUTINE libtetrabz_sort(n,key,sind)
  !
  IMPLICIT NONE
  !
  integer,INTENT(IN) :: n
  REAL(8),INTENT(inout) :: key(n)
  INTEGER,INTENT(OUT) :: sind(n)
  !
  INTEGER :: i, i0, sind0
  REAL(8) :: key0
  !
  DO i = 1, n
     sind(i) = i
  END DO
  !
  DO i = 1, n - 1
     key0 = MINVAL(key(i+1:n))
     i0   = MINLOC(key(i+1:n),1) + i
     IF(key(i) > key0) THEN
        key(i0) = key(i)
        key(i) = key0
        !
        sind0 = sind(i0)
        sind(i0) = sind(i)
        sind(i) = sind0
     END IF
  END DO
  !
END SUBROUTINE libtetrabz_sort
!
! Interpolate integration weight
!
SUBROUTINE libtetrabz_interpol_weight(nb,ngc,ngd,wc,wd)
  !
  USE libterabz_val, ONLY : nk0, indx3
  IMPLICIT NONE
  !
  integer,INTENT(IN) :: nb, ngc(3), ngd(3)
  REAL(8),INTENT(IN) :: wd(nb,nk0)
  REAL(8),INTENT(OUT) :: wc(nb,product(ngc(1:3)))
  !
  INTEGER :: i1, i2, i3, ik, nkc, nkd
  REAL(8) :: kv(3, PRODUCT(ngd(1:3)))
  !
  nkc = PRODUCT(ngc(1:3))
  nkd = PRODUCT(ngd(1:3))
  !
  ik = 0
  DO i3 = 1, ngd(3)
     DO i2 = 1, ngd(2)
        DO i1 = 1, ngd(1)
           ik = ik + 1
           kv(1:3,ik) = DBLE((/i1, i2, i3/) - 1) / DBLE(ngd(1:3))
        END DO
     END DO
  END DO
  !
  wc(1:nb,1:nkc) = 0d0
  !
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP & SHARED(nk0,nkc,nkd,nb,ngc,kv,wc,wd,indx3) &
  !$OMP PRIVATE(ik)
  !
  !$OMP DO REDUCTION(+: wc)
  DO ik = 1, nk0
     CALL libtetrabz_interpol_weight2(nkc, nb, ngc, kv(1:3,indx3(ik)), wd(1:nb,ik), wc)
  END DO
  !$OMP END DO
  !$OMP END PARALLEL
  !
END SUBROUTINE libtetrabz_interpol_weight
!
! first or third order interpolation of weights
!
SUBROUTINE libtetrabz_interpol_weight2(nk,nb,ng,ko,wi,wo)
  !
  USE libterabz_val, ONLY : ltetra, ivvec
  IMPLICIT NONE
  !
  INTEGER,INTENT(IN)  :: nk, nb, ng(3)
  REAL(8),INTENT(IN)  :: ko(3)
  REAL(8),INTENT(IN) :: wi(nb)
  REAL(8),INTENT(INOUT) :: wo(nb,nk)
  !
  INTEGER :: ikv(3), ikv1(3), ik(20), ii, it, it0, ierr
  REAL(8) :: rot(3,3), res(3), prod(3), u, x, y, z, thr = 1d-10
  !
  rot(1:3,1) = (/  2d0, - 1d0,   0d0/)
  rot(1:3,2) = (/- 1d0,   2d0, - 1d0/)
  rot(1:3,3) = (/  0d0, - 1d0,   1d0/)
  !
  ! Search nearest neighbor grid points.
  !
  res(1:3) = ko(1:3) * DBLE(ng(1:3))
  ikv(1:3) = FLOOR(res(1:3))
  res(1:3) = res(1:3) - DBLE(ikv(1:3))
  !
  DO it = 1, 6
     !
     DO ii = 1, 3
        prod(ii) = DOT_PRODUCT(DBLE(ivvec(1:3,1 + ii,it) - ivvec(1:3,1,it)), &
        &                                  res(1:3) - DBLE(ivvec(1:3,1,it))  )
     END DO
     !
     prod(1:3) = MATMUL(rot(1:3,1:3), prod(1:3))
     !
     IF(MINVAL(prod(1:3)) > - thr .AND. SUM(prod(1:3)) < 1d0 + thr) THEN
        it0 = it
        GOTO 10
     END IF
     !
  END DO
  !
  STOP "interpol"
  !
10 CONTINUE
  !
  x = prod(1)
  y = prod(2)
  z = prod(3)
  u = 1d0 - x - y - z
  !
  DO ii = 1, 20
     !
     ikv1(1:3) = ikv(1:3) + ivvec(1:3,ii,it0)
     ikv1(1:3) = MODULO(ikv1(1:3), ng(1:3))
     ik(ii) = 1 + ikv1(1) + ikv1(2) * ng(1) + ikv1(3) * ng(1) * ng(2)
     !
  END DO
  !
  IF(ltetra == 0 .OR. ltetra == 1) THEN
     !
     wo(1:nb,ik(1)) = wo(1:nb,ik(1)) + wi(1:nb) * u
     wo(1:nb,ik(2)) = wo(1:nb,ik(2)) + wi(1:nb) * x
     wo(1:nb,ik(3)) = wo(1:nb,ik(3)) + wi(1:nb) * y
     wo(1:nb,ik(4)) = wo(1:nb,ik(4)) + wi(1:nb) * z
     !
  ELSE IF(ltetra == 2) THEN
     !
     wo(1:nb,ik( 1)) = wo(1:nb,ik( 1)) + wi(1:nb) * 0.5d0 * u * (2d0 + u * (1d0 - u) + 2d0 * y * (x + z))
     wo(1:nb,ik( 2)) = wo(1:nb,ik( 2)) + wi(1:nb) * 0.5d0 * x * (2d0 + x * (1d0 - x) + 2d0 * z * (u + y))
     wo(1:nb,ik( 3)) = wo(1:nb,ik( 3)) + wi(1:nb) * 0.5d0 * y * (2d0 + y * (1d0 - y) + 2d0 * u * (x + z))
     wo(1:nb,ik( 4)) = wo(1:nb,ik( 4)) + wi(1:nb) * 0.5d0 * z * (2d0 + z * (1d0 - z) + 2d0 * x * (u + y))
     wo(1:nb,ik( 5)) = wo(1:nb,ik( 5)) + wi(1:nb) * x * u * (2d0 * y - u - 1d0) / 6d0
     wo(1:nb,ik( 6)) = wo(1:nb,ik( 6)) + wi(1:nb) * x * y * (2d0 * z - x - 1d0) / 6d0
     wo(1:nb,ik( 7)) = wo(1:nb,ik( 7)) + wi(1:nb) * y * z * (2d0 * u - y - 1d0) / 6d0
     wo(1:nb,ik( 8)) = wo(1:nb,ik( 8)) + wi(1:nb) * z * u * (2d0 * x - z - 1d0) / 6d0
     wo(1:nb,ik( 9)) = wo(1:nb,ik( 9)) + wi(1:nb) * y * u * (2d0 * y + u - 3d0) / 6d0
     wo(1:nb,ik(10)) = wo(1:nb,ik(10)) + wi(1:nb) * x * z * (2d0 * z + x - 3d0) / 6d0
     wo(1:nb,ik(11)) = wo(1:nb,ik(11)) + wi(1:nb) * y * u * (2d0 * u + y - 3d0) / 6d0
     wo(1:nb,ik(12)) = wo(1:nb,ik(12)) + wi(1:nb) * x * z * (2d0 * x + z - 3d0) / 6d0
     wo(1:nb,ik(13)) = wo(1:nb,ik(13)) + wi(1:nb) * z * u * (2d0 * y - u - 1d0) / 6d0
     wo(1:nb,ik(14)) = wo(1:nb,ik(14)) + wi(1:nb) * x * u * (2d0 * z - x - 1d0) / 6d0
     wo(1:nb,ik(15)) = wo(1:nb,ik(15)) + wi(1:nb) * x * y * (2d0 * u - y - 1d0) / 6d0
     wo(1:nb,ik(16)) = wo(1:nb,ik(16)) + wi(1:nb) * y * z * (2d0 * x - z - 1d0) / 6d0
     wo(1:nb,ik(17)) = wo(1:nb,ik(17)) + wi(1:nb) * (- x * z * u)
     wo(1:nb,ik(18)) = wo(1:nb,ik(18)) + wi(1:nb) * (- x * y * u)
     wo(1:nb,ik(19)) = wo(1:nb,ik(19)) + wi(1:nb) * (- x * y * z)
     wo(1:nb,ik(20)) = wo(1:nb,ik(20)) + wi(1:nb) * (- y * z * u)
     !
  ELSE
     !
     STOP "interpol2"
     ! 
  END IF
  !
END SUBROUTINE libtetrabz_interpol_weight2
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
END MODULE libtetrabz_routines
