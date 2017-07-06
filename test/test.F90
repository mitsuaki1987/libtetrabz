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
MODULE test_val
  !
  IMPLICIT NONE
  !
  INTEGER,SAVE :: &
  & my_proc, &
  & nb, &
  & nge(3), &
  & ngw(3), &
  & ltetra, &
  & nke, &
  & nkw
  !
  REAL(8),PARAMETER :: &
  & pi = ACOS(-1d0)
  !
  REAL(8),SAVE :: &
  & bvec(3,3), &
  & VBZ
  !
  REAL(8),ALLOCATABLE,SAVE :: &
  & eig1(:,:), &
  & eig2(:,:), &
  & mat(:)
  !
END MODULE test_val

MODULE tests
  !
  IMPLICIT NONE
  !
CONTAINS
!
!
!
SUBROUTINE test_occ
  !
#if defined(__MPI)
  USE mpi, ONLY : MPI_COMM_WORLD
#endif
  USE libtetrabz, ONLY : libtetrabz_occ
  USE test_val, ONLY : ltetra, nb, nge, ngw, nke, nkw, my_proc, &
  &                    bvec, VBZ, pi, eig1, mat
  IMPLICIT NONE
  !
  INTEGER :: ib
  REAL(8) :: val(nb), wght(nb,nkw)
  !
  eig1(1:nb,1:nke) = eig1(1:nb,1:nke) - 0.5d0
  !
#if defined(__MPI)
  CALL libtetrabz_occ(ltetra,bvec,nb,nge,eig1,ngw,wght,MPI_COMM_WORLD)
#else
  CALL libtetrabz_occ(ltetra,bvec,nb,nge,eig1,ngw,wght)
#endif
  !
  DO ib = 1, nb
     val(ib) = SUM(wght(ib,1:nkw) * mat(1:nkw))
  END DO
  IF(my_proc==0) THEN
     WRITE(*,'(a)') "# libtetrabz_occ"
     WRITE(*,'(5x,2e15.5)')         4d0 * pi / 5d0, val(1) * VBZ
     WRITE(*,'(5x,2e15.5)') pi / (5d0 * SQRT(2d0)), val(2) * VBZ
     WRITE(*,*)
  END IF
  !
  eig1(1:nb,1:nke) = eig1(1:nb,1:nke) + 0.5d0
  !
END SUBROUTINE test_occ
!
!
!
SUBROUTINE test_fermieng
  !
#if defined(__MPI)
  USE mpi, ONLY : MPI_COMM_WORLD
#endif
  USE libtetrabz, ONLY : libtetrabz_fermieng
  USE test_val, ONLY : ltetra, nb, nge, ngw, nkw, my_proc, &
  &                    bvec, VBZ, pi, eig1, mat
  IMPLICIT NONE
  !
  INTEGER :: ib
  REAL(8) :: val(nb), wght(nb,nkw), ef, nelec
  !
  nelec = (4d0 * pi / 3d0 + SQRT(2d0) * pi / 3d0) / VBZ 
#if defined(__MPI)
  CALL libtetrabz_fermieng(ltetra,bvec,nb,nge,eig1,ngw,wght,ef,nelec,MPI_COMM_WORLD)
#else
  CALL libtetrabz_fermieng(ltetra,bvec,nb,nge,eig1,ngw,wght,ef,nelec)
#endif
  !
  DO ib = 1, nb
     val(ib) = SUM(wght(ib,1:nkw) * mat(1:nkw))
  END DO
  IF(my_proc==0) THEN
     WRITE(*,'(a)') "# libtetrabz_fermieng"
     WRITE(*,'(5x,2e15.5)') 0.5d0, ef
     WRITE(*,'(5x,2e15.5)')         4d0 * pi / 5d0, val(1) * VBZ
     WRITE(*,'(5x,2e15.5)') pi / (5d0 * SQRT(2d0)), val(2) * VBZ
     WRITE(*,*)
  END IF
  !
END SUBROUTINE test_fermieng
!
!
!
SUBROUTINE test_dos
  !
#if defined(__MPI)
  USE mpi, ONLY : MPI_COMM_WORLD
#endif
  USE libtetrabz, ONLY : libtetrabz_dos
  USE test_val, ONLY : ltetra, nb, nge, ngw, nkw, my_proc, &
  &                    bvec, VBZ, pi, eig1, mat
  IMPLICIT NONE
  !
  INTEGER :: ie, ib, ne = 5
  REAL(8) :: e0(5), val(5,nb), wght(5,nb,nkw), val0(5,nb)
  !
  DO ie = 1, ne
     e0(ie) = 0.2d0 * DBLE(ie)
     val0(ie,1) = 4d0 * pi * e0(ie)**3
     IF(e0(ie) > 1d0 / SQRT(2d0)) THEN
        val0(ie,2) = SQRT(2d0) * pi * SQRT(-1d0 + 2d0 * e0(ie)**2)**3
     ELSE
        val0(ie,2) = 0d0
     END IF
  END DO
  e0(1:ne) = e0(1:ne)**2 * 0.5d0
  !
#if defined(__MPI)
  CALL libtetrabz_dos(ltetra,bvec,nb,nge,eig1,ngw,wght,ne,e0,MPI_COMM_WORLD)
#else
  CALL libtetrabz_dos(ltetra,bvec,nb,nge,eig1,ngw,wght,ne,e0)
#endif
  !
  DO ib = 1, nb
     DO ie = 1, ne
        val(ie,ib) = SUM(wght(ie,ib,1:nkw) * mat(1:nkw))
     END DO
  END DO
  !
  IF(my_proc==0) THEN
     WRITE(*,'(a)') "# libtetrabz_dos"
     DO ib = 1, nb
        DO ie = 1, ne
           WRITE(*,'(5x,2e15.5)') val0(ie,ib), val(ie,ib) * VBZ
        END DO
     END DO
     WRITE(*,*)
  END IF
  !
END SUBROUTINE test_dos
!
!
!
SUBROUTINE test_intdos
  !
#if defined(__MPI)
  USE mpi, ONLY : MPI_COMM_WORLD
#endif
  USE libtetrabz, ONLY : libtetrabz_intdos
  USE test_val, ONLY : ltetra, nb, nge, ngw, nkw, my_proc, &
  &                    bvec, VBZ, pi, eig1, mat
  IMPLICIT NONE
  !
  INTEGER :: ie, ib, ne = 5
  REAL(8) :: e0(5), val(5,nb), wght(5,nb,nkw), val0(5,nb)
  !
  DO ie = 1, ne
     e0(ie) = 0.2d0 * DBLE(ie)
     val0(ie,1) = 4d0 * pi * e0(ie)**5 / 5d0
     IF(e0(ie) > 1d0 / SQRT(2d0)) THEN
        val0(ie,2) = pi * SQRT(-1d0 + 2d0 * e0(ie)**2)**5 / (5d0 * SQRT(2d0))
     ELSE
        val0(ie,2) = 0d0
     END IF
  END DO
  e0(1:ne) = e0(1:ne)**2 * 0.5d0
  !
#if defined(__MPI)
  CALL libtetrabz_intdos(ltetra,bvec,nb,nge,eig1,ngw,wght,ne,e0,MPI_COMM_WORLD)
#else
  CALL libtetrabz_intdos(ltetra,bvec,nb,nge,eig1,ngw,wght,ne,e0)
#endif
  !
  DO ib = 1, nb
     DO ie = 1, ne
        val(ie,ib) = SUM(wght(ie,ib,1:nkw) * mat(1:nkw))
     END DO
  END DO
  !
  IF(my_proc==0) THEN
     WRITE(*,'(a)') "# libtetrabz_intdos"
     DO ib = 1, nb
        DO ie = 1, ne
           WRITE(*,'(5x,2e15.5)') val0(ie,ib), val(ie,ib) * VBZ
        END DO
     END DO
     WRITE(*,*)
  END IF
  !
END SUBROUTINE test_intdos
!
!
!
SUBROUTINE test_dblstep
  !
#if defined(__MPI)
  USE mpi, ONLY : MPI_COMM_WORLD
#endif
  USE libtetrabz, ONLY : libtetrabz_dblstep
  USE test_val, ONLY : ltetra, nb, nge, ngw, nke, nkw, my_proc, &
  &                    bvec, VBZ, pi, eig1, eig2, mat
  IMPLICIT NONE
  !
  INTEGER :: ib, jb
  REAL(8) :: val(nb,nb), wght(nb,nb,nkw)
  !
  eig1(1:nb,1:nke) = eig1(1:nb,1:nke) - 0.5d0
  eig2(1:nb,1:nke) = eig2(1:nb,1:nke) - 0.5d0
  !
#if defined(__MPI)
  CALL libtetrabz_dblstep(ltetra,bvec,nb,nge,eig1,eig2,ngw,wght,MPI_COMM_WORLD)
#else
  CALL libtetrabz_dblstep(ltetra,bvec,nb,nge,eig1,eig2,ngw,wght)
#endif
  !
  DO ib = 1, nb
     DO jb = 1, nb
        val(jb,ib) = SUM(wght(jb,ib,1:nkw) * mat(1:nkw))
     END DO
  END DO
  IF(my_proc==0) THEN
     WRITE(*,'(a)') "# libtetrabz_dblstep"
     WRITE(*,'(5x,2e15.5)') 49d0 * pi / 320d0, val(1,1) * VBZ
     WRITE(*,'(5x,2e15.5)') 0d0, val(2,1) * VBZ
     WRITE(*,'(5x,2e15.5)') pi * (512d0 * SQRT(2d0) - 319d0) / 10240d0, val(1,2) * VBZ
     WRITE(*,'(5x,2e15.5)') 0d0, val(2,2) * VBZ
     WRITE(*,*)
  END IF
  !
  eig1(1:nb,1:nke) = eig1(1:nb,1:nke) + 0.5d0
  eig2(1:nb,1:nke) = eig2(1:nb,1:nke) + 0.5d0
  !
END SUBROUTINE test_dblstep
!
!
!
SUBROUTINE test_dbldelta
  !
#if defined(__MPI)
  USE mpi, ONLY : MPI_COMM_WORLD
#endif
  USE libtetrabz, ONLY : libtetrabz_dbldelta
  USE test_val, ONLY : ltetra, nb, nge, ngw, nke, nkw, my_proc, &
  &                    bvec, VBZ, pi, eig1, eig2, mat
  IMPLICIT NONE
  !
  INTEGER :: ib, jb
  REAL(8) :: val(nb,nb), wght(nb,nb,nkw)
  !
  eig1(1:nb,1:nke) = eig1(1:nb,1:nke) - 0.5d0
  eig2(1:nb,1:nke) = eig2(1:nb,1:nke) - 0.5d0
  !
#if defined(__MPI)
  CALL libtetrabz_dbldelta(ltetra,bvec,nb,nge,eig1,eig2,ngw,wght,MPI_COMM_WORLD)
#else
  CALL libtetrabz_dbldelta(ltetra,bvec,nb,nge,eig1,eig2,ngw,wght)
#endif
  !
  DO ib = 1, nb
     DO jb = 1, nb
        val(jb,ib) = SUM(wght(jb,ib,1:nkw) * mat(1:nkw))
     END DO
  END DO
  IF(my_proc==0) THEN
     WRITE(*,'(a)') "# libtetrabz_dbldelta"
     WRITE(*,'(5x,2e15.5)') 2d0 * pi, val(1,1) * VBZ
     WRITE(*,'(5x,2e15.5)') 0d0, val(2,1) * VBZ
     WRITE(*,'(5x,2e15.5)')       pi, val(1,2) * VBZ
     WRITE(*,'(5x,2e15.5)') 0d0, val(2,2) * VBZ
     WRITE(*,*)
  END IF
  !
  eig1(1:nb,1:nke) = eig1(1:nb,1:nke) + 0.5d0
  eig2(1:nb,1:nke) = eig2(1:nb,1:nke) + 0.5d0
  !
END SUBROUTINE test_dbldelta
!
!
!
SUBROUTINE test_polstat
  !
#if defined(__MPI)
  USE mpi, ONLY : MPI_COMM_WORLD
#endif
  USE libtetrabz, ONLY : libtetrabz_polstat
  USE test_val, ONLY : ltetra, nb, nge, ngw, nke, nkw, my_proc, &
  &                    bvec, VBZ, pi, eig1, eig2, mat
  IMPLICIT NONE
  !
  INTEGER :: ib, jb
  REAL(8) :: val(nb,nb), wght(nb,nb,nkw)
  !
  eig1(1:nb,1:nke) = eig1(1:nb,1:nke) - 0.5d0
  eig2(1:nb,1:nke) = eig2(1:nb,1:nke) - 0.5d0
  !
#if defined(__MPI)
  CALL libtetrabz_polstat(ltetra,bvec,nb,nge,eig1,eig2,ngw,wght,MPI_COMM_WORLD)
#else
  CALL libtetrabz_polstat(ltetra,bvec,nb,nge,eig1,eig2,ngw,wght)
#endif
  !
  DO ib = 1, nb
     DO jb = 1, nb
        val(jb,ib) = SUM(wght(jb,ib,1:nkw) * mat(1:nkw))
     END DO
  END DO
  IF(my_proc==0) THEN
     WRITE(*,'(a)') "# libtetrabz_polstat"
     WRITE(*,'(5x,2e15.5)') pi * (68d0 + 45d0 * LOG(3d0)) / 96d0, val(1,1) * VBZ
     WRITE(*,'(5x,2e15.5)') pi * 8d0 / 5d0, val(2,1) * VBZ
     WRITE(*,'(5x,2e15.5)') pi * ( 228d0 + 22d0 * SQRT(2d0) - 96d0 * LOG(2d0) &
     &                           + 192d0 * LOG(4d0 + SQRT(2d0)) &
     &                           - 3d0 * LOG(1d0 + 2d0 * SQRT(2d0)) &
     &                           ) / 1536d0, val(1,2) * VBZ
     WRITE(*,'(5x,2e15.5)') pi * SQRT(8d0) / 5d0, val(2,2) * VBZ
     WRITE(*,*)
  END IF
  !
  eig1(1:nb,1:nke) = eig1(1:nb,1:nke) + 0.5d0
  eig2(1:nb,1:nke) = eig2(1:nb,1:nke) + 0.5d0
  !
END SUBROUTINE test_polstat
!
!
!
SUBROUTINE test_fermigr
  !
#if defined(__MPI)
  USE mpi, ONLY : MPI_COMM_WORLD
#endif
  USE libtetrabz, ONLY : libtetrabz_fermigr
  USE test_val, ONLY : ltetra, nb, nge, ngw, nke, nkw, my_proc, &
  &                    bvec, VBZ, pi, eig1, eig2, mat
  IMPLICIT NONE
  !
  INTEGER :: ne = 3, ie, ib, jb
  REAL(8) :: val(3,nb,nb), wght(3,nb,nb,nkw), val0(3,nb,nb), e0(3)
  !
  eig1(1:nb,1:nke) = eig1(1:nb,1:nke) - 0.5d0
  eig2(1:nb,1:nke) = eig2(1:nb,1:nke) - 0.5d0
  !
  e0(1) = 1d0 / 3d0
  e0(2) = 2d0 / 3d0
  e0(3) = 1d0
  val0(1,1,1) = 4d0 * pi / 9d0
  val0(2,1,1) = 1295d0 * pi / 2592d0
  val0(3,1,1) = 15d0 * pi / 32d0
  val0(1,1,2) = 5183d0 * pi / 41472d0
  val0(2,1,2) = 4559d0 * pi / 41472d0
  val0(3,1,2) = 0d0
  val0(1:3,2,1:2) = 0d0
  !
#if defined(__MPI)
  CALL libtetrabz_fermigr(ltetra,bvec,nb,nge,eig1,eig2,ngw,wght,ne,e0,MPI_COMM_WORLD)
#else
  CALL libtetrabz_fermigr(ltetra,bvec,nb,nge,eig1,eig2,ngw,wght,ne,e0)
#endif
  !
  DO ib = 1, nb
     DO jb = 1, nb
        DO ie = 1, ne
           val(ie,jb,ib) = SUM(wght(ie,jb,ib,1:nkw) * mat(1:nkw))
        END DO
     END DO
  END DO
  IF(my_proc==0) THEN
     WRITE(*,'(a)') "# libtetrabz_fermigr"
     DO ib = 1, nb
        DO jb = 1, nb
           DO ie = 1, ne
              WRITE(*,'(5x,2e15.5)') val0(ie,jb,ib), val(ie,jb,ib) * VBZ
           END DO
        END DO
     END DO
     WRITE(*,*)
  END IF
  !
  eig1(1:nb,1:nke) = eig1(1:nb,1:nke) + 0.5d0
  eig2(1:nb,1:nke) = eig2(1:nb,1:nke) + 0.5d0
  !
END SUBROUTINE test_fermigr
!
!
!
SUBROUTINE test_polcmplx
  !
#if defined(__MPI)
  USE mpi, ONLY : MPI_COMM_WORLD
#endif
  USE libtetrabz, ONLY : libtetrabz_polcmplx
  USE test_val, ONLY : ltetra, nb, nge, ngw, nke, nkw, my_proc, &
  &                    bvec, VBZ, pi, eig1, eig2, mat
  IMPLICIT NONE
  !
  INTEGER :: ne = 3, ie, ib, jb
  COMPLEX(8) :: val(3,nb,nb), wght(3,nb,nb,nkw), val0(3,nb,nb), e0(3)
  !
  eig1(1:nb,1:nke) = eig1(1:nb,1:nke) - 0.5d0
  eig2(1:nb,1:nke) = eig2(1:nb,1:nke) - 0.5d0
  !
  e0(1) = CMPLX(-2d0,    1d0, KIND(0d0))
  e0(2) = CMPLX( 0d0,    2d0, KIND(0d0))
  e0(3) = CMPLX( 1d0, -0.5d0, KIND(0d0))
  val0(1,1,1) = CMPLX(-0.838243341280338d0, - 0.734201894333234d0, KIND(0d0))
  val0(2,1,1) = CMPLX( 0.270393588876530d0, - 0.771908416949610d0, KIND(0d0))
  val0(3,1,1) = CMPLX( 0.970996830573510d0,   0.302792326476720d0, KIND(0d0))
  val0(1,1,2) = CMPLX(-0.130765724778920d0, - 0.087431218706638d0, KIND(0d0))
  val0(2,1,2) = CMPLX( 0.030121954547245d0, - 0.135354254293510d0, KIND(0d0))
  val0(3,1,2) = CMPLX( 0.178882244951203d0,   0.064232167683425d0, KIND(0d0))
  val0(1:3,2,1) = (     8d0  * pi) / (5d0 * (1d0 + 2d0 * e0(1:3)))
  val0(1:3,2,2) = (SQRT(8d0) * pi) / (5d0 * (1d0 + 4d0 * e0(1:3)))
  !
#if defined(__MPI)
  CALL libtetrabz_polcmplx(ltetra,bvec,nb,nge,eig1,eig2,ngw,wght,ne,e0,MPI_COMM_WORLD)
#else
  CALL libtetrabz_polcmplx(ltetra,bvec,nb,nge,eig1,eig2,ngw,wght,ne,e0)
#endif
  !
  DO ib = 1, nb
     DO jb = 1, nb
        DO ie = 1, ne
           val(ie,jb,ib) = SUM(wght(ie,jb,ib,1:nkw) * mat(1:nkw))
        END DO
     END DO
  END DO
  IF(my_proc==0) THEN
     WRITE(*,'(a)') "# libtetrabz_polcmplx"
     DO ib = 1, nb
        DO jb = 1, nb
           DO ie = 1, ne
              WRITE(*,'(5x,2e15.5)') DBLE( val0(ie,jb,ib)), DBLE( val(ie,jb,ib) * VBZ)
              WRITE(*,'(5x,2e15.5)') AIMAG(val0(ie,jb,ib)), AIMAG(val(ie,jb,ib) * VBZ)
           END DO
        END DO
     END DO
     WRITE(*,*)
  END IF
  !
  eig1(1:nb,1:nke) = eig1(1:nb,1:nke) + 0.5d0
  eig2(1:nb,1:nke) = eig2(1:nb,1:nke) + 0.5d0
  !
END SUBROUTINE test_polcmplx
!
END MODULE tests
!
PROGRAM test
  !
#if defined(__MPI)
  USE mpi, ONLY : MPI_COMM_WORLD
#endif
  USE tests, ONLY : test_occ, test_fermieng, test_dos, test_intdos, &
  &                 test_dblstep, test_dbldelta, test_polstat, &
  &                 test_fermigr, test_polcmplx
  USE test_val, ONLY : ltetra, nb, nge, ngw, nke, nkw, my_proc, &
  &                    bvec, VBZ, eig1, eig2, mat
  IMPLICIT NONE
  !
  INTEGER :: i1, i2, i3, ik
  REAL(8) :: kvec(3)
  !
#if defined(__MPI)
  INTEGER :: ierr
  call MPI_INIT(ierr)
  call MPI_COMM_RANK (MPI_COMM_WORLD, my_proc, ierr)
#else
  my_proc = 0
#endif
  !
  ltetra = 2
  nge(1:3) = 8
  ngw(1:3) = 8
  nke = PRODUCT(nge(1:3))
  nkw = PRODUCT(ngw(1:3))
  nb = 2
  bvec(1:3,1) = (/3d0, 0d0, 0d0/)
  bvec(1:3,2) = (/0d0, 3d0, 0d0/)
  bvec(1:3,3) = (/0d0, 0d0, 3d0/)
  VBZ = bvec(1,1) * bvec(2,2) * bvec(3,3) + bvec(1,2) * bvec(2,3) * bvec(3,1) &
  &   + bvec(1,3) * bvec(2,1) * bvec(3,2) - bvec(1,3) * bvec(2,2) * bvec(3,1) &
  &   + bvec(1,2) * bvec(2,1) * bvec(3,3) - bvec(1,1) * bvec(2,3) * bvec(3,2)
  !
  ALLOCATE(eig1(nb,nke), eig2(nb,nke), mat(nkw))
  !
  ik = 0
  DO i3 = 0, nge(3) - 1
     DO i2 = 0, nge(2) - 1
        DO i1 = 0, nge(1) - 1
           ik = ik + 1
           !
           kvec(1:3) = DBLE((/i1, i2, i3/)) / DBLE(nge(1:3))
           kvec(1:3) = kvec(1:3) - ANINT(kvec(1:3))
           kvec(1:3) = MATMUL(bvec(1:3,1:3), kvec(1:3))
           !
           eig1(1,ik) = 0.5d0 * DOT_PRODUCT(kvec(1:3), kvec(1:3))
           eig1(2,ik) = eig1(1,ik) + 0.25d0
           !
           kvec(1) = kvec(1) + 1d0
           eig2(1,ik) = 0.5d0 * DOT_PRODUCT(kvec(1:3), kvec(1:3))
           eig2(2,ik) = eig1(1,ik) + 0.5d0
           !           
        END DO
     END DO
  END DO
  !
  ik = 0
  DO i3 = 0, ngw(3) - 1
     DO i2 = 0, ngw(2) - 1
        DO i1 = 0, ngw(1) - 1
           ik = ik + 1
           !
           kvec(1:3) = DBLE((/i1, i2, i3/)) / DBLE(ngw(1:3))
           kvec(1:3) = kvec(1:3) - ANINT(kvec(1:3))
           kvec(1:3) = MATMUL(bvec(1:3,1:3), kvec(1:3))
           !
           mat(ik) = DOT_PRODUCT(kvec(1:3), kvec(1:3))
        END DO
     END DO
  END DO
  !
  IF(my_proc == 0) WRITE(*,'(a,10x,a,10x,a)') "#", "Ideal", "Result"
  !
  CALL test_occ()
  !
  CALL test_fermieng()
  !
  CALL test_dos()
  !
  CALL test_intdos()
  !
  CALL test_dblstep()
  !
  CALL test_dbldelta()
  !
  CALL test_polstat()
  !
  CALL test_fermigr()
  !
  CALL test_polcmplx()
  !
#if defined(__MPI)
  call MPI_FINALIZE(ierr)
#endif
  !
END PROGRAM test
