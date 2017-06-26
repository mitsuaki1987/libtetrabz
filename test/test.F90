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
  & mat(:,:)
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
  REAL(8) :: val, wght(nb,nkw)
  !
  IF(my_proc==0) WRITE(*,'(a)') " libtetrabz_occ"
  !
  eig1(1:nb,1:nke) = eig1(1:nb,1:nke) - 0.5d0
  !
#if defined(__MPI)
  CALL libtetrabz_occ(ltetra,bvec,nb,nge,eig1,ngw,wght,MPI_COMM_WORLD)
#else
  CALL libtetrabz_occ(ltetra,bvec,nb,nge,eig1,ngw,wght)
#endif
  !
  val = SUM(wght(1:nb,1:nkw) * mat(1:nb,1:nkw))  
  IF(my_proc==0) THEN
     WRITE(*,'(a,e15.5)') "     Ideal : ", 4d0 * pi / 5d0
     WRITE(*,'(a,e15.5)') "    Result : ", val * VBZ
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
  REAL(8) :: val, wght(nb,nkw), ef, nelec
  !
  IF(my_proc == 0) THEN
     WRITE(*,'(a)') " libtetrabz_fermieng"
  END IF
  !
  nelec = 4d0 * pi / 3d0 / VBZ
#if defined(__MPI)
  CALL libtetrabz_fermieng(ltetra,bvec,nb,nge,eig1,ngw,wght,ef,nelec,MPI_COMM_WORLD)
#else
  CALL libtetrabz_fermieng(ltetra,bvec,nb,nge,eig1,ngw,wght,ef,nelec)
#endif
  !
  val = SUM(wght(1:nb,1:nkw) * mat(1:nb,1:nkw))  
  IF(my_proc==0) THEN
     WRITE(*,'(a,2e15.5)') "     Ideal : ", 0.5d0, 4d0 * pi / 5d0
     WRITE(*,'(a,2e15.5)') "    Result : ",    ef, val * VBZ
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
  INTEGER :: ie, ne = 5
  REAL(8) :: e0(5), val(5), wght(5,nb,nkw), val0(5)
  !
  IF(my_proc==0) WRITE(*,'(a)') " libtetrabz_dos"
  !
  DO ie = 1, ne
     e0(ie) = 0.2d0 * DBLE(ie)
     val0(ie) = 4d0 * pi * e0(ie)**3
  END DO
  e0(1:ne) = e0(1:ne)**2 * 0.5d0
  !
#if defined(__MPI)
  CALL libtetrabz_dos(ltetra,bvec,nb,nge,eig1,ngw,wght,ne,e0,MPI_COMM_WORLD)
#else
  CALL libtetrabz_dos(ltetra,bvec,nb,nge,eig1,ngw,wght,ne,e0)
#endif
  !
  DO ie = 1, ne
     val(ie) = SUM(wght(ie,1:nb,1:nkw) * mat(1:nb,1:nkw))
  END DO
  !
  IF(my_proc==0) THEN
     WRITE(*,'(a,5e15.5)') "     Ideal : ", val0(1:ne)
     WRITE(*,'(a,5e15.5)') "    Result : ", val(1:ne) * VBZ
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
  INTEGER :: ie, ne = 5
  REAL(8) :: e0(5), val(5), wght(5,nb,nkw), val0(5)
  !
  IF(my_proc==0) WRITE(*,'(a)') " libtetrabz_intdos"
  !
  DO ie = 1, ne
     e0(ie) = 0.2d0 * DBLE(ie)
     val0(ie) = 4d0 * pi * e0(ie)**5 / 5d0
  END DO
  e0(1:ne) = e0(1:ne)**2 * 0.5d0
  !
#if defined(__MPI)
  CALL libtetrabz_intdos(ltetra,bvec,nb,nge,eig1,ngw,wght,ne,e0,MPI_COMM_WORLD)
#else
  CALL libtetrabz_intdos(ltetra,bvec,nb,nge,eig1,ngw,wght,ne,e0)
#endif
  !
  DO ie = 1, ne
     val(ie) = SUM(wght(ie,1:nb,1:nkw) * mat(1:nb,1:nkw))
  END DO
  !
  IF(my_proc==0) THEN
     WRITE(*,'(a,5e15.5)') "     Ideal : ", val0(1:ne)
     WRITE(*,'(a,5e15.5)') "    Result : ", val(1:ne) * VBZ
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
  REAL(8) :: val, wght(nb,nb,nkw)
  !
  IF(my_proc==0) WRITE(*,'(a)') " libtetrabz_dblstep"
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
  val = SUM(wght(1,1,1:nkw) * mat(1,1:nkw))  
  IF(my_proc==0) THEN
     WRITE(*,'(a,e15.5)') "     Ideal : ", 49d0 * pi / 320d0
     WRITE(*,'(a,e15.5)') "    Result : ", val * VBZ
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
  REAL(8) :: val, wght(nb,nb,nkw)
  !
  IF(my_proc==0) WRITE(*,'(a)') " libtetrabz_dbldelta"
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
  val = SUM(wght(1,1,1:nkw) * mat(1,1:nkw))  
  IF(my_proc==0) THEN
     WRITE(*,'(a,e15.5)') "     Ideal : ", 2d0 * pi
     WRITE(*,'(a,e15.5)') "    Result : ", val * VBZ
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
  REAL(8) :: val, wght(nb,nb,nkw)
  !
  IF(my_proc==0) WRITE(*,'(a)') " libtetrabz_polstat"
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
  val = SUM(wght(1,1,1:nkw) * mat(1,1:nkw))  
  IF(my_proc==0) THEN
     WRITE(*,'(a,e15.5)') "     Ideal : ", pi * (68d0 + 45d0 * LOG(3d0)) / 96d0
     WRITE(*,'(a,e15.5)') "    Result : ", val * VBZ
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
  INTEGER :: ne = 3, ie
  REAL(8) :: val(3), wght(3,nb,nb,nkw), val0(3), e0(3)
  !
  IF(my_proc==0) WRITE(*,'(a)') " libtetrabz_fermigr"
  !
  eig1(1:nb,1:nke) = eig1(1:nb,1:nke) - 0.5d0
  eig2(1:nb,1:nke) = eig2(1:nb,1:nke) - 0.5d0
  !
  e0(1) = 1d0 / 3d0
  e0(2) = 2d0 / 3d0
  e0(3) = 1d0
  val0(1) = 4d0 * pi / 9d0
  val0(2) = 1295d0 * pi / 2592d0
  val0(3) = 15d0 * pi / 32d0
  !
#if defined(__MPI)
  CALL libtetrabz_fermigr(ltetra,bvec,nb,nge,eig1,eig2,ngw,wght,ne,e0,MPI_COMM_WORLD)
#else
  CALL libtetrabz_fermigr(ltetra,bvec,nb,nge,eig1,eig2,ngw,wght,ne,e0)
#endif
  !
  DO ie = 1, ne
     val(ie) = SUM(wght(ie,1,1,1:nkw) * mat(1,1:nkw))
  END DO
  IF(my_proc==0) THEN
     WRITE(*,'(a,3e15.5)') "     Ideal : ", val0(1:ne)
     WRITE(*,'(a,3e15.5)') "    Result : ", val(1:ne) * VBZ
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
  &                    bvec, VBZ, eig1, eig2, mat
  IMPLICIT NONE
  !
  INTEGER :: ne = 3, ie
  COMPLEX(8) :: val(3), wght(3,nb,nb,nkw), val0(3), e0(3)
  !
  IF(my_proc==0) WRITE(*,'(a)') " libtetrabz_polcmplx"
  !
  eig1(1:nb,1:nke) = eig1(1:nb,1:nke) - 0.5d0
  eig2(1:nb,1:nke) = eig2(1:nb,1:nke) - 0.5d0
  !
  e0(1) = CMPLX(-2d0,    1d0, KIND(0d0))
  e0(2) = CMPLX( 0d0,    2d0, KIND(0d0))
  e0(3) = CMPLX( 1d0, -0.5d0, KIND(0d0))
  val0(1) = CMPLX(-0.838243341280338, - 0.734201894333234, KIND(0d0))
  val0(2) = CMPLX( 0.270393588876530, - 0.771908416949610, KIND(0d0))
  val0(3) = CMPLX( 0.970996830573510,   0.302792326476720, KIND(0d0))
  !
#if defined(__MPI)
  CALL libtetrabz_polcmplx(ltetra,bvec,nb,nge,eig1,eig2,ngw,wght,ne,e0,MPI_COMM_WORLD)
#else
  CALL libtetrabz_polcmplx(ltetra,bvec,nb,nge,eig1,eig2,ngw,wght,ne,e0)
#endif
  !
  DO ie = 1, ne
     val(ie) = SUM(wght(ie,1,1,1:nkw) * mat(1,1:nkw))
  END DO
  IF(my_proc==0) THEN
     WRITE(*,'(a,6e15.5)') "     Ideal : ", val0(1:ne)
     WRITE(*,'(a,6e15.5)') "    Result : ", val(1:ne) * VBZ
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
  nge(1:3) = 16
  ngw(1:3) = 16
  nke = PRODUCT(nge(1:3))
  nkw = PRODUCT(ngw(1:3))
  nb = 1
  bvec(1:3,1) = (/3d0, 0d0, 0d0/)
  bvec(1:3,2) = (/0d0, 3d0, 0d0/)
  bvec(1:3,3) = (/0d0, 0d0, 3d0/)
  VBZ = bvec(1,1) * bvec(2,2) * bvec(3,3) + bvec(1,2) * bvec(2,3) * bvec(3,1) &
  &   + bvec(1,3) * bvec(2,1) * bvec(3,2) - bvec(1,3) * bvec(2,2) * bvec(3,1) &
  &   + bvec(1,2) * bvec(2,1) * bvec(3,3) - bvec(1,1) * bvec(2,3) * bvec(3,2)
  !
  ALLOCATE(eig1(nb,nke), eig2(nb,nke), mat(nb,nkw))
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
           !
           kvec(1) = kvec(1) + 1d0
           eig2(1,ik) = 0.5d0 * DOT_PRODUCT(kvec(1:3), kvec(1:3))
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
           mat(1,ik) = DOT_PRODUCT(kvec(1:3), kvec(1:3))
        END DO
     END DO
  END DO
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
