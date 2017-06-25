PROGRAM test
  !
#if defined(__MPI)
  USE mpi, ONLY : 
#endif
  USE libtetrabz, ONLY : libtetrabz_fermieng
  IMPLICIT NONE
  !
  INTEGER :: ltetra, nb, nge(3), ngw(3), i1, i2, i3, ik, &
  &          nke, nkw, my_proc
  REAL(8) :: bvec(3,3), ef, kvec(3), VBZ, pi, nelec, val
  REAL(8),ALLOCATABLE :: eig(:,:), wght(:,:), mat(:,:)
  !
#if defined(__MPI)
  INTEGER :: ierr
  call MPI_INIT(ierr)
  call MPI_COMM_RANK (MPI_COMM_WORLD, my_proc, ierr)
#elif
  my_proc = 0
#endif
  !
  ltetra = 2
  nge(1:3) = 16
  ngw(1:3) = 16
  nke = PRODUCT(nge(1:3))
  nkw = PRODUCT(ngw(1:3))
  pi = ACOS(-1d0)
  nb = 1
  bvec(1:3,1) = (/3d0, 0d0, 0d0/)
  bvec(1:3,2) = (/0d0, 3d0, 0d0/)
  bvec(1:3,3) = (/0d0, 0d0, 3d0/)
  VBZ = bvec(1,1) * bvec(2,2) * bvec(3,3) + bvec(1,2) * bvec(2,3) * bvec(3,1) &
  &   + bvec(1,3) * bvec(2,1) * bvec(3,2) - bvec(1,3) * bvec(2,2) * bvec(3,1) &
  &   + bvec(1,2) * bvec(2,1) * bvec(3,3) - bvec(1,1) * bvec(2,3) * bvec(3,2)
  !
  ALLOCATE(eig(nb,nke), wght(nb,nkw), mat(nb,nkw))
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
           eig(1,ik) = 0.5d0 * DOT_PRODUCT(kvec(1:3), kvec(1:3))
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
  nelec = 4d0 * pi / 3d0 / VBZ
  CALL libtetrabz_fermieng(ltetra,bvec,nb,nge,eig,ngw,wght,ef,nelec)
  eig(1:nb,1:nke) = eig(1:nb,1:nke) - ef
  !
  IF(my_rank == 0) THEN
     WRITE(*,*) "Fermi energy"
     WRITE(*,*) "   Ideal : ", 0.5d0
     WRITE(*,*) "  Result : ", ef
  END IF
  !
  val = SUM(wght(1:nb,1:nkw) * mat(1:nb,1:nkw))  
  IF(my_rank == 0) THEN
     WRITE(*,*) "Integration"
     WRITE(*,*) "   Ideal : ", 4d0 * pi / 5d0
     WRITE(*,*) "  Result : ", val * VBZ
  END IF
  !
  write(*,*) "debug", SUM(wght(1:nb,1:nkw))
  !
#if defined{__MPI)
  call MPI_FINALIZE(ierr)
#endif
  !
END PROGRAM test
