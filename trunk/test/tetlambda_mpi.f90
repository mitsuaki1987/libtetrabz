module global
  !
  implicit none
  !
  integer,save :: &
  & petot,      &
  & my_rank,    &
  & nb,         & ! # of bands
  & ng(3),      & ! k-point (dense) grid
  & nk,         & ! # of k. ng(1) * ng(2) * ng(3)
  & nm
  !
  real(8),save :: & 
  & nelec,      & ! # of electron
  & bvec(3,3),  & ! Reciprocal lattice vector
  & dos           ! Dos(E_F)
  !
  real(8),allocatable,save :: &
  & eig1(:,:),   & ! KS enegrgy [Ry] (nb,nk)
  & eig2(:,:),   & ! KS enegrgy [Ry] (nb,nk)
  & g2(:,:),     & ! El-Ph matrix [Ry^2] (nm,nb * nb * nk)
  & omg(:),      & ! Phonon frequency [Ry] (nm)
  & lam(:)         ! El-Ph coupling constant (nm)
  !
  interface
     !
     subroutine read_elph()
     end subroutine read_elph
     !
     subroutine km_fermieng()
     end subroutine km_fermieng
     !
     function sumnelec(ef) result(wsum)
       real(8),intent(in) :: ef
       real(8) :: wsum
     end function sumnelec
     !
     subroutine calc_dos()
     end subroutine calc_dos
     !
     subroutine integrate()
     end subroutine integrate
     !
  end interface
  !
end module global
!
! Read elph file
!
subroutine read_elph()
  !
  use mpi
  use global, only : nb, ng, nk, nm, g2, omg, eig1, eig2, nelec, bvec, my_rank
  implicit none
  !
  integer :: fi = 10, ierr
  !
  if(my_rank == 0) then
     !
     open(fi, file = "elph.dat" )
     !
     read(fi,*) ng(1:3)
     nk = product(ng(1:3))
     !
     read(fi,*) nb
     !
     read(fi,*) nelec
     !
     read(fi,*) nm
     allocate(g2(nm,nb * nb * nk), omg(nm), eig1(nb,nk), eig2(nb,nk))
     !
     read(fi,*) bvec(1:3,1:3)
     !
     read(fi,*) omg(1:nm)
     !
     ! Read |g|^2
     !
     read(fi,*) g2(1:nm,1:nb * nb * nk)
     !
     ! Read Kohn-Sham energies
     !
     read(fi,*) eig1(1:nb,1:nk)
     read(fi,*) eig2(1:nb,1:nk)
     !
     close(fi)
     !
     write(*,*) "  Reciprocal lattice vector[2pi/a] : "
     write(*,'(3f10.5)') bvec(1:3,1:3)
     write(*,*) "  # of electrons : ", nelec
     write(*,*) ""
     write(*,*) "  k point grid : ", ng(1:3)
     write(*,*) "  # of bands : ", nb
     write(*,*) "  # of modes : ", nm
     !
  end if
  !
  call MPI_BCAST(nb,     1, MPI_INTEGER,          0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(ng,     3, MPI_INTEGER,          0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(nm,     1, MPI_INTEGER,          0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(nelec,  1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(bvec,   9, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  !
  nk = product(ng(1:3))
  if(my_rank /= 0) allocate(g2(nm,nb * nb * nk), omg(nm), eig1(nb,nk), eig2(nb,nk))
  !
  call MPI_BCAST(omg,                 nm, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(eig1,           nb * nk, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(eig2,           nb * nk, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(g2,   nm * nb * nb * nk, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  !
end subroutine read_elph
!
! Dos calclation
!
subroutine calc_dos()
  !
  use mpi
  use libtetrabz_mpi, only : libtetrabz_mpi_fermieng, libtetrabz_mpi_dos, libtetrabz_mpi_doubledelta
  use global, only : ng, nb, nk, nelec, bvec, eig1, eig2, g2, lam, nm, dos, my_rank
  !
  implicit none
  !
  real(8) :: ef
  real(8),allocatable :: wlam(:)
  !
  allocate(wlam(nb * nk))
  !
  ef = 0d0
  call libtetrabz_mpi_fermieng(2,MPI_COMM_WORLD, bvec,nb,ng,eig1,ng,wlam,ef,nelec)
  eig1(1:nb,1:nk) = eig1(1:nb,1:nk) - ef
  eig2(1:nb,1:nk) = eig2(1:nb,1:nk) - ef
  if(my_rank == 0) write(*,*) "  Fermi energy [Ry]", -minval(eig1(1:nb,1:nk))
  !
  call libtetrabz_mpi_dos(2,MPI_COMM_WORLD, bvec,nb,ng,eig1,ng,wlam,1,(/0d0/))
  !
  dos = sum(wlam(1:nb * nk))
  if(my_rank == 0) write(*,'(a,10e18.8)') "  DOS[/Ryd/cell] = ", dos
  !
  deallocate(wlam)
  !
  allocate(wlam(nb * nb * nk), lam(nm))
  call libtetrabz_mpi_doubledelta(2,MPI_COMM_WORLD, bvec,nb,ng,eig1,eig2,ng,wlam)
  !
  lam(1:nm) = matmul(g2(1:nm, 1:nb * nb * nk), wlam(1:nb * nb * nk))
  !
  deallocate(wlam)
  !
end subroutine calc_dos
!
! Main module
!
program main
  !
  !$ use omp_lib
  use mpi
  use global, only : read_elph, calc_dos, integrate, &
  &                  nm, omg, lam, dos, my_rank, petot
  use omp_lib
  !
  implicit none
  !
  integer :: im, ierr
  real(8) :: t1, t2
  !
  call MPI_INIT(ierr)
  call MPI_COMM_SIZE (MPI_COMM_WORLD, PETOT, ierr)
  call MPI_COMM_RANK (MPI_COMM_WORLD, my_rank, ierr)
  !
  t1 = OMP_GET_WTIME()
  !
  if(my_rank == 0) write(*,*) " # of PEs : ", petot
  !$OMP PARALLEL
  !$ if(OMP_GET_THREAD_NUM() == 0) then
  !$    if(my_rank == 0) write(*,*) " # of threads : ", OMP_GET_NUM_THREADS()
  !$ end if
  !$OMP END PARALLEL
  !
  ! Read elph.dat
  !
  call read_elph()
  !
  ! Calc. DOS
  !
  call calc_dos()
  !
  ! calc. lambda
  !
  if(my_rank == 0) write(*,*) "  mode #, frequence[Ryd], lambda : "
  do im = 1, nm
     if(omg(im) <= 0d0 ) then
        lam(im) = 0d0
     else
        lam(im) = lam(im) * 2d0 / (omg(im) * dos)
     end if
     if(my_rank == 0) write(*,*) im, omg(im), lam(im)
  end do
  !
  t2 = OMP_GET_WTIME()
  !
  if(my_rank == 0) write(*,*) ""
  if(my_rank == 0) write(*,*) "  Total time : ", t2 - t1, " sec"
  if(my_rank == 0) write(*,*) ""
  if(my_rank == 0) write(*,*) "#####  Done  #####"
  if(my_rank == 0) write(*,*) ""
  !
  call MPI_FINALIZE(ierr)
  !
end program main
