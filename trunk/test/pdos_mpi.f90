module global
  !
  implicit none
  !
  integer,save ::  &
  & petot,         &
  & my_rank,       &
  & ne,            &
  & nb,            & !
  & ng(3),         &
  & nwfc
  !
  real(8),save :: &
  & bvec(3,3)
  !
  real(8),save,allocatable :: &
  & e0(:),          & ! (ne)
  & dos(:,:),       & ! (nwfc,ne)
  & eig(:,:,:,:),   & ! (nb,ng(1),ng(2),ng(3)) 
  & wfc(:,:,:,:,:)   ! (nwfc + 1, nb,ng(1),ng(2),ng(3)) 
  !
  interface
     !
     subroutine read_dat()
     end subroutine read_dat
     !
     subroutine calc_dos()
     end subroutine calc_dos
     !
     subroutine calc_occ()
     end subroutine calc_occ
     !
     subroutine write_dos()
     end subroutine write_dos
     !
  end interface
  !
end module global
!
! Read Quantum-ESPRESSO output
!
subroutine read_dat()
  !
  use mpi
  use global, only : nb, ng, nwfc, bvec, e0, dos, ne, eig, wfc, my_rank
  !
  implicit none
  !
  integer :: fi = 10, ierr
  !
  if(my_rank == 0) then
     !
     open(fi, file = "proj.dat")
     read(fi,*) ng(1:3)
     read(fi,*) nb
     read(fi,*) nwfc
     read(fi,*) ne
     read(fi,*) bvec(1:3,1:3)
     !
     allocate(e0(ne), eig(nb, ng(1), ng(2), ng(3)), &
     &        wfc(nwfc, nb, ng(1), ng(2), ng(3)))
     !
     read(fi,*) e0(1:ne)
     read(fi,*) eig(1:nb,1:ng(1),1:ng(2),1:ng(3))
     read(fi,*) wfc(1:nwfc,1:nb,1:ng(1),1:ng(2),1:ng(3))
     close(fi)
     !
     write(*,*) "  # of bands : ", nb
     write(*,*) "  # of atomic wfc : ", nwfc - 1
     write(*,*) "  Reciprocal lattice vector[2pi/a] : "
     write(*,'(3f10.5)') bvec(1:3,1:3)
     write(*,*) "  MONKHORST_PACK_GRID : ", ng(1:3)
     write(*,*) "  Emin[Ry] : ", minval(eig(1:nb,1:ng(1),1:ng(2),1:ng(3)))
     write(*,*) "  Emax[Ry] : ", maxval(eig(1:nb,1:ng(1),1:ng(2),1:ng(3)))
     write(*,*) " Energy step [Ryd] : ", e0(2) - e0(1)
     write(*,*) " # of energy for DOS : ", ne
     !
  end if
  !
  call MPI_BCAST(nb,     1, MPI_INTEGER,          0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(nwfc,   1, MPI_INTEGER,          0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(ng,     3, MPI_INTEGER,          0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(ne,     1, MPI_INTEGER,          0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(bvec,   9, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  !
  if(my_rank /= 0) allocate(e0(ne), eig(nb, ng(1), ng(2), ng(3)), &
  &                          wfc(nwfc, nb, ng(1), ng(2), ng(3)))
  !
  call MPI_BCAST(e0,                            ne, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(eig,        nb * product(ng(1:3)), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(wfc, nwfc * nb * product(ng(1:3)), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  !
  allocate(dos(nwfc,ne))
  !
end subroutine read_dat
!
! Dos calclation
!
subroutine calc_dos()
  !
  use mpi
  use libtetrabz_mpi, only : libtetrabz_mpi_dos
  use global, only : ng, nb, bvec, eig, dos, ne, nwfc, wfc, e0
  implicit none
  !
  integer :: ie, i1, i2, i3, ib, iwfc
  real(8) :: nelec, ef
  real(8),allocatable :: wdos(:,:,:,:,:)
  !
  nelec = 0d0
  ef = 0d0
  !
  allocate(wdos(ne,nb,ng(1),ng(2),ng(3)))
  call libtetrabz_mpi_dos(2,MPI_COMM_WORLD, bvec,nb,ng,eig,ng,wdos,ne,e0)
  !
  dos(1:nwfc,1:ne) = 0d0
  !
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP & SHARED(ng,nb,ne,nwfc,dos,wfc,wdos) &
  !$OMP & PRIVATE(i1,i2,i3,ib,ie)
  !
  !$OMP DO REDUCTION(+: dos)
  do i3 = 1, ng(3)
     do i2 = 1, ng(2)
        do i1 = 1, ng(1)
           do ib = 1, nb
              do ie = 1, ne
                 dos(1:nwfc,ie) = dos(1:nwfc,ie) + wfc(1:nwfc, ib, i1, i2, i3) &
                 &                             *  wdos(  ie, ib, i1, i2, i3)
              end do
           end do
        end do
     end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL
  !
  deallocate(wdos)
  !
end subroutine calc_dos
!
!
!
subroutine calc_occ()
  !
  use mpi
  use libtetrabz_mpi, only : libtetrabz_mpi_occ
  use global, only : ng, nb, bvec, eig, nwfc, wfc, my_rank
  implicit none
  !
  integer :: ie, i1, i2, i3, ib, iwfc, iexp, n
  real(8) :: nelec, ef, occ(nwfc)
  real(8),allocatable :: wocc(:,:,:,:)
  !
  do iexp = 1, 6
     !
     n = 2**(iexp - 1)
     !
     allocate(wocc(nb,ng(1) /n,ng(2)/n,ng(3)/n))
     call libtetrabz_mpi_occ(2,MPI_COMM_WORLD, bvec,nb,ng,eig,ng / n,wocc)
     !
     occ(1:nwfc) = 0d0
     !
     !$OMP PARALLEL DEFAULT(NONE) &
     !$OMP & SHARED(ng,nb,nwfc,occ,wfc,wocc,n) &
     !$OMP & PRIVATE(i1,i2,i3,ib)
     !
     !$OMP DO REDUCTION(+: occ)
     do i3 = 1, ng(3) / n
        do i2 = 1, ng(2) / n
           do i1 = 1, ng(1) / n
              do ib = 1, nb
                 occ(1:nwfc) = occ(1:nwfc) + wfc(1:nwfc, ib, (i1-1)*n+1, (i2-1)*n+1, (i3-1)*n+1) &
                 &                         *  wocc(      ib, i1, i2, i3)
              end do
           end do
        end do
     end do
     !$OMP END DO
     !$OMP END PARALLEL
     !
     deallocate(wocc)
     !
     if(my_rank == 0) then
        write(*,*)
        write(*,*) "# WFC,  OCC "
        do iwfc = 1, nwfc
           write(*,*) iwfc, occ(iwfc)
        end do
        write(*,*)
     end if
     !
  end do
  !
end subroutine calc_occ
!
! Write dos file
!
subroutine write_dos()
  !
  use global, only : ne, nwfc, e0, dos
  !
  implicit none
  !
  integer :: fo = 20, ie
  !
  open(fo, file = "pdos.dat")
  !
  do ie = 1, ne
     write(fo,'(1000e18.8)') e0(ie), dos(1:nwfc,ie)
  end do
  !
  close(fo)
  !
end subroutine write_dos
!
! Main routine
!
program pdos
  !
  !$ use omp_lib
  use mpi
  use global, only : read_dat, calc_dos, write_dos, ng, calc_occ, my_rank, petot
  implicit none
  !
  integer :: ierr
  !
  call MPI_INIT(ierr)
  call MPI_COMM_SIZE (MPI_COMM_WORLD, PETOT, ierr)
  call MPI_COMM_RANK (MPI_COMM_WORLD, my_rank, ierr)
  !
  if(my_rank == 0) write(*,*) " # of PEs : ", petot
  !$OMP PARALLEL
  !$ if(OMP_GET_THREAD_NUM() == 0) then
  !$    if(my_rank == 0) write(*,*) " # of threads : ", OMP_GET_NUM_THREADS()
  !$ end if
  !$OMP END PARALLEL
  !
  ! Read Quantum-ESPRESSO output
  !
  if(my_rank == 0) write(*,*) ""
  if(my_rank == 0) write(*,*) "#####  Read Quantum-ESPRESSO output  #####"
  if(my_rank == 0) write(*,*) ""
  call read_dat()
  !
  if(my_rank == 0) write(*,*) ""
  if(my_rank == 0) write(*,*) "#####  Culc. Dos  #####"
  if(my_rank == 0) write(*,*) ""
  !
  call calc_dos()
  !
  call calc_occ()
  !
  if(my_rank == 0) call write_dos()
  !
  if(my_rank == 0) write(*,*) 
  if(my_rank == 0) write(*,*) "#####  Done  #####"
  if(my_rank == 0) write(*,*)
  !
  call MPI_FINALIZE(ierr)
  !
end program pdos
