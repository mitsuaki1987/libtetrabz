module global
  !
  implicit none
  !
  integer,save ::  &
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
  use global, only : nb, ng, nwfc, bvec, e0, dos, ne, eig, wfc
  !
  implicit none
  !
  integer :: fi = 10
  !
  open(fi, file = "proj.dat")
  read(fi,*) ng(1:3)
  read(fi,*) nb
  read(fi,*) nwfc
  read(fi,*) ne
  read(fi,*) bvec(1:3,1:3)
  !
  allocate(e0(ne), eig(nb, ng(1), ng(2), ng(3)), wfc(nwfc, nb, ng(1), ng(2), ng(3)))
  !
  read(fi,*) e0(1:ne)
  read(fi,*) eig(1:nb,1:ng(1),1:ng(2),1:ng(3))
  read(fi,*) wfc(1:nwfc,1:nb,1:ng(1),1:ng(2),1:ng(3))
  close(fi)
  !
  allocate(dos(nwfc,ne))
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
end subroutine read_dat
!
! Dos calclation
!
subroutine calc_dos()
  !
  use libtetra_mod, only : libtetra_dos
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
  call libtetra("dos", 2, ng, nb, ne, bvec, e0, eig, wdos)
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
  use libtetra_mod, only : libtetra_occ, libtetra_occ_int
  use global, only : ng, nb, bvec, eig, nwfc, wfc
  implicit none
  !
  integer :: ie, i1, i2, i3, ib, iwfc, iexp, n
  real(8) :: nelec, ef, occ(nwfc)
  real(8),allocatable :: wocc(:,:,:,:)
  !
  allocate(wocc(nb,ng(1),ng(2),ng(3))
  call libtetra_occ("occ", 2, ng, nb, bvec, eig, wocc)
  !
  occ(1:nwfc) = 0d0
  !
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP & SHARED(ng,nb,nwfc,occ,wfc,wocc) &
  !$OMP & PRIVATE(i1,i2,i3,ib)
  !
  !$OMP DO REDUCTION(+: occ)
  do i3 = 1, ng(3)
     do i2 = 1, ng(2)
        do i1 = 1, ng(1)
           do ib = 1, nb
              occ(1:nwfc) = occ(1:nwfc) + wfc(1:nwfc, ib, i1, i2, i3) &
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
  write(*,*)
  write(*,*) "# WFC,  OCC "
  do iwfc = 1, nwfc
     write(*,*) iwfc, occ(iwfc)
  end do
  write(*,*)
  !
  do iexp = 1, 5
     !
     n = 2**iexp
     !
     allocate(wocc(nb,ng(1) /n,ng(2)/n,ng(3)/n))
     call libtetra_occ_int("occ", 2, ng, ng / n, nb, bvec, eig, wocc)
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
     write(*,*)
     write(*,*) "# WFC,  OCC "
     do iwfc = 1, nwfc
        write(*,*) iwfc, occ(iwfc)
     end do
     write(*,*)
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
  use global, only : read_dat, calc_dos, write_dos, ng, calc_occ
  implicit none
  !
  integer :: i1, i2, i3, ik
  !
  ! Read Quantum-ESPRESSO output
  !
  write(*,*) ""
  write(*,*) "#####  Read Quantum-ESPRESSO output  #####"
  write(*,*) ""
  call read_dat()
  !
  write(*,*) ""
  write(*,*) "#####  Clc. Dos  #####"
  write(*,*) ""
  !
  call calc_dos()
  !
  call calc_occ()
  !
  call write_dos()
  !
  write(*,*) 
  write(*,*) "#####  Done  #####"
  write(*,*)
  !
end program pdos
