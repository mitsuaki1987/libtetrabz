module global
  !
  implicit none
  !
  integer,save :: &
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
  use global, only : nb, ng, nk, nm, g2, omg, eig1, eig2, nelec, bvec
  implicit none
  !
  integer :: fi = 10, fstb
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
end subroutine read_elph
!
! Dos calclation
!
subroutine calc_dos()
  !
  use libtetra, only : libtetrabz_fermieng, libtetrabz_dos, libtetrabz_doubledelta
  use global, only : ng, nb, nk, nelec, bvec, eig1, eig2, g2, lam, nm, dos
  !
  implicit none
  !
  real(8) :: ef
  real(8),allocatable :: wlam(:)
  !
  allocate(wlam(nb * nk))
  !
  ef = 0d0
  call libtetrabz_fermieng(2, ng, nb, ef, nelec, bvec, eig1, wlam)
  eig1(1:nb,1:nk) = eig1(1:nb,1:nk) - ef
  eig2(1:nb,1:nk) = eig2(1:nb,1:nk) - ef
  write(*,*) "  Fermi energy [Ry]", -minval(eig1(1:nb,1:nk))
  !
  call libtetrabz_dos(2, ng, nb, 1, bvec, (/0d0/), eig1, wlam)
  !
  dos = sum(wlam(1:nb * nk))
  write(*,'(a,10e18.8)') "  DOS[/Ryd/cell] = ", dos
  !
  deallocate(wlam)
  !
  allocate(wlam(nb * nb * nk), lam(nm))
  call libtetrabz_doubledelta(2, ng, nb, bvec, eig1, eig2, wlam)
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
  use global, only : read_elph, calc_dos, integrate, &
  &                  nm, omg, lam, dos
  use omp_lib
  !
  implicit none
  !
  integer :: im
  real(8) :: t1, t2
  !
  t1 = OMP_GET_WTIME()
  !
  !$OMP PARALLEL
  if(OMP_GET_THREAD_NUM() == 0) then
     write(*,*) " # of threads : ", OMP_GET_NUM_THREADS()
  end if
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
  write(*,*) "  mode #, frequence[Ryd], lambda : "
  do im = 1, nm
     if(omg(im) <= 0d0 ) then
        lam(im) = 0d0
     else
        lam(im) = lam(im) * 2d0 / (omg(im) * dos)
     end if
     write(*,*) im, omg(im), lam(im)
  end do
  !
  t2 = OMP_GET_WTIME()
  !
  write(*,*) ""
  write(*,*) "  Total time : ", t2 - t1, " sec"
  write(*,*) ""
  write(*,*) "#####  Done  #####"
  write(*,*) ""
  !
end program main
