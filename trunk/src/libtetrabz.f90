module libtetrabz_mod
  !
  use libtetrabz_common, only : nb, nk
  !
  interface
     !
     subroutine libtetrabz(job,ltetra_0,ngd,ngc,nb_0,ne_0,ef_0,nelec_0,bvec_0,e0,eig1,eig2,wght)
       character(*) :: job
       integer,intent(in) :: ltetra_0, ngd(3), ngc(3), nb_0, ne_0
       real(8),intent(in) :: bvec_0(3,3), nelec_0
       real(8),intent(in) :: e0(ne_0)
       real(8),intent(in) :: eig1(nb_0,product(ngd(1:3))), eig2(nb_0,product(ngd(1:3)))
       real(8),intent(inout) :: ef_0
       real(8),intent(out) :: wght(1:*)
     end subroutine libtetrabz
     !
     subroutine libtetrabz_fermieng(eig,occ)
       import nb, nk
       real(8),intent(in) :: eig(nb,nk)
       real(8),intent(out) :: occ(nb,nk)
     end subroutine libtetrabz_fermieng
     !
     subroutine libtetrabz_kgrid()
     end subroutine libtetrabz_kgrid
     !
  end interface
  !
end module libtetrabz_mod
!
! Driver routine of libtetrabz
!
subroutine libtetrabz(job,ltetra_0,ngd,ngc,nb_0,ne_0,ef_0,nelec_0,bvec_0,e0,eig1,eig2,wght)
  !
  use libtetrabz_common, only : ltetra, ng, nb, ne, ef, nk, indx1, indx2, indx3, bvec, nelec, &
  &                           libtetrabz_initialize, libtetrabz_interpol_weight, &
  &                           libtetrabz_occ, libtetrabz_dos, libtetrabz_doubledelta, &
  &                           libtetrabz_occstep, libtetrabz_polstat, libtetrabz_fermigr, libtetrabz_polimg
  !
  use libtetrabz_mod, only : libtetrabz_fermieng, libtetrabz_kgrid
  !
  implicit none
  !
  character(*) :: job
  integer,intent(in) :: ltetra_0, ngd(3), ngc(3), nb_0, ne_0
  real(8),intent(in) :: bvec_0(3,3), nelec_0
  real(8),intent(in) :: e0(ne_0)
  real(8),intent(in) :: eig1(nb_0,product(ngd(1:3))), eig2(nb_0,product(ngd(1:3)))
  real(8),intent(inout) :: ef_0
  real(8),intent(out) :: wght(1:*)
  !
  integer :: nn
  logical :: lintp
  real(8),allocatable :: wghtd(:,:)
  !
  ltetra = ltetra_0
  nb = nb_0
  ne = ne_0
  ng(1:3) = ngd(1:3)
  nk = product(ngd(1:3))
  ef = ef_0
  nelec = nelec_0
  bvec(1:3,1:3) = bvec_0(1:3,1:3)
  !
  lintp = any(ngd(1:3) /= ngc(1:3))
  !
  call libtetrabz_initialize()
  call libtetrabz_kgrid()
  !
  if(trim(job) == "fermieng") then
     !
     nn = nb
     !
     if(lintp) then
        allocate(wghtd(nn, nk))
        call libtetrabz_fermieng(eig1,wghtd)
     else
        call libtetrabz_fermieng(eig1,wght)
     end if
     !
     ef_0 = ef
     !
  else if(trim(job) == "occ") then
     !
     nn = nb
     !
     if(lintp) then
        allocate(wghtd(nn, nk))
        call libtetrabz_occ(eig1,wghtd)
     else
        call libtetrabz_occ(eig1,wght)
     end if
     !
  else if(trim(job) == "dos") then
     !
     nn = nb * ne
     !
     if(lintp) then
        allocate(wghtd(nn, nk))
        call libtetrabz_dos(eig1,e0,wghtd)
     else
        call libtetrabz_dos(eig1,e0,wght)
     end if
     !
  else if(trim(job) == "doubledelta") then
     !
     nn = nb * nb
     !
     if(lintp) then
        allocate(wghtd(nn, nk))
        call libtetrabz_doubledelta(eig1,eig2,wghtd)
     else
        call libtetrabz_doubledelta(eig1,eig2,wght)
     end if
     !
  else if(trim(job) == "occstep") then
     !
     nn = nb * nb
     !
     if(lintp) then
        allocate(wghtd(nn, nk))
        call libtetrabz_occstep(eig1,eig2,wghtd)
     else
        call libtetrabz_occstep(eig1,eig2,wght)
     end if
     !
  else if(trim(job) == "polstat") then
     !
     nn = nb * nb
     !
     if(lintp) then
        allocate(wghtd(nn, nk))
        call libtetrabz_polstat(eig1,eig2,wghtd)
     else
        call libtetrabz_polstat(eig1,eig2,wght)
     end if
     !
  else if(trim(job) == "fermigr") then
     !
     nn = ne * nb * nb
     !
     if(lintp) then
        allocate(wghtd(nn, nk))
        call libtetrabz_fermigr(eig1,eig2,e0,wghtd)
     else
        call libtetrabz_fermigr(eig1,eig2,e0,wght)
     end if
     !
  else if(trim(job) == "polimg") then
     !
     nn = 2 * ne * nb * nb
     !
     if(lintp) then
        allocate(wghtd(nn, nk))
        call libtetrabz_polimg(eig1,eig2,e0,wghtd)
     else
        call libtetrabz_polimg(eig1,eig2,e0,wght)
     end if
     !
  else
     !
     stop "Invalid job"
     !
  end if
  !
  if(lintp) then
     call libtetrabz_interpol_weight(nn,ngc,ngd,wght,wghtd)
     deallocate(wghtd)
  end if
  !
  deallocate(indx1, indx2, indx3)
  !
end subroutine libtetrabz
!
! Calculate Fermi energy
!
subroutine libtetrabz_fermieng(eig,occ)
  !
  use libtetrabz_common, only : nelec, ef, nb, nk, ef, libtetrabz_occ
  implicit none
  !
  real(8),intent(in) :: eig(nb,nk)
  real(8),intent(out) :: occ(nb,nk)
  !
  integer :: iter, maxiter = 300
  real(8) :: elw, eup, sumkmid, eps= 1.0d-10
  !
  ! find bounds for the Fermi energy.
  !
  elw = minval(eig(1:nb,1:nk))
  eup = maxval(eig(1:nb,1:nk))
  !
  !      Bisection method
  !
  do iter = 1, maxiter
     !
     ef = (eup + elw) / 2.d0
     !
     ! Calc. # of electrons 
     !
     call libtetrabz_occ(eig,occ)
     !
     sumkmid = sum(occ(1:nb,1:nk)) * 2d0
     !
     ! convergence check
     !
     if(abs(sumkmid - nelec) < eps) then
        exit
     elseif(sumkmid < nelec) then
        elw = ef
     else
        eup = ef
     endif
     !
  enddo ! iter
  !
  if(iter >= maxiter) stop "libtetrabz_omp_fermieng"
  !
end subroutine libtetrabz_fermieng
!
! Occupation
!
subroutine libtetrabz_occ(ltetra_0,ng_0,nb_0,bvec_0,eig1,wght)
  !
  use libtetrabz_common, only : ltetra, ng, nb, nk, indx1, indx2, indx3, bvec, &
  &                           libtetrabz_initialize, libtetrabz_occ1
  !
  use libtetrabz_mod, only : libtetrabz_kgrid
  !
  implicit none
  !
  integer,intent(in) :: ltetra_0, ng_0(3), nb_0
  real(8),intent(in) :: bvec_0(3,3)
  real(8),intent(in) :: eig1(nb_0,product(ng_0(1:3)))
  real(8),intent(out) :: wght(nb_0,product(ng_0(1:3)))
  !
  ltetra = ltetra_0
  nb = nb_0
  ng(1:3) = ng_0(1:3)
  nk = product(ng_0(1:3))
  bvec(1:3,1:3) = bvec_0(1:3,1:3)
  !
  call libtetrabz_initialize()
  call libtetrabz_kgrid()
  !
  call libtetrabz_occ1(eig1,wght)
  !
  deallocate(indx1, indx2, indx3)
  !
end subroutine libtetrabz_occ
!
! Occupation (interpolated)
!
subroutine libtetrabz_occ_int(ltetra_0,ngd,ngc,nb_0,bvec_0,eig1,wght)
  !
  use libtetrabz_common, only : ltetra, ng, nb, nk, indx1, indx2, indx3, bvec, &
  &                           libtetrabz_initialize, libtetrabz_interpol_weight, &
  &                           libtetrabz_occ1
  !
  use libtetrabz_mod, only : libtetrabz_kgrid
  !
  implicit none
  !
  integer,intent(in) :: ltetra_0, ngd(3), ngc(3), nb_0
  real(8),intent(in) :: bvec_0(3,3)
  real(8),intent(in) :: eig1(nb_0,product(ngd(1:3)))
  real(8),intent(out) :: wght(nb_0,product(ngd(1:3)))
  !
  real(8),allocatable :: wghtd(:,:)
  !
  ltetra = ltetra_0
  nb = nb_0
  ng(1:3) = ngd(1:3)
  nk = product(ngd(1:3))
  bvec(1:3,1:3) = bvec_0(1:3,1:3)
  !
  call libtetrabz_initialize()
  call libtetrabz_kgrid()
  !
  allocate(wghtd(nn, nk))
  call libtetrabz_occ(eig1,wghtd)
  !
  call libtetrabz_interpol_weight(nn,ngc,ngd,wght,wghtd)
  !
  deallocate(wghtd, indx1, indx2, indx3)
  !
end subroutine libtetrabz_occ
!
! Initialize grid
!
subroutine libtetrabz_kgrid()
  !
  use libtetrabz_common, only : nk, nk0, indx1, indx2, indx3, ng, ivvec, fst, lst
  !
  implicit none
  !
  integer :: it, i1, i2, i3, ii, ikv(3), nt, ik
  !
  allocate(indx1(20, 6 * nk), indx2(20, 6 * nk), indx3(20 * 6 * nk))
  !
  nt = 0
  do i3 = 1, ng(3)
     do i2  = 1, ng(2)
        do i1 = 1, ng(1)
           !
           do it = 1, 6
              !
              nt = nt + 1
              !
              do ii = 1, 20
                 !
                 ikv(1:3) = (/i1, i2, i3/) + ivvec(1:3,ii,it) - 1
                 ikv(1:3) = modulo(ikv(1:3), ng(1:3))
                 !
                 indx1(ii,nt) = 1 + ikv(1) + ng(1) * ikv(2) + ng(1) * ng(2) * ikv(3)
                 !
              end do
              !
           end do
           !
        end do
     end do
  end do
  !
  indx2(1:20,1:6 * nk) = indx1(1:20,1:6 * nk)
  indx3(1:20 * 6 * nk) = 0
  !
  nk0 = nk
  fst = 1
  lst = 6 * nk
  !
  do ik = 1, nk
     indx3(ik) = ik
  end do
  !
end subroutine libtetrabz_kgrid
