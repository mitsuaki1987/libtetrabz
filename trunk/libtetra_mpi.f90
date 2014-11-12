module libtetra_mpi_mod
  !
  use libtetra_common, only : nb, nk
  !
  integer,save :: comm
  !
  interface
     !
     subroutine libtetra_mpi(job,ltetra_0,ngd,ngc,nb_0,ne_0,ef_0,nelec_0,bvec_0,e0,eig1,eig2,wght)
       character(*) :: job
       integer,intent(in) :: ltetra_0, ngd(3), ngc(3), nb_0, ne_0
       real(8),intent(in) :: bvec_0(3,3), nelec_0
       real(8),intent(in) :: e0(ne_0)
       real(8),intent(in) :: eig1(nb_0,product(ngd(1:3))), eig2(nb_0,product(ngd(1:3)))
       real(8),intent(inout) :: ef_0
       real(8),intent(out) :: wght(1:*)
     end subroutine libtetra_mpi
     !
     subroutine libtetra_mpi_fermieng(eig,occ)
       import nb, nk
       real(8),intent(in) :: eig(nb,nk)
       real(8),intent(out) :: occ(nb,nk)
     end subroutine libtetra_mpi_fermieng
     !
     subroutine libtetra_mpi_kgrid()
     end subroutine libtetra_mpi_kgrid
     !
     subroutine libtetra_fst_and_lst()
     end subroutine libtetra_fst_and_lst
     !
  end interface
  !
end module libtetra_mpi_mod
!
! Driver routine of libtetra
!
subroutine libtetra_mpi(job,ltetra_0,comm_0,ngd,ngc,nb_0,ne_0,ef_0,nelec_0,bvec_0,e0,eig1,eig2,wght)
  !
  use mpi, only : MPI_ALLREDUCE, MPI_IN_PLACE, MPI_DOUBLE_PRECISION, MPI_SUM
  use libtetra_common, only : ltetra, ng, nb, ne, ef, nk, indx1, indx2, indx3, bvec, &
  &                           libtetra_initialize, libtetra_interpol_weight, &
  &                           libtetra_occ, libtetra_dos, libtetra_doubledelta, &
  &                           libtetra_occstep, libtetra_polstat, libtetra_fermigr, libtetra_polimg
  !
  use libtetra_mpi_mod, only : comm, libtetra_mpi_fermieng, libtetra_mpi_kgrid
  !
  character(*) :: job
  integer,intent(in) :: ltetra_0, comm_0, ngd(3), ngc(3), nb_0, ne_0
  real(8),intent(in) :: bvec_0(3,3), nelec_0
  real(8),intent(in) :: e0(ne_0)
  real(8),intent(in) :: eig1(nb_0,product(ngd(1:3))), eig2(nb_0,product(ngd(1:3)))
  real(8),intent(inout) :: ef_0
  real(8),intent(out) :: wght(1:*)
  !
  integer :: nn, ierr
  logical :: lintp
  real(8),allocatable :: wghtd(:,:)
  !
  ltetra = ltetra_0
  comm = comm_0
  nb = nb_0
  ne = ne_0
  ng(1:3) = ngd(1:3)
  nk = product(ngd(1:3))
  ef = ef_0
  nelec = nelec_0
  bvec(1:3,1:3) = bvec_0(1:3,1:3)
  !
  lintp = all(ngd(1:3) == ngc(1:3))
  !
  call libtetra_initialize()
  call libtetra_mpi_kgrid()
  !
  if(trim(job) == "fermieng") then
     !
     nn = nb
     !
     if(lintp) then
        allocate(wghtd(nn, nk))
        call libtetra_fermieng(eig1,wghtd)
     else
        call libtetra_fermieng(eig1,wght)
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
        call libtetra_occ(eig1,wghtd)
     else
        call libtetra_occ(eig1,wght)
     end if
     !
  else if(trim(job) == "dos") then
     !
     nn = nb * ne
     !
     if(lintp) then
        allocate(wghtd(nn, nk))
        call libtetra_dos(eig1,e0,wghtd)
     else
        call libtetra_dos(eig1,e0,wght)
     end if
     !
  else if(trim(job) == "doubledelta") then
     !
     nn = nb * nb
     !
     if(lintp) then
        allocate(wghtd(nn, nk))
        call libtetra_doubledelta(eig1,eig2,wghtd)
     else
        call libtetra_doubledelta(eig1,eig2,wght)
     end if
     !
  else if(trim(job) == "occstep") then
     !
     nn = nb * nb
     !
     if(lintp) then
        allocate(wghtd(nn, nk))
        call libtetra_occstep(eig1,eig2,wghtd)
     else
        call libtetra_occstep(eig1,eig2,wght)
     end if
     !
  else if(trim(job) == "polstat") then
     !
     nn = nb * nb
     !
     if(lintp) then
        allocate(wghtd(nn, nk))
        call libtetra_polstat(eig1,eig2,wghtd)
     else
        call libtetra_polstat(eig1,eig2,wght)
     end if
     !
  else if(trim(job) == "fermigr") then
     !
     nn = ne * nb * nb
     !
     if(lintp) then
        allocate(wghtd(nn, nk))
        call libtetra_fermigr(eig1,eig2,e0,wghtd)
     else
        call libtetra_fermigr(eig1,eig2,e0,wght)
     end if
     !
  else if(trim(job) == "polimg") then
     !
     nn = 2 * ne * nb * nb
     !
     if(lintp) then
        allocate(wghtd(nn, nk))
        call libtetra_polimg(eig1,eig2,e0,wghtd)
     else
        call libtetra_polimg(eig1,eig2,e0,wght)
     end if
     !
  else
     !
     stop "Invalid job"
     !
  end if
  !
  if(lintp) then
     call libtetra_interpol_weight(nn,ngc,ngd,wght,wghtd)
     deallocate(wghtd)
  end if
  !
  deallocate(indx1, indx2, indx3)
  !
  call MPI_allREDUCE(MPI_IN_PLACE, wght, nn * product(ngc(1:3)), &
  &                  MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)
  !
end subroutine libtetra_mpi
!
! Calculate Fermi energy
!
subroutine libtetra_mpi_fermieng(eig,occ)
  !
  use mpi, only : MPI_ALLREDUCE, MPI_IN_PLACE, MPI_DOUBLE_PRECISION, MPI_SUM
  use libtetra_common, only : nelec, ef, nb, nk, ef, libtetra_occ
  use libtetra_mpi_mod, only : comm
  implicit none
  !
  real(8),intent(in) :: eig(nb,nk)
  real(8),intent(out) :: occ(nb,nk)
  !
  integer :: iter, maxiter = 300, ierr
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
     call libtetra_occ(eig,occ)
     !
     sumkmid = sum(occ(1:nb,1:nk)) * 2d0
     call MPI_allREDUCE(MPI_IN_PLACE, sumkmid, 1, &
     &                  MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)
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
  if(iter >= maxiter) stop "libtetra_omp_fermieng"
  !
end subroutine libtetra_mpi_fermieng
!
! Initialize grid
!
subroutine libtetra_mpi_kgrid()
  !
  use libtetra_common, only : nk, nk0, indx1, indx2, indx3, ng, ivvec, fst, lst
  use libtetra_mpi_mod, only : libtetra_fst_and_lst
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
  indx2(1:20,1:6 * nk) = 0
  indx3(1:20 * 6 * nk) = 0
  !
  call libtetra_fst_and_lst()
  !
  nk0 = 0
  do it = fst, lst
     !
     do ii = 1, 20
        !
        do ik = 1, nk0
           !
           if(indx1(ii,it) == indx3(ik)) then
              !
              indx2(ii,it) = ik
              goto 10
              !
           end if
           !
        end do
        !
        nk0 = nk0 + 1
        indx2(ii,it) = nk0
        indx3(nk0) = indx1(ii,it)
        !
10      continue
        !
     end do
     !
  end do
  !
end subroutine libtetra_mpi_kgrid
!
! Compute cnt and dsp
!
subroutine libtetra_fst_and_lst()
  !
  use mpi, only : mpi_comm_size, mpi_comm_rank
  use libtetra_common, only : fst, lst, nk
  use libtetra_mpi_mod, only : comm
  implicit none
  !
  integer :: ii, petot, my_rank, ierr
  integer,allocatable :: cnt(:), dsp(:)
  !
  call MPI_COMM_SIZE(comm, petot, ierr)
  call MPI_COMM_RANK(comm, my_rank, ierr)
  !
  allocate(cnt(0:petot-1), dsp(0:petot-1))
  !
  cnt(0:petot-1)        = 6 * nk / petot
  cnt(0:mod(6 * nk,petot)-1) = 6 * nk / petot + 1
  dsp(0) = 0
  do ii = 1, petot - 1
     dsp(ii) = dsp(ii - 1) + cnt(ii - 1)
  end do
  !
  fst = dsp(my_rank) + 1
  lst = dsp(my_rank) + cnt(my_rank)
  !
  deallocate(cnt,dsp)
  !
end subroutine libtetra_fst_and_lst
