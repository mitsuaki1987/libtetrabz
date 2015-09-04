module global
  !
  implicit none
  !
  integer,save :: &
  & ltetra,     &
  & ivvec(3,20,6), & 
  & petot,      &
  & my_rank,    &
  & nb,         & ! # of bands
  & ngd(3),     & ! k-point (dense) grid
  & ng(3),      & ! k-point (dense) grid
  & nkd,        & ! # of k. ng(1) * ng(2) * ng(3)
  & nk,         & ! # of k. ng(1) * ng(2) * ng(3)
  & nmf
  !
  real(8),save :: & 
  & wlsm(4,20), & !
  & bvec(3,3)     ! Reciprocal lattice vector
  !
  real(8),allocatable,save :: &
  & eig1(:,:),   & ! KS enegrgy [Ry] (nb,nk)
  & eig2(:,:),   & ! KS enegrgy [Ry] (nb,nk)
  & mf(:)          ! Phonon frequency [Ry] (nm)
  !
  interface
     !
     subroutine read_elph()
     end subroutine read_elph
     !
     subroutine cnt_and_dsp(n,cnt1,dsp1)
       integer,intent(in) :: n
       integer,intent(out) :: cnt1, dsp1
     end subroutine cnt_and_dsp
     !
     subroutine fermi_fuctor()
     end subroutine fermi_fuctor
     !
     subroutine fermi_fuctor2()
     end subroutine fermi_fuctor2
     !
     subroutine tetra_type()
     end subroutine tetra_type
     !     
     subroutine tetraweight(nkd0,dgrid,indx1,indx2,wghtd)
       import nkd, nb, nmf
       integer,intent(in) :: nkd0, dgrid(3,nkd), indx1(20,6 * nkd), indx2(20,6 * nkd)
       complex(8),intent(out) :: wghtd(0:nmf,nb,nb,nkd0)
     end subroutine tetraweight
     !
     subroutine sort(n1,n2,a)
       integer,intent(in) :: n1, n2
       real(8),intent(inout) :: a(n1,n2) 
     end subroutine sort
     !
     subroutine tetra2(ei,ej,w)
       import nmf, nb
       real(8),intent(in) :: ei(4), ej(nb,4)
       real(8),intent(inout) :: w(4,2,0:nmf,nb,4)
     end subroutine tetra2
     !
     subroutine lindhard(ei,ej,w)
       import nmf
       real(8),intent(in) :: ei(4), ej(4)
       real(8),intent(inout) :: w(4,2,0:nmf,4)
     end subroutine lindhard
     !
     function lindhard_1234_0(g1,g2,g3,g4,lng1,lng2,lng3,lng4) result(w)
       real(8),intent(in) :: g1,g2,g3,g4,lng1,lng2,lng3,lng4
       real(8) :: w
     end function lindhard_1234_0
     !
     function lindhard_1234(g1,g2,g3,g4) result(w)
       real(8),intent(in) :: g1,g2,g3,g4
       real(8) :: w(2)
     end function lindhard_1234
     !
     function lindhard_1231_0(g1,g2,g3,lng1,lng2,lng3) result(w)
       real(8),intent(in) :: g1,g2,g3,lng1,lng2,lng3
       real(8) :: w
     end function lindhard_1231_0
     !
     function lindhard_1231(g1,g2,g3) result(w)
       real(8),intent(in) :: g1,g2,g3
       real(8) :: w(2)
     end function lindhard_1231
     !
     function lindhard_1233_0(g1,g2,g3,lng1,lng2,lng3) result(w)
       real(8),intent(in) :: g1,g2,g3,lng1,lng2,lng3
       real(8) :: w
     end function lindhard_1233_0
     !
     function lindhard_1233(g1,g2,g3) result(w)
       real(8),intent(in) :: g1,g2,g3
       real(8) :: w(2)
     end function lindhard_1233
     !
     function lindhard_1221_0(g1,g2,lng1,lng2) result(w)
       real(8),intent(in) :: g1,g2,lng1,lng2
       real(8) :: w
     end function lindhard_1221_0
     !
     function lindhard_1221(g1,g2) result(w)
       real(8),intent(in) :: g1,g2
       real(8) :: w(2)
     end function lindhard_1221
     !
     function lindhard_1211_0(g1,g2,lng1,lng2) result(w)
       real(8),intent(in) :: g1,g2,lng1,lng2
       real(8) :: w
     end function lindhard_1211_0
     !
     function lindhard_1211(g1,g2) result(w)
       real(8),intent(in) :: g1, g2
       real(8) :: w(2)
     end function lindhard_1211
     !
     function lindhard_1222_0(g1,g2,lng1,lng2) result(w)
       real(8),intent(in) :: g1,g2,lng1,lng2
       real(8) :: w
     end function lindhard_1222_0
     !
     function lindhard_1222(g1,g2) result(w)
       real(8),intent(in) :: g1, g2
       real(8) :: w(2)
     end function lindhard_1222
     !
     subroutine interpol_weight(nk,nb,ng,ko,wi,wo)
       integer,intent(in)  :: nk, nb, ng(3)
       real(8),intent(in)  :: ko(3)
       complex(8),intent(in) ::wi(nb)
       complex(8),intent(inout) :: wo(nb,nk)
     end subroutine interpol_weight
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
  use global, only : nb, ng, nk, ngd, nkd, nmf, mf, eig1, eig2, bvec, my_rank
  implicit none
  !
  real(8) :: nelec
  integer :: fi = 10, ierr
  real(8),allocatable :: g2(:,:)
  !
  if(my_rank == 0) then
     !
     open(fi, file = "elph.dat" )
     !
     read(fi,*) ngd(1:3)
     nkd = product(ngd(1:3))
     !
     read(fi,*) nb
     !
     read(fi,*) nelec
     !
     read(fi,*) nmf
     allocate(g2(nmf,nb * nb * nkd), mf(nmf), eig1(nb,nkd), eig2(nb,nkd))
     !
     read(fi,*) bvec(1:3,1:3)
     !
     read(fi,*) mf(1:nmf)
     !
     ! Read |g|^2
     !
     read(fi,*) g2(1:nmf,1:nb * nb * nkd)
     !
     ! Read Kohn-Sham energies
     !
     read(fi,*) eig1(1:nb,1:nkd)
     read(fi,*) eig2(1:nb,1:nkd)
     !
     close(fi)
     !
     deallocate(g2)
     !
     write(*,*) "  Reciprocal lattice vector[2pi/a] : "
     write(*,'(3f10.5)') bvec(1:3,1:3)
     write(*,*) "  # of electrons : ", nelec
     write(*,*) ""
     write(*,*) "  k point grid : ", ngd(1:3)
     write(*,*) "  # of bands : ", nb
     write(*,*) "  # of modes : ", nmf
     !
  end if
  !
  call MPI_BCAST(nb,     1, MPI_INTEGER,          0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(ngd,     3, MPI_INTEGER,          0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(nmf,    1, MPI_INTEGER,          0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(bvec,   9, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  !
  nkd = product(ngd(1:3))
  if(my_rank /= 0) allocate(mf(nmf), eig1(nb,nkd), eig2(nb,nkd))
  !
  call MPI_BCAST(eig1,           nb * nkd, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(eig2,           nb * nkd, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  !
  mf(1:3) = (/1d0, 2d0, 3d0/)
  ng(1:3) = (/8, 8, 8/)
  nk = product(ng(1:3))
  !
  eig1(1:nb,1:nkd) = eig1(1:nb,1:nkd) - 0.322946141634714
  eig2(1:nb,1:nkd) = eig2(1:nb,1:nkd) - 0.322946141634714
  !
end subroutine read_elph
!
! define shortest diagonal line & define type of tetragonal
!
subroutine tetra_type()
  !
  use mpi, only : MPI_COMM_WORLD
  USE global, ONLY : my_rank, wlsm, ivvec, ng, ltetra, bvec
  !
  implicit none
  !
  integer :: itype, i1, i2, i3, it, &
  &          divvec(4,4), ivvec0(4), ierr
  real(8) :: l(4), bvec2(3,3), bvec3(3,4)
  !
  do i1 = 1, 3
     bvec2(1:3,i1) = bvec(1:3,i1) / dble(ng(i1))
  end do
  !
  bvec3(1:3,1) = -bvec2(1:3,1) + bvec2(1:3,2) + bvec2(1:3,3)
  bvec3(1:3,2) =  bvec2(1:3,1) - bvec2(1:3,2) + bvec2(1:3,3)
  bvec3(1:3,3) =  bvec2(1:3,1) + bvec2(1:3,2) - bvec2(1:3,3)
  bvec3(1:3,4) =  bvec2(1:3,1) + bvec2(1:3,2) + bvec2(1:3,3)
  !
  ! length of delta bvec
  !
  do i1 = 1, 4
     l(i1) = dot_product(bvec3(1:3,i1),bvec3(1:3,i1))
  end do
  !
  itype = minloc(l(1:4),1)
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
  it = 0
  do i1 = 1, 3
     do i2 = 1, 3
        if(i2 == i1) cycle
        do i3 = 1, 3
           if(i3 == i1 .or. i3 == i2) cycle
           !
           it = it + 1
           !
           ivvec(1:3,1,it) = ivvec0(1:3)
           ivvec(1:3,2,it) = ivvec(1:3,1,it) + divvec(1:3,i1)
           ivvec(1:3,3,it) = ivvec(1:3,2,it) + divvec(1:3,i2)
           ivvec(1:3,4,it) = ivvec(1:3,3,it) + divvec(1:3,i3)
           !
        end do
     end do
  end do
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
  if(ltetra == 1) then
     !
     if(my_rank == 0) write(*,*) "  Linear tetrahedron method is used."
     !
     wlsm(1:4,1:20) = 0.0d0
     wlsm(1,1) = 1.0d0
     wlsm(2,2) = 1.0d0
     wlsm(3,3) = 1.0d0
     wlsm(4,4) = 1.0d0
     !
  else if(ltetra == 2) then
     !
     if(my_rank == 0) write(*,*) "  Improved tetrahedron method is used."
     !
     wlsm(1, 1: 4) = dble((/1440,    0,   30,    0/))
     wlsm(2, 1: 4) = dble((/   0, 1440,    0,   30/))
     wlsm(3, 1: 4) = dble((/  30,    0, 1440,    0/))
     wlsm(4, 1: 4) = dble((/   0,   30,    0, 1440/))
     !
     wlsm(1, 5: 8) = dble((/ -38,    7,   17,  -28/))
     wlsm(2, 5: 8) = dble((/ -28,  -38,    7,   17/))
     wlsm(3, 5: 8) = dble((/  17,  -28,  -38,    7/))
     wlsm(4, 5: 8) = dble((/   7,   17,  -28,  -38/))
     !
     wlsm(1, 9:12) = dble((/ -56,    9,  -46,    9/))
     wlsm(2, 9:12) = dble((/   9,  -56,    9,  -46/))
     wlsm(3, 9:12) = dble((/ -46,    9,  -56,    9/))
     wlsm(4, 9:12) = dble((/   9,  -46,    9,  -56/))
     !
     wlsm(1,13:16) = dble((/ -38,  -28,   17,    7/))
     wlsm(2,13:16) = dble((/   7,  -38,  -28,   17/))
     wlsm(3,13:16) = dble((/  17,    7,  -38,  -28/))
     wlsm(4,13:16) = dble((/ -28,   17,    7,  -38/))
     !
     wlsm(1,17:20) = dble((/ -18,  -18,   12,  -18/))
     wlsm(2,17:20) = dble((/ -18,  -18,  -18,   12/))
     wlsm(3,17:20) = dble((/  12,  -18,  -18,  -18/))
     wlsm(4,17:20) = dble((/ -18,   12,  -18,  -18/))
     !
     wlsm(1:4,1:20) = wlsm(1:4,1:20) / 1260d0
     !
  else
     !
     write(*,*) "Stop in tetra_type. ltetra = ", ltetra
     call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
     call MPI_FINALIZE(ierr)
     stop
     !
  end if
  !
end subroutine tetra_type
!
! Compute cnt and dsp
!
subroutine cnt_and_dsp(n,cnt1,dsp1)
  !
  use global, only : petot, my_rank
  implicit none
  !
  integer,intent(in) :: n
  integer,intent(out) :: cnt1, dsp1
  !
  integer :: ii
  integer :: cnt(0:petot-1), dsp(0:petot-1)
  !
  cnt(0:petot-1)        = n / petot
  cnt(0:mod(n,petot)-1) = n / petot + 1
  dsp(0) = 0
  do ii = 1, petot - 1
     dsp(ii) = dsp(ii - 1) + cnt(ii - 1)
  end do
  !
  cnt1 = cnt(my_rank)
  dsp1 = dsp(my_rank)
  !
end subroutine cnt_and_dsp
!
! Calculation of weight function (f(1-f')) / (e - e')
!
subroutine fermi_fuctor()
  !
  use mpi, only : MPI_IN_PLACE, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD
  use global, only : nb, nk, nkd, nmf, ng, ngd, ltetra, ivvec, my_rank, &
  &                  tetraweight, interpol_weight, tetra_type, cnt_and_dsp
  implicit none
  !
  complex(8) :: wght(0:nmf,nb,nb,nk)
  !
  integer :: cnt, dsp, nkd0, nt, it, ik, i1, i2, i3, ii, ierr, ikv(3), &
  &          dgrid(3,nkd), indx1(20, 6 * nkd), indx2(20, 6 * nkd), indx3(20 * 6 * nkd)
  real(8) :: kv(3)
  complex(8),allocatable :: wghtd(:,:,:,:)
  !
  ltetra = 2
  call tetra_type()
  !
  nt = 0
  ik = 0
  do i3 = 1, ngd(3)
     do i2  = 1, ngd(2)
        do i1 = 1, ngd(1)
           !
           ik = ik + 1
           dgrid(1:3,ik) = (/i1, i2, i3/) - 1
           !
           do it = 1, 6
              !
              nt = nt + 1
              !
              do ii = 1, 20
                 !
                 ikv(1:3) = dgrid(1:3,ik) + ivvec(1:3,ii,it)
                 ikv(1:3) = modulo(ikv(1:3), ngd(1:3))
                 !
                 indx1(ii,nt) = 1 + ikv(1) + ngd(1) * ikv(2) + ngd(1) * ngd(2) * ikv(3)
                 !
              end do
              !
           end do
           !
        end do
     end do
  end do
  !
  indx2(1:20,1:6 * nkd) = 0
  indx3(1:20 * 6 * nkd) = 0
  !
  call cnt_and_dsp(nkd * 6,cnt,dsp)
  !
  nkd0 = 0
  do it = dsp + 1, dsp + cnt
     !
     do ii = 1, 20
        !
        do ik = 1, nkd0
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
        nkd0 = nkd0 + 1
        indx2(ii,it) = nkd0
        indx3(nkd0) = indx1(ii,it)
        !
10      continue
        !
     end do
     !
  end do
  !
  allocate(wghtd(0:nmf,nb,nb,nkd0))
  !
  call tetraweight(nkd0, dgrid,indx1,indx2,wghtd)
  !
  ! Interpolation of weight
  !
  wght(0:nmf,1:nb,1:nb,1:nk) = 0d0
  do ik = 1, nkd0
     kv(1:3) = dble(dgrid(1:3,indx3(ik))) / dble(ngd(1:3))
     call interpol_weight(nk, (nmf + 1) * nb * nb, ng, kv(1:3), &
     &                     wghtd(0:nmf, 1:nb,1:nb,ik), wght)
  end do
  !
  call MPI_allREDUCE(MPI_IN_PLACE, wght, (nmf + 1) * nb * nb * nk, &
  &                  MPI_DOUBLE_COMPLEX, MPI_SUM,MPI_COMM_WORLD,ierr)
  !
  deallocate(wghtd)
  !
if(my_rank == 0) write(*,*) sum(wght(0:nmf,1:nb,1:nb,1:nk))
  if(my_rank == 0) write(51,'(1e25.15)') wght(0:nmf,1:nb,1:nb,1:nk)
  !
end subroutine fermi_fuctor
!
! Integration weight with tetrahedron method
!
subroutine tetraweight(nkd0,dgrid,indx1,indx2,wghtd)
  !
  USE global, ONLY : nkd, nb, nmf, ngd, ivvec, wlsm, eig1, eig2, &
  &                  cnt_and_dsp, sort, tetra2
  !
  implicit none
  !
  integer,intent(in) :: nkd0, dgrid(3,nkd), indx1(20,6 * nkd), indx2(20,6 * nkd)
  complex(8),intent(out) :: wghtd(0:nmf,nb,nb,nkd0)
  !
  integer :: ik, it, ib, jb, imf, ii, cnt, dsp
  real(8) :: thr = 1d-8, V
  real(8) :: e(4), a(4,4), ei(nb,4), ej(nb,4), ei2(4), ej2(nb,4), &
  &        tmp(10,0:nmf,nb,4), tmp2(10,0:nmf,nb,4), &
  &         w0(4,2,0:nmf,nb,4), w1(4,2,0:nmf,nb), w2(4,2,0:nmf,nb,4)
  !
  call cnt_and_dsp(nkd * 6,cnt,dsp)
  !
  wghtd(0:nmf,1:nb,1:nb,1:nkd0) = 0d0
  !
  w0(1:4,1:2,0:nmf,1:nb,1:4) = 0d0
  do ii = 1, 4
     w0(ii,1:2,0:nmf,1:nb,ii) = 1d0
  end do
  !
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP & SHARED(nkd,cnt,dsp,nb,dgrid,ivvec,ngd,wlsm,nmf,eig1,eig2,wghtd,thr, &
  !$OMP &        indx1,indx2,w0) &
  !$OMP & PRIVATE(ik,it,ii,ib,jb,imf, &
  !$OMP &         ei,ej,ei2,ej2,w1,w2,tmp,tmp2,a,e,V)
  !
  do it = dsp + 1, dsp + cnt
     !
     ei(1:nb, 1:4) = 0d0
     ej(1:nb, 1:4) = 0d0
     do ii = 1, 20
        !
        do ib = 1, nb
           !
           ei(ib, 1:4) = ei(ib, 1:4) + wlsm(1:4,ii) * eig1(ib, indx1(ii,it))
           ej(ib, 1:4) = ej(ib, 1:4) + wlsm(1:4,ii) * eig2(ib, indx1(ii,it))
           !
        end do
        !
     end do
     !
     !$OMP DO
     do ib = 1, nb
        !
        w1(1:4,1:2,0:nmf,1:nb) = 0d0
        !
        tmp(1:10, 0:nmf, 1:nb, 1:4) = 0d0
        tmp(   1,     0,    1, 1:4) = ei(                 ib, 1:4)
        tmp(   2,     0, 1:nb, 1:4) = ej(               1:nb, 1:4)
        tmp(3: 6, 0:nmf, 1:nb, 1:4) = w0(1:4, 1, 0:nmf, 1:nb, 1:4)
        tmp(7:10, 0:nmf, 1:nb, 1:4) = w0(1:4, 2, 0:nmf, 1:nb, 1:4)
        !
        call sort(10 * (nmf + 1) * nb, 4, tmp)
        !
        e(1:4) = tmp(1, 0, 1, 1:4)
        !
        do ii = 1, 4
           a(ii,1:4) = ( 0d0 - e(1:4) ) / (e(ii) - e(1:4))
        end do
        !
        if(e(1) <= 0d0 .and. 0d0 < e(2) ) then
           !
           ! A - 1
           !
           V = a(2,1) * a(3,1) * a(4,1)
           !
           if(V > thr) then
              !
              tmp2(1:10,0:nmf,1:nb,1) = tmp(1:10,0:nmf,1:nb,1)
              tmp2(1:10,0:nmf,1:nb,2) = tmp(1:10,0:nmf,1:nb,1) * a(1,2) &
              &                       + tmp(1:10,0:nmf,1:nb,2) * a(2,1) 
              tmp2(1:10,0:nmf,1:nb,3) = tmp(1:10,0:nmf,1:nb,1) * a(1,3) &
              &                       + tmp(1:10,0:nmf,1:nb,3) * a(3,1)
              tmp2(1:10,0:nmf,1:nb,4) = tmp(1:10,0:nmf,1:nb,1) * a(1,4) &
              &                       + tmp(1:10,0:nmf,1:nb,4) * a(4,1)
              !
              ei2(                    1:4) = tmp2(   1,     0,    1, 1:4)
              ej2(              1:nb, 1:4) = tmp2(   2,     0, 1:nb, 1:4)
              w2(1:4, 1, 0:nmf, 1:nb, 1:4) = tmp2(3: 6, 0:nmf, 1:nb, 1:4)
              w2(1:4, 2, 0:nmf, 1:nb, 1:4) = tmp2(7:10, 0:nmf, 1:nb, 1:4)
              !
              call tetra2(ei2,ej2,w2)
              !
              do ii = 1, 4
                 w1(1:4,1:2,0:nmf,1:nb) = w1(1:4,1:2,0:nmf,1:nb) &
                 &                      + w2(1:4,1:2,0:nmf,1:nb,ii) * V
              end do
              !
           end if
           !
        else if(e(2) <= 0d0 .and. 0d0 < e(3)) then
           !
           ! B - 1
           !
           V = a(3,1) * a(4,1) * a(2,4)
           !
           if(V > thr) then
              !
              tmp2(1:10,0:nmf,1:nb,1) = tmp(1:10,0:nmf,1:nb,1)
              tmp2(1:10,0:nmf,1:nb,2) = tmp(1:10,0:nmf,1:nb,1) * a(1,3) &
              &                       + tmp(1:10,0:nmf,1:nb,3) * a(3,1) 
              tmp2(1:10,0:nmf,1:nb,3) = tmp(1:10,0:nmf,1:nb,1) * a(1,4) &
              &                       + tmp(1:10,0:nmf,1:nb,4) * a(4,1) 
              tmp2(1:10,0:nmf,1:nb,4) = tmp(1:10,0:nmf,1:nb,2) * a(2,4) &
              &                       + tmp(1:10,0:nmf,1:nb,4) * a(4,2) 
              !
              ei2(                    1:4) = tmp2(   1,     0,    1, 1:4)
              ej2(              1:nb, 1:4) = tmp2(   2,     0, 1:nb, 1:4)
              w2(1:4, 1, 0:nmf, 1:nb, 1:4) = tmp2(3: 6, 0:nmf, 1:nb, 1:4)
              w2(1:4, 2, 0:nmf, 1:nb, 1:4) = tmp2(7:10, 0:nmf, 1:nb, 1:4)
              ! 
              call tetra2(ei2,ej2,w2)
              !
              do ii = 1, 4
                 w1(1:4,1:2,0:nmf,1:nb) = w1(1:4,1:2,0:nmf,1:nb) &
                 &                      + w2(1:4,1:2,0:nmf,1:nb,ii) * V
              end do
              !
           end if
           !
           ! B - 2
           !
           V = a(3,2) * a(4,2)
           !
           if(V > thr) then
              !
              tmp2(1:10,0:nmf,1:nb,1:2) = tmp(1:10,0:nmf,1:nb,1:2)
              tmp2(1:10,0:nmf,1:nb,  3) = tmp(1:10,0:nmf,1:nb,2) * a(2,3) &
              &                         + tmp(1:10,0:nmf,1:nb,3) * a(3,2) 
              tmp2(1:10,0:nmf,1:nb,  4) = tmp(1:10,0:nmf,1:nb,2) * a(2,4) &
              &                         + tmp(1:10,0:nmf,1:nb,4) * a(4,2) 
              !
              ei2(                    1:4) = tmp2(   1,     0,    1, 1:4)
              ej2(              1:nb, 1:4) = tmp2(   2,     0, 1:nb, 1:4)
              w2(1:4, 1, 0:nmf, 1:nb, 1:4) = tmp2(3: 6, 0:nmf, 1:nb, 1:4)
              w2(1:4, 2, 0:nmf, 1:nb, 1:4) = tmp2(7:10, 0:nmf, 1:nb, 1:4)
              ! 
              call tetra2(ei2,ej2,w2)
              !
              do ii = 1, 4
                 w1(1:4,1:2,0:nmf,1:nb) = w1(1:4,1:2,0:nmf,1:nb) &
                 &                      + w2(1:4,1:2,0:nmf,1:nb,ii) * V
              end do
              !
           end if
           !
           ! B - 3
           !
           V = a(2,3) * a(3,1) * a(4,2)
           !
           if(V > thr) then
              !
              tmp2(1:10,0:nmf,1:nb,1) = tmp(1:10,0:nmf,1:nb,1)
              tmp2(1:10,0:nmf,1:nb,2) = tmp(1:10,0:nmf,1:nb,1) * a(1,3) &
              &                       + tmp(1:10,0:nmf,1:nb,3) * a(3,1) 
              tmp2(1:10,0:nmf,1:nb,3) = tmp(1:10,0:nmf,1:nb,2) * a(2,3) &
              &                       + tmp(1:10,0:nmf,1:nb,3) * a(3,2) 
              tmp2(1:10,0:nmf,1:nb,4) = tmp(1:10,0:nmf,1:nb,2) * a(2,4) &
              &                       + tmp(1:10,0:nmf,1:nb,4) * a(4,2) 
              !
              ei2(                    1:4) = tmp2(   1,     0,    1, 1:4)
              ej2(              1:nb, 1:4) = tmp2(   2,     0, 1:nb, 1:4)
              w2(1:4, 1, 0:nmf, 1:nb, 1:4) = tmp2(3: 6, 0:nmf, 1:nb, 1:4)
              w2(1:4, 2, 0:nmf, 1:nb, 1:4) = tmp2(7:10, 0:nmf, 1:nb, 1:4)
              ! 
              call tetra2(ei2,ej2,w2)
              !
              do ii = 1, 4
                 w1(1:4,1:2,0:nmf,1:nb) = w1(1:4,1:2,0:nmf,1:nb) &
                 &                      + w2(1:4,1:2,0:nmf,1:nb,ii) * V
              end do
              !
           end if
           !
        else if(e(3) <= 0d0 .and. 0d0 < e(4)) then
           !
           ! C - 1
           !
           V = a(4,3)
           !
           if(V > thr) then
              !
              tmp2(1:10,0:nmf,1:nb,1:3) = tmp(1:10,0:nmf,1:nb,1:3)
              tmp2(1:10,0:nmf,1:nb,  4) = tmp(1:10,0:nmf,1:nb,3) * a(3,4) &
              &                         + tmp(1:10,0:nmf,1:nb,4) * a(4,3) 
              !
              ei2(                    1:4) = tmp2(   1,     0,    1, 1:4)
              ej2(              1:nb, 1:4) = tmp2(   2,     0, 1:nb, 1:4)
              w2(1:4, 1, 0:nmf, 1:nb, 1:4) = tmp2(3: 6, 0:nmf, 1:nb, 1:4)
              w2(1:4, 2, 0:nmf, 1:nb, 1:4) = tmp2(7:10, 0:nmf, 1:nb, 1:4)
              ! 
              call tetra2(ei2,ej2,w2)
              !
              do ii = 1, 4
                 w1(1:4,1:2,0:nmf,1:nb) = w1(1:4,1:2,0:nmf,1:nb) &
                 &                      + w2(1:4,1:2,0:nmf,1:nb,ii) * V
              end do
              !
           end if
           !
           ! C - 2
           !
           V = a(3,4) * a(4,2)
           !
           if(V > thr) then
              !
              tmp2(1:10,0:nmf,1:nb,1:2) = tmp(1:10,0:nmf,1:nb,1:2)
              tmp2(1:10,0:nmf,1:nb,  3) = tmp(1:10,0:nmf,1:nb,2) * a(2,4) &
              &                         + tmp(1:10,0:nmf,1:nb,4) * a(4,2) 
              tmp2(1:10,0:nmf,1:nb,  4) = tmp(1:10,0:nmf,1:nb,3) * a(3,4) &
              &                         + tmp(1:10,0:nmf,1:nb,4) * a(4,3) 
              !
              ei2(                    1:4) = tmp2(   1,     0,    1, 1:4)
              ej2(              1:nb, 1:4) = tmp2(   2,     0, 1:nb, 1:4)
              w2(1:4, 1, 0:nmf, 1:nb, 1:4) = tmp2(3: 6, 0:nmf, 1:nb, 1:4)
              w2(1:4, 2, 0:nmf, 1:nb, 1:4) = tmp2(7:10, 0:nmf, 1:nb, 1:4)
              ! 
              call tetra2(ei2,ej2,w2)
              !
              do ii = 1, 4
                 w1(1:4,1:2,0:nmf,1:nb) = w1(1:4,1:2,0:nmf,1:nb) &
                 &                      + w2(1:4,1:2,0:nmf,1:nb,ii) * V
              end do
              !
           end if
           !
           ! C - 3
           !
           V = a(3,4) * a(2,4) * a(4,1)
           !
           if(V > thr) then
              !
              tmp2(1:10,0:nmf,1:nb,1) = tmp(1:10,0:nmf,1:nb,1)
              tmp2(1:10,0:nmf,1:nb,2) = tmp(1:10,0:nmf,1:nb,1) * a(1,4) &
              &                       + tmp(1:10,0:nmf,1:nb,4) * a(4,1) 
              tmp2(1:10,0:nmf,1:nb,3) = tmp(1:10,0:nmf,1:nb,2) * a(2,4) &
              &                       + tmp(1:10,0:nmf,1:nb,4) * a(4,2) 
              tmp2(1:10,0:nmf,1:nb,4) = tmp(1:10,0:nmf,1:nb,3) * a(3,4) &
              &                       + tmp(1:10,0:nmf,1:nb,4) * a(4,3) 
              !
              ei2(                    1:4) = tmp2(   1,     0,    1, 1:4)
              ej2(              1:nb, 1:4) = tmp2(   2,     0, 1:nb, 1:4)
              w2(1:4, 1, 0:nmf, 1:nb, 1:4) = tmp2(3: 6, 0:nmf, 1:nb, 1:4)
              w2(1:4, 2, 0:nmf, 1:nb, 1:4) = tmp2(7:10, 0:nmf, 1:nb, 1:4)
              ! 
              call tetra2(ei2,ej2,w2)
              !
              do ii = 1, 4
                 w1(1:4,1:2,0:nmf,1:nb) = w1(1:4,1:2,0:nmf,1:nb) &
                 &                      + w2(1:4,1:2,0:nmf,1:nb,ii) * V
              end do
              !
           end if
           !
        else if(e(4) <= 0d0 ) then
           !
           ! D - 1
           !
           V = 1d0
           !             
           tmp2(1:10,0:nmf,1:nb,1:4) = tmp(1:10,0:nmf,1:nb,1:4)
           !
           ei2(                    1:4) = tmp2(   1,     0,    1, 1:4)
           ej2(              1:nb, 1:4) = tmp2(   2,     0, 1:nb, 1:4)
           w2(1:4, 1, 0:nmf, 1:nb, 1:4) = tmp2(3: 6, 0:nmf, 1:nb, 1:4)
           w2(1:4, 2, 0:nmf, 1:nb, 1:4) = tmp2(7:10, 0:nmf, 1:nb, 1:4)
           ! 
           call tetra2(ei2,ej2,w2)
           !
           do ii = 1, 4
              w1(1:4,1:2,0:nmf,1:nb) = w1(1:4,1:2,0:nmf,1:nb) &
              &                      + w2(1:4,1:2,0:nmf,1:nb,ii) * V
           end do
           !
        end if
        !
        do ii = 1, 20
           !
           do jb = 1, nb
              !
              do imf = 0, nmf
                 !
                 wghtd(imf,jb,ib,indx2(ii,it)) = wghtd(imf,jb,ib,indx2(ii,it)) &
                 &          + cmplx(sum(wlsm(1:4,ii) * w1(1:4,1,imf,jb)), &
                 &                  sum(wlsm(1:4,ii) * w1(1:4,2,imf,jb))  )
                 !
              end do ! imf = 0, nmf
              !
           end do ! jb = 1, nb
           !               
        end do ! ii = 1, 20
        !
     end do ! ib
     !$OMP END DO NOWAIT
     !
  end do ! it
  !
  !$OMP END PARALLEL
  !
  wghtd(0:nmf,1:nb,1:nb,1:nkd0) = wghtd(0:nmf,1:nb,1:nb,1:nkd0) / dble(6 * nkd)
  !
end subroutine tetraweight
!
! Simple sort
!
subroutine sort(n1,n2,a)
  !
  implicit none
  !
  integer,intent(in) :: n1, n2
  real(8),intent(inout) :: a(n1,n2) 
  !
  integer :: i, m
  real(8) :: am, atmp(n1)
  !
  do i = 1, n2 - 1
     am = minval(a(1,i+1:n2) )
     m  = minloc(a(1,i+1:n2),1) + i
     if(a(1,i) > am) then
        atmp(1:n1) = a(1:n1,m)
        a(1:n1,m) = a(1:n1,i)
        a(1:n1,i) = atmp(1:n1)
     end if
  end do
  !
end subroutine sort
!
! Tetrahedra method for 1 - f(ep)
!
subroutine tetra2(ei,ej,w)
  !
  USE global, ONLY : sort, lindhard, nb, nmf
  !
  implicit none
  !
  real(8),intent(in) :: ei(4), ej(nb,4)
  real(8),intent(inout) :: w(4,2,0:nmf,nb,4)
  !
  integer :: ii, ib
  real(8) :: V, ei2(4), ej2(4), w2(4,2,0:nmf,4), thr = 1d-8
  real(8) :: tmp(10,0:nmf,4), tmp2(10,0:nmf,4), e(4), a(4,4)
  !
  do ib = 1, nb
     !
     tmp(1:10, 0:nmf, 1:4) = 0d0
     tmp(   1,     0, 1:4) = ej(           ib,1:4)
     tmp(   2,     0, 1:4) = ei(              1:4)
     tmp(3: 6, 0:nmf, 1:4) = w(1:4,1,0:nmf,ib,1:4)
     tmp(7:10, 0:nmf, 1:4) = w(1:4,2,0:nmf,ib,1:4)
     !
     call sort(10 * (nmf + 1), 4, tmp)
     !
     e(1:4) = tmp(1, 0, 1:4)
     !
     do ii = 1, 4
        a(ii,1:4) = ( 0d0 - e(1:4) ) / (e(ii) - e(1:4) )
     end do
     !
     w(1:4,1:2,0:nmf,ib,1:4) = 0d0
     !
     if(0d0 <= e(1) ) then
        !
        ! A - 1
        !
        V = 1d0
        !
        tmp2(1:10,0:nmf,1:4) = tmp(1:10,0:nmf,1:4)
        !
        ej2(             1:4) = tmp2(   1,     0, 1:4)
        ei2(             1:4) = tmp2(   2,     0, 1:4)
        w2(1:4, 1, 0:nmf,1:4) = tmp2(3: 6, 0:nmf, 1:4)
        w2(1:4, 2, 0:nmf,1:4) = tmp2(7:10, 0:nmf, 1:4)
        !
        call lindhard(ei2,ej2,w2)
        w(1:4,1:2,0:nmf,ib,1:4) = w(1:4,1:2,0:nmf,ib,1:4) + w2(1:4,1:2,0:nmf,1:4) * V
        !
     else if((e(1) < 0d0 .and. 0d0 <= e(2)) .or. (e(1) <= 0d0 .and. 0d0 < e(2))) then
        !
        ! B - 1
        !
        V = a(1,2)
        !
        if(V > thr) then
           !
           tmp2(1:10,0:nmf,1)   = tmp(1:10,0:nmf,1) * a(1,2) + tmp(1:10,0:nmf,2) * a(2,1)
           tmp2(1:10,0:nmf,2:4) = tmp(1:10,0:nmf,2:4)
           !
           ej2(             1:4) = tmp2(   1,     0, 1:4)
           ei2(             1:4) = tmp2(   2,     0, 1:4)
           w2(1:4, 1, 0:nmf,1:4) = tmp2(3: 6, 0:nmf, 1:4)
           w2(1:4, 2, 0:nmf,1:4) = tmp2(7:10, 0:nmf, 1:4)
           !
           call lindhard(ei2,ej2,w2)
           w(1:4,1:2,0:nmf,ib,1:4) = w(1:4,1:2,0:nmf,ib,1:4) + w2(1:4,1:2,0:nmf,1:4) * V
           !
        end if
        !
        ! B - 2
        !
        V = a(1,3) * a(2,1)
        !
        if(V > thr) then
           !
           tmp2(1:10,0:nmf,1) = tmp(1:10,0:nmf,1) * a(1,2) + tmp(1:10,0:nmf,2) * a(2,1)
           tmp2(1:10,0:nmf,2) = tmp(1:10,0:nmf,1) * a(1,3) + tmp(1:10,0:nmf,3) * a(3,1)
           tmp2(1:10,0:nmf,3:4) = tmp(1:10,0:nmf,3:4)
           !
           ej2(             1:4) = tmp2(   1,     0, 1:4)
           ei2(             1:4) = tmp2(   2,     0, 1:4)
           w2(1:4, 1, 0:nmf,1:4) = tmp2(3: 6, 0:nmf, 1:4)
           w2(1:4, 2, 0:nmf,1:4) = tmp2(7:10, 0:nmf, 1:4)
           !
           call lindhard(ei2,ej2,w2)
           w(1:4,1:2,0:nmf,ib,1:4) = w(1:4,1:2,0:nmf,ib,1:4) + w2(1:4,1:2,0:nmf,1:4) * V
           !
        end if
        !
        ! B - 3
        !
        V = a(1,4) * a(2,1) * a(3,1)
        !
        if(V > thr) then
           !
           tmp2(1:10,0:nmf,1) = tmp(1:10,0:nmf,1) * a(1,2) + tmp(1:10,0:nmf,2) * a(2,1)
           tmp2(1:10,0:nmf,2) = tmp(1:10,0:nmf,1) * a(1,3) + tmp(1:10,0:nmf,3) * a(3,1)
           tmp2(1:10,0:nmf,3) = tmp(1:10,0:nmf,1) * a(1,4) + tmp(1:10,0:nmf,4) * a(4,1)
           tmp2(1:10,0:nmf,4) = tmp(1:10,0:nmf,4)
           !
           ej2(             1:4) = tmp2(   1,     0, 1:4)
           ei2(             1:4) = tmp2(   2,     0, 1:4)
           w2(1:4, 1, 0:nmf,1:4) = tmp2(3: 6, 0:nmf, 1:4)
           w2(1:4, 2, 0:nmf,1:4) = tmp2(7:10, 0:nmf, 1:4)
           !
           call lindhard(ei2,ej2,w2)
           w(1:4,1:2,0:nmf,ib,1:4) = w(1:4,1:2,0:nmf,ib,1:4) + w2(1:4,1:2,0:nmf,1:4) * V
           !
        end if
        !          
     else if((e(2) < 0d0 .and. 0d0 <= e(3)) .or. (e(2) <= 0d0 .and. 0d0 < e(3))) then
        !          
        ! C - 1
        !
        V = a(2,4) * a(1,4) * a(3,1)
        !
        if(V > thr) then
           !
           tmp2(1:10,0:nmf,1) = tmp(1:10,0:nmf,1) * a(1,3) + tmp(1:10,0:nmf,3) * a(3,1)
           tmp2(1:10,0:nmf,2) = tmp(1:10,0:nmf,1) * a(1,4) + tmp(1:10,0:nmf,4) * a(4,1)
           tmp2(1:10,0:nmf,3) = tmp(1:10,0:nmf,2) * a(2,4) + tmp(1:10,0:nmf,4) * a(4,2)
           tmp2(1:10,0:nmf,4) = tmp(1:10,0:nmf,4)
           !
           ej2(             1:4) = tmp2(   1,     0, 1:4)
           ei2(             1:4) = tmp2(   2,     0, 1:4)
           w2(1:4, 1, 0:nmf,1:4) = tmp2(3: 6, 0:nmf, 1:4)
           w2(1:4, 2, 0:nmf,1:4) = tmp2(7:10, 0:nmf, 1:4)
           !
           call lindhard(ei2,ej2,w2)
           w(1:4,1:2,0:nmf,ib,1:4) = w(1:4,1:2,0:nmf,ib,1:4) + w2(1:4,1:2,0:nmf,1:4) * V
           !
        end if
        !
        ! C - 2
        !
        V = a(1,3) * a(2,3)
        !
        if(V > thr) then
           !
           tmp2(1:10,0:nmf,1) = tmp(1:10,0:nmf,1) * a(1,3) + tmp(1:10,0:nmf,3) * a(3,1)
           tmp2(1:10,0:nmf,2) = tmp(1:10,0:nmf,2) * a(2,3) + tmp(1:10,0:nmf,3) * a(3,2)
           tmp2(1:10,0:nmf,3:4) = tmp(1:10,0:nmf,3:4)
           !
           ej2(             1:4) = tmp2(   1,     0, 1:4)
           ei2(             1:4) = tmp2(   2,     0, 1:4)
           w2(1:4, 1, 0:nmf,1:4) = tmp2(3: 6, 0:nmf, 1:4)
           w2(1:4, 2, 0:nmf,1:4) = tmp2(7:10, 0:nmf, 1:4)
           !
           call lindhard(ei2,ej2,w2)
           w(1:4,1:2,0:nmf,ib,1:4) = w(1:4,1:2,0:nmf,ib,1:4) + w2(1:4,1:2,0:nmf,1:4) * V
           !
        end if
        !
        ! C - 3
        ! 
        V = a(1,3) * a(2,4) * a(3,2)
        !
        if(V > thr) then
           !
           tmp2(1:10,0:nmf,1) = tmp(1:10,0:nmf,1) * a(1,3) + tmp(1:10,0:nmf,3) * a(3,1)
           tmp2(1:10,0:nmf,2) = tmp(1:10,0:nmf,2) * a(2,3) + tmp(1:10,0:nmf,3) * a(3,2)
           tmp2(1:10,0:nmf,3) = tmp(1:10,0:nmf,2) * a(2,4) + tmp(1:10,0:nmf,4) * a(4,2)
           tmp2(1:10,0:nmf,4) = tmp(1:10,0:nmf,4)
           !
           ej2(             1:4) = tmp2(   1,     0, 1:4)
           ei2(             1:4) = tmp2(   2,     0, 1:4)
           w2(1:4, 1, 0:nmf,1:4) = tmp2(3: 6, 0:nmf, 1:4)
           w2(1:4, 2, 0:nmf,1:4) = tmp2(7:10, 0:nmf, 1:4)
           !
           call lindhard(ei2,ej2,w2)
           w(1:4,1:2,0:nmf,ib,1:4) = w(1:4,1:2,0:nmf,ib,1:4) + w2(1:4,1:2,0:nmf,1:4) * V
           !
        end if
        !          
     else if((e(3) < 0d0 .and. 0d0 <= e(4)) .or. (e(3) <= 0d0 .and. 0d0 < e(4))) then
        !
        ! D - 1
        !
        V = a(3,4) * a(2,4) * a(1,4) 
        !          
        if(V > thr) then
           !
           tmp2(1:10,0:nmf,1) = tmp(1:10,0:nmf,1) * a(1,4) + tmp(1:10,0:nmf,4) * a(4,1)
           tmp2(1:10,0:nmf,2) = tmp(1:10,0:nmf,2) * a(2,4) + tmp(1:10,0:nmf,4) * a(4,2)
           tmp2(1:10,0:nmf,3) = tmp(1:10,0:nmf,3) * a(3,4) + tmp(1:10,0:nmf,4) * a(4,3)
           tmp2(1:10,0:nmf,4) = tmp(1:10,0:nmf,4)
           !
           ej2(             1:4) = tmp2(   1,     0, 1:4)
           ei2(             1:4) = tmp2(   2,     0, 1:4)
           w2(1:4, 1, 0:nmf,1:4) = tmp2(3: 6, 0:nmf, 1:4)
           w2(1:4, 2, 0:nmf,1:4) = tmp2(7:10, 0:nmf, 1:4)
           !
           call lindhard(ei2,ej2,w2)
           w(1:4,1:2,0:nmf,ib,1:4) = w(1:4,1:2,0:nmf,ib,1:4) + w2(1:4,1:2,0:nmf,1:4) * V
           !
        end if
        !
     end if
     !
  end do
  !
end subroutine tetra2
!
! Tetarahedra method for delta(om - ep + e)
!
subroutine lindhard(ei,ej,w)
  !
  use mpi, only : MPI_COMM_WORLD
  USE global, ONLY : nmf, mf, sort, &
  &                  lindhard_1234,   lindhard_1231,   lindhard_1233, &
  &                  lindhard_1221,   lindhard_1222,   lindhard_1211, &
  &                  lindhard_1234_0, lindhard_1231_0, lindhard_1233_0, &
  &                  lindhard_1221_0, lindhard_1222_0, lindhard_1211_0
  !
  implicit none
  !
  real(8),intent(in) :: ei(4), ej(4)
  real(8),intent(inout) :: w(4,2,0:nmf,4)
  !
  integer :: ii, imf, ierr
  real(8) :: tmp(9,0:nmf,4), w2(2,4), de(4), de0(4), lnd(4), thr, thr2
  !
  tmp(1:9, 0:nmf, 1:4) = 0d0
  tmp(  1,     0, 1:4) = ej(1:4) - ei(1:4)
  tmp(2:5, 0:nmf, 1:4) = w(1:4, 1, 0:nmf, 1:4)
  tmp(6:9, 0:nmf, 1:4) = w(1:4, 2, 0:nmf, 1:4)
  call sort(9 * (nmf + 1), 4, tmp)
  de0(          1:4) = tmp(  1,     0, 1:4)
  w(1:4,1,0:nmf,1:4) = tmp(2:5, 0:nmf, 1:4)
  w(1:4,2,0:nmf,1:4) = tmp(6:9, 0:nmf, 1:4)
  !
  thr2 = 1d-12
  !
  do ii = 1, 4
     if(de0(ii) < thr2) then
        if(ii == 3) then
           write(*,*) "Stop in lindhard. Nesting ! "
           write(*,'(9e15.5)') de0(1:4)
           call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
           call MPI_FINALIZE(ierr)
           stop
        end if
        lnd(ii) = 0d0
        de0(ii) = 0d0
     else
        lnd(ii) = log(de0(ii))
     end if
  end do
  !
  de(1:4) = de0(1:4)
  thr = maxval(de(1:4)) * 1d-3
  !
  if(abs(de(4) - de(3)) < thr ) then
     if(abs(de(4) - de(2)) < thr ) then
        if(abs(de(4) - de(1)) < thr ) then
           !
           ! de(4) = de(3) = de(2) = de(1)
           !
           w2(1,4) = 0.25d0 / de(4)
           w2(1,3) = w2(1,4)
           w2(1,2) = w2(1,4)
           w2(1,1) = w2(1,4)
           !
        else
           !
           ! de(4) = de(3) = de(2)
           !
           w2(1,4) = lindhard_1211_0(de(4),de(1),lnd(4),lnd(1))
           w2(1,3) = w2(1,4)
           w2(1,2) = w2(1,4)
           w2(1,1) = lindhard_1222_0(de(1),de(4),lnd(1),lnd(4))
           !
           if(any(w2(1,1:4) < 0d0)) then
              write(*,*) "Stop in lindhard. weighting 4=3=2, imf = ", 0
              write(*,'(100e15.5)') de(1:4)
              write(*,'(100e15.5)') w2(1,1:4)
              call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
              call MPI_FINALIZE(ierr)
              stop
           end if
           !
        end if
     else if(abs(de(2) - de(1)) < thr ) then
        !
        ! de(4) = de(3), de(2) = de(1)
        !
        w2(1,4) = lindhard_1221_0(de(4),de(2), lnd(4),lnd(2))
        w2(1,3) = w2(1,4)
        w2(1,2) = lindhard_1221_0(de(2),de(4), lnd(2),lnd(4))
        w2(1,1) = w2(1,2)
        !
        if(any(w2(1,1:4) < 0d0)) then
           write(*,*) "Stop in lindhard. weighting 4=3 2=1, imf = ", 0
           write(*,'(100e15.5)') de(1:4)
           write(*,'(100e15.5)') w2(1,1:4)
           call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
           call MPI_FINALIZE(ierr)
           stop
        end if
        !
     else
        !
        ! de(4) = de(3)
        !
        w2(1,4) = lindhard_1231_0(de(4),de(1),de(2),lnd(4),lnd(1),lnd(2))
        w2(1,3) = w2(1,4)
        w2(1,2) = lindhard_1233_0(de(2),de(1),de(4),lnd(2),lnd(1),lnd(4))
        w2(1,1) = lindhard_1233_0(de(1),de(2),de(4),lnd(1),lnd(2),lnd(4))
        !
        if(any(w2(1,1:4) < 0d0)) then
           write(*,*) "Stop in lindhard. weighting 4=3, imf = ", 0
           write(*,'(100e15.5)') de(1:4)
           write(*,'(100e15.5)') w2(1,1:4)
           call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
           call MPI_FINALIZE(ierr)
           stop
        end if
        !
     end if
  else if(abs(de(3) - de(2)) < thr) then
     if(abs(de(3) - de(1)) < thr) then
        !
        ! de(3) = de(2) = de(1)
        !
        w2(1,4) = lindhard_1222_0(de(4),de(3), lnd(4),lnd(3))
        w2(1,3) = lindhard_1211_0(de(3),de(4), lnd(3),lnd(4))
        w2(1,2) = w2(1,3)
        w2(1,1) = w2(1,3)
        !
        if(any(w2(1,1:4) < 0d0)) then
           write(*,*) "Stop in lindhard. weighting 3=2=1, imf = ", 0
           write(*,'(100e15.5)') de(1:4)
           write(*,'(100e15.5)') w2(1,1:4)
           call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
           call MPI_FINALIZE(ierr)
           stop
        end if
        !
     else
        !
        ! de(3) = de(2)
        !
        w2(1,4) = lindhard_1233_0(de(4),de(1),de(3),lnd(4),lnd(1),lnd(3))
        w2(1,3) = lindhard_1231_0(de(3),de(1),de(4),lnd(3),lnd(1),lnd(4))
        w2(1,2) = w2(1,3)
        w2(1,1) = lindhard_1233_0(de(1),de(4),de(3),lnd(1),lnd(4),lnd(3))
        !
        if(any(w2(1,1:4) < 0d0)) then
           write(*,*) "Stop in lindhard. weighting 3=2, imf = ", 0
           write(*,'(100e15.5)') de(1:4)
           write(*,'(100e15.5)') w2(1,1:4)
           call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
           call MPI_FINALIZE(ierr)
           stop
        end if
        !
     end if
  else if(abs(de(2) - de(1)) < thr) then
     !
     ! de(2) = de(1)
     !
     w2(1,4) = lindhard_1233_0(de(4),de(3),de(2),lnd(4),lnd(3),lnd(2))
     w2(1,3) = lindhard_1233_0(de(3),de(4),de(2),lnd(3),lnd(4),lnd(2))
     w2(1,2) = lindhard_1231_0(de(2),de(3),de(4),lnd(2),lnd(3),lnd(4))
     w2(1,1) = w2(1,2)
     !
     if(any(w2(1,1:4) < 0d0)) then
        write(*,*) "Stop in lindhard. weighting 2=1, imf = ", 0
        write(*,'(100e15.5)') de(1:4)
        write(*,'(100e15.5)') w2(1,1:4)
        call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
        call MPI_FINALIZE(ierr)
        stop
     end if
     !
  else
     !
     ! Different each other.
     !
     w2(1,4) = lindhard_1234_0(de(4),de(1),de(2),de(3),lnd(4),lnd(1),lnd(2),lnd(3))
     w2(1,3) = lindhard_1234_0(de(3),de(1),de(2),de(4),lnd(3),lnd(1),lnd(2),lnd(4))
     w2(1,2) = lindhard_1234_0(de(2),de(1),de(3),de(4),lnd(2),lnd(1),lnd(3),lnd(4))
     w2(1,1) = lindhard_1234_0(de(1),de(2),de(3),de(4),lnd(1),lnd(2),lnd(3),lnd(4))
     !      
     if(any(w2(1,1:4) < 0d0)) then
        write(*,*) "Stop in lindhard. weighting each other, imf = ", 0
        write(*,'(100e15.5)') de(1:4)
        write(*,'(100e15.5)') w2(1,1:4)
        call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
        call MPI_FINALIZE(ierr)
        stop
     end if
     !
  end if
  !
  do ii = 1, 4
     w(1:4,1,0,ii) = w2(1,ii) * w(1:4,1,0,ii)
     w(1:4,2,0,ii) = 0d0
  end do ! ii
  !
  ! w /= 0 part
  !
  do imf = 1, nmf
     !
     de(1:4) = de0(1:4) / mf(imf)
     !thr = maxval(de(1:4)) * 1d-3
     thr = max(1d-3,  maxval(de(1:4)) * 1d-2)
     !
     if(abs(de(4) - de(3)) < thr ) then
        if(abs(de(4) - de(2)) < thr ) then
           if(abs(de(4) - de(1)) < thr ) then
              !
              ! de(4) = de(3) = de(2) = de(1)
              !
              w2(1,4) = 0.25d0 * de(4) / ((1d0 + de(4)**2))
              w2(2,4) = 0.25d0         / ((1d0 + de(4)**2))
              w2(1:2,3) = w2(1:2,4)
              w2(1:2,2) = w2(1:2,4)
              w2(1:2,1) = w2(1:2,4)
              !
           else
              !
              ! de(4) = de(3) = de(2)
              !
              w2(1:2,4) = lindhard_1211(de(4),de(1))
              w2(1:2,3) = w2(1:2,4)
              w2(1:2,2) = w2(1:2,4)
              w2(1:2,1) = lindhard_1222(de(1),de(4))
              !
              if(any(w2(1:2,1:4) < 0d0)) then
                 write(*,*) "Stop in lindhard. weighting 4=3=2. imf = ", imf
                 write(*,'(100e15.5)') de(1:4)
                 write(*,'(2e15.5)') w2(1:2,1:4)
                 call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
                 call MPI_FINALIZE(ierr)
                 stop
              end if
              !
           end if
        else if(abs(de(2) - de(1)) < thr ) then
           !
           ! de(4) = de(3), de(2) = de(1)
           !
           w2(1:2,4) = lindhard_1221(de(4),de(2))
           w2(1:2,3) = w2(1:2,4)
           w2(1:2,2) = lindhard_1221(de(2),de(4))
           w2(1:2,1) = w2(1:2,2)
           !
           if(any(w2(1:2,1:4) < 0d0)) then
              write(*,*) "Stop in lindhard. weighting 4=3, 2=1. imf = ", imf
              write(*,'(100e15.5)') de(1:4)
              write(*,'(2e15.5)') w2(1:2,1:4)
              call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
              call MPI_FINALIZE(ierr)
              stop
           end if
           !
        else
           !
           ! de(4) = de(3)
           !
           w2(1:2,4) = lindhard_1231(de(4),de(1),de(2))
           w2(1:2,3) = w2(1:2,4)
           w2(1:2,2) = lindhard_1233(de(2),de(1),de(4))
           w2(1:2,1) = lindhard_1233(de(1),de(2),de(4))
           !
           if(any(w2(1:2,1:4) < 0d0)) then
              write(*,*) "Stop in lindhard. weighting 4=3. imf = ", imf
              write(*,'(100e15.5)') de(1:4)
              write(*,'(2e15.5)') w2(1:2,1:4)
              call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
              call MPI_FINALIZE(ierr)
              stop
           end if
           !
        end if
     else if(abs(de(3) - de(2)) < thr) then
        if(abs(de(3) - de(1)) < thr) then
           !
           ! de(3) = de(2) = de(1)
           !
           w2(1:2,4) = lindhard_1222(de(4),de(3))
           w2(1:2,3) = lindhard_1211(de(3),de(4))
           w2(1:2,2) = w2(1:2,3)
           w2(1:2,1) = w2(1:2,3)
           !
           if(any(w2(1:2,1:4) < 0d0)) then
              write(*,*) "Stop in lindhard. weighting 3=2=1. imf = ", imf
              write(*,'(100e15.5)') de(1:4)
              write(*,'(2e15.5)') w2(1:2,1:4)
              call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
              call MPI_FINALIZE(ierr)
              stop
           end if
           !
        else
           !
           ! de(3) = de(2)
           !
           w2(1:2,4) = lindhard_1233(de(4),de(1),de(3))
           w2(1:2,3) = lindhard_1231(de(3),de(1),de(4))
           w2(1:2,2) = w2(1:2,3)
           w2(1:2,1) = lindhard_1233(de(1),de(4),de(3))
           !
           if(any(w2(1:2,1:4) < 0d0)) then
              write(*,*) "Stop in lindhard. weighting 3=2. imf = ", imf
              write(*,'(100e15.5)') de(1:4)
              write(*,'(2e15.5)') w2(1:2,1:4)
              call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
              call MPI_FINALIZE(ierr)
              stop
           end if
           !
        end if
     else if(abs(de(2) - de(1)) < thr) then
        !
        ! de(2) = de(1)
        !
        w2(1:2,4) = lindhard_1233(de(4),de(3),de(2))
        w2(1:2,3) = lindhard_1233(de(3),de(4),de(2))
        w2(1:2,2) = lindhard_1231(de(2),de(3),de(4))
        w2(1:2,1) = w2(1:2,2)
        !
        if(any(w2(1:2,1:4) < 0d0)) then
           write(*,*) "Stop in lindhard. weighting 2=1. imf = ", imf
           write(*,'(100e15.5)') de(1:4)
           write(*,'(2e15.5)') w2(1:2,1:4)
           call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
           call MPI_FINALIZE(ierr)
           stop
        end if
        !
     else
        !
        ! Different each other.
        !
        w2(1:2,4) = lindhard_1234(de(4),de(1),de(2),de(3))
        w2(1:2,3) = lindhard_1234(de(3),de(1),de(2),de(4))
        w2(1:2,2) = lindhard_1234(de(2),de(1),de(3),de(4))
        w2(1:2,1) = lindhard_1234(de(1),de(2),de(3),de(4))
        !      
        if(any(w2(1:2,1:4) < 0d0)) then
           write(*,*) "Stop in lindhard. weighting each other. imf = ", imf
           write(*,'(100e15.5)') de(1:4)
           write(*,'(2e15.5)') w2(1:2,1:4)
           call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
           call MPI_FINALIZE(ierr)
           stop
        end if
        !
     end if
     !
     do ii = 1, 4
        w(1:4,1,imf,ii) = w2(1,ii) * w(1:4,1,imf,ii) /    mf(imf)
        w(1:4,2,imf,ii) = w2(2,ii) * w(1:4,2,imf,ii) / (- mf(imf))
     end do ! ii
     !
  end do ! imf
  !
end subroutine lindhard
!
! Results of Integration (1-x-y-z)/(g0+(g1-g0)x+(g2-g0)y+(g3-g0))
!  for 0<x<1, 0<y<1-x, 0<z<1-x-y
!
! 1, Different each other
!
function lindhard_1234_0(g1,g2,g3,g4,lng1,lng2,lng3,lng4) result(w)
  !
  implicit none
  !
  real(8),intent(in) :: g1,g2,g3,g4,lng1,lng2,lng3,lng4
  real(8) :: w
  !
  real(8) :: w2, w3, w4
  !
  w2 = ((lng2 - lng1)/(g2 - g1)*g2 - 1d0)*g2/(g2 - g1)
  w3 = ((lng3 - lng1)/(g3 - g1)*g3 - 1d0)*g3/(g3 - g1)
  w4 = ((lng4 - lng1)/(g4 - g1)*g4 - 1d0)*g4/(g4 - g1)
  w2 = ((w2 - w3)*g2)/(g2 - g3)
  w4 = ((w4 - w3)*g4)/(g4 - g3)
  w = (w4 - w2)/(g4 - g2)
  !
end function lindhard_1234_0
!
! 1, Different each other
!
function lindhard_1234(g1,g2,g3,g4) result(w)
  !
  implicit none
  !
  real(8),intent(in) :: g1, g2, g3, g4
  real(8) :: w(2)
  !
  real(8) :: w2, w3, w4
  !
  ! Real
  !
  w2 = 2d0*(3d0*g2**2 - 1d0)*(atan(g2) - atan(g1)) + (g2**2 - &
  &      3d0)*g2*log((1d0 + g2**2)/( 1d0 + g1**2))
  w2 = -2d0*(g2**2 - 1d0) + w2/(g2 - g1 )
  w2 = w2/(g2 - g1 )
  w3 = 2d0*(3d0*g3**2 - 1d0)*(atan(g3) - atan(g1)) + (g3**2 -  &
  &      3d0)*g3*log((1d0 + g3**2)/( 1d0 + g1**2))
  w3 = -2d0*(g3**2 - 1d0) + w3/(g3 - g1 )
  w3 = w3/(g3 - g1 )
  w4 = 2d0*(3d0*g4**2 - 1d0)*(atan(g4) - atan(g1)) + (g4**2 -  &
  &      3d0)*g4*log((1d0 + g4**2)/( 1d0 + g1**2))
  w4 = -2d0*(g4**2 - 1d0) + w4/(g4 - g1 )
  w4 = w4/(g4 - g1 )
  w2 = (w2 - w3)/(g2 - g3)
  w4 = (w4 - w3)/(g4 - g3)
  w(1) = (w4 - w2)/(2d0*(g4 - g2))
  !
  ! Imaginal
  !
  w2 = 2d0*(3d0 - g2**2)* &
  &    g2*(atan(g2) - atan(g1)) + (3d0*g2**2 - 1d0)* &
  &    log((1d0 + g2**2)/(1d0 + g1**2))
  w2 = 4d0*g2 - w2/(g2 - g1)
  w2 = w2/(g2 - g1)
  w3 = 2d0*(3d0 - g3**2)* &
  &    g3*(atan(g3) - atan(g1)) + (3d0*g3**2 - 1d0)* &
  &    log((1d0 + g3**2)/(1d0 + g1**2))
  w3 = 4d0*g3 - w3/(g3 - g1)
  w3 = w3/(g3 - g1)
  w4 = 2d0*(3d0 - g4**2)* &
  &    g4*(atan(g4) - atan(g1)) + (3d0*g4**2 - 1d0)* &
  &    log((1d0 + g4**2)/(1d0 + g1**2))
  w4 = 4d0*g4 - w4/(g4 - g1)
  w4 = w4/(g4 - g1)
  w2 = (w2 - w3)/(g2 - g3)
  w4 = (w4 - w3)/(g4 - g3)
  w(2) = (w4 - w2)/(2d0*(g4 - g2))
  !
end function lindhard_1234
!
! 2, g4 = g1
!
function lindhard_1231_0(g1,g2,g3,lng1,lng2,lng3) result(w)
  !
  implicit none
  !
  real(8),intent(in) :: g1,g2,g3,lng1,lng2,lng3
  real(8) :: w
  !
  real(8) :: w2, w3
  !
  w2 = ((lng2 - lng1)/(g2 - g1)*g2 - 1d0)*g2**2/(g2 - g1) - g1/( &
  &   2d0)
  w2 = w2/(g2 - g1)
  w3 = ((lng3 - lng1)/(g3 - g1)*g3 - 1d0)*g3**2/(g3 - g1) - g1/( &
  &   2d0)
  w3 = w3/(g3 - g1)
  w = (w3 - w2)/(g3 - g2)
  !
end function lindhard_1231_0
!
! 2, g4 = g1
!
function lindhard_1231(g1,g2,g3) result(w)
  !
  implicit none
  !
  real(8),intent(in) :: g1, g2, g3
  real(8) :: w(2)
  !
  real(8) :: w2, w3
  !
  ! Real
  !
  w2 = 2d0*(-1d0 + 3d0*g2**2)*(atan(g2) - atan(g1)) +  &
  &   g2*(-3d0 + g2**2)*log((1d0 + g2**2)/(1d0 + g1**2))
  w2 = 2d0*(1d0 - g2**2) + w2/(g2 - g1)
  w2 = -g1 + w2/(g2 - g1)
  w2 = w2/(g2 - g1)
  w3 = 2d0*(-1d0 + 3d0*g3**2)*(atan(g3) - atan(g1)) +  &
  &   g3*(-3d0 + g3**2)*log((1d0 + g3**2)/(1d0 + g1**2))
  w3 = 2d0*(1 - g3**2) + w3/(g3 - g1)
  w3 = -g1 + w3/(g3 - g1)
  w3 = w3/(g3 - g1)
  w(1) = (w3 - w2)/(2d0*(g3 - g2))
  !
  ! Imaginal
  !
  w2 = 2d0* &
  &    g2*(3d0 - g2**2)*(atan(g2) - atan(g1)) + (-1d0 + 3d0*g2**2)* &
  &    log((1d0 + g2**2)/(1d0 + g1**2))
  w2 = 4d0*g2 - w2/(g2 - g1)
  w2 = 1 + w2/(g2 - g1)
  w2 = w2/(g2 - g1)
  w3 = 2d0* &
  &    g3*(3d0 - g3**2)*(atan(g3) - atan(g1)) + (-1d0 + 3d0*g3**2)* &
  &    log((1d0 + g3**2)/(1d0 + g1**2))
  w3 = 4d0*g3 - w3/(g3 - g1)
  w3 = 1 + w3/(g3 - g1)
  w3 = w3/(g3 - g1)
  w(2) = (w3 - w2)/(2d0*(g3 - g2))
  !
end function lindhard_1231
!
! 3, g4 = g3
!
function lindhard_1233_0(g1,g2,g3,lng1,lng2,lng3) result(w)
  !
  implicit none
  !
  real(8),intent(in) :: g1,g2,g3,lng1,lng2,lng3
  real(8) :: w
  !
  real(8) :: w2, w3
  !
  w2 = (lng2 - lng1)/(g2 - g1)*g2 - 1d0
  w2 = (g2*w2)/(g2 - g1)
  w3 = (lng3 - lng1)/(g3 - g1)*g3 - 1d0
  w3 = (g3*w3)/(g3 - g1)
  w2 = (w3 - w2)/(g3 - g2)
  w3 = (lng3 - lng1)/(g3 - g1)*g3 - 1d0
  w3 = 1d0 - (2d0*w3*g1)/(g3 - g1)
  w3 = w3/(g3 - g1)
  w = (g3*w3 - g2*w2)/(g3 - g2)
  !
end function lindhard_1233_0
!
! 3, g4 = g3
!
function lindhard_1233(g1, g2, g3) result(w)
  !
  implicit none
  !
  real(8),intent(in) :: g1, g2, g3
  real(8) :: w(2)
  !
  real(8) :: w2, w3
  !
  ! Real
  !
  w2 = 2d0*(1d0 - 3d0*g2**2)*(atan(g2) - atan(g1)) +  &
  &   g2*(3d0 - g2**2)*log((1d0 + g2**2)/(1d0 + g1**2))
  w2 = 2d0*(1 - g2**2) - w2/(g2 - g1)
  w2 = w2/(g2 - g1)
  w3 = 2d0*(1d0 - 3d0*g3**2)*(atan(g3) - atan(g1)) +  &
  &   g3*(3d0 - g3**2)*log((1d0 + g3**2)/(1d0 + g1**2))
  w3 = 2d0*(1 - g3**2) - w3/(g3 - g1)
  w3 = w3/(g3 - g1)
  w2 = (w3 - w2)/(g3 - g2)
  w3 = 4d0*(1d0 - 3d0*g1*g3)*(atan(g3) - atan(g1)) + (3d0*g1 +  &
  &      3d0*g3 - 3d0*g1*g3**2 + g3**3) * log((1d0 + g3**2)/( &
  &     1d0 + g1**2))
  w3 = -4d0*(1d0 - g1**2) + w3/(g3 - g1)
  w3 = 4d0*g1 + w3/(g3 - g1)
  w3 = w3/(g3 - g1)
  w(1) = (w3 - w2)/(2d0*(g3 - g2))
  !
  ! Imaginal
  !
  w2 = 2d0* &
  &    g2*(3d0 - g2**2)*(atan(g2) - atan(g1)) + (-1d0 + 3d0*g2**2)* &
  &    log((1d0 + g2**2)/(1d0 + g1**2))
  w2 = 4d0*g2 - w2/(g2 - g1)
  w2 = w2/(g2 - g1)
  w3 = 2d0* &
  &    g3*(3d0 - g3**2)*(atan(g3) - atan(g1)) + (-1d0 + 3d0*g3**2)* &
  &    log((1d0 + g3**2)/(1d0 + g1**2))
  w3 = 4d0*g3 - w3/(g3 - g1)
  w3 = w3/(g3 - g1)
  w2 = (w3 - w2)/(g3 - g2)
  w3 = (3d0*g1 - 3d0*g1*g3**2 + 3d0*g3 + g3**3)*(atan(g3) -  &
  &      atan(g1)) + (3d0*g1*g3 - 1d0)* &
  &    log((1d0 + g3**2)/(1d0 + g1**2))
  w3 = w3/(g3 - g1) - 4d0*g1
  w3 = w3/(g3 - g1) - 2d0
  w3 = (2d0*w3)/(g3 - g1)
  w(2) = (w3 - w2)/(2d0*(g3 - g2))
  !
end function lindhard_1233
!
! 4, g4 = g1 and g3 = g2
!
function lindhard_1221_0(g1,g2,lng1,lng2) result(w)
  !
  implicit none
  !
  real(8),intent(in) :: g1, g2, lng1, lng2
  real(8) :: w
  !
  w = 1d0 - (lng2 - lng1)/(g2 - g1)*g1
  w = -1d0 + (2d0*g2*w)/(g2 - g1)
  w = -1d0 + (3d0*g2*w)/(g2 - g1)
  w = w/(2d0*(g2 - g1))
  !
end function lindhard_1221_0
!
! 4, g4 = g1 and g3 = g2
!
function lindhard_1221(g1,g2) result(w)
  !
  implicit none
  !
  real(8),intent(in) :: g1, g2
  real(8) :: w(2)
  !
  ! Real
  !
  w(1) = -2d0*(-1d0 + 2d0*g1*g2 + g2**2)*(atan(g2) -  &
  &      atan(g1)) + (g1 + 2d0*g2 - g1*g2**2)* &
  &    log((1d0 + g2**2)/(1d0 + g1**2))
  w(1) = 2d0*(-1d0 + g1**2) + w(1)/(g2 - g1)
  w(1) = 3d0*g1 + w(1)/(g2 - g1)
  w(1) = 2d0 + (3d0*w(1))/(g2 - g1)
  w(1) = w(1)/(2d0*(g2 - g1))
  !
  ! Imaginal
  !
  w(2) = 2d0*(g1 + 2d0*g2 - g1*g2**2)*(atan(g2) -  &
  &      atan(g1)) + (-1d0 + 2d0*g1*g2 + g2**2)* &
  &    log((1 + g2**2)/(1 + g1**2))
  w(2) = -4d0*g1 + w(2)/(g2 - g1)
  w(2) = -3d0 + w(2)/(g2 - g1)
  w(2) = (3d0*w(2))/(2d0*(g2 - g1)**2)
  !
end function lindhard_1221
!
! 5, g4 = g3 = g2
!
function lindhard_1222_0(g1,g2,lng1,lng2) result(w)
  !
  implicit none
  !
  real(8),intent(in) :: g1, g2, lng1, lng2
  real(8) :: w
  !
  w = (lng2 - lng1)/(g2 - g1)*g2 - 1d0
  w = (2d0*g1*w)/(g2 - g1) - 1d0
  w = (3d0*g1*w)/(g2 - g1) + 1d0
  w = w/(2d0*(g2 - g1))
  !
end function lindhard_1222_0
!
! 5, g4 = g3 = g2
!
function lindhard_1222(g1,g2) result(w)
  !
  implicit none
  !
  real(8),intent(in) :: g1, g2
  real(8) :: w(2)
  !
  ! Real
  !
  w(1) = 2d0*(-1d0 + g1**2 + 2d0*g1*g2)*(atan(g2) -  &
  &      atan(g1)) + (-2d0*g1 - g2 + g1**2*g2) * log((1d0 + g2**2)/( &
  &     1d0 + g1**2))
  w(1) = 2d0*(1d0 - g1**2) + w(1)/(g2 - g1)
  w(1) = g1 - w(1)/(g2 - g1)
  w(1) = 1d0 - (3d0*w(1))/(g2 - g1)
  w(1) = w(1)/(2d0*(g2 - g1))
  !
  ! Imaginal
  !
  w(2) = 2d0*(-2d0*g1 - g2 + g1**2*g2)*(atan(g2) - atan(g1)) + (1d0 - &
  &       g1**2 - 2d0*g1*g2) * log((1d0 + g2**2)/(1d0 + g1**2))
  w(2) = 4d0*g1 + w(2)/(g2 - g1)
  w(2) = 1d0 + w(2)/(g2 - g1)
  w(2) = (3d0*w(2))/(2d0*(g2 - g1)**2)
  !
end function lindhard_1222
!
! 6, g4 = g3 = g1
!
function lindhard_1211_0(g1,g2,lng1,lng2) result(w)
  !
  implicit none
  !
  real(8),intent(in) :: g1,g2,lng1,lng2
  real(8) :: w
  !
  w = -1d0 + (lng2 - lng1)/(g2 - g1)*g2
  w = -1d0 + (2d0*g2*w)/(g2 - g1)
  w = -1d0 + (3d0*g2*w)/(2d0*(g2 - g1))
  w = w/(3d0*(g2 - g1))
  !
end function lindhard_1211_0
!
! 6, g4 = g3 = g1
!
function lindhard_1211(g1,g2) result(w)
  !
  implicit none
  !
  real(8),intent(in) :: g1, g2
  real(8) :: w(2)
  !
  ! Real
  !
  w(1) = 2d0*(3d0*g2**2 - 1d0)*(atan(g2) - atan(g1)) +  &
  &   g2*(g2**2 - 3d0)*log((1d0 + g2**2)/(1d0 + g1**2))
  w(1) = 2d0*(1d0 - g1**2) + w(1)/(g2 - g1)
  w(1) = -5d0*g1 + w(1)/(g2 - g1)
  w(1) = -11d0 + (3d0*w(1))/(g2 - g1)
  w(1) = w(1)/(6d0*(g2 - g1))
  !
  ! Imaginal
  !
  w(2) = 2d0*g2*(-3d0 + g2**2)*(atan(g2) - atan(g1)) + (1d0 -  &
  &      3d0*g2**2)*log((1d0 + g2**2)/(1d0 + g1**2))
  w(2) = 4d0*g2 + w(2)/(g2 - g1)
  w(2) = 1d0 + w(2)/(g2 - g1)
  w(2) = w(2)/(2d0*(g2 - g1)**2)
  !
end function lindhard_1211
!
! first or third order interpolation of weights
!
subroutine interpol_weight(nk,nb,ng,ko,wi,wo)
  !
  use mpi, only : MPI_COMM_WORLD
  use global, only : ivvec, ltetra
  implicit none
  !
  integer,intent(in)  :: nk, nb, ng(3)
  real(8),intent(in)  :: ko(3)
  complex(8),intent(in) ::wi(nb)
  complex(8),intent(inout) :: wo(nb,nk)
  !
  integer :: ikv(3), ikv1(3), ik(20), ii, it, it0, ierr
  real(8) :: rot(3,3), res(3), prod(3), u, x, y, z, thr = 1d-10
  !
  rot(1:3,1) = (/  2d0, - 1d0,   0d0/)
  rot(1:3,2) = (/- 1d0,   2d0, - 1d0/)
  rot(1:3,3) = (/  0d0, - 1d0,   1d0/)
  !
  ! Search nearest neighbor grid points.
  !
  res(1:3) = ko(1:3) * dble(ng(1:3))
  ikv(1:3) = floor(res(1:3))
  res(1:3) = res(1:3) - dble(ikv(1:3))
  !
  do it = 1, 6
     !
     do ii = 1, 3
        prod(ii) = dot_product(dble(ivvec(1:3,1 + ii,it) - ivvec(1:3,1,it)), &
        &                                  res(1:3) - dble(ivvec(1:3,1,it))  )
     end do
     !
     prod(1:3) = matmul(rot(1:3,1:3), prod(1:3))
     !
     if(minval(prod(1:3)) > - thr .and. sum(prod(1:3)) < 1d0 + thr) then
        it0 = it
        goto 10
     end if
     !
  end do
  !
  write(*,*) "Stop in interpol weight."
  call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
  call MPI_FINALIZE(ierr)
  stop
  !
10 continue
  !
  x = prod(1)
  y = prod(2)
  z = prod(3)
  u = 1d0 - x - y - z
  !
  do ii = 1, 20
     !
     ikv1(1:3) = ikv(1:3) + ivvec(1:3,ii,it0)
     ikv1(1:3) = modulo(ikv1(1:3), ng(1:3))
     ik(ii) = 1 + ikv1(1) + ikv1(2) * ng(1) + ikv1(3) * ng(1) * ng(2)
     !
  end do
  !
  if(ltetra == 0 .or. ltetra == 1) then
     !
     wo(1:nb,ik(1)) = wo(1:nb,ik(1)) + wi(1:nb) * u
     wo(1:nb,ik(2)) = wo(1:nb,ik(2)) + wi(1:nb) * x
     wo(1:nb,ik(3)) = wo(1:nb,ik(3)) + wi(1:nb) * y
     wo(1:nb,ik(4)) = wo(1:nb,ik(4)) + wi(1:nb) * z
     !
  else if(ltetra == 2) then
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
  else
     !
     write(*,*) "Stop in interpol_weight. ltetra = ", ltetra
     call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
     call MPI_FINALIZE(ierr)
     stop
     ! 
  end if
  !
end subroutine interpol_weight
!
!
!
subroutine fermi_fuctor2()
  !
  use mpi
  use global, only : ltetra, bvec, nb, ngd, eig1, eig2, ng, nk, nmf, mf, my_rank
  use libtetrabz_mpi, only : libtetrabz_mpi_polstat, libtetrabz_mpi_polimg
  !
  implicit none
  !
  integer :: ib, jb, ik, imf
  real(8) :: wghts(nb,nb,nk), wghti(2,nmf,nb,nb,nk)
  complex(8) :: wghti2(0:nmf,nb,nb,nk)
  !
  call libtetrabz_mpi_polstat(ltetra,MPI_COMM_WORLD,bvec,nb,ngd,eig1,eig2,ng,wghts)
  !
  call libtetrabz_mpi_polimg( ltetra,MPI_COMM_WORLD,bvec,nb,ngd,eig1,eig2,ng,wghti,nmf,mf)
  !
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP & SHARED(wghti,wghti2,wghts,nk,nb,nmf) &
  !$OMP & PRIVATE(ik,ib,jb,imf)
  !
  !$OMP DO
  do ik = 1, nk
     do ib = 1, nb
        do jb = 1, nb
           wghti2(0,jb,ib,ik) = cmplx(wghts(jb,ib,ik), 0d0)
           do imf = 1, nmf
              wghti2(imf,jb,ib,ik) = cmplx(wghti(1,imf,jb,ib,ik), wghti(2,imf,jb,ib,ik))
           end do
        end do
     end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL
  !
if(my_rank == 0) write(*,*) sum(wghti2(0:nmf,1:nb,1:nb,1:nk))
  if(my_rank == 0) write(61,'(1e25.15)') wghti2(0:nmf,1:nb,1:nb,1:nk)
  !
end subroutine fermi_fuctor2
!
! Maing routine
!
program rpa
  !
  use mpi
  use omp_lib
  use global, only : my_rank, petot, read_elph, fermi_fuctor, fermi_fuctor2
  implicit none
  !
  integer :: ierr
  !
  ! MPI Initialize
  !
  call MPI_INIT(ierr)
  call MPI_COMM_SIZE (MPI_COMM_WORLD, PETOT, ierr)
  call MPI_COMM_RANK (MPI_COMM_WORLD, my_rank, ierr)
  !
  if(my_rank == 0) write(*,*) "  # of PEs : ", petot
  !$OMP PARALLEL
  if(my_rank == 0 .and. OMP_GET_THREAD_NUM() ==0) &
  & write(*,*) '  # of thread : ', OMP_GET_NUM_THREADS()
  !$OMP END PARALLEL
  !
  call read_elph()
  !
  call fermi_fuctor()
  !
  call fermi_fuctor2()
  !
  call MPI_FINALIZE(ierr)
  !
end program rpa
