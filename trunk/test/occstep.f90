MODULE global
  !
  implicit none
  !
  integer,save :: &
  & sttk,     &
  & lstk,     &
  & km_tetra, &
  & ivvec(3,20,6), &
  & petot,      &
  & my_rank,    &
  & nb,         & ! # of bands
  & ng(3),      & ! k-point (dense) grid
  & nk
  !
  real(8),save :: & 
  & wlsm(4,20), &
  & bvec(3,3)
  !
  integer,allocatable,save :: &
  & grid(:,:), &
  & indx(:,:,:)
  !
  real(8),allocatable,save :: &
  & beta(:,:,:), &
  & eig1(:,:),   & ! KS enegrgy [Ry] (nb,nk)
  & eig2(:,:)      ! KS enegrgy [Ry] (nb,nk)
  !
  interface
     !
     subroutine km_eig_and_grid()
     end subroutine km_eig_and_grid
     !
     subroutine km_tetra_type()
     end subroutine km_tetra_type
     !
     subroutine km_calc_beta1()
     end subroutine km_calc_beta1
     !
     subroutine km_tetra2_theta(ei,ej,w)
       import nb
       real(8),intent(in) :: ei(4), ej(nb,4)
       real(8),intent(inout) :: w(4,nb,4)
     end subroutine km_tetra2_theta
     !
     subroutine km_constweight(ei,ej,w)
       real(8),intent(in) :: ei(4), ej(4)
       real(8),intent(inout) :: w(4,4)
     end subroutine km_constweight
     !
     subroutine km_sort(n1,n2,a)
       integer,intent(in) :: n1, n2
       real(8),intent(inout) :: a(n1,n2) 
     end subroutine km_sort
     !
  end interface
  !
end MODULE global
!
! compute eig1, eig2, indx, grid
!
subroutine km_eig_and_grid()
  !
  use mpi
  USE global, ONLY : nb, ng, nk, grid, indx, eig1, eig2, my_rank, bvec, km_tetra
  !
  implicit none
  !
  INTEGER :: ik, i1, i2, i3, nm, fi = 10, ierr
  real(8) :: nelec
  real(8),allocatable :: omg(:), g2(:,:)
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
     deallocate(g2,omg)
     !
  end if
  !
  call MPI_BCAST(nb,     1, MPI_INTEGER,          0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(ng,     3, MPI_INTEGER,          0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(bvec,   9, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  !
  nk = product(ng(1:3))
  if(my_rank /= 0) allocate(eig1(nb,nk), eig2(nb,nk))
  !
  call MPI_BCAST(eig1,           nb * nk, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(eig2,           nb * nk, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  !
  eig1(1:nb,1:nk) = eig1(1:nb,1:nk) - 0.322946141634714
  eig2(1:nb,1:nk) = eig2(1:nb,1:nk) - 0.322946141634714
  !
  allocate(grid(3,nk), indx(ng(1),ng(2),ng(3)))
  !
  ! K-grid
  !
  ik = 0
  do i3 = 1, ng(3)
     do i2 = 1, ng(2)
        do i1 = 1, ng(1)
           !
           ik = ik + 1
           grid(1:3,ik) = (/ i1, i2, i3 /) - 1
           indx(i1, i2, i3) = ik
           !
        end do
     end do
  end do
  !
  km_tetra = 2
  !
end subroutine km_eig_and_grid
!
! define shortest diagonal line & define type of tetragonal
!
subroutine km_tetra_type()
  !
  USE global, ONLY : wlsm, ivvec, ng, bvec, km_tetra, my_rank
  !
  implicit none
  !
  integer :: itype, i1, i2, i3, it, &
  &          divvec(4,4), ivvec0(4)
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
  if(km_tetra == 1) then
     !
     if(my_rank ==0) write(*,*) "KM  Linear tetrahedron method is used."
     !
     wlsm(1:4,1:20) = 0.0d0
     wlsm(1,1) = 1.0d0
     wlsm(2,2) = 1.0d0
     wlsm(3,3) = 1.0d0
     wlsm(4,4) = 1.0d0
     !
  else if(km_tetra == 2) then
     !
     if(my_rank ==0) write(*,*) "KM  Improved tetrahedron method is used."
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
     stop
     !
  end if
  !
end subroutine km_tetra_type
!
! Integration weight(1st term)
!
subroutine km_calc_beta1()
  !
  USE global, ONLY : sttk, lstk, grid, ng, ivvec, wlsm, beta, indx, nb, &
  &                     eig1, eig2, km_sort, km_tetra2_theta
  !
  implicit none
  !
  integer :: ik, it, ib, jb, ikv(3), ii, ikk
  real(8) :: thr = 1d-8, V, e(4), a(4,4)
  real(8) :: ei(nb,4), ej(nb,4), ei2(4), ej2(nb,4)
  real(8) :: tmp(6,nb,4), tmp2(6,nb,4)
  real(8) :: w1(4,nb), w2(4,nb,4)
  real(8),allocatable :: w0(:,:,:)
  !
  allocate(w0(4,nb,4))
  !
  w0(1:4,1:nb,1:4) = 0d0
  do ii = 1, 4
     w0(ii,1:nb,ii) = 1d0
  end do
  tmp(1:6,1:nb,1:4) = 0d0
  !
  ! Integration
  !
  do ik = sttk, lstk
     !
     do it = 1, 6
        !
        ei(1:nb, 1:4) = 0d0
        ej(1:nb, 1:4) = 0d0
        !
        do ii = 1, 20
           !
           ikv(1:3) = grid(1:3,ik) + ivvec(1:3,ii,it)
           ikv(1:3) = modulo(ikv(1:3),ng(1:3)) + 1
           ikk = indx(ikv(1),ikv(2),ikv(3))
           !
           do ib = 1, nb
              !
              ei(ib, 1:4) = ei(ib, 1:4) + wlsm(1:4,ii) * eig1(ib, ikk)
              ej(ib, 1:4) = ej(ib, 1:4) + wlsm(1:4,ii) * eig2(ib, ikk)
              !               
           end do
           !
        end do
        !
        do ib = 1, nb
           !
           w1(1:4,1:nb) = 0d0
           !
           tmp(1,        1, 1:4) = ei(         ib, 1:4)
           tmp(2,   1:nb, 1:4) = ej(     1:nb, 1:4)
           tmp(3:6, 1:nb, 1:4) = w0(1:4, 1:nb, 1:4)
           !
           call km_sort(6 * nb, 4, tmp)
           !
           e(1:4) = tmp(1, 1, 1:4)
           !
           do ii = 1, 4
              a(ii,1:4) = ( 0d0 - e(1:4)) / (e(ii) - e(1:4))
           end do
           !
           if(e(1) <= 0d0 .and. 0d0 < e(2)) then
              !
              ! A - 1
              !
              V = a(2,1) * a(3,1) * a(4,1)
              !
              if(V > thr) then
                 !
                 tmp2(1:6,1:nb,1) = tmp(1:6,1:nb,1)
                 tmp2(1:6,1:nb,2) = tmp(1:6,1:nb,1) * a(1,2) &
                 &                  + tmp(1:6,1:nb,2) * a(2,1) 
                 tmp2(1:6,1:nb,3) = tmp(1:6,1:nb,1) * a(1,3) &
                 &                  + tmp(1:6,1:nb,3) * a(3,1)
                 tmp2(1:6,1:nb,4) = tmp(1:6,1:nb,1) * a(1,4) &
                 &                  + tmp(1:6,1:nb,4) * a(4,1)
                 !
                 ei2(            1:4) = tmp2(  1,      1, 1:4)
                 ej2(1:nb,     1:4) = tmp2(  2, 1:nb, 1:4)
                 w2(1:4, 1:nb, 1:4) = tmp2(3:6, 1:nb, 1:4)
                 ! 
                 call km_tetra2_theta(ei2,ej2,w2)
                 !
                 w1(1:4,1:nb) = w1(1:4,1:nb) &
                 &      + v * sum(w2(1:4,1:nb,1:4), 3)
                 !
              end if
              !
           else if( e(2) <= 0d0 .and. 0d0 < e(3)) then
              !
              ! B - 1
              !
              V = a(3,1) * a(4,1) * a(2,4)
              !
              if(V > thr) then
                 !
                 tmp2(1:6,1:nb,1) = tmp(1:6,1:nb,1)
                 tmp2(1:6,1:nb,2) = tmp(1:6,1:nb,1) * a(1,3) &
                 &                  + tmp(1:6,1:nb,3) * a(3,1) 
                 tmp2(1:6,1:nb,3) = tmp(1:6,1:nb,1) * a(1,4) &
                 &                  + tmp(1:6,1:nb,4) * a(4,1) 
                 tmp2(1:6,1:nb,4) = tmp(1:6,1:nb,2) * a(2,4) &
                 &                  + tmp(1:6,1:nb,4) * a(4,2) 
                 !
                 ei2(            1:4) = tmp2(  1,      1, 1:4)
                 ej2(1:nb,     1:4) = tmp2(  2, 1:nb, 1:4)
                 w2(1:4, 1:nb, 1:4) = tmp2(3:6, 1:nb, 1:4)
                 ! 
                 call km_tetra2_theta(ei2,ej2,w2)
                 !
                 w1(1:4,1:nb) = w1(1:4,1:nb) &
                 &      + V * sum(w2(1:4,1:nb,1:4), 3)
                 !
              end if
              !
              ! B - 2
              !
              V = a(3,2) * a(4,2)
              !
              if(V > thr) then
                 !
                 tmp2(1:6,1:nb,1:2) = tmp(1:6,1:nb,1:2)
                 tmp2(1:6,1:nb,3)   = tmp(1:6,1:nb,2) * a(2,3) &
                 &                    + tmp(1:6,1:nb,3) * a(3,2) 
                 tmp2(1:6,1:nb,4)   = tmp(1:6,1:nb,2) * a(2,4) &
                 &                    + tmp(1:6,1:nb,4) * a(4,2) 
                 !
                 ei2(            1:4) = tmp2(  1,      1, 1:4)
                 ej2(1:nb,     1:4) = tmp2(  2, 1:nb, 1:4)
                 w2(1:4, 1:nb, 1:4) = tmp2(3:6, 1:nb, 1:4)
                 ! 
                 call km_tetra2_theta(ei2,ej2,w2)
                 !
                 w1(1:4,1:nb) = w1(1:4,1:nb) &
                 &      + V * sum(w2(1:4,1:nb,1:4), 3)
                 !
              end if
              !
              ! B - 3
              !
              V = a(2,3) * a(3,1) * a(4,2)
              !
              if(V > thr) then
                 !
                 tmp2(1:6,1:nb,1) = tmp(1:6,1:nb,1)
                 tmp2(1:6,1:nb,2) = tmp(1:6,1:nb,1) * a(1,3) &
                 &                  + tmp(1:6,1:nb,3) * a(3,1) 
                 tmp2(1:6,1:nb,3) = tmp(1:6,1:nb,2) * a(2,3) &
                 &                  + tmp(1:6,1:nb,3) * a(3,2) 
                 tmp2(1:6,1:nb,4) = tmp(1:6,1:nb,2) * a(2,4) &
                 &                  + tmp(1:6,1:nb,4) * a(4,2) 
                 !
                 ei2(            1:4) = tmp2(  1,      1, 1:4)
                 ej2(1:nb,     1:4) = tmp2(  2, 1:nb, 1:4)
                 w2(1:4, 1:nb, 1:4) = tmp2(3:6, 1:nb, 1:4)
                 ! 
                 call km_tetra2_theta(ei2,ej2,w2)
                 !
                 w1(1:4,1:nb) = w1(1:4,1:nb) &
                 &      + V * sum(w2(1:4,1:nb,1:4), 3)
                 !
              end if
              !
           else if( e(3) <= 0d0 .and. 0d0 < e(4)) then
              !
              ! C - 1
              !
              V = a(4,3)
              !
              if(V > thr) then
                 !
                 tmp2(1:6,1:nb,1:3) = tmp(1:6,1:nb,1:3)
                 tmp2(1:6,1:nb,4)   = tmp(1:6,1:nb,3) * a(3,4) &
                 &                    + tmp(1:6,1:nb,4) * a(4,3) 
                 !
                 ei2(            1:4) = tmp2(  1,      1, 1:4)
                 ej2(1:nb,     1:4) = tmp2(  2, 1:nb, 1:4)
                 w2(1:4, 1:nb, 1:4) = tmp2(3:6, 1:nb, 1:4)
                 ! 
                 call km_tetra2_theta(ei2,ej2,w2)
                 !
                 w1(1:4,1:nb) = w1(1:4,1:nb) &
                 &      + V * sum(w2(1:4,1:nb,1:4), 3)
                 !
              end if
              !
              ! C - 2
              !
              V = a(3,4) * a(4,2)
              !
              if(V > thr) then
                 !
                 tmp2(1:6,1:nb,1:2) = tmp(1:6,1:nb,1:2)
                 tmp2(1:6,1:nb,3)   = tmp(1:6,1:nb,2) * a(2,4) &
                 &                    + tmp(1:6,1:nb,4) * a(4,2) 
                 tmp2(1:6,1:nb,4)   = tmp(1:6,1:nb,3) * a(3,4) &
                 &                    + tmp(1:6,1:nb,4) * a(4,3) 
                 !
                 ei2(            1:4) = tmp2(  1,      1, 1:4)
                 ej2(1:nb,     1:4) = tmp2(  2, 1:nb, 1:4)
                 w2(1:4, 1:nb, 1:4) = tmp2(3:6, 1:nb, 1:4)
                 ! 
                 call km_tetra2_theta(ei2,ej2,w2)
                 !
                 w1(1:4,1:nb) = w1(1:4,1:nb) &
                 &      + V * sum(w2(1:4,1:nb,1:4), 3)
                 !
              end if
              !
              ! C - 3
              !
              V = a(3,4) * a(2,4) * a(4,1)
              !
              if(V > thr) then
                 !
                 tmp2(1:6,1:nb,1) = tmp(1:6,1:nb,1)
                 tmp2(1:6,1:nb,2) = tmp(1:6,1:nb,1) * a(1,4) &
                 &                  + tmp(1:6,1:nb,4) * a(4,1) 
                 tmp2(1:6,1:nb,3) = tmp(1:6,1:nb,2) * a(2,4) &
                 &                  + tmp(1:6,1:nb,4) * a(4,2) 
                 tmp2(1:6,1:nb,4) = tmp(1:6,1:nb,3) * a(3,4) &
                 &                  + tmp(1:6,1:nb,4) * a(4,3) 
                 !
                 ei2(            1:4) = tmp2(  1,      1, 1:4)
                 ej2(1:nb,     1:4) = tmp2(  2, 1:nb, 1:4)
                 w2(1:4, 1:nb, 1:4) = tmp2(3:6, 1:nb, 1:4)
                 ! 
                 call km_tetra2_theta(ei2,ej2,w2)
                 !
                 w1(1:4,1:nb) = w1(1:4,1:nb) &
                 &      + V * sum(w2(1:4,1:nb,1:4), 3)
                 !
              end if
              !
           else if( e(4) <= 0d0 ) then
              !
              ! D - 1
              !
              V = 1d0
              !             
              tmp2(1:6,1:nb,1:4) = tmp(1:6,1:nb,1:4)
              !
              ei2(            1:4) = tmp2(  1,      1, 1:4)
              ej2(1:nb,     1:4) = tmp2(  2, 1:nb, 1:4)
              w2(1:4, 1:nb, 1:4) = tmp2(3:6, 1:nb, 1:4)
              ! 
              call km_tetra2_theta(ei2,ej2,w2)
              !
              w1(1:4,1:nb) = w1(1:4,1:nb) &
              &      + V * sum(w2(1:4,1:nb,1:4), 3)
              !
           end if
           !
           do ii = 1, 20
              !
              ikv(1:3) = grid(1:3,ik) + ivvec(1:3,ii,it)
              ikv(1:3) = modulo(ikv(1:3), ng(1:3)) + 1
              ikk = indx(ikv(1),ikv(2),ikv(3))
              !
              beta(1:nb,ib,ikk) = beta(1:nb,ib,ikk) &
              &               + matmul(wlsm(1:4,ii), w1(1:4,1:nb))
              !               
           end do
           !
        end do ! ib
        !
     end do ! it
     !
  end do ! ik
  !
  deallocate(w0)
  !
end subroutine km_calc_beta1
!
! Tetrahedra method for theta(e - ep)
!
subroutine km_tetra2_theta(ei,ej,w)
  !
  USE global, ONLY : km_sort, km_constweight, nb
  !
  implicit none
  !
  real(8),intent(in) :: ei(4), ej(nb,4)
  real(8),intent(inout) :: w(4,nb,4)
  !
  integer :: ii, ib
  real(8) :: V, ei2(4), ej2(4), w2(4,4), thr = 1d-8
  real(8) :: tmp(6,4), tmp2(6,4), e(4), a(4,4)
  !
  do ib = 1, nb
     !
     tmp(1,1:4) = ej(ib,1:4) - ei(1:4)
     tmp(2,1:4) = ei(1:4)
     tmp(3:6,1:4) = w(1:4,ib,1:4)
     !
     call km_sort(6, 4, tmp)
     !
     e(1:4) = tmp(1,1:4)
     !
     do ii = 1, 4
        a(ii,1:4) = (0d0 - e(1:4)) / (e(ii) - e(1:4))
     end do
     !
     w(1:4,ib,1:4) = 0d0
     !
     if(abs(e(1)) < thr .and. abs(e(4)) < thr) then
        !
        ! Theta(0) = 0.5
        !
        V = 0.5d0
        !
        do ii = 1, 4
           ei2(ii) = tmp(2,ii)
           ej2(ii) = tmp(1,ii)
        end do
        w2(1:4,1:4) = tmp(3:6, 1:4)
        !
        call km_constweight(ei2,ej2,w2)
        w(1:4, ib,1:4) = w( 1:4, ib, 1:4) &
        &              + w2(1:4,     1:4) * V
        !
     else if((e(1) <= 0d0 .and. 0d0 < e(2)) .or. (e(1) < 0d0 .and. 0d0 <= e(2))) then
        !
        ! A - 1
        !
        V = a(2,1) * a(3,1) * a(4,1)
        !
        if(V > thr) then
           !
           tmp2(1:6,1) = tmp(1:6,1)
           tmp2(1:6,2) = tmp(1:6,1) * a(1,2) + tmp(1:6,2) * a(2,1) 
           tmp2(1:6,3) = tmp(1:6,1) * a(1,3) + tmp(1:6,3) * a(3,1)
           tmp2(1:6,4) = tmp(1:6,1) * a(1,4) + tmp(1:6,4) * a(4,1)
           !
           do ii = 1, 4
              ei2(ii) = tmp2(2,ii)
              ej2(ii) = tmp2(1,ii)
           end do
           w2(1:4,1:4) = tmp2(3:6, 1:4)
           !
           call km_constweight(ei2,ej2,w2)
           w(1:4, ib,1:4) = w( 1:4, ib, 1:4) &
           &              + w2(1:4,     1:4) * V
           !
        end if
        !
     else if((e(2) <= 0d0 .and. 0d0 < e(3)) .or. (e(2) < 0d0 .and. 0d0 <= e(3))) then
        !
        ! B - 1
        !
        V = a(3,1) * a(4,1) * a(2,4)
        !
        if(V > thr) then
           !
           tmp2(1:6,1) = tmp(1:6,1)
           tmp2(1:6,2) = tmp(1:6,1) * a(1,3) + tmp(1:6,3) * a(3,1) 
           tmp2(1:6,3) = tmp(1:6,1) * a(1,4) + tmp(1:6,4) * a(4,1) 
           tmp2(1:6,4) = tmp(1:6,2) * a(2,4) + tmp(1:6,4) * a(4,2) 
           !
           do ii = 1, 4
              ei2(ii) = tmp2(2,ii)
              ej2(ii) = tmp2(1,ii)
           end do
           w2(1:4,1:4) = tmp2(3:6, 1:4)
           ! 
           call km_constweight(ei2,ej2,w2)
           w(1:4, ib,1:4) = w( 1:4, ib, 1:4) &
           &              + w2(1:4,     1:4) * V
           !
        end if
        !
        ! B - 2
        !
        V = a(3,2) * a(4,2)
        !
        if(V > thr) then
           !
           tmp2(1:6,1:2) = tmp(1:6,1:2)
           tmp2(1:6,3)   = tmp(1:6,2) * a(2,3) + tmp(1:6,3) * a(3,2) 
           tmp2(1:6,4)   = tmp(1:6,2) * a(2,4) + tmp(1:6,4) * a(4,2) 
           !
           do ii = 1, 4
              ei2(ii) = tmp2(2,ii)
              ej2(ii) = tmp2(1,ii)
           end do
           w2(1:4,1:4) = tmp2(3:6, 1:4)
           ! 
           call km_constweight(ei2,ej2,w2)
           w(1:4, ib,1:4) = w( 1:4, ib, 1:4) &
           &              + w2(1:4,     1:4) * V
           !
        end if
        !
        ! B - 3
        !
        V = a(2,3) * a(3,1) * a(4,2)
        !
        if(V > thr) then
           !
           tmp2(1:6,1) = tmp(1:6,1)
           tmp2(1:6,2) = tmp(1:6,1) * a(1,3) + tmp(1:6,3) * a(3,1) 
           tmp2(1:6,3) = tmp(1:6,2) * a(2,3) + tmp(1:6,3) * a(3,2) 
           tmp2(1:6,4) = tmp(1:6,2) * a(2,4) + tmp(1:6,4) * a(4,2) 
           !
           do ii = 1, 4
              ei2(ii) = tmp2(2,ii)
              ej2(ii) = tmp2(1,ii)
           end do
           w2(1:4,1:4) = tmp2(3:6, 1:4)
           ! 
           call km_constweight(ei2,ej2,w2)
           w(1:4, ib,1:4) = w( 1:4, ib, 1:4) &
           &              + w2(1:4,     1:4) * V
           !
        end if
        !
     else if((e(3) <= 0d0 .and. 0d0 < e(4)) .or. (e(3) < 0d0 .and. 0d0 <= e(4))) then
        !
        ! C - 1
        !
        V = a(4,3)
        !
        if(V > thr) then
           !
           tmp2(1:6,1:3) = tmp(1:6,1:3)
           tmp2(1:6,4)   = tmp(1:6,3) * a(3,4) + tmp(1:6,4) * a(4,3) 
           !
           do ii = 1, 4
              ei2(ii) = tmp2(2,ii)
              ej2(ii) = tmp2(1,ii)
           end do
           w2(1:4,1:4) = tmp2(3:6, 1:4)
           ! 
           call km_constweight(ei2,ej2,w2)
           w(1:4, ib,1:4) = w( 1:4, ib, 1:4) &
           &              + w2(1:4,     1:4) * V
           !
        end if
        !
        ! C - 2
        !
        V = a(3,4) * a(4,2)
        !
        if(V > thr) then
           !
           tmp2(1:6,1:2) = tmp(1:6,1:2)
           tmp2(1:6,3)   = tmp(1:6,2) * a(2,4) + tmp(1:6,4) * a(4,2) 
           tmp2(1:6,4)   = tmp(1:6,3) * a(3,4) + tmp(1:6,4) * a(4,3) 
           !
           do ii = 1, 4
              ei2(ii) = tmp2(2,ii)
              ej2(ii) = tmp2(1,ii)
           end do
           w2(1:4,1:4) = tmp2(3:6, 1:4)
           ! 
           call km_constweight(ei2,ej2,w2)
           w(1:4, ib,1:4) = w( 1:4, ib, 1:4) &
           &              + w2(1:4,     1:4) * V
           !
        end if
        !
        ! C - 3
        !
        V = a(3,4) * a(2,4) * a(4,1)
        !
        if(V > thr) then
           !
           tmp2(1:6,1) = tmp(1:6,1)
           tmp2(1:6,2) = tmp(1:6,1) * a(1,4) + tmp(1:6,4) * a(4,1) 
           tmp2(1:6,3) = tmp(1:6,2) * a(2,4) + tmp(1:6,4) * a(4,2) 
           tmp2(1:6,4) = tmp(1:6,3) * a(3,4) + tmp(1:6,4) * a(4,3) 
           !
           do ii = 1, 4
              ei2(ii) = tmp2(2,ii)
              ej2(ii) = tmp2(1,ii)
           end do
           w2(1:4,1:4) = tmp2(3:6, 1:4)
           ! 
           call km_constweight(ei2,ej2,w2)
           w(1:4, ib,1:4) = w( 1:4, ib, 1:4) &
           &              + w2(1:4,     1:4) * V
           !
        end if
        !
     else if(e(4) <= 0d0) then
        !
        ! D - 1
        !
        V = 1d0
        !             
        tmp2(1:6,1:4) = tmp(1:6,1:4)
        !
        do ii = 1, 4
           ei2(ii) = tmp2(2,ii)
           ej2(ii) = tmp2(1,ii)
        end do
        w2(1:4,1:4) = tmp2(3:6, 1:4)
        ! 
        call km_constweight(ei2,ej2,w2)
        w(1:4, ib,1:4) = w( 1:4, ib, 1:4) &
        &              + w2(1:4,     1:4) * V
        !
     end if
     !
  end do
  !
end subroutine km_tetra2_theta
!
! Tetarahedra method for delta(om - ep + e)
!
subroutine km_constweight(ei,ej,w)
  !
  implicit none
  !
  real(8),intent(in) :: ei(4), ej(4)
  real(8),intent(inout) :: w(4,4)
  !
  w(1:4,1:4) = w(1:4,1:4) / 4d0
  !
end subroutine km_constweight
!
! Simple sort
!
subroutine km_sort(n1,n2,a)
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
end subroutine km_sort
!
!
! Main routine
!
program km_beta_main
  !
  use mpi
  use global, only : my_rank, petot, nb, nk, sttk, lstk, beta, ng, bvec, km_tetra, &
  &                  eig1, eig2, km_eig_and_grid, km_tetra_type, km_calc_beta1
  use libtetrabz_mpi, only : libtetrabz_mpi_occstep
  !
  implicit none
  !
  integer :: nkpp, rest, ierr
  !
  call MPI_INIT(ierr)
  call MPI_COMM_SIZE (MPI_COMM_WORLD, PETOT, ierr)
  call MPI_COMM_RANK (MPI_COMM_WORLD, my_rank, ierr)
  !
  CALL km_eig_and_grid()
  !
  ! Work Sharing
  !
  nkpp = nk / petot
  rest = mod(nk, petot)
  if(my_rank < rest) then
     sttk = (nkpp + 1) *  my_rank + 1
     lstk = (nkpp + 1) * (my_rank + 1)
  else
     sttk = nkpp *  my_rank + 1  + rest
     lstk = nkpp * (my_rank + 1) + rest
  end if
  !
  call km_tetra_type()
  !
  allocate(beta(nb,nb,nk))
  beta(1:nb,1:nb,1:nk) = 0d0
  CALL km_calc_beta1()
  !
  call MPI_allREDUCE(MPI_IN_PLACE, beta, nb * nb * nk, &
  &                  MPI_DOUBLE_PRECISION, MPI_SUM,MPI_COMM_WORLD,ierr)
  !
  beta(1:nb,1:nb,1:nk) = beta(1:nb, 1:nb,1:nk) / dble(6 * nk)
  !
  if(my_rank == 0) write(52,'(1e25.15)') beta(1:nb,1:nb,1:nk)
  !
  call libtetrabz_mpi_occstep(km_tetra,MPI_COMM_WORLD,bvec,nb,ng,eig1,eig2,ng,beta)
  !
  if(my_rank == 0) write(62,'(1e25.15)') beta(1:nb,1:nb,1:nk)
  !
  call MPI_FINALIZE(ierr)
  !
end program km_beta_main

