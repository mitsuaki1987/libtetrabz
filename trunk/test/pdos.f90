module global
  !
  implicit none
  !
  integer,save ::  &
  & ivvec(3,20,6), &
  & ne,            &
  & nk,            & !
  & nb,            & !
  & nsym,          & !
  & ng(3),         &
  & nwfc
  !
  real(8),save :: &
  & wlsm(4,20), &
  & bvec(3,3)
  !
  integer,save,allocatable :: &
  & grid(:,:),  & ! (3, nk)
  & sym(:,:,:)    ! (3,3,nsym)
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
     subroutine rotate_k(nk0,kv,eig0,rwfc0)
       import nb, nwfc
       integer,intent(in) :: nk0
       real(8),intent(in) :: eig0(nb,nk0), rwfc0(nb,nwfc,nk0), kv(3,nk0)
     end subroutine rotate_k
     !
     subroutine tetra_type()
     end subroutine tetra_type
     !
     subroutine en_step()
     end subroutine en_step
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
  use iotk_module
  use global, only : nb, nsym, ng, nwfc, bvec, sym, e0, dos, rotate_k, ne
  !
  implicit none
  !
  integer :: fi = 10, ik, iwfc, isym, nk0, ie
  real(8) :: ef, alat, avec(3,3), emin, emax, de
  logical :: invsym
  real(8),allocatable :: eig0(:,:), rwfc0(:,:,:), kv(:,:)
  complex(8),allocatable :: cwfc(:)
  character(iotk_namlenx) :: attr, attr2
  !
  write(*,'(a,a)') "   open atomic_proj.xml"
  call iotk_open_read(fi,"atomic_proj.xml")
  !
  call iotk_scan_begin(fi,"HEADER")
  !
  call iotk_scan_dat(fi,"NUMBER_OF_BANDS",nb)
  write(*,*) "  # of bands : ", nb
  !
  call iotk_scan_dat(fi,"NUMBER_OF_K-POINTS",nk0)
  write(*,*) "  # of k : ", nk0
  !
  call iotk_scan_dat(fi,"NUMBER_OF_ATOMIC_WFC",nwfc)
  write(*,*) "  # of atomic wfc : ", nwfc
  !
  call iotk_scan_dat(fi,"FERMI_ENERGY",eF)
  ! 
  call iotk_scan_end(fi,"HEADER")
  !
  allocate(kv(3,nk0), eig0(nb,nk0), cwfc(nb), rwfc0(nb,nwfc,nk0))
  !
  ! Read k-vector
  !
  call iotk_scan_dat(fi,"K-POINTS",kv)
  !
  ! Read eigenvalue
  !    
  call iotk_scan_begin(fi,"EIGENVALUES")
  !
  do ik = 1, nk0
     !
     write(attr2,*) ik
     write(attr,'(a,a)') "K-POINT.", trim(adjustl(attr2))
     call iotk_scan_begin(fi,trim(attr))
     !
     call iotk_scan_dat(fi,"EIG", eig0(1:nb,ik))
     !
     call iotk_scan_end(fi,trim(attr))
     !
  end do
  !
  eig0(1:nb,1:nk0) = eig0(1:nb,1:nk0) - ef
  !
  call iotk_scan_end(fi,"EIGENVALUES")
  !
  ! Read projections
  !
  call iotk_scan_begin(fi,"PROJECTIONS")
  !
  do ik = 1, nk0
     !
     write(attr2,*) ik
     write(attr,'(a,a)') "K-POINT.", trim(adjustl(attr2))
     call iotk_scan_begin(fi,trim(attr))
     !
     do iwfc = 1, nwfc
        !
        write(attr2,*) iwfc
        write(attr,'(a,a)') "ATMWFC.", trim(adjustl(attr2))
        !
        call iotk_scan_dat(fi,trim(attr),cwfc(1:nb))
        !
        rwfc0(1:nb,iwfc,ik) = dble(conjg(cwfc(1:nb)) * cwfc(1:nb))
        !
     end do
     !
     write(attr2,*) ik
     write(attr,'(a,a)') "K-POINT.", trim(adjustl(attr2))
     call iotk_scan_end(fi,trim(attr))
     !
  end do
  !
  call iotk_scan_end(fi,"PROJECTIONS")
  !
  call iotk_close_read(fi)
  !
  ! Read data-file.xml
  !
  write(*,'(a,a)') "   open data-file.xml"
  call iotk_open_read(fi,"data-file.xml")
  !
  ! Read CELL PARAMETER
  !
  call iotk_scan_begin(fi,"CELL")
  !
  ! Read lattice parameter
  !
  call iotk_scan_dat(fi,"LATTICE_PARAMETER",alat)
  write(*,*) "  Lattice parameter[a.u.] : ", alat
  !
  ! Read Lattice vector
  !
  call iotk_scan_begin(fi,"DIRECT_LATTICE_VECTORS")
  call iotk_scan_dat(fi,"a1",avec(1:3,1))
  call iotk_scan_dat(fi,"a2",avec(1:3,2))
  call iotk_scan_dat(fi,"a3",avec(1:3,3))
  avec(1:3,1:3) = avec(1:3,1:3) / alat
  call iotk_scan_end(fi,"DIRECT_LATTICE_VECTORS")
  write(*,*) "  Direct lattice vector[a] : "
  write(*,'(3f10.5)') avec(1:3,1:3)
  !
  ! Read reciprocal lattice vecor
  !
  call iotk_scan_begin(fi,"RECIPROCAL_LATTICE_VECTORS")
  call iotk_scan_dat(fi,"b1",bvec(1:3,1))
  call iotk_scan_dat(fi,"b2",bvec(1:3,2))
  call iotk_scan_dat(fi,"b3",bvec(1:3,3))
  call iotk_scan_end(fi,"RECIPROCAL_LATTICE_VECTORS")
  !
  write(*,*) "  Reciprocal lattice vector[2pi/a] : "
  write(*,'(3f10.5)') bvec(1:3,1:3)
  !
  call iotk_scan_end(fi,"CELL")
  !
  ! Read # of point symmetry
  !
  call iotk_scan_begin(fi,"SYMMETRIES")
  call iotk_scan_dat(fi,"NUMBER_OF_SYMMETRIES",nsym)
  call iotk_scan_dat(fi,"INVERSION_SYMMETRY",invsym)
  if(invsym) then
     allocate(sym(3,3,nsym))
     write(*,*) "  # of BZ symmetry : ", nsym
  else
     allocate(sym(3,3,nsym * 2))
     write(*,*) "  Inversion symmetry is added."
     write(*,*) "  # of BZ symmetry : ", nsym * 2
  end if
  !
  ! Read symmmetry operators
  !
  do isym = 1, nsym
     !
     write(attr2,*) isym
     write(attr,'(a,a)') "SYMM.", trim(adjustl(attr2))
     !
     call iotk_scan_begin(fi,trim(attr))
     call iotk_scan_dat(fi,"ROTATION",sym(1:3,1:3,isym))
     call iotk_scan_end(fi,trim(attr))
     !
  end do
  !
  if(.not. invsym) then
     sym(1:3,1:3,nsym + 1:nsym + nsym) = - sym(1:3,1:3,1:nsym)
     nsym = nsym * 2
  end if
  !
  call iotk_scan_end(fi,"SYMMETRIES")
  !
  ! Read Monkhorst-Pack grid
  !
  call iotk_scan_begin(fi,"BRILLOUIN_ZONE")
  !
  attr=""
  call iotk_scan_empty(fi,"MONKHORST_PACK_GRID",attr)
  call iotk_scan_attr(attr,"nk1",ng(1))
  call iotk_scan_attr(attr,"nk2",ng(2))
  call iotk_scan_attr(attr,"nk3",ng(3))
  !
  write(*,*) "  MONKHORST_PACK_GRID : ", ng(1:3)
  call iotk_scan_end(fi,"BRILLOUIN_ZONE")
  !
  call iotk_close_read(fi)
  !
  ! Grid of K
  !
  do ik = 1, nk0
     kv(1:3,ik) = matmul(kv(1:3,ik), avec(1:3,1:3))
  end do
  !
  call rotate_k(nk0,kv,eig0,rwfc0)
  !
  deallocate(cwfc,kv,eig0,rwfc0)
  !
  nwfc = nwfc + 1
  !
end subroutine read_dat
!
! Rotation of k-points
!
subroutine rotate_k(nk0,kv,eig0,rwfc0)
  !
  use global, only : nk, nb, nwfc, ng, nsym ,sym, eig, wfc
  !
  implicit none
  !
  integer,intent(in) :: nk0
  real(8),intent(in) :: eig0(nb,nk0), rwfc0(nb,nwfc,nk0), kv(3,nk0)
  !
  integer :: isym, ik, ikv2(3), iwfc
  real(8) :: kv2(3)
  logical :: ldone(ng(1), ng(2), ng(3))
  !
  allocate(eig(nb, ng(1), ng(2), ng(3)), wfc(nwfc + 1, nb, ng(1), ng(2), ng(3)))
  !
  ldone(1:ng(1),1:ng(2),1:ng(3)) = .false.
  wfc(1,1:nb,1:ng(1),1:ng(2),1:ng(3)) = 1d0
  !
  do ik = 1, nk0
     !
     do isym = 1, nsym
        !
        kv2(1:3) = matmul(dble(sym(1:3,1:3,isym)), kv(1:3,ik)) * dble(ng(1:3))
        ikv2(1:3) = nint(kv2(1:3))
        !
        if(any(abs(kv2(1:3) - dble(ikv2(1:3))) > 1d-8)) cycle
        !
        ikv2(1:3) = modulo(ikv2(1:3), ng(1:3)) + 1
        !
        ldone(ikv2(1), ikv2(2),ikv2(3)) = .true.
        !
        eig(1:nb,ikv2(1), ikv2(2),ikv2(3)) = eig0(1:nb,ik)
        do iwfc = 1, nwfc
           wfc(iwfc + 1,1:nb,ikv2(1), ikv2(2),ikv2(3)) = rwfc0(1:nb,iwfc,ik)
        end do
        !
     end do ! End isym
     !
  end do ! End ik
  !
  ! Check
  !
  if(count(.not. ldone) /= 0) &
  &     write(*,*)  "  # of elements that are not done : ", count(.not. ldone)
  !
end subroutine rotate_k
!
! define shortest diagonal line & define type of tetragonal
!
subroutine tetra_type()
  !
  use global, only : ng, bvec, ivvec, wlsm
  implicit none
  !
  integer :: itype, i1, i2, i3, itet, divvec(4,4), ivvec0(4)
  real(8) :: l(4), bvec2(3,3), bvec3(3,4)
  !
  do i1 = 1, 3
     bvec2(1:3,i1) = bvec(1:3,i1) / real(ng(i1), 8)
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
  itet = 0
  do i1 = 1, 3
     do i2 = 1, 3
        if(i2 == i1) cycle
        do i3 = 1, 3
           if(i3 == i1 .or. i3 == i2) cycle
           !
           itet = itet + 1
           !
           ivvec(1:3,1,itet) = ivvec0(1:3)
           ivvec(1:3,2,itet) = ivvec(1:3,1,itet) + divvec(1:3,i1)
           ivvec(1:3,3,itet) = ivvec(1:3,2,itet) + divvec(1:3,i2)
           ivvec(1:3,4,itet) = ivvec(1:3,3,itet) + divvec(1:3,i3)
           !
        end do
     end do
  end do
  !
  ivvec(1:3, 5,1:6) = 2 * ivvec(1:3,1,1:6) - ivvec(1:3,2,1:6)
  ivvec(1:3, 6,1:6) = 2 * ivvec(1:3,1,1:6) - ivvec(1:3,3,1:6)
  ivvec(1:3, 7,1:6) = 2 * ivvec(1:3,1,1:6) - ivvec(1:3,4,1:6)
  !
  ivvec(1:3, 8,1:6) = 2 * ivvec(1:3,2,1:6) - ivvec(1:3,1,1:6)
  ivvec(1:3, 9,1:6) = 2 * ivvec(1:3,2,1:6) - ivvec(1:3,3,1:6)
  ivvec(1:3,10,1:6) = 2 * ivvec(1:3,2,1:6) - ivvec(1:3,4,1:6)
  !
  ivvec(1:3,11,1:6) = 2 * ivvec(1:3,3,1:6) - ivvec(1:3,1,1:6)
  ivvec(1:3,12,1:6) = 2 * ivvec(1:3,3,1:6) - ivvec(1:3,4,1:6)
  ivvec(1:3,13,1:6) = 2 * ivvec(1:3,3,1:6) - ivvec(1:3,2,1:6)
  !
  ivvec(1:3,14,1:6) = 2 * ivvec(1:3,4,1:6) - ivvec(1:3,1,1:6)
  ivvec(1:3,15,1:6) = 2 * ivvec(1:3,4,1:6) - ivvec(1:3,2,1:6)
  ivvec(1:3,16,1:6) = 2 * ivvec(1:3,4,1:6) - ivvec(1:3,3,1:6)
  !
  ivvec(1:3,17,1:6) = -ivvec(1:3,1,1:6) + ivvec(1:3,2,1:6) + ivvec(1:3,4,1:6)
  ivvec(1:3,18,1:6) =  ivvec(1:3,1,1:6) - ivvec(1:3,2,1:6) + ivvec(1:3,3,1:6)
  ivvec(1:3,19,1:6) =  ivvec(1:3,2,1:6) - ivvec(1:3,3,1:6) + ivvec(1:3,4,1:6)
  ivvec(1:3,20,1:6) =  ivvec(1:3,1,1:6) + ivvec(1:3,3,1:6) - ivvec(1:3,4,1:6)
  !
  wlsm(1, 1:10) = real((/1440,    0,   30,    0, -38, -56, -38, -28,  7,     9/), 8)
  wlsm(2, 1:10) = real((/   0, 1440,    0,   30, -28,   9,   7, -38, -38,  -56/), 8)
  wlsm(3, 1:10) = real((/  30,    0, 1440,    0,  17, -46,  17,   7, -28,    9/), 8)
  wlsm(4, 1:10) = real((/   0,   30,    0, 1440,   7,   9, -28,  17,  17,  -46/), 8)
  !
  wlsm(1,11:20) = real((/ -46,   17,   17,  -28,   9,   7, -18, -18,  12,  -18/), 8)
  wlsm(2,11:20) = real((/   9,    7,  -28,   17, -46,  17, -18, -18, -18,   12/), 8)
  wlsm(3,11:20) = real((/ -56,  -38,  -38,    7,   9, -28,  12, -18, -18,  -18/), 8)
  wlsm(4,11:20) = real((/   9,  -28,    7,  -38, -56, -38, -18,  12, -18,  -18/), 8)
  !
  wlsm(1:4,1:20) = wlsm(1:4,1:20) / 1260
  !
  !wlsm(1:4,1:20) = 0.0d0
  !wlsm(1,1:4) = real((/1, 0, 0, 0/), 8)
  !wlsm(2,1:4) = real((/0, 1, 0, 0/), 8)
  !wlsm(3,1:4) = real((/0, 0, 1, 0/), 8)
  !wlsm(4,1:4) = real((/0, 0, 0, 1/), 8)
  !
end subroutine tetra_type
!
! Define step of energy
!
subroutine en_step()
  !
  use global, only : nk, nb, eig, e0, dos, ne, ivvec, grid, wlsm, ng, nwfc
  use omp_lib
  implicit none
  !
  integer :: ik, it, ii, jj, ib, ikv(3), ie
  real(8) :: de, emin, emax, ei(4,nb)
  !
  de = 0d0
  !
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP & SHARED(nk, nb, eig, e0, dos, ne, ivvec, grid, wlsm, ng, de) &
  !$OMP & PRIVATE(ik, it, ii, jj, ib, ikv, ie, emin, emax, ei)
  !
  !$OMP DO REDUCTION(+:de)
  do ik = 1, nk
     !
     do it = 1, 6
        !
        ei(1:4,1:nb) = 0.0d0
        !
        do ii = 1, 20
           !
           ikv(1:3) = grid(1:3,ik) + ivvec(1:3, ii, it)
           ikv(1:3) = modulo(ikv(1:3), ng(1:3)) + 1
           !
           do ib = 1, nb
              ei(1:4,ib) = ei(1:4,ib) &
              &          + wlsm(1:4,ii) * eig(ib,ikv(1),ikv(2),ikv(3))
           end do
           !
        end do
        !
        do ib = 1, nb
           do ii = 1, 4
              do jj = 1, ii - 1
                 de = de + abs(ei(ii,ib) - ei(jj,ib)) 
              end do
           end do
        end do
        !
     end do ! it
     !
  end do ! ik
  !$OMP END DO
  !
  !$OMP END PARALLEL
  !
  de = de / dble(nk * nb * 6 * 6)
  !
  emin = minval(eig(1:nb, 1: ng(1), 1:ng(2), 1:ng(3)))
  emax = maxval(eig(1:nb, 1: ng(1), 1:ng(2), 1:ng(3)))
  write(*,*) "  Emin[Ry] : ", emin
  write(*,*) "  Emax[Ry] : ", emax
  !
  ne = int((emax - emin) / de) + 1
  allocate(e0(ne), dos(nwfc + 1,ne))
  !
  do ie = 1, ne
     e0(ie) = emin + de * dble(ie - 1)
  end do
  ie = minloc(abs(e0(1:ne)), 1)
  e0(1:ne) = e0(1:ne) - e0(ie)
  !
  write(*,*) " Energy step [Ryd] : ", de
  write(*,*) " # of energy for DOS : ", ne
  !
end subroutine en_step
!
! Dos calclation
!
subroutine calc_dos()
  !
  use libtetra_mod, only : libtetra
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
  call libtetra("dos", 2, ng, ng, nb, ne, ef, nelec, bvec, e0, eig, eig, wdos)
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
  use libtetra_mod, only : libtetra
  use global, only : ng, nb, bvec, eig, nwfc, wfc
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
     call libtetra("occ", 2, ng, ng / n, nb, 1, ef, nelec, bvec, (/0d0/), eig, eig, wocc)
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
  use global, only : read_dat, calc_dos, write_dos, ng, nk, grid, en_step, tetra_type, calc_occ
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
  ! Define k-point grid
  !
  nk = product(ng(1:3))
  !
  allocate(grid(3,nk))
  ik = 0
  do i3 = 1, ng(3)
     do i2 = 1, ng(2)
        do i1 = 1, ng(1)
           !
           ik = ik + 1
           grid(1:3,ik) = (/i1, i2, i3 /) - 1
           !
        end do
     end do
  end do
  !
  call tetra_type()
  !
  call en_step()
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
