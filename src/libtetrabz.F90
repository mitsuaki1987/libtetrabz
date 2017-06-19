MODULE libtetrabz_mpi_routines
  !
  IMPLICIT NONE
  !
  INTEGER,SAVE :: comm
  !
CONTAINS
!
! Initialize grid
!
SUBROUTINE libtetrabz_mpi_kgrid()
  !
  USE libtetrabz_vals, ONLY : nk, nk0, indx1, indx2, indx3, ng, ivvec, fst, lst
  !
  INTEGER :: it, i1, i2, i3, ii, ikv(3), nt, ik
  !
  ALLOCATE(indx1(20, 6 * nk), indx2(20, 6 * nk), indx3(20 * 6 * nk))
  !
  nt = 0
  DO i3 = 1, ng(3)
     DO i2  = 1, ng(2)
        DO i1 = 1, ng(1)
           !
           DO it = 1, 6
              !
              nt = nt + 1
              !
              DO ii = 1, 20
                 !
                 ikv(1:3) = (/i1, i2, i3/) + ivvec(1:3,ii,it) - 1
                 ikv(1:3) = MODULO(ikv(1:3), ng(1:3))
                 !
                 indx1(ii,nt) = 1 + ikv(1) + ng(1) * ikv(2) + ng(1) * ng(2) * ikv(3)
                 !
              END DO
              !
           END DO
           !
        END DO
     END DO
  END DO
  !
  indx2(1:20,1:6 * nk) = 0
  indx3(1:20 * 6 * nk) = 0
  !
  CALL libtetrabz_fst_and_lst()
  !
  nk0 = 0
  DO it = fst, lst
     !
     DO ii = 1, 20
        !
        DO ik = 1, nk0
           !
           IF(indx1(ii,it) == indx3(ik)) THEN
              !
              indx2(ii,it) = ik
              GOTO 10
              !
           END IF
           !
        END DO
        !
        nk0 = nk0 + 1
        indx2(ii,it) = nk0
        indx3(nk0) = indx1(ii,it)
        !
10      continue
        !
     END DO
     !
  END DO
  !
END SUBROUTINE libtetrabz_mpi_kgrid
!
! Compute cnt and dsp
!
SUBROUTINE libtetrabz_fst_and_lst()
  !
  USE mpi, ONLY : mpi_comm_size, mpi_comm_rank
  USE libtetrabz_vals, ONLY : fst, lst, nk
  !
  INTEGER :: ii, petot, my_rank, ierr
  INTEGER,ALLOCATABLE :: cnt(:), dsp(:)
  !
  CALL MPI_COMM_SIZE(comm, petot, ierr)
  CALL MPI_COMM_RANK(comm, my_rank, ierr)
  !
  ALLOCATE(cnt(0:petot-1), dsp(0:petot-1))
  !
  cnt(0:petot-1)        = 6 * nk / petot
  cnt(0:mod(6 * nk,petot)-1) = 6 * nk / petot + 1
  dsp(0) = 0
  DO ii = 1, petot - 1
     dsp(ii) = dsp(ii - 1) + cnt(ii - 1)
  END DO
  !
  fst = dsp(my_rank) + 1
  lst = dsp(my_rank) + cnt(my_rank)
  !
  DEALLOCATE(cnt,dsp)
  !
END SUBROUTINE libtetrabz_fst_and_lst
!
END MODULE libtetrabz_mpi_routines
!
!
!
MODULE libtetrabz_mpi
  !
  IMPLICIT NONE
  !
CONTAINS
!
END MODULE libtetrabz_mpi
