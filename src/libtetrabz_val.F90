MODULE libtetrabz_val
  !
  INTEGER,SAVE :: &
  & fst,          & !
  & lst,          & !
  & nk,           & !
  & nk0,          & !
  & nt,           & !
  & nb,           & !
  & ne,           & !
  & ng(3),        & !
  & ltetra,       & !
  & ivvec(3,20,6)   !
  !
  REAL(8),SAVE :: &
  & wlsm(4,20)      !
  !
  INTEGER,SAVE,ALLOCATABLE:: &
  & indx1(:,:), &
  & indx2(:,:), &
  & indx3(:)
  !
END MODULE libtetrabz_val
