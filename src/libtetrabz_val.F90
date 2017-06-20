MODULE libtetrabz_val
  !
  LOGICAL,SAVE :: &
  & lmpi, &
  & linterpol
  !
  INTEGER,SAVE :: &
  & nkBZ, &
  & comm,         & !
  & nt_local,     & !
  & nk_local,     & !
  & nb,           & !
  & ne,           & !
  & ng(3)         !
  !
  REAL(8),SAVE :: &
  & wlsm(4,20)      !
  !
  INTEGER,SAVE,ALLOCATABLE :: &
  & ik_global(:,:), &
  & ik_local(:,:)
  !
  REAL(8),SAVE,ALLOCATABLE :: &
  & kvec(:,:)
  !
END MODULE libtetrabz_val
