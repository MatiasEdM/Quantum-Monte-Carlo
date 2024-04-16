SUBROUTINE  drift(NUC,a,Rn,r,b)

  implicit none

  integer         , intent(in)  :: NUC
  double precision, intent(in)  :: a(NUC)
  double precision, intent(in)  :: Rn(3,NUC)
  double precision, intent(in)  :: r(3)
 
  double precision, external    :: MO

  double precision, intent(out) :: b(3)
 
  call GRAD_MO(NUC,a(:),Rn(:,:),r(:),b(:))

  b(:) = b(:)/MO(NUC,a(:),Rn(:,:),r(:))

  return

END SUBROUTINE drift
