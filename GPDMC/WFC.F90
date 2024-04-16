double precision FUNCTION psi(NEL,NUC,a,Rn,r)

  implicit none

  integer         , intent(in) :: NEL
  integer         , intent(in) :: NUC
  double precision, intent(in) :: a(NUC)
  double precision, intent(in) :: Rn(3,NUC)
  double precision, intent(in) :: r(3,NEL)

  double precision, external   :: MO

  integer                      :: i

  psi = 1.d0

  do i = 1,NEL,1
     psi = psi * MO(NUC,a(:),Rn(:,:),r(:,i))
  enddo

  return

END FUNCTION psi
