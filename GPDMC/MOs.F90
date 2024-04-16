double precision FUNCTION MO(NUC,a,Rn,r)

  implicit none

  integer         , intent(in) :: NUC
  double precision, intent(in) :: a(NUC)
  double precision, intent(in) :: Rn(3,NUC)
  double precision, intent(in) :: r(3)

  double precision, external   :: AO

  integer                      :: K

  MO   = 0.d0

  do K = 1,NUC,1
     MO = MO + AO(a(K),Rn(:,K),r(:))
  enddo

  return

END FUNCTION MO


double precision FUNCTION LAP_MO(NUC,a,Rn,r)

  implicit none

  integer         , intent(in) :: NUC
  double precision, intent(in) :: a(NUC)
  double precision, intent(in) :: Rn(3,NUC)
  double precision, intent(in) :: r(3)

  double precision, external   :: LAP_AO

  integer                      :: K

  LAP_MO = 0.d0

  do K = 1,NUC,1
     LAP_MO = LAP_MO + LAP_AO(a(K),Rn(:,K),r(:))
  enddo 

  return

END FUNCTION LAP_MO


SUBROUTINE GRAD_MO(NUC,a,Rn,r,grad)

  implicit none

  integer         , intent(in)  :: NUC
  double precision, intent(in)  :: a(NUC)
  double precision, intent(in)  :: Rn(3,NUC)
  double precision, intent(in)  :: r(3)

  integer                       :: K
  double precision, allocatable :: gr_ao(:)

  double precision, intent(out) :: grad(3)

  allocate( gr_ao(3) )

  grad(:) = 0.d0

  do K = 1,NUC,1
     call GRAD_AO(a(K),Rn(:,K),r(:),gr_ao(:))
     grad(:) = grad(:) + gr_ao(:)
  enddo
     
  deallocate( gr_ao )

  return

END SUBROUTINE GRAD_MO
