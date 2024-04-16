double precision FUNCTION AO(a,Rn,r)

  implicit none

  double precision, intent(in)  :: a
  double precision, intent(in)  :: Rn(3)
  double precision, intent(in)  :: r(3)

  double precision, external    :: distance

  double precision              :: dist
  double precision, allocatable :: vec(:)

  allocate( vec(3) )

  vec(:) = r(:)-Rn(:)

  dist   = distance(3,vec(:),vec(:))

  AO     = dexp(- a * dist )

  deallocate( vec )

  return

END FUNCTION AO


double precision FUNCTION LAP_AO(a,Rn,r)

  implicit none

  double precision, intent(in)  :: a
  double precision, intent(in)  :: Rn(3)
  double precision, intent(in)  :: r(3)

  double precision, external    :: distance
  double precision, external    :: AO

  double precision              :: dist
  double precision, allocatable :: vec(:)

  allocate( vec(3) )

  vec(:) = r(:) - Rn(:)

  dist   = distance(3,vec(:),vec(:))

  LAP_AO = (a**2 - (2.d0*a)/dist) * AO(a,Rn,r)

  return

END FUNCTION LAP_AO 


SUBROUTINE GRAD_AO(a,Rn,r,grad)

  implicit none

  double precision, intent(in)  :: a
  double precision, intent(in)  :: Rn(3)
  double precision, intent(in)  :: r(3)

  double precision              :: dist
  double precision, allocatable :: vec(:)

  double precision, external    :: distance
  double precision, external    :: AO

  double precision, intent(out) :: grad(3)

  allocate( vec(3) )

  vec(:)  = r(:) - Rn(:)

  dist    = distance(3,vec(:),vec(:))

  grad(:) = - a/dist * vec(:) * AO(a,Rn,r)

  deallocate( vec )

  return

END SUBROUTINE GRAD_AO
