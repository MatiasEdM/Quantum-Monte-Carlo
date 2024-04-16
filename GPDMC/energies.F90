SUBROUTINE e_loc(NEL,NUC,a,Z,Rn,r,EL,KinL,UeeL,UeNL)


  implicit none

  integer, intent(in)           :: NEL
  integer, intent(in)           :: NUC
  integer, intent(in)           :: Z(NUC)
  double precision, intent(in)  :: a(NUC)
  double precision, intent(in)  :: Rn(3,NUC)
  double precision, intent(in)  :: r(3,NEL)


  double precision, external    :: kin
  double precision, external    :: Vee
  double precision, external    :: VeN

  double precision, intent(out) :: EL
  double precision, intent(out) :: KinL
  double precision, intent(out) :: UeeL
  double precision, intent(out) :: UeNL

  KinL = kin(NEL,NUC,a(:),Rn(:,:),r(:,:))
  UeeL = Vee(NEL,r(:,:))
  UeNL = VeN(NEL,NUC,Z(:),Rn(:,:),r(:,:))
  EL   = KinL + UeeL + UeNL 

  return

END SUBROUTINE e_loc


double precision FUNCTION VeN(NEL,NUC,Z,Rn,r)

  implicit none

  integer, intent(in)           :: NEL
  integer, intent(in)           :: NUC
  integer, intent(in)           :: Z(NUC)
  double precision, intent(in)  :: Rn(3,NUC)
  double precision, intent(in)  :: r(3,NEL)

  double precision, external    :: distance

  integer                       :: i,K
  double precision              :: dist_iK
  double precision, allocatable :: vec_iK(:)


  allocate( vec_iK(3) )

  VeN = 0.d0

  do K = 1,NUC,1
     do i = 1,NEL,1
        vec_iK(:) = r(:,i) - Rn(:,K) 
        dist_iK   = distance(3,vec_iK(:),vec_iK(:))
        if (dist_iK.gt.0.d0) then 
           VeN    = VeN - Z(K) * ( 1.d0/dist_iK )
        else
           write(*,*) 'WARNING: ELECTRONS TOO CLOSE TO NUCLEUS'
           VeN    = VeN - huge(1.d0)
        endif
     enddo
  enddo

  deallocate ( vec_iK )

  return

END FUNCTION VeN

double precision FUNCTION Vee(NEL,r)

  implicit none

  integer         , intent(in)  :: NEL
  double precision, intent(in)  :: r(3,NEL)
 
  double precision, external    :: distance

  integer                       :: i, j
  double precision, allocatable :: r_ij(:)
  double precision              :: dist_ij

  allocate( r_ij(3) )

  Vee = 0.d0

  do i = 1,NEL,1
     do j = 1,NEL,1
        if (i.eq.j) cycle
        r_ij    = r(:,j) - r(:,i)
        dist_ij = distance(3,r_ij(:),r_ij(:))
        if (dist_ij.gt.0.d0) then    
           Vee  = Vee + 1.d0/dist_ij
        else
           write(*,*) 'WARNING: ELECTRONS TOO CLOSE'
           Vee = Vee + huge(1.d0)
        endif
     enddo
  enddo

  Vee = 0.5d0 * Vee

  deallocate( r_ij )

  return

END FUNCTION Vee

double precision FUNCTION VNN(NUC,Z,Rn)

  implicit none

  integer         , intent(in)  :: NUC
  integer         , intent(in)  :: Z(NUC)
  double precision, intent(in)  :: Rn(3,NUC)
 
  double precision, external    :: distance

  integer                       :: A, B
  double precision, allocatable :: R_AB(:)
  double precision              :: dist_AB

  allocate( R_AB(3) )

  VNN = 0.d0

  do A = 1,NUC,1
     do B = 1,NUC,1
        if (A.eq.B) cycle
        R_AB    = Rn(:,B) - Rn(:,A)
        dist_AB = distance(3,R_AB(:),R_AB(:))
        VNN     = VNN + (Z(A)*Z(B))/dist_AB 
     enddo
  enddo

  VNN = 0.5d0 * VNN

  deallocate( R_AB )

  return

END FUNCTION VNN


double precision FUNCTION kin(NEL,NUC,a,Rn,r)

  implicit none

  integer         , intent(in) :: NEL
  integer         , intent(in) :: NUC
  double precision, intent(in) :: a(NUC)
  double precision, intent(in) :: Rn(3,NUC)
  double precision, intent(in) :: r(3,NEL)

  double precision, external   :: MO
  double precision, external   :: LAP_MO

  integer                      :: i
 
  kin = 0.d0

  do i = 1,NEL,1
     kin = kin - 0.5d0 * LAP_MO(NUC,a(:),Rn(:,:),r(:,i)) / MO(NUC,a(:),Rn(:,:),r(:,i))
  enddo

  return

END FUNCTION kin
