SUBROUTINE read_input(Nruns, NMC, NEL, NUC, dt, tau_max, E_ref, Rn, a, Z)

  implicit none

  integer                      , intent(out) :: Nruns
  integer                      , intent(out) :: NMC
  integer                      , intent(out) :: NEL
  integer                      , intent(out) :: NUC
  double precision             , intent(out) :: dt
  double precision             , intent(out) :: tau_max
  double precision             , intent(out) :: E_ref
  double precision, allocatable, intent(out) :: Rn(:,:)
  double precision, allocatable, intent(out) :: a(:)
  integer         , allocatable, intent(out) :: Z(:)

  integer                                    :: K
  character(1)                               :: ch1
  character(10)                              :: ch2

  read(*,*) ch1
  read(*,*) ch2, Nruns
  read(*,*) ch2, NMC
  read(*,*) ch2, dt
  read(*,*) ch2, tau_max

  read(*,*)
 
  read(*,*) ch1
  read(*,*) ch2, NEL
 
  read(*,*)
 
  read(*,*) ch1
  read(*,*) ch2, NUC
  allocate( Z(NUC) , a(NUC), Rn(3,NUC))  
  do K = 1,NUC,1
     read(*,*) Z(K), a(K), Rn(:,K)
  enddo

  read(*,*)
 
  read(*,*) ch1
  read(*,*) ch2, E_ref

  return

END SUBROUTINE read_input


SUBROUTINE print_input(Nruns, NMC, NEL, NUC, dt, tau_max, E_ref, Rn, a, Z)

  implicit none

  integer         , intent(in) :: Nruns
  integer         , intent(in) :: NMC
  integer         , intent(in) :: NEL
  integer         , intent(in) :: NUC
  double precision, intent(in) :: dt
  double precision, intent(in) :: tau_max
  double precision, intent(in) :: E_ref
  double precision, intent(in) :: Rn(3,NUC)
  double precision, intent(in) :: a(NUC)
  integer         , intent(in) :: Z(NUC)

  integer                                    :: K

  write(*,*)
   
  write(*,*) '-------SIMULATION PARAMETERS----------'
  write(*,*)
  write(*,fmt=200) '  Number of runs          :', Nruns
  write(*,fmt=200) '  Number of MC iterations :', NMC
  write(*,fmt=210) '  Time step               :', dt
  write(*,fmt=210) '  Maximum projection time :', tau_max
  write(*,*) 
  write(*,*) '--------------------------------------'

  write(*,*)
 
  write(*,*) '--------------ELECTRONS---------------'
  write(*,*)
  write(*,fmt=200) '  Number of electrons    :', NEL
  write(*,*)
  write(*,*) '--------------------------------------'


  write(*,*)
 
  write(*,*) '---------------NUCLEI-----------------'
  write(*,*)
  write(*,fmt=200) '  Number of nuclei       :', NUC
  write(*,*)
  write(*,*) '  Z  |  a  |  Rx          Ry          Rz  [a.u.]'  
  write(*,*)
  do K = 1,NUC,1
     write(*,fmt=220) Z(K), a(K), Rn(:,K)
  enddo
  write(*,*)
  write(*,*) '--------------------------------------'


  write(*,*)

  write(*,*) '-------------ENERGIES-----------------'
  write(*,*)
  write(*,fmt=210) '  Reference energy       :', E_ref
  write(*,*)
  write(*,*) '--------------------------------------'


  200 format(A27,1I10)
  210 format(A27,1F10.5)
  220 format(1I6,1F6.2,3F10.3)

  return

END SUBROUTINE print_input


SUBROUTINE random_gauss(z,n)

  implicit none

  integer, intent(in)           :: n
  
  double precision              :: u(n+1)
  double precision, parameter   :: two_pi = 2.d0*dacos(-1.d0)
  integer                       :: i

  double precision, intent(out) :: z(n)

  call random_number(u)

  if (iand(n,1) == 0) then
     ! n is even
     do i=1,n,2
        z(i)   = dsqrt(-2.d0*dlog(u(i))) 
        z(i+1) = z(i) * dsin( two_pi*u(i+1) )
        z(i)   = z(i) * dcos( two_pi*u(i+1) )
     end do

  else
     ! n is odd
     do i=1,n-1,2
        z(i)   = dsqrt(-2.d0*dlog(u(i))) 
        z(i+1) = z(i) * dsin( two_pi*u(i+1) )
        z(i)   = z(i) * dcos( two_pi*u(i+1) )
     end do

     z(n)   = dsqrt(-2.d0*dlog(u(n))) 
     z(n)   = z(n) * dcos( two_pi*u(n+1) )

  end if

  return

END SUBROUTINE random_gauss


SUBROUTINE ave_error(a,m,ave,err)

  implicit none

  integer         , intent(in)  :: m
  double precision, intent(in)  :: a(m)
    
  double precision              :: variance
     
  double precision, intent(out) :: ave, err
      
  if (m.lt.1) then
     stop 'm < 1 not enough independent MC calculations'
     elseif (m.eq.1) then
           ave = a(1)
           err = 0.d0
     else
            ave      = sum(a(:))/dble(m)
            variance = sum((a(:) - ave)**2) / dble(m-1)
            err      = dsqrt(variance / dble(m) )
  endif

  return

END SUBROUTINE ave_error

double precision FUNCTION distance(n,a,b)

  implicit none

  integer, intent(in)          :: n
  double precision, intent(in) :: a(n)
  double precision, intent(in) :: b(n)

  distance = dsqrt( dot_product( a(:), b(:) ) )

  return

END FUNCTION distance
