PROGRAM QMC

  implicit none

!----------------------!
! VARIABLE DECLARATION !
!----------------------!

! PARAMETERS

  double precision, parameter   :: angs_to_au = 1.8897261254535d0 

! INPUT VARIABLES

  integer                       :: Nruns
  integer                       :: NMC
  integer                       :: NEL
  integer                       :: NUC
  double precision              :: dt
  double precision              :: tau_max
  double precision              :: E_ref
  double precision, allocatable :: Rn(:,:)
  double precision, allocatable :: a(:)
  integer         , allocatable :: Z(:)

         
! LOCAL VARIABLES

  integer                       :: iruns
  integer                       :: K
  double precision, allocatable :: E(:)
  double precision, allocatable :: accep(:)
  double precision              :: ave
  double precision              :: err
  character(1)                  :: ch1

  integer                       :: istep
  double precision, allocatable :: E_path(:)
  double precision, allocatable :: Kin_path(:)
  double precision, allocatable :: Vee_path(:)
  double precision, allocatable :: VeN_path(:)

  double precision              :: UNN

  double precision, external    :: VNN

!---------------!
! INPUT READING !
!---------------!

!  call read_input(Nruns, NMC, NEL, NUC, dt, tau_max, E_ref, Rn(:,:), a(:), Z(:))

  read(*,*) ch1
  read(*,*) Nruns
  read(*,*) NMC
  read(*,*) dt
  read(*,*) tau_max

  read(*,*)
 
  read(*,*) ch1
  read(*,*) NEL
 
  read(*,*)
 
  read(*,*) ch1
  read(*,*) NUC
  allocate( Z(NUC) , a(NUC), Rn(3,NUC))  
  do K = 1,NUC,1
     read(*,*) Z(K), a(K), Rn(:,K)
  enddo

  Rn(:,:) = Rn(:,:) * angs_to_au

  read(*,*)

  read(*,*) ch1
  read(*,*) E_ref


  allocate( E(Nruns), accep(Nruns) )
  allocate( E_path(NMC), Kin_path(NMC), Vee_path(NMC), VeN_path(NMC) )

!-------------!
! HELLO WORLD !
!-------------!

  write(*,*)
  write(*,*)'**************************************'
  write(*,*)'|          PDMC CALCULATION           |'
  write(*,*)'**************************************'
  write(*,*)

  call print_input(Nruns, NMC, NEL, NUC, dt, tau_max, E_ref, Rn(:,:), a(:), Z(:))
   
!---------------------!
! PDMC LOOP OVER RUNS !
!---------------------!

  E_path   = 0.d0
  Kin_path = 0.d0
  Vee_path = 0.d0
  VeN_path = 0.d0

  UNN      = VNN(NUC,Z(:),Rn(:,:))

  do iruns = 1,Nruns,1
     call pdmc(NMC, NEL, NUC, dt, tau_max, E_ref, Rn(:,:), a(:), Z(:), E(iruns), accep(iruns), E_path(:), Kin_path(:), &
               & Vee_path(:), VeN_path(:))
  enddo

!---------!
! RESULTS !
!---------!

  write(*,*)
  write(*,*)'**************************************'
  write(*,*)'|            PDMC RESULTS             |'
  write(*,*)'**************************************'
  write(*,*)


!----------------!
! QMC STATISTICS !
!----------------!
 
  call ave_error(E,Nruns,ave,err)
  write(*,*) '----------ELECTRONIC ENERGY-----------'
  write(*,*)
  write(*,fmt=200) '  E = ', ave, ' +/- ', err, ' [Hartree]'
  write(*,*)
  write(*,*) '-------ELECTRONIC ENERGY + VNN--------'
  write(*,*)
  write(*,fmt=200) '  E = ', ave + UNN, ' +/- ', err, ' [Hartree]'
  write(*,*)

!----------------------------!
! ACCEPTANCE RATE STATISTICS !
!----------------------------!

  call ave_error(accep,Nruns,ave,err)
  write(*,*) '----------ACCEPTANCE RATE-------------'
  write(*,*)
  write(*,fmt=210) '  A = ', ave, ' +/- ', err

!-----------------------!
! ENERGY ALONG THE PATH !
!-----------------------!

  open(unit=20,file='energy.out',access='sequential',status='replace')
 
  E_path(:)   = E_path(:)/dble(Nruns)
  Kin_path(:) = Kin_path(:)/dble(Nruns)
  Vee_path(:) = Vee_path(:)/dble(Nruns)
  VeN_path(:) = VeN_path(:)/dble(Nruns)

  do istep = 1,NMC,1
     write(20,*) istep*dt, E_path(istep), Kin_path(istep), Vee_path(istep), VeN_path(istep)
  enddo

  close(unit=20,status='keep')  



  200 format(A7,1F7.3,A5,1F7.3,A11)
  210 format(A7,1F7.3,A5,1F7.3)

  deallocate( Rn, a, Z )
  deallocate( E, accep )
  deallocate( E_path, Kin_path, Vee_path, VeN_path )

  stop

END PROGRAM QMC
