SUBROUTINE pdmc(NMC, NEL, NUC, dt, tau_max, E_ref, Rn, a, Z, energy, accep, E_path, Kin_path, Vee_path, VeN_path)

  implicit none
  
  integer         , intent(in)    :: NMC
  integer         , intent(in)    :: NUC
  integer         , intent(in)    :: NEL
  double precision, intent(in)    :: dt
  double precision, intent(in)    :: tau_max
  double precision, intent(in)    :: E_ref

  integer         , intent(in)    :: Z(NUC)
  double precision, intent(in)    :: a(NUC)
  double precision, intent(in)    :: Rn(3,NUC)

  integer                         :: i
  integer                         :: istep
  integer                         :: n_accep
  double precision                :: sq_dt
  double precision                :: psi_old, psi_new
  double precision, allocatable   :: chi(:)
  double precision, allocatable   :: r_old(:,:), r_new(:,:)
  double precision, allocatable   :: b_old(:,:), b_new(:,:)

  double precision                :: Wacc
  double precision                :: w
  double precision                :: tau
  double precision                :: EL
  double precision                :: normalization
  double precision                :: prod
  double precision                :: arg_expo
  double precision                :: wfc
  double precision                :: tr_prob
  double precision                :: ratio
  double precision                :: u

  double precision                :: KinL
  double precision                :: UeeL
  double precision                :: UeNL

  double precision                :: energy_kin
  double precision                :: energy_ee
  double precision                :: energy_eN

  double precision, external      :: psi

  double precision, intent(out)   :: energy, accep
  double precision, intent(inout) :: E_path(NMC)
  double precision, intent(inout) :: Kin_path(NMC)
  double precision, intent(inout) :: Vee_path(NMC)
  double precision, intent(inout) :: VeN_path(NMC)

!----------------!
! INITIALIZATION !
!----------------!

  allocate( chi(3) )
  allocate( r_old(3,NEL), r_new(3,NEL) )
  allocate( b_old(3,NEL), b_new(3,NEL) )

  energy        = 0.d0
  accep         = 0.d0
  Wacc          = 1.d0
  tau           = 0.d0
  normalization = 0.d0

  energy_kin    = 0.d0
  energy_ee     = 0.d0
  energy_eN     = 0.d0

  n_accep       = 0

  sq_dt         = dsqrt(dt)

  call random_seed()
  do i = 1,NEL,1
     call random_gauss(r_old(:,i),3)
  enddo
  do i = 1,NEL,1
     call drift(NUC,a(:),Rn(:,:),r_old(:,i),b_old(:,i))
  enddo
  psi_old = psi(NEL,NUC,a(:),Rn(:,:),r_old(:,:))

!--------------------------!
! MAIN QMC LOOP OVER PATHS !
!--------------------------!


  do istep = 1,NMC,1
     ![1]  Evaluate the local energy at r, El(r)

           call  e_loc(NEL,NUC,a(:),Z(:),Rn(:,:),r_old(:,:),EL,KinL,UeeL,UeNL)
     
     ![2]  Compute the contribution to the weight w(r) = exp{-dt*(El(r) - Eref)}

           w             = dexp(-dt*(EL-E_ref))

     ![3]  Update the cummulative weight W(r) = W(r)*w(r)

           Wacc          = Wacc * w

     ![4]  Accumulate the weighted energy W(r)*El(r) and the weight W(r) for the normalization

           energy          = energy + EL*Wacc

           energy_kin      = energy_kin + KinL*Wacc
           energy_ee       = energy_ee  + UeeL*Wacc
           energy_eN       = energy_eN  + UeNL*Wacc
           
           normalization   = normalization + Wacc
           
           E_path(istep)   = E_path(istep)   + energy/normalization
           Kin_path(istep) = Kin_path(istep) + energy_kin/normalization
           Vee_path(istep) = Vee_path(istep) + energy_ee/normalization
           VeN_path(istep) = VeN_path(istep) + energy_eN/normalization

     ![5]  Update the imaginary time tau = tau + dt

           tau = tau + dt

     ![6]  Check for a reset tau > tau_max

           if (tau.gt.tau_max) then
              Wacc = 1.d0
              tau  = 0.d0
           endif

     ![7]  Update the position r' = r + dt*b(r) + chi*sqrt(dt) [stochastic equation]
          
           do i = 1,NEL,1
              call random_gauss(chi,3)
              r_new(:,i) = r_old(:,i) + dt*b_old(:,i) + chi(:)*sq_dt
           enddo

     ![8]  Evaluate psi(r') and b(r') [new]

           do i = 1,NEL,1
              call drift(NUC,a(:),Rn(:,:),r_new(:,i),b_new(:,i))
           enddo
         
           psi_new  = psi(NEL,NUC,a(:),Rn(:,:),r_new(:,:))

     ![9]  Compute the ratio A = (T[r'->r]*P[r'])/(T[r->r']*P[r])

          tr_prob  = 1.d0
          do i = 1,NEL,1
             prod     = dot_product(b_new(:,i)+b_old(:,i),r_new(:,i)-r_old(:,i))
             arg_expo = prod + 0.5d0 * dt * (dot_product(b_new(:,i),b_new(:,i)) - dot_product(b_old(:,i),b_old(:,i)))
             tr_prob  = tr_prob * dexp(-arg_expo)
          enddo
          wfc      = psi_new/psi_old
          ratio    = wfc * wfc * tr_prob
                    
     ![10] Metropolis algorithm for acceptance/rejection of the move

          call random_number(u)

          if (ratio.ge.u) then
             n_accep    = n_accep + 1
             r_old(:,:) = r_new(:,:)
             b_old(:,:) = b_new(:,:)
             psi_old    = psi_new
          endif

  enddo

  energy = energy / normalization
  accep  = dble(n_accep) / dble(NMC)

  deallocate( chi )
  deallocate( r_old, r_new )
  deallocate( b_old, b_new )
  
  return

END SUBROUTINE pdmc
