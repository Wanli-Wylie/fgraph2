! =================================================================================================================================
! MODULE       : stomate_somdynamics
!
! CONTACT      : orchidee-help _at_ ipsl.jussieu.fr
!
! LICENCE      : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF       Calculate soil dynamics largely following the Century model
!!	
!!\n DESCRIPTION: None
!!
!! RECENT CHANGE(S): None
!!
!! SVN		:
!! $HeadURL: svn://forge.ipsl.jussieu.fr/orchidee/perso/albert.jornet/ORCHIDEE-CN-P/src_stomate/stomate_soilcarbon.f90 $ 
!! $Date: 2017-01-10 13:20:28 +0100 (Tue, 10 Jan 2017) $
!! $Revision: 3985 $
!! \n
!_ ================================================================================================================================

MODULE stomate_som_dynamics

  ! modules used:

  USE ioipsl_para
  USE xios_orchidee
  USE stomate_data
  USE constantes
  USE constantes_soil
  USE stomate_phosphorus
 

  IMPLICIT NONE

  ! private & public routines

  PRIVATE
  PUBLIC som_dynamics,som_dynamics_clear,nitrogen_dynamics,nitrogen_dynamics_clear, &
         control_temp_func, biologicalN2fixation

  ! Variables shared by all subroutines in this module
  
  LOGICAL, SAVE    :: firstcall_som = .TRUE.   !! Is this the first call? (true/false)
!$OMP THREADPRIVATE(firstcall_som)
  LOGICAL, SAVE    :: firstcall_nitrogen = .TRUE.   !! Is this the first call? (true/false)
!$OMP THREADPRIVATE(firstcall_nitrogen)
  ! flux fractions within carbon pools
  REAL,ALLOCATABLE,SAVE, DIMENSION(:,:,:)          :: frac_carb        !! Flux fractions between carbon pools 
                                                                              !! (second index=origin, third index=destination) 
                                                                              !! (unitless, 0-1)
!$OMP THREADPRIVATE(frac_carb)
  REAL, ALLOCATABLE,SAVE, DIMENSION(:,:)           :: frac_resp        !! Flux fractions from carbon pools to the atmosphere (respiration) (unitless, 0-1)
!$OMP THREADPRIVATE(frac_resp)

CONTAINS


!! ================================================================================================================================
!!  SUBROUTINE   : som_dynamics_clear
!!
!>\BRIEF        Set the flag ::firstcall_som to .TRUE. and as such activate sections 1.1.2 and 1.2 of the subroutine som_dynamics 
!! (see below).
!! 
!_ ================================================================================================================================
  
  SUBROUTINE som_dynamics_clear
    firstcall_som=.TRUE.
  END SUBROUTINE som_dynamics_clear


!! ================================================================================================================================
!!  SUBROUTINE   : som_dynamics
!!
!>\BRIEF        Computes the soil respiration and carbon stocks, essentially 
!! following Parton et al. (1987).
!!
!! DESCRIPTION	: The soil is divided into 3 carbon pools, with different 
!! characteristic turnover times : active (1-5 years), slow (20-40 years) 
!! and passive (200-1500 years).\n
!! There are three types of carbon transferred in the soil:\n
!! - carbon input in active and slow pools from litter decomposition,\n
!! - carbon fluxes between the three pools,\n
!! - carbon losses from the pools to the atmosphere, i.e., soil respiration.\n
!!
!! The subroutine performs the following tasks:\n
!!
!! Section 1.\n
!! The flux fractions (f) between carbon pools are defined based on Parton et 
!! al. (1987). The fractions are constants, except for the flux fraction from
!! the active pool to the slow pool, which depends on the clay content,\n
!! \latexonly
!! \input{soilcarbon_eq1.tex}
!! \endlatexonly\n
!! In addition, to each pool is assigned a constant turnover time.\n
!!
!! Section 2.\n
!! The carbon input, calculated in the stomate_litter module, is added to the 
!! carbon stock of the different pools.\n
!!
!! Section 3.\n
!! First, the outgoing carbon flux of each pool is calculated. It is 
!! proportional to the product of the carbon stock and the ratio between the 
!! iteration time step and the residence time:\n
!! \latexonly
!! \input{soilcarbon_eq2.tex}
!! \endlatexonly
!! ,\n
!! Note that in the case of crops, the additional multiplicative factor 
!! integrates the faster decomposition due to tillage (following Gervois et 
!! al. (2008)).
!! In addition, the flux from the active pool depends on the clay content:\n
!! \latexonly
!! \input{soilcarbon_eq3.tex}
!! \endlatexonly
!! ,\n
!! Each pool is then cut from the carbon amount corresponding to each outgoing
!! flux:\n
!! \latexonly
!! \input{soilcarbon_eq4.tex}
!! \endlatexonly\n
!! Second, the flux fractions lost to the atmosphere is calculated in each pool
!! by subtracting from 1 the pool-to-pool flux fractions. The soil respiration 
!! is then the summed contribution of all the pools,\n
!! \latexonly
!! \input{soilcarbon_eq5.tex}
!! \endlatexonly\n
!! Finally, each carbon pool accumulates the contribution of the other pools:
!! \latexonly
!! \input{soilcarbon_eq6.tex}
!! \endlatexonly
!!
!! Section 4.\n
!! If the flag SPINUP_ANALYTIC is set to true, the matrix A is updated following
!! Lardy (2011).
!!
!! RECENT CHANGE(S): None
!! 
!! MAIN OUTPUTS VARIABLE(S): carbon, resp_hetero_soil
!!
!! REFERENCE(S)   :
!! - Parton, W.J., D.S. Schimel, C.V. Cole, and D.S. Ojima. 1987. Analysis of 
!! factors controlling soil organic matter levels in Great Plains grasslands. 
!! Soil Sci. Soc. Am. J., 51, 1173-1179.
!! - Gervois, S., P. Ciais, N. de Noblet-Ducoudre, N. Brisson, N. Vuichard, 
!! and N. Viovy (2008), Carbon and water balance of European croplands 
!! throughout the 20th century, Global Biogeochem. Cycles, 22, GB2022, 
!! doi:10.1029/2007GB003018.
!! - Lardy, R, et al., A new method to determine soil organic carbon equilibrium,
!! Environmental Modelling & Software (2011), doi:10.1016|j.envsoft.2011.05.016
!!
!! FLOWCHART    :
!! \latexonly
!! \includegraphics[scale=0.5]{soilcarbon_flowchart.jpg}
!! \endlatexonly
!! \n
!_ ================================================================================================================================

  SUBROUTINE som_dynamics (npts, clay, silt, &
       som_input, control_temp, control_moist, &! DSGclean drainage,&
       CN_target, CP_target,som, &
       resp_hetero_soil, &
       MatrixA,n_mineralisation, p_mineralisation, & 
       CN_som_litter_longterm, CP_som_litter_longterm, tau_CN_longterm)

!! 0. Variable and parameter declaration

    !! 0.1 Input variables
    
    INTEGER, INTENT(in)                            :: npts             !! Domain size (unitless)
    REAL, DIMENSION(npts), INTENT(in)              :: clay             !! Clay fraction (unitless, 0-1) 
    REAL, DIMENSION(npts), INTENT(in)              :: silt             !! Silt fraction (unitless, 0-1)
    REAL, DIMENSION(npts,ncarb,nvm,nelements), INTENT(in) :: som_input !! Amount of Organic Matter going into the SOM pools from litter decomposition \f$(gC m^{-2} day^{-1})$\f
    REAL, DIMENSION(npts,nlevs), INTENT(in)        :: control_temp     !! Temperature control of heterotrophic respiration (unitless: 0->1)
    REAL, DIMENSION(npts,nlevs), INTENT(in)        :: control_moist    !! Moisture control of heterotrophic respiration (unitless: 0.25->1)
    !DSGclean REAL, DIMENSION(npts), INTENT(in)              :: drainage    ! water lost from the soil column by leaching  [ kg /m2 /day]
                                                                                   ! (fraction of water content)
    REAL, DIMENSION(npts,nvm,ncarb), INTENT(in)    :: CN_target        !! C to N ratio of SOM flux from one pool to another (gN m-2 dt-1)
    REAL, DIMENSION(npts,nvm,ncarb), INTENT(in)    :: CP_target        !! C to P ratio of SOM flux from one pool to another (g(C) g-1(P))
    !! 0.2 Output variables
    
    REAL, DIMENSION(npts,nvm), INTENT(out)         :: resp_hetero_soil !! Soil heterotrophic respiration \f$(gC m^{-2} (dt_sechiba one_day^{-1})^{-1})$\f

    !! 0.3 Modified variables
    
    REAL, DIMENSION(npts,ncarb,nvm,nelements), INTENT(inout) :: som             !! SOM pools: active, slow, or passive, \f$(gC m^{2})$\f
    REAL, DIMENSION(npts,nvm,nbpools,nbpools), INTENT(inout) :: MatrixA  !! Matrix containing the fluxes between the carbon pools
                                                                                !! per sechiba time step 
                                                                                !! @tex $(gC.m^2.day^{-1})$ @endtex
    REAL, DIMENSION(npts,nvm), INTENT(inout)                 :: n_mineralisation !! Release of mineral N from SOM decomposition  (gN m-2 dt-1)
    REAL, DIMENSION(npts,nvm), INTENT(inout)                 :: p_mineralisation !! Release of mineral P from SOM decomposition  (gP m-2 dt-1)
    REAL, DIMENSION(npts,nvm,nbpools), INTENT(inout)         :: CN_som_litter_longterm !! Longterm CN ratio of litter and som pools (gC/gN)
    REAL, DIMENSION(npts,nvm,nbpools), INTENT(inout)         :: CP_som_litter_longterm !! Longterm CP ratio of litter and som pools (gC/gP)
    REAL, INTENT(inout)                                      :: tau_CN_longterm        !! Counter used for calculating the longterm CN ratio of som and litter pools [seconds]

    !! 0.4 Local variables
    REAL                                           :: dt               !! Time step \f$(dt_sechiba one_day^{-1})$\f
    LOGICAL                                               :: l_error          !! Diagnostic boolean for error allocation (true/false)
    INTEGER                                        :: ier              !! Check errors in netcdf call
    REAL, SAVE, DIMENSION(ncarb)                   :: som_turn         !! Residence time in SOM pools (days)
!$OMP THREADPRIVATE(som_turn)
    REAL, DIMENSION(npts,ncarb,nelements)          :: fluxtot          !! Total flux out of carbon pools \f$(gC m^{2})$\f
    REAL, DIMENSION(npts,ncarb,ncarb,nelements)    :: flux             !! Fluxes between carbon pools \f$(gC m^{2})$\f
    REAL, DIMENSION(npts,ncarb,ncarb,nvm,nelements)    :: som_f        !!Fluxes between carbon pools \f$(gC m^{2})$\f
    CHARACTER(LEN=7), DIMENSION(ncarb)                    :: soilpools_str    !! Name of the soil pools for informative outputs (unitless)
    ! mineral nitrogen in the soil (gN/m**2)
    INTEGER                                        :: k,kk,m,j,l      !! Indices (unitless)
    CHARACTER(LEN=2), DIMENSION(nelements)                     :: element_str !! string suffix indicating element
!_ ================================================================================================================================

!! 1. Initializations
    !! printlev is the level of diagnostic information, 0 (none) to 4 (full)
    IF (printlev>=3) WRITE(numout,*) 'Entering som_dynamics' 

    dt = dt_sechiba/one_day
    IF ( firstcall_som ) THEN
        l_error = .FALSE.
        ALLOCATE (frac_carb(npts,ncarb,ncarb), stat=ier)
        l_error = l_error .OR. (ier.NE.0)
        ALLOCATE (frac_resp(npts,ncarb), stat=ier)
        l_error = l_error .OR. (ier.NE.0)
        IF (l_error) THEN
           STOP 'stomate_som_dynamics: error in memory allocation'
        ENDIF
 
        frac_carb(:,:,:)   = zero
        !! 1.1 Get soil "constants"
        !! 1.1.1 Flux fractions between carbon pools: depend on clay content, recalculated each time
        ! From active pool: depends on clay content
        frac_carb(:,iactive,ipassive) = active_to_pass_ref_frac + active_to_pass_clay_frac*clay(:)
        ! JC MOD042 not sure why use (clay+silt) here, original version is clay only
        ! frac_carb(:,iactive,islow)    = un - frac_carb(:,iactive,ipassive) - (active_to_co2_ref_frac - active_to_co2_clay_silt_frac*(clay(:)+silt(:)))
        frac_carb(:,iactive,islow)    = un - frac_carb(:,iactive,ipassive) - (active_to_co2_ref_frac - active_to_co2_clay_silt_frac*clay(:)) 
   
        ! From slow pool 
        ! JC MOD042 trunk version is 0.03 default   
        ! frac_carb(:,islow,ipassive) = slow_to_pass_ref_frac + slow_to_pass_clay_frac*clay(:)
        ! OCN doesn't use Parton 1993 formulation for frac_carb(:,islow,ipassive) 
        ! but the one of 1987 : ie = 0.03 .... 
        ! frac_carb(:,islow,iactive) = un - frac_carb(:,islow,ipassive) - slow_to_co2_ref_frac
        frac_carb(:,islow,ipassive) = frac_carb_sp
        frac_carb(:,islow,iactive) = frac_carb_sa
    
        ! From passive pool
        frac_carb(:,ipassive,iactive) = pass_to_active_ref_frac
        frac_carb(:,ipassive,islow)   = pass_to_slow_ref_frac

        ! From surface pool
        ! JC MOD042 use trunk version parameter here assuming surface turnover the same
        ! as active turnover
        ! frac_carb(:,isurface,islow) = surf_to_slow_ref_frac
        frac_carb(:,isurface,ipassive) = frac_carb(:,iactive,ipassive)
        frac_carb(:,isurface,islow) = frac_carb(:,iactive,islow)

        !! 1.1.2 Determine the respiration fraction : what's left after
        ! subtracting all the 'pool-to-pool' flux fractions
        ! Diagonal elements of frac_carb are zero
        frac_resp(:,:) = un - frac_carb(:,:,isurface) - frac_carb(:,:,iactive) - frac_carb(:,:,islow) - &
             frac_carb(:,:,ipassive) 

        !! 1.1.3 Turnover in SOM pools (in days)
        !! som_turn_ipool are the turnover (in year)
        !! It is weighted by Temp and Humidity function later
        som_turn(iactive)  = som_turn_iactive  / one_year    
        som_turn(islow)    = som_turn_islow    / one_year            
        som_turn(ipassive) = som_turn_ipassive / one_year    
        som_turn(isurface) = som_turn_isurface / one_year
        
        !! 1.2 Messages : display the residence times  
        soilpools_str(iactive)  = 'active'
        soilpools_str(islow)    = 'slow'
        soilpools_str(ipassive) = 'passive'
        soilpools_str(isurface) = 'surface'

        WRITE(numout,*) 'som_dynamics:'
        
        WRITE(numout,*) '   > minimal SOM residence time in soil pools (d):'
        DO k = 1, ncarb ! Loop over soil pools
          WRITE(numout,*) '(1, ::soilpools_str(k)):',soilpools_str(k),' : (1, ::som_turn(k)):',som_turn(k)
          WRITE(numout,*) 'FRACCARB k=',k,' ',frac_carb(1,k,isurface),frac_carb(1,k,iactive),frac_carb(1,k,islow),frac_carb(1,k,ipassive)
        ENDDO
        
        WRITE(numout,*) '   > flux fractions between soil pools: depend on clay content'
        
        firstcall_som = .FALSE.
        
    ENDIF

    !! 1.3 Set soil respiration to zero
    resp_hetero_soil(:,:) = zero
    som_f(:,:,:,:,:) = zero
!! 2. Update the SOM stocks with the different soil carbon input

    som(:,:,:,:) = som(:,:,:,:) + som_input(:,:,:,:) * dt

!! 3. Fluxes between carbon reservoirs, and to the atmosphere (respiration) \n

    !! 3.2. Calculate fluxes

    DO m = 2,nvm ! Loop over # PFTs

      !! 3.2.1. Flux out of pools

      DO k = 1, ncarb ! Loop over SOM pools from which the flux comes (active, slow, passive)
        
        DO l = 1, nelements ! Loop over elements (Carbon, Nitrogen, Phosphorus)
           ! Determine total flux out of pool
           fluxtot(:,k,l) = dt*som_turn(k) * som(:,k,m,l) * &
                control_moist(:,ibelow) * control_temp(:,ibelow) * decomp_factor(m) 

           ! Flux from active pools depends on clay content
           IF ( k .EQ. iactive ) THEN
              fluxtot(:,k,l) = fluxtot(:,k,l) * ( un - som_turn_iactive_clay_frac * clay(:) )
           ENDIF
     
           ! Update the loss in each carbon pool
           som(:,k,m,l) = som(:,k,m,l) - fluxtot(:,k,l)
        ENDDO

        ! Fluxes towards the other pools (k -> kk)
        DO kk = 1, ncarb ! Loop over the SOM pools where the flux goes
          ! Carbon flux
          flux(:,k,kk,icarbon)   = frac_carb(:,k,kk) * fluxtot(:,k,icarbon)
          ! Nitrogen flux - Function of the C stock of the 'departure' pool 
          ! and of the C to N target ratio of the 'arrival' pool
          flux(:,k,kk,initrogen) = frac_carb(:,k,kk) * fluxtot(:,k,icarbon) / & 
               CN_target(:,m,kk)
          ! Phosphorus flux - Function of the N stock of the 'departure' pool 
          ! and of the rNto P target ratio of the 'arrival' pool
          flux(:,k,kk,iphosphorus) = frac_carb(:,k,kk) * fluxtot(:,k,icarbon) / & 
                CP_target(:,m,kk)
        ENDDO
        
     ENDDO ! End of loop over SOM pools

      
      !! 3.2.2 respiration
      !BE CAREFUL: Here resp_hetero_soil is divided by dt to have a value which corresponds to
      ! the sechiba time step but then in stomate.f90 resp_hetero_soil is multiplied by dt.
      ! Perhaps it could be simplified. Moreover, we must totally adapt the routines to the dtradia/one_day
      ! time step and avoid some constructions that could create bug during future developments.
      !
      
      resp_hetero_soil(:,m) = &
         ( frac_resp(:,iactive)  * fluxtot(:,iactive ,icarbon) + &
           frac_resp(:,islow)    * fluxtot(:,islow   ,icarbon) + &
           frac_resp(:,ipassive) * fluxtot(:,ipassive,icarbon) + &
           frac_resp(:,isurface) * fluxtot(:,isurface,icarbon)  ) / dt
      !! save fluxes for output 
      som_f(:,:,:,m,:) = flux(:,:,:,:) 
      !! 3.2.3 add fluxes to active, slow, and passive pools
      
      DO k = 1, ncarb ! Loop over SOM pools

         DO l = 1, nelements ! Loop over elements (Carbon, Nitrogen, Phosphorus)
            som(:,k,m,l)   = som(:,k,m,l)          + &
                              flux(:,iactive ,k,l) + &
                              flux(:,ipassive,k,l) + &
                              flux(:,islow   ,k,l) + &
                              flux(:,isurface,k,l)
         ENDDO

      !! 3.2.4 compute mineralization
         n_mineralisation(:,m) = n_mineralisation(:,m) + fluxtot(:,k,initrogen) - &
               (flux(:,k,iactive,initrogen)+ flux(:,k,ipassive,initrogen)       + &
                flux(:,k,islow,initrogen)  + flux(:,k,isurface,initrogen))

         p_mineralisation(:,m) = p_mineralisation(:,m) + fluxtot(:,k,iphosphorus) - &
                (flux(:,k,iactive,iphosphorus) + flux(:,k,ipassive,iphosphorus)   + &
                 flux(:,k,islow,iphosphorus)   + flux(:,k,isurface,iphosphorus))

     ENDDO ! Loop over SOM pools
      
   ENDDO ! End loop over PFTs

   ! For Yilong
   ! resp_hetero_soil has been divided by dt above
    CALL xios_orchidee_send_field("HETERO_RESP_SOIL",resp_hetero_soil(:,:))

    DO l=1,nelements
       IF     (l == icarbon) THEN
          element_str(l) = '_c'
       ELSEIF (l == initrogen) THEN
          element_str(l) = '_n'
       ELSEIF (l == iphosphorus) THEN
          element_str(l) = '_p'
       ELSE
          STOP 'Define element_str'
       ENDIF
       CALL xios_orchidee_send_field("SOMF_ACTIVE_PASSIVE"//TRIM(element_str(l)),som_f(:,iactive,ipassive,:,l)/dt)
       CALL xios_orchidee_send_field("SOMF_ACTIVE_SLOW"//TRIM(element_str(l)),som_f(:,iactive,islow,:,l)/dt)
       CALL xios_orchidee_send_field("SOMF_SLOW_PASSIVE"//TRIM(element_str(l)),som_f(:,islow,ipassive,:,l)/dt)
       CALL xios_orchidee_send_field("SOMF_SLOW_ACTIVE"//TRIM(element_str(l)),som_f(:,islow,iactive,:,l)/dt)
       CALL xios_orchidee_send_field("SOMF_PASSIVE_ACTIVE"//TRIM(element_str(l)),som_f(:,ipassive,iactive,:,l)/dt)
       CALL xios_orchidee_send_field("SOMF_PASSIVE_SLOW"//TRIM(element_str(l)),som_f(:,ipassive,islow,:,l)/dt)
       CALL xios_orchidee_send_field("SOMF_SURF_PASSIVE"//TRIM(element_str(l)),som_f(:,isurface,ipassive,:,l)/dt)
       CALL xios_orchidee_send_field("SOMF_SURF_SLOW"//TRIM(element_str(l)),som_f(:,isurface,islow,:,l)/dt)
    ENDDO
    
 !! 4. (Quasi-)Analytical Spin-up
    
    !! 4.1.1 Finish to fill MatrixA with fluxes between soil pools
    
    IF (spinup_analytic) THEN

       DO m = 2,nvm 

          ! flux leaving the active pool
          MatrixA(:,m,iactive_pool,iactive_pool) = moins_un * &
               dt*som_turn(iactive) * &
               control_moist(:,ibelow) * control_temp(:,ibelow) * &
               ( 1. - som_turn_iactive_clay_frac * clay(:)) * decomp_factor(m)

          ! flux received by the active pool from the slow pool
          MatrixA(:,m,iactive_pool,islow_pool) =  frac_carb(:,islow,iactive)*dt*som_turn(islow) * &
               control_moist(:,ibelow) * control_temp(:,ibelow) * decomp_factor(m)

          ! flux received by the active pool from the passive pool
          MatrixA(:,m,iactive_pool,ipassive_pool) =  frac_carb(:,ipassive,iactive)*dt*som_turn(ipassive) * &
               control_moist(:,ibelow) * control_temp(:,ibelow) * decomp_factor(m)

          ! flux leaving the slow pool
          MatrixA(:,m,islow_pool,islow_pool) = moins_un * &
               dt*som_turn(islow) * &
               control_moist(:,ibelow) * control_temp(:,ibelow) * decomp_factor(m)

          ! flux received by the slow pool from the active pool
          MatrixA(:,m,islow_pool,iactive_pool) =  frac_carb(:,iactive,islow) *&
               dt*som_turn(iactive) * &
               control_moist(:,ibelow) * control_temp(:,ibelow) * &
               ( 1. - som_turn_iactive_clay_frac * clay(:) ) * decomp_factor(m)

          ! flux received by the slow pool from the surface pool
          MatrixA(:,m,islow_pool,isurface_pool) =  frac_carb(:,isurface,islow) *&
               dt*som_turn(isurface) * &
               control_moist(:,ibelow) * control_temp(:,ibelow) * decomp_factor(m)

          ! flux leaving the passive pool
          MatrixA(:,m,ipassive_pool,ipassive_pool) =  moins_un * &
               dt*som_turn(ipassive) * &
               control_moist(:,ibelow) * control_temp(:,ibelow) * decomp_factor(m)      

          ! flux received by the passive pool from the active pool
          MatrixA(:,m,ipassive_pool,iactive_pool) =  frac_carb(:,iactive,ipassive)* &
               dt*som_turn(iactive) * &
               control_moist(:,ibelow) * control_temp(:,ibelow) *&
               ( 1. - som_turn_iactive_clay_frac * clay(:) ) * decomp_factor(m)

          ! flux received by the passive pool from the slow pool
          MatrixA(:,m,ipassive_pool,islow_pool) =  frac_carb(:,islow,ipassive) * &
               dt*som_turn(islow) * &
               control_moist(:,ibelow) * control_temp(:,ibelow) * decomp_factor(m)

          ! flux leaving the surface pool
          MatrixA(:,m,isurface_pool,isurface_pool) =  moins_un * &
               dt*som_turn(isurface) * &
               control_moist(:,ibelow) * control_temp(:,ibelow) * decomp_factor(m)      
          
          ! Long term stoichiometry
          ! C-to-N 
          ! see comments in stomate_litter for the rationale behin these changes


          !WHERE (som(:,isurface,m,initrogen) .GT. min_stomate)
          WHERE (som(:,isurface,m,icarbon) .GT. min_stomate)
             CN_som_litter_longterm(:,m,isurface_pool) = ( CN_som_litter_longterm(:,m,isurface_pool) * (tau_CN_longterm-dt) &
                  + som(:,isurface,m,icarbon)/som(:,isurface,m,initrogen) * dt)/ (tau_CN_longterm)
             CP_som_litter_longterm(:,m,isurface_pool) = ( CP_som_litter_longterm(:,m,isurface_pool) * (tau_CN_longterm-dt) &
                  + som(:,isurface,m,icarbon)/som(:,isurface,m,iphosphorus) * dt)/ (tau_CN_longterm)
          ELSEWHERE
             CN_som_litter_longterm(:,m,isurface_pool) = zero
             CP_som_litter_longterm(:,m,isurface_pool) = zero
          ENDWHERE

          !WHERE (som(:,iactive,m,initrogen) .GT. min_stomate)
          WHERE (som(:,iactive,m,icarbon) .GT. min_stomate)
             CN_som_litter_longterm(:,m,iactive_pool)  = ( CN_som_litter_longterm(:,m,iactive_pool) * (tau_CN_longterm-dt) &
                  + som(:,iactive,m,icarbon)/som(:,iactive,m,initrogen) * dt)/ (tau_CN_longterm)
             CP_som_litter_longterm(:,m,iactive_pool)  = ( CP_som_litter_longterm(:,m,iactive_pool) * (tau_CN_longterm-dt) &
                  + som(:,iactive,m,icarbon)/som(:,iactive,m,iphosphorus) * dt)/ (tau_CN_longterm)
          ELSEWHERE
             CN_som_litter_longterm(:,m,iactive_pool)  = zero
             CP_som_litter_longterm(:,m,iactive_pool)  = zero
          ENDWHERE

          !WHERE(som(:,islow,m,initrogen) .GT. min_stomate)
          WHERE(som(:,islow,m,icarbon) .GT. min_stomate)
             CN_som_litter_longterm(:,m,islow_pool)    = ( CN_som_litter_longterm(:,m,islow_pool) * (tau_CN_longterm-dt) &
               + som(:,islow,m,icarbon)/som(:,islow,m,initrogen) * dt)/ (tau_CN_longterm)
             CP_som_litter_longterm(:,m,islow_pool)    = ( CP_som_litter_longterm(:,m,islow_pool) * (tau_CN_longterm-dt) &
               + som(:,islow,m,icarbon)/som(:,islow,m,iphosphorus) * dt)/ (tau_CN_longterm)
          ELSEWHERE
             CN_som_litter_longterm(:,m,islow_pool)    =  zero
             CP_som_litter_longterm(:,m,islow_pool)    =  zero
          ENDWHERE

          !WHERE(som(:,ipassive,m,initrogen) .GT. min_stomate)
          WHERE(som(:,ipassive,m,icarbon) .GT. min_stomate)
             CN_som_litter_longterm(:,m,ipassive_pool) =  ( CN_som_litter_longterm(:,m,ipassive_pool) * (tau_CN_longterm-dt) &
                  + som(:,ipassive,m,icarbon)/som(:,ipassive,m,initrogen) * dt)/ (tau_CN_longterm)
             CP_som_litter_longterm(:,m,ipassive_pool) =  ( CP_som_litter_longterm(:,m,ipassive_pool) * (tau_CN_longterm-dt) &
                  + som(:,ipassive,m,icarbon)/som(:,ipassive,m,iphosphorus) * dt)/ (tau_CN_longterm)
          ELSEWHERE
             CN_som_litter_longterm(:,m,ipassive_pool) =   zero
             CP_som_litter_longterm(:,m,ipassive_pool) =   zero
          ENDWHERE
          ! C-to-P 

          IF (printlev>=4) WRITE(numout,*)'Finish to fill MatrixA'

       ENDDO ! Loop over # PFTS

       ! 4.2 Add Identity for each submatrix(7,7) 

       DO j = 1,nbpools
          MatrixA(:,:,j,j) = MatrixA(:,:,j,j) + un 
       ENDDO

    ENDIF ! (spinup_analytic)

    IF (printlev>=4) WRITE(numout,*) 'Leaving som_dynamics'
    
  END SUBROUTINE som_dynamics


!! ================================================================================================================================
!!  SUBROUTINE   : nitrogen_dynamics_clear
!!
!>\BRIEF        Set the flag ::firstcall to .TRUE. 
!! 
!! 
!_ ================================================================================================================================
  
  SUBROUTINE nitrogen_dynamics_clear
    firstcall_nitrogen=.TRUE.
  END SUBROUTINE nitrogen_dynamics_clear



  ! This is essentially inspired by DNDC, but very simplified to avoid having
  ! to calculate microbe growth for the moment. Builds on the physico-chemical
  ! reactions with some fixed assumptions.
  SUBROUTINE nitrogen_dynamics(npts, contfrac, clay, sand, &
                               temp_sol, tmc_pft, drainage_pft, swc_pft,veget_max, resp_sol, &
                               som, biomass, n_input, pH, &
                               npp_week,                 &
                               mineralisation, pb, max_eau_var, plant_uptake, Bd,soil_n_min,p_O2,bact, &
                               grm_nfert, agri_nfert, fclover, BNF_clover)
    !
    ! 0 declarations
    !

    ! 0.1 input

    INTEGER, INTENT(in)                                   :: npts     !! Domain size
    REAL,DIMENSION (npts), INTENT (in)                    :: contfrac !! Fraction of continent in the grid cell (unitless)
    REAL, DIMENSION(npts), INTENT(in)                     :: clay     !! clay fraction (between 0 and 1)
    REAL, DIMENSION(npts), INTENT(in)                     :: sand     !! sand fraction (between 0 and 1)
    REAL, DIMENSION(npts), INTENT(in)                     :: temp_sol !! soil temperature (degC)
    REAL, DIMENSION(npts,nvm), INTENT(in)                 :: swc_pft  !! Relative soil water content per PFT
                                                                             
    REAL, DIMENSION(npts,nvm), INTENT(in)                 :: drainage_pft    !!  water lost from the soil column by leaching (drainage & runoff)
                                                                             !! [kg/m2/timestep] 
    REAL, DIMENSION(npts,nvm), INTENT(in)                 :: tmc_pft  !!  total water in the soil column 
                                                                             !! [kg/m2]
    REAL, DIMENSION(npts,nvm), INTENT(in)                 :: veget_max!! fraction of a vegetation (0-1)
    REAL, DIMENSION(npts,nvm), INTENT(in)                 :: resp_sol !! carbon respired from below ground 
                                                                             !! (hetero+autotrophic) (gC/m**2/day)
    REAL, DIMENSION(npts,ncarb,nvm,nelements), INTENT(in) :: som      !! SOM (gC(or N)/m**2)
    REAL, DIMENSION(npts,nvm,nparts,nelements), INTENT(in):: biomass  !! Biomass pools (gC(or N)/m**2)
    ! ninput=4
    REAL, DIMENSION(npts,nvm,ninput), INTENT(inout)       :: n_input  !! nitrogen inputs into the soil  (gN/m**2/day)
                                                                             !!  NH4 and NOX from the atmosphere, NH4 from BNF,
                                                                             !!   agricultural fertiliser as NH4/NO3 
    REAL, DIMENSION(npts), INTENT(in)                     :: pH       !! soil pH
    REAL, DIMENSION(npts,nvm)  ,INTENT(in)                :: npp_week !! "weekly" NPP (gC/day/(m**2 covered)
    REAL, DIMENSION(npts,nvm), INTENT(in)                 :: mineralisation !! net nitrogen mineralisation of decomposing SOM
                                                                                   !! (gN/m**2/day), supposed to be NH4
    REAL, DIMENSION(npts), INTENT(in)                     :: pb             !! Air pressure (hPa)
    REAL, DIMENSION(npts), INTENT(in)                     :: Bd             !! Bulk density (kg/m**3)
    REAL, DIMENSION(npts), INTENT(in)                     :: max_eau_var    !! Maximum water content of the soil   
    REAL, DIMENSION(npts,nvm,ninput), INTENT(in)          :: grm_nfert      !! nfert from GRM module
    REAL, DIMENSION(npts,nvm,ninput), INTENT(in)          :: agri_nfert     !! nfert from agricultural application (read from input maps)
    REAL, DIMENSION(npts,nvm), INTENT(in)                 :: fclover        !! fraction of clover in grass/legume mixture

    REAL, DIMENSION(npts,nvm,nnspec),INTENT(inout)        :: soil_n_min     !! mineral nitrogen in the soil (gN/m**2) 
                                                                                   !! (first index=npts, second index=nvm, third index=nnspec)  
    REAL, DIMENSION(npts,nvm),INTENT(inout)               :: p_O2           !! partial pressure of oxigen in the soil (hPa)
                                                                                   !! (first index=npts, second index=nvm)
                      
    REAL, DIMENSION(npts,nvm),INTENT(inout)               :: bact           !! denitrifier biomass (gC/m**2)
                                                                                   !! (first index=npts, second index=nvm)


    ! 0.2 output

    REAL, DIMENSION(npts,nvm,nionspec), INTENT(out)       :: plant_uptake   !! Uptake of soil N by platns 
                                                                                         !! (gN/m**2/timestep) 
    ! BNF by clover (legume) nodules (assuming directly into plant rather than
    ! soil), it is different from the bnf fixed by soil bacteria
    REAL, DIMENSION(npts,nvm) , INTENT(out)               :: BNF_clover     !! see above

    ! 0.3 local
    REAL                                                  :: dt             !! Time step \f$(dt_sechiba one_day^{-1})$\f
    LOGICAL                                                      :: l_error        !! Diagnostic boolean for error allocation (true/false)
    INTEGER                                               :: ier            !! Check errors in netcdf call
    REAL, DIMENSION(npts,nvm,nionspec)                    :: leaching       !! mineral nitrogen leached from the soil 
    REAL, DIMENSION(npts,nvm,nionspec)                    :: root_uptake    !! Uptake of soil N per root biomass
    REAL, DIMENSION(npts,nvm)                             :: immob          !! N immobilized (gN/m**2/day)

    REAL, DIMENSION(npts)                                 :: afps_max       !! maximum pore-volume of the soil
                                                                                   !! Table 2, Li et al., 2000
                                                                                   !! (fraction)                       
    REAL, DIMENSION(npts,nvm)                             :: afps           !! pore-volume of the soil
                                                                                   !! Table 2, Li et al., 2000
                                                                                   !! (fraction)   
    REAL, DIMENSION(npts,nvm)                             :: D_s            !! Oxygen diffusion in soil (m**2/day)
    REAL, DIMENSION(npts,nvm)                             :: mol_O2_resp    !! Density of Moles of O2 related to respiration term 
                                                                                   !! ((molesO2 m-3) (gC m-2)-1)
    REAL, DIMENSION(npts,nvm)                             :: p_O2_resp      !! O2 partial pressure related to the respiration term 
                                                                                   !! ((hPa) (gC m-2)-1)
    REAL, DIMENSION(npts)                                 :: p_O2air        !! Oxygen partial pressure in air (hPa)
    REAL, DIMENSION(npts)                                 :: g_O2           !! gradient of O2 partial pressure (hPa m-1)
    REAL, DIMENSION(npts)                                 :: d_O2           !! change in O2 partial pressure (hPa day-1)
    REAL, DIMENSION(npts,nvm)                             :: anvf           !! Volumetric fraction of anaerobic microsites 
                                                                                   !! (fraction of pore-volume)
    REAL, DIMENSION(npts)                                 :: FixNH4         !! Fraction of adsorbed NH4+ (-)
    REAL, DIMENSION(npts,nvm)                             :: n_adsorbed     !! Ammonium adsorpted (gN/m**2)
                                                                                   !! based on Li et al. 1992, JGR, Table 4
    REAL, DIMENSION(npts,nvm)                             :: fw             !! Effect of soil moisture on nitrification (-)
                                                                                   !! Zhang et al. 2002, Ecological Modelling, appendix A, page 101
   REAL, DIMENSION(npts)                 :: var_temp_sol  !! Temperature function used for calc. ft_nit (-)
                                                                 !! Zhang et al. 2002, Ecological Modelling, appendix A, page 101
    REAL, DIMENSION(npts)                :: ft_nit        !! Effect of temperature on nitrification (-)
                                                                 !! Zhang et al. 2002, Ecological Modelling, appendix A, page 101
    REAL, DIMENSION(npts)                :: fph           !! Effect of pH on nitrification (-)
                                                                 !! Zhang et al. 2002, Ecological Modelling, appendix A, page 101
    REAL, DIMENSION(npts)                :: ftv           !! Effect of temperature on NO2 or NO production 
                                                                 !! during nitrification (-)
                                                                 !! Zhang et al. 2002, Ecological Modelling, appendix A, page 102
    REAL, DIMENSION(npts,nvm,3)          :: nitrification !! N-compounds production (NO3, N2O, NO) related 
                                                                 !! to nitrification process (gN/m**2/tstep)
    REAL, DIMENSION(npts)                :: ft_denit      !! Temperature response of relative growth rate of 
                                                                 !! total denitrifiers (-)
                                                                 !! Eq. 2 Table 4 of Li et al., 2000
    REAL, DIMENSION(npts)                :: fph_no3       !! Soil pH response of relative growth rate of 
                                                                 !! NO3 denitrifiers (-)
                                                                 !! Eq. 2 Table 4 of Li et al., 2000
    REAL, DIMENSION(npts)                :: fph_no        !! Soil pH response of relative growth rate of 
                                                                 !! NO denitrifiers (-)
                                                                 !! Eq. 2 Table 4 of Li et al., 2000
    REAL, DIMENSION(npts)                :: fph_n2o       !! Soil pH response of relative growth rate of 
                                                                 !! N2O denitrifiers (-)
                                                                 !! Eq. 2 Table 4 of Li et al., 2000

    REAL, DIMENSION(npts)                :: fph_no3_ref   !! Soil pH response of NO3 denitrification at reference pH
    REAL, DIMENSION(npts)                :: fph_no_ref    !! Soil pH response of NO denitrification at reference pH
    REAL, DIMENSION(npts)                :: fph_n2o_ref   !! Soil pH response of N2O denitrification at reference pH
                                                                 
                                                                 

    REAL, DIMENSION(npts)                :: Kn_conv         !! Conversion from kgN/m3 to gN/m3(pore space) !DSG
    REAL, DIMENSION(npts)                :: mu_NO3          !! Relative growth rate of NO3 denitrifiers (hour**-1)
                                                                   !! Eq.1 Table 4 Li et al., 2000 
    REAL, DIMENSION(npts)                :: mu_N2O          !! Relative growth rate of N2O denitrifiers (hour**-1)
                                                                   !! Eq.1 Table 4 Li et al., 2000
    REAL, DIMENSION(npts)                :: mu_NO           !! Relative growth rate of NO denitrifiers (hour**-1)
                                                                   !! Eq.1 Table 4 Li et al., 2000  
    REAL, DIMENSION(npts)                :: sum_n           !! sum of all N species in the soil (gN/m**2)
    REAL, DIMENSION(npts,nvm,3)          :: denitrification !! N-compounds consumption (NO3, N2O, NO) related 
                                                                   !! to denitrificaiton process (gN/m**2/tstep)
    REAL, DIMENSION(npts)                :: dn_bact         !! denitrifier biomass change (kgC/m**3/timestep)
    REAL, DIMENSION(npts)                :: bact_min        !! minimum denitrifier biomass (kgC/m**3)
    REAL, DIMENSION(npts)                :: bact_max        !! minimum denitrifier biomass (kgC/m**3)

    REAL, DIMENSION(npts)                :: ft_uptake       !! Temperature response of N uptake by plants (-)
    REAL                                 :: conv_fac_vmax   !! Conversion from (umol (gDW)-1 h-1) to (gN (gC)-1 timestep-1)
    REAL, DIMENSION(nvm)                 :: nc_leaf_min     !! Minimal NC ratio of leaf (gN / gC)
    REAL, DIMENSION(nvm)                 :: nc_leaf_max     !! Maximal NC ratio of leaf (gN / gC)
    REAL, DIMENSION(npts)                :: lab_n           !! Labile nitrogen in plants (gN/m**2) 
    REAL, DIMENSION(npts)                :: lab_c           !! Labile carbon in plants (gC/m**2) 
    REAL, DIMENSION(npts)                :: NCplant         !! NC ratio of the plant (gN / gC)
                                                                   !! Eq. (9) p. 3 of SM of Zaehle & Friend, 2010
    REAL, DIMENSION(npts,nvm)            :: f_NCplant       !!  Response of Nitrogen uptake by plants 
                                                                   !! to N/C ratio of the labile pool

    REAL, DIMENSION(npts,nvm)            :: f_PNplant       !!  Response of BNF 
                                                                   !! to P/N ratio of the labile pool
    REAL, DIMENSION(npts,nvm)            :: f_Nmin          !! Response of BNF 
                                                                   !! to avail. soilmineral N
                                                                  

    REAL                                 :: conv_fac_concent!! Conversion factor from (umol per litter) to (gN m-2) 
    REAL, DIMENSION(npts)                :: frac_nh3        !! dissociation of [NH3] to [NH4+] (-)
    REAL, DIMENSION(npts)                :: emm_fac         !! Factor for reducing NH3 emissions (-)
    REAL, DIMENSION(npts)                :: F_clay          !! Response of N-emissions to clay fraction (-)
    REAL, DIMENSION(npts,nvm,nnspec)     :: emission        !! volatile losses of nitrogen 
                                                                   !! (gN/m**2/timestep)
    INTEGER                              :: m               !! Index to loop over the nvm PFT's
    INTEGER                              :: n               !! Index to loop over the npts grid boxes
    INTEGER                              :: ji              !! Index to loop over the npts grid boxes

    REAL, DIMENSION(npts)                :: qc              !! scalar for denitrification rates

    REAL, DIMENSION(npts,nvm)            :: f_drain         !! fraction of tmc which has been drained in the last time step 

    ! DSG mass conservation
    REAL, DIMENSION(npts,nvm)    :: mass_before             !! Temporary variable
    REAL, DIMENSION(npts,nvm)    :: mass_change             !! Temporary variable
    REAL, DIMENSION(npts,nvm)    :: mass_after              !! Temporary variable
    REAL, DIMENSION(npts)        :: tempdiff                !! Temporary variable

    LOGICAL, DIMENSION(npts,nvm)        :: N_added                 !! Nitrogen added during spinup_analytic to avoid overdepletion [true/false]
    REAL, DIMENSION(npts,nvm)    :: N_added_cons            !! Nitrogen added during spinup_cnp to avoid overdepletion [gN/m2/dt]
    REAL, PARAMETER              :: pH_ref=6.5              !! reference pH for pH response of dentrification ( Zaehle et al.,2011)


    REAL, PARAMETER              :: BNF_1bnfclovermin = 0.15!! clover N uptake reduction
    REAL, DIMENSION(npts)        :: vartmp                  !! Temporary variable 
!===============================================================================

!===============================================================================
    

    !! bavard is the level of diagnostic information, 0 (none) to 4 (full)
    IF (printlev>=3) WRITE(numout,*) 'Entering nitrogen_dynamics' 
    IF(printlev>=4)THEN
       WRITE(numout,*) 'CHECK values in nitrogen dynamics'
       WRITE(numout,*) 'soil_n_min ',soil_n_min(test_grid,test_pft,:)
    ENDIF

    dt = dt_sechiba/one_day

    IF ( firstcall_nitrogen ) THEN

        firstcall_nitrogen = .FALSE.

    ENDIF


    !DSG mass conservation ========================================
     mass_before(:,:) = SUM(soil_n_min(:,:,:),DIM=3)

    ! Initializations
    !DSGseq: to easily adjust the sequence of N consuming processes, we need to
    !initialize them with zero
    leaching(:,:,:)      = zero 
    nitrification(:,:,:) = zero
    denitrification(:,:,:) = zero
    plant_uptake(:,:,:)  = zero
    BNF_clover(:,:) = zero
    !DSGseq:

!=======================================
    ! 
    ! 1. Immobilisation of mass from decomposition
    !
    ! immobilisation has absolute priority to avoid mass conservation problems
    ! the code in litter and soilcarbon has to make sure that soil ammonium can never be 
    ! more than exhausted completely by immobilisation!!!
    ! 
    immob(:,:) = zero
    N_added_cons(:,:) = zero
    WHERE(mineralisation(:,:).LT.0.)
       immob(:,:) = - mineralisation(:,:)
       ! In case, the N related to immobilisation is higher that the N
       ! available in the [NH4+] pool, we take the remaining N from the 
       ! [NO3-] pool
       ! JC Debug there is possibility that mineralisation is negative
       ! and its abs value larger than available soil_n_min
       ! nitrate pool can be depleted
       soil_n_min(:,:,initrate) = soil_n_min(:,:,initrate) - &
            MAX(0.0,immob(:,:)-(soil_n_min(:,:,iammonium)-min_stomate))

       soil_n_min(:,:,iammonium) = soil_n_min(:,:,iammonium) - &
            MIN(immob(:,:),soil_n_min(:,:,iammonium)-min_stomate)
    ENDWHERE
    !! give warning for negative mineralisation that can not be deal as immb
    !IF (ANY(soil_n_min(:,:,initrate) .LT. zero)) THEN
    !  WRITE (numout,*) 'WARNING nitrogen_dynamics',&
    !       'negative mineralisation that can not be deal as immb, &
    !                   due to already low mineral N. &
    !        We will add N to avoid overdepletion, only when spinup_cnp or &
    !        impose_cn'
    !ENDIF
    !DSG: we need to avoid overdepletion of mineral N 
    ! (a) during the analytic spinup
    ! (b) when we impose CN ratios
    IF (spinup_cnp.OR.impose_cn) THEN 
      IF(printlev>=4)THEN
          WRITE(numout,*) 'during the spinup period we add N in case of overdepletion'
      ENDIF

      N_added(:,:) = .FALSE.
      DO m=1,nvm 
         WHERE ((soil_n_min(:,m,initrate) .LT. zero)) 
            N_added(:,m)             = .TRUE.
            N_added_cons(:,m) = - soil_n_min(:,m,initrate)
            soil_n_min(:,m,initrate) = zero
         ENDWHERE
      ENDDO

      IF(printlev>=4)THEN
          IF (ANY(N_added)) THEN 
              WRITE(numout,*) 'we added N to soil_n_min pools to avoid overdepletion'
              WRITE(numout,*) 'in : ', COUNT(N_added)
              WRITE(numout,*) 'of (nvm*npts): ', npts*nvm
          ENDIF
      ENDIF

    ENDIF

    IF(printlev>=4)THEN
       WRITE(numout,*) 'CHECK values after mineralisation'
       WRITE(numout,*) 'soil_n_min ',soil_n_min(test_grid,test_pft,:)
    ENDIF



!=======================================
    ! 
    ! 2. Plant uptake of available ionic forms of nitrogen 
    !
    ! DSG: To not artifically penalyze plant relative to microbe by the sequence of
    ! consumption processes we put them directly after immobilization
    !
    ! Comment from OCN : Temperature function of uptake is similar to SOM decomposition
    ! to avoid N accumulation at low temperatures 
    ! In addition in the SM of Zaehle & Friend, 2010 P. 3 it is mentioned
    ! "The uptake rate is observed to be sensitive to root temperature which is included as f(T), 
    ! thereby following the temperature sensitivity of net N mineralization"
    tempdiff = temp_sol(:)+tp_00
    ft_uptake(:) = control_temp_func (npts, tempdiff)

    !We suppress this T control on root uptake due to lack of evidence of such a
    !control:
    ft_uptake(:) = un
  
  
    
    ! Vmax of nitrogen uptake (in umol (g DryWeight_root)-1 h-1)

    !! calculate the scaling function according to eq. 10 , p. 3 of SM of Zaehle & Friend, 2010
    !! here we calculate the response function for the N/C ratio of labile plant biomass
    CALL calc_f_XYplant(npts,nvm,  &
                        biomass,   &
                        initrogen, & ! = X
                        icarbon,   & ! = Y
                        f_NCplant  &
                         )

    !JC MOD040 for crops no effect of plant NC on N uptake
    DO m = 2, nvm
      IF (.NOT. natural(m)) THEN
        f_NCplant(:,m) = un
      ENDIF
    ENDDO

    !! calculate root uptake of ammonium & nitrate
    DO m =iammonium,initrate
       CALL root_conductivity(npts,nvm,                      &  
                              soil_n_min(:,:,m),max_eau_var(:),     & 
                              vmax_N_uptake(m), K_N_min(m),            &
                              low_K_N_min(m), 14.,                &
                              dt,                            &
                              root_uptake(:,:,m))
       !DSG no effect of plant NC on uptake 
       !WRITE (numout,*) 'no effect of plant NC on N uptake'
       !f_NCplant(:,:) = un

       ! JCADD 20170427 BNF for grassland
       IF (GRM_allow_BNF) THEN
         ! If account for BNF from PaSim
         ! The mineral N uptake of clover is reduced by 15% (bnfclover = 0.15)
         ! as compared to the uptake by the grass. (Schmid et al., 2001 eq 12)
         plant_uptake(:,:,m) = root_uptake(:,:,m)  * biomass(:,:,iroot,icarbon) * &
            SPREAD(ft_uptake(:),NCOPIES=nvm,DIM=2 )* f_NCplant(:,:) * &
            (un - fclover * BNF_1bnfclovermin)

       ELSE
         ! multiply plant mass with the uptake capacity of roots and scale with N demand and temperature
         plant_uptake(:,:,m) = root_uptake(:,:,m)  * biomass(:,:,iroot,icarbon) *SPREAD(ft_uptake(:),NCOPIES=nvm,DIM=2 )* f_NCplant(:,:)
       ENDIF
       ! END JCADD

       ! DSG starve them 
       ! plant_uptake(:,:,m) = zero
       ! WRITE(6,*) 'we starve the plants'
    ENDDO

    IF(printlev>=4) THEN
       WRITE(numout,*) 'f_NCplant=',f_NCplant(test_grid,test_pft)
       WRITE(numout,*) 'root_uptake=',root_uptake(test_grid,test_pft,:)
       WRITE(numout,*) 'low_K_N_min',low_K_N_min
       WRITE(numout,*) 'K_N_min',K_N_min
       WRITE(numout,*) 'vmax_N_uptake',vmax_N_uptake
       WRITE(numout,*) 'max_eau_var(test_grid)',max_eau_var(test_grid)
       WRITE(numout,*) 'biomass(test_grid,test_pft,iroot,icarbon)',biomass(test_grid,test_pft,iroot,icarbon)
       WRITE(numout,*) 'plant_uptake(1)=',plant_uptake(test_grid,test_pft,:)
       WRITE(numout,*) 'denitrification(test_grid,test_pft,1)',denitrification(test_grid,test_pft,1)
       WRITE(numout,*) 'nitrification(test_grid,test_pft,1)',nitrification(test_grid,test_pft,1)
       WRITE(numout,*) 'leaching(test_grid,test_pft,initrate)',leaching(test_grid,test_pft,:)
    ENDIF

    !! to ensure mass conservation when changing the sequence of mineral N
    !! consuming processes I wrote account for all other consuming processes
    !! even if there are not yet computed here:
    plant_uptake(:,:,iammonium) = MAX(MIN(soil_n_min(:,:,iammonium) - & 
          (leaching(:,:,iammonium) + nitrification(:,:,i_nh4_to_no3) &
            + nitrification(:,:,i_nh4_to_no) + nitrification(:,:,i_nh4_to_n2o)), &  
           plant_uptake(:,:,iammonium)),zero)
    plant_uptake(:,:,initrate) = MAX(MIN(soil_n_min(:,:,initrate) - & 
          leaching(:,:,initrate) - denitrification(:,:,i_no3_to_n2) &
             +  nitrification(:,:,i_nh4_to_no3), plant_uptake(:,:,initrate)),zero)

    ! JCADD 20170427 BNF for grassland
    ! after calculation of plant_uptake
    CALL calc_BNF_clover(npts,nvm,dt, &
                soil_n_min,fclover,plant_uptake, &
                temp_sol,biomass,BNF_clover)

    ! 
    ! 3. Adsorption of ammonium
    !
    ! Ammonium adsorption, reduce ammonium concentrations by this ammount
    ! before computing the not yet computed N consuming processes.
    ! based on Li et al. 1992, JGR, Table 4
    ! Bd : Bulk density kg m-3
    ! FixNH4 : Fraction of adsorbed NH4+
    ! FixNH4=[0.41-0.47 log(NH4) clay/clay_max]
    ! NH4+ concentration in the soil liquid, gN kg-1 soil (p. 9774 of Li et al., 1992)
    ! but Zhang et al. (2000) use the same equation but NH4+ is defined as
    ! NH4+ in a soil layer in kgN ha-1
    !DSG! In OCN, NH4 seems defined as kgN kg-1 or m-3 of water 
    ! In OCN, NH4 is defined as mol(NH4) m-3 of water 
    !DSG! 
    ! Comment of N. Vuichard: I don't know which definition is the good one ...
    !  ... DSG: to be cosistent with root uptake, I take here N per volume
    !  of soil pores; 
    ! DSGminN001: OCN uses molar concentration (18 for NH4)
    !             I do the same
    
    ! a_FixNH4 = 0.41
    ! b_FixNH4 = -0.47
    ! clay_max = 0.63
    n_adsorbed(:,:) = zero
    DO m = 1, nvm
       WHERE((soil_n_min(:,m,iammonium).GT.min_stomate).AND.(tmc_pft(:,m).GT.min_stomate))
          FixNH4(:) = a_FixNH4 + b_FixNH4 * &
               !DSG log10(soil_n_min(:,m,iammonium) / ( zmaxh * Bd(:) ) ) &
               !DSGminN001 log10(soil_n_min(:,m,iammonium) / ( zmaxh * max_eau_var ) ) & !DSG: max_eau_var is not yet accounted for in zmaxh at this point of the code!!
               log10(soil_n_min(:,m,iammonium) / ( tmc_pft(:,m)/1000. * 18. ) ) & !DSG: like done in OCN, we use here molar concentration in soil water [m3/m2]
               !DSGminN001
               * MIN(clay(:),clay_max) / clay_max
       ELSEWHERE
          FixNH4(:) = zero
       ENDWHERE
       !DSGdebug: avoid negative values:
       FixNH4(:)=MAX(zero,FixNH4(:)) 

       ! In OCN, we do not multiply by FixNH4 but by FixNH4/(1+FixNH4)
       ! Comment of N. Vuichard: It is not clear if we should keep the formulation of OCN or not
       ! So far, we keep the original formulation
       ! DSG: I go with the OCN formulation 
       n_adsorbed(:,m) = MIN(soil_n_min(:,m,iammonium),soil_n_min(:,m,iammonium) * (FixNH4(:)/(un+FixNH4)))
       ! correction for low temperatures (from OCN)
       WHERE (temp_sol.LT.zero) 
          n_adsorbed(:,m) = n_adsorbed(:,m) * 0.5
       ENDWHERE
    ENDDO
    
    ! update soil mineral ammonium pool
    ! JC comment I think the current update is to make sure adsorbed ammonium
    ! will not go out of the system through leaching and emission
    ! it will be added back at the end of the subroutine
    soil_n_min(:,:,iammonium) = soil_n_min(:,:,iammonium) - n_adsorbed(:,:) 

    IF(printlev>=4)THEN
       WRITE(numout,*) 'CHECK values after adsorption'
       WRITE(numout,*) 'soil_n_min ',soil_n_min(test_grid,test_pft,:)
       WRITE(numout,*) 'Bd(test_grid)',Bd(test_grid)
       WRITE(numout,*) 'max_eau_var(test_grid)',max_eau_var(test_grid)
       WRITE(numout,*) 'tmc_pft(test_grid,test_pft)',tmc_pft(test_grid,test_pft)
       WRITE(numout,*) 'clay(test_grid)',clay(test_grid)
       WRITE(numout,*) 'clay_max', clay_max
    ENDIF


!=======================================
    ! 
    ! 4. Losses of nitrogen by leaching
    !
    WHERE((tmc_pft(:,:)-drainage_pft(:,:)) .NE. zero) 
        f_drain(:,:) = drainage_pft(:,:)/(tmc_pft(:,:)-drainage_pft(:,:))
    ELSEWHERE
        f_drain(:,:) = zero
    ENDWHERE

    DO m = 1, nvm
          leaching(:,m,iammonium) = MAX(MIN(soil_n_min(:,m,iammonium)* f_drain(:,m), &
                                             soil_n_min(:,m,iammonium) - plant_uptake(:,m,iammonium)), zero)
          leaching(:,m,initrate)  = MAX(MIN(soil_n_min(:,m,initrate)* f_drain(:,m), &
                                             soil_n_min(:,m,initrate) - plant_uptake(:,m,initrate) ),zero)
    ENDDO
    
!=======================================
    !  5.0 ANAEROBIC MICROSITE: this information is need to have simultaneous denitrification and
    ! nitrification
 
    ! APPROACH: calculate the anaerobic microsites as in OCN (Zaehle et al. 2011) 
    ! this is simplified  version of the DNDC approach which increases stability of the caluculation

    ! DSG_MICT_MERGE: oxygen diffusion from MICT could be used here to define the
    ! anaerobic fraction of the soil.

    DO m = 1, nvm
       ! function made up by SZ assuming a critical moisture threshold of 0.8 to have
       ! strongly varying anox behaviour at air entry
       anvf(:,m) = MIN(.9,MAX(0.1,(un-exp(-(swc_pft(:,m)/scal_anvf)**(scal_anvf*10)))))

       ! correction for low temperatures
       WHERE(temp_sol(:).LE.-0.001)
            anvf(:,m) = MIN(.98,MAX(zero,(un-exp(-(MAX(zero,swc_pft(:,m)-0.5)/scal_anvf)**(scal_anvf*10)))))
       ENDWHERE

       ! oxygen diffusion coefficient only used to calculate gasoues loss rates (emissions)
       ! not for micorsites (like NV)
       D_s(:,m) =  MAX(0.005,( un - swc_pft(:,m) )) / ( dt / z_decomp )

       ! correction for low temperatures
       WHERE(temp_sol(:).LE.-46.02+min_stomate)
          D_s(:,m)= zero 
       ENDWHERE

    ENDDO

!=======================================
       ! 5.1 Nitrification of NH4 to NO3 in the oxygenated part of the soil  
       ! OCN: simplified scheme based on Zhang et al. 2002 using thresholds
       ! DSG: makes more sense than keeping it zero (see above)
    WHERE(swc_pft(:,:).GT.0.02)
          fw(:,:)=un
    ELSEWHERE 
          fw(:,:)=zero
    ENDWHERE
  
    ! Effect of temperature on nitrification
    ! OCN: simplified scheme based on Zhang et al. 2002, with subzero degree
    ! nitrification rates based on Xi Ru  
    ! DSG: smooth curve over broad range of temperatures:

!DSG_hotpot
!DSG: the optimum function is not defined for very high or low temperatures:
!     thus we set the function to zero outside the range -40:75 C
!     the limits are given by the equation 
    WHERE(temp_sol(:).LE.-39.0) ! too cold
          ft_nit(:) = zero
    ELSEWHERE(temp_sol(:).GT.75.0) ! too hot
          ft_nit(:) = zero
    ELSEWHERE ! in this range the function is defined
          ft_nit(:) =( ((70. - temp_sol(:))/(70.-38.) )**12 ) * exp(12.*(temp_sol(:)-38.0)/(70.-38.))
    ENDWHERE
!DSG_hotpot
  
    ! Effect of pH on nitrification
    ! Zhang et al. 2002, Ecological Modelling, appendix A, page 101
    ! fph_0 = -1.2314
    ! fph_1 = 0.7347
    ! fph_2 = -0.0604
    !    fph(:) = fph_0 + fph_1 * pH(:) + fph_2 * ph(:)**2
    !JC revise a scalling factor was added in OCN
    fph(:) = fph_0 + fph_1 * (pH(:)+scal_ph) + fph_2 * (pH(:)+scal_ph)**2
    fph(:) = MAX(0.0, fph(:) )
  
    ! Effect of temperature on NO2 or NO production during nitrification
    ! Zhang et al. 2002, Ecological Modelling, appendix A, page 102
    ! ftv_0 = 2.72
    ! ftv_1 = 34.6
    ! ftv_2 = 9615.
!DSG_hotpot
    ! DSG:  exponential increase with T: here we need to set limits as we
    ! cannot expect the process to increase forever with warming.
    ! du to lack of information I set it to 35C
    WHERE(temp_sol(:).LT.35.) 
        ftv(:) = ftv_0 **(ftv_1 - ftv_2 /(temp_sol(:) + tp_00) )
    ELSEWHERE
        ftv(:) = ftv_0 **(ftv_1 - ftv_2 /(35. + tp_00) )
    ENDWHERE
!DSG_hotpot
    ftv(:) = MAX(0.0, ftv(:) )
  
    DO m = 1, nvm
       ! i_nh4_to_no3 = 1
       ! i_nh4_to_n2o = 3
       ! i_nh4_to_no  = 2
       !
       ! 5.1.2 Actual nitrification rate - NH4 to NO3
       ! 
       ! Equation for the nitrification rate probably from Schmid et al., 2001, Nutr. Cycl. Agro (eq.1)
       ! but the environmental factors used are from Zhang (2002) who used an other equation
       ! I don't know how this can be mixed together ?
       ! In addition, the formulation mixed the one from Schmid with the use of the anaerobic balloon defined by Li et al., 2000. I don't know how this can be mixed together ?
       ! Last, in OCN, the default fraction of N-NH4 which is converted to N-NO3 appears equal to 2 day-1 (a factor "2" in the equation) - In Schmid et al., the nitrification rate at 20 C and field capacity (knitrif,20) was set to 0.2 d−1 (Ten time less...)
       ! k_nitrif = 0.2
!JC revise this function follows Zaehle et al., 2011 with combined temperature
!and ph factors, while no water factor applied k_nitrif=0.2
! in OCN k_nitrif = 1.2 is chosen for achieve the 0.1 per day at 20degC at
! intermediate soil ph and moisture levels
       nitrification(:,m,i_nh4_to_no3) = MIN(fw(:,m) * fph(:) * ft_nit(:) * k_nitrif * dt * & 
            soil_n_min(:,m,iammonium) * (1.0 - anvf(:,m)) ,soil_n_min(:,m,iammonium)- leaching(:,m,iammonium)- plant_uptake(:,m,iammonium))

       IF(printlev>=4 .AND. m == test_pft)THEN
          WRITE(numout,*) 'CHECK values after nitrification nh4 to no3'
          WRITE(numout,*) 'PFT=',m
          WRITE(numout,*) 'nitrification(:,m,i_nh4_to_no3) ',nitrification(test_grid,m,i_nh4_to_no3)
          WRITE(numout,*) 'fw=',fw(test_grid,test_pft)
          WRITE(numout,*) 'fph=',fph(test_grid)
          WRITE(numout,*) 'ft_nit=',ft_nit(test_grid)
          WRITE(numout,*) 'anvf=',anvf(test_grid,m)
       ENDIF
      
       !
       ! 5.1.3 Emission of N2O during nitrification - NH4 to N2O
       ! 
       ! From Zhang et al., 2002 - Appendix A p. 102
       ! Reference n2o production per N-NO3 produced g N-N2O  (g N-NO3)-1
       ! n2o_nitrif_p = 0.0006
 
       !JC revise using Zhang et al., 2002  
       ! revise again using XuRi's emission factor RN2ON < 0.1-0.2% day-1
       !         nitrification(:,m,i_nh4_to_n2o) = RN2ON * nitrification(:,m,i_nh4_to_no3)
       ! 0.0008 is chosen by OCN to achieve 0.2% at ~30degC?
       nitrification(:,m,i_nh4_to_n2o) = ftv(:) * 0.0008 * nitrification(:,m,i_nh4_to_no3)

       IF(printlev>=4 .AND. m==test_pft)THEN
          WRITE(numout,*) 'CHECK values after nitrification nh4 to n2o'
          WRITE(numout,*) 'nitrification(:,m,i_nh4_to_n2o) ',nitrification(test_grid,m,i_nh4_to_n2o)
          WRITE(numout,*) 'ftv=',ftv(test_grid)
          WRITE(numout,*) 'swc_pft=',swc_pft(test_grid,test_pft)
       ENDIF

       nitrification(:,m,i_nh4_to_n2o) = MIN(nitrification(:,m,i_nh4_to_no3),nitrification(:,m,i_nh4_to_n2o))  
 
       !
       ! 2.2.4 Production of NO during nitrification - NH4 to NO
       ! NO production during nitrification - Zhang et al., 2002 - Appendix p.102
       ! Reference NO production per N-NO3 produced g N-NO  (g N-NO3)-1
       ! no_nitrif_p = 0.0025
       !JC revise keep it as Zhang et al., 2002
       nitrification(:,m,i_nh4_to_no) = ftv(:) * no_nitrif_p * nitrification(:,m,i_nh4_to_no3)
       !revise again to XuRi's value RNON 0.1%-4% day-1 2% mean is used by XuRi
       ! nitrification(:,m,i_nh4_to_no) = RNON * nitrification(:,m,i_nh4_to_no3)

       IF(printlev>=4 .AND. m == test_pft)THEN
          WRITE(numout,*) 'CHECK values after nitrification nh4 to no'
          WRITE(numout,*) 'nitrification(:,m,i_nh4_to_no) ',nitrification(test_grid,m,i_nh4_to_no)
       ENDIF
       ! NO production from chemodenitrification
       ! based on Kesik et al., 2005, Biogeosciences
       ! BUT Kesik et al. used NO2 concentration and not NO3 production as it is done in OCN
       ! and modification of a multiplicative constant from 300 (Kesik) to 30 (OCN)
       ! without clear motivation - I don't how this is reliable 
       ! chemo_t0  = -31494.
       ! R  = 8.3144
       ! chemo_ph0 = -1.62
       ! chemo_0   = 30.
       ! chemo_1   = 16565.
       ! DSG: it might be worth a try to  use the original value of 300 (Kesik) here; it could work out (not tested  yet)

       !JC revise a scalling factor was added in OCN
       nitrification(:,m,i_nh4_to_no) = nitrification(:,m,i_nh4_to_no) + &
            ( chemo_0 * chemo_1 * exp(chemo_ph0 * (pH(:)+scal_ph) ) * &
            exp(chemo_t0/((temp_sol(:)+tp_00)*RR))) * nitrification(:,m,i_nh4_to_no3)

       nitrification(:,m,i_nh4_to_no) = MIN(nitrification(:,m,i_nh4_to_no3)-nitrification(:,m,i_nh4_to_n2o), &
            nitrification(:,m,i_nh4_to_no))

       IF(printlev>=4 .AND. m == test_pft)THEN
          WRITE(numout,*) 'CHECK values after nitrification nh4 to no chemodenitrification'
          WRITE(numout,*) 'nitrification(:,m,i_nh4_to_no) ',nitrification(test_grid,m,i_nh4_to_no)
          WRITE(numout,*) 'pH(:)=',pH(test_grid)
          WRITE(numout,*) 'temp_sol(:)=',temp_sol(test_grid)
       ENDIF

       ! In OCN, NO production and N2O production is deduced from the NO3 production as calculated
       ! in 2.2.2 (see below) -
        nitrification(:,m,i_nh4_to_no3) = nitrification(:,m,i_nh4_to_no3) - &
            (nitrification(:,m,i_nh4_to_no) + nitrification(:,m,i_nh4_to_n2o))
    ENDDO

!=======================================
    ! 2.3 Denitrification processes
    ! denitrification as in OCN (Zaehle et al. 2011); its a simplified version of the approach of Li et al. (2001)

    ! A temperature response of denitrification
    WHERE (temp_sol(:).GT.-46.01)
       ft_denit(:) = exp ( 308.56 * ( un/68.02 - un/ (temp_sol(:) + 46.02 ) ) )
    ELSEWHERE
       ft_denit(:) = zero
    ENDWHERE
    !JC revise if we want to follow Li et al., 2000 
    !! Temperature response of relative growth rate of total denitrifiers
    !! Eq. 2 Table 4 of Li et al., 2000
    !! ft_denit_0 = 2.
    !! ft_denit_1 = 22.5
    !! ft_denit_2 = 10.
    !ft_denit(:) = 2.**((temp_sol(:)-22.5)/10.)


    !JC revise fph from Li et al., 2000 But it is used for calculate relative growth
    !rate of total denitrifiers       
    ! Soil pH response of relative growth rate of total denitrifiers
    ! Eq. 2 Table 4 of Li et al., 2000
    ! See also comment about parenthesis' position in Thesis of Vincent
    ! Prieur (page 50)
    ! fph_no3_0  = 4.25
    ! fph_no3_1  = 0.5
    ! fph_no_0  = 5.25
    ! fph_no_1  = 1.
    ! fph_n2o_0  = 6.25
    ! fph_n2o_1  = 1.5
    ! B soil pH response of denitrification for nitrate, nox and nitrous oxide;
    ! this is simplified from Li et al (2000) 
    !JC revise a scalling factor was added in OCN
    ! fph_no3(:) = 1.0 - 1.0 / ( 1.0 + EXP((pH(:)-fph_no3_0)/fph_no3_1)) !nparam_scal  not done yet
    ! fph_no(:)  = 1.0 - 1.0 / ( 1.0 + EXP((pH(:)-fph_no_0) /fph_no_1))  !nparam_scal  not done yet
    ! fph_n2o(:) = 1.0 - 1.0 / ( 1.0 + EXP((pH(:)-fph_n2o_0)/fph_n2o_1)) !nparam_scal: not done yet
    fph_no3(:) = 1.0 - 1.0 / ( 1.0 + EXP((pH(:)+scal_ph-fph_no3_0)/fph_no3_1))
    fph_no(:)  = 1.0 - 1.0 / ( 1.0 + EXP((pH(:)+scal_ph-fph_no_0) /fph_no_1))
    fph_n2o(:) = 1.0 - 1.0 / ( 1.0 + EXP((pH(:)+scal_ph-fph_n2o_0)/fph_n2o_1))

    ! reference
    fph_no3_ref(:) = 1.0 - 1.0 / ( 1.0 + EXP((pH_ref-fph_no3_0)/fph_no3_1)) 
    fph_no_ref(:)  = 1.0 - 1.0 / ( 1.0 + EXP((pH_ref-fph_no_0) /fph_no_1)) 
    fph_n2o_ref(:) = 1.0 - 1.0 / ( 1.0 + EXP((pH_ref-fph_n2o_0)/fph_n2o_1))

    !JC revise From Li et al., 1992 it should be 0.083 kg N m-3
    ! if 1m soil depth, it should be 83 g N m-2
    ! if 0.01m soil depth, it is 0.83 g N m-2 But the value of 0.83 seems more
    ! reasonable if we see the normal nitrate concentration in soil
    ! C Half Saturation of N oxydes (gN/m2); value from OCN (Zaehle et al. 2011)
    Kn_conv(:) = .83 ! value chosen to get a nice curve for soil minberal nitrate 1-2g/m2

    DO m = 1, nvm
       ! D stocks of all N species (NO3, NO, N2O) (gN/m**2)

       ! E activity response to N concentration (Zaehle et al. 2011)
       mu_no3(:) = soil_n_min(:,m,initrate) / (soil_n_min(:,m,initrate) + Kn_conv(:))
       mu_no(:)  = soil_n_min(:,m,inox)     / (soil_n_min(:,m,inox)     + Kn_conv(:))
       mu_n2o(:) = soil_n_min(:,m,initrous) / (soil_n_min(:,m,initrous) + Kn_conv(:))
       
       !JC revise the qc is calculated from assumed soil respiration
       ! which is calculated using turnover time of active soil pool 7.3
       ! 17 g C m-3 is from Li et al., 1992 value 0.017 kg C m-3
       ! F activity of wetted soil active SOM; ignore temperature dependecy
       ! to avoid double accounting; includes unit conversion to gC and
       ! turnover time; from OCN; parameter values unclear (askSZ)
       qc(:) = ( som(:,iactive,m,icarbon) * un/(un/7.3*one_year) ) / &
                 ( som(:,iactive,m,icarbon) * un/(un/7.3 * one_year) + 17. )
          
       ! G compute actual denitrification rates from OCN (Zaehle et al., 2011)
       ! G1:
       denitrification(:,m,i_no3_to_n2) = MIN(MAX(soil_n_min(:,m,initrate)-leaching(:,m,initrate)-plant_uptake(:,m,initrate),zero), & 
                                                anvf(:,m) * soil_n_min(:,m,initrate) * &
                                                qc(:) * mu_no3(:) * ft_denit(:) * fph_no3(:) * dt) !DSG: *24?

       !JC revise
       ! here I think SZ used a relative increase of n2o production
       ! as writen in NMIP protocol
       ! 0.001 * 10 might indicate a maximum 1% can be as N2O increase
       ! should denitrification(:,m,i_nox_to_n2o) be increase by maximum 1% of the
       ! total denitrification, and increase also the part that will be goes to N2?
       ! in that case, denitrification(:,m,i_n2o_to_n2) should be added to 
       ! denitrification(:,m,i_nox_to_n2o)
       ! OR the relative fph has account for this?
       ! the 0.001 and 0.04 is tuned to get reported maximum at 25degC and constants at
       ! 6.5 PH
       ! G2: parameters are tuned (by SZ) to get reported maximum at 25degC and constants at pH6.5
       denitrification(:,m,i_no3_to_nox) = MIN(soil_n_min(:,m,inox), &
                                            0.001 * MIN(10.,MAX(0.,(fph_no3(:)/fph_n2o(:))/(fph_no3_ref(:)/fph_n2o_ref(:)))) &
                                                                 * ft_denit(:) * denitrification(:,m,i_no3_to_n2) )
       !JC revise
       !from XuRi & Prentice 2008 N2 increase should be Ndenit - NOdenit - N2Odenit
       !here 0.04 * 10 might indicate a maximum 40% of denitrification can be goes to
       !N2 
       !I still do not understand these two calculation
       denitrification(:,m,i_no3_to_n2o)  = MIN(soil_n_min(:,m,initrous), &
                                            0.04  * MIN(10.,MAX(0.,(fph_no(:)/fph_n2o(:))/(fph_no_ref(:)/fph_n2o_ref(:)))) &
                                                                 * ft_denit(:) * denitrification(:,m,i_no3_to_n2) )

       denitrification(:,m,i_no3_to_n2) = MIN(soil_n_min(:,m,initrate), &
              denitrification(:,m,i_no3_to_n2) - denitrification(:,m,i_no3_to_nox) - &
              denitrification(:,m,i_no3_to_n2o))

    ENDDO

!=======================================
  
    !
    ! 4. Loss of ionic forms of N through drainage and gaseous forms through
    ! volatilisation (the emission part should eventually be calculated in Sechiba's diffuco routines
    !
    ! DSG: very close to OCN, but a with a difference in the frac_nh3. I think
    !      OCN has it indeed wrong, thus I stick with NV version unless SZ can
    !      convince me his formulation makes sense.
    
    DO m = 1,nvm
  
       ! 4.1 Loss of NH4 due to volatilisation of NH3
       ! See Li et al. 1992, JGR, Table 4
  
       ! Current dissociation of [NH3] to [NH4+]
       ! See Table 4 of Li et al. 1992 and Appendix A of Zhang et al. 2002
       ! log(K_NH4) - log(K_H20) = log(NH4/NH3) + pH
       ! See also formula in "DISSOCIATION CONSTANTS OF INORGANIC ACIDS AND BASES" pdf file
       ! pK_H2O = -log(K_H2O) = 14
       ! pk_NH4 = -log(K_NH4) = 9.25
       ! pK_H2O - pK_NH4 = log([NH4+]/[NH3]) + pH
       ! [NH4+]/[NH3] = 1O^(pK_H2O - pK_NH4 - pH) = 10^(4.75-pH)
       
       ! In OCN, one makes use of frac_nh3. That should be the NH3/NH4 ratio. It is defined as:
       ! frac_nh3(:) = 10.0**(4.25-pH(:)) / (1. + 10.0**(4.25-pH(:)))
       
       ! CommentS of N. Vuichard : I have several interrogations about the equation and value used
       ! in OCN. But I have also concerns about the formulation of Li et al. ... 
       ! 1/ 
       ! The formulation of Li et al. doesn't match with the formulas in 
       ! "DISSOCIATION CONSTANTS OF INORGANIC ACIDS AND BASES" or with the formulas
       ! of http://www.onlinebiochemistry.com/obj-512/Chap4-StudNotes.html
       ! To my opinion, one should replace log(K_NH4+) by log(K_NH3) in the formulation of Li et al.
       ! with the relationship pKa + pKb = pKwater (where a and b are acid and base)
       ! this leads to [NH4+]/[NH3] = 10^(pK_NH4 - pH) = 10^(9.25 - pH)
       ! 2/
       ! Wether I'm right or not about 1/, I don't understand the value used in OCN (4.25). 
       ! It should be either 4.75 or 9.25 but 4.25 looks strange
       ! 3/ 
       ! OCN used a formulation for [NH3]/[NH4+] of the type: X/(1+X) with X=[NH3]/[NH4+]
       ! This leads to X/(1+X)=([NH3]/[NH4+])/(([NH4+]/[NH4+])+([NH3]/[NH4+]))
       ! or X/(1+X) = [NH3] / ( [NH3] + [NH4+] )
       ! This means that the value stored in soil_n_min(:,:,iammonium) corresponds to the total
       ! N of both [NH4+] and [NH3]. This makes sense to my opinion. But I wonder if one should not
       ! account for this partitioning between [NH4+] and [NH3] in other processes. To check.
       ! 4/ 
       ! The X value should relate to [NH3]/[NH4+] but to my opinion the X value used 
       ! in the equation in OCN corresponds to [NH4+]/[NH3]. Is this a bug ? 
   
       ! In conclusion, I propose the formulation
       frac_nh3(:) = 10.0**(pH(:)-pk_NH4) / (1. + 10.0**(pH(:)-pk_NH4))
       !DSG: I agree with NV about the NH4/NH3 vs NH3/NH4 
       !     also about the pk_NH4; it should be 9.25 25degreeC
       !     I sent an email to SZ 2016/08/08 asking for clarification
  
       ! This seems a patch added to OCN in order to (as mentioned in OCN) 
       ! reduced emissions at low concentration (problem of only one soil layer)
       ! high conentrations are usually associated with fertiliser events -> top layer
       ! and thus increased emission
       ! NOT ACTIVATE HERE 
       ! emm_fac(:) = 0.01
       ! WHERE(soil_n_min(:,m,iammonium).GT.0.01.AND.soil_n_min(:,m,iammonium).LE.4.)
       !      emm_fac(:) = MAX(0.01,1-exp(-(soil_n_min(:,m,iammonium)/1.75)**8)) 
       ! ENDWHERE       
       ! WHERE(soil_n_min(:,m,iammonium).GT.4.)
       !      emm_fac(:)=1.0
       ! ENDWHERE
       emm_fac(:)=1.0
      !DSG: not activated yet, as we have too much N in the soils
  
       ! 4.2 Volatilisation of gasous species, Table 4, Li et al. 2000 I,
       !     using diffusivity of oxigen in air as a surrogate (from OCN)
       !     assumes no effect of air concentration on diffusion
       !     takes standard depth as reference
       
       ! Clay limitation
       ! F_clay_0 = 0.13
       ! F_clay_1 = -0.079
       F_clay(:) = F_clay_0 + F_clay_1 * clay(:)
  
       ! In OCN, one used formulation of the type
       ! emission(:,m,inox-1) = & 
       !       MIN( d_ox(:) * soil_n_min(:,m,inox) * (0.13-0.079*clay(:)) * dt / z_decomp
       ! It should have the unit d_ox * soil_n_min * dt / z_decomp
       !                         m2 day-1 * gN m-2 * day / m
       !                         gN / m which is not homogeneous with the unit expected (gn m-2)
       ! To my opinion, one should not use soil_n_min (gN m-2) but a volumetric concentration (gN m-3)
       ! But I'm not clear what is the appropriate volume to consider (volume of soil ?)
       ! DSG: volume of water:
  
       ! NH4 emission (gN m-2 per time step)
       WHERE(tmc_pft(:,m).GT.zero)
          emission(:,m,iammonium) = D_s(:,m) * emm_fac(:) * frac_nh3(:) * soil_n_min(:,m,iammonium)  / (tmc_pft(:,m)/1000.) &
               * F_clay(:) / z_decomp * dt
               
          ! NO NO3 emission 
          emission(:,m,initrate) = 0.
  
          ! NO2 emission (gN m-2 per time step)
          emission(:,m,inox) = D_s(:,m) * soil_n_min(:,m,inox)        /(tmc_pft(:,m)/1000.)  * F_clay(:) / z_decomp * dt
  
          ! N2O emission (gN m-2 per time step)
          emission(:,m,initrous) = D_s(:,m) * soil_n_min(:,m,initrous)/(tmc_pft(:,m)/1000.) * F_clay(:) / z_decomp * dt
                      
          ! N2 emission (gN m-2 per time step)
          emission(:,m,idinitro) = D_s(:,m) * soil_n_min(:,m,idinitro)/(tmc_pft(:,m)/1000.) * F_clay(:) / z_decomp * dt
       ELSEWHERE
          emission(:,m,iammonium) = zero
          emission(:,m,initrate)  = zero
          emission(:,m,inox)      = zero
          emission(:,m,initrous)  = zero
          emission(:,m,idinitro)  = zero
       ENDWHERE 
    ENDDO
    ! don't do it like in OCN; OCN seems to be buggy in respect to frac_nh3,
    ! rest is the same: askSZ 





    ! ensure mass conservation
    emission(:,:,iammonium) = MIN(emission(:,:,iammonium), &
         soil_n_min(:,:,iammonium) - nitrification(:,:,i_nh4_to_no3) &
         - nitrification(:,:,i_nh4_to_no) - nitrification(:,:,i_nh4_to_n2o) &
         - leaching(:,:,iammonium) - plant_uptake(:,:,iammonium))  

    emission(:,:,inox) = MIN(emission(:,:,inox), &
         soil_n_min(:,:,inox) + nitrification(:,:,i_nh4_to_no) & 
          + denitrification(:,:,i_no3_to_nox))

    emission(:,:,initrous) = MIN(emission(:,:,initrous), &
         soil_n_min(:,:,initrous) + nitrification(:,:,i_nh4_to_n2o) & 
           + denitrification(:,:,i_no3_to_n2o)) 

    emission(:,:,idinitro) = MIN(emission(:,:,idinitro), &
          soil_n_min(:,:,idinitro)  + denitrification(:,:,i_no3_to_n2)) 
    !
    ! 5. Update pools
    !
  
    ! 5.1 update pools of nitrogen in the soil from nitrification and 
    ! denitrification, plant uptake, leaching and volatile emissions,
    ! desorption from clay, and net mineralisation
    !
    ! To my opinion, for a better consistency, I would recommand to consider the leaching separately 
    ! when calculating the ammonium and nitrate budget. Leaching is calculated from sechiba
    ! Might be better to remove leaching at the top of the routine especially due to the later calculation of
    ! NH4+ and NO3- concentration that will vary with the soil water content
    ! Let's imagine that from one time step to another, the change in soil water content is only due to leaching
    ! We don't want that the NH4+ and NO3- concentration vary from one time step to the other. The best way to avoid
    ! this is to remove first the leaching from the NH4+ and NO3- pools
    ! THIS IS NOT DONE YET

    IF(printlev>=4)THEN
       WRITE(numout,*) 'CHECK values before update'
       WRITE(numout,*) 'nitrification ',nitrification(test_grid,test_pft,:)
       WRITE(numout,*) 'denitrification ',denitrification(test_grid,test_pft,:)
       WRITE(numout,*) 'leaching ',leaching(test_grid,test_pft,:)
       WRITE(numout,*) 'emission ',emission(test_grid,test_pft,:)
       WRITE(numout,*) 'plant_uptake ',plant_uptake(test_grid,test_pft,:)
       WRITE(numout,*) 'mineralisation ',mineralisation(test_grid,test_pft)
       WRITE(numout,*) 'immob ',immob(test_grid,test_pft)
    ENDIF

    ! DSG: The problem with all the following bilancing is that we sometimes end up with
    ! very minor negative values for soil_n_min due to machine precision which will mess up many
    ! calculations which depend on soil_n_min having positive value. 
    ! I hate to do so, but I introduce max(x,zero) function, which can easily result in violation of mass
    ! conservation; thus it is crucial to ensure mass conservation by the use of
    ! checks like the ones I use. The root of all evil is the structure of the routine nitrogen_dynamics 
    ! (it is not the only routine of ORCHIDEE) which makes the code prone to mass violation.

    soil_n_min(:,:,iammonium) = MAX(soil_n_min(:,:,iammonium) + n_adsorbed(:,:) & 
         - nitrification(:,:,i_nh4_to_no3) - nitrification(:,:,i_nh4_to_no) - nitrification(:,:,i_nh4_to_n2o) &
         - leaching(:,:,iammonium) - emission(:,:,iammonium) &
         + mineralisation(:,:) + immob(:,:) - plant_uptake(:,:,iammonium), zero)
  
    soil_n_min(:,:,initrate) = MAX(soil_n_min(:,:,initrate) + nitrification(:,:,i_nh4_to_no3) & 
          - denitrification(:,:,i_no3_to_n2) - denitrification(:,:,i_no3_to_nox) &
          - denitrification(:,:,i_no3_to_n2o) - leaching(:,:,initrate) & 
         - plant_uptake(:,:,initrate), zero)
  
    soil_n_min(:,:,inox) =  MAX(soil_n_min(:,:,inox) + nitrification(:,:,i_nh4_to_no) & 
         + denitrification(:,:,i_no3_to_nox) - emission(:,:,inox) , zero)
  
    soil_n_min(:,:,initrous) = MAX(soil_n_min(:,:,initrous) + nitrification(:,:,i_nh4_to_n2o) & 
          + denitrification(:,:,i_no3_to_n2o) - emission(:,:,initrous),zero)

     soil_n_min(:,:,idinitro) = MAX(soil_n_min(:,:,idinitro)  + denitrification(:,:,i_no3_to_n2) & 
         - emission(:,:,idinitro), zero)


    IF(printlev>=4)THEN
       WRITE(numout,*) 'CHECK values after update'
       WRITE(numout,*) 'soil_n_min ',soil_n_min(test_grid,test_pft,:)
    ENDIF

    mass_change(:,:) = zero
    DO m=1,nvm
       ! Deposition of NHx and NOy
       WHERE(veget_max(:,m).GT.min_stomate.AND.som(:,iactive,m,icarbon).GT.min_stomate) 
          soil_n_min(:,m,iammonium) = soil_n_min(:,m,iammonium) & 
               + n_input(:,m,iammonium)*dt + &
            ! JCADD 20170421 for CNP merge
               grm_nfert(:,m,iammonium)*dt + &
               agri_nfert(:,m,iammonium)*dt
            ! END JCADD
          soil_n_min(:,m,initrate) = soil_n_min(:,m,initrate) & 
               + n_input(:,m,initrate)*dt + &
            ! JCADD 20170421 for CNP merge
               grm_nfert(:,m,initrate)*dt + &
               agri_nfert(:,m,initrate)*dt
            ! END JCADD
          !mass_cons check:     
          mass_change(:,m) = mass_change(:,m) + &
                  n_input(:,m,iammonium)*dt + n_input(:,m,initrate)*dt + &
                  grm_nfert(:,m,iammonium)*dt + grm_nfert(:,m,initrate)*dt  + &
                  agri_nfert(:,m,iammonium)*dt + agri_nfert(:,m,initrate)*dt
       ENDWHERE

       WHERE(veget_max(:,m).GT.min_stomate.AND.som(:,iactive,m,icarbon).LE.min_stomate) 
          leaching(:,m,iammonium)=leaching(:,m,iammonium) + &
              n_input(:,m,iammonium)*dt+grm_nfert(:,m,iammonium)*dt + &
              agri_nfert(:,m,iammonium)*dt
          leaching(:,m,initrate)=leaching(:,m,initrate) + &
              n_input(:,m,initrate)*dt+grm_nfert(:,m,initrate)*dt + &
              agri_nfert(:,m,initrate)*dt
          !mass_cons check:     
          ! JC comments why there will be leaching of all n_input?
          ! soil can not have n if there is no som
          ! In this case, N input to baresoils is leached, no N added to
          ! soil_n_min. However, in the mass_change calculation at the end of
          ! this subroutine, it minus all leaching, while this direct leaching
          ! part is acturally not from soil. So I first added this part to avoid
          ! over leaching and keep it mass conserved.
          mass_change(:,m) = mass_change(:,m) + &!+ n_input(:,m,initrate) + n_input(:,m,iammonium)
              (n_input(:,m,iammonium)*dt+grm_nfert(:,m,iammonium)*dt + &
              agri_nfert(:,m,iammonium)*dt) + &
              (n_input(:,m,initrate)*dt+grm_nfert(:,m,initrate)*dt + &
              agri_nfert(:,m,initrate)*dt)
       ENDWHERE
       ! BNF

       ! calculate the response function of BNF to P/N ratio of labile pools
       ! this scaling function is an extension for the Cleveland NPP-BNF
       ! correlation model to reduce BNF when plants dont need it (=N/P is high)
       ! JC Debug NOTE the calc_f_XYplant is applied for every PFT
       ! it can be moved out of the DO nvm cycle, but it will not provide any output
       ! since f_PNplant has been set to 0
       IF(ok_pcycle) THEN
          CALL calc_f_XYplant(npts,nvm,  &
                              biomass,   &
                              iphosphorus, & ! = X
                              initrogen,   & ! = Y
                              f_PNplant  &
                               )
       ELSE
           f_PNplant = un
       ENDIF

       ! sigmoidal scaling of BNF (using a product inhibitation assumption)
       ! the two coefficient are chosen that f_Nmin is close to 1 for soil_n_min
       ! and 0 for soil_n_min >1. This might need some further calibration
       f_Nmin(:,m) = un - un/(un +exp((-(                                      &
                            soil_n_min(:,m,initrate)+soil_n_min(:,m,iammonium) &
                                             )*10.)+7.))

       ! set to zero; deactivate BNF scaling via plant tissue N/P by setting it to zero
       f_PNplant=zero

       IF ( natural(m) ) THEN 
          IF (read_bnf) THEN 
             ! in the presence of a organic soil component assume that there is also BNF
             ! as long as plant available nitrogen is not too high
             WHERE(.NOT.(som(:,iactive,m,icarbon).GT.min_stomate.AND. & 
               soil_n_min(:,m,iammonium)+soil_n_min(:,m,initrate).LT.max_soil_n_bnf(m)))

                 n_input(:,m,ibnf) = zero
             ENDWHERE
          ELSE

             ! JCMOD060 29Mar2018
             ! biologicalN2fixation will not be applied for grassland if
             ! GRM_allow_BNF is activated
                IF (is_tree(m) .OR. .NOT. GRM_allow_BNF) THEN             
!               !! Calculate Biological N2 fixation according to Cleveland et al (1999) following CLM & JSBACH; modified

                CALL biologicalN2fixation (npts       &
                                           ,npp_week(:,m)  &
   !demand driven:
              !                             ,f_NCplant(:,m) & 
              !                             ,f_PNplant(:,m) & 
   !supply driven:
                                           ,f_Nmin(:,m)       &
                                           ,f_PNplant(:,m)    & ! is set to zero
   !supply driven:
                                           ,n_input(:,m,ibnf)  &  ! in g(N)/m2/day
                                           )
  
              IF(printlev>=4)THEN
                   WRITE(numout,*)'f_PNplant',f_PNplant(test_grid,m)
                   WRITE(numout,*)'npp_week',npp_week(test_grid,m)
                   WRITE(numout,*)'dt',dt
              ENDIF
                 ENDIF ! for tree only or GRM_allow_BNF is not activated
          ENDIF ! read_bnf


          soil_n_min(:,m,iammonium) = soil_n_min(:,m,iammonium) + n_input(:,m,ibnf)*dt
          !mass_cons check:     
          mass_change(:,m) = mass_change(:,m) + n_input(:,m,ibnf)*dt

       ENDIF
    ENDDO

    IF(printlev>=4)THEN
       WRITE(numout,*) 'CHECK values after N input'
       WRITE(numout,*) 'soil_n_min ',soil_n_min(test_grid,test_pft,:)
       WRITE(numout,*) 'end of nitrogen_dynamics'
    ENDIF

    !DSG 
    IF(ANY(soil_n_min(:,:,:).LT.zero)) THEN
       DO m=1,nvm
         DO n=1,npts
           IF (ANY(soil_n_min(n,m,:).LT.zero)) THEN
               WRITE(numout,*) 'box',n
               WRITE(numout,*) 'PFT',m
               WRITE(numout,*) 'soil_n_min',soil_n_min(n,m,:)
               WRITE(numout,*) 'n_input',n_input(n,:,:)
               WRITE(numout,*) 'nitrification(n,m,i_nh4_to_no3)', nitrification(n,m,i_nh4_to_no3)
               WRITE(numout,*) 'nitrification(n,m,i_nh4_to_no)' , nitrification(n,m,i_nh4_to_no)
               WRITE(numout,*) 'nitrification(n,m,i_nh4_to_n2o)', nitrification(n,m,i_nh4_to_n2o)
               WRITE(numout,*) 'nitrification(n,m,i_no3_to_nox)', nitrification(n,m,i_no3_to_nox)
               WRITE(numout,*) 'nitrification(n,m,i_nox_to_n2o)', nitrification(n,m,i_no3_to_nox)
               WRITE(numout,*) 'nitrification(n,m,i_no3_to_n2o)' , nitrification(n,m,i_no3_to_n2o)
               WRITE(numout,*) 'leaching(n,m,:)',leaching(n,m,:)
               WRITE(numout,*) 'emission(n,m,:)',emission(n,m,:)
               WRITE(numout,*) 'mineralisation(n,m)',mineralisation(n,m)
               WRITE(numout,*) 'immob(n,m)',immob(n,m)
               WRITE(numout,*) 'plant_uptake(n,m,:)',plant_uptake(n,m,:)
           ENDIF
         ENDDO 
       ENDDO
       STOP

    ENDIF

    !DSG mass conservation ========================================
    mass_change(:,:)   = mass_change(:,:)                           &  ! already includes input
                          ! JC Debug, not all mineralisation has been applied,
                          ! for positive mineralisation, all has been added to soil_n_min
                          ! for negative ones, only used_immob has been reduced from soil_n_min
                          ! + mineralisation(:,:)                     &
                          + mineralisation(:,:) + N_added_cons(:,:) &
                          - SUM(leaching(:,:,:),DIM=3)                      &
                          - SUM(emission(:,:,:),DIM=3)                      &
                          - SUM(plant_uptake(:,:,:),DIM=3)
    mass_after = SUM(soil_n_min(:,:,:),DIM=3) 
    CALL cons_mass(mass_before(:,:),              &  ! mass before
          mass_after(:,:),               &  ! mass after
          mass_change(:,:),              &  ! net of fluxes
          'end of nitrogen_dynamics' )

    CALL histwrite_p (hist_id_stomate, 'BNF_CLOVER', itime, &
         BNF_clover(:,:)/dt, npts*nvm, horipft_index) 
    CALL histwrite_p (hist_id_stomate, 'N_UPTAKE_NH4', itime, &
         plant_uptake(:,:,iammonium)/dt, npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'N_UPTAKE_NO3', itime, &
         plant_uptake(:,:,initrate)/dt, npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'N_MINERALISATION', itime, &
         mineralisation(:,:)/dt, npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'SOIL_NH4', itime, &
         soil_n_min(:,:,iammonium), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'SOIL_NO3', itime, &
         soil_n_min(:,:,initrate), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'SOIL_NOX', itime, &
         soil_n_min(:,:,inox), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'SOIL_N2O', itime, &
         soil_n_min(:,:,initrous), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'SOIL_N2', itime, &
         soil_n_min(:,:,idinitro), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'SOIL_P_OX', itime, &
         p_o2(:,:), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'NH3_EMISSION', itime, &
         emission(:,:,iammonium)/dt, npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'NOX_EMISSION', itime, &
         emission(:,:,inox)/dt, npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'N2O_EMISSION', itime, &
         emission(:,:,initrous)/dt, npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'N2_EMISSION', itime, &
         emission(:,:,idinitro)/dt, npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'NO3_LEACHING', itime, &
         leaching(:,:,initrate)/dt, npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'NH4_LEACHING', itime, &
         leaching(:,:,iammonium)/dt, npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'NITRIFICATION', itime, &
         nitrification(:,:,i_nh4_to_no3)/dt, npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'DENITRIFICATION', itime, &
         denitrification(:,:,i_no3_to_n2o)/dt, npts*nvm, horipft_index)

! DSG: they are in day-1
    CALL histwrite_p (hist_id_stomate, 'NHX_DEPOSITION', itime, &
         n_input(:,1,iatm_ammo), npts, hori_index)
    CALL histwrite_p (hist_id_stomate, 'NOY_DEPOSITION', itime, &
         n_input(:,1,iatm_nitr), npts, hori_index)
    CALL histwrite_p (hist_id_stomate, 'BNF', itime, &
         n_input(:,:,ibnf), npts*nvm, horipft_index)

!DSG_IFERT    CALL histwrite_p (hist_id_stomate, 'N_FERTILISER', itime, &
!DSG_IFERT         n_input(:,1,ifert), npts, hori_index)  
    CALL histwrite_p (hist_id_stomate, 'AGRI_FERT_NH4', itime, &
         agri_nfert(:,:,iammonium), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'AGRI_FERT_NO3', itime, &
         agri_nfert(:,:,initrate), npts*nvm, horipft_index)

    CALL xios_orchidee_send_field("BNF_CLOVER",BNF_clover(:,:)/dt)
    CALL xios_orchidee_send_field("N_UPTAKE_NH4",plant_uptake(:,:,iammonium)/dt)
    CALL xios_orchidee_send_field("N_UPTAKE_NO3",plant_uptake(:,:,initrate)/dt)
    CALL xios_orchidee_send_field("N_MINERALISATION",mineralisation(:,:)/dt)
    CALL xios_orchidee_send_field("N_IMMOBILISATION",immob(:,:)/dt)
    CALL xios_orchidee_send_field("N_ADDED_CONS",N_added_cons(:,:)/dt)    
    CALL xios_orchidee_send_field("SOIL_NH4",soil_n_min(:,:,iammonium))
    CALL xios_orchidee_send_field("SOIL_NO3",soil_n_min(:,:,initrate))
    CALL xios_orchidee_send_field("SOIL_NOX",soil_n_min(:,:,inox))
    CALL xios_orchidee_send_field("SOIL_N2O",soil_n_min(:,:,initrous))
    CALL xios_orchidee_send_field("SOIL_N2",soil_n_min(:,:,idinitro))
    CALL xios_orchidee_send_field("SOIL_P_OX",p_o2(:,:))
    CALL xios_orchidee_send_field("NH3_EMISSION",emission(:,:,iammonium)/dt)
    CALL xios_orchidee_send_field("NOX_EMISSION",emission(:,:,inox)/dt)
    CALL xios_orchidee_send_field("N2O_EMISSION",emission(:,:,initrous)/dt)
    CALL xios_orchidee_send_field("N2_EMISSION",emission(:,:,idinitro)/dt)
    CALL xios_orchidee_send_field("NH4_LEACHING",leaching(:,:,iammonium)/dt)
    CALL xios_orchidee_send_field("NO3_LEACHING",leaching(:,:,initrate)/dt)
    CALL xios_orchidee_send_field("NITRIFICATION",nitrification(:,:,i_nh4_to_no3)/dt)
    CALL xios_orchidee_send_field("DENITRIFICATION",denitrification(:,:,i_no3_to_n2)/dt)
    CALL xios_orchidee_send_field("NITRIF_NO3",nitrification(:,:,i_nh4_to_no3)/dt)
    CALL xios_orchidee_send_field("DENITRIF_NOX",denitrification(:,:,i_no3_to_nox)/dt)
    CALL xios_orchidee_send_field("NITRIF_N2O",nitrification(:,:,i_nh4_to_n2o)/dt)
    CALL xios_orchidee_send_field("DENITRIF_N2O",denitrification(:,:,i_no3_to_n2o)/dt)
    CALL xios_orchidee_send_field("NITRIF_NO",nitrification(:,:,i_nh4_to_no)/dt)
    CALL xios_orchidee_send_field("DENITRIF_N2",denitrification(:,:,i_no3_to_n2)/dt)

    CALL xios_orchidee_send_field("NHX_DEPOSITION",n_input(:,1,iatm_ammo))
    CALL xios_orchidee_send_field("NOY_DEPOSITION",n_input(:,1,iatm_nitr))
    CALL xios_orchidee_send_field("BNF",n_input(:,:,ibnf))
!DSG_IFERT    CALL xios_orchidee_send_field("N_FERTILISER",n_input(:,1,ifert))  
    CALL xios_orchidee_send_field("AGRI_FERT_NH4",agri_nfert(:,:,iammonium))
    CALL xios_orchidee_send_field("AGRI_FERT_NO3",agri_nfert(:,:,initrate))

    CALL xios_orchidee_send_field("fN2O",SUM(emission(:,:,initrous)*veget_max,dim=2)/dt/1e3/one_day*contfrac)
    vartmp(:)=zero
    vartmp(:)=n_input(:,1,iatm_ammo)+n_input(:,1,iatm_nitr)
    CALL xios_orchidee_send_field("fNdep",vartmp(:)/1e3/one_day*contfrac)
    CALL xios_orchidee_send_field("fBNF",SUM(n_input(:,:,ibnf)*veget_max,dim=2)/1e3/one_day*contfrac)
    vartmp(:)=SUM((plant_uptake(:,:,iammonium)+plant_uptake(:,:,initrate))*veget_max,dim=2)
    CALL xios_orchidee_send_field("fNup",vartmp(:)/dt/1e3/one_day*contfrac)
    CALL xios_orchidee_send_field("fNnetmin",SUM(mineralisation(:,:)*veget_max,dim=2)/dt/1e3/one_day*contfrac)
    CALL xios_orchidee_send_field("fNit",SUM(nitrification(:,:,i_nh4_to_no3)*veget_max,dim=2)/dt/1e3/one_day*contfrac)
    CALL xios_orchidee_send_field("fDenit",SUM(denitrification(:,:,i_no3_to_n2)*veget_max,dim=2)/dt/1e3/one_day*contfrac)
    vartmp(:)=SUM((leaching(:,:,iammonium)+leaching(:,:,initrate))*veget_max,dim=2)
    CALL xios_orchidee_send_field("fNleach",vartmp(:)/dt/1e3/one_day*contfrac)
    CALL xios_orchidee_send_field("fNH3",SUM(emission(:,:,iammonium)*veget_max,dim=2)/dt/1e3/one_day*contfrac)
    CALL xios_orchidee_send_field("fNOX",SUM(emission(:,:,inox)*veget_max,dim=2)/dt/1e3/one_day*contfrac)
    CALL xios_orchidee_send_field("fN2",SUM(emission(:,:,idinitro)*veget_max,dim=2)/dt/1e3/one_day*contfrac)


 
  END SUBROUTINE nitrogen_dynamics

!! ================================================================================================================================
!! SUBROUTINE     : biologicalN2fixation
!>\BRIEF          Calculate the N fixed from the atmosphere and made available to bioa
!                 1. following Cleveland et al. (1999) using NPP as driver instead of ET (Thornton et al., 2007,Parida et al, 2010, Goll et al.,2012)
!                 2. the maximum rate is adjusted to the highest average rate for non-tropical biomes found in table 13 of Cleveland et al. 1999. 
!                    as the tropical fixation rates in Cleveland are flawed (see Sullivan et al. 2014, PNAS)


  SUBROUTINE biologicalN2fixation (npts       &
                                   ,npp_week  &
                                   ,f_NCplant & 
                                   ,f_PNplant & 
                                   ,bnf       &
                                   )

    INTEGER, INTENT(in)                            :: npts             !! Domain size (unitless)
    ! driver
    REAL, DIMENSION(npts)  ,INTENT(in)         :: npp_week         !! "weekly" NPP (gC/day/(m**2 covered)
    REAL, DIMENSION(npts)  ,INTENT(in)         :: f_NCplant        !! We scale BNF using the same equation as was used to scale N uptake
    REAL, DIMENSION(npts)  ,INTENT(in)         :: f_PNplant        !! We scale BNF using the same equation as was used to scale N uptake

    !output
    REAL, DIMENSION(npts)  ,INTENT(out)        :: bnf              !! atmospheric nitrogen fixed (gN/m2/day) ! checked

    ! BNF scaling using NPP: (ab)use of the empricial correlation between BNF and
    ! evpotranspiration from Cleveland et al (1999) following NCAR & JSBACH

    ! conversion from daily NPP to annual NPP & annual BNF to daily BNF
    bnf(:) = (BNF_scal * (un - exp(BNF_coef * npp_week/one_day*one_year))) /one_year * one_day

    ! modification of original model to account for downregulation of BNF when
    ! N is in ample supply: (DSG: might it be better to make it function of
    ! mineralN?)
    bnf(:) = bnf(:) * f_NCplant(:) * (un - f_PNplant(:))

    bnf(:) = max(zero,bnf(:)) !! To prevent negative BNF in case of negative NPP   


  END SUBROUTINE biologicalN2fixation

!!
!================================================================================================================================
!! SUBROUTINE     : calc_BNF_clover
!>\BRIEF          Calculate the BNF the atmosphere by legumes in managed
!grassland
!                 Method default following PaSim: Vuichard Thesis p56 eq34 35
!                 ! Schmid 2001 equation (11)
!                 new Method following Lazzarotto et al. 2009

  SUBROUTINE calc_BNF_clover(npts,nvm,dt, &
                soil_n_min,fclover,plant_uptake, &
                temp_sol,biomass,BNF_clover)

    INTEGER, INTENT(in)                               :: npts        !! Domain size
    INTEGER, INTENT(in)                               :: nvm         !! Domain size (unitless
    REAL, INTENT(in)                                  :: dt          !! time step
    REAL, DIMENSION(npts,nvm,nnspec),INTENT(in)       :: soil_n_min  !! mineral nitrogen in the soil (gN/m**2)
                                                                            !! (first index=npts, second index=nvm, third index=nnspec)
    REAL, DIMENSION(npts,nvm), INTENT(in)             :: fclover     !! legume fraction
    REAL, DIMENSION(npts,nvm,nionspec), INTENT(in)    :: plant_uptake!! Uptake of soil N by platns
                                                                            !! (gN/m**2/timestep)
    REAL, DIMENSION(npts), INTENT(in)                 :: temp_sol    !! soil temperature (degC)
    REAL, DIMENSION(npts,nvm,nparts,nelements), INTENT(in):: biomass !! Biomass pools (gC(or N)/m**2)
    !! BNF by clover (legume) nodules (assuming directly into plant rather than
    !! soil), it is different from the bnf fixed by soil bacteria
    REAL, DIMENSION(npts,nvm) , INTENT(out)           :: BNF_clover  !! see above

    ! local variables and parameters
    REAL, DIMENSION(npts)        :: Un_1neff                !! BNF parameter
    REAL, PARAMETER              :: Un_1kneff = 1.0         !! BNF parameter (gN/m2)
    REAL, PARAMETER              :: BNF_1ecostbnf = 0.5     !! parameter from PaSim
    REAL, PARAMETER              :: BNF_1bnfclovermin = 0.15!! clover N uptake reduction
    ! parameters subroutine :: fixation_legume_2
    ! According to Lazzarotto et al. 2009
    REAL, PARAMETER              :: T_BNF1 = 282.15         !! parameter defining the temperature dependence of BNF (K)
    REAL, PARAMETER              :: T_BNF2 = 286.15         !! parameter defining the temperature dependence of BNF (K)
    REAL, PARAMETER              :: T_BNF3 = 299.15         !! SAVE defining the temperature dependence of BNF (K)
    REAL, PARAMETER              :: T_BNF4 = 303.15         !! SAVE defining the temperature dependence of BNF (K)
    REAL, PARAMETER              :: BNF_max= 0.0024         !! maximal BNF rate by legume (kg N (kg root DM)-1 d-1)
    REAL, PARAMETER              :: fUN_leg= 0.95           !! maximal fraction of N uptake by legume from the soil
    REAL, PARAMETER              :: K_BNF  = 6              !! substrate concentration at which BNF is half maximal (g N m-2)
    REAL                         :: fT_BNF                  !! Temperature factor for BNF
    REAL, DIMENSION(npts)        :: vartmp                  !! Temporary variable

    INTEGER                      :: m,n                     !! Index to loop over 


    BNF_clover(:,:) = zero
    DO m = 2, nvm
       ! only for grassland PFT
       IF (natural(m) .AND. .NOT. is_tree(m) .AND. GRM_allow_BNF) THEN
          IF (GRM_BNF_newmethod .EQV. .FALSE.) THEN
          ! BNF calculation derived from PaSim
          ! for grassland only natural and not tree
          ! Vuichard Thesis p56 eq34 35
          ! Schmid 2001 equation (11)
          ! this one is based on the plant_uptake calculated as PaSim
          Un_1neff(:) = (soil_n_min(:,m,iammonium) + soil_n_min(:,m,initrate))/&
             (soil_n_min(:,m,iammonium) + soil_n_min(:,m,initrate)+Un_1kneff)
          BNF_clover(:,m) = MAX(zero, &
             BNF_1ecostbnf * fclover(:,m) * &
             (plant_uptake(:,m,iammonium)+plant_uptake(:,m,initrate)) * &
             (un-Un_1neff(:)*(un-BNF_1bnfclovermin)) / &
             (un-fclover(:,m)*BNF_1bnfclovermin)/&
             Un_1neff(:))
!!          ! following original calculation from Schwinning and Parsons, The
!!          ! clover fraction can gain from symbiotic N2 fixation 60% of the
!!          ! difference between the mineral uptake and potential N uptake.
!!          BNF_clover(:,m) = MAX(zero, &
!!             BNF_1ecostbnf * fclover(:,m) * &
!!             (root_uptake(:,m,iammonium) + root_uptake(:,m,initrate))  * &
!!             biomass(:,m,iroot,icarbon) *ft_uptake(:) * &
!!             (un-f_NCplant(:,m)) )

          ! above the BNF_clover is g N/m2/dt the same as plant_uptake
!         WRITE(numout,*) 'BNF_clover',BNF_clover(1,m)
          ELSE
          ! Lazzarotto et al. 2009
            DO n=1, npts
              IF (temp_sol(n)+tp_00 .LE. T_BNF1) THEN
                fT_BNF = 0.0
              ELSE IF ( (temp_sol(n)+tp_00 .GT. T_BNF1) .AND. (temp_sol(n)+tp_00 .LE. T_BNF2) ) THEN
                fT_BNF = (temp_sol(n)+tp_00 - T_BNF1)/(T_BNF2-T_BNF1)
              ELSE IF ( (temp_sol(n)+tp_00 .GT. T_BNF2) .AND. (temp_sol(n)+tp_00 .LE. T_BNF3) ) THEN
                fT_BNF = 1.0
              ELSE IF ( (temp_sol(n)+tp_00 .GT. T_BNF3) .AND. (temp_sol(n)+tp_00 .LE. T_BNF4) ) THEN
                fT_BNF = 1.0 - (temp_sol(n)+tp_00-T_BNF3)/(T_BNF4-T_BNF3)
              ELSE
                fT_BNF = 0.0
              END IF
              ! calculation gives g N/m2/day need to *dt
              ! here BNF_max = 0.0024 g g-1DM day-1
              BNF_clover(n,m) = MAX(zero, fT_BNF * BNF_1ecostbnf * BNF_max * fclover(n,m) * &
                  (biomass(n,m,iroot,icarbon)/0.45) * & ! root biomass convert to gDM m-2
                  (un - fUN_leg * (soil_n_min(n,m,iammonium) + soil_n_min(n,m,initrate))/&
                  (soil_n_min(n,m,iammonium) + soil_n_min(n,m,initrate)+K_BNF)) )
              BNF_clover(n,m) =  BNF_clover(n,m) * dt
            END DO ! npts
          END IF ! method selection
       END IF
    END DO

  END SUBROUTINE calc_BNF_clover


!! ================================================================================================================================
!! FUNCTION     : control_temp_func
!!
!>\BRIEF        Calculate temperature control for litter and soild C decomposition
!!
!! DESCRIPTION  : Calculate temperature control factor applied
!! to litter decomposition and to soil carbon decomposition in
!! stomate_soilcarbon.f90 using the following equation: \n
!! \latexonly
!! \input{control_temp_func1.tex}
!! \endlatexonly
!! \n
!! with T the temperature control factor, temp the temperature in Kelvin of 
!! the air (for aboveground litter) or of the soil (for belowground litter 
!! and soil)
!! Then, the function is limited in its maximal range to 1:\n
!! \latexonly
!! \input{control_temp_func2.tex}
!! \endlatexonly
!! \n
!! RECENT CHANGE(S) : None
!!
!! RETURN VALUE: ::tempfunc_result
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART : None
!! \n
!_ ================================================================================================================================

  FUNCTION control_temp_func (npts, temp_in) RESULT (tempfunc_result)

  !! 0. Variable and parameter declaration
    
    !! 0.1 Input variables
    INTEGER, INTENT(in)                 :: npts            !! Domain size - number of land pixels (unitless)
    REAL, DIMENSION(npts), INTENT(in)   :: temp_in         !! Temperature (K)

    !! 0.2 Output variables
    REAL, DIMENSION(npts)               :: tempfunc_result !! Temperature control factor (0-1, unitless)

    !! 0.3 Modified variables

    !! 0.4 Local variables

!_ ================================================================================================================================

    tempfunc_result(:) = exp( soil_Q10 * ( temp_in(:) - (ZeroCelsius+tsoil_ref)) / Q10 )
    tempfunc_result(:) = MIN( un, tempfunc_result(:) )

  END FUNCTION control_temp_func


END MODULE stomate_som_dynamics
