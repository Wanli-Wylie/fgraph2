! =================================================================================================================================
! MODULE       : stomate_phosphorus
!
! CONTACT      : orchidee-help _at_ ipsl.jussieu.fr
!
! LICENCE      : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF       Gathers the main elements for the cycling of phosphorus
!!	
!!\n DESCRIPTION: None
!!
!! RECENT CHANGE(S): None
!!
!! SVN :
!! $HeadURL: $
!! $Date: $
!! $Revision: $
!! \n
!================================================================================================================================


MODULE stomate_phosphorus

  ! modules used:
  USE xios_orchidee
  USE constantes
  USE stomate_data
  USE ioipsl_para
  USE ioipsl
  USE constantes_soil
  USE ieee_arithmetic 

  IMPLICIT NONE

  ! private & public routines
  PRIVATE

  PUBLIC check_mass, cons_mass, &
         phosphorus_dynamics_clear,phosphorus_dynamics, &
         phosphorus_weathering, phosphorus_init, &
         root_conductivity, calc_f_XYplant, &
         root_surface_concentration, enhanced_weathering, &
         phosphorus_freundlich

  INTERFACE cons_mass
    MODULE PROCEDURE cons_mass_r2d, cons_mass_r3d
  END INTERFACE

  ! variables shared by all subroutines in this module
  LOGICAL, SAVE    :: firstcall_phosphorus = .TRUE.   !! Is this the first call? (true/false)

!$OMP THREADPRIVATE(firstcall_phosphorus)

  REAL, ALLOCATABLE, SAVE, DIMENSION (:)   :: Ea       !! Activation energy for Arrhenius term 
!$OMP THREADPRIVATE(Ea)
  REAL, ALLOCATABLE, SAVE, DIMENSION (:)   :: bS_emp   !! empirical factor to relate 
                                                              !! P release by silicate weathering to runoff
!$OMP THREADPRIVATE(bS_emp)
  REAL, ALLOCATABLE, SAVE, DIMENSION (:)   :: bC_emp   !! empirical factor to relate 
                                                              !! P release by carbonate weathering to runoff
!$OMP THREADPRIVATE(bC_emp)
      
  REAL, ALLOCATABLE, SAVE, DIMENSION (:) :: K_freund         !! Freundlich-coefficient for P soprtion [ L(water)kg-1(soil)]
!$OMP THREADPRIVATE(K_freund)
  REAL, ALLOCATABLE, SAVE, DIMENSION (:) :: n_freund         !! Freundlich exponent for P sorption [ ] 
!$OMP THREADPRIVATE(n_freund)

 !DSG_tau_sorb_usda
  REAL, ALLOCATABLE, SAVE, DIMENSION (:) :: tau_sorb       !! 
!$OMP THREADPRIVATE(tau_sorb)
 !DSG_tau_sorb_usda


CONTAINS 

!! ================================================================================================================================
!!  SUBROUTINE   : phosphorus_dynamics_clear
!!
!>\BRIEF        Set the flag ::firstcall to .TRUE. 
!! 
!! 
!_ ================================================================================================================================
  
  SUBROUTINE phosphorus_dynamics_clear
    firstcall_phosphorus=.TRUE.
  END SUBROUTINE phosphorus_dynamics_clear

!! ================================================================================================================================
!! SUBROUTINE    : phosphorus_dynamics
!!
!>\BRIEF        computes the dynamisc of mineral phosphorus pools 
!! tbc...



  SUBROUTINE  phosphorus_dynamics(npts,dt,                          &
                                  soil_orders,                      &
                                  clay,silt, bulk,                  &
                                  drainage_pft,tmc_pft,swc_pft,     &
                                  max_eau_var,                      &
                                  veget_max,                        &
                                  temp_sol,                         &
                                  p_input_dt,                       &
                                  p_mineralisation,                 &
                                  soil_p_min,                       &
                                  f_Pdissolved,                     &
                                  som,                              &
                                  plant_p_uptake, strong_sorption, leaching, &
                                  biomass,                          &
! JCADD merge CNP GRM
                                  grm_pfert,agri_pfert)
! END JCADD

    INTEGER, INTENT(in)                                         :: npts              !! Domain size (unitless)
    REAL, INTENT(in)                                            :: dt                !! Time step of Stomate [fraction of a day]
    INTEGER, DIMENSION(npts),INTENT(in)                         :: soil_orders       !! dominant USDA soil order of the grid box

    REAL,DIMENSION(npts),INTENT(in)                             :: clay              !! Clay fraction of soil (0-1, unitless)
    REAL,DIMENSION(npts),INTENT(in)                             :: silt              !! Silt fraction of soil (0-1, unitless) 
    REAL,DIMENSION(npts),INTENT(in)                             :: bulk              !! Bulk density (kg/m**3) 
    
    ! 0.1 input - climatic drivers
    REAL, DIMENSION(npts,nvm), INTENT(in)                       :: drainage_pft      !! soil water lost by drainage     [kg/m2/timestep]
    REAL, DIMENSION(npts,nvm), INTENT(in)                       :: tmc_pft           !! total soilwater contetn [kg/m2]
    REAL, DIMENSION(npts,nvm), INTENT(in)                       :: swc_pft           !! relative soilwater contetn [m3/m3]
    REAL,DIMENSION (npts)    , INTENT(in)                       :: max_eau_var       !! Maximum water content of the soil   
                                                                                            !! @tex ($kg m^{-2}$) @endtex 
    REAL, DIMENSION(npts,nvm), INTENT(in)                       :: veget_max         !! fraction of a vegetation (0-1)
    REAL, DIMENSION(npts)    , INTENT(in)                       :: temp_sol          !! soil temperature [K]

    ! 0.1 input - phosphorus state variables
    REAL, DIMENSION(npts,nvm), INTENT(in)                       :: p_mineralisation ! net phosphorus mineralisation of decomposing SOM
    REAL, DIMENSION(npts,nvm), INTENT(in)                       :: p_input_dt        !! P released from chemical weathering of rocks (gP/m**2/timestep)
                                                                                            !! P from deposition (gP/m**2/timestep)
                                                                                            !! P from fertilizer (gP/m**2/timestep)
    REAL, DIMENSION(npts,nvm,nparts,nelements), INTENT(in)      :: biomass           !! Biomass pools (gC(or N or P)/m**2)

    ! 0.1 input grass land managemant
! JCADD merge CNP GRM
    REAL, DIMENSION(npts,nvm),INTENT(in)       :: grm_pfert                          !! DESCRIPTION MISSING
    REAL, DIMENSION(npts,nvm),INTENT(in)       :: agri_pfert                         !! DESCRIPTION MISSING
! END JCADD

    ! 0.2 output - phosphorus state variables 
    REAL, DIMENSION(npts,nvm,npspec), INTENT(inout)             :: soil_p_min        !! mineral phosphorus in the soil (gP/m**2) 
                                                                                            !!(first index=npts, second index=nvm, third index=npspec
    REAL, DIMENSION(npts,nvm)       , INTENT(inout)             :: f_Pdissolved      !! fraction of dissolved P in contact with roots [ ]
                                                                                            !! (first index=npts, second index=nvm ) 
    REAL, DIMENSION(npts,ncarb,nvm,nelements), INTENT(inout)    :: som               ! SOM (gC(or N or P)/m**2)

    REAL, DIMENSION(npts,nvm), INTENT(out)                      :: plant_p_uptake    !! Uptake of soil P by plants 
    REAL, DIMENSION(npts,nvm), INTENT(out)                      :: strong_sorption   !! P losses to strong soprtion
    REAL, DIMENSION(npts,nvm), INTENT(out)                      :: leaching          !! mineral phosphorus in solution leached from the soil 


    REAL, DIMENSION(npts,nvm,ncarb)                            :: BC_mineralisation  !! mineralisation flux via phosphatases (gP/m**2/timestep)


    ! 0.3 locals

    INTEGER                              :: m                           !! Index to loop over the nvm PFT's
    INTEGER                              :: n                           !! Index to loop over the ncarb SOM 

    REAL, DIMENSION(npts,nvm)            :: root_p_uptake               !! Uptake of soil P per root biomass
    REAL, DIMENSION(npts,nvm)            :: p_immobilisation            !! P immobilized (gP/m**2/day)

    REAL, DIMENSION(npts)                :: ft_uptake                   !! Temperature response of N uptake by plants (-)


    REAL, DIMENSION(npts)                :: lab_p                       !! Labile phosphorus in plants (gP/m**2) 
    REAL, DIMENSION(npts)                :: lab_n                       !! Labile nitrogen in plants   (gN/m**2) 

    REAL, DIMENSION(nvm)                 :: pn_leaf_min                 !! Minimal PN ratio of leaf (gP / gN)
    REAL, DIMENSION(nvm)                 :: pn_leaf_max                 !! Maximal PN ratio of leaf (gP / gN)

    REAL, DIMENSION(npts,ncarb)          :: f_PNsom                     !! Response of biochemical mineralisation to P content of SOM (deactivated)
    REAL, DIMENSION(npts,ncarb)          :: f_Pmin                      !! Response of biochemical mineralisation to dissolved P concentration
    
    REAL, DIMENSION(npts,ncarb)          :: PNsom                       !! PN ratio of som             (gP / gN)
                                                                               
    REAL, DIMENSION(npts,nvm)            :: f_PNplant                   !!  Response function for phosphorus uptake by plants 
                                                                               !!  to P/N ratio of the labile pool
    REAL, DIMENSION(npts,nvm)            :: f_PCplant                   !!  Response function for biochemical mineralization 
                                                                               !!  to P/C ratio of the labile pool
                                      
    
    REAL, DIMENSION(npts,nvm)            :: soil_p_conc                 !! Concentration of P in soil solution

    REAL                                 :: conv_fac_vmax               !! Conversion from (umol (gDW)-1 h-1) to (g(P) (gC)-1 timestep-1)
    REAL, DIMENSION(npts)                :: conv_fac_concent            !! Conversion factor from (umol per litter) to (g(P) m-2) 

    !
    REAL, DIMENSION(npts,nvm,nelements)  :: mass_before                 !! Biomass pools (gC(or N or P)/m**2)
    REAL, DIMENSION(npts,nvm,nelements)  :: mass_change                 !! Biomass pools (gC(or N or P)/m**2)
    !BCM
    REAL, DIMENSION(ncarb)         :: bcm_turn                    !! Residence time in SOM pools (days^-1)

    REAL, DIMENSION(npts,nvm)      :: labile_turn                 !! Turnover of dissolved P [days^{-1}] as defined by Helfenstein et al. 2018:
                                                                               !! only considering fluxes between dissolved and sorbed

    REAL, DIMENSION(npts,nvm)            :: f_drain                     !! fraction of tmc which has been drained in the last time step 

    REAL, DIMENSION(npts,nvm)                      :: tmpvar ! Temporary var 
    REAL, DIMENSION(npts,nvm)                      :: tmppool ! Temporary var for pool

    REAL, DIMENSION(npts,nvm)                      :: losses ! Temporary var
    REAL, DIMENSION(npts,nvm)                      :: mass_after        !! Temporary variable 
    REAL, DIMENSION(npts,nvm)                      :: P_added_cons !! Temporary variable for mass conservation


    !! printlev is the level of diagnostic information, 0 (none) to 4 (full)
    IF (printlev.GE.3) WRITE(numout,*) 'Entering phosphorus_dynamics' , printlev.GE.3
    IF (printlev.GE.3) WRITE(numout,*) 'firstcall_phosphorus' , firstcall_phosphorus

    !!  ===============================
    !!  0.1 Initialization of variables
    !!  ===============================
    tmpvar            = zero
    tmppool           = zero
    losses            = zero
    mass_after        = zero
    f_drain           = zero
    bcm_turn          = zero
    mass_change       = zero
    mass_before       = zero
    conv_fac_concent  = zero
    conv_fac_vmax     = zero
    soil_p_conc       = zero
    f_PCplant         = zero
    f_PNplant         = zero
    PNsom             = zero
    f_Pmin            = zero
    f_PNsom           = zero
    pn_leaf_max       = zero
    pn_leaf_min       = zero
    lab_n             = zero
    lab_p             = zero
    ft_uptake         = zero
    p_immobilisation  = zero
    root_p_uptake     = zero
    BC_mineralisation = zero
    P_added_cons      = zero
    leaching          = zero
    strong_sorption   = zero
    plant_p_uptake    = zero
    ! set these fluxes to zero, to be able to easily shift the sequence of P
    ! consuming processes 
    BC_mineralisation(:,:,:) = zero
    leaching(:,:)            = zero
    plant_p_uptake(:,:)      = zero

    IF ( firstcall_phosphorus) THEN

       WRITE(numout,*) 'initialize phosphorus variables'
       ! Get parameter values:
       CALL phosphorus_init

       ! Check if the mineral P pools were NOT read in from restart    
       ! if so initialize them assuming high P availability; 
       ! this has to be done using the sorption routine
       ! to ensure both pools are in chemical equilibrium 
       IF ((ALL(soil_p_min(:,:,ipdissolved) == undef)) .OR. &
          ( ALL(soil_p_min(:,:,ipsorbed)    == undef ))) THEN

          CALL phosphorus_freundlich(npts, dt,                    &
                                     soil_orders,                 &
                                     bulk,                        &
                                     max_eau_var,                 &
                                     soil_p_min(:,:,ipdissolved), &
                                     soil_p_min(:,:,ipsorbed),    &
                                     strong_sorption(:,:),        &
                                     P_added_cons(:,:))

          WRITE(numout,*) ' values for desorbed(labile) and sorbed P are initialized from scratch'
       END IF
       firstcall_phosphorus = .FALSE.
    ENDIF

    ! === DSG mass conservation ======================
    mass_before(:,:,1) = SUM(soil_p_min(:,:,1:2),DIM=3)

    !!  ================================
    !!  1.0 Immobilisation of phosphorus
    !!  ================================
    ! this follows the consideration of N immobilisation (Zaehle et al. 2010):
    ! immobilisation has absolute priority to avoid mass conservation problems
    ! the code in litter and soilcarbon has to make sure that dissolved P can never be 
    ! more than exhausted completely by immobilisation!!! 

     IF(printlev.GE.3)THEN
        WRITE(numout,*) 'CHECK values before mineralisation'
        WRITE(numout,*) 'P dissolved ',soil_p_min(test_grid,test_pft,ipdissolved)
        WRITE(numout,*) 'mineralisation ',p_mineralisation(test_grid,:)
     ENDIF

     WHERE(p_mineralisation(:,:).LT.zero)
        p_immobilisation(:,:) = - p_mineralisation(:,:)
     ELSEWHERE
        p_immobilisation(:,:) = zero
     ENDWHERE


     IF(printlev.GE.3)THEN
        WRITE(numout,*) 'CHECK values after mineralisation'
        WRITE(numout,*) 'P dissolved ',soil_p_min(test_grid,test_pft,ipdissolved)
        WRITE(numout,*) 'immob',p_immobilisation(test_grid,:)
     ENDIF
   

    !!  ===========================================
    !!  2.0 Leaching of dissolved labile phosphorus
    !!  ===========================================
    ! this follows the consideration of N leaching (Zaehle et al. 2010):
    !! DSG: There is a comment in stomate_soilcarbon.f90 which states leaching of
    !! nitrogen should be handled in sechiba (section 5. Update pools); 
    !! P leaching should be handled consistent with nitrogen leaching 
    !! DSGseq: as long as we do not have a soil resistance, I leave leaching
    !before plant uptake

    WHERE((tmc_pft(:,:)-drainage_pft(:,:)) .NE. zero) 
        f_drain(:,:) = drainage_pft(:,:)/(tmc_pft(:,:)-drainage_pft(:,:))
    ELSEWHERE
        f_drain(:,:) = zero
    ENDWHERE

    DO m = 1, nvm
        leaching(:,m) = MIN((soil_p_min(:,m,ipdissolved) -plant_p_uptake(:,m))&
                       * f_drain(:,m), soil_p_min(:,m,ipdissolved) )
    ENDDO

    !!  ================================================
    !!  3.0 Plant uptake & Biochemical mineralisation(*)
    !!  ================================================
    ! (*)biochemical mineralization is the hyrdolization of phosphodiester-bonded
    ! P via phosphatses.


    ! 3.1 Scaling functions:
    !!  ====================

    ! 3.1.1 Temperature response function:
    ! Comment from OCN : Temperature function of uptake is similar to SOM decomposition
    ! to avoid P accumulation at low temperatures:
    ! DSG: we also use this for the biochemical mineralization for consistency

    ft_uptake(:) = MAX(control_temp_func (npts, temp_sol(:)), zero )
    
    !We suppress this T control on root uptake due to lack of evidence of such a
    !control:
    ft_uptake(:) = un

    ! 3.1.2 Nutrition status response function
    ! Biochemical mineralisation (BCM) is energetically expensive and usually BCM rates are close
    ! to zero under sufficient P availability

    ! Call the routine which calculate the nutrition status response function
    CALL calc_f_XYplant(npts,nvm,  &
                        biomass,   &
                      !DSG xtrawide  iphosphorus,& ! = X
                      !DSG xtrawide  initrogen, &  ! = Y
                        initrogen, &  ! = Y
                        iphosphorus,& ! = X
                      !DSG 
                        f_PNplant  &
                         )


    DO m = 2, nvm

       
       ! ===============================
       ! 3.2 Biochemical mineralisation 
       ! ===============================
       ! this is more or less following Goll et al. (2012) assuming no BCM in
       ! case vegetation has enough P and maximal in case N/P is max.

       ! the turnover parameters have to be calibrated later checking the
       ! latitudinal pattern of SOM N:P ratios
       ! this is completely made up here: 
        bcm_turn(isurface) = bcm_turn_isurface / one_year    
        bcm_turn(iactive)  = bcm_turn_iactive  / one_year    
        bcm_turn(islow)    = bcm_turn_islow    / one_year            
        bcm_turn(ipassive) = bcm_turn_ipassive / one_year    

       DO n = 1,ncarb ! from iactive to ipassive calculate the BCmineralisation flux
                  ! DSG: should we omitt the passive pool here?
           BC_mineralisation(:,m,n) = som(:,n,m,iphosphorus)*bcm_turn(n)*dt &
           ! the labile P function doesn't really work due to sorpotion etc
           !   * MAX(f_Pmin(:,n),f_PNplant(:,m)) * ft_uptake(:) ! scaling functions 
              * f_PNplant(:,m) * ft_uptake(:) ! scaling functions 
           !  

          ! update som stocks
          som(:,n,m,iphosphorus) = som(:,n,m,iphosphorus) - BC_mineralisation(:,m,n)    
       ENDDO
    ENDDO

    ! ======================================
    ! 3.3 Dissolved labile P uptake by roots
    ! ======================================
    ! Here we assume that concentration of dissolved P in contact to the root
    ! surface cannot be replenished fast enought to sustain a homogenous P
    ! concentration in soil solution

    ! 3.3.1 Calculate the dissolved labile P concentration at the root surface
    CALL root_surface_concentration(npts,nvm,dt,                                            &  
                                    soil_orders,                                            &
                                    clay,silt, bulk,                                        &
                                    biomass(:,:,iroot,icarbon),soil_p_min(:,:,ipdissolved), &
                                    tmc_pft(:,:), swc_pft(:,:), f_Pdissolved(:,:))


    ! one should use totwater instead of max_eau_var(:) in the next subroutine:
    ! The problem is that his leads to huge uptake rates when totwater-> zero,
    ! because we do not account for the water mass flow needed to uptake P
    ! thus I use soil volume instead.

    ! tmpvar for comput. efficiency:
    tmpvar = f_Pdissolved*soil_p_min(:,:,ipdissolved)
    
    ! 3.3.2 Calculate root uptake capacity
    CALL root_conductivity(npts,nvm,                                   &  
                           tmpvar,max_eau_var(:), & 
                           vmax_P_uptake, K_P_min,                     &
                           low_K_P_min, 31.,                           &
                           dt,                                         &
                           root_p_uptake(:,:))

    ! 3.3.3 Calculate the plant nutrition scaling function for P uptake
    CALL calc_f_XYplant(npts,nvm,  &
                        biomass,   &
                      !DSG xtrawide  iphosphorus,& ! = X
                      !DSG xtrawide  initrogen, &  ! = Y
                        initrogen, &  ! = Y
                        iphosphorus,& ! = X
                      !DSG 
                        f_PNplant  &
                         )

    ! 3.3.4 Calculate actual root P uptake
    ! multiply plant mass with the uptake capacity of roots and scale with P demand and temperature
    ! DSG remark: we should NOT scale it with PC instead of PN because 
    ! of the two step procedure of stoichiometric flexibility: first we adjust
    ! C:nutrient, in a second step we adjust N:P 
    plant_p_uptake(:,:) = root_p_uptake(:,:)  * biomass(:,:,iroot,icarbon) *SPREAD(ft_uptake(:),NCOPIES=nvm,DIM=2 )* f_PNplant(:,:)

    !DSG: starve  for debugging:
    ! plant_p_uptake(:,:) = zero
    !DSG: starve 

    ! ensure soil P min doesn't get depleted
    plant_p_uptake(:,:) = MIN(soil_p_min(:,:,ipdissolved)-leaching(:,:),plant_p_uptake(:,:))
    
    ! 3.3.5 update fraction of dissolved P in root contact
    WHERE (soil_p_min(:,:,ipdissolved).GT.zero) 
      f_Pdissolved(:,:)   =  (soil_p_min(:,:,ipdissolved)*f_Pdissolved(:,:)-plant_p_uptake(:,:)) &
                               /soil_p_min(:,:,ipdissolved)
    ELSEWHERE
      f_Pdissolved(:,:)   =  un
    ENDWHERE

    IF(printlev.GE.3)THEN
       WRITE(numout,*) 'CHECK values in phosphorus dynamics'
       WRITE(numout,*) 'root_p_uptake  ',root_p_uptake(test_grid,test_pft)
       WRITE(numout,*) 'vmax_P_uptake  ',vmax_P_uptake
       WRITE(numout,*) 'low_K_P_min',low_K_P_min
       WRITE(numout,*) 'soil_p_min(test_grid,test_pft,idissolved)',soil_p_min(test_grid,test_pft,ipdissolved)
       WRITE(numout,*) 'biomass(root)   ',biomass(test_grid,test_pft,iroot,icarbon)
       WRITE(numout,*) 'ft_uptake       ',ft_uptake(test_grid)
       WRITE(numout,*) 'f_PNplant       ',f_PNplant(test_grid,test_pft)
       WRITE(numout,*) 'BC_mineralisation(test_grid,test_pft,SUM)',SUM(BC_mineralisation(test_grid,test_pft,:),DIM=1)

    ENDIF

    IF(printlev.GE.3)THEN
        WRITE(numout,*) 'before update pools'
        WRITE(numout,*) 'soil_p_min(test_grid,test_pft,idissolved)',soil_p_min(test_grid,test_pft,ipdissolved)
        WRITE(numout,*) 'soil_p_min(test_grid,test_pft,isorbed)',soil_p_min(test_grid,test_pft,ipsorbed)
        WRITE(numout,*) 'p_input_dt(test_grid)',p_input_dt(test_grid,test_pft)
        WRITE(numout,*) 'p_mineralisation(test_grid,test_pft)',p_mineralisation(test_grid,test_pft)
        WRITE(numout,*) 'plant_p_uptake(test_grid,test_pft)',-plant_p_uptake(test_grid,test_pft)
        WRITE(numout,*) 'leaching(test_grid,test_pft)',  -leaching(test_grid,test_pft)
    ENDIF

    ! ======================================
    ! 4.0 Update labile P pools
    ! ======================================
    
    ! gains 
    tmpvar = p_input_dt + p_mineralisation                    & ! gains continue
                        + SUM(BC_mineralisation(:,:,:),DIM=3) & ! gains continue
                        + grm_pfert(:,:)*dt + agri_pfert(:,:)*dt
    ! losses
    losses = leaching + plant_p_uptake

    ! we diagnose the turnover of the labile P pool 
    !(only exchange between sorbed and dissolved is considered)
    ! part 1: remember sorbed pool
    tmppool = soil_p_min(:,:,ipsorbed)

    CALL phosphorus_freundlich(npts,dt,                             &
                             soil_orders,                           &
                             bulk,                                  &
                             max_eau_var,                           &
                             soil_p_min(:,:,ipdissolved),           & ! desorbed P
                             soil_p_min(:,:,ipsorbed),              & ! adsorbed P
                             strong_sorption(:,:),                  & ! losses to strong sorption
                             P_added_cons(:,:),                     &
                             tmpvar,                                & ! gains
                             losses ) ! losses


    ! part 2a: get the change in sorbed pool (w/o considering strongly sorption)
    labile_turn(:,:) = ABS(tmppool - (soil_p_min(:,:,ipsorbed) + strong_sorption(:,:)))
    ! part 2b: divide the dissolved pool by the flux between dissolved and
    ! sorbed to get the turnover as define by Helfenstein et al. 2018
    WHERE ( labile_turn(:,:) .GT. zero )
         labile_turn(:,:) = soil_p_min(:,:,ipdissolved) / labile_turn(:,:)
    ELSEWHERE
         ! still to be confirmed it is 'undef' here to get a missing value in netCDF output:
         labile_turn(:,:) = undef
    ENDWHERE

    ! ======================================
    ! 5.0 Write Output
    ! ======================================

    CALL histwrite_p (hist_id_stomate, 'LABILE_TURN', itime, &
         labile_turn(:,:)/dt,       npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'P_LEACHING', itime, &
         leaching(:,:)/dt,       npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'P_UPTAKE', itime, &
         plant_p_uptake(:,:)/dt, npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'P_MINERALISATION', itime, &
          p_mineralisation(:,:)/dt, npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'BC_MINERALISATION', itime, &
          SUM(BC_mineralisation(:,:,:)/dt,DIM=3), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'P_SSORPTION', itime, &
          strong_sorption(:,:)/dt, npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'SOIL_P_DISSOLVED', itime, &
         soil_p_min(:,:,ipdissolved), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'SOIL_P_SORBED', itime, &
         soil_p_min(:,:,ipsorbed), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'AGRI_FERT_P', itime, &
         agri_pfert(:,:), npts*nvm, horipft_index)

    CALL xios_orchidee_send_field('LABILE_TURN',labile_turn(:,:)/dt)
    CALL xios_orchidee_send_field('P_LEACHING',leaching(:,:)/dt)
    CALL xios_orchidee_send_field('P_UPTAKE',plant_p_uptake(:,:)/dt)
    CALL xios_orchidee_send_field('P_MINERALISATION',p_mineralisation(:,:)/dt)
    CALL xios_orchidee_send_field('BC_MINERALISATION',SUM(BC_mineralisation(:,:,:)/dt,DIM=3))
    CALL xios_orchidee_send_field('P_SSORPTION',strong_sorption(:,:)/dt)
    CALL xios_orchidee_send_field('SOIL_P_DISSOLVED',soil_p_min(:,:,ipdissolved))
    CALL xios_orchidee_send_field('SOIL_P_SORBED',soil_p_min(:,:,ipsorbed))
    CALL xios_orchidee_send_field("AGRI_FERT_P",agri_pfert(:,:))

    ! ======================================
    ! 6.0 Final Mass conservation check
    ! =================================

    ! === DSG mass conservation ======
    IF(printlev.GE.3)THEN
       WRITE(numout,*) 'END of phosphorus_dynamics'
       WRITE(numout,*) 'soil_p_min(test_grid,test_pft,idissolved)',soil_p_min(test_grid,test_pft,ipdissolved)
       WRITE(numout,*) 'soil_p_min(test_grid,test_pft,isorbed)',soil_p_min(test_grid,test_pft,ipsorbed)
       WRITE(numout,*) 'p_input_dt(test_grid)',p_input_dt(test_grid,test_pft)
       WRITE(numout,*) 'p_mineralisation(test_grid,test_pft)',p_mineralisation(test_grid,test_pft)
       WRITE(numout,*) 'plant_p_uptake(test_grid,test_pft)',-plant_p_uptake(test_grid,test_pft)
       WRITE(numout,*) 'leaching(test_grid,test_pft)',  -leaching(test_grid,test_pft)
       WRITE(numout,*) 'strong_sorption(test_grid,test_pft)',  -strong_sorption(test_grid,test_pft)
    ENDIF
    ! =================================================
    mass_change(:,:,1)      = p_input_dt(:,:)  &
                               + p_mineralisation(:,:)+P_added_cons(:,:)  &
                               !BCM
                               + SUM(BC_mineralisation(:,:,:),DIM=3)  & 
                               + grm_pfert(:,:)*dt + agri_pfert(:,:)*dt & 
                               !BCM
                               - plant_p_uptake(:,:) - leaching(:,:)  &
                               - strong_sorption(:,:)

    IF (.NOT.(spinup_analytic)) THEN
       mass_after = SUM(soil_p_min(:,:,1:2),DIM=3) ! POSSIBLE ERROR HERE?
       CALL cons_mass(mass_before(:,:,1),             &  ! mass before
             mass_after(:,:),              &  ! mass after
             mass_change(:,:,1)              &  ! net of fluxes
            )
    ENDIF
    ! =================================================

  END SUBROUTINE phosphorus_dynamics   

  !! -- Update "mineral" Pools using Freundlich Isotherme
  !     Most sophisticated approach for labile P sorption which is currently implemented 
  !     
  !     
  !     

  SUBROUTINE phosphorus_freundlich(npts, dt,      &
                                 soil_orders,     &
                                 bulk,            &
                                 max_eau_var,     &
                                 desorbed,        &
                                 adsorbed,        &
                                 strong_sorption, &
                                 P_added_cons,    &
                                 gains,           &
                                 losses, NoStrongSorption)
    !
    INTEGER, INTENT(in)                                :: npts              !! Domain size (unitless)
    REAL, INTENT(in)                                   :: dt                !! Time step [fraction of a day]


!DSGrefine
    REAL,DIMENSION(npts),INTENT(in)                    :: bulk              !! Bulk density (kg/m**3)  
    REAL, DIMENSION(npts),INTENT(in)                   :: max_eau_var       !! maximum water holding capacity of soils  [kg m-2]
!DSGrefine
    
    ! input
    REAL, DIMENSION(npts,nvm), INTENT(in) , OPTIONAL   :: gains             !! gain flux     [gP/m2/timestep]
    REAL, DIMENSION(npts,nvm), INTENT(in) , OPTIONAL   :: losses            !! loss flux     [gP/m2/timestep]
    REAL, DIMENSION(npts,nvm), INTENT(in) , OPTIONAL   :: NoStrongSorption  !! we compute no strong sorption losses 
                                                                                   !! as this is the 2nd  call of this routine within the same timestep 
                                                                                   

    INTEGER, DIMENSION(npts), INTENT(in)               :: soil_orders       !! dominant USDA soil order

    ! output
    REAL, DIMENSION(npts,nvm), INTENT(inout)           :: desorbed         !! mineral nitrogen in desorbed form (gP/m2) 
                                                                                  !!(first index=npts, second index=nvm, third index=nnspec
    REAL, DIMENSION(npts,nvm), INTENT(inout)           :: adsorbed         !! mineral nitrogen in adsorbed form (gP/m2) 
    REAL, DIMENSION(npts,nvm), INTENT(out)             :: strong_sorption  !! the flux of sorbed into strongyl sorbed form (gP/m2)
                                                                                  !!(first index=npts, second index=nvm, third index=nnspec
    REAL, DIMENSION(npts,nvm), INTENT(out)             :: P_added_cons
    ! locals
    INTEGER                                            :: l                !! Index to loop over the npts grid boxes
    INTEGER                                            :: m                !! Index to loop over the nvm PFT's
    INTEGER                                            :: n                !! Index to loop over the p_nscm soil orders
    
    REAL, DIMENSION(npts,nvm)                          :: distributor      !! ratio which distributes losses and gains among dissolved and sorbed P
    REAL, DIMENSION(npts,nvm)                          :: deficit          !! the deficit of the dissolved pool, in case losses fluxes exceed supply 
    REAL, DIMENSION(npts,nvm)                          :: labileMinP       !! the net change of dissolved and adsorbed labile P
    LOGICAL                                                   :: init_pools       !! in case =TRUE: this routine handles the initialization of pools 
                                                                                  !! when not read in from restart files
   
    LOGICAL, DIMENSION(npts,nvm)                              :: P_added          !! we added P to min_p_soil pools to avoid overdepletion 
                                                                                  !! during analytic spinup [true/false]

    REAL, DIMENSION(npts)                              :: conv_des         !! conversion factor [m2 L-1 ] 
    REAL, DIMENSION(npts)                              :: conv_sorb        !! conversion factor [kg1 m-2] 
    REAL, DIMENSION(npts,nvm)                          :: K_freund_grid    !! Freundlich coefficient; units corrected [ ]
    REAL, DIMENSION(npts,nvm)                          :: n_freund_grid    !! Freundlich exponent [ ]
 !DSG_tau_sorb_usda
    REAL, DIMENSION(npts,nvm)                          :: tau_sorb_grid    !! 
 !DSG_tau_sorb_usda

    REAL, DIMENSION(npts,nvm)                          :: original_P !! temporary variables for record zero of negative mineral P
!
    LOGICAL                                                   ::  ssorpt_from_diss ! strongly sorption losses are a function of dissolved pool not sorbed pool
    ssorpt_from_diss=.FALSE.
!

    ! check if this routine handles (1) changes in pools or (2) initialization of pools 
    
    IF (PRESENT(gains) .OR. PRESENT(losses)) THEN
            init_pools=.FALSE.
            IF (.NOT.(PRESENT(gains).AND.PRESENT(losses))) THEN
               CALL ipslerr_p(3, 'phosphorus_freundlich', 'We will STOP now.', &
                  & 'The combination of input varialbes to the subroutine is &
                  invalid.','')
            END IF

    ELSE ! not present gains and losses:
            init_pools=.TRUE.
    END IF

    ! Freundlich Isotherm:
    ! adsorbed = K_F * dissolved**(1/n)

    ! differentiated form:
    ! d adsorbed    dissolved**(1/n-1)  d dissolved
    ! ---------- =  ------------------ ------------
    !     d t               n             d t

    ! DSG: check Goll et al. (2012) SI to see how to derive the dynamics of
    ! dissolved and adsorbed: just substitute the Langmuir equation with the
    ! Freundlich equation

    ! to have it correct one should convert desorbed into mg/l(H2O) and adsorbed
    ! into mg/kg(soil). however currently we lack a proper discretisation of the soil
    ! volume thus use pore space to substitute water volume throughout the
    ! nutrient code. 
    !DSG: the unit conversion factor doesnt work at the moment:
    ! g m-2 -> mg / l    using max_eau_varwhich is in [kg m-2]: / max_eau_var(n)
    conv_des(:)=(un/max_eau_var(:))*1.e3
    !conv_des(:)=un
    ! g m-2 -> mg / kg   using bulk soil density which is in [kg m-3] and 0.3 m depth according to the measurement data :  / (bulk(n)*un)
    conv_sorb(:)=((bulk(:)*0.3))/1.e3
    !conv_sorb(:)=un
    ! initialize P_added_cons
    P_added_cons(:,:) = zero

    ! correct K_freund for the units 
    DO l=1,npts
!DSGtuning: 
      K_freund_grid(l,:)=SPREAD(K_freund(soil_orders(l))*conv_des(l)*conv_sorb(l),NCOPIES=nvm,DIM=1)*sorb_tune
!DSGtuning
      n_freund_grid(l,:)=SPREAD(n_freund(soil_orders(l)),NCOPIES=nvm,DIM=1)
      ! 
      tau_sorb_grid(l,:)=SPREAD(tau_sorb(soil_orders(l)),NCOPIES=nvm,DIM=1)
    ENDDO
    
    IF (init_pools) THEN ! initialize pools in chemical equilbrium
      WRITE(numout,*) ' value of desorbed(dissolved) P is set to [g/m2]= ',dissolved_init
      strong_sorption(:,:) = zero
      desorbed(:,:) = dissolved_init  ! initialize with default value 
      adsorbed(:,:) = (K_freund_grid(:,:) * (desorbed(:,:))**(un/n_freund_grid(:,:))) 
    ELSE ! handle changes in pools
      ! calculate the dissolution/sorbtion ratio
      !  distributor(l,:) = (un + A )**(-1)
      !                 A = K_freund(soil_orders(l)) * (C/n_freund(soil_orders(l)))
      !                 C = desorbed(l,:)**((un/n_freund(soil_orders(l)))-un)

      distributor(:,:) = zero
      ! in the case of desorbed has a very small value
      ! we calculate distributor using a minimum value to avoid numerical
      ! overflow/divided-by-0 problems
      WHERE (desorbed(:,:) .GE. min_stomate)  
        distributor(:,:) = ((un +  K_freund_grid(:,:) *                     &
                              ((desorbed(:,:))**((un/n_freund_grid(:,:))-un) &
                              /n_freund_grid(:,:)))**(-1))  
      ! the distributor can be changed by chelating agents and acids [not
      ! implemented yet]
      ! ====== (see Buendia et al (2014) for possible parametrization)
      ELSEWHERE
        distributor(:,:) = ((un +  K_freund_grid(:,:) *                     &
                              ((min_stomate)**((un/n_freund_grid(:,:))-un) &
                              /n_freund_grid(:,:)))**(-1))
      ENDWHERE

      IF (PRESENT(NoStrongSorption)) THEN
        strong_sorption(:,:) = zero
      ELSEIF (ssorpt_from_diss) THEN
        ! new formulation
        !strong_sorption(:,:) = (un/(tau_sorbed))*dt*desorbed(:,:)
        strong_sorption(:,:) = (un/(tau_sorb_grid(:,:)))*dt*desorbed(:,:)
      ELSE
        ! original formulation:
        WHERE (desorbed(:,:) .GE. min_stomate)  
         ! strong_sorption(:,:) = (un/(tau_sorbed))*dt                         &
         strong_sorption(:,:) = (un/(tau_sorb_grid(:,:)))*dt                         &
                           * K_freund_grid(:,:)  * (desorbed(:,:))**(un/n_freund_grid(:,:))
        ELSEWHERE
          strong_sorption(:,:) = 1e-10
        ENDWHERE
      ENDIF

      ! net change in labile P:
      labileMinP(:,:) =  gains(:,:)  - losses(:,:)     &
                         - strong_sorption(:,:)

      ! update desorbed labile P:
      desorbed(:,:) = desorbed(:,:) + distributor(:,:)* labileMinP(:,:)

      ! update adsorbed labile P:
      adsorbed(:,:) = adsorbed(:,:) + (un - distributor(:,:)) * labileMinP(:,:)

      ! mass conservation has to be ensured:
      ! immobilisation can theoretically exceed supply:
      deficit(:,:) = MIN(desorbed(:,:),zero)                   
      deficit(:,:) = deficit(:,:) + MIN(adsorbed(:,:),zero)                   
    ENDIF 

    IF (printlev>=3) THEN
       WRITE(6,*) 'desorbed(test_grid,test_pft)',desorbed(test_grid,test_pft)
       WRITE(6,*) 'adsorbed(test_grid,test_pft)',adsorbed(test_grid,test_pft)
       WRITE(6,*) 'distributor(test_grid,test_pft)',distributor(test_grid,test_pft)
     ENDIF

    ! During a spinup we can have overdepletion of the mineral P : immob > desorbed
    P_added_cons(:,:) = zero
    IF (spinup_cnp.OR.impose_cn.OR.impose_np) THEN 
      P_added(:,:) = .FALSE.

      DO m=1,nvm 
        WHERE ((desorbed(:,m) .LE. zero)) 
            P_added_cons(:,m) = - desorbed(:,m) + 1.e-10
            desorbed (:,m) = 1.e-10        ! we set the pools to a value close to zero (zero wouldnt work with the distributor calc.)
            WHERE (adsorbed (:,m) .LE. zero)
              P_added_cons(:,m) = P_added_cons(:,m) - adsorbed (:,m) + &
                      (K_freund_grid(:,m) * (desorbed(:,m))**(un/n_freund_grid(:,m))) 
              adsorbed (:,m) = (K_freund_grid(:,m) * (desorbed(:,m))**(un/n_freund_grid(:,m))) 
            ELSEWHERE
              P_added_cons(:,m) = P_added_cons(:,m) - adsorbed (:,m) + &
                      (K_freund_grid(:,m) * (desorbed(:,m))**(un/n_freund_grid(:,m)))
              adsorbed (:,m) = (K_freund_grid(:,m) * (desorbed(:,m))**(un/n_freund_grid(:,m)))
            ENDWHERE
            P_added  (:,m) = .TRUE. ! mark the box
        ELSEWHERE ((adsorbed(:,m) .LT. zero)) !DSGdebug:
            ! this case happens only when desorbed is still positive
            ! we should stop in this case?
            P_added_cons(:,m) = - desorbed(:,m) + 1.e-10
            desorbed (:,m) = 1.e-10
            P_added_cons(:,m) = P_added_cons(:,m) - adsorbed(:,m) + &
                   (K_freund_grid(:,m) * (desorbed(:,m))**(un/n_freund_grid(:,m)))
            adsorbed (:,m) = (K_freund_grid(:,m) * (desorbed(:,m))**(un/n_freund_grid(:,m)))
            P_added  (:,m) = .TRUE. ! mark the box
        ENDWHERE 
      ENDDO

     !DSG IF(printlev>=4)THEN
          IF (ANY(P_added)) THEN 
              WRITE(numout,*) 'we added P to soil_p_min pools to avoid overdepletion'
              WRITE(numout,*) 'in : ', COUNT(P_added)
              WRITE(numout,*) 'of (nvm*npts): ', npts*nvm
          ENDIF
     !DSG ENDIF
    ENDIF

    ! debug
    IF (ANY(desorbed < zero) .OR. ANY(adsorbed < zero)) THEN
        WRITE(6,*) 'FATAL ERROR '
        WRITE(6,*) 'Freundlich equation: there is not enough P to satisfy P demand:'
        DO m=1,npts
           IF(ANY(desorbed(m,:).LT.zero).OR.ANY(adsorbed(m,:).LT.zero)) THEN
             WRITE(6,*) 'm'        ,m
             WRITE(6,*) 'deficit(m,:)'        ,deficit        (m,:)
             WRITE(6,*) 'gains(m,:)'          ,gains          (m,:)
             WRITE(6,*) 'losses(m,:)'         ,losses         (m,:)
             WRITE(6,*) 'strong_sorption(m,:)',strong_sorption(m,:)
             WRITE(6,*) 'distributor(m,:)'    ,distributor    (m,:)
             WRITE(6,*) 'desorbed(m,:)'       ,desorbed       (m,:)
             WRITE(6,*) 'adsorbed(m,:)'       ,adsorbed       (m,:)
           ENDIF
        ENDDO
        STOP 'Freundlich equation: there is not enough P to satisfy P demand'
    END IF

  END SUBROUTINE phosphorus_freundlich

!! ================================================================================================================================
!! SUBROUTINE    : phosphorus_weathering
!!
!>\BRIEF        computes the release of phosphorus from minerals by 
!!              weathering
!!
!! DESCRIPTION   : Here the flux of phosphorus entering the biosphere from
!!                 the chemical weathering of P-bearing minerals is computed. 
!!                 The flux is computed following Goll et al (2014); the model
!!                 is descriped in Hartmann et al. (2014).
!!                 The stock of mineral P is assumed to stay constant; this
!!                 limits the applicabilty of this model to time scales < 10,000 yr.
!!                 The model is driven by air temperature, not by soil
!!                 temperature as air temperature was used to parametrize the
!!                 empirical model. 
!!
!!                 P release (Pr) is calculated:                
!!                 Pr = b * q *f_t * s
!!                 where
!!                 q: annual runoff [mm/a]
!!                 s: factor to correct for soil shielding
!!                 b: empirical factor 
!!                 temperature response term f_t is given by:
!!                 f_t = exp(-Ea/R*(1/T-1/Tref))
!!                 with
!!                 Ea:  apparent activation energy (depends on lithology)
!!                 R: gas constant [J/mol /K] 
!!                 T: annual mean 2m air temperature in [K]
!!                 Tref: 284.15K
!!
!!                 ORCHIDEE specific implementation: the long-term (3yr) average of runoff and temperature are used
!!                 to get smooth changes; original time step of model is annual (Goll et al.,2014)
!!
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): :: p_release
!!
!! REFERENCE(S)   : 
!! - Hartmann, J., N. Moosdorf, R. Sauerwald, M. Hinderer and J. West. 2014.
!! Global chemical weathering and associated P-release — The role of lithology, 
!! temperature and soil properties", Chemical Geology, 363, 145-163
!! - Goll , D.S, Moosdorf, N., Hartmann, J., and Brovkin, V. 2013.
!! Climate-driven changes in chemical weathering and associated phosphorus
!! release since 1850 : Implications for the land carbon balance,
!! DOI: 10.1002/2014GL059471
!!
!! FLOWCHART    : 
!! \latexonly 
!! \includegraphics[scale=0.5]{clearcutflow.jpg}
!! \endlatexonly
!! \n
!_ ================================================================================================================================



  SUBROUTINE  phosphorus_weathering(npts,                                               &
                                    lith_frac, soil_shield,                             &
                                    drainage_longterm,runoff_longterm,t2m_longterm,     &
                                    p_release)


    INTEGER, INTENT(in)                            :: npts              !! Domain size (unitless)
    ! lithological & soil information
    REAL,  DIMENSION(npts,nlcm), INTENT(in)        :: lith_frac         !! area fractions of lithologies [0-1]
    REAL,  DIMENSION(npts)     , INTENT(in)        :: soil_shield       !! factor to correct P weathering for soil shielding effects [0-1]

    ! climatic drivers
    REAL, DIMENSION(npts), INTENT(in)              :: drainage_longterm !! annual sum drainage        (long term average)     (kg/m2/year) 
    REAL, DIMENSION(npts), INTENT(in)              :: runoff_longterm   !! annual sum drainage        (long term average)     (kg/m2/year)
    REAL, DIMENSION(npts), INTENT(in)              :: t2m_longterm      !! average 2m air temperature (long term average)     (K)

    ! output
    REAL, DIMENSION(npts), INTENT(out)             :: p_release         !! P release from chemical weathering of rocks (gP/m**2/day)

    ! locals

    REAL, PARAMETER                 :: Tref=284.15       !! refernce temperature (mean annual temperatuer Japan) [K]
    REAL, DIMENSION(npts,nlcm)      :: Ftemp             !! temperature response function 
    REAL, DIMENSION(npts,nlcm)      :: prelease         !! P release from chemical weathering of rocks  of each lithological fraction(gP/m**2/timestep)

    REAL, DIMENSION(npts)           :: tmp_shield
    INTEGER                         :: dsg1

    IF (printlev>=3) THEN
       WRITE(numout,*) 'Entering phosphorus_weathering'
       WRITE(numout,*) 'Ea(:);',Ea(:)
       WRITE(numout,*) 'drainage_longterm(:);',drainage_longterm(test_grid)
       WRITE(numout,*) 'runoff_longterm(:);',runoff_longterm(test_grid)
       WRITE(numout,*) 't2m_longterm(:)',t2m_longterm(test_grid)
       WRITE(numout,*) 'lith_frac(:);',lith_frac(test_grid,:)
       WRITE(numout,*) 'soil_shield(:);',soil_shield(test_grid)
    ENDIF

    ! just to make sure
    prelease(:,:) = zero

   ! 1. temperature dependence of chemical weathering; Arrhenius
    DO dsg1=1,npts
        WHERE (Ea(:) .NE. undef)  
            Ftemp(dsg1,:) = exp(-Ea(:)/RR*(1/t2m_longterm(dsg1)-1/Tref))
        ELSEWHERE
            Ftemp(dsg1,:) = un ! no temperature effect
        ENDWHERE
    ENDDO

    DO dsg1=1,npts
        ! if there is water running
        IF (drainage_longterm(dsg1)+runoff_longterm(dsg1).GT. min_stomate) THEN  
           ! 2.1 release of P from silicate weathering [g/m2/timestep]:
           WHERE (bS_emp(:).NE.undef)
               prelease(dsg1,:) = MAX(zero,                                          &
                                  Ftemp(dsg1,:) * bS_emp(:) * lith_frac(dsg1,:)      & 
                                  * (drainage_longterm(dsg1) + runoff_longterm(dsg1)))
           ELSEWHERE
               prelease(dsg1,:) = zero   ! silicate-free lithology
           ENDWHERE

           ! 2.2 release of P from carbonate weathering:
           WHERE (bC_emp(:) .NE. undef)
               prelease(dsg1,:) =  prelease(dsg1,:) +                                & ! add P release from carbonate weathering 
                                  MAX(zero,                                          &
                                  bC_emp(:) * lith_frac(dsg1,:)                      & ! no temperature dependence  (Hartmann et al. (2014))
                                  * (drainage_longterm(dsg1) + runoff_longterm(dsg1)))
           ENDWHERE
       ENDIF
    ENDDO

    !  3. Correct for soil shielding effects and sum weathering rates up to the grid box
!DSG modify increase minimum soil_shield to 0.3
! for tropical weathering See mail
    tmp_shield(:) = soil_shield(:)
    ! not needed: tmp_shield(:) = MAX(0.3,tmp_shield(:)) 
    p_release(:) = tmp_shield(:)*SUM(prelease(:,:),DIM=2)/365.

    IF (printlev>=3) WRITE(numout,*) 'Leaving phosphorus_weathering'

  END SUBROUTINE phosphorus_weathering


!! ================================================================================================================================
!! SUBROUTINE    : enhanced weathering
!!
!>\BRIEF        computes the weathering of added minerals according to their size and
!!              type  - 
!!              
!!
!! DESCRIPTION   :  We use the weathering model designed by Strefler et al and
!!                  adapted it to simulate P release besides CO2 consumption
!!                 
!!                 
!!                
!! MAIN OUTPUT VARIABLE(S): :: 
!!
!! REFERENCE(S)   :  Strefler, Jessica, et al. "Potential and costs of Carbon
!!                   Dioxide Removal by Enhanced Weathering of rocks." Environmental Research
!!                   Letters 13.3 (2018): 034010.
!!                   Bandstra J Z and Brantley S L 2008 Data fitting techniques
!!                   with applications to mineral dissolution kinetics Kinetics of Water-Rock
!!                   Interaction ed S L Brantley, J D Kubicki and A F White (New York: Springer) pp
!!                   211–57
!!                   Goll et al (unplished)



  SUBROUTINE  enhanced_weathering(npts,dt,           & !
                                  temp_sol,pH,      & ! drivers
                                  EW_grain, EW_flux ) ! output


    INTEGER, INTENT(in)                                   :: npts              !! Domain size (unitless)
    REAL   , INTENT(in)                                   :: dt                !! time step [fraction of one day]
    ! drivers: pH, temperatures
    REAL, DIMENSION(npts)                , INTENT(in)     :: temp_sol          !! soil temperature [K]
    REAL, DIMENSION(npts)                , INTENT(in)     :: pH                !! soil pH

    ! state variables: grain_size, grain_mass
    REAL, DIMENSION(npts,nvm,nminerals,2), INTENT(inout)  :: EW_grain   !! Mass of minerals [g m-2] (1)/ average diameter of mineral grain [10-6 m] (2)
    REAL, DIMENSION(npts,nvm,nECfluxes)  , INTENT(out)    :: EW_flux    !! The fluxes deriving from weahtering (=1=ico2 consumption, =2=iPr release)

    ! Locals
    REAL, DIMENSION(npts,nvm,nminerals)                   :: SSA           !!  Specific surface area of mineral grains [m2 g-1]
    REAL, DIMENSION(npts,nvm,nminerals)                   :: diss_rate     !!  Dissolution rate of mineral grains [g m-2 yr-1]
    
    REAL, DIMENSION(npts,nminerals)                       :: f_temp        !! Temperature response function [ ]
    REAL, DIMENSION(npts,nminerals)                       :: weath_rate    !! Weathering rate [mol m-2 s-1]

    REAL, DIMENSION(npts,nvm,nminerals)                   :: grain_mass    !! dummy for Mass of minerals [g m-2]
    REAL, DIMENSION(npts,nvm,nminerals)                   :: grain_size    !! dummy for Size of mineral grain [g m-2]
    REAL, DIMENSION(npts,nvm,nminerals,nECfluxes)         :: tmp_EC        !! dummy 

    REAL, PARAMETER                                       :: sec_per_day=86400.

    ! parameters (should be externalized):
    ! # depends on scenario: # depends on mineral: grain_std_atm_weight, Pcontent, CO2cons
    REAL, DIMENSION(nminerals)                            :: saw           !!    standard atomic weight of mineral [ g mol-1]
    REAL, DIMENSION(nminerals)                            :: Ea_EC         !!    activation energy of mineral  [J mol-1]
    REAL, DIMENSION(nminerals)                            :: pH_offset     !!    
    REAL, DIMENSION(nminerals)                            :: k_H           !!    [mol m-2 s-1]
    REAL, DIMENSION(nminerals)                            :: n_H           !!    [  ]
    REAL, DIMENSION(nminerals)                            :: k_OH          !!    [mol m-2 s-1]
    REAL, DIMENSION(nminerals)                            :: n_OH          !!    [  ]
    REAL, DIMENSION(nminerals)                            :: P_content     !! P content of mineral [%]
    REAL, DIMENSION(nminerals)                            :: CO2_cons      !! CO2 consumption of wathering [g(c) g(rock)]
    LOGICAL    , DIMENSION(nminerals)                           :: temp_effect    !! 

    INTEGER  :: m

    IF (printlev>=3) THEN
       WRITE(numout,*) 'Entering enhanced_weathering'
    ENDIF

    ! fill the local variables for grain size and mass
    grain_mass(:,:,:) = EW_grain(:,:,:,1)
    grain_size(:,:,:) = EW_grain(:,:,:,2)
    ! initialize
    EW_flux   = zero

!=================================================================== 
!DSG: this could go  somewhere else as it is static in time
    ! Basalt :
    saw        (ibasalt)=124.
    Ea_EC      (ibasalt)=47500.
    pH_offset  (ibasalt)=14.
    k_H        (ibasalt)=588.
    n_H        (ibasalt)=-1.16
    k_OH       (ibasalt)=0.0822
    n_OH       (ibasalt)=0.16
    temp_effect(ibasalt)=.TRUE.
    ! Value from database  Earthchem DB
    P_content  (ibasalt)=0.00151
    ! in brackets is conversion from CO2 to C
    CO2_cons   (ibasalt)=0.3*(12.0096/(12.0096+2*15.999))

    ! Forsterite : 
    saw        (iforsterite)=140.7
    Ea_EC      (iforsterite)=undef
    pH_offset  (iforsterite)=undef
    k_H        (iforsterite)=5.55e-8
    n_H        (iforsterite)=-0.372 
    k_OH       (iforsterite)=zero
    n_OH       (iforsterite)=undef
    temp_effect(iforsterite)=.FALSE.
    P_content  (iforsterite)=0.0 
    ! in brackets is conversion from CO2 to C
    CO2_cons   (iforsterite)=1.25*(12.0096/(12.0096+2*15.999))

    IF (printlev>=3) THEN
       WRITE(numout,*) 'Enhanced_weathering: calculated specific surface area'
    ENDIF
     
    WHERE (grain_mass(:,:,:).GT.min_stomate)
       ! Specific surface area of grain (Strefler et al 2018):
       SSA(:,:,:) = 69.18 * grain_size(:,:,:)**(-1.24)
    ELSEWHERE
       SSA(:,:,:) = undef
    ENDWHERE
!DSG
!=================================================================== 

   ! Weathering rate (Strefler et al 2018) based on Brandstra & Brantley (2008)
   DO m = 1, npts ! loop over land points
       ! A. temperature dependence of weathering; Arrhenius (RR is gas constant)
       
       ! f_temp range around 10^-9; thus k_H is factor 10^8 smaller than for
       ! minerals w/o temp effect ; not sure this is the right thing to do
       WHERE(.NOT.(Ea_EC(:)==undef))
           f_temp(m,:) = exp(-Ea_EC(:)/(8.314472*temp_sol(m)))
       ELSEWHERE 
           f_temp(m,:) = un
       ENDWHERE

       ! Weathering rate [ mol m-2 s-1]
       WHERE(k_OH(:).GT. zero)
          weath_rate(m,:)   = (k_H(:) * 10**(n_H(:)*pH(m))  + k_OH(:) * 10**(n_OH(:)*(pH(m)-pH_offset(:))) ) * f_temp(m,:)
       ELSEWHERE 
          weath_rate(m,:)   = k_H(:)  * 10**(n_H(:)*pH(m)) * f_temp(m,:) 
       ENDWHERE
   ENDDO

!!DGS tested ! works!
!   WRITE(numout,*) 'Enhanced_weathering: calculated weathering rate1'
!   ! log10(          (a      * 10.^(-b     *pH)  + c       * 10.^(d*(pH-14)))       * exp(-e.      /(f *g          ))                 )
!   tmp_EC(1,1,1,1) = ( 588.* 10**(-1.16*4.)  + 0.0822 *   10**(0.16*(4.-14)) ) * exp(-47500/(8.314472*(298.15)))
!   WRITE(numout,*) 'Enhanced_weathering: calculated weathering rate2'
!   tmp_EC(2,1,1,1) = (k_H(2) * 10**(n_H(2)*7.)  + k_OH(2) * 10**(n_OH(2)*(7.-pH_offset(2))) ) * exp(-Ea_EC(2)/(RR*(25.+273.15)))
!   WRITE(numout,*) 'Enhanced_weathering: calculated weathering rate3'
!   tmp_EC(3,1,1,1) = (k_H(2) * 10**(n_H(2)*9.)  + k_OH(2) * 10**(n_OH(2)*(9.-pH_offset(2))) ) * exp(-Ea_EC(2)/(RR*(25.+273.15)))
!   WRITE(numout,*) 'pH=4: -10.13'
!   WRITE(numout,*) 'pH=4', log10(tmp_EC(1,1,1,1))
!   WRITE(numout,*) 'pH=7: -10.53'
!   WRITE(numout,*) 'pH=7', log10(tmp_EC(2,1,1,1))
!   WRITE(numout,*) 'pH=9: -10.21'
!   WRITE(numout,*) 'pH=9', log10(tmp_EC(3,1,1,1))
!!!DGS test
!!!DSGdbg


   ! Dissolution rate [g m-2 s-1] (Strefler et al. 2018)
   DO m = 1, nvm
       WHERE (grain_mass(:,m,:).GT.min_stomate) 
           diss_rate(:,m,:) = SSA(:,m,:) * weath_rate(:,:) * SPREAD(saw(:),DIM=1,NCOPIES=npts) * grain_mass(:,m,:)
       ELSEWHERE
           diss_rate(:,m,:) = zero
       ENDWHERE
   END DO
  ! WRITE(numout,*) 'diss_rate(:,:,:)',diss_rate(:,:,:)

   ! from g m-2 s-1 to g m-2 day-1
   diss_rate(:,:,:) = diss_rate(:,:,:)*86400.

   ! Update grain mass: g m-2 day-1  to g m-2 timestep-1
   grain_mass(:,:,:) =  grain_mass(:,:,:) - diss_rate(:,:,:)*dt
   
   ! calculate P release and CO2 consumption
   DO m = 1,nminerals
       tmp_EC(:,:,m,ico2cons)  = diss_rate(:,:,m)*CO2_cons(m)
       tmp_EC(:,:,m,iPrelease) = diss_rate(:,:,m)*P_content(m)
   ENDDO
   ! sum over mineral types 
   EW_flux(:,:,iPrelease)  = SUM(tmp_EC(:,:,:,iPrelease),DIM=3)
   EW_flux(:,:,ico2cons)   = SUM(tmp_EC(:,:,:,ico2cons) ,DIM=3)

   ! For now we do not account for changes in grain size, 
   ! we only update grain mass:
   EW_grain(:,:,:,1) = grain_mass(:,:,:)

  ! WRITE(numout,*) 'EW_grain(:,:,:,:)',EW_grain(:,:,:,:)
   

   CALL xios_orchidee_send_field('EW_CO2C'      ,EW_flux(:,:,1))
   CALL xios_orchidee_send_field('EW_PR'        ,EW_flux(:,:,2))
   CALL xios_orchidee_send_field('EW_BASALT'    ,EW_grain(:,:,ibasalt,1))

   IF (printlev>=3) THEN
      WRITE(numout,*) 'Leaving enhanced_weathering'
   ENDIF

  END SUBROUTINE enhanced_weathering

!! ================================================================================================================================
!! SUBROUTINE   : phosphorus_init
!!
!>\BRIEF        Initializations and memory allocation   
!!
!! DESCRIPTION  :
!! - 1 Some initializations
!! - 2 
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S) : 
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
!!_ phosphorus_init

  SUBROUTINE phosphorus_init

    ! interface description
    INTEGER                                     :: ier                   !! Error code



!_ ================================================================================================================================

    ! 1.0 initialisation
    IF (firstcall_phosphorus) THEN 
        firstcall_phosphorus=.FALSE.
    ELSE 
       WRITE (numout,*) ' firstcall_phosphorus false . we stop '
       STOP 'phosphorus_init'
    ENDIF
    
    ! 2.0 allocate memory 
    ! 2.1 weathering parameters:
    ALLOCATE (Ea(nlcm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in Ea allocation. We stop. We need nlcm words = ',nlcm
       STOP 'phosphorus_init'
    END IF

    ALLOCATE (bS_emp(nlcm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in bS_emp allocation. We stop. We need nlcm words = ',nlcm
       STOP 'phosphorus_init'
    END IF
    ALLOCATE (bC_emp(nlcm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in bC_emp allocation. We stop. We need nlcm words = ',nlcm
       STOP 'phosphorus_init'
    END IF
    
    ! Attention: the nscm_usda here correspond to the 12 USDA soil order not to
    ! 12 texture classes; like elsewhere in the code. This inconsistency could
    ! lead to soil order specific parameters which are not consistent with each
    ! other. The handling of the hydroligical parameters should use USDA soil
    ! order and not texture classes to avoid problems. 

    ALLOCATE (K_freund(p_nscm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in K_freund allocation. We stop. We need nlcm words = p_nscm ',p_nscm
       STOP 'phosphorus_init'
    END IF
    ALLOCATE (n_freund(p_nscm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in n_freund allocation. We stop. We need nlcm words = p_nscm ',p_nscm
       STOP 'phosphorus_init'
    END IF
 !DSG_tau_sorb_usda
    ALLOCATE (tau_sorb(p_nscm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in tau_sorb allocation. We stop. We need nlcm words = p_nscm ',p_nscm
       STOP 'phosphorus_init'
    END IF
 !DSG_tau_sorb_usda

    ! 3.0 assign values to the parameter according to the classes
    ! 3.1 Lithological parameter; only nclm=16 supported (GLiM)
    SELECTCASE (nlcm)
    CASE (16)
       Ea(:)     = Ea_GLiM(:)       
       bS_emp(:) = bS_emp_GLiM(:)       
       bC_emp(:) = bC_emp_GLiM(:)       

    CASE DEFAULT
       WRITE (numout,*) 'Unsupported lithological type classification.Currently, only GLiM classification is supported'
       STOP 'phosphorus_init'
    ENDSELECT

    ! 3.2 soil parameter; only nscm=12 supported (USDA)
    SELECTCASE (p_nscm)
    CASE (12)
       K_freund(:) = K_freund_usda(:)
       n_freund(:) = n_freund_usda(:)
       tau_sorb(:) = tau_sorb_usda(:)
    CASE DEFAULT
       WRITE (numout,*) 'Unsupported soil type classification.Currently, only USDA classification is supported for phosphorus'
       WRITE (numout,*) 'p_nscm',p_nscm
       STOP 'phosphorus_init'
    ENDSELECT

    END SUBROUTINE phosphorus_init

! ===================================
!  ROUTINES USED FOR N AND P CYCLE
! ===================================

  SUBROUTINE calc_f_XYplant(npts,nvm, &
                            biomass,  &
                            xelement, &
                            yelement, &
                            f_XYplant &
                            )

    ! this subroutine computes the scaling functions f_NCplant, f_PNplant, and
    ! f_PCplant which are used to scale processes according to the stoichiometric status of
    ! vegetation; it is equation (10) p. 3 of SM of Zaehle & Friend, 2010 

    INTEGER, INTENT(in)                                        :: npts          !! Domain size (unitless)
    INTEGER, INTENT(in)                                        :: nvm           !! Domain size (unitless)

    REAL, DIMENSION(npts,nvm,nparts,nelements),INTENT(in)      :: biomass       !! Biomass pools [g m-2]
    INTEGER, INTENT(in)                                        :: xelement      !! index for the element X
    INTEGER, INTENT(in)                                        :: yelement      !! index for the element Y

    ! out
    REAL, DIMENSION(npts,nvm),INTENT(out)                      :: f_XYplant    !! scaling function using the ratio X to Y [0-1]

    ! locals
    REAL, DIMENSION(npts,nvm)   :: lab_x      !! labile mass of element X in plants
    REAL, DIMENSION(npts,nvm)   :: lab_y      !! labile mass of element Y in plants
    REAL, DIMENSION(npts)       :: XYplant    !! ratio of X to Y in plants

    ! parameters to be set
    REAL, DIMENSION(nvm)         :: fyx_root
    REAL, DIMENSION(nvm)         :: yx_leaf_prescribed
    REAL, DIMENSION(nvm)         :: xy_leaf_max
    REAL, DIMENSION(nvm)         :: xy_leaf_min

    INTEGER (i_std)                     :: m

    ! check with which element combination this routine was called
    ! and set the need parameters according to the combination
    IF ((xelement.EQ.2) .AND. (yelement.EQ.1)) THEN  ! N to C
        fyx_root(:)           = fcn_root(:)
        yx_leaf_prescribed(:) = cn_leaf_init(:) 

       ! minimal and maximal XY ratio of leaf (gX / gY)
        xy_leaf_max(1)        = zero
        xy_leaf_max(2:nvm)    = un/cn_leaf_min(2:nvm)       

        xy_leaf_min(1)        = -un/min_stomate
        xy_leaf_min(2:nvm)    =  un/cn_leaf_max(2:nvm)       

    ELSEIF ((xelement.EQ.3) .AND. (yelement.EQ.2)) THEN ! P to N
        fyx_root(:)           = fnp_root(:)
        yx_leaf_prescribed(:) = np_leaf_init(:) 

       ! minimal and maximal XY ratio of leaf (gX / gY)
        xy_leaf_max(1)        = zero
        xy_leaf_max(2:nvm)    = un/np_leaf_min(2:nvm)       

        xy_leaf_min(1)        = -un/min_stomate
        xy_leaf_min(2:nvm)    =  un/np_leaf_max(2:nvm)       


    ELSEIF ((xelement.EQ.2) .AND. (yelement.EQ.3)) THEN ! N to P
        fyx_root(:)           = fnp_root(:)
        yx_leaf_prescribed(:) = un/np_leaf_init(:) 
!       ! minimal and maximal XY ratio of leaf (gX / gY)
        ! we don't need to set them as we use a different scaling function 
!DSGnarrow03
!        IF (ANY(ABS(np_leaf_max(2:nvm)-60.).GT.min_stomate).OR. &
!            ANY(ABS(np_leaf_min(2:nvm)-5.).GT.min_stomate))THEN
!       IF (ANY(ABS(np_leaf_max(2:nvm)-35.).GT.min_stomate).OR. &
!           ANY(ABS(np_leaf_min(2:nvm)-12.).GT.min_stomate))THEN
       IF (ANY(ABS(np_leaf_max(2:nvm)-30.).GT.min_stomate).OR. &
           ANY(ABS(np_leaf_min(2:nvm)-12.).GT.min_stomate))THEN
!        
!DSGnarrow03
              CALL ipslerr_p ( 3, 'calc_f_XYplant: sigmoid scaling is not adjusted to this np_leaf range','model STOP', &
                                                 &         '','')
        ENDIF

!DSG   !DSG C:P ratio is not used for scaling in ORCHIDEE:
!    ELSEIF ((xelement.EQ.3) .AND. (yelement.EQ.1)) THEN ! P to C
!        fyx_root(:)           = fnp_root(:)*fcn_root(:)
!        yx_leaf_prescribed(:) = np_leaf_init(:)* cn_leaf_init 
!
!       ! minimal and maximal XY ratio of leaf (gX / gY)
!        xy_leaf_max(1)        = zero
!        xy_leaf_max(2:nvm)    = un/(np_leaf_min(2:nvm)*cn_leaf_min(2:nvm))    
!
!        xy_leaf_min(1)        = -un/min_stomate
!        xy_leaf_min(2:nvm)    =  un/(np_leaf_max(2:nvm)*cn_leaf_max(2:nvm))      

    ELSE
       STOP 'subroutine f_XYplant: subroutine does not support the element combination'
    ENDIF

    ! compute the XY ratio of labile compounds in plant (gX / gY)
    ! adapted from equation (9) p. 3 of SM of Zaehle & Friend, 2010
    lab_x(:,:) = biomass(:,:,ilabile,xelement) + &
         biomass(:,:,iroot,xelement) + biomass(:,:,ileaf,xelement)
    lab_y(:,:) = biomass(:,:,ilabile,yelement) + &
         biomass(:,:,iroot,yelement) + biomass(:,:,ileaf,yelement)

    ! set to zero
    f_XYplant(:,:) = zero
    XYplant = zero

    DO m=1,nvm

        IF(.NOT.((xelement.EQ.2) .AND. (yelement.EQ.3))) THEN
           ! LINEAR SCALING BETWEEN MAX AND MIN:
           ! set the X to Y status of plants 
           WHERE(lab_y(:,m).GT.min_stomate)
              ! in case we have biomass we set the actual ratio
              XYplant(:) = lab_x(:,m) / lab_y(:,m)
           ELSEWHERE
              ! for the case we have no biomass we use the prescribed ratio of roots
              XYplant(:) = fyx_root(m)/yx_leaf_prescribed(m) 
           ENDWHERE
           !WRITE(6,*) 'XYplant',XYplant(test_grid)

           ! Originally this was designed to:
           ! "Nitrogen demand responds to N/C ratio of the labile pool for each PFT
           ! zero uptake at a NC in the labile pool corresponding to maximal leaf NC 
           ! and 1. for a NC corresponding to a minimal leaf NC
           ! See equation (10) p. 3 of SM of Zaehle & Friend, 2010" 
           ! DSG: linear decrease of f with increase in xy
           ! we make it work for other combination of elements
           f_XYplant(:,m) = max( ( XYplant(:) - xy_leaf_max(m) ) / ( xy_leaf_min(m) - xy_leaf_max(m) ), zero )
           f_XYplant(:,m) = min( f_XYplant(:,m), un )
       ELSE
           ! SIGMOIDAL SCALING FOR WIDE N:P RANGE
           ! set the X to Y status of plants 
           WHERE(lab_y(:,m).GT.min_stomate)
              ! in case we have biomass we set the actual ratio
              XYplant(:) = lab_x(:,m) / lab_y(:,m)
           ELSEWHERE
              ! for the case we have no biomass we use the prescribed ratio of roots
              XYplant(:) = fyx_root(m)/yx_leaf_prescribed(m) 
           ENDWHERE
           ! the parameter here are fitted to the extrawide range of 5 to 60 N:P
           ! to get a strong increase in P aquiring process around N:P of 20
!DSGnarrow03
!           f_XYplant(:,m) = MAX(MIN(un/(un +exp((-( XYplant(:) )/2.)+10.)), un),zero)
!NHsink
!          f_XYplant(:,m) = MAX(MIN(un/(un +exp((-( XYplant(:) )/1.3)+15.)), un),zero)
           f_XYplant(:,m) = MAX(MIN(un/(un +exp((-( XYplant(:) )/1.6)+12.)), un),zero)
!NHsink
!DSGnarrow03
       ENDIF    

    ENDDO

  END SUBROUTINE calc_f_XYplant

  SUBROUTINE root_conductivity(npts,nvm,                      &  
                               soil_nutrient,max_eau_var,     & 
                               vmax_uptake, K_min,            &
                               low_K_min, saw,                &
                               dt,                            &
                               uptake_capacity_max)

    ! this subroutine calculates the maximum uptake of nutrient from the soil solution per g(C) root
    ! it is adapted from the uptake of nitrate and ammonia by Zaehle & Post 2010
  
    INTEGER, INTENT(in)                       :: npts          !! Domain size (unitless)
    INTEGER, INTENT(in)                       :: nvm           !! Domain size (unitless)

    REAL, DIMENSION(npts,nvm),INTENT(in)      :: soil_nutrient       !! nutrient which is available for uptake [g m-2]
    REAL, DIMENSION(npts),INTENT(in)          :: max_eau_var         !! maximum water holding capacity of soils [kg m-2]

    REAL, INTENT(in)                          :: vmax_uptake         !! [umol(N)/g(DW)/h
    REAL, INTENT(in)                          :: K_min               !! 
    REAL, INTENT(in)                          :: low_K_min           !! 
    REAL, INTENT(in)                          :: saw                 !! standard atomic weight of nutrient [g(nutrient) mol-1]
    REAL, INTENT(in)                          :: dt                  !! time step

    REAL, DIMENSION(npts,nvm),INTENT(out)     :: uptake_capacity_max !! maximum uptake rate per g(C) in root [g(nutrient) g (root C)-1 timestep-1]
    ! locals
    REAL                                      :: conv_fac_vmax       !! conversion factor
    REAL,DIMENSION(npts)                      :: conv_fac_concent    !! conversion factor
    INTEGER                                   :: m                   !! index    

     IF(printlev.GE.3)THEN
        WRITE(numout,*) 'Entering root_conductivity'
     ENDIF
   
    ! calculate converstation factors following Zaehle & Post (2010): 
    ! 1. Conversion from [umol(nutrient) g(DW)-1 h-1] to [g(nutrient) g(C)-1 timestep-1]
    conv_fac_vmax= 24. * dt * 2. * 1.e-6 * saw

    ! 2. Conversion factor from [umol(nutrient) per litter] to [g(nutrient) m-2] 
    ! saw  : molar mass for P (gP mol-1)
    ! 10^3 : conversion factor (dm3 to m3)
    ! 10-6 : conversion factor (ug to g)
    ! max_eau_var   (kg/m2) : water holding capacity as approximation of pore space
    ! (Smith et al. 2014, "Implications of incorporating N cycling and N
    ! limitations on primary production in an individual-based dynamic
    ! vegetation model")
    conv_fac_concent(:) = saw * 1.e3 * 1.e-6 * max_eau_var(:)/1.e3      

    DO m=1,nvm
       WHERE(max_eau_var.GT.zero)
         uptake_capacity_max(:,m) =  vmax_uptake * conv_fac_vmax * soil_nutrient(:,m) &
                                     *( low_K_min / conv_fac_concent(:) &
                                     + un / ( soil_nutrient(:,m) + K_min * conv_fac_concent(:) ) ) 
       ELSEWHERE
         uptake_capacity_max(:,m) =  zero
       ENDWHERE
    ENDDO
     IF(printlev.GE.3)THEN
        WRITE(numout,*) 'Leaving root_conductivity'
     ENDIF

  END SUBROUTINE root_conductivity


  SUBROUTINE root_surface_concentration(npts,nvm,dt,                   &  
                                        soil_orders,                   &
                                        clay,silt, bulk,               &
                                        root_mass,Pdissolved,          & 
                                        tmc_pft, swc_pft,f_Pdissolved)

    ! This routine computes the fraction of dissolved labile P which is in contact with the
    ! root surface. We use Fick's Law to simulate the diffusion of phosphate from the sourounding to the root zone, omitting advection. 
    ! Currently, we do not use information about the root profile in soils, but
    ! assume all fine roots (all root biomass is assume to be fine root biomass) are evenly distributed in the first 1m of soils.
   
 
    INTEGER, INTENT(in)                       :: npts          !! Domain size (unitless)
    INTEGER, INTENT(in)                       :: nvm           !! Domain size (unitless)
    REAL, INTENT(in)                          :: dt            !! time step
   
    INTEGER, DIMENSION(npts), INTENT(in)      :: soil_orders   !! dominant USDA soil order of grid 
    REAL,DIMENSION(npts),INTENT(in)           :: clay          !! Clay fraction of soil (0-1, unitless) 
    REAL,DIMENSION(npts),INTENT(in)           :: silt          !! Silt fraction of soil (0-1, unitless) 
    REAL,DIMENSION(npts),INTENT(in)           :: bulk          !! Bulk density (kg/m**3) 

    REAL, DIMENSION(npts,nvm), INTENT(in)     :: root_mass     !! root carbon [g(C) m-2]
    REAL, DIMENSION(npts,nvm), INTENT(in)     :: Pdissolved    !! dissolved labile P in soil solution [g m-2]
 
    REAL, DIMENSION(npts,nvm), INTENT(in)     :: tmc_pft       !! totalsoil water content per PFT [kg/m2]
    REAL, DIMENSION(npts,nvm), INTENT(in)     :: swc_pft       !! relative water content per PF [ ]
 
    REAL, DIMENSION(npts,nvm), INTENT(inout)  :: f_Pdissolved  !! fraction of dissolved labile P in soil solution in contact with root surface [ ] 
 
    ! locals 
    REAL, DIMENSION(npts)                    :: RLD            !! root length density [m m-3]
    REAL, DIMENSION(npts)                    :: r_diff         !! average distance P has to diffuse [m]
    REAL, DIMENSION(npts,nvm)                :: F_diff         !! diffusion flux to the root surface [ g m-2 timestep-1]

    REAL, DIMENSION(npts)                    :: swc_th         !! threshold soil water content at which solute diffusivity 
                                                                      !! (ratio of diffusion coefficients in soil and free water)
                                                                      !! approaches zero
                                                                      !! [m3 m-3]
    REAL, DIMENSION(npts,nvm)                :: f_tortuosity   !! tortuosity factor based on soil water content [ ] 
    REAL, DIMENSION(npts,nvm)                :: coeff_D        !! diffusion coefficient [cm day-1]
    REAL, DIMENSION(npts,nvm)                :: deltaPconc     !! concentration difference between root surface and boundary [mg mL-1]]  ; DSG: mL=cm3
    
    INTEGER                                      :: i,j
    
    ! conversion factors
    REAL                                        ::  conv_D_ref  !! for D_ref
    REAL                                        ::  conv_bulk   !! for bulk
    
    ! parameters ; should all be externalized
    ! root paremeters: 
    REAL,  DIMENSION(nvm)                       :: r_d                     !! specific root density [g(biomass) m-3(root)]; in FRED: (F00709) 'root tissue density (RTD)'
    REAL,  DIMENSION(nvm)                       :: r_r                     !! fine root radius [m]; in FRED: (F00679) 'root diameter'

     REAL,  DIMENSION(p_nscm)                      ::  f1
     REAL,  DIMENSION(p_nscm)                      ::  f2
     REAL,  DIMENSION(p_nscm)                      ::  swc_is

    ! diffusion constant for P in water
    REAL,  PARAMETER                            :: D_ref=0.759             !! Diffusion coefficient of phosphate in free water at 25degreeC [cm-2 day] (Mollier et al 2008))

! Parameters relating clay fraction, silt fraction, and soil bulk density to diffusion coefficient D_ref (Olesen et al. (2001); eq5 page 93)
    REAL, PARAMETER                                    ::  OlesenC1=0.81        !! empirical parameter relating  clay fraction
    REAL, PARAMETER                                    ::  OlesenC2=-0.90       !! empirical parameter relating  clay fraction **2
    REAL, PARAMETER                                    ::  OlesenS =-0.07       !! empirical parameter relating  silt fraction   
    REAL, PARAMETER                                    ::  OlesenB1=-0.60       !! empirical parameter relating  bulk density 
    REAL, PARAMETER                                    ::  OlesenB2=0.22        !! empirical parameter relating  bulk desnity**2
    REAL, PARAMETER                                    ::  OlesenA =0.42        !! empirical parameter offset 



     IF(printlev.GE.3)THEN
        WRITE(numout,*) 'Entering root_surface_concentration'
     ENDIF

!    convert from cm2 day-1 -> m2 timestep-1
     conv_D_ref=1.e4*dt
!    bulk density conversion from kg/m3 (ORCHIDEE) to Mg/m3 (Olesen et al. (2001) eq5)
     conv_bulk=1.e-3

    ! from Bonan et al (2014) Table3; at the moment globally uniform, as soon as
    ! FRED (http://roots.ornl.gov/) beccomes availabel these parameters can be
    ! refined and externalized:
    r_d(:) = 0.31e6 
    r_r(:) = 0.29e-3


!===================================================================================================    
! DSG: globally uniform parametrization
!DSG EXTERNALIZE  PARAMETERS
       ! from Mollier et al (2008) (original Barraclough and Tinker, 1981)
       ! Marko & Bruno plan to get soil order specific data 
       ! till then I use these parameters globally
       f1(:)    =1.58
       f2(:)    =-0.17
       swc_is(:)=0.12
!DSG EXTERNALIZE PARAMETERS

!    WRITE(6,*) 'swc(test_grid,test_pft)',swc(test_grid)
  
    ! 1. compute the diffusion coefficent [cm2 day-1] following Mollier et al. (2008) eq.10&11 (orig. Barraclough and Tinker (1981))
    ! 1.1 first compute the tortuosity factor based on soil water content (eq11) using a broken line function

    DO i=1,npts
      DO j=1,nvm
!DSG: not needed anymore as swc_pft is now really relative soil moisture content
!instead of whatever (tmc-tmcr)/(tmcs-tmcr)
!DSG        IF (swc(i,j).LT. swc_is(soil_orders(i))) THEN
!DSG           ! the original equation can get negative for very small values:
          f_tortuosity(i,j) = (swc_pft(i,j) *(f1(soil_orders(i))*swc_pft(i,j)+f2(soil_orders(i))))/(swc_is(soil_orders(i))) 
!DSG           !  I keep the diffusitivity factor fixed to the level of 0.12:
!DSG           f_tortuosity(i,j)= f1(soil_orders(i))*0.12+f2(soil_orders(i))
!DSG        ELSE
!DSG           f_tortuosity(i,j)= f1(soil_orders(i))*swc_pft(i,j)+f2(soil_orders(i))
!DSG        ENDIF
      ENDDO
    ENDDO
    ! 1.2 then compute the diffusion coefficient and convert it to ORCHIDEE
    ! units [m2 timestep-1]:
    coeff_D(:,:) = D_ref*conv_D_ref * swc_pft(:,:) * f_tortuosity(:,:)

!===================================================================================================    
! DSG: globally uniform parametrization
!DSG: accouting for soil heterogenity:
!DSG: doesnt work as often the swc in ORCHIDEE is well below the threshold and
!     thus diffusion is zero killing everything
!     Overall, I found that the Mollier equation covers the variation in
!     tortuosity factors not too bad -> so Mollier seems to  be at least for now
!     (2-layer hydro) thei better choice
! DSG: I keep this code; might be usefull later

!DSGdeactivated 
!DSGdeactivated     ! 1. compute the diffusion coefficent [m2 s-1] following Barraclough and Tinker (1981) 
!DSGdeactivated     !    taking into account heterogeinity in soil properties following Olesen et al. (2001)
!DSGdeactivated 
!DSGdeactivated     ! 1.1 first compute from silt, clay and soil bulk density  the threshold soil water content following Olesen et al. (2001) (eq5):
!DSGdeactivated     swc_th(:) = OlesenC1 * clay(:) + OlesenC2 * clay(:)**2 + OlesenS*silt(:) + OlesenB1*(conv_bulk*bulk) + OlesenB2*(conv_bulk*bulk(:))**2 + OlesenA
!DSGdeactivated     ! 1.2 compute the tortuosity factor following Olesen et al. (2001) (eq7a&b):
!DSGdeactivated     WHERE (swc(:).GT.swc_th(:)) 
!DSGdeactivated         f_tortuosity(:) = 1.1 * (swc(:)-swc_th(:))
!DSGdeactivated     ELSEWHERE
!DSGdeactivated         f_tortuosity(:) = zero
!DSGdeactivated     ENDWHERE
!DSGdeactivated     ! 1.3 then compute the diffusion coefficient following Barraclough and Tinker (1981) and convert to ORCHIDEE units
!DSGdeactivated     coeff_D(:) = D_ref*conv_D_ref * swc(:) * f_tortuosity(:)
!DSGdeactivated 
!DSGdeactivated !DSGdebug
!DSGdeactivated     WRITE(6,*) 'swc_th(:);',swc_th(:)
!DSGdeactivated     WRITE(6,*) 'swc(:);',swc(:)
!DSGdeactivated     WRITE(6,*) 'f_tortuosity(:) new method',f_tortuosity(:)
!DSGdeactivated !DSGdebug

!    WRITE(numout,*) 'root_mass(:,:)',root_mass(:,:)
!    WRITE(numout,*) 'r_r(:)',   r_r(:)
!    WRITE(numout,*) 'r_d(:)',r_d(:)
!    WRITE(numout,*) 'zmaxh', zmaxh
!    WRITE(numout,*) 'pi',pi
!    WRITE(numout,*) 'Pdissolved(:,:)',Pdissolved(:,:)
!    WRITE(numout,*) 'f_Pdissolved(:,:)',f_Pdissolved(:,:)
!    WRITE(numout,*) 'tmc_pft(:,:)',tmc_pft(:,:)
!    WRITE(numout,*) 'coeff_D(:,j)',coeff_D(:,:)

    DO j=1,nvm
      WHERE((root_mass(:,j).GT.min_stomate/1.e3).AND.(Pdissolved(:,j).GT.min_stomate/1.e3)) 
        ! 2.1 compute the root length density (eq.A30) [m m-3] ; 
        ! convert mass carbon to biomass with factor 2; 
        
        !DSG no soil layering yet RLD(:,j)    = (root_mass(:,j)*2./soil_layer_thickness)/(r_d(j) * pi * r_r(j)**2) 
        RLD(:)    = (2.*root_mass(:,j)/zmaxh)/(r_d(j) * pi * r_r(j)**2) 

        ! 2.2 compute half of the average distance between roots [m] (text between eq.A23/A24 in Bonan et al (2014))    
        r_diff(:) = (pi*RLD(:))**(-0.5)

        ! 2.3 restrict the diffusion path to 10 cm; this is needed for deciduous trees 
        ! (in a distance of 10 cm and more from the root surface, the concentration is not affected by plant uptake even under extreme conditions (Li et al.(1991))
        r_diff(:) = MIN(0.1,r_diff(:))

        ! 3. compute the difference in dissolved P concentration between root surface and diffusion boundary with Fick Law
        !    assuming soil water is homogenously distributed in the soil volume
        !    [ convert totwater [kg/m2] to [g/m2]]

        ! 3.1 concentration difference [g m-3]
        WHERE (tmc_pft(:,j).GT.zero)
            deltaPconc(:,j) = (Pdissolved(:,j)*(f_Pdissolved(:,j)-un))/(tmc_pft(:,j)*1E-3) 
        ELSEWHERE
            deltaPconc(:,j) = un
        ENDWHERE
 
        ! 3.2 the diffusion flux [g cm-2 timestep-1] following Ficks Law;
        F_diff(:,j) = - coeff_D(:,j)*(deltaPconc(:,j))/r_diff(:)

        ! 4. update the fraction of P concentration at root surface with the F_diff
        f_Pdissolved(:,j) = MAX(MIN((f_Pdissolved(:,j)*Pdissolved(:,j) + F_diff(:,j))/Pdissolved(:,j),un),zero)

      ELSEWHERE ! no roots no depletion zone
        f_Pdissolved(:,j) = un
      ENDWHERE
    ENDDO

!    WRITE(6,*) 'r_r(test_pft)',r_r(test_pft) 
!    WRITE(6,*) 'r_d(test_pft)',r_d(test_pft) 
!    WRITE(6,*) 'pi',pi
!    WRITE(6,*) 'root_mass(test_grid,test_pft)',root_mass(test_grid,test_pft)
!    WRITE(6,*) 'RLD(test_grid,test_pft)',RLD(test_grid,test_pft)
!    WRITE(6,*) 'half average distance for diffusion [m]:'
!    WRITE(6,*) 'r_diff(test_grid,test_pft)',r_diff(test_grid,test_pft)
!!
!!
!    WRITE(6,*) 'Diffusion flux [g m2 timestep-1]:'
!    WRITE(6,*) 'F_diff(test_grid,test_pft)',F_diff(test_grid,test_pft)
!!    WRITE(6,*) 'old fraction of dissolved P at root:'
!    WRITE(6,*) 'Pdissolved(test_grid,test_pft)',Pdissolved(test_grid,test_pft)
!    WRITE(6,*) 'f_Pdissolved(test_grid,test_pft)',f_Pdissolved(test_grid,test_pft)
    

   CALL histwrite_p (hist_id_stomate, 'f_PDISSOLVED', itime, &
        f_Pdissolved(:,:),       npts*nvm, horipft_index)
   CALL xios_orchidee_send_field('f_PDISSOLVED',f_Pdissolved(:,:))
 
     IF(printlev.GE.3)THEN
        WRITE(numout,*) 'Leaving root_surface_concentration'
     ENDIF

  END SUBROUTINE root_surface_concentration


! ===================================
! DEBUG TOOLS========================
! ===================================

  SUBROUTINE check_mass(npts,biomass,identifier)
! this routines checks 
! if any pool has negative values
! if biomass pools have a CNP ratio which is in the prescribed range


  INTEGER, INTENT(in)                                               :: npts              !! Domain size (unitless)
  ! biomass pools to be checked
  !REAL, DIMENSION(npts,nvm,ncirc,nparts,nelements), INTENT(in)      :: circ_class_biomass   !! Biomass components of the model tree  
  REAL, DIMENSION(npts,nvm,nparts,nelements),INTENT(in)             :: biomass   !! Biomass components of the model tree  

  
  CHARACTER(LEN=*), INTENT(in),OPTIONAL    ::  identifier                   !! string with to identify where this routine was called form
  

  ! locals
  INTEGER                                                                  :: dsg0,dsg1,dsg1a,dsg2,dsg3
  LOGICAL                                                                  :: crash

  CHARACTER(LEN=47)                              :: ident_routine

  REAL,DIMENSION(nparts)                 :: carbon 
  REAL,DIMENSION(nparts)                 :: nitrogen
  REAL,DIMENSION(nparts)                 :: phosphorus

  REAL,DIMENSION(nparts)                 :: cn_status 
  REAL,DIMENSION(nparts)                 :: np_status 

  REAL                                   :: cn_min
  REAL                                   :: np_min
  REAL                                   :: cn_max
  REAL                                   :: np_max

  LOGICAL    , DIMENSION(npts,nvm,nparts,2)     :: inrange
  LOGICAL                                       :: no_ratio
  LOGICAL                                       :: check_som

  REAL, PARAMETER                               :: CNP_accuracy=1.e-3


 ! IF (PRESENT(som)) THEN
 !     check_som=.TRUE.
 ! ELSE
     check_som=.FALSE.
 ! ENDIF 
 IF (PRESENT(identifier)) THEN
     ident_routine=identifier
 ELSE
      ident_routine=' no information'
 ENDIF


    ! 1. check for negative biomass pools:
    !WRITE (6,*) '=== BIOMASS CHECK ==='
    crash=.FALSE.
   IF (ANY(biomass(:,:,:,:).LT.zero))    THEN
       DO dsg0=1,npts
         DO dsg1=1,nvm
           DO dsg2=1,nparts
             DO dsg3=1,nelements
               IF ((biomass(dsg0,dsg1,dsg2,dsg3) .LT. zero).OR. & ! negative
                   (biomass(dsg0,dsg1,dsg2,dsg3) .GT. large_value).OR. & ! infinity 
                   (isnan(biomass(dsg0,dsg1,dsg2,dsg3))))   THEN ! NaN
                        WRITE (numout,*) 'FATAL negative biomass detected'
                        WRITE (numout,*) 'routine:', ident_routine
                        WRITE (numout,*) 'biomass = ', biomass(dsg0,dsg1,dsg2,dsg3)
                        WRITE (numout,*) 'box     = ',dsg0
                        WRITE (numout,*) 'PFT     = ',dsg1
                        WRITE (numout,*) 'npart   = ',dsg2
                        WRITE (numout,*) '1=leaf,2=sapabov,3=sapbelo,4=heartabov,5=heatbelow,6=root,7=fruit,8=carbres,9=ilabile'
                        WRITE (numout,*) 'element = ',dsg3
                ENDIF
             ENDDO
           ENDDO
         ENDDO
       ENDDO

       crash=.TRUE.
    END IF

    IF (crash) THEN
          WRITE (6,*) " FATAL: negative biomass"
           CALL ipslerr_p ( 3, 'A negative biomass pool was detected:','model STOP', &
                               &         '','')
    !ELSE
    !  WRITE (6,*) '     passed          '
    ENDIF

    ! 2. check for stoichiometry of biomass
    ! WRITE (6,*) '=== STOICHIOMETRY CHECK BIOMASS ==='
    inrange(:,:,:,:) =.FALSE.

    no_ratio  =.FALSE.
    np_status = -999
    cn_status = -999

    DO dsg0=1,npts
      DO dsg1=1,nvm 
        
        ! no need to check bare soil:
        IF (( dsg1 == 1 )) THEN
            inrange(dsg0,dsg1,:,:) = .TRUE.
        ELSE 
                 
          carbon    (:) = biomass(dsg0,dsg1,:,icarbon)
          nitrogen  (:) = biomass(dsg0,dsg1,:,initrogen)
          phosphorus(:) = biomass(dsg0,dsg1,:,iphosphorus)
                 
          ! if there is any biomass 
          IF ((ANY(carbon(:) .GT. min_stomate) )          &
             .OR. (ANY(nitrogen(:) .GT. min_stomate ))   &
             .OR. (ANY(phosphorus(:) .GT. min_stomate )) ) THEN

            WHERE (nitrogen.GT.zero) 
                cn_status(:) = carbon(:)  /nitrogen(:)
            ELSEWHERE
                cn_status(:) = cn_leaf_init(dsg1)
            ENDWHERE

            WHERE (phosphorus.GT.zero) 
               np_status(:) = nitrogen(:)/phosphorus(:)
            ELSEWHERE
               np_status(:) = np_leaf_init(dsg1)
            ENDWHERE

            ! loop over biomass compartments
            DO dsg2=1,nparts
              IF (dsg2 .EQ. ileaf) THEN
                    !  WRITE(6,*) 'leafi(1)?',dsg2
                    IF(impose_cn) THEN
                        cn_max = cn_leaf_init(dsg1) 
                        cn_min = cn_leaf_init(dsg1) 
                    ELSE
                        cn_max = cn_leaf_max(dsg1) 
                        cn_min = cn_leaf_min(dsg1) 
                    ENDIF
                    np_max = np_leaf_max(dsg1) 
                    np_min = np_leaf_min(dsg1) 
                    no_ratio  =.FALSE.
              ELSEIF (dsg2 .LT. iroot) THEN
                    !   Wood
                    IF (is_tree(dsg1)) THEN
                       ! since okrigidwood=TRUE we keep wood CNP ratio constant
                       cn_max = (cn_leaf_min(dsg1) / fcn_wood(dsg1))
                       cn_min = (cn_leaf_min(dsg1) / fcn_wood(dsg1))
                       np_max = (np_leaf_min(dsg1) / fnp_wood(dsg1))
                       np_min = (np_leaf_min(dsg1) / fnp_wood(dsg1))
                    ELSE
                       cn_max = cn_leaf_max(dsg1) / fcn_wood(dsg1)
                       cn_min = cn_leaf_min(dsg1) / fcn_wood(dsg1)
                       np_max = np_leaf_max(dsg1) / fnp_wood(dsg1)
                       np_min = np_leaf_min(dsg1) / fnp_wood(dsg1)
                    ENDIF
                    no_ratio  =.FALSE.
              ! We don't check fruits as their stoichiometry can deviate
              ! slightly from the prescribed range due to missing code in
              ! stomate_growth_fun transloc_N / transloc_P calcultion at the end
              ! of stomate_growth_fun:
              ELSEIF (dsg2 .EQ. iroot) THEN
              !
                    !   Roots:
                    cn_max = cn_leaf_max(dsg1) / fcn_root(dsg1)
                    cn_min = cn_leaf_min(dsg1) / fcn_root(dsg1)
                    np_max = np_leaf_max(dsg1) / fnp_root(dsg1)
                    np_min = np_leaf_min(dsg1) / fnp_root(dsg1)
                    no_ratio  =.FALSE.
              ELSE
                    ! Other pools have no fixed ratio
                    no_ratio=.TRUE.
              END IF

              ! start checking
              IF ((.NOT.no_ratio).AND.(biomass(dsg0,dsg1,dsg2,icarbon) .GT. min_stomate)) THEN

                ! check CN ratio
               !DSG debug IF((cn_status(dsg2).LT.cn_min)      &
               !DSG debug    .OR. (cn_status(dsg2).GT.cn_max) &
                IF((cn_status(dsg2)-cn_min.LT.-CNP_accuracy)      &
                   .OR. (cn_status(dsg2)-cn_max.GT.CNP_accuracy) &
                   ) THEN
                   inrange(dsg0,dsg1,dsg2,1) = .FALSE.
                    !  WRITE(6,*) 'out of range ',dsg2
                    WRITE(6,*) 'biomass(dsg0,dsg1,:,icarbon)'  ,biomass(dsg0,dsg1,:,icarbon)
                    WRITE(6,*) 'biomass(dsg0,dsg1,:,initrogen)',biomass(dsg0,dsg1,:,initrogen)
                ELSE
                   inrange(dsg0,dsg1,dsg2,1) = .TRUE.
                ENDIF

                ! check NP ratio
                !DSG debug IF(     (np_status(dsg2).LT.np_min) &
                !DSG debug    .OR. (np_status(dsg2).GT.np_max) &
                IF((np_status(dsg2)-np_min.LT.-CNP_accuracy)      &
                   .OR. (np_status(dsg2)-np_max.GT.CNP_accuracy) &
                   ) THEN
                      inrange(dsg0,dsg1,dsg2,2) = .FALSE.
                      WRITE(6,*) 'out of range ',dsg2
                ELSE
                   inrange(dsg0,dsg1,dsg2,2) = .TRUE.
                ENDIF

                ! write information to output
                IF (.NOT.(ALL(inrange(dsg0,dsg1,dsg2,:))))  THEN
                  WRITE (numout,*) 'routine:', ident_routine
                  WRITE (numout,*) '=== PFT ',dsg1, ' ==='
                  WRITE (numout,*) '=== box ',dsg0, ' ==='
                  WRITE (numout,*) 'plant part ', dsg2
                  WRITE (numout,*) '1=leaf,2=sapabov,3=sapbelo,4=heartabov,5=heatbelow,6=root,7=fruit,8=carbres,9=ilabile'
                  WRITE (numout,*) 'STATUS=FATAL: stoichemtry out of range' 

                  IF (.NOT.inrange(dsg0,dsg1,dsg2,1)) THEN
                    WRITE (numout,*) 'CN= '          ,cn_status
                    WRITE (numout,*) 'allowed range ',cn_min, ' to ',cn_max
                    WRITE (numout,*) 'accuarcy (+/-)', CNP_accuracy
                  ELSE IF (.NOT.inrange(dsg0,dsg1,dsg2,2)) THEN
                    WRITE (numout,*) 'NP= '          ,np_status
                    WRITE (numout,*) 'allowed range ',np_min, ' to ',np_max
                    WRITE (numout,*) 'accuarcy (+/-)', CNP_accuracy
                  END IF
                ENDIF

              ELSE
                inrange(dsg0,dsg1,dsg2,:) = .TRUE.
              ENDIF

            ENDDO
          ELSE
              ! no biomass no range
              inrange(dsg0,dsg1,:,:) = .TRUE.
          ENDIF
        ENDIF   ! bares/ soil crops
      ENDDO
    ENDDO

    IF (.NOT.ok_pcycle) inrange(:,:,:,2)=.TRUE.
    IF (.NOT.ok_ncycle) inrange(:,:,:,1)=.TRUE.

    ! crash model if out of range
    IF (ANY(.NOT.inrange(:,:,:,:))) THEN
           CALL ipslerr_p ( 3, 'stoichiometry is outside boundaries','model STOP', &
                               &         '','')
    !ELSE
    !     WRITE (6,*) '     passed          '
    ENDIF

  END SUBROUTINE check_mass

  SUBROUTINE cons_mass_r3d(mass_before,      &
                          mass_after,       &
                          mass_change,      &
                          identifier)

  ! biomass pools to be checked
  REAL, DIMENSION(:,:,:),INTENT(in)             :: mass_before   !! Biomass old
  REAL, DIMENSION(:,:,:),INTENT(in)             :: mass_after    !! Biomass new (npts,nvm,nelements)
  REAL, DIMENSION(:,:,:),INTENT(in)             :: mass_change   !! Biomass change

  CHARACTER(LEN=*), INTENT(in),OPTIONAL    ::  identifier                   !! string with to identify where this routine was called form


  ! locals
  INTEGER                                                           :: dsg1,dsg1a,dsg2,dsg3, ierr
  INTEGER                                                           :: nelem
  LOGICAL,ALLOCATABLE, DIMENSION(:,:,:)                             :: crash
  REAL, ALLOCATABLE, DIMENSION(:,:,:)                        :: mass_conserv           !! conservation 
  CHARACTER(LEN=47)                                                 :: ident_routine
  ! this subroutine checks for mass conservation
  INTEGER                                     :: npts          !! Domain size (unitless)
  INTEGER                                     :: nvm           !! Domain size (unitless)
  INTEGER                                     :: nelements     !! Domain size (unitless)
  
! ===================================================================================================
  npts = SIZE(mass_before, DIM=1)
  nvm = SIZE(mass_before, DIM=2)
  nelements = SIZE(mass_before, DIM=3)

  IF (PRESENT(identifier)) THEN
      ident_routine=identifier
  ELSE
       ident_routine=' no information'
  ENDIF

  ALLOCATE(crash(npts, nvm, nelements), stat=ierr)
  IF (ierr /= 0) CALL ipslerr_p(3, 'cons_mass_r3d', 'Memory allocation problem with ', 'crash variable', '')

  ALLOCATE(mass_conserv(npts, nvm, nelements), stat=ierr)
  IF (ierr /= 0) CALL ipslerr_p(3, 'cons_mass_r3d', 'Memory allocation problem with ', 'mass_conserv variable', '')

  crash(:,:,:)=.FALSE.
  ! bookkeeping:
  mass_conserv(:,:,:)          = mass_before(:,:,:) - mass_after(:,:,:) + mass_change(:,:,:)

  IF (ok_pcycle) THEN
      nelem=nelements
  ELSE
      nelem=nelements-1
  ENDIF

  IF (ok_ncycle) THEN
      nelem=nelem-1
  ENDIF

  ! check the bookkeeping
  DO dsg3=1,nelem
    DO dsg1=1,npts
      DO dsg2=1,nvm
         IF (ABS(mass_conserv(dsg1,dsg2,dsg3)).GT. (min_stomate)) THEN
                WRITE (numout,*) 'FATAL mass conservation failed (positive value = leak)'
                WRITE (numout,*) 'routine:', ident_routine
                WRITE (numout,*) ' limit: ',min_stomate
                WRITE (numout,*) ' limit: ',min_stomate
                WRITE (numout,*) ' for element :', dsg3, ' of nelements: ',nelements
                WRITE (numout,*) ' gridpoint: ',dsg1 , ' of ngrids: ',npts
                WRITE (numout,*) ' PFT: ',dsg2 , ' of npfts: ',nvm
                WRITE (numout,*) ' mismatch =', mass_conserv(dsg1,dsg2,dsg3)
                WRITE (numout,*) ' mass(before) =', mass_before(dsg1,dsg2,dsg3)
                WRITE (numout,*) ' mass(after) =', mass_after(dsg1,dsg2,dsg3)
                WRITE (numout,*) ' test_grid =', test_grid
                WRITE (numout,*) ' test_pft =', test_pft
                crash(dsg1,dsg2,dsg3)=.TRUE.
         ENDIF  
      ENDDO
    ENDDO
  ENDDO

  IF ((mass_conservation).AND.(ANY(crash(:,:,:)))) THEN
     CALL ipslerr_p ( 3, 'mass_conservation: FAILED','model STOP', &
          &         '','')
  ENDIF

  DEALLOCATE(mass_conserv)
  DEALLOCATE(crash)

  END SUBROUTINE cons_mass_r3d

  SUBROUTINE cons_mass_r2d(mass_before,      &
                          mass_after,       &
                          mass_change,      &
                          identifier)


  ! biomass pools to be checked
  REAL, DIMENSION(:,:),INTENT(in)             :: mass_before   !! Biomass old
  REAL, DIMENSION(:,:),INTENT(in)             :: mass_after    !! Biomass new (npts,nvm)
  REAL, DIMENSION(:,:),INTENT(in)             :: mass_change   !! Biomass change

  CHARACTER(LEN=*), INTENT(in),OPTIONAL    ::  identifier                   !! string with to identify where this routine was called form


  ! locals
  INTEGER                                                           :: dsg1,dsg1a,dsg2, ierr
  LOGICAL,ALLOCATABLE,DIMENSION(:,:)                             :: crash
  REAL,ALLOCATABLE, DIMENSION(:,:)                        :: mass_conserv           !! conservation 
  CHARACTER(LEN=47)                                                 :: ident_routine
  ! this subroutine checks for mass conservation
  INTEGER                                     :: npts          !! Domain size (unitless)
  INTEGER                                     :: nvm           !! Domain size (unitless)
! ===================================================================================================
  npts = SIZE(mass_before, DIM=1)
  nvm = SIZE(mass_before, DIM=2)

  ALLOCATE(crash(npts, nvm), stat=ierr)
  IF (ierr /= 0) CALL ipslerr_p(3, 'cons_mass_r2d', 'Memory allocation problem with ', 'crash variable', '')

  ALLOCATE(mass_conserv(npts, nvm), stat=ierr)
  IF (ierr /= 0) CALL ipslerr_p(3, 'cons_mass_r2d', 'Memory allocation problem with ', 'mass_conserv variable', '')

  IF (PRESENT(identifier)) THEN
      ident_routine=identifier
  ELSE
       ident_routine=' no information'
  ENDIF


  crash(:,:)=.FALSE.
  ! bookkeeping:
  mass_conserv(:,:)          = mass_before(:,:) - mass_after(:,:) + mass_change(:,:)

  ! check the bookkeeping
  DO dsg1=1,npts
    DO dsg2=1,nvm
       IF (ABS(mass_conserv(dsg1,dsg2)).GT. (min_stomate)) THEN
              WRITE (numout,*) 'FATAL mass conservation failed (positive value = leak)'
              WRITE (numout,*) 'routine:', ident_routine
              WRITE (numout,*) ' limit: ',min_stomate
              WRITE (numout,*) ' limit: ',min_stomate
              WRITE (numout,*) ' gridpoint: ',dsg1 , ' of ngrids: ',npts
              WRITE (numout,*) ' PFT: ',dsg2 , ' of npfts: ',nvm
              WRITE (numout,*) ' mismatch =', mass_conserv(dsg1,dsg2)
              WRITE (numout,*) ' mass(before) =', mass_before(dsg1,dsg2)
              WRITE (numout,*) ' mass(after) =', mass_after(dsg1,dsg2)
              WRITE (numout,*) ' test_grid =', test_grid
              WRITE (numout,*) ' test_pft =', test_pft
              crash(dsg1,dsg2)=.TRUE.
       ENDIF  
    ENDDO
  ENDDO

  IF ((mass_conservation).AND.(ANY(crash(:,:)))) THEN
     CALL ipslerr_p ( 3, 'mass_conservation: FAILED','model STOP', &
          &         '','')
  ENDIF

  DEALLOCATE(mass_conserv)
  DEALLOCATE(crash)

  END SUBROUTINE cons_mass_r2d


 !=========================
  ! the following is copied from stomate_soilcarbon.f90
  ! the module of nitrogen phosphorus and carbon should be restructured to avoid
  ! such things

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

END MODULE stomate_phosphorus
