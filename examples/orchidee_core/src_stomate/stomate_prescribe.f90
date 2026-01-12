! =================================================================================================================================
! MODULE       : stomate_prescribe
!
! CONTACT      : orchidee-help _at_ listes.ipsl.fr
!
! LICENCE      : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF         Initialize and update density, crown area.
!!
!!\n DESCRIPTION: None
!!
!! RECENT CHANGE(S): None
!!
!! REFERENCE(S)	:
!!
!! SVN          :
!! $HeadURL: svn://forge.ipsl.jussieu.fr/orchidee/branches/ORCHIDEE-CN-P/ORCHIDEE/src_stomate/stomate_prescribe.f90 $
!! $Date: 2018-10-09 10:10:10 +0200 (二, 2018-10-09) $
!! $Revision: 5473 $
!! \n
!_ ================================================================================================================================

MODULE stomate_prescribe

  ! modules used:

  USE ioipsl_para
  USE stomate_data
  USE pft_parameters
  USE constantes
  USE function_library,    ONLY: calculate_c0_alloc

  IMPLICIT NONE

  ! private & public routines

  PRIVATE
  PUBLIC prescribe,prescribe_clear

    ! first call
    LOGICAL, SAVE                                              :: firstcall_prescribe = .TRUE.
!$OMP THREADPRIVATE(firstcall_prescribe)

CONTAINS

! =================================================================================================================================
!! SUBROUTINE   : prescribe_clear
!!
!>\BRIEF        : Set the firstcall_prescribe flag back to .TRUE. to prepare for the next simulation.
!_=================================================================================================================================

  SUBROUTINE prescribe_clear
    firstcall_prescribe=.TRUE.
  END SUBROUTINE prescribe_clear


!! ================================================================================================================================
!! SUBROUTINE   : prescribe
!!
!>\BRIEF         Works only with static vegetation and agricultural PFT. Initialize biomass,
!!               density, presence in the first call and update them in the following.
!!
!! DESCRIPTION (functional, design, flags): \n
!! This module works only with static vegetation and agricultural PFT.
!! In the first call, initialize density of individuals, biomass,
!! and leaf age distribution to some reasonable value. In the following calls,
!! these variables are updated.
!!
!! To fulfill these purposes, pipe model are used:
!! \latexonly 
!!     \input{prescribe1.tex}
!!     \input{prescribe2.tex}
!! \endlatexonly
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLES(S): ::ind, ::cn_ind, ::leaf_frac
!!
!! REFERENCES   :
!! - Krinner, G., N. Viovy, et al. (2005). "A dynamic global vegetation model 
!!   for studies of the coupled atmosphere-biosphere system." Global 
!!   Biogeochemical Cycles 19: GB1015, doi:1010.1029/2003GB002199.
!! - Sitch, S., B. Smith, et al. (2003), Evaluation of ecosystem dynamics,
!!   plant geography and terrestrial carbon cycling in the LPJ dynamic 
!!   global vegetation model, Global Change Biology, 9, 161-185.
!! - McDowell, N., Barnard, H., Bond, B.J., Hinckley, T., Hubbard, R.M., Ishii, 
!!   H., Köstner, B., Magnani, F. Marshall, J.D., Meinzer, F.C., Phillips, N., 
!!   Ryan, M.G., Whitehead D. 2002. The relationship between tree height and leaf 
!!   area: sapwood area ratio. Oecologia, 132:12–20.
!! - Novick, K., Oren, R., Stoy, P., Juang, F.-Y., Siqueira, M., Katul, G. 2009. 
!!   The relationship between reference canopy conductance and simplified hydraulic 
!!   architecture. Advances in water resources 32, 809-819. 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

 SUBROUTINE prescribe (npts, veget_cov_max, dt, PFTpresent, &
               everywhere, when_growthinit, biomass, leaf_frac, &
               ind, co2_to_bm,n_to_bm,p_to_bm, &
               KF, &
               senescence, age, npp_longterm, lm_lastyearmax, k_latosa_adapt)

!! 0. Parameters and variables declaration

   !! 0.1 Input variables

    INTEGER, INTENT(in)                        :: npts              !! Domain size (unitless)
    REAL, INTENT(in)                           :: dt                !! time step (dt_days)   
    REAL, DIMENSION(:,:), INTENT(in)           :: veget_cov_max     !! "maximal" coverage fraction of a PFT 
                                                                           !! (LAI -> infinity) on ground. May sum to
                                                                           !! less than unity if the pixel has
                                                                           !! nobio area. (unitless; 0-1)

    !! 0.2 Output variables
 
   
    !! 0.3 Modified variables

    LOGICAL, DIMENSION(:,:), INTENT(inout)            :: PFTpresent         !! PFT present (0 or 1)
    LOGICAL, DIMENSION(:,:), INTENT(inout)            :: senescence         !! Flag for setting senescence stage (only 
                                                                            !! for deciduous trees)
    REAL, DIMENSION(:,:), INTENT(inout)        :: everywhere         !! is the PFT everywhere in the grid box or 
                                                                            !! very localized (after its introduction) (?)
    REAL, DIMENSION(:,:), INTENT(inout)        :: when_growthinit    !! how many days ago was the beginning of 
                                                                            !! the growing season (days)
    REAL, DIMENSION(:,:), INTENT(inout)        :: ind                !! Density of individuals at the stand level 
                                                                            !! @tex $(m^{-2})$ @endtex
    REAL, DIMENSION(:,:), INTENT(inout)        :: npp_longterm       !! "long term" net primary productivity
                                                                            !! @tex ($gC m^{-2} year^{-1}$) @endtex
    REAL, DIMENSION(:,:), INTENT(inout)        :: age                !! mean age (years)
    REAL, DIMENSION(:,:), INTENT(inout)        :: lm_lastyearmax     !! last year's maximum leaf mass for each PFT 
                                                                            !! @tex ($gC m^{-2}$) @endtex
    REAL, DIMENSION(:,:,:,:), INTENT(inout)    :: biomass            !! Stand level biomass 
                                                                            !! tex $(gC.m^{-2})$ @endtex

    REAL, DIMENSION(:,:), INTENT(inout)        :: co2_to_bm          !! CO2 taken from the atmosphere to get C
                                                                            !! to create the seedlings 
                                                                            !!  @tex (gC.m^{-2}dt^{-1})$ @endtex  
    REAL, DIMENSION(npts,nvm), INTENT(inout)   :: n_to_bm            !! N taken to create the seedlings
    REAL, DIMENSION(npts,nvm), INTENT(inout)   :: p_to_bm            !! P taken to create the seedlings
    REAL, DIMENSION(:,:,:,:), INTENT(inout)    :: leaf_frac          !! fraction of leaves in leaf age 
                                                                            !! class (unitless;0-1)
    REAL, DIMENSION(:,:), INTENT(inout)        :: KF                 !! Scaling factor to convert sapwood mass
                                                                            !! into leaf mass (m) - this variable is 
                                                                            !! passed to other routines
    REAL, DIMENSION(:,:), INTENT(inout)        :: k_latosa_adapt     !! Leaf to sapwood area adapted for long 
                                                                            !! term water stress (m)
  

    !! 0.4 Local variables

    REAL, DIMENSION(nvm)                        :: c0_alloc           !! Root to sapwood tradeoff parameter
    INTEGER                                     :: ipts, ivm, ipar, k !! index (unitless)
    INTEGER                                     :: iele, imbc   !! index (unitless)
    INTEGER                                     :: deb,fin, imaxt     !! index (unitless)
    REAL, DIMENSION(npts,nvm,nparts,nelements)  :: sync_biomass       !! Temporary stand level biomass 
                                                                             !! @tex $(gC.m^{-2})$ @endtex
    REAL, DIMENSION(npts,nvm)                   :: sync_ind           !! Temporary density of individuals at the 
                                                                             !! stand level @tex $(m^{-2})$ @endtex
    REAL, DIMENSION(npts,nvm)                   :: k_latosa           !! Height dependent base value to calculate 
                                                                             !! KF (-) 
!    REAL, DIMENSION(nvm,nparts,nelements) :: bm_sapl            !! Sapling biomass for the functional 
                                                                             !! allocation with a dimension for the 
                                                                             !! circumference classes
                                                                             !! @tex $(gC.ind^{-1})$ @endtex
    REAL, DIMENSION(npts,nvm)                   :: LF                 !! Scaling factor to convert sapwood mass
                                                                             !! into root mass (unitless)
    REAL, DIMENSION(npts,nvm)                   :: lstress_fac        !! Light stress factor, based on total
                                                                             !! transmitted light (unitless, 0-1)
    REAL, DIMENSION(nparts)                     :: bm_init            !! Biomass needed to initiate the next 
                                                                             !! planting @tex $(gC m^{-2})$ @endtex 
    REAL                                        :: nb_trees_i         !! Number of trees in each twentith 
                                                                             !! circumference quantile of the 
                                                                             !! distribution (ind) 
    REAL                                        :: excedent           !! Number of trees after truncation to be 
                                                                             !! reallocated to smaller quantiles of the 
                                                                             !! distribution (ind) 
    REAL                                        :: ave_tree_height    !! The height of the ideal tree in each
                                                                             !! circumference class of the 
                                                                             !! distribution (ind)...not saved since
                                                                             !! it should only be used for prescribing
    REAL, DIMENSION(npts,nvm,nmbcomp,nelements) :: check_intern       !! Contains the components of the internal
                                                                             !! mass balance chech for this routine
                                                                             !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    REAL, DIMENSION(npts,nvm,nelements)         :: closure_intern     !! Check closure of internal mass balance
                                                                             !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    REAL, DIMENSION(npts,nvm,nelements)         :: pool_start         !! Start and end pool of this routine 
                                                                             !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    REAL, DIMENSION(npts,nvm,nelements)         :: pool_end           !! Start and end pool of this routine 
                                                                             !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    REAL                                        :: delta_KF           !! Difference between new and old estimate 
                                                                             !! of KF while iterating
    REAL                                        :: min_height_init           
    REAL                                        :: max_height_init           
    REAL                                        :: sapwood_density
    REAL                                        :: dia_init           
    REAL                                        :: wood_init                    
    INTEGER                                     :: iloop

    !++++ temp ++++
!!$    REAL, DIMENSION(npts,nvm,ncirc)             :: height,dia,cn       !! Tree height calculated from allometric 
!!$                                                                              !! relationships (m) 
    REAL                                        :: cn_leaf,cn_root,cn_wood  
                                                                                 !! CN ratio of leaves, root and wood pool
                                                                                 !! (gC/gN)        
    REAL                                        :: np_leaf,np_root,np_wood  
                                                                                 !! NP ratio of leaves, root and wood pool
                                                                                 !! (gN/gP)        
    REAL                                        :: fcn_wood_act           !! Actual wood CN to leaf CN ratio; 
    REAL                                        :: fnp_wood_act           !! Actual wood NP to leaf NP ratio; 

!_ ================================================================================================================================

 !! 1. Initialize biomass at first call

    !! 1.1.1 Initialize check for mass balance closure
    !  The mass balance is calculated at the end of this routine
    !  in section 4
    !  Initial biomass pool
    pool_start(:,:,:) = zero
    DO ipar = 1,nparts
       DO iele = 1,nelements
          pool_start(:,:,iele) = pool_start(:,:,iele) + &
               (biomass(:,:,ipar,iele) * veget_cov_max(:,:))
       ENDDO
    ENDDO

   ! co2_to_bm is has as intent inout, the variable accumulates
   ! carbon over the course of a day. Use the difference between
   ! start and the end of this routine
   check_intern(:,:,iatm2land,icarbon) = - un * co2_to_bm(:,:) * veget_cov_max(:,:) * dt

   ! Prescribe was taken out of the first call loop. Because of the firstcall here, 
   ! nothing was done when the vegetation was killed in lpj_gap and the vegetation
   ! thus stayed dead forever.  We want it to regrow.  So instead, let's loop over 
   ! every point and PFT.  If there is supposed to be vegetation here (according to 
   ! veget_cov_max) but there isn't (according to ind), then we grow it.
   DO ivm = 2,nvm ! Loop over PFTs


      ! Preparations: Calculate current tissue stoichiometries
      cn_leaf=cn_leaf_init(ivm)
      cn_root=cn_leaf/fcn_root(ivm)
      np_leaf=np_leaf_init(ivm)
      np_root=np_leaf/fnp_root(ivm)
      IF ( is_tree(ivm) ) THEN
          ! for rigid wood CNP ratios we need to update fcnp_wood_act:
          fcn_wood_act = fcn_wood(ivm) * cn_leaf/cn_leaf_min(ivm)
          fnp_wood_act = fnp_wood(ivm) * np_leaf/np_leaf_min(ivm)
      ELSE
          ! flexible wood CNP ratios we dont need to
          fcn_wood_act = fcn_wood(ivm)
          fnp_wood_act = fnp_wood(ivm)
      ENDIF
      ! using updated fcnp_wood_act we can get also the wood CNP:
      cn_wood=cn_leaf/fcn_wood_act
      np_wood=np_leaf/fnp_wood_act

      DO ipts = 1, npts ! Loop over pixels


         ! We are supposed to have vegetation, but we don't have any: prescribe some.
         IF((veget_cov_max(ipts,ivm) .GT. min_stomate) .AND. (ind(ipts,ivm) .LE. min_stomate))THEN

            ! Initilaize tree level variables
            ind(ipts,ivm) = zero
            biomass(ipts,ivm,:,:) = zero

            ! Assume that there is no memory for k_latosa between
            ! different generations (this is probably true when trees
            ! are planted but seems less likely for natural regeneration
            ! Allowing for memory has caused problems with the LAI from
            ! the second generation onwards 
            k_latosa_adapt(ipts,ivm) = k_latosa_min(ivm)

            ! Write output
            IF (printlev >= 4 ) THEN
               WRITE(numout,*) 'Prescribe (stomate_prescribe.f90):'
               WRITE(numout,*) '   > Imposing initial biomass for prescribed trees, '// &
                  'initial reserve mass for prescribed grasses.'
               WRITE(numout,*) '   > Declaring prescribed PFTs present.'
               WRITE(numout,*) '   > Grid point: ',ipts,' PFT Type: ',ivm
               WRITE(numout,*) 'k_latosa_adapt: ',k_latosa_adapt(ipts,ivm)
            ENDIF
            !! 2. Calculate the vegetation characteristics of a newly established vegetation

            !! 2.1 Stress factors

            ! Note that light and water stress have an opposite effect on KF. Waterstress 
            ! will decrease KF because more C should be allocated to the roots. Light 
            ! stress will increase KF because more C should be allocated to the leaves.

            ! We will need the c0_alloc, it only depends on PFT and the effective 
            ! longiveties
            c0_alloc(ivm) = calculate_c0_alloc(ipts, ivm, tau_root(ivm), &
                 tau_sap(ivm))

            ! Lightstress varies from 0 to 1 and is calculated from the canopy structure (veget)
            ! Given that there is no vegetation at this point, the lightstress cannot be calculated 
            ! and is set to 1. However as soon as we grow a canopy it will experience light stress
            ! and so KF should be adjusted. If this is not done, we can't prescribe tall vegetation
            ! which is a useful feature to speed up optimisation and testing. We will have iterate
            ! over light stress and KF to prescribe vegetation that is in balance with its 
            ! light environment. 
            lstress_fac(ipts,ivm) = un

            ! Initial vegetation has to be prescribed only when the vegetation is static
            IF ( ( .NOT. ok_dgvm ) .AND. &
                 ( veget_cov_max(ipts,ivm) .GT. min_stomate ) )  THEN

               !! 2.2 Initialize woody PFT's 
               !      Use veget_cov_max to check whether the PFT is present. If the PFT is present but it has no
               !      biomass, prescribe its biomass (gC m-{2}).
               IF ( is_tree(ivm) .AND. &
                    ( veget_cov_max(ipts,ivm) .GT. min_stomate ) .AND. &
                    ( SUM( biomass(ipts,ivm,:,icarbon) ) .LE. min_stomate ) ) THEN

                  ! PFT is present so it needs to be initialized 
                  IF (veget_cov_max(ipts,ivm) .GT. min_stomate) THEN

                     ! The forest stand starts with the number of individuals as prescribed in constantes_mtc.f90
                     ! This number is defined per hectare whereas these calculations are per m2
                     ind(ipts,ivm) = nmaxtrees(ivm) * ha_to_m2   

                     !! 2.3 Initialize height distribution
                     !  Initializing height distribution...I've taken this from the circumference initilization of Valentin
                     !  we want ncirc bins between height_init_min and height_init_max
                     min_height_init = height_init_min(ivm)
                     max_height_init = height_init_max(ivm)

                     !! 2.4 Calculate height distribution

                     !! 2.4.1 Height distribution of a high stand
                     !  Deterministic initial distribution following a truncated exponential law 
                     !  if they are coppices, we handle this differently


                     ! this array is going to be used to create the sapling for this class so that we can create 
                     ! the total amount of biomass
                     ! in this class...after this we will always calculate the height/circ/diam/whatever 
                     ! from the biomass in the class
                     ave_tree_height=(0.5_r_std)*&
                          (max_height_init+min_height_init)

                     ! I am changing this so that our lowest class is min_height_init and our biggest class
                     ! is max_height_init.  With this change, we can use the same code the redistribute
                     ! the trees if the biomass in one of our height classes is equal to zero later on due to mortality or
                     ! thinning.
                     !+++ This caused some strange behavior, so change it back for now.
                     !                        ave_tree_height(icir)=min_height_init+REAL(icir-1,r_std)/REAL(ncirc-1,r_std)*&
                     !                             (max_height_init-min_height_init)


                     IF (printlev>=4) THEN
                        WRITE(numout,*)"min initial height ",min_height_init,'max initial height ',max_height_init
                     ENDIF

!                     ENDIF ! check FM type

                     !! 2.4.2 Circumference distribution of a coppice
                     ! I'm not see why this should be any different than a standard prescribe.
                     ! Ideally, the plantings should be stems 20-25 cm long stuck into the soil,
                     ! but the allocation scheme will have the same issues with this that it does
                     ! with planting a small sapling, i.e. it won't like it.  The stems do not
                     ! follow standard allocation rules.  Then again, neither do coppices.

!!$                     IF (forest_managed(ipts,ivm) == ifm_src) THEN
!!$                        CALL ipslerr(3,'stomate_prescribe.f90','Problem with SRC','','')
!!$                     ENDIF

                     !! 2.5 Allocation factors
                     !  Sapwood to root ratio
                     !  Following Magnani et al. 2000 "In order to decreases hydraulic resistance,
                     !  the investment of carbon in fine roots or sapwood yields to the plant very
                     !  different returns, both because of different hydraulic conductivities and 
                     !  because of the strong impact of plant height on shoot resistance. On the 
                     !  other hand, fine roots and sapwood have markedly different longevities and 
                     !  the cost of production, discounted for turnover, will differ accordingly.
                     !  Optimal growth under hydraulic constraints requires that the ratio of marginal 
                     !  hydraulic returns to marginal annual cost for carbon investment in either roots 
                     !  or sapwood be the same (Bloom, Chapin & Mooney 1985; Case & Fair 1989). This
                     !  is formalized in equation (13) and further derived to obtain equation (17) 
                     !  in Magnani et al 2000. The latter is implemented here.
                     !  Pipe_density is given in gC/m-3, convert to kg/m-3. And apply equation (17) 
                     !  in Magnani et al 2000. Note that c0_alloc was calculated at the start of this 
                     !  routine. The calculation itself is done in function_library
                     
                     ! Calculate leaf area to sapwood area
                     ! To be consistent with the hydraulic limitations and pipe theory,
                     ! k_latosa is calculated from equation (18) in Magnani et al.
                     ! To do so, total hydraulic resistance and tree height need to known. This
                     ! poses a problem as the resistance depends on the leaf area and the leaf 
                     ! area on the resistance. There is no independent equation and equations 12
                     ! and 18 depend on each other and substitution would be circular. Hence 
                     ! prescribed k_latosa values were obtained from observational records 
                     ! and are given in mtc_parameters.f90. 

                     ! The relationship between height and k_latosa as reported in McDowell 
                     ! et al 2002 and Novick et al 2009 is implemented to adjust k_latosa for 
                     ! the height of the stand. This did NOT result in a realistic model behavior 
                     !!$ k_latosa(ipts,ivm) = wstress_fac(ipts,ivm) * &
                     !!$    (k_latosa_max(ivm) - (latosa_height(ivm) * &
                     !!$    (SUM( nb_trees_i(:) * ave_tree_height(:) ) / SUM( nb_trees_i(:) ))))

                     ! Alternatively, k_latosa is also reported to be a function of diameter 
                     ! (i.e. stand thinning, Simonin et al 2006, Tree Physiology, 26:493-503).
                     ! Here the relationship with thinning was interpreted as a realtionship with
                     ! light stress. Note that light stress cannot be calculated at this time in
                     ! the model because there is no canopy (that's why we are in prescribe!) and
                     ! so there is no lightstress. lstress was therefore set to one (see above).
                     ! We prefered this redundancy in the code because it makes it clear that
                     ! k_latosa is calculated in the same way in prescribe.f90 and growth_fun_all.f90
                     ! +++CHECK+++
                     ! How do we want to deal with waterstress
!!$                     k_latosa(ipts,ivm) = k_latosa_min(ivm) + &
!!$                          (wstress_fac(ipts,ivm) * lstress_fac(ipts,ivm) * &
!!$                          (k_latosa_max(ivm)-k_latosa_min(ivm)))
!!$                     k_latosa(ipts,ivm) = wstress_fac(ipts,ivm) * (k_latosa_min(ivm) + &
!!$                          (lstress_fac(ipts,ivm) * &
!!$                          (k_latosa_max(ivm)-k_latosa_min(ivm))))
                     k_latosa(ipts,ivm) = (k_latosa_adapt(ipts,ivm) + &
                            (lstress_fac(ipts,ivm) * &
                          (k_latosa_max(ivm)-k_latosa_min(ivm))))
                     ! ++++++++++++
                     

                     ! Also k_latosa has been reported to be a function of CO2 concentration 
                     ! (Atwell et al. 2003, Tree Physiology, 23:13-21 and Pakati et al. 2000, 
                     ! Global Change Biology, 6:889-897). This effect is not accounted for in 
                     ! the current code

                     ! Calculate conversion coefficient for sapwood area to leaf area 
                     ! (1) The scaling parameter between leaf and sapwood mass is derived from
                     ! LA_ind = k_latosa * SA_ind, where LA_ind = leaf area of an individual, SA_ind is the 
                     ! sapwood area of an individual and k_latosa a pipe-model parameter
                     ! (2) LA_ind = Cl * sla
                     ! (3) Cs = SA_ind * height * wooddensity * tree_ff
                     ! Substitute (2) and (3) in (1)
                     ! Cl = Cs * k1 / (wooddensity * sla * tree_ff * height)
                     ! Cl = Cs*KF/height, where KF is in (m)
                     ! KF is passed to the allocation routine and it is saved in the restart file.
                     KF(ipts,ivm) = k_latosa(ipts,ivm) / ( sla(ivm) * pipe_density(ivm) * tree_ff(ivm))
                     
                     ! Initialize delta_KF to get the DO WHILE started
                     delta_KF = un

                     iloop=0
                     DO WHILE (delta_KF .GT. max_delta_KF)

                        !  If there is a WHILE loop, there always needs to be a check on the number
                        ! of loops to make sure we don't get stuck in an infinite loop.  This
                        ! number is completely arbitrary.
                        iloop=iloop+1
                        IF(iloop > 1000)THEN
                           WRITE(numout,*) 'Taking too long to converge the delta_KF loop in prescribe!'
                           WRITE(numout,*) 'iloop,delta_KF,max_delta_KF,ivm,ipts: ',&
                                iloop,delta_KF,max_delta_KF,ivm,ipts
                           CALL ipslerr_p (3,'stomate_prescribe',&
                                'Taking too long to converge the delta_KF','','')
                                
                        ENDIF

                        !! 2.6 Create saplings for each height class
                        !  Now we create the saplings for each class, based on the height
                           
                        ! The assumption we make is that we plant trees of 2 to 3 years old rather than 
                        ! growing trees from seeds. The allometric relationship between height and 
                        ! diameter is derived from mature tree and likely unrealistic for saplings. 
                        ! The height of the saplings is prescribed and determines the reserves which are 
                        ! especially important for deciduous species which need to survive on their 
                        ! reserves for the first year (new phenology scheme requires annual mean values
                        ! to get started) 
                        dia_init = ( ave_tree_height / pipe_tune2(ivm) ) ** ( 1. / pipe_tune3(ivm) )
                        wood_init = ( ave_tree_height * pi / 4. * (dia_init) ** 2. ) * &
                             pipe_density(ivm) * tree_ff(ivm)
                        
                        ! The woody biomass is contained in four components. Thus, wood_init = isapabove + 
                        ! isapbelow + iheartabove + iheartbelow. Given that isapbelow = isapbelow and 
                        ! iheartabove = iheartbelow =  bm_sapl_heartabove*isapabove. If bm_sapl_heartbelow = 
                        ! bm_sapl_heartabove = 0.2, then isapabove = wood_init/2.4
                        bm_sapl(ivm,isapabove,icarbon) = wood_init / (2. + bm_sapl_heartabove + bm_sapl_heartbelow) 
                        bm_sapl(ivm,isapbelow,icarbon) = bm_sapl(ivm,isapabove,icarbon)
                        bm_sapl(ivm,iheartabove,icarbon) =  bm_sapl_heartabove * bm_sapl(ivm,isapabove,icarbon)
                        bm_sapl(ivm,iheartbelow,icarbon) =  bm_sapl_heartbelow * bm_sapl(ivm,isapbelow,icarbon)          
                        
                        ! Use the allometric relationships to calculate initial leaf and root mass 
                        bm_sapl(ivm,ileaf,icarbon) = ( bm_sapl(ivm,isapabove,icarbon) + &
                             bm_sapl(ivm,isapbelow,icarbon) ) * KF(ipts,ivm) /  ave_tree_height 
                        !+++CHECK+++
                        !How do we want to deal with water stress? wstress is accounted for through c0
!!$                           bm_sapl(ivm,iroot,icarbon) = bm_sapl(ivm,ileaf,icarbon) / ( KF(ipts,ivm) * c0_alloc(ivm) )
                        bm_sapl(ivm,iroot,icarbon) = bm_sapl(ivm,ileaf,icarbon) / &
                             ( KF(ipts,ivm) * c0_alloc(ivm) )
                        !+++++++++++
                        
                        ! Write initial values
                        IF (printlev>=4) THEN 
                           WRITE(numout,*) ' PFT type: ',ivm
                           WRITE(numout,*) '       root to sapwood tradeoff p :', c0_alloc(ivm)
                           WRITE(numout,*) 'height_init, dia_init, wood_init, ', &
                                ave_tree_height , dia_init, wood_init
                           WRITE(numout,*) 'pipe_density, ',pipe_density(ivm)
                        ENDIF
                        
                        !++++++ CHECK ++++++
                        ! The carbohydrate reserves do not seem to be set before this line.  This
                        ! is a problem since it then uses an unitililized value.  Therefore, I
                        ! will initilize it.
                        ! Should this really be zero? If so the code below is nonesensical and icarbres
                        ! could be omitted. nonesensical code = 2 * (bm_sapl(ivm,icir,icarbres,icarbon) + &
                        !        bm_sapl(ivm,icir,ileaf,icarbon) + bm_sapl(ivm,icir,iroot,icarbon))
                        bm_sapl(ivm,icarbres,icarbon)=zero
                        !++++++++++++++++++
                        
                        ! Pools that are defined in the same way for trees and grasses      
                        bm_sapl(ivm,ifruit,icarbon) = zero
                        
                        !+++CHECK+++
                        ! There is an inconsistency in the calculation - most pools are in gN 
                        ! but leaves is in gC. The correction is proposed, that implies that 
                        ! the parameter labile_reserve will need to be tuned
!!$                           bm_sapl(ivm,ilabile,icarbon) = labile_to_total * &
!!$                                (bm_sapl(ivm,ileaf,icarbon) / cn_leaf_prescribed(ivm)  + &
!!$                                fcn_root(ivm) * bm_sapl(ivm,iroot,icarbon) + fcn_wood(ivm) * &
!!$                                (bm_sapl(ivm,isapabove,icarbon) + bm_sapl(ivm,isapbelow,icarbon) + &
!!$                                bm_sapl(ivm,icarbres,icarbon)))

                        bm_sapl(ivm,ilabile,icarbon) = labile_to_total * (bm_sapl(ivm,ileaf,icarbon)  + &
                             fcn_root(ivm) * bm_sapl(ivm,iroot,icarbon) + fcn_wood_act * &
                             (bm_sapl(ivm,isapabove,icarbon) + bm_sapl(ivm,isapbelow,icarbon) + &
                             bm_sapl(ivm,icarbres,icarbon)))
                        !+++++++++++
                        
                        ! Avoid deciduous PFTs to have leaves out at establishment
                        ! Whether the saplings have leaves or don't have leaves the first year doesn't really matter
                        ! Either the approach is correct in the northern hemisphere or in the southern hemisphere. Note
                        ! that the resource limitation approach starts with the sapling having leaves.
                        ! Anyhow, a spin-up is needed to avoid issues with the initial conditions
                        IF ( pheno_type(ivm) .NE. 1 ) THEN
                           
                           ! Not evergreen. Deciduous PFTs now need to survive an extra year before bud burst. To 
                           ! ensure survival there are several options: (a) either the height of the initial 
                           ! vegetation is increased (this results in more reserves) or (b) the reserves could be
                           ! increased. The second option may result in numerical issues further down
                           ! the code as the optimal reserve level is calculated from the other biomass pools.
                           ! Also some initial tests showed that higher results simply resulted in more respiration.
                           ! Use taller trees to start with.
                           bm_sapl(ivm,icarbres,icarbon) = (bm_sapl(ivm,icarbres,icarbon) + &
                                bm_sapl(ivm,ileaf,icarbon) + bm_sapl(ivm,iroot,icarbon))
                           bm_sapl(ivm,ileaf,icarbon) = zero
                           bm_sapl(ivm,iroot,icarbon) = zero

                           ! When deciduous trees have no leaves they are senescent
                           senescence(ipts,ivm) = .TRUE.
                           
                        ELSE
                           
                           ! Initilize carbohydrate reserves for evergreen PFTs
                           bm_sapl(ivm,icarbres,icarbon) = zero

                           ! Evergreen trees never go into senescence
                           senescence(ipts,ivm) = .FALSE.
                           
                        ENDIF
                        
                        IF (printlev>=4) THEN
                           WRITE(numout,*) '       sapling biomass (gC):',ivm,ipts
                           WRITE(numout,*) '         leaves: (::bm_sapl(ivm,ileaf,icarbon))',&
                                bm_sapl(ivm,ileaf,icarbon)
                           WRITE(numout,*) '         sap above ground: (::bm_sapl(ivm,ispabove,icarbon)):',&
                                bm_sapl(ivm,isapabove,icarbon)
                           WRITE(numout,*) '         sap below ground: (::bm_sapl(ivm,isapbelow,icarbon))',&
                                bm_sapl(ivm,isapbelow,icarbon)
                           WRITE(numout,*) '         heartwood above ground: (::bm_sapl(ivm,iheartabove,icarbon))',&
                                bm_sapl(ivm,iheartabove,icarbon)
                           WRITE(numout,*) '         heartwood below ground: (::bm_sapl(ivm,iheartbelow,icarbon))',&
                                bm_sapl(ivm,iheartbelow,icarbon)
                           WRITE(numout,*) '         roots: (::bm_sapl(ivm,iroot,icarbon))',&
                                bm_sapl(ivm,iroot,icarbon)
                           WRITE(numout,*) '         fruits: (::bm_sapl(ivm,ifruit,icarbon))',&
                                bm_sapl(ivm,ifruit,icarbon)
                           WRITE(numout,*) '         carbohydrate reserve: (::bm_sapl(ivm,icarbres,icarbon))',&
                                bm_sapl(ivm,icarbres,icarbon)
                           WRITE(numout,*) '         labile reserve: (::bm_sapl(ivm,ilabile,icarbon))',&
                                bm_sapl(ivm,ilabile,icarbon)
                        ENDIF

                        DO k=1,nparts
                          IF (k.EQ.ileaf) THEN 
                             bm_sapl(ivm,k,initrogen)   = bm_sapl(ivm,k,icarbon)  / cn_leaf 
                             bm_sapl(ivm,k,iphosphorus) = bm_sapl(ivm,k,initrogen) / np_leaf 
                          ELSE IF (k.LT.iroot) THEN 
                             bm_sapl(ivm,k,initrogen)   = bm_sapl(ivm,k,icarbon)  / cn_wood
                             bm_sapl(ivm,k,iphosphorus) = bm_sapl(ivm,k,initrogen) / np_wood
                          ELSE 
                             bm_sapl(ivm,k,initrogen)   = bm_sapl(ivm,k,icarbon)  / cn_root
                             bm_sapl(ivm,k,iphosphorus) = bm_sapl(ivm,k,initrogen) / np_root
                          ENDIF
                        ENDDO
                    
                        
                        IF (printlev>=4) THEN
                           WRITE(numout,*)'Initial distribution, method 2',ivm
                           WRITE(numout,*)'Average trees height (m): ',ave_tree_height
                        ENDIF
                        
                        !! 2.7 Determine the biomass
                        !  I do this based on the biomass in each sapling and the number of trees in each
                        !  circumference class...we need the biomass in an average tree
                        biomass(ipts,ivm,:,:)= &
                             bm_sapl(ivm,:,:)*ind(ipts,ivm)

                        ! The light stress should be calculated making use of Pgap so it accounts for LAI
                        ! crown dimensions and tree distribution. However, this would be computationally 
                        ! expensive so we just use a first order estimate based on light attenuation model 
                        ! by Lambert-Beer. When LAI is low, a lot of light reaches the forest floor and so
                        ! KF should increase to make use of the available light by growing leaves
                        lstress_fac = exp(-biomass(ipts,ivm,ileaf,icarbon) * sla(ivm) * 0.5)
                        ! Causing large differences between first and second prescribe
                        delta_KF = ABS (KF(ipts,ivm) - ((k_latosa_adapt(ipts,ivm) + &
                          (lstress_fac(ipts,ivm) * (k_latosa_max(ivm)-k_latosa_min(ivm))))) / &
                          ( sla(ivm) * tree_ff(ivm) * pipe_density(ivm) ))
                        KF(ipts,ivm) = ((k_latosa_adapt(ipts,ivm) + &
                          (lstress_fac(ipts,ivm) * &
                          (k_latosa_max(ivm)-k_latosa_min(ivm)))) / &
                          ( sla(ivm) * tree_ff(ivm) * pipe_density(ivm) ))

                        IF(printlev>=4)THEN
                           WRITE(numout,*) 'prescribe delta_KF, ', delta_KF
                           WRITE(numout,*) 'prescribe lstress, ', lstress_fac
                        ENDIF
                     END DO

                     IF (printlev>=4) THEN
                        WRITE(numout,*)'Initial biomass distribution'
                        WRITE(numout,*)'End initial biomass distribution',ind(ipts,ivm)
                        WRITE(numout,*)"biomass(ipts,ivm,:,icarbon)",biomass(ipts,ivm,:,icarbon)
                        WRITE(numout,*)"biomass(ipts,ivm,:,initrogen)",biomass(ipts,ivm,:,initrogen)
                     END IF

                     IF (ivm .EQ. test_pft .AND. printlev>=4) THEN
                        WRITE(numout,*) 'Check prescribe'
                        WRITE(numout,*) 'stomate_prescribe::init ind  ',&
                             ipts,ivm,ind(ipts,ivm)
                     ENDIF

                  ! PFT is not present
                  ELSE

                     ! At the stand level
                     biomass(ipts,ivm,:,:) = zero
                     ind(ipts,ivm) = zero

                  ENDIF

                  ! Set leaf age classes, all leaves are current year leaves
                  leaf_frac(ipts,ivm,:,:) = zero
                  leaf_frac(ipts,ivm,:,1) = un

                  !+++CHECK+++
                  ! Set time since last beginning of growing season but only
                  ! for the first day of the whole simulation. When the model
                  ! is initialized when_growthinit is set to undef. In subsequent
                  ! time steps it should have a value. For trees without phenology
                  ! the growing season starts at the moment the PFT is prescribed
                  IF (when_growthinit(ipts,ivm) .EQ. undef) THEN
                     when_growthinit(ipts,ivm) = 200
                  ENDIF
                  !+++++++++++

                  ! Seasonal trees have no leaves at beginning
                  ! Saplings of evergreen trees have a leaf mass on day 1 and mass in the other components
                  ! saplings of deciduous trees have no leaves on day 1 but mass in the other components.
                  IF ( pheno_model(ivm) .NE. 'none' ) THEN

                     ! Add the carbon from the leaves to the reserve pool
                     biomass(ipts,ivm,icarbres,:) = biomass(ipts,ivm,icarbres,:) + biomass(ipts,ivm,ileaf,:)
                     biomass(ipts,ivm,ileaf,:) = zero
                     leaf_frac(ipts,ivm,:,1) = zero

                     !+++CHECK+++
                     ! Set time since last beginning of growing season but only
                     ! for the first day of the whole simulation. When the model
                     ! is initialized when_growthinit is set to undef. In subsequent
                     ! time steps it should have a value. The phenology module 
                     ! prevents leaf onset soon after senescence, by setting 
                     ! ::when_growthinit to a value, leaf offset 
                     ! will occur at the first opportunity 
!!$                     when_growthinit(ipts,ivm) = large_value
                     IF (when_growthinit(ipts,ivm) .EQ. undef) THEN
                        when_growthinit(ipts,ivm) = 200
                     ENDIF
                     !++++++++++++

                     ! Redundant, flag has already been set. When there are no leaves, the tree is in senescence
                     senescence(ipts,ivm) = .TRUE.

                  ENDIF ! pheno_model(ivm)

                  ! The biomass to build the saplings is taken from the atmosphere, keep track of
                  ! amount to calculate the C-balance closure
                  co2_to_bm(ipts,ivm) = co2_to_bm(ipts,ivm) + ( SUM(biomass(ipts,ivm,:,icarbon))  / dt )          
                  n_to_bm(ipts,ivm) = n_to_bm(ipts,ivm) + ( SUM(biomass(ipts,ivm,:,initrogen))  / dt )
                  p_to_bm(ipts,ivm) = p_to_bm(ipts,ivm) + ( SUM(biomass(ipts,ivm,:,iphosphorus))  / dt )

               ENDIF ! tree(ivm)

               !! 2.8 Initialize grassy PFTs
               !! Use veget_cov_max to check whether the PFT is present. If the PFT is present but it 
               !! has no biomass, prescribe its biomass (gC m-{2}). It is assumed that at day 1 
               !! grasses have all their biomass in the reserve pool. The criteria exclude crops.
               !! Crops are no longer prescribed but planted the day that begin_leaves is true.    

               IF ( ( .NOT. is_tree(ivm) ) .AND. &
                    ( veget_cov_max(ipts,ivm) .GT. min_stomate ) .AND. &
                    ( SUM( biomass(ipts,ivm,:,icarbon) ) .LE. min_stomate ) ) THEN

                  !+++TEMP+++
                  IF(printlev>=4 .AND. test_pft == ivm)THEN
                     WRITE(numout,*) 'We will prescribe a new vegetation, the old one died'
                  ENDIF
                  !++++++++++

                  ! For grasses we assume that an individual grass is 1 m2 of grass. This is set
                  ! in nmaxtrees in pft_parameters.f90. It could be set to another value but 
                  ! with the current code this should not have any meaning. The grassland
                  ! does not necessarily covers the whole 1 m2 so adjust for the canopy cover
                  ind(ipts,ivm) = nmaxtrees(ivm) * ha_to_m2 * canopy_cover(ivm)

                  ! now we generate the size of a single grass sapling...this was all taken from
                  ! stomate_data.f90...we do not deal with circumference classes for grasses
                  ! and crops, but we want to keep the arrays the same as for the trees so
                  ! we put the sapling information into the first circumference class
                
                  !+++CHECK+++
                  ! Calculate the sapwood to leaf mass in a similar way as has been done for trees.
                  ! For trees this approach had been justified by observations. For grasses such
                  ! justification is not supported by observations but we didn't try to find it.
                  ! Needs more work by someone interested in grasses. There might be a more elegant
                  ! solution making use of a well observed parameter. 
!!$                  k_latosa(ipts,ivm) = k_latosa_min(ivm) + &
!!$                       (wstress_fac(ipts,ivm) * lstress_fac(ipts,ivm) * &
!!$                       (k_latosa_max(ivm)-k_latosa_min(ivm)))
!!$                  k_latosa(ipts,ivm) = wstress_fac(ipts,ivm) * (k_latosa_min(ivm) + &
!!$                       (lstress_fac(ipts,ivm) * &
!!$                       (k_latosa_max(ivm)-k_latosa_min(ivm))))
                  k_latosa(ipts,ivm) = (k_latosa_adapt(ipts,ivm) + &
!                       (lstress_fac(ipts,ivm) * &
                       (k_latosa_max(ivm)-k_latosa_min(ivm)))

                  ! The mass of the structural carbon relates to the mass of the leaves through
                  ! a prescribed parameter ::k_latosa
                  KF(ipts,ivm) = k_latosa(ipts,ivm)
                  !+++++++++++

                  ! Calculate leaf to root area    
                  LF(ipts,ivm) = c0_alloc(ivm) * KF(ipts,ivm)

                  !---TEMP---
                  IF(printlev>=4.AND. test_pft == ivm)THEN
                     WRITE(numout,*) 'KF, ', k_latosa(ipts,ivm), KF(ipts,ivm)
                     WRITE(numout,*) 'LF, c0_alloc, ', c0_alloc(ivm) * KF(ipts,ivm), &
                       c0_alloc(ivm)
                  ENDIF
                  !----------

                  ! initialize everything to make sure there are not random values floating around
                  ! for ncirc != 1
                  bm_sapl(ivm,:,:) = val_exp

                  ! Similar as for trees, the initial height of the vegetation was defined 
                  bm_sapl(ivm,ileaf,icarbon) = height_init_min(ivm) / lai_to_height(ivm) / sla(ivm)

                  ! Use allometric relationships to define the root mass based on leaf mass. Some 
                  ! sapwood mass is needed to store the reserves. An arbitrairy fraction of 5% was
                  ! used 
                  bm_sapl(ivm,iroot,icarbon) = bm_sapl(ivm,ileaf,icarbon) / LF(ipts,ivm)
                  bm_sapl(ivm,isapabove,icarbon) = bm_sapl(ivm,ileaf,icarbon) / KF(ipts,ivm)

                  ! Some of the biomass components that exist for trees are undefined for grasses
                  bm_sapl(ivm,isapbelow,icarbon) = zero
                  bm_sapl(ivm,iheartabove,icarbon) = zero
                  bm_sapl(ivm,iheartbelow,icarbon) = zero
                  bm_sapl(ivm,ifruit,icarbon) = zero
                  bm_sapl(ivm,icarbres,icarbon) = zero

                  ! Pools that are defined in the same way for trees and grasses      
                  bm_sapl(ivm,ilabile,icarbon) = labile_to_total * (bm_sapl(ivm,ileaf,icarbon) + &
                       fcn_root(ivm) * bm_sapl(ivm,iroot,icarbon) + fcn_wood_act * &
                       (bm_sapl(ivm,isapabove,icarbon) + bm_sapl(ivm,isapbelow,icarbon) + &
                       bm_sapl(ivm,icarbres,icarbon)))

                  ! Avoid deciduous PFTs to have leaves out at establishment
                  ! Whether the saplings have leaves or don't have leaves the first year doesn't really matter
                  ! Either the approach is correct in the northern hemisphere or in the southern hemisphere. Note
                  ! that the resource limitation approach starts with the sapling having leaves.
                  ! Anyhow, a spin-up is needed to avoid issues with the initial conditions
                  IF ( pheno_type(ivm) .NE. 1 ) THEN

                     ! Grasses/crops. Deciduous PFTs now need to survive an extra year before bud burst. To 
                     ! ensure survival there are several options: (a) either the height of the initial 
                     ! vegetation is increased (this results in more reserves) or (b) the reserves could be
                     ! increased. The second option may result in result in numerical issues further down
                     ! the code as the optimal reserve level is calculated from the other biomass pools
                     bm_sapl(ivm,icarbres,icarbon) = bm_sapl(ivm,icarbres,icarbon) + bm_sapl(ivm,ileaf,icarbon) + &
                          bm_sapl(ivm,iroot,icarbon) + bm_sapl(ivm,isapabove,icarbon)
                     ! there should be no tissues other than reserves
                     bm_sapl(ivm,ileaf,icarbon) = zero
                     bm_sapl(ivm,iroot,icarbon) = zero
                     bm_sapl(ivm,isapabove,icarbon) = zero

                     ! When there are no leaves, the crop/grass is in senescence
                     senescence(ipts,ivm) = .TRUE.

                  ELSE

                     ! Initilize carbohydrate reserves for evergreen PFTs
                     bm_sapl(ivm,icarbres,icarbon) = zero 

                     ! Evergreen plants never go into senescence
                     senescence(ipts,ivm) = .FALSE.

                  ENDIF
   

                  IF (printlev>=4) THEN
                     WRITE(numout,*) '       sapling biomass (gC):',1,ivm,ipts
                     WRITE(numout,*) '         leaves: (::bm_sapl(ivm,ileaf,icarbon))',&
                          bm_sapl(ivm,ileaf,icarbon)
                     WRITE(numout,*) '         sap above ground: (::bm_sapl(ivm,ispabove,icarbon)):',&
                          bm_sapl(ivm,isapabove,icarbon)
                     WRITE(numout,*) '         sap below ground: (::bm_sapl(ivm,isapbelow,icarbon))',&
                          bm_sapl(ivm,isapbelow,icarbon)
                     WRITE(numout,*) '         heartwood above ground: (::bm_sapl(ivm,iheartabove,icarbon))',&
                          bm_sapl(ivm,iheartabove,icarbon)
                     WRITE(numout,*) '         heartwood below ground: (::bm_sapl(ivm,iheartbelow,icarbon))',&
                          bm_sapl(ivm,iheartbelow,icarbon)
                     WRITE(numout,*) '         roots: (::bm_sapl(ivm,iroot,icarbon))',&
                          bm_sapl(ivm,iroot,icarbon)
                     WRITE(numout,*) '         fruits: (::bm_sapl(ivm,ifruit,icarbon))',&
                          bm_sapl(ivm,ifruit,icarbon)
                     WRITE(numout,*) '         carbohydrate reserve: (::bm_sapl(ivm,icarbres,icarbon))',&
                          bm_sapl(ivm,icarbres,icarbon)
                     WRITE(numout,*) '         labile reserve: (::bm_sapl(ivm,ilabile,icarbon))',&
                          bm_sapl(ivm,ilabile,icarbon)
                  ENDIF

                  DO k=1,nparts 
                     IF (k.EQ.ileaf) THEN 
                        bm_sapl(ivm,k,initrogen)   = bm_sapl(ivm,k,icarbon)  / cn_leaf 
                        bm_sapl(ivm,k,iphosphorus) = bm_sapl(ivm,k,initrogen) / np_leaf 
                     ELSE IF (k.LT.iroot) THEN 
                        bm_sapl(ivm,k,initrogen)   = bm_sapl(ivm,k,icarbon)  / cn_wood 
                        bm_sapl(ivm,k,iphosphorus) = bm_sapl(ivm,k,initrogen) / np_wood
                     ELSE 
                        bm_sapl(ivm,k,initrogen)   = bm_sapl(ivm,k,icarbon)  / cn_root 
                        bm_sapl(ivm,k,iphosphorus) = bm_sapl(ivm,k,initrogen) / np_root 
                     ENDIF
                  ENDDO

                  ! Write initial values
                  IF (printlev>=4) THEN
                     WRITE(numout,*) '       root to sapwood tradeoff (LF) : ', c0_alloc(ivm)
                     WRITE(numout,*) '       grass sapling biomass: ',bm_sapl(ivm,:,icarbon)
                  ENDIF

                  ! Initial biomass (g C m-2)
                  biomass(ipts,ivm,:,:) = bm_sapl(ivm,:,:) * ind(ipts,ivm)
                  IF (printlev>=4) THEN 
                     WRITE(numout,*) 'bm_grass, prescribe, ', biomass(ipts,ivm,:,icarbon)
                  ENDIF


                  ! Set leaf age classes -> all leaves will be current year leaves
                  leaf_frac(ipts,ivm,:,:) = zero
                  leaf_frac(ipts,ivm,:,1) = un

                  ! Set time since last beginning of growing season but only
                  ! for the first day of the whole simulation. When the model
                  ! is initialized when_growthinit is set to undef. In subsequent
                  ! time steps it should have a value.
!!$                  when_growthinit(ipts,ivm) = large_value
!                  IF (when_growthinit(ipts,ivm) .EQ. undef) THEN
!                     when_growthinit(ipts,ivm) = 200
!                  ENDIF
                  when_growthinit(ipts,ivm) = large_value
                  ! The biomass to build the saplings is taken from the atmosphere, keep track of
                  ! amount to calculate the C-balance closure
                  co2_to_bm(ipts,ivm) = co2_to_bm(ipts,ivm) + ( SUM(biomass(ipts,ivm,:,icarbon))  / dt )
                  n_to_bm(ipts,ivm) = n_to_bm(ipts,ivm) + (SUM(biomass(ipts,ivm,:,initrogen))  / dt )
                  p_to_bm(ipts,ivm) = p_to_bm(ipts,ivm) + (SUM(biomass(ipts,ivm,:,iphosphorus))  / dt )

               ENDIF ! .NOT. tree(ivm)

               !! 2.3 Declare PFT present 
               !! Now that the PFT has biomass it should be declared 'present'
               !! everywhere in that grid box. Assign some additional properties
               PFTpresent(ipts,ivm) = .TRUE.
               everywhere(ipts,ivm) = un
               age(ipts,ivm) = zero
               npp_longterm(ipts,ivm) = npp_longterm_init
               lm_lastyearmax(ipts,ivm) = zero

            ENDIF   ! not ok_dgvm  or agricultural

         ENDIF ! IF (veget_cov_max .GT. zero .AND. ind .EQ. zero)

      ENDDO ! loop over pixels

   ENDDO ! loop over PFTs

 !! 4. Calculate components of the mass balance

   !! 4.1 Calculate final biomass
   pool_end(:,:,:) = zero 
   DO ipar = 1,nparts
      DO iele = 1,nelements
         pool_end(:,:,iele) = pool_end(:,:,iele) + &
              (biomass(:,:,ipar,iele) * veget_cov_max(:,:))
      ENDDO
   ENDDO

   !! 4.2 Calculate mass balance
   check_intern(:,:,iatm2land,icarbon) = check_intern(:,:,iatm2land,icarbon) + &
        co2_to_bm(:,:) * veget_cov_max(:,:) * dt
   check_intern(:,:,iland2atm,icarbon) = -un * zero
   check_intern(:,:,ilat2out,icarbon) = zero
   check_intern(:,:,ilat2in,icarbon) = -un * zero
   check_intern(:,:,ipoolchange,icarbon) = -un * (pool_end(:,:,icarbon) - pool_start(:,:,icarbon))
   closure_intern = zero
   DO imbc = 1,nmbcomp
      closure_intern(:,:,icarbon) = closure_intern(:,:,icarbon) + check_intern(:,:,imbc,icarbon)
   ENDDO

   !! 4.3 Write outcome
   DO ipts=1,npts
      DO ivm=1,nvm
         IF(printlev>=4 .AND. test_pft == ivm .AND. ipts==test_grid)THEN
             WRITE(numout,*) ' ENDprescribe - total biomass(C), ',SUM(biomass(ipts,ivm,:,icarbon))
             WRITE(numout,*) ' ENDprescribe - total biomass(N), ',SUM(biomass(ipts,ivm,:,initrogen))
             WRITE(numout,*) ' ENDprescribe - total biomass(P), ',SUM(biomass(ipts,ivm,:,iphosphorus))
         ENDIF
         IF(ABS(closure_intern(ipts,ivm,icarbon)) .LE. min_stomate)THEN
            IF (ld_massbal) WRITE(numout,*) 'Mass balance closure in prescribe_prognostic'
         ELSE
            WRITE(numout,*) 'Error: mass balance is not closed in prescribe_prognostic'
            WRITE(numout,*) '   ipts,ivm; ', ipts,ivm
            WRITE(numout,*) '   Difference is, ', closure_intern(ipts,ivm,icarbon)
            WRITE(numout,*) '   pool_end,pool_start: ', pool_end(ipts,ivm,icarbon), pool_start(ipts,ivm,icarbon)
            WRITE(numout,*) '   check_intern,co2_to_bm,pool_end,veget_cov_max: ', &
                 check_intern(ipts,ivm,iatm2land,icarbon),co2_to_bm(ipts,ivm), veget_cov_max(ipts,ivm)
            IF(ld_stop)THEN
               CALL ipslerr_p (3,'prescribe_prognostic', 'Mass balance error.','','')
            ENDIF
         ENDIF
      ENDDO
   ENDDO

   firstcall_prescribe = .FALSE.   

   END SUBROUTINE prescribe


END MODULE stomate_prescribe
