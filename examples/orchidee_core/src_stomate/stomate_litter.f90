! =================================================================================================================================
! MODULE       : stomate_litter
!
! CONTACT      : orchidee-help _at_ listes.ipsl.fr
!
! LICENCE      : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF       Update litter and lignine content after litter fall and 
!! calculating litter decomposition.      
!!
!!\n DESCRIPTION: None
!!
!! RECENT CHANGE(S): None
!!
!! REFERENCE(S)	: None
!!
!! SVN          :
!! $HeadURL: svn://forge.ipsl.jussieu.fr/orchidee/branches/ORCHIDEE-CN-P/ORCHIDEE/src_stomate/stomate_litter.f90 $
!! $Date: 2021-04-28 15:30:16 +0200 (三, 2021-04-28) $
!! $Revision: 7166 $
!! \n
!_ ================================================================================================================================

MODULE stomate_litter

  ! modules used:

  USE ioipsl_para
  USE xios_orchidee
  USE stomate_data
  USE constantes
  USE constantes_soil
  USE pft_parameters

  IMPLICIT NONE

  ! private & public routines

  PRIVATE
  PUBLIC littercalc,littercalc_clear, deadleaf

  LOGICAL, SAVE                        :: firstcall_litter = .TRUE.       !! first call
!$OMP THREADPRIVATE(firstcall_litter)

CONTAINS

!! ================================================================================================================================
!! SUBROUTINE   : littercalc_calc
!!
!!\BRIEF        Set the flag ::firstcall_litter to .TRUE. and as such activate section
!! 1.1 of the subroutine littercalc (see below).
!!
!! DESCRIPTION  : None
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S) : None
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE littercalc_clear
    firstcall_litter =.TRUE.
  END SUBROUTINE littercalc_clear


!! ================================================================================================================================
!! SUBROUTINE   : littercalc
!!
!!\BRIEF        Calculation of the litter decomposition and therefore of the 
!! heterotrophic respiration from litter following Parton et al. (1987).
!!
!! DESCRIPTION  : The littercal routine splits the litter in 4 pools: 
!! aboveground metaboblic, aboveground structural, belowground metabolic and 
!! belowground structural. the fraction (F) of plant material going to metabolic 
!! and structural is defined following Parton et al. (1987)
!! \latexonly
!! \input{littercalc1.tex}
!! \endlatexonly
!! \n
!! where L is the lignin content of the plant carbon pools considered and CN 
!! its CN ratio. L and CN are fixed parameters for each plant carbon pools,
!! therefore it is the ratio between each plant carbon pool within a PFT, which
!! controlled the part of the total litter, that will be considered as
!! recalcitrant (i.e. structural litter) or labile (i.e. metabolic litter).\n  
!! 
!! The routine calculates the fraction of aboveground litter which is metabolic
!! or structural (the litterpart variable) which is then used in lpj_fire.f90.\n 
!! 
!! In the section 2, the routine calculate the new plant material entering the
!! litter pools by phenological death of plants organs (corresponding to the
!! variable turnover) and by fire, herbivory and others non phenological causes
!! (variable bm_to_litter). This calculation is first done for each PFT and then
!! the values calculated for each PFT are added up. Following the same approach
!! the lignin content of the total structural litter is calculated and will be
!! then used as a factor control of the decomposition of the structural litter
!! (lignin_struc) in the section 5.1.2. A test is performed to avoid that we add
!! more lignin than structural litter. Finally, the variable litterpart is
!! updated.\n
!! 
!! In the section 3 and 4 the temperature and the moisture controlling the
!! decomposition are calculated for above and belowground. For aboveground
!! litter, air temperature and litter moisture are calculated in sechiba and used 
!! directly. For belowground, soil temperature and moisture are also calculated 
!! in sechiba but are modulated as a function of the soil depth. The modulation 
!! is a multiplying factor exponentially distributed between 0 (in depth) and 1
!! in surface.\n  
!! 
!! Then, in the section 5, the routine calculates the structural litter decomposition 
!! (C) following first order kinetics following Parton et al. (1987).
!! \latexonly
!! \input{littercalc2.tex}
!! \endlatexonly
!! \n
!! with k the decomposition rate of the structural litter. 
!! k corresponds to
!! \latexonly
!! \input{littercalc3.tex}
!! \endlatexonly
!! \n
!! with littertau the turnover rate, T a function of the temperature and M a function of
!! the moisture described below.\n
!!  
!! Then, the fraction of dead leaves (DL) composed by aboveground structural litter is
!! calculated as following
!! \latexonly
!! \input{littercalc4.tex}
!! \endlatexonly
!! \n
!! with k the decomposition rate of the structural litter previously
!! described.\n
!!
!! In the section 5.1, the fraction of decomposed structural litter
!! incorporated to the soil (Input) and its associated heterotrophic respiration are
!! calculated. For structural litter, the C decomposed could go in the active
!! soil carbon pool or in the slow carbon, as described in 
!! stomate_soilcarbon.f90.\n
!! \latexonly
!! \input{littercalc5.tex}
!! \endlatexonly
!! \n
!! with f a parameter describing the fraction of structural litter incorporated
!! into the considered soil carbon pool, C the amount of litter decomposed and L 
!! the amount of lignin in the litter. The litter decomposed which is not
!! incorporated into the soil is respired.\n
!!
!! In the section 5.2, the fraction of decomposed metabolic litter
!! incorporated to the soil and its associated heterotrophic respiration are
!! calculated with the same approaches presented for 5.1 but no control factor
!! depending on the lignin content are used.\n
!! 
!! In the section 6 the dead leaf cover is calculated through a call to the 
!! deadleaf subroutine presented below.\n
!!
!! In the section 7, if the flag SPINUP_ANALYTIC is set to true, we fill MatrixA
!! and VectorB following Lardy(2011).
!!
!! MAIN OUTPUT VARIABLES: ::deadleaf_cover, ::resp_hetero_litter, ::soilcarbon_input, 
!! ::control_temp, ::control_moist
!!
!! REFERENCES:
!! - Parton, WJ, Schimel, DS, Cole, CV, and Ojima, DS. 1987. Analysis
!! of factors controlling soil organic matter levels in Great Plains
!! grasslands. Soil Science Society of America journal (USA)
!! (51):1173-1179.
!! - Lardy, R, et al., A new method to determine soil organic carbon equilibrium,
!! Environmental Modelling & Software (2011), doi:10.1016|j.envsoft.2011.05.016
!!
!! FLOWCHART    :
!! \latexonly
!! \includegraphics(scale=0.5){littercalcflow.jpg}
!! \endlatexonly
!! \n
!_ ================================================================================================================================

  SUBROUTINE littercalc (npts, &
       turnover, bm_to_litter, &
       veget_cov_max, tsurf, tsoil, soilhum, litterhum, soil_n_min, &
       litter, &
       litter_avail, litter_not_avail, litter_avail_frac, &
       dead_leaves, lignin_struc, &
       lignin_wood, n_mineralisation, p_mineralisation, &
       deadleaf_cover, resp_hetero_litter, &
       som_input, control_temp, control_moist, &
       MatrixA, VectorB, CN_target, CN_som_litter_longterm, tau_CN_longterm, &
       CP_target, CP_som_litter_longterm, &
       sla_calc,do_slow)

    !! 0. Variable and parameter declaration
    
    !! 0.1 Input variables

    INTEGER, INTENT(in)                                  :: npts               !! Domain size - number of grid pixels
    REAL, DIMENSION(npts,nvm,nparts,nelements), INTENT(in) :: turnover         !! Turnover rates of plant biomass 
                                                                                      !! @tex $(gC m^{-2} dt\_slow^{-1})$ @endtex 
    REAL, DIMENSION(npts,nvm,nparts,nelements), INTENT(in) :: bm_to_litter     !! Conversion of biomass to litter 
                                                                                      !! @tex $(gC m^{-2} dt\_slow^{-1})$ @endtex 
    REAL,DIMENSION(npts,nvm),INTENT(in)                  :: veget_cov_max      !! PFT "Maximal" coverage fraction of a PFT 
                                                                                      !! defined in the input vegetation map 
                                                                                      !! @tex $(m^2 m^{-2})$ @endtex 
    REAL, DIMENSION(npts), INTENT(in)                    :: tsurf              !! Temperature (K) at the surface
    REAL, DIMENSION(npts,nslm), INTENT(in)               :: tsoil              !! Soil temperature (K)
    REAL, DIMENSION(npts,nslm), INTENT(in)               :: soilhum            !! Daily soil humidity of each soil layer 
                                                                                      !! (unitless)
    REAL, DIMENSION(npts), INTENT(in)                    :: litterhum          !! Daily litter humidity (unitless)
    REAL, DIMENSION(npts,nvm,nnspec),INTENT(in)          :: soil_n_min         !! mineral nitrogen in the soil (gN/m**2)  
                                                                                      !! (first index=npts, second index=nvm, third index=nnspec)  
    REAL, DIMENSION(npts,nlitt,nvm), INTENT(in)          :: litter_not_avail   !! litter not edible for animal
    REAL,DIMENSION(npts,nvm),INTENT(in)                  :: sla_calc           !! leaf age-related SLA
    LOGICAL,INTENT(in)                                          :: do_slow            !! called only once per day

    !! 0.2 Output variables
    
    REAL, DIMENSION(npts), INTENT(out)                   :: deadleaf_cover     !! Fraction of soil covered by dead leaves 
                                                                                      !! over all PFTs (0-1, unitless)
    REAL, DIMENSION(npts,nvm), INTENT(out)               :: resp_hetero_litter !! Litter heterotrophic respiration. The unit
                                                                                      !! is given by m^2 of ground.  
                                                                                      !! @tex $(gC dt_sechiba one\_day^{-1}) m^{-2})$ @endtex 
    REAL, DIMENSION(npts,ncarb,nvm,nelements), INTENT(out) :: som_input        !! Quantity of Carbon (or Nitrogen) going into SOM pools
                                                                                      !! from litter decomposition. The unit is  
                                                                                      !! given by m^2 of ground 
                                                                                      !! @tex $(gC(orN) m^{-2} dt\_slow^{-1})$ @endtex 
    REAL, DIMENSION(npts,nlevs), INTENT(out)             :: control_temp       !! Temperature control of heterotrophic 
                                                                                      !! respiration, above and below (0-1, 
                                                                                      !! unitless)
    REAL, DIMENSION(npts,nlevs), INTENT(out)             :: control_moist      !! Moisture control of heterotrophic 
                                                                                      !! respiration (0.25-1, unitless)
    REAL, DIMENSION(npts,nvm,nbpools,nbpools), INTENT(out) :: MatrixA          !! Matrix containing the fluxes between the
                                                                                      !! carbon pools per sechiba time step 
                                                                                      !! @tex $(gC.m^2.day^{-1})$ @endtex
    REAL, DIMENSION(npts,nvm,nbpools), INTENT(out)         :: VectorB          !! Vector containing the litter increase per
                                                                                      !! sechiba time step
                                                                                      !! @tex $(gC m^{-2})$ @endtex
    REAL, DIMENSION(npts,nvm,ncarb), INTENT(out)           :: CN_target        !! C to N ratio of SOM flux from one pool to another (gC gN-1)
    REAL, DIMENSION(npts,nvm,ncarb), INTENT(out)           :: CP_target        !! C to P ratio of SOM flux from one pool to another (gC gP-1)
    REAL, DIMENSION(npts,nlitt,nvm), INTENT(out)           :: litter_avail_frac!! litter fraction that is edible for animal
   
    !! 0.3 Modified variables
    
    REAL, DIMENSION(npts,nlitt,nvm,nlevs,nelements), INTENT(inout) :: litter   !! Metabolic and structural litter,above and
                                                                                      !! below ground. The unit is given by m^2 of 
                                                                                      !! ground @tex $(gC m^{-2})$ @endtex
    REAL, DIMENSION(npts,nlitt,nvm), INTENT(inout)       :: litter_avail       !! edible litter for animal 
    REAL, DIMENSION(npts,nvm,nlitt), INTENT(inout)       :: dead_leaves        !! Dead leaves per ground unit area, per PFT,
                                                                                      !! metabolic and structural in 
                                                                                      !! @tex $(gC m^{-2})$ @endtex
    REAL, DIMENSION(npts,nvm,nlevs), INTENT(inout)       :: lignin_struc       !! Ratio Lignin content in structural litter,
                                                                                      !! above and below ground, (0-1, unitless)
    REAL, DIMENSION(npts,nvm,nlevs), INTENT(inout)       :: lignin_wood        !! Ratio Lignin/Carbon in woody litter,	
                                                                                      !! above and below ground, (0-1, unitless)
    REAL, DIMENSION(npts,nvm), INTENT(inout)             :: n_mineralisation   !! mineralisation of N (gN m-2 timestep-1)
    REAL, DIMENSION(npts,nvm), INTENT(inout)             :: p_mineralisation   !! mineralisation of P (gP m-2 timestep-1)

    REAL, DIMENSION(npts,nvm,nbpools), INTENT(inout)     :: CN_som_litter_longterm !! Longterm CN ratio of litter and som pools (gC/gN)
    REAL, INTENT(inout)                                  :: tau_CN_longterm        !! Counter used for calculating the longterm CN_ratio of som and litter pools [seconds]
    REAL, DIMENSION(npts,nvm,nbpools), INTENT(inout)     :: CP_som_litter_longterm !! Longterm NP ratio of litter and som pools (gN/gP)

    !! 0.4 Local variables
 
    REAL                                                 :: dt                 !! Number of sechiba(fast processes) time-step per day 
    REAL, ALLOCATABLE, SAVE, DIMENSION(:,:,:,:)         :: litterfrac         !! The fraction of leaves, wood, etc. that 
                                                                                      !! goes into metabolic and structural 
                                                                                      !! litterpools (0-1, unitless)
!$OMP THREADPRIVATE(litterfrac)
    REAL, SAVE, ALLOCATABLE, DIMENSION(:)                :: z_soil             !! Soil levels (m)
!$OMP THREADPRIVATE(z_soil)
    REAL, DIMENSION(npts)                                :: rpc                !! Integration constant for vertical root 
                                                                                      !! profiles (unitless)
    REAL, SAVE, DIMENSION(nlitt)                         :: litter_turn        !! Turnover time in litter pools (days)
!$OMP THREADPRIVATE(litter_turn)
    REAL, SAVE, DIMENSION(nlitt,ncarb,nlevs)             :: frac_soil          !! Fraction of litter that goes into soil 
                                                                                      !! (litter -> carbon, above and below). The
                                                                                      !! remaining part goes to the atmosphere
!$OMP THREADPRIVATE(frac_soil)
    REAL, DIMENSION(npts)                                :: tsoil_decomp       !! Temperature used for decompostition in 
                                                                                      !! soil (K)
    REAL, DIMENSION(npts)                                :: soilhum_decomp     !! Humidity used for decompostition in soil
                                                                                      !! (unitless)
    REAL, DIMENSION(npts)                                :: fd                 !! Fraction of structural or metabolic litter
                                                                                      !! decomposed (unitless)
    REAL, DIMENSION(npts,nelements)                      :: qd                 !! Quantity of structural or metabolic litter
                                                                                      !! decomposed @tex $(gC m^{-2})$ @endtex
    REAL, DIMENSION(npts,nvm,nlevs)                      :: old_struc          !! Old structural litter, above and below 
                                                                                      !! @tex $(gC m^{-2})$ @endtex
    REAL, DIMENSION(npts,nvm,nlevs)                      :: old_woody          !! Old woody litter, above and below 
                                                                                      !! @tex $(gC m^{-2})$ @endtex
    REAL, DIMENSION(npts,nlitt,nvm,nlevs,nelements)      :: litter_inc         !! Increase of metabolic and structural 
                                                                                      !! litter, above and below ground. The unit 
                                                                                      !! is given by m^2 of ground. 
                                                                                      !! @tex $(gC m^{-2})$ @endtex
    REAL, DIMENSION(npts,nvm,nlevs)                      :: lignin_struc_inc   !! Lignin increase in structural litter, 
                                                                                      !! above and below ground. The unit is given 
                                                                                      !! by m^2 of ground. 
                                                                                      !! @tex $(gC m^{-2})$ @endtex                                             
    REAL, DIMENSION(npts,nvm,nlevs)                      :: lignin_wood_inc    !! Lignin increase in woody litter, 
                                                                                      !! above and below ground. The unit is given 
                                                                                      !! by m^2 of ground. 
                                                                                      !! @tex $(gC m^{-2})$ @endtex
    REAL, DIMENSION(npts)                                :: zdiff_min          !! Intermediate field for looking for minimum
                                                                                      !! of what? this is not used in the code. 
                                                                                      !! [??CHECK] could we delete it?
    CHARACTER(LEN=10), DIMENSION(nlitt)                         :: litter_str         !! Messages to write output information about
                                                                                      !! the litter
    CHARACTER(LEN=22), DIMENSION(nparts)                        :: part_str           !! Messages to write output information about
                                                                                      !! the plant
    CHARACTER(LEN=7), DIMENSION(ncarb)                          :: carbon_str         !! Messages to write output information about
                                                                                      !! the soil carbon
    CHARACTER(LEN=5), DIMENSION(nlevs)                          :: level_str          !! Messages to write output information about
                                                                                      !! the level (aboveground or belowground litter)
    INTEGER                                              :: i,j,k,l,m          !! Indices (unitless)
    REAL, DIMENSION(npts,nvm)                            :: f_soil_n_min       !! soil_n_min function response used for defing C to N target ratios (-) 
    REAL, DIMENSION(npts,nvm)                            :: f_pnc              !! plant nitrogen concentration function response used for defining C to N target ratios (-)  
    INTEGER                                              :: itarget            !! target som pool 
    REAL, DIMENSION(npts,nvm,nparts)                     :: CN                 !! CN ratio of the litter pools 
    REAL, DIMENSION(npts,nelements)                      :: summ  
    INTEGER                                              :: ier                !! Error handling
    CHARACTER(LEN=2), DIMENSION(nelements)                     :: element_str          !! string suffix indicating element
    REAL, DIMENSION(npts,nlitt,ncarb,nvm,nelements)         :: som_input_comp  !! fluxes of som_input (litter source -> som target)
!_ ================================================================================================================================

    IF (printlev>=3) WRITE(numout,*) 'Entering littercalc'

  !! 1. Initialisations of the different fields during the first call of the routine 
    dt = dt_sechiba/one_day
    IF ( firstcall_litter ) THEN


       IF ( .NOT. ALLOCATED(litterfrac) ) THEN
          ALLOCATE(litterfrac(npts,nvm,nparts,nlitt), stat=ier)
          IF (ier /= 0) CALL ipslerr_p(3, 'littercalc', 'ALLOCATE litterfrac', &
                    'Not enough memory', 'npts * nvm * nparts * nlitt')
       ENDIF


       !! 1.1.4 residence times in litter pools (days)
       litter_turn(imetabolic)  = turn_metabolic / one_year      ! .5 years
       litter_turn(istructural) = turn_struct / one_year     ! 3 years
       litter_turn(iwoody)      = turn_woody / one_year    !!!!???? 30 years (2.45)
      
       !! 1.1.5 decomposition flux fraction that goes into soil
       !       (litter -> carbon, above and below)
       !       1-frac_soil goes into atmosphere
       frac_soil(:,:,:) = zero

       ! structural litter: lignin fraction goes into slow pool + respiration,
       !                    rest into active pool + respiration
       frac_soil(istructural,isurface,iabove) = frac_soil_struct_sua
       frac_soil(istructural,iactive,ibelow) = frac_soil_struct_ab
       frac_soil(istructural,islow,iabove) = frac_soil_struct_sa
       frac_soil(istructural,islow,ibelow) = frac_soil_struct_sb

       ! metabolic litter: all goes into active pool + respiration.
       !   Nothing into slow or passive pool.
       frac_soil(imetabolic,isurface,iabove) = frac_soil_metab_sua
       frac_soil(imetabolic,iactive,ibelow) = frac_soil_metab_ab

       
       !! 1.2 soil levels
       ALLOCATE(z_soil(0:nslm), stat=ier)
       IF ( ier /= 0 ) CALL ipslerr_p(3,'littercalc','Pb in allocate of z_soil','','')

       z_soil(0) = zero
       z_soil(1:nslm) = diaglev(1:nslm)

       
       !! 1.3 messages
       litter_str(imetabolic) = 'metabolic'
       litter_str(istructural) = 'structural'
       litter_str(iwoody) = 'woody'

       carbon_str(iactive) = 'active'
       carbon_str(isurface) = 'surface'
       carbon_str(islow) = 'slow'
       carbon_str(ipassive) = 'passive'

       level_str(iabove) = 'above'
       level_str(ibelow) = 'below'

       part_str(ileaf) = 'leaves'
       part_str(isapabove) = 'sap above ground'
       part_str(isapbelow) = 'sap below ground'
       part_str(iheartabove) = 'heartwood above ground'
       part_str(iheartbelow) = 'heartwood below ground'
       part_str(iroot) = 'roots'
       part_str(ifruit) = 'fruits'
       part_str(icarbres) = 'carbohydrate reserve'
       part_str(ilabile) = 'labile reserve'

       IF (printlev >=2) THEN
          WRITE(numout,*) 'litter:'

          WRITE(numout,*) '   > C/N ratios: '
          DO k = 1, nparts
             WRITE(numout,*) '       ', part_str(k), ': ',CN_fix(k)
          ENDDO

          WRITE(numout,*) '   > Lignine/C ratios: '
          DO k = 1, nparts
             WRITE(numout,*) '       ', part_str(k), ': ',LC(:,k)
          ENDDO

          WRITE(numout,*) '   > fraction of compartment that goes into litter: '
          ALLOCATE(litterfrac(npts,nvm,nparts,nlitt), stat=ier)
            DO k = 1, nparts
              DO m = 1, nlitt
                WRITE(numout,*) '       ', part_str(k), '-> ',litter_str(m), ':',litterfrac(:,:,k,m)
              ENDDO
          ENDDO

          WRITE(numout,*) '   > scaling depth for decomposition (m): ',z_decomp

          WRITE(numout,*) '   > litter decomposition flux fraction that really goes '
          WRITE(numout,*) '     into carbon pools (rest into the atmosphere):'
          DO m = 1, nlitt
             DO l = 1, nlevs
                DO k = 1, ncarb
                   WRITE(numout,*) '       ',litter_str(m),' ',level_str(l),' -> ',&
                        carbon_str(k),':', frac_soil(m,k,l)
                ENDDO
             ENDDO
          ENDDO
       END IF ! printlev >=2
       firstcall_litter = .FALSE.

    ENDIF


    !! 1.1.3 litter fractions:
    !!   what fraction of leaves, wood, etc. goes into metabolic and structural litterpools

    ! PFT 1 needs to be initialised for when everything is printed below
    litterfrac(:,1,:,:) = zero
    CN(:,:,:)           = zero

    DO k = 1, nparts
       WHERE((bm_to_litter(:,:,k,initrogen)+turnover(:,:,k,initrogen)).GT.min_stomate) 
          CN(:,:,k) = (bm_to_litter(:,:,k,icarbon)+turnover(:,:,k,icarbon))/ & 
               (bm_to_litter(:,:,k,initrogen)+turnover(:,:,k,initrogen))
       ELSEWHERE
          CN(:,:,k) = CN_fix(k)
       ENDWHERE
       DO j = 2,nvm

          IF ( ((k == isapabove) .OR. (k == isapbelow) .OR. (k == iheartabove) .OR. (k == iheartbelow)) &
               .AND. is_tree(j) ) THEN
             litterfrac(:,j,k,iwoody) = 1.
             litterfrac(:,j,k,imetabolic) = zero
             litterfrac(:,j,k,istructural) = zero
          ELSE 
             IF( (k == icarbres) .OR. (k == ilabile)) THEN
                litterfrac(:,j,k,imetabolic) = 1.0
             ELSE
                litterfrac(:,j,k,imetabolic) = MAX(metabolic_ref_frac - metabolic_LN_ratio  * LC(j,k) * CN(:,j,k),zero)
             ENDIF
             litterfrac(:,j,k,istructural) = 1. - litterfrac(:,j,k,imetabolic)
             litterfrac(:,j,k,iwoody) = 0.0
             
          ENDIF
          
       END DO
       
    END DO

    CN_target(:,:,:)=0.0

    ! C to N target ratios of differnt pools
    f_soil_n_min = MIN((soil_n_min(:,:,iammonium) + soil_n_min(:,:,initrate)), 2.)

    !DSG no flexibility
    f_soil_n_min(:,:) = SPREAD(fixed_fsoiln,NCOPIES=npts,DIM=1)

    ! CN ratio ranges between 3 and 15 for active pool
    ! Figure 4 of Parton et al. (1993) show lower CN values for active pool of 2 rather than 3
    CN_target(:,:,iactive)= CN_target_iactive_ref + CN_target_iactive_Nmin * f_soil_n_min(:,:)
    ! CN ratio ranges between 12 and 20 for slow pool
    CN_target(:,:,islow)= CN_target_islow_ref + CN_target_islow_Nmin * f_soil_n_min(:,:)
    ! CN ratio ranges between 7 and 10 for passive pool
    ! OCN uses a fixed value of 9. (don't know why)
    ! Figure 4 of Parton et al. (1993) show lower CN values for passive pool of 3 rather than 7
    !DSG: we keep the passive CN fixed (
    CN_target(:,:,ipassive)= CN_target_ipassive_ref + CN_target_ipassive_Nmin * f_soil_n_min(:,:)


    ! plant nitrogen content (%)
    f_pnc(:,:)=zero
    WHERE(litter(:,imetabolic,:,iabove,icarbon)+litter(:,istructural,:,iabove,icarbon).GT. min_stomate)
       f_pnc(:,:)=MIN((litter(:,imetabolic,:,iabove,initrogen)+litter(:,istructural,:,iabove,initrogen)) / &
            (2.*(litter(:,imetabolic,:,iabove,icarbon)+litter(:,istructural,:,iabove,icarbon))) * 100., 2.)
    ELSEWHERE
       f_pnc(:,:)=0.
    ENDWHERE

    CN_target(:,:,isurface) = CN_target_isurface_ref + CN_target_isurface_pnc * f_pnc(:,:)


    ! 1.1.4. CP stoichiometry 
    ! CURRENT STATUS: CP ratio is not flexible
    CP_target(:,:,:)=0.0

    ! C to P target ratios of different SOM pools

    ! we keep the target NP ratio of SOM fixed (biochemical mineralization will cause the actual NP ratio to vary):
     CP_target(:,:,iactive)  = CN_target(:,:,iactive  ) * NP_target_iactive_ref
     CP_target(:,:,islow)    = CN_target(:,:,islow    ) * NP_target_islow_ref 
     CP_target(:,:,ipassive) = CN_target(:,:,ipassive ) * NP_target_ipassive_ref 
     CP_target(:,:,isurface) = CN_target(:,:,isurface ) * NP_target_islow_ref 
    
    !! 1.4 set output to zero
    deadleaf_cover(:)       = zero
    resp_hetero_litter(:,:) = zero
    som_input(:,:,:,:)      = zero
    som_input_comp(:,:,:,:,:)      = zero
    
  !! 2. Add biomass to different litterpools (per m^2 of ground)
    
    !! 2.1 first, save old structural litter (needed for lignin fractions).
    !     above/below
    DO l = 1, nlevs !Loop over litter levels (above and below ground)
       DO m = 1,nvm !Loop over PFTs

          old_struc(:,m,l) = litter(:,istructural,m,l,icarbon)
          old_woody(:,m,l) = litter(:,iwoody,m,l,icarbon)
       ENDDO
    ENDDO

    
    !! 2.2 update litter, dead leaves, and lignin content in structural litter
    litter_inc(:,:,:,:,:) = zero
    lignin_struc_inc(:,:,:) = zero
    lignin_wood_inc(:,:,:) = zero

    DO j = 1,nvm !Loop over PFTs

       !! 2.2.1 litter
       DO k = 1, nlitt    !Loop over litter pools (metabolic and structural)

          DO l = 1, nelements  ! Loop over element pools (carbon and nitrogen) 
          !! 2.2.2 calculate litter increase (per m^2 of ground).
          !       Only a given fracion of fruit turnover is directly coverted into litter.
          !       Litter increase for each PFT, structural and metabolic, above/below
          litter_inc(:,k,j,iabove,l) =                                           &
               litterfrac(:,j,ileaf      ,k) * bm_to_litter(:,j,ileaf      ,l) + &
               litterfrac(:,j,isapabove  ,k) * bm_to_litter(:,j,isapabove  ,l) + &
               litterfrac(:,j,iheartabove,k) * bm_to_litter(:,j,iheartabove,l) + &
               litterfrac(:,j,ifruit     ,k) * bm_to_litter(:,j,ifruit     ,l) + &
               litterfrac(:,j,icarbres   ,k) * bm_to_litter(:,j,icarbres   ,l) + &
               litterfrac(:,j,ilabile    ,k) * bm_to_litter(:,j,ilabile    ,l) + &
               litterfrac(:,j,ileaf      ,k) * turnover(:,j,ileaf      ,l)     + &
               litterfrac(:,j,isapabove  ,k) * turnover(:,j,isapabove  ,l)     + &
               litterfrac(:,j,iheartabove,k) * turnover(:,j,iheartabove,l)     + &
               litterfrac(:,j,ifruit     ,k) * turnover(:,j,ifruit     ,l)     + &
               litterfrac(:,j,icarbres   ,k) * turnover(:,j,icarbres   ,l)     + & 
               litterfrac(:,j,ilabile    ,k) * turnover(:,j,ilabile    ,l)

          litter_inc(:,k,j,ibelow,l) = &
               litterfrac(:,j,isapbelow,k)   * bm_to_litter(:,j,isapbelow  ,l) + &
               litterfrac(:,j,iheartbelow,k) * bm_to_litter(:,j,iheartbelow,l) + &
               litterfrac(:,j,iroot,k)       * bm_to_litter(:,j,iroot      ,l) + &
               litterfrac(:,j,isapbelow,k)   * turnover(:,j,isapbelow  ,l)     + &
               litterfrac(:,j,iheartbelow,k) * turnover(:,j,iheartbelow,l)     + &
               litterfrac(:,j,iroot,k)       * turnover(:,j,iroot      ,l)    
             
          ENDDO
          !! 2.2.3 dead leaves, for soil cover.
          dead_leaves(:,j,k) = &
               dead_leaves(:,j,k) + &
               litterfrac(:,j,ileaf,k) * ( bm_to_litter(:,j,ileaf,icarbon) + turnover(:,j,ileaf,icarbon) )
       ENDDO

       
       !! 2.2.4 lignin increase in structural litter

       lignin_struc_inc(:,j,iabove) = &
            LC(j,ileaf)       * bm_to_litter(:,j,ileaf,icarbon) + &
            LC(j,isapabove)   * bm_to_litter(:,j,isapabove,icarbon) + &
            LC(j,iheartabove) * bm_to_litter(:,j,iheartabove,icarbon) + &
            LC(j,ifruit)      * bm_to_litter(:,j,ifruit,icarbon) + &
            LC(j,icarbres)    * bm_to_litter(:,j,icarbres,icarbon) + &
            LC(j,ilabile)     * bm_to_litter(:,j,ilabile,icarbon) + &
            LC(j,ileaf)       * turnover(:,j,ileaf,icarbon) + &
            LC(j,isapabove)   * turnover(:,j,isapabove,icarbon) + &
            LC(j,iheartabove) * turnover(:,j,iheartabove,icarbon) + &
            LC(j,ifruit)      * turnover(:,j,ifruit,icarbon) + &
            LC(j,icarbres)    * turnover(:,j,icarbres,icarbon) + &
            LC(j,ilabile)     * turnover(:,j,ilabile,icarbon)


       lignin_struc_inc(:,j,ibelow) = &
            LC(j,isapbelow)   * bm_to_litter(:,j,isapbelow,icarbon) + &
            LC(j,iheartbelow) * bm_to_litter(:,j,iheartbelow,icarbon) + &
            LC(j,iroot)       * bm_to_litter(:,j,iroot,icarbon) + &
            LC(j,isapbelow)   * turnover(:,j,isapbelow,icarbon) + &
            LC(j,iheartbelow) * turnover(:,j,iheartbelow,icarbon) + &
            LC(j,iroot)       * turnover(:,j,iroot,icarbon)

       !! 2.2.4 lignin increase in woody litter
       
       lignin_wood_inc(:,j,iabove) = &
            LC(j,isapabove)   * bm_to_litter(:,j,isapabove,icarbon) + &
            LC(j,iheartabove) * bm_to_litter(:,j,iheartabove,icarbon) + &
            LC(j,isapabove)   * turnover(:,j,isapabove,icarbon) + &
            LC(j,iheartabove) * turnover(:,j,iheartabove,icarbon) 
       
       lignin_wood_inc(:,j,ibelow) = &
            LC(j,isapbelow)   * bm_to_litter(:,j,isapbelow,icarbon) + &
            LC(j,iheartbelow) * bm_to_litter(:,j,iheartbelow,icarbon) + &
            LC(j,isapbelow)   * turnover(:,j,isapbelow,icarbon) + &
            LC(j,iheartbelow) * turnover(:,j,iheartbelow,icarbon) 
       
    ENDDO


    !! 2.2.5 add new litter (struct/met, above/below)
    litter(:,:,:,:,:) = litter(:,:,:,:,:) + litter_inc(:,:,:,:,:)
    !gmjc for grazing litter
    IF (do_slow) THEN
      DO j=2,nvm
        !! grazed grassland or natural grassland with wild animal
        IF ((is_grassland_manag(j) .AND. is_grassland_grazed(j)) .OR. &
           ((.NOT. is_tree(j)) .AND. natural(j) .AND. &
             & (.NOT. is_grassland_cut(j)) .AND. (.NOT.is_grassland_grazed(j)))) THEN
          WHERE (litter(:,:,j,iabove,icarbon) .GE. litter_not_avail(:,:,j) &
             .AND. litter(:,:,j,iabove,icarbon) .GT. 0.0)
            litter_avail_frac(:,:,j) = (litter(:,:,j,iabove,icarbon) - litter_not_avail(:,:,j)) &
                  & / litter(:,:,j,iabove,icarbon)
          ! if litter not available equal to or larger than litter
          ELSEWHERE
            litter_avail_frac(:,:,j) = 0.0
          ENDWHERE
          IF (ANY(litter_avail_frac .LT. 0.0)) THEN
            WRITE (numout,*) 'frac error',litter(:,:,j,iabove,icarbon),litter_not_avail(:,:,j),litter_avail_frac(:,:,j)
            STOP 'error fraction litter available for grazing < 0'
          ENDIF
        ELSE
          litter_avail_frac(:,:,j) = 1.0
        ENDIF
      ENDDO
    ENDIF
    !end gmjc for grazing litter

    !! 2.2.6 for security: can't add more lignin than structural litter (above/below)
    DO l = 1, nlevs !Loop over litter levels (above and below ground)
       DO m = 1,nvm !Lopp over PFTs

          lignin_struc_inc(:,m,l) = &
               MIN( lignin_struc_inc(:,m,l), litter_inc(:,istructural,m,l,icarbon) )
          lignin_wood_inc(:,m,l) = &
               MIN( lignin_wood_inc(:,m,l), litter_inc(:,iwoody,m,l,icarbon) )

       ENDDO
    ENDDO


    !! 2.2.7 new lignin content: add old lignin and lignin increase, divide by 
    !!       total structural litter (above/below)

    WHERE ( litter(:,istructural,:,:,icarbon) .GT. min_stomate )
       lignin_struc(:,:,:) =  &
            ( lignin_struc(:,:,:)*old_struc(:,:,:) + lignin_struc_inc(:,:,:) ) / &
            litter(:,istructural,:,:,icarbon)
    ELSEWHERE
       lignin_struc(:,:,:) = zero
    ENDWHERE

    WHERE ( litter(:,iwoody,:,:,icarbon) .GT. min_stomate )
       lignin_wood(:,:,:) =  &
            ( lignin_wood(:,:,:)*old_woody(:,:,:) + lignin_wood_inc(:,:,:) ) / &
            litter(:,iwoody,:,:,icarbon)
    ELSEWHERE
       lignin_wood(:,:,:) = zero
    ENDWHERE
    

    !! 3. Temperature control on decay: Factor between 0 and 1

    !! 3.1 above: surface temperature
    control_temp(:,iabove) = control_temp_func (npts, tsurf)

    
    !! 3.2 below: convolution of temperature and decomposer profiles
    !!            (exponential decomposer profile supposed)
   
    !! 3.2.1 rpc is an integration constant such that the integral of the root profile is 1.
    rpc(:) = un / ( un - EXP( -z_soil(nslm) / z_decomp ) )

    !! 3.2.2 integrate over the nslm levels
    tsoil_decomp(:) = zero

    DO l = 1, nslm

       tsoil_decomp(:) = &
            tsoil_decomp(:) + tsoil(:,l) * rpc(:) * &
            ( EXP( -z_soil(l-1)/z_decomp ) - EXP( -z_soil(l)/z_decomp ) )

    ENDDO

    control_temp(:,ibelow) = control_temp_func (npts, tsoil_decomp)

  !! 4. Moisture control. Factor between 0 and 1
    
    !! 4.1 above the ground: litter humidity
    control_moist(:,iabove) = control_moist_func (npts, litterhum)

    !
    !! 4.2 below: convolution of humidity and decomposer profiles
    !            (exponential decomposer profile supposed)

    !! 4.2.1 rpc is an integration constant such that the integral of the root profile is 1.
    rpc(:) = un / ( un - EXP( -z_soil(nslm) / z_decomp ) )

    !! 4.2.2 integrate over the nslm levels
    soilhum_decomp(:) = zero

    DO l = 1, nslm !Loop over soil levels

       soilhum_decomp(:) = &
            soilhum_decomp(:) + soilhum(:,l) * rpc(:) * &
            ( EXP( -z_soil(l-1)/z_decomp ) - EXP( -z_soil(l)/z_decomp ) )

    ENDDO

    control_moist(:,ibelow) = control_moist_func (npts, soilhum_decomp)

  !! 5. fluxes from litter to carbon pools and respiration

    DO l = 1, nlevs !Loop over litter levels (above and below ground)
       DO m = 1,nvm !Loop over PFTs

          IF (l.EQ.iabove)THEN 
             itarget=isurface 
          ELSE 
             itarget=iactive 
          ENDIF
          !! 5.1 structural litter: goes into active and slow carbon pools + respiration

          !! 5.1.1 total quantity of structural litter which is decomposed
          fd(:) = dt*litter_turn(istructural) * &
               control_temp(:,l) * control_moist(:,l) * exp( -litter_struct_coef * lignin_struc(:,m,l) )

          DO k = 1,nelements
             
             qd(:,k) = litter(:,istructural,m,l,k) * fd(:)

          END DO

          litter(:,istructural,m,l,:) = litter(:,istructural,m,l,:) - qd(:,:)
          n_mineralisation(:,m) = n_mineralisation(:,m) + qd(:,initrogen) 
          p_mineralisation(:,m) = p_mineralisation(:,m) + qd(:,iphosphorus)

          !! 5.1.2 decompose same fraction of structural part of dead leaves. Not exact
          !!       as lignine content is not the same as that of the total structural litter.
          ! to avoid a multiple (for ibelow and iabove) modification of dead_leaves,
          ! we do this test to do this calcul only ones in 1,nlev loop
          if (l == iabove)  dead_leaves(:,m,istructural) = dead_leaves(:,m,istructural) * ( un - fd(:) )

          !! 5.1.3 non-lignin fraction of structural litter goes into
          !!       active (or surface) carbon pool + respiration 
          som_input(:,itarget,m,icarbon) = som_input(:,itarget,m,icarbon) + & 
               frac_soil(istructural,itarget,l) * qd(:,icarbon) * ( 1. - lignin_struc(:,m,l) ) / dt 

          som_input_comp(:,istructural,itarget,m,icarbon) = som_input_comp(:,istructural,itarget,m,icarbon) + &
               frac_soil(istructural,itarget,l) * & 
               qd(:,icarbon) * ( 1. - lignin_struc(:,m,l) ) / dt
          som_input_comp(:,istructural,itarget,m,initrogen) = som_input_comp(:,istructural,itarget,m,initrogen) + &
               frac_soil(istructural,itarget,l) * & 
               qd(:,initrogen) * ( 1. - lignin_struc(:,m,l) ) / dt
          som_input_comp(:,istructural,itarget,m,iphosphorus) = som_input_comp(:,istructural,itarget,m,iphosphorus) + &
               frac_soil(istructural,itarget,l) * &
               qd(:,iphosphorus) * ( 1. - lignin_struc(:,m,l) ) / dt
          
          !BE CAREFUL: Here resp_hetero_litter is divided by dt to have a value which corresponds to 
          ! the sechiba time step but then in stomate.f90 resp_hetero_litter is multiplied by dt. 
          ! Perhaps it could be simplified. Moreover, we must totally adapt the routines to the dtradia/one_day 
          ! time step and avoid some constructions that could create bug during future developments. 
          resp_hetero_litter(:,m) = resp_hetero_litter(:,m) + &
               ( 1. - frac_soil(istructural,itarget,l) ) * qd(:,icarbon) * &
               ( 1. - lignin_struc(:,m,l) ) / dt

          !! 5.1.4 lignin fraction of structural litter goes into
          !!       slow carbon pool + respiration
          som_input(:,islow,m,icarbon) = som_input(:,islow,m,icarbon) + &
               frac_soil(istructural,islow,l) * qd(:,icarbon) * lignin_struc(:,m,l) / dt

          som_input_comp(:,istructural,islow,m,icarbon) = som_input_comp(:,istructural,islow,m,icarbon) + & 
               frac_soil(istructural,islow,l) * &
               qd(:,icarbon) * lignin_struc(:,m,l) / dt
          som_input_comp(:,istructural,islow,m,initrogen) = som_input_comp(:,istructural,islow,m,initrogen) + &
          frac_soil(istructural,islow,l) * &
               qd(:,initrogen) * lignin_struc(:,m,l) / dt
          som_input_comp(:,istructural,islow,m,iphosphorus) = som_input_comp(:,istructural,islow,m,iphosphorus) + &
          frac_soil(istructural,islow,l) * &
               qd(:,iphosphorus) * lignin_struc(:,m,l) / dt

      !BE CAREFUL: Here resp_hetero_litter is divided by dt to have a value which corresponds to
      ! the sechiba time step but then in stomate.f90 resp_hetero_litter is multiplied by dt.
      ! Perhaps it could be simplified. Moreover, we must totally adapt the routines to the dt_sechiba/one_day
      ! time step and avoid some constructions that could create bug during future developments.
          resp_hetero_litter(:,m) = resp_hetero_litter(:,m) + &
               ( 1. - frac_soil(istructural,islow,l) ) * qd(:,icarbon) * lignin_struc(:,m,l) / dt

          
          !! 5.2 metabolic litter goes into active carbon pool + respiration
         
          !! 5.2.1 total quantity of metabolic litter that is decomposed
          fd(:) = dt*litter_turn(imetabolic) * control_temp(:,l) * control_moist(:,l)

          DO k = 1,nelements
          
             qd(:,k) = litter(:,imetabolic,m,l,k) * fd(:)

          END DO

          litter(:,imetabolic,m,l,:) = litter(:,imetabolic,m,l,:) - qd(:,:)

          n_mineralisation(:,m) = n_mineralisation(:,m) + qd(:,initrogen) 
          p_mineralisation(:,m) = p_mineralisation(:,m) + qd(:,iphosphorus)
          !! 5.2.2 decompose same fraction of metabolic part of dead leaves.
          !  to avoid a multiple (for ibelow and iabove) modification of dead_leaves,
          !  we do this test to do this calcul only ones in 1,nlev loop
          if (l == iabove)  dead_leaves(:,m,imetabolic) = dead_leaves(:,m,imetabolic) * ( 1. - fd(:) )

          !! 5.2.3 put decomposed litter into active (or surface) pool + respiration 
          som_input(:,itarget,m,icarbon) = som_input(:,itarget,m,icarbon) + & 
               frac_soil(imetabolic,itarget,l) * qd(:,icarbon) / dt 

          som_input_comp(:,imetabolic,itarget,m,icarbon) = som_input_comp(:,imetabolic,itarget,m,icarbon) + &
               frac_soil(imetabolic,itarget,l) * &
               qd(:,icarbon) / dt
          som_input_comp(:,imetabolic,itarget,m,initrogen) = som_input_comp(:,imetabolic,itarget,m,initrogen) + &
               frac_soil(imetabolic,itarget,l) * &
               qd(:,initrogen) / dt
          som_input_comp(:,imetabolic,itarget,m,iphosphorus) = som_input_comp(:,imetabolic,itarget,m,iphosphorus) + &
          frac_soil(imetabolic,itarget,l) * &
               qd(:,iphosphorus) / dt
          
          !BE CAREFUL: Here resp_hetero_litter is divided by dt to have a value which corresponds to 
          ! the sechiba time step but then in stomate.f90 resp_hetero_litter is multiplied by dt. 
          ! Perhaps it could be simplified. Moreover, we must totally adapt the routines to the dtradia/one_day 
          ! time step and avoid some constructions that could create bug during future developments.
          resp_hetero_litter(:,m) = resp_hetero_litter(:,m) + &
               ( 1. - frac_soil(imetabolic,itarget,l) ) * qd(:,icarbon) / dt

          !! 5.3 woody litter: goes into active and slow carbon pools + respiration
          
          !! 5.3.1 total quantity of woody litter which is decomposed
    
          fd(:) = dt*litter_turn(iwoody) * &
               control_temp(:,l) * control_moist(:,l) * EXP( -3. * lignin_wood(:,m,l) )
          
          DO k = 1,nelements 
             
             qd(:,k) = litter(:,iwoody,m,l,k) * fd(:)

          END DO
          
          litter(:,iwoody,m,l,:) = litter(:,iwoody,m,l,:) - qd(:,:)
          n_mineralisation(:,m) = n_mineralisation(:,m) + qd(:,initrogen)
          p_mineralisation(:,m) = p_mineralisation(:,m) + qd(:,iphosphorus)
          
          !! 5.3.2 non-lignin fraction of woody litter goes into
          !!       active/structural carbon pool + respiration (per time unit)
          
          som_input(:,itarget,m,icarbon) = som_input(:,itarget,m,icarbon) + &
               frac_soil(istructural,itarget,l) * qd(:,icarbon) * ( 1. - lignin_wood(:,m,l) ) / dt

          som_input_comp(:,iwoody,itarget,m,icarbon) = som_input_comp(:,iwoody,itarget,m,icarbon) + &
               frac_soil(istructural,itarget,l) * &
               qd(:,icarbon) * ( 1. - lignin_wood(:,m,l) ) / dt
          som_input_comp(:,iwoody,itarget,m,initrogen) = som_input_comp(:,iwoody,itarget,m,initrogen) + &
               frac_soil(istructural,itarget,l) * &
               qd(:,initrogen) * ( 1. - lignin_wood(:,m,l) ) / dt
          som_input_comp(:,iwoody,itarget,m,iphosphorus) = som_input_comp(:,iwoody,itarget,m,iphosphorus) + & 
               frac_soil(istructural,itarget,l) * &
               qd(:,iphosphorus) * ( 1. - lignin_wood(:,m,l) ) / dt

    
          resp_hetero_litter(:,m) = resp_hetero_litter(:,m) + &
               ( 1. - frac_soil(istructural,itarget,l) ) * qd(:,icarbon) * &
               ( 1. - lignin_wood(:,m,l) ) / dt
    
          !! 5.3.3 lignin fraction of woody litter goes into
          !!       slow carbon pool + respiration (per time unit)
          
          som_input(:,islow,m,icarbon) = som_input(:,islow,m,icarbon) + &
               frac_soil(istructural,islow,l) * qd(:,icarbon) * lignin_wood(:,m,l) / dt

          som_input_comp(:,iwoody,islow,m,icarbon) = som_input_comp(:,iwoody,islow,m,icarbon) + frac_soil(istructural,islow,l) * &
               qd(:,icarbon) * lignin_wood(:,m,l) / dt
          som_input_comp(:,iwoody,islow,m,initrogen) = som_input_comp(:,iwoody,islow,m,initrogen) + frac_soil(istructural,islow,l) * &
               qd(:,initrogen) * lignin_wood(:,m,l) / dt
          som_input_comp(:,iwoody,islow,m,iphosphorus) = som_input_comp(:,iwoody,islow,m,iphosphorus) + frac_soil(istructural,islow,l) * &
               qd(:,iphosphorus) * lignin_wood(:,m,l) / dt
    
          resp_hetero_litter(:,m) = resp_hetero_litter(:,m) + &
               ( 1. - frac_soil(istructural,islow,l) ) * qd(:,icarbon) * lignin_wood(:,m,l) / dt


       ENDDO
    ENDDO

    ! Nitrogen flux from litter to SOM 
    som_input(:,iactive ,:,initrogen) = som_input(:,iactive ,:,icarbon) /CN_target(:,:,iactive) 
    som_input(:,islow   ,:,initrogen) = som_input(:,islow   ,:,icarbon) /CN_target(:,:,islow) 
    som_input(:,isurface,:,initrogen) = som_input(:,isurface,:,icarbon) /CN_target(:,:,isurface) 

    ! Nitrogen mineralization / immobilization :
    n_mineralisation(:,:) = n_mineralisation(:,:) - & 
         ( som_input(:,iactive,:,initrogen) + & 
           som_input(:,islow,:,initrogen)   + & 
           som_input(:,isurface,:,initrogen))*dt !! multiply by dt !! 

    ! Phosphorus flux from litter to SOM 
    som_input(:,iactive ,:,iphosphorus) = som_input(:,iactive ,:,icarbon) / CP_target(:,:,iactive)
    som_input(:,islow   ,:,iphosphorus) = som_input(:,islow   ,:,icarbon) / CP_target(:,:,islow)
    som_input(:,isurface,:,iphosphorus) = som_input(:,isurface,:,icarbon) / CP_target(:,:,isurface)

    ! Phosphorus mineralization / immobilization :
    p_mineralisation(:,:) = p_mineralisation(:,:)      - &
               ( som_input(:,iactive ,:,iphosphorus)   + &
                 som_input(:,islow   ,:,iphosphorus)   + &
                 som_input(:,isurface,:,iphosphorus))*dt !! multiply by dt !!

    
    DO m=1,nvm 
       summ(:,:)=zero 
       DO l = 1, nlevs 
          DO k = 1,nlitt 
             summ(:,:)=summ(:,:)+litter(:,k,m,l,:) 
          ENDDO
       ENDDO
       WHERE(summ(:,icarbon).LE.min_stomate) 
       
          resp_hetero_litter(:,m) = resp_hetero_litter(:,m) + & 
               summ(:,icarbon) 
          n_mineralisation(:,m) = n_mineralisation(:,m) + & 
               summ(:,initrogen) 
          p_mineralisation(:,m)   = p_mineralisation(:,m) + &
              summ(:,iphosphorus)

          litter(:,imetabolic ,m,iabove,icarbon)     = zero
          litter(:,istructural,m,iabove,icarbon)     = zero
          litter(:,iwoody     ,m,iabove,icarbon)     = zero
          litter(:,imetabolic ,m,ibelow,icarbon)     = zero
          litter(:,istructural,m,ibelow,icarbon)     = zero
          litter(:,iwoody     ,m,ibelow,icarbon)     = zero

          litter(:,imetabolic ,m,iabove,initrogen)   = zero
          litter(:,istructural,m,iabove,initrogen)   = zero
          litter(:,iwoody     ,m,iabove,initrogen)   = zero
          litter(:,imetabolic ,m,ibelow,initrogen)   = zero
          litter(:,istructural,m,ibelow,initrogen)   = zero
          litter(:,iwoody     ,m,ibelow,initrogen)   = zero

          litter(:,imetabolic ,m,iabove,iphosphorus) = zero
          litter(:,istructural,m,iabove,iphosphorus) = zero
          litter(:,iwoody     ,m,iabove,iphosphorus) = zero
          litter(:,imetabolic ,m,ibelow,iphosphorus) = zero
          litter(:,istructural,m,ibelow,iphosphorus) = zero
          litter(:,iwoody     ,m,ibelow,iphosphorus) = zero
       ENDWHERE
    ENDDO

    ! For YANN
    CALL xios_orchidee_send_field("N_MINERALISATION_LIT",n_mineralisation(:,:)/dt)
    CALL xios_orchidee_send_field("P_MINERALISATION_LIT",p_mineralisation(:,:)/dt)
    ! resp_hetero_litter has been divided by dt above
    CALL xios_orchidee_send_field("HETERO_RESP_LIT",resp_hetero_litter(:,:))

    !! JC ADD for Yilong Wang 15N
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

       ! bm_to_litter into different litter pools
       CALL xios_orchidee_send_field("BM_LEAF_MET_AB"//TRIM(element_str(l)),litterfrac(:,:,ileaf,imetabolic) * & 
            bm_to_litter(:,:,ileaf,l)/dt)
       CALL xios_orchidee_send_field("BM_LEAF_STR_AB"//TRIM(element_str(l)),litterfrac(:,:,ileaf,istructural) * & 
            bm_to_litter(:,:,ileaf,l)/dt)
       CALL xios_orchidee_send_field("BM_SAP_AB_MET_AB"//TRIM(element_str(l)),litterfrac(:,:,isapabove,imetabolic) * & 
            bm_to_litter(:,:,isapabove,l)/dt)
       CALL xios_orchidee_send_field("BM_SAP_AB_STR_AB"//TRIM(element_str(l)),litterfrac(:,:,isapabove,istructural) * & 
            bm_to_litter(:,:,isapabove,l)/dt)
       CALL xios_orchidee_send_field("BM_SAP_AB_WOD_AB"//TRIM(element_str(l)),litterfrac(:,:,isapabove,iwoody) * & 
            bm_to_litter(:,:,isapabove,l)/dt)
       CALL xios_orchidee_send_field("BM_HEART_AB_MET_AB"//TRIM(element_str(l)),litterfrac(:,:,iheartabove,imetabolic) * & 
            bm_to_litter(:,:,iheartabove,l)/dt)
       CALL xios_orchidee_send_field("BM_HEART_AB_STR_AB"//TRIM(element_str(l)),litterfrac(:,:,iheartabove,istructural) * & 
            bm_to_litter(:,:,iheartabove,l)/dt)
       CALL xios_orchidee_send_field("BM_HEART_AB_WOD_AB"//TRIM(element_str(l)),litterfrac(:,:,iheartabove,iwoody) * & 
            bm_to_litter(:,:,iheartabove,l)/dt)

       CALL xios_orchidee_send_field("BM_SAP_BE_MET_BE"//TRIM(element_str(l)),litterfrac(:,:,isapbelow,imetabolic) * & 
            bm_to_litter(:,:,isapbelow,l)/dt)
       CALL xios_orchidee_send_field("BM_SAP_BE_STR_BE"//TRIM(element_str(l)),litterfrac(:,:,isapbelow,istructural) * & 
            bm_to_litter(:,:,isapbelow,l)/dt)
       CALL xios_orchidee_send_field("BM_SAP_BE_WOD_BE"//TRIM(element_str(l)),litterfrac(:,:,isapbelow,iwoody) * & 
            bm_to_litter(:,:,isapbelow,l)/dt)
       CALL xios_orchidee_send_field("BM_HEART_BE_MET_BE"//TRIM(element_str(l)),litterfrac(:,:,iheartbelow,imetabolic) * & 
            bm_to_litter(:,:,iheartbelow,l)/dt)
       CALL xios_orchidee_send_field("BM_HEART_BE_STR_BE"//TRIM(element_str(l)),litterfrac(:,:,iheartbelow,istructural) * & 
            bm_to_litter(:,:,iheartbelow,l)/dt)
       CALL xios_orchidee_send_field("BM_HEART_BE_WOD_BE"//TRIM(element_str(l)),litterfrac(:,:,iheartbelow,iwoody) * & 
            bm_to_litter(:,:,iheartbelow,l)/dt)

       CALL xios_orchidee_send_field("BM_FRUIT_MET_AB"//TRIM(element_str(l)),litterfrac(:,:,ifruit,imetabolic) * & 
            bm_to_litter(:,:,ifruit,l)/dt)
       CALL xios_orchidee_send_field("BM_FRUIT_STR_AB"//TRIM(element_str(l)),litterfrac(:,:,ifruit,istructural) * & 
            bm_to_litter(:,:,ifruit,l)/dt)
       CALL xios_orchidee_send_field("BM_RES_MET_AB"//TRIM(element_str(l)),litterfrac(:,:,icarbres,imetabolic) * & 
            bm_to_litter(:,:,icarbres,l)/dt)
       CALL xios_orchidee_send_field("BM_LAB_MET_AB"//TRIM(element_str(l)),litterfrac(:,:,ilabile,imetabolic) * & 
            bm_to_litter(:,:,ilabile,l)/dt)
       CALL xios_orchidee_send_field("BM_ROOT_MET_BE"//TRIM(element_str(l)),litterfrac(:,:,iroot,imetabolic) * & 
            bm_to_litter(:,:,iroot,l)/dt) 
       CALL xios_orchidee_send_field("BM_ROOT_STR_BE"//TRIM(element_str(l)),litterfrac(:,:,iroot,istructural) * & 
            bm_to_litter(:,:,iroot,l)/dt)

       ! turnover into different litter pools
       CALL xios_orchidee_send_field("TURN_LEAF_MET_AB"//TRIM(element_str(l)),litterfrac(:,:,ileaf,imetabolic) * & 
            turnover(:,:,ileaf,l)/dt)
       CALL xios_orchidee_send_field("TURN_LEAF_STR_AB"//TRIM(element_str(l)),litterfrac(:,:,ileaf,istructural) * & 
            turnover(:,:,ileaf,l)/dt)
       CALL xios_orchidee_send_field("TURN_SAP_AB_MET_AB"//TRIM(element_str(l)),litterfrac(:,:,isapabove,imetabolic) * & 
            turnover(:,:,isapabove,l)/dt)
       CALL xios_orchidee_send_field("TURN_SAP_AB_STR_AB"//TRIM(element_str(l)),litterfrac(:,:,isapabove,istructural) * & 
            turnover(:,:,isapabove,l)/dt)
       CALL xios_orchidee_send_field("TURN_SAP_AB_WOD_AB"//TRIM(element_str(l)),litterfrac(:,:,isapabove,iwoody) * & 
            turnover(:,:,isapabove,l)/dt)
       CALL xios_orchidee_send_field("TURN_HEART_AB_MET_AB"//TRIM(element_str(l)),litterfrac(:,:,iheartabove,imetabolic) * & 
            turnover(:,:,iheartabove,l)/dt)
       CALL xios_orchidee_send_field("TURN_HEART_AB_STR_AB"//TRIM(element_str(l)),litterfrac(:,:,iheartabove,istructural) * & 
            turnover(:,:,iheartabove,l)/dt)
       CALL xios_orchidee_send_field("TURN_HEART_AB_WOD_AB"//TRIM(element_str(l)),litterfrac(:,:,iheartabove,iwoody) * & 
            turnover(:,:,iheartabove,l)/dt)
       CALL xios_orchidee_send_field("TURN_SAP_BE_MET_BE"//TRIM(element_str(l)),litterfrac(:,:,isapbelow,imetabolic) * & 
            turnover(:,:,isapbelow,l)/dt)
       CALL xios_orchidee_send_field("TURN_SAP_BE_STR_BE"//TRIM(element_str(l)),litterfrac(:,:,isapbelow,istructural) * & 
            turnover(:,:,isapbelow,l)/dt)
       CALL xios_orchidee_send_field("TURN_SAP_BE_WOD_BE"//TRIM(element_str(l)),litterfrac(:,:,isapbelow,iwoody) * & 
            turnover(:,:,isapbelow,l)/dt)
       CALL xios_orchidee_send_field("TURN_HEART_BE_MET_BE"//TRIM(element_str(l)),litterfrac(:,:,iheartbelow,imetabolic) * & 
            turnover(:,:,iheartbelow,l)/dt)
       CALL xios_orchidee_send_field("TURN_HEART_BE_STR_BE"//TRIM(element_str(l)),litterfrac(:,:,iheartbelow,istructural) * & 
            turnover(:,:,iheartbelow,l)/dt)
       CALL xios_orchidee_send_field("TURN_HEART_BE_WOD_BE"//TRIM(element_str(l)),litterfrac(:,:,iheartbelow,iwoody) * & 
            turnover(:,:,iheartbelow,l)/dt)

       CALL xios_orchidee_send_field("TURN_FRUIT_MET_AB"//TRIM(element_str(l)),litterfrac(:,:,ifruit,imetabolic) * & 
            turnover(:,:,ifruit,l)/dt)
       CALL xios_orchidee_send_field("TURN_FRUIT_STR_AB"//TRIM(element_str(l)),litterfrac(:,:,ifruit,istructural) * & 
            turnover(:,:,ifruit,l)/dt)
       CALL xios_orchidee_send_field("TURN_RES_MET_AB"//TRIM(element_str(l)),litterfrac(:,:,icarbres,imetabolic) * & 
            turnover(:,:,icarbres,l)/dt)
       CALL xios_orchidee_send_field("TURN_LAB_MET_AB"//TRIM(element_str(l)),litterfrac(:,:,ilabile,imetabolic) * & 
            turnover(:,:,ilabile,l)/dt)
       CALL xios_orchidee_send_field("TURN_ROOT_MET_BE"//TRIM(element_str(l)),litterfrac(:,:,iroot,imetabolic) * & 
            turnover(:,:,iroot,l)/dt)
       CALL xios_orchidee_send_field("TURN_ROOT_STR_BE"//TRIM(element_str(l)),litterfrac(:,:,iroot,istructural) * & 
            turnover(:,:,iroot,l)/dt)
       ! fluxes from litter to soil pools is not possible to output directly
       ! they are rescaled to meet the target CN ratio
       ! Here we used the fluxes of carbon 
       CALL xios_orchidee_send_field("LIT_STR_ACTIVE"//TRIM(element_str(l)),som_input_comp(:,istructural,iactive,:,l))
       CALL xios_orchidee_send_field("LIT_STR_SURF"//TRIM(element_str(l)),som_input_comp(:,istructural,isurface,:,l))
       CALL xios_orchidee_send_field("LIT_STR_SLOW"//TRIM(element_str(l)),som_input_comp(:,istructural,islow,:,l))
       CALL xios_orchidee_send_field("LIT_MET_ACTIVE"//TRIM(element_str(l)),som_input_comp(:,imetabolic,iactive,:,l))
       CALL xios_orchidee_send_field("LIT_MET_SURF"//TRIM(element_str(l)),som_input_comp(:,imetabolic,isurface,:,l))
       CALL xios_orchidee_send_field("LIT_WOD_ACTIVE"//TRIM(element_str(l)),som_input_comp(:,iwoody,iactive,:,l))
       CALL xios_orchidee_send_field("LIT_WOD_SURF"//TRIM(element_str(l)),som_input_comp(:,iwoody,isurface,:,l))
       CALL xios_orchidee_send_field("LIT_WOD_SLOW"//TRIM(element_str(l)),som_input_comp(:,iwoody,islow,:,l))
       ! fluxes from liter to soil pools target pools only
       CALL xios_orchidee_send_field("SOM_INPUT_ACTIVE"//TRIM(element_str(l)),som_input(:,iactive,:,l))
       CALL xios_orchidee_send_field("SOM_INPUT_SLOW"//TRIM(element_str(l)),som_input(:,islow,:,l))
       CALL xios_orchidee_send_field("SOM_INPUT_PASSIVE"//TRIM(element_str(l)),som_input(:,ipassive,:,l))
       CALL xios_orchidee_send_field("SOM_INPUT_SURF"//TRIM(element_str(l)),som_input(:,isurface,:,l))
    ENDDO
  !! 6. calculate fraction of total soil covered by dead leaves

    CALL deadleaf (npts, veget_cov_max, dead_leaves, deadleaf_cover, &
                   sla_calc)

  !! 7. (Quasi-)Analytical Spin-up : Start filling MatrixA

    IF (spinup_analytic) THEN

       MatrixA(:,:,:,:) = zero
       VectorB(:,:,:) = zero
       
       
       DO m = 1,nvm

          !- MatrixA : carbon fluxes leaving the litter
          
          MatrixA(:,m,istructural_above,istructural_above)= - dt*litter_turn(istructural) * &
               control_temp(:,iabove) * control_moist(:,iabove) * exp( -litter_struct_coef * lignin_struc(:,m,iabove) )
          
          MatrixA(:,m,istructural_below,istructural_below) = - dt*litter_turn(istructural) * &
               control_temp(:,ibelow) * control_moist(:,ibelow) * exp( -litter_struct_coef * lignin_struc(:,m,ibelow) )
          
          MatrixA(:,m,imetabolic_above,imetabolic_above) = - dt*litter_turn(imetabolic) * & 
               control_temp(:,iabove) * control_moist(:,iabove)
          
          MatrixA(:,m,imetabolic_below,imetabolic_below) = - dt*litter_turn(imetabolic) * & 
               control_temp(:,ibelow) * control_moist(:,ibelow)
          
          ! Flux leaving the woody above litter pool :
          MatrixA(:, m, iwoody_above, iwoody_above) = - dt * litter_turn(iwoody) * control_temp(:,iabove) * &
               control_moist(:,iabove) * exp( -3. * lignin_wood(:,m,iabove) )

          ! Flux leaving the woody below litter pool :
          MatrixA(:, m, iwoody_below, iwoody_below) = - dt * litter_turn(iwoody) * control_temp(:,ibelow) * &
               control_moist(:,ibelow) * exp( -3. * lignin_wood(:,m,ibelow))

          ! Flux received by the carbon surface from the woody above litter pool :
          MatrixA(:, m, isurface_pool, iwoody_above) = frac_soil(istructural, isurface, iabove) * & 
               dt *litter_turn(iwoody) * & 
               control_temp(:,iabove) * &
               control_moist(:,iabove) *  &
               exp( -3. * lignin_wood(:,m,iabove) ) * ( 1. -  lignin_wood(:,m,iabove) ) 

          ! Flux received by the carbon active from the woody below litter pool :
          MatrixA(:, m, iactive_pool, iwoody_below) = frac_soil(istructural, iactive, ibelow) * & 
               dt *litter_turn(iwoody) * &
               control_temp(:,ibelow) * &
               control_moist(:,ibelow) * &
               exp( -3. * lignin_wood(:,m,ibelow) ) * ( 1. -  lignin_wood(:,m,ibelow) ) 

          ! Flux received by the carbon slow from the woody above litter pool :
          MatrixA(:, m, islow_pool, iwoody_above) = frac_soil(istructural, islow, iabove) * &
               dt *litter_turn(iwoody) * &
               control_temp(:,iabove) * &
               control_moist(:,iabove) * &
               exp( -3. * lignin_wood(:,m,iabove) ) * lignin_wood(:,m,iabove)

          ! Flux received by the carbon slow from the woody below litter pool :
          MatrixA(:, m, islow_pool, iwoody_below) =  frac_soil(istructural, islow, ibelow) * & 
               dt *litter_turn(iwoody) * & 
               control_temp(:,ibelow) * &
               control_moist(:,ibelow) * &
               exp( -3. * lignin_wood(:,m,ibelow) ) * lignin_wood(:,m,ibelow)

                    
          !- MatrixA : carbon fluxes between the litter and the pools (the rest of the matrix is filled in stomate_soilcarbon.f90)
          MatrixA(:,m,isurface_pool,istructural_above) = frac_soil(istructural,isurface,iabove) * &
               dt*litter_turn(istructural) * &                    
               control_temp(:,iabove) * control_moist(:,iabove) * & 
               exp( -litter_struct_coef * lignin_struc(:,m,iabove) ) * &
               ( 1. - lignin_struc(:,m,iabove) ) 
                         

          MatrixA(:,m,iactive_pool,istructural_below) = frac_soil(istructural,iactive,ibelow) * &
               dt*litter_turn(istructural) * &                     
               control_temp(:,ibelow) * control_moist(:,ibelow) * & 
               exp( -litter_struct_coef * lignin_struc(:,m,ibelow) ) * &
               ( 1. - lignin_struc(:,m,ibelow) ) 
          
          MatrixA(:,m,isurface_pool,imetabolic_above) =  frac_soil(imetabolic,isurface,iabove) * &
               dt*litter_turn(imetabolic) * control_temp(:,iabove) * control_moist(:,iabove) 
           
          MatrixA(:,m,iactive_pool,imetabolic_below) =  frac_soil(imetabolic,iactive,ibelow) * &
               dt*litter_turn(imetabolic) * control_temp(:,ibelow) * control_moist(:,ibelow)          
                    
          MatrixA(:,m,islow_pool,istructural_above) = frac_soil(istructural,islow,iabove) * &
               dt*litter_turn(istructural) * &                   
               control_temp(:,iabove) * control_moist(:,iabove) * &
               exp( -litter_struct_coef * lignin_struc(:,m,iabove) )* &
               lignin_struc(:,m,iabove)
          
          
          MatrixA(:,m,islow_pool,istructural_below) = frac_soil(istructural,islow,ibelow) * &
               dt*litter_turn(istructural) * &   
               control_temp(:,ibelow) * control_moist(:,ibelow) *  &
                  exp( -litter_struct_coef * lignin_struc(:,m,ibelow) )* &
                  lignin_struc(:,m,ibelow) 
          
          
          !- VectorB : carbon input -
          
          VectorB(:,m,istructural_above) = litter_inc(:,istructural,m,iabove,icarbon)
          VectorB(:,m,istructural_below) = litter_inc(:,istructural,m,ibelow,icarbon)
          VectorB(:,m,imetabolic_above)  = litter_inc(:,imetabolic,m,iabove,icarbon)
          VectorB(:,m,imetabolic_below)  = litter_inc(:,imetabolic,m,ibelow,icarbon)
          VectorB(:,m,iwoody_above)      = litter_inc(:,iwoody,m,iabove,icarbon)
          VectorB(:,m,iwoody_below)      = litter_inc(:,iwoody,m,ibelow,icarbon)

          IF (printlev>=4) WRITE(numout,*) 'We filled MatrixA and VectorB' 
         

          ! DSG: there WHERE statements are problematic ; we need to check the
          ! icarbon stocks rather the nutrient stocks as the nutrients stocks << carbon stocks
          ! and we need to have nutrient when we have carbon
          !DSG exchange initrogen 
          !WHERE(litter(:,istructural,m,iabove,initrogen) .GT. min_stomate)
          WHERE(litter(:,istructural,m,iabove,icarbon) .GT. min_stomate)
             CN_som_litter_longterm(:,m,istructural_above) = ( CN_som_litter_longterm(:,m,istructural_above) * (tau_CN_longterm-dt) &
                  + litter(:,istructural,m,iabove,icarbon)/litter(:,istructural,m,iabove,initrogen) * dt)/ (tau_CN_longterm)

             CP_som_litter_longterm(:,m,istructural_above) = ( CP_som_litter_longterm(:,m,istructural_above) * (tau_CN_longterm-dt) &
                  + litter(:,istructural,m,iabove,icarbon)/litter(:,istructural,m,iabove,iphosphorus) * dt)/ (tau_CN_longterm)
          ELSEWHERE
             CN_som_litter_longterm(:,m,istructural_above) = zero
             CP_som_litter_longterm(:,m,istructural_above) = zero
          ENDWHERE
          
          !WHERE(litter(:,istructural,m,ibelow,initrogen) .GT. min_stomate)
          WHERE(litter(:,istructural,m,ibelow,icarbon) .GT. min_stomate)
             CN_som_litter_longterm(:,m,istructural_below) = ( CN_som_litter_longterm(:,m,istructural_below) * (tau_CN_longterm-dt) &
                  + litter(:,istructural,m,ibelow,icarbon)/litter(:,istructural,m,ibelow,initrogen) * dt)/ (tau_CN_longterm)

             CP_som_litter_longterm(:,m,istructural_below) = ( CP_som_litter_longterm(:,m,istructural_below) * (tau_CN_longterm-dt) &
                  + litter(:,istructural,m,ibelow,icarbon)/litter(:,istructural,m,ibelow,iphosphorus) * dt)/ (tau_CN_longterm)
          ELSEWHERE
             CN_som_litter_longterm(:,m,istructural_below) = zero
             CP_som_litter_longterm(:,m,istructural_below) = zero
          ENDWHERE

          !WHERE(litter(:,imetabolic,m,iabove,initrogen) .GT. min_stomate)
          WHERE(litter(:,imetabolic,m,iabove,icarbon) .GT. min_stomate)
             CN_som_litter_longterm(:,m,imetabolic_above) = ( CN_som_litter_longterm(:,m,imetabolic_above) * (tau_CN_longterm-dt) &
                  + litter(:,imetabolic,m,iabove,icarbon)/litter(:,imetabolic,m,iabove,initrogen) * dt)/ (tau_CN_longterm)

             CP_som_litter_longterm(:,m,imetabolic_above) = ( CP_som_litter_longterm(:,m,imetabolic_above) * (tau_CN_longterm-dt) &
                  + litter(:,imetabolic,m,iabove,icarbon)/litter(:,imetabolic,m,iabove,iphosphorus) * dt)/ (tau_CN_longterm)
          ELSEWHERE
             CN_som_litter_longterm(:,m,imetabolic_above) = zero
             CP_som_litter_longterm(:,m,imetabolic_above) = zero
          ENDWHERE

          !WHERE(litter(:,imetabolic,m,ibelow,initrogen) .GT. min_stomate)
          WHERE(litter(:,imetabolic,m,ibelow,icarbon) .GT. min_stomate)
             CN_som_litter_longterm(:,m,imetabolic_below) = ( CN_som_litter_longterm(:,m,imetabolic_below) * (tau_CN_longterm-dt) &
                  + litter(:,imetabolic,m,ibelow,icarbon)/litter(:,imetabolic,m,ibelow,initrogen) * dt)/ (tau_CN_longterm)

             CP_som_litter_longterm(:,m,imetabolic_below) = ( CP_som_litter_longterm(:,m,imetabolic_below) * (tau_CN_longterm-dt) &
                  + litter(:,imetabolic,m,ibelow,icarbon)/litter(:,imetabolic,m,ibelow,iphosphorus) * dt)/ (tau_CN_longterm)
          ELSEWHERE
             CN_som_litter_longterm(:,m,imetabolic_below) = zero
             CP_som_litter_longterm(:,m,imetabolic_below) = zero
          ENDWHERE

          !WHERE(litter(:,iwoody,m,iabove,initrogen) .GT. min_stomate)
          WHERE(litter(:,iwoody,m,iabove,icarbon) .GT. min_stomate)
             CN_som_litter_longterm(:,m,iwoody_above) = ( CN_som_litter_longterm(:,m,iwoody_above) * (tau_CN_longterm-dt) &
                  + litter(:,iwoody,m,iabove,icarbon)/litter(:,iwoody,m,iabove,initrogen) * dt)/ (tau_CN_longterm)

             CP_som_litter_longterm(:,m,iwoody_above) = ( CP_som_litter_longterm(:,m,iwoody_above) * (tau_CN_longterm-dt) &
                  + litter(:,iwoody,m,iabove,icarbon)/litter(:,iwoody,m,iabove,iphosphorus) * dt)/ (tau_CN_longterm)
          ELSEWHERE
             CN_som_litter_longterm(:,m,iwoody_above) = zero
             CP_som_litter_longterm(:,m,iwoody_above) = zero
          ENDWHERE
          
          !WHERE(litter(:,iwoody,m,ibelow,initrogen) .GT. min_stomate)
          WHERE(litter(:,iwoody,m,ibelow,icarbon) .GT. min_stomate)
             CN_som_litter_longterm(:,m,iwoody_below) = ( CN_som_litter_longterm(:,m,iwoody_below) * (tau_CN_longterm-dt) &
                  + litter(:,iwoody,m,ibelow,icarbon)/litter(:,iwoody,m,ibelow,initrogen) * dt)/ (tau_CN_longterm)

             CP_som_litter_longterm(:,m,iwoody_below) = ( CP_som_litter_longterm(:,m,iwoody_below) * (tau_CN_longterm-dt) &
                  + litter(:,iwoody,m,ibelow,icarbon)/litter(:,iwoody,m,ibelow,iphosphorus) * dt)/ (tau_CN_longterm)
          ELSEWHERE
             CN_som_litter_longterm(:,m,iwoody_below) = zero
             CP_som_litter_longterm(:,m,iwoody_below) = zero
          ENDWHERE

       ENDDO ! Loop over # PFTs
 
          
    ENDIF ! spinup analytic 

    IF (printlev>=4) WRITE(numout,*) 'Leaving littercalc'

  END SUBROUTINE littercalc


!! ==============================================================================================================================\n
!! SUBROUTINE   : deadleaf
!!
!>\BRIEF        This routine calculates the deadleafcover. 
!!
!! DESCRIPTION  : It first calculates the lai corresponding to the dead leaves (LAI) using 
!! the dead leaves carbon content (DL) the specific leaf area (sla) and the 
!! maximal coverage fraction of a PFT (vegetmax) using the following equations:
!! \latexonly
!! \input{deadleaf1.tex}
!! \endlatexonly
!! \n
!! Then, the dead leaf cover (DLC) is calculated as following:\n
!! \latexonly
!! \input{deadleaf2.tex}
!! \endlatexonly
!! \n
!! 
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE: ::deadleaf_cover
!! 
!! REFERENCE(S) : None
!!
!! FLOWCHART : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE deadleaf (npts, veget_cov_max, dead_leaves, deadleaf_cover, &
                       sla_calc)

  !! 0. Variable and parameter declaration
    
    !! 0.1 Input variables

    INTEGER, INTENT(in)                          :: npts           !! Domain size - number of grid pixels (unitless)
    REAL, DIMENSION(npts,nvm,nlitt), INTENT(in)  :: dead_leaves    !! Dead leaves per ground unit area, per PFT, 
                                                                          !! metabolic and structural  
                                                                          !! @tex $(gC m^{-2})$ @endtex
    REAL,DIMENSION(npts,nvm),INTENT(in)          :: veget_cov_max  !! PFT "Maximal" coverage fraction of a PFT defined in 
                                                                          !! the input vegetation map 
                                                                          !! @tex $(m^2 m^{-2})$ @endtex 
    REAL,DIMENSION(npts,nvm),INTENT(in)          :: sla_calc       !! leaf age-related SLA
    
    !! 0.2 Output variables
    
    REAL, DIMENSION(npts), INTENT(out)           :: deadleaf_cover !! Fraction of soil covered by dead leaves over all PFTs
                                                                          !! (0-1, unitless)

    !! 0.3 Modified variables

    !! 0.4 Local variables

    REAL, DIMENSION(npts)                        :: dead_lai       !! LAI of dead leaves @tex $(m^2 m^{-2})$ @endtex
    INTEGER                                      :: j              !! Index (unitless)
!_ ================================================================================================================================
    
  !! 1. LAI of dead leaves
  
    dead_lai(:) = zero

    DO j = 1,nvm !Loop over PFTs
       dead_lai(:) = dead_lai(:) + ( dead_leaves(:,j,imetabolic) + &
                     dead_leaves(:,j,istructural) ) * sla_calc(:,j) &
            * veget_cov_max(:,j)
    ENDDO

  !! 2. fraction of soil covered by dead leaves

    deadleaf_cover(:) = un - exp( - 0.5 * dead_lai(:) )

    IF (printlev>=4) WRITE(numout,*) 'Leaving deadleaf'

  END SUBROUTINE deadleaf


!! ================================================================================================================================
!! FUNCTION     : control_moist_func
!!
!>\BRIEF        Calculate moisture control for litter and soil C decomposition
!!
!! DESCRIPTION  : Calculate moisture control factor applied
!! to litter decomposition and to soil carbon decomposition in
!! stomate_soilcarbon.f90 using the following equation: \n
!! \latexonly
!! \input{control_moist_func1.tex}
!! \endlatexonly
!! \n
!! with M the moisture control factor and soilmoisutre, the soil moisture 
!! calculated in sechiba.
!! Then, the function is ranged between Moistcont_min and 1:\n
!! \latexonly
!! \input{control_moist_func2.tex}
!! \endlatexonly
!! \n
!! RECENT CHANGE(S) : None
!!
!! RETURN VALUE : ::moistfunc_result
!! 
!! REFERENCE(S) : None
!!
!! FLOWCHART : None
!! \n
!_ ================================================================================================================================
  
  FUNCTION control_moist_func (npts, moist_in) RESULT (moistfunc_result)

  !! 0. Variable and parameter declaration
    
    !! 0.1 Input variables
          
    INTEGER, INTENT(in)               :: npts                !! Domain size - number of grid pixel (unitless)
    REAL, DIMENSION(npts), INTENT(in) :: moist_in            !! relative humidity (unitless)

    !! 0.2 Output variables
   
    REAL, DIMENSION(npts)             :: moistfunc_result    !! Moisture control factor (0.25-1, unitless)

    !! 0.3 Modified variables

    !! 0.4 Local variables

!_ ================================================================================================================================

    moistfunc_result(:) = -moist_coeff(1) * moist_in(:) * moist_in(:) + moist_coeff(2)* moist_in(:) - moist_coeff(3)
    moistfunc_result(:) = MAX( moistcont_min, MIN( un, moistfunc_result(:) ) )

  END FUNCTION control_moist_func


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

END MODULE stomate_litter
