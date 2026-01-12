
! =================================================================================================================================
! MODULE       : stomate_io
!
! CONTACT      : orchidee-help _at_ listes.ipsl.fr
!
! LICENCE      : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF       Module for read and write of restart files for all stomate modules.
!!
!!\n DESCRIPTION : This module contains the subroutines readstart and writerestart. All variables that will be read or written
!!                 are passed as argument to the subroutines. The subroutine readstart is called from stomate_initialize and 
!!                 writerestart is called from stomate_finalize.
!!                 Note: Not all variables saved in the start files are absolutely necessary. However, Sechiba's and Stomate's 
!!                 PFTs are not necessarily identical, and for that case this information needs to be saved.
!!
!!
!! RECENT CHANGE(S) : None
!!
!! REFERENCE(S)	: None
!!
!! SVN :
!! $HeadURL: svn://forge.ipsl.jussieu.fr/orchidee/branches/ORCHIDEE-CN-P/ORCHIDEE/src_stomate/stomate_io.f90 $
!! $Date: 2021-04-29 14:43:02 +0200 (四, 2021-04-29) $
!! $Revision: 7172 $
!! \n
!_ ================================================================================================================================
MODULE stomate_io
  USE stomate_data
  USE constantes
  USE constantes_soil
  USE mod_orchidee_para
  USE ioipsl_para 
  !-
  IMPLICIT NONE
  !-

  PRIVATE
  PUBLIC readstart, writerestart
  PUBLIC get_deposition_clear, readep
  PUBLIC get_lithology_clear, readlithology
  PUBLIC get_soil_orders
  !-
  ! first call?
  !-
  LOGICAL,SAVE :: firstcall_io = .TRUE.
!$OMP THREADPRIVATE(firstcall_io)
  LOGICAL, SAVE                                :: firstcall_lith = .TRUE. !! first call?
!$OMP THREADPRIVATE(firstcall_lith)
  LOGICAL, SAVE                                :: firstcall_sorder = .TRUE. !! first call?
!$OMP THREADPRIVATE(firstcall_sorder)
  LOGICAL, SAVE                                :: firstcall_dep = .TRUE. !! first call?
!$OMP THREADPRIVATE(firstcall_dep)
  !-
  ! reference temperature (K)
  !-
  REAL,ALLOCATABLE,DIMENSION(:),SAVE :: trefe
!$OMP THREADPRIVATE(trefe)
  REAL, ALLOCATABLE, DIMENSION(:,:), SAVE :: depo               !! inputs  of N&P by deposition [g/m2/day]
!$OMP THREADPRIVATE(depo)
  REAL, ALLOCATABLE, DIMENSION(:,:), SAVE :: lith               !! area fractions of lithologies [0-1]
!$OMP THREADPRIVATE(lith)
  REAL, ALLOCATABLE, DIMENSION(:)  , SAVE :: soilshield         !! factor to correct P weathering for soil shielding effects [0-1]
!$OMP THREADPRIVATE(soilshield)
 INTEGER, ALLOCATABLE, DIMENSION(:)  , SAVE :: soilorder          !! soil order 
!$OMP THREADPRIVATE(soilorder)

  !! Computational interpolation constants 
  INTEGER :: INTERPOL_AVG = 0                           !! Apply average fraction on interpolation
  INTEGER :: INTERPOL_BAR = 1                           !! Select biggest area found during interpolation,

  !-
CONTAINS

!! ================================================================================================================================
!! SUBROUTINE   : readstart
!!
!>\BRIEF        Read all variables for stomate from restart file. 
!!
!! DESCRIPTION  : Read all variables for stomate from restart file. 
!!                Initialize the variables if they were not found in the restart file or if there was no restart file.
!!                
!! \n
!_ ================================================================================================================================

  SUBROUTINE readstart &
       & (npts, index, lalo, resolution, t2m, dt_days, date_loc, &
       &  ind, adapted, regenerate, moiavail_daily, max_eau_var_daily, gdd_init_date, litterhum_daily, &
       &  t2m_daily, t2m_min_daily, tsurf_daily, tsoil_daily, &
       &  soilhum_daily, &
       &  drainage_longterm, runoff_longterm, &
       &  precip_daily, &
       &  gpp_daily, npp_daily, turnover_daily, &
       &  moiavail_month, moiavail_week, moiavail_growingseason, t2m_longterm, tau_longterm, &
       &  t2m_month, t2m_week, tsoil_month, soilhum_month, &
       &  fireindex, firelitter, &
       &  maxmoiavail_lastyear, maxmoiavail_thisyear, &
       &  minmoiavail_lastyear, minmoiavail_thisyear, &
       &  maxgppweek_lastyear, maxgppweek_thisyear, &
       &  gdd0_lastyear, gdd0_thisyear, precip_lastyear, precip_thisyear, &
       &  gdd_m5_dormance,  gdd_from_growthinit, gdd_midwinter, ncd_dormance, ngd_minus5, &
       &  PFTpresent, npp_longterm, lm_lastyearmax, lm_thisyearmax, &
       &  maxfpc_lastyear, maxfpc_thisyear, &
       &  turnover_longterm, gpp_week, npp_week, biomass, resp_maint_part, &
       &  leaf_age, leaf_frac, senescence, when_growthinit, age, &
       &  resp_hetero, resp_maint, resp_growth, co2_fire, co2_to_bm_dgvm, &
       &  n_to_bm, p_to_bm,veget_lastlight, everywhere, need_adjacent, RIP_time, &
       &  time_hum_min, hum_min_dormance, &
       &  litter, dead_leaves, &
       &  som, lignin_struc, lignin_wood, turnover_time, &
       &  prod10,prod100,flux10, flux100, &
       &  convflux, cflux_prod10, cflux_prod100, &
       &  prod10_harvest,prod100_harvest,flux10_harvest, flux100_harvest, &
       &  convflux_harvest, cflux_prod10_harvest, cflux_prod100_harvest, &
       &  woodharvestpft, bm_to_litter, carb_mass_total, &
       &  Tseason, Tseason_length, Tseason_tmp, & 
       &  Tmin_spring_time, begin_leaves, onset_date, &
       &  global_years, ok_equilibrium, nbp_accu, nbp_flux, &
       &  MatrixV, VectorU, previous_stock, current_stock, assim_param, & 
       &  CN_som_litter_longterm, tau_CN_longterm, KF, k_latosa_adapt, &
       &  CP_som_litter_longterm, np_leaf_avg_season, soil_p_min, & ! P variables
       &  f_Pdissolved, EW_grain,&
       &  rue_longterm, cn_leaf_avg_season, nstress_season, pstress_season, &
       &  lai_target_longterm, &
       &  soil_n_min, p_O2,bact, &
       &  record_reserve_grass, &
       &  wshtotsum, sr_ugb, sla_calc, nb_ani, grazed_frac, &
       &  import_yield, t2m_14, litter_not_avail, nb_grazingdays, &
       &  after_snow, after_wet, wet1day, wet2day, GRM_devstage, &
       &  days_senescence)
    !---------------------------------------------------------------------
    !- read start file
    !---------------------------------------------------------------------
    !-
    ! 0 declarations
    !-
    ! 0.1 input
    !-
    ! Domain size
    INTEGER,INTENT(in) :: npts
    ! Indices of the points on the map
    INTEGER,DIMENSION(npts),INTENT(in) :: index
    ! Geogr. coordinates (latitude,longitude) (degrees)
    REAL,DIMENSION(npts,2),INTENT(in) :: lalo
    ! size in x an y of the grid (m)
    REAL,DIMENSION(npts,2),INTENT(in) :: resolution
    REAL,DIMENSION(npts),INTENT(in)   :: t2m                !! 2 m air temperature from forcing file or coupled model (K)
    !-
    ! 0.2 output
    !-
    ! time step of STOMATE in days
    REAL,INTENT(out) :: dt_days
    ! date_loc (d)
    INTEGER,INTENT(out) :: date_loc
    ! density of individuals (1/m**2)
    REAL,DIMENSION(npts,nvm),INTENT(out) :: ind
    ! Winter too cold? between 0 and 1
    REAL,DIMENSION(npts,nvm),INTENT(out) :: adapted
    ! Winter sufficiently cold? between 0 and 1
    REAL,DIMENSION(npts,nvm),INTENT(out) :: regenerate
    ! daily moisture availability
    REAL,DIMENSION(npts,nvm),INTENT(out) :: moiavail_daily
    ! date for beginning of gdd count
    REAL,DIMENSION(npts,2),INTENT(out) :: gdd_init_date
    ! daily litter humidity
    REAL,DIMENSION(npts),INTENT(out)      :: litterhum_daily
    ! daily 2 meter temperatures (K)
    REAL,DIMENSION(npts),INTENT(out)      :: t2m_daily
    ! daily minimum 2 meter temperatures (K)
    REAL,DIMENSION(npts),INTENT(out)      :: t2m_min_daily
    ! daily surface temperatures (K)
    REAL,DIMENSION(npts),INTENT(out)      :: tsurf_daily
    ! daily soil temperatures (K)
    REAL,DIMENSION(npts,nslm),INTENT(out) :: tsoil_daily
    ! daily soil humidity
    REAL,DIMENSION(npts,nslm),INTENT(out) :: soilhum_daily
    ! daily precipitations (mm/day) (for phenology)
    ! daily 
    REAL,DIMENSION(npts),INTENT(out)      :: max_eau_var_daily       !!
    !
    REAL, DIMENSION(npts), INTENT(out)    :: drainage_longterm       !! "long term" annual sum of drainage 
                                                                            !! (kgH2O)/m**2/year)
    REAL, DIMENSION(npts), INTENT(out)    :: runoff_longterm         !! "long term" annual sum of runoff 
                                                                                                   !! (kgH2O)/m**2/year)
    REAL, DIMENSION(npts,nvm), INTENT(out)              :: lai_target_longterm    !! long-term target lai

    REAL,DIMENSION(npts),INTENT(out)      :: precip_daily
    ! daily gross primary productivity (gC/m**2/day)
    REAL,DIMENSION(npts,nvm),INTENT(out) :: gpp_daily
    ! daily net primary productivity (gC/m**2/day)
    REAL,DIMENSION(npts,nvm),INTENT(out) :: npp_daily
    ! daily turnover rates (gC/m**2/day)
    REAL,DIMENSION(npts,nvm,nparts,nelements),INTENT(out) :: turnover_daily
    ! "monthly" moisture availability
    REAL,DIMENSION(npts,nvm),INTENT(out) :: moiavail_month
    ! "weekly" moisture availability
    REAL,DIMENSION(npts,nvm),INTENT(out) :: moiavail_week
    ! mean growing season moisture availability (used for allocation response)
    REAL,DIMENSION(npts,nvm),INTENT(out) :: moiavail_growingseason
    ! "long term" 2 meter temperatures (K)
    REAL,DIMENSION(npts),INTENT(out)      :: t2m_longterm
    ! "tau_longterm"
    REAL, INTENT(out)        :: tau_longterm
    ! "monthly" 2 meter temperatures (K)
    REAL,DIMENSION(npts),INTENT(out)      :: t2m_month
    ! "seasonal" 2 meter temperatures (K) 
    REAL,DIMENSION(npts),INTENT(out)      :: Tseason
    ! temporary variable to calculate Tseason
    REAL,DIMENSION(npts),INTENT(out)      :: Tseason_length
    ! temporary variable to calculate Tseason
    REAL,DIMENSION(npts),INTENT(out)      :: Tseason_tmp
    REAL,DIMENSION(npts,nvm),INTENT(out)  :: Tmin_spring_time
    REAL,DIMENSION(npts,nvm),INTENT(out)  :: onset_date
    LOGICAL,DIMENSION(npts,nvm),INTENT(out)      :: begin_leaves

    ! "weekly" 2 meter temperatures (K)
    REAL,DIMENSION(npts),INTENT(out)      :: t2m_week
    ! "monthly" soil temperatures (K)
    REAL,DIMENSION(npts,nslm),INTENT(out) :: tsoil_month
    ! "monthly" soil humidity
    REAL,DIMENSION(npts,nslm),INTENT(out) :: soilhum_month
    ! Probability of fire
    REAL,DIMENSION(npts,nvm),INTENT(out) :: fireindex
    ! Longer term total litter above the ground, gC/m**2 of ground
    REAL,DIMENSION(npts,nvm),INTENT(out) :: firelitter
    ! last year's maximum moisture availability
    REAL,DIMENSION(npts,nvm),INTENT(out) :: maxmoiavail_lastyear
    ! this year's maximum moisture availability
    REAL,DIMENSION(npts,nvm),INTENT(out) :: maxmoiavail_thisyear
    ! last year's minimum moisture availability
    REAL,DIMENSION(npts,nvm),INTENT(out) :: minmoiavail_lastyear
    ! this year's minimum moisture availability
    REAL,DIMENSION(npts,nvm),INTENT(out) :: minmoiavail_thisyear
    ! last year's maximum weekly GPP
    REAL,DIMENSION(npts,nvm),INTENT(out) :: maxgppweek_lastyear
    ! this year's maximum weekly GPP
    REAL,DIMENSION(npts,nvm),INTENT(out) :: maxgppweek_thisyear
    ! last year's annual GDD0
    REAL,DIMENSION(npts),INTENT(out)      :: gdd0_lastyear
    ! this year's annual GDD0
    REAL,DIMENSION(npts),INTENT(out)      :: gdd0_thisyear
    ! last year's annual precipitation (mm/year)
    REAL,DIMENSION(npts),INTENT(out)      :: precip_lastyear
    ! this year's annual precipitation (mm/year)
    REAL,DIMENSION(npts),INTENT(out)      :: precip_thisyear
    ! growing degree days, threshold -5 deg C (for phenology)
    REAL,DIMENSION(npts,nvm),INTENT(out) :: gdd_m5_dormance
    ! growing degree days, from begin of season
    REAL,DIMENSION(npts,nvm),INTENT(out) :: gdd_from_growthinit
    ! growing degree days since midwinter (for phenology)
    REAL,DIMENSION(npts,nvm),INTENT(out) :: gdd_midwinter
    ! number of chilling days since leaves were lost (for phenology)
    REAL,DIMENSION(npts,nvm),INTENT(out) :: ncd_dormance
    ! number of growing days, threshold -5 deg C (for phenology)
    REAL,DIMENSION(npts,nvm),INTENT(out) :: ngd_minus5
    ! PFT exists (equivalent to fpc_max > 0 for natural PFTs)
    LOGICAL,DIMENSION(npts,nvm),INTENT(out)    :: PFTpresent
    ! "long term" net primary productivity (gC/m**2/year)
    REAL,DIMENSION(npts,nvm),INTENT(out) :: npp_longterm
    ! last year's maximum leaf mass, for each PFT (gC/m**2)
    REAL,DIMENSION(npts,nvm),INTENT(out) :: lm_lastyearmax
    ! this year's maximum leaf mass, for each PFT (gC/m**2)
    REAL,DIMENSION(npts,nvm),INTENT(out) :: lm_thisyearmax
    ! last year's maximum fpc for each natural PFT, on ground
    REAL,DIMENSION(npts,nvm),INTENT(out) :: maxfpc_lastyear
    ! this year's maximum fpc for each PFT,
    ! on *total* ground (see stomate_season)
    REAL,DIMENSION(npts,nvm),INTENT(out) :: maxfpc_thisyear
    ! "long term" turnover rate (gC/m**2/year)
    REAL,DIMENSION(npts,nvm,nparts,nelements),INTENT(out) :: turnover_longterm
    ! "weekly" GPP (gC/day/(m**2 covered)
    REAL,DIMENSION(npts,nvm),INTENT(out) :: gpp_week
    ! biomass (gC/m**2)
    REAL,DIMENSION(npts,nvm),INTENT(out) :: npp_week
    ! biomass (gC/m**2)
    REAL,DIMENSION(npts,nvm,nparts,nelements),INTENT(out) :: biomass
    ! maintenance resp (gC/m**2)
    REAL,DIMENSION(npts,nvm,nparts),INTENT(out) :: resp_maint_part
    ! leaf age (days)
    REAL,DIMENSION(npts,nvm,nparts,nleafages),INTENT(out) :: leaf_age
    ! fraction of leaves in leaf age class
    REAL,DIMENSION(npts,nvm,nparts,nleafages),INTENT(out) :: leaf_frac
    ! is the plant senescent ? 
    !(only for deciduous trees - carbohydrate reserve)
    LOGICAL,DIMENSION(npts,nvm),INTENT(out) :: senescence
    ! how many days ago was the beginning of the growing season
    REAL,DIMENSION(npts,nvm),INTENT(out) :: when_growthinit
    ! mean age (years)
    REAL,DIMENSION(npts,nvm),INTENT(out) :: age
    ! heterotrophic respiration (gC/day/m**2)
    REAL,DIMENSION(npts,nvm),INTENT(out) :: resp_hetero
    ! maintenance respiration (gC/day/m**2)
    REAL,DIMENSION(npts,nvm),INTENT(out) :: resp_maint
    ! growth respiration (gC/day/m**2)
    REAL,DIMENSION(npts,nvm),INTENT(out) :: resp_growth
    ! carbon emitted into the atmosphere by fire (living and dead biomass)
    ! (in gC/m**2/time step)
    REAL,DIMENSION(npts,nvm),INTENT(out) :: co2_fire
    ! biomass uptaken (gC/(m**2 of total ground)/day)
    REAL,DIMENSION(npts,nvm),INTENT(out) :: co2_to_bm_dgvm
    ! biomass uptaken (gN/(m**2 of total ground)/day)
    REAL,DIMENSION(npts,nvm),INTENT(out) :: n_to_bm
    REAL,DIMENSION(npts,nvm),INTENT(out) :: p_to_bm
    ! vegetation fractions (on ground) after last light competition
    REAL,DIMENSION(npts,nvm),INTENT(out) :: veget_lastlight
    ! is the PFT everywhere in the grid box or very localized
    ! (after its introduction)
    REAL,DIMENSION(npts,nvm),INTENT(out) :: everywhere
    ! in order for this PFT to be introduced,
    ! does it have to be present in an adjacent grid box?
    LOGICAL,DIMENSION(npts,nvm),INTENT(out) :: need_adjacent
    ! How much time ago was the PFT eliminated for the last time (y)
    REAL,DIMENSION(npts,nvm),INTENT(out) :: RIP_time
    ! time elapsed since strongest moisture availability (d)
    REAL,DIMENSION(npts,nvm),INTENT(out) :: time_hum_min
    ! minimum moisture during dormance
    REAL,DIMENSION(npts,nvm),INTENT(out) :: hum_min_dormance
    ! fraction of litter above the ground belonging to different PFTs
    ! separated for natural and agricultural PFTs.
    REAL,DIMENSION(npts,nlitt,nvm,nlevs,nelements),INTENT(out):: litter
    ! dead leaves on ground, per PFT, metabolic and structural,
    ! in gC/(m**2 of ground)
    REAL,DIMENSION(npts,nvm,nlitt),INTENT(out) :: dead_leaves
    ! Soil Organic Matter pool: active, slow, or passive, (gC (or N)/m**2)
    REAL,DIMENSION(npts,ncarb,nvm,nelements),INTENT(out) :: som
    ! ratio Lignine/Carbon in structural litter, above and below ground,(gC/m**2)
    REAL,DIMENSION(npts,nvm,nlevs),INTENT(out) :: lignin_struc
    ! ratio Lignine/Carbon in woody litter, above and below ground,(gC/m**2)
    REAL,DIMENSION(npts,nvm,nlevs),INTENT(out) :: lignin_wood
    REAL,DIMENSION(npts,nvm),INTENT(out) :: turnover_time

    ! For Spinup matrix resolution
    INTEGER, INTENT(out) :: global_years   
    LOGICAL, DIMENSION(npts), INTENT(out) :: ok_equilibrium
  !DSG  LOGICAL,  INTENT(out)                 :: ok_spunup !anal_spin
    REAL, DIMENSION(npts), INTENT(out) :: nbp_accu  !! Accumulated Net Biospheric Production over the year
    REAL, DIMENSION(npts), INTENT(out) :: nbp_flux  !! Net Biospheric Production over the year
    !-
    REAL, DIMENSION(npts,nvm,nbpools,nbpools), INTENT(out) :: MatrixV
    REAL, DIMENSION(npts,nvm,nbpools), INTENT(out) :: VectorU
    REAL, DIMENSION(npts,nvm,nbpools), INTENT(out) :: previous_stock
    REAL, DIMENSION(npts,nvm,nbpools), INTENT(out) :: current_stock
    REAL, DIMENSION(npts,nvm,nbpools), INTENT(out) :: CN_som_litter_longterm
    REAL, INTENT(out)                              :: tau_CN_longterm  !! Counter used for calculating the longterm CN ratio of SOM and litter pools (seconds)
    REAL, DIMENSION(npts,nvm,npco2),   INTENT(out) :: assim_param

    REAL, DIMENSION(npts,nvm), INTENT(out)           :: record_reserve_grass !! record maximum grass reserve pool

    REAL, DIMENSION(npts,nvm), INTENT(out)           :: KF             !! Scaling factor to convert sapwood mass                     
                                                                              !! into leaf mass (m)
    REAL, DIMENSION(npts,nvm), INTENT(out)           :: k_latosa_adapt !! Leaf to sapwood area adapted for water       
                                                                              !! stress. Adaptation takes place at the 
                                                                              !! end of the year (m)
    REAL, DIMENSION(npts,nvm), INTENT(out)                       :: rue_longterm            !! longterm radiation use efficiency
    REAL, DIMENSION(npts,nvm), INTENT(out)                       :: cn_leaf_avg_season      !! Seasonal average CN ratio of leaves 
    REAL, DIMENSION(npts,nvm), INTENT(out)                       :: nstress_season          !! N-related seasonal stress (used for allocation) 
    REAL, DIMENSION(npts,nvm), INTENT(out)                       :: pstress_season          !! P-related seasonal stress (used for allocation) 
   
    REAL, DIMENSION(npts,nvm,nnspec), INTENT(out)                :: soil_n_min              !! mineral nitrogen in the soil (gN/m**2)  
                                                                                                   !! (first index=npts, second index=nvm, third index=nnspec) 
    REAL, DIMENSION(npts,nvm), INTENT(out)                       :: p_O2                    !! partial pressure of oxigen in the soil (hPa)
                                                                                                   !! (first index=npts, second index=nvm)
                      
    REAL, DIMENSION(npts,nvm), INTENT(out)                       :: bact                    !! denitrifier biomass (gC/m**2)
                                                                                                   !! (first index=npts, second index=nvm)
    REAL, DIMENSION(npts,nvm,nbpools), INTENT(out)               :: CP_som_litter_longterm
    REAL, DIMENSION(npts,nvm), INTENT(out)                       :: np_leaf_avg_season      !! Seasonal average CP ratio of leaves 
    REAL, DIMENSION(npts,nvm,npspec), INTENT(out)                :: soil_p_min              !! mineral phosphorus in the soil (gP/m**2)  
                                                                                                   !! (first index=npts, second index=nvm, third index=npspec) 
    REAL, DIMENSION(npts,nvm), INTENT(out)                       :: f_Pdissolved            !! fraction of dissolved P in contact with roots [ ]
                                                                                                   !! (first index=npts, second index=nvm ) 
    REAL, DIMENSION(npts,nvm,nminerals,2), INTENT(out)           :: EW_grain               !! Mass and grains size of crushed rocks 

    REAL, DIMENSION(npts,nvm),INTENT(out)   :: sla_calc      !! leaf age-related SLA
    REAL, DIMENSION(npts,nvm),INTENT(out)   :: wshtotsum     !! accumulated harvested grass biomass
    REAL, DIMENSION(npts,nvm),INTENT(out)   :: sr_ugb        !! grazing stocking rate of grazing livestock
    REAL, DIMENSION(npts,nvm),INTENT(out)   :: nb_ani        !! optimal livestock density
    REAL, DIMENSION(npts,nvm),INTENT(out)   :: grazed_frac   !! optimal fraction for grazing (in contract to mowning)
    REAL, DIMENSION(npts,nvm), INTENT(out)  :: import_yield  !! annual total harvested grass biomass yield
    REAL, DIMENSION(npts),INTENT(out)       :: t2m_14        !! ''14 days'' 2 meter air temperature
    REAL, DIMENSION(npts,nlitt,nvm), INTENT(out)  :: litter_not_avail !! litter not edible for animals
    REAL, DIMENSION(npts,nvm), INTENT(out)  :: nb_grazingdays!!annual total number of days animal grazing in field 
    REAL, DIMENSION(npts),INTENT(out)       :: after_snow    !! day counter after snow melt to prevent grazing
    REAL, DIMENSION(npts),INTENT(out)       :: after_wet     !! day counter after wet soil is detected to prevent grazing
    REAL, DIMENSION(npts),INTENT(out)       :: wet1day       !! accumulated days with wet soil
    REAL, DIMENSION(npts),INTENT(out)       :: wet2day       !! accumulated days with wet soil
    REAL, DIMENSION(npts,nvm), INTENT(out)  :: GRM_devstage  !! development stage of grasses
    REAL, DIMENSION(npts,nvm), INTENT(out)  :: days_senescence!! day counter after senescence is detected
 
    ! 0.4 local
    !-
    ! date, real
    REAL :: date_real
    ! PFT exists (equivalent to fpc_max > 0 for natural PFTs), real
    REAL,DIMENSION(npts,nvm) :: PFTpresent_real
    ! is the plant senescent ?
    ! (only for deciduous trees - carbohydrate reserve), real
    REAL,DIMENSION(npts,nvm) :: senescence_real
    REAL,DIMENSION(npts,nvm) :: begin_leaves_real
    ! in order for this PFT to be introduced,
    ! does it have to be present in an adjacent grid box? - real
    REAL,DIMENSION(npts,nvm) :: need_adjacent_real
    REAL, DIMENSION(1) :: vartmp  !! temporary variable because restget/restput needs an array and not a scalar
    ! To store variables names for I/O
    CHARACTER(LEN=80) :: var_name
    ! string suffix indicating an index
    CHARACTER(LEN=10) :: part_str
    ! string suffix indicating litter type
    CHARACTER(LEN=4),DIMENSION(nlitt) :: litter_str
    ! string suffix indicating level
    CHARACTER(LEN=2),DIMENSION(nlevs) :: level_str
    CHARACTER(LEN=4),DIMENSION(nleafages) :: age_str
    ! temporary storage
    REAL,DIMENSION(1) :: xtmp
    ! index
    INTEGER :: j,k,l,m
    ! reference temperature (K)

    CHARACTER(LEN=2),DIMENSION(nelements) :: element_str   !! string suffix indicating element
    REAL, DIMENSION(1) :: temp_global_years
    CHARACTER(LEN=6), DIMENSION(nbpools) :: pools_str
    REAL, DIMENSION(npts) :: ok_equilibrium_real    
    ! land cover change variables 
    ! products remaining in the 10/100 year-turnover pool after the annual release for each compartment
    ! (10 or 100 + 1 : input from year of land cover change)
    REAL,DIMENSION(npts,0:10),INTENT(out)                         :: prod10
    REAL,DIMENSION(npts,0:100),INTENT(out)                        :: prod100
    ! annual release from the 10/100 year-turnover pool compartments
    REAL,DIMENSION(npts,10),INTENT(out)                           :: flux10
    REAL,DIMENSION(npts,100),INTENT(out)                          :: flux100
    REAL, DIMENSION(npts), INTENT(out)                            :: convflux
    REAL, DIMENSION(npts), INTENT(out)                            :: cflux_prod10
    REAL, DIMENSION(npts), INTENT(out)                            :: cflux_prod100
    ! wood harvest variables 
    ! products remaining in the 10/100 year-turnover pool after the annual release for each compartment
    ! (10 or 100 + 1 : input from year of land cover change)
    REAL,DIMENSION(npts,0:10),INTENT(out)                           :: prod10_harvest
    REAL,DIMENSION(npts,0:100),INTENT(out)                          :: prod100_harvest
    ! annual release from the 10/100 year-turnover pool compartments
    REAL,DIMENSION(npts,10),INTENT(out)                           :: flux10_harvest
    REAL,DIMENSION(npts,100),INTENT(out)                          :: flux100_harvest
    REAL, DIMENSION(npts), INTENT(out)                            :: convflux_harvest
    REAL, DIMENSION(npts), INTENT(out)                            :: cflux_prod10_harvest
    REAL, DIMENSION(npts), INTENT(out)                            :: cflux_prod100_harvest
    REAL, DIMENSION(npts,nvm), INTENT(out)                        :: woodharvestpft
    REAL,DIMENSION(npts,nvm,nparts,nelements),INTENT(out)         :: bm_to_litter
    REAL,DIMENSION(npts),INTENT(out)                              :: carb_mass_total
    REAL,DIMENSION(npts,nvm)                                      :: vcmax_tmp
    !---------------------------------------------------------------------
    IF (printlev >= 3) WRITE(numout,*) 'Entering readstart'
    !-
    ! 1 string definitions
    !-
    DO l=1,nlitt
       IF     (l == imetabolic) THEN
          litter_str(l) = 'met'
       ELSEIF (l == istructural) THEN
          litter_str(l) = 'str'
       ELSEIF (l == iwoody) THEN
          litter_str(l) = 'wood'
       ELSE
          CALL ipslerr_p(3,'stomate_io readstart', 'Define litter_str','','')
       ENDIF
    ENDDO
    !-
    DO l=1,nlevs
       IF     (l == iabove) THEN
          level_str(l) = 'ab'
       ELSEIF (l == ibelow) THEN
          level_str(l) = 'be'
       ELSE
          CALL ipslerr_p(3,'stomate_io readstart','Define level_str','','')
       ENDIF
    ENDDO

    !JC separate tissue age
    !-
    DO l=1,nleafages
       IF     (l == 1) THEN
          age_str(l) = 'age1'
       ELSEIF (l == 2) THEN
          age_str(l) = 'age2'
       ELSEIF (l == 3) THEN
          age_str(l) = 'age3'
       ELSEIF (l == 4) THEN
          age_str(l) = 'age4'
       ELSE
          CALL ipslerr_p(3,'stomate_io readstart','Define age_str','','')
       ENDIF
    ENDDO

    pools_str(1:nbpools) =(/'str_ab ','str_be ','met_ab ','met_be ','wood_ab','wood_be',& 
         & 'actif  ','slow   ','passif ','surface'/) 
    !-
    DO l=1,nelements
       IF     (l == icarbon) THEN
          element_str(l) = ''
       ELSEIF (l == initrogen) THEN
          element_str(l) = '_n'
       ELSEIF (l == iphosphorus) THEN
          element_str(l) = '_p'
       ELSE
          CALL ipslerr_p(3,'stomate_io readstart','Define element_str','','')
       ENDIF
    ENDDO
    !-
    ! 2 run control
    !-
    ! 2.2 time step of STOMATE in days
    !-
    IF (is_root_prc) THEN
       var_name = 'dt_days'
       CALL restget (rest_id_stomate, var_name, 1   , 1     , 1, itime, &
            &                 .TRUE., xtmp)
       dt_days = xtmp(1)
       IF (dt_days == val_exp) dt_days = un
    ENDIF
    CALL bcast(dt_days)
    !-
    ! 2.3 date
    !-
    IF (is_root_prc) THEN
       var_name = 'date'
       CALL restget (rest_id_stomate, var_name, 1   , 1     , 1, itime, &
            &                 .TRUE., xtmp)
       date_real = xtmp(1)
       IF (date_real == val_exp) date_real = zero
       date_loc = NINT(date_real)
    ENDIF
    CALL bcast(date_loc)
    !-
    ! 3 daily meteorological variables
    !-
    moiavail_daily(:,:) = val_exp
    var_name = 'moiavail_daily'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &              .TRUE., moiavail_daily, 'gather', nbp_glo, index_g)
    IF (ALL(moiavail_daily(:,:) == val_exp)) moiavail_daily(:,:) = zero
    !-
    max_eau_var_daily(:) = val_exp
    var_name = 'max_eau_var_daily'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo,   1  , 1, itime, &
         &              .TRUE., max_eau_var_daily, 'gather', nbp_glo, index_g)
    IF (ALL(max_eau_var_daily(:) == val_exp)) max_eau_var_daily(:) = zero
    !-
    gdd_init_date(:,:) = val_exp
    var_name = 'gdd_init_date'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, 2 , 1, itime, &
         &              .TRUE., gdd_init_date, 'gather', nbp_glo, index_g)
    ! Keep val_exp as initial value for gdd_init_date(:,2)
    IF (ALL(gdd_init_date(:,1) == val_exp)) gdd_init_date(:,1) = 365.

    !-
    litterhum_daily(:) = val_exp
    var_name = 'litterhum_daily'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, 1     , 1, itime, &
         &              .TRUE., litterhum_daily, 'gather', nbp_glo, index_g)
    IF (ALL(litterhum_daily(:) == val_exp)) litterhum_daily(:) = zero
    !-
    t2m_daily(:) = val_exp
    var_name = 't2m_daily'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, 1     , 1, itime, &
         &                .TRUE., t2m_daily, 'gather', nbp_glo, index_g)
    IF (ALL(t2m_daily(:) == val_exp)) t2m_daily(:) = zero
    !-
    t2m_min_daily(:) = val_exp
    var_name = 't2m_min_daily'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, 1     , 1, itime, &
         &                .TRUE., t2m_min_daily, 'gather', nbp_glo, index_g)
    IF (ALL(t2m_min_daily(:) == val_exp)) t2m_min_daily(:) = large_value
    !-
    tsurf_daily(:) = val_exp
    var_name = 'tsurf_daily'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, 1     , 1, itime, &
         &                .TRUE., tsurf_daily, 'gather', nbp_glo, index_g)
    ! The initial value is set to the current temperature at 2m
    IF (ALL(tsurf_daily(:) == val_exp)) tsurf_daily(:) = t2m(:)
    !-
    tsoil_daily(:,:) = val_exp
    var_name = 'tsoil_daily'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo,   nslm, 1, itime, &
         &                .TRUE., tsoil_daily, 'gather', nbp_glo, index_g)
    IF (ALL(tsoil_daily(:,:) == val_exp)) tsoil_daily(:,:) = zero
    !-
    soilhum_daily(:,:) = val_exp
    var_name = 'soilhum_daily'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo,   nslm, 1, itime, &
         &                .TRUE., soilhum_daily, 'gather', nbp_glo, index_g)
    IF (ALL(soilhum_daily(:,:) == val_exp)) soilhum_daily(:,:) = zero

    !-
    drainage_longterm(:) = val_exp
    var_name = 'drainage_longterm'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo,   1  , 1, itime, &
         &              .TRUE., drainage_longterm, 'gather', nbp_glo, index_g)
    IF (ALL(drainage_longterm(:) == val_exp)) drainage_longterm(:) = zero
    !-
    runoff_longterm(:) = val_exp
    var_name = 'runoff_longterm'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo,   1  , 1, itime, &
         &              .TRUE., runoff_longterm, 'gather', nbp_glo, index_g)
    IF (ALL(runoff_longterm(:) == val_exp)) runoff_longterm(:) = zero
    !-
    lai_target_longterm(:,:) = val_exp
    var_name = 'lai_target_longterm'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo,   nvm, 1, itime, &
         &              .TRUE., lai_target_longterm, 'gather', nbp_glo, index_g)
    IF (ALL(lai_target_longterm(:,:) == val_exp)) lai_target_longterm(:,:) = zero
    !-
    precip_daily(:) = val_exp
    var_name = 'precip_daily'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, 1     , 1, itime, &
         &                .TRUE., precip_daily, 'gather', nbp_glo, index_g)
    IF (ALL(precip_daily(:) == val_exp)) precip_daily(:) = zero
    !-
    ! 4 productivities
    !-
    gpp_daily(:,:) = val_exp
    var_name = 'gpp_daily'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &              .TRUE., gpp_daily, 'gather', nbp_glo, index_g)
    IF (ALL(gpp_daily(:,:) == val_exp)) gpp_daily(:,:) = zero
    !-
    npp_daily(:,:) = val_exp
    var_name = 'npp_daily'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &              .TRUE., npp_daily, 'gather', nbp_glo, index_g)
    IF (ALL(npp_daily(:,:) == val_exp)) npp_daily(:,:) = zero
    !-
    turnover_daily(:,:,:,:) = val_exp
    DO l = 1,nelements
       DO k = 1,nparts
          WRITE(part_str,'(I2)') k
          IF (k < 10) part_str(1:1) = '0'
          var_name = 'turnover_daily_'//part_str(1:LEN_TRIM(part_str))//element_str(l)
          CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
               &                .TRUE., turnover_daily(:,:,k,l), 'gather', nbp_glo, index_g)
          IF (ALL(turnover_daily(:,:,k,l) == val_exp)) &
               &       turnover_daily(:,:,k,l) = zero
       ENDDO
    END DO
    !-
    ! 5 monthly meteorological variables
    !-
    moiavail_month(:,:) = val_exp
    var_name = 'moiavail_month'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &              .TRUE., moiavail_month, 'gather', nbp_glo, index_g)
    IF (ALL(moiavail_month(:,:) == val_exp)) moiavail_month(:,:) = zero
    !-
    moiavail_week(:,:) = val_exp
    var_name = 'moiavail_week'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &              .TRUE., moiavail_week, 'gather', nbp_glo, index_g)
    IF (ALL(moiavail_week(:,:) == val_exp)) moiavail_week(:,:) = zero

    moiavail_growingseason(:,:) = val_exp
    var_name = 'moiavail_grow'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &              .TRUE., moiavail_growingseason, 'gather', nbp_glo, index_g)
    IF (ALL(moiavail_growingseason(:,:) == val_exp)) moiavail_growingseason(:,:) = un
    

    !
    ! Longterm temperature at 2m
    !
    var_name = 't2m_longterm'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, 1     , 1, itime, &
         &              .TRUE., t2m_longterm, 'gather', nbp_glo, index_g)

    IF (ALL(t2m_longterm(:) == val_exp)) THEN
       ! t2m_longterm is not in restart file
       ! The initial value for the reference temperature is set to the current temperature
       t2m_longterm(:)=t2m(:)
       ! Set the counter to 2 time steps
       tau_longterm=2
    ELSE
       ! t2m_longterm was in the restart file
       ! Now read tau_longterm
       ! tau_longterm is a scalar, therefor only master process read this value
       IF (is_root_prc) THEN
          CALL restget (rest_id_stomate, 'tau_longterm', 1 ,1  , 1, itime, &
               .TRUE., vartmp)
          IF (vartmp(1) == val_exp) THEN
             ! tau_longterm is not found in restart file. 
             ! This is not normal as t2m_longterm was in restart file. Write a warning and initialize it to tau_longterm_max
             CALL ipslerr(2, 'stomate_io readstart','tau_longterm was not in restart file',&
                  'But t2m_longterm was in restart file','')
             tau_longterm = tau_longterm_max
          ELSE
             tau_longterm = vartmp(1)
          END IF
       ENDIF
       CALL bcast(tau_longterm)

    END IF
    !-
    t2m_month(:) = val_exp
    var_name = 't2m_month'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, 1     , 1, itime, &
         &              .TRUE., t2m_month, 'gather', nbp_glo, index_g)
    IF (ALL(t2m_month(:) == val_exp)) t2m_month(:) = t2m(:)
    
    CALL restget_p (rest_id_stomate, 'Tseason', nbp_glo, 1     , 1, itime, &
         .TRUE., Tseason, 'gather', nbp_glo, index_g)
    IF (ALL(Tseason(:) == val_exp)) Tseason(:) = t2m(:)
    
    CALL restget_p (rest_id_stomate,'Tseason_length', nbp_glo, 1     , 1, itime, &
         .TRUE., Tseason_length, 'gather', nbp_glo, index_g)
    IF (ALL(Tseason_length(:) == val_exp)) Tseason_length(:) = zero
    
    CALL restget_p (rest_id_stomate, 'Tseason_tmp', nbp_glo, 1     , 1, itime, &
         .TRUE., Tseason_tmp, 'gather', nbp_glo, index_g)
    IF (ALL(Tseason_tmp(:) == val_exp)) Tseason_tmp(:) = zero

    CALL restget_p (rest_id_stomate, 'Tmin_spring_time', nbp_glo, nvm, 1, itime, &
         .TRUE., Tmin_spring_time, 'gather', nbp_glo, index_g)
    IF (ALL(Tmin_spring_time(:,:) == val_exp)) Tmin_spring_time(:,:) = zero
    
    CALL restget_p (rest_id_stomate, 'onset_date', nbp_glo, nvm  , 1, itime, &
         .TRUE., onset_date(:,:), 'gather', nbp_glo, index_g)
    IF (ALL(onset_date(:,:) == val_exp)) onset_date(:,:) = zero

    t2m_week(:) = val_exp
    var_name = 't2m_week'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, 1     , 1, itime, &
         &              .TRUE., t2m_week, 'gather', nbp_glo, index_g)
    ! The initial value is set to the current temperature
    IF (ALL(t2m_week(:) == val_exp)) t2m_week(:) = t2m(:)
    
    tsoil_month(:,:) = val_exp
    var_name = 'tsoil_month'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo,   nslm, 1, itime, &
         &              .TRUE., tsoil_month, 'gather', nbp_glo, index_g)

    ! The initial value is set to the current temperature
    IF (ALL(tsoil_month(:,:) == val_exp)) THEN
       DO l=1,nslm
          tsoil_month(:,l) = t2m(:)
       ENDDO
    ENDIF
    !-
    soilhum_month(:,:) = val_exp
    var_name = 'soilhum_month'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo,   nslm, 1, itime, &
         &              .TRUE., soilhum_month, 'gather', nbp_glo, index_g)
    IF (ALL(soilhum_month(:,:) == val_exp)) soilhum_month(:,:) = zero
    !-
    ! 6 fire probability
    !-
    fireindex(:,:) = val_exp
    var_name = 'fireindex'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &              .TRUE., fireindex, 'gather', nbp_glo, index_g)
    IF (ALL(fireindex(:,:) == val_exp)) fireindex(:,:) = zero
    !-
    firelitter(:,:) = val_exp
    var_name = 'firelitter'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &              .TRUE., firelitter, 'gather', nbp_glo, index_g)
    IF (ALL(firelitter(:,:) == val_exp)) firelitter(:,:) = zero
    !-
    ! 7 maximum and minimum moisture availabilities for tropic phenology
    !-
    maxmoiavail_lastyear(:,:) = val_exp
    var_name = 'maxmoistr_last'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &              .TRUE., maxmoiavail_lastyear, 'gather', nbp_glo, index_g)
    IF (ALL(maxmoiavail_lastyear(:,:) == val_exp)) &
         &     maxmoiavail_lastyear(:,:) = zero
    !-
    maxmoiavail_thisyear(:,:) = val_exp
    var_name = 'maxmoistr_this'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &              .TRUE., maxmoiavail_thisyear, 'gather', nbp_glo, index_g)
    IF (ALL(maxmoiavail_thisyear(:,:) == val_exp)) &
         &     maxmoiavail_thisyear(:,:) = zero
    !-
    minmoiavail_lastyear(:,:) = val_exp
    var_name = 'minmoistr_last'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &              .TRUE., minmoiavail_lastyear, 'gather', nbp_glo, index_g)
    IF (ALL(minmoiavail_lastyear(:,:) == val_exp)) &
         &     minmoiavail_lastyear(:,:) = un
    !-
    minmoiavail_thisyear(:,:) = val_exp
    var_name = 'minmoistr_this'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &              .TRUE., minmoiavail_thisyear, 'gather', nbp_glo, index_g)
    IF (ALL( minmoiavail_thisyear(:,:) == val_exp)) &
         &     minmoiavail_thisyear(:,:) = un
    !-
    ! 8 maximum "weekly" GPP
    !-
    maxgppweek_lastyear(:,:) = val_exp
    var_name = 'maxgppweek_lastyear'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &              .TRUE., maxgppweek_lastyear, 'gather', nbp_glo, index_g)
    IF (ALL(maxgppweek_lastyear(:,:) == val_exp)) &
         &     maxgppweek_lastyear(:,:) = zero
    !-
    maxgppweek_thisyear(:,:) = val_exp
    var_name = 'maxgppweek_thisyear'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &              .TRUE., maxgppweek_thisyear, 'gather', nbp_glo, index_g)
    IF (ALL(maxgppweek_thisyear(:,:) == val_exp)) &
         &     maxgppweek_thisyear(:,:) = zero
    !-
    ! 9 annual GDD0
    !-
    gdd0_thisyear(:) = val_exp
    var_name = 'gdd0_thisyear'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, 1     , 1, itime, &
         &              .TRUE., gdd0_thisyear, 'gather', nbp_glo, index_g)
    IF (ALL(gdd0_thisyear(:) == val_exp)) gdd0_thisyear(:) = zero
    !-
    gdd0_lastyear(:) = val_exp
    var_name = 'gdd0_lastyear'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, 1     , 1, itime, &
         &              .TRUE., gdd0_lastyear, 'gather', nbp_glo, index_g)
    IF (ALL(gdd0_lastyear(:) == val_exp)) gdd0_lastyear(:) = gdd_crit_estab
    !-
    ! 10 annual precipitation
    !-
    precip_thisyear(:) = val_exp
    var_name = 'precip_thisyear'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, 1     , 1, itime, &
         &              .TRUE., precip_thisyear, 'gather', nbp_glo, index_g)
    IF (ALL(precip_thisyear(:) == val_exp)) precip_thisyear(:) = zero
    !-
    precip_lastyear(:) = val_exp
    var_name = 'precip_lastyear'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, 1     , 1, itime, &
         &              .TRUE., precip_lastyear, 'gather', nbp_glo, index_g)
    IF (ALL(precip_lastyear(:) == val_exp)) &
         &     precip_lastyear(:) = precip_crit
    !-
    ! 11 derived "biometeorological" variables
    !-
    gdd_m5_dormance(:,:) = val_exp
    var_name = 'gdd_m5_dormance'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &              .TRUE., gdd_m5_dormance, 'gather', nbp_glo, index_g)
    IF (ALL(gdd_m5_dormance(:,:) == val_exp)) &
         &     gdd_m5_dormance(:,:) = undef
    !-
    gdd_from_growthinit(:,:) = val_exp
    var_name = 'gdd_from_growthinit'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &              .TRUE., gdd_from_growthinit, 'gather', nbp_glo, index_g)
    IF (ALL(gdd_from_growthinit(:,:) == val_exp)) &
         &     gdd_from_growthinit(:,:) = zero
    !-
    gdd_midwinter(:,:) = val_exp
    var_name = 'gdd_midwinter'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &              .TRUE., gdd_midwinter, 'gather', nbp_glo, index_g)
    IF (ALL(gdd_midwinter(:,:) == val_exp)) gdd_midwinter(:,:) = undef
    !-
    ncd_dormance(:,:) = val_exp
    var_name = 'ncd_dormance'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &              .TRUE., ncd_dormance, 'gather', nbp_glo, index_g)
    IF (ALL(ncd_dormance(:,:) == val_exp)) ncd_dormance(:,:) = undef
    !-
    ngd_minus5(:,:) = val_exp
    var_name = 'ngd_minus5'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &              .TRUE., ngd_minus5, 'gather', nbp_glo, index_g)
    IF (ALL(ngd_minus5(:,:) == val_exp)) ngd_minus5(:,:) = zero
    !-
    time_hum_min(:,:) = val_exp
    var_name = 'time_hum_min'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &              .TRUE., time_hum_min, 'gather', nbp_glo, index_g)
    IF (ALL(time_hum_min(:,:) == val_exp)) time_hum_min(:,:) = undef
    !-
    hum_min_dormance(:,:) = val_exp
    var_name = 'hum_min_dormance'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &              .TRUE., hum_min_dormance, 'gather', nbp_glo, index_g)
    IF (ALL(hum_min_dormance(:,:) == val_exp)) &
         &     hum_min_dormance(:,:) = undef
    !-
    ! 12 Plant status
    !-
    PFTpresent_real(:,:) = val_exp
    var_name = 'PFTpresent'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &              .TRUE., PFTpresent_real, 'gather', nbp_glo, index_g)
    IF (ALL(PFTpresent_real(:,:) == val_exp)) PFTpresent_real(:,:) = zero
    WHERE (PFTpresent_real(:,:) >= .5)
       PFTpresent = .TRUE.
    ELSEWHERE
       PFTpresent = .FALSE.
    ENDWHERE
    !-
    ind(:,:) = val_exp
    var_name = 'ind'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &              .TRUE., ind, 'gather', nbp_glo, index_g)
    IF (ALL(ind(:,:) == val_exp)) ind(:,:) = zero
    !-
    adapted(:,:) = val_exp
    var_name = 'adapted'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &              .TRUE., adapted, 'gather', nbp_glo, index_g)
    IF (ALL(adapted(:,:) == val_exp)) adapted(:,:) = zero
    !-
    regenerate(:,:) = val_exp
    var_name = 'regenerate'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &              .TRUE., regenerate, 'gather', nbp_glo, index_g)
    IF (ALL(regenerate(:,:) == val_exp)) regenerate(:,:) = zero
    !-
    npp_longterm(:,:) = val_exp
    var_name = 'npp_longterm'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &              .TRUE., npp_longterm, 'gather', nbp_glo, index_g)
    IF (ALL(npp_longterm(:,:) == val_exp)) npp_longterm(:,:) = zero
    !-
    lm_lastyearmax(:,:) = val_exp
    var_name = 'lm_lastyearmax'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &              .TRUE., lm_lastyearmax, 'gather', nbp_glo, index_g)
    IF (ALL(lm_lastyearmax(:,:) == val_exp)) lm_lastyearmax(:,:) = zero
    !-
    lm_thisyearmax(:,:) = val_exp
    var_name = 'lm_thisyearmax'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &              .TRUE., lm_thisyearmax, 'gather', nbp_glo, index_g)
    IF (ALL(lm_thisyearmax(:,:) == val_exp)) lm_thisyearmax(:,:) = zero
    !-
    maxfpc_lastyear(:,:) = val_exp
    var_name = 'maxfpc_lastyear'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &              .TRUE., maxfpc_lastyear, 'gather', nbp_glo, index_g)
    IF (ALL(maxfpc_lastyear(:,:) == val_exp)) maxfpc_lastyear(:,:) = zero
    !-
    maxfpc_thisyear(:,:) = val_exp
    var_name = 'maxfpc_thisyear'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &              .TRUE., maxfpc_thisyear, 'gather', nbp_glo, index_g)
    IF (ALL(maxfpc_thisyear(:,:) == val_exp)) maxfpc_thisyear(:,:) = zero
    !-
    turnover_time(:,:) = val_exp
    var_name = 'turnover_time'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &              .TRUE., turnover_time, 'gather', nbp_glo, index_g)
    IF ( ALL( turnover_time(:,:) == val_exp)) turnover_time(:,:) = 100.
    !-
    turnover_longterm(:,:,:,:) = val_exp
    DO l = 1,nelements
       DO k = 1,nparts
          WRITE(part_str,'(I2)') k
          IF ( k < 10 ) part_str(1:1) = '0'
          var_name = 'turnover_loterm_'//part_str(1:LEN_TRIM(part_str))//element_str(l)
          CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
               &              .TRUE., turnover_longterm(:,:,k,l), 'gather', nbp_glo, index_g)
          IF (ALL(turnover_longterm(:,:,k,l) == val_exp)) &
               &       turnover_longterm(:,:,k,l) = zero
       ENDDO
    END DO
    !-
    gpp_week(:,:) = val_exp
    var_name = 'gpp_week'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &              .TRUE., gpp_week, 'gather', nbp_glo, index_g)
    IF (ALL(gpp_week(:,:) == val_exp)) gpp_week(:,:) = zero
    !-
    npp_week(:,:) = val_exp
    var_name = 'npp_week'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &              .TRUE., npp_week, 'gather', nbp_glo, index_g)
    IF (ALL(npp_week(:,:) == val_exp)) npp_week(:,:) = zero
    !-
    biomass(:,:,:,:) = val_exp
    DO l = 1,nelements
       DO k = 1,nparts
          WRITE(part_str,'(I2)') k
          IF ( k < 10 ) part_str(1:1) = '0'
          var_name = 'biomass_'//part_str(1:LEN_TRIM(part_str))//element_str(l)
          CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
               &                   .TRUE., biomass(:,:,k,l), 'gather', nbp_glo, index_g)
          IF (ALL(biomass(:,:,k,l) == val_exp)) biomass(:,:,k,l) = zero
       ENDDO
    END DO
    !-
    resp_maint_part(:,:,:) = val_exp
    DO k=1,nparts
       WRITE(part_str,'(I2)') k
       IF ( k < 10 ) part_str(1:1) = '0'
       var_name = 'maint_resp_'//part_str(1:LEN_TRIM(part_str))
       CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
            &                   .TRUE., resp_maint_part(:,:,k), 'gather', nbp_glo, index_g)
       IF (ALL(resp_maint_part(:,:,k) == val_exp)) resp_maint_part(:,:,k) = zero
    ENDDO
    !JC separate tissue age
    !! old leaf_age leaf_frac
    !    !-
    !    leaf_age(:,:,:) = val_exp
    !    DO m=1,nleafages
    !       WRITE (part_str,'(I2)') m
    !       IF ( m < 10 ) part_str(1:1) = '0'
    !       var_name = 'leaf_age_'//part_str(1:LEN_TRIM(part_str))
    !       CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
    !            &                   .TRUE., leaf_age(:,:,m), 'gather', nbp_glo, index_g)
    !       IF (ALL(leaf_age(:,:,m) == val_exp)) leaf_age(:,:,m) = zero
    !    ENDDO
    !    !-
    !    leaf_frac(:,:,:) = val_exp
    !    DO m=1,nleafages
    !       WRITE(part_str,'(I2)') m
    !       IF ( m < 10 ) part_str(1:1) = '0'
    !       var_name = 'leaf_frac_'//part_str(1:LEN_TRIM(part_str))
    !       CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
    !            &                  .TRUE., leaf_frac(:,:,m), 'gather', nbp_glo, index_g)
    !       IF (ALL(leaf_frac(:,:,m) == val_exp)) leaf_frac(:,:,m) = zero
    !    ENDDO
    !-
    leaf_age(:,:,:,:) = val_exp
    DO l = 1,nleafages
       DO k = 1,nparts
          WRITE(part_str,'(I2)') k
          IF ( k < 10 ) part_str(1:1) = '0'
          var_name = 'leaf_age_'//part_str(1:LEN_TRIM(part_str))//'_'//TRIM(age_str(l))
          CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
               &                   .TRUE., leaf_age(:,:,k,l), 'gather', nbp_glo, index_g)
          IF (ALL(leaf_age(:,:,k,l) == val_exp)) leaf_age(:,:,k,l) = zero
       ENDDO
    END DO
    !-
    leaf_frac(:,:,:,:) = val_exp
    DO l = 1,nleafages
       DO k = 1,nparts
          WRITE(part_str,'(I2)') k
          IF ( k < 10 ) part_str(1:1) = '0'
          var_name = 'leaf_frac_'//part_str(1:LEN_TRIM(part_str))//'_'//TRIM(age_str(l))
          CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
               &                   .TRUE., leaf_frac(:,:,k,l), 'gather', nbp_glo, index_g)
          IF (ALL(leaf_frac(:,:,k,l) == val_exp)) leaf_frac(:,:,k,l) = zero
       ENDDO
    END DO
    !End JC separate tissue age
    !-
    senescence_real(:,:) = val_exp
    var_name = 'senescence'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &                .TRUE., senescence_real, 'gather', nbp_glo, index_g)
    IF (ALL(senescence_real(:,:) == val_exp)) senescence_real(:,:) = zero
    WHERE ( senescence_real(:,:) >= .5 )
       senescence = .TRUE.
    ELSEWHERE
       senescence = .FALSE.
    ENDWHERE


    ! Read real value for begin_leaves
    CALL restget_p (rest_id_stomate, 'begin_leaves', nbp_glo, nvm  , 1, itime, &
         .TRUE., begin_leaves_real, 'gather', nbp_glo, index_g)
    IF (ALL(begin_leaves_real(:,:) == val_exp)) begin_leaves_real(:,:) = zero

    ! Transform into logical needed by the modele
    WHERE ( begin_leaves_real(:,:) >= 0.5 )
       begin_leaves = .TRUE.
    ELSEWHERE
       begin_leaves = .FALSE.
    ENDWHERE


    when_growthinit(:,:) = val_exp
    var_name = 'when_growthinit'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &                .TRUE., when_growthinit, 'gather', nbp_glo, index_g)
    IF (ALL(when_growthinit(:,:) == val_exp)) &
         &     when_growthinit(:,:) = zero
    !-
    age(:,:) = val_exp
    var_name = 'age'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &                .TRUE., age, 'gather', nbp_glo, index_g)
    IF (ALL(age(:,:) == val_exp)) age(:,:) = zero
    !-
    ! 13 CO2
    !-
    resp_hetero(:,:) = val_exp
    var_name = 'resp_hetero'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                .TRUE., resp_hetero, 'gather', nbp_glo, index_g)
    IF (ALL(resp_hetero(:,:) == val_exp)) resp_hetero(:,:) = zero
    !-
    resp_maint(:,:) = val_exp
    var_name = 'resp_maint'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &                .TRUE., resp_maint, 'gather', nbp_glo, index_g)
    IF (ALL(resp_maint(:,:) == val_exp)) resp_maint(:,:) = zero
    !-
    resp_growth(:,:) = val_exp
    var_name = 'resp_growth'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &                .TRUE., resp_growth, 'gather', nbp_glo, index_g)
    IF (ALL(resp_growth(:,:) == val_exp)) resp_growth(:,:) = zero
    !-
    co2_fire(:,:) = val_exp
    var_name = 'co2_fire'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &                .TRUE., co2_fire, 'gather', nbp_glo, index_g)
    IF (ALL(co2_fire(:,:) == val_exp)) co2_fire(:,:) = zero
    !-
    co2_to_bm_dgvm(:,:) = val_exp
    var_name = 'co2_to_bm_dgvm'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &                .TRUE., co2_to_bm_dgvm, 'gather', nbp_glo, index_g)
    IF (ALL(co2_to_bm_dgvm(:,:) == val_exp)) co2_to_bm_dgvm(:,:) = zero
    !-
    n_to_bm(:,:) = val_exp
    var_name = 'n_to_bm'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &                .TRUE., n_to_bm, 'gather', nbp_glo, index_g)
    IF (ALL(n_to_bm(:,:) == val_exp)) n_to_bm(:,:) = zero
    !-
    p_to_bm(:,:) = val_exp
    var_name = 'p_to_bm'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &                .TRUE., p_to_bm, 'gather', nbp_glo, index_g)
    IF (ALL(p_to_bm(:,:) == val_exp)) p_to_bm(:,:) = zero
    !-
    ! 14 vegetation distribution after last light competition
    !-
    veget_lastlight(:,:) = val_exp
    var_name = 'veget_lastlight'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &                .TRUE., veget_lastlight, 'gather', nbp_glo, index_g)
    IF (ALL(veget_lastlight(:,:) == val_exp)) veget_lastlight(:,:) = zero
    !-
    ! 15 establishment criteria
    !-
    everywhere(:,:) = val_exp
    var_name = 'everywhere'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &                .TRUE., everywhere, 'gather', nbp_glo, index_g)
    IF (ALL(everywhere(:,:) == val_exp)) everywhere(:,:) = zero
    !-
    need_adjacent_real(:,:) = val_exp
    var_name = 'need_adjacent'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &                .TRUE., need_adjacent_real, 'gather', nbp_glo, index_g)
    IF (ALL(need_adjacent_real(:,:) == val_exp)) &
         &     need_adjacent_real(:,:) = zero
    WHERE ( need_adjacent_real(:,:) >= .5 )
       need_adjacent = .TRUE.
    ELSEWHERE
       need_adjacent = .FALSE.
    ENDWHERE
    !-
    RIP_time(:,:) = val_exp
    var_name = 'RIP_time'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &                .TRUE., RIP_time, 'gather', nbp_glo, index_g)
    IF (ALL(RIP_time(:,:) == val_exp)) RIP_time(:,:) = large_value
    !-
    ! 17 litter
    !-
    litter(:,:,:,:,:) = val_exp
    DO k = 1,nelements
       DO l = 1,nlevs
          DO m = 1,nvm
             WRITE (part_str, '(I2)') m
             IF (m<10) part_str(1:1)='0'
             var_name = 'litter_'//part_str(1:LEN_TRIM(part_str))//'_'//level_str(l)//element_str(k)
             CALL restget_p (rest_id_stomate, var_name, nbp_glo, nlitt , 1, itime, &
                  &                     .TRUE., litter(:,:,m,l,k), 'gather', nbp_glo, index_g)
             IF (ALL(litter(:,:,m,l,k) == val_exp)) litter(:,:,m,l,k) = zero
          ENDDO
       ENDDO
    END DO
    !-
    dead_leaves(:,:,:) = val_exp
    DO l=1,nlitt
       var_name = 'dead_leaves_'//TRIM(litter_str(l))
       CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
            &                   .TRUE., dead_leaves(:,:,l), 'gather', nbp_glo, index_g)
       IF (ALL(dead_leaves(:,:,l) == val_exp)) dead_leaves(:,:,l) = zero
    ENDDO
    !-
    som(:,:,:,:) = val_exp
    DO m=1,nvm
       WRITE (part_str, '(I2)') m
       IF (m<10) part_str(1:1)='0'
       var_name = 'carbon_'//part_str(1:LEN_TRIM(part_str))
       CALL restget_p (rest_id_stomate, var_name, nbp_glo, ncarb , 1, itime, &
            &                   .TRUE., som(:,:,m,icarbon), 'gather', nbp_glo, index_g) 
       IF (ALL(som(:,:,m,icarbon) == val_exp)) THEN 
        !  som(:,iactive ,m,icarbon) = 1000. 
        !  som(:,isurface,m,icarbon) = 1000. 
        !  som(:,islow   ,m,icarbon) = 3000. 
        !  som(:,ipassive,m,icarbon) = 5000.
          som(:,iactive ,m,icarbon) = 1. 
          som(:,isurface,m,icarbon) = 1. 
          som(:,islow   ,m,icarbon) = 3. 
          som(:,ipassive,m,icarbon) = 5.
       ENDIF
    ENDDO

    DO m=1,nvm 
       WRITE (part_str, '(I2)') m 
       IF (m<10) part_str(1:1)='0' 
       var_name = 'nitrogen_'//part_str(1:LEN_TRIM(part_str)) 
       CALL restget_p (rest_id_stomate, var_name, nbp_glo, ncarb , 1, itime, & 
            &                   .TRUE., som(:,:,m,initrogen), 'gather', nbp_glo, index_g) 
       IF (ALL(som(:,:,m,initrogen) == val_exp)) THEN 
          som(:,iactive,m,initrogen)  = som(:,iactive,m,icarbon) / CN_target_iactive_ref 
          som(:,isurface,m,initrogen) = som(:,isurface,m,icarbon) / CN_target_isurface_ref 
          som(:,islow,m,initrogen)    = som(:,islow,m,icarbon) / CN_target_islow_ref 
          som(:,ipassive,m,initrogen) =  som(:,ipassive,m,icarbon) / CN_target_ipassive_ref 
       ENDIF
    ENDDO

    DO m=1,nvm 
       WRITE (part_str, '(I2)') m 
       IF (m<10) part_str(1:1)='0' 
       var_name = 'phosphorus_'//part_str(1:LEN_TRIM(part_str)) 
       CALL restget_p (rest_id_stomate, var_name, nbp_glo, ncarb , 1, itime, & 
            &                   .TRUE., som(:,:,m,iphosphorus), 'gather', nbp_glo, index_g) 
       IF (ALL(som(:,:,m,iphosphorus) == val_exp)) THEN 
          som(:,iactive ,m,iphosphorus) = som(:,iactive ,m,icarbon) /  &
                                          (CN_target_iactive_ref  * NP_target_iactive_ref)
          som(:,isurface,m,iphosphorus) = som(:,isurface,m,icarbon) /  &
                                          (CN_target_isurface_ref * NP_target_isurface_ref)
          som(:,islow   ,m,iphosphorus) = som(:,islow   ,m,icarbon) /  &
                                           (CN_target_islow_ref    * NP_target_islow_ref)
          som(:,ipassive,m,iphosphorus) = som(:,ipassive,m,icarbon) /  &
                                           (CN_target_ipassive_ref * NP_target_ipassive_ref) 
       ENDIF
    ENDDO


    !-
    lignin_struc(:,:,:) = val_exp
    DO l=1,nlevs
       var_name = 'lignin_struc_'//level_str(l)
       CALL restget_p &
            &    (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
            &     .TRUE., lignin_struc(:,:,l), 'gather', nbp_glo, index_g)
       IF (ALL(lignin_struc(:,:,l) == val_exp)) lignin_struc(:,:,l) = zero
    ENDDO
    !-
    lignin_wood(:,:,:) = val_exp
    DO l=1,nlevs
       var_name = 'lignin_wood_'//level_str(l)
       CALL restget_p &
            &    (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
            &     .TRUE., lignin_wood(:,:,l), 'gather', nbp_glo, index_g)
       IF (ALL(lignin_wood(:,:,l) == val_exp)) lignin_wood(:,:,l) = zero
    ENDDO
    !-
    ! 18 land cover change
    !-
    prod10(:,:) = val_exp
    var_name = 'prod10'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, 11     , 1, itime, &
         &                .TRUE., prod10, 'gather', nbp_glo, index_g)
    IF (ALL(prod10(:,:) == val_exp)) prod10(:,:) = zero

    prod100(:,:) = val_exp
    var_name = 'prod100'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, 101     , 1, itime, &
         &                .TRUE., prod100, 'gather', nbp_glo, index_g)
    IF (ALL(prod100(:,:) == val_exp)) prod100(:,:) = zero


    flux10(:,:) = val_exp
    var_name = 'flux10'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, 10     , 1, itime, &
         &                .TRUE., flux10, 'gather', nbp_glo, index_g)
    IF (ALL(flux10(:,:) == val_exp)) flux10(:,:) = zero

    flux100(:,:) = val_exp
    var_name = 'flux100'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, 100     , 1, itime, &
         &                .TRUE., flux100, 'gather', nbp_glo, index_g)
    IF (ALL(flux100(:,:) == val_exp)) flux100(:,:) = zero

    convflux(:) = val_exp
    var_name = 'convflux'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, 1     , 1, itime, &
         &              .TRUE., convflux, 'gather', nbp_glo, index_g)
    IF (ALL(convflux(:) == val_exp)) convflux(:) = zero

    cflux_prod10(:) = val_exp
    var_name = 'cflux_prod10'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, 1     , 1, itime, &
         &              .TRUE., cflux_prod10, 'gather', nbp_glo, index_g)
    IF (ALL(cflux_prod10(:) == val_exp)) cflux_prod10(:) = zero

    cflux_prod100(:) = val_exp
    var_name = 'cflux_prod100'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, 1     , 1, itime, &
         &              .TRUE., cflux_prod100, 'gather', nbp_glo, index_g)
    IF (ALL(cflux_prod100(:) == val_exp)) cflux_prod100(:) = zero

    !-
    ! 18-bis wood harvest
    !-
    IF (do_wood_harvest) THEN
       prod10_harvest(:,:) = val_exp
       var_name = 'prod10_harvest'
       CALL restget_p (rest_id_stomate, var_name, nbp_glo, 11     , 1, itime, &
            .TRUE., prod10_harvest, 'gather', nbp_glo, index_g)
       IF (ALL(prod10_harvest(:,:) == val_exp)) prod10_harvest(:,:) = zero
       
       prod100_harvest(:,:) = val_exp
       var_name = 'prod100_harvest'
       CALL restget_p (rest_id_stomate, var_name, nbp_glo, 101     , 1, itime, &
            .TRUE., prod100_harvest, 'gather', nbp_glo, index_g)
       IF (ALL(prod100_harvest(:,:) == val_exp)) prod100_harvest(:,:) = zero
       
       flux10_harvest(:,:) = val_exp
       var_name = 'flux10_harvest'
       CALL restget_p (rest_id_stomate, var_name, nbp_glo, 10     , 1, itime, &
            .TRUE., flux10_harvest, 'gather', nbp_glo, index_g)
       IF (ALL(flux10_harvest(:,:) == val_exp)) flux10_harvest(:,:) = zero
       
       flux100_harvest(:,:) = val_exp
       var_name = 'flux100_harvest'
       CALL restget_p (rest_id_stomate, var_name, nbp_glo, 100     , 1, itime, &
            .TRUE., flux100_harvest, 'gather', nbp_glo, index_g)
       IF (ALL(flux100_harvest(:,:) == val_exp)) flux100_harvest(:,:) = zero
       
       convflux_harvest(:) = val_exp
       var_name = 'convflux_harvest'
       CALL restget_p (rest_id_stomate, var_name, nbp_glo, 1     , 1, itime, &
            .TRUE., convflux_harvest, 'gather', nbp_glo, index_g)
       IF (ALL(convflux_harvest(:) == val_exp)) convflux_harvest(:) = zero
       
       cflux_prod10_harvest(:) = val_exp
       var_name = 'cflux_prod10_harvest'
       CALL restget_p (rest_id_stomate, var_name, nbp_glo, 1     , 1, itime, &
            .TRUE., cflux_prod10_harvest, 'gather', nbp_glo, index_g)
       IF (ALL(cflux_prod10_harvest(:) == val_exp)) cflux_prod10_harvest(:) = zero
       
       cflux_prod100_harvest(:) = val_exp
       var_name = 'cfluxprod100_harvest'
       CALL restget_p (rest_id_stomate, var_name, nbp_glo, 1     , 1, itime, &
            .TRUE., cflux_prod100_harvest, 'gather', nbp_glo, index_g)
       IF (ALL(cflux_prod100_harvest(:) == val_exp)) cflux_prod100_harvest(:) = zero
       
       woodharvestpft(:,:) = val_exp
       var_name = 'woodharvestpft'
       CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
            .TRUE., woodharvestpft, 'gather', nbp_glo, index_g)
       IF (ALL(woodharvestpft(:,:) == val_exp)) woodharvestpft(:,:) = zero
    END IF


    bm_to_litter(:,:,:,:) = val_exp
    DO l = 1,nelements
       DO k = 1,nparts
          WRITE(part_str,'(I2)') k
          IF ( k < 10 ) part_str(1:1) = '0'
          var_name = 'bm_to_litter_'//part_str(1:LEN_TRIM(part_str))//element_str(l)
          CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
               &                .TRUE., bm_to_litter(:,:,k,l), & 
               &'gather', nbp_glo, index_g)
          IF (ALL(bm_to_litter(:,:,k,l) == val_exp)) bm_to_litter(:,:,k,l) = zero
       ENDDO
    END DO

    carb_mass_total(:) = val_exp
    var_name = 'carb_mass_total'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, 1     , 1, itime, &
         &              .TRUE., carb_mass_total, 'gather', nbp_glo, index_g)
    IF (ALL(carb_mass_total(:) == val_exp)) carb_mass_total(:) = zero

    ! gmjc
    wshtotsum(:,:) = val_exp
    var_name = 'wshtotsum'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
   &              .TRUE., wshtotsum, 'gather', nbp_glo, index_g)
      IF (ALL(wshtotsum(:,:) == val_exp)) wshtotsum(:,:) = zero
  !-
    sr_ugb(:,:) = val_exp
    var_name = 'sr_ugb'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
   &              .TRUE., sr_ugb, 'gather', nbp_glo, index_g)
      IF (ALL(sr_ugb(:,:) == val_exp)) sr_ugb(:,:) = zero
  !-
    sla_calc(:,:) = val_exp
    var_name = 'sla_calc'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
   &              .TRUE., sla_calc, 'gather', nbp_glo, index_g)
  !    IF (ALL(sla_calc(:,:) == val_exp)) sla_calc(:,:) = zero
    IF (ALL(sla_calc(:,:) == val_exp)) THEN
       DO j=1,nvm
          sla_calc(:,j) = sla(j)
       END DO
    END IF
  !-
    nb_ani(:,:) = val_exp
    var_name = 'nb_ani'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
   &              .TRUE., nb_ani, 'gather', nbp_glo, index_g)
      IF (ALL(nb_ani(:,:) == val_exp)) nb_ani(:,:) = zero
  !-
    grazed_frac(:,:) = val_exp
    var_name = 'grazed_frac'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
   &              .TRUE., grazed_frac, 'gather', nbp_glo, index_g)
      IF (ALL(grazed_frac(:,:) == val_exp)) grazed_frac(:,:) = zero
  !-
    import_yield(:,:) = val_exp

    var_name = 'import_yield'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
   &              .TRUE., import_yield, 'gather', nbp_glo, index_g)
      IF (ALL(import_yield(:,:) == val_exp)) import_yield(:,:) = zero
  
  !-
     t2m_14(:) = val_exp
    var_name = 't2m_14'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, 1     , 1, itime, &
   &              .TRUE., t2m_14, 'gather', nbp_glo, index_g)
    IF (ALL(t2m_14(:) == val_exp)) t2m_14(:) = t2m(:)
  !
  
  ! in current trunk version there is no 3D restart
      litter_not_avail(:,:,:) = val_exp
      DO l=1,nlitt
         var_name = 'litter_not_avail_'//TRIM(litter_str(l))
         CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
              &                   .TRUE., litter_not_avail(:,l,:), 'gather', nbp_glo, index_g)
         IF (ALL(litter_not_avail(:,l,:) == val_exp)) litter_not_avail(:,l,:) = zero
      ENDDO
  !-
    nb_grazingdays(:,:) = val_exp
    var_name = 'nb_grazingdays'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
   &              .TRUE., nb_grazingdays, 'gather', nbp_glo, index_g)
      IF (ALL(nb_grazingdays(:,:) == val_exp)) nb_grazingdays(:,:) = zero
  
  !-
    after_snow(:) = val_exp
    var_name = 'after_snow'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, 1     , 1, itime, &
   &              .TRUE., after_snow, 'gather', nbp_glo, index_g)
      IF (ALL(after_snow(:) == val_exp)) after_snow(:) = zero
  !-
    after_wet(:) = val_exp
    var_name = 'after_wet'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, 1     , 1, itime, &
   &              .TRUE., after_wet, 'gather', nbp_glo, index_g)
      IF (ALL(after_wet(:) == val_exp)) after_wet(:) = zero
  !-
    wet1day(:) = val_exp
    var_name = 'wet1day'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, 1     , 1, itime, &
   &              .TRUE., wet1day, 'gather', nbp_glo, index_g)
      IF (ALL(wet1day(:) == val_exp)) wet1day(:) = 6
  !-
    wet2day(:) = val_exp
    var_name = 'wet2day'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, 1     , 1, itime, &
   &              .TRUE., wet2day, 'gather', nbp_glo, index_g)
      IF (ALL(wet2day(:) == val_exp)) wet2day(:) = 6
    GRM_devstage(:,:) = val_exp
    var_name = 'GRM_devstage'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
   &              .TRUE., GRM_devstage, 'gather', nbp_glo, index_g)
      IF (ALL(GRM_devstage(:,:) == val_exp)) GRM_devstage(:,:) = 2.0
    days_senescence(:,:) = val_exp
    var_name = 'days_senescence'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
   &              .TRUE., days_senescence, 'gather', nbp_glo, index_g)
      IF (ALL(days_senescence(:,:) == val_exp)) days_senescence(:,:) = zero
  !-
  !end gmjc


    !-
    ! 19. Spinup
    !-
    IF (spinup_analytic) THEN

       IF (is_root_prc) THEN
          temp_global_years(1) = val_exp
          var_name = 'Global_years'
          CALL restget (rest_id_stomate, var_name, 1 ,1  , 1, itime, &
               &                .TRUE., temp_global_years)
          IF(temp_global_years(1) == val_exp) temp_global_years(1) = zero
          global_years = INT(temp_global_years(1))
       ENDIF
       CALL bcast(global_years)

       nbp_accu(:) = val_exp
       var_name = 'nbp_sum'
       CALL restget_p (rest_id_stomate, var_name, nbp_glo, 1     , 1, itime, &
            &              .TRUE., nbp_accu, 'gather', nbp_glo, index_g)
       IF (ALL(nbp_accu(:) == val_exp)) nbp_accu(:) = zero    

       nbp_flux(:) = val_exp
       var_name = 'nbp_flux'
       CALL restget_p (rest_id_stomate, var_name, nbp_glo, 1     , 1, itime, &
            &              .TRUE., nbp_flux, 'gather', nbp_glo, index_g)
       IF (ALL(nbp_flux(:) == val_exp)) nbp_flux(:) = zero     

       !-
       ok_equilibrium_real(:) = val_exp
       var_name = 'ok_equilibrium'
       CALL restget_p (rest_id_stomate, var_name, nbp_glo , 1  , 1, itime, &
            &                .TRUE., ok_equilibrium_real,'gather', nbp_glo, index_g)
       IF (ALL(ok_equilibrium_real(:) == val_exp)) ok_equilibrium_real(:) = zero
       WHERE(ok_equilibrium_real(:) >= 0.5) 
          ok_equilibrium = .TRUE.
       ELSEWHERE
          ok_equilibrium = .FALSE.
       ENDWHERE

       MatrixV(:,:,:,:) = val_exp
       DO k = 1,nbpools
          DO j = 1,nbpools
             WRITE(part_str,'(I2)') k
             IF (k < 10) part_str(1:1) = '0'             
             var_name = 'MatrixV_'//part_str(1:LEN_TRIM(part_str))//'_'//TRIM(pools_str(j))
             CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm , 1, itime, &
                  &                     .TRUE., MatrixV(:,:,k,j), 'gather', nbp_glo, index_g)
          ENDDO
       ENDDO
       ! If nothing is found in the restart file, we initialize each submatrix by identity
       IF (ALL(MatrixV(:,:,:,:) == val_exp))  THEN 
          MatrixV(:,:,:,:) = zero
          DO l = 1,nbpools
             MatrixV(:,:,l,l) = un
          END DO
       END IF

       VectorU(:,:,:)  = val_exp
       DO k= 1,nbpools
          WRITE(part_str,'(I2)') k
          IF (k < 10) part_str(1:1) = '0' 
          var_name = 'Vector_U_'//part_str(1:LEN_TRIM(part_str))
          CALL restget_p &
               &    (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
               &     .TRUE., VectorU(:,:,k), 'gather', nbp_glo, index_g)
          IF (ALL(VectorU(:,:,k) == val_exp))  VectorU(:,:,k) = zero
       ENDDO
       
       previous_stock(:,:,:)  = val_exp
       DO k= 1,nbpools
          WRITE(part_str,'(I2)') k
          IF (k < 10) part_str(1:1) = '0' 
          var_name = 'previous_stock_'//part_str(1:LEN_TRIM(part_str))
          CALL restget_p &
               &    (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
               &     .TRUE., previous_stock(:,:,k), 'gather', nbp_glo, index_g)
          IF (ALL(previous_stock(:,:,k) == val_exp))  previous_stock(:,:,k) = undef_sechiba
       ENDDO
       
       current_stock(:,:,:)  = val_exp
       DO k= 1,nbpools
          WRITE(part_str,'(I2)') k
          IF (k < 10) part_str(1:1) = '0' 
          var_name = 'current_stock_'//part_str(1:LEN_TRIM(part_str))
          CALL restget_p &
               &    (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
               &     .TRUE., current_stock(:,:,k), 'gather', nbp_glo, index_g)
          IF (ALL(current_stock(:,:,k) == val_exp))  current_stock(:,:,k) = zero
       ENDDO

       CN_som_litter_longterm(:,:,:)  = val_exp
       DO k= 1,nbpools
          WRITE(part_str,'(I2)') k
          IF (k < 10) part_str(1:1) = '0' 
          var_name = 'CN_longterm_'//part_str(1:LEN_TRIM(part_str))
          CALL restget_p &
               &    (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
               &     .TRUE., CN_som_litter_longterm(:,:,k), 'gather', nbp_glo, index_g)
          IF (ALL(CN_som_litter_longterm(:,:,k) == val_exp))  CN_som_litter_longterm(:,:,k) = zero
       ENDDO

       vartmp(1) = val_exp
       var_name = 'tau_CN_longterm'
       IF (is_root_prc) THEN
          CALL restget (rest_id_stomate, var_name, 1 ,1  , 1, itime, &
               .TRUE., vartmp)
          IF (vartmp(1) == val_exp)  THEN
              tau_CN_longterm = dt_sechiba/one_day
          ELSE
              tau_CN_longterm = vartmp(1)
          ENDIF
       ENDIF

       CP_som_litter_longterm(:,:,:)  = val_exp
       DO k= 1,nbpools
          WRITE(part_str,'(I2)') k
          IF (k < 10) part_str(1:1) = '0' 
          var_name = 'CP_longterm_'//part_str(1:LEN_TRIM(part_str))
          CALL restget_p &
               &    (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
               &     .TRUE., CP_som_litter_longterm(:,:,k), 'gather', nbp_glo, index_g)
          IF (ALL(CP_som_litter_longterm(:,:,k) == val_exp))  CP_som_litter_longterm(:,:,k) = zero
       ENDDO

       CALL bcast(tau_CN_longterm)

 
         
    ENDIF ! spinup_matrix_method

    record_reserve_grass(:,:) = val_exp
    var_name = 'record_reserve_grass'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &                .TRUE., record_reserve_grass, 'gather', nbp_glo, index_g)
    IF (ALL(record_reserve_grass(:,:) == val_exp)) record_reserve_grass(:,:) = 200.

    KF(:,:) = val_exp
    var_name = 'KF'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &                .TRUE., KF, 'gather', nbp_glo, index_g)
    ! I don't want to set it equal to zero, since this is a problem if these
    ! values are not here!  Better it blows up later on
    !IF (ALL(KF(:,:) == val_exp)) KF(:,:) = zero

    k_latosa_adapt(:,:) = val_exp
    var_name = 'k_latosa_adapt'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &                .TRUE., k_latosa_adapt, 'gather', nbp_glo, index_g)
    DO m = 1,nvm
       IF (ALL(k_latosa_adapt(:,m) == val_exp)) k_latosa_adapt(:,m) = k_latosa_min(m)
    ENDDO

    rue_longterm(:,:) = val_exp
    var_name = 'rue_longterm'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm , 1, itime, &
         &        .TRUE., rue_longterm(:,:), 'gather', nbp_glo, index_g)
    IF (ALL(rue_longterm(:,:) == val_exp)) rue_longterm(:,:) = 1.

    cn_leaf_avg_season(:,:) = val_exp 
    var_name = 'cn_leaf_avg_season' 
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, & 
         &              .TRUE., cn_leaf_avg_season, 'gather', nbp_glo, index_g) 
    IF ( ALL(cn_leaf_avg_season(:,:) == val_exp) ) THEN 
       DO m=1,nvm 
          cn_leaf_avg_season(:,m) = cn_leaf_min(m)
       ENDDO
    ENDIF
    
    nstress_season(:,:) = val_exp 
    var_name = 'nstress_season' 
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, & 
         &              .TRUE., nstress_season, 'gather', nbp_glo, index_g) 
    IF ( ALL(nstress_season(:,:) == val_exp) ) nstress_season(:,:)=1.0 

    pstress_season(:,:) = val_exp 
    var_name = 'pstress_season' 
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, & 
         &              .TRUE., pstress_season, 'gather', nbp_glo, index_g) 
    IF ( ALL(pstress_season(:,:) == val_exp) ) pstress_season(:,:)=1.0 

    
    soil_n_min(:,:,:) = val_exp 
    DO k=1,nnspec 
       WRITE(part_str,'(I1)') k 
       var_name = 'soil_n_min_'//part_str(1:LEN_TRIM(part_str)) 
       CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, & 
            &              .TRUE., soil_n_min(:,:,k), 'gather', nbp_glo, index_g) 
       IF ( ALL(soil_n_min(:,:,k) == val_exp) ) soil_n_min(:,:,k)=2.
    ENDDO

    np_leaf_avg_season(:,:) = val_exp 
    var_name = 'np_leaf_avg_season' 
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, & 
         &              .TRUE., np_leaf_avg_season, 'gather', nbp_glo, index_g) 
    IF ( ALL(np_leaf_avg_season(:,:) == val_exp) ) THEN 
       DO m=1,nvm 
          !DSGlab03 np_leaf_avg_season(:,m) = np_leaf_init(m) 
          np_leaf_avg_season(:,m) = np_leaf_min(m)
          !DSGlab03 
       ENDDO
    ENDIF
    
    soil_p_min(:,:,:) = val_exp 
    DO k=1,npspec 
       WRITE(part_str,'(I1)') k 
       var_name = 'soil_p_min_'//part_str(1:LEN_TRIM(part_str)) 
       CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, & 
            &              .TRUE., soil_p_min(:,:,k), 'gather', nbp_glo, index_g) 
       !DSG: never inititialze any value here, then "undef"
       ! if you wan to give an initial amount to mineral P 
       ! use the variable "labile_init" in src_parameters/constantes_var.f90
       IF ( ALL(soil_p_min(:,:,k) == val_exp) ) soil_p_min(:,:,k)=undef
       ! 
    ENDDO

    f_Pdissolved(:,:) = val_exp 
    var_name = 'f_Pdissolved'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, & 
         &              .TRUE., f_Pdissolved, 'gather', nbp_glo, index_g) 
    IF ( ALL(f_Pdissolved(:,:) == val_exp) ) THEN 
      f_Pdissolved(:,:) = un
    ENDIF

    !-
    EW_grain(:,:,:,:) = val_exp
    DO l = 1,nminerals
       DO k = 1,2
          WRITE(part_str,'(I2)') k
          IF ( k < 10 ) part_str(1:1) = '0'
          var_name = 'EW_grain_'//part_str(1:LEN_TRIM(part_str))//element_str(l)
          CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
               &                   .TRUE., EW_grain(:,:,l,k), 'gather', nbp_glo, index_g)
       ENDDO
       ! If not found in restart file:
       IF (ALL(EW_grain(:,:,l,1) == val_exp)) EW_grain(:,:,l,1) = zero  ! 
       IF (ALL(EW_grain(:,:,l,2) == val_exp)) EW_grain(:,:,l,2) = undef ! 
       ! If use the the force stock mode, we set the size:
       IF (EW_FORCE_STOCK) EW_grain(:,:,l,2) = EW_TargetSize(l)
    END DO

    IF (printlev>=1) THEN
        WRITE (numout,*) 'After: reading restart:'
        WRITE (numout,*) 'EW_target_stock:',EW_targetStock
        WRITE (numout,*) 'MAXVAL(EW_grain(:,:,:,1))',  MAXVAL(EW_grain(:,:,:,1))
        WRITE (numout,*) 'EW_target_size:' ,EW_targetSize
        WRITE (numout,*) 'MAXVAL(EW_grain(:,:,:,2))',  MAXVAL(EW_grain(:,:,:,2))
    ENDIF

    p_O2(:,:) = val_exp 
    var_name = 'pO2'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, & 
         &              .TRUE., p_O2(:,:), 'gather', nbp_glo, index_g) 
    IF ( ALL(p_O2(:,:) == val_exp) ) p_O2(:,:)=200 
 
    bact(:,:) = val_exp 
    var_name = 'bact'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, & 
         &              .TRUE., bact(:,:), 'gather', nbp_glo, index_g) 
    IF ( ALL(bact(:,:) == val_exp) ) bact(:,:)=10 
 


    ! Read assim_param from restart file. The initialization of assim_param will 
    ! be done in stomate_var_init if the variable is not in the restart file.
    assim_param(:,:,:)  = val_exp
    DO k= 1,npco2
       WRITE(part_str,'(I2)') k
       IF (k < 10) part_str(1:1) = '0' 
       var_name = 'assim_param_'//part_str(1:LEN_TRIM(part_str))
       CALL restget_p &
            &    (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
            &     .TRUE., assim_param(:,:,k), 'gather', nbp_glo, index_g)
    END DO
 
    IF (printlev >= 4) WRITE(numout,*) 'Leaving readstart'
    !-----------------------
  END SUBROUTINE readstart

!! ================================================================================================================================
!! SUBROUTINE   : writerestart
!!
!>\BRIEF        Write all variables for stomate from restart file. 
!!
!! DESCRIPTION  : Write all variables for stomate from restart file. 
!!                
!! \n
!_ ================================================================================================================================

  SUBROUTINE writerestart &
       & (npts, index, dt_days, date_loc, &
       &  ind, adapted, regenerate, moiavail_daily, max_eau_var_daily, &
       &  gdd_init_date, litterhum_daily, t2m_daily, t2m_min_daily, tsurf_daily, tsoil_daily, &
       &  soilhum_daily, &
       &  drainage_longterm, runoff_longterm, &
       &  precip_daily, gpp_daily, npp_daily, &
       &  turnover_daily, moiavail_month, moiavail_week, moiavail_growingseason,&
       &  t2m_longterm, tau_longterm, t2m_month, t2m_week, &
       &  tsoil_month, soilhum_month, fireindex, firelitter, &
       &  maxmoiavail_lastyear, maxmoiavail_thisyear, &
       &  minmoiavail_lastyear, minmoiavail_thisyear, &
       &  maxgppweek_lastyear, maxgppweek_thisyear, &
       &  gdd0_lastyear, gdd0_thisyear, precip_lastyear, precip_thisyear, &
       &  gdd_m5_dormance, gdd_from_growthinit, gdd_midwinter, ncd_dormance, ngd_minus5, &
       &  PFTpresent, npp_longterm, lm_lastyearmax, lm_thisyearmax, &
       &  maxfpc_lastyear, maxfpc_thisyear, &
       &  turnover_longterm, gpp_week, npp_week, biomass, resp_maint_part, &
       &  leaf_age, leaf_frac, senescence, when_growthinit, age, &
       &  resp_hetero, resp_maint, resp_growth, co2_fire, co2_to_bm_dgvm, &
       &  n_to_bm, p_to_bm, veget_lastlight, everywhere, need_adjacent, RIP_time, &
       &  time_hum_min, hum_min_dormance, &
       &  litter, dead_leaves, &
       &  som, lignin_struc, lignin_wood, turnover_time, &
       &  prod10,prod100 ,flux10, flux100, &
       &  convflux, cflux_prod10, cflux_prod100, & 
       &  prod10_harvest,prod100_harvest ,flux10_harvest, flux100_harvest, &
       &  convflux_harvest, cflux_prod10_harvest, cflux_prod100_harvest, &
       &  woodharvestpft, bm_to_litter, carb_mass_total, &
       &  Tseason, Tseason_length, Tseason_tmp, & 
       &  Tmin_spring_time, begin_leaves, onset_date, &
       &  global_years, ok_equilibrium, nbp_accu, nbp_flux, &
       &  MatrixV, VectorU, previous_stock, current_stock, assim_param, &
       &  CN_som_litter_longterm, tau_CN_longterm, KF, k_latosa_adapt, &
       &  CP_som_litter_longterm, np_leaf_avg_season, soil_p_min, & !DSG: P  variables
       &  f_Pdissolved, EW_grain, &
       &  rue_longterm, cn_leaf_avg_season, nstress_season, pstress_season, &
       &  lai_target_longterm, &
       &  soil_n_min, p_O2, bact, &
       &  record_reserve_grass, &
       &  wshtotsum, sr_ugb, sla_calc, nb_ani, grazed_frac, &
       &  import_yield, t2m_14, litter_not_avail, nb_grazingdays, &
       &  after_snow, after_wet, wet1day, wet2day, GRM_devstage, &
       &  days_senescence)
    !---------------------------------------------------------------------
    !- write restart file
    !---------------------------------------------------------------------
    !-
    ! 0 declarations
    !-
    ! 0.1 input
    !-
    ! Domain size
    INTEGER,INTENT(in) :: npts
    ! Indices of the points on the map
    INTEGER,DIMENSION(npts),INTENT(in) :: index
    ! time step of STOMATE in days
    REAL,INTENT(in) :: dt_days
    ! date_loc (d)
    INTEGER,INTENT(in) :: date_loc
    ! density of individuals (1/m**2)
    REAL,DIMENSION(npts,nvm),INTENT(in) :: ind
    ! Winter too cold? between 0 and 1
    REAL,DIMENSION(npts,nvm),INTENT(in) :: adapted
    ! Winter sufficiently cold? between 0 and 1
    REAL,DIMENSION(npts,nvm),INTENT(in) :: regenerate
    ! daily moisture availability
    REAL,DIMENSION(npts,nvm),INTENT(in) :: moiavail_daily
    ! gdd init date
    REAL,DIMENSION(npts,2),INTENT(in) :: gdd_init_date
    ! daily litter humidity
    REAL,DIMENSION(npts),INTENT(in) :: litterhum_daily
    ! daily 2 meter temperatures (K)
    REAL,DIMENSION(npts),INTENT(in) :: t2m_daily
    ! daily minimum 2 meter temperatures (K)
    REAL,DIMENSION(npts),INTENT(in) :: t2m_min_daily
    ! daily surface temperatures (K)
    REAL,DIMENSION(npts),INTENT(in) :: tsurf_daily
    ! daily soil temperatures (K)
    REAL,DIMENSION(npts,nslm),INTENT(in) :: tsoil_daily
    ! daily soil humidity
    REAL,DIMENSION(npts,nslm),INTENT(in) :: soilhum_daily
    ! daily precipitations (mm/day) (for phenology)
    REAL,DIMENSION(npts),INTENT(in) :: max_eau_var_daily
    ! daily precipitations (mm/day) (for phenology)

    REAL, DIMENSION(npts), INTENT(in)                           :: drainage_longterm       !! "long term" annual sum of drainage (kgH2O/m/yr)
    REAL, DIMENSION(npts), INTENT(in)                           :: runoff_longterm         !! "long term" annual sum of runoff   (kgH2O/m/yr)
    REAL, DIMENSION(npts,nvm), INTENT(in)                      :: lai_target_longterm     !! "long term" lai target [m2/m2]

    REAL,DIMENSION(npts),INTENT(in) :: precip_daily
    ! daily gross primary productivity (gC/m**2/day)
    REAL,DIMENSION(npts,nvm),INTENT(in) :: gpp_daily
    ! daily net primary productivity (gC/m**2/day)
    REAL,DIMENSION(npts,nvm),INTENT(in) :: npp_daily
    ! daily turnover rates (gC/m**2/day)
    REAL,DIMENSION(npts,nvm,nparts,nelements),INTENT(in) :: turnover_daily
    ! "monthly" moisture availability
    REAL,DIMENSION(npts,nvm),INTENT(in) :: moiavail_month
    ! "weekly" moisture availability
    REAL,DIMENSION(npts,nvm),INTENT(in) :: moiavail_week
    ! mean growing season moisture availability 
    REAL,DIMENSION(npts,nvm),INTENT(in) :: moiavail_growingseason
    ! "long term" 2 meter temperatures (K)
    REAL,DIMENSION(npts),INTENT(in) :: t2m_longterm
    ! "tau_longterm"
    REAL, INTENT(IN)             :: tau_longterm
    ! "monthly" 2 meter temperatures (K)
    REAL,DIMENSION(npts),INTENT(in) :: t2m_month
    ! "seasonal" 2 meter temperatures (K) 
    REAL,DIMENSION(npts),INTENT(in)      :: Tseason
    ! temporary variable to calculate Tseason
    REAL,DIMENSION(npts),INTENT(in)      :: Tseason_length
    ! temporary variable to calculate Tseason
    REAL,DIMENSION(npts),INTENT(in)      :: Tseason_tmp
    REAL,DIMENSION(npts,nvm),INTENT(in)  :: Tmin_spring_time
    REAL,DIMENSION(npts,nvm),INTENT(in)  :: onset_date
    LOGICAL,DIMENSION(npts,nvm),INTENT(in)      :: begin_leaves

    ! "weekly" 2 meter temperatures (K)
    REAL,DIMENSION(npts),INTENT(in) :: t2m_week
    ! "monthly" soil temperatures (K)
    REAL,DIMENSION(npts,nslm),INTENT(in) :: tsoil_month
    ! "monthly" soil humidity
    REAL,DIMENSION(npts,nslm),INTENT(in) :: soilhum_month
    ! Probability of fire
    REAL,DIMENSION(npts,nvm),INTENT(in) :: fireindex
    ! Longer term total litter above the ground, gC/m**2 of ground
    REAL,DIMENSION(npts,nvm),INTENT(in) :: firelitter
    ! last year's maximum moisture availability
    REAL,DIMENSION(npts,nvm),INTENT(in) :: maxmoiavail_lastyear
    ! this year's maximum moisture availability
    REAL,DIMENSION(npts,nvm),INTENT(in) :: maxmoiavail_thisyear
    ! last year's minimum moisture availability
    REAL,DIMENSION(npts,nvm),INTENT(in) :: minmoiavail_lastyear
    ! this year's minimum moisture availability
    REAL,DIMENSION(npts,nvm),INTENT(in) :: minmoiavail_thisyear
    ! last year's maximum weekly GPP
    REAL,DIMENSION(npts,nvm),INTENT(in) :: maxgppweek_lastyear
    ! this year's maximum weekly GPP
    REAL,DIMENSION(npts,nvm),INTENT(in) :: maxgppweek_thisyear
    ! last year's annual GDD0
    REAL,DIMENSION(npts),INTENT(in) :: gdd0_lastyear
    ! this year's annual GDD0
    REAL,DIMENSION(npts),INTENT(in) :: gdd0_thisyear
    ! last year's annual precipitation (mm/year)
    REAL,DIMENSION(npts),INTENT(in) :: precip_lastyear
    ! this year's annual precipitation (mm/year)
    REAL,DIMENSION(npts),INTENT(in) :: precip_thisyear
    ! growing degree days, threshold -5 deg C (for phenology)
    REAL,DIMENSION(npts,nvm),INTENT(in) :: gdd_m5_dormance
    ! growing degree days, from begin of season (crops)
    REAL,DIMENSION(npts,nvm),INTENT(in) :: gdd_from_growthinit
    ! growing degree days since midwinter (for phenology)
    REAL,DIMENSION(npts,nvm),INTENT(in) :: gdd_midwinter
    ! number of chilling days since leaves were lost (for phenology)
    REAL,DIMENSION(npts,nvm),INTENT(in) :: ncd_dormance
    ! number of growing days, threshold -5 deg C (for phenology)
    REAL,DIMENSION(npts,nvm),INTENT(in) :: ngd_minus5
    ! PFT exists (equivalent to fpc_max > 0 for natural PFTs)
    LOGICAL,DIMENSION(npts,nvm),INTENT(in) :: PFTpresent
    ! "long term" net primary productivity (gC/m**2/year)
    REAL,DIMENSION(npts,nvm),INTENT(in) :: npp_longterm
    ! last year's maximum leaf mass, for each PFT (gC/m**2)
    REAL,DIMENSION(npts,nvm),INTENT(in) :: lm_lastyearmax
    ! this year's maximum leaf mass, for each PFT (gC/m**2)
    REAL,DIMENSION(npts,nvm),INTENT(in) :: lm_thisyearmax
    ! last year's maximum fpc for each natural PFT, on ground
    REAL,DIMENSION(npts,nvm),INTENT(in) :: maxfpc_lastyear
    ! this year's maximum fpc for each PFT,
    ! on *total* ground (see stomate_season)
    REAL,DIMENSION(npts,nvm),INTENT(in) :: maxfpc_thisyear
    ! "long term" turnover rate (gC/m**2/year)
    REAL,DIMENSION(npts,nvm,nparts,nelements),INTENT(in) :: turnover_longterm
    ! "weekly" GPP (gC/day/(m**2 covered)
    REAL,DIMENSION(npts,nvm),INTENT(in) :: gpp_week
    ! biomass (gC/m**2)
    REAL,DIMENSION(npts,nvm),INTENT(in) :: npp_week
    ! biomass (gC/m**2)
    REAL,DIMENSION(npts,nvm,nparts,nelements),INTENT(in) :: biomass
    ! maintenance respiration (gC/m**2)
    REAL,DIMENSION(npts,nvm,nparts),INTENT(in) :: resp_maint_part
    ! leaf age (days)
    REAL,DIMENSION(npts,nvm,nparts,nleafages),INTENT(in) :: leaf_age
    ! fraction of leaves in leaf age class
    REAL,DIMENSION(npts,nvm,nparts,nleafages),INTENT(in) :: leaf_frac
    ! is the plant senescent ?
    ! (only for deciduous trees - carbohydrate reserve)
    LOGICAL,DIMENSION(npts,nvm),INTENT(in) :: senescence
    ! how many days ago was the beginning of the growing season
    REAL,DIMENSION(npts,nvm),INTENT(in) :: when_growthinit
    ! mean age (years)
    REAL,DIMENSION(npts,nvm),INTENT(in) :: age
    ! heterotrophic respiration (gC/day/m**2)
    REAL,DIMENSION(npts,nvm),INTENT(in) :: resp_hetero
    ! maintenance respiration (gC/day/m**2)
    REAL,DIMENSION(npts,nvm),INTENT(in) :: resp_maint
    ! growth respiration (gC/day/m**2)
    REAL,DIMENSION(npts,nvm),INTENT(in) :: resp_growth
    ! carbon emitted into the atmosphere by fire (living and dead biomass)
    ! (in gC/m**2/time step)
    REAL,DIMENSION(npts,nvm),INTENT(in) :: co2_fire
    ! biomass uptaken (gC/(m**2 of total ground)/day)
    REAL,DIMENSION(npts,nvm),INTENT(in) :: co2_to_bm_dgvm
    ! biomass uptaken (gN/(m**2 of total ground)/day)
    REAL,DIMENSION(npts,nvm),INTENT(in) :: n_to_bm
    REAL,DIMENSION(npts,nvm),INTENT(in) :: p_to_bm
    ! vegetation fractions (on ground) after last light competition
    REAL,DIMENSION(npts,nvm),INTENT(in) :: veget_lastlight
    ! is the PFT everywhere in the grid box or very localized
    ! (after its introduction)
    REAL,DIMENSION(npts,nvm),INTENT(in) :: everywhere
    ! in order for this PFT to be introduced,
    ! does it have to be present in an adjacent grid box?
    LOGICAL,DIMENSION(npts,nvm),INTENT(in) :: need_adjacent
    ! How much time ago was the PFT eliminated for the last time (y)
    REAL,DIMENSION(npts,nvm),INTENT(in) :: RIP_time
    ! time elapsed since strongest moisture availability (d)
    REAL,DIMENSION(npts,nvm),INTENT(in) :: time_hum_min
    ! minimum moisture during dormance
    REAL,DIMENSION(npts,nvm),INTENT(in) :: hum_min_dormance
    ! fraction of litter above the ground belonging to different PFTs
    REAL,DIMENSION(npts,nlitt,nvm,nlevs,nelements),INTENT(in) :: litter
    ! dead leaves on ground, per PFT, metabolic and structural,
    ! in gC/(m**2 of ground)
    REAL,DIMENSION(npts,nvm,nlitt),INTENT(in) :: dead_leaves
    ! SOM pool: active, slow, or passive, (gC(orN)/m**2)
    REAL,DIMENSION(npts,ncarb,nvm,nelements),INTENT(in) :: som
   ! ratio Lignine/Carbon in structural litter, above and below ground, (gC/m**2)
    REAL,DIMENSION(npts,nvm,nlevs),INTENT(in) :: lignin_struc
    ! turnover_time of leaves
    ! ratio Lignine/Carbon in woody litter, above and below ground, (gC/m**2)
    REAL,DIMENSION(npts,nvm,nlevs),INTENT(in) :: lignin_wood
    ! turnover_time of leaves

    REAL,DIMENSION(npts,nvm),INTENT(in) :: turnover_time

    ! For Spinup matrix resolution
    INTEGER, INTENT(in) :: global_years   
    LOGICAL, DIMENSION(npts), INTENT(in) :: ok_equilibrium
 !DSG   LOGICAL, INTENT(in) :: ok_spunup ! anal_spin
    REAL, DIMENSION(npts), INTENT(in) :: nbp_accu  !! Accumulated Net Biospheric Production over the year 
    REAL, DIMENSION(npts), INTENT(in) :: nbp_flux  !! Net Biospheric Production over the year 
    !-
    REAL, DIMENSION(npts,nvm,nbpools,nbpools), INTENT(in) :: MatrixV
    REAL, DIMENSION(npts,nvm,nbpools), INTENT(in) :: VectorU
    REAL, DIMENSION(npts,nvm,nbpools), INTENT(in) :: previous_stock
    REAL, DIMENSION(npts,nvm,nbpools), INTENT(in) :: current_stock
    REAL, DIMENSION(npts,nvm,nbpools), INTENT(in) :: CN_som_litter_longterm
    REAL, INTENT(in)                              :: tau_CN_longterm  !! Counter used for calculating the longterm CN ratio of SOM and litter pools (seconds)
    REAL, DIMENSION(npts,nvm,npco2),   INTENT(in) :: assim_param

    REAL, DIMENSION(npts,nvm), INTENT(in)         :: record_reserve_grass !! record maximum grass reserve pool

    REAL, DIMENSION(npts,nvm), INTENT(in)         :: KF                !! Scaling factor to convert sapwood mass
                                                                              !! into leaf mass (m)
    REAL, DIMENSION(npts,nvm), INTENT(in)         :: k_latosa_adapt    !! Leaf to sapwood area adapted for water  
                                                                              !! stress. Adaptation takes place at the  
                                                                              !! end of the year (m)
    REAL, DIMENSION(npts,nvm), INTENT(in)         :: rue_longterm      !! Longterm radiation use efficiency

    REAL, DIMENSION(npts,nvm), INTENT(in)         :: cn_leaf_avg_season    !! Seasonal CN ratio of leaves 
    REAL, DIMENSION(npts,nvm), INTENT(in)         :: nstress_season    !! N-related seasonal stress (used for allocation) 
    REAL, DIMENSION(npts,nvm), INTENT(in)         :: pstress_season    !! P-related seasonal stress (used for allocation) 
    
    REAL, DIMENSION(npts,nvm,nnspec), INTENT(in)  :: soil_n_min        !! mineral nitrogen in the soil (gN/m**2)  
                                                                              !! (first index=npts, second index=nvm, third index=nnspec)   
    REAL, DIMENSION(npts,nvm), INTENT(in)        :: p_O2               !! partial pressure of oxigen in the soil (hPa)
                                                                              !! (first index=npts, second index=nvm)
                                                        
    REAL, DIMENSION(npts,nvm), INTENT(in)        :: bact               !! denitrifier biomass (gC/m**2)
                                                                                                   !! (first index=npts, second index=nvm)
    REAL, DIMENSION(npts,nvm,nbpools), INTENT(in) :: CP_som_litter_longterm
    REAL, DIMENSION(npts,nvm), INTENT(in)         :: np_leaf_avg_season     !! Seasonal NP ratio of leaves 
    REAL, DIMENSION(npts,nvm,npspec), INTENT(in)  :: soil_p_min             !! mineral phosphorus in the soil (gP/m**2)  
                                                                                   !! (first index=npts, second index=nvm, third index=npspec)   
    REAL, DIMENSION(npts,nvm), INTENT(in)         :: f_Pdissolved            !! fraction of dissolved P in contact with roots [ ]
                                                                                    !! (first index=npts, second index=nvm ) 
    REAL,DIMENSION(npts,nvm,nminerals,2),INTENT(in) :: EW_grain             !!!! Mass and grains size of crushed rocks  
    REAL, DIMENSION(npts,nvm),INTENT(in)   :: sla_calc      !! leaf age-related SLA
    REAL, DIMENSION(npts,nvm),INTENT(in)   :: wshtotsum     !! accumulated harvested grass biomass
    REAL, DIMENSION(npts,nvm),INTENT(in)   :: sr_ugb        !! grazing stocking rate of grazing livestock
    REAL, DIMENSION(npts,nvm),INTENT(in)   :: nb_ani        !! optimal livestock density
    REAL, DIMENSION(npts,nvm),INTENT(in)   :: grazed_frac   !! optimal fraction for grazing (in contract to mowning)
    REAL, DIMENSION(npts,nvm), INTENT(in)  :: import_yield  !! annual total harvested grass biomass yield
    REAL, DIMENSION(npts),INTENT(in)       :: t2m_14        !! ''14 days'' 2 meter air temperature
    REAL, DIMENSION(npts,nlitt,nvm), INTENT(in)  :: litter_not_avail !! litter not edible for animals
    REAL, DIMENSION(npts,nvm), INTENT(in)  :: nb_grazingdays!!annual total number of days animal grazing in field
    REAL, DIMENSION(npts),INTENT(in)       :: after_snow    !! day counter after snow melt to prevent grazing
    REAL, DIMENSION(npts),INTENT(in)       :: after_wet     !! day counter after wet soil is detected to prevent grazing
    REAL, DIMENSION(npts),INTENT(in)       :: wet1day       !! accumulated days with wet soil
    REAL, DIMENSION(npts),INTENT(in)       :: wet2day       !! accumulated days with wet soil
    REAL, DIMENSION(npts,nvm), INTENT(in)  :: GRM_devstage  !! development stage of grasses
    REAL, DIMENSION(npts,nvm), INTENT(in)  :: days_senescence!! day counter after senescence is detected

    !-
    ! 0.2 local
    !-
    ! date, real
    REAL :: date_real
    ! PFT exists (equivalent to fpc_max > 0 for natural PFTs), real
    REAL,DIMENSION(npts,nvm) :: PFTpresent_real
    ! is the plant senescent ?
    ! (only for deciduous trees - carbohydrate reserve), real
    REAL,DIMENSION(npts,nvm) :: senescence_real
    REAL,DIMENSION(npts,nvm) :: begin_leaves_real

    ! in order for this PFT to be introduced,
    ! does it have to be present in an adjacent grid box? - real
    REAL,DIMENSION(npts,nvm) :: need_adjacent_real
    ! To store variables names for I/O
    CHARACTER(LEN=80) :: var_name
    ! string suffix indicating an index
    CHARACTER(LEN=10) :: part_str
    ! string suffix indicating litter type
    CHARACTER(LEN=4),DIMENSION(nlitt) :: litter_str
    ! string suffix indicating level
    CHARACTER(LEN=2),DIMENSION(nlevs) :: level_str
!JC separate tissue age
    CHARACTER(LEN=4),DIMENSION(nleafages) :: age_str
    ! temporary storage
    REAL,DIMENSION(1) :: xtmp
    REAL, DIMENSION(1) :: vartmp  !! temporary variable because restget/restput needs a variable with DIMESION(:)
    ! index
    INTEGER :: j,k,l,m
    CHARACTER(LEN=2),DIMENSION(nelements) :: element_str  !! string suffix indicating element
    REAL, DIMENSION(1) :: temp_global_years
    CHARACTER(LEN=6),DIMENSION(nbpools) :: pools_str
    REAL, DIMENSION(npts) :: ok_equilibrium_real    

    ! land cover change variables 
    ! products remaining in the 10/100 year-turnover pool after the annual release for each compartment
    ! (10 or 100 + 1 : input from year of land cover change)
    REAL,DIMENSION(npts,0:10),INTENT(in)                           :: prod10
    REAL,DIMENSION(npts,0:100),INTENT(in)                          :: prod100
    ! annual release from the 10/100 year-turnover pool compartments
    REAL,DIMENSION(npts,10),INTENT(in)                           :: flux10
    REAL,DIMENSION(npts,100),INTENT(in)                          :: flux100
    REAL, DIMENSION(npts), INTENT(in)                            :: convflux
    REAL, DIMENSION(npts), INTENT(in)                            :: cflux_prod10
    REAL, DIMENSION(npts), INTENT(in)                            :: cflux_prod100

    ! wood harvest variables 
    ! products remaining in the 10/100 year-turnover pool after the annual release for each compartment
    ! (10 or 100 + 1 : input from year of land cover change)
    REAL,DIMENSION(npts,0:10),INTENT(in)                           :: prod10_harvest
    REAL,DIMENSION(npts,0:100),INTENT(in)                          :: prod100_harvest
    ! annual release from the 10/100 year-turnover pool compartments
    REAL,DIMENSION(npts,10),INTENT(in)                           :: flux10_harvest
    REAL,DIMENSION(npts,100),INTENT(in)                          :: flux100_harvest
    REAL, DIMENSION(npts), INTENT(in)                            :: convflux_harvest
    REAL, DIMENSION(npts), INTENT(in)                            :: cflux_prod10_harvest
    REAL, DIMENSION(npts), INTENT(in)                            :: cflux_prod100_harvest
    REAL, DIMENSION(npts,nvm), INTENT(in)                        :: woodharvestpft
    REAL,DIMENSION(npts,nvm,nparts,nelements),INTENT(in)         :: bm_to_litter
    REAL,DIMENSION(npts),INTENT(in)                              :: carb_mass_total
    !---------------------------------------------------------------------
    IF (printlev >= 3) WRITE(numout,*) 'Entering writerestart'

    !-
    ! 1 string definitions
    !-
    DO l=1,nlitt
       IF     (l == imetabolic) THEN
          litter_str(l) = 'met'
       ELSEIF (l == istructural) THEN
          litter_str(l) = 'str'
       ELSEIF (l == iwoody) THEN
          litter_str(l) = 'wood'
       ELSE
          CALL ipslerr_p(3,'stomate_io writerestart','Define litter_str','','')
       ENDIF
    ENDDO
    !-
    DO l=1,nlevs
       IF     (l == iabove) THEN
          level_str(l) = 'ab'
       ELSEIF (l == ibelow) THEN
          level_str(l) = 'be'
       ELSE
          CALL ipslerr_p(3,'stomate_io writerestart','Define level_str','','')
       ENDIF
    ENDDO
    !JC separate tissue age
    !-
    DO l=1,nleafages
       IF     (l == 1) THEN
          age_str(l) = 'age1'
       ELSEIF (l == 2) THEN
          age_str(l) = 'age2'
       ELSEIF (l == 3) THEN
          age_str(l) = 'age3'
       ELSEIF (l == 4) THEN
          age_str(l) = 'age4'
       ELSE
          CALL ipslerr_p(3,'stomate_io readstart','Define age_str','','')
       ENDIF
    ENDDO
    !-
    DO l=1,nelements
       IF     (l == icarbon) THEN
          element_str(l) = ''
       ELSEIF (l == initrogen) THEN
          element_str(l) = '_n'
       ELSEIF (l == iphosphorus) THEN
          element_str(l) = '_p'
       ELSE
          CALL ipslerr_p(3,'stomate_io writerestart','Define element_str','','')
       ENDIF
    ENDDO
    !-
    pools_str(1:nbpools) =(/'str_ab ','str_be ','met_ab ','met_be ','wood_ab','wood_be', & 
         & 'actif  ','slow   ','passif ','surface'/) 

    !-
    IF (is_root_prc) THEN
       CALL ioconf_setatt_p ('UNITS','-')
       CALL ioconf_setatt_p ('LONG_NAME',' ')
    ENDIF
    !-
    ! 2 run control
    !-
    ! 2.2 time step of STOMATE in days
    !-
    IF (is_root_prc) THEN
       var_name = 'dt_days'
       xtmp(1) = dt_days
       CALL restput (rest_id_stomate, var_name, 1, 1, 1, itime, xtmp)
    ENDIF
    !-
    ! 2.3 date
    !-
    IF (is_root_prc) THEN
       var_name = 'date'
       date_real = REAL(date_loc,r_std)
       xtmp(1) = date_real
       CALL restput (rest_id_stomate, var_name, 1, 1, 1, itime, xtmp)
    ENDIF
    !-
    ! 3 daily meteorological variables
    !-
    var_name = 'moiavail_daily'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                moiavail_daily, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'gdd_init_date'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo,    2, 1, itime, &
         &              gdd_init_date, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'litterhum_daily'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo,    1, 1, itime, &
         &                litterhum_daily, 'scatter', nbp_glo, index_g)
    !-
    var_name = 't2m_daily'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo,    1, 1, itime, &
         &                t2m_daily, 'scatter', nbp_glo, index_g)
    !-
    var_name = 't2m_min_daily'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo,    1, 1, itime, &
         &                t2m_min_daily, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'tsurf_daily'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo,    1, 1, itime, &
         &                tsurf_daily, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'tsoil_daily'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nslm, 1, itime, &
         &                tsoil_daily, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'soilhum_daily'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nslm, 1, itime, &
         &                soilhum_daily, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'precip_daily'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo,    1, 1, itime, &
         &                precip_daily, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'max_eau_var_daily'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, 1, 1, itime, &
         &                max_eau_var_daily, 'scatter', nbp_glo, index_g)
    !-
    !-
    ! 4 productivities
    !-
    var_name = 'gpp_daily'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                gpp_daily, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'npp_daily'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                npp_daily, 'scatter', nbp_glo, index_g)
    !-
    DO l = 1,nelements
       DO k = 1,nparts
          WRITE(part_str,'(I2)') k
          IF (k < 10) part_str(1:1) = '0'
          var_name = 'turnover_daily_'//part_str(1:LEN_TRIM(part_str))//element_str(l)
          CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
               &                   turnover_daily(:,:,k,l), 'scatter', nbp_glo, index_g)
       ENDDO
    END DO
    !-
    ! 5 monthly meteorological variables
    !-
    var_name = 'moiavail_month'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                moiavail_month, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'moiavail_week'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                moiavail_week, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'moiavail_grow'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                moiavail_growingseason, 'scatter', nbp_glo, index_g)
    !-
    var_name = 't2m_longterm'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo,    1, 1, itime, &
         &                t2m_longterm, 'scatter', nbp_glo, index_g)
    
    IF (is_root_prc) THEN
       var_name='tau_longterm'
       vartmp(1)=tau_longterm
       CALL restput (rest_id_stomate, var_name, 1, 1, 1, itime, vartmp)
    ENDIF
       

    var_name = 't2m_month'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo,    1, 1, itime, &
                         t2m_month, 'scatter', nbp_glo, index_g)
    

    CALL restput_p (rest_id_stomate, 'Tseason', nbp_glo,    1, 1, itime, &
         Tseason, 'scatter', nbp_glo, index_g)
    
    CALL restput_p (rest_id_stomate, 'Tseason_length', nbp_glo,    1, 1, itime, &
         Tseason_length, 'scatter', nbp_glo, index_g)
    
    CALL restput_p (rest_id_stomate, 'Tseason_tmp', nbp_glo,    1, 1, itime, &
         Tseason_tmp, 'scatter', nbp_glo, index_g)
    
    CALL restput_p (rest_id_stomate, 'Tmin_spring_time', nbp_glo, nvm, 1, itime, &
         Tmin_spring_time, 'scatter', nbp_glo, index_g)
    
    CALL restput_p (rest_id_stomate, 'onset_date', nbp_glo, nvm, 1, itime, &
         onset_date(:,:), 'scatter', nbp_glo, index_g)
    
    var_name = 't2m_week'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo,    1, 1, itime, &
         &                t2m_week, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'tsoil_month'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nslm, 1, itime, &
         &                tsoil_month, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'soilhum_month'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nslm, 1, itime, &
         &                soilhum_month, 'scatter', nbp_glo, index_g)
    !-
    ! 6 fire probability
    !-
    var_name = 'fireindex'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                fireindex, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'firelitter'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                firelitter, 'scatter', nbp_glo, index_g)
    !-
    ! 7 maximum and minimum moisture availabilities for tropic phenology
    !-
    var_name = 'maxmoistr_last'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                maxmoiavail_lastyear, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'maxmoistr_this'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                maxmoiavail_thisyear, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'minmoistr_last'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                minmoiavail_lastyear, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'minmoistr_this'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                minmoiavail_thisyear, 'scatter', nbp_glo, index_g)

    !-
    ! X long term soil water loss fluxes
    !-

    var_name = 'drainage_longterm'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo,   1, 1, itime, &
         &                drainage_longterm, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'runoff_longterm'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo,   1, 1, itime, &
         &                runoff_longterm, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'lai_target_longterm'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                lai_target_longterm, 'scatter', nbp_glo, index_g)
    !-
    ! 8 maximum "weekly" GPP
    !-
    var_name = 'maxgppweek_lastyear'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                maxgppweek_lastyear, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'maxgppweek_thisyear'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                maxgppweek_thisyear, 'scatter', nbp_glo, index_g)
    !-
    ! 9 annual GDD0
    !-
    var_name = 'gdd0_thisyear'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, 1, 1, itime, &
         &                gdd0_thisyear, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'gdd0_lastyear'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, 1, 1, itime, &
         &                gdd0_lastyear, 'scatter', nbp_glo, index_g)
    !-
    ! 10 annual precipitation
    !-
    var_name = 'precip_thisyear'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, 1, 1, itime, &
         &                precip_thisyear, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'precip_lastyear'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, 1, 1, itime, &
         &                precip_lastyear, 'scatter', nbp_glo, index_g)
    !-
    ! 11 derived "biometeorological" variables
    !-
    var_name = 'gdd_m5_dormance'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                gdd_m5_dormance, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'gdd_from_growthinit'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &              gdd_from_growthinit, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'gdd_midwinter'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                gdd_midwinter, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'ncd_dormance'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                ncd_dormance, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'ngd_minus5'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                ngd_minus5, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'time_hum_min'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                time_hum_min, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'hum_min_dormance'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                hum_min_dormance, 'scatter', nbp_glo, index_g)
    !-
    ! 12 Plant status
    !-
    var_name = 'PFTpresent'
    WHERE ( PFTpresent(:,:) )
       PFTpresent_real = un
    ELSEWHERE
       PFTpresent_real = zero
    ENDWHERE
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                PFTpresent_real, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'ind'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                ind, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'turnover_time'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                turnover_time, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'adapted'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                adapted, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'regenerate'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                regenerate, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'npp_longterm'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                npp_longterm, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'lm_lastyearmax'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                lm_lastyearmax, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'lm_thisyearmax'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                lm_thisyearmax, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'maxfpc_lastyear'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                maxfpc_lastyear, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'maxfpc_thisyear'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                maxfpc_thisyear, 'scatter', nbp_glo, index_g)
    !-
    DO l = 1,nelements
       DO k = 1,nparts
          WRITE(part_str,'(I2)') k
          IF (k < 10) part_str(1:1) = '0'
          var_name = 'turnover_loterm_'//part_str(1:LEN_TRIM(part_str))//element_str(l)
          CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
               &                   turnover_longterm(:,:,k,l), 'scatter', nbp_glo, index_g)
       ENDDO
    END DO
    !-
    var_name = 'gpp_week'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                gpp_week, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'npp_week'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                npp_week, 'scatter', nbp_glo, index_g)
    !-
    DO l = 1,nelements
       DO k = 1,nparts
          WRITE(part_str,'(I2)') k
          IF (k < 10) part_str(1:1) = '0'
          var_name = 'biomass_'//part_str(1:LEN_TRIM(part_str))//element_str(l)
          CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
               &                   biomass(:,:,k,l), 'scatter', nbp_glo, index_g)
       ENDDO
    END DO

    !-
    DO l = 1,nminerals
       DO k = 1,2
          WRITE(part_str,'(I2)') k
          IF (k < 10) part_str(1:1) = '0'
          var_name = 'EW_grain_'//part_str(1:LEN_TRIM(part_str))//element_str(l)
          CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
               &                   EW_grain(:,:,l,k), 'scatter', nbp_glo, index_g)
       ENDDO
    END DO

    !-
    DO k=1,nparts
       WRITE(part_str,'(I2)') k
       IF (k < 10) part_str(1:1) = '0'
       var_name = 'maint_resp_'//part_str(1:LEN_TRIM(part_str))
       CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
            &                   resp_maint_part(:,:,k), 'scatter', nbp_glo, index_g)
    ENDDO
    !JC separate tissue age
    !! old leaf_age leaf_frac
    !-
    !    DO m=1,nleafages
    !       WRITE(part_str,'(I2)') m
    !       IF (m < 10) part_str(1:1) = '0'
    !       var_name = 'leaf_age_'//part_str(1:LEN_TRIM(part_str))
    !       CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
    !            &                  leaf_age(:,:,m), 'scatter', nbp_glo, index_g)
    !    ENDDO
    !    !-
    !    DO m=1,nleafages
    !       WRITE(part_str,'(I2)') m
    !       IF (m < 10) part_str(1:1) = '0'
    !       var_name = 'leaf_frac_'//part_str(1:LEN_TRIM(part_str))
    !       CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
    !            &                   leaf_frac(:,:,m), 'scatter', nbp_glo, index_g)
    !    ENDDO
    !-
    DO l = 1,nleafages
       DO k = 1,nparts
          WRITE(part_str,'(I2)') k
          IF (k < 10) part_str(1:1) = '0'
          var_name = 'leaf_age_'//part_str(1:LEN_TRIM(part_str))//'_'//TRIM(age_str(l))
          CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
               &                   leaf_age(:,:,k,l), 'scatter', nbp_glo, index_g)
       ENDDO
    END DO
    !-
    DO l = 1,nleafages
       DO k = 1,nparts
          WRITE(part_str,'(I2)') k
          IF (k < 10) part_str(1:1) = '0'
          var_name = 'leaf_frac_'//part_str(1:LEN_TRIM(part_str))//'_'//TRIM(age_str(l))
          CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
               &                   leaf_frac(:,:,k,l), 'scatter', nbp_glo, index_g)
       ENDDO
    END DO
    !End JC separate tissue age
    !-
    var_name = 'senescence'
    WHERE ( senescence(:,:) )
       senescence_real = un
    ELSEWHERE
       senescence_real = zero
    ENDWHERE
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                senescence_real, 'scatter', nbp_glo, index_g)
 
    ! Transform the logical variable begin_leaves to real before writing to restart file
    WHERE ( begin_leaves(:,:) )
       begin_leaves_real = un
    ELSEWHERE
       begin_leaves_real = zero
    ENDWHERE
    CALL restput_p (rest_id_stomate, 'begin_leaves', nbp_glo, nvm, 1, itime, &
         begin_leaves_real, 'scatter', nbp_glo, index_g)


    var_name = 'when_growthinit'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                when_growthinit, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'age'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &                age, 'scatter', nbp_glo, index_g)
    !-
    ! 13 CO2
    !-
    var_name = 'resp_hetero'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                resp_hetero, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'resp_maint'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                resp_maint, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'resp_growth'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                resp_growth, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'co2_fire'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo,  nvm, 1, itime, &
         &                co2_fire, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'co2_to_bm_dgvm'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                co2_to_bm_dgvm, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'n_to_bm'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                n_to_bm, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'p_to_bm'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                p_to_bm, 'scatter', nbp_glo, index_g)
    !-
    ! 14 vegetation distribution after last light competition
    !-
    var_name = 'veget_lastlight'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                veget_lastlight, 'scatter', nbp_glo, index_g)
    !-
    ! 15 establishment criteria
    !-
    var_name = 'everywhere'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                everywhere, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'need_adjacent'
    WHERE (need_adjacent(:,:))
       need_adjacent_real = un
    ELSEWHERE
       need_adjacent_real = zero
    ENDWHERE
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                need_adjacent_real, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'RIP_time'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                RIP_time, 'scatter', nbp_glo, index_g)
    !-
    ! 17 litter
    !-
    DO k = 1,nelements
       DO l = 1,nlevs
          DO m = 1,nvm
             WRITE (part_str, '(I2)') m
             IF (m<10) part_str(1:1)='0'
             var_name = 'litter_'//part_str(1:LEN_TRIM(part_str))//'_'//level_str(l)//element_str(k)
             CALL restput_p (rest_id_stomate, var_name, nbp_glo, nlitt, 1, itime, &
                  &                     litter(:,:,m,l,k), 'scatter', nbp_glo, index_g)
          ENDDO
       ENDDO
    END DO
    !-
    DO l=1,nlitt
       var_name = 'dead_leaves_'//TRIM(litter_str(l))
       CALL restput_p (rest_id_stomate, var_name, nbp_glo,  nvm, 1, itime, &
            &                   dead_leaves(:,:,l), 'scatter', nbp_glo, index_g)
    ENDDO
    !-
    DO m=1,nvm
       WRITE (part_str, '(I2)') m
       IF (m<10) part_str(1:1)='0'
       var_name = 'carbon_'//part_str(1:LEN_TRIM(part_str))
       CALL restput_p (rest_id_stomate, var_name, nbp_glo, ncarb, 1, itime, &
            &                   som(:,:,m,icarbon), 'scatter', nbp_glo, index_g)
    ENDDO
    !-
    DO m=1,nvm 
       WRITE (part_str, '(I2)') m 
       IF (m<10) part_str(1:1)='0' 
       var_name = 'nitrogen_'//part_str(1:LEN_TRIM(part_str)) 
       CALL restput_p (rest_id_stomate, var_name, nbp_glo, ncarb, 1, itime, & 
            &                   som(:,:,m,initrogen), 'scatter', nbp_glo, index_g) 
    ENDDO
    !-
    DO m=1,nvm 
       WRITE (part_str, '(I2)') m 
       IF (m<10) part_str(1:1)='0' 
       var_name = 'phosphorus_'//part_str(1:LEN_TRIM(part_str)) 
       CALL restput_p (rest_id_stomate, var_name, nbp_glo, ncarb, 1, itime, & 
            &                   som(:,:,m,iphosphorus), 'scatter', nbp_glo, index_g) 
    ENDDO
    !- 
    DO l=1,nlevs
       var_name = 'lignin_struc_'//level_str(l)
       CALL restput_p &
            &      (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
            &       lignin_struc(:,:,l), 'scatter', nbp_glo, index_g)
    ENDDO
    !- 
    DO l=1,nlevs
       var_name = 'lignin_wood_'//level_str(l)
       CALL restput_p &
            &      (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
            &       lignin_wood(:,:,l), 'scatter', nbp_glo, index_g)
    ENDDO
    !-
    ! 18 land cover change
    !-
    var_name = 'prod10'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, 11, 1, itime, &
         &                prod10, 'scatter', nbp_glo, index_g)
    var_name = 'prod100'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, 101, 1, itime, &
         &                prod100, 'scatter', nbp_glo, index_g)
    var_name = 'flux10'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, 10, 1, itime, &
         &                flux10, 'scatter', nbp_glo, index_g)
    var_name = 'flux100'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, 100, 1, itime, &
         &                flux100, 'scatter', nbp_glo, index_g)

    var_name = 'convflux'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, 1, 1, itime, &
         &              convflux, 'scatter', nbp_glo, index_g)
    var_name = 'cflux_prod10'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, 1, 1, itime, &
         &              cflux_prod10, 'scatter', nbp_glo, index_g)
    var_name = 'cflux_prod100'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, 1, 1, itime, &
         &              cflux_prod100, 'scatter', nbp_glo, index_g)

    !-
    ! 18-bis wood harvest
    !-
    IF (do_wood_harvest) THEN
       var_name = 'prod10_harvest'
       CALL restput_p (rest_id_stomate, var_name, nbp_glo, 11, 1, itime, &
            prod10_harvest, 'scatter', nbp_glo, index_g)
       var_name = 'prod100_harvest'
       CALL restput_p (rest_id_stomate, var_name, nbp_glo, 101, 1, itime, &
            prod100_harvest, 'scatter', nbp_glo, index_g)
       var_name = 'flux10_harvest'
       CALL restput_p (rest_id_stomate, var_name, nbp_glo, 10, 1, itime, &
            flux10_harvest, 'scatter', nbp_glo, index_g)
       var_name = 'flux100_harvest'
       CALL restput_p (rest_id_stomate, var_name, nbp_glo, 100, 1, itime, &
            flux100_harvest, 'scatter', nbp_glo, index_g)
       var_name = 'convflux_harvest'
       CALL restput_p (rest_id_stomate, var_name, nbp_glo, 1, 1, itime, &
            convflux_harvest, 'scatter', nbp_glo, index_g)
       var_name = 'cflux_prod10_harvest'
       CALL restput_p (rest_id_stomate, var_name, nbp_glo, 1, 1, itime, &
            cflux_prod10_harvest, 'scatter', nbp_glo, index_g)
       var_name = 'cfluxprod100_harvest'
       CALL restput_p (rest_id_stomate, var_name, nbp_glo, 1, 1, itime, &
            cflux_prod100_harvest, 'scatter', nbp_glo, index_g)
       var_name = 'woodharvestpft'
       CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
            woodharvestpft, 'scatter', nbp_glo, index_g)
    END IF

    DO l = 1,nelements
       DO k = 1,nparts
          WRITE(part_str,'(I2)') k
          IF (k < 10) part_str(1:1) = '0'
          var_name = 'bm_to_litter_'//part_str(1:LEN_TRIM(part_str))//element_str(l)
          CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
               &                bm_to_litter(:,:,k,l), 'scatter', nbp_glo, index_g)
       ENDDO
    END DO

    var_name = 'carb_mass_total'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, 1, 1, itime, &
         &              carb_mass_total, 'scatter', nbp_glo, index_g)
    !gmjc
      var_name = 'wshtotsum'
      CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
     &              wshtotsum, 'scatter', nbp_glo, index_g)
    !-
      var_name = 'sr_ugb'
      CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
     &              sr_ugb, 'scatter', nbp_glo, index_g)
    !-
      var_name = 'sla_calc'
      CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
     &              sla_calc, 'scatter', nbp_glo, index_g)
    !-
      var_name = 'nb_ani'
      CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
     &              nb_ani, 'scatter', nbp_glo, index_g)
    !-
      var_name = 'grazed_frac'
      CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
     &              grazed_frac, 'scatter', nbp_glo, index_g)
    !-
      var_name = 'import_yield'
      CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
     &              import_yield, 'scatter', nbp_glo, index_g)
    !-
      var_name = 't2m_14'
      CALL restput_p (rest_id_stomate, var_name, nbp_glo,    1, 1, itime, &
     &              t2m_14, 'scatter', nbp_glo, index_g)
    
    ! in current trunk version there is no 3D restart 
    ! only str and met not wood
        DO l=1,2 !nlitt
           var_name = 'litter_not_avail_'//TRIM(litter_str(l))
           CALL restput_p (rest_id_stomate, var_name, nbp_glo,  nvm, 1, itime, &
                &      litter_not_avail(:,l,:), 'scatter', nbp_glo,index_g)
        ENDDO
    !-
    !  CALL restput_p (rest_id_stomate, 'litter_not_avail', nbp_glo,  nlitt, nvm, itime, &
    ! &                   litter_not_avail, 'scatter', nbp_glo, index_g)
    !-
      var_name = 'nb_grazingdays'
      CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
     &              nb_grazingdays, 'scatter', nbp_glo, index_g)
    !-
    !-
      var_name = 'after_snow'
      CALL restput_p (rest_id_stomate, var_name, nbp_glo,    1, 1, itime, &
     &              after_snow, 'scatter', nbp_glo, index_g)
    !-
      var_name = 'after_wet'
      CALL restput_p (rest_id_stomate, var_name, nbp_glo,    1, 1, itime, &
     &              after_wet, 'scatter', nbp_glo, index_g)
    !-
      var_name = 'wet1day'
      CALL restput_p (rest_id_stomate, var_name, nbp_glo,    1, 1, itime, &
     &              wet1day, 'scatter', nbp_glo, index_g)
    !-
      var_name = 'wet2day'
      CALL restput_p (rest_id_stomate, var_name, nbp_glo,    1, 1, itime, &
     &              wet2day, 'scatter', nbp_glo, index_g)
    !-
      var_name = 'GRM_devstage'
      CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
     &              GRM_devstage, 'scatter', nbp_glo, index_g)
      var_name = 'days_senescence'
      CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
     &              days_senescence, 'scatter', nbp_glo, index_g)
    !end gmjc

    !-
    ! 19. Spinup
    !-
    IF (spinup_analytic) THEN

       IF (is_root_prc) THEN
          temp_global_years(1) = REAL(global_years)
          var_name='Global_years'
          CALL restput (rest_id_stomate, var_name, 1, 1, 1, itime, temp_global_years)
       ENDIF
       
       var_name = 'nbp_sum'
       CALL restput_p (rest_id_stomate, var_name, nbp_glo, 1, 1, itime, &
            &              nbp_accu, 'scatter', nbp_glo, index_g)

       var_name = 'nbp_flux'
       CALL restput_p (rest_id_stomate, var_name, nbp_glo, 1, 1, itime, &
            &              nbp_flux, 'scatter', nbp_glo, index_g)

       var_name = 'ok_equilibrium'
       WHERE(ok_equilibrium(:))
          ok_equilibrium_real = un
       ELSEWHERE
          ok_equilibrium_real = zero
       ENDWHERE
       CALL restput_p (rest_id_stomate, var_name, nbp_glo, 1, 1, itime, &
            &               ok_equilibrium_real, 'scatter', nbp_glo, index_g)
       
       DO k = 1,nbpools
          DO j = 1,nbpools
             WRITE(part_str,'(I2)') k
             IF (k < 10) part_str(1:1) = '0'             
             var_name = 'MatrixV_'//part_str(1:LEN_TRIM(part_str))//'_'//TRIM(pools_str(j))
             CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
                  &                MatrixV(:,:,k,j), 'scatter', nbp_glo, index_g)
          ENDDO
       ENDDO
          
       DO k = 1,nbpools
          WRITE(part_str,'(I2)') k
          IF (k < 10) part_str(1:1) = '0' 
          var_name = 'Vector_U_'//part_str(1:LEN_TRIM(part_str))
          CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
               &                VectorU(:,:,k), 'scatter', nbp_glo, index_g)
       ENDDO
          
       DO k = 1,nbpools
          WRITE(part_str,'(I2)') k
          IF (k < 10) part_str(1:1) = '0' 
          var_name = 'previous_stock_'//part_str(1:LEN_TRIM(part_str))
          CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
               &                previous_stock(:,:,k), 'scatter', nbp_glo, index_g)
       ENDDO
          
       DO k = 1,nbpools
          WRITE(part_str,'(I2)') k
          IF (k < 10) part_str(1:1) = '0' 
          var_name = 'current_stock_'//part_str(1:LEN_TRIM(part_str))
          CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
               &                current_stock(:,:,k), 'scatter', nbp_glo, index_g)
       ENDDO

       DO k = 1,nbpools
          WRITE(part_str,'(I2)') k
          IF (k < 10) part_str(1:1) = '0' 
          var_name = 'CN_longterm_'//part_str(1:LEN_TRIM(part_str))
          CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
               &                CN_som_litter_longterm(:,:,k), 'scatter', nbp_glo, index_g)
       ENDDO

       IF (is_root_prc) THEN
          var_name='tau_CN_longterm'
          vartmp(1)=tau_CN_longterm
          CALL restput (rest_id_stomate, var_name, 1, 1, 1, itime, vartmp)
       ENDIF

       DO k = 1,nbpools
          WRITE(part_str,'(I2)') k
          IF (k < 10) part_str(1:1) = '0' 
          var_name = 'CP_longterm_'//part_str(1:LEN_TRIM(part_str))
          CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
               &                CP_som_litter_longterm(:,:,k), 'scatter', nbp_glo, index_g)
       ENDDO


            

    ENDIF !(spinup_analytic)

    var_name = 'record_reserve_grass'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
          &              record_reserve_grass(:,:), 'scatter', nbp_glo, index_g)

    !-
    var_name = 'KF'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
          &              KF(:,:), 'scatter', nbp_glo, index_g)
    !-
    var_name = 'k_latosa_adapt'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
          &              k_latosa_adapt(:,:), 'scatter', nbp_glo, index_g)
    !-

    var_name = 'rue_longterm'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
          &              rue_longterm(:,:), 'scatter', nbp_glo, index_g)
    !-
    var_name = 'cn_leaf_avg_season' 
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, & 
         &              cn_leaf_avg_season(:,:), 'scatter', nbp_glo, index_g) 
    !- 
    var_name = 'nstress_season' 
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, & 
         &              nstress_season(:,:), 'scatter', nbp_glo, index_g) 
    !- 
    var_name = 'pstress_season' 
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, & 
         &              pstress_season(:,:), 'scatter', nbp_glo, index_g) 


    !- 
    DO k=1,nnspec 
       WRITE(part_str,'(I1)') k 
       var_name = 'soil_n_min_'//part_str(1:LEN_TRIM(part_str)) 
       CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, & 
            &              soil_n_min(:,:,k), 'scatter', nbp_glo, index_g) 
    ENDDO
    !-
    var_name = 'np_leaf_avg_season' 
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, & 
         &              np_leaf_avg_season(:,:), 'scatter', nbp_glo, index_g) 
    !- 
    DO k=1,npspec 
       WRITE(part_str,'(I1)') k 
       var_name = 'soil_p_min_'//part_str(1:LEN_TRIM(part_str)) 
       CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, & 
            &              soil_p_min(:,:,k), 'scatter', nbp_glo, index_g) 
    ENDDO
    ! -
    var_name = 'f_Pdissolved' 
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, & 
         &              f_Pdissolved(:,:), 'scatter', nbp_glo, index_g) 
    !- 
    var_name = 'p_O2' 
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, & 
         &              p_O2(:,:), 'scatter', nbp_glo, index_g) 
    !- 
    !- 
    var_name = 'bact' 
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, & 
         &              bact(:,:), 'scatter', nbp_glo, index_g) 
    !- 
    !- 

    DO k = 1,npco2
       WRITE(part_str,'(I2)') k
       IF (k < 10) part_str(1:1) = '0' 
       var_name = 'assim_param_'//part_str(1:LEN_TRIM(part_str))
       CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
            &                assim_param(:,:,k), 'scatter', nbp_glo, index_g)
    ENDDO
       

    IF (printlev >= 4) WRITE(numout,*) 'Leaving writerestart'
    !--------------------------
  END SUBROUTINE writerestart
  !-
  !===
  !-
!! ================================================================================================================================
!! SUBROUTINE   : readep
!!
!>\BRIEF         
!!
!! DESCRIPTION  : None
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): deposition
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE readep (npts, lalo, resolution, deposition)

    IMPLICIT NONE

    !! 0. Parameters and variables declaration
    
    !! 0.1 Input variables  
    
    INTEGER, INTENT(in)                  :: npts        !! Domain size
    REAL, DIMENSION(npts,2), INTENT(in)  :: lalo        !! Geogr. coordinates (latitude,longitude) (degrees)
    REAL, DIMENSION(npts,2), INTENT(in)  :: resolution  !! size in x an y of the grid (m)

    !! 0.3 Modified variables

    REAL, DIMENSION(npts,3), INTENT(inout) :: deposition  !! N&P inputs by deposition [g/m2/day] (NOY=1,NHX=2,P=3)

!_ ================================================================================================================================


    IF ( ok_pcycle.OR.ok_ncycle ) THEN
       CALL get_deposition (npts, lalo, resolution, deposition)
    ENDIF

  END SUBROUTINE readep 

!=====================================================================================
!! SUBROUTINE   : readlitholgoy
!!
!>\BRIEF         
!!
!! DESCRIPTION  : reads in input fields of lithological and soil order characteristics
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): lith_frac, soil_shield, soil_orders
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
  SUBROUTINE readlithology (npts, lalo, resolution, soil_order_imposed, lith_frac, soil_shield, soil_orders)

    IMPLICIT NONE

    !! 0. Parameters and variables declaration
    
    !! 0.1 Input variables  
    
    INTEGER, INTENT(in)                  :: npts        !! Domain size
    REAL, DIMENSION(npts,2), INTENT(in)  :: lalo        !! Geogr. coordinates (latitude,longitude) (degrees)
    REAL, DIMENSION(npts,2), INTENT(in)  :: resolution  !! size in x an y of the grid (m)
    REAL, DIMENSION(npts)     , INTENT(in)   :: soil_order_imposed  !! soil orders we impose globally if impsoils=true

    !! 0.3 Modified variables

    REAL,    DIMENSION(npts,nlcm), INTENT(inout) :: lith_frac    !! area fraction for lithologies [0-1]
    REAL,    DIMENSION(npts)     , INTENT(inout) :: soil_shield  !! factor to correct P weathering for soil shielding effects [0-1]
    INTEGER, DIMENSION(npts)     , INTENT(inout) :: soil_orders  !! soil orders (USDA)

!_ ================================================================================================================================

    IF ( ok_pcycle) THEN
       CALL get_lithology (npts, lalo, resolution, lith_frac, soil_shield)
       CALL get_soil_orders (npts, lalo, resolution, soil_order_imposed,soil_orders)
    ENDIF

  END SUBROUTINE readlithology

!! ================================================================================================================================
!! SUBROUTINE   :  get_lithology_clear
!!
!>\BRIEF         
!!
!! DESCRIPTION  : None
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S):
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE get_lithology_clear

    firstcall_lith = .TRUE.
    IF (ALLOCATED (lith)) DEALLOCATE( lith )
    IF (ALLOCATED (soilshield)) DEALLOCATE( soilshield)

  END SUBROUTINE get_lithology_clear




!! ================================================================================================================================
!! SUBROUTINE   :  get_deposition_clear
!!
!>\BRIEF         
!!
!! DESCRIPTION  : None
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S):
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE get_deposition_clear

    firstcall_dep = .TRUE.
    IF (ALLOCATED (depo)) DEALLOCATE( depo )

  END SUBROUTINE get_deposition_clear
  SUBROUTINE get_lithology (npts, lalo, resolution, lith_out, soilshield_out)

    IMPLICIT NONE

    !! 0. Parameters and variables declaration
    
    !! 0.1 Input variables 

    INTEGER, INTENT(in)                 :: npts                 !! Domain size
    REAL, DIMENSION(npts,2), INTENT(in) :: lalo                 !! Geogr. coordinates (latitude,longitude) (degrees)
    REAL, DIMENSION(npts,2), INTENT(in) :: resolution           !! size in x an y of the grid (m)

    !! 0.2 Output variables

    REAL, DIMENSION(npts,nlcm), INTENT(out)  :: lith_out        !! fractions of gridarea underlied by respective lithological type

    REAL, DIMENSION(npts), INTENT(out)  :: soilshield_out       !! correction factor for P release for soil shielding (Hartmann et al. 2014)

    !! 0.4 Local variables

    INTEGER, PARAMETER                  :: nbvmax = 200         !!
    CHARACTER(LEN=80)                          :: filename             !!
    INTEGER                             :: iml, jml             !!
    INTEGER                             :: lml                  !!
    INTEGER                             :: tml                  !!
    INTEGER                             :: fid                  !! 
    INTEGER                             :: ib, ip, jp           !!
    INTEGER                             :: fopt                 !!
    INTEGER                             :: ilf                  !!
    INTEGER                             :: lastjp               !! 
    REAL                                :: date                 !!
    REAL                                :: dt                   !!
    REAL                                :: coslat               !!
    INTEGER                             :: itau(1)              !!
    REAL, ALLOCATABLE, DIMENSION(:)     :: lev                  !!
    REAL, ALLOCATABLE, DIMENSION(:,:)   :: lat_rel, lon_rel     !!   
    REAL, ALLOCATABLE, DIMENSION(:,:)   :: lat_ful, lon_ful     !! 
    !DSG: these two could be combined
    REAL, ALLOCATABLE, DIMENSION(:,:,:) :: lithology_file       !! 
    REAL, ALLOCATABLE, DIMENSION(:,:)   :: lithology_file2      !! 
    !DSG: these should be combined
    REAL, ALLOCATABLE, DIMENSION(:,:)   :: loup_rel, lolow_rel  !!
    REAL, ALLOCATABLE, DIMENSION(:,:)   :: laup_rel, lalow_rel  !!
    REAL                                :: lon_up, lon_low      !! 
    REAL                                :: lat_up, lat_low      !!
    REAL                                :: ax, ay, sgn          !!
    REAL, DIMENSION(nbvmax)             :: area                 !!
    !DSG: these two could be combined
    REAL, DIMENSION(nbvmax,nlcm)        :: tt                   !!
    REAL, DIMENSION(nbvmax)             :: tt2                  !!
    !DSG: these two could be combined
    REAL                                :: resx, resy           !!
    LOGICAL                                    :: do_again             !!


!_ ================================================================================================================================

    !! 1 If this is the first call, calculate the lithological fractions
    !!   and keep it in memory

    IF (printlev >= 4) WRITE(numout,*) 'Entering get_lithology'


    IF (firstcall_lith) THEN
       !---
       !-- 1.1 only do this once
       !---
       firstcall_lith = .FALSE.
       !---
       !-- 1.2 allocate the field
       !---
       ALLOCATE( lith(npts,nlcm) )
       ALLOCATE( soilshield(npts) )
       !---
       !-- 1.3 read and interpolate the lithology file
       !---
       !-- Needs to be a configurable variable
       !---
       !Config Key   = LITHOLOGY_FILE
       !Config Desc  = Name of file from which the lithology is read
       !Config If    = OK_STOMATE
       !Config Def   = lithology.nc
       !Config Help  = The name of the file to be opened to read
       !Config         the lithological fractions.
       !Config         The data from this file is then interpolated
       !Config         to the grid of of the model. 
       !Config Units = [FILE]
       !---
        filename = 'lithology.nc'
        CALL getin_p('LITHOLOGY_FILE',filename)
       !---
       IF (is_root_prc) CALL flininfo(filename,iml, jml, lml, tml, fid)
       CALL bcast(iml)
       CALL bcast(jml)
       CALL bcast(lml)
       CALL bcast(tml)
       !---
       ALLOCATE (lat_rel(iml,jml))
       ALLOCATE (lon_rel(iml,jml))
       ALLOCATE (lev(lml))
       ALLOCATE (laup_rel(iml,jml))
       ALLOCATE (loup_rel(iml,jml))
       ALLOCATE (lalow_rel(iml,jml))
       ALLOCATE (lolow_rel(iml,jml))
       ALLOCATE (lat_ful(iml+2,jml+2))
       ALLOCATE (lon_ful(iml+2,jml+2))
       ALLOCATE (lithology_file(iml,jml,nlcm))
       ALLOCATE (lithology_file2(iml,jml))
       !---

       IF (is_root_prc) CALL flinopen (filename, .FALSE., iml, jml, lml, &
            &                                   lon_rel, lat_rel, lev, tml, itau, date, dt, fid)


       CALL bcast(lon_rel)
       CALL bcast(lat_rel)
       CALL bcast(itau)
       CALL bcast(date)
       CALL bcast(dt)

       !---

       IF (is_root_prc) CALL flinget (fid, 'lith_frac', iml, jml, lml, tml, &
            &                                  2, 1, lithology_file(:,:,:))

       IF (is_root_prc) CALL flinget (fid, 'soil_shield', iml, jml, lml, tml, &
            &                                  2, 1, lithology_file2(:,:))
       CALL bcast(lithology_file)
       CALL bcast(lithology_file2)

       !---
       IF (is_root_prc) CALL flinclo (fid)
       !---
       !-- Duplicate the border assuming we have a global grid
       !-- going from west to east
       !---
       lon_ful(2:iml+1,2:jml+1) = lon_rel(1:iml,1:jml)
       lat_ful(2:iml+1,2:jml+1) = lat_rel(1:iml,1:jml)
       !---
       IF ( lon_rel(iml,1) < lon_ful(2,2)) THEN
          lon_ful(1,2:jml+1) = lon_rel(iml,1:jml)
          lat_ful(1,2:jml+1) = lat_rel(iml,1:jml)
       ELSE
          lon_ful(1,2:jml+1) = lon_rel(iml,1:jml)-360
          lat_ful(1,2:jml+1) = lat_rel(iml,1:jml)
       ENDIF
       !---
       IF ( lon_rel(1,1) > lon_ful(iml+1,2)) THEN
          lon_ful(iml+2,2:jml+1) = lon_rel(1,1:jml)
          lat_ful(iml+2,2:jml+1) = lat_rel(1,1:jml)
       ELSE
          lon_ful(iml+2,2:jml+1) = lon_rel(1,1:jml)+360
          lat_ful(iml+2,2:jml+1) = lat_rel(1,1:jml)
       ENDIF
       !---
       sgn = lat_rel(1,1)/ABS(lat_rel(1,1))
       lat_ful(2:iml+1,1) = sgn*180 - lat_rel(1:iml,1)
       sgn = lat_rel(1,jml)/ABS(lat_rel(1,jml))
       lat_ful(2:iml+1,jml+2) = sgn*180 - lat_rel(1:iml,jml)
       lat_ful(1,1) = lat_ful(iml+1,1)
       lat_ful(iml+2,1) = lat_ful(2,1)
       lat_ful(1,jml+2) = lat_ful(iml+1,jml+2)
       lat_ful(iml+2,jml+2) = lat_ful(2,jml+2)
       !---
       !-- Add the longitude lines to the top and bottom
       !---
       lon_ful(:,1) = lon_ful(:,2)
       lon_ful(:,jml+2) = lon_ful(:,jml+1)
       !---
       !-- Get the upper and lower limits of each grid box
       !---
       DO ip=1,iml
          DO jp=1,jml
             loup_rel(ip,jp) = &
                  &        MAX(0.5*(lon_ful(ip,jp+1)+lon_ful(ip+1,jp+1)), &
                  &            0.5*(lon_ful(ip+1,jp+1)+lon_ful(ip+2,jp+1)))
             lolow_rel(ip,jp) = &
                  &        MIN(0.5*(lon_ful(ip,jp+1)+lon_ful(ip+1,jp+1)), &
                  &            0.5*(lon_ful(ip+1,jp+1)+lon_ful(ip+2,jp+1)))
             laup_rel(ip,jp) = &
                  &        MAX(0.5*(lat_ful(ip+1,jp)+lat_ful(ip+1,jp+1)), &
                  &            0.5*(lat_ful(ip+1,jp+1)+lat_ful(ip+1,jp+2)))
             lalow_rel(ip,jp) = &
                  &        MIN(0.5*(lat_ful(ip+1,jp)+lat_ful(ip+1,jp+1)), &
                  &            0.5*(lat_ful(ip+1,jp+1)+lat_ful(ip+1,jp+2)))
          ENDDO
       ENDDO
       !---
       !-- Now we take each grid point and find out which values
       !-- from the forcing we need to average
       !---

       !
       tt(:,:) = val_exp
       tt2(:)  = val_exp

       DO ib=1,npts
          
          resx = resolution(ib,1)
          resy = resolution(ib,2)
          !-----
          do_again = .TRUE.
          !-----
          DO WHILE (do_again)
             !-----
             do_again = .FALSE.
             !-------
             !------ We find the 4 limits of the grid-box.
             !------ As we transform the resolution of the model into longitudes
             !------ and latitudes we do not have the problem of periodicity.
             !------ coslat is a help variable here !
             !-------
             coslat = MAX(COS(lalo(ib,1)*pi/180.),mincos)*pi/180.*R_Earth
             !-------
             lon_up  = lalo(ib,2)+resx/(2.0*coslat)
             lon_low = lalo(ib,2)-resx/(2.0*coslat)
             !-------
             coslat  = pi/180.*R_Earth
             !-------
             lat_up  = lalo(ib,1)+resy/(2.0*coslat)
             lat_low = lalo(ib,1)-resy/(2.0*coslat)
             !-------
             !------ Find the grid boxes from the data that go into
             !------ the model's boxes.
             !------ We still work as if we had a regular grid !
             !------ Well it needs to be localy regular so that
             !------ the longitude at the latitude of the last found point
             !------ is close to the one of the next point.
             !-------
             fopt = 0
             lastjp = 1
             DO ip=1,iml
                !---------
                !-------- Either the center of the data grid point is in the interval
                !-------- of the model grid or the East and West limits of the data
                !-------- grid point are on either sides of the border of the data grid
                !---------
                IF (      lon_rel(ip,lastjp) > lon_low &
                     &            .AND. lon_rel(ip,lastjp) < lon_up &
                     &             .OR. lolow_rel(ip,lastjp) < lon_low &
                     &            .AND. loup_rel(ip,lastjp) > lon_low &
                     &             .OR. lolow_rel(ip,lastjp) < lon_up &
                     &            .AND. loup_rel(ip,lastjp) > lon_up ) THEN
                   DO jp=1,jml
                      !-------------
                      !------------ Now that we have the longitude let us find the latitude
                      !-------------
                      IF (      lat_rel(ip,jp) > lat_low &
                           &                 .AND. lat_rel(ip,jp) < lat_up &
                           &                  .OR. lalow_rel(ip,jp) < lat_low &
                           &                 .AND. laup_rel(ip,jp) > lat_low &
                           &                  .OR. lalow_rel(ip,jp) < lat_up &
                           &                 .AND. laup_rel(ip,jp) > lat_up) THEN
                         lastjp = jp
                         !---------------
                         fopt = fopt + 1
                         IF ( fopt > nbvmax) THEN
                            WRITE(numout,*) &
                                 &                       'Please increase nbvmax in subroutine get_lithology',ib
                            STOP
                         ELSE
                            !-----------------
                            !---------------- Get the area of the fine grid in the model grid
                            !-----------------
                            coslat = MAX(COS(lat_rel(ip,jp)*pi/180.),mincos)
                            ax =  ( MIN(lon_up,loup_rel(ip,jp)) &
                                 &                       -MAX(lon_low,lolow_rel(ip,jp))) &
                                 &                     *pi/180.*R_Earth*coslat
                            ay =  ( MIN(lat_up,laup_rel(ip,jp)) &
                                 &                       -MAX(lat_low,lalow_rel(ip,jp))) &
                                 &                     *pi/180.*R_Earth
                            area(fopt) = ax*ay
                            tt(fopt,:) = lithology_file(ip,jp,:)
                            tt2(fopt)  = lithology_file2(ip,jp)
                         ENDIF
                      ENDIF
                   ENDDO
                ENDIF
             ENDDO
            
             !!
             ! set an average lithological composition to all boxes
             ! which are not covered (.GT.un) in the lithology file
             WHERE(tt(1:fopt,1) .GT. un)
                   tt(:,1)  = 0.0019021
                   tt(:,2)  = 0.27252
                   tt(:,3)  = 0.10161
                   tt(:,4)  = 0.00069925
                   tt(:,5)  = 0.048499
                   tt(:,6)  = 0.0081311
                   tt(:,7)  = 0.0035592
                   tt(:,8)  = 0.0044056
                   tt(:,9)  = 0.059839
                   tt(:,10) = 0.13715
                   tt(:,11) = 0.13197
                   tt(:,12) = 0.16555
                   tt(:,13) = 0.0085954
                   tt(:,14) = 0.032713
                   tt(:,15) = 0.016390
                   tt(:,16) = 0.0064659
             ENDWHERE

             ! soil shielding
             WHERE( tt2(1:fopt) .GT. un)
                   tt2(:) = 1.
             ENDWHERE

             !!

             !-------
             !------ Check that we found some points
             !-------
             lith(ib,:)     = zero
             soilshield(ib) = zero
             !-------
             IF (fopt == 0) THEN
                do_again = .TRUE.
                !-------
                !------ increase search radius
                !-------
                resx = resx*2.
                resy = resy*2.
                IF ( resx > 2.*pi*R_Earth .OR. resy > pi*R_Earth ) THEN
                   STOP 'get_lithology: found no point'
                ENDIF
             ELSE
                sgn = zero
                !-------
                !------ Compute the dominat lithology
                !-------
                DO ilf=1,fopt
                   lith(ib,:)     = lith(ib,:)     + tt(ilf,:) * area(ilf)
                   soilshield(ib) = soilshield(ib) + tt2(ilf)  * area(ilf)
                   sgn = sgn + area(ilf)
                ENDDO
                !-------
                !------ Normalize the surface
                !-------
                IF (sgn < min_sechiba) THEN
                   do_again = .TRUE.
                   !---------
                   !-------- increase search radius
                   !---------
                   resx = resx * 2.
                   resy = resy * 2.
                   IF ( resx > 2.*pi*R_Earth .OR. resy > pi*R_Earth ) THEN
                      STOP 'get_lithology: found no point'
                   ENDIF
                ELSE
                   lith(ib,:)     = lith(ib,:) / sgn
                   soilshield(ib) = soilshield(ib) / sgn
                ENDIF
             ENDIF
          ENDDO
       ENDDO
       !-
       ! deallocate
       !-
       DEALLOCATE (lat_rel)
       DEALLOCATE (lon_rel)
       DEALLOCATE (laup_rel)
       DEALLOCATE (loup_rel)
       DEALLOCATE (lalow_rel)
       DEALLOCATE (lolow_rel)
       DEALLOCATE (lat_ful)
       DEALLOCATE (lon_ful)
       DEALLOCATE (lithology_file)
       DEALLOCATE (lithology_file2)
    ENDIF
    !-
    ! 2 test if the fractions still add up to 1 
    !-

    !DSG this test doesn't work because of missing land points in the lithology file
    IF (ANY(ABS(un - SUM(lith(1:npts,:),DIM=2)).GT. min_stomate*1e3)) THEN
        WRITE (6,*) 'maxval',MAXVAL(SUM(lith(1:npts,:),DIM=2))
        WRITE (6,*) 'minval',MINVAL(SUM(lith(1:npts,:),DIM=2))
        STOP 'get_lithology: the interpolation lead to fractions which dont add up to one'
    ENDIF

    !-
    ! 3 output the lithological fractions & soil shield factor
    !-
    lith_out(:,:) = lith(:,:)
    soilshield_out(:) = soilshield(:)

    IF (printlev >= 4) WRITE(numout,*) 'Leaving get_lithology'

  END SUBROUTINE get_lithology

!! ================================================================================================================================
!! SUBROUTINE   : get_deposition
!!
!>\BRIEF         This subroutine reads the deposition of N & P from a boundary condition file. Its an adaption of get_reftemp.
!!
!! DESCRIPTION  : This subroutine reads the deposition of N & P from a boundary condition file. No time dimension, yet. 
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): None
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE get_deposition (npts, lalo, resolution, depo_out)

    IMPLICIT NONE

    !! 0. Parameters and variables declaration
    
    !! 0.1 Input variables 

    INTEGER, INTENT(in)                 :: npts                 !! Domain size
    REAL, DIMENSION(npts,2), INTENT(in) :: lalo                 !! Geogr. coordinates (latitude,longitude) (degrees)
    REAL, DIMENSION(npts,2), INTENT(in) :: resolution           !! size in x an y of the grid (m)

    !! 0.2 Output variables

    REAL, DIMENSION(npts,3), INTENT(out)  :: depo_out           !! input by deposition of either N (NOY=1,NHX=2) or P(=3) deposition

    !! 0.4 Local variables

    INTEGER, PARAMETER                  :: nbvmax = 200         !!
    CHARACTER(LEN=80)                          :: filename             !!
    INTEGER                             :: iml, jml             !!
    INTEGER                             :: lml                  !!
    INTEGER                             :: tml                  !!
    INTEGER                             :: fid                  !! 
    INTEGER                             :: ib, ip, jp           !!
    INTEGER                             :: fopt                 !!
    INTEGER                             :: ilf                  !!
    INTEGER                             :: lastjp               !! 
    REAL, DIMENSION(1)                  :: lev                  !!
    REAL                                :: date                 !!
    REAL                                :: dt                   !!
    REAL                                :: coslat               !!
    INTEGER                             :: itau(1)              !!
    REAL, ALLOCATABLE, DIMENSION(:,:)   :: lat_rel, lon_rel     !!   
    REAL, ALLOCATABLE, DIMENSION(:,:)   :: lat_ful, lon_ful     !! 
    REAL, ALLOCATABLE, DIMENSION(:,:,:) :: deposition_file      !! 
    REAL, ALLOCATABLE, DIMENSION(:,:)   :: loup_rel, lolow_rel  !!
    REAL, ALLOCATABLE, DIMENSION(:,:)   :: laup_rel, lalow_rel  !!
    REAL                                :: lon_up, lon_low      !! 
    REAL                                :: lat_up, lat_low      !!
    REAL                                :: ax, ay, sgn          !!
    REAL, DIMENSION(nbvmax)             :: area                 !!
    REAL, DIMENSION(nbvmax,3)           :: tt                   !!
    REAL                                :: resx, resy           !!
    LOGICAL                                    :: do_again             !!

!_ ================================================================================================================================

    !! 1 If this is the first call, calculate the deposition
    !!   and keep it in memory

    IF (printlev >= 4) WRITE(numout,*) 'Entering get_deposition'

    IF (firstcall_dep) THEN
       !---
       !-- 1.1 only do this once
       !---
       firstcall_dep = .FALSE.
       !---
       !-- 1.2 allocate the field
       !---
       ALLOCATE( depo(npts,3) )
       !---
       !-- 1.3 read and interpolate the deposition file
       !---
       !-- Needs to be a configurable variable
       !---
       !Config Key   = DEPOSITION_FILE
       !Config Desc  = Name of file from which the deposition is read
       !Config If    = OK_STOMATE
       !Config Def   = deposition.nc
       !Config Help  = The name of the file to be opened to read
       !Config         the reference surface temperature.
       !Config         The data from this file is then interpolated
       !Config         to the grid of of the model.
       !Config Units = [FILE]
       !---
        filename = 'deposition.nc'
       CALL getin_p('DEPOSITION_FILE',filename)

       !---
       IF (is_root_prc) CALL flininfo(filename,iml, jml, lml, tml, fid)
       CALL bcast(iml)
       CALL bcast(jml)
       CALL bcast(lml)
       CALL bcast(tml)
       !---
       ALLOCATE (lat_rel(iml,jml))
       ALLOCATE (lon_rel(iml,jml))
       ALLOCATE (laup_rel(iml,jml))
       ALLOCATE (loup_rel(iml,jml))
       ALLOCATE (lalow_rel(iml,jml))
       ALLOCATE (lolow_rel(iml,jml))
       ALLOCATE (lat_ful(iml+2,jml+2))
       ALLOCATE (lon_ful(iml+2,jml+2))
       ALLOCATE (deposition_file(iml,jml,3))
       !---
       IF (is_root_prc) CALL flinopen (filename, .FALSE., iml, jml, lml, &
            &                                   lon_rel, lat_rel, lev, tml, itau, date, dt, fid)
       CALL bcast(lon_rel)
       CALL bcast(lat_rel)
       CALL bcast(itau)
       CALL bcast(date)
       CALL bcast(dt)

       !---
     !DSG: since the last merge N deposition is already read in slowproc. 
           ! but I will continue using this; long-term, this must be changed
       IF (is_root_prc) CALL flinget (fid, 'NOY', iml, jml, lml, tml, &
            &                                  1, 1, deposition_file(:,:,1))
       IF (is_root_prc) CALL flinget (fid, 'NHX', iml, jml, lml, tml, &
            &                                  1, 1, deposition_file(:,:,2))

       IF (is_root_prc) CALL flinget (fid, 'p_depo', iml, jml, lml, tml, &
            &                                  1, 1, deposition_file(:,:,3))
       CALL bcast(deposition_file)
       !---
       IF (is_root_prc) CALL flinclo (fid)
       !---
       !-- Duplicate the border assuming we have a global grid
       !-- going from west to east
       !---
       lon_ful(2:iml+1,2:jml+1) = lon_rel(1:iml,1:jml)
       lat_ful(2:iml+1,2:jml+1) = lat_rel(1:iml,1:jml)
       !---
       IF ( lon_rel(iml,1) < lon_ful(2,2)) THEN
          lon_ful(1,2:jml+1) = lon_rel(iml,1:jml)
          lat_ful(1,2:jml+1) = lat_rel(iml,1:jml)
       ELSE
          lon_ful(1,2:jml+1) = lon_rel(iml,1:jml)-360
          lat_ful(1,2:jml+1) = lat_rel(iml,1:jml)
       ENDIF
       !---
       IF ( lon_rel(1,1) > lon_ful(iml+1,2)) THEN
          lon_ful(iml+2,2:jml+1) = lon_rel(1,1:jml)
          lat_ful(iml+2,2:jml+1) = lat_rel(1,1:jml)
       ELSE
          lon_ful(iml+2,2:jml+1) = lon_rel(1,1:jml)+360
          lat_ful(iml+2,2:jml+1) = lat_rel(1,1:jml)
       ENDIF
       !---
       sgn = lat_rel(1,1)/ABS(lat_rel(1,1))
       lat_ful(2:iml+1,1) = sgn*180 - lat_rel(1:iml,1)
       sgn = lat_rel(1,jml)/ABS(lat_rel(1,jml))
       lat_ful(2:iml+1,jml+2) = sgn*180 - lat_rel(1:iml,jml)
       lat_ful(1,1) = lat_ful(iml+1,1)
       lat_ful(iml+2,1) = lat_ful(2,1)
       lat_ful(1,jml+2) = lat_ful(iml+1,jml+2)
       lat_ful(iml+2,jml+2) = lat_ful(2,jml+2)
       !---
       !-- Add the longitude lines to the top and bottom
       !---
       lon_ful(:,1) = lon_ful(:,2)
       lon_ful(:,jml+2) = lon_ful(:,jml+1)
       !---
       !-- Get the upper and lower limits of each grid box
       !---
       DO ip=1,iml
          DO jp=1,jml
             loup_rel(ip,jp) = &
                  &        MAX(0.5*(lon_ful(ip,jp+1)+lon_ful(ip+1,jp+1)), &
                  &            0.5*(lon_ful(ip+1,jp+1)+lon_ful(ip+2,jp+1)))
             lolow_rel(ip,jp) = &
                  &        MIN(0.5*(lon_ful(ip,jp+1)+lon_ful(ip+1,jp+1)), &
                  &            0.5*(lon_ful(ip+1,jp+1)+lon_ful(ip+2,jp+1)))
             laup_rel(ip,jp) = &
                  &        MAX(0.5*(lat_ful(ip+1,jp)+lat_ful(ip+1,jp+1)), &
                  &            0.5*(lat_ful(ip+1,jp+1)+lat_ful(ip+1,jp+2)))
             lalow_rel(ip,jp) = &
                  &        MIN(0.5*(lat_ful(ip+1,jp)+lat_ful(ip+1,jp+1)), &
                  &            0.5*(lat_ful(ip+1,jp+1)+lat_ful(ip+1,jp+2)))
          ENDDO
       ENDDO
       !---
       !-- Now we take each grid point and find out which values
       !-- from the forcing we need to average
       !---
       DO ib=1,npts
          !-----
          resx = resolution(ib,1)
          resy = resolution(ib,2)
          !-----
          do_again = .TRUE.
          !-----
          DO WHILE (do_again)
             !-----
             do_again = .FALSE.
             !-------
             !------ We find the 4 limits of the grid-box.
             !------ As we transform the resolution of the model into longitudes
             !------ and latitudes we do not have the problem of periodicity.
             !------ coslat is a help variable here !
             !-------
             coslat = MAX(COS(lalo(ib,1)*pi/180.),mincos)*pi/180.*R_Earth
             !-------
             lon_up  = lalo(ib,2)+resx/(2.0*coslat)
             lon_low = lalo(ib,2)-resx/(2.0*coslat)
             !-------
             coslat  = pi/180.*R_Earth
             !-------
             lat_up  = lalo(ib,1)+resy/(2.0*coslat)
             lat_low = lalo(ib,1)-resy/(2.0*coslat)
             !-------
             !------ Find the grid boxes from the data that go into
             !------ the model's boxes.
             !------ We still work as if we had a regular grid !
             !------ Well it needs to be localy regular so that
             !------ the longitude at the latitude of the last found point
             !------ is close to the one of the next point.
             !-------
             fopt = 0
             lastjp = 1
             DO ip=1,iml
                !---------
                !-------- Either the center of the data grid point is in the interval
                !-------- of the model grid or the East and West limits of the data
                !-------- grid point are on either sides of the border of the data grid
                !---------
                IF (      lon_rel(ip,lastjp) > lon_low &
                     &            .AND. lon_rel(ip,lastjp) < lon_up &
                     &             .OR. lolow_rel(ip,lastjp) < lon_low &
                     &            .AND. loup_rel(ip,lastjp) > lon_low &
                     &             .OR. lolow_rel(ip,lastjp) < lon_up &
                     &            .AND. loup_rel(ip,lastjp) > lon_up ) THEN
                   DO jp=1,jml
                      !-------------
                      !------------ Now that we have the longitude let us find the latitude
                      !-------------
                      IF (      lat_rel(ip,jp) > lat_low &
                           &                 .AND. lat_rel(ip,jp) < lat_up &
                           &                  .OR. lalow_rel(ip,jp) < lat_low &
                           &                 .AND. laup_rel(ip,jp) > lat_low &
                           &                  .OR. lalow_rel(ip,jp) < lat_up &
                           &                 .AND. laup_rel(ip,jp) > lat_up) THEN
                         lastjp = jp
                         !---------------
                         fopt = fopt + 1
                         IF ( fopt > nbvmax) THEN
                            WRITE(numout,*) &
                                 &                       'Please increase nbvmax in subroutine get_deposition',ib
                            STOP
                         ELSE
                            !-----------------
                            !---------------- Get the area of the fine grid in the model grid
                            !-----------------
                            coslat = MAX(COS(lat_rel(ip,jp)*pi/180.),mincos)
                            ax =  ( MIN(lon_up,loup_rel(ip,jp)) &
                                 &                       -MAX(lon_low,lolow_rel(ip,jp))) &
                                 &                     *pi/180.*R_Earth*coslat
                            ay =  ( MIN(lat_up,laup_rel(ip,jp)) &
                                 &                       -MAX(lat_low,lalow_rel(ip,jp))) &
                                 &                     *pi/180.*R_Earth
                            area(fopt) = ax*ay
                            tt(fopt,:) = deposition_file(ip,jp,:)
                         ENDIF
                      ENDIF
                   ENDDO
                ENDIF
             ENDDO
             !-------
             !------ Check that we found some points
             !-------
             depo(ib,:) = zero
             !-------
             IF (fopt == 0) THEN
                do_again = .TRUE.
                !-------
                !------ increase search radius
                !-------
                resx = resx*2.
                resy = resy*2.
                IF ( resx > 2.*pi*R_Earth .OR. resy > pi*R_Earth ) THEN
                   STOP 'get_deposition: found no point'
                ENDIF
             ELSE
                sgn = zero
                !-------
                !------ Compute the average deposition
                !-------
                DO ilf=1,fopt
                   depo(ib,1) = depo(ib,1) + tt(ilf,1) * area(ilf)
                   depo(ib,2) = depo(ib,2) + tt(ilf,2) * area(ilf)
                   depo(ib,3) = depo(ib,3) + tt(ilf,3) * area(ilf)

                   sgn = sgn + area(ilf)
                ENDDO
                !-------
                !------ Normalize the surface
                !-------
                IF (sgn < min_sechiba) THEN
                   do_again = .TRUE.
                   !---------
                   !-------- increase search radius
                   !---------
                   resx = resx * 2.
                   resy = resy * 2.
                   IF ( resx > 2.*pi*R_Earth .OR. resy > pi*R_Earth ) THEN
                      STOP 'get_deposition: found no point'
                   ENDIF
                ELSE
                   depo(ib,1) = depo(ib,1) / sgn
                   depo(ib,2) = depo(ib,2) / sgn
                   depo(ib,3) = depo(ib,3) / sgn
                ENDIF
             ENDIF
          ENDDO
       ENDDO
       !-
       ! deallocate
       !-
       DEALLOCATE (lat_rel)
       DEALLOCATE (lon_rel)
       DEALLOCATE (laup_rel)
       DEALLOCATE (loup_rel)
       DEALLOCATE (lalow_rel)
       DEALLOCATE (lolow_rel)
       DEALLOCATE (lat_ful)
       DEALLOCATE (lon_ful)
       DEALLOCATE (deposition_file)
    ENDIF
    !-
    ! 2 output the deposition
    !-
    ! convert units g m-2 s-1 to g m-2 day-1
    depo_out(:,:) = depo(:,:)*(60.*60.*24.)

    IF (printlev >= 4) WRITE(numout,*) 'Leaving get_deposition'

  END SUBROUTINE get_deposition


!
!=
!
!! ================================================================================================================================
!! SUBROUTINE   : get_soil_orders
!!
!>\BRIEF        : This subroutine reads the soil order from a boundary condition file. Its an adaption of get_reftemp.
!!
!! DESCRIPTION  : This subroutine reads the soil order from a boundary condition file. Its an adaption of get_reftemp.
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): None
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE get_soil_orders (npts, lalo, resolution, soil_order_imposed,soilorder_out)

    IMPLICIT NONE

    !! 0. Parameters and variables declaration
    
    !! 0.1 Input variables 

    INTEGER, INTENT(in)                 :: npts                 !! Domain size
    REAL, DIMENSION(npts,2), INTENT(in) :: lalo                 !! Geogr. coordinates (latitude,longitude) (degrees)
    REAL, DIMENSION(npts,2), INTENT(in) :: resolution           !! size in x an y of the grid (m)
    REAL, DIMENSION(npts), INTENT(in)   :: soil_order_imposed    !! soil order USDA when imposed

    !! 0.2 Output variables

    INTEGER, DIMENSION(npts), INTENT(out)  :: soilorder_out     !! soil order USDA 

    !! 0.4 Local variables

    INTEGER, PARAMETER                  :: nbvmax = 200         !!
    CHARACTER(LEN=80)                          :: filename             !!
    INTEGER                             :: iml, jml             !!
    INTEGER                             :: lml                  !!
    INTEGER                             :: tml                  !!
    INTEGER                             :: fid                  !! 
    INTEGER                             :: ib, ip, jp           !!
    INTEGER                             :: fopt                 !!
    INTEGER                             :: lastjp               !! 
    REAL, DIMENSION(1)                  :: lev                  !!
    REAL                                :: date                 !!
    REAL                                :: dt                   !!
    REAL                                :: coslat               !!
    INTEGER                             :: itau(1)              !!
    REAL, ALLOCATABLE, DIMENSION(:,:)   :: lat_rel, lon_rel     !!   
    REAL, ALLOCATABLE, DIMENSION(:,:)   :: lat_ful, lon_ful     !! 
    REAL, ALLOCATABLE, DIMENSION(:,:)   :: soil_order_file      !! 
    REAL, ALLOCATABLE, DIMENSION(:,:)   :: loup_rel, lolow_rel  !!
    REAL, ALLOCATABLE, DIMENSION(:,:)   :: laup_rel, lalow_rel  !!
    REAL                                :: lon_up, lon_low      !! 
    REAL                                :: lat_up, lat_low      !!
    REAL                                :: ax, ay, sgn          !!
    REAL, DIMENSION(nbvmax)             :: area                 !!
    REAL, DIMENSION(nbvmax)             :: tt                   !!
    REAL                                :: resx, resy           !!
    LOGICAL                                    :: do_again             !!
    INTEGER, DIMENSION(1)               :: idx                  !!

!_ ================================================================================================================================

    !! 1 If this is the first call, read in sorder
    !!   and keep it in memory

    IF (printlev >= 4) WRITE(numout,*) 'Entering get_soil_orders'

    IF (.NOT.firstcall_sorder) THEN
        CALL ipslerr_p(3, 'get_soil_order', 'We will STOP now. This subroutine', &
         & 'is already initialized', '' )
    ENDIF
       !---
       !-- 1.1 only do this once
       !---
       firstcall_sorder = .FALSE.
       !---
       !-- 1.2 allocate the field
       !---
       ALLOCATE( soilorder(npts) )
       idx = -1
       !---
       !-- 1.3 read and interpolate the Soil order file
       !---
       !-- Needs to be a configurable variable
       !---
       !Config Key   = SOILORDER_FILE
       !Config Desc  = Name of file from which the soil order is read
       !Config If    = OK_STOMATE
       !Config Def   = USDA_SoilSuborder.nc
       !Config Help  = The name of the file to be opened to read
       !Config         the reference surface temperature.
       !Config         The data from this file is then interpolated
       !Config         to the grid of of the model.
       !Config Units = [FILE]
       !---
        filename = 'USDA_SoilSuborder.nc'
        CALL getin_p('SOILORDER_FILE',filename)
       !---
     


       IF (is_root_prc) CALL flininfo(filename,iml, jml, lml, tml, fid)
       CALL bcast(iml)
       CALL bcast(jml)
       CALL bcast(lml)
       CALL bcast(tml)
       !---
       ALLOCATE (lat_rel(iml,jml))
       ALLOCATE (lon_rel(iml,jml))
       ALLOCATE (laup_rel(iml,jml))
       ALLOCATE (loup_rel(iml,jml))
       ALLOCATE (lalow_rel(iml,jml))
       ALLOCATE (lolow_rel(iml,jml))
       ALLOCATE (lat_ful(iml+2,jml+2))
       ALLOCATE (lon_ful(iml+2,jml+2))
       ALLOCATE (soil_order_file(iml,jml))
       !---
       IF (is_root_prc) CALL flinopen (filename, .FALSE., iml, jml, lml, &
            &                                   lon_rel, lat_rel, lev, tml, itau, date, dt, fid)
       CALL bcast(lon_rel)
       CALL bcast(lat_rel)
       CALL bcast(itau)
       CALL bcast(date)
       CALL bcast(dt)

       !---
       IF (is_root_prc) CALL flinget (fid, 'USDA_soil_map', iml, jml, lml, tml, &
            &                                  1, 1, soil_order_file(:,:))

       CALL bcast(soil_order_file)
       !---
       IF (is_root_prc) CALL flinclo (fid)
       !---
       !-- Duplicate the border assuming we have a global grid
       !-- going from west to east
       !---
       lon_ful(2:iml+1,2:jml+1) = lon_rel(1:iml,1:jml)
       lat_ful(2:iml+1,2:jml+1) = lat_rel(1:iml,1:jml)
       !---


       IF ( lon_rel(iml,1) < lon_ful(2,2)) THEN
          lon_ful(1,2:jml+1) = lon_rel(iml,1:jml)
          lat_ful(1,2:jml+1) = lat_rel(iml,1:jml)
       ELSE
          lon_ful(1,2:jml+1) = lon_rel(iml,1:jml)-360
          lat_ful(1,2:jml+1) = lat_rel(iml,1:jml)
       ENDIF
       !---
       IF ( lon_rel(1,1) > lon_ful(iml+1,2)) THEN
          lon_ful(iml+2,2:jml+1) = lon_rel(1,1:jml)
          lat_ful(iml+2,2:jml+1) = lat_rel(1,1:jml)
       ELSE
          lon_ful(iml+2,2:jml+1) = lon_rel(1,1:jml)+360
          lat_ful(iml+2,2:jml+1) = lat_rel(1,1:jml)
       ENDIF
       !---
       sgn = lat_rel(1,1)/ABS(lat_rel(1,1))
       lat_ful(2:iml+1,1) = sgn*180 - lat_rel(1:iml,1)
       sgn = lat_rel(1,jml)/ABS(lat_rel(1,jml))
       lat_ful(2:iml+1,jml+2) = sgn*180 - lat_rel(1:iml,jml)
       lat_ful(1,1) = lat_ful(iml+1,1)
       lat_ful(iml+2,1) = lat_ful(2,1)
       lat_ful(1,jml+2) = lat_ful(iml+1,jml+2)
       lat_ful(iml+2,jml+2) = lat_ful(2,jml+2)
       !---
       !-- Add the longitude lines to the top and bottom
       !---
       lon_ful(:,1) = lon_ful(:,2)
       lon_ful(:,jml+2) = lon_ful(:,jml+1)
       !---
       !-- Get the upper and lower limits of each grid box
       !---
       DO ip=1,iml
          DO jp=1,jml
             loup_rel(ip,jp) = &
                  &        MAX(0.5*(lon_ful(ip,jp+1)+lon_ful(ip+1,jp+1)), &
                  &            0.5*(lon_ful(ip+1,jp+1)+lon_ful(ip+2,jp+1)))
             lolow_rel(ip,jp) = &
                  &        MIN(0.5*(lon_ful(ip,jp+1)+lon_ful(ip+1,jp+1)), &
                  &            0.5*(lon_ful(ip+1,jp+1)+lon_ful(ip+2,jp+1)))
             laup_rel(ip,jp) = &
                  &        MAX(0.5*(lat_ful(ip+1,jp)+lat_ful(ip+1,jp+1)), &
                  &            0.5*(lat_ful(ip+1,jp+1)+lat_ful(ip+1,jp+2)))
             lalow_rel(ip,jp) = &
                  &        MIN(0.5*(lat_ful(ip+1,jp)+lat_ful(ip+1,jp+1)), &
                  &            0.5*(lat_ful(ip+1,jp+1)+lat_ful(ip+1,jp+2)))
          ENDDO
       ENDDO
       !---
       !-- Now we take each grid point and find out which values
       !-- from the forcing we need to average
       !---
       DO ib=1,npts
          !-----
          resx = resolution(ib,1)
          resy = resolution(ib,2)
          !-----
          do_again = .TRUE.
          !-----
          DO WHILE (do_again)
             !-----
             do_again = .FALSE.
             !-------
             !------ We find the 4 limits of the grid-box.
             !------ As we transform the resolution of the model into longitudes
             !------ and latitudes we do not have the problem of periodicity.
             !------ coslat is a help variable here !
             !-------
             coslat = MAX(COS(lalo(ib,1)*pi/180.),mincos)*pi/180.*R_Earth
             !-------
             lon_up  = lalo(ib,2)+resx/(2.0*coslat)
             lon_low = lalo(ib,2)-resx/(2.0*coslat)
             !-------
             coslat  = pi/180.*R_Earth
             !-------
             lat_up  = lalo(ib,1)+resy/(2.0*coslat)
             lat_low = lalo(ib,1)-resy/(2.0*coslat)
             !-------
             !------ Find the grid boxes from the data that go into
             !------ the model's boxes.
             !------ We still work as if we had a regular grid !
             !------ Well it needs to be localy regular so that
             !------ the longitude at the latitude of the last found point
             !------ is close to the one of the next point.
             !-------
             fopt = zero
             !DSGfix_nbvmax
             area(:)= zero
             !DSGfix_nbvmax
             lastjp = un
             DO ip=1,iml
                !---------
                !-------- Either the center of the data grid point is in the interval
                !-------- of the model grid or the East and West limits of the data
                !-------- grid point are on either sides of the border of the data grid
                !---------
                IF (      lon_rel(ip,lastjp) > lon_low &
                     &            .AND. lon_rel(ip,lastjp) < lon_up &
                     &             .OR. lolow_rel(ip,lastjp) < lon_low &
                     &            .AND. loup_rel(ip,lastjp) > lon_low &
                     &             .OR. lolow_rel(ip,lastjp) < lon_up &
                     &            .AND. loup_rel(ip,lastjp) > lon_up ) THEN
                   DO jp=1,jml
                      !-------------
                      !------------ Now that we have the longitude let us find the latitude
                      !-------------
                      IF (      lat_rel(ip,jp) > lat_low &
                           &                 .AND. lat_rel(ip,jp) < lat_up &
                           &                  .OR. lalow_rel(ip,jp) < lat_low &
                           &                 .AND. laup_rel(ip,jp) > lat_low &
                           &                  .OR. lalow_rel(ip,jp) < lat_up &
                           &                 .AND. laup_rel(ip,jp) > lat_up) THEN
                         lastjp = jp


                         !---------------
                         !DSG: FM NN update: start
                    !DSGnope     IF(soil_order_file(ip,jp) .LE. 12) THEN   !@FM
                            fopt = fopt + 1
                            IF ( fopt > nbvmax) THEN
                                WRITE(numout,*)  'Please increase nbvmax in subroutine get_soil_orders',ib
                               STOP
                            ELSE
                               !-----------------
                               !---------------- Get the area of the fine grid in the model grid
                               !-----------------
                               coslat = MAX(COS(lat_rel(ip,jp)*pi/180.),mincos)
                               ax =  ( MIN(lon_up,loup_rel(ip,jp)) &
                                    &                       -MAX(lon_low,lolow_rel(ip,jp))) &
                                    &                     *pi/180.*R_Earth*coslat
                               ay =  ( MIN(lat_up,laup_rel(ip,jp)) &
                                    &                       -MAX(lat_low,lalow_rel(ip,jp))) &
                                    &                     *pi/180.*R_Earth
                               area(fopt) = ax*ay
                          !    assign soil order 
                                IF(soil_order_file(ip,jp) .GT. 12) THEN ! no soil order information
                                   WRITE(numout,*) &
                                        &                   'no soil order found in file for grid point',ib
                                   WRITE(numout,*) &
                                     &                      'we assign an aridisol:', 3
                                   tt(fopt)  = 3.
                                ELSE
                                  tt(fopt) = soil_order_file(ip,jp)
                                ENDIF
                               !DSG: FM NN update: end
                            ENDIF
                        !DSGnope ENDIF ! soil_order_file .LE. 12


                      ENDIF
                   ENDDO
                ENDIF
             ENDDO
             !-------
             !------ Check that we found some points
             !-------
             soilorder(ib) = zero
             !-------
             IF (fopt == 0) THEN
                do_again = .TRUE.
                !-------
                !------ increase search radius
                !-------
                resx = resx*2.
                resy = resy*2.
                IF ( resx > 2.*pi*R_Earth .OR. resy > pi*R_Earth ) THEN
                   STOP 'get_soil_order: found no point'
                ENDIF
             ELSE
                sgn = zero
                !-------
                !------ Select the soil order with biggest area
                !-------
                idx = MAXLOC(area)
                soilorder(ib) = tt(idx(1))
                sgn = SUM(area)
                !-------
                !------ Normalize the surface
                !-------
                IF (sgn < min_sechiba) THEN
                   do_again = .TRUE.
                   !---------
                   !-------- increase search radius
                   !---------
                   resx = resx * 2.
                   resy = resy * 2.
                   IF ( resx > 2.*pi*R_Earth .OR. resy > pi*R_Earth ) THEN
                      STOP 'get_soil_order: found no point'
                   ENDIF
                ENDIF
             ENDIF
          ENDDO
       ENDDO
       !-
       ! deallocate
       !-
       DEALLOCATE (lat_rel)
       DEALLOCATE (lon_rel)
       DEALLOCATE (laup_rel)
       DEALLOCATE (loup_rel)
       DEALLOCATE (lalow_rel)
       DEALLOCATE (lolow_rel)
       DEALLOCATE (lat_ful)
       DEALLOCATE (lon_ful)
       DEALLOCATE (soil_order_file)
    !-
    ! 2 output the soil order as read from file
    !-
    WRITE(numout,*) 'IMPOSE_SOILS=', impsoils
    WRITE(numout,*) 'We impose soil order=', soilorder(:)

    IF (.NOT.impsoils) THEN
      soilorder_out(:) = INT(soilorder(:))
    ELSE 
      WRITE(numout,*) 'WE impose the soil order to all grids to='
      soilorder_out(:) = INT(soil_order_imposed(:))
    ENDIF

    IF (printlev >= 4) WRITE(numout,*) 'Leaving get_soil_orders'

  END SUBROUTINE get_soil_orders
!
!=
!
END MODULE stomate_io
