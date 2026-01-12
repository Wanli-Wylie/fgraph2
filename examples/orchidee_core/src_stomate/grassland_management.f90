! =================================================================================================================================
! MODULE       : grassland_management
!
! CONTACT      : orchidee-help _at_ ipsl.jussieu.fr
!
! LICENCE      : IPSL (2006)
! This software is governed by the CeCILL licence see
! ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF       Groups the subroutines that: (1) initialize all variables in
!! grassland management, (2) call subroutines of major grassland management 
!! modules (cut/harvest, grazing, fertilization) , (3) read external maps 
!! such as mangement, livestock density, fertilization, wild animal density
!! (4) calculate plants status such as developement stage and tgrowth
!!
!!\n DESCRIPTION : None
!!
!! RECENT CHANGE(S) : None
!!
!! REFERENCE(S) : None
!!
!! \n
!_
!================================================================================================================================

MODULE grassland_management
  ! this module include grassland management from PaSim
  ! graze - cut - fertilisation 
  ! with auto management or user-defined management

  USE grassland_constantes
  USE constantes
  USE grassland_fonctions
  USE grassland_grazing
  !USE applic_plant
  USE grassland_cutting
  USE grassland_fertilisation
  USE ioipsl
  USE ioipsl_para
  USE mod_orchidee_para
  USE xios_orchidee
  USE netcdf
  USE defprec
  USE grid
  USE matrix_resolution
  USE interpol_help  ! necessary for management map input
!  USE parallel

  IMPLICIT NONE

  PUBLIC grassmanag_clear,slowproc_GRM_ManageInput

  LOGICAL, SAVE :: first_call_grassland_manag = .TRUE. 
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: intake
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: intakemax
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: intake_litter
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: intake_animal_litter
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: litter_avail_totDM
  INTEGER, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: grazing_litter
  ! shoot dry matter afer cutting Kg/m2
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: wshtotcutinit
  ! lai after cutting (m**2 leaf/m**2)
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: lcutinit     
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: devstage
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: faecesc
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: faecesn
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: urinen
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: urinec
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: nel
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: nanimaltot
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: tgrowth               
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: wsh
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: wshtotinit
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: wshtot
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: wr
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: wrtot
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: wanimal
  ! concentration totale en N (kg n/kg)
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: ntot   
  ! concentration en C du substrat de la plante (kg C/kg)
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: c      
  ! concentration en N du substrat de la plante (kg N/kg)
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: n      
  ! n in structral mass kgN/kg
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: fn     
  ! n concentration of apoplast (kgN/kg)
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: napo   
  ! n concentration of symplast (kgN/kg)
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: nsym   
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: wnapo
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: wnsym
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: wn
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: nanimal
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: tanimal
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: danimal             
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: tcut  
  ! day of fertilisation (management) (d)           
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: tfert   
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: Nliquidmanure    
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: nslurry           
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: Nsolidmanure      
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: legume_fraction
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:)    :: soil_fertility
  ! Threshold shoot dry matter, under which animals are moved out (kg/m2)
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: Animalwgrazingmin
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: AnimalkintakeM
  ! parameter for calculation of vegetation compartement selection by animals (-)
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: AnimalDiscremineQualite
  ! Valeurs associ�es � la croissance a�rienne entre 2 �v�nements de fertilisation (autogestion fauches)
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: controle_azote 
  ! Carbon flux from Organic fertilization to metabolic SOM pool (kg C m-2 day-1)
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: fcOrganicFertmetabolicsum    
  ! Carbon flux from Organic fertilization to strcutural SOM pool (kg C m-2 day-1)
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: fcOrganicFertstructsum       
  ! Nitrogen flux from Organic fertilization to strcutural SOM pool (kg N m-2 day-1)
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: fnOrganicFertmetabolicsum    
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: fnOrganicFertstructsum
  ! Nitrogen flux coming from slurry and liquid manure (k N.m-2)
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: fnOrganicFerturinesum        
  ! Nitrogen deposition (kg N m-2 year-1)
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: fnatmsum                     
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: controle_azote_sum
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   :: nfertamm
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   :: nfertnit
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: intakesum
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: intakensum
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: intake_animal
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: intake_animalsum
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: PIYcow
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: PIMcow
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: BCSYcow
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: BCSMcow
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: PICcow
  ! Age of dairy primi cow
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: AGE_cow_P           
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: AGE_cow_M
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: Autogestion_out
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: Forage_quantity
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: tcut_modif
  INTEGER, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: countschedule
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: mux
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: mugmean
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: sigx
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: sigy
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: gmeanslope
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: gzero
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: gcor      
  INTEGER (i_std)  , ALLOCATABLE, SAVE, DIMENSION(:,:)   :: cuttingend
  LOGICAL   , ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: tcut_verif
  INTEGER   , ALLOCATABLE, SAVE, DIMENSION(:,:)   :: regcount
  INTEGER                                       :: tcutmodel
!  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: wshcutinit      
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: gmean           
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: tgmean           
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: wc_frac              
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: wgn             
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: tasum           
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: loss            
  ! perte en C lors de la fauche
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: lossc                 
  ! perte en N lors de la fauche
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: lossn                 
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: tlossstart         
  INTEGER , ALLOCATABLE, SAVE, DIMENSION(:,:)   :: flag_fertilisation
  INTEGER, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: fertcount
  ! C/N dans le r�servoir stucturel  = 150
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: c2nratiostruct        
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: nfertammtot
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: nfertnittot
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: nfertammtotyear
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: nfertnittotyear
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: nfertammtotprevyear
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: nfertnittotprevyear
  ! metabolic C in slurry and manure (kg C/m**2/d)
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: fcOrganicFertmetabolic 
  ! structural C in slurry and manure (kg C/m**2/d)
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: fcOrganicFertstruct   
  ! urine N in slurry and manure (kg N/m**2/d)
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: fnOrganicFerturine     
  ! metabolic N in slurry and manure (kg N/m**2/d)
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: fnOrganicFertmetabolic            
  ! struct N in slurry and manure (kg N/m**2/d)
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: fnOrganicFertstruct
  ! urine P in slurry and manure (kg N/m**2/d)
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: fpOrganicFerturine
  ! metabolic P in slurry and manure (kg N/m**2/d)
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: fpOrganicFertmetabolic
  ! struct P in slurry and manure (kg N/m**2/d)
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: fpOrganicFertstruct

  ! variables pour l'auto gestion de nicolas
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: nsatur_somerror_temp
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: nsatur_somerror
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: tfert_modif
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: nnonlimit_SOMerror
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: nnonlimit_SOMerrormax
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: controle_azote_sum_mem
  INTEGER , ALLOCATABLE, SAVE, DIMENSION(:,:)   :: n_auto
  INTEGER , ALLOCATABLE, SAVE, DIMENSION(:,:)   :: stoplimitant
  INTEGER , ALLOCATABLE, SAVE, DIMENSION(:,:)   :: fertcount_start
  INTEGER , ALLOCATABLE, SAVE, DIMENSION(:,:)   :: fertcount_current
  LOGICAL   , ALLOCATABLE, SAVE, DIMENSION(:,:)   :: fertil_year
  REAL, ALLOCATABLE, SAVE, DIMENSION(:)   :: toto
  REAL,ALLOCATABLE, SAVE, DIMENSION(:,:)   :: wshtotsumprevyear
  REAL, ALLOCATABLE, SAVE, DIMENSION(:)   :: tmp_sr_ugb_C3
  REAL, ALLOCATABLE, SAVE, DIMENSION(:)   :: tmp_nb_ani_C3
  REAL, ALLOCATABLE, SAVE, DIMENSION(:)   :: tmp_grazed_frac_C3
  REAL, ALLOCATABLE, SAVE, DIMENSION(:)   :: tmp_import_yield_C3
  REAL, ALLOCATABLE, SAVE, DIMENSION(:)   :: tmp_wshtotsum_C3
  REAL, ALLOCATABLE, SAVE, DIMENSION(:)   :: tmp_sr_ugb_C4
  REAL, ALLOCATABLE, SAVE, DIMENSION(:)   :: tmp_nb_ani_C4
  REAL, ALLOCATABLE, SAVE, DIMENSION(:)   :: tmp_grazed_frac_C4
  REAL, ALLOCATABLE, SAVE, DIMENSION(:)   :: tmp_import_yield_C4
  REAL, ALLOCATABLE, SAVE, DIMENSION(:)   :: tmp_wshtotsum_C4

  REAL,ALLOCATABLE, SAVE, DIMENSION(:,:)   :: DM_cutyearly
  REAL,ALLOCATABLE, SAVE, DIMENSION(:,:)   :: C_cutyearly
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: YIELD_RETURN
  REAL, ALLOCATABLE, SAVE, DIMENSION(:)   :: sr_ugb_init

  INTEGER                  , SAVE                 :: cut_year = 10
  INTEGER                  , SAVE                 :: compt_fert
  INTEGER                  , SAVE                 :: min_fert
  INTEGER                  , SAVE                 :: fert_max
  INTEGER                  , SAVE                 :: i_compt
  REAL                     , SAVE                 :: deltatt             ! = 0.02
  ! couter number of years of simulation     
  INTEGER                  , SAVE                 :: count_year            
  INTEGER                  , SAVE                 :: year_count1
  ! couter number of years of simulation
  INTEGER                  , SAVE                 :: year_count2
  ! couter number of years of simulation

  CHARACTER(len=500), ALLOCATABLE,SAVE,DIMENSION (:)     :: file_management
  CHARACTER(len=500)              , SAVE                 :: file_param_init
  CHARACTER(len=500)              , SAVE                 :: file_import_yield

  INTEGER                  , SAVE                 :: Type_animal
  INTEGER                  , SAVE                 :: mcut_C3
  INTEGER                  , SAVE                 :: mauto_C3
  INTEGER                  , SAVE                 :: mcut_C4
  INTEGER                  , SAVE                 :: mauto_C4
  ! yearly total azote by fertilization
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: apport_azote
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: trampling

  ! new variables for get map of management
  INTEGER                  , SAVE                 :: f_management_map 
  CHARACTER(len=500)              , SAVE                 :: management_map
  CHARACTER(len=500)              , SAVE                 :: fertility_map

  INTEGER                  , SAVE                 :: f_deposition_map
  CHARACTER(len=500)              , SAVE                 :: deposition_map
  INTEGER                  , SAVE                 :: f_grazing_map
  CHARACTER(len=500)              , SAVE                 :: grazing_map
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: ndeposition
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: N_fert_total
!  REAL                     , SAVE                 :: N_effect
  REAL,ALLOCATABLE,    SAVE , DIMENSION(:,:) :: compt_cut
  REAL,ALLOCATABLE,    SAVE , DIMENSION(:,:) :: frequency_cut
  REAL,ALLOCATABLE,    SAVE , DIMENSION(:,:) :: sr_wild
  INTEGER,ALLOCATABLE,    SAVE , DIMENSION(:,:) :: flag_cutting
  ! from applic_plant
  REAL, DIMENSION(:), ALLOCATABLE, SAVE :: tamean1
  REAL, DIMENSION(:), ALLOCATABLE, SAVE :: tamean2
  REAL, DIMENSION(:), ALLOCATABLE, SAVE :: tamean3
  REAL, DIMENSION(:), ALLOCATABLE, SAVE :: tamean4
  REAL, DIMENSION(:), ALLOCATABLE, SAVE :: tamean5
  REAL, DIMENSION(:), ALLOCATABLE, SAVE :: tamean6
  REAL, DIMENSION(:), ALLOCATABLE, SAVE :: tameand
  REAL, DIMENSION(:), ALLOCATABLE, SAVE :: tameanw
  REAL, DIMENSION(:), ALLOCATABLE, SAVE :: tacumm
  REAL, DIMENSION(:), ALLOCATABLE, SAVE :: tacummprev
  REAL, DIMENSION(:), ALLOCATABLE, SAVE :: tsoilcumm
  REAL, DIMENSION(:), ALLOCATABLE, SAVE :: tsoilcummprev
  REAL, DIMENSION(:), ALLOCATABLE, SAVE :: tsoilmeand
  REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: tcut0

  ! flag for verify that without twice fertilisation at same time
  LOGICAL     , ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: tfert_verif
  LOGICAL     , ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: tfert_verif2  ! idem
  LOGICAL     , ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: tfert_verif3  ! idem

  REAL,ALLOCATABLE,    SAVE , DIMENSION(:,:) :: Fert_sn
  REAL,ALLOCATABLE,    SAVE , DIMENSION(:,:) :: Fert_on
  REAL,ALLOCATABLE,    SAVE , DIMENSION(:,:) :: Fert_PRP
  ! JCADD 20170421 for CNP merge
  REAL,ALLOCATABLE,    SAVE , DIMENSION(:,:,:,:) :: losscnp
  REAL,ALLOCATABLE,    SAVE , DIMENSION(:,:,:) :: grassharvestcnp
  REAL,ALLOCATABLE,    SAVE , DIMENSION(:,:,:) :: intakecnp
  REAL,ALLOCATABLE,    SAVE , DIMENSION(:,:,:,:,:) :: litter_intakecnp
  REAL,ALLOCATABLE,    SAVE , DIMENSION(:,:,:,:) :: tramplingcnp
  REAL,ALLOCATABLE,    SAVE , DIMENSION(:,:,:) :: faecescnp
  REAL,ALLOCATABLE,    SAVE , DIMENSION(:,:,:) :: urinecnp
  REAL,ALLOCATABLE,    SAVE , DIMENSION(:,:,:) :: milkcnp
  REAL,ALLOCATABLE,    SAVE , DIMENSION(:,:,:) :: liveweightcnp
  REAL,ALLOCATABLE,    SAVE , DIMENSION(:,:) :: ranimal_c
  REAL,ALLOCATABLE,    SAVE , DIMENSION(:,:) :: EFch4_c
  REAL,ALLOCATABLE,    SAVE , DIMENSION(:,:,:) :: grassorganicfertcnp
  REAL,ALLOCATABLE,    SAVE , DIMENSION(:,:,:) :: pfert_mineral
  REAL, SAVE :: flam_left ! fraction of leaf in the residual biomass after harvest
  ! END JCADD

CONTAINS

  SUBROUTINE main_grassland_management(&
     npts ,lalo, neighbours, resolution, contfrac, &
     dt, tjulian, t2m_daily, t2m_min_daily, &
     t2m_14, tsoil, snowfall_daily, biomass, bm_to_litter, &
     litter, litter_avail, litter_not_avail , &
!     !spitfire
!     fuel_1hr, fuel_10hr, &
!     fuel_100hr, fuel_1000hr, &
!     !end spitfire
     new_day, new_year, when_growthinit_cut, nb_grazingdays,&
     lai, sla_calc, leaf_age, leaf_frac, &
     wshtotsum, sr_ugb, compt_ugb, &
     nb_ani, grazed_frac, import_yield, &
!gmjc top 5 layer grassland soil moisture for grazing
     moiavail_daily, tmc_topgrass_daily,fc_grazing, snowmass_daily,&
     after_snow, after_wet, wet1day, wet2day, &
!end gmjc 
     harvest_gm, ranimal_gm, ch4_pft_gm, cinput_gm, n2o_pft_gm, &
     ! JCADD 20170421 for CNP merge
     grm_nfert, grm_pfert,&
     ! variable for mass balance calculation
     grm_masschangecnp,GRM_devstage,when_growthinit)
     ! END JCADD

    INTEGER                                , INTENT(in)   :: npts   
    INTEGER,DIMENSION(npts,8),INTENT(in) :: neighbours        !!Neighoring grid points if land for the DGVM
                                                                         !!(unitless)
    REAL,DIMENSION(npts,2),INTENT(in)    :: lalo              !!Geographical coordinates (latitude,longitude)
    REAL,DIMENSION(npts,2),INTENT(in)    :: resolution        !! Size in x an y of the grid (m) - surface area of
                                                                         !! the gridbox
    REAL,DIMENSION (npts), INTENT (in)   :: contfrac          !!Fraction of continent in the grid cell (unitless)
    REAL                             , INTENT(in)   :: dt          
    INTEGER                             , INTENT(in)   :: tjulian 
    ! julien day
    REAL, DIMENSION(npts)            , INTENT(in)   :: t2m_daily      
    ! air temperature 
    REAL, DIMENSION(npts)            , INTENT(in)   ::  t2m_min_daily
    ! daily minimum temperature
    REAL, DIMENSION(npts)            , INTENT(in)   ::  t2m_14   
    ! 14 days mean temperature
    REAL, DIMENSION(npts)            , INTENT(in)   :: tsoil     
    ! soil surface t (k)
    REAL, DIMENSION(npts)            , INTENT(in)   :: snowfall_daily       
    ! snow fall (mm/d)
    REAL, DIMENSION(npts)            , INTENT(in)   :: snowmass_daily
    ! snow mass (kg/m2)
    REAL, DIMENSION(npts,nvm,nparts,nelements), INTENT(inout):: biomass       
    REAL, DIMENSION(npts,nvm,nparts,nelements), INTENT(inout):: bm_to_litter 
    ! conv of biomass to litter (gC/(m**2/agri ground)) / day
    REAL, DIMENSION(npts,nlitt,nvm,nlevs,nelements), INTENT(inout):: litter
    REAL, DIMENSION(npts,nlitt,nvm), INTENT(inout):: litter_avail
    REAL, DIMENSION(npts,nlitt,nvm) , INTENT(inout):: litter_not_avail
    LOGICAL                                , INTENT(in)   :: new_day   
    ! flag indicate new day
    LOGICAL                                , INTENT(in)   :: new_year   
    ! flag indicate new year
    REAL, DIMENSION(npts,nvm), INTENT(inout)       :: when_growthinit_cut
    ! how many days ago was the beginning of the last cut
    REAL, DIMENSION(npts,nvm), INTENT(inout)       :: nb_grazingdays
    REAL, DIMENSION(npts,nvm)       , INTENT(inout)  :: lai        
    ! leaf area index OF AN INDIVIDUAL PLANT
!    !spitfire
!    REAL, DIMENSION(npts,nvm,nlitt),INTENT(inout)        :: fuel_1hr
!    REAL, DIMENSION(npts,nvm,nlitt),INTENT(inout)        :: fuel_10hr
!    REAL, DIMENSION(npts,nvm,nlitt),INTENT(inout)        :: fuel_100hr
!    REAL, DIMENSION(npts,nvm,nlitt),INTENT(inout)        :: fuel_1000hr
!    !end spitfire
    REAL, DIMENSION(npts,nvm)       , INTENT(in)  :: sla_calc
    REAL, DIMENSION(npts,nvm,nparts,nleafages)       , INTENT(inout)  :: leaf_frac
    REAL, DIMENSION(npts,nvm,nparts,nleafages)       , INTENT(inout)  :: leaf_age
    REAL, DIMENSION(npts,nvm), INTENT(inout)                 :: wshtotsum
    REAL, DIMENSION(npts,nvm), INTENT(inout)  ::  sr_ugb
    REAL, DIMENSION(npts,nvm), INTENT(inout)  ::  compt_ugb
    REAL, DIMENSION(npts,nvm), INTENT(inout)  ::  nb_ani
    REAL, DIMENSION(npts,nvm), INTENT(inout)  ::  grazed_frac
    REAL, DIMENSION(npts,nvm), INTENT(inout)  ::  import_yield
!gmjc top 5 layer grassland soil moisture for grazing
    REAL, DIMENSION(npts,nvm), INTENT(in)           :: moiavail_daily
    REAL,DIMENSION (npts), INTENT(in)       :: tmc_topgrass_daily
    REAL,DIMENSION (npts), INTENT(in)       :: fc_grazing
    REAL,DIMENSION (npts), INTENT(inout)    :: after_snow
    REAL,DIMENSION (npts), INTENT(inout)    :: after_wet
    REAL,DIMENSION (npts), INTENT(inout)    :: wet1day
    REAL,DIMENSION (npts), INTENT(inout)    :: wet2day
!end gmjc
    REAL, DIMENSION(npts,nvm), INTENT(inout)  :: harvest_gm
    REAL, DIMENSION(npts,nvm), INTENT(inout)  :: ranimal_gm
    REAL, DIMENSION(npts,nvm), INTENT(inout)  :: ch4_pft_gm
    REAL, DIMENSION(npts,nvm), INTENT(inout)  :: cinput_gm
    REAL, DIMENSION(npts,nvm), INTENT(inout)  :: n2o_pft_gm
    ! JCADD 20170421 for CNP merge
    REAL, DIMENSION(npts,nvm,ninput), INTENT(out)  :: grm_nfert 
    REAL, DIMENSION(npts,nvm), INTENT(out)  :: grm_pfert
    REAL, DIMENSION(npts,nvm,nelements), INTENT(out)  :: grm_masschangecnp
    REAL, DIMENSION(npts,nvm), INTENT(inout)  :: GRM_devstage
    REAL, DIMENSION(npts,nvm), INTENT(in)  :: when_growthinit
    ! END JCADD

    LOGICAL :: l_error = .FALSE.
    INTEGER :: ier, i, j, k,h, m,l
    REAL, DIMENSION(npts)        :: xtmp_npts
    REAL, DIMENSION(npts,ngmean) :: xtmp_npts_3d
    REAL, DIMENSION(npts,nvm)        :: regcount_real
    REAL, DIMENSION(npts,nvm)        :: fertcount_real
    INTEGER :: fertcount_next
    REAL :: intakemax_t
    REAL :: wanimal_t
    REAL, DIMENSION(ncut)        ::wshtotcutinit_t
    REAL, DIMENSION(npts,nvm)        :: lm_before
    REAL, DIMENSION(npts,nvm)        :: lm_after
    REAL, DIMENSION(npts,nvm)        :: bm_cut

    REAL, PARAMETER       :: n2o_EF1 = 0.01
    REAL, PARAMETER       :: n2o_EF2 = 0.015
    REAL, PARAMETER       :: n2o_EF3 = 0.01
    REAL, PARAMETER       :: n2o_EF4 = 0.0075
    REAL, PARAMETER       :: n2o_FracGASF = 0.10
    REAL, PARAMETER       :: n2o_FracGASM = 0.20
    REAL, PARAMETER       :: n2o_FracLEACH_H = 0.30

    CHARACTER(LEN=2), DIMENSION(nelements)                     :: element_str
!! string suffix indicating element

!    REAL, DIMENSION(npts,nvm)        :: N_fert_total
    REAL :: fertility_legume_t

    ! 1. initialisations
    init_grassland : IF (first_call_grassland_manag) THEN

      first_call_grassland_manag = .FALSE. 

      ! 1.1 allocate variables

      ALLOCATE (intake                (npts,nvm)          , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'intake', 'Not enough memory', '')
      ALLOCATE (intakemax             (npts,nvm)          , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'intakemax', 'Not enough memory', '')
      ALLOCATE (intake_litter         (npts,nvm)          , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'intake_litter', 'Not enough memory', '')
      ALLOCATE (intake_animal_litter  (npts,nvm)          , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'intake_animal_litter', 'Not enough memory', '')
      ALLOCATE (grazing_litter        (npts,nvm)          , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'grazing_litter', 'Not enough memory', '')
      ALLOCATE (litter_avail_totDM    (npts,nvm)          , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'litter_avail_totDM', 'Not enough memory', '')
      ALLOCATE (wshtotcutinit         (npts,nvm,ncut)     , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'wshtotcutinit', 'Not enough memory', '')
      ALLOCATE (lcutinit              (npts,nvm,ncut)     , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'lcutinit', 'Not enough memory', '')
      ALLOCATE (devstage              (npts,nvm)          , stat=ier)   
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'devstage', 'Not enough memory', '')
      ALLOCATE (faecesc               (npts,nvm)          , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'faecesc', 'Not enough memory', '')
      ALLOCATE (faecesn               (npts,nvm)          , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'faecesn', 'Not enough memory', '')
      ALLOCATE (urinen                (npts,nvm)          , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'urinen', 'Not enough memory', '')
      ALLOCATE (urinec                (npts,nvm)          , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'urinec', 'Not enough memory', '')
      ALLOCATE (nel                   (npts,nvm)          , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'nel', 'Not enough memory', '')
      ALLOCATE (nanimaltot            (npts,nvm)          , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'nanimaltot', 'Not enough memory', '')
      ALLOCATE (tgrowth               (npts,nvm)          , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'tgrowth', 'Not enough memory', '')
      ALLOCATE (wsh                   (npts,nvm)          , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'wsh', 'Not enough memory', '')
      ALLOCATE (wshtot                (npts,nvm)          , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'wshtot', 'Not enough memory', '')
      ALLOCATE (wshtotinit            (npts,nvm)          , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'wshtotinit', 'Not enough memory', '')
      ALLOCATE (wr                    (npts,nvm)          , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'wr', 'Not enough memory', '')
      ALLOCATE (wrtot                 (npts,nvm)          , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'wrtot', 'Not enough memory', '')
      ALLOCATE (wanimal               (npts,nvm)          , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'wanimal', 'Not enough memory', '')
      ALLOCATE (ntot                  (npts,nvm)          , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'ntot', 'Not enough memory', '')
      ALLOCATE (c                     (npts,nvm)          , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'c', 'Not enough memory', '')
      ALLOCATE (n                     (npts,nvm)          , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'n', 'Not enough memory', '')
      ALLOCATE (fn                    (npts,nvm)          , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'fn', 'Not enough memory', '')
      ALLOCATE (napo                  (npts,nvm)          , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'napo', 'Not enough memory', '')
      ALLOCATE (nsym                  (npts,nvm)          , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'nsym', 'Not enough memory', '')
      ALLOCATE (wnapo                 (npts,nvm)          , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'wnapo', 'Not enough memory', '')
      ALLOCATE (wnsym                 (npts,nvm)          , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'wnsym', 'Not enough memory', '')
      ALLOCATE (wn                    (npts,nvm)          , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'wn', 'Not enough memory', '')
      ALLOCATE (nanimal               (npts,nvm,nstocking), stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'nanimal', 'Not enough memory', '')
      ALLOCATE (tanimal               (npts,nvm,nstocking), stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'tanimal', 'Not enough memory', '')
      ALLOCATE (danimal               (npts,nvm,nstocking), stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'danimal', 'Not enough memory', '')
      ALLOCATE (tcut                  (npts,nvm,nstocking), stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'tcut', 'Not enough memory', '')
      ALLOCATE (tfert                 (npts,nvm,nstocking), stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'tfert', 'Not enough memory', '')
      ALLOCATE (Nliquidmanure         (npts,nvm,nstocking), stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'Nliquidmanure', 'Not enough memory', '')
      ALLOCATE (nslurry               (npts,nvm,nstocking), stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'nslurry', 'Not enough memory', '')
      ALLOCATE (Nsolidmanure          (npts,nvm,nstocking), stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'Nsolidmanure', 'Not enough memory', '')
      ALLOCATE (legume_fraction       (npts,nvm)          , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'legume_fraction', 'Not enough memory', '')
      ALLOCATE (soil_fertility        (npts,nvm)          , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'soil_fertility', 'Not enough memory', '')
      ALLOCATE (Animalwgrazingmin     (npts,nvm)          , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'Animalwgrazingmin', 'Not enough memory', '')
      ALLOCATE (AnimalkintakeM        (npts,nvm)          , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'AnimalkintakeM', 'Not enough memory', '')
      ALLOCATE (AnimalDiscremineQualite (npts,nvm)        , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'AnimalDiscremineQualite', 'Not enough memory', '')
      ALLOCATE (controle_azote        (npts,nvm,nstocking), stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'controle_azote', 'Not enough memory', '')
      ALLOCATE (fcOrganicFertmetabolicsum (npts,nvm)      , stat=ier)  
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'fcOrganicFertmetabolicsum', 'Not enough memory', '')
      ALLOCATE (fcOrganicFertstructsum (npts,nvm)         , stat=ier)  
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'fcOrganicFertstructsum', 'Not enough memory', '')
      ALLOCATE (fnOrganicFertmetabolicsum (npts,nvm)      , stat=ier)  
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'fnOrganicFertmetabolicsum', 'Not enough memory', '')
      ALLOCATE (fnOrganicFertstructsum (npts,nvm)         , stat=ier)  
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'fnOrganicFertstructsum', 'Not enough memory', '')
      ALLOCATE (fnOrganicFerturinesum (npts,nvm)          , stat=ier)  
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'fnOrganicFerturinesum', 'Not enough memory', '')
      ALLOCATE (fnatmsum              (npts,nvm)          , stat=ier)  
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'fnatmsum', 'Not enough memory', '')
      ALLOCATE (controle_azote_sum    (npts,nvm)          , stat=ier)  
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'controle_azote_sum', 'Not enough memory', '')
      ALLOCATE (nfertamm              (npts,nvm,nstocking), stat=ier) 
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'nfertamm', 'Not enough memory', '')
      ALLOCATE (nfertnit              (npts,nvm,nstocking), stat=ier) 
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'nfertnit', 'Not enough memory', '')
      ALLOCATE (intakesum             (npts,nvm)          , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'intakesum', 'Not enough memory', '')
      ALLOCATE (intakensum            (npts,nvm)          , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'intakensum', 'Not enough memory', '')
      ALLOCATE (intake_animal         (npts,nvm)          , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'intake_animal', 'Not enough memory', '')
      ALLOCATE (intake_animalsum      (npts,nvm)          , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'intake_animalsum', 'Not enough memory', '')
      ALLOCATE (PIYcow                (npts,nvm,nstocking), stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'PIYcow', 'Not enough memory', '')
      ALLOCATE (PIMcow                (npts,nvm,nstocking), stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'PIMcow', 'Not enough memory', '')
      ALLOCATE (BCSYcow               (npts,nvm,nstocking), stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'BCSYcow', 'Not enough memory', '')
      ALLOCATE (BCSMcow               (npts,nvm,nstocking), stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'BCSMcow', 'Not enough memory', '')
      ALLOCATE (PICcow                (npts,nvm,nstocking), stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'PICcow', 'Not enough memory', '')
      ALLOCATE (AGE_cow_P             (npts,nvm,nstocking), stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'AGE_cow_P', 'Not enough memory', '')
      ALLOCATE (AGE_cow_M             (npts,nvm,nstocking), stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'AGE_cow_M', 'Not enough memory', '')
      ALLOCATE (Autogestion_out       (npts,nvm,n_out)    , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'Autogestion_out', 'Not enough memory', '')
      ALLOCATE (Forage_quantity       (npts,nvm,nstocking), stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'Forage_quantity', 'Not enough memory', '')
      ALLOCATE (tcut_modif            (npts,nvm,nstocking), stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'tcut_modif', 'Not enough memory', '')
      ALLOCATE (countschedule         (npts,nvm)          , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'countschedule', 'Not enough memory', '')
      ALLOCATE (mux                   (npts,nvm)          , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'mux', 'Not enough memory', '')
      ALLOCATE (mugmean               (npts,nvm)          , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'mugmean', 'Not enough memory', '')
      ALLOCATE (sigx                  (npts,nvm)          , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'sigx', 'Not enough memory', '')
      ALLOCATE (sigy                  (npts,nvm)          , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'sigy', 'Not enough memory', '')
      ALLOCATE (gmeanslope            (npts,nvm)          , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'gmeanslope', 'Not enough memory', '')
      ALLOCATE (gzero                 (npts,nvm)          , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'gzero', 'Not enough memory', '')
      ALLOCATE (gcor                  (npts,nvm)          , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'gcor', 'Not enough memory', '')
      ALLOCATE (cuttingend            (npts,nvm)          , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'cuttingend', 'Not enough memory', '')
      ALLOCATE (tcut_verif            (npts,nvm,nstocking), stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'tcut_verif', 'Not enough memory', '')
      ALLOCATE (tfert_verif           (npts,nvm,nstocking), stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'tfert_verif', 'Not enough memory', '')
      ALLOCATE (tfert_verif2          (npts,nvm,nstocking), stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'tfert_verif2', 'Not enough memory', '')
      ALLOCATE (tfert_verif3          (npts,nvm,nstocking), stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'tfert_verif3', 'Not enough memory', '')
      ALLOCATE (regcount              (npts,nvm)          , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'regcount', 'Not enough memory', '')
!      ALLOCATE (wshcutinit            (npts,nvm)          , stat=ier)
!      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'wshcutinit', 'Not enough memory', '')
      ALLOCATE (gmean                 (npts,nvm,ngmean)   , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'gmean', 'Not enough memory', '')
      ALLOCATE (tgmean                (npts,nvm,ngmean)   , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'tgmean', 'Not enough memory', '')
      ALLOCATE (wc_frac               (npts,nvm)          , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'wc_frac', 'Not enough memory', '')
      ALLOCATE (wgn                   (npts,nvm)          , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'wgn', 'Not enough memory', '')
      ALLOCATE (tasum                 (npts,nvm)          , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'tasum', 'Not enough memory', '')
      ALLOCATE (loss                  (npts,nvm)          , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'loss', 'Not enough memory', '')
      ALLOCATE (lossc                 (npts,nvm)          , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'lossc', 'Not enough memory', '')
      ALLOCATE (lossn                 (npts,nvm)          , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'lossn', 'Not enough memory', '')
      ALLOCATE (tlossstart            (npts,nvm)          , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'tlossstart', 'Not enough memory', '')
      ALLOCATE (flag_fertilisation    (npts,nvm)          , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'flag_fertilisation', 'Not enough memory', '')
      ALLOCATE (fertcount             (npts,nvm)          , stat=ier) 
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'fertcount', 'Not enough memory', '')
      ALLOCATE (c2nratiostruct        (npts,nvm)          , stat=ier)  
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'c2nratiostruct', 'Not enough memory', '')
      ALLOCATE (nfertammtot           (npts,nvm)          , stat=ier)  
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'nfertammtot', 'Not enough memory', '')
      ALLOCATE (nfertnittot           (npts,nvm)          , stat=ier)  
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'nfertnittot', 'Not enough memory', '')
      ALLOCATE (nfertammtotyear       (npts,nvm)          , stat=ier)   
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'nfertammtotyear', 'Not enough memory', '')
      ALLOCATE (nfertnittotyear       (npts,nvm)          , stat=ier)   
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'nfertnittotyear', 'Not enough memory', '')
      ALLOCATE (nfertammtotprevyear   (npts,nvm)          , stat=ier)   
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'nfertammtotprevyear', 'Not enough memory', '')
      ALLOCATE (nfertnittotprevyear   (npts,nvm)          , stat=ier)   
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'nfertnittotprevyear', 'Not enough memory', '')
      ALLOCATE (fcOrganicFertmetabolic (npts,nvm)         , stat=ier)  
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'fcOrganicFertmetabolic', 'Not enough memory', '')
      ALLOCATE (fcOrganicFertstruct   (npts,nvm)          , stat=ier)  
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'fcOrganicFertstruct', 'Not enough memory', '')
      ALLOCATE (fnOrganicFerturine    (npts,nvm)          , stat=ier)  
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'fnOrganicFerturine', 'Not enough memory', '')
      ALLOCATE (fnOrganicFertstruct   (npts,nvm)          , stat=ier)  
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'fnOrganicFertstruct', 'Not enough memory', '')
      ALLOCATE (fnOrganicFertmetabolic (npts,nvm)         , stat=ier)  
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'fnOrganicFertmetabolic', 'Not enough memory', '')
      ALLOCATE (fpOrganicFerturine    (npts,nvm)          , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'fpOrganicFerturine', 'Not enough memory', '')
      ALLOCATE (fpOrganicFertstruct   (npts,nvm)          , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'fpOrganicFertstruct', 'Not enough memory', '')
      ALLOCATE (fpOrganicFertmetabolic (npts,nvm)         , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'fpOrganicFertmetabolic', 'Not enough memory', '')
      ALLOCATE (nsatur_somerror_temp  (npts,nvm)          , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'nsatur_somerror_temp', 'Not enough memory', '')
      ALLOCATE (nsatur_somerror       (npts,nvm)          , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'nsatur_somerror', 'Not enough memory', '')
      ALLOCATE (tfert_modif           (npts,nvm,nstocking), stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'tfert_modif', 'Not enough memory', '')
      ALLOCATE (nnonlimit_SOMerror    (npts,nvm)          , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'nnonlimit_SOMerror', 'Not enough memory', '')
      ALLOCATE (nnonlimit_SOMerrormax (npts,nvm)          , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'nnonlimit_SOMerrormax', 'Not enough memory', '')
      ALLOCATE (controle_azote_sum_mem (npts,nvm)         , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'controle_azote_sum_mem', 'Not enough memory', '')
      ALLOCATE (n_auto                (npts,nvm)          , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'n_auto', 'Not enough memory', '')
      ALLOCATE (stoplimitant          (npts,nvm)          , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'stoplimitant', 'Not enough memory', '')
      ALLOCATE (fertcount_start       (npts,nvm)          , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'fertcount_start', 'Not enough memory', '')
      ALLOCATE (fertcount_current     (npts,nvm)          , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'fertcount_current', 'Not enough memory', '')
      ALLOCATE (wshtotsumprev         (npts,nvm)          , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'wshtotsumprev', 'Not enough memory', '')
      ALLOCATE (fertil_year           (npts,nvm)          , stat=ier)  
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'fertil_year', 'Not enough memory', '')
      ALLOCATE (toto                  (npts)              , stat=ier)  
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'toto', 'Not enough memory', '')
      ALLOCATE (apport_azote          (npts,nvm)          , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'apport_azote', 'Not enough memory', '')
      ALLOCATE (trampling             (npts,nvm)          , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'trampling', 'Not enough memory', '')
      ALLOCATE (wshtotsumprevyear     (npts,nvm)          , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'wshtotsumprevyear', 'Not enough memory', '')
      ALLOCATE(file_management        (nvm)               , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'file_management', 'Not enough memory', '')
      ALLOCATE (tmp_sr_ugb_C3         (npts)              , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'tmp_sr_ugb_C3', 'Not enough memory', '')
      ALLOCATE (tmp_nb_ani_C3         (npts)              , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'tmp_nb_ani_C3', 'Not enough memory', '')
      ALLOCATE (tmp_grazed_frac_C3    (npts)              , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'tmp_grazed_frac_C3', 'Not enough memory', '')
      ALLOCATE (tmp_import_yield_C3   (npts)              , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'tmp_import_yield_C3', 'Not enough memory', '')
      ALLOCATE (tmp_wshtotsum_C3      (npts)              , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'tmp_wshtotsum_C3', 'Not enough memory', '')
      ALLOCATE (tmp_sr_ugb_C4         (npts)              , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'tmp_sr_ugb_C4', 'Not enough memory', '')
      ALLOCATE (tmp_nb_ani_C4         (npts)              , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'tmp_nb_ani_C4', 'Not enough memory', '')
      ALLOCATE (tmp_grazed_frac_C4    (npts)              , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'tmp_grazed_frac_C4', 'Not enough memory', '')
      ALLOCATE (tmp_import_yield_C4   (npts)              , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'tmp_import_yield_C4', 'Not enough memory', '')
      ALLOCATE (tmp_wshtotsum_C4      (npts)              , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'tmp_wshtotsum_C4', 'Not enough memory', '')
      ALLOCATE (DM_cutyearly          (npts,nvm)          , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'DM_cutyearly', 'Not enough memory', '')
      ALLOCATE (C_cutyearly           (npts,nvm)          , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'C_cutyearly', 'Not enough memory', '')
      ALLOCATE (YIELD_RETURN          (npts,nvm)          , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'YIELD_RETURN', 'Not enough memory', '')
      ALLOCATE (sr_ugb_init           (npts)              , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'sr_ugb_init', 'Not enough memory', '')
      ALLOCATE (N_fert_total          (npts,nvm)          , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'N_fert_total', 'Not enough memory', '')
      ALLOCATE (ndeposition           (npts,nvm)          , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'ndeposition', 'Not enough memory', '')
      ALLOCATE (compt_cut             (npts,nvm)          , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'compt_cut', 'Not enough memory', '')
      ALLOCATE (frequency_cut         (npts,nvm)          , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'frequency_cut', 'Not enough memory', '')
      ALLOCATE (sr_wild               (npts,nvm)          , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'sr_wild', 'Not enough memory', '')
      ALLOCATE (flag_cutting          (npts,nvm)          , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'flag_cutting', 'Not enough memory', '')


      ! JCADD 20170421 for CNP merge
      ALLOCATE (losscnp          (npts,nvm,nparts,nelements)          , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'losscnp', 'Not enough memory', '')
      ALLOCATE (grassharvestcnp          (npts,nvm,nelements)          , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'grassharvestcnp', 'Not enough memory', '')
      ALLOCATE (intakecnp          (npts,nvm,nelements)          , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'intakecnp', 'Not enough memory', '')
      ALLOCATE (litter_intakecnp          (npts,nlitt,nvm,nlevs,nelements)   , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'litter_intakecnp', 'Not enough memory', '')
      ALLOCATE (tramplingcnp          (npts,nvm,nparts,nelements)          , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'tramplingcnp', 'Not enough memory', '')
      ALLOCATE (faecescnp          (npts,nvm,nelements)          , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'faecescnp', 'Not enough memory', '')
      ALLOCATE (urinecnp          (npts,nvm,nelements)          , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'urinecnp', 'Not enough memory', '')
      ALLOCATE (milkcnp          (npts,nvm,nelements)          , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'milkcnp', 'Not enough memory', '')
      ALLOCATE (liveweightcnp          (npts,nvm,nelements)          , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'liveweightcnp', 'Not enough memory', '')
      ALLOCATE (ranimal_c         (npts,nvm)          , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'ranimal_c', 'Not enough memory', '')
      ALLOCATE (EFch4_c          (npts,nvm)          , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'EFch4_c', 'Not enough memory', '')
      ALLOCATE (grassorganicfertcnp          (npts,nvm,nelements)          , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'grassorganicfertcnp', 'Not enough memory', '')
      ALLOCATE (pfert_mineral          (npts,nvm,nstocking)          , stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_GM', 'pfert_mineral', 'Not enough memory', '')
      ! END JCADD

      ! from applic_plant
      ALLOCATE (tamean1        (npts), stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_appl_pre_animal', 'tamean1', 'Not enough memory', '')
      ALLOCATE (tamean2        (npts), stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_appl_pre_animal', 'tamean2', 'Not enough memory', '')
      ALLOCATE (tamean3        (npts), stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_appl_pre_animal', 'tamean3', 'Not enough memory', '')
      ALLOCATE (tamean4        (npts), stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_appl_pre_animal', 'tamean4', 'Not enough memory', '')
      ALLOCATE (tamean5        (npts), stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_appl_pre_animal', 'tamean5', 'Not enough memory', '')
      ALLOCATE (tamean6        (npts), stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_appl_pre_animal', 'tamean6', 'Not enough memory', '')
      ALLOCATE (tameand        (npts), stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_appl_pre_animal', 'tameand', 'Not enough memory', '')
      ALLOCATE (tameanw        (npts), stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_appl_pre_animal', 'tameanw', 'Not enough memory', '')
      ALLOCATE (tacumm         (npts), stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_appl_pre_animal', 'tacumm', 'Not enough memory', '')
      ALLOCATE (tacummprev     (npts), stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_appl_pre_animal', 'tacummprev', 'Not enough memory', '')
      ALLOCATE (tsoilcumm      (npts), stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_appl_pre_animal', 'tsoilcumm', 'Not enough memory', '')
      ALLOCATE (tsoilcummprev  (npts), stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_appl_pre_animal', 'tsoilcummprev', 'Not enough memory', '')
      ALLOCATE (tsoilmeand     (npts), stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_appl_pre_animal', 'tsoilmeand', 'Not enough memory', '')
      ALLOCATE (tcut0          (npts,nvm), stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_appl_pre_animal', 'tcut0', 'Not enough memory', '')
      ALLOCATE (Fert_sn          (npts,nvm), stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_appl_pre_animal', 'Fert_sn', 'Not enough memory', '')
      ALLOCATE (Fert_on          (npts,nvm), stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_appl_pre_animal', 'Fert_on', 'Not enough memory', '')
      ALLOCATE (Fert_PRP          (npts,nvm), stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'Main_appl_pre_animal', 'Fert_PRP', 'Not enough memory', '')

      ! 1.2 set flags and variables need to read in Pasim

      flam_left = 0.25
      CALL getin_p('FLAM_LEFT',flam_left)
      WRITE(numout,*)  'FLAM_LEFT',flam_left
 
      ! saturant N supply
      f_saturant = 0
      CALL getin_p('GRM_F_SATURANT',f_saturant)
      ! N fertilization without limitation
      f_nonlimitant = 0
      CALL getin_p('GRM_F_NONLIMITANT',f_nonlimitant)
      ! f_autogestion = 1-5
      ! 1: auto cut for PFT m_auto
      ! 2: auto graze for PFT m_auto
      ! 3: auto cut and graze for PFT m_cut and m_grazed with increasing sr_ugb
      ! 4: auto cut and graze for PFT m_cut and m_grazed with constant sr_ugb
      ! 5: auto graze for PFT m_grazed with grazing litter during winter for LGM period
      f_autogestion = 0
      CALL getin_p('GRM_F_AUTOGESTION',f_autogestion)
      WRITE(numout,*)  'GRM_F_AUTOGESTION',f_autogestion
      ! whether animal is fed by extra feedstuffs
      f_complementation = 0
      CALL getin_p('GRM_F_COMPLEMENTATION',f_complementation)
      ! whether apply fertilizer
      f_fertilization = 1         
      CALL getin_p('GRM_F_FERTILIZATION',f_fertilization)
      ! JCCOMMENT 10April2015 already set and read at src_parameter
      !      ! number of management year cycled
      !      nb_year_management(:) = 0
      !      CALL getin_p('NB_YEAR_MANAGEMENT',nb_year_management)
      !      WRITE(numout,*) 'NB_YEAR_MANAGEMENT',nb_year_management
      ! f_postauto = 0-5
      ! 1: after f_autogestion=2 with varied sr_ugb and nb_ani
      ! 2: after f_postauto=1 with varied sr_ugb and nb_ani
      ! 3: simulation with constant sr_ugb and grazed_frac
      ! 4: simulation with increasing sr_ugb and constant grazed_frac
      ! 5: global simulation with prescribed sr_ugb from external file
      f_postauto = 0
      CALL getin_p('GRM_F_POSTAUTO',f_postauto)
      WRITE(numout,*)  'GRM_F_POSTAUTO',f_postauto
      ! the maximum impact to vcmax due to N fertilization
      ! N_effect =0.0 - 1.0
      N_effect=0.6
      CALL getin_p('GRM_N_EFFECT',N_effect)
      IF (N_effect .LT. 0.0 .OR. N_effect .GT. 1.0) THEN
        N_effect =0.6
      ENDIF

      ! 1.3 READ INITIAL CONDITION FILE FOR OLD/NEW ANIMAL MODULE

      file_param_init='/ccc/work/cont003/dsm/p529chan/input_gm/laq-int.init_cond.par'

      CALL getin_p('GRM_FILE_PARAM_INIT',file_param_init)
      WRITE (numout,*) 'GRM_FILE_PARAM_INIT',file_param_init
      OPEN(unit=61, file = file_param_init)

      READ(61, *, iostat = ier) toto(:)
      READ(61, *, iostat = ier) (wshtotcutinit_t(h), h=1,ncut)
      READ(61, *, iostat = ier) toto(:)
      READ(61, *, iostat = ier) toto(:)
      READ(61, *, iostat = ier) toto(:)

      READ(61, *, iostat = ier) toto(:)
      READ(61, *, iostat = ier) toto(:)
      READ(61, *, iostat = ier) toto(:)
      READ(61, *, iostat = ier) toto(:)
      READ(61, *, iostat = ier) toto(:)

      READ(61, *, iostat = ier) toto(:)
      READ(61, *, iostat = ier) toto(:)
      READ(61, *, iostat = ier) toto(:)
      READ(61, *, iostat = ier) toto(:)
      READ(61, *, iostat = ier) toto(:)

      READ(61, *, iostat = ier) toto(:)
      READ(61, *, iostat = ier) toto(:)
      READ(61, *, iostat = ier) toto(:)
      READ(61, *, iostat = ier) toto(:)
      READ(61, *, iostat = ier) toto(:)

      READ(61, *, iostat = ier) toto(:)
      READ(61, *, iostat = ier) toto(:)
      READ(61, *, iostat = ier) toto(:)
      READ(61, *, iostat = ier) toto(:)
      READ(61, *, iostat = ier) toto(:)

      READ(61, *, iostat = ier) toto(:)
      READ(61, *, iostat = ier) toto(:)
      READ(61, *, iostat = ier) toto(:)
      READ(61, *, iostat = ier) tcutmodel
      READ(61, *, iostat = ier) intakemax_t

      READ(61, *, iostat = ier) wanimal_t
      READ(61, *, iostat = ier) Type_animal
      READ(61, *, iostat = ier) toto(:)
      READ(61, *, iostat = ier) toto(:)
      READ(61, *, iostat = ier) toto(:)

      READ(61, *, iostat = ier) toto(:)
      READ(61, *, iostat = ier) toto(:)
      READ(61, *, iostat = ier) toto(:)
      READ(61, *, iostat = ier) toto(:)
      READ(61, *, iostat = ier) toto(:)

      READ(61, *, iostat = ier) toto(:)
      READ(61, *, iostat = ier) toto(:)
      READ(61, *, iostat = ier) toto(:)
      READ(61, *, iostat = ier) toto(:)
      READ(61, *, iostat = ier) toto(:)

      READ(61, *, iostat = ier) toto(:)
      READ(61, *, iostat = ier) toto(:)

      intakemax(:,:)=intakemax_t
      wanimal(:,:)=wanimal_t
      DO h=1,ncut
        wshtotcutinit(:,:,h)=wshtotcutinit_t(h)
      END DO
      CLOSE (61)
      WRITE(numout,*) 'wshtotcutinit_t',wshtotcutinit_t
      WRITE(numout,*) 'Type_animal',Type_animal

      ! 1.4 set constantes and initialise variables allocated above

      ! 1.4.1 initialisation the variables lied on animals cattle or sheep ?
      ! Type_Animal = 1,2,3 : Dairy cows, suckler cows, cows in old module 
      ! Type_Animal = 4,5 : Dairy heifers, suckler heifers
      ! Type_Animal = 6 : sheep in old module
      IF (Type_Animal==3)  THEN ! old module
      !090810 AIG changement du seuil de sortie des animaux Animalwgrazingmin trop faible
      ! changement de AnimalkintakeM pour garder que l'ingere garde la meme pente
      ! en fonction de la biomasse disponible et pour eviter un artefact de calcul
      ! Animalwgrazingmin = 0.03 ! Threshold shoot dry matter, under which animals are moved out for old module (kg.m-2)  
        Animalwgrazingmin(:,:)        = 0.03  ! N. Vuichard
        !AnimalkintakeM           = 0.18 ! AI Graux
        AnimalkintakeM(:,:)           = 0.1   ! N. Vuichard
        AnimalDiscremineQualite(:,:)  = 2    ! AI Graux  
      ELSEIF (Type_Animal .EQ. 6)THEN ! Sheep
        Animalwgrazingmin(:,:)        = 0.015
        AnimalkintakeM(:,:)           = 0.045
        AnimalDiscremineQualite(:,:)  = 3
      ELSE !new module
        Animalwgrazingmin(:,:)        = 0.03
        AnimalkintakeM(:,:)           = 0.1
        !Animalwgrazingmin        = 0.11 ! AI Graux ! unsued in the new module
        !AnimalkintakeM           = 0.18 ! AI Graux ! unsued in the new module
        AnimalDiscremineQualite(:,:)  = 2    ! AI Graux  
      ENDIF ! Type_Animal

      ! 1.4.2 concentrations : mean value
      c(:,:)                  = 0.0365122     !  4.22e-02
      n(:,:)                  = 0.00732556    !  8.17e-03
      napo(:,:)               = 0.000542054   !  6.39e-04
      nsym(:,:)               = 0.0108071     !  6.15e-03
      fn(:,:)                 = 0.0316223     !  4.15e-02   ! 2.64e-02
      ntot(:,:)               = 0.03471895    !  2.89e-02

      ! 1.4.3 initialisations of variables allocated above
      intake(:,:)                = 0.0
      intake_litter(:,:)         = 0.0
      intake_animal_litter(:,:)  = 0.0
      grazing_litter(:,:)        = 2
      litter_avail_totDM(:,:)    = 0.0
      devstage(:,:)              = GRM_devstage
      faecesc(:,:)               = 0.0
      faecesn(:,:)               = 0.0
      urinen(:,:)                = 0.0
      urinec(:,:)                = 0.0
      nel(:,:)                   = 0.0
      nanimaltot(:,:)            = 0.0
      !grazingcstruct(:,:)        = 0.0
      !grazingnstruct(:,:)        = 0.0
      tgrowth(:,:)               = 0.0
      wshtot(:,:) = (biomass(:,:,ileaf,icarbon) + biomass(:,:,isapabove,icarbon) + &
                    & biomass(:,:,ifruit,icarbon))/(1000*CtoDM) ! Unit: kgDM/m2
      wsh(:,:) = wshtot(:,:) / (1.0 + (mc /12.0)*c(:,:) + (mn /14.0)*n(:,:) )
      wshtotinit(:,:)            = wshtot(:,:)         
      wrtot(:,:) = (biomass(:,:,iroot,icarbon) + biomass(:,:,isapbelow,icarbon))/ &
                   & (1000*CtoDM)   ! Unit: kg/m2
      wr(:,:) = wrtot(:,:) / (1.0 + (mc /12.0)*c(:,:) + (mn /14.0)*n(:,:) )
      wnapo(:,:)                 = 0.0
      wnsym(:,:)                 = 0.0
      mux(:,:)                   = 0.0
      mugmean(:,:)               = 0.0
      sigx(:,:)                  = 0.0
      sigy(:,:)                  = 0.0
      gmeanslope(:,:)            = 0.0
      gzero(:,:)                 = 0.0
      gcor(:,:)                  = 0.0
      countschedule(:,:)         = 0
      cuttingend(:,:)            = 0
      regcount(:,:)              = 1
      gmean(:,:,:)               = 0.0
      tgmean(:,:,:)              = 0.0
      wc_frac(:,:)               = 0.0
      wgn(:,:)                   = 0.0
      tasum(:,:)                 = 0.0
      loss(:,:)                  = 0.0
      ! JCADD 20170421 for CNP merge
      losscnp(:,:,:,:) = 0.0
      grassharvestcnp(:,:,:) = 0.0
      grm_nfert(:,:,:) = 0.0
      grm_pfert(:,:) = 0.0
      intakecnp(:,:,:) = 0.0
      litter_intakecnp(:,:,:,:,:) = 0.0
      tramplingcnp(:,:,:,:) = 0.0
      faecescnp(:,:,:) = 0.0
      urinecnp(:,:,:) = 0.0
      milkcnp(:,:,:) = 0.0
      liveweightcnp(:,:,:) = 0.0
      ranimal_c(:,:) = 0.0
      EFch4_c(:,:) = 0.0
      grassorganicfertcnp(:,:,:) = 0.0
      pfert_mineral(:,:,:) = 0.0
      ! END JCADD
      lossc(:,:)                 = 0.0
      lossn(:,:)                 = 0.0
      tlossstart(:,:)            = 0.0
!      wshcutinit(:,:)            = 0.0
      deltatt                     = dt
      fertcount(:,:)             = 0.0
      c2nratiostruct(:,:)        = 150.0
      nfertammtot(:,:)           = 0.0
      nfertnittot(:,:)           = 0.0
      nfertammtotyear(:,:)       = 0.0
      nfertnittotyear(:,:)       = 0.0
      nfertammtotprevyear(:,:)   = 0.0
      nfertnittotprevyear(:,:)   = 0.0
      fcOrganicFertmetabolic(:,:)      = 0.0
      fcOrganicFertstruct(:,:)         = 0.0
      fnOrganicFertmetabolic(:,:)      = 0.0
      fnOrganicFertstruct(:,:)         = 0.0
      fnOrganicFerturine(:,:)          = 0.0
      fpOrganicFertmetabolic(:,:)      = 0.0
      fpOrganicFertstruct(:,:)         = 0.0
      fpOrganicFerturine(:,:)          = 0.0
      flag_fertilisation(:,:)    = 0
      fertil_year(:,:)           = .TRUE.
      tcut_verif(:,:,:)          = .FALSE. 
      tfert_verif(:,:,:)         = .FALSE. 
      tfert_verif2(:,:,:)        = .FALSE.
      tfert_verif3(:,:,:)        = .FALSE.
      nsatur_somerror_temp(:,:)          = 0.0
      nsatur_somerror(:,:)               = 0.0
      stoplimitant(:,:)                  = 0
      fertcount_start(:,:)               = 0
      fertcount_current(:,:)             = 0
      nnonlimit_SOMerror(:,:)            = 0.0
      nnonlimit_SOMerrormax(:,:)         = 0.5
      controle_azote_sum_mem(:,:)        = 0.0
      n_auto(:,:)                        = 4
      flag_fertilisation(:,:)            = 0
      YIELD_RETURN(:,:) = 0.0
      sr_ugb_init(:) = 0.0
      compt_cut(:,:) =0.0
      frequency_cut(:,:) =0.0
      sr_wild(:,:) = 0.0
      flag_cutting(:,:) = 0
      tamean1(:)         = 273.0
      tamean2(:)         = 273.0
      tamean3(:)         = 273.0
      tamean4(:)         = 273.0
      tamean5(:)         = 273.0
      tamean6(:)         = 273.0
      tameand(:)         = 273.0
      tameanw(:)         = 0.0
      tacumm(:)          = 0.0
      tacummprev(:)      = 0.0
      tsoilcumm(:)       = 0.0
      tsoilcummprev(:)   = 0.0
      tsoilmeand(:)      = 273.0
      tcut0(:,:)         = 0.0
      N_fert_total(:,:) = 0.0
      Fert_on(:,:) = 0.0
      Fert_sn(:,:) = 0.0
      Fert_PRP(:,:) = 0.0
      ! 1.5 Define PFT that used for optimization, cutting, and grazing
      DO j=2,nvm
        IF (is_grassland_cut(j) .AND. (.NOT. is_grassland_grazed(j)) .AND. &
           (.NOT. is_c4(j)) .AND. (.NOT.is_tree(j))) THEN
          mcut_C3=j
        END IF
        IF (is_grassland_manag(j) .AND. (.NOT. is_grassland_cut(j)) .AND. &
           (.NOT.is_grassland_grazed(j)) .AND. (.NOT. is_c4(j)) .AND. &
           (.NOT. is_tree(j))) THEN
          mauto_C3=j
        END IF
        IF (is_grassland_manag(j) .AND. (is_grassland_grazed(j)) .AND. &
           (.NOT. is_grassland_cut(j)) .AND. (.NOT. is_c4(j)) .AND. &
           (.NOT. is_tree(j))) THEN
          mgraze_C3=j
        END IF
        IF (is_grassland_cut(j) .AND. (.NOT. is_grassland_grazed(j)) .AND. &
           (is_c4(j)) .AND. (.NOT. is_tree(j))) THEN
          mcut_C4=j
        END IF
        IF ( is_grassland_manag(j) .AND. (.NOT. is_grassland_cut(j)) .AND. &
           (.NOT. is_grassland_grazed(j)) .AND. (is_c4(j)) .AND. &
           (.NOT. is_tree(j))) THEN
          mauto_C4=j
        END IF
        IF ( is_grassland_manag(j) .AND. (is_grassland_grazed(j)) .AND. &
           (.NOT. is_grassland_cut(j)) .AND. (is_c4(j)) .AND. &
           (.NOT. is_tree(j))) THEN
          mgraze_C4=j
        END IF
        IF ((.NOT. is_grassland_manag(j)) .AND. (.NOT. is_grassland_grazed(j)) .AND. &
           (.NOT. is_grassland_cut(j)) .AND. (.NOT. is_c4(j)) .AND. &
           (.NOT. is_tree(j)) .AND. natural(j)) THEN
          mnatural_C3=j
        END IF
        IF ((.NOT. is_grassland_manag(j)) .AND. (.NOT. is_grassland_grazed(j)).AND. &
          (.NOT. is_grassland_cut(j)) .AND. (is_c4(j)) .AND. &
          (.NOT. is_tree(j)) .AND. natural(j)) THEN
          mnatural_C4=j
        END IF
      END DO ! nvm
      !WRITE(numout,*) 'PFT_M',mauto_C3,mcut_C3,mgraze_C3,mauto_C4,mcut_C4,mgraze_C4,mnatural_C3,mnatural_C4 
      ! avoid PFT = 0
      IF (mauto_C4 .EQ. 0) THEN 
        mauto_C4=1
      ENDIF
      IF (mcut_C4 .EQ. 0) THEN
        mcut_C4=1
      ENDIF
      IF (mgraze_C4 .EQ. 0) THEN
        mgraze_C4=1
      ENDIF
      IF (mauto_C3 .EQ. 0) THEN
        mauto_C3=1
      ENDIF
      IF (mcut_C3 .EQ. 0) THEN
        mcut_C3=1
      ENDIF
      IF (mgraze_C3 .EQ. 0) THEN
        mgraze_C3=1
      ENDIF
      IF (mnatural_C4 .EQ. 0) THEN
        mnatural_C4=1
      ENDIF
      IF (mnatural_C3 .EQ. 0) THEN
        mnatural_C3=1
      ENDIF
      WRITE(numout,*) 'PFT_M2',mauto_C3,mcut_C3,mgraze_C3,mauto_C4,mcut_C4,mgraze_C4,mnatural_C3,mnatural_C4
  
      ! 1.6 Initialization of management related parameters
      ! for each management option
      IF (f_postauto .EQ. 0) THEN
        IF (f_autogestion .EQ. 1) THEN
      ! 1: auto cut for PFT m_auto
          ! keep wshtotsum for mauto_C3 and C4  
          sr_ugb         = 1e-5
          compt_ugb      = 0.0
          nb_ani         = 5e-6
          grazed_frac         = 0.50
        ELSE IF (f_autogestion .EQ. 2) THEN
      ! 2: auto graze for PFT m_auto
          ! keep wshtotsum for each year calculation of import_yield
          cut_year = 10 
          CALL getin_p('GRM_NB_CUT_YEAR',cut_year)
          WHERE ( wshtotsum(:,mauto_C3) .GE. 0.0)
            import_yield(:,mauto_C3) = wshtotsum(:,mauto_C3) / cut_year
          ENDWHERE
          WHERE ( wshtotsum(:,mauto_C4) .GE. 0.0)
            import_yield(:,mauto_C4) = wshtotsum(:,mauto_C4) / cut_year
          ENDWHERE
          ! keep sr_ugb, compt_ugb, nb_ani, grazed_frac 
          ! infact it could be keep just the value read from restart
          compt_ugb      = 0.0
        ELSE IF ((f_autogestion .GE. 3) .AND. &
                (f_autogestion .LE. 5)) THEN 
      ! 3: auto cut and graze for PFT m_cut and m_grazed with increasing sr_ugb
      ! 4: auto cut and graze for PFT m_cut and m_grazed with constant sr_ugb
      ! 5: auto graze for PFT m_grazed with grazing litter during winter for LGM
          ! keep the grazing variables from restart
          compt_ugb      = 0.0
          wshtotsum (:,:) = 0.0
          import_yield (:,:) = 0.0
        ELSE ! f_postauto = 0 and f_autogestion > 5
          sr_ugb         = 1e-5
          compt_ugb      = 0.0
          nb_ani         = 5e-6
          grazed_frac         = 0.50
          wshtotsum (:,:) = 0.0
          import_yield (:,:) = 0.0
        ENDIF ! f_autogestion 
      ELSE IF (f_postauto .EQ. 1) THEN
      ! 1: after f_autogestion=2 with varied sr_ugb and nb_ani
      ! ONLY run for one years
        tmp_sr_ugb_C3(:)=sr_ugb(:,mauto_C3)
        tmp_sr_ugb_C4(:)=sr_ugb(:,mauto_C4)
        sr_ugb         = 1e-5
        sr_ugb(:,mgraze_C3)      = tmp_sr_ugb_C3(:)
        sr_ugb(:,mgraze_C4)      = tmp_sr_ugb_C4(:)
        tmp_nb_ani_C3(:)=nb_ani(:,mauto_C3)
        tmp_nb_ani_C4(:)=nb_ani(:,mauto_C4)
        nb_ani         = 5e-6
        nb_ani(:,mgraze_C3)         = tmp_nb_ani_C3(:)
        nb_ani(:,mgraze_C4)         = tmp_nb_ani_C4(:)
        tmp_grazed_frac_C3(:)=grazed_frac(:,mauto_C3)
        tmp_grazed_frac_C4(:)=grazed_frac(:,mauto_C4)
        grazed_frac         = 0.50
        grazed_frac(:,mgraze_C3)         = tmp_grazed_frac_C3(:)
        grazed_frac(:,mgraze_C4)         = tmp_grazed_frac_C4(:)
        WHERE (sr_ugb(:,mgraze_C3) .GT. 0.0)
          grazed_frac(:,mgraze_C3)  = nb_ani(:,mgraze_C3)/sr_ugb(:,mgraze_C3)
        ELSEWHERE
          grazed_frac(:,mgraze_C3)  = tmp_grazed_frac_C3(:)
        ENDWHERE
        WHERE (sr_ugb(:,mgraze_C4) .GT. 0.0)
          grazed_frac(:,mgraze_C4)  = nb_ani(:,mgraze_C4)/sr_ugb(:,mgraze_C4)
        ELSEWHERE
          grazed_frac(:,mgraze_C4)  = tmp_grazed_frac_C4(:)
        ENDWHERE
        compt_ugb      = 0.0
        wshtotsum (:,:) = 0.0
        tmp_import_yield_C3(:) = import_yield(:,mauto_C3)
        tmp_import_yield_C4(:) = import_yield(:,mauto_C4)
        import_yield = 0.0
        import_yield (:,mgraze_C3) = tmp_import_yield_C3(:)
        import_yield (:,mgraze_C4) = tmp_import_yield_C4(:)
      ELSE IF (f_postauto .GE. 2) THEN
      ! 2: after f_postauto=1 with varied sr_ugb and nb_ani
      ! 3: simulation with constant sr_ugb and grazed_frac
      ! 4: simulation with increasing sr_ugb and constant grazed_frac
      ! 5: global simulation with prescribed sr_ugb from external file
          ! keep the grazing variables from restart
        compt_ugb      = 0.0
        wshtotsum (:,:) = 0.0
        ! keep import_yield from restart for mean value saving
        ! import_yield = 0.0
      ELSE
        sr_ugb         = 1e-5
        compt_ugb      = 0.0
        nb_ani         = 5e-6
        grazed_frac         = 0.50
        compt_ugb      = 0.0
        wshtotsum (:,:) = 0.0
        import_yield = 0.0
      ENDIF ! f_autogestion or f_postauto
 
      wshtotsumprevyear(:,:) = 0.0 
      DM_cutyearly(:,:)=0.0
      C_cutyearly(:,:) =0.0
      wshtotsumprev   (:,:) = 0.0
      controle_azote(:,:,:)       = 0.0
      controle_azote_sum(:,:)        = 0.0
      trampling(:,:)              = 0.0
      count_year            = 1
      year_count1 = 0
      year_count2 = 0
      tcut(:,:,:) = 500.0
      tfert(:,:,:) = 500.0
      nfertamm(:,:,:) = 0.0
      nfertnit(:,:,:) = 0.0
      nanimal(:,:,:) = 0.0
      tanimal(:,:,:) = 500.0
      danimal(:,:,:) = 0.0
      nliquidmanure(:,:,:) = 0.0
      nslurry(:,:,:) = 0.0
      nsolidmanure(:,:,:) = 0.0
      legume_fraction(:,:) =0.0
      soil_fertility(:,:) = 1.0
      ndeposition(:,:) = 0.0
      PIYcow(:,:,:) = 0.0
      PIMcow(:,:,:) = 0.0
      BCSYcow(:,:,:) = 0.0
      BCSMcow(:,:,:) = 0.0
      PICcow(:,:,:) = 0.0
      AGE_cow_P(:,:,:) = 36.0
      AGE_cow_M(:,:,:) = 54.0
      Forage_quantity(:,:,:) = 0
  
      IF (blabla_pasim) PRINT *, 'PASIM : end memory allocation'
 
      ! 1.7 read management maps/files 
      ! get_map of 1 spatial .nc file or 0 old txt/dat file  
      f_management_map = 0
      CALL getin_p ('GRM_F_MANAGEMENT_MAP',f_management_map)
      WRITE(numout,*)  'GRM_F_MANAGEMENT_MAP',f_management_map
      f_deposition_map = 0
      CALL getin_p ('GRM_F_DEPOSITION_MAP',f_deposition_map)
      WRITE(numout,*)  'GRM_F_DEPOSITION_MAP',f_deposition_map  
      f_grazing_map = 0
      CALL getin_p ('GRM_F_GRAZING_MAP',f_grazing_map)
      WRITE(numout,*)  'GRM_F_GRAZING_MAP',f_grazing_map
   
      IF (f_management_map .EQ. 1) THEN
        management_map='GRM_input.nc'
        !'/ccc/work/cont003/dsm/p529chan/data/eur_management_interpolated.nc'
        CALL getin_p('GRM_MANAGEMENT_MAP',management_map)
        WRITE(numout,*) 'GRM_MANAGEMENT_MAP',management_map
        fertility_map='/ccc/work/cont003/dsm/p529chan/data/eur_fertility.nc'
        CALL getin_p('GRM_FERTILITY_MAP',fertility_map)
        WRITE(numout,*) 'GRM_FERTILITY_MAP',fertility_map
  
        deposition_map='GRM_input.nc'
        !'/ccc/work/cont003/dsm/p529chan/data/eur_Ndeposition_NCAR.nc'
        CALL getin_p('GRM_DEPOSITION_MAP',deposition_map)
        WRITE(numout,*) 'GRM_DEPOSITION_MAP',deposition_map  
        grazing_map='GRM_input.nc'
        !'/ccc/scratch/cont003/dsm/p529chan/glbdata/glb_sr_ugb_1961_2010_adjusted.nc'
        CALL getin_p('GRM_GRAZING_MAP',grazing_map)
        WRITE(numout,*) 'GRM_GRAZING_MAP',grazing_map
  
        ! read management map    
        CALL reading_map_manag(&
               npts,lalo, neighbours, resolution, contfrac, & 
               count_year, nb_year_management,& 
               management_intensity,&
               management_start,&
               tcut, tfert, nfertamm, nfertnit, pfert_mineral, &
               nanimal, tanimal, danimal,&
               nliquidmanure, nslurry, nsolidmanure,&
               legume_fraction,soil_fertility,&
               deposition_start,ndeposition,sr_ugb,sr_wild)
!        ! calculate effect of N fertilizer to vcmax
!        CALL calc_N_limfert(&
!               npts,nfertamm, nfertnit,&
!               nliquidmanure, nslurry, nsolidmanure,&
!               legume_fraction,soil_fertility,ndeposition,&
!               N_fert_total,N_limfert) 
      ELSE
        ! re-initial management variables
        tcut(:,:,:) = 500.0
        tfert(:,:,:) = 500.0
        nfertamm(:,:,:) = 0.0
        nfertnit(:,:,:) = 0.0
        pfert_mineral(:,:,:) = 0.0
        nanimal(:,:,:) = 0.0
        tanimal(:,:,:) = 500.0
        danimal(:,:,:) = 0.0
        nliquidmanure(:,:,:) = 0.0
        nslurry(:,:,:) = 0.0
        nsolidmanure(:,:,:) = 0.0
        ndeposition(:,:) = 0.0
  !!! delete FIRE_MANAGEMENT READ: not used in LGM
        CALL getin_p('GRM_FILE_MANAGEMENT',file_management)
        WRITE(numout,*)  'GRM_FILE_MANAGEMENT',file_management
        IF (blabla_pasim) PRINT *, 'PASIM : reading management conditions'
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!! READ NEW MANAGEMENT TXT DAT FILE gmjc
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        CALL reading_new_animal(&
             npts           , &
             nb_year_management , &
             tcutmodel      , &
             tcut           , &
             tfert          , &
             nfertamm       , &
             nfertnit       , &
             nanimal        , &
             tanimal        , &
             danimal        , &
             nliquidmanure  , &
             nslurry        , &
             nsolidmanure   , &
             PIYcow         , &
             PIMcow         , &
             BCSYcow        , &
             BCSMcow        , &
             PICcow         , &
             AGE_cow_P      , &
             AGE_cow_M      , &
             Forage_quantity)
  !      CALL getin('SOIL_FERTILITY',fertility_legume_t)
  !      soil_fertility(:,:)=fertility_legume_t
  !      CALL getin('LEGUME_FRACTION',fertility_legume_t)
  !      legume_fraction(:,:)=fertility_legume_t
  !
  !      CALL calc_N_limfert(&
  !             npts,nfertamm, nfertnit,&
  !             nliquidmanure, nslurry, nsolidmanure,&
  !             legume_fraction,soil_fertility,ndeposition,&
  !             N_fert_total,N_limfert)
  
      ENDIF ! f_management_map

      DO k=1,nstocking
        WHERE (tfert(:,:,k) .NE. 500) 
          apport_azote(:,:) = apport_azote(:,:) + nfertamm(:,:,k) + nfertnit(:,:,k)    
        END WHERE  
      END DO
        !************************************************
        !************************************************
        ! modifs Nico 20/07/2004
        !************************************************
        !************************************************
        ! MODIF INN
      IF (f_nonlimitant .EQ. 1) THEN
        IF (f_autogestion .NE. 2) THEN
          WHERE (tcut(:,:,1) .EQ. 500.0)
            stoplimitant(:,:) = 1
          END WHERE
        ENDIF
        DO j=2,nvm
          DO i=1,npts
            IF (tfert(i,j,1) .EQ. 500.0) THEN
              stoplimitant(i,j) = 1
            ELSE
              compt_fert = 1
              min_fert   = 1
              DO WHILE (tfert(i,j,compt_fert) .NE. 500.0)
!                 print *, compt_fert, min_fert
!                 print *, controle_azote(i,j,compt_fert)
!                 print *, controle_azote(i,j,min_fert)
                IF (controle_azote(i,j,compt_fert) .GT. controle_azote(i,j,min_fert)) THEN
                  min_fert = compt_fert
                ENDIF
                  compt_fert = compt_fert + 1
              END DO
              fert_max = compt_fert - 1
              IF ((min_fert - 1) .EQ. 0) THEN
                fertcount_start(i,j) = fert_max
              ELSE
                fertcount_start(i,j) = min_fert - 1
              ENDIF
                i_compt = min_fert + 1
              DO WHILE ( tfert(i,j,i_compt) .NE. 500.0 )
                controle_azote(i,j,i_compt) = controle_azote(i,j,i_compt - 1)+&
                  controle_azote(i,j,i_compt)
                i_compt = i_compt + 1
              END DO
              IF ( min_fert .NE. 1. ) THEN
                controle_azote(i,j,1) = controle_azote(i,j,1) + controle_azote(i,j,fert_max)
                i_compt = 2
                DO WHILE (i_compt .NE. min_fert)
                  controle_azote(i,j,i_compt) = controle_azote(i,j,i_compt-1)+&
                    controle_azote(i,j,i_compt)
                  i_compt = i_compt + 1
                END DO
              ENDIF
            ENDIF
          END DO ! i
        END DO !j
          fertcount_current(:,:) = fertcount_start(:,:)
      ENDIF
        ! fin initialisation auto gestion nicolas
    END IF init_grassland

    ! 2 updating variables each day (new_day) 
    ! update the root/shoot dry matter variables
    wshtot(:,:) = (biomass(:,:,ileaf,icarbon) + biomass(:,:,isapabove,icarbon) + &
                 & biomass(:,:,ifruit,icarbon))/(1000*CtoDM) ! Unit: kgDM/m2
    wsh(:,:) = wshtot(:,:) / (1.0 + (mc /12.0)*c(:,:) + (mn /14.0)*n(:,:) )
    wrtot(:,:) = (biomass(:,:,iroot,icarbon) + biomass(:,:,isapbelow,icarbon))/ &
                 & (1000*CtoDM)   ! Unit: kg/m2
    wr(:,:) = wrtot(:,:) / (1.0 + (mc /12.0)*c(:,:) + (mn /14.0)*n(:,:) )

    WRITE (numout,*)  'tjulian', tjulian
    n_day : IF (new_day) THEN

      ! GMEAN
      ! Taux de croissance moyen de la repousse
      h  = 1

      DO WHILE (h  .LT. ngmean)
        gmean(:,:,h ) = gmean(:,:,h +1)
        h  = h  + 1
      END DO

      DO j=2,nvm  
        DO i=1,npts
          IF ((tgrowth(i,j) .GT. 0.0) .AND. (devstage(i,j) .GE. 2.0)) THEN
            gmean(i,j,ngmean) = MAX (0.0, (wshtot(i,j) - wshtotcutinit(i,j,regcount(i,j)))/tgrowth(i,j))
          ELSEIF ((tgrowth(i,j) .GT. 0.0) .AND. (devstage(i,j) .GT. 0.0) .AND. &
            & (regcount(i,j) .GT. 1)) THEN
            gmean(i,j,ngmean) = MAX (0.0, (wshtot(i,j) - wshtotcutinit(i,j,regcount(i,j)))/tgrowth(i,j))
          ELSEIF ((tgrowth(i,j) .GT. 0.0) .AND. (devstage(i,j) .GT. 0.0) .AND. &
            & (regcount(i,j) .EQ. 1)) THEN
            gmean(i,j,ngmean) = MAX (0.0, (wshtot(i,j)  - wshtotinit(i,j))/tgrowth(i,j))
          ELSE
            gmean(i,j,ngmean) = 0.0
          ENDIF
        END DO
      ENDDO
 
      h = 1
      DO WHILE (h .LE. ngmean) 
        tgmean(:,:,h) = h
        h = h + 1
      END DO
      ! for each new day the daily value should be reset to 0
      ! otherwise, it will be repeatly cumulated
      Fert_on(:,:) = 0.0
      Fert_sn(:,:) = 0.0
      Fert_PRP(:,:) = 0.0      

      ! 1.4.2 concentrations : mean value
      c(:,:)                  = 0.0365122     !  4.22e-02
      n(:,:)                  = 0.00732556    !  8.17e-03
      napo(:,:)               = 0.000542054   !  6.39e-04
      nsym(:,:)               = 0.0108071     !  6.15e-03
      fn(:,:)                 = 0.0316223     !  4.15e-02   ! 2.64e-02
      ntot(:,:)               = 0.03471895    !  2.89e-02
      intake(:,:)                = 0.0
      intake_litter(:,:)         = 0.0
      intake_animal_litter(:,:)  = 0.0
      grazing_litter(:,:)        = 2
      litter_avail_totDM(:,:)    = 0.0
      faecesc(:,:)               = 0.0
      faecesn(:,:)               = 0.0
      urinen(:,:)                = 0.0
      urinec(:,:)                = 0.0
      nel(:,:)                   = 0.0
      nanimaltot(:,:)            = 0.0
      loss = 0.0 
      ! JCADD 20170421 for CNP merge
      ! set daily variable to 0
      losscnp(:,:,:,:) = 0.0
      grassharvestcnp(:,:,:) = 0.0
      grm_nfert(:,:,:) = 0.0
      grm_pfert(:,:) = 0.0
      intakecnp(:,:,:) = 0.0
      litter_intakecnp(:,:,:,:,:) = 0.0
      tramplingcnp(:,:,:,:) = 0.0
      faecescnp(:,:,:) = 0.0
      urinecnp(:,:,:) = 0.0
      milkcnp(:,:,:) = 0.0
      liveweightcnp(:,:,:) = 0.0
      ranimal_c(:,:) = 0.0
      EFch4_c(:,:) = 0.0
      grassorganicfertcnp(:,:,:) = 0.0
      ! END JCADD

!JC MOD024 
      devstage(:,:) = GRM_devstage(:,:)
    END IF n_day

    ! 3 updating variables at the end of the year (last day = new_year =
    n_year : IF (new_year) THEN
      WRITE(numout,*) 'LastTsYear gm'
      tcut_verif(:,:,:)         = .FALSE. 
      fertil_year(:,:)          = .TRUE. 
      tasum(:,:)                = 0.0
      regcount(:,:)             = 1
      nfertammtotprevyear(:,:)  = nfertammtot 
      nfertnittotprevyear(:,:)  = nfertnittot 
      fertcount(:,:)            = 0
      nfertammtotyear(:,:)      = 0.0
      nfertnittotyear(:,:)      = 0.0
      fnatmsum(:,:)             = 0.0
      tfert_verif(:,:,:)        = .FALSE.
      tfert_verif2(:,:,:)       = .FALSE.
      tfert_verif3(:,:,:)       = .FALSE.
      fcOrganicFertmetabolicsum(:,:) = 0.0
      fcOrganicFertstructsum(:,:)    = 0.0
      fnOrganicFertmetabolicsum(:,:) = 0.0
      fnOrganicFertstructsum(:,:)    = 0.0
      fnOrganicFerturinesum(:,:)     = 0.0
!      devstage(:,:)             = 0.0
      fertcount(:,:)            = 0
      tgrowth (:,:)             = 0.0
      tfert_modif(:,:,:)        = 500.0
      frequency_cut(:,:) = compt_cut(:,:)
      compt_cut(:,:) = 0.0

      IF (f_saturant .EQ. 1) THEN
         nfertamm(:,:,:)  = 0.025
         nfertnit(:,:,:)  = 0.025
         pfert_mineral(:,:,:)  = 0.005
         nsatur_somerror(:,:)      = 0.0
         nsatur_somerror_temp(:,:) = 0.0
      END IF
      ! calculate annual grass forage production
      DM_cutyearly(:,:)= wshtotsum(:,:)-wshtotsumprevyear(:,:)
      C_cutyearly(:,:) = DM_cutyearly(:,:) * 1000 * CtoDM
      ! should be after calculating the import_yield
      !wshtotsumprevyear(:,:) = wshtotsum(:,:)

      ! calculate import_yield saved to restart for output
      ! and for updating grazing variables in grazing subroutine
      IF ((f_postauto .GE. 1) .OR. (f_autogestion  .EQ. 4) .OR. &
         !!!! JCMODIF 290714 for postaut = 5
         (f_autogestion  .EQ. 3) ) THEN
        import_yield(:,mgraze_C3) = wshtotsum(:,mcut_C3)-wshtotsumprevyear(:,mcut_C3)
        wshtotsumprevyear(:,mcut_C3) = wshtotsum(:,mcut_C3)
        import_yield(:,mgraze_C4) = wshtotsum(:,mcut_C4)-wshtotsumprevyear(:,mcut_C4)
        wshtotsumprevyear(:,mcut_C4) = wshtotsum(:,mcut_C4)
        ! if as trunk that restart and initiailize every year
        ! save wshtotsumprevyear is useless
        ! but if as NV driver run for many years
        ! save wshtotsumprevyear is necessary because wshtotsum will keep
        ! inscresing
      END IF

      wshtotsumprevyear(:,:) = wshtotsum(:,:)
      wshtotsumprev(:,:)          = 0.0
      c(:,:)                  = 0.0365122     !  4.22e-02
      n(:,:)                  = 0.00732556    !  8.17e-03
      napo(:,:)               = 0.000542054   !  6.39e-04
      nsym(:,:)               = 0.0108071     !  6.15e-03
      fn(:,:)                 = 0.0316223     !  4.15e-02   ! 2.64e-02
      ntot(:,:)               = 0.03471895    !  2.89e-02  

      ! count_year is useless for trunk driver
      ! only necessary for NV driver run for many years
      ! unless save it to restart in the future
      count_year = count_year + 1
      IF (count_year .LT. 30) THEN
        year_count1 = count_year-1
        year_count2 = 0
      ELSEIF (count_year .GE. 30) THEN
        year_count1 = 29
        year_count2 = count_year - 29
      ELSE 
        year_count1 = 29
        year_count2 = 21
      ENDIF

!JCCOMMENT There is no need to read again in standard trunk driver
      ! read management map at the end of year
      ! is only useful for NV driver and multi-year maps
      ! for trunk driver, the annual file will be changed every year
      ! and read when initialize

      ! get_map of spatial .nc file or old txt/dat file  
      IF (f_management_map .EQ. 1) THEN
        ! re-initial management variables
        tcut(:,:,:) = 500.0
        tfert(:,:,:) = 500.0
        nfertamm(:,:,:) = 0.0
        nfertnit(:,:,:) = 0.0
        pfert_mineral(:,:,:) = 0.0
        nanimal(:,:,:) = 0.0
        tanimal(:,:,:) = 500.0
        danimal(:,:,:) = 0.0
        nliquidmanure(:,:,:) = 0.0
        nslurry(:,:,:) = 0.0
        nsolidmanure(:,:,:) = 0.0
        ndeposition(:,:) = 0.0
        CALL reading_map_manag(&  
               npts, lalo, neighbours, resolution, contfrac, &
               count_year, nb_year_management,&
               management_intensity,&
               management_start,& 
               tcut, tfert, nfertamm, nfertnit,pfert_mineral,&
               nanimal, tanimal, danimal,&
               nliquidmanure, nslurry, nsolidmanure,&
               legume_fraction,soil_fertility,&
               deposition_start,ndeposition,sr_ugb,sr_wild)
!!!  
!!!        CALL calc_N_limfert(&
!!!               npts,nfertamm, nfertnit,&
!!!               nliquidmanure, nslurry, nsolidmanure,&
!!!               legume_fraction,soil_fertility,ndeposition,&
!!!               N_fert_total,N_limfert)
!!!  
      ELSE
        IF (ANY(nb_year_management(:) .GT. 1)) THEN
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!! READ NEW MANAGEMENT FILE gmjc
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! delete FIRE_MANAGEMENT READ: not used in LGM
            CALL reading_new_animal(&
               npts           , &
               nb_year_management , &
               tcutmodel      , &
               tcut           , &
               tfert          , &
               nfertamm       , &
               nfertnit       , &
               nanimal        , &
               tanimal        , &
               danimal        , &
               nliquidmanure  , &
               nslurry        , &
               nsolidmanure   , &
               PIYcow         , &
               PIMcow         , &
               BCSYcow        , &
               BCSMcow        , &
               PICcow         , &
               AGE_cow_P      , &
               AGE_cow_M      , &
               Forage_quantity)
!!!          ! re-initial management variables
!!!          tcut(:,:,:) = 500.0
!!!          tfert(:,:,:) = 500.0
!!!          nfertamm(:,:,:) = 0.0
!!!          nfertnit(:,:,:) = 0.0
!!!          nanimal(:,:,:) = 0.0
!!!          tanimal(:,:,:) = 500.0
!!!          danimal(:,:,:) = 0.0
!!!          nliquidmanure(:,:,:) = 0.0
!!!          nslurry(:,:,:) = 0.0
!!!          nsolidmanure(:,:,:) = 0.0
!!!          ndeposition(:,:) = 0.0
!!!          CALL calc_N_limfert(&
!!!                 npts,nfertamm, nfertnit,&
!!!                 nliquidmanure, nslurry, nsolidmanure,&
!!!                 legume_fraction,soil_fertility,ndeposition,&
!!!                 N_fert_total,N_limfert)
!!!    
        END IF ! nb_year_management
!!!
!!!        DO k=1,nstocking
!!!          WHERE (tfert(:,:,k) .NE. 500) 
!!!            apport_azote(:,:) = apport_azote(:,:)  + nfertamm(:,:,k) + nfertnit(:,:,k)
!!!          END WHERE
!!!        END DO
!!!
      END IF ! f_management_map

    END IF n_year

    ! 4 fertilization
    ! Fertilisation from PaSim 2011
    ! 4.1 ****** RUN USERS OR RUN SATURANT *****
    users_or_saturant_fert : IF ((tcutmodel .EQ. 0) .AND. (f_saturant .EQ. 0)) THEN

      ! flag_fertilisation : flag for spatialization of cutting

      DO k=1,nstocking
        flag_fertilisation(:,:) = 0
        DO j=2 ,nvm
          DO i=1,npts
            IF (tfert_verif(i,j,k) .EQV. .FALSE. ) THEN
              IF ((tfert(i,j,k)-0.5 .LE. REAL(tjulian)) .AND. &
                    (tfert(i,j,k)+0.9 .GE. REAL(tjulian)) .AND. &
                    (tfert_verif(i,j,k) .EQV. .FALSE.)) THEN
                tfert_verif(i,j,k) = .TRUE.
                flag_fertilisation(i,j) = 1 
                Fert_sn(i,j) = Fert_sn(i,j) + nfertamm(i,j,k) + nfertnit(i,j,k)
                !counter for fertilizer application
                fertcount(i,j)  = fertcount(i,j)  + 1
                !mineral fertilization
                nfertammtot(i,j)  = nfertammtot(i,j)  + nfertamm(i,j,fertcount(i,j))
                nfertnittot(i,j)  = nfertnittot(i,j)  + nfertnit(i,j,fertcount(i,j))
      
                nfertammtotyear(i,j)  = nfertammtot(i,j)  - nfertammtotprevyear(i,j)
                nfertnittotyear(i,j)  = nfertnittot(i,j)  - nfertnittotprevyear(i,j)
                ! JCADD 20170421 for CNP merge
                ! nfertamm input is kg m-2
                grm_nfert(i,j,iammonium) = nfertamm(i,j,fertcount(i,j))*1e3
                grm_nfert(i,j,initrate) = nfertnit(i,j,fertcount(i,j))*1e3
                grm_pfert(i,j) = pfert_mineral(i,j,fertcount(i,j))*1e3
                ! END JCADD
                IF (f_nonlimitant .EQ. 1.) THEN
                  controle_azote_sum_mem(i,j) = controle_azote_sum(i,j)
                  IF (.NOT. (tfert_verif3(i,j,fertcount_start(i,j))) ) THEN
                    tfert_verif3(i,j,fertcount_start(i,j))= .TRUE.
                    controle_azote_sum(i,j) = 0.
                  END IF
                END IF
              END IF
            END IF
          END DO ! i
        END DO ! j
      END DO ! k
      !*****************************************
      ! MODIFS NICO AUTO MANAGEMENT DE PASIM
      !*****************************************

    ELSE IF (f_saturant .EQ. 1) THEN   !***** RUN SATURANT *******

      flag_fertilisation(:,:) = 0
      DO j=2 ,nvm
        DO i=1,npts
          IF (( tjulian .GE. tfert(i,j,fertcount(i,j) + 1)) .AND. &
             ( tjulian .LT. (tfert(i,j,fertcount(i,j) + 2) - 1))) THEN
             !JCmodif 110523  with problem
             !above means tjulian between two tfert 
             !undaily(i) uptake n daily always=0
             !thetas volumetric water content in soil layer h
             !thetasfc water field capacity
             !!!!! For we did not consider undaily , there will be no point need to fert???                  
             ! IF ((undaily(i) .GT. 0.0) .AND. (thetas(i,1) .LE. thetasfc(i,1))) THEN
             flag_fertilisation(i,j) = 1
              !counter for fertilizer application
              fertcount(i,j)  = fertcount(i,j)  + 1
              !mineral fertilization
              nfertammtot(i,j)  = nfertammtot(i,j)  + nfertamm(i,j,fertcount(i,j))
              nfertnittot(i,j)  = nfertnittot(i,j)  + nfertnit(i,j,fertcount(i,j))

              nfertammtotyear(i,j)  = nfertammtot(i,j)  - nfertammtotprevyear(i,j)
              nfertnittotyear(i,j)  = nfertnittot(i,j)  - nfertnittotprevyear(i,j)
                ! JCADD 20170421 for CNP merge
                ! nfertamm input is kg m-2
                grm_nfert(i,j,iammonium) = nfertamm(i,j,fertcount(i,j))*1e3
                grm_nfert(i,j,initrate) = nfertnit(i,j,fertcount(i,j))*1e3
                grm_pfert(i,j) = pfert_mineral(i,j,fertcount(i,j))*1e3
                ! END JCADD
              IF (f_nonlimitant .EQ. 1.) THEN
                controle_azote_sum_mem(i,j) = controle_azote_sum(i,j)
                IF (.NOT. (tfert_verif3(i,j,fertcount_start(i,j))) ) THEN
                  tfert_verif3(i,j,fertcount_start(i,j))= .TRUE.
                  controle_azote_sum(i,j) = 0.
                END IF
              END IF

              IF (controle_azote_sum(i,j) .GT. 0.) THEN
                 nsatur_somerror_temp(i,j) = &
                   ABS(controle_azote(i,j,fertcount(i,j)) - controle_azote_sum(i,j)) / &
                   controle_azote_sum(i,j)
              ENDIF
              IF (nsatur_somerror_temp(i,j) .GT. nsatur_somerror(i,j)) THEN
                nsatur_somerror(i,j) = nsatur_somerror_temp(i,j)
              ENDIF
              controle_azote(i,j,fertcount(i,j) ) = controle_azote_sum(i,j)
              tfert_modif(i,j,fertcount(i,j) )    = tjulian
          ELSE
              flag_fertilisation(i,j) = 2

            IF ((tjulian .GE. (tfert(i,j,fertcount(i,j)+2)-1))  .AND. &
              (tjulian .LE. (tfert(i,j,fertcount(i,j)+2)-0.1)) .AND. &
              (.NOT.(tfert_verif2(i,j,fertcount(i,j)+2)) )) THEN
                flag_fertilisation(i,j) = 1
                tfert_verif2(i,j,fertcount(i,j)+2) = .TRUE.
                nfertamm(i,j,fertcount(i,j) + 1) = 0.
                nfertnit(i,j,fertcount(i,j) + 1) = 0.
                tfert_modif(i,j,fertcount(i,j) + 1) = 500.0

                !counter for fertilizer application
                fertcount(i,j)  = fertcount(i,j)  + 1
                !mineral fertilization
                nfertammtot(i,j)  = nfertammtot(i,j)  + nfertamm(i,j,fertcount(i,j))
                nfertnittot(i,j)  = nfertnittot(i,j)  + nfertnit(i,j,fertcount(i,j))

                nfertammtotyear(i,j)  = nfertammtot(i,j)  - nfertammtotprevyear(i,j)
                nfertnittotyear(i,j)  = nfertnittot(i,j)  - nfertnittotprevyear(i,j)
                ! JCADD 20170421 for CNP merge
                ! nfertamm input is kg m-2
                grm_nfert(i,j,iammonium) = nfertamm(i,j,fertcount(i,j))*1e3
                grm_nfert(i,j,initrate) = nfertnit(i,j,fertcount(i,j))*1e3
                grm_pfert(i,j) = pfert_mineral(i,j,fertcount(i,j))*1e3
                ! END JCADD
                IF (f_nonlimitant .EQ. 1.) THEN
                  controle_azote_sum_mem(i,j) = controle_azote_sum(i,j)
                  IF (.NOT. (tfert_verif3(i,j,fertcount_start(i,j))) ) THEN
                    tfert_verif3(i,j,fertcount_start(i,j))= .TRUE.
                    controle_azote_sum(i,j) = 0.
                  END IF
                END IF

            END IF

          END IF
        END DO
      END DO

! JC comment not sure what is it for
!      WHERE (flag_fertilisation(:,:) .NE. 2)
!        flag_fertilisation(:,:) = 0
!      END WHERE
    END IF users_or_saturant_fert


    ! 4.2 ***** RUN NONLIMITANT *****
    ! recherche des erreurs pour l'�quilibre
    ! recherche de stoplimitant (fin du run)
    run_nonlimitant : IF ((f_nonlimitant .EQ. 1) .AND. (ANY(stoplimitant(:,:) .EQ. 0))) THEN   ! any ?

      DO j=2,nvm
        DO i=1,npts
          ! search the last time of fertilization
          IF (tfert(i,j,fertcount_current(i,j) + 1) .EQ. 500) THEN
              fertcount_next = 1
          ELSE
              fertcount_next = fertcount_current(i,j) + 1
          ENDIF

          ! if tjulian correspond to next time of fertilization
          IF ((tjulian .GE. tfert(i,j,fertcount_next)) .AND. &
             (tjulian .LE. tfert(i,j,fertcount_next)+0.9) .AND. &
             (tfert_verif2(i,j,fertcount_next) .EQV. .FALSE.)) THEN

              tfert_verif2(i,j,fertcount_next) = .TRUE.

              ! calcul de somerror
              IF(controle_azote(i,j,fertcount_next) .GT. 0.) THEN
                  nnonlimit_SOMerror(i,j) = &
                     (controle_azote(i,j,fertcount_next) - controle_azote_sum_mem(i,j))/ &
                     controle_azote(i,j,fertcount_next)
              ELSE
                  nnonlimit_SOMerror(i,j) = 0.
              ENDIF
              ! on regarde si on ne d�passe pas l'erreur max voulue
              ! puis on r�ajuste cette erreur max suivant dans quel cas
              ! nous sommes
              IF (nnonlimit_SOMerror(i,j) .GT. nnonlimit_SOMerrormax(i,j)) THEN

                nfertamm(i,j,fertcount_current(i,j)) = nfertamm(i,j,fertcount_current(i,j)) + 0.00125
                nfertnit(i,j,fertcount_current(i,j)) = nfertnit(i,j,fertcount_current(i,j)) + 0.00125
                ! JCADD 20170421 for CNP merge
                ! nfertamm input is kg m-2
                grm_nfert(i,j,iammonium) = nfertamm(i,j,fertcount(i,j))*1e3
                grm_nfert(i,j,initrate) = nfertnit(i,j,fertcount(i,j))*1e3
                grm_pfert(i,j) = pfert_mineral(i,j,fertcount(i,j))*1e3
                ! END JCADD
                PRINT *, '!!! apport en azote !!! pour fertcount_current = ', fertcount_current(i,j) &
                              ,' nfertamm= ',nfertamm(i,j,fertcount_current(i,j))
              ELSE
                  fertcount_current(i,j) = fertcount_current(i,j) + 1

                  IF(tfert(i,j,fertcount_current(i,j)) .EQ. 500.) THEN
                      fertcount_current(i,j) = 1
                  ENDIF

                  IF (fertcount_current(i,j) .EQ. fertcount_start(i,j)) THEN
                      nnonlimit_SOMerrormax(i,j) = nnonlimit_SOMerrormax(i,j) - n_auto(i,j)*0.05
                      n_auto(i,j) = n_auto(i,j) - 1.
                      IF ( nnonlimit_SOMerrormax(i,j) .LE. 0.) THEN
                         nnonlimit_SOMerrormax(i,j)=0.025
                      ENDIF
                      IF (n_auto(i,j) .LT. 0.) THEN
                          stoplimitant(i,j) = 1
                          print *,'*********************************'
                          print *,'stoplimitant =1 '
                          print *,'********************************'
                      ENDIF ! n_auto
                  ENDIF ! fertcount_current
              ENDIF ! nnonlimit
          ENDIF ! tjulian
        END DO ! npts
      END DO ! nvm
    END IF run_nonlimitant

    ! 4.3 run spatialize fertilization
    ! calculating organic C input into soil
    CALL fertilisation_pas_temps(&
       npts                           , &
       fertcount                      , &
       dt                             , &
       tjulian                        , &
       deltatt                         , &
       tfert                          , &
       Nliquidmanure                  , &
       nslurry                        , &
       Nsolidmanure                   , &
       fcOrganicFertmetabolic         , &
       fcOrganicFertstruct            , &
       fnOrganicFerturine             , &
       fnOrganicFertmetabolic         , &
       fnOrganicFertstruct            , &
       fpOrganicFerturine             , &
       fpOrganicFertmetabolic         , &
       fpOrganicFertstruct            , &
       c2nratiostruct                 , &
       Fert_on)

    DO j=2,nvm
      CALL Euler_funct(npts, dt, MAX(0.0,(t2m_daily - 278.15)), tasum(:,j))
    END DO

    ! 5 calculate variables that not included in ORCHIDEE
    ! liste : 
    ! * devstage              
    ! * tgrowth               
    CALL Main_appl_pre_animal(&
       npts                  , &
       dt                    , &
       tjulian               , &
       t2m_daily                    , &
       tsoil                 , &
       new_day               , &
       new_year              , &
       regcount              , &
       tcut                  , &
       devstage              , &
       tgrowth              )

    ! 6 start grazing practice
    ! 6.1 updating available litter for wild animal grazing
    !gmjc prepare litter_avail for grazing litter
    ! kg DM/m^2
!    litter_avail_totDM(:,:) = (litter_avail(:,istructural,:) + &
!       & litter_avail(:,imetabolic,:)) / (1000. * CtoDM) 
! JCADD 20170421 for CNP merge
    ! for CNP version, all non-woody litter can be eaten
    ! currently, I did not differenciate avail and non_avail in CNP version
    litter_avail_totDM(:,:) = (litter(:,imetabolic,:,iabove,icarbon) + &
                              litter(:,istructural,:,iabove,icarbon) ) / (1000.* CtoDM)
! END JCADD
    !end gmjc
    ! 6.2 grazing
!!    IF ((Type_animal.EQ.3).OR.(Type_animal.EQ.6)) THEN ! old animal module
      CALL Animaux_main(&
       npts, dt, devstage, wsh, intakemax, &
       snowfall_daily, wshtot, Animalwgrazingmin, &
       AnimalkintakeM, nel, wanimal, nanimaltot, &
       ntot, intake, urinen, faecesn, urinec, faecesc, &
       tgrowth, new_year, new_day, &
       nanimal, tanimal, danimal, &
       tcutmodel, tjulian, import_yield, &
       intakesum, intakensum, fn, c, n, leaf_frac, &
       intake_animal, intake_animalsum, &
       biomass, trampling, sr_ugb,sr_wild, &
       compt_ugb, nb_ani, grazed_frac,AnimalDiscremineQualite, &
       YIELD_RETURN,sr_ugb_init,year_count1,year_count2, & 
       grazing_litter, litter_avail_totDM, &
       intake_animal_litter, intake_litter,nb_grazingdays, &
!gmjc top 5 layer grassland soil moisture for grazing
       moiavail_daily, tmc_topgrass_daily,fc_grazing, &
       after_snow, after_wet, wet1day, wet2day, &
       snowmass_daily, t2m_daily, &
!end gmjc
       ranimal_gm, ch4_pft_gm, Fert_PRP, &
! JCADD 20170421 for CNP merge
       intakecnp,litter_intakecnp,tramplingcnp,faecescnp,urinecnp, &
       milkcnp,liveweightcnp,ranimal_c,EFch4_c,litter,when_growthinit)
! END JCADD
!!    ELSE ! new animal module
!!
!!      CALL Animaux_main_dynamic(&
!!        npts, dt, devstage                  , &
!!        intakemax, snowfall_daily, wshtot, wsh        , &
!!        nel, nanimaltot                     , &
!!        intake                              , &
!!        import_yield                        , &
!!        new_year, new_day                   , &
!!        nanimal, tanimal, danimal           , &
!!        PIYcow, PIMcow, BCSYcow             , &
!!        BCSMcow, PICcow, AGE_cow_P          , &
!!        AGE_cow_M, tcutmodel, tjulian       , &
!!        intakesum                           , &
!!        intakensum, fn, ntot, c, n,leaf_frac, &
!!        intake_animal, intake_animalsum     , &
!!        t2m_min_daily, type_animal          , &
!!        t2m_daily, intakemax, Autogestion_out      , &
!!        Forage_quantity,t2m_14              , &
!!        intake_tolerance                    , &
!!        q_max_complement                    , &
!!        biomass, urinen, faecesn, urinec,faecesc, &
!!        file_param_init,trampling,sr_ugb,sr_wild, &
!!        compt_ugb, nb_ani,grazed_frac,AnimalDiscremineQualite, &
!!        grazing_litter, nb_grazingdays) 
!!
!!    ENDIF ! Type_Animal

    ! 7 CUTTING
    ! Cutting Management: auto_fauche and user_fauche
    flag_cutting(:,:) = 0
    ! 7.1 user defined cut 
    user_fauche : IF ((f_autogestion .EQ. 0) .AND. (f_postauto .EQ. 0)) THEN

      flag_cutting(:,:) = 0
      when_growthinit_cut(:,:) = when_growthinit_cut(:,:) + dt
      lm_before(:,:)= biomass(:,:,ileaf,icarbon)

      DO k=1,nstocking

        flag_cutting(:,:) = 0

        DO j=2,nvm
          IF (is_grassland_manag(j) )THEN

            IF (ANY(tcut_verif(:,j,k) .EQV. .FALSE.)) THEN
              WHERE ((tjulian .GE. tcut(:,j,k)) .AND. (tjulian .LE. tcut(:,j,k)+0.9) .AND. &
                (tcut_verif(:,j,k) .EQV. .FALSE.))

                tcut_verif(:,j,k) = .TRUE. 
                flag_cutting(:,j) = 1
                compt_cut(:,j) = compt_cut(:,j) + 1
                when_growthinit_cut(:,j) = 0.0                  
              END WHERE
            END IF
          END IF
        END DO              
        IF (ANY (flag_cutting(:,:) .EQ. 1)) THEN 

          IF (blabla_pasim)  PRINT *, 'cutting users', tjulian

!!          CALL cutting_spa(&
!!                   npts              , &
!!                   tjulian           , &
!!                   flag_cutting      , &
!!                   wshtotcutinit     , &
!!                   lcutinit          , &
!!                   wsh               , &
!!                   wshtot            , &
!!                   wr                , &
!!                   c                 , &
!!                   n                 , &
!!                   napo              , &
!!                   nsym              , &
!!                   fn                , &
!!                   tjulian           , &
!!                   nel               , &
!!                   biomass           , &
!!                   devstage          , &
!!                   regcount          , &
!!                   gmean             , &
!!                   wc_frac                , &
!!                   wnapo             , &
!!                   wnsym             , &
!!                   wgn               , &
!!                   tasum             , &
!!                   tgrowth           , &
!!                   loss              , &
!!                   lossc             , &
!!                   lossn             , &
!!                   tlossstart        , &
!!                   lai               , &
!!                   tcut              , &
!!                   tcut_modif        , &
!!                   wshtotsum         , &
!!                   controle_azote_sum)
          CALL cutting_spa_cnp (&
            npts,tjulian,flag_cutting,wshtotcutinit, &
            biomass,loss,losscnp,grassharvestcnp, &
            devstage,tasum,tgrowth,gmean,regcount, &
            tlossstart,wshtotsum,wshtotsumprev,tcut,tcut_modif,flam_left)      


          WHERE ((wsh + wr .GT. 0.0).AND. (flag_cutting .EQ. 1)) 
            c = wc_frac / (wsh + wr)
            n = (wnapo + wnsym) / (wsh + wr) 
            fn = wgn / (wr + wsh)
            napo = wnapo / (wsh + wr)
            nsym = wnsym / (wsh + wr)
          END WHERE
               
          WHERE (wshtot + wrtot .GT. 0.0)
            ntot = (wnapo + wnsym + wgn) / (wshtot + wrtot)
          END WHERE
        END IF

      END DO !nstocking

    END IF user_fauche

    ! 7.2 auto cut
    n_day_autofauche :  IF (new_day) THEN
      DO  j=2,nvm
        CALL linreg_pasim (&
           npts          , &
           ngmean        , &
           tgmean(:,j,:)        , &
           gmean(:,j,:)         , &
           ngmean        , &
           misval        , &
           mux(:,j)           , &
           mugmean(:,j)       , &
           sigx(:,j)          , &
           sigy(:,j)          , &
           gmeanslope(:,j)    , &
           gzero(:,j)         , &
           gcor(:,j))
      END DO
      countschedule(:,:)  = 1

      auto_fauche : IF (f_autogestion .EQ. 1) THEN ! for optimalize sr_ugb and nb_ani

        flag_cutting(:,:) = 0
        when_growthinit_cut(:,:) = when_growthinit_cut(:,:) + dt
        DO j=2,nvm
          IF (is_grassland_manag(j) .AND. (.NOT. is_grassland_cut(j)) .AND. &
            (.NOT. is_grassland_grazed(j)))THEN

        ! FIRST test for automanagement > 45 days 
          WHERE((nanimal(:,j,1) .EQ. 0.0) .AND. (cuttingend(:,j) .EQ. 0) .AND. &
            (countschedule(:,j) .EQ. 1) .AND. (((tgrowth(:,j) .GE. tgrowthmin) .AND. &
            (gmean(:,j,ngmean) .GT. 0.0) .AND.(lai(:,j) .GE. 2.5)  .AND. &
            (devstage(:,j) .GT. devstagemin ) .AND. &
            (gmeanslope(:,j) .LT. gmeansloperel * mugmean(:,j)))))

            flag_cutting(:,j)  = 1
            countschedule(:,j) = countschedule(:,j)  + 1             
            compt_cut(:,j) = compt_cut(:,j) + 1
          END WHERE 
          END IF
        END DO !nvm

        ! If there is at least one point concerned (flag_cutting = 1)
        IF (ANY (flag_cutting(:,:) .EQ. 1)) THEN 
!          IF (blabla_pasim) PRINT *, 'FAUCHE AVEC METHODE NV ', tjulian
            ! There will be one fertilization the day after cutting
            ! A COURT-CIRCUITER si couplage autogestion ferti avec INN
            ! AIG 06/10/2009

            IF (f_fertilization.NE.1) THEN
              DO j=2,nvm   
                DO i=1,npts
                  IF (flag_cutting(i,j) .EQ. 1) THEN
                    tfert(i,j,regcount(i,j) + 1) = tjulian + 1
!                    print*, 'FERTILISATION AVEC METHODE NV', tjulian
                  END IF
                END DO
              END DO  
            END IF

!!            CALL cutting_spa(&
!!               npts              , &
!!               tjulian           , &
!!               flag_cutting      , &
!!               wshtotcutinit     , &
!!               lcutinit          , &
!!               wsh               , &
!!               wshtot            , &
!!               wr                , &
!!               c                 , &
!!               n                 , &
!!               napo              , &
!!               nsym              , &
!!               fn                , &
!!               tjulian           , &
!!               nel               , &
!!               biomass           , &
!!               devstage          , &
!!               regcount          , &
!!               gmean             , &
!!               wc_frac                , &
!!               wnapo             , &
!!               wnsym             , &
!!               wgn               , &
!!               tasum             , &
!!               tgrowth           , &
!!               loss              , &
!!               lossc             , &
!!               lossn             , &
!!               tlossstart        , &
!!               lai               , &
!!               tcut              , &
!!               tcut_modif        , &
!!               wshtotsum         , &
!!               controle_azote_sum)

          CALL cutting_spa_cnp (&
            npts,tjulian,flag_cutting,wshtotcutinit, &
            biomass,loss,losscnp,grassharvestcnp, &
            devstage,tasum,tgrowth,gmean,regcount, &
            tlossstart,wshtotsum,wshtotsumprev,tcut,tcut_modif,flam_left)

            ! ******************************************************
            ! update plant c n concentrations
            ! ******************************************************

            WHERE ((wsh + wr .GT. 0.0)  .AND. (flag_cutting .EQ. 1) )
                c = wc_frac / (wsh + wr)
                n = (wnapo + wnsym) / (wsh + wr) 
                fn = wgn / (wr + wsh)
                napo = wnapo / (wsh + wr)
                nsym = wnsym / (wsh + wr)
            END WHERE

            WHERE (wshtot + wrtot .GT. 0.0)
                ntot = (wnapo + wnsym + wgn) / (wshtot + wrtot)
            END WHERE

          END IF

          WHERE (flag_cutting(:,:) .EQ. 1)
            when_growthinit_cut(:,:) = 0.0
          END WHERE

        ! SECOND test for automanagement lai & accumulated temperature over shreshold
        flag_cutting(:,:) = 0

        DO j=2,nvm
          IF (is_grassland_manag(j) .AND. (.NOT. is_grassland_cut(j)) .AND. &
            (.NOT. is_grassland_grazed(j))) THEN
        
            WHERE ((countschedule(:,j) .EQ. 1) .AND. (nanimal(:,j,1) .EQ. 0.0) .AND. &
              (devstage(:,j) .LT. 2.0) .AND. (tasum(:,j) .GE. tasumrep ) .AND. &
              (lai(:,j) .GE. 2.5))

              flag_cutting(:,j) = 1

              countschedule(:,j) = countschedule(:,j)  + 1
                compt_cut(:,j) = compt_cut(:,j) + 1
            END WHERE
          END IF
        END DO !nvm      

        ! If there is at least one point concerned
        IF (ANY (flag_cutting(:,:) .EQ. 1)) THEN 
!            IF (blabla_pasim) PRINT *, 'FAUCHE AVEC METHODE NV', tjulian

            ! There will be one fertilization the day after cutting
            ! MODIF INN
            !courciruiter le calcul de tfert si f_fertiliZation = 0
            IF (f_fertilization.NE.1) THEN
              DO j=2,nvm
                DO i=1,npts
                  IF (flag_cutting(i,j) .EQ. 1) THEN
                       tfert(i,j,regcount(i,j) + 1) = tjulian + 1
!                      print*, 'FERTILISATION AVEC METHODE NV', tjulian
                   END IF
                END DO
              END DO
            END IF

!!            CALL cutting_spa(&
!!               npts              , &
!!               tjulian           , &
!!               flag_cutting      , &
!!               wshtotcutinit     , &
!!               lcutinit          , &
!!               wsh               , &
!!               wshtot            , &
!!               wr                , &
!!               c                 , &
!!               n                 , &
!!               napo              , &
!!               nsym              , &
!!               fn                , &
!!               tjulian           , &
!!               nel               , &
!!               biomass           , &
!!               devstage          , &
!!               regcount          , &
!!               gmean             , &
!!               wc_frac                , &
!!               wnapo             , &
!!               wnsym             , &
!!               wgn               , &
!!               tasum             , &
!!               tgrowth           , &
!!               loss              , &
!!               lossc             , &
!!               lossn             , &
!!               tlossstart        , &
!!               lai               , &
!!               tcut              , &
!!               tcut_modif        , &
!!               wshtotsum         , &
!!               controle_azote_sum)
          CALL cutting_spa_cnp (&
            npts,tjulian,flag_cutting,wshtotcutinit, &
            biomass,loss,losscnp,grassharvestcnp, &
            devstage,tasum,tgrowth,gmean,regcount, &
            tlossstart,wshtotsum,wshtotsumprev,tcut,tcut_modif,flam_left)

            ! ******************************************************
            ! update plant c n concentrations 
            ! ******************************************************
            WHERE ((wsh + wr .GT. 0.0) .AND. (flag_cutting .EQ. 1))
                c = wc_frac / (wsh + wr)
                n = (wnapo + wnsym) / (wsh + wr) 
                fn = wgn / (wr + wsh)
                napo = wnapo / (wsh + wr)
                nsym = wnsym / (wsh + wr)
            END WHERE

            WHERE (wshtot + wrtot .GT. 0.0)
                ntot = (wnapo + wnsym + wgn) / (wshtot + wrtot)
            END WHERE

        END IF

      WHERE (flag_cutting(:,:) .EQ. 1)
        when_growthinit_cut(:,:) = 0.0
      END WHERE

      !If there are ncut cutting, it's finish
    
      WHERE (regcount(:,:) .EQ. ncut ) 
        cuttingend(:,:) = 1
      END WHERE
     
      !end of the cutting season by snow fall

    ELSE IF ((f_postauto .EQ.1 ) .OR. (f_autogestion .EQ. 3) .OR. &
      (f_autogestion .EQ. 4) &
      !! JCMODIF 290714 for postauto 5
      .OR. (f_postauto .GE. 2)) THEN

        flag_cutting(:,:) = 0
        when_growthinit_cut(:,:) = when_growthinit_cut(:,:) + dt
! gmjc 07082016 reset lossc to zero for history writing
! NOTE:the flag_cutting will be determined twice a day, thus we cannot
! reset them in the cutting subroutine
    loss(:,:) = 0.0
    lossc(:,:) = 0.0
    lossn(:,:) = 0.0
    tlossstart(:,:) = 500.0
        ! FIRST test for automanagement
        WHERE((nanimal(:,mcut_C3,1) .EQ. 0.0) .AND. (cuttingend(:,mcut_C3) .EQ. 0) .AND. &
          (countschedule(:,mcut_C3) .EQ. 1) .AND. (((tgrowth(:,mcut_C3).GE.tgrowthmin) .AND. &
          (gmean(:,mcut_C3,ngmean).GT. 0.0) .AND. &
          (lai(:,mcut_C3) .GE. 2.5)  .AND. (devstage(:,mcut_C3) .GT. devstagemin ) .AND. &
          (gmeanslope(:,mcut_C3) .LT.gmeansloperel * mugmean(:,mcut_C3)))))

            flag_cutting(:,mcut_C3)  = 1
            countschedule(:,mcut_C3) = countschedule(:,mcut_C3)  + 1
            compt_cut(:,mcut_C3) = compt_cut(:,mcut_C3) + 1
        END WHERE

        WHERE((nanimal(:,mcut_C4,1) .EQ. 0.0) .AND. (cuttingend(:,mcut_C4) .EQ. 0).AND. &
          (countschedule(:,mcut_C4) .EQ. 1) .AND. (((tgrowth(:,mcut_C4).GE.tgrowthmin).AND. &
          (gmean(:,mcut_C4,ngmean).GT. 0.0) .AND. &
          (lai(:,mcut_C4) .GE. 2.5)  .AND. (devstage(:,mcut_C4) .GT. devstagemin ).AND. &
          (gmeanslope(:,mcut_C4) .LT.gmeansloperel * mugmean(:,mcut_C4)))))

            flag_cutting(:,mcut_C4)  = 1
            countschedule(:,mcut_C4) = countschedule(:,mcut_C4)  + 1
            compt_cut(:,mcut_C4) = compt_cut(:,mcut_C4) + 1
        END WHERE

        ! If there is at least one point concerned (flag_cutting = 1)

        IF ((ANY(flag_cutting(:,mcut_C3) .EQ. 1)) .OR. &
            (ANY(flag_cutting(:,mcut_C4) .EQ. 1))) THEN
!            IF (blabla_pasim) PRINT *, 'FAUCHE AVEC METHODE NV ', tjulian
                ! There will be one fertilization the day after cutting
                ! A COURT-CIRCUITER si couplage autogestion ferti avec INN
                ! AIG 06/10/2009

            IF (f_fertilization.NE.1) THEN
             DO j=2,nvm
               DO i=1,npts
                  IF (flag_cutting(i,j) .EQ. 1) THEN
                      tfert(i,j,regcount(i,j) + 1) = tjulian + 1
!                      print*, 'FERTILISATION AVEC METHODE NV', tjulian
                  END IF
                END DO
              END DO
            END IF
!!           CALL cutting_spa(&
!!               npts              , &
!!               tjulian           , &
!!               flag_cutting      , &
!!               wshtotcutinit     , &
!!               lcutinit          , &
!!               wsh               , &
!!               wshtot            , &
!!               wr                , &
!!               c                 , &
!!               n                 , &
!!               napo              , &
!!               nsym              , &
!!               fn                , &
!!               tjulian           , &
!!               nel               , &
!!               biomass           , &
!!               devstage          , &
!!               regcount          , &
!!               gmean             , &
!!               wc_frac                , &
!!               wnapo             , &
!!               wnsym             , &
!!               wgn               , &
!!               tasum             , &
!!               tgrowth           , &
!!               loss              , &
!!               lossc             , &
!!               lossn             , &
!!               tlossstart        , &
!!               lai               , &
!!               tcut              , &
!!               tcut_modif        , &
!!               wshtotsum         , &
!!               controle_azote_sum)
          CALL cutting_spa_cnp (&
            npts,tjulian,flag_cutting,wshtotcutinit, &
            biomass,loss,losscnp,grassharvestcnp, &
            devstage,tasum,tgrowth,gmean,regcount, &
            tlossstart,wshtotsum,wshtotsumprev,tcut,tcut_modif,flam_left)

            ! ******************************************************
            ! mise � jour des concentrations
            ! ******************************************************

            WHERE ((wsh + wr .GT. 0.0)  .AND. (flag_cutting .EQ. 1) )
                c = wc_frac / (wsh + wr)
                n = (wnapo + wnsym) / (wsh + wr)
                fn = wgn / (wr + wsh)
                napo = wnapo / (wsh + wr)
                nsym = wnsym / (wsh + wr)
            END WHERE

            WHERE (wshtot + wrtot .GT. 0.0)
                ntot = (wnapo + wnsym + wgn) / (wshtot + wrtot)
            END WHERE

        END IF

    WHERE (flag_cutting(:,mcut_C3) .EQ. 1)
        when_growthinit_cut(:,mcut_C3) = 0.0
    END WHERE
    WHERE (flag_cutting(:,mcut_C4) .EQ. 1)
        when_growthinit_cut(:,mcut_C4) = 0.0
    END WHERE

        ! SECOND test for automanagement
        flag_cutting(:,mcut_C3) = 0
        flag_cutting(:,mcut_C4) = 0

        WHERE ((countschedule(:,mcut_C3) .EQ. 1) .AND. (nanimal(:,mcut_C3,1) .EQ.0.0) .AND. &
          (devstage(:,mcut_C3) .LT. 2.0) .AND. (tasum(:,mcut_C3) .GE. tasumrep ) .AND. &
          (lai(:,mcut_C3) .GE. 2.5))

            flag_cutting(:,mcut_C3) = 1
            countschedule(:,mcut_C3) = countschedule(:,mcut_C3)  + 1
            compt_cut(:,mcut_C3) = compt_cut(:,mcut_C3) + 1
        END WHERE

        WHERE ((countschedule(:,mcut_C4) .EQ. 1) .AND. (nanimal(:,mcut_C4,1).EQ.0.0) .AND. &
          (devstage(:,mcut_C4) .LT. 2.0) .AND. (tasum(:,mcut_C4) .GE. tasumrep ).AND. &
          (lai(:,mcut_C4) .GE. 2.5))

            flag_cutting(:,mcut_C4) = 1
            countschedule(:,mcut_C4) = countschedule(:,mcut_C4)  + 1
            compt_cut(:,mcut_C4) = compt_cut(:,mcut_C4) + 1
        END WHERE

       ! If there is at least one point concerned
        IF ((ANY(flag_cutting(:,mcut_C3) .EQ. 1)) .OR. &
            (ANY(flag_cutting(:,mcut_C4) .EQ. 1))) THEN

!            IF (blabla_pasim) PRINT *, 'FAUCHE AVEC METHODE NV', tjulian

            ! There will be one fertilization the day after cutting
            ! MODIF INN
            !courciruiter le calcul de tfert si f_fertiliZation = 0
            IF (f_fertilization.NE.1) THEN
              DO j=2,nvm
                DO i=1,npts
                  IF (flag_cutting(i,j) .EQ. 1) THEN
                       tfert(i,j,regcount(i,j) + 1) = tjulian + 1
!                      print*, 'FERTILISATION AVEC METHODE NV', tjulian
                   END IF
                END DO
              END DO
            END IF

!!           CALL cutting_spa(&
!!               npts              , &
!!               tjulian           , &
!!               flag_cutting      , &
!!               wshtotcutinit     , &
!!               lcutinit          , &
!!               wsh               , &
!!               wshtot            , &
!!               wr                , &
!!               c                 , &
!!               n                 , &
!!               napo              , &
!!               nsym              , &
!!               fn                , &
!!               tjulian           , &
!!               nel               , &
!!               biomass           , &
!!               devstage          , &
!!               regcount          , &
!!               gmean             , &
!!               wc_frac                , &
!!               wnapo             , &
!!               wnsym             , &
!!               wgn               , &
!!               tasum             , &
!!               tgrowth           , &
!!               loss              , &
!!               lossc             , &
!!               lossn             , &
!!               tlossstart        , &
!!               lai               , &
!!               tcut              , &
!!               tcut_modif        , &
!!               wshtotsum         , &
!!               controle_azote_sum)
          CALL cutting_spa_cnp (&
            npts,tjulian,flag_cutting,wshtotcutinit, &
            biomass,loss,losscnp,grassharvestcnp, &
            devstage,tasum,tgrowth,gmean,regcount, &
            tlossstart,wshtotsum,wshtotsumprev,tcut,tcut_modif,flam_left)

            ! ******************************************************
            ! update plant c n concentrations
            ! ******************************************************
           WHERE ((wsh + wr .GT. 0.0) .AND. (flag_cutting .EQ. 1))
                c = wc_frac / (wsh + wr)
                n = (wnapo + wnsym) / (wsh + wr)
                fn = wgn / (wr + wsh)
                napo = wnapo / (wsh + wr)
                nsym = wnsym / (wsh + wr)
            END WHERE

            WHERE (wshtot + wrtot .GT. 0.0)
                ntot = (wnapo + wnsym + wgn) / (wshtot + wrtot)
            END WHERE

        END IF

      WHERE (flag_cutting(:,mcut_C3) .EQ. 1)
        when_growthinit_cut(:,mcut_C3) = 0.0
      END WHERE

      !If there are ncut cutting, it's finish
      WHERE (regcount(:,mcut_C3) .EQ. ncut )
        cuttingend(:,mcut_C3) = 1
      END WHERE

      WHERE (flag_cutting(:,mcut_C4) .EQ. 1)
        when_growthinit_cut(:,mcut_C4) = 0.0
      END WHERE

        !If there are ncut cutting, it's finish

        WHERE (regcount(:,mcut_C4) .EQ. ncut )
            cuttingend(:,mcut_C4) = 1
        END WHERE

        !end of the cutting season by snow fall

      END IF auto_fauche
    END IF n_day_autofauche

    ! 8 updating plant and soil variables after management practice
    ! maintenant nous allons regarder les changements que le management apporte � Orchidee. 
    ! Dans un premier temps uniquement ceux sur le Carbone vu qu'Orchidee n'a pas d'azote.
    ! ******************************************************
    ! 8.1 updating soil status 
    ! ******************************************************
!!    CALL chg_sol_bio(&
!!       npts                     , &
!!       tjulian                  , &
!!       bm_to_litter             , &
!!       litter                   , &
!!       litter_avail             , &
!!       litter_not_avail         , &
!!!       !spitfire
!!!       fuel_1hr, &
!!!       fuel_10hr, &
!!!       fuel_100hr, &
!!!       fuel_1000hr, &
!!!       !end spitfire
!!       litter_avail_totDM         , &
!!       intake_litter            , &
!!       biomass                  , &
!!       faecesc                  , &
!!       urinec                   , &
!!       fcOrganicFertmetabolic   , &
!!       fcOrganicFertstruct      , &
!!       fnOrganicFerturine       , &
!!       fnOrganicFertstruct      , &
!!       fnOrganicFertmetabolic    , &
!!       trampling                , &
!!       YIELD_RETURN             , &
!!       harvest_gm, cinput_gm)

  CALL chg_sol_bio_cnp(&
     npts                     , &
     tjulian                  , &
     bm_to_litter             , &
     litter                   , &
     litter_intakecnp         , &
     tramplingcnp             , &
     faecescnp                , &
     urinecnp                 , &
     fcOrganicFertmetabolic   , &
     fcOrganicFertstruct      , &
     fnOrganicFerturine       , &
     fnOrganicFertmetabolic   , &
     fnOrganicFertstruct      , &
     fpOrganicFerturine       , &
     fpOrganicFertmetabolic   , &
     fpOrganicFertstruct      , &
     grm_nfert                , &
     grm_pfert                , &
     grassorganicfertcnp      , &
     harvest_gm               , &
     cinput_gm)

    ! JC NOTED
    ! the grassland management has included many fluxes from or to
    ! environment (such as fertilization input, harvested biomass output,
    ! milk output, liveweight output, animal respiration output,
    ! and enteric fermentation ch4 emission)
    ! HOWEVER, it should be noted that grm_nfert excluding N in urinecnp
    !  are N fluxes
    ! that will directly goes into soil in the next time step, thus they
    ! are not included in the mass balance calculation
    ! BUT N in urinecnp comes from biomass or litter while exported to soil
    ! so should be seen as output 
     grm_masschangecnp(:,:,icarbon) = grassorganicfertcnp(:,:,icarbon) - &
          milkcnp(:,:,icarbon) - liveweightcnp(:,:,icarbon) - &
          ranimal_c(:,:) - EFch4_c(:,:) - &
          grassharvestcnp(:,:,icarbon)
     grm_masschangecnp(:,:,initrogen)   = grassorganicfertcnp(:,:,initrogen) - &
          milkcnp(:,:,initrogen) - liveweightcnp(:,:,initrogen) - &
          grassharvestcnp(:,:,initrogen) - urinecnp(:,:,initrogen)
     grm_masschangecnp(:,:,iphosphorus) = grassorganicfertcnp(:,:,iphosphorus) - &
          milkcnp(:,:,iphosphorus) - liveweightcnp(:,:,iphosphorus) - &
          grassharvestcnp(:,:,iphosphorus) - urinecnp(:,:,iphosphorus)

!jcadd calculate N2O emission
  ! Fert_sn mineral fertilizer N kg N m-2 d-1
  ! Fert_on organic fertilizer N kg N m-2 d-1
  ! Fert_PRP N in grazing excreta kg N m-2 d-1
  ! ndeposition N deposition kg N ha yr-1
  n2o_pft_gm = &
                ! Direct emission
                ((Fert_sn+Fert_on + ndeposition/1e4/365.) * n2o_EF1 + &
                Fert_PRP * n2o_EF2) + &
                ! Volatilization
                ((Fert_sn + ndeposition/1e4/365.) * n2o_FracGASF + &
                (Fert_on + Fert_PRP) * n2o_FracGASM) * n2o_EF3 + &
                ! Leaching
                (Fert_sn+Fert_on + ndeposition/1e4/365. + Fert_PRP) * &
                n2o_FracLEACH_H * n2o_EF4
!end jcadd

    lai(:,:) = biomass(:,:,ileaf,icarbon)*sla_calc(:,:)

! JC MOD024 devstage now calculated in stomate_lpj
    GRM_devstage(:,:) = devstage
    ! 8.2 write history 
    ! HISTWRITE
    CALL xios_orchidee_send_field("FERT_SN",Fert_sn)
    CALL xios_orchidee_send_field("FERT_ON",Fert_on)
    CALL xios_orchidee_send_field("FERT_PRP",Fert_PRP)
    CALL xios_orchidee_send_field("WSHTOT",wshtot)
    CALL xios_orchidee_send_field("WRTOT",wrtot)
    CALL xios_orchidee_send_field("WSHTOTSUM",wshtotsum)
    CALL xios_orchidee_send_field("SR_UGB",sr_ugb)
    CALL xios_orchidee_send_field("FCORGFERTMET",fcOrganicFertmetabolic)
    CALL xios_orchidee_send_field("FCORGFERTSTR",fcOrganicFertstruct)
    CALL xios_orchidee_send_field("LOSS",loss)
    CALL xios_orchidee_send_field("C_CUTYEARLY",C_cutyearly)
    CALL xios_orchidee_send_field("FREQUENCY_CUT",frequency_cut)
    CALL xios_orchidee_send_field("NFERT_TOTAL",N_fert_total)
!    CALL xios_orchidee_send_field("NDEP",ndeposition)
    CALL xios_orchidee_send_field("TMCGRASS_DAILY",tmc_topgrass_daily)
    CALL xios_orchidee_send_field("FC_GRAZING",fc_grazing)
    CALL xios_orchidee_send_field("CINPUT_GM",cinput_gm)
    CALL xios_orchidee_send_field("HARVEST_GM",harvest_gm)
    ! JCADD 20170421 for CNP merge
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
       CALL xios_orchidee_send_field("GRM_LOSS"//TRIM(element_str(l)),SUM(losscnp(:,:,:,l),DIM=3))
       CALL histwrite_p (hist_id_stomate, 'GRM_LOSS'//TRIM(element_str(l)), itime, &
            SUM(losscnp(:,:,:,l),DIM=3),npts*nvm, horipft_index)
       CALL xios_orchidee_send_field("GRM_HARVEST"//TRIM(element_str(l)),grassharvestcnp(:,:,l))
       CALL histwrite_p (hist_id_stomate, 'GRM_HARVEST'//TRIM(element_str(l)), itime, &
            grassharvestcnp(:,:,l),npts*nvm, horipft_index)
       CALL xios_orchidee_send_field("GRM_INTAKE"//TRIM(element_str(l)),intakecnp(:,:,l))
       CALL histwrite_p (hist_id_stomate, 'GRM_INTAKE'//TRIM(element_str(l)),itime, &
            intakecnp(:,:,l),npts*nvm, horipft_index)
       CALL xios_orchidee_send_field("GRM_LIT_INTAKE"//TRIM(element_str(l)),SUM(SUM(litter_intakecnp(:,:,:,:,l),DIM=4),DIM=2))
       CALL histwrite_p (hist_id_stomate, 'GRM_LIT_INTAKE'//TRIM(element_str(l)),itime, &
            SUM(SUM(litter_intakecnp(:,:,:,:,l),DIM=4),DIM=2),npts*nvm, horipft_index)
       CALL xios_orchidee_send_field("GRM_TRAMPLING"//TRIM(element_str(l)),SUM(tramplingcnp(:,:,:,l),DIM=3))
       CALL histwrite_p (hist_id_stomate, 'GRM_TRAMPLING'//TRIM(element_str(l)),itime, &
            SUM(tramplingcnp(:,:,:,l),DIM=3),npts*nvm, horipft_index)
       CALL xios_orchidee_send_field("GRM_MILK"//TRIM(element_str(l)),milkcnp(:,:,l))
       CALL histwrite_p (hist_id_stomate, 'GRM_MILK'//TRIM(element_str(l)),itime, &
            milkcnp(:,:,l),npts*nvm, horipft_index)
       CALL xios_orchidee_send_field("GRM_LIVEWEIGHT"//TRIM(element_str(l)),liveweightcnp(:,:,l))
       CALL histwrite_p (hist_id_stomate, 'GRM_LIVEWEIGHT'//TRIM(element_str(l)), itime, &
            liveweightcnp(:,:,l),npts*nvm, horipft_index)
       CALL xios_orchidee_send_field("GRM_FAECES"//TRIM(element_str(l)),faecescnp(:,:,l))
       CALL histwrite_p (hist_id_stomate, 'GRM_FAECES'//TRIM(element_str(l)),itime, &
            faecescnp(:,:,l),npts*nvm, horipft_index)
       CALL xios_orchidee_send_field("GRM_URINE"//TRIM(element_str(l)),urinecnp(:,:,l))
       CALL histwrite_p (hist_id_stomate, 'GRM_URINE'//TRIM(element_str(l)),itime, &
            urinecnp(:,:,l),npts*nvm, horipft_index)
       CALL xios_orchidee_send_field("GRM_ORGFERT"//TRIM(element_str(l)),grassorganicfertcnp(:,:,l))
       CALL histwrite_p (hist_id_stomate, 'GRM_ORGFERT'//TRIM(element_str(l)),itime, &
            grassorganicfertcnp(:,:,l),npts*nvm, horipft_index)

    END DO
    CALL xios_orchidee_send_field("GRM_RANIMAL_C",ranimal_c)
    CALL histwrite_p(hist_id_stomate,'GRM_RANIMAL_C',itime,ranimal_c,npts*nvm, horipft_index)
    CALL xios_orchidee_send_field("GRM_EFCH4_C",EFch4_c)
    CALL histwrite_p(hist_id_stomate,'GRM_EFCH4_C',itime,EFch4_c,npts*nvm, horipft_index)        
    CALL xios_orchidee_send_field("GRM_NFERT_AMM",grm_nfert(:,:,1))
    CALL histwrite_p(hist_id_stomate,'GRM_NFERT_AMM',itime,grm_nfert(:,:,1),npts*nvm, horipft_index)
    CALL xios_orchidee_send_field("GRM_NFERT_NIT",grm_nfert(:,:,2))
    CALL histwrite_p(hist_id_stomate,'GRM_NFERT_NIT',itime,grm_nfert(:,:,2),npts*nvm, horipft_index)
    CALL xios_orchidee_send_field("GRM_PFERT",grm_pfert(:,:))
    CALL histwrite_p(hist_id_stomate,'GRM_PFERT',itime,grm_pfert(:,:),npts*nvm, horipft_index)

    ! END JCADD
    regcount_real  = regcount
    fertcount_real = fertcount
!    CALL histwrite_p(hist_id_stomate ,'YIELD_RETURN',itime,YIELD_RETURN,npts*nvm, horipft_index)
!    CALL histwrite_p(hist_id_stomate ,'REGCOUNT' ,itime ,regcount_real , npts*nvm, horipft_index)
!    CALL histwrite_p(hist_id_stomate ,'FERTCOUNT',itime ,fertcount_real, npts*nvm, horipft_index)
!    CALL histwrite_p(hist_id_stomate ,'GMEAN1',itime ,gmean(:,:,1) ,npts*nvm, horipft_index)
!    CALL histwrite_p(hist_id_stomate ,'GMEAN2',itime ,gmean(:,:,2) ,npts*nvm, horipft_index)
!    CALL histwrite_p(hist_id_stomate ,'GMEAN3',itime ,gmean(:,:,3) ,npts*nvm, horipft_index)
!    CALL histwrite_p(hist_id_stomate ,'GMEAN4',itime ,gmean(:,:,4) ,npts*nvm, horipft_index)
!    CALL histwrite_p(hist_id_stomate ,'GMEAN5',itime ,gmean(:,:,5) ,npts*nvm, horipft_index)
!    CALL histwrite_p(hist_id_stomate ,'GMEAN6',itime ,gmean(:,:,6) ,npts*nvm, horipft_index)
!    CALL histwrite_p(hist_id_stomate ,'GMEAN7',itime ,gmean(:,:,7) ,npts*nvm, horipft_index)
!    CALL histwrite_p(hist_id_stomate ,'GMEAN8',itime ,gmean(:,:,8) ,npts*nvm, horipft_index)
!    CALL histwrite_p(hist_id_stomate ,'GMEAN9',itime ,gmean(:,:,9) ,npts*nvm, horipft_index)
!    CALL histwrite_p(hist_id_stomate ,'GMEAN0',itime ,gmean(:,:,10) ,npts*nvm, horipft_index)
    CALL histwrite_p(hist_id_stomate ,'WSH'   ,itime , wsh   , npts*nvm, horipft_index)
    CALL histwrite_p(hist_id_stomate ,'WSHTOT',itime , wshtot, npts*nvm, horipft_index)
    CALL histwrite_p(hist_id_stomate ,'WR',    itime , wr,     npts*nvm, horipft_index)
    CALL histwrite_p(hist_id_stomate ,'WRTOT', itime , wrtot,  npts*nvm, horipft_index)
    CALL histwrite_p(hist_id_stomate ,'WSHTOTSUM', itime , wshtotsum,  npts*nvm, horipft_index)
    CALL histwrite_p(hist_id_stomate ,'SR_UGB', itime , sr_ugb,  npts*nvm,horipft_index)
    ! HISTWRITE POUR LA FERTIILSATION
    CALL histwrite_p(hist_id_stomate ,'FCORGFERTMET',itime , fcOrganicFertmetabolic,npts*nvm, horipft_index)
    CALL histwrite_p(hist_id_stomate ,'FCORGFERTSTR'   ,itime , fcOrganicFertstruct   ,npts*nvm, horipft_index)
    !CALL histwrite_p(hist_id_stomate ,'FNORGFERTURINE'    ,itime , fnOrganicFerturine    ,npts*nvm, horipft_index)
    !CALL histwrite_p(hist_id_stomate ,'FNORGFERTSTR'   ,itime , fnOrganicFertstruct   ,npts*nvm, horipft_index)
    !CALL histwrite_p(hist_id_stomate ,'FNORGFERTMET',itime , fnOrganicFertmetabolic,npts*nvm, horipft_index)
    !CALL histwrite_p(hist_id_stomate ,'NFERTNITTOT'          ,itime , nfertnit(:,:,1)    ,npts*nvm, horipft_index)
    !CALL histwrite_p(hist_id_stomate ,'NFERTAMMTOT'          ,itime , nfertamm(:,:,1)    ,npts*nvm, horipft_index)
    ! HISTWRITE POUR LA FAUCHE
    CALL histwrite_p(hist_id_stomate ,'LOSS' ,itime ,loss  ,npts*nvm, horipft_index)
!    CALL histwrite_p(hist_id_stomate ,'LOSSC',itime ,lossc ,npts*nvm, horipft_index)
!    CALL histwrite_p(hist_id_stomate ,'LOSSN',itime ,lossn ,npts*nvm, horipft_index)
    CALL histwrite_p(hist_id_stomate ,'DM_CUTYEARLY',itime ,DM_cutyearly ,npts*nvm, horipft_index)
    CALL histwrite_p(hist_id_stomate ,'C_CUTYEARLY',itime ,C_cutyearly ,npts*nvm, horipft_index)
    CALL histwrite_p(hist_id_stomate ,'COMPT_CUT' ,itime ,compt_cut  ,npts*nvm, horipft_index)
    CALL histwrite_p(hist_id_stomate ,'FREQUENCY_CUT',itime ,frequency_cut ,npts*nvm, horipft_index)
    CALL histwrite_p(hist_id_stomate ,'NFERT_TOTAL',itime ,N_fert_total ,npts*nvm, horipft_index)
!    CALL histwrite_p(hist_id_stomate ,'NDEP',itime ,ndeposition ,npts*nvm,horipft_index)
!    CALL histwrite_p(hist_id_stomate ,'LEGUME_FRACTION',itime ,legume_fraction ,npts*nvm, horipft_index)
!    CALL histwrite_p(hist_id_stomate ,'SOIL_FERTILITY',itime ,soil_fertility ,npts*nvm, horipft_index)
!    CALL histwrite_p(hist_id_stomate ,'C'       ,itime, c       , npts*nvm, horipft_index)
!    CALL histwrite_p(hist_id_stomate ,'N'       ,itime, n       , npts*nvm, horipft_index)
!    CALL histwrite_p(hist_id_stomate ,'FN'      ,itime, fn      , npts*nvm, horipft_index)
!    CALL histwrite_p(hist_id_stomate ,'NTOT'    ,itime, ntot    , npts*nvm, horipft_index)
!    CALL histwrite_p(hist_id_stomate ,'NAPO'    ,itime, napo    , npts*nvm, horipft_index)
!    CALL histwrite_p(hist_id_stomate ,'NSYM'    ,itime, nsym    , npts*nvm, horipft_index)
!    CALL histwrite_p(hist_id_stomate ,'DEVSTAGE',itime, devstage, npts*nvm, horipft_index)
    CALL histwrite_p(hist_id_stomate ,'TGROWTH' ,itime, tgrowth , npts*nvm, horipft_index)
!    CALL histwrite_p(hist_id_stomate ,'GRAZINGCSTRUCT',itime, grazingcstruct      , npts*nvm, horipft_index)
!    CALL histwrite_p(hist_id_stomate ,'GRAZINGNSTRUCT',itime, grazingnstruct      , npts*nvm, horipft_index)
!    CALL histwrite_p(hist_id_stomate ,'GRAZINGWN'     ,itime, Substrate_grazingwn, npts*nvm, horipft_index)
!    CALL histwrite_p(hist_id_stomate ,'GRAZINGWC'     ,itime, Substrate_grazingwc, npts*nvm, horipft_index)
!gmjc top 5 layer grassland soil moisture for grazing
    CALL histwrite_p(hist_id_stomate ,'TMCGRASS_DAILY',itime,tmc_topgrass_daily, npts, hori_index)
    CALL histwrite_p(hist_id_stomate ,'FC_GRAZING',itime,fc_grazing, npts, hori_index)
!end gmjc
  END SUBROUTINE main_grassland_management

  ! modules calculating devstage and tgrowth for grazing
  ! liste of functions calculated
  ! - devstage
  ! - tgrowth
  ! - dndfi
  SUBROUTINE Main_appl_pre_animal(&
     npts                  , &
     dt                    , &
     tjulian               , &
     t2m_daily                    , &
     tsoil                 , &
     new_day               , &
     new_year              , &
     regcount              , &
     tcut                  , &
     devstage              , &
     tgrowth               )

    INTEGER (i_std)                      , INTENT(in)  :: npts
    LOGICAL                              , INTENT(in)  :: new_day
    LOGICAL                              , INTENT(in)  :: new_year
    REAL                          , INTENT(in)  :: dt
    INTEGER                       , INTENT(in)  :: tjulian
    REAL, DIMENSION(npts)         , INTENT(in)  :: t2m_daily
    ! air temperature (K)
    REAL, DIMENSION(npts)          , INTENT(in)  :: tsoil
    ! soil surface temperature
    REAL, DIMENSION(npts,nvm,nstocking), INTENT(in)  :: tcut
    INTEGER   , DIMENSION(npts,nvm) , INTENT(in)  :: regcount
    REAL, DIMENSION(npts,nvm)       , INTENT(out) :: devstage
    ! state of developpement of growth
    REAL, DIMENSION(npts,nvm)       , INTENT(out) :: tgrowth
    ! time from last cut (d)

    INTEGER :: ier
    REAL, DIMENSION(npts)      :: xtmp_npts
    IF (new_year) THEN
        tcut0(:,:) = 0.0
    END IF

    CALL cal_devstage(npts, dt, t2m_daily, tsoil, new_day, &
            new_year, regcount, devstage)
    CALL cal_tgrowth(npts, dt, devstage, tjulian, new_day, &
            new_year, regcount, tcut, tgrowth)


  END SUBROUTINE Main_appl_pre_animal

  ! module calculating devstage
  SUBROUTINE cal_devstage(&
                npts,dt,t2m_daily,tsoil,new_day, &
                new_year, regcount, devstage)

    INTEGER (i_std)                   , INTENT(in)  :: npts
    REAL                 , INTENT(in)  :: dt
    REAL, DIMENSION(npts), INTENT(in)  :: t2m_daily
    REAL, DIMENSION(npts), INTENT(in)  :: tsoil
    LOGICAL                    , INTENT(in)  :: new_day
    LOGICAL                    , INTENT(in)  :: new_year
    INTEGER   , DIMENSION(npts,nvm), INTENT(in)  :: regcount
    REAL, DIMENSION(npts,nvm), INTENT(inout) :: devstage

    INTEGER :: i,j

    CALL Euler_funct(npts,dt,t2m_daily,tacumm)
    CALL Euler_funct(npts,dt,tsoil,tsoilcumm)

!    CALL histwrite_p(hist_id_stomate ,'TSOILCUMM',itime,tsoilcumm, npts, hori_index)


    IF (new_day) THEN

      tamean1(:)  = tamean2(:)
      tamean2(:)  = tamean3(:)
      tamean3(:)  = tamean4(:)
      tamean4(:)  = tamean5(:)
      tamean5(:)  = tamean6(:)
      tamean6(:)  = tameand(:)

      tameand(:) = tacumm(:) - tacummprev(:)
      tacummprev(:) = tacumm(:)

      tameanw(:)  = (&
         tamean1(:) + &
         tamean2(:) + &
         tamean3(:) + &
         tamean4(:) + &
         tamean5(:) + &
         tamean6(:) + &
         tameand(:))/7.0

      tsoilmeand(:) = tsoilcumm(:) - tsoilcummprev(:)
      tsoilcummprev(:) = tsoilcumm(:)

! JC MOD024 it is now calculated in stomate_lpj
!!      DO j=2,nvm
!!        DO i=1,npts
!!
!!          IF ((devstage(i,j) .LE. 0.0) .AND. ( (tameanw(i) .GT. trep) .OR. &
!!             (regcount(i,j) .EQ. 2) ) ) THEN
!!
!!              devstage(i,j) = MAX(0.0, tameand(i) - tbase)/tasumrep
!!
!!          ELSEIF ((devstage(i,j) .GT. 0.0) .AND. &
!!                 (tsoilmeand(i) .GT. tbase) .AND. &
!!                 (devstage(i,j) .LT. 2.0) ) THEN
!!
!!              devstage(i,j) = devstage(i,j) + MAX(0.0, tameand(i) - &
!!                              tbase)/tasumrep
!!
!!          ELSE
!!              devstage(i,j) = devstage(i,j)
!!
!!          ENDIF
!!        END DO ! npts
!!      END DO ! nvm
    END IF

!!    IF (new_year) THEN
!!
!!      devstage(:,:) = 0.0
!!
!!    END IF

  END SUBROUTINE cal_devstage

  ! module calculating tgrowth
  SUBROUTINE cal_tgrowth(&
                npts, dt, devstage, tjulian, new_day, &
                new_year, regcount, tcut, tgrowth)

    INTEGER                        , INTENT(in)  :: npts
    REAL                           , INTENT(in)  :: dt
    INTEGER                        , INTENT(in)  :: tjulian    ! julien day (d)
    LOGICAL                               , INTENT(in)  :: new_day
    LOGICAL                               , INTENT(in)  :: new_year
    REAL, DIMENSION(npts,nvm)      , INTENT(in)  :: devstage   ! state of developpement
    INTEGER   , DIMENSION(npts,nvm), INTENT(in)  :: regcount   ! number of cut
    REAL, DIMENSION(npts,nvm,nstocking), INTENT(in)  :: tcut   ! cut date
    REAL, DIMENSION(npts,nvm)      , INTENT(out) :: tgrowth    ! regrowth time after last cut (d)

    INTEGER :: i,j

    IF (new_day) THEN

      ! TGROWTH
      !(robson, m. j. et al., 1988)
      DO j=2,nvm
        WHERE ((devstage(:,j) .GT. 0.0) .AND. (tcut0(:,j) .LE. 0.0))

          tcut0(:,j)  = FLOAT(tjulian)

        END WHERE
      END DO

      DO j=2,nvm
        DO i=1,npts
          IF ((regcount(i,j) .EQ. 1) .AND. (tcut0(i,j) .LE. 0.0)) THEN

            tgrowth(i,j)  = 0.0
          ELSEIF (regcount(i,j) .EQ. 1) THEN

            tgrowth(i,j) = tjulian  - tcut0(i,j)

          ELSE

            tgrowth(i,j) = tjulian - tcut(i,j,regcount(i,j)-1)

          ENDIF
        END DO ! npts
      END DO ! nvm
    END IF

    IF (new_year) THEN
      tgrowth(:,:) = 0.0
    END IF
  END SUBROUTINE cal_tgrowth

  ! module updating soil status
  SUBROUTINE chg_sol_bio(&
     npts                     , &
     tjulian                  , &
     bm_to_litter             , &
     litter                   , &
     litter_avail             , &
     litter_not_avail         , &
!     !spitfire
!     fuel_1hr, &
!     fuel_10hr, &
!     fuel_100hr, &
!     fuel_1000hr, &
!     !end spitfire
     litter_avail_totDM         , &
     intake_litter            , &
     biomass                  , &
     faecesc                  , &
     urinec                   , &
     fcOrganicFertmetabolic    , &       
     fcOrganicFertstruct       , &
     fnOrganicFerturine        , &
     fnOrganicFertstruct       , &
     fnOrganicFertmetabolic    , &
     trampling                 , &
     YIELD_RETURN              , &
     harvest_gm, cinput_gm)

    INTEGER                                , INTENT(in)   :: npts
    INTEGER                             , INTENT(in)   :: tjulian                 ! jour julien    
    REAL, DIMENSION(npts,nvm,nparts,nelements), INTENT(inout):: bm_to_litter 
    REAL, DIMENSION(npts,nlitt,nvm,nlevs,nelements), INTENT(inout):: litter
    REAL, DIMENSION(npts,nlitt,nvm), INTENT(inout):: litter_avail 
    REAL, DIMENSION(npts,nlitt,nvm), INTENT(inout):: litter_not_avail 
!    !spitfire
!    REAL, DIMENSION(npts,nvm,nlitt),INTENT(inout)        :: fuel_1hr
!    REAL, DIMENSION(npts,nvm,nlitt),INTENT(inout)        :: fuel_10hr
!    REAL, DIMENSION(npts,nvm,nlitt),INTENT(inout)        :: fuel_100hr
!    REAL, DIMENSION(npts,nvm,nlitt),INTENT(inout)        :: fuel_1000hr
!    !end spitfire
    REAL, DIMENSION(npts,nvm), INTENT(inout):: litter_avail_totDM
    REAL, DIMENSION(npts,nvm), INTENT(in):: intake_litter
    REAL, DIMENSION(npts,nvm)            , INTENT(in)   :: faecesc
    REAL, DIMENSION(npts,nvm)            , INTENT(in)   :: urinec
    REAL, DIMENSION(npts,nvm,nparts,nelements), INTENT(in)   :: biomass           
    ! totalit� de masse s�che du shoot(kg/m**2)
    REAL, DIMENSION(npts,nvm)            , INTENT(in)   :: fcOrganicFertmetabolic
    ! metabolic C in slurry and manure (kg C/m**2/d)
    REAL, DIMENSION(npts,nvm)            , INTENT(in)   :: fcOrganicFertstruct 
    ! structural C in slurry and manure (kg C/m**2/d)
    REAL, DIMENSION(npts,nvm)            , INTENT(in)   :: fnOrganicFerturine    
    ! urine N in slurry and manure (kg N/m**2/d)
    REAL, DIMENSION(npts,nvm)            , INTENT(in)   :: fnOrganicFertstruct   
    ! structural N in slurry and manure (kg N/m**2/d)
    REAL, DIMENSION(npts,nvm)            , INTENT(in)   :: fnOrganicFertmetabolic  
    ! metabolic N in slurry and manure (kg N/m**2/d)           
    REAL, DIMENSION(npts,nvm)            , INTENT(in)   :: trampling
    REAL, DIMENSION(npts,nvm)            , INTENT(inout)   :: YIELD_RETURN

    REAL, DIMENSION(npts,nvm), INTENT(inout)  :: harvest_gm
    REAL, DIMENSION(npts,nvm), INTENT(inout)  :: cinput_gm

    REAL, DIMENSION(npts,nvm) :: litter_avail_totDM_old
    REAL, DIMENSION(npts,nvm) :: fcloss  
    REAL, DIMENSION(npts,nvm) :: fnloss
    REAL, DIMENSION(npts,nvm) :: floss
    REAL, DIMENSION(npts,nvm) :: fcplantsoil
    REAL, DIMENSION(npts,nvm) :: fnplantsoil
    REAL, DIMENSION(npts,nvm) :: fplantsoil
    REAL, DIMENSION(npts,nvm) :: l2nratio
    REAL, DIMENSION(npts,nvm) :: fmetabolic
    REAL, DIMENSION(npts,nvm) :: manure_barn
    INTEGER :: j
    REAL, PARAMETER       :: yieldloss    = 0.05
!    !spitfire
!    REAL, DIMENSION(npts,nvm,nlitt)       :: fuel_all_type
!    REAL, DIMENSION(npts,nvm,nlitt,4)     :: fuel_type_frac
!    !end spitfire

    IF (blabla_pasim) PRINT *, 'PASIM main grassland : call chg_sol_bio'
    fmetabolic = 0.0
    
!    !spitfire
!    fuel_type_frac(:,:,:,:) = zero
!    fuel_all_type(:,:,:) = fuel_1hr(:,:,:) + fuel_10hr(:,:,:) + &
!                           fuel_100hr(:,:,:) + fuel_1000hr(:,:,:)
!    WHERE (fuel_all_type(:,:,:) .GT. min_stomate)
!      fuel_type_frac(:,:,:,1) = fuel_1hr(:,:,:)/fuel_all_type(:,:,:)
!      fuel_type_frac(:,:,:,2) = fuel_10hr(:,:,:)/fuel_all_type(:,:,:)
!      fuel_type_frac(:,:,:,3) = fuel_100hr(:,:,:)/fuel_all_type(:,:,:)
!      fuel_type_frac(:,:,:,4) = fuel_1000hr(:,:,:)/fuel_all_type(:,:,:)
!    ENDWHERE
!    !end spitfire   

    DO j=2,nvm
      ! tlossstart will be set to tjulian when cutting is trigered at the point
      ! deltatt = 1.000000
      WHERE ((tjulian .GE. tlossstart(:,j)) .AND. &
            (tjulian .LT. (tlossstart(:,j) + deltatt/2.0))) 

        fcloss(:,j) = lossc(:,j)/deltatt 
        fnloss(:,j) = lossn(:,j)/deltatt 
        floss(:,j)  = loss(:,j) /deltatt 
        ! loss is 5% of total harvest, 95% is exported out of ecosystem
        ! kgC m-2 day-1 -> g C m-2 day-1
        harvest_gm(:,j) = fcloss(:,j)/yieldloss*(1-yieldloss)*1e3

      ELSEWHERE

        fcloss(:,j) = zero
        fnloss(:,j) = zero
        floss(:,j)  = zero
        harvest_gm(:,j) = zero
      END WHERE

      cinput_gm(:,j) = (fcOrganicFertstruct(:,j)+fcOrganicFertmetabolic(:,j))*1e3 

      fcplantsoil(:,j) = fcloss(:,j)
      fnplantsoil(:,j) = fnloss(:,j)
      fplantsoil(:,j)  = floss(:,j) 

      WHERE (fnplantsoil(:,j) .GT. 0.0) 

        l2nratio(:,j) = fligninresidue * fplantsoil(:,j)/ fnplantsoil(:,j)

      ELSEWHERE

        l2nratio(:,j) = 0.0

      END WHERE

      IF (is_grassland_cut(j).AND.(.NOT.is_grassland_grazed(j)))THEN
       
      ! Manure produced at barn
      ! 0.05 yieldloss 0.95 import_yield/harvest 0.85 loss during trasportation 
      ! 0.12 fraction of manure spread to field in total intake dry matter at barn (0.85*0.95harvest) 
      !JCcomment for accounting for Manurefert only not return 
      !        WHERE (YIELD_RETURN(:,:) .GT. 0.0)
      !          manure_barn(:,j) = fcplantsoil(:,j) / 0.05 * 0.95 * 0.85 *0.12 +&
      !            YIELD_RETURN(:,:) * CtoDM 
      !          YIELD_RETURN(:,:) = 0.0
      !        ELSEWHERE
      !          manure_barn(:,j) = fcplantsoil(:,j) / 0.05 * 0.95 * 0.85 *0.12
      !        ENDWHERE  
      manure_barn(:,j) = 0.0

      ELSE 
      manure_barn(:,j) = 0.0
      END IF
 
      fmetabolic(:,j) = MAX(0.625,MIN(0.85 - 0.018 * l2nratio(:,j), 1.0 - fligninresidue))
      bm_to_litter(:,j,ileaf,icarbon) = bm_to_litter(:,j,ileaf,icarbon) + &
         & fmetabolic(:,j) * (fcplantsoil(:,j) * 1000.0 + trampling(:,j))

      bm_to_litter(:,j,isapabove,icarbon) = bm_to_litter(:,j,isapabove,icarbon) + & 
         & (1.0 - fmetabolic(:,j)) * (fcplantsoil(:,j) * 1000.0 + trampling(:,j)) 
      litter_avail_totDM_old(:,j) = litter_avail_totDM(:,j)
      ! new litter available tot DM after intake litter
      litter_avail_totDM(:,j) = litter_avail_totDM(:,j) - intake_litter(:,j)
      IF (ANY(litter_avail_totDM(:,j) .LT. -0.01 ) ) THEN
        WRITE(numout,*) 'zd ','litter avail', j, litter_avail_totDM_old(:,j)
        WRITE(numout,*) 'zd ','intake litter', j, intake_litter(:,j)
        STOP 'available litter is not enough for grazing'

      ENDIF
      ! litter available C left is recalculated 
      ! assuming the same structural and metabolic fraction    
      WHERE (litter_avail_totDM_old(:,j) .GT. 0.0 )
      litter_avail(:,istructural,j) = litter_avail(:,istructural,j) * &
            & (litter_avail_totDM(:,j)/litter_avail_totDM_old(:,j))
      litter_avail(:,imetabolic,j) = litter_avail(:,imetabolic,j) * &
            & (litter_avail_totDM(:,j)/litter_avail_totDM_old(:,j))
      ELSEWHERE
      litter_avail(:,istructural,j) = litter_avail(:,istructural,j)
      litter_avail(:,imetabolic,j) = litter_avail(:,imetabolic,j)
      ENDWHERE
      ! new litter not available after manure/urine

      litter_not_avail(:,istructural,j) = litter_not_avail(:,istructural,j) + &
            & (faecesc(:,j) + urinec(:,j) + manure_barn(:,j) ) * 1000.0 * (1.0 - fmetabolic(:,j)) + &
            &  fcOrganicFertstruct(:,j) * 1000.0

      litter_not_avail(:,imetabolic,j) = litter_not_avail(:,imetabolic,j) + &
            & (faecesc(:,j) + urinec(:,j) + manure_barn(:,j) ) * 1000.0 * fmetabolic(:,j) + &
            &  fcOrganicFertmetabolic(:,j) * 1000.0
      ! update litter
      litter(:,:,j,iabove,icarbon) = litter_avail(:,:,j) + litter_not_avail(:,:,j)
!      !spitfire
!      fuel_1hr(:,j,:) = litter(:,:,j,iabove,icarbon) * fuel_type_frac(:,j,:,1) 
!      fuel_10hr(:,j,:) = litter(:,:,j,iabove,icarbon) * fuel_type_frac(:,j,:,2) 
!      fuel_100hr(:,j,:) = litter(:,:,j,iabove,icarbon) * fuel_type_frac(:,j,:,3) 
!      fuel_1000hr(:,j,:) = litter(:,:,j,iabove,icarbon) * fuel_type_frac(:,j,:,4) 
!      !endspit

    END DO
  END SUBROUTINE chg_sol_bio

  ! JCADD 20170421 for CNP merge
  SUBROUTINE chg_sol_bio_cnp(&
     npts                     , &
     tjulian                  , &
     bm_to_litter             , &
     litter                   , &
     litter_intakecnp         , &
     tramplingcnp             , &
     faecescnp                , &
     urinecnp                 , &
     fcOrganicFertmetabolic   , &
     fcOrganicFertstruct      , &
     fnOrganicFerturine       , &
     fnOrganicFertmetabolic   , &
     fnOrganicFertstruct      , &
     fpOrganicFerturine       , &
     fpOrganicFertmetabolic   , &
     fpOrganicFertstruct      , &
     grm_nfert                , &
     grm_pfert                , &
     grassorganicfertcnp      , &
     harvest_gm               , &
     cinput_gm)

    INTEGER                                , INTENT(in)   :: npts
    INTEGER                         , INTENT(in)   :: tjulian ! jour julien     
    REAL, DIMENSION(npts,nvm,nparts,nelements), INTENT(inout) :: bm_to_litter 
    REAL, DIMENSION(npts,nlitt,nvm,nlevs,nelements), INTENT(inout) :: litter
    REAL, DIMENSION(npts,nlitt,nvm,nlevs,nelements), INTENT(in) :: litter_intakecnp
    REAL, DIMENSION(npts,nvm,nparts,nelements), INTENT(in)   :: tramplingcnp
    REAL, DIMENSION(npts,nvm,nelements), INTENT(in)   :: faecescnp
    REAL, DIMENSION(npts,nvm,nelements), INTENT(in)   :: urinecnp
    REAL, DIMENSION(npts,nvm)            , INTENT(in)   :: fcOrganicFertmetabolic
    ! metabolic C in slurry and manure (kg C/m**2/d)
    REAL, DIMENSION(npts,nvm)            , INTENT(in)   :: fcOrganicFertstruct
    ! structural C in slurry and manure (kg C/m**2/d)
    REAL, DIMENSION(npts,nvm)            , INTENT(in)   :: fnOrganicFerturine
    ! urine N in slurry and manure (kg N/m**2/d)
    REAL, DIMENSION(npts,nvm)            , INTENT(in)   :: fnOrganicFertstruct
    ! structural N in slurry and manure (kg N/m**2/d)
    REAL, DIMENSION(npts,nvm)            , INTENT(in)   :: fnOrganicFertmetabolic
    ! metabolic N in slurry and manure (kg N/m**2/d)
    REAL, DIMENSION(npts,nvm)            , INTENT(in)   :: fpOrganicFerturine
    ! urine P  in slurry and manure (kg N/m**2/d)
    REAL, DIMENSION(npts,nvm)            , INTENT(in)   :: fpOrganicFertstruct
    ! structural P in slurry and manure (kg N/m**2/d)
    REAL, DIMENSION(npts,nvm)            , INTENT(in)   :: fpOrganicFertmetabolic
    ! metabolic P in slurry and manure (kg N/m**2/d)
    REAL, DIMENSION(npts,nvm,ninput), INTENT(inout)  :: grm_nfert
    REAL, DIMENSION(npts,nvm), INTENT(inout)  :: grm_pfert
    REAL, DIMENSION(npts,nvm,nelements), INTENT(inout) :: grassorganicfertcnp
    REAL, DIMENSION(npts,nvm), INTENT(inout)  :: harvest_gm
    REAL, DIMENSION(npts,nvm), INTENT(inout)  :: cinput_gm

    INTEGER :: j
    REAL, PARAMETER       :: yieldloss    = 0.05
    REAL, DIMENSION(npts,nvm) :: fmetabolic
    REAL, DIMENSION(npts) :: l2nratio

    l2nratio = 0.0
    fmetabolic = 0.0

    DO j=2,nvm
      ! tlossstart will be set to tjulian when cutting is trigered at the point
      ! deltatt = 1.000000
      WHERE ((tjulian .GE. tlossstart(:,j)) .AND. &
            (tjulian .LT. (tlossstart(:,j) + deltatt/2.0)))
        ! grassharvestcnp is directly calculated from biomass gCm-2
        harvest_gm(:,j) = grassharvestcnp(:,j,icarbon)
      ELSEWHERE
        harvest_gm(:,j) = zero
      ENDWHERE

      ! l2nratio is calculated from n in dry matter bm_to_litter
      ! i.e., Lignin:nitrogen ratio of plant residues
      ! from PaSim, which is different from that of ORCHIDEE in stomate_litter
      ! if choosing LC_sapabove - 0.09 and CN = 20, the fmetabolic will always
      ! be very close to 0.85 metabolic_ref_frac 
      WHERE (bm_to_litter(:,j,isapabove,initrogen) .GT. 0.0)
        l2nratio(:) = fligninresidue * bm_to_litter(:,j,isapabove,icarbon)/0.45 &
                      / bm_to_litter(:,j,isapabove,initrogen)
      ELSEWHERE
        l2nratio(:) = 0.0
      ENDWHERE

      ! fligninresidue = 0.25
      ! fmetabolic is between 0.625 and 0.75
      fmetabolic(:,j) = MAX(0.625,MIN(0.85 - 0.018 * l2nratio(:), 1.0 - fligninresidue))

      ! update bm_to_litter from losscnp and tramplingcnp
      bm_to_litter(:,j,:,:) = bm_to_litter(:,j,:,:) + losscnp(:,j,:,:) + &
                              tramplingcnp(:,j,:,:)
      ! update grm_nfert from fnOrganicFerturine and urinecnp(:,:,initrogen)
      ! they goes into soil ammoniacal N directly (derived from PaSim)
      ! urinecnp g CNP m-2 d-1, fnOrganicFerturine kg N m-2 d-1
      grm_nfert(:,j,iammonium) = grm_nfert(:,j,iammonium) + &
                urinecnp(:,j,initrogen) + fnOrganicFerturine(:,j)*1e3

      grm_pfert(:,j) = grm_pfert(:,j) + &
                urinecnp(:,j,iphosphorus) + fpOrganicFerturine(:,j)*1e3
      ! calculate organic fertilizer that will be put as litter
      ! which should be accounted for in mass balance calculation
      grassorganicfertcnp(:,j,icarbon) = fcOrganicFertmetabolic(:,j)*1e3 + &
           fcOrganicFertstruct(:,j)*1e3
      grassorganicfertcnp(:,j,initrogen) = fnOrganicFertmetabolic(:,j)*1e3 + &
           fnOrganicFertstruct(:,j)*1e3 
      grassorganicfertcnp(:,j,iphosphorus) = fpOrganicFertmetabolic(:,j)*1e3 + &
           fpOrganicFertstruct(:,j)*1e3
      ! update litter from 
      ! litter_intakecnp
      litter(:,:,:,:,:) = litter(:,:,:,:,:) - litter_intakecnp(:,:,:,:,:)

      ! Organic fertiliser (metabolic struct) c n p?
      ! faecescnp, urinec and p
      ! litter C metabolic from Organic fert, faeces metabolic part and urine 
      litter(:,imetabolic,j,iabove,icarbon) = litter(:,imetabolic,j,iabove,icarbon) + &
             fcOrganicFertmetabolic(:,j)*1e3 + fmetabolic(:,j) * faecescnp(:,j,icarbon) + &
             urinecnp(:,j,icarbon)
      ! litter N metabolic from Organic fert, faeces metabolic part
      ! urinen directly as amm in grm_nfert
      litter(:,imetabolic,j,iabove,initrogen) = litter(:,imetabolic,j,iabove,initrogen) + &
             fnOrganicFertmetabolic(:,j)*1e3 + fmetabolic(:,j) * faecescnp(:,j,initrogen)
      ! litter P metabolic from faeces metabolic part
      litter(:,imetabolic,j,iabove,iphosphorus) = litter(:,imetabolic,j,iabove,iphosphorus) + &
             fpOrganicFertmetabolic(:,j)*1e3 + fmetabolic(:,j) * faecescnp(:,j,iphosphorus)      

      ! litter C struct from Organic fert, faeces struct part
      litter(:,istructural,j,iabove,icarbon) = litter(:,istructural,j,iabove,icarbon) + &
             fcOrganicFertstruct(:,j)*1e3 + (un - fmetabolic(:,j)) * faecescnp(:,j,icarbon)
      ! litter N struct from Organic fert, faeces struct part
      litter(:,istructural,j,iabove,initrogen) = litter(:,istructural,j,iabove,initrogen) + &
             fnOrganicFertstruct(:,j)*1e3 + (un - fmetabolic(:,j)) * faecescnp(:,j,initrogen)
      ! litter P struct from faeces struct part
      litter(:,istructural,j,iabove,iphosphorus) = litter(:,istructural,j,iabove,iphosphorus) + &
             fpOrganicFertstruct(:,j)*1e3 + (un - fmetabolic(:,j)) * faecescnp(:,j,iphosphorus)


    END DO ! nvm
  END SUBROUTINE chg_sol_bio_cnp

  ! clear memory used by grassland management module
  SUBROUTINE grassmanag_clear
    IF (ALLOCATED(intake)) DEALLOCATE(intake)
    IF (ALLOCATED(intakemax)) DEALLOCATE(intakemax)
    IF (ALLOCATED(intake_litter)) DEALLOCATE(intake_litter)
    IF (ALLOCATED(intake_animal_litter)) DEALLOCATE(intake_animal_litter)
    IF (ALLOCATED(grazing_litter)) DEALLOCATE(grazing_litter)
    IF (ALLOCATED(litter_avail_totDM)) DEALLOCATE(litter_avail_totDM)
    IF (ALLOCATED(wshtotcutinit)) DEALLOCATE(wshtotcutinit)
    IF (ALLOCATED(lcutinit)) DEALLOCATE(lcutinit)
    IF (ALLOCATED(devstage)) DEALLOCATE(devstage)
    IF (ALLOCATED(faecesc)) DEALLOCATE(faecesc)
    IF (ALLOCATED(faecesn)) DEALLOCATE(faecesn)
    IF (ALLOCATED(urinen)) DEALLOCATE(urinen)
    IF (ALLOCATED(urinec)) DEALLOCATE(urinec)
    IF (ALLOCATED(nel)) DEALLOCATE(nel)
    IF (ALLOCATED(nanimaltot)) DEALLOCATE(nanimaltot)
    IF (ALLOCATED(tgrowth)) DEALLOCATE(tgrowth)
    IF (ALLOCATED(wsh)) DEALLOCATE(wsh)
    IF (ALLOCATED(wshtot)) DEALLOCATE(wshtot)
    IF (ALLOCATED(wshtotinit)) DEALLOCATE(wshtotinit)
    IF (ALLOCATED(wr)) DEALLOCATE(wr)
    IF (ALLOCATED(wrtot)) DEALLOCATE(wrtot)
    IF (ALLOCATED(wanimal)) DEALLOCATE(wanimal)
    IF (ALLOCATED(ntot)) DEALLOCATE(ntot)
    IF (ALLOCATED(c)) DEALLOCATE(c)
    IF (ALLOCATED(n)) DEALLOCATE(n)
    IF (ALLOCATED(fn)) DEALLOCATE(fn)
    IF (ALLOCATED(napo)) DEALLOCATE(napo)
    IF (ALLOCATED(nsym)) DEALLOCATE(nsym)
    IF (ALLOCATED(wnapo)) DEALLOCATE(wnapo)
    IF (ALLOCATED(wnsym)) DEALLOCATE(wnsym)
    IF (ALLOCATED(wn)) DEALLOCATE(wn)
    IF (ALLOCATED(nanimal)) DEALLOCATE(nanimal)
    IF (ALLOCATED(tanimal)) DEALLOCATE(tanimal)
    IF (ALLOCATED(danimal)) DEALLOCATE(danimal)
    IF (ALLOCATED(tcut)) DEALLOCATE(tcut)
    IF (ALLOCATED(tfert)) DEALLOCATE(tfert)
    IF (ALLOCATED(Nliquidmanure)) DEALLOCATE(Nliquidmanure)
    IF (ALLOCATED(nslurry)) DEALLOCATE(nslurry)
    IF (ALLOCATED(Nsolidmanure)) DEALLOCATE(Nsolidmanure)
    IF (ALLOCATED(legume_fraction)) DEALLOCATE(legume_fraction)
    IF (ALLOCATED(soil_fertility)) DEALLOCATE(soil_fertility)
    IF (ALLOCATED(Animalwgrazingmin)) DEALLOCATE(Animalwgrazingmin)
    IF (ALLOCATED(AnimalkintakeM)) DEALLOCATE(AnimalkintakeM)
    IF (ALLOCATED(AnimalDiscremineQualite)) DEALLOCATE(AnimalDiscremineQualite)
    IF (ALLOCATED(controle_azote)) DEALLOCATE(controle_azote)
    IF (ALLOCATED(fcOrganicFertmetabolicsum)) DEALLOCATE(fcOrganicFertmetabolicsum)
    IF (ALLOCATED(fcOrganicFertstructsum)) DEALLOCATE(fcOrganicFertstructsum)
    IF (ALLOCATED(fnOrganicFertmetabolicsum)) DEALLOCATE(fnOrganicFertmetabolicsum)
    IF (ALLOCATED(fnOrganicFertstructsum)) DEALLOCATE(fnOrganicFertstructsum)
    IF (ALLOCATED(fnOrganicFerturinesum)) DEALLOCATE(fnOrganicFerturinesum)
    IF (ALLOCATED(fnatmsum)) DEALLOCATE(fnatmsum)
    IF (ALLOCATED(controle_azote_sum)) DEALLOCATE(controle_azote_sum)
    IF (ALLOCATED(nfertamm)) DEALLOCATE(nfertamm)
    IF (ALLOCATED(nfertnit)) DEALLOCATE(nfertnit)
    IF (ALLOCATED(intakesum)) DEALLOCATE(intakesum)
    IF (ALLOCATED(intakensum)) DEALLOCATE(intakensum)
    IF (ALLOCATED(intake_animal)) DEALLOCATE(intake_animal)
    IF (ALLOCATED(intake_animalsum)) DEALLOCATE(intake_animalsum)
    IF (ALLOCATED(PIYcow)) DEALLOCATE(PIYcow)
    IF (ALLOCATED(PIMcow)) DEALLOCATE(PIMcow)
    IF (ALLOCATED(BCSYcow)) DEALLOCATE(BCSYcow)
    IF (ALLOCATED(BCSMcow)) DEALLOCATE(BCSMcow)
    IF (ALLOCATED(PICcow)) DEALLOCATE(PICcow)
    IF (ALLOCATED(AGE_cow_P)) DEALLOCATE(AGE_cow_P)
    IF (ALLOCATED(AGE_cow_M)) DEALLOCATE(AGE_cow_M)
    IF (ALLOCATED(Autogestion_out)) DEALLOCATE(Autogestion_out)
    IF (ALLOCATED(Forage_quantity)) DEALLOCATE(Forage_quantity)
    IF (ALLOCATED(tcut_modif)) DEALLOCATE(tcut_modif)
    IF (ALLOCATED(countschedule)) DEALLOCATE(countschedule)
    IF (ALLOCATED(mux)) DEALLOCATE(mux)
    IF (ALLOCATED(mugmean)) DEALLOCATE(mugmean)
    IF (ALLOCATED(sigx)) DEALLOCATE(sigx)
    IF (ALLOCATED(sigy)) DEALLOCATE(sigy)
    IF (ALLOCATED(gmeanslope)) DEALLOCATE(gmeanslope)
    IF (ALLOCATED(gzero)) DEALLOCATE(gzero)
    IF (ALLOCATED(gcor)) DEALLOCATE(gcor)
    IF (ALLOCATED(cuttingend)) DEALLOCATE(cuttingend)
    IF (ALLOCATED(tcut_verif)) DEALLOCATE(tcut_verif)
    IF (ALLOCATED(tfert_verif)) DEALLOCATE(tfert_verif)
    IF (ALLOCATED(tfert_verif2)) DEALLOCATE(tfert_verif2)
    IF (ALLOCATED(tfert_verif3)) DEALLOCATE(tfert_verif3)
    IF (ALLOCATED(regcount)) DEALLOCATE(regcount)
!    IF (ALLOCATED(wshcutinit)) DEALLOCATE(wshcutinit)
    IF (ALLOCATED(gmean)) DEALLOCATE(gmean)
    IF (ALLOCATED(tgmean)) DEALLOCATE(tgmean)
    IF (ALLOCATED(wc_frac)) DEALLOCATE(wc_frac)
    IF (ALLOCATED(wgn)) DEALLOCATE(wgn)
    IF (ALLOCATED(tasum)) DEALLOCATE(tasum)
    IF (ALLOCATED(loss)) DEALLOCATE(loss)
    IF (ALLOCATED(lossc)) DEALLOCATE(lossc)
    IF (ALLOCATED(lossn)) DEALLOCATE(lossn)
    IF (ALLOCATED(tlossstart)) DEALLOCATE(tlossstart)
    IF (ALLOCATED(flag_fertilisation)) DEALLOCATE(flag_fertilisation)
    IF (ALLOCATED(fertcount)) DEALLOCATE(fertcount)
    IF (ALLOCATED(c2nratiostruct)) DEALLOCATE(c2nratiostruct)
    IF (ALLOCATED(nfertammtot)) DEALLOCATE(nfertammtot)
    IF (ALLOCATED(nfertnittot)) DEALLOCATE(nfertnittot)
    IF (ALLOCATED(nfertammtotyear)) DEALLOCATE(nfertammtotyear)
    IF (ALLOCATED(nfertnittotyear)) DEALLOCATE(nfertnittotyear)
    IF (ALLOCATED(nfertammtotprevyear)) DEALLOCATE(nfertammtotprevyear)
    IF (ALLOCATED(nfertnittotprevyear)) DEALLOCATE(nfertnittotprevyear)
    IF (ALLOCATED(fcOrganicFertmetabolic)) DEALLOCATE(fcOrganicFertmetabolic)
    IF (ALLOCATED(fcOrganicFertstruct)) DEALLOCATE(fcOrganicFertstruct)
    IF (ALLOCATED(fnOrganicFerturine)) DEALLOCATE(fnOrganicFerturine)
    IF (ALLOCATED(fnOrganicFertstruct)) DEALLOCATE(fnOrganicFertstruct)
    IF (ALLOCATED(fnOrganicFertmetabolic)) DEALLOCATE(fnOrganicFertmetabolic)
    IF (ALLOCATED(nsatur_somerror_temp)) DEALLOCATE(nsatur_somerror_temp)
    IF (ALLOCATED(nsatur_somerror)) DEALLOCATE(nsatur_somerror)
    IF (ALLOCATED(tfert_modif)) DEALLOCATE(tfert_modif)
    IF (ALLOCATED(nnonlimit_SOMerror)) DEALLOCATE(nnonlimit_SOMerror)
    IF (ALLOCATED(nnonlimit_SOMerrormax)) DEALLOCATE(nnonlimit_SOMerrormax)
    IF (ALLOCATED(controle_azote_sum_mem)) DEALLOCATE(controle_azote_sum_mem)
    IF (ALLOCATED(n_auto)) DEALLOCATE(n_auto)
    IF (ALLOCATED(stoplimitant)) DEALLOCATE(stoplimitant)
    IF (ALLOCATED(fertcount_start)) DEALLOCATE(fertcount_start)
    IF (ALLOCATED(fertcount_current)) DEALLOCATE(fertcount_current)
    IF (ALLOCATED(wshtotsumprev)) DEALLOCATE(wshtotsumprev)
    IF (ALLOCATED(fertil_year)) DEALLOCATE(fertil_year)
    IF (ALLOCATED(toto)) DEALLOCATE(toto)
    IF (ALLOCATED(apport_azote)) DEALLOCATE(apport_azote)
    IF (ALLOCATED(trampling)) DEALLOCATE(trampling)
    IF (ALLOCATED(wshtotsumprevyear)) DEALLOCATE(wshtotsumprevyear)
    IF (ALLOCATED(file_management)) DEALLOCATE(file_management)
    IF (ALLOCATED(tmp_sr_ugb_C3)) DEALLOCATE(tmp_sr_ugb_C3)
    IF (ALLOCATED(tmp_nb_ani_C3)) DEALLOCATE(tmp_nb_ani_C3)
    IF (ALLOCATED(tmp_grazed_frac_C3)) DEALLOCATE(tmp_grazed_frac_C3)
    IF (ALLOCATED(tmp_import_yield_C3)) DEALLOCATE(tmp_import_yield_C3)
    IF (ALLOCATED(tmp_wshtotsum_C3)) DEALLOCATE(tmp_wshtotsum_C3)
    IF (ALLOCATED(tmp_sr_ugb_C4)) DEALLOCATE(tmp_sr_ugb_C4)
    IF (ALLOCATED(tmp_nb_ani_C4)) DEALLOCATE(tmp_nb_ani_C4)
    IF (ALLOCATED(tmp_grazed_frac_C4)) DEALLOCATE(tmp_grazed_frac_C4)
    IF (ALLOCATED(tmp_import_yield_C4)) DEALLOCATE(tmp_import_yield_C4)
    IF (ALLOCATED(tmp_wshtotsum_C4)) DEALLOCATE(tmp_wshtotsum_C4)
    IF (ALLOCATED(DM_cutyearly)) DEALLOCATE(DM_cutyearly)
    IF (ALLOCATED(C_cutyearly)) DEALLOCATE(C_cutyearly)
    IF (ALLOCATED(YIELD_RETURN)) DEALLOCATE(YIELD_RETURN)
    IF (ALLOCATED(sr_ugb_init)) DEALLOCATE(sr_ugb_init)
    IF (ALLOCATED(N_fert_total)) DEALLOCATE(N_fert_total)
    IF (ALLOCATED(ndeposition)) DEALLOCATE(ndeposition)
    IF (ALLOCATED(compt_cut)) DEALLOCATE(compt_cut)
    IF (ALLOCATED(frequency_cut)) DEALLOCATE(frequency_cut)
    IF (ALLOCATED(sr_wild)) DEALLOCATE(sr_wild)
    IF (ALLOCATED(flag_cutting)) DEALLOCATE(flag_cutting)
    ! from applic_plant
    IF (ALLOCATED(tamean1)) DEALLOCATE(tamean1)
    IF (ALLOCATED(tamean2)) DEALLOCATE(tamean2)
    IF (ALLOCATED(tamean3)) DEALLOCATE(tamean3)
    IF (ALLOCATED(tamean4)) DEALLOCATE(tamean4)
    IF (ALLOCATED(tamean5)) DEALLOCATE(tamean5)
    IF (ALLOCATED(tamean6)) DEALLOCATE(tamean6)
    IF (ALLOCATED(tameand)) DEALLOCATE(tameand)
    IF (ALLOCATED(tameanw)) DEALLOCATE(tameanw)
    IF (ALLOCATED(tacumm)) DEALLOCATE(tacumm)
    IF (ALLOCATED(tacummprev)) DEALLOCATE(tacummprev)
    IF (ALLOCATED(tsoilcumm)) DEALLOCATE(tsoilcumm)
    IF (ALLOCATED(tsoilcummprev)) DEALLOCATE(tsoilcummprev)
    IF (ALLOCATED(tsoilmeand)) DEALLOCATE(tsoilmeand)
    IF (ALLOCATED(tcut0)) DEALLOCATE(tcut0)
    IF (ALLOCATED(Fert_sn)) DEALLOCATE(Fert_sn)
    IF (ALLOCATED(Fert_on)) DEALLOCATE(Fert_on)
    IF (ALLOCATED(Fert_PRP)) DEALLOCATE(Fert_PRP)
    ! JCADD 20170421 for CNP merge
    IF (ALLOCATED(losscnp)) DEALLOCATE(losscnp)
    IF (ALLOCATED(grassharvestcnp)) DEALLOCATE(grassharvestcnp)
    IF (ALLOCATED(intakecnp)) DEALLOCATE(intakecnp)
    IF (ALLOCATED(litter_intakecnp)) DEALLOCATE(litter_intakecnp)
    IF (ALLOCATED(tramplingcnp)) DEALLOCATE(tramplingcnp)
    IF (ALLOCATED(faecescnp)) DEALLOCATE(faecescnp)
    IF (ALLOCATED(urinecnp)) DEALLOCATE(urinecnp)
    IF (ALLOCATED(milkcnp)) DEALLOCATE(milkcnp)
    IF (ALLOCATED(liveweightcnp)) DEALLOCATE(liveweightcnp)
    IF (ALLOCATED(ranimal_c)) DEALLOCATE(ranimal_c)
    IF (ALLOCATED(EFch4_c)) DEALLOCATE(EFch4_c)
    IF (ALLOCATED(grassorganicfertcnp)) DEALLOCATE(grassorganicfertcnp)
    IF (ALLOCATED(fpOrganicFerturine)) DEALLOCATE(fpOrganicFerturine)
    IF (ALLOCATED(fpOrganicFertstruct)) DEALLOCATE(fpOrganicFertstruct)
    IF (ALLOCATED(fpOrganicFertmetabolic)) DEALLOCATE(fpOrganicFertmetabolic)
    IF (ALLOCATED(pfert_mineral)) DEALLOCATE(pfert_mineral)
    ! END JCADD

    CALL animal_clear

  END SUBROUTINE grassmanag_clear

  ! ________________________________________________________________
  ! Functions read management.dat text file (not used for non-site simulation
  ! ________________________________________________________________

  SUBROUTINE reading_new_animal(&
           npts           , &
           nb_year_management, &
           tcutmodel      , &
           tcut           , &
           tfert          , &
           nfertamm       , &
           nfertnit       , &
           nanimal        , &
           tanimal        , &
           danimal        , &
           nliquidmanure  , &
           nslurry        , &
           nsolidmanure   , &
           PIYcow         , &
           PIMcow         , &
           BCSYcow        , &
           BCSMcow        , &
           PICcow         , &
           AGE_cow_P      , &
           AGE_cow_M      , &
           Forage_quantity)

    INTEGER (i_std)                             , INTENT(in)  :: npts
    INTEGER                              , INTENT(in) :: tcutmodel
    INTEGER,DIMENSION(nvm)             , INTENT(in) :: nb_year_management
    REAL, DIMENSION(npts,nvm,nstocking), INTENT(out) :: tcut
    REAL, DIMENSION(npts,nvm,nstocking), INTENT(out) :: tfert
    REAL, DIMENSION(npts,nvm,nstocking), INTENT(out) :: nfertamm
    REAL, DIMENSION(npts,nvm,nstocking), INTENT(out) :: nfertnit
    REAL, DIMENSION(npts,nvm,nstocking), INTENT(out) :: nanimal
    REAL, DIMENSION(npts,nvm,nstocking), INTENT(out) :: tanimal
    REAL, DIMENSION(npts,nvm,nstocking), INTENT(out) :: danimal
    REAL, DIMENSION(npts,nvm,nstocking), INTENT(out) :: nliquidmanure
    REAL, DIMENSION(npts,nvm,nstocking), INTENT(out) :: nslurry
    REAL, DIMENSION(npts,nvm,nstocking), INTENT(out) :: nsolidmanure
    REAL, DIMENSION(npts,nvm,nstocking), INTENT(out) :: PIYcow
    REAL, DIMENSION(npts,nvm,nstocking), INTENT(out) :: PIMcow
    REAL, DIMENSION(npts,nvm,nstocking), INTENT(out) :: BCSYcow
    REAL, DIMENSION(npts,nvm,nstocking), INTENT(out) :: BCSMcow
    REAL, DIMENSION(npts,nvm,nstocking), INTENT(out) :: PICcow
    REAL, DIMENSION(npts,nvm,nstocking), INTENT(out) :: AGE_cow_P
    REAL, DIMENSION(npts,nvm,nstocking), INTENT(out) :: AGE_cow_M
    REAL, DIMENSION(npts,nvm,nstocking), INTENT(out) :: Forage_quantity

    REAL, DIMENSION(nstocking)          :: nanimal_t
    REAL, DIMENSION(nstocking)          :: tanimal_t
    REAL, DIMENSION(nstocking)          :: danimal_t
    REAL, DIMENSION(nstocking)          :: tcut_t
    REAL, DIMENSION(nstocking)          :: tfert_t
    REAL, DIMENSION(nstocking)          :: nfertamm_t
    REAL, DIMENSION(nstocking)          :: nfertnit_t
    REAL, DIMENSION(nstocking)          :: nliquidmanure_t
    REAL, DIMENSION(nstocking)          :: nslurry_t
    REAL, DIMENSION(nstocking)          :: nsolidmanure_t

    REAL, DIMENSION(nstocking)          :: PIYcow_t
    REAL, DIMENSION(nstocking)          :: PIMcow_t
    REAL, DIMENSION(nstocking)          :: BCSYcow_t
    REAL, DIMENSION(nstocking)          :: BCSMcow_t
    REAL, DIMENSION(nstocking)          :: PICcow_t
    REAL, DIMENSION(nstocking)          :: AGE_cow_P_t
    REAL, DIMENSION(nstocking)          :: AGE_cow_M_t
    REAL, DIMENSION(nstocking)          :: Forage_quantity_t

    INTEGER            :: ier, i, year, fin,j
    CHARACTER(len=200) :: description

    DO j=2,nvm
      IF (is_grassland_manag(j)) THEN

      OPEN(unit=60, file = file_management(j))

      READ(60, *   , iostat = ier) description
      read_management : IF (tcutmodel .EQ. 0) THEN
        IF (blabla_pasim) PRINT *, 'USERS MANAGEMENT'

        IF (nb_year_management(j) .LT. 1 ) STOP 'error with the nb_year_management'

        IF (MOD(count_year,nb_year_management(j))  .EQ. 0) THEN
            fin = nb_year_management(j)
        ELSE
            fin = MOD(count_year,nb_year_management(j))
        END IF

        DO year = 1, fin
            READ(60, *, iostat=ier) tcut_t(:)
            READ(60, *, iostat=ier) tfert_t(:)
            READ(60, *, iostat=ier) nfertamm_t(:)
            READ(60, *, iostat=ier) nfertnit_t(:)
            READ(60, *, iostat=ier) nanimal_t(:)
            READ(60, *, iostat=ier) tanimal_t(:)
            READ(60, *, iostat=ier) danimal_t(:)
            READ(60, *, iostat=ier) nliquidmanure_t(:)
            READ(60, *, iostat=ier) nslurry_t(:)
            READ(60, *, iostat=ier) nsolidmanure_t(:)

            READ(60, *, iostat=ier) PIYcow_t(:)
            READ(60, *, iostat=ier) PIMcow_t(:)
            READ(60, *, iostat=ier) BCSYcow_t(:)
            READ(60, *, iostat=ier) BCSMcow_t(:)
            READ(60, *, iostat=ier) PICcow_t(:)
            READ(60, *, iostat=ier) AGE_cow_P_t(:)
            READ(60, *, iostat=ier) AGE_cow_M_t(:)
            READ(60, *, iostat=ier) Forage_quantity_t(:)
          DO i=1,npts
            nanimal(i,j,:)=nanimal_t(:)
            tanimal(i,j,:)=tanimal_t(:)
            danimal(i,j,:)=danimal_t(:)
            tcut(i,j,:)=tcut_t(:)
            tfert(i,j,:)=tfert_t(:)
            ! for management.dat file 
            ! fertilizations are all kgN m-2 
            nfertamm(i,j,:)=nfertamm_t(:)
            nfertnit(i,j,:)=nfertnit_t(:)
            nliquidmanure(i,j,:)=nliquidmanure_t(:)
            nslurry(i,j,:)=nslurry_t(:)
            nsolidmanure(i,j,:)=nsolidmanure_t(:)

            PIYcow(i,j,:)=PIYcow_t(:)
            PIMcow(i,j,:)=PIMcow_t(:)
            BCSYcow(i,j,:)=BCSYcow_t(:)
            BCSMcow(i,j,:)=BCSMcow_t(:)
            PICcow(i,j,:)=PICcow_t(:)
            AGE_cow_P(i,j,:)=AGE_cow_P_t(:)
            AGE_cow_M(i,j,:)=AGE_cow_M_t(:)
            Forage_quantity(i,j,:)=Forage_quantity_t(:)

          END DO
        END DO
        
      WRITE(numout,*) 'fin',fin
      WRITE(numout,*) 'tcut_t',tcut_t(:)
      WRITE(numout,*) 'tfert_t',tfert_t(:)
      WRITE(numout,*) 'nfertamm_t',nfertamm_t(:)
      WRITE(numout,*) 'nanimal_t',nanimal_t(:)
      WRITE(numout,*) 'tanimal_t',tanimal_t(:)
      WRITE(numout,*) 'danimal_t',danimal_t(:)
      WRITE(numout,*) 'nliquidmanure_t',nliquidmanure_t(:)
      WRITE(numout,*) 'nslurry_t',nslurry_t(:)
      WRITE(numout,*) 'nsolidmanure_t',nsolidmanure_t(:)

      ELSE IF (tcutmodel .EQ. 1) THEN

        PRINT *, 'AUTO MANAGEMENT'
        READ(61, *, iostat = ier) toto(:)
        READ(61, *, iostat = ier) toto(:)
        READ(61, *, iostat = ier) toto(:)
        READ(61, *, iostat = ier) toto(:)

        READ(60,     *, iostat=ier)    nanimal_t
        DO i=1,npts
          nanimal(i,j,:)=nanimal_t(:)
        END DO
      ELSE

        STOP 'PASIM ERROR :: tcutmodel must be 0 or 1'

      END IF read_management
      CLOSE(60)

    ENDIF ! is_grassland_manag
    END DO !nvm
  END SUBROUTINE reading_new_animal

  ! subroutine for reading management from map nc file

  SUBROUTINE reading_map_manag(&
           npts, lalo, neighbours, resolution, contfrac, &
           count_year     , &
           nb_year_management, &
           management_intensity, &
           management_start, &
           tcut           , &
           tfert          , &
           nfertamm       , &
           nfertnit       , &
           pfert_mineral  , &
           nanimal        , &
           tanimal        , &
           danimal        , &
           nliquidmanure  , &
           nslurry        , &
           nsolidmanure   , &
           legume_fraction, &
           soil_fertility , &
           deposition_start, &
           ndeposition, &
           sr_ugb, &
           sr_wild)

    INTEGER (i_std)                             , INTENT(in)  :: npts
    INTEGER,DIMENSION(npts,8),INTENT(in) :: neighbours        !!Neighoring grid points if land for the DGVM
                                                                         !!(unitless)
    REAL,DIMENSION(npts,2),INTENT(in)    :: lalo              !!Geographical coordinates (latitude,longitude)
                                                                         !! for pixels (degrees)
    REAL,DIMENSION(npts,2),INTENT(in)    :: resolution        !! Size in x an y of the grid (m) - surface area of
                                                                         !! the gridbox
    REAL,DIMENSION (npts), INTENT (in)   :: contfrac          !! Fraction of continent in the grid cell (unitless)
    INTEGER (i_std)                             , INTENT(in)  :: count_year
    INTEGER,DIMENSION(nvm)             , INTENT(in) :: nb_year_management
    INTEGER,DIMENSION(nvm)             , INTENT(in) :: management_intensity
    INTEGER,DIMENSION(nvm)             , INTENT(in) :: management_start
    REAL, DIMENSION(npts,nvm,nstocking), INTENT(out) :: tcut
    REAL, DIMENSION(npts,nvm,nstocking), INTENT(out) :: tfert
    REAL, DIMENSION(npts,nvm,nstocking), INTENT(out) :: nfertamm
    REAL, DIMENSION(npts,nvm,nstocking), INTENT(out) :: nfertnit
    REAL, DIMENSION(npts,nvm,nstocking), INTENT(out) :: pfert_mineral
    REAL, DIMENSION(npts,nvm,nstocking), INTENT(out) :: nanimal
    REAL, DIMENSION(npts,nvm,nstocking), INTENT(out) :: tanimal
    REAL, DIMENSION(npts,nvm,nstocking), INTENT(out) :: danimal
    REAL, DIMENSION(npts,nvm,nstocking), INTENT(out) :: nliquidmanure
    REAL, DIMENSION(npts,nvm,nstocking), INTENT(out) :: nslurry
    REAL, DIMENSION(npts,nvm,nstocking), INTENT(out) :: nsolidmanure
    REAL, DIMENSION(npts,nvm)          , INTENT(out) :: legume_fraction
    REAL, DIMENSION(npts,nvm)          , INTENT(out) :: soil_fertility
    INTEGER,DIMENSION(nvm)             , INTENT(in) :: deposition_start
    REAL, DIMENSION(npts,nvm)          , INTENT(out) :: ndeposition
    REAL, DIMENSION(npts,nvm)          , INTENT(inout) :: sr_ugb
    REAL, DIMENSION(npts,nvm)          , INTENT(inout) :: sr_wild
    ! new variables for get map of management
    REAL, DIMENSION(npts)        :: nfert_temp
    REAL, DIMENSION(npts)        :: nmanure_temp
    REAL, DIMENSION(npts)        :: nanimal_temp
    REAL, DIMENSION(npts)        :: tcut_temp
    REAL, DIMENSION(npts)        :: grazing_temp
    REAL, DIMENSION(npts)        :: wild_temp
    INTEGER                      :: management_year
    INTEGER                      :: deposition_year
    REAL, ALLOCATABLE,DIMENSION(:,:,:)        :: manage1
    REAL, ALLOCATABLE,DIMENSION(:,:,:)        :: manage2
    REAL, ALLOCATABLE,DIMENSION(:,:,:)        :: manage3
    REAL, ALLOCATABLE,DIMENSION(:,:,:)        :: manage4
    REAL, ALLOCATABLE,DIMENSION(:,:,:)        :: manage5
    REAL, ALLOCATABLE,DIMENSION(:,:,:)        :: manage6
    INTEGER            :: j
    ! CROP spec
    INTEGER                                        :: yrlen
    CHARACTER(LEN=30)                                     :: strManage
    CHARACTER(LEN=30)                                     :: strVar1
    CHARACTER(LEN=30)                                     :: strVar2
    CHARACTER(LEN=30)                                     :: strVar3
    CHARACTER(LEN=30)                                     :: strVar4
    CHARACTER(LEN=30)                                     :: strVar5
    CHARACTER(LEN=30)                                     :: strVar6
      !initialize variables
      tcut(:,:,:) = 500.0
      tfert(:,:,:) = 500.0
      nfertamm(:,:,:) = 0.0
      nfertnit(:,:,:) = 0.0
      pfert_mineral(:,:,:) = 0.0
      nanimal(:,:,:) = 0.0
      tanimal(:,:,:) = 500.0
      danimal(:,:,:) = 0.0
      nliquidmanure(:,:,:) = 0.0
      nslurry(:,:,:) = 0.0 
      nsolidmanure(:,:,:) = 0.0

      legume_fraction(:,:) =0.0
      soil_fertility(:,:) = 1.0
      ndeposition(:,:) = 0.0 

      nfert_temp(:) =0.0
      nmanure_temp(:) =0.0
      nanimal_temp(:) = 0.0
      grazing_temp(:) =0.0
      wild_temp(:) =0.0
     
      !JCMODIF to avoid read nc file many times (cost lot of CPU time)
      ! modify the processes to read only once the file and save all variables
      ! then put the temporary variables to PFTs that need it
      ! though in this case the ManagInput module is only for postauto = 5
      ! rather than general reading
      yrlen=1
      ! read file
      ! fixed variable name
      ! if run non-global simulation, should still present all 6 variables
      strManage = "GRM_MANAGEMENT_MAP"
      strVar1 = "Ndep"
      strVar2 = "Nmanure"
      strVar3 = "Nmineral"
      strVar4 = "Tfert"
      strVar5 = "sr_ugb"
      strVar6 = "sr_wild"

      CALL slowproc_GRM_ManageInput(npts,lalo,neighbours,resolution,contfrac, &
               strManage,strVar1,manage1,strVar2,manage2,strVar3,manage3,&
                         strVar4,manage4,strVar5,manage5,strVar6,manage6,yrlen)
        ! gmjc add this for grids fail to read grid (fopt=0)
        WHERE (manage1 .EQ. val_exp)
          manage1 = 0.0
        ENDWHERE
        WHERE (manage2 .EQ. val_exp)
          manage2 = 0.0
        ENDWHERE
        WHERE (manage3 .EQ. val_exp)
          manage3 = 0.0
        ENDWHERE
        ! Tfert default is 500 (no Tfert)
        WHERE (manage4 .EQ. val_exp)
          manage4 = 500.0
        ENDWHERE
        WHERE (manage5 .EQ. val_exp)
          manage5 = 0.0
        ENDWHERE
        WHERE (manage6 .EQ. val_exp)
          manage6 = 0.0
        ENDWHERE

      DO j=2,nvm
        ! NOT necessary in reading 2D management since they are all 1 year per
        ! file
        ! IF (nb_year_management(j) .LT. 1 ) STOP 'error with the nb_year_management'
        ! get which year of management should be read 
        IF (MOD(count_year,nb_year_management(j))  .EQ. 0) THEN
            management_year = nb_year_management(j) + management_start(j)-1
            deposition_year = nb_year_management(j) + deposition_start(j)-1
        ELSE
            management_year = MOD(count_year,nb_year_management(j)) + management_start(j)-1
            deposition_year = MOD(count_year,nb_year_management(j)) + deposition_start(j)-1
        END IF
        WRITE(numout,*)  management_year,deposition_year
        !!!! read deposition global file for all grassland including nature
        IF ( (.NOT. is_tree(j)) .AND. natural(j) .AND. (f_deposition_map .EQ. 1)) THEN
          ndeposition(:,j)=manage1(:,1,1) 
        ELSE
          ndeposition(:,j)=0.0
        ENDIF
        !!!! read fertilization global file
        IF (management_intensity(j) .EQ. 4) THEN
          nslurry(:,j,1)=manage2(:,1,1)/10000.
          nfertamm(:,j,1)=0.5*manage3(:,1,1)/10000.         
          nfertnit(:,j,1)=0.5*manage3(:,1,1)/10000. 
          ! tfert at global scale is not defined, set to 1st April
          tfert(:,j,1) = manage4(:,1,1)    
!JC comment for now the tfert is only available for grid with mineral N input
          tfert(:,j,1)=90
        ENDIF
        !!!! read sr_ugb global file
        IF (f_postauto .EQ. 5 .AND. f_grazing_map .EQ. 1) THEN
          sr_ugb(:,mgraze_C3)=manage5(:,1,1)/10000.
          sr_ugb(:,mgraze_C4)=manage5(:,1,1)/10000.
          WHERE (sr_ugb(:,mgraze_C3) .GT. 0.001)
            sr_ugb(:,mgraze_C3) = 0.001
            sr_ugb(:,mgraze_C4) = 0.001
          END WHERE
          if (ANY(sr_ugb(:,mgraze_C3) .EQ. 0.001)) then 
          print *, 'error sr_ugb',sr_ugb(:,mgraze_C3)
          endif
          !!!! read sr_wild global file wild animal density 
          !!!! only natural grassland will be grazed by wild animal
          !!!! only when f_autogestion = 5 or f_postauto = 5
          IF ((.NOT. is_tree(j)) .AND. natural(j) .AND. &
             & (.NOT. is_grassland_cut(j)) .AND. (.NOT.is_grassland_grazed(j))) THEN
            sr_wild(:,j)=manage6(:,1,1)/10000.
          ENDIF
        ENDIF

        IF (management_intensity(j) .EQ. 1) THEN
        ! low intensity of management in Europe ;  NOT used anymore
        ELSEIF (management_intensity(j) .EQ. 2) THEN
        ! middle intensity of management in Europe ; for Leip et al., data only
          nslurry(:,j,1)=manage2(:,1,1)/10000.
          nfertamm(:,j,1)=0.5*manage3(:,1,1)/10000.
          nfertnit(:,j,1)=0.5*manage3(:,1,1)/10000.
          tfert(:,j,1)=90!manage4(:,1,1)
        ELSEIF (management_intensity(j) .EQ. 3) THEN
        ! high intensity of management in Europe ;  NOT used anymore
        ENDIF

      END DO ! nvm
    END SUBROUTINE reading_map_manag
!!
!!    ! subrouting calculate Nitrogen effect to vcmax
!!    SUBROUTINE calc_N_limfert(&
!!             npts,nfertamm, nfertnit,&
!!             nliquidmanure, nslurry, nsolidmanure,&
!!             legume_fraction,soil_fertility,ndeposition,&
!!             N_fert_total,N_limfert)
!!
!!    INTEGER (i_std)                             , INTENT(in)  :: npts
!!    REAL, DIMENSION(npts,nvm,nstocking), INTENT(in) :: nfertamm
!!    REAL, DIMENSION(npts,nvm,nstocking), INTENT(in) :: nfertnit
!!    REAL, DIMENSION(npts,nvm,nstocking), INTENT(in) :: nliquidmanure
!!    REAL, DIMENSION(npts,nvm,nstocking), INTENT(in) :: nslurry
!!    REAL, DIMENSION(npts,nvm,nstocking), INTENT(in) :: nsolidmanure
!!    REAL, DIMENSION(npts,nvm)          , INTENT(in) :: legume_fraction
!!    REAL, DIMENSION(npts,nvm)          , INTENT(in) :: soil_fertility
!!    REAL, DIMENSION(npts,nvm)          , INTENT(out) :: N_fert_total
!!    REAL, DIMENSION(npts,nvm)          , INTENT(out) :: N_limfert
!!    REAL, DIMENSION(npts,nvm)          , INTENT(in) :: ndeposition
!!
!!    INTEGER :: k,j
!!
!!      N_fert_total(:,:) = 0.0 
!!      DO k=1,nstocking
!!        N_fert_total(:,:) = (N_fert_total(:,:) + nfertamm(:,:,k) + &
!!                            nfertnit(:,:,k) + nliquidmanure(:,:,k) + &
!!                            nslurry(:,:,k) + nsolidmanure(:,:,k))
!!      ENDDO
!!      N_fert_total(:,:) = N_fert_total(:,:) * 10000. + ndeposition(:,:)
!!      N_fert_total(:,1) = 0.0
!!!      DO j=2,nvm
!!!        IF ((management_intensity(j) .EQ. 2).AND. (.NOT. is_c4(j))) THEN
!!!          N_fert_total(:,mcut_C3)=N_fert_total(:,j)
!!!          N_fert_total(:,mgraze_C3)=N_fert_total(:,j)
!!!        ENDIF
!!!        IF ((management_intensity(j) .EQ. 2).AND. (is_c4(j))) THEN
!!!          N_fert_total(:,mcut_C4)=N_fert_total(:,j)
!!!          N_fert_total(:,mgraze_C4)=N_fert_total(:,j)
!!!        ENDIF
!!!
!!!      ENDDO
!!      !gmjc new fertilization effect
!!      ! linear
!!      !N_limfert(:,:) = 1.0 + (1.60-1.0)/320 * N_fert_total(:,:) 
!!      ! index
!!      N_limfert(:,:) = 1. + N_effect - N_effect * (0.75 ** (N_fert_total(:,:)/30))
!!
!!      WHERE (N_limfert(:,:) .LT. 1.0) 
!!        N_limfert(:,:) = 1.0
!!      ELSEWHERE (N_limfert(:,:) .GT. 2.0)
!!        N_limfert(:,:) = 1.+N_effect
!!      ENDWHERE
!!
!!  END SUBROUTINE calc_N_limfert

! Author: Xuhui Wang
! Date: Oct. 18th, 2010
! Interpolate (extract) Planting Date information
! for a specific crop type
! Modified by Jinfeng Chang
! Date: Dec. 1st, 2014 
! General management map reading for grassland management module
! Modified by Jinfeng Chang
! Date: Apr. 21st, 2016
! to speed up the running, only open management file once 
! reading all variables
  SUBROUTINE slowproc_GRM_ManageInput(npts,lalo,neighbours,resolution,contfrac, &
               strIn,varname1,manage1,varname2,manage2,varname3,manage3,&
                         varname4,manage4,varname5,manage5,varname6,manage6,yrlen)

!    INTEGER, parameter :: i_std = 4
!    REAL, parameter :: r_std = 8
    !
    ! 0.1 INPUT
    !
    INTEGER, INTENT(in)  :: npts         ! Number of points for which the data needs to be interpolated (extracted)
    REAL, INTENT(in) :: lalo(npts,2)     ! Vector of latitude and longtitude
    INTEGER, INTENT(in)  :: neighbours(npts,8)   ! Vectors of neighbours for each grid point
    ! (1=N, 2=NE, 3=E, 4=SE, 5=S, 6=SW, 7=W, 8=NW)
    REAL,INTENT(in)  :: resolution(npts,2)   ! The size in km of each grid box in lat and lon
    !REAL             :: resolution(npts,2)   ! The size in km of each
    !grid box in lat and lon
    REAL,INTENT(in)  :: contfrac(npts)   ! The fraction of land in each grid box
    CHARACTER(LEN=30),INTENT(in) :: strIn       ! getin parameter and Call Sign of the management data
    CHARACTER(LEN=30),INTENT(in) :: varname1     ! variable name in the nc file
    CHARACTER(LEN=30),INTENT(in) :: varname2     ! variable name in the nc file
    CHARACTER(LEN=30),INTENT(in) :: varname3     ! variable name in the nc file
    CHARACTER(LEN=30),INTENT(in) :: varname4     ! variable name in the nc file
    CHARACTER(LEN=30),INTENT(in) :: varname5     ! variable name in the nc file
    CHARACTER(LEN=30),INTENT(in) :: varname6     ! variable name in the nc file

    !
    ! 0.2 OUTPUT
    !
    REAL,ALLOCATABLE, INTENT(out)    :: manage1(:,:,:)    ! The planting date of the crop: npts, veg, year
    REAL,ALLOCATABLE, INTENT(out)    :: manage2(:,:,:)
    REAL,ALLOCATABLE, INTENT(out)    :: manage3(:,:,:)
    REAL,ALLOCATABLE, INTENT(out)    :: manage4(:,:,:)
    REAL,ALLOCATABLE, INTENT(out)    :: manage5(:,:,:)
    REAL,ALLOCATABLE, INTENT(out)    :: manage6(:,:,:)
    ! nvm is the number of PFTs, there may not be planting date for all the PFTs
    INTEGER, INTENT(out)             :: yrlen            ! year length of the output matrix
    !
    ! 0.3 LOCAL
    !
    INTEGER      :: nbvmax       ! a parameter for interpolation
    REAL         :: myres(npts,2)
    CHARACTER(LEN=80)       :: filename
    INTEGER      :: iml, jml, lml, tml, fid, fid1
    INTEGER      :: ip, jp, ib, ilf, fopt, it ! for-loop variable
    INTEGER      :: nbexp
    REAL         :: lev(1), date, dt
    REAL         :: missing_val
    INTEGER      :: itau(1)

    INTEGER      :: nb_dim
    INTEGER,DIMENSION(flio_max_var_dims) :: l_d_w
    LOGICAL         :: l_ex

    REAL, ALLOCATABLE, DIMENSION(:,:)    :: lat_rel, lon_rel
    REAL, ALLOCATABLE, DIMENSION(:,:,:,:)    :: manage_mat1 ! LON LAT VEGET, Time
    REAL, ALLOCATABLE, DIMENSION(:,:,:,:)    :: manage_mat2 
    REAL, ALLOCATABLE, DIMENSION(:,:,:,:)    :: manage_mat3
    REAL, ALLOCATABLE, DIMENSION(:,:,:,:)    :: manage_mat4
    REAL, ALLOCATABLE, DIMENSION(:,:,:,:)    :: manage_mat5
    REAL, ALLOCATABLE, DIMENSION(:,:,:,:)    :: manage_mat6
! JC for loop variables
    REAL, ALLOCATABLE, DIMENSION(:,:,:,:,:)   :: manage_mat_all
    REAL, ALLOCATABLE, DIMENSION(:,:,:,:)    :: manage_all
    INTEGER      :: nb_var_manag
    INTEGER      :: iv_manag
! end gmjc
    INTEGER, ALLOCATABLE, DIMENSION(:,:) :: mask
    REAL, ALLOCATABLE, DIMENSION(:,:)    :: temp_data
    REAL, ALLOCATABLE, DIMENSION(:,:)    :: sub_area
    INTEGER, ALLOCATABLE, DIMENSION(:,:,:)   :: sub_index

    REAL :: sgn, sum_float
    INTEGER :: ivgt   ! , icyc, pltcyc
    CHARACTER(LEN=30) :: callsign
    LOGICAL :: ok_interpol
    INTEGER :: ALLOC_ERR
    LOGICAL :: mydebug = .true.

!   ! croptype = TRIM(croptype) !if croptype is a string
!   ! else a switch expression is needed
!   filename = "/work/cont003/p529tan/WXH/plt_date_modif.nc" ! default input
!   file
!   ! String operation needed
    filename = "PlantingDate.nc"
    CALL getin_p(strIn,filename)

    IF (is_root_prc) THEN
    ! ? what does is_root_prc mean?
        CALL flininfo(filename, iml, jml, lml, tml, fid)
    ENDIF
    CALL bcast(iml)
    CALL bcast(jml)
    ! CALL bcast(lml)
    ! CALL bcast(tml)

    ! Printing information for debugging
    IF (mydebug) THEN
        WRITE(numout, *) "Xuhui's debug info for slowproc_ManageInput #1:"
        WRITE(numout, *) "string in: ", strIn
        WRITE(numout, *) "variable name: ", varname1,varname2,varname3,&
                            varname4,varname5,varname6
        WRITE(numout, *) "filename is: ", filename
        WRITE(numout, *) "Dimension 1, lon, iml:", iml
        WRITE(numout, *) "Dimension 2, lat, jml:", jml
        WRITE(numout, *) "Dimension 3, veget, lml:", lml
        WRITE(numout, *) "Dimension 4, time, tml:", tml
    ENDIF
    ! apparently, flinget function is not designed to take veget but levels to
    ! be the
    ! 3rd dimension, modification to lml is needed

!JG all flio calls must be done by is_root_prc
    IF (is_root_prc) THEN
       CALL flioopfd(filename,fid1)
       ! JC here only use dimension of the first variable
       CALL flioinqv(fid1,v_n=varname1, l_ex = l_ex, nb_dims = nb_dim, len_dims =l_d_w)
       IF (lml == 0) THEN
          ! CALL
          ! flioinqv(fid1,v_n="PLNTDT",l_ex=l_ex,nb_dims=nb_dim,len_dims=l_d_w)
          lml=l_d_w(3)
          IF (mydebug) THEN
              WRITE(numout, *) "len_dims: ", l_d_w
              WRITE(numout, *) "lml AFTER revision"
              WRITE(numout, *) "lml: ", lml
          ENDIF
       ENDIF
       IF (mydebug) THEN
           WRITE(numout,*) "nb_dim: ", nb_dim
           WRITE(numout,*) "resolution: ", resolution(1,:)
       ENDIF

       IF (nb_dim .NE. 4) THEN
          WRITE(numout,*) "dimension not supported for ", nb_dim
       ENDIF
       tml = l_d_w(4)
       !yrlen = tml
    END IF
    IF (mydebug) THEN
        WRITE(numout, *) "Now the tml is, :", tml
        WRITE(numout, *) "Now the lml is:", lml
    ENDIF

!JG REMVOVE    CALL flioclo(fid1)
    CALL bcast(lml)
    CALL bcast(tml)
    CALL bcast(nb_dim)
    ! CALL bcast(plantcyc)

    ! JG yrlen must not be done after bcast(tml)
    yrlen = tml
    nb_var_manag = INT(6)
    
    ALLOC_ERR=-1
    ALLOCATE(manage_all(nb_var_manag,npts,lml,tml),STAT=ALLOC_ERR)
    IF (ALLOC_ERR/=0) THEN
        WRITE(numout,*) "ERROR IN ALLOCATION OF manage: ", ALLOC_ERR
    ENDIF
    WRITE(numout,*) "manage ALLOCATED"
    !CALL bcast(manage)
    ALLOC_ERR=-1
    ALLOCATE(manage1(npts,lml,tml),STAT=ALLOC_ERR)
    IF (ALLOC_ERR/=0) THEN
        WRITE(numout,*) "ERROR IN ALLOCATION OF manage: ", ALLOC_ERR
    ENDIF
    WRITE(numout,*) "manage ALLOCATED"
    !CALL bcast(manage)
    ALLOC_ERR=-1
    ALLOCATE(manage2(npts,lml,tml),STAT=ALLOC_ERR)
    IF (ALLOC_ERR/=0) THEN
        WRITE(numout,*) "ERROR IN ALLOCATION OF manage: ", ALLOC_ERR
    ENDIF
    WRITE(numout,*) "manage ALLOCATED"
    !CALL bcast(manage)
    ALLOC_ERR=-1
    ALLOCATE(manage3(npts,lml,tml),STAT=ALLOC_ERR)
    IF (ALLOC_ERR/=0) THEN
        WRITE(numout,*) "ERROR IN ALLOCATION OF manage: ", ALLOC_ERR
    ENDIF
    WRITE(numout,*) "manage ALLOCATED"
    !CALL bcast(manage)
    ALLOC_ERR=-1
    ALLOCATE(manage4(npts,lml,tml),STAT=ALLOC_ERR)
    IF (ALLOC_ERR/=0) THEN
        WRITE(numout,*) "ERROR IN ALLOCATION OF manage: ", ALLOC_ERR
    ENDIF
    WRITE(numout,*) "manage ALLOCATED"
    !CALL bcast(manage)
    ALLOC_ERR=-1
    ALLOCATE(manage5(npts,lml,tml),STAT=ALLOC_ERR)
    IF (ALLOC_ERR/=0) THEN
        WRITE(numout,*) "ERROR IN ALLOCATION OF manage: ", ALLOC_ERR
    ENDIF
    WRITE(numout,*) "manage ALLOCATED"
    !CALL bcast(manage)
    ALLOC_ERR=-1
    ALLOCATE(manage6(npts,lml,tml),STAT=ALLOC_ERR)
    IF (ALLOC_ERR/=0) THEN
        WRITE(numout,*) "ERROR IN ALLOCATION OF manage: ", ALLOC_ERR
    ENDIF
    WRITE(numout,*) "manage ALLOCATED"
    !CALL bcast(manage)
    !
    ALLOC_ERR=-1
    ALLOCATE(lat_rel(iml,jml), STAT=ALLOC_ERR)
      IF (ALLOC_ERR/=0) THEN
        WRITE(numout,*) "ERROR IN ALLOCATION of lat_rel : ",ALLOC_ERR
        STOP
    ENDIF
!    CALL bcast(lat_rel)

    ALLOC_ERR=-1
    ALLOCATE(lon_rel(iml,jml), STAT=ALLOC_ERR)
   IF (ALLOC_ERR/=0) THEN
        WRITE(numout,*) "ERROR IN ALLOCATION of lon_rel : ",ALLOC_ERR
        STOP
    ENDIF
!    CALL bcast(lon_rel)

    ALLOC_ERR=-1
    ALLOCATE(mask(iml,jml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR/=0) THEN
        WRITE(numout,*) "ERROR IN ALLOCATION of mask : ",ALLOC_ERR
        STOP
    ENDIF
!    CALL bcast(mask)


    ALLOC_ERR=-1
    ALLOCATE(manage_mat_all(nb_var_manag,iml,jml,lml,tml), STAT=ALLOC_ERR)
    ! !lml is supposed to be nvm (number of PFTs), if not ,change it
    IF (ALLOC_ERR/=0) THEN
        WRITE(numout,*) "ERROR IN ALLOCATION of manage_mat : ",ALLOC_ERR
        STOP
    ENDIF

    ALLOC_ERR=-1
    ALLOCATE(manage_mat1(iml,jml,lml,tml), STAT=ALLOC_ERR)
    ! !lml is supposed to be nvm (number of PFTs), if not ,change it
    IF (ALLOC_ERR/=0) THEN
        WRITE(numout,*) "ERROR IN ALLOCATION of manage_mat : ",ALLOC_ERR
        STOP
    ENDIF
!    CALL bcast(manage_mat)
    ALLOC_ERR=-1
    ALLOCATE(manage_mat2(iml,jml,lml,tml), STAT=ALLOC_ERR)
    ! !lml is supposed to be nvm (number of PFTs), if not ,change it
    IF (ALLOC_ERR/=0) THEN
        WRITE(numout,*) "ERROR IN ALLOCATION of manage_mat : ",ALLOC_ERR
        STOP
    ENDIF
!    CALL bcast(manage_mat)
    ALLOC_ERR=-1
    ALLOCATE(manage_mat3(iml,jml,lml,tml), STAT=ALLOC_ERR)
    ! !lml is supposed to be nvm (number of PFTs), if not ,change it
    IF (ALLOC_ERR/=0) THEN
        WRITE(numout,*) "ERROR IN ALLOCATION of manage_mat : ",ALLOC_ERR
        STOP
    ENDIF
!    CALL bcast(manage_mat)
    ALLOC_ERR=-1
    ALLOCATE(manage_mat4(iml,jml,lml,tml), STAT=ALLOC_ERR)
    ! !lml is supposed to be nvm (number of PFTs), if not ,change it
    IF (ALLOC_ERR/=0) THEN
        WRITE(numout,*) "ERROR IN ALLOCATION of manage_mat : ",ALLOC_ERR
        STOP
    ENDIF
!    CALL bcast(manage_mat)
    ALLOC_ERR=-1
    ALLOCATE(manage_mat5(iml,jml,lml,tml), STAT=ALLOC_ERR)
    ! !lml is supposed to be nvm (number of PFTs), if not ,change it
    IF (ALLOC_ERR/=0) THEN
        WRITE(numout,*) "ERROR IN ALLOCATION of manage_mat : ",ALLOC_ERR
        STOP
    ENDIF
!    CALL bcast(manage_mat)
    ALLOC_ERR=-1
    ALLOCATE(manage_mat6(iml,jml,lml,tml), STAT=ALLOC_ERR)
    ! !lml is supposed to be nvm (number of PFTs), if not ,change it
    IF (ALLOC_ERR/=0) THEN
        WRITE(numout,*) "ERROR IN ALLOCATION of manage_mat : ",ALLOC_ERR
        STOP
    ENDIF
!    CALL bcast(manage_mat)
!    WRITE (numout,*) 'bcast manage_mat'

    ! input of some attributes
    IF (is_root_prc) THEN
! JG with the flioclo, done before this was not ok. Now ok
        CALL flinget(fid, 'LON', iml, jml, lml, tml, 1, 1, lon_rel)
        CALL flinget(fid, 'LAT', iml, jml, lml, tml, 1, 1, lat_rel)
    ENDIF
    CALL bcast(lon_rel)
    CALL bcast(lat_rel)
    WRITE (numout,*) 'lon_rel size: ', SIZE(lon_rel)
    WRITE (numout,*) 'lat_rel size: ', SIZE(lat_rel)


    ! input of the matrix
    IF (is_root_prc) THEN
        ! CALL flinget(fid, 'PLNTDT', iml, jml, lml, tml, 1, 1, plntdt_mat)
! JG remove CALL flioopfd: already done
!       CALL flioopfd(filename,fid1)
        CALL fliogetv(fid1,trim(varname1),manage_mat1,start=(/1,1,1,1/),count=(/iml,jml,lml,tml/))
        CALL fliogetv(fid1,trim(varname2),manage_mat2,start=(/1,1,1,1/),count=(/iml,jml,lml,tml/))
        CALL fliogetv(fid1,trim(varname3),manage_mat3,start=(/1,1,1,1/),count=(/iml,jml,lml,tml/))
        CALL fliogetv(fid1,trim(varname4),manage_mat4,start=(/1,1,1,1/),count=(/iml,jml,lml,tml/))
        CALL fliogetv(fid1,trim(varname5),manage_mat5,start=(/1,1,1,1/),count=(/iml,jml,lml,tml/))
        CALL fliogetv(fid1,trim(varname6),manage_mat6,start=(/1,1,1,1/),count=(/iml,jml,lml,tml/))
        ! get missing_val
        CALL fliogeta(fid1,varname1,'missing_value',missing_val)

        CALL flioclo(fid1)
    ENDIF
    CALL bcast(manage_mat1)
    CALL bcast(manage_mat2)
    CALL bcast(manage_mat3)
    CALL bcast(manage_mat4)
    CALL bcast(manage_mat5)
    CALL bcast(manage_mat6)
    WRITE (numout,*) 'bcast manage_mat'

    ! JC combine matrix for loop
    manage_mat_all(1,:,:,:,:) = manage_mat1
    manage_mat_all(2,:,:,:,:) = manage_mat2
    manage_mat_all(3,:,:,:,:) = manage_mat3
    manage_mat_all(4,:,:,:,:) = manage_mat4
    manage_mat_all(5,:,:,:,:) = manage_mat5
    manage_mat_all(6,:,:,:,:) = manage_mat6   

    ! WRITE(numout,*) 'manage_mat size: ',SIZE(manage_mat)
    ! WRITE(numout,*) 'missing value: ', missing_val
    ! WRITE(numout,*) 'lat(361,284): ',lat_rel(361,284)
    ! WRITE(numout,*) 'lon(361,284): ',lon_rel(361,284)
    ! WRITE(numout,*) 'plntdt(361,284,1,1): ',plntdt_mat(361,284,1,1)

    IF (is_root_prc) CALL flinclo(fid)

    manage1(:,:,:) = zero ! npts veget year
    manage2(:,:,:) = zero ! npts veget year
    manage3(:,:,:) = zero ! npts veget year
    manage4(:,:,:) = zero ! npts veget year
    manage5(:,:,:) = zero ! npts veget year
    manage6(:,:,:) = zero ! npts veget year
    manage_all(:,:,:,:) = zero

DO iv_manag = 1,nb_var_manag
    DO it = 1,tml
        DO ivgt = 1,lml ! ? We can suppose PFTs less than 10 are natural veg without planting date, but not now
!            IF (.NOT. natural(ivgt)) THEN
                WRITE(numout,*) "variable, veget, time: ",iv_manag, ivgt,it
                nbexp = 0
                ! the number of exceptions

!JCCOMMENT GRM_input.nc for every grid value >=0
! thus mask = un
                ! mask of available value
!                mask(:,:) = zero;  ! Defined in constante.f90
             IF (iv_manag .EQ. 1) THEN
                mask(:,:) = un
!                DO ip = 1,iml
!                    DO jp = 1,jml
!                        IF ((manage_mat_all(iv_manag,ip,jp,ivgt,it) .GT. min_sechiba) .AND. &
!                        (manage_mat_all(iv_manag,ip,jp,ivgt,it) /= missing_val)) THEN
!                            mask(ip,jp) = un;  ! Defined in constante.f90
!                            ! here we assumed that for each plant cycle at each
!                            ! there might be missing data at different grid
!                            ! in this case, mask has to be generated each plant
!                            ! cycle each time step
!                        ENDIF
!                    ENDDO
!                ENDDO

                ! Interpolation started
                nbvmax = 200
                ! the maximum amount of fine grids that one coarse grid may have

                callsign = strIn

                ok_interpol = .FALSE.

                DO WHILE ( .NOT. ok_interpol )
                    WRITE(numout,*) "Pojection arrays for ", callsign, ":"
                    WRITE(numout,*) "nbvmax = ", nbvmax

                    ALLOC_ERR = -1
                    ALLOCATE(temp_data(nbvmax,lml), STAT=ALLOC_ERR)
                    IF (ALLOC_ERR /=0) THEN
                        WRITE(numout,*) "ERROR IN ALLOCATION OF temp_data :", ALLOC_ERR
                        STOP
                    ENDIF
                    ALLOC_ERR = -1
                    ALLOCATE(sub_index(npts,nbvmax,2), STAT=ALLOC_ERR)
                    IF (ALLOC_ERR /=0) THEN
                        WRITE(numout,*) "ERROR IN ALLOCATION OF sub_index :", ALLOC_ERR
                        STOP
                    ENDIF
                    sub_index(:,:,:) = zero
                    ALLOC_ERR = -1
                    ALLOCATE(sub_area(npts, nbvmax), STAT=ALLOC_ERR)
                    IF (ALLOC_ERR /=0) THEN
                        WRITE(numout,*) "ERROR IN ALLOCATION OF sub_area :",ALLOC_ERR
                        STOP
                    ENDIF
                    sub_area(:,:) = zero
                    myres(:,:) = resolution(:,:)/1000  !m -> km
                    !write(numout,*) "resolution updated: ", myres(1,:), " km"
                    !CALL bcast(myres)
!                    CALL bcast(myres)

                    write(*,*) "calling aggregate_p? "
                   CALL aggregate_p(npts, lalo, neighbours, myres, contfrac, &
                    &                iml, jml, lon_rel, lat_rel, mask, callsign, &
                    &                nbvmax, sub_index, sub_area, ok_interpol)
                    write(numout,*) "wu: we finished aggregate_p:) "

                    IF ( .NOT. ok_interpol ) THEN
                        DEALLOCATE(temp_data)
                        DEALLOCATE(sub_index)
                        DEALLOCATE(sub_area)
                        nbvmax = nbvmax * 2
                    ENDIF
                ENDDO

                WRITE(numout,*) "called aggregate_p"
             ENDIF ! only call aggregate once
                ! assign the values to plantdate
                ! values should be given to all PFTs
                DO ib = 1, npts
                    ! examing all sub_point we found
                    fopt = COUNT(sub_area(ib,:)>zero)

                    ! confirm that we found some points
                    IF ( fopt .EQ. 0) THEN
                        nbexp = nbexp + 1
                        manage_all(iv_manag,ib,ivgt,it) = val_exp
                    ELSE
                        DO ilf = 1,fopt
                            ! !Not to get lat and lon in wrong order
                            temp_data(ilf,ivgt) = manage_mat_all(iv_manag,sub_index(ib,ilf,1),sub_index(ib,ilf,2),ivgt,it)
                        ENDDO

                        sgn = zero
                        sum_float = zero
                        DO ilf = 1,fopt
                            ! average the data weighted by area ! better to
                            ! multiply
                            ! PFT HERE
                            ! need to add management specific judgem
                                sum_float = sum_float + temp_data(ilf,ivgt)*sub_area(ib,ilf)
                                sgn = sgn + sub_area(ib,ilf)
                        ENDDO

                        ! Normalize the surface
                        ! sgn can be a scaler, however, to prepare it for future
                        ! incorporation of fraction
                        ! I make it a vector with nvm values which are equal to
                        ! each
                        ! other
                        IF ( sgn .LT. min_sechiba) THEN
                            nbexp = nbexp + 1
                            manage_all(iv_manag,ib,ivgt,it) = val_exp ! plantdate_default(ivgt)
                        ELSE
                        ! ANINT is used for plant date integer
                        ! BUT not for grassland management input
                            !manage(ib,ivgt,it) = ANINT(sum_float/sgn)
                            manage_all(iv_manag,ib,ivgt,it) = sum_float/sgn
                        ENDIF

                    ENDIF

                ENDDO ! ib
                WRITE(numout,*) 'fopt subarea',fopt!,sub_area

                IF ( nbexp .GT. 0) THEN
                    WRITE(numout,*) 'slowproc_ManageInput : exp_val was applied in', nbexp, 'grid(s)'
                    WRITE(numout,*) 'slowproc_ManageInput : These are either coastal points or having missing data'
                ENDIF
!JC keep the sub_area sub_index till all variables read
!                DEALLOCATE (sub_area)
!                DEALLOCATE (sub_index)
!                DEALLOCATE (temp_data)
                ! WRITE(numout,*) 'Planting Date of Site 1 veget ',ivgt,' :
                ! ',plantdate(1,ivgt,icyc)
!            ENDIF
        ENDDO
        ! End of Veget cycle
    ENDDO
    ! End of Time Axis cycle
ENDDO
! End of variables
! gmjc 
                DEALLOCATE (sub_area)
                DEALLOCATE (sub_index)
                DEALLOCATE (temp_data)

    DEALLOCATE (lat_rel)
    DEALLOCATE (lon_rel)
    DEALLOCATE (mask)
! gmjc multivariables
    manage1(:,:,:) = manage_all(1,:,:,:) 
    manage2(:,:,:) = manage_all(2,:,:,:)
    manage3(:,:,:) = manage_all(3,:,:,:)
    manage4(:,:,:) = manage_all(4,:,:,:)
    manage5(:,:,:) = manage_all(5,:,:,:)
    manage6(:,:,:) = manage_all(6,:,:,:)
    DEALLOCATE (manage_mat1)
    DEALLOCATE (manage_mat2)
    DEALLOCATE (manage_mat3)
    DEALLOCATE (manage_mat4)
    DEALLOCATE (manage_mat5)
    DEALLOCATE (manage_mat6)
    DEALLOCATE (manage_mat_all)
    DEALLOCATE (manage_all)

    WRITE (numout,*) 'Output Management Date:'
    WRITE (numout,*) 'time_step 1:'
    WRITE (numout,*) manage1(1,:,1), manage2(1,:,1), manage3(1,:,1), &
                     manage4(1,:,1), manage5(1,:,1), manage6(1,:,1)
    IF (tml>1) THEN
        WRITE (numout,*) 'time_step 2:'
        WRITE (numout,*) manage1(1,:,2)
    ENDIF
    WRITE (numout,*) '***END of DEBUG INFO slowproc_ManageInput***'
    RETURN

  END SUBROUTINE slowproc_GRM_ManageInput
! End of Edition by Xuhui, Mar. 16th 2011

END MODULE grassland_management
