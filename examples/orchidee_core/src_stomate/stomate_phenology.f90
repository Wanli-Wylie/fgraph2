! =================================================================================================================================
! MODULE        : stomate_phenology
!
! CONTACT       : orchidee-help _at_ listes.ipsl.fr
!
! LICENCE       : IPSL (2006). This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF        This module manages the beginning of the growing season (leaf onset).
!!      
!!\n DESCRIPTION: None
!!
!! RECENT CHANGE(S): None
!!
!! SVN          :
!! $HeadURL: svn://forge.ipsl.jussieu.fr/orchidee/branches/ORCHIDEE-CN-P/ORCHIDEE/src_stomate/stomate_phenology.f90 $ 
!! $Date: 2018-09-14 22:28:33 +0200 (五, 2018-09-14) $
!! $Revision: 5411 $
!! \n
!_ =================================================================================================================================

MODULE stomate_phenology

  ! modules used:
  USE xios_orchidee
  USE ioipsl_para
  USE stomate_data
  USE constantes
  USE pft_parameters
  USE function_library,    ONLY: calculate_c0_alloc, wood_to_ba

  IMPLICIT NONE

  ! private & public routines

  PRIVATE
  PUBLIC phenology,phenology_clear

  ! first call
  LOGICAL, SAVE                                              :: firstcall_all_phenology = .TRUE.
!$OMP THREADPRIVATE(firstcall_all_phenology)
  LOGICAL, SAVE                                              :: firstcall_hum = .TRUE.
!$OMP THREADPRIVATE(firstcall_hum)
  LOGICAL, SAVE                                              :: firstcall_moi = .TRUE.
!$OMP THREADPRIVATE(firstcall_moi)
  LOGICAL, SAVE                                              :: firstcall_humgdd = .TRUE.
!$OMP THREADPRIVATE(firstcall_humgdd)
  LOGICAL, SAVE                                              :: firstcall_moigdd = .TRUE.
!$OMP THREADPRIVATE(firstcall_moigdd)

CONTAINS


!! ================================================================================================================================
!! SUBROUTINE   : phenology_clear
!!
!>\BRIEF          Flags setting   
!!
!! DESCRIPTION  : This subroutine sets flags 
!!                ::firstcall_all_phenology, ::firstcall_hum, ::firstcall_moi, ::firstcall_humgdd, 
!!                ::firstcall_moigdd to .TRUE., and therefore activates section 1.1 of each 
!!                subroutine which writes messages to the output. \n
!!                This subroutine is called at the beginning of ::stomateLpj_clear in the 
!!                ::stomate_lpj module.
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): ::firstcall_all_phenology, ::firstcall_hum, ::firstcall_moi, ::firstcall_humgdd, 
!!                ::firstcall_moigdd
!!
!! REFERENCE(S)  : None
!!
!! FLOWCHART     : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE phenology_clear
    firstcall_all_phenology=.TRUE.
    firstcall_hum=.TRUE.
    firstcall_moi = .TRUE.
    firstcall_humgdd = .TRUE.
    firstcall_moigdd = .TRUE.
  END SUBROUTINE phenology_clear


!! ================================================================================================================================
!! SUBROUTINE   : phenology
!!
!>\BRIEF          This subroutine controls the detection of the beginning of the growing season 
!!                (if dormancy has been long enough), leaf onset, given favourable biometeorological 
!!                conditions, and leaf growth and biomass allocation when leaf biomass is low (i.e. 
!!                at the start of the growing season.
!!
!! DESCRIPTION  : This subroutine is called by the module ::stomate_lpj and deals with the beginning of the  
!!                growing season. First it is established whether the beginning of the growing season is
!!                allowed. This occurs if the dormance period has been long enough (i.e. greater 
!!                than a minimum PFT-dependent threshold, specified by ::lowgpp_time), 
!!                AND if the last beginning of the growing season was a sufficiently long time ago 
!!                (i.e. when the growing season length is greater than a minimum threshold, specified
!!                by ::min_growthinit_time, which is defined in this module to be 300 days. \n
!!                The dormancy time-length is represented by the variable 
!!                ::time_lowgpp, which is calculated in ::stomate_season. It is increased by 
!!                the stomate time step when the weekly GPP is lower than a threshold. Otherwise
!!                it is set to zero. \n
!!                ::lowgpp_time is set for each PFT in ::stomate_data from a table of all
!!                PFT values (::lowgpp_time_tab), which is defined in ::stomate_constants. \n
!!                The growing season length is given by ::when_growthinit, which increases
!!                by the stomate time-step at each call to this phenology module, except for when
!!                leaf onset is detected, when it is set to 0. \n
!!                If these two conditions are met, leaf onset occurs if the biometeorological 
!!                conditions are also met. This is determined by the leaf onset models, which are
!!                biome-specific. Each PFT is looped over (ignoring bare soil).
!!                The onset phenology model is selected, (according to the parameter 
!!                ::pheno_model, which is initialised in stomate_data), and called. \n
!!                There are six leaf onset phenology models currently being used by ORCHIDEE. 
!!                These are: 'hum' and 'moi', which are based exclusively on moisture conditions,
!!                'humgdd' and 'moigdd', which are based on both temperature and moisture conditions,
!!                'ncdgdd', which is based on a "chilling" requirement for leaf onset, and 
!!                'ngd', which is based on the number of growing days since the temperature was 
!!                above a certain threshold, to account for the end of soil frost.
!!                Those models which are based mostly on temperature conditions are used for
!!                temperate and boreal biomes, and those which include a moisture condition are used
!!                for tropical biomes. More detail on the biometeorological conditions is provided
!!                in the sections on the individual onset models. \n
!!                The moisture conditions are based on the concept of plant "moisture availability".
!!                This is based on the soil humidity (relative soil moisture), but is moderated by
!!                the root density profile, as per the equation:
!!                \latexonly
!!                \input{phenology_moiavail_eqn1.tex}
!!                \endlatexonly
!!                \n
!!                Although some studies have shown that the length of the photoperiod is important
!!                in determining onset (and senescence) dates, this is not considered in the current
!!                versions of the onset models (Krinner et al., 2005). \n
!!                If conditions are favourable, leaf onset occurs (::begin_leaves is set to TRUE), 
!!                ::when_growthinit is set to 0.0, and the growing season has begun. \n
!!                Following the detection of leaf onset, biomass is allocated from the carbohydrate 
!!                reserves equally to the leaves and roots IF the leaf biomass is lower than a minimum
!!                threshold, which is calculated in this subroutine from the parameter
!!                ::lai_initmin, divided by the specific leaf area (both of which are
!!                PFT-dependent and set in ::stomate_constants). \n
!!                Finally, if biomass is required to be allocated from the carbohydrate reserve 
!!                because the leaf biomass is too low, the leaf age and leaf age distribution is 
!!                re-set. In this case the youngest age class fraction is set to 1 and all other   
!!                leaf age class fractions are set to 0. All leaf ages are set to 0. If there is 
!!                no biomass in the carbohydrate reserve, leaf onset will not occur and the PFT
!!                will disappear from the grid cell (Krinner et al., 2005). \n
!!                This subrouting is called in ::stomate_lpj.
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): ::biomass, 
!!                        ::when_growthinit,
!!                        ::leaf age distribution
!!                        ::leaf fraction
!!
!! REFERENCE(S) :
!! - Krinner, G., N. Viovy, N. de Noblet-Ducoudre, J. Ogee, J. Polcher, P. 
!! Friedlingstein, P. Ciais, S. Sitch and I.C. Prentice (2005), A dynamic global
!! vegetation model for studies of the coupled atmosphere-biosphere system, Global
!! Biogeochemical Cycles, 19, doi:10.1029/2003GB002199.
!!
!! FLOWCHART    : 
!! \latexonly
!! \includegraphics[scale = 1]{phenology_flowchart.png}
!! \endlatexonly
!! \n
!_ ================================================================================================================================

  SUBROUTINE phenology (npts, dt, PFTpresent, &
       veget_cov_max, &
       t2m_longterm, t2m_month, t2m_week, gpp, &
       maxmoiavail_lastyear, minmoiavail_lastyear, &
       moiavail_month, moiavail_week, &
       gdd_m5_dormance, gdd_midwinter, ncd_dormance, ngd_minus5, &
       senescence, time_hum_min, &
       biomass, som, leaf_frac, leaf_age, &
       when_growthinit, co2_to_bm, n_to_bm, p_to_bm,&
       begin_leaves, ind, KF, N_support, P_support, &
       cn_leaf_avg_season, np_leaf_avg_season, &
       sla_calc)

    !
    !! 0. Variable and parameter declaration
    !

    !
    !! 0.1 Input variables
    !
    INTEGER, INTENT(in)                                          :: npts                 !! Domain size - number of grid 
                                                                                                !! cells (unitless) 
    REAL, INTENT(in)                                             :: dt                   !! time step (dt_days)
    LOGICAL, DIMENSION(npts,nvm), INTENT(in)                            :: PFTpresent           !! PFT exists (true/false)
    REAL, DIMENSION(npts,nvm), INTENT(in)                        :: veget_cov_max        !! "maximal" coverage fraction of a 
                                                                                                !! PFT (LAI -> infinity) on ground 
                                                                                                !! (0-1, unitless)
    REAL, DIMENSION(npts), INTENT(in)                            :: t2m_longterm         !! "long term" 2 meter reference 
                                                                                                !! temperatures (K) 
    REAL, DIMENSION(npts), INTENT(in)                            :: t2m_month            !! "monthly" 2-meter temperatures 
                                                                                                !! (K) 
    REAL, DIMENSION(npts), INTENT(in)                            :: t2m_week             !! "weekly" 2-meter temperatures (K)
    REAL, DIMENSION(npts,nvm), INTENT(in)                        :: gpp                  !! daily gross primary productivity 
                                                                                                !! @tex ($gC m^{-2} of 
                                                                                                !! ground/day$) @endtex 
    REAL, DIMENSION(npts,nvm), INTENT(in)                        :: maxmoiavail_lastyear !! last year's maximum moisture 
                                                                                                !! availability (0-1, unitless)
    REAL, DIMENSION(npts,nvm), INTENT(in)                        :: minmoiavail_lastyear !! last year's minimum moisture 
                                                                                                !! availability (0-1, unitless)
    REAL, DIMENSION(npts,nvm), INTENT(in)                        :: moiavail_month       !! "monthly" moisture availability 
                                                                                                !! (0-1, unitless)
    REAL, DIMENSION(npts,nvm), INTENT(in)                        :: moiavail_week        !! "weekly" moisture availability 
                                                                                                !! (0-1, unitless)
    REAL, DIMENSION(npts,nvm), INTENT(in)                        :: gdd_m5_dormance      !! growing degree days above a 
                                                                                                !! threshold of -5 deg C (C) 
    REAL, DIMENSION(npts,nvm), INTENT(in)                        :: ncd_dormance         !! number of chilling days since 
                                                                                                !! leaves were lost (days) 
    REAL, DIMENSION(npts,nvm), INTENT(in)                        :: ngd_minus5           !! number of growing days above a 
                                                                                                !! threshold of -5 deg C (days) 
    LOGICAL, DIMENSION(npts,nvm), INTENT(in)                            :: senescence           !! is the plant senescent? (only 
                                                                                                !! for deciduous trees - 
                                                                                                !! carbohydrate reserve) 
                                                                                                !! (true/false) 
    REAL, DIMENSION(npts,nvm), INTENT(in)                        :: time_hum_min         !! time elapsed since strongest 
                                                                                                !! moisture availability (days) 
    REAL, DIMENSION(npts,nvm), INTENT(in)                        :: ind                  !! Stand level number of individuals
                                                                                                !! @tex $(ind m^{-2})$ @endtex 
    REAL, DIMENSION(npts,nvm), INTENT(in)                        :: KF                   !! Scaling factor to convert sapwood mass
                                                                                                !! into leaf mass (m)
    REAL, DIMENSION(npts,nvm), INTENT(in)                        :: sla_calc             !! leaf age-related SLA
    !
    !! 0.2 Ouput variables 
    !

    !
    !! 0.3 Modified variables
    !
    REAL, DIMENSION(npts,nvm,nparts,nelements), INTENT(inout)    :: biomass              !! biomass @tex ($gC m^{-2} of 
                                                                                                !! ground$) @endtex
    REAL,DIMENSION(npts,ncarb,nvm,nelements), INTENT(inout)      :: som   !! SOM pool: active, slow, or passive
  !! @tex ($g(C or N) m^{-2}$) @endtex
    REAL, DIMENSION(npts,nvm,nparts,nleafages), INTENT(inout)    :: leaf_frac            !! fraction of leaves in leaf age 
                                                                                                !! class (0-1, unitless)
    REAL, DIMENSION(npts,nvm,nparts,nleafages), INTENT(inout)    :: leaf_age             !! leaf age (days)
    REAL, DIMENSION(npts,nvm), INTENT(inout)                     :: when_growthinit      !! how many days since the 
                                                                                                !! beginning of the growing season 
                                                                                                !! (days) 
    REAL, DIMENSION(npts,nvm), INTENT(inout)                     :: co2_to_bm            !! co2 taken up by carbohydrate 
                                                                                                !! reserve at the beginning of the 
                                                                                                !! growing season @tex ($gC m^{-2} 
                                                                                                !! of total ground/day$) @endtex 
                                                                                                ! NV passge 2D
    REAL, DIMENSION(npts,nvm), INTENT(inout)                     :: n_to_bm              !! N taken up from ?? when 
                                                                                                !! introducing a new PFT (introduced for 
                                                                                                !! nitrogen balance closure) 
                                                                                                !! @tex $(gN m^{-2} dtslow^{-1})$ @endtex 
    REAL, DIMENSION(npts,nvm), INTENT(inout)                     :: p_to_bm              !! P taken up from ?? when 
                                                                                                !! introducing a new PFT (introduced for 
                                                                                                !! nitrogen balance closure) 
                                                                                                !! @tex $(gP m^{-2} dtslow^{-1})$ @endtex 
    REAL, DIMENSION(npts,nvm), INTENT(inout)                     :: gdd_midwinter        !! growing degree days, since midwinter
    REAL, DIMENSION(:,:), INTENT(inout)                          :: cn_leaf_avg_season   !! Average leaf nitrogen concentration (C:N) of the growing season
    REAL, DIMENSION(:,:), INTENT(inout)                          :: np_leaf_avg_season   !! Average leaf phosphorus concentration (N:P) of the growing season
    REAL, DIMENSION(npts,nvm), INTENT(inout)                     :: N_support            !! Nitrogen which is added to the ecosystem to support vegetation growth (only used when IMPOSE_CN .EQ. true) (gN/m2/day)
    REAL, DIMENSION(npts,nvm), INTENT(inout)                     :: P_support            !! Phosphorus which is added to the ecosystem to support vegetation growth (only used when IMPOSE_CN .EQ. true) (gP/m2/day)
    LOGICAL, DIMENSION(npts,nvm), INTENT(out)                           :: begin_leaves         !! signal to start putting leaves on (true/false) 

    !
    !! 0.4 Local variables
    !
    LOGICAL, DIMENSION(npts,nvm)                                        :: allow_initpheno      !! are we allowed to decalre the 
                                                                                                !! beginning of the growing 
                                                                                                !! season? (true/false) 
    REAL                                                         :: bm_wanted            !! biomass we would like to have 
                                                                                                !! @tex ($gC m^{-2} of ground$) 
                                                                                                !! @endtex 
    REAL                                                         :: bm_use               !! biomass we use (from 
                                                                                                !! carbohydrate reserve or from 
                                                                                                !! atmosphere) @tex ($gC m^{-2} of 
                                                                                                !! ground$) @endtex

    REAL                                                         :: bm_wanted_n          !! biomass we would like to have 
                                                                                                !! @tex ($gN m^{-2} of ground$) 
                                                                                                !! @endtex 
    REAL                                                         :: n_use                !! nitrogen biomass we use (from  
                                                                                                !! labile or reserve pool @tex ($gN m^{-2} of  
                                                                                                !! ground$) @endtex 
    REAL                                                         :: bm_wanted_p          !! biomass we would like to have 
                                                                                                !! @tex ($gP m^{-2} of ground$) 
                                                                                                !! @endtex 
    REAL                                                         :: p_use                !! phosphorus biomass we use (from  
                                                                                                !! labile or reserve pool @tex ($gP m^{-2} of  
                                                                                                !! ground$) @endtex 
    REAL                                                         :: n_avail              !! nitrogen in
                                                                                                !! labile or reserve pool @tex ($gN m^{-2} of  
                                                                                                !! ground$) @endtex 
    REAL                                                         :: p_avail              !! phosphorus in
                                                                                                !! labile or reserve pool @tex ($gP m^{-2} of  
                                                                                                !! ground$) @endtex 
    REAL                                                         :: deficit              !! Carbon that needs to be taken from the
                                                                                                !! labile and or reserve pools
                                                                                                !! @tex $(gC m^{-2})$ @endtex     

    REAL                                                        :: circ_class_ba         !! basal area of the model tree in each
                                                                                                !! circ class @tex $(m^{2} m^{-2})$ @endtex 
    REAL                                                        :: Cl_tree               !! Individual plant, leaf compartment
                                                                                                !! @tex $(gC tree^{-1})$ @endtex
    REAL                                                        :: Cr_tree               !! Individual plant, root compartment
                                                                                                !! @tex $(gC tree^{-1})$ @endtex
    REAL                                                        :: Cs_tree               !! Individual plant, sapwood compartment
                                                                                                !! @tex $(gC. tree^{-1})$ @endtex
    REAL                                                        :: Cs_grass             !! Individual plant, sapwood compartment
                                                                                               !! @tex $(gC. ind^{-1})$ @endtex
    REAL                                                        :: Cl_init              !! Initial leaf carbon required to start
                                                                                               !! the growing season
                                                                                               !! @tex $(gC tree^{-1})$ @endtex
    REAL                                                        :: Cr_init              !! Initial root carbon required to start
                                                                                               !! the growing season
                                                                                               !! @tex $(gC tree^{-1})$ @endtex
    REAL                                                        :: height               !! Tree height calculated from allometric 
                                                                                               !! relationships (m)
    REAL                                                        :: lm_min               !! minimum leaf mass @tex ($gC 
                                                                                               !! m^{-2} of ground$) @endtex 
    LOGICAL(r_std), DIMENSION(npts)                                    :: age_reset            !! does the leaf age distribution 
                                                                                               !! have to be reset? (true/false) 
    INTEGER                                                     :: i,j,m,ipar           !! indices (unitless)
    REAL, DIMENSION(npts,nvm)                                   :: histvar              !! controls the history output 
                                                                                               !! level - 0: nothing is written; 
                                                                                               !! 10: everything is written 
                                                                                               !! (0-10, unitless) 

    REAL                                                        :: c0_alloc             !! Root to sapwood tradeoff parameter
    REAL                                                        :: LF                   !! Scaling factor to convert sapwood mass
                                                                                               !! into root mass (unitless)
    REAL                                                        :: cn_leaf_use          !! CN ratio used for allocationg biomass (gC/gN)
    REAL                                                        :: cp_leaf_use          !! CP ratio used for allocationg biomass (gC/gP)
    REAL                                                        :: np_leaf_use          !! NP ratio used for allocationg biomass (gN/gP)
    REAL                                                        :: temp_som             !! value record som
    REAL                                                        :: need_som             !! value record how much som is needed for grass initialization
    REAL                                                        :: after_som            !! value record the rest som
    INTEGER                                                     :: k
!_ ================================================================================================================================

    IF (printlev>=3) WRITE(numout,*) 'Entering phenology'

    !
    !! 1. first call - output message giving the setting of the ::always_init
    !!    and ::min_growthinit_time parameters.
    !

    IF ( firstcall_all_phenology ) THEN

       IF (printlev >= 2) THEN
          WRITE(numout,*) 'phenology:'

          WRITE(numout,*) '   > take carbon from atmosphere if carbohydrate' // &
               ' reserve too small (::always_init): ', always_init
          
          WRITE(numout,*) '   > minimum time since last beginning of a growing' // &
               ' season (d) (::min_growthinit_time): ', min_growthinit_time
       END IF
       firstcall_all_phenology = .FALSE.

    ENDIF

    ! initialisation
    Cr_init=zero
    Cl_init=zero
    temp_som=zero
    need_som=zero
    after_som=zero
    !
    !! 2. Detection of the beginning of the growing season.
    !

    !
    !! 2.1 allow detection of the beginning of the growing season if dormance was
    !!     long enough (i.e. when ::time_lowgpp, which is calculated in ::stomate_season, 
    !!     is above a certain PFT-dependent threshold, ::lowgpp_time, 
    !!     which is given in ::stomate_constants),
    !!     AND the last beginning of growing season was a sufficiently long time ago 
    !!     (i.e. when ::when_growthinit, which is calculated in this module, 
    !!     is greater than ::min_growthinit_time, which is declared at the beginning of this module).
    !!     If these conditions are met, allow_initpheno is set to TRUE. Each PFT is looped over.
    !

    allow_initpheno(:,1) = .FALSE.
    DO j = 2,nvm

       WHERE ( when_growthinit(:,j) .GT. min_growthinit_time )
          allow_initpheno(:,j) = .TRUE.
          ! senescence (:,j) = .FALSE.
       ELSEWHERE
          allow_initpheno(:,j) = .FALSE.
       ENDWHERE

    ENDDO

    WHERE(allow_initpheno)
       histvar=un
    ELSEWHERE
       histvar=zero
    ENDWHERE

    CALL xios_orchidee_send_field("ALLOW_INITPHENO",histvar)
    CALL histwrite_p (hist_id_stomate, 'ALLOW_INITPHENO', itime, histvar, npts*nvm, horipft_index)

    !
    !! 2.2 increase the ::when_growthinit counter, which gives the number of days since the beginning of the growing season.
    !!     Needed for allocation and for the detection of the beginning of the growing season.
    !

    when_growthinit(:,:) = when_growthinit(:,:) + dt

    !
    !! 3. Leaf onset.
    !!    Check biometeorological conditions using the onset phenological models, 
    !!    which are different for each PFT group (i.e. grass versus tropical etc. 
    !!    See below for more detail on the different models and which PFTs use each model).
    !

    !! - By default: phenology does not start (::begin_leaves set to FALSE).
    begin_leaves(:,:) = .FALSE.

    !! - The onset phenology model is selected, (according to the parameter ::pheno_model, 
    !! which is initialised in stomate_data), and called.
    !! Each PFT is looped over (ignoring bare soil). 
    !! If conditions are favourable, begin_leaves is set to TRUE.
    
    ! parameter used in all the differents models of phenology 
    t_always = ZeroCelsius + t_always_add

    DO j = 2,nvm ! Loop over # PFTs

       SELECT CASE ( pheno_model(j) )

       CASE ( 'hum' )

          CALL pheno_hum (npts, j, PFTpresent, allow_initpheno, &
               moiavail_month, moiavail_week, &
               maxmoiavail_lastyear, minmoiavail_lastyear, &
               begin_leaves)

       CASE ( 'moi' )

          CALL pheno_moi (npts, j, PFTpresent, allow_initpheno, &
               time_hum_min, &
               moiavail_month, moiavail_week, &
               begin_leaves)


       CASE ( 'ncdgdd' )

          CALL pheno_ncdgdd (npts, j, PFTpresent, allow_initpheno, &
               ncd_dormance, gdd_midwinter, &
               t2m_month, t2m_week, begin_leaves)

       CASE ( 'ngd' )

          CALL pheno_ngd (npts, j, PFTpresent, allow_initpheno, ngd_minus5, &
               t2m_month, t2m_week, begin_leaves)

       CASE ( 'humgdd' )

          CALL pheno_humgdd (npts, j, PFTpresent, allow_initpheno, gdd_m5_dormance, &
               maxmoiavail_lastyear, minmoiavail_lastyear, &
               t2m_longterm, t2m_month, t2m_week, &
               moiavail_week, moiavail_month, &
               begin_leaves)

       CASE ( 'moigdd' )

          CALL pheno_moigdd (npts, j, PFTpresent, allow_initpheno, gdd_m5_dormance, &
               time_hum_min, &
               t2m_longterm, t2m_month, t2m_week, &
               moiavail_week, moiavail_month, &
               begin_leaves)

       CASE ( 'none' )

          ! no action

       CASE default

          WRITE(numout,*) 'phenology: don''t know how to treat this PFT.'
          WRITE(numout,*) '  number: (::j)',j
          WRITE(numout,*) '  phenology model (::pheno_model(j)) : ',pheno_model(j)
          CALL ipslerr_p(3,'stomate phenology','Cannot treat this PFT','','')

       END SELECT

    ENDDO

    WHERE(begin_leaves)
       histvar=un
    ELSEWHERE
       histvar=zero
    ENDWHERE

    CALL xios_orchidee_send_field("BEGIN_LEAVES",histvar)

    CALL histwrite_p (hist_id_stomate, 'BEGIN_LEAVES', itime, histvar, npts*nvm, horipft_index)

    !
    !! 4. Leaf growth and biomass allocation when leaf biomass is low.
    !!   Leaves start to grow if biometeorological conditions are favourable (::begin_leaves == TRUE) and if
    !!   leaf growth is allowed (::allow_initpheno == TRUE).
    !!   PFTs and then grid cells are looped over.
    !

    DO j = 2,nvm ! Loop over # PFTs

       age_reset(:) = .FALSE.
       
       DO i = 1,npts

          ! We might need the c0_alloc factor, so let's calculate it.
          c0_alloc = calculate_c0_alloc(i, j, tau_root(j), &
               tau_sap(j))

          ! initialize temporary variables that may be used
          temp_som=zero
          need_som=zero
          after_som=zero

          IF (begin_leaves(i,j) .AND. (SUM(biomass(i,j,:,icarbon)) .GE. min_stomate  &
               .OR. always_init .OR. senescence_type(j) .EQ. 'crop') ) THEN

             !! 4.1 First minimum biomass is calculated using the following equation:
             !! \latexonly
             !! \input{phenology_lm_min_eqn2.tex}
             !! \endlatexonly
             !! \n  
             lm_min = lai_initmin(j) / sla_calc(i,j)
  
             ! The minimum leaf biomass is prescribed by ::lm_min which in turn
             ! is basically prescribed through ::lai_initmin. However, lm_min could exceed
             ! the leaf mass that is required to respect the allometric relationships
             ! therefore this leaf biomass should be calculated
             
             ! Calculate the allocation factors see stomate_growth_fun_alloc.f90 and
             ! stomate_prescribe.f90 for more details 
             LF = c0_alloc * KF(i,j)



             IF ( is_tree(j) ) THEN

                ! Calculate the available biomass is roots, sapwood and leaves (gC tree-1)
                Cl_tree = biomass(i,j,ileaf,icarbon) 
                Cs_tree = ( biomass(i,j,isapabove,icarbon) + &
                     biomass(i,j,isapbelow,icarbon) )
                Cr_tree = biomass(i,j,iroot,icarbon)

                ! Calculate the stand structure (gC m-2)
                circ_class_ba = wood_to_ba(biomass(i,j,:,icarbon),j) 
                height = pipe_tune2(j)*(4/pi*circ_class_ba)**(pipe_tune3(j)/2)

                ! Calculate the leaves and roots that need to be grown see 
                ! stomate_growth_fun_alloc.f90 for more details 
                ! Note that the units for Cl_init and Cr_init are gC m-2
                Cl_init = MIN( lm_min, ( KF(i,j) * Cs_tree / height ) )
                Cr_init = Cl_init / LF

             ELSE  ! grass/crop
                !JC MOD033 when begin_leaves the KF and LF should be the maximum
                LF = c0_alloc * k_latosa_max(j) 
                Cl_init = lm_min
                Cr_init = Cl_init / LF

             ENDIF

             ! ==================================
             ! Calculate the C to nutrient ratios:  
             ! ==================================
             IF(biomass(i,j,ileaf,icarbon) .GT. min_stomate) THEN
                cn_leaf_use=biomass(i,j,ileaf,icarbon)/biomass(i,j,ileaf,initrogen)
                cp_leaf_use=biomass(i,j,ileaf,icarbon)/biomass(i,j,ileaf,iphosphorus)
                np_leaf_use=biomass(i,j,ileaf,initrogen)/biomass(i,j,ileaf,iphosphorus)
             ELSE
                ! during budburst (no leaves) we try to have a high nutrient content; to allow adaptation of
                ! ecosystem to long-term environmental changes, we use the average nutrient
                ! content which was achieved last growing season (highest doesn't
                ! make any sense):
                cn_leaf_use=cn_leaf_avg_season(i,j)
                np_leaf_use=np_leaf_avg_season(i,j)
                cp_leaf_use=cn_leaf_use*np_leaf_use
             ENDIF
             ! if we impose nutrient concentration we use _init (from run.def)
             ! or mean ratio from restart file (cn_leaf_avg_season &
             ! np_leaf_avg_season)
             IF (impose_cn) THEN
                IF (impose_cn_restart) THEN
                  cn_leaf_use=cn_leaf_avg_season(i,j)
                ELSE
                  cn_leaf_use=cn_leaf_init(j)
                ENDIF
             ENDIF
             IF (impose_np) THEN
                IF (impose_np_restart) THEN
                  np_leaf_use=np_leaf_avg_season(i,j)
                ELSE
                  np_leaf_use=np_leaf_init(j)
                ENDIF
             ENDIF
             ! restrict ratios to the range prescribed; it's needed because
             ! machine precission error can lead to ratio slightly outside the
             ! range
             cn_leaf_use = MAX(MIN(cn_leaf_use,cn_leaf_max(j)),cn_leaf_min(j))
             np_leaf_use = MAX(MIN(np_leaf_use,np_leaf_max(j)),np_leaf_min(j))
             ! cp_leaf_use should be calculated
             cp_leaf_use=cn_leaf_use*np_leaf_use

             ! ====================================
             ! Calculate nutrient availability:
             ! ====================================
             ! it is assumed that 90% of the available nutrients can be
             ! allocated to new leaves in a time step
             ! nitrogen available:
             n_avail = MAX(0.9*(biomass(i,j,ilabile,initrogen)   + biomass(i,j,icarbres,initrogen)),zero) 
             ! phosphorus available:
             p_avail = MAX(0.9*(biomass(i,j,ilabile,iphosphorus) + biomass(i,j,icarbres,iphosphorus)),zero) 

             ! ========================================
             ! Calculate if we need to grow more leaves 
             ! and how much leaves can be produced
             ! ========================================
             ! If leaf biomass is lower than the minimum biomass then biomass must 
             ! be allocated from the labile pool to leaves and roots.
             IF ( biomass(i,j,ileaf,icarbon) .LT. Cl_init ) THEN

                IF ( printlev>=4 .AND. j == test_pft .AND. i == test_grid ) THEN
                  WRITE(numout,*) 'We need to grow some leaves'
                ENDIF

                ! First calculate how much biomass is wanted/required 
                ! (::bm_wanted = 2 x the minimum leaf biomass).
                bm_wanted = (Cl_init -  biomass(i,j,ileaf,icarbon)) &
                     + MAX( Cr_init -  biomass(i,j,iroot,icarbon), zero)

                ! =====================================
                ! Specific setting for forced phenology
                ! =====================================
                ! If the biomass in the carbohydrate reserves is less than the required biomass
                ! take the required amount of carbon from the atmosphere and put it into the
                ! labile pool. This only occurs if the parameter ::always_init 
                ! (set at beginning of this ::subroutine) is TRUE. Default is FALSE. 
                !JC MOD038 for new cnp crop strategy, we have to always initial crops
                !JC MOD047 initialize grass from som to avoid soil depletion
                !DSG: this is not really nice
                ! this initialized grass from seed like crops:
                IF ( always_init .OR. senescence_type(j) .EQ. 'crop') THEN! .OR. &
                    !(.NOT. is_tree(j) .AND. natural(j)) ) THEN

                   bm_wanted_n = MAX( Cl_init/cn_leaf_use - biomass(i,j,ileaf,initrogen), zero) &
                        + MAX( Cr_init *fcn_root(j)/cn_leaf_use - biomass(i,j,iroot,initrogen), zero)

                   bm_wanted_p = MAX( Cl_init/cp_leaf_use - biomass(i,j,ileaf,iphosphorus), zero) &
                        + MAX( Cr_init *(fcn_root(j)*fnp_root(j))/cp_leaf_use - biomass(i,j,iroot,iphosphorus), zero)

                   IF ( (biomass(i,j,ilabile,icarbon)+biomass(i,j,icarbres,icarbon)) .LT. bm_wanted ) THEN
                      IF (always_init .OR. senescence_type(j) .EQ. 'crop') THEN
                      co2_to_bm(i,j) = co2_to_bm(i,j) + ( bm_wanted - (biomass(i,j,ilabile,icarbon)+biomass(i,j,icarbres,icarbon) )) / dt
                      biomass(i,j,ilabile,icarbon) = biomass(i,j,ilabile,icarbon) + &
                           ( bm_wanted -  (biomass(i,j,ilabile,icarbon)+biomass(i,j,icarbres,icarbon)) )
                      ELSE IF (.NOT. is_tree(j) .AND. natural(j)) THEN
                        ! for grasses, assuming seed from som
                        ! However, it is very DANGEROUS!!!!!
                        temp_som=zero
                        need_som=zero
                        after_som=zero
                        temp_som = SUM(som(i,:,j,icarbon))
                        need_som = bm_wanted - (biomass(i,j,ilabile,icarbon)+biomass(i,j,icarbres,icarbon) )
                        IF ( temp_som .GT. min_stomate .AND. &
                            temp_som .GT. need_som) THEN                          
                          after_som = temp_som - need_som
                          DO k = 1, ncarb
                            som(i,k,j,icarbon) = som(i,k,j,icarbon) * after_som / temp_som 
                          ENDDO
                          biomass(i,j,ilabile,icarbon) = biomass(i,j,ilabile,icarbon) + &
                               need_som
                        ENDIF
                      ENDIF
                   ENDIF

                   ! same for N
                   IF ( (biomass(i,j,ilabile,initrogen)+biomass(i,j,icarbres,initrogen)) .LT. bm_wanted_n ) THEN
                      IF (always_init .OR. senescence_type(j) .EQ. 'crop') THEN
                      !JC MOD038 for new cnp crop strategy, seed C:N 6-17 is about half of tissue 15-50 (Gan et
                      !al., 2011); seed N:P 2.6-7.4 is about 0.9 of tissue 3-8 (Be´langeret al.,
                      !2012) Here maize only, and for test, further detail adjustment is required.
                      n_to_bm(i,j) = n_to_bm(i,j) + 3.0 *( bm_wanted_n - (biomass(i,j,ilabile,initrogen)+biomass(i,j,icarbres,initrogen)) ) / dt
                      biomass(i,j,ilabile,initrogen) = biomass(i,j,ilabile,initrogen) + &
                           3.0 * ( bm_wanted_n -  (biomass(i,j,ilabile,initrogen)+biomass(i,j,icarbres,initrogen)) )
                      ! update nitrogen availability:
                      n_avail = MAX(0.9*(biomass(i,j,ilabile,initrogen)   + biomass(i,j,icarbres,initrogen)),zero) 
                      ELSE IF (.NOT. is_tree(j) .AND. natural(j)) THEN
                        ! for grasses, assuming seed from som
                        ! However, it is very DANGEROUS!!!!!
                        temp_som=zero
                        need_som=zero
                        after_som=zero
                        temp_som = SUM(som(i,:,j,initrogen))
                        need_som = 3.0 * (bm_wanted_n - (biomass(i,j,ilabile,initrogen)+biomass(i,j,icarbres,initrogen) ))
                        IF ( temp_som .GT. min_stomate .AND. &
                            temp_som .GT. need_som) THEN
                          after_som = temp_som - need_som
                          DO k = 1, ncarb
                            som(i,k,j,initrogen) = som(i,k,j,initrogen) * after_som / temp_som
                          ENDDO
                          biomass(i,j,ilabile,initrogen) = biomass(i,j,ilabile,initrogen) + &
                               need_som
                          ! update nitrogen availability:
                          n_avail = MAX(0.9*(biomass(i,j,ilabile,initrogen)   + biomass(i,j,icarbres,initrogen)),zero)
                        ENDIF
                      ENDIF

                   ENDIF
                   ! same for P
                   IF ( (biomass(i,j,ilabile,iphosphorus)+biomass(i,j,icarbres,iphosphorus)) .LT. bm_wanted_p ) THEN

                      IF (always_init .OR. senescence_type(j) .EQ. 'crop') THEN
                      p_to_bm(i,j) = p_to_bm(i,j) + 3.0 * 1.1 *( bm_wanted_p - (biomass(i,j,ilabile,iphosphorus)+biomass(i,j,icarbres,iphosphorus)) ) / dt
                      biomass(i,j,ilabile,iphosphorus) = biomass(i,j,ilabile,iphosphorus) + &
                           3.0 * 1.1 *( bm_wanted_p -  (biomass(i,j,ilabile,iphosphorus)+biomass(i,j,icarbres,iphosphorus)) )
                      ! update phosphorus availability:
                      p_avail = MAX(0.9*(biomass(i,j,ilabile,iphosphorus) + biomass(i,j,icarbres,iphosphorus)),zero) 
                      ELSE IF (.NOT. is_tree(j) .AND. natural(j)) THEN
                        ! for grasses, assuming seed from som
                        ! However, it is very DANGEROUS!!!!!
                        temp_som=zero
                        need_som=zero
                        after_som=zero
                        temp_som = SUM(som(i,:,j,iphosphorus))
                        need_som = 3.0 * 1.1 * (bm_wanted_p - (biomass(i,j,ilabile,iphosphorus)+biomass(i,j,icarbres,iphosphorus) ))
                        IF ( temp_som .GT. min_stomate .AND. &
                            temp_som .GT. need_som) THEN
                          after_som = temp_som - need_som
                          DO k = 1, ncarb
                            som(i,k,j,iphosphorus) = som(i,k,j,iphosphorus) * after_som / temp_som
                          ENDDO
                          biomass(i,j,ilabile,iphosphorus) = biomass(i,j,ilabile,iphosphorus) + &
                               need_som
                          ! update nitrogen availability:
                          p_avail = MAX(0.9*(biomass(i,j,ilabile,iphosphorus)   + biomass(i,j,icarbres,iphosphorus)),zero)
                        ENDIF
                      ENDIF
                   ENDIF
                ! to avoid error message for not initialized
                ELSE
                  bm_wanted_n = val_exp
                  bm_wanted_p = val_exp 
                ENDIF ! always_init

                !=======================================================
                ! Calculate how much biomass we need to allocate to have the
                ! mass of leaves and root we want:
                !=======================================================
                ! The biomass available for new leaves is set to be the minimum of the biomass of 
                ! the labile pool (if we do not take carbon the atmosphere) and the wanted biomass.
                bm_use = MIN( biomass(i,j,ilabile,icarbon) + biomass(i,j,icarbres,icarbon), &
                     bm_wanted )

                ! bm_use may differ from Cr_init and Cl_init, the final allocation will depend 
                ! on bm_use. Distrinute bm_use over leaves and roots following allometric relationships
                ! (i) bm_use = Cl_incp + Cr_incp
                ! (ii) Cr_incp = (Cl_incp+Cl)/LF - Cr
                ! Substitue (ii) in (i) and solve for Cl_inc
                ! <=> Cl_incp = (LF*(b_incp+Cr)-Cl)/(1+LF)
                ! Because this is the start of the growing season, Cr = Cl = 0
                !JC Debug MOD012
                ! For grasses, sometimes the Cr and Cl != 0
                ! Thus the Cl_init should not determined only by bm_use
                ! But the sum of bm_use and original Cl Cr
                ! otherwise the calculation of n_use and p_use will be wrong
                !                Cl_init = LF * bm_use / (un + LF )
                !JC Debug MOD032
                ! only if both ileaf < Cl_init and iroot < Cr_init
                ! we apply the following initialization for both ileaf and iroot
                IF (biomass(i,j,iroot,icarbon) .LT. Cr_init) THEN                

                  IF ( printlev>=4 .AND. j == test_pft .AND. i == test_grid ) THEN
                      WRITE(numout,*) 'We need to grow some roots'
                  ENDIF

                  Cl_init = LF * (bm_use+biomass(i,j,ileaf,icarbon)+biomass(i,j,iroot,icarbon)) / (un + LF)
                  Cr_init = Cl_init / LF
                ELSE ! root biomass >= Cr_init, bm_use is for Cl_init only
                  Cl_init = bm_use + biomass(i,j,ileaf,icarbon)
                  ! and we set new Cr_init as current root biomass
                  Cr_init = biomass(i,j,iroot,icarbon)
                  ! update LF
                  LF = Cl_init/Cr_init
                ENDIF
                !End JC Debug

                !=============================================================
                ! Calculate how much nitrogen we allocated, potentially reduce
                ! new biomass if not enough nitrogen is available
                !============================================================
                ! First: derive the amount of N needed for the carbon to be
                ! allocated assuming no change in C:N ratio
                !JC Debug MOD012
                ! in case Cr and Cl != 0, especially for the first year restarted from Conly run
                ! the C:N ratio of leaf and root could be inconsistent.
                ! thus there need to assume a translocation of N between leaf and root
                ! thus the calculation of n_use and p_use should be modified
                ! BUT it is possible this strategy does not fit for decideous forest, 
                ! need to discuss?
                !DSG: this can lead to problems as you might have not enough nutrients to
                !achieve this, or?
                n_use = MAX( (Cl_init+Cr_init *fcn_root(j))/cn_leaf_use - &
                              (biomass(i,j,ileaf,initrogen)+biomass(i,j,iroot,initrogen)), zero)

                IF ( printlev>=4 .AND. j == test_pft .AND. i == test_grid ) THEN
                   WRITE(numout,*) 'bm_use =', bm_use
                   WRITE(numout,*) 'Cl_init =', Cl_init
                   WRITE(numout,*) 'Cr_init =', Cr_init
                   WRITE(numout,*) 'cn_leaf_use  =', cn_leaf_use
                   WRITE(numout,*) 'fcn_root  =', fcn_root(j)
                   WRITE(numout,*) 'cp_leaf_use  =', cp_leaf_use
                   WRITE(numout,*) 'n_avail=', n_avail
                   WRITE(numout,*) 'n_use  =', n_use
                ENDIF

                ! in case there is not enough nutrients available we 
                ! have to make sure the following 3 points stay true
                !       cn_leaf_use < cn_leaf_max
                !       cp_leaf_use < cn_leaf_max*np_leaf_max
                !       np_leaf_min < np_leaf_use <np_leaf_max
                ! so ...

                ! ... in case there is not enough N available we ...
                IF (n_use .GT. n_avail) THEN

                   IF ( printlev>=4 .AND. j == test_pft .AND. i == test_grid ) THEN
                       WRITE(numout,*) 'not enough N for bm_use and current cn_leaf_use'
                   ENDIF

                   IF (impose_cn) THEN 
                    ! ... either add it magically in case we impose CN  ...
                       N_support(i,j) = N_support(i,j)+ (n_use &
                             - n_avail )

                       biomass(i,j,ilabile,initrogen) = biomass(i,j,ilabile,initrogen) + (n_use &
                             - n_avail )
                   ELSE
                      !  ... or try to use the leaf C:N which would be possible ...
                      cn_leaf_use=(Cl_init + Cr_init *fcn_root(j))/ &
                               (n_avail +                           &

                                   biomass(i,j,ileaf,initrogen)+biomass(i,j,iroot,initrogen))
                      cp_leaf_use = cn_leaf_use * np_leaf_use

                      ! ... if the CN gets too high than we need to reduce it
                      ! to the maximum allowed
                      IF(cn_leaf_use .GT. cn_leaf_max(j)) THEN
                         cn_leaf_use = cn_leaf_max(j)
                         cp_leaf_use = cn_leaf_use * np_leaf_use

                         IF ( printlev>=4 .AND. j == test_pft .AND. i == test_grid ) THEN
                             WRITE(numout,*) 'we use max CN ratio and reduce bm_use'
                         ENDIF

                         !JC Debug MOD032
                         ! here, we use new_LF instead LF in case there is overwinter root biomass
                         ! it should not affect other PFTs since new_LF = LF for trees
                         ! we update the carbon biomass which can be allocated
                         bm_use = (1+LF)/(LF+fcn_root(j))* &
                                      (( n_avail +         &
                                  biomass(i,j,ileaf,initrogen)+biomass(i,j,iroot,initrogen)) * cn_leaf_use) - &
                                  (biomass(i,j,ileaf,icarbon)+biomass(i,j,iroot,icarbon))
                         ! and also update Cl_init and Cr_init - as bm_use, Cl_init and Cr_init always have to be consistent 
                         ! due to the messy way of calculating 
                         Cl_init = LF *(bm_use+biomass(i,j,ileaf,icarbon)+biomass(i,j,iroot,icarbon)) / (un + LF)
                         Cr_init = Cl_init / LF
                      ELSE
                         IF ( printlev>=4 .AND. j == test_pft .AND. i == test_grid ) THEN
                             WRITE(numout,*) 'we use increased CN ratio, no need to reduce bm_use'
                         ENDIF

                      ENDIF

                   ENDIF ! impose_cn
                ENDIF ! (n_use .GT. n_avail)
                
                ! ..  now we check about P availability

                ! derive the amount of P needed to for the carbon to be
                ! allocated assuming no change in N:P ratio, but taking into
                ! account a potential increased C:N ratio
                p_use = MAX( (Cl_init+Cr_init*(fcn_root(j)*fnp_root(j)))/(cp_leaf_use) - &
                              (biomass(i,j,ileaf,iphosphorus)+biomass(i,j,iroot,iphosphorus)), zero)

                IF ( printlev>=4 .AND. j == test_pft .AND. i == test_grid ) THEN
                   WRITE(numout,*) 'p_avail=', p_avail
                   WRITE(numout,*) 'p_use  =', p_use
                ENDIF
                
                ! ... in case there is not enough P available we ...
                IF (p_use.GT.p_avail) THEN

                   IF (impose_np) THEN 
                      ! calculate the amount of P to be added
                       P_support(i,j) = P_support(i,j) + (p_use &
                                  - p_avail )

                      ! add P needed for the new leaf ... 
                       biomass(i,j,ilabile,iphosphorus) = biomass(i,j,ilabile,iphosphorus) + (p_use &
                      ! minus all P which can be supported by internal storage     
                                  - p_avail )
                   ELSE
                      !  ... or try to use the leaf N:P which would be possible ...
                      cp_leaf_use = (Cl_init + Cr_init*fcn_root(j)*fnp_root(j))/ &
                                               ( p_avail  +  &
                                   biomass(i,j,ileaf,iphosphorus)+biomass(i,j,iroot,iphosphorus))

                      ! if C:P ratio too high
                      ! bm_use should be limited by P
                      ! cn_leaf_use has been adjusted to the maximum that can be
                      ! supported by N
                      ! Thus if the following situation happen, it means that
                      ! current P can not support the bm_use calculated after N
                      ! limiatation
                      IF (cp_leaf_use .GT. (cn_leaf_use * np_leaf_max(j))) THEN
                         cp_leaf_use = cn_leaf_use * np_leaf_max(j)
                         np_leaf_use = np_leaf_max(j)
                        !JC Debug MOD032
                        ! here, we use new_LF instead LF in case there is overwinter root biomass
                        ! it should not affect other PFTs since new_LF = LF for trees

                        ! Now we now the possible CP stoichiometry so we can
                        ! recalculate bm_use
                        bm_use = (1+LF)/(LF+fcn_root(j)*fnp_root(j))* &
                                  (( p_avail + &
                                  biomass(i,j,ileaf,iphosphorus)+biomass(i,j,iroot,iphosphorus)) * cp_leaf_use) - &
                                  (biomass(i,j,ileaf,icarbon)+biomass(i,j,iroot,icarbon))
                        !DSG warning. In case BM_USE is reduced and we only
                        !want to grow only leaves OR only roots, we will have
                        !a instant transformation of the other tissue to the one
                        !we want to grow:
                        Cl_init = LF *(bm_use+biomass(i,j,ileaf,icarbon)+biomass(i,j,iroot,icarbon)) / (un + LF)
                        Cr_init = Cl_init / LF
                      ENDIF
                   ENDIF
                ENDIF

                ! Finally we know how much maximally could be supported by 
                ! nutrients

                ! first, we recalcuate N and P use

                n_use = MAX( (Cl_init+Cr_init *fcn_root(j))/cn_leaf_use - &
                              (biomass(i,j,ileaf,initrogen)+biomass(i,j,iroot,initrogen)), zero) 

                p_use = MAX( (Cl_init+Cr_init *(fcn_root(j)*fnp_root(j)))/(cp_leaf_use) - &
                              (biomass(i,j,ileaf,iphosphorus)+biomass(i,j,iroot,iphosphorus)), zero)
                 
                IF ( printlev>=4 .AND. j == test_pft .AND. i == test_grid ) THEN
                   WRITE(numout,*) 'CNP use as consistent with availabilities:'
                   WRITE(numout,*) 'bm_use=', bm_use
                   WRITE(numout,*) 'n_use  =', n_use
                   WRITE(numout,*) 'p_use  =', p_use
                   WRITE(numout,*) 'Cl_init =', Cl_init
                   WRITE(numout,*) 'Cr_init =', Cr_init
                   WRITE(numout,*) 'cn_leaf_use  =', cn_leaf_use
                   WRITE(numout,*) 'cp_leaf_use  =', cp_leaf_use
                ENDIF

                ! second, we remove the carbon to be allocated from the labile/reserve C pool 
                !JC MOD027 it is not suitable to delete this initialization, because when start
                !from scratch, grasses still need initialization (i.e., leaf/root biomass <
                !min_stomate
                ! In another case, if root >= min_stomate but leaf < min_stomate
                ! the initialization is needed too for leaf, to allow use_reserve in allocation
                !So, I change MOD023 a little, when there is absolutely no grasses or no leaf, do
                !initialization
                !JC MOD023 delete grasses initial because I change overwinter turnover
                IF (biomass(i,j,ileaf,icarbon) .LT. Cl_init) THEN

                  IF ( printlev>=4 .AND. j == test_pft .AND. i == test_grid ) THEN
                      WRITE(numout,*) 'LEAF MASS LESS THAN CL_INIT'
                  ENDIF

                  IF ( is_tree(j) .OR. &
                  (.NOT. is_tree(j) .AND. .NOT. natural(j))) THEN
                    ! Decrease labile pool biomass by the amount that's been allocated 
                    ! to the leaves and roots. If the labile pool is depleted use carbon
                    ! from the reserve pool.
                    deficit = bm_use - biomass(i,j,ilabile,icarbon)
   
                    ! There is enough carbon in the labile pool
                    IF (deficit .LT. zero) THEN
   
                       biomass(i,j,ilabile,icarbon) = biomass(i,j,ilabile,icarbon) - bm_use
   
                    ! Deplete the labile pool, use the reserve pool
                    ELSE
   
                       biomass(i,j,ilabile,icarbon) = zero
                       biomass(i,j,icarbres,icarbon) = MAX(biomass(i,j,icarbres,icarbon) - deficit, zero)                   
   
                    ENDIF

                  ELSE ! grasses
                    ! for grasses, only when no leaf,
                    ! initialization is done,
                    ! otherwise, it will be initialized by allocation
                    IF ((biomass(i,j,ileaf,icarbon) .LE. min_stomate)) THEN
                      ! Decrease labile pool biomass by the amount that's been
                      ! allocated
                      ! to the leaves and roots. If the labile pool is depleted
                      ! use carbon
                      ! from the reserve pool.
                      deficit = bm_use - biomass(i,j,ilabile,icarbon)

                      ! There is enough carbon in the labile pool
                      IF (deficit .LT. zero) THEN

                         biomass(i,j,ilabile,icarbon) = biomass(i,j,ilabile,icarbon) - bm_use

                      ! Deplete the labile pool, use the reserve pool
                      ELSE

                         biomass(i,j,ilabile,icarbon) = zero
                         biomass(i,j,icarbres,icarbon) = MAX(biomass(i,j,icarbres,icarbon) - deficit, zero)

                      ENDIF  
                    ENDIF ! grasses more conditions
                 ! third, we fill leaf and root pool with C,N and P
                 !        amd remove the N,P from the labile pools
               
                  ENDIF ! trees/crops and grasses

                 ! third, we fill leaf and root pool with C,N and P
                 !        amd remove the N,P from the labile pools
  !                 IF ( is_tree(j) ) THEN
                    IF ( is_tree(j) .OR.                                       & ! IF tree OR
                         (.NOT. is_tree(j) .AND. natural(j).AND.               & ! grass ..
                           (biomass(i,j,ileaf,icarbon) .LE. min_stomate))) THEN  !.. with no leaves
   
                      ! Distribute the biomass over the leaves and roots (gC tree-1)
                      biomass(i,j,ileaf,icarbon) = Cl_init
                      biomass(i,j,iroot,icarbon) = Cr_init 
                      ! in case we have increased the CN ratio  of leaves and root we might need to retranslocate
                      ! N from leaf & root to labile pool. Thus we record it ...
                      deficit = MAX( biomass(i,j,ileaf,initrogen) + biomass(i,j,iroot,initrogen) &
                                     - ( biomass(i,j,ileaf,icarbon)/cn_leaf_use                  &
                                         + biomass(i,j,iroot,icarbon)*fcn_root(j)/cn_leaf_use + n_use), zero) 
                      ! ... before we set new N concentration in leaf and roots:
                      biomass(i,j,ileaf,initrogen) = biomass(i,j,ileaf,icarbon)/cn_leaf_use
                      biomass(i,j,iroot,initrogen) = biomass(i,j,iroot,icarbon)*fcn_root(j)/cn_leaf_use
                      ! ... and also in reserves and labile pool: 
                      biomass(i,j,icarbres,initrogen)=biomass(i,j,icarbres,initrogen)+MIN(biomass(i,j,ilabile,initrogen)-n_use,zero)
                      biomass(i,j,ilabile,initrogen)= MAX( biomass(i,j,ilabile,initrogen) - n_use + deficit,zero)

                      ! in case we have increased the NP ratio  of leaves and root we might need to retranslocate
                      ! P from leaf & root to labile pool. Thus we record it ...
                      deficit = MAX(biomass(i,j,ileaf,iphosphorus) + biomass(i,j,iroot,iphosphorus)     &
                                    - (biomass(i,j,ileaf,icarbon)/cp_leaf_use                           &
                                       + biomass(i,j,iroot,icarbon)* (fcn_root(j)*fnp_root(j))/cp_leaf_use) + p_use,zero)
                      ! ... before we set new P concentration in leaf and roots:
                      biomass(i,j,ileaf,iphosphorus) = biomass(i,j,ileaf,icarbon)/cp_leaf_use
                      biomass(i,j,iroot,iphosphorus) = biomass(i,j,iroot,icarbon)* (fcn_root(j)*fnp_root(j))/cp_leaf_use
                      ! ... and also in reserves and labile pool: 
                      biomass(i,j,icarbres,iphosphorus)=biomass(i,j,icarbres,iphosphorus)+MIN(biomass(i,j,ilabile,iphosphorus)-p_use,zero)
                      biomass(i,j,ilabile,iphosphorus) = MAX( biomass(i,j,ilabile,iphosphorus) - p_use + deficit,zero)

                   ! crops
                   ! JC MOD023 only crop will initial lai, grasses have more conditions
                   ELSE IF (.NOT. natural(j) ) THEN
                  
                      ! Cl_init and Cr_init were recalculated to properly account for bm_use (the available C)
                      ! Distribute the biomass over the leaves and roots (gC tree-1)
                      biomass(i,j,ileaf,icarbon) = Cl_init
                      biomass(i,j,iroot,icarbon) = Cr_init 
                      ! in case we have increased the CN ratio  of leaves and root we might need to retranslocate
                      ! N from leaf & root to labile pool. Thus we record it ...
                      deficit = MAX( biomass(i,j,ileaf,initrogen) + biomass(i,j,iroot,initrogen) &
                                     - ( biomass(i,j,ileaf,icarbon)/cn_leaf_use                  &
                                         + biomass(i,j,iroot,icarbon)*fcn_root(j)/cn_leaf_use + n_use), zero) 
                      ! ... before we set new N concentration in leaf and roots:
                      biomass(i,j,ileaf,initrogen) = biomass(i,j,ileaf,icarbon)/cn_leaf_use
                      biomass(i,j,iroot,initrogen) = biomass(i,j,iroot,icarbon)*fcn_root(j)/cn_leaf_use
                      ! ... and also in reserves and labile pool: 
                      biomass(i,j,icarbres,initrogen)=biomass(i,j,icarbres,initrogen)+MIN(biomass(i,j,ilabile,initrogen)-n_use,zero)
                      biomass(i,j,ilabile,initrogen)= MAX( biomass(i,j,ilabile,initrogen) - n_use + deficit,zero)

                      ! in case we have increased the NP ratio  of leaves and root we might need to retranslocate
                      ! P from leaf & root to labile pool. Thus we record it ...
                      deficit = MAX(biomass(i,j,ileaf,iphosphorus) + biomass(i,j,iroot,iphosphorus)     &
                                    - (biomass(i,j,ileaf,icarbon)/cp_leaf_use                           &
                                       + biomass(i,j,iroot,icarbon)* (fcn_root(j)*fnp_root(j))/cp_leaf_use) + p_use,zero)
                      ! ... before we set new P concentration in leaf and roots:
                      biomass(i,j,ileaf,iphosphorus) = biomass(i,j,ileaf,icarbon)/cp_leaf_use
                      biomass(i,j,iroot,iphosphorus) = biomass(i,j,iroot,icarbon)* (fcn_root(j)*fnp_root(j))/cp_leaf_use
                      ! ... and also in reserves and labile pool: 
                      biomass(i,j,icarbres,iphosphorus)=biomass(i,j,icarbres,iphosphorus)+MIN(biomass(i,j,ilabile,iphosphorus)-p_use,zero)
                      biomass(i,j,ilabile,iphosphorus) = MAX( biomass(i,j,ilabile,iphosphorus) - p_use + deficit,zero)
                      
  !DSG: the following bit has still the bug and is anyway the same as for tree,
  !thus I merged this case with (is_tree)                    
                   ! for grasses, only when absolutely not leaf and root
                   ! DSG: ^^ this comment is not accurate, you do it when there is roots:
  !                 ELSE IF (.NOT. is_tree(j) .AND. natural(j) .AND. &
  !                         (biomass(i,j,ileaf,icarbon) .LE. min_stomate)) THEN
  !                    ! Cl_init and Cr_init were recalculated to properly
  !                    ! account for bm_use (the available C)
  !                    biomass(i,j,ileaf,icarbon) = Cl_init

 !!                     biomass(i,j,iroot,icarbon) = MAX( biomass(i,j,iroot,icarbon), Cr_init )
  !                    biomass(i,j,iroot,icarbon) =  Cr_init

  !                    biomass(i,j,ileaf,initrogen) = biomass(i,j,ileaf,icarbon)/cn_leaf_use
  !                    biomass(i,j,iroot,initrogen) = biomass(i,j,iroot,icarbon)*fcn_root(j)/cn_leaf_use

  !                    biomass(i,j,icarbres,initrogen)=biomass(i,j,icarbres,initrogen)+MIN(biomass(i,j,ilabile,initrogen)-n_use,0.)
  !                    biomass(i,j,ilabile,initrogen)= MAX( biomass(i,j,ilabile,initrogen) - n_use,0.)

  !                    biomass(i,j,ileaf,iphosphorus) = biomass(i,j,ileaf,icarbon)/cp_leaf_use
  !                    biomass(i,j,iroot,iphosphorus) = biomass(i,j,iroot,icarbon)*(fcn_root(j)*fnp_root(j))/cp_leaf_use

  !                    biomass(i,j,icarbres,iphosphorus)=biomass(i,j,icarbres,iphosphorus)+MIN(biomass(i,j,ilabile,iphosphorus)-p_use,0.)
  !                    biomass(i,j,ilabile,iphosphorus)= MAX( biomass(i,j,ilabile,iphosphorus) - p_use,0.)
                   ENDIF ! tree/crops/grasses
                   ! set reset leaf age distribution (::age_reset) flag. Default is TRUE.
                   ! done later for better vectorization)
                   age_reset(i) = .TRUE.
                ! JC seems not need this
                ENDIF ! leaf mass is very low (biomass(i,j,ileaf,icarbon) .LT. Cl_init)

             ENDIF  ! leaf mass is very low

             !! 4.3 reset when_growthinit counter: start of the growing season
             when_growthinit(i,j) = zero

          ENDIF    ! start of the growing season
          
       ENDDO      ! loop over grid points


       !! 4.4 reset leaf age distribution where necessary (i.e. when age_reset is TRUE)
       !! simply say that everything is in the youngest age class

       ! Set the youngest age class fraction to 1 and all other leaf age class fractions to 0.
       DO ipar = 1,nparts
         IF (ipar .NE. iroot) THEN
           WHERE ( age_reset(:) )
               leaf_frac(:,j,ipar,1) = un
           ENDWHERE
  
           DO m = 2, nleafages
  
              WHERE ( age_reset(:) )
                leaf_frac(:,j,ipar,m) = zero
              ENDWHERE
  
           ENDDO ! nleafages
  
           ! Ages - set all leaf ages to 0.
           DO m = 1, nleafages
  
              WHERE ( age_reset(:) )
                  
                 leaf_age(:,j,ipar,m) = zero
  
              ENDWHERE
  
           ENDDO ! nleafages
         ENDIF
       ENDDO ! npart
       
    ENDDO ! loop over # PFTs



    IF (printlev>=3) WRITE(numout,*) 'Leaving phenology'

  END SUBROUTINE phenology


!! ================================================================================================================================
!! SUBROUTINE   : pheno_hum 
!!
!>\BRIEF          The 'hum' onset model initiate leaf onset based exclusively on moisture 
!!                availability criteria. 
!!                Currently no PFTs are assigned to this onset model.
!!
!! DESCRIPTION  : This model is for tropical biomes, where temperatures are high but moisture
!!                might be a limiting factor on growth. It is based on leaf onset model 4a in 
!!                Botta et al. (2000), which adopts the approach of Le Roux (1995). \n
!!                Leaf onset occurs if the monthly moisture availability is still quite
!!                low (i.e. lower than the weekly availability), but the weekly availability is 
!!                higher than the critical threshold ::availability_crit (as it reacts faster), 
!!                which indicates the weekly moisture availability is increasing.
!!                OR if the monthly moisture availability is high enough (i.e. above the 
!!                threshold value ::moiavail_always), leaf onset is initiated if this has not 
!!                already happened. This allows vegetation in arid areas to respond to rapidly
!!                changing soil moisture conditions (Krinner et al., 2005). \n
!!                The critical weekly moisture availability threshold (::availability_crit), is
!!                calculated in this subroutine, and is a function of last year's maximum and
!!                minimum moisture availability and the PFT-dependent parameter
!!                ::hum_frac, which specifies how much of last year's available 
!!                moisture is required for leaf onset, as per the equation:
!!                \latexonly
!!                \input{phenology_moi_availcrit_eqn3.tex}
!!                \endlatexonly
!!                \n
!!                ::hum_frac is set for each PFT in ::stomate_data from a table
!!                which contains all the PFT values (::hum_frac_tab) in ::stomate_constants. \n
!!                Last year's maximum and minimum moisture availability and the monthly and 
!!                weekly moisture availability are  
!!                The ::pheno_hum subroutine is called in the subroutine ::phenology. 
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): ::begin_leaves - specifies whether leaf growth can start.
!!
!! REFERENCE(S) : 
!! - Botta, A., N. Viovy, P. Ciais, P. Friedlingstein and P. Monfray (2000), 
!! A global prognostic scheme of leaf onset using satellite data,
!! Global Change Biology, 207, 337-347.
!! - Le Roux, X. (1995), Etude et modelisation des echanges d'eau et d'energie
!! sol-vegetation-atmosphere dans une savane humide, PhD Thesis, University
!! Pierre et Marie Curie, Paris, France.
!! - Krinner, G., N. Viovy, N. de Noblet-Ducoudre, J. Ogee, J. Polcher, P. 
!! Friedlingstein, P. Ciais, S. Sitch and I.C. Prentice (2005), A dynamic global
!! vegetation model for studies of the coupled atmosphere-biosphere system, Global
!! Biogeochemical Cycles, 19, doi:10.1029/2003GB002199.
!!
!! FLOWCHART    : 
!! \latexonly
!! \includegraphics[scale = 1]{pheno_hum.png}
!! \endlatexonly
!! \n             
!_ ================================================================================================================================

  SUBROUTINE pheno_hum (npts, j, PFTpresent, allow_initpheno, &
       moiavail_month, moiavail_week, &
       maxmoiavail_lastyear, minmoiavail_lastyear, &
       begin_leaves)

    !
    !! 0. Variable and parameter declarations
    !

    !
    !! 0.1 Input variables
    !
    INTEGER, INTENT(in)                                             :: npts                  !! Domain size - number of 
                                                                                                    !! grid cells (unitless) 
    INTEGER, INTENT(in)                                             :: j                     !! PFT index (unitless)
    LOGICAL, DIMENSION(npts,nvm), INTENT(in)                               :: PFTpresent            !! PFT exists (true/false)
    LOGICAL, DIMENSION(npts,nvm), INTENT(in)                               :: allow_initpheno       !! are we allowed to 
                                                                                                    !! declare the beginning of 
                                                                                                    !! the growing season? 
                                                                                                    !! (true/false) 
    REAL, DIMENSION(npts,nvm), INTENT(in)                           :: moiavail_month        !! "monthly" moisture 
                                                                                                    !! availability (0-1, unitless)
    REAL, DIMENSION(npts,nvm), INTENT(in)                           :: moiavail_week         !! "weekly" moisture 
                                                                                                    !! availability (0-1, unitless)
    REAL, DIMENSION(npts,nvm), INTENT(in)                           :: maxmoiavail_lastyear  !! last year's maximum 
                                                                                                    !! moisture availability 
                                                                                                    !! (0-1, unitless)
    REAL, DIMENSION(npts,nvm), INTENT(in)                           :: minmoiavail_lastyear  !! last year's minimum 
                                                                                                    !! moisture availability 
                                                                                                    !! (0-1, unitless)

    !
    !! 0.2 Output variables
    !

    !
    !! 0.3 Modified variables
    !
    LOGICAL, DIMENSION(npts,nvm), INTENT(inout)                            :: begin_leaves          !! signal to start putting 
                                                                                                    !! leaves on (true/false) 

    !
    !! 0.4 Local variables
    !
    REAL                                                            :: moiavail_always       !! critical monthly 
                                                                                                    !! moisture availability - set 
                                                                                                    !! for tree or grass 
                                                                                                    !! (0-1, unitless)
    REAL, DIMENSION(npts)                                           :: availability_crit     !! critical weekly moisture 
                                                                                                    !! availability (0-1, unitless)
    INTEGER                                                         :: i                     !! index (unitless)

!_ ================================================================================================================================

    IF (printlev>=3) WRITE(numout,*) 'Entering hum'

    !
    !! 1. Initializations
    !

    !
    !! 1.1 first call - outputs the name of onset model and the moisture availability 
    !!     parameters for tree and grass
    !

    IF ( firstcall_hum ) THEN

       IF (printlev >= 2) THEN
          WRITE(numout,*) 'pheno_hum:'
          WRITE(numout,*) '   > moisture availability above which moisture tendency doesn''t matter: '
          WRITE(numout,*) '         trees (::moiavail_always_tree): ', moiavail_always_tree
          WRITE(numout,*) '         grasses (::moiavail_always_grass):', moiavail_always_grass
       END IF
       firstcall_hum = .FALSE.

    ENDIF

    !
    !! 1.2 initialize output
    !

    begin_leaves(:,j) = .FALSE.

    !
    !! 1.3 check the critical value ::hum_frac is defined. If not, stop.
    !

    IF ( hum_frac(j) .EQ. undef ) THEN

       WRITE(numout,*) 'hum: hum_frac is undefined for PFT (::j)',j
       CALL ipslerr_p(3,'stomate phenology','hum_frac is undefined for this PFT','','')

    ENDIF

    !
    !! 1.4 set the critical monthly moisture availability above which we always detect the beginning of the
    !!     growing season - set as the moisture availability for trees or grass.
    !

    IF ( is_tree(j) ) THEN
       moiavail_always = moiavail_always_tree
    ELSE
       moiavail_always = moiavail_always_grass
    ENDIF

    !
    !! 2. Check if biometeorological conditions are favourable for leaf growth.
    !! The PFT has to be there and start of growing season must be allowed
    !

    DO i = 1, npts

       IF ( PFTpresent(i,j) .AND. allow_initpheno(i,j) ) THEN

          !! 2.1 Calculate the critical weekly moisture availability: depends linearly on the last year 
          !! minimum and maximum moisture availabilities, and on the parameter ::hum_frac.

          availability_crit(i) = minmoiavail_lastyear(i,j) + hum_frac(j) * &
               ( maxmoiavail_lastyear(i,j) - minmoiavail_lastyear(i,j) )

          !! 2.2 Determine if growing season should start (if so, ::begin_leaves set to TRUE).
          !!     Leaf onset occurs if the monthly moisture availability is still quite
          !!     low (i.e. lower than the weekly availability), but the weekly availability is 
          !!     already higher than the critical threshold ::availability_crit (as it reacts faster), 
          !!     which indicates the weekly moisture availability is increasing.
          !!     OR if the monthly moisture availability is high enough (i.e. above the threshold value 
          !!     ::moiavail_always), leaf onset is initiated if this has not already happened.

          IF ( ( ( moiavail_week(i,j)  .GE. availability_crit(i) ) .AND. &
               ( moiavail_month(i,j) .LT. moiavail_week(i,j) )   ) .OR. &
               ( moiavail_month(i,j) .GE. moiavail_always )                ) THEN
             begin_leaves(i,j) = .TRUE.
          ENDIF

       ENDIF        ! PFT there and start of growing season allowed

    ENDDO ! end loop over grid points

    IF (printlev>=4) WRITE(numout,*) 'Leaving hum'

  END SUBROUTINE pheno_hum


!! ================================================================================================================================
!! SUBROUTINE   : pheno_moi
!!
!>\BRIEF          The 'moi' onset model (::pheno_moi) initiates leaf onset based exclusively 
!!                on moisture availability criteria. 
!!                It is very similar to the 'hum' onset model but instead of the weekly moisture 
!!                availability being higher than a constant threshold, the condition is that the 
!!                moisture minimum happened a sufficiently long time ago. 
!!                Currently PFT 3 (Tropical Broad-leaved Raingreen) is assigned to this model.
!!
!! DESCRIPTION  : This model is for tropical biomes, where temperatures are high but moisture
!!                might be a limiting factor on growth. It is based on leaf onset model 4b in 
!!                Botta et al. (2000).
!!                Leaf onset begins if the plant moisture availability minimum was a sufficiently  
!!                time ago, as specified by the PFT-dependent parameter ::hum_min_time 
!!                AND if the "monthly" moisture availability is lower than the "weekly"
!!                availability (indicating that soil moisture is increasing).
!!                OR if the monthly moisture availability is high enough (i.e. above the threshold 
!!                value ::moiavail_always), leaf onset is initiated if this has not already 
!!                happened. \n
!!                ::hum_min_time is set for each PFT in ::stomate_data, and is 
!!                defined in the table ::hum_min_time_tab in ::stomate_constants. \n
!!                ::moiavail_always is defined for both tree and grass in this subroutine 
!!                (set to 1. and 0.6 respectively). \n
!!                The ::pheno_moi subroutine is called in the subroutine ::phenology. 
!!
!! RECENT CHANGE(S): None
!!        
!! MAIN OUTPUT VARIABLE(S): ::begin_leaves - specifies whether leaf growth can start.
!!
!! REFERENCE(S) : 
!! - Botta, A., N. Viovy, P. Ciais, P. Friedlingstein and P. Monfray (2000), 
!! A global prognostic scheme of leaf onset using satellite data,
!! Global Change Biology, 207, 337-347.
!! - Krinner, G., N. Viovy, N. de Noblet-Ducoudre, J. Ogee, J. Polcher, P. 
!! Friedlingstein, P. Ciais, S. Sitch and I.C. Prentice (2005), A dynamic global
!! vegetation model for studies of the coupled atmosphere-biosphere system, Global
!! Biogeochemical Cycles, 19, doi:10.1029/2003GB002199.
!!
!! FLOWCHART    : 
!! \latexonly
!! \includegraphics[scale = 1]{pheno_moi.png}
!! \endlatexonly
!! \n
!_ ================================================================================================================================

  SUBROUTINE pheno_moi (npts, j, PFTpresent, allow_initpheno, &
       time_hum_min, &
       moiavail_month, moiavail_week, &
       begin_leaves)

    !
    !! 0. Variable and parameter declaration
    !

    !
    !! 0.1 Input variables
    !
    INTEGER, INTENT(in)                               :: npts            !! Domain size - number of grid cells (unitless)
    INTEGER, INTENT(in)                               :: j               !! PFT index (unitless)
    LOGICAL, DIMENSION(npts,nvm), INTENT(in)                 :: PFTpresent      !! PFT exists (true/false)
    LOGICAL, DIMENSION(npts,nvm), INTENT(in)                 :: allow_initpheno !! are we allowed to declare the beginning of the 
                                                                                !! growing season? (true/false) 
    REAL, DIMENSION(npts,nvm), INTENT(in)             :: time_hum_min    !! time elapsed since strongest moisture 
                                                                                !! availability (days) 
    REAL, DIMENSION(npts,nvm), INTENT(in)             :: moiavail_month  !! "monthly" moisture availability (0-1, unitless)
    REAL, DIMENSION(npts,nvm), INTENT(in)             :: moiavail_week   !! "weekly" moisture availability (0-1, unitless)

    !
    !! 0.2 Output variables
    !
    LOGICAL, DIMENSION(npts,nvm), INTENT(inout)              :: begin_leaves    !! signal to start putting leaves on (true/false)

    !
    !! 0.3 Modified variables
    !

    !
    !! 0.4 Local variables
    !
    REAL                                              :: moiavail_always                 !! critical moisture availability - 
                                                                                                !! set for tree or grass 
                                                                                                !! (0-1, unitless)
    INTEGER                                           :: i                               !! index (unitless)

!_ ================================================================================================================================

    IF (printlev>=3) WRITE(numout,*) 'Entering moi'

    !
    !! 1. Initializations
    !

    !
    !! 1.1 first call - outputs the name of onset model and the moisture availability 
    !!     parameters for tree and grass
    !

    IF ( firstcall_moi ) THEN
       IF (printlev >= 2) THEN
          WRITE(numout,*) 'pheno_moi:'
          WRITE(numout,*) '   > moisture availability above which moisture tendency doesn''t matter: '
          WRITE(numout,*) '         trees (::moiavail_always_tree):', moiavail_always_tree
          WRITE(numout,*) '         grasses (::moiavail_always_grass):', moiavail_always_grass
       END IF
       firstcall_moi = .FALSE.

    ENDIF

    !
    !! 1.2 initialize output
    !

    begin_leaves(:,j) = .FALSE.

    !
    !! 1.3 check the critical value ::hum_min_time is definded. If not, stop
    !

    IF ( hum_min_time(j) .EQ. undef ) THEN

       WRITE(numout,*) 'moi: hum_min_time is undefined for PFT (::j) ',j
       CALL ipslerr_p(3,'stomate phenology','hum_min_time is undefined for this PFT','','')

    ENDIF

    !
    !! 1.4 set the critical monthly moisture availability above which we always detect the beginning of the
    !!     growing season - set as the moisture availability for trees or grass.
    !

    IF ( is_tree(j) ) THEN
       moiavail_always = moiavail_always_tree
    ELSE
       moiavail_always = moiavail_always_grass
    ENDIF

    !
    !! 2. Check if biometeorological conditions are favourable for leaf growth.
    !! The PFT has to be there and start of growing season must be allowed.
    !

    DO i = 1, npts

       IF ( PFTpresent(i,j) .AND. allow_initpheno(i,j) ) THEN
          
          !! 2.1 Determine if growing season should start (if so, ::begin_leaves set to TRUE).
          !!     The favorable season starts if the moisture minimum (::time_hum_min) was a sufficiently long 
          !!     time ago, i.e. greater than the threshold specified by the parameter ::hum_min_time 
          !!     and if the "monthly" moisture availability is lower than the "weekly"
          !!     availability (indicating that soil moisture is increasing).
          !!     OR if the monthly moisture availability is high enough (i.e. above the threshold value 
          !!     ::moiavail_always), initiate the growing season if this has not happened yet.

          IF  ( ( ( moiavail_week(i,j) .GT. moiavail_month(i,j) ) .AND. &
               ( time_hum_min(i,j) .GT. hum_min_time(j) )    ) .OR. &
               ( moiavail_month(i,j) .GE. moiavail_always )                     ) THEN
             begin_leaves(i,j) = .TRUE.
          ENDIF

       ENDIF        ! PFT there and start of growing season allowed

    ENDDO ! end loop over grid points

    IF (printlev>=4) WRITE(numout,*) 'Leaving moi'

  END SUBROUTINE pheno_moi


!! ================================================================================================================================
!! SUBROUTINE   : pheno_humgdd
!!
!>\BRIEF          The 'humgdd' onset model initiates leaf onset based on mixed conditions of 
!!                temperature and moisture availability criteria. 
!!                Currently no PFTs are assigned to this onset model. 
!!
!! DESCRIPTION  : In this model the Growing Degree Day (GDD) model (Chuine, 2000) is combined 
!!                with the 'hum' onset model (::pheno_hum), which has previously been described,
!!                in order to account for dependence on both temperature and moisture conditions 
!!                in warmer climates. \n. 
!!                The GDD model specifies that daily temperatures above a threshold of -5  
!!                degrees C are summed, minus this threshold, giving the GDD, starting from 
!!                the beginning of the dormancy period (::time_lowgpp>0), i.e. since the leaves 
!!                were lost. \n.
!!                The dormancy time-length is represented by the variable 
!!                ::time_lowgpp, which is calculated in ::stomate_season. It is increased by 
!!                the stomate time step when the weekly GPP is lower than a threshold. Otherwise
!!                it is set to zero. \n
!!                Leaf onset begins when the a PFT-dependent GDD-threshold is reached.
!!                In addition there are temperature and moisture conditions.
!!                The temperature condition specifies that the monthly temperature has to be 
!!                higher than a constant threshold (::t_always) OR
!!                the weekly temperature is higher than the monthly temperature.
!!                There has to be at least some moisture. The moisture condition 
!!                is exactly the same as the 'hum' onset model (::pheno_hum), which has already
!!                been described. \n
!!                The GDD (::gdd_m5_dormance) is calculated in ::stomate_season. GDD is set to 
!!                undef if beginning of the growing season detected, i.e. when there is GPP 
!!                (::time_lowgpp>0).
!!                The parameter ::t_always is defined as 10 degrees C in this subroutine, 
!!                as are the parameters ::moisture_avail_tree and ::moisture_avail_grass 
!!                (set to 1 and 0.6 respectively), which are used in the moisture condition 
!!                (see ::pheno_moi onset model description). \n
!!                The PFT-dependent GDD threshold (::gdd_crit) is calculated as in the onset 
!!                model ::pheno_humgdd, using the equation:
!!                \latexonly
!!                \input{phenology_hummoigdd_gddcrit_eqn4.tex}
!!                \endlatexonly
!!                \n
!!                The three GDDcrit parameters (::gdd(j,*)) are set for each PFT in 
!!                ::stomate_data, and three tables defining each of the three critical GDD
!!                parameters for each PFT is given in ::gdd_crit1_tab, ::gdd_crit2_tab and 
!!                ::gdd_crit3_tab in ::stomate_constants. \n
!!                The ::pheno_humgdd subroutine is called in the subroutine ::phenology. 
!!
!! RECENT CHANGES: None
!!                
!! MAIN OUTPUT VARIABLES: ::begin_leaves - specifies whether leaf growth can start
!!
!! REFERENCE(S) : 
!! - Botta, A., N. Viovy, P. Ciais, P. Friedlingstein and P. Monfray (2000), 
!! A global prognostic scheme of leaf onset using satellite data,
!! Global Change Biology, 207, 337-347.
!! - Chuine, I (2000), A unified model for the budburst of trees, Journal of 
!! Theoretical Biology, 207, 337-347.
!! - Krinner, G., N. Viovy, N. de Noblet-Ducoudre, J. Ogee, J. Polcher, P. 
!! Friedlingstein, P. Ciais, S. Sitch and I.C. Prentice (2005), A dynamic global
!! vegetation model for studies of the coupled atmosphere-biosphere system, Global
!! Biogeochemical Cycles, 19, doi:10.1029/2003GB002199.
!!
!! FLOWCHART    : 
!! \latexonly
!! \includegraphics[scale = 1]{pheno_humgdd.png}
!! \endlatexonly
!! \n             
!_ ================================================================================================================================

  SUBROUTINE pheno_humgdd (npts, j, PFTpresent, allow_initpheno, gdd, &
       maxmoiavail_lastyear, minmoiavail_lastyear, &
       t2m_longterm, t2m_month, t2m_week, &
       moiavail_week, moiavail_month, &
       begin_leaves)

    !
    !! 0. Variable and parameter declaration
    !

    !
    !! 0.1 Input variables
    !
    INTEGER, INTENT(in)                               :: npts                    !! Domain size - number of grid cells 
                                                                                        !! (unitless) 
    INTEGER, INTENT(in)                               :: j                       !! PFT index (unitless)
    LOGICAL, DIMENSION(npts,nvm), INTENT(in)                 :: PFTpresent              !! PFT exists (true/false)
    LOGICAL, DIMENSION(npts,nvm), INTENT(in)                 :: allow_initpheno         !! are we allowed to declare the beginning 
                                                                                        !! of the growing season? (true/false) 
    REAL, DIMENSION(npts,nvm), INTENT(in)             :: gdd                     !! growing degree days, calculated since 
                                                                                        !! leaves have fallen (C) 
    REAL, DIMENSION(npts,nvm), INTENT(in)             :: maxmoiavail_lastyear    !! last year's maximum moisture 
                                                                                        !! availability (0-1, unitless)
    REAL, DIMENSION(npts,nvm), INTENT(in)             :: minmoiavail_lastyear    !! last year's minimum moisture 
                                                                                        !! availability (0-1, unitless)
    REAL, DIMENSION(npts), INTENT(in)                 :: t2m_longterm            !! "long term" 2 meter temperatures (K)
    REAL, DIMENSION(npts), INTENT(in)                 :: t2m_month               !! "monthly" 2-meter temperatures (K)
    REAL, DIMENSION(npts), INTENT(in)                 :: t2m_week                !! "weekly" 2-meter temperatures (K)
    REAL, DIMENSION(npts,nvm), INTENT(in)             :: moiavail_week           !! "weekly" moisture availability 
                                                                                        !! (0-1, unitless)
    REAL, DIMENSION(npts,nvm), INTENT(in)             :: moiavail_month          !! "monthly" moisture availability 
                                                                                        !! (0-1, unitless)

    !
    !! 0.2 Output variables
    !
    LOGICAL, DIMENSION(npts,nvm), INTENT(inout)              :: begin_leaves            !! signal to start putting leaves on 
                                                                                        !! (true/false) 

    !
    !! 0.3 Modified variables
    !

    !
    !! 0.4 Local variables
    !
    REAL                                              :: moiavail_always                 !! critical moisture availability - 
                                                                                                !! set for tree or grass 
                                                                                                !! (0-1, unitless)
    REAL, DIMENSION(npts)                             :: moiavail_crit                   !! critical moisture availability 
                                                                                                !! (0-1, unitless)
    REAL, DIMENSION(npts)                             :: tl                              !! long term temperature (C)
    REAL, DIMENSION(npts)                             :: gdd_crit                        !! critical GDD (C)
    INTEGER                                           :: i                               !! index (unitless)

!_ ================================================================================================================================

    IF (printlev>=3) WRITE(numout,*) 'Entering humgdd'

    !
    !! 1. Initializations
    !

    !
    !! 1.1 first call - outputs the name of the onset model, the values of the  
    !!     moisture availability parameters for tree and grass, and the value of the 
    !!     critical monthly temperature.
    !

    IF ( firstcall_humgdd ) THEN

       IF (printlev >= 2) THEN
          WRITE(numout,*) 'pheno_humgdd:'
          WRITE(numout,*) '   > moisture availability above which moisture tendency doesn''t matter: '
          WRITE(numout,*) '         trees (::moiavail_always_tree): ', moiavail_always_tree
          WRITE(numout,*) '         grasses (::moiavail_always_grass): ', moiavail_always_grass
          WRITE(numout,*) '   > monthly temp. above which temp. tendency doesn''t matter: ', &
               t_always
       END IF

       firstcall_humgdd = .FALSE.

    ENDIF

    !
    !! 1.2 initialize output
    !

    begin_leaves(:,j) = .FALSE.

    !
    !! 1.3 check the critical values ::gdd and ::pheno_crit_hum_frac are defined.
    !!     If not, stop.
    !

    IF ( ANY(pheno_gdd_crit(j,:) .EQ. undef) ) THEN

       WRITE(numout,*) 'humgdd: pheno_gdd_crit is undefined for PFT (::j) ',j
       CALL ipslerr_p(3,'stomate phenology','pheno_gdd_crit is undefined for this PFT','','')

    ENDIF

    IF ( hum_frac(j) .EQ. undef ) THEN

       WRITE(numout,*) 'humgdd: hum_frac is undefined for PFT (::j) ',j
       CALL ipslerr_p(3,'stomate phenology','hum_frac is undefined for this PFT','','')

    ENDIF

    !
    !! 1.4 set the critical moisture availability above which we always detect the beginning of the
    !!     growing season - set as the moisture availability for trees or grass.
    !

    IF ( is_tree(j) ) THEN
       moiavail_always = moiavail_always_tree
    ELSE
       moiavail_always = moiavail_always_grass
    ENDIF

    !
    !! 2. Check if biometeorological conditions are favourable for leaf growth.
    !!   The PFT has to be there, start of growing season must be allowed, 
    !!   and GDD has to be defined.
    !

    DO i = 1, npts

       IF ( PFTpresent(i,j) .AND. allow_initpheno(i,j) .AND. &
            ( gdd(i,j) .NE. undef )                           ) THEN

          !! 2.1 Calculate the critical weekly moisture availability: depends linearly on the last year 
          !! minimum and maximum moisture availabilities, and on the parameter ::hum_frac.,
          !! (as in the ::pheno_hum model), as per the equation:

          moiavail_crit(i) = minmoiavail_lastyear(i,j) + hum_frac(j) * &
               ( maxmoiavail_lastyear(i,j) - minmoiavail_lastyear(i,j) )

          !! 2.2 Calculate the critical GDD (::gdd_crit), which is a function of the PFT-dependent 
          !!     critical GDD and the "long term" 2 meter air temperatures.  

          tl(i) =  t2m_longterm(i) - ZeroCelsius
          gdd_crit(i) = pheno_gdd_crit(j,1) + tl(i)*pheno_gdd_crit(j,2) + &
               tl(i)*tl(i)*pheno_gdd_crit(j,3)
          
          !! 2.3 Determine if the growing season should start (if so, ::begin_leaves set to TRUE).
          !!     - Has the critical gdd been reached and is the temperature increasing?
          !!     - Is there at least some humidity/moisture availability?
          !!     This occurs if the critical gdd (::gdd_crit) has been reached 
          !!     AND that is temperature increasing, which is true either if the monthly
          !!     temperature being higher than the threshold ::t_always, OR if the weekly
          !!     temperature is higher than the monthly, 
          !!     AND finally that there is sufficient moisture availability, which is 
          !!     the same condition as for the ::pheno_hum onset model.

          IF ( ( gdd(i,j) .GE. gdd_crit(i) ) .AND. &
               ( ( t2m_week(i) .GT. t2m_month(i) ) .OR. &
               ( t2m_month(i) .GT. t_always )          ) .AND. &
               ( ( ( moiavail_week(i,j)  .GE. moiavail_crit(i) ) .AND. &
               ( moiavail_month(i,j) .LT. moiavail_crit(i) )        ) .OR. &
               ( moiavail_month(i,j) .GE. moiavail_always )                   ) )  THEN
             begin_leaves(i,j) = .TRUE.
          ENDIF

       ENDIF        ! PFT there and start of growing season allowed

    ENDDO ! End loop over grid points

    IF (printlev>=4) WRITE(numout,*) 'Leaving humgdd'

  END SUBROUTINE pheno_humgdd


!! ================================================================================================================================
!! SUBROUTINE   : pheno_moigdd
!!
!>\BRIEF          The 'moigdd' onset model initiates leaf onset based on mixed temperature 
!!                and moisture availability criteria.
!!                Currently PFTs 10 - 13 (C3 and C4 grass, and C3 and C4 agriculture) 
!!                are assigned to this model. 
!!
!! DESCRIPTION  : This onset model combines the GDD model (Chuine, 2000), as described for 
!!                the 'humgdd' onset model (::pheno_humgdd), and the 'moi' model, in order 
!!                to account for dependence on both temperature and moisture conditions in
!!                warmer climates. \n
!!                Leaf onset begins when the a PFT-dependent GDD threshold is reached.
!!                In addition there are temperature and moisture conditions.
!!                The temperature condition specifies that the monthly temperature has to be 
!!                higher than a constant threshold (::t_always) OR
!!                the weekly temperature is higher than the monthly temperature.
!!                There has to be at least some moisture. The moisture condition 
!!                is exactly the same as the 'moi' onset model (::pheno_moi), which has
!!                already been described. \n
!!                GDD is set to undef if beginning of the growing season detected.
!!                As in the ::pheno_humgdd model, the parameter ::t_always is defined as 
!!                10 degrees C in this subroutine, as are the parameters ::moisture_avail_tree
!!                and ::moisture_avail_grass (set to 1 and 0.6 respectively), which are used
!!                in the moisture condition (see ::pheno_moi onset model description). \n
!!                The PFT-dependent GDD threshold (::gdd_crit) is calculated as in the onset 
!!                model ::pheno_humgdd, using the equation:
!!                \latexonly
!!                \input{phenology_hummoigdd_gddcrit_eqn4.tex}
!!                \endlatexonly
!!                \n
!!                where i and j are the grid cell and PFT respectively.
!!                The three GDDcrit parameters (::gdd(j,*)) are set for each PFT in 
!!                ::stomate_data, and three tables defining each of the three critical GDD
!!                parameters for each PFT is given in ::gdd_crit1_tab, ::gdd_crit2_tab and 
!!                ::gdd_crit3_tab in ::stomate_constants. \n
!!                The ::pheno_moigdd subroutine is called in the subroutine ::phenology. 
!!
!! RECENT CHANGE(S): Added temperature threshold for C4 grass (pheno_moigdd_t_crit), Dan Zhu april 2015
!!                
!! MAIN OUTPUT VARIABLE(S): ::begin_leaves - specifies whether leaf growth can start
!!
!! REFERENCE(S) : 
!! - Botta, A., N. Viovy, P. Ciais, P. Friedlingstein and P. Monfray (2000), 
!! A global prognostic scheme of leaf onset using satellite data,
!! Global Change Biology, 207, 337-347.
!! - Chuine, I (2000), A unified model for the budburst of trees, Journal of 
!! Theoretical Biology, 207, 337-347.
!! - Krinner, G., N. Viovy, N. de Noblet-Ducoudre, J. Ogee, J. Polcher, P. 
!! Friedlingstein, P. Ciais, S. Sitch and I.C. Prentice (2005), A dynamic global
!! vegetation model for studies of the coupled atmosphere-biosphere system, Global
!! Biogeochemical Cycles, 19, doi:10.1029/2003GB002199.
!! - Still et al., Global distribution of C3 and C4 vegetation: Carbon cycle implications, 
!! 2003, Global Biogeochemmical Cycles, DOI: 10.1029/2001GB001807. 
!!
!! FLOWCHART    : 
!! \latexonly
!! \includegraphics[scale = 1]{pheno_moigdd.png}
!! \endlatexonly
!! \n
!_ ================================================================================================================================

  SUBROUTINE pheno_moigdd (npts, j, PFTpresent, allow_initpheno, gdd, &
       time_hum_min, &
       t2m_longterm, t2m_month, t2m_week, &
       moiavail_week, moiavail_month, &
       begin_leaves)

    !
    !! 0. Variable and parameter declaration
    !

    !
    !! 0.1 Input variables
    !
    INTEGER, INTENT(in)                               :: npts            !! Domain size - number of grid cells (unitless)
    INTEGER, INTENT(in)                               :: j               !! PFT index (unitless)
    LOGICAL, DIMENSION(npts,nvm), INTENT(in)                 :: PFTpresent      !! PFT exists (true/false)
    LOGICAL, DIMENSION(npts,nvm), INTENT(in)                 :: allow_initpheno !! are we allowed to decalre the beginning of the 
                                                                                !! growing season? (true/false) 
    REAL, DIMENSION(npts,nvm), INTENT(in)             :: gdd             !! growing degree days, calculated since leaves 
                                                                                !! have fallen (C) 
    REAL, DIMENSION(npts,nvm), INTENT(in)             :: time_hum_min    !! time elapsed since strongest moisture 
                                                                                !! availability (days) 
    REAL, DIMENSION(npts), INTENT(in)                 :: t2m_longterm    !! "long term" 2 meter temperatures (K)
    REAL, DIMENSION(npts), INTENT(in)                 :: t2m_month       !! "monthly" 2-meter temperatures (K)
    REAL, DIMENSION(npts), INTENT(in)                 :: t2m_week        !! "weekly" 2-meter temperatures (K)
    REAL, DIMENSION(npts,nvm), INTENT(in)             :: moiavail_week   !! "weekly" moisture availability (0-1, unitless)
    REAL, DIMENSION(npts,nvm), INTENT(in)             :: moiavail_month  !! "monthly" moisture availability (0-1, unitless)

    !
    !! 0.2 Output variables
    !

    !
    !! 0.3 Modified variables
    !
    LOGICAL, DIMENSION(npts,nvm), INTENT(inout)              :: begin_leaves    !! signal to start putting leaves on (true/false)

    !
    !! 0.4 Local variables
    !
    REAL                                              :: moiavail_always !! critical moisture availability - 
                                                                                !! set for tree or grass 
                                                                                !! (0-1, unitless)
    REAL, DIMENSION(npts)                             :: tl              !! long term temperature (C)
    REAL, DIMENSION(npts)                             :: gdd_crit        !! critical GDD (C)
    INTEGER                                           :: i               !! index (unitless)

!_ ================================================================================================================================

    IF (printlev>=3) WRITE(numout,*) 'Entering moigdd'

    !
    !! 1. Initializations
    !

    !
    !! 1.1 first call - outputs the name of the onset model, the values of the  
    !!     moisture availability parameters for tree and grass, and the value of the 
    !!     critical monthly temperature.
    !

    IF ( firstcall_moigdd ) THEN

       IF (printlev >= 2) THEN
          WRITE(numout,*) 'pheno_moigdd:'
          WRITE(numout,*) '   > moisture availability above which moisture tendency doesn''t matter: '
          WRITE(numout,*) '         trees (::moiavail_always_tree) :', moiavail_always_tree
          WRITE(numout,*) '         grasses (::moiavail_always_grass) :', moiavail_always_grass
          WRITE(numout,*) '   > monthly temp. above which temp. tendency doesn''t matter (::t_always): ', &
               t_always
       END IF
       firstcall_moigdd = .FALSE.

    ENDIF

    !
    !! 1.2 initialize output
    !

    begin_leaves(:,j) = .FALSE.

    !
    !! 1.3 check the critical values ::gdd and ::pheno_crit_hum_min_time are defined.
    !!     If not, stop.
    !

    IF ( ANY(pheno_gdd_crit(j,:) .EQ. undef) ) THEN

       WRITE(numout,*) 'moigdd: pheno_gdd_crit is undefined for PFT',j
       CALL ipslerr_p(3,'stomate phenology','pheno_gdd is undefined for this PFT','','')

    ENDIF

    IF ( hum_min_time(j) .EQ. undef ) THEN

       WRITE(numout,*) 'moigdd: hum_min_time is undefined for PFT',j
       CALL ipslerr_p(3,'stomate phenology','hum_min is undefined for this PFT','','')

    ENDIF

    !
    !! 1.4 set the critical moisture availability above which we always detect the beginning of the
    !!     growing season - set as the moisture availability for trees or grass.
    !

    IF ( is_tree(j) ) THEN
       moiavail_always = moiavail_always_tree
    ELSE
       moiavail_always = moiavail_always_grass
    ENDIF

    !
    !! 2. Check if biometeorological conditions are favourable for leaf growth.
    !!    The PFT has to be there, the start of growing season must be allowed, 
    !!    and GDD has to be defined.
    !

    DO i = 1, npts

       IF ( PFTpresent(i,j) .AND. allow_initpheno(i,j) .AND. &
            ( gdd(i,j) .NE. undef )                           ) THEN
          
          !! 2.1 Calculate the critical GDD (::gdd_crit), which is a function of the PFT-dependent 
          !!     critical GDD and the "long term" 2 meter air temperatures 
          
          tl(i) = t2m_longterm(i) - ZeroCelsius
          gdd_crit(i) = pheno_gdd_crit(j,1) + tl(i)*pheno_gdd_crit(j,2) + &
               tl(i)*tl(i)*pheno_gdd_crit(j,3)

          !! 2.2 Determine if the growing season should start (if so, ::begin_leaves set to TRUE).
          !!     This occurs if the critical gdd (::gdd_crit) has been reached 
          !!     AND that is temperature increasing, which is true either if the monthly
          !!     temperature being higher than the threshold ::t_always, OR if the weekly
          !!     temperature is higher than the monthly, 
          !!     AND finally that there is sufficient moisture availability, which is 
          !!     the same condition as for the ::pheno_moi onset model.
          !!     AND when pheno_moigdd_t_crit is set(for C4 grass), if the average temperature threshold is reached

          IF ( ( gdd(i,j) .GE. gdd_crit(i) ) .AND. &
               ( ( t2m_week(i) .GT. t2m_month(i) ) .OR. &
                 ( t2m_month(i) .GT. t_always )  ) .AND. &
               ( ( ( time_hum_min(i,j) .GT. hum_min_time(j) ) .AND. &
                 ( moiavail_week(i,j) .GT. moiavail_month(i,j) ) ) .OR. &
                 ( moiavail_month(i,j) .GE. moiavail_always )  ) .AND. &
               ( ( pheno_moigdd_t_crit(j) == undef ) .OR. &
                 (t2m_month(i) .GT. (ZeroCelsius + pheno_moigdd_t_crit(j))) ) ) THEN

             begin_leaves(i,j) = .TRUE.
             
          ENDIF

       ENDIF        ! PFT there and start of growing season allowed

    ENDDO

    IF (printlev>=4) WRITE(numout,*) 'Leaving moigdd'

  END SUBROUTINE pheno_moigdd


!! ================================================================================================================================
!! SUBROUTINE   : pheno_ncdgdd
!!
!>\BRIEF          The Number of Chilling Days - Growing Degree Day (NCD-GDD) model initiates 
!!                leaf onset if a certain relationship between the number of chilling days (NCD) 
!!                since leaves were lost, and the growing degree days (GDD) since midwinter, is 
!!                fulfilled. 
!!                Currently PFT 6 (Temperate Broad-leaved Summergreen) and PFT 8 (Boreal Broad-
!!                leaved Summergreen) are assigned to this model.
!! 
!! DESCRIPTION  : Experiments have shown that some
!!                species have a "chilling" requirement, i.e. their physiology needs cold 
!!                temperatures to trigger the mechanism that will allow the following budburst 
!!                (e.g. Orlandi et al., 2004). 
!!                An increase in chilling days, defined as a day with a daily mean air
!!                temperature below a PFT-dependent threshold, reduces a plant's GDD demand 
!!                (Cannell and Smith, 1986; Murray et al., (1989); Botta et al., 2000).
!!                The GDD threshold therefore decreases as NCD 
!!                increases, using the following empirical negative explonential law:
!!                \latexonly
!!                \input{phenology_ncdgdd_gddmin_eqn5.tex}
!!                \endlatexonly
!!                \n
!!                The constants used have been calibrated against data CHECK FOR REFERENCE OR PERSON WHO DID UPDATE.
!!                Leaf onset begins if the GDD is higher than the calculated minimum GDD
!!                (dependent upon NCD) AND if the weekly temperature is higher than the monthly 
!!                temperature. This is to ensure the temperature is increasing. \n
!!                The dormancy time-length is represented by the variable 
!!                ::time_lowgpp, which is calculated in ::stomate_season. It is increased by 
!!                the stomate time step when the weekly GPP is lower than a threshold. Otherwise
!!                it is set to zero. \n
!!                The NCD (::ncd_dormance) is calculated in ::stomate_season as  
!!                the number of days with a temperature below a PFT-dependent constant threshold
!!                (::ncdgdd_temp), starting from the beginning of the dormancy period
!!                (::time_lowgpp>0), i.e. since the leaves were lost. \n
!!                The growing degree day sum of the temperatures higher than 
!!                ::ncdgdd_temp (GDD) since midwinter (::gdd_midwinter) 
!!                is also calculated in ::stomate_season.
!!                Midwinter is detected if the monthly temperature is lower than the weekly
!!                temperature AND  the monthly temperature is lower than the long-term
!!                temperature. ::gdd_minter is therefore set to 0 at the beginning of midwinter
!!                and increased with each temperature greater than the PFT-dependent threshold.
!!                When midsummer is detected (the opposite of the above conditions), 
!!                ::gdd_midwinter is set to undef.
!!                CHECK! WHEN TO START OF DORMANCY BEEN MODIFIED FROM BOTTA- ADD IN?
!!                The ::pheno_ncdgdd subroutine is called in the subroutine ::phenology. 
!!
!! RECENT CHANGE(S): None
!!                
!! MAIN OUTPUT VARIABLE(S): ::begin_leaves - specifies whether leaf growth can start
!!
!! REFERENCE(S) :
!! - Botta, A., N. Viovy, P. Ciais, P. Friedlingstein and P. Monfray (2000), 
!! A global prognostic scheme of leaf onset using satellite data,
!! Global Change Biology, 207, 337-347.
!! - Cannell, M.J.R. and R.I. Smith (1986), Climatic warming, spring budburst and
!! frost damage on trees, Journal of Applied Ecology, 23, 177-191.
!! - Krinner, G., N. Viovy, N. de Noblet-Ducoudre, J. Ogee, J. Polcher, P. 
!! Friedlingstein, P. Ciais, S. Sitch and I.C. Prentice (2005), A dynamic global
!! vegetation model for studies of the coupled atmosphere-biosphere system, Global
!! Biogeochemical Cycles, 19, doi:10.1029/2003GB002199.
!! - Murray, M.B., G.R. Cannell and R.I. Smith (1989), Date of budburst of fifteen 
!! tree species in Britain following climatic warming, Journal of Applied Ecology,
!! 26, 693-700.
!! - Orlandi, F., H. Garcia-Mozo, L.V. Ezquerra, B. Romano, E. Dominquez, C. Galan,
!! and M. Fornaciari (2004), Phenological olive chilling requirements in Umbria
!! (Italy) and Andalusia (Spain), Plant Biosystems, 138, 111-116. 
!!
!! FLOWCHART    : 
!! \latexonly
!! \includegraphics[scale = 1]{pheno_ncdgdd.png}
!! \endlatexonly
!! \n
!_ ================================================================================================================================

  SUBROUTINE pheno_ncdgdd (npts, j, PFTpresent, allow_initpheno, &
       ncd_dormance, gdd_midwinter, &
       t2m_month, t2m_week, begin_leaves)

    !
    !! 0. Variable and parameter declaration
    !

    !
    !! 0.1 Input variables
    !
    INTEGER, INTENT(in)                               :: npts            !! Domain size - number of grid cells (unitless)
    INTEGER, INTENT(in)                               :: j               !! PFT index (unitless)
    LOGICAL, DIMENSION(npts,nvm), INTENT(in)                 :: PFTpresent      !! PFT exists (true/false)
    LOGICAL, DIMENSION(npts,nvm), INTENT(in)                 :: allow_initpheno !! are we allowed to declare the beginning of the 
                                                                                !! growing season? (true/false) 
    REAL, DIMENSION(npts,nvm), INTENT(in)             :: ncd_dormance    !! number of chilling days since leaves were lost 
                                                                                !! (days) 
    REAL, DIMENSION(npts,nvm), INTENT(inout)          :: gdd_midwinter   !! growing degree days since midwinter (C)
    REAL, DIMENSION(npts), INTENT(in)                 :: t2m_month       !! "monthly" 2-meter temperatures (K)
    REAL, DIMENSION(npts), INTENT(in)                 :: t2m_week        !! "weekly" 2-meter temperatures (K)

    !
    !! 0.2 Output variables
    !

    !
    !! 0.3 Modified variables
    !
    LOGICAL, DIMENSION(npts,nvm), INTENT(inout)              :: begin_leaves    !! signal to start putting leaves on (true/false)

    !
    !! 0.4 Local variables
    !
    INTEGER                                           :: i               !! index (unitless)
    REAL                                              :: gdd_min         !! critical gdd (C)

!_ ================================================================================================================================

    IF (printlev>=3) WRITE(numout,*) 'Entering ncdgdd'

    !
    !! 1. Initializations
    !

    !
    !! 1.1 initialize output
    !

    begin_leaves(:,j) = .FALSE.

    !
    !! 1.2 check the critical value ::ncdgdd_temp is defined.
    !!     If not, stop.
    !

    IF ( ncdgdd_temp(j) .EQ. undef ) THEN

       WRITE(numout,*) 'ncdgdd: ncdgdd_temp is undefined for PFT (::j) ',j
       CALL ipslerr_p(3,'stomate phenology','ncdgdd_temp this PFT','','')

    ENDIF

    !
    !! 2. Check if biometeorological conditions are favourable for leaf growth.    
    !!    PFT has to be there and start of growing season must be allowed.
    !

    DO i = 1, npts ! loop over grid points

       IF ( PFTpresent(i,j) .AND. allow_initpheno(i,j) .AND. &
            ( gdd_midwinter(i,j) .NE. undef ) .AND. &
            ( ncd_dormance(i,j) .NE. undef )                  ) THEN

          !! 2.1 Calculate the critical gdd, which is related to ::ncd_dormance
          !!     using an empirical negative exponential law as described above.           

          gdd_min = ( gddncd_ref / exp(gddncd_curve*ncd_dormance(i,j)) - gddncd_offset )

          !! 2.2 Determine if the growing season should start (if so, ::begin_leaves set to TRUE).
          !!     This occurs if the critical GDD been reached AND the temperatures are increasing.
          !!     If the growing season has started, ::gdd_midwinter is set to "undef". 

          IF ( ( gdd_midwinter(i,j) .GE. gdd_min ) .AND. &
               ( t2m_week(i) .GT. t2m_month(i) ) ) THEN
             begin_leaves(i,j) = .TRUE.
             gdd_midwinter(i,j)=undef
          ENDIF

       ENDIF        ! PFT there and start of growing season allowed

    ENDDO ! end loop over grid points

    IF (printlev>=4) WRITE(numout,*) 'Leaving ncdgdd'

  END SUBROUTINE pheno_ncdgdd


!! ================================================================================================================================
!! SUBROUTINE   : pheno_ngd
!!
!>\BRIEF          The Number of Growing Days (NGD) leaf onset model initiates leaf onset if the NGD, 
!!                defined as the number of days with temperature above a constant threshold, 
!!                exceeds a critical value.
!!                Currently PFT 9 (Boreal Leedleleaf Summergreen) is assigned to this model. 
!!
!! DESCRIPTION    The NGD model is a variant of the GDD model. The model was proposed by Botta et
!!                al. (2000) for boreal and arctic biomes, and is designed to estimate 
!!                leaf onset after the end of soil frost. 
!!                The NDG (::ngd_minus5) is the number of days with a daily mean air 
!!                temperature of greater than -5 degrees C, 
!!                starting from the beginning of the dormancy period (i.e. time since the leaves 
!!                were lost/GPP below a certain threshold).
!!                Leaf onset begins if the NGD is higher than the PFT-dependent constant threshold, 
!!                ::ngd,  AND if the weekly temperature is higher than the monthly 
!!                temperature. \n
!!                The dormancy time-length is represented by the variable 
!!                ::time_lowgpp, which is calculated in ::stomate_season. It is increased by 
!!                the stomate time step when the weekly GPP is lower than a threshold. Otherwise
!!                it is set to zero. \n
!!                ::ngd_minus5 is also calculated in ::stomate_season. It is initialised at the
!!                beginning of the dormancy period (::time_lowgpp>0), and increased by the 
!!                stomate time step when the temperature > -5 degrees C. \n
!!                ::ngd is set for each PFT in ::stomate_data, and a 
!!                table defining the minimum NGD for each PFT is given in ::ngd_crit_tab
!!                in ::stomate_constants. \n
!!                The ::pheno_ngd subroutine is called in the subroutine ::phenology.      
!!
!! RECENT CHANGE(S): None
!!                
!! MAIN OUTPUT VARIABLE(S): ::begin_leaves - specifies whether leaf growth can start
!!
!! REFERENCE(S) : 
!! - Botta, A., N. Viovy, P. Ciais, P. Friedlingstein and P. Monfray (2000), 
!! A global prognostic scheme of leaf onset using satellite data,
!! Global Change Biology, 207, 337-347. 
!! - Krinner, G., N. Viovy, N. de Noblet-Ducoudre, J. Ogee, J. Polcher, P. 
!! Friedlingstein, P. Ciais, S. Sitch and I.C. Prentice (2005), A dynamic global
!! vegetation model for studies of the coupled atmosphere-biosphere system, Global
!! Biogeochemical Cycles, 19, doi:10.1029/2003GB002199.
!!
!! FLOWCHART    : 
!! \latexonly
!! \includegraphics[scale = 1]{pheno_ngd.png}
!! \endlatexonly
!! \n
!_ ================================================================================================================================

  SUBROUTINE pheno_ngd (npts, j, PFTpresent, allow_initpheno, ngd, &
       t2m_month, t2m_week, begin_leaves)

    !
    !! 0. Variable and parameter declaration
    !

    !
    !! 0.1 Input variables
    !
    INTEGER, INTENT(in)                               :: npts            !! Domain size - number of grid cells (unitless)
    INTEGER, INTENT(in)                               :: j               !! PFT index (unitless)
    LOGICAL, DIMENSION(npts,nvm), INTENT(in)                 :: PFTpresent      !! PFT exists (true/false)
    LOGICAL, DIMENSION(npts,nvm), INTENT(in)                 :: allow_initpheno !! are we allowed to declare the beginning of the 
                                                                                !! growing season? (true/false) 
    REAL, DIMENSION(npts,nvm), INTENT(in)             :: ngd             !! growing degree days (C)
    REAL, DIMENSION(npts), INTENT(in)                 :: t2m_month       !! "monthly" 2-meter temperatures (K)
    REAL, DIMENSION(npts), INTENT(in)                 :: t2m_week        !! "weekly" 2-meter temperatures (K)

    !
    !! 0.2 Output variables
    !

    !
    !! 0.3 Modified variables
    !
    LOGICAL, DIMENSION(npts,nvm), INTENT(inout)              :: begin_leaves    !! signal to start putting leaves on (true/false)

    !
    !! 0.4 Local variables
    !
    INTEGER                                           :: i               !! index (unitless)

    !! =========================================================================

    IF (printlev>=3) WRITE(numout,*) 'Entering ngd'

    !
    !! 1. Initializations
    !

    !
    !! 1.1 initialize output
    !

    begin_leaves(:,j) = .FALSE.

    !
    !! 1.2 check the critical value ::ngd_crit is defined.
    !!     If not, stop.
    !

    IF ( ngd_crit(j) .EQ. undef ) THEN

       WRITE(numout,*) 'ngd: ngd_crit is undefined for PFT (::j) ',j
       CALL ipslerr_p(3,'stomate phenology','ngd_crit is undefined for this PFT','','')

    ENDIF

    !
    !! 2. Check if biometeorological conditions are favourable for leaf growth.
    !!    PFT has to be there and start of growing season must be allowed.
    !

    DO i = 1, npts

       IF ( PFTpresent(i,j) .AND. allow_initpheno(i,j) ) THEN

          !! 2.1 Determine if the growing season should start (if so, ::begin_leaves set to TRUE).
          !!     This occurs if the critical NGD has been reached AND are temperatures increasing.

          IF ( ( ngd(i,j) .GE. ngd_crit(j) ) .AND. &
               ( t2m_week(i) .GT. t2m_month(i) )        ) THEN
             begin_leaves(i,j) = .TRUE.
          ENDIF

       ENDIF        ! PFT there and start of growing season allowed

    ENDDO ! end loop over grid points

    IF (printlev>=4) WRITE(numout,*) 'Leaving ngd'

  END SUBROUTINE pheno_ngd

END MODULE stomate_phenology
