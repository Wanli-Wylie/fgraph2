! ================================================================================================================================
! MODULE       : stomate_lpj
!
! CONTACT      : orchidee-help _at_ listes.ipsl.fr
!
! LICENCE      : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF       Main entry point for daily processes in STOMATE and LPJ (phenology, 
!! allocation, kill, turn, light, establish, crown, cover, lcchange)
!!
!!\n DESCRIPTION: None
!!
!! RECENT CHANGE(S): None
!!
!! REFERENCE(S) : None
!!
!! SVN          :
!! $HeadURL: svn://forge.ipsl.jussieu.fr/orchidee/branches/ORCHIDEE-CN-P/ORCHIDEE/src_stomate/stomate_lpj.f90 $
!! $Date: 2021-05-06 11:05:55 +0200 (四, 2021-05-06) $
!! $Revision: 7177 $
!! \n
!_ ================================================================================================================================

MODULE stomate_lpj

  ! modules used:

  USE ioipsl_para
  USE xios_orchidee
  USE grid
  USE stomate_data
  USE constantes
  USE constantes_soil
  USE pft_parameters
  USE lpj_constraints
  USE lpj_pftinout
  USE lpj_kill
  USE lpj_crown
  USE lpj_fire
  USE lpj_gap
  USE lpj_light
  USE lpj_establish
  USE lpj_cover
  USE stomate_prescribe
  USE stomate_phenology
  USE stomate_growth_fun_all
  USE stomate_turnover
  USE stomate_litter
  USE stomate_som_dynamics
  USE stomate_vmax
  USE stomate_lcchange
  USE stomate_phosphorus
!  USE Orch_Write_field_p
  USE stomate_woodharvest
  USE grassland_constantes
  USE grassland_management
  USE time

  IMPLICIT NONE

  ! private & public routines

  PRIVATE
  PUBLIC StomateLpj,StomateLpj_clear

CONTAINS


!! ================================================================================================================================
!! SUBROUTINE   : StomateLpj_clear
!!
!>\BRIEF        Re-initialisation of variable
!!
!! DESCRIPTION  : This subroutine reinitializes variables. To be used if we want to relaunch 
!! ORCHIDEE but the routine is not used in current version.
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S): None
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE StomateLpj_clear

    CALL prescribe_clear
    CALL phenology_clear
    CALL turn_clear
    CALL som_dynamics_clear
    CALL constraints_clear
    CALL establish_clear
    CALL fire_clear
    CALL gap_clear
    CALL light_clear
    CALL pftinout_clear
  END SUBROUTINE StomateLpj_clear


!! ================================================================================================================================
!! SUBROUTINE   : StomateLPJ
!!
!>\BRIEF        Main entry point for daily processes in STOMATE and LPJ, structures the call sequence 
!!              to the different processes such as dispersion, establishment, competition and mortality of PFT's.
!! 
!! DESCRIPTION  : This routine is the main entry point to all processes calculated on a 
!! daily time step. Is mainly devoted to call the different STOMATE and LPJ routines 
!! depending of the ok_dgvm (is dynamic veg used) and lpj_constant_mortality (is background mortality used).
!! It also prepares the cumulative 
!! fluxes or pools (e.g TOTAL_M TOTAL_BM_LITTER etc...)
!!
!! This routine makes frequent use of "weekly", "monthly" and "long term" variables. Quotion is used because
!! by default "weekly" denotes 7 days, by default "monthly" denotes 20 days and by default "Long term" denotes
!! 3 years. dtslow refers to 24 hours (1 day).
!!
!!
!! RECENT CHANGE(S) : None
!! 
!! MAIN OUTPUT VARIABLE(S): All variables related to stomate and required for LPJ dynamic vegetation mode.
!!
!! REFERENCE(S) : 
!! - Krinner, G., N. Viovy, N. de Noblet-Ducoudré, J. Ogeé, J. Polcher, P. Friedlingstein, P. Ciais, S. Sitch, 
!! and I. C. Prentice. 2005. A dynamic global vegetation model for studies of the coupled atmosphere-biosphere 
!! system. Global Biogeochemical Cycles 19:GB1015, doi:1010.1029/2003GB002199.
!! - Sitch, S., B. Smith, I. C. Prentice, A. Arneth, A. Bondeau, W. Cramer, J. O. Kaplan, S. Levis, W. Lucht, 
!! M. T. Sykes, K. Thonicke, and S. Venevsky. 2003. Evaluation of ecosystem dynamics, plant geography and 
!! terrestrial carbon cycling in the LPJ dynamic global vegetation model. Global Change Biology 9:161-185.
!!
!! FLOWCHART    : Update with existing flowchart from N Viovy (Jan 19, 2012)
!! \n
!_ ================================================================================================================================
 
  SUBROUTINE StomateLpj (npts, dt_days, &
       neighbours, resolution, &
       herbivores, &
       soil_orders,max_eau_var, bulk, &
       tsurf_daily, tsoil_daily, t2m_daily, t2m_min_daily, &
       litterhum_daily, soilhum_daily, &
       maxmoiavail_lastyear, minmoiavail_lastyear, &
       gdd0_lastyear, precip_lastyear, &
       moiavail_month, moiavail_week, t2m_longterm, t2m_month, t2m_week, &
       tsoil_month, soilhum_month, &
       gdd_m5_dormance, gdd_from_growthinit, gdd_midwinter, ncd_dormance, ngd_minus5, &
       turnover_longterm, gpp_daily, &
!DSGlab_fac
       gpp_week, &
!DSG
       time_hum_min, maxfpc_lastyear, resp_maint_part, &
       PFTpresent, age, fireindex, firelitter, &
       leaf_age, leaf_frac, biomass, ind, adapted, regenerate, &
       senescence, when_growthinit, &
       litter, &
       litter_avail, litter_not_avail, litter_avail_frac, &
       dead_leaves, som, lignin_struc, &
       lignin_wood, veget_cov_max,veget_cov, veget_cov_max_new, woodharvest ,npp_longterm, lm_lastyearmax, veget_lastlight, &
       everywhere, need_adjacent, RIP_time, &
       lai, rprof,npp_daily, turnover_daily, turnover_time,&
       control_moist, control_temp, som_input, &
       co2_to_bm, n_to_bm, p_to_bm,co2_fire, resp_hetero, resp_maint, resp_growth, &
       height, deadleaf_cover, vcmax, &
       jmax,                          &
       nue,bm_to_litter, &
       prod10,prod100,flux10, flux100, &
       convflux,cflux_prod10,cflux_prod100, &
       prod10_harvest,prod100_harvest,flux10_harvest, flux100_harvest, &
       convflux_harvest,cflux_prod10_harvest,cflux_prod100_harvest, woodharvestpft, &
       harvest_above, harvest_bioenergy,carb_mass_total, fpc_max, MatrixA, &
       Tseason, Tmin_spring_time, begin_leaves, onset_date, KF, k_latosa_adapt, &
       cn_leaf_avg_season, np_leaf_avg_season, nstress_season, &
       pstress_season, &
       lai_target, lai_target_longterm, &
       moiavail_growingseason,soil_n_min, soil_p_min,rue_longterm, &
       n_uptake_daily, p_uptake_daily, &
       record_reserve_grass, &
       wshtotsum, sr_ugb, compt_ugb, nb_ani, grazed_frac, &
       import_yield, sla_age1, t2m_14, sla_calc, snowfall_daily, day_of_year, &
       when_growthinit_cut, nb_grazingdays, &
       moiavail_daily,tmc_topgrass_daily,fc_grazing,snowmass_daily,&
       after_snow, after_wet, wet1day, wet2day, &
       grm_nfert,grm_pfert,BNF_clover_daily, &
       GRM_devstage,agri_fert_save,agri_nfert,&
       agri_pfert, rest_Ngrazing, grazed_above,&
       days_senescence)
    
  !! 0. Variable and parameter declaration

    !! 0.1 input

    INTEGER, INTENT(in)                                 :: npts                 !! Domain size (unitless)
    REAL, INTENT(in)                                    :: dt_days              !! Time step of Stomate (days)
    INTEGER, DIMENSION(npts,NbNeighb), INTENT(in)       :: neighbours           !! Indices of the 8 neighbours of each grid 
                                                                                       !! point [1=North and then clockwise] 
    REAL, DIMENSION(npts,2), INTENT(in)                 :: resolution           !! Resolution at each grid point (m)  
                                                                                       !! [1=E-W, 2=N-S] 
    REAL, DIMENSION(npts,nvm), INTENT(in)               :: herbivores           !! Time constant of probability of a leaf to 
                                                                                       !! be eaten by a herbivore (days) 
    INTEGER, DIMENSION(npts),INTENT(in)                 :: soil_orders       !! dominant USDA soil order of the grid box
    REAL,DIMENSION (npts)    , INTENT(in)               :: max_eau_var       !! Maximum water content of the soil   
                                                                                            !! @tex ($kg m^{-2}$) @endtex 
    REAL,DIMENSION(npts),INTENT(in)                     :: bulk              !! Bulk density (kg/m**3) 


    REAL, DIMENSION(npts), INTENT(in)                   :: tsurf_daily          !! Daily surface temperatures (K)
    REAL, DIMENSION(npts,nslm), INTENT(in)              :: tsoil_daily          !! Daily soil temperatures (K)
    REAL, DIMENSION(npts), INTENT(in)                   :: t2m_daily            !! Daily 2 meter temperatures (K)
    REAL, DIMENSION(npts), INTENT(in)                   :: t2m_min_daily        !! Daily minimum 2 meter temperatures (K)
    REAL, DIMENSION(npts), INTENT(in)                   :: litterhum_daily      !! Daily litter humidity (0 to 1, unitless)
    REAL, DIMENSION(npts,nslm), INTENT(in)              :: soilhum_daily        !! Daily soil humidity (0 to 1, unitless)
    REAL, DIMENSION(npts,nvm), INTENT(in)               :: maxmoiavail_lastyear !! Last year's maximum moisture availability 
                                                                                       !! (0 to 1, unitless) 
    REAL, DIMENSION(npts,nvm), INTENT(in)               :: minmoiavail_lastyear !! Last year's minimum moisture availability 
                                                                                       !! (0 to 1, unitless) 
    REAL, DIMENSION(npts), INTENT(in)                   :: gdd0_lastyear        !! Last year's GDD0 (K)
    REAL, DIMENSION(npts), INTENT(in)                   :: precip_lastyear      !! Lastyear's precipitation 
                                                                                       !! @tex $(mm year^{-1})$ @endtex
                                                                                       !! to determine if establishment possible
    REAL, DIMENSION(npts,nvm), INTENT(in)               :: moiavail_month       !! "Monthly" moisture availability (0 to 1, 
                                                                                       !! unitless) 
    REAL, DIMENSION(npts,nvm), INTENT(in)               :: moiavail_week        !! "Weekly" moisture availability 
                                                                                       !! (0 to 1, unitless)
    REAL, DIMENSION(npts), INTENT(in)                   :: t2m_longterm         !! "Long term" 2 meter reference 
                                                                                       !! temperatures (K) 
    REAL, DIMENSION(npts), INTENT(in)                   :: t2m_month            !! "Monthly" 2-meter temperatures (K)
    REAL, DIMENSION(npts), INTENT(in)                   :: t2m_week             !! "Weekly" 2-meter temperatures (K)
    REAL, DIMENSION(npts), INTENT(in)                   :: Tseason              !! "seasonal" 2-meter temperatures (K)
    REAL, DIMENSION(npts,nvm), INTENT(in)               :: Tmin_spring_time     !! Number of days after begin_leaves (leaf onset) 
    REAL, DIMENSION(npts,nvm), INTENT(in)               :: onset_date           !! Date in the year at when the leaves started to grow(begin_leaves)
    REAL, DIMENSION(npts,nslm), INTENT(in)              :: tsoil_month          !! "Monthly" soil temperatures (K)
    REAL, DIMENSION(npts,nslm), INTENT(in)              :: soilhum_month        !! "Monthly" soil humidity
                                                                                       !! (0 to 1, unitless) 
    REAL, DIMENSION(npts,nvm), INTENT(in)               :: gdd_m5_dormance      !! Growing degree days (K), threshold -5 deg 
                                                                                       !! C (for phenology) 
    REAL, DIMENSION(npts,nvm), INTENT(in)               :: gdd_from_growthinit  !! growing degree days, since growthinit for crops
    REAL, DIMENSION(npts,nvm), INTENT(inout)            :: gdd_midwinter        !! Growing degree days (K), since midwinter 
                                                                                       !! (for phenology) - this is written to the history files 
    REAL, DIMENSION(npts,nvm), INTENT(in)               :: ncd_dormance         !! Number of chilling days (days), since 
                                                                                       !! leaves were lost (for phenology) 
    REAL, DIMENSION(npts,nvm), INTENT(in)               :: ngd_minus5           !! Number of growing days (days), threshold 
                                                                                       !! -5 deg C (for phenology) 
    REAL, DIMENSION(npts,nvm,nparts), INTENT(in)        :: turnover_longterm    !! "Long term" turnover rate  
                                                                                       !! @tex $(gC m^{-2} year^{-1})$ @endtex 
    REAL, DIMENSION(npts,nvm), INTENT(inout)            :: gpp_daily            !! Daily gross primary productivity  
                                                                                       !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL, DIMENSION(npts,nvm), INTENT(in)               :: gpp_week             !! Weekly gross primary productivity   
                                                                                       !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL, DIMENSION(npts,nvm), INTENT(in)               :: time_hum_min         !! Time elapsed since strongest moisture 
                                                                                       !! availability (days) 
    REAL, DIMENSION(npts,nvm), INTENT(in)               :: maxfpc_lastyear      !! Last year's maximum foliage projected
                                                                                       !! coverage for each natural PFT,
                                                                                       !! @tex $(m^2 m^{-2})$ @endtex 
    REAL, DIMENSION(npts,nvm,nparts), INTENT(in)        :: resp_maint_part      !! Maintenance respiration of different 
                                                                                       !! plant parts  
                                                                                       !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL, DIMENSION(npts,nvm), INTENT(in)               :: fpc_max              !! "Maximal" coverage fraction of a PFT (LAI 
                                                                                       !! -> infinity) on ground  
                                                                                       !! @tex $(m^2 m^{-2})$ @endtex 
    REAL, DIMENSION(npts,nvm),INTENT(in)                :: veget_cov_max_new    !! New "maximal" coverage fraction of a PFT 
    REAL, DIMENSION(npts,nvm), INTENT(inout)            :: cn_leaf_avg_season   !! Seasonal average CN ratio of leaves 
    REAL, DIMENSION(npts,nvm), INTENT(inout)            :: np_leaf_avg_season   !! Seasonal average NP ratio of leaves 
    REAL, DIMENSION(npts,nvm), INTENT(in)               :: nstress_season       !! N-related seasonal stress (used for allocation)    
    REAL, DIMENSION(npts,nvm), INTENT(in)               :: pstress_season       !! P-related seasonal stress (used for allocation)    
    REAL, DIMENSION(npts,nvm), INTENT(in)               :: moiavail_growingseason !! mean growing season moisture availability 

     REAL, DIMENSION(npts),INTENT(in)                   :: woodharvest          !! Harvested wood biomass (gC m-2 yr-1)

    REAL, DIMENSION(npts,nvm), INTENT(in)               :: lai_target_longterm  !! seasonal lai target 
    REAL, DIMENSION(npts,nvm), INTENT(out)              :: lai_target           !!  lai target 

  !! 0.2 Output variables
    
    REAL, DIMENSION(npts,nvm), INTENT(out)              :: npp_daily            !! Net primary productivity 
                                                                                       !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL, DIMENSION(npts,nvm,nparts,nelements), INTENT(out) :: turnover_daily       !! Turnover rates 
                                                                                       !! @tex $(gC m^{-2} year^{-1})$ @endtex 
    REAL, DIMENSION(npts,nvm), INTENT(out)              :: co2_to_bm            !! CO2 taken up from atmosphere when 
                                                                                       !! introducing a new PFT (introduced for 
                                                                                       !! carbon balance closure) 
                                                                                       !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL, DIMENSION(npts,nvm), INTENT(out)              :: n_to_bm              !! N taken up from ?? when 
                                                                                       !! introducing a new PFT (introduced for 
                                                                                       !! carbon balance closure) 
                                                                                       !! @tex $(gN m^{-2} dtslow^{-1})$ @endtex 
    REAL, DIMENSION(npts,nvm), INTENT(out)              :: p_to_bm              !! P taken up from ?? when 
                                                                                       !! introducing a new PFT (introduced for 
                                                                                       !! carbon balance closure) 
                                                                                       !! @tex $(gP m^{-2} dtslow^{-1})$ @endtex 
    REAL, DIMENSION(npts,nvm), INTENT(out)              :: co2_fire             !! Carbon emitted into the atmosphere by 
                                                                                       !! fire (living and dead biomass)  
                                                                                       !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL, DIMENSION(npts,nvm), INTENT(inout)            :: resp_hetero          !! Heterotrophic respiration
                                                                                       !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL, DIMENSION(npts,nvm), INTENT(out)              :: resp_maint           !! Maintenance respiration  
                                                                                       !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL, DIMENSION(npts,nvm), INTENT(out)              :: resp_growth          !! Growth respiration  
                                                                                       !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    
    REAL, DIMENSION(npts), INTENT(inout)                :: deadleaf_cover       !! Fraction of soil covered by dead leaves 
                                                                                       !! (0 to 1, unitless) 
    REAL, DIMENSION(npts,nvm,2), INTENT(out)              :: vcmax                !! Maximum rate of carboxylation ; third dimension 1=actual, 2=potential
    REAL, DIMENSION(npts,nvm,2), INTENT(out)              :: jmax                 !! maximum rate of electron transport; third dimension 1=actual, 2=potential

    REAL, DIMENSION(npts,nvm,nparts,nelements), INTENT(out):: bm_to_litter      !! Conversion of biomass to litter 
                                                                                       !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    LOGICAL, DIMENSION(npts,nvm), INTENT(out)                  :: begin_leaves         !! signal to start putting leaves on (true/false)
    REAL,DIMENSION(npts,nvm), INTENT(out)               :: nue                  !! Nitrogen use Efficiency with impact of leaf age (umol CO2 (gN)-1 s-1) 

    !! 0.3 Modified variables
    
    REAL, DIMENSION(npts,nvm), INTENT(inout)            :: height               !! Height of vegetation (m) 
    REAL, DIMENSION(npts,nlevs), INTENT(inout)          :: control_moist        !! Moisture control of heterotrophic 
                                                                                       !! respiration (0 to 1, unitless) 
    REAL, DIMENSION(npts,nlevs), INTENT(inout)          :: control_temp         !! Temperature control of heterotrophic 
                                                                                       !! respiration, above and below 
                                                                                       !! (0 to 1, unitless) 
    REAL, DIMENSION(npts,ncarb,nvm,nelements), INTENT(inout) :: som_input       !! Quantity of carbon going into carbon 
                                                                                       !! pools from litter decomposition  
                                                                                       !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL, DIMENSION(npts,nvm), INTENT(inout)            :: lai                  !! Leaf area index OF AN INDIVIDUAL PLANT,
                                                                                       !! where a PFT contains n indentical plants
                                                                                       !! i.e., using the mean individual approach 
                                                                                       !! @tex $(m^2 m^{-2})$ @endtex 
    REAL, DIMENSION(npts,nvm), INTENT(inout)            :: rprof                !! Prescribed root depth (m) 
    LOGICAL, DIMENSION(npts,nvm), INTENT(inout)                :: PFTpresent           !! Tab indicating which PFTs are present in 
                                                                                       !! each pixel 
    REAL, DIMENSION(npts,nvm), INTENT(inout)            :: age                  !! Age (years)    
    REAL, DIMENSION(npts,nvm), INTENT(inout)            :: fireindex            !! Probability of fire (0 to 1, unitless)
    REAL, DIMENSION(npts,nvm), INTENT(inout)            :: firelitter           !! Longer term litter above the ground that 
                                                                                       !! can be burned, @tex $(gC m^{-2})$ @endtex 
    REAL, DIMENSION(npts,nvm,nparts,nleafages), INTENT(inout)  :: leaf_age      !! Leaf age (days)
    REAL, DIMENSION(npts,nvm,nparts,nleafages), INTENT(inout)  :: leaf_frac     !! Fraction of leaves in leaf age class, 
                                                                                       !! (0 to 1, unitless)
    REAL, DIMENSION(npts,nvm,nparts,nelements), INTENT(inout) :: biomass        !! Biomass @tex $(gC m^{-2})$ @endtex
    REAL, DIMENSION(npts,nvm), INTENT(inout)            :: ind                  !! Density of individuals 
                                                                                       !! @tex $(m^{-2})$ @endtex 
    REAL, DIMENSION(npts,nvm), INTENT(inout)            :: adapted              !! Adaptation of PFT (killed if too cold) 
                                                                                       !! (0 to 1, unitless) 
    REAL, DIMENSION(npts,nvm), INTENT(inout)            :: regenerate           !! "Fitness": Winter sufficiently cold for 
                                                                                       !! PFT regeneration ? (0 to 1, unitless) 
    LOGICAL, DIMENSION(npts,nvm), INTENT(inout)                :: senescence           !! Flag for setting senescence stage (only 
                                                                                       !! for deciduous trees) 
    REAL, DIMENSION(npts,nvm), INTENT(inout)            :: when_growthinit      !! How many days ago was the beginning of 
                                                                                       !! the growing season (days) 
    REAL, DIMENSION(npts,nlitt,nvm,nlevs,nelements), INTENT(inout):: litter     !! Metabolic and structural litter, above 
                                                                                       !! and below ground 
                                                                                       !! @tex $(gC m^{-2})$ @endtex 
    REAL, DIMENSION(npts,nvm,nlitt), INTENT(inout)      :: dead_leaves          !! Dead leaves on ground, per PFT, metabolic 
                                                                                       !! and structural,  
                                                                                       !! @tex $(gC m^{-2})$ @endtex 
    REAL, DIMENSION(npts,ncarb,nvm,nelements), INTENT(inout)      :: som        !! Carbon pool: active, slow, or passive, 
                                                                                       !! @tex $(gC m^{-2})$ @endtex  
    REAL, DIMENSION(npts,nvm,nlevs), INTENT(inout)      :: lignin_struc         !! Ratio of Lignin/Carbon in structural 
                                                                                       !! litter, above and below ground,  
                                                                                       !! @tex $(gC m^{-2})$ @endtex 
    REAL, DIMENSION(npts,nvm,nlevs), INTENT(inout)      :: lignin_wood          !! Ratio of Lignin/Carbon in woody
                                                                                       !! @tex $(gC m^{-2})$ @endtex 
    REAL, DIMENSION(npts,nvm), INTENT(inout)            :: veget_cov_max        !! "Maximal" coverage fraction of a PFT (LAI 
                                                                                       !! -> infinity) on ground 
    REAL,DIMENSION(npts,nvm), INTENT(in)                :: veget_cov            !! Fractional coverage: actually share of the pixel 
                                                                                       !! covered by a PFT (fraction of ground area), 
                                                                                       !! taking into account LAI ??(= grid scale fpc)?? 
    REAL, DIMENSION(npts,nvm), INTENT(inout)            :: npp_longterm         !! "Long term" mean yearly primary 
                                                                                       !! productivity 
                                                                                       !! @tex $(gC m^{-2} year^{-1})$ @endtex 
    REAL, DIMENSION(npts,nvm), INTENT(inout)            :: lm_lastyearmax       !! Last year's maximum leaf mass, for each 
                                                                                       !! PFT @tex $(gC m^{-2})$ @endtex 
    REAL, DIMENSION(npts,nvm), INTENT(inout)            :: veget_lastlight      !! Vegetation fractions (on ground) after 
                                                                                       !! last light competition  
                                                                                       !! @tex $(m^2 m^{-2})$ @endtex 
    REAL, DIMENSION(npts,nvm), INTENT(inout)            :: everywhere           !! Is the PFT everywhere in the grid box or 
                                                                                       !! very localized (after its introduction) 
                                                                                       !! (unitless) 
    LOGICAL, DIMENSION(npts,nvm), INTENT(inout)                :: need_adjacent        !! In order for this PFT to be introduced, 
                                                                                       !! does it have to be present in an 
                                                                                       !! adjacent grid box? 
    REAL, DIMENSION(npts,nvm), INTENT(inout)            :: RIP_time             !! How much time ago was the PFT eliminated 
                                                                                       !! for the last time (y) 
    REAL, DIMENSION(npts,nvm), INTENT(inout)            :: turnover_time        !! Turnover_time of leaves for grasses 
                                                                                       !! (days)
    REAL,DIMENSION(npts,0:10), INTENT(inout)            :: prod10               !! Products remaining in the 10
                                                                                       !! year-turnover pool after the annual 
                                                                                       !! release for each compartment (10
                                                                                       !! + 1 : input from year of land cover 
                                                                                       !! change) @tex $(gC m^{-2})$ @endtex 
    REAL,DIMENSION(npts,0:100), INTENT(inout)           :: prod100              !! Products remaining in the 100 
                                                                                       !! year-turnover pool after the annual 
                                                                                       !! release for each compartment (100 
                                                                                       !! + 1 : input from year of land cover 
                                                                                       !! change) @tex $(gC m^{-2})$ @endtex 
    REAL,DIMENSION(npts,10), INTENT(inout)              :: flux10               !! Annual release from the 10
                                                                                       !! year-turnover pool compartments  
                                                                                       !! @tex $(gC m^{-2})$ @endtex 
    REAL,DIMENSION(npts,100), INTENT(inout)             :: flux100              !! Annual release from the 100 
                                                                                       !! year-turnover pool compartments  
                                                                                       !! @tex $(gC m^{-2})$ @endtex 
    REAL,DIMENSION(npts), INTENT(inout)                 :: convflux             !! Release during first year following land 
                                                                                       !! cover change @tex $(gC m^{-2})$ @endtex 
    REAL,DIMENSION(npts), INTENT(inout)                 :: cflux_prod10         !! Total annual release from the 10 
                                                                                       !! year-turnover pool 
                                                                                       !! @tex $(gC m^{-2})$ @endtex 
    REAL,DIMENSION(npts), INTENT(inout)                 :: cflux_prod100        !! Total annual release from the 100 
                                                                                       !! year-turnover pool 
                                                                                       !! @tex $(gC m^{-2})$ @endtex 
    REAL,DIMENSION(npts,0:10), INTENT(inout)            :: prod10_harvest       !! Products remaining in the 10
                                                                                       !! year-turnover pool after the annual 
                                                                                       !! release for each compartment (10
                                                                                       !! + 1 : input from year of wood harvest)
                                                                                       !! @tex $(gC m^{-2})$ @endtex 
    REAL,DIMENSION(npts,0:100), INTENT(inout)           :: prod100_harvest      !! Products remaining in the 100 
                                                                                       !! year-turnover pool after the annual 
                                                                                       !! release for each compartment (100 
                                                                                       !! + 1 : input from year of wood harvest)
                                                                                       !! @tex $(gC m^{-2})$ @endtex 
    REAL,DIMENSION(npts,10), INTENT(inout)              :: flux10_harvest       !! Annual release from the 10
                                                                                       !! year-turnover pool compartments  
                                                                                       !! @tex $(gC m^{-2})$ @endtex 
    REAL,DIMENSION(npts,100), INTENT(inout)             :: flux100_harvest      !! Annual release from the 100 
                                                                                       !! year-turnover pool compartments  
                                                                                       !! @tex $(gC m^{-2})$ @endtex 
    REAL,DIMENSION(npts), INTENT(inout)                 :: convflux_harvest     !! Release during first year following wood 
                                                                                       !! harvest @tex $(gC m^{-2})$ @endtex 
    REAL,DIMENSION(npts), INTENT(inout)                 :: cflux_prod10_harvest !! Total annual release from the 10 
                                                                                       !! year-turnover pool 
                                                                                       !! @tex $(gC m^{-2})$ @endtex 
    REAL,DIMENSION(npts), INTENT(inout)                 :: cflux_prod100_harvest!! Total annual release from the 100 
                                                                                       !! year-turnover pool 
                                                                                       !! @tex $(gC m^{-2})$ @endtex 
    REAL,DIMENSION(npts,nvm), INTENT(inout)             :: woodharvestpft       !! Harvested wood biomass (gC m-2 dt_stomate-1)
    REAL, DIMENSION(npts), INTENT(inout)                :: harvest_above        !! Harvest above ground biomass for 
                                                                                       !! agriculture @tex $(gC m^{-2})$ @endtex 
   REAL, DIMENSION(npts, nvm), INTENT(inout)           :: harvest_bioenergy    !! Harvest biomass for bioenergy
                                                                                       !! @tex $(gC m^{-2})$ @endtex
    REAL, DIMENSION(npts), INTENT(inout)                :: carb_mass_total      !! Carbon Mass total (soil, litter, veg) 
                                                                                       !! @tex $(gC m^{-2})$ @endtex  
    REAL, DIMENSION(npts,nvm,nbpools,nbpools), INTENT(inout) :: MatrixA         !! Matrix containing the fluxes  
                                                                                       !! between the carbon pools
                                                                                       !! per sechiba time step 
                                                                                       !! @tex $(gC.m^2.day^{-1})$ @endtex
    REAL, DIMENSION(:,:), INTENT(inout)                 :: KF                   !! Scaling factor to convert sapwood mass
                                                                                       !! into leaf mass (m)
    REAL, DIMENSION(:,:), INTENT(inout)                 :: k_latosa_adapt       !! Leaf to sapwood area adapted for long 
                                                                                       !! term water stress (m)
    REAL, DIMENSION(npts,nvm,nionspec), INTENT(inout)   :: n_uptake_daily       !! Uptake of soil N by plants  
    REAL, DIMENSION(npts,nvm), INTENT(inout)            :: p_uptake_daily       !! Uptake of soil P by plants  
       !! (gN/m**2/day)  
    REAL, DIMENSION(npts,nvm,nnspec), INTENT(inout)     :: soil_n_min           !! mineral nitrogen in the soil (gN/m**2)  
    REAL, DIMENSION(npts,nvm,npspec), INTENT(inout)     :: soil_p_min           !! mineral phosphorus in the soil (gP/m**2)  
    REAL, DIMENSION(:,:), INTENT(inout)                 :: rue_longterm         !! Longterm radiation use efficiency 
    REAL, DIMENSION(npts,nvm), INTENT(inout)            :: record_reserve_grass !! record maximum grass reserve pool
                                                                                       !! (??units??) 
    REAL, DIMENSION(npts,nlitt,nvm), INTENT(out)        :: litter_avail         !! edible litter for animals
    REAL, DIMENSION(npts,nlitt,nvm) , INTENT(out)       :: litter_not_avail     !! litter not edible for animals
    REAL, DIMENSION(npts,nlitt,nvm), INTENT(in)         :: litter_avail_frac    !! edible litter fraction
    REAL, DIMENSION(npts,nvm), INTENT(inout)            :: sla_calc             !! leaf age-related SLA
    REAL, DIMENSION(npts), INTENT(in)                   :: t2m_14               !! "14days" 2 meter temperatures (K)
    REAL, DIMENSION(npts,nvm), INTENT(inout)            :: wshtotsum            !! accumulated harvested grass biomass
    REAL, DIMENSION(npts,nvm), INTENT(inout)            :: sr_ugb               !! Optimised stocking density (animal m-2)
    REAL, DIMENSION(npts,nvm), INTENT(inout)            :: compt_ugb            !! Counter of grazing days (d)
    REAL, DIMENSION(npts,nvm), INTENT(inout)            :: nb_ani               !! optimal livestock density
    REAL, DIMENSION(npts,nvm), INTENT(inout)            :: grazed_frac          !! optimal fraction for grazing (in contract to mowning)
    REAL, DIMENSION(npts,nvm), INTENT(inout)            :: import_yield         !! annual total harvested grass biomass yield
    REAL, DIMENSION(npts,nvm), INTENT(inout)            :: sla_age1             !! maximum SLA for youngest leaves
    INTEGER, INTENT(in)                                 :: day_of_year          !! DOY
    REAL, DIMENSION(npts,nvm), INTENT(inout)            :: when_growthinit_cut  !! growing days after grass harvest
    REAL, DIMENSION(npts,nvm), INTENT(inout)            :: nb_grazingdays       !! annual total number of days animal grazing in field
    REAL, DIMENSION(npts), INTENT(in)                   :: snowfall_daily       !! daily snow fall
    REAL, DIMENSION(npts), INTENT(in)                   :: snowmass_daily       !! daily snow mass

    REAL, DIMENSION(npts,nvm), INTENT(in)               :: moiavail_daily       !! Daily moisture availability (0-1, unitless)

    REAL,DIMENSION (npts), INTENT(in)                   :: tmc_topgrass_daily   !! Daily top 5 layer soil moisture (m^3 m^-3)
    REAL,DIMENSION (npts), INTENT(in)                   :: fc_grazing           !! soil moisture threshold for detecting wet soil
    REAL,DIMENSION (npts), INTENT(inout)                :: after_snow           !! day counter after snow melt to prevent grazing
    REAL,DIMENSION (npts), INTENT(inout)                :: after_wet            !! day counter after wet soil is detected to prevent grazing
    REAL,DIMENSION (npts), INTENT(inout)                :: wet1day              !! accumulated days with wet soil
    REAL,DIMENSION (npts), INTENT(inout)                :: wet2day              !! accumulated days with wet soil
    REAL, DIMENSION(npts,nvm,ninput), INTENT(out)       :: grm_nfert            !! N fertilizer applied to grassland in GRM module
    REAL, DIMENSION(npts,nvm), INTENT(out)              :: grm_pfert            !! P fertilizer applied to grassland in GRM module
    REAL, DIMENSION(npts,nvm), INTENT(in)               :: BNF_clover_daily     !! daily clover BNF in grass-legume mixtures
    REAL, DIMENSION(npts,nvm), INTENT(inout)            :: GRM_devstage         !! development stage of grasses
    REAL, DIMENSION(npts,nvm,5), INTENT(in)             :: agri_fert_save       !! saved fertilization map from agricultural practices
    REAL, DIMENSION(npts,nvm,ninput), INTENT(out)       :: agri_nfert           !! N fertilization by agricultural practices
    REAL, DIMENSION(npts,nvm), INTENT(out)              :: agri_pfert           !! P fertilization by agricultural practices
    REAL,DIMENSION (npts), INTENT(inout)                :: grazed_above         !! total grazed biomass in simplified grazing module
    REAL, DIMENSION(npts,nvm), INTENT(inout)            :: rest_Ngrazing        !! the residual N applied to pasture/rangeland that did not output from biomass grazed
    REAL, DIMENSION(npts,nvm), INTENT(inout)            :: days_senescence      !! day counter after senescence is detected

    !! 0.4 Local variables

    REAL, DIMENSION(npts,nvm,nparts,nelements)           :: recycling_daily     !! Recycled nutrients from tissue prior shedding

    REAL, DIMENSION(npts,nvm,nelements)                  :: tot_bm_to_litter    !! Total conversion of biomass to litter 
                                                                                       !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL, DIMENSION(npts,nvm,nelements)                  :: tot_live_biomass    !! Total living biomass  
                                                                                       !! @tex $(gC m{-2})$ @endtex 
    REAL, DIMENSION(npts,nvm,nparts,nelements)           :: bm_alloc            !! Biomass increase, i.e. NPP per plant part 
                                                                                       !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL, DIMENSION(npts,nvm,nelements)                  :: tot_turnover        !! Total turnover rate  
                                                                                       !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL, DIMENSION(npts,nvm)                            :: tot_litter_soil_carb!! Total soil and litter carbon  
                                                                                       !! @tex $(gC m^{-2})$ @endtex 
    REAL, DIMENSION(npts,nvm)                            :: tot_litter_carb     !! Total litter carbon 
                                                                                       !! @tex $(gC m^{-2})$ @endtex 
    REAL, DIMENSION(npts,nvm)                            :: tot_soil_carb       !! Total soil carbon  
                                                                                       !! @tex $(gC m^{-2})$ @endtex 
    REAL, DIMENSION(npts,nvm)                            :: tot_litter_nitr     !! Total litter nitrogen
                                                                                       !! @tex $(gN m^{-2})$ @endtex
    REAL, DIMENSION(npts,nvm)                            :: tot_soil_nitr       !! Total soil nitrogen
                                                                                       !! @tex $(gN m^{-2})$ @endtex
    REAL, DIMENSION(npts,nvm)                            :: tot_litter_phos     !! Total litter phosphorus
                                                                                       !! @tex $(gP m^{-2})$ @endtex
    REAL, DIMENSION(npts,nvm)                            :: tot_soil_phos       !! Total soil phosphorus
                                                                                       !! @tex $(gP m^{-2})$ @endtex
    REAL, DIMENSION(npts)                                :: carb_mass_variation !! Carbon Mass variation  
                                                                                       !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL, DIMENSION(npts,nvm)                            :: cn_ind              !! Crown area of individuals 
                                                                                       !! @tex $(m^{2})$ @endtex 
    REAL, DIMENSION(npts,nvm)                            :: woodmass_ind        !! Woodmass of individuals (gC) 
    REAL, DIMENSION(npts,nvm,nparts)                     :: f_alloc             !! Fraction that goes into plant part 
                                                                                       !! (0 to 1, unitless) 
    REAL, DIMENSION(npts)                                :: avail_tree          !! Space availability for trees 
                                                                                       !! (0 to 1, unitless) 
    REAL, DIMENSION(npts)                                :: avail_grass         !! Space availability for grasses 
                                                                                       !! (0 to 1, unitless) 
    INTEGER                                                     :: j,k,l
    REAL,DIMENSION(npts)                                 :: prod10_total        !! Total products remaining in the pool 
                                                                                       !! after the annual release 
                                                                                       !! @tex $(gC m^{-2})$ @endtex 
    REAL,DIMENSION(npts)                                 :: prod100_total       !! Total products remaining in the pool 
                                                                                       !! after the annual release 
                                                                                       !! @tex $(gC m^{-2})$ @endtex 
    REAL,DIMENSION(npts)                                 :: cflux_prod_total    !! Total flux from conflux and the 10/100 
                                                                                       !! year-turnover pool 
                                                                                       !! @tex $(gC m^{-2} year^{-1})$ @endtex 
    REAL,DIMENSION(npts)                                 :: nflux_prod_total    !! Total flux from land cover change                                    
       !! @tex $(gN m^{-2} year^{-1})$ @endtex 
    REAL,DIMENSION(npts)                                 :: pflux_prod_total    !! Total flux from land cover change
       !! @tex $(gN m^{-2} year^{-1})$ @endtex 
    REAL,DIMENSION(npts)                                 :: prod10_harvest_total!! Total products remaining in the pool 
                                                                                       !! after the annual release 
                                                                                       !! @tex $(gC m^{-2})$ @endtex 
    REAL,DIMENSION(npts)                                 :: prod100_harvest_total!! Total products remaining in the pool 
                                                                                       !! after the annual release 
                                                                                       !! @tex $(gC m^{-2})$ @endtex 
    REAL,DIMENSION(npts)                                 :: cflux_prod_harvest_total!! Total flux from conflux and the 10/100 
                                                                                       !! year-turnover pool 
                                                                                       !! @tex $(gC m^{-2} year^{-1})$ @endtex 
    REAL,DIMENSION(npts,nvm)                             :: veget_cov_max_tmp   !! "Maximal" coverage fraction of a PFT  
                                                                                       !! (LAI-> infinity) on ground (unitless) 
    REAL, DIMENSION(npts,nvm)                            :: mortality           !! Fraction of individual dying this time 
                                                                                       !! step (0 to 1, unitless) 
    REAL, DIMENSION(npts)                                :: vartmp              !! Temporary variable used to add history
    REAL, DIMENSION(npts,nvm)                            :: histvar             !! History variables
    REAL, DIMENSION(npts,nvm)                            :: use_reserve         !! Mass taken from carbohydrate reserve 
                                                                                       !! @tex $(gC m^{-2})$ @endtex
    CHARACTER(LEN=2), DIMENSION(nelements)                     :: element_str          !! string suffix indicating element 
    REAL, DIMENSION(npts,nvm)                           :: vcmax_new            !! Vmax leafN relationship from Kattge et al. 2009 

    REAL, DIMENSION(npts,nvm)                           :: N_support            !! Nitrogen which is added to the ecosystem 
                                                                                       !! to support vegetation growth (only used when IMPOSE_CN .EQ. true) (gC/m2/day)
    REAL, DIMENSION(npts,nvm)                           :: P_support            !! Phosphorus which is added to the ecosystem 
                                                                                       !! to support vegetation growth (only used when IMPOSE_CN .EQ. true) (gP/m2/day)

    
    INTEGER                                                     ::    dsg1,dsg2,dsg3
    REAL, DIMENSION(npts,nvm,nelements)                :: mass_before
    REAL, DIMENSION(npts,nvm,nelements)                :: mass_before_2
    REAL, DIMENSION(npts,nvm,nelements)                :: mass_change

    REAL, PARAMETER :: sla_conv=.5                                              !! conversion of sla m2/g(C) ->m2/g(DW)

    !gmjc lcchange of managed grassland
    ! "maximal" coverage fraction of a PFT (LAI -> infinity) on ground
    INTEGER                                 :: ier
    LOGICAL                                        :: l_error =.FALSE.
    ! variables to get closure carbon cycle for nbp
    REAL, DIMENSION(npts,nvm)               :: harvest_gm         !! GRM harvested C
    REAL, DIMENSION(npts,nvm)               :: ranimal_gm         !! GRM animal respiration C
    REAL, DIMENSION(npts,nvm)               :: ch4_pft_gm         !! GRM PFT specific CH4-C emission
    REAL, DIMENSION(npts)                   :: ch4_gm             !! GRM CH4-C emission
    REAL, DIMENSION(npts,nvm)               :: cinput_gm          !! GRM C input through fertilization
    REAL, DIMENSION(npts)                   :: co2_gm             !! GRM CO2 taken
    REAL,DIMENSION(npts,nvm)                :: veget_cov_max_gm   !! GRM temporary vegetation cover
    REAL, DIMENSION(npts)                   :: veget_exist_gm     !! GRM temporaty vegetation existance
    REAL,DIMENSION(npts,nvm)                :: n2o_pft_gm         !! GRM PFT specific N2O-N emission
    REAL, DIMENSION(npts)                   :: n2o_gm             !! GRM N2O-N emission
    REAL,DIMENSION(npts,nvm,nelements)      :: grm_masschangecnp  !! GRM mass changes of CNP
    REAL,DIMENSION(npts,nvm,nelements)      :: fert_masschangecnp !! Fertilization mass changes of CNP
    REAL,DIMENSION(npts,nvm,nelements)      :: fire_masschangecnp !! FIRE mass changes of CNP
    REAL,DIMENSION(npts,nvm,nelements)      :: harvest_masschangecnp !! CROP harvest mass changes of CNP
    REAL,DIMENSION(npts,nvm,nelements)      :: harvest_BE_masschangecnp !! Bioenergy CROP harvest mass changes of CNP
    !end gmjc
    REAL, DIMENSION(npts,nvm,nelements)     :: mass_after        !! Temporary variable 
!_ ================================================================================================================================

    mass_after = zero
    ! init mass check variables
    mass_before(:,:,:) = zero
    mass_change(:,:,:) = zero
    grm_masschangecnp(:,:,:) = zero
    fert_masschangecnp(:,:,:) = zero
    fire_masschangecnp(:,:,:) = zero
    harvest_masschangecnp(:,:,:) = zero
    harvest_BE_masschangecnp(:,:,:) = zero
    IF (printlev>=3) WRITE(numout,*) 'Entering stomate_lpj'
 
  !! 1. Initializations
    
    !! 1.1 Initialize variables to zero
    co2_to_bm(:,:) = zero
    n_to_bm(:,:) = zero
    p_to_bm(:,:) = zero
    co2_fire(:,:) = zero
    npp_daily(:,:) = zero
    resp_maint(:,:) = zero
    resp_growth(:,:)      = zero
    harvest_above(:)      = zero
    harvest_bioenergy(:,:)= zero
    bm_to_litter(:,:,:,:) = zero
    cn_ind(:,:) = zero
    woodmass_ind(:,:) = zero
    turnover_daily(:,:,:,:) = zero
    recycling_daily(:,:,:,:) = zero
    use_reserve(:,:) = zero
    !! 1.2  Initialize variables to veget_cov_max
    veget_cov_max_tmp(:,:) = veget_cov_max(:,:)
    N_support(:,:) = zero
    P_support(:,:) = zero

    !! Initialize GRM variables for nbp to zero
    harvest_gm(:,:) = zero
    ranimal_gm(:,:) = zero
    ch4_pft_gm(:,:) = zero
    cinput_gm(:,:) = zero
    co2_gm(:) = zero
    ch4_gm(:) = zero
    n2o_gm(:) = zero
    n2o_pft_gm(:,:) = zero
    veget_cov_max_gm(:,:) = zero
    veget_exist_gm(:) = zero
    grazed_above(:) = zero

    !DSG mass conservation ========================================
    mass_before_2(:,:,:) = SUM(som(:,:,:,:),DIM=2) + SUM(SUM(litter(:,:,:,:,:),DIM=4),DIM=2) +  SUM(biomass(:,:,:,:),DIM=3)
    mass_before_2(:,:,iphosphorus) = mass_before_2(:,:,iphosphorus) +  SUM(soil_p_min(:,:,:),DIM=3)
    mass_before_2(:,:,initrogen)   = mass_before_2(:,:,initrogen)   +  SUM(soil_n_min(:,:,:),DIM=3)

    !! 1.3 Calculate some vegetation characteristics
    
    !! 1.3.1 Calculate some vegetation characteristics 
    !        Calculate cn_ind (individual crown mass) and individual height from
    !        state variables if running DGVM or dynamic mortality in static cover mode
    !??        Explain (maybe in the header once) why you mulitply with veget_cov_max in the DGVM
    !??        and why you don't multiply with veget_cov_max in stomate.
    IF ( ok_dgvm .OR. .NOT.lpj_gap_const_mort) THEN
       IF(ok_dgvm) THEN
          WHERE (ind(:,:).GT.min_stomate)
             woodmass_ind(:,:) = &
                  ((biomass(:,:,isapabove,icarbon)+biomass(:,:,isapbelow,icarbon) &
                  +biomass(:,:,iheartabove,icarbon)+biomass(:,:,iheartbelow,icarbon)) & 
                  *veget_cov_max(:,:))/ind(:,:)
          ENDWHERE
       ELSE
          WHERE (ind(:,:).GT.min_stomate)
             woodmass_ind(:,:) = &
                  (biomass(:,:,isapabove,icarbon)+biomass(:,:,isapbelow,icarbon) &
                  +biomass(:,:,iheartabove,icarbon)+biomass(:,:,iheartbelow,icarbon))/ind(:,:)
          ENDWHERE
       ENDIF

       ! DSG: does not modify biomass
       CALL crown (npts,  PFTpresent, &
            ind, biomass, woodmass_ind, &
            veget_cov_max, cn_ind, height)
    ENDIF

    !! 1.3.2 Prescribe characteristics if the vegetation is not dynamic
    !        IF the DGVM is not activated, the density of individuals and their crown
    !        areas don't matter, but they should be defined for the case we switch on
    !        the DGVM afterwards. At the first call, if the DGVM is not activated, 
    !        impose a minimum biomass for prescribed PFTs and declare them present.
    IF (printlev>=4) WRITE(numout,*) 'before prescribe'
    IF (printlev>=4) WRITE(numout,*) 'ind(test_grid,test_pft)=',ind(test_grid,test_pft)

    !DSG mass conservation ========================================
    mass_before(:,:,:) = SUM(biomass(:,:,:,:),DIM=3)

    CALL prescribe (npts, &
         veget_cov_max, dt_days, PFTpresent, everywhere, when_growthinit, &
         biomass, leaf_frac, ind, co2_to_bm,n_to_bm,p_to_bm, &
         KF,senescence,age,npp_longterm,&
         lm_lastyearmax,k_latosa_adapt)

    IF(dsg_debug) THEN
       !DSG mass conservation ============================================
       CALL  check_mass(npts,biomass(:,:,:,:),'lpj: after prescribe')

       mass_after(:,:,:) = SUM(biomass(:,:,:,:),DIM=3)
       mass_change(:,:,icarbon)     = co2_to_bm
       mass_change(:,:,initrogen)   = n_to_bm
       mass_change(:,:,iphosphorus) = p_to_bm

       CALL cons_mass( mass_before(:,:,:),     &  ! mass before
              mass_after,                      &  ! mass after
              mass_change(:,:,:),              &  ! net of fluxes
              'lpj: after prescribe' )
    ENDIF

    IF (printlev>=4) WRITE(numout,*) 'Leaving prescribe'
    IF (printlev>=4) WRITE(numout,*) 'ind(test_grid,test_pft)=',ind(test_grid,test_pft)
  !! 2. Climatic constraints for PFT presence and regenerativeness

    !   Call this even when DGVM is not activated so that "adapted" and "regenerate"
    !   are kept up to date for the moment when the DGVM is activated.
    CALL constraints (npts, dt_days, &
         t2m_month, t2m_min_daily,when_growthinit, Tseason, &
         adapted, regenerate)

    
  !! 3. Determine introduction and elimination of PTS based on climate criteria
 
    IF ( ok_dgvm ) THEN

       !DSG mass conservation ========================================
       mass_before(:,:,:) = SUM(biomass(:,:,:,:),DIM=3)
      
       !! 3.1 Calculate introduction and elimination
       CALL pftinout (npts, dt_days, adapted, regenerate, &
            neighbours, veget_cov_max, &
            biomass, ind, cn_ind, age, leaf_frac, npp_longterm, lm_lastyearmax, senescence, &
            PFTpresent, everywhere, when_growthinit, need_adjacent, RIP_time, &
            co2_to_bm, n_to_bm, p_to_bm, &
            avail_tree, avail_grass, &
            sla_calc)
 
       IF (dsg_debug) THEN
          CALL  check_mass(npts,biomass(:,:,:,:), 'lpj: after pftinout')
          !DSG mass conservation ============================================
          WRITE (numout,*) 'MASS CONSERVATON: pftinout is designed to violate against mass conservation'
          WRITE (numout,*) 'but it should be not detectable according to comments'
          mass_change(:,:,icarbon)     = zero
          mass_change(:,:,initrogen)   = zero
          mass_change(:,:,iphosphorus) = zero
          mass_after = SUM(biomass(:,:,:,:),DIM=3)

          CALL cons_mass( mass_before(:,:,:),              &  ! mass before
                 mass_after,                      &  ! mass after
                 mass_change(:,:,:),              &  ! net of fluxes
                 'lpj: after pftinout' )
       ENDIF

       !! 3.2 Reset attributes for eliminated PFTs.
       !     This also kills PFTs that had 0 leafmass during the last year. The message
       !     "... after pftinout" is misleading in this case.


       !DSG mass conservation ========================================
       mass_before(:,:,:) = SUM(biomass(:,:,:,:),DIM=3)

       CALL kill (npts, 'pftinout  ', lm_lastyearmax, &
            ind, PFTpresent, cn_ind, biomass, senescence, RIP_time, &
            lai, age, leaf_age, leaf_frac, npp_longterm, &
            when_growthinit, everywhere, veget_cov_max, bm_to_litter)

       IF(dsg_debug) THEN
          CALL  check_mass(npts,biomass(:,:,:,:), 'lpj: after kill')
          !DSG mass conservation ============================================
          mass_change(:,:,icarbon)     = SUM(bm_to_litter(:,:,:,icarbon),DIM=3)
          mass_change(:,:,initrogen)   = SUM(bm_to_litter(:,:,:,initrogen),DIM=3)
          mass_change(:,:,iphosphorus) = SUM(bm_to_litter(:,:,:,iphosphorus),DIM=3)
          mass_after = SUM(biomass(:,:,:,:),DIM=3)
          CALL cons_mass( mass_before(:,:,:),               &  ! mass before
                 mass_after,                       &  ! mass after
                 mass_change(:,:,:),               &  ! net of fluxes
                 'lpj: after kill' )
       ENDIF
       
       !! 3.3 Calculate woodmass of individual tree
       IF(ok_dgvm) THEN
          WHERE ((ind(:,:).GT.min_stomate))
             woodmass_ind(:,:) = &
                  ((biomass(:,:,isapabove,icarbon) + biomass(:,:,isapbelow,icarbon) &
                  + biomass(:,:,iheartabove,icarbon) + biomass(:,:,iheartbelow,icarbon))*veget_cov_max(:,:))/ind(:,:)
          ENDWHERE
       ELSE
          WHERE ((ind(:,:).GT.min_stomate))
             woodmass_ind(:,:) =(biomass(:,:,isapabove,icarbon) + biomass(:,:,isapbelow,icarbon) &
                  + biomass(:,:,iheartabove,icarbon) + biomass(:,:,iheartbelow,icarbon))/ind(:,:)
          ENDWHERE
       ENDIF
       
       ! Calculate crown area and diameter for all PFTs (including the newly established)
       CALL crown (npts, PFTpresent, &
            ind, biomass, woodmass_ind, &
            veget_cov_max, cn_ind, height)

    ENDIF
    
  !! 4. Phenology

    !! 4.1 Write values to history file
    !      Current values for ::when_growthinit 
    CALL xios_orchidee_send_field("WHEN_GROWTHINIT",when_growthinit)

    CALL histwrite_p (hist_id_stomate, 'WHEN_GROWTHINIT', itime, when_growthinit, npts*nvm, horipft_index)

    ! Set and write values for ::PFTpresent
    WHERE(PFTpresent)
       histvar=un
    ELSEWHERE
       histvar=zero
    ENDWHERE

    CALL xios_orchidee_send_field("PFTPRESENT",histvar)

    CALL histwrite_p (hist_id_stomate, 'PFTPRESENT', itime, histvar, npts*nvm, horipft_index)

    ! Set and write values for gdd_midwinter
    WHERE(gdd_midwinter.EQ.undef)
       histvar=val_exp
    ELSEWHERE
       histvar=gdd_midwinter
    ENDWHERE

    CALL xios_orchidee_send_field("GDD_MIDWINTER",histvar)

    CALL histwrite_p (hist_id_stomate, 'GDD_MIDWINTER', itime, histvar, npts*nvm, horipft_index)

    ! Set and write values for gdd_m5_dormance
    WHERE(gdd_m5_dormance.EQ.undef)
       histvar=val_exp
    ELSEWHERE
       histvar=gdd_m5_dormance
    ENDWHERE
    
    CALL xios_orchidee_send_field('GDD_M5_DORMANCE',histvar)
    CALL histwrite_p (hist_id_stomate, 'GDD_M5_DORMANCE', itime, histvar, npts*nvm, horipft_index)

    ! Set and write values for ncd_dormance
    WHERE(ncd_dormance.EQ.undef)
       histvar=val_exp
    ELSEWHERE
       histvar=ncd_dormance
    ENDWHERE

    CALL xios_orchidee_send_field("NCD_DORMANCE",histvar)

    CALL histwrite_p (hist_id_stomate, 'NCD_DORMANCE', itime, histvar, npts*nvm, horipft_index)

    !DSG mass conservation ========================================
    mass_before(:,:,:) = SUM(biomass(:,:,:,:),DIM=3)
    ! we might have already added co2,n or p to bm in prescribe thus we need to
    ! substract here again:
    mass_before(:,:,icarbon)     = mass_before(:,:,icarbon)     - co2_to_bm(:,:)
    mass_before(:,:,initrogen)   = mass_before(:,:,initrogen)   - n_to_bm(:,:)
    mass_before(:,:,iphosphorus) = mass_before(:,:,iphosphorus) - p_to_bm(:,:)

    IF (printlev>=3) THEN
       WRITE(numout,*) 'stomateLPJ: before phenology_prognostic'
       WRITE(numout,*) 'biomass(test_grid, test_pft,:,icarbon)',biomass(test_grid, test_pft,:,icarbon)
       WRITE(numout,*) 'biomass(test_grid, test_pft,:,initrogen)',biomass(test_grid, test_pft,:,initrogen)
       WRITE(numout,*) 'biomass(test_grid, test_pft,:,iphosphorus)',biomass(test_grid, test_pft,:,iphosphorus)

       WRITE(numout,*) 'co2_to_bm(test_grid, test_pft)',co2_to_bm(test_grid, test_pft)
       WRITE(numout,*) 'n_to_bm(test_grid, test_pft)',n_to_bm(test_grid, test_pft)
       WRITE(numout,*) 'p_to_bm(test_grid, test_pft)',p_to_bm(test_grid, test_pft)
    ENDIF
      

    !! 4.2 Calculate phenology
    CALL phenology (npts, dt_days, PFTpresent, &
         veget_cov_max, &
         t2m_longterm, t2m_month, t2m_week, gpp_daily, &
         maxmoiavail_lastyear, minmoiavail_lastyear, &
         moiavail_month, moiavail_week, &
         gdd_m5_dormance, gdd_midwinter, ncd_dormance, ngd_minus5, &
         senescence, time_hum_min, &
         biomass,som, leaf_frac, leaf_age, &
         when_growthinit, co2_to_bm, n_to_bm, p_to_bm, &
         begin_leaves,ind,KF,N_support, P_support, &
         cn_leaf_avg_season, np_leaf_avg_season, &
         sla_calc)
    !DSG: the variable senescene is not used in phenology
    
    IF (dsg_debug) THEN
       !DSG debug: ===================================================
       WRITE(6,*) 'location in code: '
       WRITE(6,*) 'stomateLPJ: after phenology_prognostic'

       IF (printlev>=3) THEN
          WRITE(numout,*) 'stomateLPJ: after phenology_prognostic'
          WRITE(numout,*) 'biomass(test_grid, test_pft,:,icarbon)',biomass(test_grid, test_pft,:,icarbon)
          WRITE(numout,*) 'biomass(test_grid, test_pft,:,initrogen)',biomass(test_grid, test_pft,:,initrogen)
          WRITE(numout,*) 'biomass(test_grid, test_pft,:,iphosphorus)',biomass(test_grid, test_pft,:,iphosphorus)

          WRITE(numout,*) 'co2_to_bm(test_grid, test_pft)',co2_to_bm(test_grid, test_pft)
          WRITE(numout,*) 'n_to_bm(test_grid, test_pft)',n_to_bm(test_grid, test_pft)
          WRITE(numout,*) 'p_to_bm(test_grid, test_pft)',p_to_bm(test_grid, test_pft)
       ENDIF
      
       CALL  check_mass(npts,biomass(:,:,:,:), 'lpj: after phenology_prognostic')

      !DSG mass conservation ============================================
      mass_after(:,:,:) = SUM(biomass(:,:,:,:),DIM=3)
      mass_change(:,:,icarbon)     = co2_to_bm(:,:)
      mass_change(:,:,initrogen)   = n_to_bm(:,:)
      mass_change(:,:,iphosphorus) = p_to_bm(:,:)

      CALL cons_mass( mass_before(:,:,:),               &  ! mass before
            mass_after,                       &  ! mass afte
            mass_change(:,:,:),               &  ! net of fluxes
            'lpj: after phenology'  )

!DSG: check being added       
! ======================================


       IF(printlev>=3)THEN
           WRITE(numout,*)    'biomass(test_grid,test_pft,:,icarbon)',biomass(test_grid,test_pft,:,icarbon)
           WRITE(numout,*)    'biomass(test_grid,test_pft,ilabile,initrogen)',biomass(test_grid,test_pft,ilabile,initrogen)
           WRITE(numout,*)    'gpp_daily(test_grid,test_pft)' ,gpp_daily(test_grid,test_pft)
           WRITE(numout,*)    'npp_daily(test_grid,test_pft)' ,npp_daily(test_grid,test_pft)
           WRITE(numout,*)    'resp_maint(test_grid,test_pft)',resp_maint(test_grid,test_pft)
           WRITE(numout,*)    'resp_growth(test_grid,test_pft)',resp_growth(test_grid,test_pft)
           WRITE(numout,*)    'dt_days',dt_days
        ENDIF

     !DSG mass conservation ========================================


     mass_before(:,:,:)           = SUM(biomass(:,:,:,:),DIM=3)
     mass_before(:,:,initrogen)   = mass_before(:,:,initrogen) - N_support
     mass_before(:,:,iphosphorus) = mass_before(:,:,iphosphorus) - P_support
    ENDIF
 
    ! JC MOD024 calculate GRM_devstage
    CALL calc_GRM_devstage(npts, t2m_daily, t2m_week, tsurf_daily, begin_leaves, GRM_devstage)
    ! End JC MOD024

    CALL growth_fun_all (npts, dt_days, veget_cov_max, veget_cov, PFTpresent, &
        senescence, when_growthinit, t2m_week, &
        nstress_season, pstress_season,moiavail_growingseason, &
        gpp_daily, resp_maint_part, resp_maint, &
        gpp_week, &
        resp_growth, npp_daily, biomass, age, &
        leaf_age, leaf_frac, use_reserve, &
        ind, rue_longterm, KF, k_latosa_adapt, &
        cn_leaf_avg_season, n_uptake_daily, N_support, &
        np_leaf_avg_season, p_uptake_daily, P_support, &
        lai_target, lai_target_longterm, &
        record_reserve_grass, &
        sla_calc, sla_age1, BNF_clover_daily, &
        GRM_devstage, when_growthinit_cut, days_senescence)

     IF(dsg_debug) THEN 
        !DSG debug: ===================================================
        IF(PRINTLEV>=3) THEN
             WRITE(numout,*) 'location in code: '
             WRITE(numout,*) 'stomateLPJ: after growth_fun'
             WRITE(numout,*)    'SUM(Biomass(test_grid,test_pft,:,icarbon))',  SUM(biomass(test_grid,test_pft,:,icarbon))
             WRITE(numout,*)    'SUM(Biomass(test_grid,test_pft,:,initrogen))',SUM(biomass(test_grid,test_pft,:,initrogen))
             WRITE(numout,*)    'biomass(test_grid,test_pft,:,icarbon)',biomass(test_grid,test_pft,:,icarbon)
             WRITE(numout,*)    'biomass(test_grid,test_pft,:,initrogen)',biomass(test_grid,test_pft,:,initrogen)
             WRITE(numout,*)    'gpp_daily(test_grid,test_pft)' ,gpp_daily(test_grid,test_pft)
             WRITE(numout,*)    'npp_daily(test_grid,test_pft)' ,npp_daily(test_grid,test_pft)
             WRITE(numout,*)    'resp_maint(test_grid,test_pft)',resp_maint(test_grid,test_pft)
             WRITE(numout,*)    'resp_growth(test_grid,test_pft)',resp_growth(test_grid,test_pft)
             WRITE(numout,*)    'SUM(n_uptake_daily(test_grid,test_pft,:))',SUM(n_uptake_daily(test_grid,test_pft,:))
             WRITE(numout,*)    'n_support(test_grid,test_pft)',N_support(test_grid,test_pft)
             WRITE(numout,*)    'BNF_clover_daily(test_grid,test_pft)',BNF_clover_daily(test_grid,test_pft)
             WRITE(numout,*)    'dt_days',dt_days
        ENDIF 
        CALL check_mass(npts,biomass(:,:,:,:), 'lpj: after growth_fun')

        !DSG mass conservation ============================================
        mass_change(:,:,icarbon)     = npp_daily(:,:)*dt_days
        IF (impose_cn) THEN
           mass_change(:,:,initrogen)   = SUM(n_uptake_daily(:,:,:),DIM=3) + N_support(:,:) + &
              BNF_clover_daily(:,:)
        ELSE
           mass_change(:,:,initrogen)   = SUM(n_uptake_daily(:,:,:),DIM=3) + &
               BNF_clover_daily(:,:)
        ENDIF

        IF (impose_cn) THEN
            mass_change(:,:,iphosphorus) = p_uptake_daily(:,:)  + P_support(:,:)
        ELSE
            mass_change(:,:,iphosphorus) = p_uptake_daily(:,:) 
        ENDIF

        mass_after = SUM(biomass(:,:,:,:),DIM=3)
        CALL cons_mass( mass_before(:,:,:),               &  ! mass before
              mass_after,                       &  ! mass afte
              mass_change(:,:,:),               &  ! net of fluxes
              'lpj: after growth_fun'  )
     ENDIF

  !! 6. NPP, maintenance and growth respiration

    !! 6.1 Calculate NPP and respiration terms
!    CALL npp_calc (npts, dt_days, &
!         PFTpresent, &
!         t2m_daily, tsoil_daily, lai, rprof, &
!         gpp_daily, f_alloc, bm_alloc, resp_maint_part,&
!         biomass, leaf_age, leaf_frac, age, &
!         resp_maint, resp_growth, npp_daily)

    !! 6.2 Kill slow growing PFTs in DGVM or STOMATE with constant mortality
    IF ( ok_dgvm .OR. .NOT.lpj_gap_const_mort) THEN

       CALL kill (npts, 'npp       ', lm_lastyearmax,  &
            ind, PFTpresent, cn_ind, biomass, senescence, RIP_time, &
            lai, age, leaf_age, leaf_frac, npp_longterm, &
            when_growthinit, everywhere, veget_cov_max, bm_to_litter)

       !! 6.2.1 Update wood biomass      
       !        For the DGVM
       IF(ok_dgvm) THEN
          WHERE (ind(:,:).GT.min_stomate)
             woodmass_ind(:,:) = &
                  ((biomass(:,:,isapabove,icarbon) + biomass(:,:,isapbelow,icarbon) &
                  + biomass(:,:,iheartabove,icarbon) + biomass(:,:,iheartbelow,icarbon)) & 
                  *veget_cov_max(:,:))/ind(:,:)
          ENDWHERE

       ! For all pixels with individuals
       ELSE
          WHERE (ind(:,:).GT.min_stomate)
             woodmass_ind(:,:) = &
                  (biomass(:,:,isapabove,icarbon) + biomass(:,:,isapbelow,icarbon) &
                  + biomass(:,:,iheartabove,icarbon) + biomass(:,:,iheartbelow,icarbon))/ind(:,:)
          ENDWHERE
       ENDIF ! ok_dgvm

       !! 6.2.2 New crown area and maximum vegetation cover after growth
       CALL crown (npts, PFTpresent, &
            ind, biomass, woodmass_ind,&
            veget_cov_max, cn_ind, height)

    ENDIF ! ok_dgvm

  !! 7. fire
    !DSG mass conservation ========================================
    mass_before(:,:,:) = SUM(biomass(:,:,:,:),DIM=3) + &
                               SUM(bm_to_litter(:,:,:,:),DIM=3) + &
                               SUM(SUM(litter(:,:,:,:,:),DIM=4),DIM=2)

    !! 7.1. Burn PFTs
    CALL fire (npts, dt_days, &
         litterhum_daily, t2m_daily, lignin_struc, lignin_wood, veget_cov_max, &
         fireindex, firelitter, biomass, ind, &
         litter, dead_leaves, bm_to_litter, &
         co2_fire, fire_masschangecnp, MatrixA)
    !gmjc
    ! after fire burning
    litter_avail(:,:,:) = litter(:,:,:,iabove,icarbon) * &
            litter_avail_frac(:,:,:)
    litter_not_avail(:,:,:) = litter(:,:,:,iabove,icarbon) * &
            (1.0 - litter_avail_frac(:,:,:))
    !end gmjc

    !! 7.2 Kill PFTs in DGVM
    IF ( ok_dgvm ) THEN
       ! reset attributes for eliminated PFTs
       CALL kill (npts, 'fire      ', lm_lastyearmax, &
            ind, PFTpresent, cn_ind, biomass, senescence, RIP_time, &
            lai, age, leaf_age, leaf_frac, npp_longterm, &
            when_growthinit, everywhere, veget_cov_max, bm_to_litter)
    ENDIF ! ok_dgvm

    IF (dsg_debug) THEN
      !DSG debug: ===================================================
      CALL  check_mass(npts,biomass(:,:,:,:), 'stomateLPJ: after fire')
      !DSG mass conservation ============================================
      ! JC NOTED
      ! the fire has included many fluxes to
      ! environment (such as emission of gaseous N and aerosol NP)
      mass_change(:,:,icarbon)     = -fire_masschangecnp(:,:,icarbon)
      mass_change(:,:,initrogen)   = -fire_masschangecnp(:,:,initrogen)
      mass_change(:,:,iphosphorus) = -fire_masschangecnp(:,:,iphosphorus)

      mass_after = SUM(biomass(:,:,:,:),DIM=3) + SUM(bm_to_litter(:,:,:,:),DIM=3) + &
            SUM(SUM(litter(:,:,:,:,:),DIM=4),DIM=2)

      CALL cons_mass( mass_before(:,:,:),     &  ! mass before
            mass_after(:,:,:),                & ! mass after
            mass_change(:,:,:),               &  ! net of fluxes
            'lpj: after fire')
    ENDIF

  !! 8. Tree mortality
    !DSG mass conservation ========================================
    mass_before(:,:,:) = SUM(biomass(:,:,:,:),DIM=3) + SUM(bm_to_litter(:,:,:,:),DIM=3)

    !IF (dsg_debug) THEN
    !    WRITE(6,*) 'location in code: '
    !    WRITE(6,*) 'stomateLPJ: before gap'
    !    IF(PRINTLEV>=3) THEN
    !    !    WRITE(6,*) 'bm_to_litter(test_grid,test_pft,:,icarbon)',bm_to_litter(test_grid,test_pft,:,icarbon)
    !        WRITE(6,*)'mass_before(test_grid,test_pft,icarbon)',mass_before(test_grid,test_pft,icarbon)
    !        WRITE(6,*)'mass_before(test_grid,test_pft,initrogen)',mass_before(test_grid,test_pft,initrogen)
    !    !    WRITE(6,*) 'bm_to_litter(test_grid,test_pft,:,initrogen)',bm_to_litter(test_grid,test_pft,:,initrogen)
    !    !    WRITE(6,*)'biomass(test_grid,test_pft,:,initrogen)',biomass(test_grid,test_pft,:,initrogen)
    !    !    WRITE(6,*) 'bm_to_litter(test_grid,test_pft,:,iphosphorus)',bm_to_litter(test_grid,test_pft,:,iphosphorus)
    !    !    WRITE(6,*)'biomass(test_grid,test_pft,:,iphosphorus',biomass(test_grid,test_pft,:,iphosphorus)
    !    ENDIF
    !ENDIF

    ! Does not depend on age, therefore does not change crown area.
    CALL gap (npts, dt_days, &
         npp_longterm, turnover_longterm, lm_lastyearmax, &
         PFTpresent, t2m_min_daily, Tmin_spring_time, &
         biomass, ind, bm_to_litter, mortality, &
         sla_calc)

    IF ( ok_dgvm ) THEN

       ! reset attributes for eliminated PFTs
       CALL kill (npts, 'gap       ', lm_lastyearmax, &
            ind, PFTpresent, cn_ind, biomass, senescence, RIP_time, &
            lai, age, leaf_age, leaf_frac, npp_longterm, &
            when_growthinit, everywhere, veget_cov_max, bm_to_litter)

    ENDIF

    IF(dsg_debug) THEN
       !DSG debug: ===================================================
       IF(PRINTLEV>=3) THEN
            WRITE(numout,*) 'location in code: '
            WRITE(numout,*) 'stomateLPJ: after gap n kill (1)'
            WRITE(numout,*) 'test_pft',test_pft
            WRITE(numout,*) 'test_grid',test_grid
            WRITE(numout,*) 'biomass(test_grid,test_pft,:,icarbon)',biomass(test_grid,test_pft,:,icarbon)
            WRITE(numout,*) 'biomass(test_grid,test_pft,:,initrogen)',biomass(test_grid,test_pft,:,initrogen)
            WRITE(numout,*) 'bm_to_litter(test_grid,test_pft,:,icarbon)',bm_to_litter(test_grid,test_pft,:,icarbon)
            WRITE(numout,*) 'bm_to_litter(test_grid,test_pft,:,initrogen)',bm_to_litter(test_grid,test_pft,:,initrogen)
       ENDIF
       CALL check_mass(npts,biomass(:,:,:,:), 'lpj: after gap n kill')

       !DSG mass conservation ============================================
       mass_change(:,:,icarbon)     = -SUM(bm_to_litter(:,:,:,icarbon),DIM=3)
       mass_change(:,:,initrogen)   = -SUM(bm_to_litter(:,:,:,initrogen),DIM=3)
       mass_change(:,:,iphosphorus) = -SUM(bm_to_litter(:,:,:,iphosphorus),DIM=3)

       mass_after = SUM(biomass(:,:,:,:),DIM=3)
       CALL cons_mass( mass_before(:,:,:),               &  ! mass before
             mass_after,                       &  ! mass after
             mass_change(:,:,:),               &  ! net of fluxes
             'lpj: after gap n kill')
     ENDIF


!! =====================================
!! DSG: mass conservation issue: START

  !! 10. Leaf senescence, new lai and other turnover processes
    !DSG mass conservation ========================================
    mass_before(:,:,:) = SUM(biomass(:,:,:,:),DIM=3)  + SUM(bm_to_litter(:,:,:,:),DIM=3)
       WRITE(6,*) 'location in code: '

    CALL turn (npts, dt_days, PFTpresent, &
         herbivores, &
         maxmoiavail_lastyear, minmoiavail_lastyear, &
         moiavail_week,  moiavail_month,t2m_longterm, t2m_month, t2m_week, veget_cov_max, &
         gdd_from_growthinit, leaf_age, leaf_frac, age, lai, biomass, &
         turnover_daily, recycling_daily, senescence,turnover_time, &
         sla_calc,GRM_enable_grazing,when_growthinit,days_senescence)

    !! 10. Light competition
    
    !! If not using constant mortality then kill with light competition
!    IF ( ok_dgvm .OR. .NOT.(lpj_gap_const_mort) ) THEN
    IF ( ok_dgvm ) THEN
 
       !! 10.1 Light competition
       CALL light (npts, dt_days, &
            veget_cov_max, fpc_max, PFTpresent, cn_ind, lai, maxfpc_lastyear, &
            lm_lastyearmax, ind, biomass, veget_lastlight, bm_to_litter, mortality, &
            sla_calc)
       
       !! 10.2 Reset attributes for eliminated PFTs
       CALL kill (npts, 'light     ', lm_lastyearmax, &
            ind, PFTpresent, cn_ind, biomass, senescence, RIP_time, &
            lai, age, leaf_age, leaf_frac, npp_longterm, &
            when_growthinit, everywhere, veget_cov_max, bm_to_litter)

    ENDIF
    IF(dsg_debug) THEN
       !DSG debug: ===================================================
       IF(PRINTLEV>=3) THEN
          WRITE(numout,*) 'location in code: '
          WRITE(numout,*) 'stomateLPJ: after turn'
          WRITE(numout,*) 'bm_to_litter(test_grid,test_pft,:,icarbon)',bm_to_litter(test_grid,test_pft,:,icarbon)
          WRITE(numout,*)'biomass(test_grid,test_pft,:,icarbon)',biomass(test_grid,test_pft,:,icarbon)
          WRITE(numout,*) 'bm_to_litter(test_grid,test_pft,:,initrogen)',bm_to_litter(test_grid,test_pft,:,initrogen)
          WRITE(numout,*)'biomass(test_grid,test_pft,:,initrogen)',biomass(test_grid,test_pft,:,initrogen)
          WRITE(numout,*) 'bm_to_litter(test_grid,test_pft,:,iphosphorus)',bm_to_litter(test_grid,test_pft,:,iphosphorus)
          WRITE(numout,*)'biomass(test_grid,test_pft,:,iphosphorus',biomass(test_grid,test_pft,:,iphosphorus)
       ENDIF
       CALL  check_mass(npts,biomass(:,:,:,:), 'lpj: after turn')

       !DSG mass conservation ============================================
       mass_change(:,:,icarbon)     = -SUM(turnover_daily(:,:,:,icarbon),DIM=3) &
                                      -SUM(bm_to_litter(:,:,:,icarbon),DIM=3)  ! in case ok_dvgm=true
       mass_change(:,:,initrogen)   = -SUM(turnover_daily(:,:,:,initrogen),DIM=3) &
                                      -SUM(bm_to_litter(:,:,:,initrogen),DIM=3)
       mass_change(:,:,iphosphorus) = -SUM(turnover_daily(:,:,:,iphosphorus),DIM=3) &
                                      -SUM(bm_to_litter(:,:,:,iphosphorus),DIM=3)

       mass_after = SUM(biomass(:,:,:,:),DIM=3)
       CALL cons_mass( mass_before(:,:,:),               &  ! mass before
             mass_after,                       &  ! mass afte
             mass_change(:,:,:),               &  ! net of fluxes
             'lpj: after turn')
    ENDIF

!! DSG: mass conservation issue: END
!! =====================================
    
  !! 11. Establishment of saplings
    
    IF ( ok_dgvm .OR. .NOT.lpj_gap_const_mort ) THEN
       !DSG mass conservation ========================================
       mass_before(:,:,:) = SUM(biomass(:,:,:,:),DIM=3) + SUM(bm_to_litter(:,:,:,:),DIM=3)

       !! 11.1 Establish new plants
       CALL establish (npts, dt_days, PFTpresent, regenerate, &
            neighbours, resolution, need_adjacent, herbivores, &
            soil_orders, max_eau_var, bulk,                     &
            precip_lastyear, gdd0_lastyear, lm_lastyearmax, &
            cn_ind, lai, avail_tree, avail_grass, npp_longterm, &
            leaf_age, leaf_frac, &
            ind, biomass, age, everywhere, co2_to_bm, &
            soil_n_min, n_uptake_daily, soil_p_min, p_uptake_daily, &
            nstress_season, pstress_season,veget_cov_max, woodmass_ind, &
            mortality, bm_to_litter, &
            sla_calc)

       IF (printlev>=3) WRITE (numout,*) 'after establish soil_n_min(test_grid,test_pft,:):',soil_n_min(test_grid,test_pft,:)
       !! 12.2 Calculate new crown area (and maximum vegetation cover)
       CALL crown (npts, PFTpresent, &
            ind, biomass, woodmass_ind, &
            veget_cov_max, cn_ind, height)

       IF(dsg_debug) THEN
          !DSG debug: ===================================================
          CALL  check_mass(npts,biomass(:,:,:,:), 'lpj: after establish')

          !DSG mass conservation ============================================
          mass_change(:,:,icarbon)     = -SUM(bm_to_litter(:,:,:,icarbon),DIM=3)
          mass_change(:,:,initrogen)   = -SUM(bm_to_litter(:,:,:,initrogen),DIM=3)
          mass_change(:,:,iphosphorus) = -SUM(bm_to_litter(:,:,:,iphosphorus),DIM=3)

          mass_after =  SUM(biomass(:,:,:,:),DIM=3)
          CALL cons_mass( mass_before(:,:,:),               &  ! mass before
                mass_after,                       &  ! mass afte
                mass_change(:,:,:),               &  ! net of fluxes
                'lpj: after establish')
       ENDIF

    ENDIF

    !JC add NPinput
    IF (allow_agri_fert) THEN
     mass_before(:,:,:) = SUM(som(:,:,:,:),DIM=2) + SUM(SUM(litter(:,:,:,:,:),DIM=4),DIM=2)
       
      CALL apply_agri_fert(npts,dt_days,day_of_year,veget_cov_max,&
             agri_fert_save,litter,som,agri_nfert,agri_pfert,fert_masschangecnp)  
       IF(dsg_debug) THEN
          !DSG debug: ===================================================

          CALL  check_mass(npts,biomass(:,:,:,:), 'stomateLPJ: after apply_agri_fert')

          !DSG mass conservation ============================================
          !JC NOTE: only solid organic fertilizer is applied to litter
          ! mineral fertilizer and the urine N in organic fertilizer
          ! is output through agri_nfert and agri_pfert
          ! and will be applied directly to soil mineral N P pools next day
          mass_change(:,:,icarbon)     = fert_masschangecnp(:,:,icarbon)
          mass_change(:,:,initrogen)   = fert_masschangecnp(:,:,initrogen)
          mass_change(:,:,iphosphorus) = fert_masschangecnp(:,:,iphosphorus)

          mass_after = SUM(som(:,:,:,:),DIM=2) + SUM(SUM(litter(:,:,:,:,:),DIM=4),DIM=2)
          CALL cons_mass( mass_before(:,:,:),               &  ! mass before
                mass_after,                       &  ! mass afte
                mass_change(:,:,:),                &  ! net of fluxes
                'lpj: after apply_agri_fert')
       ENDIF

    ENDIF 
    !End JC add NPinput

    !gmjc Grassland_management
     mass_before(:,:,:) = SUM(biomass(:,:,:,:),DIM=3) + &
                                SUM(bm_to_litter(:,:,:,:),DIM=3) + &
                                SUM(SUM(litter(:,:,:,:,:),DIM=4),DIM=2)
    !
    ! 13 calculate grazing by animals or cutting for forage
    !
    IF (GRM_enable_grazing) THEN
        CALL main_grassland_management(&
           npts, lalo, neighbours, resolution, contfrac, &
           dt_days        , &
           day_of_year    , &
           t2m_daily      , &
           t2m_min_daily  , &
           t2m_14         , &
           tsurf_daily    , &
           snowfall_daily , &
           biomass        , &
           bm_to_litter   , &
           litter         , &
           litter_avail   , &
           litter_not_avail , &
!           !spitfire
!           fuel_1hr(:,:,:,icarbon), &
!           fuel_10hr(:,:,:,icarbon), &
!           fuel_100hr(:,:,:,icarbon), &
!           fuel_1000hr(:,:,:,icarbon), &
!           !end spitfire
           .TRUE., LastTsYear, when_growthinit_cut, nb_grazingdays, &
           lai,sla_calc,leaf_age,leaf_frac, &
           wshtotsum,sr_ugb,compt_ugb, &
           nb_ani,grazed_frac,import_yield, &
           moiavail_daily,tmc_topgrass_daily,fc_grazing, snowmass_daily, &
           after_snow, after_wet, wet1day, wet2day, &
           harvest_gm, ranimal_gm, ch4_pft_gm, cinput_gm, n2o_pft_gm, &
           grm_nfert, grm_pfert, &
           grm_masschangecnp,GRM_devstage,when_growthinit)

    ENDIF
       IF(dsg_debug) THEN
          !DSG debug: ===================================================
          CALL  check_mass(npts,biomass(:,:,:,:), 'stomateLPJ: after main_grassland_management')

          !DSG mass conservation ============================================
          ! JC NOTED 
          ! the grassland management has included many fluxes from or to
          ! environment (such as fertilization input, harvested biomass output,
          ! milk output, liveweight output, animal respiration output,
          ! and enteric fermentation ch4 emission)
          ! HOWEVER, it should be noted that grm_nfert (including N in urinecnp) 
          ! and grm_ndep are N fluxes
          ! grm_pfert (including N in urinecnp) are P fluxes
          ! that will directly goes into soil in the next time step, thus they
          ! are not included in the mass balance calculation
          mass_change(:,:,icarbon)     = grm_masschangecnp(:,:,icarbon)
          mass_change(:,:,initrogen)   = grm_masschangecnp(:,:,initrogen)
          mass_change(:,:,iphosphorus) = grm_masschangecnp(:,:,iphosphorus)

          mass_after = SUM(biomass(:,:,:,:),DIM=3) + SUM(bm_to_litter(:,:,:,:),DIM=3) + &
                SUM(SUM(litter(:,:,:,:,:),DIM=4),DIM=2)

          CALL cons_mass( mass_before(:,:,:),     &  ! mass before
                mass_after(:,:,:),                & ! mass after
                mass_change(:,:,:),               &  ! net of fluxes
                'lpj: after main_grassland_management')
       ENDIF
       !end gmjc


  !! 13. Calculate final LAI and vegetation cover
    
    CALL cover (npts, cn_ind, ind, biomass, &
         veget_cov_max, veget_cov_max_tmp, lai, &
         litter, &
         litter_avail, litter_not_avail, &
         som, turnover_daily, bm_to_litter, &
         co2_to_bm, co2_fire, resp_hetero, resp_maint, resp_growth, gpp_daily, &
         lignin_struc, lignin_wood, soil_n_min, soil_p_min) ! DSG: here the ratio of idissolved and isorbed is not changed, thus manipulation of soil_p_min is okay.


    IF (printlev>=3) WRITE (numout,*) 'after cover soil_n_min(test_grid,test_pft,:):',soil_n_min(test_grid,test_pft,:)
  !! 14. Update litter pools to account for harvest
 
    ! the whole litter stuff:
    !    litter update, lignin content, PFT parts, litter decay, 
    !    litter heterotrophic respiration, dead leaf soil cover.
    !    No vertical discretisation in the soil for litter decay.\n
    ! added by shilong for harvest
    IF(harvest_agri) THEN
       CALL harvest(npts, dt_days, veget_cov_max, &
            bm_to_litter, turnover_daily, &
            harvest_above, harvest_masschangecnp)
    ENDIF

    IF(harvest_BE) THEN
       IF (printlev>=3) WRITE (numout,*) 'we call harvest_biograss'
       CALL harvest_biograss(npts, dt_days, veget_cov_max, &
            bm_to_litter, turnover_daily, &
            harvest_bioenergy, harvest_BE_masschangecnp)
    ENDIF

    !! XX. Update turnover_daily to account for grazing
    ! this simplified subroutine is used to balance the
    ! assumed N input through grazing excreta deposition
    IF (allow_agri_fert) THEN
       CALL simplified_grazing(npts, dt_days, day_of_year, veget_cov_max, &
            agri_fert_save, turnover_daily, rest_Ngrazing, grazed_above)
    ENDIF

    !! 14. Land cover change if it is time to do so
    !! The flag do_now_lcchange is set in slowproc_main at the same time as the vegetation is read from file.
    !! The vegetation fractions are not updated yet and will be updated in the end of sechiba_main.
    IF (do_now_stomate_woodharvest) THEN 
       CALL woodharvest_main(npts, dt_days, veget_cov_max, &
            biomass, &
            flux10_harvest,flux100_harvest, prod10_harvest,prod100_harvest,&
            convflux_harvest,cflux_prod10_harvest,cflux_prod100_harvest,&
            woodharvest,woodharvestpft)
       do_now_stomate_woodharvest=.FALSE.
    ENDIF

    IF(dsg_debug) THEN
       !DSG check if mass gets negative:
       CALL  check_mass(npts,biomass(:,:,:,:),'lpj: after wood_harvest')
    ENDIF

    IF (do_now_stomate_lcchange) THEN
       CALL lcchange_main (npts, dt_days, veget_cov_max, veget_cov_max_new, &
            biomass, ind, age, PFTpresent, senescence, when_growthinit, everywhere, &
            co2_to_bm, n_to_bm, p_to_bm, &
            bm_to_litter, turnover_daily, bm_sapl, cn_ind,flux10,flux100, &
            prod10,prod100,convflux,cflux_prod10,cflux_prod100,nflux_prod_total,pflux_prod_total,leaf_frac,&
            flux10_harvest,flux100_harvest, prod10_harvest,prod100_harvest,&
            convflux_harvest,cflux_prod10_harvest,cflux_prod100_harvest,&        
            npp_longterm, lm_lastyearmax, litter, &
            litter_avail, litter_not_avail, &
            som, soil_n_min, soil_p_min, KF, k_latosa_adapt,rue_longterm, &  ! DSG: here the ratio of idissolved and isorbed is not changed, thus manipulation of soil_p_min is okay.
            lignin_struc, lignin_wood,&
            woodharvest, woodharvestpft)  ! do_wood_harvest 
            
       do_now_stomate_lcchange=.FALSE.

       ! Set the flag done_stomate_lcchange to be used in the end of sechiba_main to update the fractions.
       done_stomate_lcchange=.TRUE.
    ENDIF

    IF(dsg_debug) THEN
       !DSG check if mass gets negative:
       CALL  check_mass(npts,biomass(:,:,:,:),'lpj: after lcchange_main')
    ENDIF

     IF (printlev>=3) WRITE (numout,*) 'after lcchange soil_n_min(test_grid,test_pft,:):',soil_n_min(test_grid,test_pft,:)

    !! 15. Calculate vcmax 
    CALL vmax (npts, dt_days, &
         biomass,             &
         leaf_age, leaf_frac, &
         vcmax, jmax, nue)



    !MM déplacement pour initialisation correcte des grandeurs cumulées :
    cflux_prod_total(:) = convflux(:) + cflux_prod10(:) + cflux_prod100(:)
    prod10_total(:)=SUM(prod10,dim=2)
    prod100_total(:)=SUM(prod100,dim=2)

    cflux_prod_harvest_total(:) = convflux_harvest(:) + cflux_prod10_harvest(:) + cflux_prod100_harvest(:)
    prod10_harvest_total(:)=SUM(prod10_harvest,dim=2)
    prod100_harvest_total(:)=SUM(prod100_harvest,dim=2)
    
  !! 16. Total heterotrophic respiration

    tot_soil_carb(:,:) = zero
    tot_litter_carb(:,:) = zero
    tot_soil_nitr(:,:) = zero
    tot_litter_nitr(:,:) = zero
    tot_soil_phos(:,:) = zero
    tot_litter_phos(:,:) = zero

    DO j=2,nvm

       tot_litter_carb(:,j) = tot_litter_carb(:,j) + (litter(:,istructural,j,iabove,icarbon) + &
            &          litter(:,imetabolic,j,iabove,icarbon) + litter(:,iwoody,j,iabove,icarbon) + &
            &          litter(:,istructural,j,ibelow,icarbon) + litter(:,imetabolic,j,ibelow,icarbon)+ &
            &          litter(:,iwoody,j,ibelow,icarbon))

       tot_soil_carb(:,j) = tot_soil_carb(:,j) + (som(:,iactive,j,icarbon) + &
            &          som(:,islow,j,icarbon)+  som(:,ipassive,j,icarbon)+  som(:,isurface,j,icarbon))

       tot_litter_nitr(:,j) = tot_litter_nitr(:,j) + (litter(:,istructural,j,iabove,initrogen) + &
            &          litter(:,imetabolic,j,iabove,initrogen) + litter(:,iwoody,j,iabove,initrogen) + &
            &          litter(:,istructural,j,ibelow,initrogen) + litter(:,imetabolic,j,ibelow,initrogen)+ &
            &          litter(:,iwoody,j,ibelow,initrogen))

       tot_soil_nitr(:,j) = tot_soil_nitr(:,j) + (som(:,iactive,j,initrogen) + &
            &          som(:,islow,j,initrogen)+  som(:,ipassive,j,initrogen)+ som(:,isurface,j,initrogen))

       tot_litter_phos(:,j) = tot_litter_phos(:,j) + (litter(:,istructural,j,iabove,iphosphorus) + &
            &          litter(:,imetabolic,j,iabove,iphosphorus) + litter(:,iwoody,j,iabove,iphosphorus) + &
            &          litter(:,istructural,j,ibelow,iphosphorus) + litter(:,imetabolic,j,ibelow,iphosphorus)+ &
            &          litter(:,iwoody,j,ibelow,iphosphorus))

       tot_soil_phos(:,j) = tot_soil_phos(:,j) + (som(:,iactive,j,iphosphorus) + &
            &          som(:,islow,j,iphosphorus)+  som(:,ipassive,j,iphosphorus)+ som(:,isurface,j,iphosphorus))
    ENDDO
    tot_litter_soil_carb(:,:) = tot_litter_carb(:,:) + tot_soil_carb(:,:)

    tot_live_biomass(:,:,:) = biomass(:,:,ileaf,:) + biomass(:,:,isapabove,:) + biomass(:,:,isapbelow,:) +&
             &                    biomass(:,:,iheartabove,:) + biomass(:,:,iheartbelow,:) + &
             &                    biomass(:,:,iroot,:) + biomass(:,:,ifruit,:) + biomass(:,:,icarbres,:) + biomass(:,:,ilabile,:)


    tot_turnover(:,:,:) = turnover_daily(:,:,ileaf,:) + turnover_daily(:,:,isapabove,:) + &
         &         turnover_daily(:,:,isapbelow,:) + turnover_daily(:,:,iheartabove,:) + &
         &         turnover_daily(:,:,iheartbelow,:) + turnover_daily(:,:,iroot,:) + &
         &         turnover_daily(:,:,ifruit,:) + turnover_daily(:,:,icarbres,:) + turnover_daily(:,:,ilabile,:) 

    tot_bm_to_litter(:,:,:) = bm_to_litter(:,:,ileaf,:) + bm_to_litter(:,:,isapabove,:) +&
         &             bm_to_litter(:,:,isapbelow,:) + bm_to_litter(:,:,iheartbelow,:) +&
         &             bm_to_litter(:,:,iheartabove,:) + bm_to_litter(:,:,iroot,:) + &
         &             bm_to_litter(:,:,ifruit,:) + bm_to_litter(:,:,icarbres,:) + bm_to_litter(:,:,ilabile,:)

    carb_mass_variation(:)=-carb_mass_total(:)
    carb_mass_total(:)=SUM((tot_live_biomass(:,:,icarbon)+tot_litter_carb+tot_soil_carb)*veget_cov_max,dim=2) + &
         &                 (prod10_total + prod100_total) +  (prod10_harvest_total + prod100_harvest_total)
    carb_mass_variation(:)=carb_mass_total(:)+carb_mass_variation(:)

    !DSG mass conservation ========================================
    mass_after(:,:,:) = SUM(som(:,:,:,:),DIM=2) + SUM(SUM(litter(:,:,:,:,:),DIM=4),DIM=2) +  SUM(biomass(:,:,:,:),DIM=3)
    mass_after(:,:,iphosphorus) = mass_after(:,:,iphosphorus) +  SUM(soil_p_min(:,:,:),DIM=3)
    mass_after(:,:,initrogen)   = mass_after(:,:,initrogen)   +  SUM(soil_n_min(:,:,:),DIM=3)

    ! MISSING: crop harvest per PFT. I had to disable mass conservation checks
    ! for crops in stomate_phosphorus. This needs to be reverted as soon harvest
    ! flux is available to be added here:
    mass_change(:,:,:)     =  - SUM(bm_to_litter(:,:,:,:),DIM=3)   &
                              - SUM(turnover_daily(:,:,:,:),DIM=3) &
                              - harvest_masschangecnp(:,:,:)       & ! food  crop harvest export (loss)
                              - harvest_BE_masschangecnp(:,:,:)    & ! fiber crop harvest export (loss)
                              - fire_masschangecnp(:,:,:)          & ! fire emissions (loss)
                              + grm_masschangecnp(:,:,:)           & ! input - output (gain)
                              + fert_masschangecnp(:,:,:)            ! input (gain)

    mass_change(:,:,icarbon)      =   mass_change(:,:,icarbon)     + npp_daily(:,:)*dt_days  + co2_to_bm(:,:)

    mass_change(:,:,initrogen)    =   mass_change(:,:,initrogen)   + SUM(n_uptake_daily(:,:,:),DIM=3) + N_support(:,:) + BNF_clover_daily(:,:) + n_to_bm(:,:)
    mass_change(:,:,iphosphorus)  =   mass_change(:,:,iphosphorus) + p_uptake_daily(:,:) + P_support(:,:) + p_to_bm(:,:)

    CALL cons_mass( mass_before_2(:,:,:),   &  ! mass before
           mass_after,                      &  ! mass after
           mass_change(:,:,:),              &  ! net of fluxes
           'lpj: all routines' )
    
  !! 17. Write history

    CALL xios_orchidee_send_field("RESOLUTION_X",resolution(:,1))
    CALL xios_orchidee_send_field("RESOLUTION_Y",resolution(:,2))
    CALL xios_orchidee_send_field("CONTFRAC_STOMATE",contfrac(:))
    CALL xios_orchidee_send_field("T2M_MONTH",t2m_month)
    CALL xios_orchidee_send_field("T2M_WEEK",t2m_week)
    CALL xios_orchidee_send_field("TSEASON",Tseason)
    CALL xios_orchidee_send_field("TMIN_SPRING_TIME", Tmin_spring_time)
    CALL xios_orchidee_send_field("ONSET_DATE",onset_date)
    CALL xios_orchidee_send_field("FPC_MAX",fpc_max)
    CALL xios_orchidee_send_field("MAXFPC_LASTYEAR",maxfpc_lastyear)
    CALL xios_orchidee_send_field("HET_RESP",resp_hetero(:,:))
    CALL xios_orchidee_send_field("CO2_FIRE",co2_fire)
    CALL xios_orchidee_send_field("CO2_TAKEN",co2_to_bm)
    IF(ok_ncycle) CALL xios_orchidee_send_field("N_TAKEN",n_to_bm)
    IF(ok_pcycle) CALL xios_orchidee_send_field("P_TAKEN",p_to_bm)
    CALL xios_orchidee_send_field("LAI",lai)
    CALL xios_orchidee_send_field("VEGET_COV_MAX",veget_cov_max)
    CALL xios_orchidee_send_field("NPP_STOMATE",npp_daily)
    CALL xios_orchidee_send_field("GPP",gpp_daily)
    CALL xios_orchidee_send_field("IND",ind)
    CALL xios_orchidee_send_field("CN_IND",cn_ind)
    CALL xios_orchidee_send_field("WOODMASS_IND",woodmass_ind)
    CALL xios_orchidee_send_field("MOISTRESS",moiavail_week)
    CALL xios_orchidee_send_field("MAINT_RESP",resp_maint)
    CALL xios_orchidee_send_field("GROWTH_RESP",resp_growth)

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
      
       CALL xios_orchidee_send_field("TOTAL_M"//TRIM(element_str(l)),tot_live_biomass(:,:,l))
       CALL xios_orchidee_send_field("LEAF_M"//TRIM(element_str(l)),biomass(:,:,ileaf,l))
       CALL xios_orchidee_send_field("SAP_M_AB"//TRIM(element_str(l)),biomass(:,:,isapabove,l))
       CALL xios_orchidee_send_field("SAP_M_BE"//TRIM(element_str(l)),biomass(:,:,isapbelow,l))
       CALL xios_orchidee_send_field("HEART_M_AB"//TRIM(element_str(l)),biomass(:,:,iheartabove,l))
       CALL xios_orchidee_send_field("HEART_M_BE"//TRIM(element_str(l)),biomass(:,:,iheartbelow,l))
       CALL xios_orchidee_send_field("ROOT_M"//TRIM(element_str(l)),biomass(:,:,iroot,l))
       CALL xios_orchidee_send_field("FRUIT_M"//TRIM(element_str(l)),biomass(:,:,ifruit,l))
       CALL xios_orchidee_send_field("LABILE_M"//TRIM(element_str(l)),biomass(:,:,ilabile,l))
       CALL xios_orchidee_send_field("RESERVE_M"//TRIM(element_str(l)),biomass(:,:,icarbres,l))
       CALL xios_orchidee_send_field("TOTAL_TURN"//TRIM(element_str(l)),tot_turnover(:,:,l))
       CALL xios_orchidee_send_field("LEAF_TURN"//TRIM(element_str(l)),turnover_daily(:,:,ileaf,l))
       CALL xios_orchidee_send_field("SAP_AB_TURN"//TRIM(element_str(l)),turnover_daily(:,:,isapabove,l))
       CALL xios_orchidee_send_field("ROOT_TURN"//TRIM(element_str(l)),turnover_daily(:,:,iroot,l))
       CALL xios_orchidee_send_field("FRUIT_TURN"//TRIM(element_str(l)),turnover_daily(:,:,ifruit,l))
       CALL xios_orchidee_send_field("RESERVE_TURN"//TRIM(element_str(l)),turnover_daily(:,:,icarbres,l))
       CALL xios_orchidee_send_field("LABILE_TURN"//TRIM(element_str(l)),turnover_daily(:,:,ilabile,l))
       IF  ( (l == initrogen).OR. (l == iphosphorus) ) THEN
          CALL xios_orchidee_send_field("LEAF_RECYCLE"//TRIM(element_str(l)),recycling_daily(:,:,ileaf,l))
          CALL xios_orchidee_send_field("ROOT_RECYCLE"//TRIM(element_str(l)),recycling_daily(:,:,iroot,l))
          CALL xios_orchidee_send_field("SAP_AB_RECYCLE"//TRIM(element_str(l)),recycling_daily(:,:,isapabove,l))
       ENDIF
!DSG: this is not per day, but per timestep!!!
       CALL xios_orchidee_send_field("TOTAL_BM_LITTER"//TRIM(element_str(l)),tot_bm_to_litter(:,:,l))
       CALL xios_orchidee_send_field("LEAF_BM_LITTER"//TRIM(element_str(l)),bm_to_litter(:,:,ileaf,l))
       CALL xios_orchidee_send_field("SAP_AB_BM_LITTER"//TRIM(element_str(l)),bm_to_litter(:,:,isapabove,l))
       CALL xios_orchidee_send_field("SAP_BE_BM_LITTER"//TRIM(element_str(l)),bm_to_litter(:,:,isapbelow,l))
       CALL xios_orchidee_send_field("HEART_AB_BM_LITTER"//TRIM(element_str(l)),bm_to_litter(:,:,iheartabove,l))
       CALL xios_orchidee_send_field("HEART_BE_BM_LITTER"//TRIM(element_str(l)),bm_to_litter(:,:,iheartbelow,l))
       CALL xios_orchidee_send_field("ROOT_BM_LITTER"//TRIM(element_str(l)),bm_to_litter(:,:,iroot,l))
       CALL xios_orchidee_send_field("FRUIT_BM_LITTER"//TRIM(element_str(l)),bm_to_litter(:,:,ifruit,l))
       CALL xios_orchidee_send_field("LABILE_BM_LITTER"//TRIM(element_str(l)),bm_to_litter(:,:,ilabile,l))
       CALL xios_orchidee_send_field("RESERVE_BM_LITTER"//TRIM(element_str(l)),bm_to_litter(:,:,icarbres,l))
!DSG: this is not per day!!!
       CALL xios_orchidee_send_field("LITTER_STR_AB"//TRIM(element_str(l)),litter(:,istructural,:,iabove,l))
       CALL xios_orchidee_send_field("LITTER_MET_AB"//TRIM(element_str(l)),litter(:,imetabolic,:,iabove,l))
       CALL xios_orchidee_send_field("LITTER_WOD_AB"//TRIM(element_str(l)),litter(:,iwoody,:,iabove,l))
       CALL xios_orchidee_send_field("LITTER_STR_BE"//TRIM(element_str(l)),litter(:,istructural,:,ibelow,l))
       CALL xios_orchidee_send_field("LITTER_MET_BE"//TRIM(element_str(l)),litter(:,imetabolic,:,ibelow,l))
       CALL xios_orchidee_send_field("LITTER_WOD_BE"//TRIM(element_str(l)),litter(:,iwoody,:,ibelow,l))
       CALL xios_orchidee_send_field("SOIL_ACTIVE"//TRIM(element_str(l)),som(:,iactive,:,l))
       CALL xios_orchidee_send_field("SOIL_SLOW"//TRIM(element_str(l)),som(:,islow,:,l))
       CALL xios_orchidee_send_field("SOIL_PASSIVE"//TRIM(element_str(l)),som(:,ipassive,:,l))
       CALL xios_orchidee_send_field("SOIL_SURF"//TRIM(element_str(l)),som(:,isurface,:,l))
       ! Write in stomat_litter.f90
       !CALL xios_orchidee_send_field("SOM_INPUT_ACTIVE"//TRIM(element_str(l)),som_input(:,iactive,:,l))
       !CALL xios_orchidee_send_field("SOM_INPUT_SLOW"//TRIM(element_str(l)),som_input(:,islow,:,l))
       !CALL xios_orchidee_send_field("SOM_INPUT_PASSIVE"//TRIM(element_str(l)),som_input(:,ipassive,:,l))
       !CALL xios_orchidee_send_field("SOM_INPUT_SURF"//TRIM(element_str(l)),som_input(:,isurface,:,l))

    ENDDO

    CALL xios_orchidee_send_field("DEADLEAF_COVER",deadleaf_cover)
    CALL xios_orchidee_send_field("TOTAL_SOIL_CARB",tot_litter_soil_carb)

    CALL xios_orchidee_send_field("LITTERHUM",litterhum_daily)
    CALL xios_orchidee_send_field("TURNOVER_TIME",turnover_time)
    CALL xios_orchidee_send_field("PROD10",prod10)
    CALL xios_orchidee_send_field("FLUX10",flux10)
    CALL xios_orchidee_send_field("PROD100",prod100)
    CALL xios_orchidee_send_field("FLUX100",flux100)
    CALL xios_orchidee_send_field("CONVFLUX",convflux)
    CALL xios_orchidee_send_field("CFLUX_PROD10",cflux_prod10)
    CALL xios_orchidee_send_field("CFLUX_PROD100",cflux_prod100)
    CALL xios_orchidee_send_field("PROD10_HARVEST",prod10_harvest)
    CALL xios_orchidee_send_field("FLUX10_HARVEST",flux10_harvest)
    CALL xios_orchidee_send_field("PROD100_HARVEST",prod100_harvest)
    CALL xios_orchidee_send_field("FLUX100_HARVEST",flux100_harvest)
    CALL xios_orchidee_send_field("CONVFLUX_HARVEST",convflux_harvest)
    CALL xios_orchidee_send_field("CFLUX_PROD10_HARVEST",cflux_prod10_harvest)
    CALL xios_orchidee_send_field("CFLUX_PROD100_HARVEST",cflux_prod100_harvest)
    CALL xios_orchidee_send_field("WOOD_HARVEST",woodharvest/one_year*dt_days)
    CALL xios_orchidee_send_field("WOOD_HARVEST_PFT",woodharvestpft)
    CALL xios_orchidee_send_field("HARVEST_ABOVE",harvest_above)
    CALL xios_orchidee_send_field("HARVEST_BIOENERGY",harvest_bioenergy)
    CALL xios_orchidee_send_field("GRAZED_ABOVE",grazed_above)
    CALL xios_orchidee_send_field("GRM_DEVSTAGE",GRM_devstage)

    vcmax_new=zero
    ! DSG: vcmax new is the nitrogen derived vcmax (which is used instead of
    ! vcmax25)
    ! DSG: it is based on Kattge et al (2009)
    DO j=2,nvm
           !DSG avoid division by zero; not quite sure how LAI>zero and leaf C =zero can occur
           ! - but it does. This might indicate something went wrong earlier.
       WHERE((lai(:,j) .GT. min_stomate) &
          !DSGadded
             .AND.(biomass(:,j,ileaf,icarbon).GT. zero) )
          !DSGadded
           ! use the top leaf N concentration
           vcmax_new(:,j) = nue(:,j)*biomass(:,j,ileaf,initrogen) * ext_coeff_N(j)  &
                          / ( 1. - exp(-ext_coeff_N(j) * sla_calc(:,j)*biomass(:,j,ileaf,icarbon)) ) 
       ELSEWHERE
           vcmax_new(:,j) = zero
       ENDWHERE
    ENDDO

    CALL xios_orchidee_send_field("VCMAX_ELLSWORTH",vcmax(:,:,1))
    CALL xios_orchidee_send_field("VCMAX_KATTGE"   ,vcmax_new)
    CALL xios_orchidee_send_field("AGE",age)
    CALL xios_orchidee_send_field("HEIGHT",height)
    CALL xios_orchidee_send_field("FIREINDEX",fireindex(:,:))

    !gmjc
    CALL xios_orchidee_send_field("T2M_14",t2m_14(:))
    CALL xios_orchidee_send_field("LITTER_STR_AVAIL",litter_avail(:,istructural,:))
    CALL xios_orchidee_send_field("LITTER_MET_AVAIL",litter_avail(:,imetabolic,:))
    CALL xios_orchidee_send_field("LITTER_STR_NAVAIL",litter_not_avail(:,istructural,:))
    CALL xios_orchidee_send_field("LITTER_MET_NAVAIL",litter_not_avail(:,imetabolic,:))
    CALL xios_orchidee_send_field("LITTER_STR_AVAILF",litter_avail_frac(:,istructural,:))
    CALL xios_orchidee_send_field("LITTER_MET_AVAILF",litter_avail_frac(:,imetabolic,:))
    !calculate grassland co2 fluxes
    DO j=2,nvm
      IF ((.NOT. is_tree(j)) .AND. natural(j)) THEN
        veget_cov_max_gm(:,j) = veget_cov_max(:,j)
      ENDIF
    END DO ! nvm
    veget_exist_gm(:) = SUM(veget_cov_max_gm,dim=2)
    WHERE (veget_exist_gm(:) .GT. 0.0)
      co2_gm(:) = SUM((gpp_daily-(resp_maint+resp_growth+resp_hetero)-co2_fire &
                    -harvest_gm-ranimal_gm-ch4_pft_gm+cinput_gm) &
                    *veget_cov_max_gm,dim=2)/veget_exist_gm
      ch4_gm(:) = SUM(ch4_pft_gm*veget_cov_max_gm,dim=2)/veget_exist_gm
      n2o_gm(:) = SUM(n2o_pft_gm*veget_cov_max_gm,dim=2)/veget_exist_gm
    ELSEWHERE
      co2_gm(:) = zero
      ch4_gm(:) = zero
      n2o_gm(:) = zero
    ENDWHERE
    CALL xios_orchidee_send_field("CO2_GM",co2_gm)
    CALL xios_orchidee_send_field("CH4_GM",ch4_gm)
    CALL xios_orchidee_send_field("N2O_GM",n2o_gm)
    CALL xios_orchidee_send_field("N2O_PFT_GM",n2o_pft_gm)
    !end gmjc

    CALL xios_orchidee_send_field("nVeg",SUM(tot_live_biomass(:,:,initrogen)*veget_cov_max,dim=2)/1e3*contfrac)
    CALL xios_orchidee_send_field("nLitter",SUM(tot_litter_nitr*veget_cov_max,dim=2)/1e3*contfrac)
    CALL xios_orchidee_send_field("nSoil",SUM(tot_soil_nitr*veget_cov_max,dim=2)/1e3*contfrac)
    CALL xios_orchidee_send_field("pVeg",SUM(tot_live_biomass(:,:,iphosphorus)*veget_cov_max,dim=2)/1e3*contfrac)
    CALL xios_orchidee_send_field("pLitter",SUM(tot_litter_phos*veget_cov_max,dim=2)/1e3*contfrac)
    CALL xios_orchidee_send_field("pSoil",SUM(tot_soil_phos*veget_cov_max,dim=2)/1e3*contfrac)


! ipcc history
     CALL xios_orchidee_send_field("cVeg",SUM(tot_live_biomass(:,:,icarbon)*veget_cov_max,dim=2)/1e3*contfrac)
     CALL xios_orchidee_send_field("cLitter",SUM(tot_litter_carb*veget_cov_max,dim=2)/1e3*contfrac)
     CALL xios_orchidee_send_field("cSoil",SUM(tot_soil_carb*veget_cov_max,dim=2)/1e3*contfrac)
     CALL xios_orchidee_send_field("cProduct",(prod10_total + prod100_total)/1e3)
     CALL xios_orchidee_send_field("cMassVariation",carb_mass_variation/1e3/one_day*contfrac)
     CALL xios_orchidee_send_field("lai_ipcc",SUM(lai*veget_cov_max,dim=2)*contfrac)
     CALL xios_orchidee_send_field("gpp_ipcc",SUM(gpp_daily*veget_cov_max,dim=2)/1e3/one_day*contfrac)
     CALL xios_orchidee_send_field("ra",SUM((resp_maint+resp_growth)*veget_cov_max,dim=2)/1e3/one_day*contfrac)
     CALL xios_orchidee_send_field("npp_ipcc",SUM(npp_daily*veget_cov_max,dim=2)/1e3/one_day*contfrac)
     CALL xios_orchidee_send_field("rh",SUM(resp_hetero*veget_cov_max,dim=2)/1e3/one_day*contfrac)
     CALL xios_orchidee_send_field("fFire",SUM(co2_fire*veget_cov_max,dim=2)/1e3/one_day*contfrac)
     CALL xios_orchidee_send_field("fHarvest",(harvest_above+SUM(harvest_bioenergy*veget_cov_max,dim=2))/1e3/one_day*contfrac)
     CALL xios_orchidee_send_field("fLuc",cflux_prod_total/1e3/one_day*contfrac)
     IF (allow_agri_fert) THEN
     ! include the input and output of organic carbon
     ! input is the organic part of agricultural fertiliser over cropland
     ! (manure), pasture (manure + grazing deposition), rangeland (grazing
     ! deposition). output is the harvsted crop and grazed grass
     CALL xios_orchidee_send_field("nbp",(SUM((gpp_daily-(resp_maint+resp_growth+resp_hetero)-co2_fire &
                                               + fert_masschangecnp(:,:,icarbon)- harvest_bioenergy) &
          & *veget_cov_max,dim=2)-cflux_prod_total-harvest_above-grazed_above)/1e3/one_day*contfrac)
     ELSE 
     CALL xios_orchidee_send_field("nbp",(SUM((gpp_daily-(resp_maint+resp_growth+resp_hetero) &
                                               -co2_fire-harvest_bioenergy) &
          &        *veget_cov_max,dim=2)-cflux_prod_total-harvest_above)/1e3/one_day*contfrac)
     ENDIF
     CALL xios_orchidee_send_field("fVegLitter",SUM((tot_bm_to_litter(:,:,icarbon) + tot_turnover(:,:,icarbon))*&
          veget_cov_max,dim=2)/1e3/one_day*contfrac)
     CALL xios_orchidee_send_field("fLitterSoil",SUM(SUM(som_input(:,:,:,icarbon),dim=2)*veget_cov_max,dim=2)/1e3/one_day)
     CALL xios_orchidee_send_field("cLeaf",SUM(biomass(:,:,ileaf,icarbon)*veget_cov_max,dim=2)/1e3*contfrac)
     CALL xios_orchidee_send_field("cWood",SUM((biomass(:,:,isapabove,icarbon)+biomass(:,:,iheartabove,icarbon))*&
          veget_cov_max,dim=2)/1e3*contfrac)
     CALL xios_orchidee_send_field("cRoot",SUM(( biomass(:,:,iroot,icarbon) + biomass(:,:,isapbelow,icarbon) + &
          biomass(:,:,iheartbelow,icarbon) )*veget_cov_max,dim=2)/1e3*contfrac)
     CALL xios_orchidee_send_field("cMisc",SUM(( biomass(:,:,icarbres,icarbon) + biomass(:,:,ifruit,icarbon))*&
          veget_cov_max,dim=2)/1e3*contfrac)
     CALL xios_orchidee_send_field("cLitterAbove",SUM((litter(:,istructural,:,iabove,icarbon)+&
          litter(:,imetabolic,:,iabove,icarbon))*veget_cov_max,dim=2)/1e3*contfrac)
     CALL xios_orchidee_send_field("cLitterBelow",SUM((litter(:,istructural,:,ibelow,icarbon)+&
          litter(:,imetabolic,:,ibelow,icarbon))*veget_cov_max,dim=2)/1e3*contfrac)
     CALL xios_orchidee_send_field("cSoilFast",SUM((som(:,iactive,:,icarbon)+som(:,isurface,:,icarbon))*veget_cov_max,dim=2)/1e3)
     CALL xios_orchidee_send_field("cSoilMedium",SUM(som(:,islow,:,icarbon)*veget_cov_max,dim=2)/1e3)
     CALL xios_orchidee_send_field("cSoilSlow",SUM(som(:,ipassive,:,icarbon)*veget_cov_max,dim=2)/1e3)
       DO j=1,nvm
          histvar(:,j)=veget_cov_max(:,j)*contfrac(:)*100
       ENDDO
     CALL xios_orchidee_send_field("landCoverFrac",histvar)
       vartmp(:)=zero
       DO j = 2,nvm
          IF (is_deciduous(j)) THEN
             vartmp(:) = vartmp(:) + veget_cov_max(:,j)*contfrac*100
          ENDIF
       ENDDO
     CALL xios_orchidee_send_field("treeFracPrimDec",vartmp)
       vartmp(:)=zero
       DO j = 2,nvm
          IF (is_evergreen(j)) THEN
             vartmp(:) = vartmp(:) + veget_cov_max(:,j)*contfrac*100
          ENDIF
       ENDDO
     CALL xios_orchidee_send_field("treeFracPrimEver",vartmp)
       vartmp(:)=zero
       DO j = 2,nvm
          IF ( .NOT.(is_c4(j)) ) THEN
             vartmp(:) = vartmp(:) + veget_cov_max(:,j)*contfrac*100
          ENDIF
       ENDDO
     CALL xios_orchidee_send_field("c3PftFrac",vartmp)
       vartmp(:)=zero
       DO j = 2,nvm
          IF ( is_c4(j) ) THEN
             vartmp(:) = vartmp(:) + veget_cov_max(:,j)*contfrac*100
          ENDIF
       ENDDO
     CALL xios_orchidee_send_field("c4PftFrac",vartmp)
     CALL xios_orchidee_send_field("rGrowth",SUM(resp_growth*veget_cov_max,dim=2)/1e3/one_day*contfrac)
     CALL xios_orchidee_send_field("rMaint",SUM(resp_maint*veget_cov_max,dim=2)/1e3/one_day*contfrac)

    CALL histwrite_p (hist_id_stomate, 'RESOLUTION_X', itime, &
         resolution(:,1), npts, hori_index)
    CALL histwrite_p (hist_id_stomate, 'RESOLUTION_Y', itime, &
         resolution(:,2), npts, hori_index)
    CALL histwrite_p (hist_id_stomate, 'CONTFRAC', itime, &
         contfrac(:), npts, hori_index)

    CALL histwrite_p (hist_id_stomate, 'DEADLEAF_COVER', itime, &
         deadleaf_cover, npts, hori_index)

    CALL histwrite_p (hist_id_stomate, 'TOTAL_SOIL_CARB', itime, &
         tot_litter_soil_carb, npts*nvm, horipft_index)

    CALL histwrite_p (hist_id_stomate, 'T2M_MONTH', itime, &
         t2m_month, npts, hori_index)
    CALL histwrite_p (hist_id_stomate, 'T2M_WEEK', itime, &
         t2m_week, npts, hori_index)
    CALL histwrite_p (hist_id_stomate, 'TSEASON', itime, &
         Tseason, npts, hori_index)
    CALL histwrite_p (hist_id_stomate, 'TMIN_SPRING_TIME', itime, &
         Tmin_spring_time, npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'ONSET_DATE', itime, &
         onset_date(:,:), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'FPC_MAX', itime, &
         fpc_max, npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'MAXFPC_LASTYEAR', itime, &
         maxfpc_lastyear, npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'HET_RESP', itime, &
         resp_hetero(:,:), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'FIREINDEX', itime, &
         fireindex(:,:), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'LITTERHUM', itime, &
         litterhum_daily, npts, hori_index)
    CALL histwrite_p (hist_id_stomate, 'CO2_FIRE', itime, &
         co2_fire, npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'CO2_TAKEN', itime, &
         co2_to_bm, npts*nvm, horipft_index)
    IF(ok_ncycle) CALL histwrite_p (hist_id_stomate, 'N_TAKEN', itime, &
         n_to_bm, npts*nvm, horipft_index)
    IF(ok_pcycle) CALL histwrite_p (hist_id_stomate, 'P_TAKEN', itime, &
         p_to_bm, npts*nvm, horipft_index)
    ! land cover change
    CALL histwrite_p (hist_id_stomate, 'CONVFLUX', itime, &
         convflux, npts, hori_index)
    CALL histwrite_p (hist_id_stomate, 'CFLUX_PROD10', itime, &
         cflux_prod10, npts, hori_index)
    CALL histwrite_p (hist_id_stomate, 'CFLUX_PROD100', itime, &
         cflux_prod100, npts, hori_index)
    CALL histwrite_p (hist_id_stomate, 'CONVFLUX_HARVEST', itime, &
         convflux_harvest, npts, hori_index)
    CALL histwrite_p (hist_id_stomate, 'CFLUX_PROD10_HARVEST', itime, &
         cflux_prod10_harvest, npts, hori_index)
    CALL histwrite_p (hist_id_stomate, 'CFLUX_PROD100_HARVES', itime, &
         cflux_prod100_harvest, npts, hori_index)
    CALL histwrite_p (hist_id_stomate, 'WOOD_HARVEST', itime, &
         woodharvest/one_year*dt_days, npts, hori_index)
    CALL histwrite_p (hist_id_stomate, 'WOOD_HARVEST_PFT', itime, &
         woodharvestpft, npts*nvm, horipft_index)

    CALL histwrite_p (hist_id_stomate, 'HARVEST_ABOVE', itime, &
         harvest_above, npts, hori_index)
    CALL histwrite_p (hist_id_stomate, 'GRM_DEVSTAGE', itime, &
         GRM_devstage, npts*nvm, horipft_index)

    CALL histwrite_p (hist_id_stomate, 'LAI', itime, &
         lai, npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'VEGET_COV_MAX', itime, &
         veget_cov_max, npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'NPP', itime, &
         npp_daily, npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'GPP', itime, &
         gpp_daily, npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'IND', itime, &
         ind, npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'CN_IND', itime, &
         cn_ind, npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'WOODMASS_IND', itime, &
         woodmass_ind, npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'MAINT_RESP', itime, &
         resp_maint, npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'GROWTH_RESP', itime, &
         resp_growth, npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'AGE', itime, &
         age, npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'HEIGHT', itime, &
         height, npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'MOISTRESS', itime, &
         moiavail_week, npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'VCMAX_ELLSWORTH', itime, &
         vcmax(:,:,1), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'VCMAX_KATTGE', itime, &
         vcmax_new, npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'TURNOVER_TIME', itime, &
         turnover_time, npts*nvm, horipft_index)
    ! land cover change
    CALL histwrite_p (hist_id_stomate, 'PROD10', itime, &
         prod10, npts*11, horip11_index)
    CALL histwrite_p (hist_id_stomate, 'PROD100', itime, &
         prod100, npts*101, horip101_index)
    CALL histwrite_p (hist_id_stomate, 'FLUX10', itime, &
         flux10, npts*10, horip10_index)
    CALL histwrite_p (hist_id_stomate, 'FLUX100', itime, &
         flux100, npts*100, horip100_index)
    ! gmjc
    CALL histwrite_p(hist_id_stomate ,'T2M_14'   ,itime, &
         t2m_14, npts, hori_index)
    CALL histwrite_p (hist_id_stomate, 'LITTER_STR_AVAIL', itime, &
         litter_avail(:,istructural,:), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'LITTER_MET_AVAIL', itime, &
         litter_avail(:,imetabolic,:), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'LITTER_STR_NAVAIL', itime, &
         litter_not_avail(:,istructural,:), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'LITTER_MET_NAVAIL', itime, &
         litter_not_avail(:,imetabolic,:), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'LITTER_STR_AVAILF', itime, &
         litter_avail_frac(:,istructural,:), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'LITTER_MET_AVAILF', itime, &
         litter_avail_frac(:,imetabolic,:), npts*nvm, horipft_index)
    ! end gmjc

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
      
       CALL histwrite_p (hist_id_stomate, 'LITTER_STR_AB'//TRIM(element_str(l)), itime, & 
            litter(:,istructural,:,iabove,l), npts*nvm, horipft_index) 
       CALL histwrite_p (hist_id_stomate, 'LITTER_MET_AB'//TRIM(element_str(l)), itime, & 
            litter(:,imetabolic,:,iabove,l), npts*nvm, horipft_index) 
       CALL histwrite_p (hist_id_stomate, 'LITTER_STR_BE'//TRIM(element_str(l)), itime, & 
            litter(:,istructural,:,ibelow,l), npts*nvm, horipft_index) 
       CALL histwrite_p (hist_id_stomate, 'LITTER_MET_BE'//TRIM(element_str(l)), itime, & 
            litter(:,imetabolic,:,ibelow,l), npts*nvm, horipft_index) 
       CALL histwrite_p (hist_id_stomate, 'LITTER_WOD_AB'//TRIM(element_str(l)), itime, & 
            litter(:,iwoody,:,iabove,l), npts*nvm, horipft_index) 
       CALL histwrite_p (hist_id_stomate, 'LITTER_WOD_BE'//TRIM(element_str(l)), itime, & 
            litter(:,iwoody,:,ibelow,l), npts*nvm, horipft_index) 
       CALL histwrite_p (hist_id_stomate, 'SOIL_ACTIVE'//TRIM(element_str(l)), itime, & 
            som(:,iactive,:,l), npts*nvm, horipft_index) 
       CALL histwrite_p (hist_id_stomate, 'SOIL_SLOW'//TRIM(element_str(l)), itime, & 
            som(:,islow,:,l), npts*nvm, horipft_index) 
       CALL histwrite_p (hist_id_stomate, 'SOIL_PASSIVE'//TRIM(element_str(l)), itime, & 
            som(:,ipassive,:,l), npts*nvm, horipft_index) 
       CALL histwrite_p (hist_id_stomate, 'SOIL_SURF'//TRIM(element_str(l)), itime, & 
            som(:,isurface,:,l), npts*nvm, horipft_index) 
       CALL histwrite_p (hist_id_stomate, 'TOTAL_M'//TRIM(element_str(l)), itime, & 
            tot_live_biomass(:,:,l), npts*nvm, horipft_index) 
       CALL histwrite_p (hist_id_stomate, 'LEAF_M'//TRIM(element_str(l)), itime, & 
            biomass(:,:,ileaf,l), npts*nvm, horipft_index) 
       CALL histwrite_p (hist_id_stomate, 'SAP_M_AB'//TRIM(element_str(l)), itime, & 
            biomass(:,:,isapabove,l), npts*nvm, horipft_index) 
       CALL histwrite_p (hist_id_stomate, 'SAP_M_BE'//TRIM(element_str(l)), itime, & 
            biomass(:,:,isapbelow,l), npts*nvm, horipft_index) 
       CALL histwrite_p (hist_id_stomate, 'HEART_M_AB'//TRIM(element_str(l)), itime, & 
            biomass(:,:,iheartabove,l), npts*nvm, horipft_index) 
       CALL histwrite_p (hist_id_stomate, 'HEART_M_BE'//TRIM(element_str(l)), itime, & 
            biomass(:,:,iheartbelow,l), npts*nvm, horipft_index) 
       CALL histwrite_p (hist_id_stomate, 'ROOT_M'//TRIM(element_str(l)), itime, & 
            biomass(:,:,iroot,l), npts*nvm, horipft_index) 
       CALL histwrite_p (hist_id_stomate, 'FRUIT_M'//TRIM(element_str(l)), itime, & 
            biomass(:,:,ifruit,l), npts*nvm, horipft_index) 
       CALL histwrite_p (hist_id_stomate, 'RESERVE_M'//TRIM(element_str(l)), itime, & 
            biomass(:,:,icarbres,l), npts*nvm, horipft_index) 
       CALL histwrite_p (hist_id_stomate, 'LABILE_M'//TRIM(element_str(l)), itime, & 
            biomass(:,:,ilabile,l), npts*nvm, horipft_index) 
       CALL histwrite_p (hist_id_stomate, 'TOTAL_TURN'//TRIM(element_str(l)), itime, & 
            tot_turnover(:,:,l), npts*nvm, horipft_index) 
       CALL histwrite_p (hist_id_stomate, 'LEAF_TURN'//TRIM(element_str(l)), itime, & 
            turnover_daily(:,:,ileaf,l), npts*nvm, horipft_index) 
       CALL histwrite_p (hist_id_stomate, 'SAP_AB_TURN'//TRIM(element_str(l)), itime, & 
            turnover_daily(:,:,isapabove,l), npts*nvm, horipft_index) 
       CALL histwrite_p (hist_id_stomate, 'ROOT_TURN'//TRIM(element_str(l)), itime, & 
            turnover_daily(:,:,iroot,l), npts*nvm, horipft_index) 
       CALL histwrite_p (hist_id_stomate, 'FRUIT_TURN'//TRIM(element_str(l)), itime, & 
            turnover_daily(:,:,ifruit,l), npts*nvm, horipft_index) 

       CALL histwrite_p (hist_id_stomate, 'LEAF_RECYCLE'//TRIM(element_str(l)), itime, & 
            recycling_daily(:,:,ileaf,l), npts*nvm, horipft_index) 
       CALL histwrite_p (hist_id_stomate, 'ROOT_RECYCLE'//TRIM(element_str(l)), itime, & 
            recycling_daily(:,:,iroot,l), npts*nvm, horipft_index) 
       CALL histwrite_p (hist_id_stomate, 'TOTAL_BM_LITTER'//TRIM(element_str(l)), itime, & 
            tot_bm_to_litter(:,:,l), npts*nvm, horipft_index) 
       CALL histwrite_p (hist_id_stomate, 'LEAF_BM_LITTER'//TRIM(element_str(l)), itime, & 
            bm_to_litter(:,:,ileaf,l), npts*nvm, horipft_index) 
       CALL histwrite_p (hist_id_stomate, 'SAP_AB_BM_LITTER'//TRIM(element_str(l)), itime, & 
            bm_to_litter(:,:,isapabove,l), npts*nvm, horipft_index) 
       CALL histwrite_p (hist_id_stomate, 'SAP_BE_BM_LITTER'//TRIM(element_str(l)), itime, & 
            bm_to_litter(:,:,isapbelow,l), npts*nvm, horipft_index) 
       CALL histwrite_p (hist_id_stomate, 'HEART_AB_BM_LITTER'//TRIM(element_str(l)), itime, & 
            bm_to_litter(:,:,iheartabove,l), npts*nvm, horipft_index) 
       CALL histwrite_p (hist_id_stomate, 'HEART_BE_BM_LITTER'//TRIM(element_str(l)), itime, & 
            bm_to_litter(:,:,iheartbelow,l), npts*nvm, horipft_index) 
       CALL histwrite_p (hist_id_stomate, 'ROOT_BM_LITTER'//TRIM(element_str(l)), itime, & 
            bm_to_litter(:,:,iroot,l), npts*nvm, horipft_index) 
       CALL histwrite_p (hist_id_stomate, 'FRUIT_BM_LITTER'//TRIM(element_str(l)), itime, & 
            bm_to_litter(:,:,ifruit,l), npts*nvm, horipft_index) 
       CALL histwrite_p (hist_id_stomate, 'RESERVE_BM_LITTER'//TRIM(element_str(l)), itime, & 
            bm_to_litter(:,:,icarbres,l), npts*nvm, horipft_index) 
       CALL histwrite_p (hist_id_stomate, 'LABILE_BM_LITTER'//TRIM(element_str(l)), itime, &
            bm_to_litter(:,:,ilabile,l), npts*nvm, horipft_index)
    ENDDO

    IF (impose_cn .AND. ok_ncycle) THEN
      CALL xios_orchidee_send_field('N_SUPPORT',N_support(:,:))
      CALL histwrite_p (hist_id_stomate, 'N_SUPPORT', itime, &
         N_support(:,:), npts*nvm, horipft_index)
    ENDIF

    IF (impose_np .AND. ok_pcycle) THEN
      CALL xios_orchidee_send_field('P_SUPPORT',P_support(:,:))
      CALL histwrite_p (hist_id_stomate, 'P_SUPPORT', itime, &
         P_support(:,:), npts*nvm, horipft_index)
    ENDIF

    IF ( hist_id_stomate_IPCC > 0 ) THEN
       vartmp(:)=SUM(tot_live_biomass(:,:,icarbon)*veget_cov_max,dim=2)/1e3
       CALL histwrite_p (hist_id_stomate_IPCC, "cVeg", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM(tot_litter_carb*veget_cov_max,dim=2)/1e3
       CALL histwrite_p (hist_id_stomate_IPCC, "cLitter", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM(tot_soil_carb*veget_cov_max,dim=2)/1e3
       CALL histwrite_p (hist_id_stomate_IPCC, "cSoil", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=(prod10_total + prod100_total + prod10_harvest_total + prod100_harvest_total)/1e3
       CALL histwrite_p (hist_id_stomate_IPCC, "cProduct", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=carb_mass_variation/1e3/one_day
       CALL histwrite_p (hist_id_stomate_IPCC, "cMassVariation", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM(lai*veget_cov_max,dim=2)
       CALL histwrite_p (hist_id_stomate_IPCC, "lai", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM(gpp_daily*veget_cov_max,dim=2)/1e3/one_day
       CALL histwrite_p (hist_id_stomate_IPCC, "gpp", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM((resp_maint+resp_growth)*veget_cov_max,dim=2)/1e3/one_day
       CALL histwrite_p (hist_id_stomate_IPCC, "ra", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM(npp_daily*veget_cov_max,dim=2)/1e3/one_day
       CALL histwrite_p (hist_id_stomate_IPCC, "npp", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM(resp_hetero*veget_cov_max,dim=2)/1e3/one_day
       CALL histwrite_p (hist_id_stomate_IPCC, "rh", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM(co2_fire*veget_cov_max,dim=2)/1e3/one_day
       CALL histwrite_p (hist_id_stomate_IPCC, "fFire", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=harvest_above/1e3/one_day
       CALL histwrite_p (hist_id_stomate_IPCC, "fHarvest", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=cflux_prod_total/1e3/one_day
       CALL histwrite_p (hist_id_stomate_IPCC, "fLuc", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=cflux_prod_harvest_total/1e3/one_day
       CALL histwrite_p (hist_id_stomate_IPCC, "fWoodharvest", itime, &
            vartmp, npts, hori_index)

       IF (allow_agri_fert) THEN
       vartmp(:)=(SUM((gpp_daily+co2_to_bm-(resp_maint+resp_growth+resp_hetero)-co2_fire + &
            &        fert_masschangecnp(:,:,icarbon)) &
            &        *veget_cov_max,dim=2)-cflux_prod_total-cflux_prod_harvest_total- &
            &         harvest_above - grazed_above)/1e3/one_day
       ELSE
       vartmp(:)=(SUM((gpp_daily+co2_to_bm-(resp_maint+resp_growth+resp_hetero)-co2_fire) &
            &        *veget_cov_max,dim=2)-cflux_prod_total-cflux_prod_harvest_total-harvest_above)/1e3/one_day
       ENDIF

       CALL histwrite_p (hist_id_stomate_IPCC, "nbp", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM((tot_bm_to_litter(:,:,icarbon) + tot_turnover(:,:,icarbon))*veget_cov_max,dim=2)/1e3/one_day
       CALL histwrite_p (hist_id_stomate_IPCC, "fVegLitter", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM(SUM(som_input(:,:,:,icarbon),dim=2)*veget_cov_max,dim=2)/1e3/one_day
       CALL histwrite_p (hist_id_stomate_IPCC, "fLitterSoil", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM(biomass(:,:,ileaf,icarbon)*veget_cov_max,dim=2)/1e3
       CALL histwrite_p (hist_id_stomate_IPCC, "cLeaf", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM((biomass(:,:,isapabove,icarbon)+biomass(:,:,iheartabove,icarbon))*veget_cov_max,dim=2)/1e3
       CALL histwrite_p (hist_id_stomate_IPCC, "cWood", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM(( biomass(:,:,iroot,icarbon) + biomass(:,:,isapbelow,icarbon) + biomass(:,:,iheartbelow,icarbon) ) &
            &        *veget_cov_max,dim=2)/1e3
       CALL histwrite_p (hist_id_stomate_IPCC, "cRoot", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM(( biomass(:,:,icarbres,icarbon) +  biomass(:,:,ilabile,icarbon) + biomass(:,:,ifruit,icarbon))*veget_cov_max,dim=2)/1e3
       CALL histwrite_p (hist_id_stomate_IPCC, "cMisc", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM((litter(:,istructural,:,iabove,icarbon)+litter(:,imetabolic,:,iabove,icarbon))*&
            veget_cov_max,dim=2)/1e3
       CALL histwrite_p (hist_id_stomate_IPCC, "cLitterAbove", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM((litter(:,istructural,:,ibelow,icarbon)+litter(:,imetabolic,:,ibelow,icarbon))*&
            veget_cov_max,dim=2)/1e3
       CALL histwrite_p (hist_id_stomate_IPCC, "cLitterBelow", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM((som(:,iactive,:,icarbon)+som(:,isurface,:,icarbon))*veget_cov_max,dim=2)/1e3
       CALL histwrite_p (hist_id_stomate_IPCC, "cSoilFast", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM(som(:,islow,:,icarbon)*veget_cov_max,dim=2)/1e3
       CALL histwrite_p (hist_id_stomate_IPCC, "cSoilMedium", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM(som(:,ipassive,:,icarbon)*veget_cov_max,dim=2)/1e3
       CALL histwrite_p (hist_id_stomate_IPCC, "cSoilSlow", itime, &
            vartmp, npts, hori_index)
       DO j=1,nvm
          histvar(:,j)=veget_cov_max(:,j)*100
       ENDDO
       CALL histwrite_p (hist_id_stomate_IPCC, "landCoverFrac", itime, &
            histvar, npts*nvm, horipft_index)
       !-
       vartmp(:)=zero
       DO j = 2,nvm
          IF (is_deciduous(j)) THEN
             vartmp(:) = vartmp(:) + veget_cov_max(:,j)*100
          ENDIF
       ENDDO
       CALL histwrite_p (hist_id_stomate_IPCC, "treeFracPrimDec", itime, &
            vartmp, npts, hori_index)
       !-
       vartmp(:)=zero
       DO j = 2,nvm
          IF (is_evergreen(j)) THEN
             vartmp(:) = vartmp(:) + veget_cov_max(:,j)*100
          ENDIF
       ENDDO
       CALL histwrite_p (hist_id_stomate_IPCC, "treeFracPrimEver", itime, &
            vartmp, npts, hori_index)
       !-
       vartmp(:)=zero
       DO j = 2,nvm
          IF ( .NOT.(is_c4(j)) ) THEN
             vartmp(:) = vartmp(:) + veget_cov_max(:,j)*100
          ENDIF
       ENDDO
       CALL histwrite_p (hist_id_stomate_IPCC, "c3PftFrac", itime, &
            vartmp, npts, hori_index)
       !-
       vartmp(:)=zero
       DO j = 2,nvm
          IF ( is_c4(j) ) THEN
             vartmp(:) = vartmp(:) + veget_cov_max(:,j)*100
          ENDIF
       ENDDO
       CALL histwrite_p (hist_id_stomate_IPCC, "c4PftFrac", itime, &
            vartmp, npts, hori_index)
       !-
       vartmp(:)=SUM(resp_growth*veget_cov_max,dim=2)/1e3/one_day
       CALL histwrite_p (hist_id_stomate_IPCC, "rGrowth", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM(resp_maint*veget_cov_max,dim=2)/1e3/one_day
       CALL histwrite_p (hist_id_stomate_IPCC, "rMaint", itime, &
            vartmp, npts, hori_index)
       !DSG vartmp(:)=SUM(bm_alloc(:,:,ileaf,icarbon)*veget_cov_max,dim=2)/1e3/one_day
       !DSG CALL histwrite_p (hist_id_stomate_IPCC, "nppLeaf", itime, &
       !DSG      vartmp, npts, hori_index)
       !DSG vartmp(:)=SUM(bm_alloc(:,:,isapabove,icarbon)*veget_cov_max,dim=2)/1e3/one_day
       !DSG CALL histwrite_p (hist_id_stomate_IPCC, "nppWood", itime, &
       !DSG      vartmp, npts, hori_index)
       !DSG vartmp(:)=SUM(( bm_alloc(:,:,isapbelow,icarbon) + bm_alloc(:,:,iroot,icarbon) )*veget_cov_max,dim=2)/1e3/one_day
       !DSG CALL histwrite_p (hist_id_stomate_IPCC, "nppRoot", itime, &
       !DSG      vartmp, npts, hori_index)

       CALL histwrite_p (hist_id_stomate_IPCC, 'RESOLUTION_X', itime, &
            resolution(:,1), npts, hori_index)
       CALL histwrite_p (hist_id_stomate_IPCC, 'RESOLUTION_Y', itime, &
            resolution(:,2), npts, hori_index)
       CALL histwrite_p (hist_id_stomate_IPCC, 'CONTFRAC', itime, &
            contfrac(:), npts, hori_index)

    ENDIF

    IF (printlev>=3) WRITE (numout,*) 'soil_n_min(test_grid,test_pft,:):',soil_n_min(test_grid,test_pft,:)
    IF (printlev>=4) WRITE(numout,*) 'Leaving stomate_lpj'

  END SUBROUTINE StomateLpj

!! ================================================================================================================================
!! SUBROUTINE   : harvest
!!
!>\BRIEF        Harvest of croplands
!!
!! DESCRIPTION  : To take into account biomass harvest from crop (mainly to take 
!! into account for the reduced litter input and then decreased soil carbon. it is a 
!! constant (40\%) fraction of above ground biomass. Comparable concept is used
!! harvest biomass from bioenergy grass plantation
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S): ::harvest_above the harvested biomass
!!
!! REFERENCE(S) :
!! - Piao, S., P. Ciais, P. Friedlingstein, N. de Noblet-Ducoudre, P. Cadule, N. Viovy, and T. Wang. 2009. 
!!   Spatiotemporal patterns of terrestrial carbon cycle during the 20th century. Global Biogeochemical 
!!   Cycles 23:doi:10.1029/2008GB003339.
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE harvest(npts, dt_days, veget_cov_max, &
       bm_to_litter, turnover_daily, &
       harvest_above, harvest_masschangecnp)

  !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER, INTENT(in)                                    :: npts            !! Domain size (unitless) 
    REAL, INTENT(in)                                :: dt_days         !! Time step (days)                               
    REAL, DIMENSION(npts,nvm), INTENT(in)           :: veget_cov_max       !! new "maximal" coverage fraction of a PFT (LAI -> 
                                                                              !! infinity) on ground @tex $(m^2 m^{-2})$ @endtex 
    
   !! 0.2 Output variables
   
   !! 0.3 Modified variables

    REAL, DIMENSION(npts,nvm,nparts,nelements), INTENT(inout) :: bm_to_litter !! [DISPENSABLE] conversion of biomass to litter 
                                                                                     !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL, DIMENSION(npts,nvm,nparts,nelements), INTENT(inout) :: turnover_daily   !! Turnover rates 
                                                                               !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL, DIMENSION(npts), INTENT(inout)            :: harvest_above    !! harvest above ground biomass for agriculture 
                                                                               !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL, DIMENSION(npts,nvm,nelements), INTENT(inout) :: harvest_masschangecnp   !! crop harvest mass change 
    !! 0.4 Local variables

    INTEGER                                         :: i, j, k, l, m    !! indices                       
    REAL                                            :: above_old        !! biomass of previous time step 
                                                                               !! @tex $(gC m^{-2})$ @endtex 
    CHARACTER(LEN=2), DIMENSION(nelements)                     :: element_str         !! string suffix indicating element
!_ ================================================================================================================================

  !! 1. Yearly initialisation

    above_old             = zero
    harvest_above         = zero
    harvest_masschangecnp = zero

    DO i = 1, npts
       DO j = 1,nvm
       !DSG: why is below ground biomass harvested? 
       !DSG: make sure bioenergy tiles are not harvested:
          IF ((.NOT. natural(j)).AND.(.NOT.bioenergy(j))) THEN
             above_old = turnover_daily(i,j,ileaf,icarbon) + turnover_daily(i,j,isapabove,icarbon) + &
                  &       turnover_daily(i,j,iheartabove,icarbon) + turnover_daily(i,j,ifruit,icarbon) + &
                  &       turnover_daily(i,j,icarbres,icarbon) + turnover_daily(i,j,ilabile,icarbon) + turnover_daily(i,j,isapbelow,icarbon) + &
                  &       turnover_daily(i,j,iheartbelow,icarbon) + turnover_daily(i,j,iroot,icarbon)

             !JC MOD037 harvest for icarbon, initrogen and iphosphorus at the same time
             ! and record harvested cnp
             harvest_masschangecnp(i,j,:) = (un - frac_turnover_daily) * &
                 SUM(turnover_daily(i,j,:,:),DIM=1)

             turnover_daily(i,j,ileaf,:)       = turnover_daily(i,j,ileaf,:)*frac_turnover_daily
             turnover_daily(i,j,isapabove,:)   = turnover_daily(i,j,isapabove,:)*frac_turnover_daily
             turnover_daily(i,j,isapbelow,:)   = turnover_daily(i,j,isapbelow,:)*frac_turnover_daily
             turnover_daily(i,j,iheartabove,:) = turnover_daily(i,j,iheartabove,:)*frac_turnover_daily
             turnover_daily(i,j,iheartbelow,:) = turnover_daily(i,j,iheartbelow,:)*frac_turnover_daily
             turnover_daily(i,j,iroot,:)       = turnover_daily(i,j,iroot,:)*frac_turnover_daily
             turnover_daily(i,j,ifruit,:)      = turnover_daily(i,j,ifruit,:)*frac_turnover_daily
             turnover_daily(i,j,icarbres,:)    = turnover_daily(i,j,icarbres,:)*frac_turnover_daily
             turnover_daily(i,j,ilabile,:)     = turnover_daily(i,j,ilabile,:)*frac_turnover_daily
             harvest_above(i)  = harvest_above(i) + veget_cov_max(i,j) * above_old *(un - frac_turnover_daily)

          ENDIF
       ENDDO
    ENDDO


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
       CALL xios_orchidee_send_field("CROP_HARVEST"//TRIM(element_str(l)),harvest_masschangecnp(:,:,l))
       CALL histwrite_p (hist_id_stomate, 'CROP_HARVEST'//TRIM(element_str(l)), itime, &
            harvest_masschangecnp(:,:,l),npts*nvm, horipft_index)
    END DO

!!$    harvest_above = harvest_above
  END SUBROUTINE harvest


!! ================================================================================================================================
!! SUBROUTINE   : harvest_bioenergygrass
!!
!>\BRIEF        Harvest of bioenergy grass crops
!!
!! DESCRIPTION  : To take into account biomass harvest from crop (mainly to take 
!! into account for the reduced litter input and then decreased soil carbon. it is a 
!! constant  fraction of above ground biomass without storage (which is assumed belowground).
!!
!! RECENT CHANGE(S) : Bioenergy grass harvest
!!
!! MAIN OUTPUT VARIABLE(S): ::harvest_above the harvested biomass
!!
!! REFERENCE(S) :
!! - Wei Li
!!   
!!   
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE harvest_biograss(npts, dt_days, veget_cov_max, &
       bm_to_litter, turnover_daily, &
       harvest_bioenergy, harvest_BE_masschangecnp)

  !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER, INTENT(in)                                    :: npts            !! Domain size (unitless) 
    REAL, INTENT(in)                                :: dt_days         !! Time step (days)                               
    REAL, DIMENSION(npts,nvm), INTENT(in)           :: veget_cov_max       !! new "maximal" coverage fraction of a PFT (LAI -> 
                                                                              !! infinity) on ground @tex $(m^2 m^{-2})$ @endtex 
    
   !! 0.2 Output variables
   
   !! 0.3 Modified variables

    REAL, DIMENSION(npts,nvm,nparts,nelements), INTENT(inout) :: bm_to_litter !! [DISPENSABLE] conversion of biomass to litter 
                                                                                     !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL, DIMENSION(npts,nvm,nparts,nelements), INTENT(inout) :: turnover_daily   !! Turnover rates 
                                                                               !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL, DIMENSION(npts,nvm), INTENT(inout)            :: harvest_bioenergy !! harvest above ground biomass for bioenergy
                                                                               !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL, DIMENSION(npts,nvm,nelements), INTENT(inout) :: harvest_BE_masschangecnp   !! bioenergy harvest mass change 
    !! 0.4 Local variables

    INTEGER                                         :: i, j, k, l, m    !! indices                       
    REAL                                            :: turnover_above_old        !! dummy
                                                                               !! @tex $(gC m^{-2})$ @endtex 
    CHARACTER(LEN=2), DIMENSION(nelements)                     :: element_str         !! string suffix indicating element
!_ ================================================================================================================================

  !! 1. Yearly initialisation

    turnover_above_old    = zero
    harvest_BE_masschangecnp = zero
    ! DSG: not sure about setting it here to zero:
    harvest_bioenergy     = zero 

    DO i = 1, npts
       DO j = 1,nvm
          IF ( bioenergy(j) .AND. (.NOT. is_tree(j) ) ) THEN
             ! Harvest only AGB w/o storage which is assumed belowground (e.g. rhizomes)
             turnover_above_old = turnover_daily(i,j,ileaf,icarbon) + turnover_daily(i,j,isapabove,icarbon) + &
                  &       turnover_daily(i,j,iheartabove,icarbon) + turnover_daily(i,j,ifruit,icarbon) 
             ! compute harvested carbon only mass:
             harvest_bioenergy(i,j)  = harvest_bioenergy(i,j) + veget_cov_max(i,j) * turnover_above_old * frac_bioenergy_harvest(j)

             ! compute turnover of CNP biomass due to bioenergy harvest:
             turnover_daily(i,j,ileaf,:)       = turnover_daily(i,j,ileaf,:)      * (un - frac_bioenergy_harvest(j))
             turnover_daily(i,j,isapabove,:)   = turnover_daily(i,j,isapabove,:)  * (un - frac_bioenergy_harvest(j))
             turnover_daily(i,j,iheartabove,:) = turnover_daily(i,j,iheartabove,:)* (un - frac_bioenergy_harvest(j))
             turnover_daily(i,j,ifruit,:)      = turnover_daily(i,j,ifruit,:)     * (un - frac_bioenergy_harvest(j))
             ! sum CNP biomass turnover up for diagnostic output
             harvest_BE_masschangecnp(i,j,:) = turnover_daily(i,j,ileaf,:) + turnover_daily(i,j,isapabove,:) &
                                               + turnover_daily(i,j,iheartabove,:) + turnover_daily(i,j,ifruit,:)
          ENDIF
       ENDDO
    ENDDO

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
       ! CALL xios_orchidee_send_field("GRASS_BE_HARVEST"//TRIM(element_str(l)),harvest_BE_masschangecnp(:,:,l))
       ! CALL histwrite_p (hist_id_stomate, 'GRASS_BE_HARVEST'//TRIM(element_str(l)), itime, &
       !      harvest_BE_masschangecnp(:,:,l),npts*nvm, horipft_index)
    END DO

!!$    harvest_above = harvest_above
  END SUBROUTINE harvest_biograss



  SUBROUTINE simplified_grazing(npts, dt_days, day_of_year, veget_cov_max, &
       agri_fert_save, turnover_daily, rest_Ngrazing, grazed_above)

    INTEGER, INTENT(in)                                    :: npts
    REAL, INTENT(in)                                :: dt_days
    INTEGER                    , INTENT(in)   :: day_of_year
    REAL, DIMENSION(npts,nvm), INTENT(in)           :: veget_cov_max
    REAL, DIMENSION(npts,nvm,5), INTENT(in)   :: agri_fert_save
    REAL, DIMENSION(npts,nvm,nparts,nelements), INTENT(inout) :: turnover_daily
    REAL, DIMENSION(npts,nvm), INTENT(inout) :: rest_Ngrazing
    REAL, DIMENSION(npts), INTENT(inout) :: grazed_above

    INTEGER                                         :: i, j, k, l, m
    CHARACTER(LEN=2), DIMENSION(nelements)                     :: element_str
    REAL                                 :: above_old
    REAL                                 :: above_new
    REAL, DIMENSION(npts,nvm)        :: tmp_Ngrazing
    REAL, DIMENSION(npts,nvm)        :: record_deficit
    REAL, DIMENSION(npts,nvm,nelements)  :: grazedcnp

    above_old = zero
    above_new = zero
    grazed_above = zero
    grazedcnp(:,:,:) = zero
    tmp_Ngrazing(:,:)  = agri_fert_save(:,:,4)
    record_deficit(:,:) = zero

    IF (REAL(day_of_year) .GT. 0.5 .AND. REAL(day_of_year) .LT. 1.5) THEN
      rest_Ngrazing(:,:) = tmp_Ngrazing(:,:)
      WHERE (veget_cov_max .LE. min_stomate)
        rest_Ngrazing(:,:) = zero
      ENDWHERE
    ENDIF

    DO i = 1, npts
       DO j = 1,nvm
          IF (is_agri_fert(j) .AND. (pasture(j) .OR. rangeland(j))) THEN
            ! only when there is still N to export as grazing manure      
            IF (veget_cov_max(i,j) .GT. min_stomate .AND. &
                tmp_Ngrazing(i,j) .GT. zero .AND. rest_Ngrazing(i,j) .GT. zero) THEN
              ! only above ground biomass can be grazed
              above_old = turnover_daily(i,j,ileaf,initrogen) + turnover_daily(i,j,isapabove,initrogen)
              IF (above_old .GT. min_stomate) THEN
                ! rest_Ngrazing can be reduced to zero
                rest_Ngrazing(i,j) = rest_Ngrazing(i,j) - MIN(above_old,rest_Ngrazing(i,j))
                above_new = above_old - MIN(above_old,rest_Ngrazing(i,j))
                DO l=1,2 ! only for carbon and nitrogen at current stage
                  grazedcnp(i,j,l) = turnover_daily(i,j,ileaf,l) * &
                     (un - above_new/above_old)
                  turnover_daily(i,j,ileaf,l) = turnover_daily(i,j,ileaf,l)*above_new/above_old
                  turnover_daily(i,j,isapabove,l) = turnover_daily(i,j,isapabove,l)*above_new/above_old  
                ENDDO
                grazed_above(i) = grazed_above(i) + veget_cov_max(i,j) * grazedcnp(i,j,icarbon)
              ENDIF
            ENDIF
          ENDIF
       ENDDO
    ENDDO
    ! last day of a year record the deficit N
    IF (REAL(day_of_year) .GT. 363.5 .AND. REAL(day_of_year) .LT. 364.5) THEN
      record_deficit(:,:) = rest_Ngrazing(:,:)
    ENDIF

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
       CALL xios_orchidee_send_field("GRASS_GRAZED"//TRIM(element_str(l)),grazedcnp(:,:,l))
!       CALL histwrite_p (hist_id_stomate, 'GRASS_GRAZED'//TRIM(element_str(l)), itime, &
!            grazedcnp(:,:,l),npts*nvm, horipft_index)
       CALL xios_orchidee_send_field("GNDEPOT_DEFICIT",record_deficit(:,:))
    END DO


  END SUBROUTINE simplified_grazing

  SUBROUTINE calc_GRM_devstage(npts, &
               t2m_daily, t2m_week, &
               tsurf_daily, begin_leaves, GRM_devstage)

    INTEGER (i_std)                   , INTENT(in)  :: npts
    REAL, DIMENSION(npts), INTENT(in)  :: t2m_daily
    REAL, DIMENSION(npts), INTENT(in)  :: t2m_week
    REAL, DIMENSION(npts), INTENT(in)  :: tsurf_daily
    LOGICAL, DIMENSION(npts,nvm), INTENT(in) :: begin_leaves
    REAL, DIMENSION(npts,nvm), INTENT(inout) :: GRM_devstage

    INTEGER                                         :: i, j

    DO j=2,nvm
      DO i=1,npts
        ! reset devstage at the beginning of growing season
        IF (begin_leaves(i,j) .EQV. .TRUE.) THEN
          GRM_devstage(i,j) = MAX(0.0, t2m_daily(i) - tbase)/tasumrep

        ! when tsoil > tbase = 5oC 
        ELSE IF ((GRM_devstage(i,j) .GE. 0.0) .AND. &
                 (tsurf_daily(i) .GT. tbase) .AND. &
                 (GRM_devstage(i,j) .LT. 2.0))  THEN

          GRM_devstage(i,j) = GRM_devstage(i,j) + MAX(0.0, t2m_daily(i) - &
                              tbase)/tasumrep
        ELSE 
          GRM_devstage(i,j) = GRM_devstage(i,j)
        ENDIF

      ENDDO
    ENDDO 

    WHERE (GRM_devstage .GT. 2.0) 
      GRM_devstage(:,:) = 2.0
    ENDWHERE

  END SUBROUTINE calc_GRM_devstage

  SUBROUTINE apply_agri_fert(npts,dt_days,day_of_year,veget_cov_max, &
               agri_fert_save,litter,som,agri_nfert,agri_pfert,agri_masschangecnp)

    INTEGER (i_std)                   , INTENT(in)   :: npts
    REAL                       , INTENT(in)   :: dt_days
    INTEGER                    , INTENT(in)   :: day_of_year
    REAL, DIMENSION(npts,nvm), INTENT(in)  :: veget_cov_max
    REAL, DIMENSION(npts,nvm,5), INTENT(in)   :: agri_fert_save
    REAL, DIMENSION(npts,nlitt,nvm,nlevs,nelements), INTENT(inout)  :: litter
    REAL, DIMENSION(npts,ncarb,nvm,nelements), INTENT(inout)      :: som
    REAL, DIMENSION(npts,nvm,ninput), INTENT(out)  :: agri_nfert
    REAL, DIMENSION(npts,nvm), INTENT(out)  :: agri_pfert
    REAL, DIMENSION(npts,nvm,nelements), INTENT(out)  :: agri_masschangecnp

    INTEGER                                         :: i, j
    REAL, DIMENSION(npts,nvm)        :: tmp_Nmineral
    REAL, DIMENSION(npts,nvm)        :: tmp_Nmanure
    REAL, DIMENSION(npts,nvm)        :: tmp_Pmineral
    REAL, DIMENSION(npts,nvm)        :: tmp_Ngrazing
    REAL, DIMENSION(npts,nvm)        :: tmp_Tfert 
    REAL, DIMENSION(npts,nvm)        :: real_Nmineral
    REAL, DIMENSION(npts,nvm)        :: real_Nmanure
    REAL, DIMENSION(npts,nvm)        :: real_Ngrazing
    REAL, DIMENSION(npts,nvm)        :: real_Pmineral
    REAL, PARAMETER       :: org_furine = 0.9
    REAL, PARAMETER       :: org_fmetabolic = 0.9
    REAL, PARAMETER       :: org_cnratio = 10.0
    REAL, PARAMETER       :: org_pnratio = 0.2
    REAL, PARAMETER       :: min_famm = 0.5
    REAL, PARAMETER       :: fresh_fmetabolic = 0.3
    REAL, PARAMETER       :: fresh_fisurface = 0.3431
    REAL, PARAMETER       :: fresh_furine = 0.5
    REAL, PARAMETER       :: fresh_depdays = 90.
    REAL, PARAMETER       :: rcnrr = 100.
    REAL, PARAMETER       :: rcnh = 10.
    REAL, PARAMETER       :: rcnfeces = 16.
    REAL :: tmp_clitter
    REAL :: tmp_csoil
    REAL :: tmp_nlitter
    REAL :: tmp_nsoil
    REAL :: tmp_plitter
    REAL :: tmp_psoil
    ! save fertilizer input to temporary variable
    tmp_Nmineral(:,:) = agri_fert_save(:,:,1)
    tmp_Nmanure(:,:)  = agri_fert_save(:,:,2)
    tmp_Pmineral(:,:) = agri_fert_save(:,:,3)
    tmp_Ngrazing(:,:)  = agri_fert_save(:,:,4)
    tmp_Tfert(:,:)    = agri_fert_save(:,:,5)
    ! initialize output variables
    agri_nfert(:,:,:) = zero
    agri_pfert(:,:) = zero
    agri_masschangecnp(:,:,:) = zero    
    real_Nmineral(:,:) = zero
    real_Nmanure(:,:) = zero
    real_Ngrazing(:,:) = zero
    real_Pmineral(:,:) = zero
    tmp_clitter = zero
    tmp_csoil = zero
    tmp_nlitter = zero
    tmp_nsoil = zero
    tmp_plitter = zero
    tmp_psoil = zero

    DO j=2,nvm
      DO i=1,npts
        ! for grazing deposition of N
        ! I assumed a even deposition for 90 days rather than once
        ! to mimic the grazing activity, though it should only be applied during
        ! growing/grazing season
        ! Tfert is set to DOY 120 180 240, thus Tfert+90 will still be within
        ! the same year
          tmp_clitter = zero
          tmp_csoil = zero
          tmp_nlitter = zero
          tmp_nsoil = zero
          tmp_plitter = zero
          tmp_psoil = zero
        IF (is_agri_fert(j) .AND. (pasture(j) .OR. rangeland(j)) .AND. &
            veget_cov_max(i,j) .GT. min_stomate .AND. tmp_Ngrazing(i,j) .GT. zero .AND. &
            (REAL(day_of_year) .GT. tmp_Tfert(i,j)-0.5) .AND. &
            (REAL(day_of_year) .LT. tmp_Tfert(i,j)+fresh_depdays-0.5)) THEN 
          ! record real deposition
          real_Ngrazing(i,j) = tmp_Ngrazing(i,j)/fresh_depdays
          ! 50% of fresh excreta as urine N
          ! urine C P is negligible
          agri_nfert(i,j,iammonium) = agri_nfert(i,j,iammonium) + &
             real_Ngrazing(i,j) * fresh_furine
          ! the rest (50%) of fresh excreta as waste feces
          ! and it will goes to litter and soil organic pool
          ! following Appendix 1 Eq.1 of Li et al., 2012
          ! Here, whit rcnfeces = 16
          ! tmp_nsoil is 84/90*fecesN tmp_nlitter is 6/90*fecesN
          ! N:P ratio is 10:3 see below
          tmp_nsoil = (rcnrr - rcnfeces)*(un-fresh_furine)*real_Ngrazing(i,j) / &
                       (rcnrr-rcnh)
          tmp_nlitter = (un - (rcnrr - rcnfeces)/(rcnrr-rcnh)) * &
                       (un-fresh_furine)*real_Ngrazing(i,j)
          tmp_csoil = rcnh*tmp_nsoil
          tmp_clitter = rcnrr*tmp_nlitter
          ! At this stage, I don't add P for manure depositon, since all the P comes from
          ! grazed biomass. To avoid over deposition (NP ratio in manure is lower than in
          ! biomass), I assumed grazed/depoted P is the same and do not change any
          ! turnover. But for N, it accelated the turnover (from biomass to urine and
          ! feces) and change the CN ratio to 8:1 (C respired by animal). In this case, we should
          ! be able to partly capture the grassland management effect on nutrient cycles,
          ! before the GRM module become fully functioning.
!          tmp_psoil = tmp_nsoil*org_pnratio
!          tmp_plitter = tmp_nlitter*org_pnratio
          ! record the organic manure deposition
          agri_masschangecnp(i,j,icarbon) = agri_masschangecnp(i,j,icarbon) + &
             tmp_csoil+tmp_clitter
          agri_masschangecnp(i,j,initrogen) = agri_masschangecnp(i,j,initrogen) + &
             tmp_nsoil+tmp_nlitter
          agri_masschangecnp(i,j,iphosphorus) = agri_masschangecnp(i,j,iphosphorus) + &
             tmp_psoil+tmp_plitter
          ! then put the feces part to litter and soil organic pool
          ! for litter the fresh_fmetabolic = 0.3 is chosen to reproduce the
          ! litter decomposition rate of 0.02/day in Li et al., 2012 Table 2
          litter(i,imetabolic,j,iabove,icarbon) = litter(i,imetabolic,j,iabove,icarbon) + &
             tmp_clitter * fresh_fmetabolic
          litter(i,imetabolic,j,iabove,initrogen) = litter(i,imetabolic,j,iabove,initrogen) + &
             tmp_nlitter * fresh_fmetabolic
          litter(i,imetabolic,j,iabove,iphosphorus) = litter(i,imetabolic,j,iabove,iphosphorus) + &
             tmp_plitter * fresh_fmetabolic
          litter(i,istructural,j,iabove,icarbon) = litter(i,istructural,j,iabove,icarbon) + &
             tmp_clitter * (un-fresh_fmetabolic)
          litter(i,istructural,j,iabove,initrogen) = litter(i,istructural,j,iabove,initrogen) + &
             tmp_nlitter * (un-fresh_fmetabolic)
          litter(i,istructural,j,iabove,iphosphorus) = litter(i,istructural,j,iabove,iphosphorus) + &
             tmp_plitter * (un-fresh_fmetabolic)
          ! for soil, the fresh_fisurface = 0.3431 is chosen to reproduce the
          ! soil decomposition rate of 0.006/day in Li et al., 2012 Table 2
          som(i,isurface,j,icarbon) = som(i,isurface,j,icarbon) + &
             tmp_csoil * fresh_fisurface
          som(i,isurface,j,initrogen) = som(i,isurface,j,initrogen) + &
             tmp_nsoil * fresh_fisurface
          som(i,isurface,j,iphosphorus) = som(i,isurface,j,iphosphorus) + &
             tmp_psoil * fresh_fisurface
          som(i,islow,j,icarbon) = som(i,islow,j,icarbon) + &
             tmp_csoil * (un-fresh_fisurface)
          som(i,islow,j,initrogen) = som(i,islow,j,initrogen) + &
             tmp_nsoil * (un-fresh_fisurface)
          som(i,islow,j,iphosphorus) = som(i,islow,j,iphosphorus) + &
             tmp_psoil * (un-fresh_fisurface)

        ENDIF

        ! it is time to synthetic and organic fertilization
        IF (veget_cov_max(i,j) .GT. min_stomate .AND. &
            (tmp_Tfert(i,j)-0.5 .LT. REAL(day_of_year)) .AND. &
            (tmp_Tfert(i,j)+0.5 .GT. REAL(day_of_year))) THEN 
          ! only for the day with fertilization, record it for history write
          real_Nmineral(i,j) = tmp_Nmineral(i,j)
          real_Nmanure(i,j) = tmp_Nmanure(i,j)
          real_Pmineral(i,j) = tmp_Pmineral(i,j)
          ! first apply mineral N and P
          agri_nfert(i,j,iammonium) = agri_nfert(i,j,iammonium) + real_Nmineral(i,j) * min_famm
          agri_nfert(i,j,initrate)  = agri_nfert(i,j,initrate) + real_Nmineral(i,j) * (un - min_famm)
          agri_pfert(i,j)           = agri_pfert(i,j) + real_Pmineral(i,j)

          ! then apply manure C, N and P to litter
          ! considering typical C:N:P ratio in manure
          ! Here, assumptions were made for the C:N:P ratios, 
          ! C:N ratio of slurry 10:1 was used considering its easy application
          ! Soussana and Lemaire, 2014
          ! N:P ratio we chose 10:3 as mean ratio
          ! Phosphorus excreted via urine by ruminants is usually negligible
          ! because a possible surplus of digested P is recirculated to the
          ! rumen and from there excreted in faeces instead of in urine.
          ! http://www.fao.org/WAIRDOCS/LEAD/X6113E/x6113e06.htm
          ! 10:3 of N:P ratio is assumed
          ! according to US states data
          ! https://www.epa.gov/nutrient-policy-data/estimated-animal-agriculture-nitrogen-and-phosphorus-manure

          ! NEED to revise considering the C:N:P ratio in plant and animal
          ! for metabolic ratio, I now used fvgcmetabolic=0.8 in GRM

          ! NOTE: the slurry contain only part of solid manure N (4-10%) 
          ! Here, we assumed a slurry application for N
          ! use a N in urine as 0.9 = fkanurine 
          ! for the solid part 0.9 = fkacmetabolic in GRM
          ! the urine N part is applied directly to the soil mineral N pool as
          ! mineral ammonium fertilizer

          ! NOTE: DLEM used a 7:1 ratio of NH4+:NO3- for manure application
          ! the fractions used here actually generally close to that in Li et
          ! al., 2012, after treatment in farm, the manure contains much more
          ! inorganic N than fresh manure, and during manure storage and
          ! treatment, NH4+ could be the major part of slurry.
          agri_nfert(i,j,iammonium) = agri_nfert(i,j,iammonium) + &
             real_Nmanure(i,j) * org_furine
          ! Urine C and P is negligible
 
          litter(i,imetabolic,j,iabove,icarbon) = litter(i,imetabolic,j,iabove,icarbon) + &
             real_Nmanure(i,j) * (un - org_furine) * org_fmetabolic * org_cnratio
          litter(i,imetabolic,j,iabove,initrogen) = litter(i,imetabolic,j,iabove,initrogen) + &
             real_Nmanure(i,j) * (un - org_furine) * org_fmetabolic
          litter(i,imetabolic,j,iabove,iphosphorus) = litter(i,imetabolic,j,iabove,iphosphorus) + &
             real_Nmanure(i,j) * org_fmetabolic * org_pnratio

          litter(i,istructural,j,iabove,icarbon) = litter(i,istructural,j,iabove,icarbon) + &
             real_Nmanure(i,j) * (un - org_furine) * (un-org_fmetabolic) * org_cnratio
          litter(i,istructural,j,iabove,initrogen) = litter(i,istructural,j,iabove,initrogen) + &
             real_Nmanure(i,j) * (un - org_furine) * (un-org_fmetabolic)
          litter(i,istructural,j,iabove,iphosphorus) = litter(i,istructural,j,iabove,iphosphorus) + &
             real_Nmanure(i,j) * (un-org_fmetabolic) * org_pnratio

          agri_masschangecnp(i,j,icarbon) = agri_masschangecnp(i,j,icarbon) + &
             real_Nmanure(i,j) * (un - org_furine) * org_cnratio
          agri_masschangecnp(i,j,initrogen) = agri_masschangecnp(i,j,initrogen) + &
             real_Nmanure(i,j) * (un - org_furine)
          agri_masschangecnp(i,j,iphosphorus) = agri_masschangecnp(i,j,iphosphorus) + &
             real_Nmanure(i,j) * org_pnratio
        ENDIF
      ENDDO ! i npts
    ENDDO ! j nvm

      ! write history only organic fert to litter and 
      ! + grazing manure deposition to litter and soil is write here
      ! the mineral part and the urine part is write in nitrogen_dynamics,
      ! and phosphorus_dynamics 
      CALL xios_orchidee_send_field('AGRI_ORGFERT_c',agri_masschangecnp(:,:,icarbon))
      CALL histwrite_p (hist_id_stomate, 'AGRI_ORGFERT_c', itime, &
         agri_masschangecnp(:,:,icarbon), npts*nvm, horipft_index)
      CALL xios_orchidee_send_field('AGRI_ORGFERT_n',agri_masschangecnp(:,:,initrogen))
      CALL histwrite_p (hist_id_stomate, 'AGRI_ORGFERT_n', itime, &
         agri_masschangecnp(:,:,initrogen), npts*nvm, horipft_index)
      CALL xios_orchidee_send_field('AGRI_ORGFERT_p',agri_masschangecnp(:,:,iphosphorus))
      CALL histwrite_p (hist_id_stomate, 'AGRI_ORGFERT_p', itime, &
         agri_masschangecnp(:,:,iphosphorus), npts*nvm, horipft_index)

    CALL xios_orchidee_send_field("fNfer",SUM(real_Nmineral(:,:)*veget_cov_max,dim=2)/1e3/one_day*contfrac)
    CALL xios_orchidee_send_field("fMnr",SUM(real_Nmanure(:,:)*veget_cov_max,dim=2)/1e3/one_day*contfrac)
    CALL xios_orchidee_send_field("fMgr",SUM(real_Ngrazing(:,:)*veget_cov_max,dim=2)/1e3/one_day*contfrac)

  END SUBROUTINE apply_agri_fert

END MODULE stomate_lpj
