! =================================================================================================================================
! MODULE       : stomate
!
! CONTACT      : orchidee-help _at_ listes.ipsl.fr
!
! LICENCE      : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF       Groups the subroutines that: (1) initialize all variables in 
!! stomate, (2) read and write forcing files of stomate and the soil component,
!! (3) aggregates and convert variables to handle the different time steps 
!! between sechiba and stomate, (4) call subroutines that govern major stomate
!! processes (litter, soil, and vegetation dynamics) and (5) structures these tasks 
!! in stomate_main
!!
!!\n DESCRIPTION : None
!!
!! RECENT CHANGE(S) : None
!!
!! REFERENCE(S)	: None
!!
!! SVN :
!! $HeadURL: svn://forge.ipsl.jussieu.fr/orchidee/branches/ORCHIDEE-CN-P/ORCHIDEE/src_stomate/stomate.f90 $
!! $Date: 2021-05-06 11:05:55 +0200 (四, 2021-05-06) $
!! $Revision: 7177 $
!! \n
!_ ================================================================================================================================

MODULE stomate

  ! Modules used:
  USE netcdf
  USE defprec
  USE grid
  USE time, ONLY : one_day, one_year, dt_sechiba, dt_stomate, LastTsYear, LastTsMonth
  USE time, ONLY : year_end, month_end, day_end, sec_end
  USE constantes
  USE constantes_soil
  USE pft_parameters
  USE stomate_io
  USE stomate_data
  USE stomate_season
  USE stomate_lpj
  USE stomate_litter
  USE stomate_phosphorus
  USE stomate_vmax
  USE stomate_som_dynamics
  USE stomate_resp
  USE mod_orchidee_para
  USE ioipsl_para 
  USE xios_orchidee

  USE matrix_resolution
  
  IMPLICIT NONE

  ! Private & public routines

  PRIVATE
  PUBLIC stomate_main,stomate_clear,init_forcing, stomate_forcing_read, stomate_initialize, stomate_finalize

  INTERFACE stomate_accu
     MODULE PROCEDURE stomate_accu_r1d, stomate_accu_r2d, stomate_accu_r3d, stomate_accu_r4d
  END INTERFACE

  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:,:,:):: biomass              !! Biomass per ground area @tex $(gC m^{-2})$ @endtex 
!$OMP THREADPRIVATE(biomass)
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:)    :: veget_cov_max        !! Maximal fractional coverage: maximum share of a pixel
                                                                         !! taken by a PFT 
!$OMP THREADPRIVATE(veget_cov_max)
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:)    :: ind                  !! Vegetation density, number of individuals per unit 
                                                                         !! ground area @tex $(m^{-2})$ @endtex 
!$OMP THREADPRIVATE(ind)
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:)    :: age                  !! Age of PFT it normalized by biomass - can increase and
                                                                         !! decrease - (years)
!$OMP THREADPRIVATE(age)
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:)    :: adapted              !! Winter too cold for PFT to survive (0-1, unitless)
!$OMP THREADPRIVATE(adapted)
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:)    :: regenerate           !! Winter sufficiently cold to produce viable seeds 
                                                                         !! (0-1, unitless)
!$OMP THREADPRIVATE(regenerate)
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:)    :: everywhere           !! Is the PFT everywhere in the grid box or very localized 
                                                                         !! (after its intoduction)
!$OMP THREADPRIVATE(everywhere)
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:)    :: fireindex            !! Probability of fire (unitless)
!$OMP THREADPRIVATE(fireindex)
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:)    :: veget_lastlight      !! Vegetation fractions (on ground) after last light 
                                                                         !! competition (unitless) 
!$OMP THREADPRIVATE(veget_lastlight)
  REAL, ALLOCATABLE,SAVE,DIMENSION(:,:)   :: fpc_max              !! "maximal" coverage fraction of a grid box (LAI -> 
                                                                         !! infinity) on ground. [??CHECK??] It's set to zero here, 
                                                                         !! and then is used once in lpj_light.f90 to test if 
                                                                         !! fpc_nat is greater than it. Something seems missing
!$OMP THREADPRIVATE(fpc_max)
  LOGICAL,ALLOCATABLE,SAVE,DIMENSION(:,:)        :: PFTpresent           !! PFT exists (equivalent to veget > 0 for natural PFTs)
!$OMP THREADPRIVATE(PFTpresent)
  LOGICAL,ALLOCATABLE,SAVE,DIMENSION(:,:)        :: senescence           !! The PFT is senescent
!$OMP THREADPRIVATE(senescence)
  LOGICAL,ALLOCATABLE,SAVE,DIMENSION(:,:)        :: begin_leaves         !! Signal to start putting leaves on (true/false)
!$OMP THREADPRIVATE(begin_leaves)
  LOGICAL,ALLOCATABLE,SAVE,DIMENSION(:,:)        :: need_adjacent        !! This PFT needs to be in present in an adjacent gridbox 
                                                                         !! if it is to be introduced in a new gridbox
!$OMP THREADPRIVATE(need_adjacent)
!--
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:)    :: humrel_daily         !! Daily plant available water -root profile weighted 
                                                                         !! (0-1, unitless)
!$OMP THREADPRIVATE(humrel_daily)
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:)    :: humrel_week          !! "Weekly" plant available water -root profile weighted
                                                                         !! (0-1, unitless)
!$OMP THREADPRIVATE(humrel_week)
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:)    :: humrel_month         !! "Monthly" plant available water -root profile weighted
                                                                         !! (0-1, unitless)
!$OMP THREADPRIVATE(humrel_month)
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:)    :: maxhumrel_lastyear   !! Last year's max plant available water -root profile 
                                                                         !! weighted (0-1, unitless)
!$OMP THREADPRIVATE(maxhumrel_lastyear)
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:)    :: maxhumrel_thisyear   !! This year's max plant available water -root profile 
                                                                         !! weighted (0-1, unitless) 
!$OMP THREADPRIVATE(maxhumrel_thisyear)
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:)    :: minhumrel_lastyear   !! Last year's min plant available water -root profile 
                                                                         !! weighted (0-1, unitless)  
!$OMP THREADPRIVATE(minhumrel_lastyear)
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:)    :: minhumrel_thisyear   !! This year's minimum plant available water -root profile
                                                                         !! weighted (0-1, unitless)
!$OMP THREADPRIVATE(minhumrel_thisyear)
!---  
  REAL,ALLOCATABLE,SAVE,DIMENSION(:)      :: t2m_daily            !! Daily air temperature at 2 meter (K)
!$OMP THREADPRIVATE(t2m_daily)

  REAL,ALLOCATABLE,SAVE,DIMENSION(:)      :: Tseason              !! "seasonal" 2 meter temperatures (K)
!$OMP THREADPRIVATE(Tseason)
  REAL,ALLOCATABLE,SAVE,DIMENSION(:)      :: Tseason_length       !! temporary variable to calculate Tseason
!$OMP THREADPRIVATE(Tseason_length)
  REAL,ALLOCATABLE,SAVE,DIMENSION(:)      :: Tseason_tmp          !! temporary variable to calculate Tseason
!$OMP THREADPRIVATE(Tseason_tmp)
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:)    :: Tmin_spring_time     !! Number of days after begin_leaves (leaf onset) 
!$OMP THREADPRIVATE(Tmin_spring_time)
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:)    :: onset_date           !! Date in the year at when the leaves started to grow(begin_leaves), only for diagnostics.
!$OMP THREADPRIVATE(onset_date)
  REAL,ALLOCATABLE,SAVE,DIMENSION(:)      :: t2m_week             !! Mean "weekly" (default 7 days) air temperature at 2 
                                                                         !! meter (K)  
!$OMP THREADPRIVATE(t2m_week)
  REAL,ALLOCATABLE,SAVE,DIMENSION(:)      :: t2m_month            !! Mean "monthly" (default 20 days) air temperature at 2 
                                                                         !! meter (K)
!$OMP THREADPRIVATE(t2m_month)
  REAL,ALLOCATABLE,SAVE,DIMENSION(:)      :: t2m_longterm         !! Mean "Long term" (default 3 years) air temperature at 
                                                                         !! 2 meter (K) 
!$OMP THREADPRIVATE(t2m_longterm)
  REAL,ALLOCATABLE,SAVE,DIMENSION(:)      :: t2m_min_daily        !! Daily minimum air temperature at 2 meter (K)
!$OMP THREADPRIVATE(t2m_min_daily)
  REAL,ALLOCATABLE,SAVE,DIMENSION(:)      :: tsurf_daily          !! Daily surface temperatures (K)
!$OMP THREADPRIVATE(tsurf_daily)
!---
  REAL,ALLOCATABLE,SAVE,DIMENSION(:)      :: precip_daily         !! Daily precipitations sum @tex $(mm day^{-1})$ @endtex
!$OMP THREADPRIVATE(precip_daily)
  REAL,ALLOCATABLE,SAVE,DIMENSION(:)      :: precip_lastyear      !! Last year's annual precipitation sum 
                                                                         !! @tex $??(mm year^{-1})$ @endtex
!$OMP THREADPRIVATE(precip_lastyear)
  REAL,ALLOCATABLE,SAVE,DIMENSION(:)      :: precip_thisyear      !! This year's annual precipitation sum 
                                                                         !! @tex $??(mm year^{-1})$ @endtex 
!$OMP THREADPRIVATE(precip_thisyear)
!---
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:)    :: soilhum_daily        !! Daily soil humidity (0-1, unitless)
!$OMP THREADPRIVATE(soilhum_daily)
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:)    :: soilhum_month        !! Soil humidity - integrated over a month (0-1, unitless) 
!$OMP THREADPRIVATE(soilhum_month)
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:)    :: tsoil_daily          !! Daily soil temperatures (K)
!$OMP THREADPRIVATE(tsoil_daily)
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:)    :: tsoil_month          !! Soil temperatures at each soil layer integrated over a
                                                                         !! month (K) 
!$OMP THREADPRIVATE(tsoil_month)
!--- 
  REAL,ALLOCATABLE,SAVE,DIMENSION(:)      :: litterhum_daily      !! Daily litter humidity (0-1, unitless)
!$OMP THREADPRIVATE(litterhum_daily)
!---
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:,:)  :: control_moist        !! Moisture control of heterotrophic respiration 
                                                                         !! (0-1, unitless)
!$OMP THREADPRIVATE(control_moist)
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:,:)    :: drainage           !! Water lost from the soil column by leaching (kg/m^2/timestep) dt_sechiba?
!$OMP THREADPRIVATE(drainage)
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:)      :: drainage_daily     !! Water lost from the soil column by leaching daily (kg/m^2/day)
!$OMP THREADPRIVATE(drainage_daily)
  REAL,ALLOCATABLE,SAVE,DIMENSION(:)        :: drainage_longterm  !! Mean "Long term" (default 3 years) annual sum of drainage (kg/m^2/yr) 
!$OMP THREADPRIVATE(drainage_longterm)
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:,:)    :: runoff             !! Water lost from the soil column by leaching (kg/m^2/timestep) 
!$OMP THREADPRIVATE(runoff)
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:)      :: runoff_daily       !! Water lost from the soil column by leaching daily (kg/m^2/day)
!$OMP THREADPRIVATE(runoff_daily)
  REAL,ALLOCATABLE,SAVE,DIMENSION(:)      :: runoff_longterm      !! Mean "Long term" (default 3 years) annual sum of runoff (kg/m^2/yr)
!$OMP THREADPRIVATE(runoff_longterm)


  REAL,ALLOCATABLE,SAVE, DIMENSION(:,:)   :: n_mineralisation_d   !! net nitrogen mineralisation of decomposing SOM
                                                                         !!   (gN/m**2/day), supposed to be NH4 
  REAL,ALLOCATABLE,SAVE, DIMENSION(:,:,:) :: n_uptake_daily       !! Uptake of soil N by plants 
                                                                         !! (gN/m**2/day) 

  REAL,ALLOCATABLE,SAVE, DIMENSION(:,:) :: p_mineralisation_d     !! net phosphorus mineralisation of decomposing SOM
                                                                         !!   (g(P)/m**2/day), supposed to labile (organic & inorganic)
  REAL,ALLOCATABLE,SAVE, DIMENSION(:,:)  :: p_uptake_daily        !! Uptake of soil P by plants 
                                                                         !! (gP/m**2/day) 
  REAL,ALLOCATABLE,SAVE, DIMENSION(:)    :: p_release_daily       !! release of P from chemical weathering of rocks
                                                                         !!   (gP/m**2/day), supposed to be PO4 
  REAL,ALLOCATABLE,SAVE, DIMENSION(:,:)    :: p_ssorption_daily   !! strongly sorption of sorbed labile P
                                                                         !!   (gP/m**2/day), supposed to be PO4 
  REAL,ALLOCATABLE,SAVE, DIMENSION(:,:)    :: p_leaching_daily    !! leaching of dissolved labile P
                                                                         !!   (gP/m**2/day), supposed to be PO4 
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:,:) :: soil_p_min            !! mineral phosphorus in the soil (gP/m**2) 
                                                                         !! (first index=npts, second index=nvm, third index=npspec)  
                                                                         !! DSG: it also includes organic P in labile form, 
                                                                         !! which amount is << organic P in passive, slow & fast pools)


  REAL,ALLOCATABLE,SAVE, DIMENSION(:,:,:,:) :: EW_grain           !! Enhanced Weathering: grain mass(1) [g m-2] and grain size(2) [mu m]
  REAL,ALLOCATABLE,SAVE, DIMENSION(:,:,:)   :: EW_flux            !! Enhanced weathering: CO2 cons. flux (iCO2cons), P release (iPrelease)

!$OMP THREADPRIVATE(soil_p_min)

  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:) :: f_Pdissolved            !! fraction of mineral phosphorus in dissolved in soil solution in contact with roots [ ]
!$OMP THREADPRIVATE(f_Pdissolved)



  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:)  :: deposition           !!  nutrient inputs by atmospheric deposition [g/m2/day] 
!$OMP THREADPRIVATE(deposition)

  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:)  :: lith_frac            !!  lithological fractions [0-1] 
!$OMP THREADPRIVATE(lith_frac)

  REAL,ALLOCATABLE, SAVE,DIMENSION(:)   :: soil_shield             !! factor to correct P weathering for soil shielding effects [0-1]
!$OMP THREADPRIVATE(soil_shield)
INTEGER,ALLOCATABLE, SAVE,DIMENSION(:)   :: soil_orders          !! dominant USDA soil order of gridbox
!$OMP THREADPRIVATE(soil_soilorders)

  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:,:)  :: control_temp         !! Temperature control of heterotrophic respiration at the
                                                                         !! different soil levels (0-1, unitless)
!$OMP THREADPRIVATE(control_temp)
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:)    :: control_moist_daily  !! Moisture control of heterotrophic respiration daily 
                                                                         !! (0-1, unitless)
!$OMP THREADPRIVATE(control_moist_daily)
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:)    :: control_temp_daily   !! Temperature control of heterotrophic respiration, above
                                                                         !! and below daily (0-1, unitless)
!$OMP THREADPRIVATE(control_temp_daily)
!---
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:)    :: gdd_init_date        !! inital date for gdd count 
!$OMP THREADPRIVATE(gdd_init_date)

  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:)    :: gdd_from_growthinit  !! gdd from beginning of season (C)
!$OMP THREADPRIVATE(gdd_from_growthinit)
  REAL,ALLOCATABLE,SAVE,DIMENSION(:)      :: gdd0_lastyear        !! Last year's annual Growing Degree Days,
                                                                         !! threshold 0 deg C (K) 
!$OMP THREADPRIVATE(gdd0_lastyear)
  REAL,ALLOCATABLE,SAVE,DIMENSION(:)      :: gdd0_thisyear        !! This year's annual Growing Degree Days,
                                                                         !! threshold 0 deg C (K)
!$OMP THREADPRIVATE(gdd0_thisyear)
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:)    :: gdd_m5_dormance      !! Growing degree days for onset of growing season, 
                                                                         !! threshold -5 deg C (K)
!$OMP THREADPRIVATE(gdd_m5_dormance)
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:)    :: gdd_midwinter        !! Growing degree days for onset of growing season, 
                                                                         !! since midwinter (K)
!$OMP THREADPRIVATE(gdd_midwinter)
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:)    :: ncd_dormance         !! Number of chilling days since leaves were lost (days) 
!$OMP THREADPRIVATE(ncd_dormance)
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:)    :: ngd_minus5           !! Number of growing days, threshold -5 deg C (days)
!$OMP THREADPRIVATE(ngd_minus5)
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:)    :: hum_min_dormance     !! Minimum moisture during dormance (0-1, unitless) 
!$OMP THREADPRIVATE(hum_min_dormance)
!---
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:)    :: gpp_daily            !! Daily gross primary productivity per ground area 
                                                                         !! @tex $(gC m^{-2} dt_stomate^{-1})$ @endtex
!$OMP THREADPRIVATE(gpp_daily) 
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:)    :: gpp_week             !! Mean "weekly" (default 7 days) GPP  
                                                                         !! @tex $(gC m^{-2} dt_stomate^{-1})$ @endtex
!$OMP THREADPRIVATE(gpp_week)
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:)    :: npp_week             !! Mean "weekly" (default 7 days) NPP  
                                                                         !! @tex $(gC m^{-2} dt_slow^{-1})$ @endtex
!$OMP THREADPRIVATE(npp_week)
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:)    :: maxgppweek_lastyear  !! Last year's maximum "weekly" GPP  
                                                                         !! @tex $(gC m^{-2} dt_stomate^{-1})$ @endtex 
!$OMP THREADPRIVATE(maxgppweek_lastyear)
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:)    :: maxgppweek_thisyear  !! This year's maximum "weekly" GPP  
                                                                         !! @tex $(gC m^{-2} dt_stomate^{-1})$ @endtex  
!$OMP THREADPRIVATE(maxgppweek_thisyear)
!---
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:)    :: npp_daily            !! Daily net primary productivity per ground area 
                                                                         !! @tex $(gC m^{-2} dt_stomate^{-1})$ @endtex 
!$OMP THREADPRIVATE(npp_daily)
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:)    :: npp_longterm         !! "Long term" (default 3 years) net primary productivity 
                                                                         !! per ground area  
                                                                         !! @tex $(gC m^{-2} year^{-1})$ @endtex   
!$OMP THREADPRIVATE(npp_longterm)
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:)    :: npp_equil            !! Equilibrium NPP written to forcesoil 
                                                                         !! @tex $(gC m^{-2} year^{-1})$ @endtex
!$OMP THREADPRIVATE(npp_equil)
  REAL,ALLOCATABLE,SAVE,DIMENSION(:)      :: npp_tot              !! Total NPP written to forcesoil 
                                                                         !! @tex $(gC m^{-2} year^{-1})$ @endtex 
!$OMP THREADPRIVATE(npp_tot)
!---
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:,:)  :: resp_maint_part_radia!! Maintenance respiration of different plant parts per 
                                                                         !! total ground area at Sechiba time step  
                                                                         !! @tex $(gC m^{-2} dt_sechiba^{-1})$ @endtex
!$OMP THREADPRIVATE(resp_maint_part_radia)
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:,:)  :: resp_maint_part      !! Maintenance respiration of different plant parts per
                                                                         !! total ground area at Stomate time step 
                                                                         !! @tex $(gC m^{-2} dt_stomate^{-1})$ @endtex
!$OMP THREADPRIVATE(resp_maint_part)
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:)    :: resp_maint_radia     !! Maintenance respiration per ground area at Sechiba time
                                                                         !! step   
                                                                         !! @tex $(gC m^{-2} dt_sechiba^{-1})$ @endtex
!$OMP THREADPRIVATE(resp_maint_radia)
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:)    :: resp_maint_d         !! Maintenance respiration per ground area at Stomate time 
                                                                         !! step  
                                                                         !! @tex $(gC m^{-2} dt_stomate^{-1})$ @endtex
!$OMP THREADPRIVATE(resp_maint_d)
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:)    :: resp_growth_d        !! Growth respiration per ground area 
                                                                         !! @tex $(gC m^{-2} dt_stomate^{-1})$ @endtex
!$OMP THREADPRIVATE(resp_growth_d)
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:)    :: resp_hetero_d        !! Heterotrophic respiration per ground area 
                                                                         !! @tex $(gC m^{-2} dt_stomate^{-1})$ @endtex
!$OMP THREADPRIVATE(resp_hetero_d)
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:)    :: resp_hetero_radia    !! Heterothrophic respiration per ground area at Sechiba
                                                                         !! time step 
                                                                         !! @tex $(gC m^{-2} dt_sechiba^{-1})$ @endtex 
!$OMP THREADPRIVATE(resp_hetero_radia)
!---
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:)     :: turnover_time       !! Turnover time of grasses 
                                                                         !! @tex $(dt_stomate^{-1})$ @endtex 
!$OMP THREADPRIVATE(turnover_time)
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:,:,:) :: turnover_daily      !! Senescence-driven turnover (better: mortality) of 
                                                                         !! leaves and roots  
                                                                         !! @tex $(gC m^{-2} dt_stomate^{-1})$ @endtex
!$OMP THREADPRIVATE(turnover_daily)
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:,:,:) :: turnover_littercalc !! Senescence-driven turnover (better: mortality) of 
                                                                         !! leaves and roots at Sechiba time step 
                                                                         !! @tex $(gC m^{-2} dt_sechiba^{-1})$ @endtex 
!$OMP THREADPRIVATE(turnover_littercalc)
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:,:,:) :: turnover_longterm   !! "Long term" (default 3 years) senescence-driven 
                                                                         !! turnover (better: mortality) of leaves and roots 
                                                                         !! @tex $(gC m^{-2} year^{-1})$ @endtex
!$OMP THREADPRIVATE(turnover_longterm)
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:,:,:) :: bm_to_litter        !! Background (not senescence-driven) mortality of biomass
                                                                         !! @tex $(gC m^{-2} dt_stomate^{-1})$ @endtex
!$OMP THREADPRIVATE(bm_to_litter)
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:,:,:) :: bm_to_littercalc    !! conversion of biomass to litter per ground area at 
                                                                         !! Sechiba time step 
                                                                         !! @tex $(gC m^{-2} dt_sechiba^{-1})$ @endtex 
!$OMP THREADPRIVATE(bm_to_littercalc)
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:,:)  :: dead_leaves          !! Metabolic and structural pools of dead leaves on ground
                                                                         !! per PFT @tex $(gC m^{-2})$ @endtex 
!$OMP THREADPRIVATE(dead_leaves)
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:,:,:,:):: litter             !! Above and below ground metabolic and structural litter 
                                                                         !! per ground area 
                                                                         !! @tex $(gC m^{-2})$ @endtex 
!$OMP THREADPRIVATE(litter)
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:)    :: firelitter           !! Total litter above the ground that could potentially 
                                                                         !! burn @tex $(gC m^{-2})$ @endtex 
!$OMP THREADPRIVATE(firelitter)
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:,:,:):: carbon_input     !! Quantity of carbon going into carbon pools from litter
                                                                         !! decomposition per ground area  at Sechiba time step 
                                                                         !! @tex $(gC m^{-2} dt_sechiba^{-1})$ @endtex 
!$OMP THREADPRIVATE(carbon_input)
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:,:,:):: nitrogen_input       !! Quantity of nitrogen going into nitrogen pools from litter 
                                                                         !! decomposition per ground area  at Sechiba time step  
                                                                         !! @tex $(gC m^{-2} dtradia^{-1})$ @endtex  
!$OMP THREADPRIVATE(nitrogen_input) 
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:,:,:)  :: som_input_daily !! Daily quantity of carbon going into carbon pools from
                                                                           !! litter decomposition per ground area 
                                                                           !! @tex $(gC m^{-2} day^{-1})$ @endtex 
!$OMP THREADPRIVATE(som_input_daily)
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:,:,:)  :: som                !! Soil organic matter  pools per ground area: active, slow, or 
                                                                         !! passive, @tex $(gC or N m^{-2})$ @endtex 
!$OMP THREADPRIVATE(som)
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:,:)  :: lignin_struc         !! Ratio Lignine/Carbon in structural litter for above and
                                                                         !! below ground compartments (unitless)
!$OMP THREADPRIVATE(lignin_struc)
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:,:)  :: lignin_wood          !! Ratio Lignine/Carbon in woody litter for above and
                                                                         !! below ground compartments (unitless)	
!$OMP THREADPRIVATE(lignin_wood)
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:)    :: lm_lastyearmax       !! Last year's maximum leaf mass per ground area for each
                                                                         !! PFT @tex $(gC m^{-2})$ @endtex  
!$OMP THREADPRIVATE(lm_lastyearmax)
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:)    :: lm_thisyearmax       !! This year's maximum leaf mass per ground area for each
                                                                         !! PFT @tex $(gC m^{-2})$ @endtex  
!$OMP THREADPRIVATE(lm_thisyearmax)
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:)    :: maxfpc_lastyear      !! Last year's maximum fpc for each natural PFT, on ground
                                                                         !! [??CHECK] fpc but this ones look ok (computed in 
                                                                         !! season, used in light)?? 
!$OMP THREADPRIVATE(maxfpc_lastyear)
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:)    :: maxfpc_thisyear      !! This year's maximum fpc for each PFT, on ground (see 
                                                                         !! stomate_season), [??CHECK] fpc but this ones look ok 
                                                                         !! (computed in season, used in light)??
!$OMP THREADPRIVATE(maxfpc_thisyear)
!---
  REAL, ALLOCATABLE,SAVE,DIMENSION(:,:,:,:)  :: leaf_age             !! Age of different leaf classes (days)
!$OMP THREADPRIVATE(leaf_age)
  REAL, ALLOCATABLE,SAVE,DIMENSION(:,:,:,:)  :: leaf_frac            !! PFT fraction of leaf mass in leaf age class (0-1, 
                                                                         !! unitless) 
!$OMP THREADPRIVATE(leaf_frac)
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:)    :: when_growthinit      !! Days since beginning of growing season (days)
!$OMP THREADPRIVATE(when_growthinit)
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:)    :: herbivores           !! Time constant of probability of a leaf to be eaten by a
                                                                         !! herbivore (days)
!$OMP THREADPRIVATE(herbivores)
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:)    :: RIP_time             !! How much time ago was the PFT eliminated for the last 
                                                                         !! time (year)
!$OMP THREADPRIVATE(RIP_time)
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:)    :: time_hum_min         !! Time elapsed since strongest moisture limitation (days) 
!$OMP THREADPRIVATE(time_hum_min)
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:)    :: drain_daily          !! daily fraction of water lost from the soil column by leaching (-)  ! DSG: is it actual the fraction??!
!$OMP THREADPRIVATE(drain_daily)

!---
  REAL,ALLOCATABLE,DIMENSION(:,:)         :: clay_fm              !! Soil clay content (0-1, unitless), parallel computing
  REAL,ALLOCATABLE,DIMENSION(:,:)         :: clay_fm_g            !! Soil clay content (0-1, unitless), parallel computing
  REAL,ALLOCATABLE,DIMENSION(:,:)         :: silt_fm              !! Soil silt content (0-1, unitless), parallel computing 
  REAL,ALLOCATABLE,DIMENSION(:,:)         :: silt_fm_g            !! Soil silt content (0-1, unitless), parallel computing 
  REAL,ALLOCATABLE,DIMENSION(:,:)         :: bulk_fm              !! Bulk density (kg/m**3), parallel computing    
  REAL,ALLOCATABLE,DIMENSION(:,:)         :: bulk_fm_g            !! Bulk density (kg/m**3), parallel computing 

  REAL,ALLOCATABLE,DIMENSION(:,:)         :: precip_fm            !! Daily precipitations sum @tex $(mm day^{-1})$ @endtex,
                                                                         !! parallel computing
  REAL,ALLOCATABLE,DIMENSION(:,:)         :: precip_fm_g          !! Daily precipitations sum @tex $(mm day^{-1})$ @endtex,
                                                                         !! parallel computing
  REAL,ALLOCATABLE,DIMENSION(:,:)         :: litterhum_daily_fm   !! Daily relative humidity of litter (0-1, unitless), 
                                                                         !! parallel computing
  REAL,ALLOCATABLE,DIMENSION(:,:)         :: litterhum_daily_fm_g !! Daily relative humidity of litter (0-1, unitless), 
                                                                         !! parallel computing
  REAL,ALLOCATABLE,DIMENSION(:,:)         :: t2m_daily_fm         !! Daily air temperature at 2 meter (K), parallel 
                                                                         !! computing
  REAL,ALLOCATABLE,DIMENSION(:,:)         :: t2m_daily_fm_g       !! Daily air temperature at 2 meter (K), parallel 
                                                                         !! computing
  REAL,ALLOCATABLE,DIMENSION(:,:)         :: t2m_min_daily_fm     !! Daily minimum air temperature at 2 meter (K), 
                                                                         !! parallel computing
  REAL,ALLOCATABLE,DIMENSION(:,:)         :: t2m_min_daily_fm_g   !! Daily minimum air temperature at 2 meter (K), 
                                                                         !! parallel computing
  REAL,ALLOCATABLE,DIMENSION(:,:)         :: tsurf_daily_fm       !! Daily surface temperatures (K), parallel 
                                                                         !! computing
  REAL,ALLOCATABLE,DIMENSION(:,:)         :: tsurf_daily_fm_g     !! Daily surface temperatures (K), parallel 
                                                                         !! computing 
  REAL,ALLOCATABLE,DIMENSION(:,:,:)       :: tsoil_daily_fm       !! Daily soil temperatures (K), parallel computing 
  REAL,ALLOCATABLE,DIMENSION(:,:,:)       :: tsoil_daily_fm_g     !! Daily soil temperatures (K), parallel computing
  REAL,ALLOCATABLE,DIMENSION(:,:,:)       :: soilhum_daily_fm     !! Daily soil humidity (0-1, unitless), parallel computing
  REAL,ALLOCATABLE,DIMENSION(:,:,:)       :: soilhum_daily_fm_g   !! Daily soil humidity (0-1, unitless), parallel computing
  REAL,ALLOCATABLE,DIMENSION(:,:)       :: max_eau_var_daily_fm !! Daily Maximum water content of the soil @tex $(kg m^{-2})$ @endtex
                                                                         !! parallel computing
  REAL,ALLOCATABLE,DIMENSION(:,:)       :: max_eau_var_daily_fm_g  !! Daily Maximum water content of the soil @tex $(kg m^{-2})$ @endtex
                                                                         !! parallel computing
  REAL,ALLOCATABLE,DIMENSION(:,:,:)       :: humrel_daily_fm      !! Daily relative humidity of atmosphere (0-1, unitless), 
                                                                         !! parallel computing
  REAL,ALLOCATABLE,DIMENSION(:,:,:)       :: humrel_daily_fm_g    !! Daily relative humidity of atmosphere (0-1, unitless), 
                                                                         !! parallel computing
  REAL,ALLOCATABLE,DIMENSION(:,:,:)       :: gpp_daily_fm         !! Daily gross primary productivity per ground area 
                                                                         !! @tex $(gC m^{-2} dt_stomate^{-1})$ @endtex, 
                                                                         !! parallel computing
  REAL,ALLOCATABLE,DIMENSION(:,:,:)       :: gpp_daily_fm_g       !! Daily gross primary productivity per ground area 
                                                                         !! @tex $(gC m^{-2} day^{-1})$ @endtex, 
                                                                         !! parallel computing
  REAL,ALLOCATABLE,DIMENSION(:,:,:)       :: veget_fm             !! Vegetation coverage taking into account non-biological
                                                                         !! coverage (unitless), parallel computing
  REAL,ALLOCATABLE,DIMENSION(:,:,:)       :: veget_fm_g           !! Vegetation coverage taking into account non-biological
                                                                         !! coverage (unitless), parallel computing
  REAL,ALLOCATABLE,DIMENSION(:,:,:)       :: veget_max_fm         !! Maximum vegetation coverage taking into account 
                                                                         !! non-biological coverage (unitless), parallel computing
  REAL,ALLOCATABLE,DIMENSION(:,:,:)       :: veget_max_fm_g       !! Maximum vegetation coverage taking into account none 
                                                                         !! biological coverage (unitless), parallel computing
  REAL,ALLOCATABLE,DIMENSION(:,:,:)       :: lai_fm               !! Leaf area index @tex $@tex $(m^2 m^{-2})$ @endtex$ @endtex, 
                                                                         !! parallel computing
  REAL,ALLOCATABLE,DIMENSION(:,:,:)       :: lai_fm_g             !! Leaf area index @tex $@tex $(m^2 m^{-2})$ @endtex$ @endtex, 
                                                                         !! parallel computing
  REAL,ALLOCATABLE,DIMENSION(:,:,:)        :: drainage_fm         !! daily fraction of water lost from the soil column by leaching (-)                                                                           ! DSG: is it actual the fraction? or kg/m2/time?
  REAL,ALLOCATABLE,DIMENSION(:,:,:)        :: drainage_fm_g       !! drain 
                                                                         !! parallel computing 
  REAL,ALLOCATABLE,DIMENSION(:,:,:)        :: runoff_fm           !! daily sum of water lost from the soil column by leaching (kg/m2/day) 
                                                                        
  REAL,ALLOCATABLE,DIMENSION(:,:,:)        :: runoff_fm_g         !! drain
                                                                         !! parallel computing

 REAL, ALLOCATABLE,SAVE,DIMENSION(:,:)   :: cn_leaf_avg_season    !! Seasonal average CN ratio of leaves 
!$OMP THREADPRIVATE(cn_leaf_avg_season)
 REAL, ALLOCATABLE,SAVE,DIMENSION(:,:)   :: np_leaf_avg_season    !! Seasonal average NP ratio of leaves 
!$OMP THREADPRIVATE(np_leaf_avg_season)
 REAL, ALLOCATABLE,SAVE,DIMENSION(:,:)   :: nstress_season        !! N-related seasonal stress (used for allocation) 
!$OMP THREADPRIVATE(nstress_season)
 REAL, ALLOCATABLE,SAVE,DIMENSION(:,:)   :: pstress_season        !! P-related seasonal stress (used for allocation) 
!$OMP THREADPRIVATE(pstress_season)
 REAL,ALLOCATABLE,SAVE,DIMENSION(:,:)    :: moiavail_growingseason!! mean growing season moisture availability (used for allocation response)
!$OMP THREADPRIVATE(moiavail_growingseason)
 REAL,ALLOCATABLE,SAVE,DIMENSION(:,:)    :: lai_target            !! lai target 
 REAL,ALLOCATABLE,SAVE,DIMENSION(:,:)    :: lai_target_longterm   !! lai target seasonal
!$OMP THREADPRIVATE(lai_target_longterm)

 REAL,ALLOCATABLE,SAVE,DIMENSION(:,:,:)  :: soil_n_min            !! mineral nitrogen in the soil (gN/m**2)   
                                                                         !! (first index=npts, second index=nvm, third index=nnspec)   
!$OMP THREADPRIVATE(soil_n_min)
 REAL, ALLOCATABLE, SAVE, DIMENSION(:,:)           :: p_O2             !! partial pressure of oxigen in the soil (hPa)
                                                                              !! (first index=npts, second index=nvm)
                      
 REAL, ALLOCATABLE, SAVE, DIMENSION(:,:)           :: bact             !! denitrifier biomass (gC/m**2)
                                                                              !! (first index=npts, second index=nvm)
!---
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:)    :: co2_fire             !! Carbon emitted to the atmosphere by burning living 
                                                                         !! and dead biomass 
                                                                         !! @tex $(gC m^{-2} dt_stomate^{-1})$ @endtex 
!$OMP THREADPRIVATE(co2_fire)
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:)    :: co2_to_bm_dgvm       !! Psuedo-photosynthesis,C used to provide seedlings with
                                                                         !! an initial biomass, arbitrarily removed from the 
                                                                         !! atmosphere  
                                                                         !! @tex $(gC m^{-2} dt_stomate^{-1})$ @endtex 
!$OMP THREADPRIVATE(co2_to_bm_dgvm)
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:)    :: n_to_bm              !! N taken from ?? to provide seedlings with
                                                                         !! an initial N biomass
                                                                         !! @tex $(gN m^{-2} dt_stomate^{-1})$ @endtex 
!$OMP THREADPRIVATE(n_to_bm)
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:)    :: p_to_bm              !! P taken from ?? to provide seedlings with
                                                                         !! an initial P biomass
                                                                         !! @tex $(gP m^{-2} dt_stomate^{-1})$ @endtex 
!$OMP THREADPRIVATE(p_to_bm)
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:)    :: co2_flux_daily       !! Daily net CO2 flux between atmosphere and biosphere 
!$OMP THREADPRIVATE(co2_flux_daily)
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:)    :: nep_daily            !! Daily net CO2 flux (positive from atmosphere to land)
                                                                         !! @tex $(gC m^{-2} dt_stomate^{-1})$ @endtex
!$OMP THREADPRIVATE(nep_daily)
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:)    :: nep_monthly          !! Monthly net CO2 flux (positive from atmosphere to land) 
                                                                         !! @tex $(gC m^{-2} dt_stomate^{-1})$ @endtex
!$OMP THREADPRIVATE(nep_monthly)
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:)    :: prod10               !! Wood products remaining in the 10 year-turnover pool 
                                                                         !! after the annual release for each compartment 
                                                                         !! @tex $(gC m^{-2})$ @endtex    
                                                                         !! (0:10 input from year of land cover change),
                                                                         !! dimension(#pixels,0:10 years
!$OMP THREADPRIVATE(prod10)
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:)    :: prod100              !! Wood products remaining in the 100 year-turnover pool
                                                                         !! after the annual release for each compartment
                                                                         !! @tex $(gC m^{-2})$ @endtex  
                                                                         !! (0:100 input from year of land cover change), 
                                                                         !! dimension(#pixels,0:100 years)
!$OMP THREADPRIVATE(prod100)
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:)    :: flux10               !! Wood decomposition from the 10 year-turnover pool 
                                                                         !! compartments 
                                                                         !! @tex $(gC m^{-2} year^{-1})$ @endtex 
                                                                         !! dimension(#pixels,0:10)  
!$OMP THREADPRIVATE(flux10)
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:)    :: flux100              !! Wood decomposition from the 100 year-turnover pool 
                                                                         !! compartments 
                                                                         !! @tex $(gC m^{-2} year^{-1})$ @endtex
                                                                         !! dimension(#pixels,0:100)
!$OMP THREADPRIVATE(flux100)
  REAL,ALLOCATABLE,SAVE,DIMENSION(:)      :: convflux             !! Release during first year following land cover change 
                                                                         !! (paper, burned, etc...) 
                                                                         !! @tex $(gC m^{-2} year^{-1})$ @endtex  
!$OMP THREADPRIVATE(convflux)
  REAL,ALLOCATABLE,SAVE,DIMENSION(:)      :: cflux_prod10         !! Total annual release from the 10 year-turnover pool
                                                                         !! sum of flux10  
                                                                         !! @tex $(gC m^{-2} year^{-1})$ @endtex
!$OMP THREADPRIVATE(cflux_prod10)
  REAL,ALLOCATABLE,SAVE,DIMENSION(:)      :: cflux_prod100        !! Total annual release from the 100 year-turnover pool 
                                                                         !! sum of flux100 
                                                                         !! @tex $(gC m^{-2} year^{-1})$ @endtex
!$OMP THREADPRIVATE(cflux_prod100)
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:)    :: prod10_harvest       !! Wood products remaining in the 10 year-turnover pool 
                                                                         !! after the annual release for each compartment 
                                                                         !! @tex $(gC m^{-2})$ @endtex    
                                                                         !! (0:10 input from year of wood harvest),
                                                                         !! dimension(#pixels,0:10 years
!$OMP THREADPRIVATE(prod10_harvest)
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:)    :: prod100_harvest      !! Wood products remaining in the 100 year-turnover pool
                                                                         !! after the annual release for each compartment
                                                                         !! @tex $(gC m^{-2})$ @endtex  
                                                                         !! (0:100 input from year of wood harvest), 
                                                                         !! dimension(#pixels,0:100 years)
!$OMP THREADPRIVATE(prod100_harvest)
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:)    :: flux10_harvest       !! Wood decomposition from the 10 year-turnover pool 
                                                                         !! compartments 
                                                                         !! @tex $(gC m^{-2} year^{-1})$ @endtex 
                                                                         !! dimension(#pixels,0:10)  
!$OMP THREADPRIVATE(flux10_harvest)
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:)    :: flux100_harvest      !! Wood decomposition from the 100 year-turnover pool 
                                                                         !! compartments 
                                                                         !! @tex $(gC m^{-2} year^{-1})$ @endtex
                                                                         !! dimension(#pixels,0:100)
!$OMP THREADPRIVATE(flux100_harvest)
  REAL,ALLOCATABLE,SAVE,DIMENSION(:)      :: convflux_harvest     !! Release during first year following wood harvest 
                                                                         !! (paper, burned, etc...) 
                                                                         !! @tex $(gC m^{-2} year^{-1})$ @endtex  
!$OMP THREADPRIVATE(convflux_harvest)
  REAL,ALLOCATABLE,SAVE,DIMENSION(:)      :: cflux_prod10_harvest !! Total annual release from the 10 year-turnover pool
                                                                         !! sum of flux10  
                                                                         !! @tex $(gC m^{-2} year^{-1})$ @endtex
!$OMP THREADPRIVATE(cflux_prod10_harvest)
  REAL,ALLOCATABLE,SAVE,DIMENSION(:)      :: cflux_prod100_harvest!! Total annual release from the 100 year-turnover pool 
                                                                         !! sum of flux100 
                                                                         !! @tex $(gC m^{-2} year^{-1})$ @endtex
!$OMP THREADPRIVATE(cflux_prod100_harvest)
  REAL, ALLOCATABLE, SAVE, DIMENSION (:,:):: woodharvestpft       !! New year wood harvest per  PFT
!$OMP THREADPRIVATE(woodharvestpft)
  REAL,ALLOCATABLE,SAVE,DIMENSION(:)      :: harvest_above        !! Harvest of above ground biomass for agriculture -not 
                                                                         !! just from land use change 
                                                                         !! @tex $(??gC m^{-2} dt_stomate^{-1})$ @endtex
!$OMP THREADPRIVATE(harvest_above)
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:)      :: harvest_bioenergy    !! Harvest of above ground biomass for bioenergy
                                                                         !! 
                                                                         !! @tex $(??gC m^{-2} dt_stomate^{-1})$ @endtex
!$OMP THREADPRIVATE(harvest_bioenergy)
  REAL,ALLOCATABLE,SAVE,DIMENSION(:)      :: carb_mass_total      !! Total on-site and off-site C pool 
                                                                         !! @tex $(??gC m^{-2})$ @endtex                        
!$OMP THREADPRIVATE(carb_mass_total)
!---
  REAL, SAVE                              :: tau_longterm
!$OMP THREADPRIVATE(tau_longterm)
  REAL,SAVE                               :: dt_days=zero         !! Time step of STOMATE (days) 
!$OMP THREADPRIVATE(dt_days)
  INTEGER,SAVE                            :: days_since_beg=0     !! Number of full days done since the start of the simulation
!$OMP THREADPRIVATE(days_since_beg)
  INTEGER,ALLOCATABLE,SAVE,DIMENSION(:)   :: nforce               !! Number of states calculated for the soil forcing 
                                                                         !! variables (unitless), dimension(::nparan*::nbyear) both 
                                                                         !! given in the run definition file    
!$OMP THREADPRIVATE(nforce)
  INTEGER,ALLOCATABLE,SAVE,DIMENSION(:)   :: isf                  !! Index for number of time steps that can be stored in 
                                                                         !! memory (unitless), dimension (#nsfm)
!$OMP THREADPRIVATE(isf)
  INTEGER,ALLOCATABLE,SAVE,DIMENSION(:)   :: nf_cumul             !! Number of years over which the average is calculated in
                                                                         !! forcesoil when cumul flag is set, dimension (#nsft)
                                                                         !! [??CHECK] definition the dimension is number of 
                                                                         !! timesteps in a year?
!$OMP THREADPRIVATE(nf_cumul)
  INTEGER, SAVE                           :: spinup_period        !! Period of years used to calculate the resolution of the system for spinup analytic. 
                                                                         !! This period correspond in most cases to the period of years of forcing data used

  INTEGER,PARAMETER                              :: r_typ = nf90_real4   !! Specify data format (server dependent)
  LOGICAL,ALLOCATABLE,SAVE,DIMENSION(:)          :: nf_written           !! Flag indicating whether the forcing data have been 
                                                                         !! written
!$OMP THREADPRIVATE(nf_written)
!---
  LOGICAL, SAVE                                  :: do_slow=.FALSE.      !! Flag that determines whether stomate_accu calculates
                                                                         !! the sum(do_slow=.FALSE.) or the mean 
                                                                         !! (do_slow=.TRUE.)
!$OMP THREADPRIVATE(do_slow)
  LOGICAL, SAVE                                  :: l_first_stomate = .TRUE.!! Is this the first call of stomate?
!$OMP THREADPRIVATE(l_first_stomate)
  LOGICAL, SAVE                                  :: cumul_forcing=.FALSE.!! flag for cumul of forcing if teststomate
!$OMP THREADPRIVATE(cumul_forcing)
  LOGICAL, SAVE                                  :: cumul_Cforcing=.FALSE.  !! Flag, if internal parameter cumul_Cforcing is 
                                                                            !! TRUE then ::nbyear (defined in run definition 
                                                                            !! file will be forced to 1 later in this module. If 
                                                                            !! FALSE the mean over ::nbyear is written in forcesoil
!$OMP THREADPRIVATE(cumul_Cforcing)
!---   
  REAL,ALLOCATABLE,SAVE,DIMENSION(:)      :: harvest_above_monthly   !! [??CHECK] post-processing - should be removed?
!$OMP THREADPRIVATE(harvest_above_monthly)
  REAL,ALLOCATABLE,SAVE,DIMENSION(:)      :: cflux_prod_monthly      !! [??CHECK] post-processing - should be removed?
!$OMP THREADPRIVATE(cflux_prod_monthly)
!---
  INTEGER, SAVE                               :: global_years        !! Global counter of years (year)
!$OMP THREADPRIVATE(global_years)
  LOGICAL, ALLOCATABLE, SAVE, DIMENSION(:)           :: ok_equilibrium      !! Logical array marking the points where the resolution is ok 
                                                                            !! (true/false)
!$OMP THREADPRIVATE(ok_equilibrium)
  LOGICAL, ALLOCATABLE, SAVE, DIMENSION(:)           :: carbon_eq           !! Logical array to mark the carbon pools at equilibrium ? 
                                                                            !! If true, the job stops. (true/false)
!$OMP THREADPRIVATE(carbon_eq)
  REAL, ALLOCATABLE, SAVE, DIMENSION(:)       :: nbp_accu            !! Accumulated Net Biospheric Production over the year (gC.m^2 )
!$OMP THREADPRIVATE(nbp_accu)
  REAL, ALLOCATABLE, SAVE, DIMENSION(:)       :: nbp_flux            !! Net Biospheric Production (gC.m^2.day^{-1})
!$OMP THREADPRIVATE(nbp_flux)
  REAL, ALLOCATABLE, DIMENSION(:,:,:,:)       :: matrixA             !! Matrix containing the fluxes between the carbon pools
                                                                            !! per sechiba time step 
                                                                            !! @tex $(gC.m^2.day^{-1})$ @endtex
!$OMP THREADPRIVATE(matrixA)
  REAL, ALLOCATABLE, DIMENSION(:,:,:)         :: vectorB             !! Vector containing the litter increase per sechiba time step
                                                                            !! @tex $(gC m^{-2})$ @endtex
!$OMP THREADPRIVATE(vectorB)
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:,:,:) :: MatrixV             !! Matrix containing the accumulated values of matrixA 
!$OMP THREADPRIVATE(MatrixV)
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   :: VectorU             !! Matrix containing the accumulated values of VectorB
!$OMP THREADPRIVATE(VectorU)
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:,:,:) :: MatrixW             !! Matrix containing the opposite of matrixA
!$OMP THREADPRIVATE(MatrixW)
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   :: previous_stock      !! Array containing the carbon stock calculated by the analytical
                                                                            !! method in the previous resolution
!$OMP THREADPRIVATE(previous_stock)
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   :: current_stock       !! Array containing the carbon stock calculated by the analytical
                                                                            !! method in the current resolution 
!$OMP THREADPRIVATE(current_stock)
  REAL, SAVE                                  :: eps_carbon          !! Stopping criterion for carbon pools (unitless,0-1)
!$OMP THREADPRIVATE(eps_carbon)
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  :: CN_som_litter_longterm !! Longterm CN ratio of litter and som pools (gC/gN)
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  :: CP_som_litter_longterm !! Longterm CP ratio of litter and som pools (gC/gP)

  REAL, SAVE                                 :: tau_CN_longterm      !! Counter used for calculating the longterm CN ratio of SOM and litter pools (seconds)    

  ! Functional Allocation
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:)     :: KF               !! Scaling factor to convert sapwood mass
                                                                         !! into leaf mass (m). The initial value is calculated
                                                                         !! in prescribe and updated during allocation
!$OMP THREADPRIVATE(KF)

!JC Debug MOD005
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:)     :: record_reserve_grass !!record maximum grass reserve pool
!End JC Debug MOD005
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:)     :: k_latosa_adapt   !! Leaf to sapwood area adapted for waterstress. 
                                                                         !! Adaptation takes place at the end of the year
                                                                         !! (m)
!$OMP THREADPRIVATE(k_latosa_adapt)
  REAL, ALLOCATABLE,SAVE,DIMENSION(:,:)   :: rue_longterm         !! Longterm radiation use efficiency (??units??)
!$OMP THREADPRIVATE(rue_longterm)
  REAL,SAVE                                   :: dt_forcesoil        !! Time step of soil forcing file (days)
!$OMP THREADPRIVATE(dt_forcesoil)



  INTEGER,PARAMETER                           :: nparanmax=366       !! Maximum number of time steps per year for forcesoil
  INTEGER,SAVE                                :: nparan              !! Number of time steps per year for forcesoil read from run definition (unitless) 
!$OMP THREADPRIVATE(nparan)
  INTEGER,SAVE                                :: nbyear=1            !! Number of years saved for forcesoil (unitless) 
!$OMP THREADPRIVATE(nbyear)
  INTEGER,SAVE                                :: iatt                !! Time step of forcing of soil processes (iatt = 1 to ::nparan*::nbyear) 
!$OMP THREADPRIVATE(iatt)
  INTEGER,SAVE                                :: iatt_old=1          !! Previous ::iatt
!$OMP THREADPRIVATE(iatt_old)
  INTEGER,SAVE                                :: nsfm                !! Number of time steps that can be stored in memory (unitless) 
!$OMP THREADPRIVATE(nsfm)
  INTEGER,SAVE                                :: nsft                !! Number of time steps in a year (unitless)
!$OMP THREADPRIVATE(nsft)
  INTEGER,SAVE                                :: iisf                !! Current pointer for teststomate (unitless)
!$OMP THREADPRIVATE(iisf)
  CHARACTER(LEN=100), SAVE                           :: forcing_name        !! Name of forcing file 1
!$OMP THREADPRIVATE(forcing_name)
  CHARACTER(LEN=100), SAVE                           :: Cforcing_name       !! Name of forcing file 2
!$OMP THREADPRIVATE(Cforcing_name)
  INTEGER,SAVE                                :: Cforcing_id         !! File identifer of file 2
!$OMP THREADPRIVATE(Cforcing_id)    
  INTEGER,PARAMETER                           :: ndm = 10            !! Maximum number of dimensions (unitless)

  REAL,ALLOCATABLE,DIMENSION(:)               :: max_eau_var_daily   !! Daily Maximum water content of the soil @tex $(kg m^{-2})$ @endtex
!$OMP THREADPRIVATE(Cforcing_id)    

  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:,:)    :: litter_avail          !! edible litter for animals
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:,:)    :: litter_not_avail      !! litter not edible for animals
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:,:)    :: litter_avail_frac     !! edible litter fraction
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:)      :: sla_calc              !! leaf age-related SLA
  REAL,ALLOCATABLE,SAVE,DIMENSION(:)        :: t2m_14                !! "14days" 2 meter temperatures (K)
  REAL,ALLOCATABLE, SAVE, DIMENSION(:,:)   :: wshtotsum             !! accumulated harvested grass biomass
  REAL,ALLOCATABLE, SAVE, DIMENSION(:,:)    :: sr_ugb                !! Optimised stocking density (animal m-2)
  REAL,ALLOCATABLE, SAVE, DIMENSION(:,:)    :: compt_ugb             !! Counter of grazing days (d)
  REAL,ALLOCATABLE, SAVE, DIMENSION(:,:)    :: nb_ani                !! optimal livestock density
  REAL,ALLOCATABLE, SAVE, DIMENSION(:,:)    :: grazed_frac           !! optimal fraction for grazing (in contract to mowning)
  REAL,ALLOCATABLE, SAVE, DIMENSION(:,:)    :: import_yield          !! annual total harvested grass biomass yield
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:)      :: sla_age1              !! maximum SLA for youngest leaves
  REAL,DIMENSION(:,:), ALLOCATABLE, SAVE    :: when_growthinit_cut   !! growing days after grass harvest
  REAL,ALLOCATABLE, SAVE, DIMENSION(:,:)    :: nb_grazingdays        !! annual total number of days animal grazing in field
  REAL,ALLOCATABLE,SAVE,DIMENSION(:)        :: snowfall_daily        !! daily snow fall
  REAL,ALLOCATABLE,SAVE,DIMENSION(:)        :: snowmass_daily        !! daily snow mass
  REAL,ALLOCATABLE,SAVE,DIMENSION(:)        :: tmc_topgrass_daily    !! Daily top 5 layer soil moisture (m^3 m^-3)
  REAL,ALLOCATABLE,SAVE,DIMENSION(:)        :: after_snow            !! day counter after snow melt to prevent grazing
  REAL,ALLOCATABLE,SAVE,DIMENSION(:)        :: after_wet             !! day counter after wet soil is detected to prevent grazing
  REAL,ALLOCATABLE,SAVE,DIMENSION(:)        :: wet1day               !! accumulated days with wet soil
  REAL,ALLOCATABLE,SAVE,DIMENSION(:)        :: wet2day               !! accumulated days with wet soil
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:,:)    :: grm_nfert             !! N fertilizer applied to grassland in GRM module
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:)      :: grm_pfert             !! P fertilizer applied to grassland in GRM module
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:)      :: BNF_clover            !! clover BNF in grass-legume mixtures
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:)      :: BNF_clover_daily      !! daily clover BNF in grass-legume mixtures
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:)      :: GRM_devstage          !! development stage of grasses
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:)      :: p_input_dt            !! P input (including weathering, deposition etc.) per time step
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:,:)    :: agri_nfert            !! N fertilization by agricultural practices
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:)      :: agri_pfert            !! P fertilization by agricultural practices
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:)      :: rest_Ngrazing         !! the residual N applied to pasture/rangeland that did not output from biomass grazed
  REAL,ALLOCATABLE,SAVE,DIMENSION(:)        :: grazed_above          !! total grazed biomass in simplified grazing module
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:)      :: days_senescence       !! day counter after senescence is detected
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:)      :: fclover               !! fraction of clover in grass/legume mixture

PUBLIC clay_fm, silt_fm, bulk_fm, humrel_daily_fm, litterhum_daily_fm, t2m_daily_fm, &
   & t2m_min_daily_fm, tsurf_daily_fm, tsoil_daily_fm, soilhum_daily_fm, &
   & precip_fm, gpp_daily_fm, veget_fm, veget_max_fm, lai_fm, max_eau_var_daily_fm
! PUBLIC  dt_days, date, do_slow, ok_equilibrium
   PUBLIC  dt_days, do_slow, ok_equilibrium
PUBLIC isf, days_since_beg, nf_written

CONTAINS
  

!! ================================================================================================================================
!! SUBROUTINE 	: stomate_initialize
!!
!>\BRIEF        Initialization routine for stomate module. 
!!
!! DESCRIPTION  : Initialization routine for stomate module. Read options from parameter file, allocate variables, read variables 
!!                from restart file and initialize variables if necessary. 
!!                
!! \n
!_ ================================================================================================================================

SUBROUTINE stomate_initialize &
        (kjit,           kjpij,             kjpindex,                        &
         rest_id_stom,   hist_id_stom,      hist_id_stom_IPCC,               &
         index,          lalo,              neighbours,   resolution,        &
         contfrac,       totfrac_nobio,     clay,         silt,              &
         bulk,           t2m,                                                &
         lai,            veget,             veget_max, soil_order_imposed,   &
         co2_flux,       fco2_lu,           deadleaf_cover,  assim_param, temp_growth)

    IMPLICIT NONE
    !! 0. Variable and parameter declaration
    !! 0.1 Input variables
    INTEGER,INTENT(in)                       :: kjit              !! Time step number (unitless)
    INTEGER,INTENT(in)                       :: kjpij             !! Total size of the un-compressed grid (unitless)
    INTEGER,INTENT(in)                       :: kjpindex          !! Domain size - terrestrial pixels only (unitless)
    INTEGER,INTENT(in)                       :: rest_id_stom      !! STOMATE's _Restart_ file identifier (unitless)
    INTEGER,INTENT(in)                       :: hist_id_stom      !! STOMATE's _history_ file identifier (unitless)
    INTEGER,INTENT(in)                       :: hist_id_stom_IPCC !! STOMATE's IPCC _history_ file identifier(unitless) 
    INTEGER,DIMENSION(kjpindex),INTENT(in)   :: index             !! The indices of the terrestrial pixels only (unitless) 
    REAL,DIMENSION(kjpindex,2),INTENT(in)    :: lalo              !! Geographical coordinates (latitude,longitude) for pixels (degrees) 
    INTEGER,DIMENSION(kjpindex,NbNeighb),INTENT(in) :: neighbours !! Neighoring grid points if land for the DGVM (unitless) 
    REAL,DIMENSION(kjpindex,2),INTENT(in)    :: resolution        !! Size in x an y of the grid (m) - surface area of the gridbox 
    REAL,DIMENSION (kjpindex), INTENT (in)   :: contfrac          !! Fraction of continent in the grid cell (unitless)
    REAL,DIMENSION(kjpindex),INTENT(in)      :: totfrac_nobio     !! Fraction of grid cell covered by lakes, land ice, cities, ... (unitless) 
    REAL,DIMENSION(kjpindex),INTENT(in)      :: clay              !! Clay fraction of soil (0-1, unitless)
    REAL,DIMENSION(kjpindex),INTENT(in)      :: silt              !! Silt fraction of soil (0-1, unitless) 
    REAL,DIMENSION(kjpindex),INTENT(in)      :: bulk              !! Bulk density (kg/m**3) 
    REAL,DIMENSION(kjpindex),INTENT(in)      :: t2m               !! 2 m air temperature (K)
    REAL,DIMENSION(kjpindex,nvm),INTENT(in)  :: lai               !! Leaf area inex @tex $(m^2 m^{-2})$ @endtex
    REAL,DIMENSION(kjpindex,nvm),INTENT(in)  :: veget             !! Fraction of vegetation type including 
                                                                         !! non-biological fraction (unitless) 
    REAL,DIMENSION(kjpindex,nvm),INTENT(in)  :: veget_max         !! Maximum fraction of vegetation type including 
                                                                         !! non-biological fraction (unitless) 

    REAL,DIMENSION(kjpindex),INTENT(in)      :: soil_order_imposed !! USDA soil order we impose globally

    !! 0.2 Output variables
    REAL,DIMENSION(kjpindex,nvm),INTENT(out) :: co2_flux          !! CO2 flux between atmosphere and biosphere
    REAL,DIMENSION(kjpindex),INTENT(out)     :: fco2_lu           !! CO2 flux between atmosphere and biosphere from land-use (without forest management)  
    REAL,DIMENSION(kjpindex),INTENT(out)     :: deadleaf_cover    !! Fraction of soil covered by dead leaves (unitless)
    REAL,DIMENSION(kjpindex,nvm,npco2),INTENT(out) :: assim_param !! min+max+opt temperatures (K) & vmax for photosynthesis  
                                                                         !! @tex $(\mu mol m^{-2}s^{-1})$ @endtex  
    REAL,DIMENSION(kjpindex),INTENT(out)     :: temp_growth       !! Growth temperature (Â°C)  
                                                                         !! Is equal to t2m_month 

    !! 0.3 Local variables
    REAL                                   :: dt_days_read             !! STOMATE time step read in restart file (days)
    INTEGER                                :: l,k,ji, jv, i, j, m      !! indices    
    REAL,PARAMETER                         :: max_dt_days = 5.         !! Maximum STOMATE time step (days)
    REAL,DIMENSION(kjpindex,nvm)           :: rprof                    !! Coefficient of the exponential functions that 
                                                                              !! relates root density to soil depth (unitless) 
    REAL,DIMENSION(kjpindex,nvm)           :: gpp_daily_x              !! "Daily" gpp for teststomate  
                                                                              !! @tex $(??gC m^{-2} dt_stomate^{-1})$ @endtex 
    REAL,DIMENSION(kjpindex,nvm)           :: veget_cov                !! Fractional coverage: actually share of the pixel 
                                                                              !! covered by a PFT (fraction of ground area), 
                                                                              !! taking into account LAI ??(= grid scale fpc)?? 
    INTEGER                                :: ier                      !! Check errors in netcdf call (unitless)

    INTEGER                                :: max_totsize              !! Memory management - maximum memory size (Mb)
    INTEGER                                :: totsize_1step            !! Memory management - memory required to store one 
                                                                              !! time step on one processor (Mb) 
    INTEGER                                :: totsize_tmp              !! Memory management - memory required to store one 
                                                                              !! time step on all processors(Mb) 
    INTEGER                                :: vid                      !! Variable identifer of netCDF (unitless)
    INTEGER                                :: nneigh                   !! Number of neighbouring pixels
    INTEGER                                :: direct                   !! 
    INTEGER,DIMENSION(ndm)                 :: d_id                     !! 

    !JC add for reading NPinput
    REAL, ALLOCATABLE,DIMENSION(:,:,:)        :: manage1
    REAL, ALLOCATABLE,DIMENSION(:,:,:)        :: manage2
    REAL, ALLOCATABLE,DIMENSION(:,:,:)        :: manage3
    REAL, ALLOCATABLE,DIMENSION(:,:,:)        :: manage4
    REAL, ALLOCATABLE,DIMENSION(:,:,:)        :: manage5
    REAL, ALLOCATABLE,DIMENSION(:,:,:)        :: manage6
    INTEGER                                        :: yrlen
    CHARACTER(LEN=30)                                     :: strManage
    CHARACTER(LEN=30)                                     :: strVar1
    CHARACTER(LEN=30)                                     :: strVar2
    CHARACTER(LEN=30)                                     :: strVar3
    CHARACTER(LEN=30)                                     :: strVar4
    CHARACTER(LEN=30)                                     :: strVar5
    CHARACTER(LEN=30)                                     :: strVar6
    !End JC add NPinput
!_ ================================================================================================================================
    
    !! 1. Initialize variable
    !! Update flag
    l_first_stomate = .FALSE.
    
    !! 1.1 Store current time step in a common variable
    itime = kjit
    
    !![DISPENSABLE] 1.2 Copy the depth of the different soil layers from diaglev specified in slow_proc
    
    !! 1.3 PFT rooting depth across pixels, humescte is pre-defined 
    ! (constantes_veg.f90). It is defined as the coefficient of an exponential 
    ! function relating root density to depth 
    DO j=1,nvm
       rprof(:,j) = 1./humcste(j)
    ENDDO
    
    !! 1.4.0 Parameters for spinup
    !
    eps_carbon = 0.01
    !Config Key   = EPS_CARBON
    !Config Desc  = Allowed error on carbon stock
    !Config If    = SPINUP_ANALYTIC
    !Config Def   = 0.01
    !Config Help  = 
    !Config Units = [%]   
    CALL getin_p('EPS_CARBON',eps_carbon)       
    
    
    !Config Key   = SPINUP_PERIOD
    !Config Desc  = Period to calulcate equilibrium during spinup analytic
    !Config If    = SPINUP_ANALYTIC
    !Config Def   = -1
    !Config Help  = Period corresponds in most cases to the number of years of forcing data used in the spinup.
    !Config Units = [years]   
    spinup_period = -1
    CALL getin_p('SPINUP_PERIOD',spinup_period)       

    ! Check spinup_period values. 
    ! For periods uptil 6 years, to obtain equilibrium, a bigger period have to be used 
    ! and therefore spinup_period is adjusted to 10 years. 
    IF (spinup_analytic) THEN
       IF (spinup_period <= 0) THEN
          WRITE(numout,*) 'Error in parameter spinup_period. This parameter must be > 0 : spinup_period=',spinup_period
          CALL ipslerr_p (3,'stomate_initialize', &
               'Parameter spinup_period must be set to a positive integer.', &
               'Set this parameter to the number of years of forcing data used for the spinup.', &
               '')
       ELSE IF (spinup_period <= 6) THEN
          ! Adjust to bigger period. The period must be a multiple of the original period.
          WRITE(numout,*) 'Initial spinup_period =',spinup_period,' will be not adjusted.'
       END IF
       IF (printlev >=1) WRITE(numout,*) 'Spinup analytic is activated using eps_carbon=',&
            eps_carbon, ' and spinup_period=',spinup_period
    END IF
    
    !! 1.4.0 Initialization of PFT specific parameters 
    ! Initialization of PFT specific parameters that have no value
    ! for the bare soil PFT i.e. fire resistance, flamability, maximum lai,
    ! settings for growing degree days (GDD), settings for senescence, 
    ! respiration coefficients, photosynthesis, etc.
    ! [DISPENSABLE]
    
    !! 1.4.1 Allocate memory for all variables in stomate
    ! Allocate memory for all variables in stomate, build new index
    ! tables accounting for the PFTs, read and check flags and set file
    ! identifier for restart and history files.
    CALL stomate_init (kjpij, kjpindex, index, lalo, &
         rest_id_stom, hist_id_stom, hist_id_stom_IPCC)
    
    !! 1.4.2 Initialization of PFT specific parameters
    ! Initialization of PFT specific parameters i.e. sla from leaf life, 
    ! sapling characteristics (biomass), migration speed, critical diameter,
    ! coldest tolerable temperature, critical values for phenology, maximum
    ! life time of leaves, respiration coefficients and photosynthesis.
    ! The subroutine also communicates settings read by stomate_constant_init.
    CALL data (kjpindex, lalo)
    
    !! 1.4.3 Initial conditions
    
    !! 1.4.3.1 Read initial values for STOMATE's variables from the _restart_ file
    ! ??Shouldn't this be included in stomate_init?? Looks like an initialization!
    co2_flux(:,:) = zero
    fco2_lu(:) = zero
    
    ! Get values from _restart_ file. Note that only ::kjpindex, ::index, ::lalo 
    ! and ::resolution are input variables, all others are output variables.
    CALL readstart &
         (kjpindex, index, lalo, resolution, t2m, &
         dt_days_read, days_since_beg, &
         ind, adapted, regenerate, &
         humrel_daily, max_eau_var_daily, gdd_init_date, litterhum_daily, &
         t2m_daily, t2m_min_daily, tsurf_daily, tsoil_daily, &
         soilhum_daily, &
         drainage_longterm, runoff_longterm,                     &
         precip_daily, &
         gpp_daily, npp_daily, turnover_daily, &
         humrel_month, humrel_week, moiavail_growingseason,&
         t2m_longterm, tau_longterm, t2m_month, t2m_week, &
         tsoil_month, soilhum_month, fireindex, firelitter, &
         maxhumrel_lastyear, maxhumrel_thisyear, &
         minhumrel_lastyear, minhumrel_thisyear, &
         maxgppweek_lastyear, maxgppweek_thisyear, &
         gdd0_lastyear, gdd0_thisyear, &
         precip_lastyear, precip_thisyear, &
         gdd_m5_dormance,  gdd_from_growthinit, gdd_midwinter, ncd_dormance, ngd_minus5, &
         PFTpresent, npp_longterm, lm_lastyearmax, lm_thisyearmax, &
         maxfpc_lastyear, maxfpc_thisyear, &
         turnover_longterm, gpp_week, npp_week, biomass, resp_maint_part, &
         leaf_age, leaf_frac, &
         senescence, when_growthinit, age, &
         resp_hetero_d, resp_maint_d, resp_growth_d, co2_fire, co2_to_bm_dgvm, &
         n_to_bm, p_to_bm,veget_lastlight, everywhere, need_adjacent, RIP_time, &
         time_hum_min, hum_min_dormance, &
         litter, dead_leaves, &
         som, lignin_struc, lignin_wood,turnover_time,&
         prod10,prod100,flux10, flux100, &
         convflux, cflux_prod10, cflux_prod100, &
         prod10_harvest,prod100_harvest,flux10_harvest, flux100_harvest, &
         convflux_harvest, cflux_prod10_harvest, cflux_prod100_harvest, &
         woodharvestpft, bm_to_litter, carb_mass_total, &
         Tseason, Tseason_length, Tseason_tmp, &
         Tmin_spring_time, begin_leaves, onset_date, &
         global_years, ok_equilibrium, nbp_accu, nbp_flux, &
         MatrixV, VectorU, previous_stock, current_stock, assim_param, CN_som_litter_longterm,tau_CN_longterm, KF, k_latosa_adapt, &
         CP_som_litter_longterm, np_leaf_avg_season, soil_p_min, &
         f_Pdissolved, EW_grain, &
         rue_longterm,cn_leaf_avg_season,nstress_season, pstress_season, &
         lai_target_longterm, &
         soil_n_min,p_O2,bact, &
         record_reserve_grass, &
         wshtotsum, sr_ugb, sla_calc, nb_ani, grazed_frac, &
         import_yield, t2m_14, litter_not_avail, nb_grazingdays, &
         after_snow, after_wet, wet1day, wet2day, GRM_devstage , &
         days_senescence)
         
    !! 1.4.4.1 If phosphorus cycle is active, lithology is read in
    IF (ok_pcycle) THEN
       ALLOCATE(lith_frac(kjpindex,nlcm))
       ALLOCATE(soil_shield(kjpindex))
       ALLOCATE(soil_orders(kjpindex))
       CALL readlithology (kjpindex, lalo, resolution, soil_order_imposed,lith_frac, soil_shield, soil_orders)

    ENDIF

    !! 1.4.5 Check time step
       
    !! 1.4.5.1 Allow STOMATE's time step to change although this is dangerous
    IF (dt_days /= dt_days_read) THEN
       WRITE(numout,*) 'slow_processes: STOMATE time step changes:', &
            & dt_days_read,' -> ',dt_days
    ENDIF
    
    !! 1.4.5.2 Time step has to be a multiple of a full day
    IF ( ( dt_days-REAL(NINT(dt_days),r_std) ) > min_stomate ) THEN
       WRITE(numout,*) 'slow_processes: STOMATE time step is not a mutiple of a full day:', &
            & dt_days,' days.'
       STOP
    ENDIF
    
    !! 1.4.5.3 upper limit to STOMATE's time step
    IF ( dt_days > max_dt_days ) THEN
       WRITE(numout,*) 'slow_processes: STOMATE time step exceeds the maximum value:', &
            & dt_days,' days > ', max_dt_days, ' days.'  
       STOP
    ENDIF
    
    !! 1.4.5.4 STOMATE time step must not be less than the forcing time step
    IF ( dt_sechiba > dt_days*one_day ) THEN
       WRITE(numout,*) &
            & 'slow_processes: STOMATE time step ::dt_days smaller than forcing time step ::dt_sechiba'
       STOP
    ENDIF
    
    !! 1.4.5.6 Final message on time step
    IF (printlev >=2) WRITE(numout,*) 'Slow_processes, STOMATE time step (days): ', dt_days
    
    !! 1.4.6 Write forcing file for teststomate
    IF (ok_co2 .AND. allow_forcing_write) THEN
       
       !Config Key   = STOMATE_FORCING_NAME
       !Config Desc  = Name of STOMATE's forcing file
       !Config If    = OK_STOMATE
       !Config Def   = NONE
       !Config Help  = Name that will be given
       !Config         to STOMATE's offline forcing file
       !Config         Compatible with Nicolas Viovy's driver
       !Config Units = [FILE]
       forcing_name = stomate_forcing_name
       CALL getin_p('STOMATE_FORCING_NAME',forcing_name)
       
       IF ( TRIM(forcing_name) /= 'NONE' ) THEN
          
          !! 1.4.6.1 Calculate steps that can be stored in memory
          ! Action for the root processor only (parallel computing)  
          IF (is_root_prc) CALL SYSTEM ('rm -f '//TRIM(forcing_name))
          IF (printlev>=2) WRITE(numout,*) 'writing a forcing file for STOMATE.'
          
          !Config Key   = STOMATE_FORCING_MEMSIZE
          !Config Desc  = Size of STOMATE forcing data in memory 
          !Config If    = OK_STOMATE
          !Config Def   = 50
          !Config Help  = This variable determines how many
          !Config         forcing states will be kept in memory.
          !Config         Must be a compromise between memory
          !Config         use and frequeny of disk access.
          !Config Units = [MegaBytes]
          max_totsize = 50
          CALL getin_p('STOMATE_FORCING_MEMSIZE', max_totsize)      
          max_totsize = max_totsize*1000000
          
          totsize_1step = &
               SIZE(clay)*KIND(clay) &
               +SIZE(silt)*KIND(silt) &
               +SIZE(bulk)*KIND(bulk) &
               +SIZE(humrel_daily)*KIND(humrel_daily) &
               +SIZE(litterhum_daily)*KIND(litterhum_daily) &
               +SIZE(t2m_daily)*KIND(t2m_daily) &
               +SIZE(t2m_min_daily)*KIND(t2m_min_daily) &
               +SIZE(tsurf_daily)*KIND(tsurf_daily) &
               +SIZE(tsoil_daily)*KIND(tsoil_daily) &
               +SIZE(soilhum_daily)*KIND(soilhum_daily) &
               +SIZE(precip_daily)*KIND(precip_daily) &
               +SIZE(gpp_daily_x)*KIND(gpp_daily_x) &
               +SIZE(veget)*KIND(veget) &
               +SIZE(veget_max)*KIND(veget_max) &
               +SIZE(lai)*KIND(lai) &    
               +SIZE(runoff_daily)*KIND(runoff_daily) &
               +SIZE(drainage_daily)*KIND(drainage_daily) 
          
          ! Totsize_1step is the size on a single processor, sum
          ! all processors and send to all processors
          CALL reduce_sum(totsize_1step,totsize_tmp)
          CALL bcast(totsize_tmp)
          totsize_1step=totsize_tmp
          
          ! Total number of forcing steps
          nsft = INT(one_year/(dt_stomate/one_day))
          
          ! Number of forcing steps in memory
          nsfm = MIN(nsft, &
               MAX(1,NINT( REAL(max_totsize,r_std) &
               /REAL(totsize_1step,r_std))))
            
             
          !! 1.6.4.2 Allocate memory for variables containing forcing data  
          ! and initialize variables (set to zero).
          CALL init_forcing (kjpindex,nsfm,nsft)
          
          ! Indexing for writing forcing file
          isf(:) = (/ (i,i=1,nsfm) /)
          nf_written(:) = .FALSE.
          nf_cumul(:) = 0
          iisf = 0
          
          !! 1.6.4.3 Create netcdf file
          ! Create, define and populate a netcdf file containing the forcing data.
          ! For the root processor only (parallel computing). NF90_ are functions
          ! from and external library.  
          IF (is_root_prc) THEN
             
             ! Create new netCDF dataset
             ier = NF90_CREATE (TRIM(forcing_name),NF90_SHARE,forcing_id)
             
             ! Add variable attribute
             ! Note ::iim_g and ::jjm_g are dimensions of the global field and 
             ! ::nbp_glo is the number of global continental points
             ier = NF90_PUT_ATT (forcing_id,NF90_GLOBAL,'dt_sechiba',dt_sechiba)
             ier = NF90_PUT_ATT (forcing_id,NF90_GLOBAL,'dt_stomate',dt_stomate)
             ier = NF90_PUT_ATT (forcing_id,NF90_GLOBAL, &
                  'nsft',REAL(nsft,r_std))
             ier = NF90_PUT_ATT (forcing_id,NF90_GLOBAL, &
                  'kjpij',REAL(iim_g*jjm_g,r_std))
             ier = NF90_PUT_ATT (forcing_id,NF90_GLOBAL, &
                  'kjpindex',REAL(nbp_glo,r_std))
             
             ! Add new dimension
             ier = NF90_DEF_DIM (forcing_id,'points',nbp_glo,d_id(1))
             ier = NF90_DEF_DIM (forcing_id,'layers',nslm,d_id(2))
             ier = NF90_DEF_DIM (forcing_id,'pft',nvm,d_id(3))
             direct=2
             ier = NF90_DEF_DIM (forcing_id,'direction',direct,d_id(4))
             nneigh=8
             ier = NF90_DEF_DIM (forcing_id,'nneigh',nneigh,d_id(5))
             ier = NF90_DEF_DIM (forcing_id,'time',NF90_UNLIMITED,d_id(6))
             ier = NF90_DEF_DIM (forcing_id,'nbparts',nparts,d_id(7))
             
             ! Add new variable
             ier = NF90_DEF_VAR (forcing_id,'points',    r_typ,d_id(1),vid)
             ier = NF90_DEF_VAR (forcing_id,'layers',    r_typ,d_id(2),vid)
             ier = NF90_DEF_VAR (forcing_id,'pft',       r_typ,d_id(3),vid)
             ier = NF90_DEF_VAR (forcing_id,'direction', r_typ,d_id(4),vid)
             ier = NF90_DEF_VAR (forcing_id,'nneigh',    r_typ,d_id(5),vid)
             ier = NF90_DEF_VAR (forcing_id,'time',      r_typ,d_id(6),vid)
             ier = NF90_DEF_VAR (forcing_id,'nbparts',   r_typ,d_id(7),vid)
             ier = NF90_DEF_VAR (forcing_id,'index',     r_typ,d_id(1),vid)
             ier = NF90_DEF_VAR (forcing_id,'contfrac',  r_typ,d_id(1),vid) 
             ier = NF90_DEF_VAR (forcing_id,'lalo', &
                  r_typ,(/ d_id(1),d_id(4) /),vid)
             ier = NF90_DEF_VAR (forcing_id,'neighbours', &
                  r_typ,(/ d_id(1),d_id(5) /),vid)
             ier = NF90_DEF_VAR (forcing_id,'resolution', &
                  r_typ,(/ d_id(1),d_id(4) /),vid)
             ier = NF90_DEF_VAR (forcing_id,'clay', &
                  r_typ,(/ d_id(1),d_id(6) /),vid)
             ier = NF90_DEF_VAR (forcing_id,'silt', & 
                  r_typ,(/ d_id(1),d_id(6) /),vid) 
             ier = NF90_DEF_VAR (forcing_id,'bulk', & 
                  r_typ,(/ d_id(1),d_id(6) /),vid) 
             ier = NF90_DEF_VAR (forcing_id,'humrel', &
                  r_typ,(/ d_id(1),d_id(3),d_id(6) /),vid)
             ier = NF90_DEF_VAR (forcing_id,'litterhum', &
                  r_typ,(/ d_id(1),d_id(6) /),vid)
             ier = NF90_DEF_VAR (forcing_id,'t2m', &
                  r_typ,(/ d_id(1),d_id(6) /),vid)
             ier = NF90_DEF_VAR (forcing_id,'t2m_min', &
                  r_typ,(/ d_id(1),d_id(6) /),vid)
             ier = NF90_DEF_VAR (forcing_id,'tsurf', &
                  r_typ,(/ d_id(1),d_id(6) /),vid)
             ier = NF90_DEF_VAR (forcing_id,'tsoil', &
                  r_typ,(/ d_id(1),d_id(2),d_id(6) /),vid)
             ier = NF90_DEF_VAR (forcing_id,'soilhum', &
                  r_typ,(/ d_id(1),d_id(2),d_id(6) /),vid)
             ier = NF90_DEF_VAR (forcing_id,'precip', &
                  r_typ,(/ d_id(1),d_id(6) /),vid)
             ier = NF90_DEF_VAR (forcing_id,'gpp', &
                  r_typ,(/ d_id(1),d_id(3),d_id(6) /),vid)
             ier = NF90_DEF_VAR (forcing_id,'veget', &
                  r_typ,(/ d_id(1),d_id(3),d_id(6) /),vid)
             ier = NF90_DEF_VAR (forcing_id,'veget_max', &
                  r_typ,(/ d_id(1),d_id(3),d_id(6) /),vid)
             ier = NF90_DEF_VAR (forcing_id,'lai', &
                  r_typ,(/ d_id(1),d_id(3),d_id(6) /),vid)
             ier = NF90_DEF_VAR (forcing_id,'drainage', & 
                  r_typ,(/ d_id(1),d_id(3),d_id(6) /),vid) 
             ier = NF90_DEF_VAR (forcing_id,'runoff', & 
                  r_typ,(/ d_id(1),d_id(3),d_id(6) /),vid) 
             ier = NF90_ENDDEF (forcing_id)
             
             ! Given the name of a varaible, nf90_inq_varid finds the variable 
             ! ID (::vid). Put data value(s) into variable ::vid
             ier = NF90_INQ_VARID (forcing_id,'points',vid)
             ier = NF90_PUT_VAR (forcing_id,vid, &
                  (/(REAL(i,r_std),i=1,nbp_glo) /))
             ier = NF90_INQ_VARID (forcing_id,'layers',vid)
             ier = NF90_PUT_VAR (forcing_id,vid,(/(REAL(i,r_std),i=1,nslm)/))
             ier = NF90_INQ_VARID (forcing_id,'pft',vid)
             ier = NF90_PUT_VAR (forcing_id,vid,(/(REAL(i,r_std),i=1,nvm)/))
             ier = NF90_INQ_VARID (forcing_id,'direction',vid)
             ier = NF90_PUT_VAR (forcing_id,vid,(/(REAL(i,r_std),i=1,2)/))
             ier = NF90_INQ_VARID (forcing_id,'nneigh',vid)
             ier = NF90_PUT_VAR (forcing_id,vid,(/(REAL(i,r_std),i=1,8)/))
             ier = NF90_INQ_VARID (forcing_id,'time',vid)
             ier = NF90_PUT_VAR (forcing_id,vid,(/(REAL(i,r_std),i=1,nsft)/))
             ier = NF90_INQ_VARID (forcing_id,'nbparts',vid)
             ier = NF90_PUT_VAR (forcing_id,vid,(/(REAL(i,r_std),i=1,nparts)/))
             ier = NF90_INQ_VARID (forcing_id,'index',vid)  
             ier = NF90_PUT_VAR (forcing_id,vid,REAL(index_g,r_std))
             ier = NF90_INQ_VARID (forcing_id,'contfrac',vid)
             ier = NF90_PUT_VAR (forcing_id,vid,REAL(contfrac_g,r_std))
             ier = NF90_INQ_VARID (forcing_id,'lalo',vid)
             ier = NF90_PUT_VAR (forcing_id,vid,lalo_g)
             !ym attention a neighbours, a modifier plus tard      
             ier = NF90_INQ_VARID (forcing_id,'neighbours',vid)
             ier = NF90_PUT_VAR (forcing_id,vid,REAL(neighbours_g,r_std))
             ier = NF90_INQ_VARID (forcing_id,'resolution',vid)
             ier = NF90_PUT_VAR (forcing_id,vid,resolution_g)
          ENDIF ! is_root_prc
       ENDIF ! (forcing_name) /= 'NONE'
    ENDIF ! ok_co2 =.TRUE.
    
    !! 1.4.7 write forcing file for forcesoil
    !! 1.4.7.1 Initialize
    !Config Key   = STOMATE_CFORCING_NAME
    !Config Desc  = Name of STOMATE's carbon forcing file
    !Config If    = OK_STOMATE
    !Config Def   = NONE
    !Config Help  = Name that will be given to STOMATE's carbon
    !Config         offline forcing file
    !Config         Compatible with Nicolas Viovy's driver
    !Config Units = [FILE]
    Cforcing_name = stomate_Cforcing_name
    CALL getin_p('STOMATE_CFORCING_NAME',Cforcing_name)
    
    IF ( TRIM(Cforcing_name) /= 'NONE' ) THEN
       
       ! Time step of forcesoil
       !Config Key   = FORCESOIL_STEP_PER_YEAR
       !Config Desc  = Number of time steps per year for carbon spinup.
       !Config If    = OK_STOMATE
       !Config Def   = 365
       !Config Help  = Number of time steps per year for carbon spinup.
       !Config Units = [days, months, year]
       nparan = 365
       CALL getin_p('FORCESOIL_STEP_PER_YEAR', nparan)
       
       ! Correct if setting is out of bounds 
       IF ( nparan < 1 ) nparan = 1
       
       !Config Key   = FORCESOIL_NB_YEAR
       !Config Desc  = Number of years saved for carbon spinup.
       !Config If    = OK_STOMATE
       !Config Def   = 1
       !Config Help  = Number of years saved for carbon spinup. If internal parameter cumul_Cforcing is TRUE in stomate.f90
       !Config         Then this parameter is forced to one.
       !Config Units = [years]
       CALL getin_p('FORCESOIL_NB_YEAR', nbyear)
       
       ! Set ::nbyear to 1. if ::cumul_Cforcing=.TRUE.
       IF ( cumul_Cforcing ) THEN
          CALL ipslerr_p (1,'stomate', &
               'Internal parameter cumul_Cforcing is TRUE in stomate.f90', &
               'Parameter FORCESOIL_NB_YEAR is therefore forced to 1.', &
               '::nbyear is thus set to 1.')
          nbyear=1
       ENDIF
       
       ! Make use of ::nparan to calculate ::dt_forcesoil
       dt_forcesoil = zero
       nparan = nparan+1
       DO WHILE ( dt_forcesoil < dt_stomate/one_day )
          nparan = nparan-1
          IF ( nparan < 1 ) THEN
             STOP 'Problem with number of soil forcing time steps ::nparan < 1.'
          ENDIF
          dt_forcesoil = one_year/REAL(nparan,r_std)
       ENDDO
       IF ( nparan > nparanmax ) THEN
          STOP 'Problem with number of soil forcing time steps ::nparan > ::nparanmax'
       ENDIF
       IF (printlev>=2) WRITE(numout,*) 'Time step of soil forcing (d): ',dt_forcesoil
       
       ! Allocate memory for the forcing variables of soil dynamics
       ALLOCATE( nforce(nparan*nbyear))
       nforce(:) = 0
       ALLOCATE(control_moist(kjpindex,nlevs,nparan*nbyear))
       ALLOCATE(npp_equil(kjpindex,nparan*nbyear))
       ALLOCATE(npp_tot(kjpindex))
       ALLOCATE(control_temp(kjpindex,nlevs,nparan*nbyear))
       ALLOCATE(carbon_input(kjpindex,ncarb,nvm,nparan*nbyear)) 
       ALLOCATE(nitrogen_input(kjpindex,ncarb,nvm,nparan*nbyear)) 
       ALLOCATE(drainage(kjpindex,nvm,nparan*nbyear))         
       ALLOCATE(runoff  (kjpindex,nvm,nparan*nbyear))         

       ! Initialize variables, set to zero
       control_moist(:,:,:) = zero
       npp_equil(:,:) = zero
       npp_tot(:) = zero
       control_temp(:,:,:) = zero
       carbon_input(:,:,:,:) = zero
       nitrogen_input(:,:,:,:) = zero
       drainage(:,:,:) = zero 
       runoff  (:,:,:) = zero 
    ENDIF ! Cforcing_name) /= 'NONE'
    
    !! 1.4.8 Calculate STOMATE's vegetation fractions from veget, veget_max
    DO j=1,nvm
       WHERE ((1.-totfrac_nobio(:)) > min_sechiba)       
          ! Pixels with vegetation
          veget_cov(:,j) = veget(:,j)/( 1.-totfrac_nobio(:) )
          veget_cov_max(:,j) = veget_max(:,j)/( 1.-totfrac_nobio(:) )
       ELSEWHERE
          ! Pixels without vegetation
          veget_cov(:,j) = zero
          veget_cov_max(:,j) = zero
       ENDWHERE
    ENDDO ! Loop over PFTs

    !! 1.4.9 Initialize non-zero variables
    CALL stomate_var_init &
         (kjpindex, veget_cov_max, leaf_age, leaf_frac, &
         dead_leaves, &
         veget, lai, deadleaf_cover, assim_param)
    
    ! Initialize land cover change variable
    ! ??Should be integrated in the subroutine?? 
    harvest_above(:) = zero
    harvest_bioenergy(:,:) = zero
    grazed_above(:) = zero
 
    ! Initialize temp_growth
    temp_growth(:)=t2m_month(:)-tp_00 

      
  END SUBROUTINE stomate_initialize
  

!! ================================================================================================================================
!! SUBROUTINE 	: stomate_main
!!
!>\BRIEF        Manages variable initialisation, reading and writing forcing 
!! files, aggregating data at stomate's time step (dt_stomate), aggregating data
!! at longer time scale (i.e. for phenology) and uses these forcing to calculate
!! CO2 fluxes (NPP and respirations) and C-pools (litter, soil, biomass, ...)
!!
!! DESCRIPTION  : The subroutine manages 
!! divers tasks:
!! (1) Initializing all variables of stomate (first call)
!! (2) Reading and writing forcing data (last call)
!! (3) Adding CO2 fluxes to the IPCC history files
!! (4) Converting the time steps of variables to maintain consistency between
!! sechiba and stomate
!! (5) Use these variables to call stomate_lpj, maint_respiration, littercalc,
!! som. The called subroutines handle: climate constraints 
!! for PFTs, PFT dynamics, Phenology, Allocation, NPP (based on GPP and
!! authothropic respiration), fire, mortality, vmax, assimilation temperatures,
!! all turnover processes, light competition, sapling establishment, lai,  
!! land cover change and litter and soil dynamics.
!! (6) Use the spin-up method developed by Lardy (2011)(only if SPINUP_ANALYTIC 
!! is set to TRUE).
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S): deadleaf_cover, assim_param, lai, height, veget, 
!! veget_max, resp_maint, 
!! resp_hetero,resp_growth, co2_flux, fco2_lu.
!!
!! REFERENCES	: 
!! - Lardy, R, et al., A new method to determine soil organic carbon equilibrium,
!! Environmental Modelling & Software (2011), doi:10.1016|j.envsoft.2011.05.016
!!
!! FLOWCHART    : 
!! \latexonly 
!! \includegraphics[scale=0.5]{stomatemainflow.png}
!! \endlatexonly
!! \n
!_ ================================================================================================================================
  
SUBROUTINE stomate_main &
       & (kjit, kjpij, kjpindex, &
       &  index, lalo, neighbours, resolution, contfrac, totfrac_nobio, clay, &
       &  silt, bulk, t2m, temp_sol, stempdiag, &
       &  humrel, shumdiag, litterhumdiag, precip_rain, precip_snow, &
       &  tmc_pft, drainage_pft, runoff_pft, swc_pft, &
       &  gpp, deadleaf_cover, assim_param, &
       &  max_eau_var, &
       &  lai, frac_age, height, veget, veget_max, &
       &  veget_max_new, woodharvest, totfrac_nobio_new, &
       &  hist_id, hist2_id, rest_id_stom, hist_id_stom, hist_id_stom_IPCC, &
       &  co2_flux, fco2_lu, resp_maint, resp_hetero, resp_growth, temp_growth, soil_pH, pb, n_input,p_input, EW_grain_add, &
!gmjc
       &  tmc_topgrass,fc_grazing, snow, agri_fert_save)
!end gmjc  
    
    IMPLICIT NONE

    
  !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER,INTENT(in)                       :: kjit              !! Time step number (unitless)
    INTEGER,INTENT(in)                       :: kjpindex          !! Domain size - terrestrial pixels only (unitless)
    INTEGER,INTENT(in)                       :: kjpij             !! Total size of the un-compressed grid (unitless)
    INTEGER,INTENT(in)                       :: rest_id_stom      !! STOMATE's _Restart_ file identifier (unitless)
    INTEGER,INTENT(in)                       :: hist_id_stom      !! STOMATE's _history_ file identifier (unitless)
    INTEGER,INTENT(in)                       :: hist_id_stom_IPCC !! STOMATE's IPCC _history_ file identifier 
                                                                         !! (unitless) 
    INTEGER,DIMENSION(kjpindex),INTENT(in)   :: index             !! Indices of the pixels on the map. Stomate uses a 
                                                                         !! reduced grid excluding oceans. ::index contains 
                                                                         !! the indices of the terrestrial pixels only 
                                                                         !! (unitless) 
    INTEGER,DIMENSION(kjpindex,NbNeighb),INTENT(in) :: neighbours !! Neighoring grid points if land for the DGVM 
                                                                         !! (unitless) 
    REAL,DIMENSION(kjpindex,2),INTENT(in)    :: lalo              !! Geographical coordinates (latitude,longitude) 
                                                                         !! for pixels (degrees) 
    REAL,DIMENSION(kjpindex,2),INTENT(in)    :: resolution        !! Size in x an y of the grid (m) - surface area of 
                                                                         !! the gridbox 
    REAL,DIMENSION (kjpindex), INTENT (in)   :: contfrac          !! Fraction of continent in the grid cell (unitless)
    REAL,DIMENSION(kjpindex),INTENT(in)      :: totfrac_nobio     !! Fraction of grid cell covered by lakes, land 
                                                                         !! ice, cities, ... (unitless) 
    REAL,DIMENSION(kjpindex),INTENT(in)      :: clay              !! Clay fraction of soil (0-1, unitless)
    REAL,DIMENSION(kjpindex),INTENT(in)      :: silt              !! Silt fraction of soil (0-1, unitless) 
    REAL,DIMENSION(kjpindex),INTENT(in)      :: bulk              !! Bulk density (kg/m**3) 
    REAL,DIMENSION(kjpindex,nvm),INTENT(in)  :: humrel            !! Relative humidity ("moisture availability") 
                                                                         !! (0-1, unitless) 
    REAL,DIMENSION(kjpindex),INTENT(in)      :: t2m               !! 2 m air temperature (K)
    REAL,DIMENSION(kjpindex),INTENT(in)      :: temp_sol          !! Surface temperature (K)
    REAL,DIMENSION(kjpindex,nslm),INTENT(in) :: stempdiag         !! Soil temperature (K)
    REAL,DIMENSION(kjpindex,nslm),INTENT(in) :: shumdiag          !! Relative soil moisture (0-1, unitless)
    REAL,DIMENSION(kjpindex),INTENT(in)      :: litterhumdiag     !! Litter humidity (0-1, unitless)
    REAL,DIMENSION(kjpindex),INTENT(in)      :: precip_rain       !! Rain precipitation  
                                                                         !! @tex $(mm dt_stomate^{-1})$ @endtex 
    REAL,DIMENSION(kjpindex),INTENT(in)      :: precip_snow       !! Snow precipitation  
                                                                         !! @tex $(mm dt_stomate^{-1})$ @endtex 
!============ DSG: the units are mm dt_sechiba-1
    REAL, DIMENSION (kjpindex,nvm), INTENT(in)   :: tmc_pft             !! Total soil water per PFT (mm/m2) 
    REAL, DIMENSION (kjpindex,nvm), INTENT(in)   :: drainage_pft        !! Drainage per PFT (mm/m2)   
    REAL, DIMENSION (kjpindex,nvm), INTENT(in)   :: runoff_pft          !! Runoff per PFT (mm/m2)   
    REAL, DIMENSION (kjpindex,nvm), INTENT(in)   :: swc_pft             !! Relative Soil water content per pft (-)     
!============ 

    REAL, DIMENSION (kjpindex), INTENT(in)   :: max_eau_var       !! Maximum water content of the soil @tex $(kg m^{-2})$ @endtex

    REAL,DIMENSION(kjpindex,nvm),INTENT(in)  :: gpp               !! GPP of total ground area  
                                                                         !! @tex $(gC m^{-2} time step^{-1})$ @endtex 
                                                                         !! Calculated in sechiba, account for vegetation 
                                                                         !! cover and effective time step to obtain ::gpp_d 
    REAL,DIMENSION(kjpindex,nvm),INTENT(in)  :: veget_max_new     !! New "maximal" coverage fraction of a PFT: only if 
                                                                         !! vegetation is updated in slowproc
    REAL,DIMENSION(kjpindex),INTENT(in)      :: woodharvest       !! Harvested wood biomass (gC m-2 yr-1)
    REAL,DIMENSION(kjpindex),INTENT(in)      :: totfrac_nobio_new !! New fraction of nobio per gridcell

    REAL,DIMENSION(kjpindex),INTENT(in)       :: soil_pH          !! soil pH 	 
    REAL,DIMENSION(kjpindex), INTENT(in)      :: pb               !! Air pressure (hPa) 

    INTEGER,INTENT(in)                       :: hist_id           !! ?? [DISPENSABLE] SECHIBA's _history_ file 
                                                                         !! identifier 
    INTEGER,INTENT(in)                       :: hist2_id          !! ?? [DISPENSABLE] SECHIBA's _history_ file 2 
                                                                         !! identifier 
    REAL,DIMENSION(kjpindex,nvm,5),INTENT(in) :: agri_fert_save   !! saved fertilization maps of agricultural practices
    REAL,DIMENSION (kjpindex), INTENT(in)     :: tmc_topgrass     !! daily top 5 layer mean soil moisture for grazing
    REAL,DIMENSION (kjpindex), INTENT(in)     :: fc_grazing       !! soil moisture threshold for grazing
    REAL,DIMENSION(kjpindex), INTENT(in)      :: snow             !! Snow mass [Kg/m^2]

    !! 0.2 Output variables

    REAL,DIMENSION(kjpindex,nvm),INTENT(out) :: co2_flux          !! CO2 flux between atmosphere and biosphere per 
                                                                         !! average ground area 
                                                                         !! @tex $(gC m^{-2} dt_stomate^{-1})$ @endtex  
                                                                         !! [??CHECK] sign convention? 
    REAL,DIMENSION(kjpindex),INTENT(out)     :: fco2_lu           !! CO2 flux between atmosphere and biosphere from 
                                                                         !! land-use (without forest management)  
                                                                         !! @tex $(gC m^{-2} dt_stomate^{-1})$ @endtex  
                                                                         !! [??CHECK] sign convention? 
    REAL,DIMENSION(kjpindex,nvm),INTENT(out) :: resp_maint        !! Maitenance component of autotrophic respiration in 
                                                                         !! @tex $(gC m^{-2} dt_stomate^{-1})$ @endtex 
    REAL,DIMENSION(kjpindex,nvm),INTENT(out) :: resp_growth       !! Growth component of autotrophic respiration in 
                                                                         !! @tex ($gC m^{-2} dt_stomate^{-1}$) @endtex
    REAL,DIMENSION(kjpindex,nvm),INTENT(out) :: resp_hetero       !! Heterotrophic respiration in  
                                                                         !! @tex $(gC m^{-2} dt_stomate^{-1})$ @endtex  
    REAL,DIMENSION(kjpindex),INTENT(out)     :: temp_growth       !! Growth temperature (Â°C)  
                                                                         !! Is equal to t2m_month 

    !! 0.3 Modified
   
    REAL,DIMENSION(kjpindex,nvm),INTENT(inout)       :: lai            !! Leaf area inex @tex $(m^2 m^{-2})$ @endtex
    REAL,DIMENSION(kjpindex,nvm),INTENT(in)          :: veget          !! Fraction of vegetation type including 
                                                                              !! non-biological fraction (unitless) 
    REAL,DIMENSION(kjpindex,nvm),INTENT(inout)       :: veget_max      !! Maximum fraction of vegetation type including 
                                                                              !! non-biological fraction (unitless) 
    REAL,DIMENSION(kjpindex,nvm),INTENT(inout)       :: height         !! Height of vegetation (m)
    REAL,DIMENSION(kjpindex,nvm,npco2),INTENT(inout) :: assim_param    !! min+max+opt temperatures (K) & vmax, nue and leaf N for 
                                                                              !! photosynthesis  
                                                                              !! @tex $(\mu mol m^{-2}s^{-1})$ @endtex  
    REAL,DIMENSION(kjpindex),INTENT(inout)           :: deadleaf_cover !! Fraction of soil covered by dead leaves 
                                                                              !! (unitless) 
    REAL,DIMENSION(kjpindex,nvm,nleafages),INTENT(inout):: frac_age    !! Age efficacity from STOMATE 
    REAL,DIMENSION(kjpindex,nvm,ninput),INTENT(inout) :: n_input       !! Nitrogen inputs into the soil  (gN/m**2/day)

    REAL,DIMENSION(kjpindex,nvm,pinput),INTENT(inout) :: p_input       !! Phosphorus inputs into the soil  (gP/m**2/day) 

    REAL,DIMENSION(kjpindex,nvm,nminerals,2),INTENT(inout) :: EW_grain_add !! Enhanced Weathering addition: grain mass(1) [g m-2] and grain size(2) [mu m]
    !! 0.4 local variables
    
    REAL                                   :: dt_days_read             !! STOMATE time step read in restart file (days)
    INTEGER                                :: l,k,ji, jv, i, j, m      !! indices    
    REAL,PARAMETER                         :: max_dt_days = 5.         !! Maximum STOMATE time step (days)
    REAL                                   :: hist_days                !! Writing frequency for history file (days)
    REAL,DIMENSION(0:nslm)                 :: z_soil                   !! Variable to store depth of the different soil 
                                                                              !! layers (m) 
    REAL,DIMENSION(kjpindex,nvm)           :: rprof                    !! Coefficient of the exponential functions that 
                                                                              !! relates root density to soil depth (unitless) 
    REAL,DIMENSION(kjpindex)               :: cvegtot                  !! Total "vegetation" cover (unitless)
    REAL,DIMENSION(kjpindex)               :: precip                   !! Total liquid and solid precipitation  
                                                                              !! @tex $(??mm dt_stomate^{-1})$ @endtex 
    REAL,DIMENSION(kjpindex,nvm)           :: gpp_d                    !! Gross primary productivity per ground area 
                                                                              !! @tex $(??gC m^{-2} dt_stomate^{-1})$ @endtex  
    REAL,DIMENSION(kjpindex,nvm)           :: gpp_daily_x              !! "Daily" gpp for teststomate  
                                                                              !! @tex $(??gC m^{-2} dt_stomate^{-1})$ @endtex 
    REAL,DIMENSION(kjpindex,nvm)           :: resp_hetero_litter       !! Litter heterotrophic respiration per ground area 
                                                                              !! @tex $(gC m^{-2} day^{-1})$ @endtex  
                                                                              !! ??Same variable is also used to 
                                                                              !! store heterotrophic respiration per ground area 
                                                                              !! over ::dt_sechiba?? 
    REAL,DIMENSION(kjpindex,nvm)           :: resp_hetero_soil         !! soil heterotrophic respiration  
                                                                              !! @tex $(gC m^{-2} day^{-1})$ @endtex
    REAL,DIMENSION(kjpindex,nvm)           :: veget_cov                !! Fractional coverage: actually share of the pixel 
                                                                              !! covered by a PFT (fraction of ground area), 
                                                                              !! taking into account LAI ??(= grid scale fpc)?? 
    REAL,DIMENSION(kjpindex,nvm)           :: veget_cov_max_new        !! New value for maximal fractional coverage (unitless)
    REAL,DIMENSION(kjpindex,nvm,2)         :: vcmax                    !! Maximum rate of carboxylation
                                                                              !! third dimension: actual Vcmax=1 ; potential Vcmax=2
    REAL,DIMENSION(kjpindex,nvm,2)         :: jmax                     !! maximum rate of electron transport
                                                                              !! third dimension: actual Jmax=1 ; potential Jmax=2
    REAL,DIMENSION(kjpindex,nvm)           :: nue                      !! Nitrogen use Efficiency with impact of leaf age (umol CO2 (gN)-1 s-1)                                                                               !! @tex $(\mumol m^{-2} s^{-1})$ @endtex
    REAL,DIMENSION(kjpindex,nlevs)         :: control_moist_inst       !! Moisture control of heterotrophic respiration 
                                                                              !! (0-1, unitless) 
    REAL,DIMENSION(kjpindex,nlevs)         :: control_temp_inst        !! Temperature control of heterotrophic 
                                                                              !! respiration, above and below (0-1, unitless) 
    REAL,DIMENSION(kjpindex,ncarb,nvm,nelements)     :: som_input_inst    !! Quantity of carbon going into carbon pools from 
                                                                              !! litter decomposition 
                                                                              !! @tex $(gC m^{-2} day^{-1})$ @endtex 
    
    REAL,DIMENSION(kjpindex)               :: stemp30cm                   !! Average Soil temperature of first 8 layers [K]

    INTEGER                                :: ier                      !! Check errors in netcdf call (unitless)
    REAL                                   :: sf_time                  !! Intermediate variable to calculate current time 
                                                                              !! step 
    INTEGER                                :: max_totsize              !! Memory management - maximum memory size (Mb)
    INTEGER                                :: totsize_1step            !! Memory management - memory required to store one 
                                                                              !! time step on one processor (Mb) 
    INTEGER                                :: totsize_tmp              !! Memory management - memory required to store one 
                                                                              !! time step on all processors(Mb) 
    REAL                                   :: xn                       !! How many times have we treated in this forcing 
    REAL, DIMENSION(kjpindex)              :: vartmp                   !! Temporary variable
    INTEGER                                :: vid                      !! Variable identifer of netCDF (unitless)
    INTEGER                                :: nneigh                   !! Number of neighbouring pixels
    INTEGER                                :: direct                   !! ??
    INTEGER,DIMENSION(ndm)                 :: d_id                     !! ??
    REAL                                   :: net_nep_monthly          !! Integrated nep_monthly over all grid-cells on local domain
    REAL                                   :: net_nep_monthly_sum      !! Integrated nep_monthly over all grid-cells on total domain(global)
    REAL                                   :: net_cflux_prod_monthly_sum    !! AR5 output?? gC m2 month-1 (one variable for 
                                                                                   !! reduce_sum and one for bcast??), parallel 
                                                                                   !! computing 
    REAL                                   :: net_cflux_prod_monthly_tot    !! AR5 output?? gC m2 month-1 (one variable for 
                                                                                   !! reduce_sum and one for bcast??), parallel 
                                                                                   !! computing 
    REAL                                   :: net_harvest_above_monthly_sum !! AR5 output?? gC m2 month-1 (one variable for 
                                                                                   !! reduce_sum and one for bcast??), parallel 
                                                                                   !! computing 
    REAL                                   :: net_harvest_above_monthly_tot !! AR5 output?? gC m2 month-1 (one variable for 
                                                                                   !! reduce_sum and one for bcast??), parallel 
                                                                                   !! computing 
    REAL                                   :: net_biosp_prod_monthly_sum    !! AR5 output?? gC m2 month-1 (one variable for 
                                                                                   !! reduce_sum and one for bcast??), parallel 
                                                                                   !! computing 
    REAL                                   :: net_biosp_prod_monthly_tot    !! AR5 output?? gC m2 month-1 (one variable for 
                                                                                   !! reduce_sum and one for bcast??), parallel 
                                                                                   !! computing 
    REAL, DIMENSION(kjpindex,nvm,nbpools)  :: carbon_stock                  !! Array containing the carbon stock for each pool
                                                                                   !! used by ORCHIDEE
   REAL, DIMENSION(kjpindex,nvm,nionspec) :: n_uptake                       !! Uptake of soil N by plants  
   !! (gN/m**2/timestep)  
   REAL, DIMENSION(kjpindex,nvm)          :: n_mineralisation               !! net nitrogen mineralisation of decomposing SOM 
   !!   (gN/m**2/day), supposed to be NH4 
    REAL, DIMENSION(kjpindex,nvm)          :: p_uptake                      !! Uptake of soil P by plants 
                                                                                   !! (gP/m**2/timestep) 
    REAL, DIMENSION(kjpindex,nvm)          :: p_leaching                    !! Leaching of dissolved P 
                                                                                   !! (gP/m**2/timestep) 
    REAL, DIMENSION(kjpindex,nvm)          :: p_ssorption                   !! Strongly sorption of sorbed P
                                                                                   !! (gP/m**2/timestep) 

    REAL, DIMENSION(kjpindex,nvm)          :: p_mineralisation              !! net phosphorus mineralisation/immobilization of decomposing SOM
                                                                                   !!   (gP/m**2/day)

   REAL,DIMENSION(kjpindex,nvm)       :: resp_total_soil                    !! soil heterotrophic respiration (gC/day/m**2 by PFT) 
   REAL, DIMENSION(kjpindex,nvm,ncarb)     :: CN_target                      !! C to N ratio of SOM flux from one pool to another (gN m-2 dt-1) 
    REAL, DIMENSION(kjpindex,nvm,ncarb)    :: CP_target                 !! N to P ratio of SOM flux from one pool to another (g(N) g-1(P))

!EXTERNALIZE: it is also used in stomate_resp.f90
    REAL, PARAMETER                  :: SNC=0.004            !! Structural nitrogen   [gN g-1C] based on C:N of dead wood (White et al., 2000)
    REAL, PARAMETER                  :: SPC=0.00036          !! Structural phosphorus [gN g-1C] based on C:N of dead wood (White et al., 2000)
                                                                    !!                       and N:P of wood from Sardans et al (2015) ~11
                                                                    !!                       and N:P of dead wood from Weedon et al (2009) ~11
!EXTERNALIZE
    REAL , DIMENSION(kjpindex,nvm,nparts)          :: cn_test
    REAL , DIMENSION(kjpindex,nvm,nparts)          :: np_test

    REAL , DIMENSION(kjpindex,nvm,nelements)          :: mass_before
    REAL , DIMENSION(kjpindex,nvm,nelements)          :: mass_before_2
    REAL , DIMENSION(kjpindex,nvm,nelements)          :: mass_after
    REAL , DIMENSION(kjpindex,nvm,nelements)          :: mass_change
    REAL , DIMENSION(kjpindex,nvm,nelements)          :: flux_before

!    REAL , DIMENSION(kjpindex,nvm)          :: totPchange

    REAL, DIMENSION (kjpindex,nvm)                 :: tmpvar          !! Runoff per PFT (mm/m2)   

    REAL, DIMENSION(kjpindex)          :: sand_arg     ! clay fraction (between 0 and 1)
    REAL, DIMENSION(kjpindex)          :: tempsol_arg ! soil temperature (degC)
    REAL, DIMENSION(kjpindex,nvm)      :: drainagepft_arg


!_ ================================================================================================================================
    
    !! 1. Initialize variables

    sand_arg        = zero
    tempsol_arg     = zero
    drainagepft_arg = zero

    !! 1.1 Store current time step in a common variable
    itime = kjit
    
    !![DISPENSABLE] 1.2 Copy the depth of the different soil layers from diaglev specified in slow_proc
    !! 1.3 PFT rooting depth across pixels, humescte is pre-defined 
    ! (constantes_veg.f90). It is defined as the coefficient of an exponential 
    ! function relating root density to depth 
    DO j=1,nvm
       rprof(:,j) = 1./humcste(j)
    ENDDO
    
    !! 1.4 Initialize first call
    ! Set growth respiration to zero
    resp_growth=zero

    ! Check that initialization is done
    IF (l_first_stomate) CALL ipslerr_p(3,'stomate_main','Initialization not yet done.','','')
    
    IF (printlev >= 4) THEN
       WRITE(numout,*) 'stomate_main: date=',days_since_beg,' ymds=', year_end, month_end, day_end, sec_end, &
            ' itime=', itime, ' do_slow=',do_slow
       WRITE(numout,*) 'latlon of test_grid:',lalo(test_grid,:)
    ENDIF

    IF (printlev>=3) WRITE (numout,*) 'BEGIN stomate_main soil_n_min(test_grid,test_pft,:):',soil_n_min(test_grid,test_pft,:)

!! 3. Special treatment for some input arrays.
    
    !! 3.1 Sum of liquid and solid precipitation
    precip(:) = ( precip_rain(:) + precip_snow(:) )*one_day/dt_sechiba
    
    !! 3.2 Calculate STOMATE's vegetation fractions from veget and veget_max
    DO j=1,nvm 
       WHERE ((1.-totfrac_nobio(:)) > min_sechiba)
          ! Pixels with vegetation
          veget_cov(:,j) = veget(:,j)/( 1.-totfrac_nobio(:) )
          veget_cov_max(:,j) = veget_max(:,j)/( 1.-totfrac_nobio(:) )
       ELSEWHERE
          ! Pixels without vegetation
          veget_cov(:,j) = zero
          veget_cov_max(:,j) = zero
       ENDWHERE
    ENDDO

    IF ( do_now_stomate_lcchange ) THEN
       DO j=1,nvm
          WHERE ((1.-totfrac_nobio_new(:)) > min_sechiba)
             ! Pixels with vegetation
             veget_cov_max_new(:,j) = veget_max_new(:,j)/( 1.-totfrac_nobio_new(:) )
          ELSEWHERE
             ! Pixels without vegetation
             veget_cov_max_new(:,j) = zero
          ENDWHERE
       ENDDO
    ENDIF

    !! 3.3 Adjust time step of GPP 
    ! No GPP for bare soil
    gpp_d(:,1) = zero
    ! GPP per PFT
    DO j = 2,nvm   
       WHERE (veget_cov_max(:,j) > min_stomate)
          ! The PFT is available on the pixel
          gpp_d(:,j) =  gpp(:,j)/ veget_cov_max(:,j)* one_day/dt_sechiba  
       ELSEWHERE
          ! The PFT is absent on the pixel
          gpp_d(:,j) = zero
       ENDWHERE
    ENDDO

  !! 4. Calculate variables for dt_stomate (i.e. "daily")

    ! Note: If dt_days /= 1, then variables 'xx_daily' (eg. half-daily or bi-daily) are by definition
    ! not expressed on a daily basis. This is not a problem but could be
    ! confusing

    !! 4.1 Accumulate instantaneous variables (do_slow=.FALSE.) 
    ! Accumulate instantaneous variables (do_slow=.FALSE.) and eventually 
    ! calculate daily mean value (do_slow=.TRUE.) 
    CALL stomate_accu (do_slow, humrel,        humrel_daily)
    CALL stomate_accu (do_slow, max_eau_var,   max_eau_var_daily)
    CALL stomate_accu (do_slow, litterhumdiag, litterhum_daily)
    CALL stomate_accu (do_slow, t2m,           t2m_daily)
    CALL stomate_accu (do_slow, temp_sol,      tsurf_daily)
    CALL stomate_accu (do_slow, stempdiag,     tsoil_daily)
    CALL stomate_accu (do_slow, shumdiag,      soilhum_daily)
    CALL stomate_accu (do_slow, precip,        precip_daily)
    CALL stomate_accu (do_slow, gpp_d,         gpp_daily)
    CALL stomate_accu (do_slow, drainage_pft*one_day/dt_sechiba, drainage_daily)
    CALL stomate_accu (do_slow, runoff_pft*one_day/dt_sechiba,   runoff_daily) 

    CALL stomate_accu (do_slow, precip_snow, snowfall_daily)
    CALL stomate_accu (do_slow, snow,        snowmass_daily)
    CALL stomate_accu (do_slow, tmc_topgrass, tmc_topgrass_daily)

    !! 4.2 Daily minimum temperature
    t2m_min_daily(:) = MIN( t2m(:), t2m_min_daily(:) )

    !! 4.3 Calculate maintenance respiration
    ! Note: lai is passed as output argument to overcome previous problems with 
    ! natural and agricultural vegetation types. 
    CALL maint_respiration &
         & (kjpindex,t2m,stempdiag,  &
         & rprof,biomass,resp_maint_part_radia)
    
    ! Aggregate maintenance respiration across the different plant parts 
    resp_maint_radia(:,:) = zero
    DO j=2,nvm
       DO k= 1, nparts
          resp_maint_radia(:,j) = resp_maint_radia(:,j) &
               & + resp_maint_part_radia(:,j,k)
       ENDDO
    ENDDO
    
    ! Maintenance respiration separated by plant parts
    resp_maint_part(:,:,:) = resp_maint_part(:,:,:) &
         & + resp_maint_part_radia(:,:,:)
    
    !! 4.4 Litter dynamics and litter heterothropic respiration 
    ! Including: litter update, lignin content, PFT parts, litter decay,
    ! litter heterotrophic respiration, dead leaf soil cover.
    ! Note: there is no vertical discretisation in the soil for litter decay.
    turnover_littercalc(:,:,:,:) = turnover_daily(:,:,:,:) * dt_sechiba/one_day
    bm_to_littercalc   (:,:,:,:) = bm_to_litter  (:,:,:,:) * dt_sechiba/one_day   

    n_mineralisation(:,:) = zero 
    p_mineralisation(:,:) = zero 


    IF(dsg_debug) THEN 
        ! Record the total P mass in the beginning of stomate iteration:
        mass_before_2(:,:,:) = SUM(som(:,:,:,:),DIM=2) +  SUM(SUM(litter(:,:,:,:,:),DIM=4),DIM=2) + SUM(biomass(:,:,:,:),DIM=3) 
        ! need to correct for the fluxes not yet added to litter but already removed from vegetation:
        mass_before_2(:,:,:) = mass_before_2(:,:,:) + SUM(turnover_daily(:,:,:,:),DIM=3) + SUM(bm_to_litter(:,:,:,:),DIM=3)
        mass_before_2(:,:,iphosphorus) = mass_before_2(:,:,iphosphorus) + SUM(soil_p_min(:,:,:),DIM=3)
        mass_before_2(:,:,initrogen) = mass_before_2(:,:,initrogen) + SUM(soil_n_min(:,:,:),DIM=3)
        WRITE (numout,*) 'BEGIN stomate_main total P=', mass_before_2(test_grid,test_pft,iphosphorus)
    ENDIF

    ! TEST P BALANCE

    
    ! === DSG mass conservation ======================
    mass_before(:,:,:) =SUM(som(:,:,:,:),DIM=2) +  SUM(SUM(litter(:,:,:,:,:),DIM=4),DIM=2) 

    CALL littercalc (kjpindex, &
         turnover_littercalc, bm_to_littercalc, &
         veget_cov_max, temp_sol, stempdiag, shumdiag, litterhumdiag, &
         soil_n_min, litter, &
         litter_avail, litter_not_avail, litter_avail_frac, &
         dead_leaves, lignin_struc, &
         lignin_wood, n_mineralisation, p_mineralisation, deadleaf_cover, resp_hetero_litter, &
         som_input_inst, control_temp_inst, control_moist_inst, &
         matrixA, vectorB, CN_target,CN_som_litter_longterm,tau_CN_longterm, &
         CP_target,CP_som_litter_longterm, &
         sla_calc,do_slow)

    ! Heterothropic litter respiration during time step ::dt_sechiba @tex $(gC m^{-2})$ @endtex
    resp_hetero_litter(:,:) = resp_hetero_litter(:,:) * dt_sechiba/one_day
    
    !! 4.5 Soil carbon dynamics and soil heterotrophic respiration
    ! Note: there is no vertical discretisation in the soil for litter decay.

    CALL som_dynamics (kjpindex, clay, silt, &
         som_input_inst, control_temp_inst, control_moist_inst, &
         CN_target, CP_target, som, resp_hetero_soil, matrixA, &
         n_mineralisation, p_mineralisation, &
         CN_som_litter_longterm, CP_som_litter_longterm, tau_CN_longterm)

    IF(dsg_debug) THEN 
       ! === DSG mass conservation ======================
       !WRITE(6,*) 'location in code: '
       !WRITE(6,*) 'END of som_dynamics (including litter)'
       
       ! Disable call as its in a loop and thus inflating the out stream
       !CALL    check_mass(kjpindex,biomass(:,:,:,:), 'END of som_dynamics(including litter)')
       ! =================================================
       mass_change(:,:,:)      = SUM(bm_to_littercalc(:,:,:,:),DIM=3)   
       mass_change(:,:,:)      = mass_change(:,:,:) + SUM(turnover_littercalc(:,:,:,:),DIM=3)  
       !DSG: this is fucking insane !!!! resp_hetero_soil and resp_hetero_litter
       !have different timesteps here !!!!
       mass_change(:,:,1)      = mass_change(:,:,1) - (resp_hetero_soil(:,:)*dt_sechiba/one_day + resp_hetero_litter(:,:) ) 
       mass_change(:,:,2)      = mass_change(:,:,2) - n_mineralisation(:,:) 
       mass_change(:,:,3)      = mass_change(:,:,3) - p_mineralisation(:,:) 
   
       mass_after(:,:,:)  = SUM(som(:,:,:,:),DIM=2) +  SUM(SUM(litter(:,:,:,:,:),DIM=4),DIM=2)
   
       CALL cons_mass(mass_before(:,:,:),               &  ! mass before
             mass_after(:,:,:),                         &  ! mass after
             mass_change(:,:,:),                        &  ! net of fluxes
             'END of som_dynamics(including litter)') 
    ENDIF
    
    ! Heterothropic soil respiration during time step ::dt_sechiba @tex $(gC m^{-2})$ @endtex 
    resp_hetero_soil(:,:) = resp_hetero_soil(:,:) * dt_sechiba/one_day

    ! Total heterothrophic respiration during time step ::dt_sechiba @tex $(gC m^{-2})$ @endtex
    resp_hetero_radia(:,:) = resp_hetero_litter(:,:) + resp_hetero_soil(:,:)
    resp_hetero_d(:,:) = resp_hetero_d(:,:) + resp_hetero_radia(:,:)

    
    ! Sum heterotrophic and autotrophic respiration in soil 
    resp_total_soil(:,:) = resp_hetero_radia(:,:) + & 
         resp_maint_part_radia(:,:,isapbelow) + resp_maint_part_radia(:,:,iroot)
    
    IF (printlev>=3) WRITE (numout,*) '4.5'
    IF (printlev>=3) WRITE (numout,*) 'resp_hetero_litter(test_grid,test_pft):',resp_hetero_litter(test_grid,test_pft)
    IF (printlev>=3) WRITE (numout,*) 'resp_hetero_soil(test_grid,test_pft):',resp_hetero_soil(test_grid,test_pft)
    IF (printlev>=3) WRITE (numout,*) 'resp_maint_part_radia(test_grid,test_pft,isapbelow):',resp_maint_part_radia(test_grid,test_pft,isapbelow)
    IF (printlev>=3) WRITE (numout,*) 'resp_maint_part_radia(test_grid,test_pft,iroot):',resp_maint_part_radia(test_grid,test_pft,iroot)


    IF(ok_ncycle) THEN
       
       n_uptake(:,:,:) = zero

       ! tempereature for root uptake and min N processes:
       ! OCN assumes N dynamics happen first 30cm of soil;
       ! first 8 layers correspond to first ~37cm of soil
       ! so we use the average temperature of the first 8 layers:
      
       !stemp30cm(:) = (stempdiag(:,1:8)) ! * dlt(1:8))/SUM(dlt(1:8))
       stemp30cm(:) = SUM(stempdiag(:,1:8) * SPREAD(dlt(1:8),NCOPIES=kjpindex,DIM=1),DIM=2) &
                      & /SPREAD(SUM(dlt(1:8),DIM=1),NCOPIES=kjpindex,DIM=1)
        
       sand_arg = MAX(zero, un - silt - clay)
       tempsol_arg = stemp30cm-tp_00
       drainagepft_arg = (drainage_pft+runoff_pft)

       CALL nitrogen_dynamics(kjpindex, contfrac, clay, sand_arg, & 
!DSG: hydrology layering: not an issue
!DSG: temperature layering: we take average of first 30cm
!            stempdiag-tp_00, tmc_pft, (drainage_pft+runoff_pft),swc_pft, veget_cov_max, resp_total_soil, som, &  
            tempsol_arg, tmc_pft, drainagepft_arg, swc_pft, veget_cov_max, resp_total_soil, som, &  
            
!DSG: layering
            biomass, n_input, soil_ph,           &
            npp_week,                            &
            n_mineralisation, pb,                &  
            max_eau_var,                         &
            n_uptake, bulk,soil_n_min,p_O2,bact, &
            grm_nfert, agri_nfert, fclover, BNF_clover) 

        n_uptake_daily(:,:,:)   = n_uptake_daily(:,:,:)   + n_uptake(:,:,:)
        n_mineralisation_d(:,:) = n_mineralisation_d(:,:) + n_mineralisation(:,:)
    ELSE
        soil_n_min(:,:,:)       = un
        n_uptake_daily(:,:,:)   = zero
        n_uptake(:,:,:)         = zero
        n_mineralisation_d(:,:) = zero
        n_mineralisation(:,:)   = zero
        n_uptake_daily(:,:,:)   = n_uptake_daily(:,:,:)   + n_uptake(:,:,:)  
        n_mineralisation_d(:,:) = n_mineralisation_d(:,:) + n_mineralisation(:,:)  
    ENDIF

    ! JCADD 20170427 BNF for grassland
    BNF_clover_daily(:,:) = BNF_clover_daily(:,:) + BNF_clover(:,:)
    ! END JCADD
    IF(ok_pcycle )THEN
       ! We need to calculate the fast variable p_release every time step
      
       p_input_dt(:,:) = SUM(p_input(:,:,:),DIM=3)*dt_sechiba/one_day

       tmpvar = drainage_pft+runoff_pft
       CALL phosphorus_dynamics(kjpindex,dt_sechiba/one_day,                          &
                                  soil_orders,                                        &
                                  clay,silt, bulk,                                    &
                                  tmpvar ,tmc_pft,swc_pft,                            & 
                                  max_eau_var,                                        &
                                  veget_max,                                          &
!DSG: temperature layering: we take average of first 30cm
!                                  stempdiag,                                          &
                                  stemp30cm,                                          &
!DSG: layering
                                  p_input_dt,                                         &
                                  p_mineralisation,                                   &
                                  soil_p_min,                                         &
                                  f_Pdissolved,                                       &
                                  som,                                                &
                                  p_uptake, p_ssorption, p_leaching,                  &
                                  biomass,                                            &
                                  grm_pfert,agri_pfert)

       p_uptake_daily(:,:)     = p_uptake_daily(:,:)    + p_uptake(:,:) 
       p_ssorption_daily(:,:)  = p_ssorption_daily(:,:) + p_ssorption(:,:) 
       p_leaching_daily(:,:)   = p_leaching_daily(:,:)  + p_leaching(:,:) 
       p_mineralisation_d(:,:) = p_mineralisation_d(:,:) + p_mineralisation(:,:) 


       ! write daily value?
       CALL xios_orchidee_send_field('P_DEPOSITION',p_input(:,:,iatm_p)) 
       CALL xios_orchidee_send_field("P_RELEASE"   ,p_input(:,:,iweat_p))

!       CALL histwrite_p (hist_id_stomate, 'P_DEPOSITION', itime, &
!            p_input(:,:,iatm_p), kjpindex, hori_index) 


!       CALL histwrite_p (hist_id_stomate, 'P_RELEASE', itime, &
!            p_input(:,:,iweat_p), kjpindex, hori_index)

       CALL xios_orchidee_send_field('SOIL_ORDERS',REAL(soil_orders(:)))
       CALL xios_orchidee_send_field('LITH_FRAC',lith_frac(:,12))
       CALL xios_orchidee_send_field('SOIL_SHIELD',soil_shield(:))
    ELSE ! no P cycle

       soil_p_min(:,:,:)            = un  
       f_Pdissolved(:,:)            = un
       p_uptake_daily(:,:)          = zero
       p_leaching_daily(:,:)        = zero
       p_ssorption_daily(:,:)       = zero
       p_mineralisation_d(:,:)      = zero

    ENDIF
    !DSG: mass conservation check is done in the phosphorus_dynamics routine


    !! 4.6 Accumulate instantaneous variables (do_slow=.FALSE.) 
    ! Accumulate instantaneous variables (do_slow=.FALSE.) and eventually 
    ! calculate daily mean value (do_slow=.TRUE.) 
    CALL stomate_accu (do_slow, control_moist_inst, control_moist_daily)
    CALL stomate_accu (do_slow, control_temp_inst,  control_temp_daily)
    CALL stomate_accu (do_slow, som_input_inst, som_input_daily)

!! 5. Daily processes - performed at the end of the day
    IF (do_slow) THEN

    !   IF(dsg_debug) THEN
    !       ! TEST P BALANCE
    !       mass_after(:,:,iphosphorus) =SUM(som(:,:,:,iphosphorus),DIM=2) +  SUM(SUM(litter(:,:,:,:,iphosphorus),DIM=4),DIM=2) + SUM(biomass(:,:,:,iphosphorus),DIM=3) + SUM(soil_p_min(:,:,:),DIM=3)
    !       ! in the meanwhile we have added all the litter from yesterday, but removed plant uptkae  
    !       mass_after(:,:,iphosphorus) = mass_before_2(:,:,iphosphorus)  + p_uptake_daily(:,:)
    !       WRITE (numout,*) 'AFTER fast processes: total P              =', mass_after(test_grid,test_pft,iphosphorus)
    !       WRITE (numout,*) 'AFTER fast processes: change in P stocks   =', mass_after(test_grid,test_pft,iphosphorus) - mass_before_2(test_grid,test_pft,iphosphorus)
    !   ENDIF

       !! 5.1 Update lai
       ! Use lai from stomate
       ! ?? check if this is the only time ok_pheno is used??
       ! ?? Looks like it is the only time. But this variables probably is defined 
       ! in stomate_constants or something, in which case, it is difficult to track.
       IF (ok_pheno) THEN
          !! 5.1.1 Update LAI 
          ! Set lai of bare soil to zero
          lai(:,ibare_sechiba) = zero
          ! lai for all PFTs
          DO j = 2, nvm
             lai(:,j) = biomass(:,j,ileaf,icarbon)*sla(j)
          ENDDO
          frac_age(:,:,:) = leaf_frac(:,:,ileaf,:)
       ELSE 
          ! 5.1.2 Use a prescribed lai
          ! WARNING: code in setlai is identical to the lines above
          ! Update subroutine if LAI has to be forced 
          CALL  setlai(kjpindex,lai) 
          frac_age(:,:,:) = zero
       ENDIF

       !! 5.2 Calculate long-term "meteorological" and biological parameters
       ! mainly in support of calculating phenology. If LastTsYear=.TRUE.
       ! annual values are update (i.e. xx_lastyear).
       CALL season &
            &          (kjpindex, dt_days, &
            &           veget_cov, veget_cov_max, &
            &           humrel_daily, t2m_daily, tsoil_daily, soilhum_daily, lalo, &
            &           precip_daily, npp_daily, biomass, &
            &           turnover_daily, gpp_daily, when_growthinit, &
            &           maxhumrel_lastyear, maxhumrel_thisyear, &
            &           minhumrel_lastyear, minhumrel_thisyear, &
            &           maxgppweek_lastyear, maxgppweek_thisyear, &
            &           gdd0_lastyear, gdd0_thisyear, &
            &           precip_lastyear, precip_thisyear, &
            &           lm_lastyearmax, lm_thisyearmax, &
            &           maxfpc_lastyear, maxfpc_thisyear, &
            &           humrel_month, humrel_week, t2m_longterm, tau_longterm, &
            &           t2m_month, t2m_week, tsoil_month, soilhum_month, &
            &           drainage_daily, runoff_daily, &
            &           drainage_longterm, runoff_longterm, &
            &           npp_longterm, turnover_longterm, gpp_week, npp_week, &
            &           gdd_m5_dormance, gdd_midwinter, ncd_dormance, ngd_minus5, &
            &           time_hum_min, hum_min_dormance, gdd_init_date, &
            &           gdd_from_growthinit, herbivores, &
            &           Tseason, Tseason_length, Tseason_tmp, &
            &           Tmin_spring_time, t2m_min_daily, begin_leaves, onset_date, &
            &           cn_leaf_avg_season, np_leaf_avg_season, nstress_season, &
            &           pstress_season, & 
            &           lai_target, lai_target_longterm, &
            &           moiavail_growingseason,rue_longterm, &
            &           t2m_14, sla_calc)
        ! DSG: this routine (stomate_season) doesn't manipulate mass, so no mass conservation
        ! check needed here

       !! 5.DSG - update P release 
        IF( ok_pcycle )THEN
           IF (.NOT. read_pweat) THEN
               CALL phosphorus_weathering(kjpindex,                                          &
                                         lith_frac, soil_shield,                            &
                                         drainage_longterm, runoff_longterm, t2m_longterm,  & ! in goes annual hydrol flux (long term average)
                                         p_release_daily)                                     ! out: daily flux DSG: should be made consisten
               IF(enhanced_weat) THEN
                   CALL enhanced_weathering(kjpindex, dt_days, &
                                            temp_sol,soil_pH,  &
                                            EW_grain, EW_flux)

                  !DSG_remove IF (printlev>=1) THEN
                  !DSG_remove      WRITE (numout,*) 'EW:'
                  !DSG_remove      WRITE (numout,*) 'MAXVAL(EW_flux(:,:,ico2cons))',  MAXVAL(EW_flux(:,:,ico2cons))
                  !DSG_remove      WRITE (numout,*) 'MINVAL(EW_flux(:,:,ico2cons))',  MINVAL(EW_flux(:,:,ico2cons))
                  !DSG_remove      WRITE (numout,*) 'MAXVAL(EW_flux(:,:,iPrelease))',  MAXVAL(EW_flux(:,:,iPrelease))
                  !DSG_remove      WRITE (numout,*) 'MINVAL(EW_flux(:,:,iPrelease))',  MINVAL(EW_flux(:,:,iPrelease))
                  !DSG_remove ENDIF

                  IF (LastTsYear) THEN
                  ! DSG: we add the EW application at the end of the year
                      !EW_forcStock=.TRUE.
                      !EW_targetStock = 3.kg
                      !EW_targetSize = 20 um  

!<<<<<<< .mine
                      !DSG_remove IF (printlev>=1) THEN
                      !DSG_remove     WRITE (numout,*) 'Before:'
                      !DSG_remove     WRITE (numout,*) 'EW_target_stock:',EW_targetStock
                      !DSG_remove     WRITE (numout,*) 'MAXVAL(EW_grain(:,:,:,1))',  MAXVAL(EW_grain(:,:,:,1))
                      !DSG_remove     WRITE (numout,*) 'MINVAL(EW_grain(:,:,:,1))',  MINVAL(EW_grain(:,:,:,1))
                      !DSG_remove     WRITE (numout,*) 'EW_target_size:' ,EW_targetSize
                      !DSG_remove     WRITE (numout,*) 'MAXVAL(EW_grain(:,:,:,2))',  MAXVAL(EW_grain(:,:,:,2))
                      !DSG_remove ENDIF

                      IF (EW_force_stock)  THEN
                      IF (printlev>=1) WRITE (numout,*) 'Last TimeStep: we add rock powder'

                         ! We reapply every year to keep a certain soil stock
                         DO i= 1,nminerals
                             WHERE (EW_grain(:,:,i,1) .LT. EW_targetStock(i))
                                 ! record how much we must add:
                                 EW_grain_add(:,:,i,1) = EW_targetStock(i) - EW_grain(:,:,i,1)
                                 ! then set the stock back to the target value:
                                 EW_grain(:,:,i,1)     = EW_targetStock(i)
                                 
                                 EW_grain    (:,:,i,2) = EW_targetSize(i)
                                 EW_grain_add(:,:,i,2) = EW_targetSize(i)
                             ENDWHERE
                         ENDDO

                      ELSE ! default: we read the application rate from forcing file:
                         ! if no minerals have been added so far .OR. the size of added minerals
                         ! equals  the size of already applied minerals, we add new minerals  
                         WHERE ((EW_grain(:,:,:,1) .EQ. zero).OR.      &
                                (EW_grain(:,:,:,2) .EQ. undef).OR.   &
                                (EW_grain(:,:,:,2) .EQ. EW_grain_add(:,:,:,2)))
                             EW_grain(:,:,:,1) = EW_grain(:,:,:,1) + EW_grain_add(:,:,:,1)
                             EW_grain(:,:,:,2) = EW_grain_add(:,:,:,2)
                         ENDWHERE
                         IF (ANY(EW_grain(:,:,:,2) .NE. EW_grain_add(:,:,:,2))) THEN
                             CALL ipslerr_p (1,'stomate', &
                                  'Enhanced weathering: the size of the minerals to be added', &
                                  ' is different from size of already present minerals', &
                                  'MODEL STOP; check input files')
                         ENDIF
!=======
!                      IF (EW_forcStock)  THEN
!                         ! We reapply every year to keep a certain soil stock
!                         DO i= 1,nminerals
!                             WHERE (EW_grain(:,:,i,1) .LT. EW_targetStock(i))
!                                 ! record how much we must add:
!                                 EW_grain_add(:,:,i,1) = EW_targetStock(i) - EW_grain(:,:,i,1)
!                                 EW_grain_add(:,:,i,2) = EW_targetSize(i)
!
!                                 EW_grain(:,:,i,1) = EW_targetStock(i)
!                                 EW_grain(:,:,i,2) = EW_targetSize(i)
!                             ENDWHERE
!
!                         ENDDO
!                      ELSE
!                      ! DEFAULT MODE ( application rate read from forcing file:
!                         ! if no minerals have been added so far .OR. the size of added minerals
!                         ! equals the size of already applied minerals (we cannot mix different sizes), we add new minerals  
!                         WHERE ((EW_grain(:,:,:,1) .EQ. zero).OR.      &
!                                (EW_grain(:,:,:,2) .EQ. undef).OR.   &
!                                (EW_grain(:,:,:,2) .EQ. EW_grain_add(:,:,:,2)))
!                             EW_grain(:,:,:,1) = EW_grain(:,:,:,1) + EW_grain_add(:,:,:,1)
!                             EW_grain(:,:,:,2) = EW_grain_add(:,:,:,2)
!                         ENDWHERE
!                         IF (ANY(EW_grain(:,:,:,2) .NE. EW_grain_add(:,:,:,2))) THEN
!                             CALL ipslerr_p (1,'stomate', &
!                                  'ENhanced weathering: the size of the minerals to be added', &
!                                  ' is different from size of already present minerals', &
!                                  'MODEL STOP; check input files')
!                         ENDIF
!>>>>>>> .r5976
                      ENDIF
!<<<<<<< .mine

                      !DSG_remove IF (printlev>=1) THEN
                      !DSG_remove     WRITE (numout,*) 'After:'
                      !DSG_remove     WRITE (numout,*) 'EW_target_stock:',EW_targetStock
                      !DSG_remove     WRITE (numout,*) 'MAXVAL(EW_grain(:,:,:,1))',  MAXVAL(EW_grain(:,:,:,1))
                      !DSG_remove     WRITE (numout,*) 'MINVAL(EW_grain(:,:,:,1))',  MINVAL(EW_grain(:,:,:,1))
                      !DSG_remove     WRITE (numout,*) 'EW_target_size:' ,EW_targetSize
                      !DSG_remove     WRITE (numout,*) 'MAXVAL(EW_grain(:,:,:,2))',  MAXVAL(EW_grain(:,:,:,2))
                      !DSG_remove ENDIF

                  ! send EW_grain_add(:,:,:,1) to xios to record how much we have added
                  CALL xios_orchidee_send_field('EW_BASALT_ADDED',EW_grain_add(:,:,ibasalt,1)) ! divided by sec per year

!=======
!
!>>>>>>> .r5976
                  ENDIF

               ELSE
                  EW_flux = zero
                  EW_grain (:,:,:,1) = zero
                  EW_grain (:,:,:,2) = undef
               ENDIF
               p_input(:,:,iweat_p) = SPREAD(p_release_daily(:),NCOPIES=nvm,DIM=2)  &
                                      + EW_flux(:,:,iPrelease)
           ENDIF                                                    
        ENDIF
        ! DSG: this routine (stomate_weathering) doesn't manipulate mass, so no mass conservation
        ! check needed here

       
       !! 5.3 Use all processes included in stomate

       !! 5.3.1  Activate stomate processes 
       ! Activate stomate processes (the complete list of processes depends 
       ! on whether the DGVM is used or not). Processes include: climate constraints 
       ! for PFTs, PFT dynamics, Phenology, Allocation, NPP (based on GPP and
       ! authothropic respiration), fire, mortality, vmax, assimilation temperatures,
       ! all turnover processes, light competition, sapling establishment, lai and 
       ! land cover change.

       IF (printlev>=3) WRITE (numout,*) 'BEFORE StomateLPJ soil_n_min(test_grid,test_pft,:):',soil_n_min(test_grid,test_pft,:)

       CALL StomateLpj &
            &            (kjpindex, dt_days, &
            &             neighbours, resolution, &
            &             herbivores, &
            &             soil_orders,max_eau_var, bulk, &
            &             tsurf_daily, tsoil_daily, t2m_daily, t2m_min_daily, &
            &             litterhum_daily, soilhum_daily, &
            &             maxhumrel_lastyear, minhumrel_lastyear, &
            &             gdd0_lastyear, precip_lastyear, &
            &             humrel_month, humrel_week, t2m_longterm, t2m_month, t2m_week, &
            &             tsoil_month, soilhum_month, &
            &             gdd_m5_dormance,  gdd_from_growthinit, gdd_midwinter, ncd_dormance, ngd_minus5, &
            &             turnover_longterm, gpp_daily, &
            &             gpp_week, &
            &             time_hum_min, maxfpc_lastyear, resp_maint_part,&
            &             PFTpresent, age, fireindex, firelitter, &
            &             leaf_age, leaf_frac, biomass, ind, adapted, regenerate, &
            &             senescence, when_growthinit, litter, &
            &             litter_avail, litter_not_avail, litter_avail_frac, &
            &             dead_leaves, som, lignin_struc, lignin_wood, &
            &             veget_cov_max,veget_cov, veget_cov_max_new, woodharvest, npp_longterm, lm_lastyearmax, &
            &             veget_lastlight, everywhere, need_adjacent, RIP_time, &
            &             lai, rprof,npp_daily, turnover_daily, turnover_time,&
            &             control_moist_inst, control_temp_inst, som_input_inst, &
            &             co2_to_bm_dgvm, n_to_bm, p_to_bm,co2_fire, &
            &             resp_hetero_d, resp_maint_d, resp_growth_d, &
            &             height, deadleaf_cover, vcmax, jmax, nue, &
            &             bm_to_litter,&
            &             prod10, prod100, flux10, flux100, &
            &             convflux, cflux_prod10, cflux_prod100, &
            &             prod10_harvest, prod100_harvest, flux10_harvest, flux100_harvest, &
            &             convflux_harvest, cflux_prod10_harvest, cflux_prod100_harvest, woodharvestpft, & 
            &             harvest_above, harvest_bioenergy,carb_mass_total, &
            &             fpc_max, matrixA, &
            &             Tseason, Tmin_spring_time, begin_leaves, onset_date, KF, k_latosa_adapt,&
            &             cn_leaf_avg_season, np_leaf_avg_season, nstress_season, &
            &             pstress_season, &
            &             lai_target, lai_target_longterm, &
            &             moiavail_growingseason, soil_n_min, soil_p_min, rue_longterm, &
            &             n_uptake_daily, p_uptake_daily, &
            &             record_reserve_grass, &
            &             wshtotsum, sr_ugb, compt_ugb, nb_ani, grazed_frac, &
            &             import_yield, sla_age1, t2m_14, sla_calc,snowfall_daily,MODULO(days_since_beg,365),&
            &             when_growthinit_cut, nb_grazingdays, &
            &             humrel_daily,tmc_topgrass_daily,fc_grazing,snowmass_daily, &
            &             after_snow, after_wet, wet1day, wet2day, &
            &             grm_nfert,grm_pfert,BNF_clover_daily, &
            &             GRM_devstage,agri_fert_save,agri_nfert, &
            &             agri_pfert, rest_Ngrazing, grazed_above, &
            &             days_senescence)

       DO j = 2,nvm
          WHERE ( k_latosa_adapt(:,j) .GE. k_latosa_max(j) )
                
             k_latosa_adapt(:,j) = k_latosa_max(j)
             !!JC Debug MOD018
             !! whether it is appropriate to limit the k_latosa_adapt between max and min?
             ! NOTE: during a long spinup, if we want to change k_latosa during
             ! the simulation (e.g., by reset k_latosa_max k_latosa_min, we have to activate it.
          ELSEWHERE (k_latosa_adapt(:,j) .LE. k_latosa_min(j) )
             
             k_latosa_adapt(:,j) = k_latosa_min(j)
             !!End JC Debug
          ENDWHERE

       ENDDO
       !! Mass conservation check is done in stomatelpj; thus we don't need to
       !do it here

       IF (printlev>=3) WRITE (numout,*) 'AFTER StomateLPJ soil_n_min(test_grid,test_pft,:):',soil_n_min(test_grid,test_pft,:)

       CALL histwrite_p (hist_id_stomate, 'K_LATOSA_ADAPT', itime, &
            k_latosa_adapt(:,:), kjpindex*nvm, horipft_index)
       CALL xios_orchidee_send_field('K_LATOSA_ADAPT',k_latosa_adapt(:,:))
       
       !! 5.3.2 Calculate the total CO2 flux from land use change
       fco2_lu(:) = convflux(:) &
            &             + cflux_prod10(:)  &
            &             + cflux_prod100(:) &
            &             + harvest_above(:) &
            &             + convflux_harvest(:) &
            &             + cflux_prod10_harvest(:)  &
            &             + cflux_prod100_harvest(:)
       
       !! 5.4 Calculate veget and veget_max
       veget_max(:,:) = zero 
       DO j = 1, nvm
          veget_max(:,j) = veget_max(:,j) + &
               & veget_cov_max(:,j) * ( 1.-totfrac_nobio(:) )
       ENDDO
       
       !! 5.5 Photosynthesis parameters
       assim_param(:,:,ivcmax) = zero
       assim_param(:,:,ijmax)  = zero
       assim_param(:,:,inue)   = zero
       assim_param(:,:,ileafn) = zero
!DSGpot    
       assim_param(:,:,ivcmax_pot) = zero
       assim_param(:,:,ijmax_pot)  = zero
!DSGpot    
       DO j = 2,nvm
          assim_param(:,j,ivcmax) = vcmax(:,j,1)
          assim_param(:,j,ijmax)  = jmax(:,j,1)
          assim_param(:,j,inue)   = nue(:,j)
          assim_param(:,j,ileafn) = biomass(:,j,ileaf,initrogen) & 
                                    - (SNC*biomass(:,j,ileaf,icarbon)) ! account for structural leaf N 
!DSGpot    
          assim_param(:,j,ivcmax_pot) = vcmax(:,j,2)
          assim_param(:,j,ijmax_pot)  = jmax(:,j,2)
!DSGpot    
       ENDDO
       
       !! 5.6 Update forcing variables for soil carbon
       IF (TRIM(Cforcing_name) /= 'NONE') THEN
          npp_tot(:) = 0
          DO j=2,nvm
             npp_tot(:) = npp_tot(:) + npp_daily(:,j)
          ENDDO
          ! ::nbyear Number of years saved for carbon spinup
          sf_time = MODULO(REAL(days_since_beg,r_std)-1,one_year*REAL(nbyear,r_std))
          iatt=FLOOR(sf_time/dt_forcesoil) + 1
          IF (iatt == 0) iatt = iatt_old + 1
          IF ((iatt<iatt_old) .and. (.not. cumul_Cforcing)) THEN
             nforce(:)=0
             carbon_input(:,:,:,:) = zero
             nitrogen_input(:,:,:,:) = zero
             control_moist(:,:,:) = zero
             control_temp(:,:,:) = zero
             npp_equil(:,:) = zero
             drainage(:,:,:) = zero
             runoff(:,:,:) = zero
          ENDIF
          iatt_old = iatt
          ! Update forcing
          nforce(iatt) = nforce(iatt)+1
          carbon_input(:,:,:,iatt) = carbon_input(:,:,:,iatt) + som_input_daily(:,:,:,icarbon)
          nitrogen_input(:,:,:,iatt) = nitrogen_input(:,:,:,iatt) + som_input_daily(:,:,:,initrogen)
          control_moist(:,:,iatt) = control_moist(:,:,iatt) + control_moist_daily(:,:)
          control_temp(:,:,iatt) = control_temp(:,:,iatt) + control_temp_daily(:,:)
          npp_equil(:,iatt)      = npp_equil(:,iatt) + npp_tot(:)
          drainage(:,:,iatt)     = drainage(:,:,iatt) + drainage_daily(:,:)
          runoff  (:,:,iatt)     = runoff  (:,:,iatt) + runoff_daily  (:,:)
       ENDIF
       
       !! 5.8 Write forcing file if ::ok_co2=.TRUE.
       ! Note: if STOMATE is run in coupled mode the forcing file is written
       ! If run in stand-alone mode, the forcing file is read!
       IF ( ok_co2 .AND. TRIM(forcing_name) /= 'NONE' ) THEN
          
          !! 5.8.1 Convert GPP to sechiba time steps
          ! GPP is multiplied by coverage to obtain forcing @tex $(gC m^{-2} dt_stomate^{-1})$\f \end@tex $(m^2 m^{-2})$ @endtexonly
          ! @tex$ m^{-2}$ @endtex remains in the units because ::veget_cov_max is a fraction, not a 
          ! surface area. In sechiba values are ponderated by surface and frac_no_bio. 
          ! At the beginning of stomate, the units are converted. 
          ! When we use forcesoil we call sechiba_main and so we need the have the same units as in sechiba.
          gpp_daily_x(:,:) = zero
          DO j = 2, nvm             
             gpp_daily_x(:,j) = gpp_daily_x(:,j) + &
              & gpp_daily(:,j) * dt_stomate / one_day * veget_cov_max(:,j)
          ENDDO
          
          ! Bare soil moisture availability has not been treated
          ! in STOMATE, update it here
          humrel_daily(:,ibare_sechiba) = humrel(:,ibare_sechiba)   

          ! Update index to store the next forcing step in memory
          iisf = iisf+1

          ! How many times have we treated this forcing state
          xn = REAL(nf_cumul(isf(iisf)),r_std)
          
          !! 5.8.2 Cumulate forcing variables
          ! Cumulate forcing variables (calculate average)
          ! Note: precipitation is multiplied by dt_stomate/one_day to be consistent with 
          ! the units in sechiba
          IF (cumul_forcing) THEN
             clay_fm(:,iisf) = (xn*clay_fm(:,iisf)+clay(:))/(xn+1.)
             silt_fm(:,iisf) = (xn*silt_fm(:,iisf)+silt(:))/(xn+1.) 
             bulk_fm(:,iisf) = (xn*bulk_fm(:,iisf)+bulk(:))/(xn+1.) 
             humrel_daily_fm(:,:,iisf) = &
                  & (xn*humrel_daily_fm(:,:,iisf) + humrel_daily(:,:))/(xn+1.)
             litterhum_daily_fm(:,iisf) = &
                  & (xn*litterhum_daily_fm(:,iisf)+litterhum_daily(:))/(xn+1.)
             t2m_daily_fm(:,iisf) = &
                  & (xn*t2m_daily_fm(:,iisf)+t2m_daily(:))/(xn+1.)
             t2m_min_daily_fm(:,iisf) = &
                  & (xn*t2m_min_daily_fm(:,iisf)+t2m_min_daily(:))/(xn+1.)
             tsurf_daily_fm(:,iisf) = &
                  & (xn*tsurf_daily_fm(:,iisf)+tsurf_daily(:))/(xn+1.)
             tsoil_daily_fm(:,:,iisf) = &
                  & (xn*tsoil_daily_fm(:,:,iisf)+tsoil_daily(:,:))/(xn+1.)
             soilhum_daily_fm(:,:,iisf) = &
                  & (xn*soilhum_daily_fm(:,:,iisf)+soilhum_daily(:,:))/(xn+1.)
             precip_fm(:,iisf) = &
                  & (xn*precip_fm(:,iisf)+precip_daily(:)*dt_stomate/one_day)/(xn+1.)
             gpp_daily_fm(:,:,iisf) = &
                  & (xn*gpp_daily_fm(:,:,iisf) + gpp_daily_x(:,:))/(xn+1.)
             veget_fm(:,:,iisf) = &
                  & (xn*veget_fm(:,:,iisf) + veget(:,:) )/(xn+1.)
             veget_max_fm(:,:,iisf) = &
                  & (xn*veget_max_fm(:,:,iisf) + veget_max(:,:) )/(xn+1.)
             lai_fm(:,:,iisf) = &
                  & (xn*lai_fm(:,:,iisf) + lai(:,:) )/(xn+1.)
             drainage_fm(:,:,iisf) = & 
                  & (xn*drainage_fm(:,:,iisf) + drainage_daily(:,:) )/(xn+1.) 
             runoff_fm(:,:,iisf) = & 
                  & (xn*runoff_fm(:,:,iisf) + runoff_daily(:,:) )/(xn+1.) 
             max_eau_var_daily_fm(:,iisf) = &
                  & (xn*max_eau_var_daily_fm(:,iisf) + max_eau_var_daily(:))/(xn+1.)
          ELSE
             ! Here we just calculate the values
             clay_fm(:,iisf) = clay(:)
             silt_fm(:,iisf) = silt(:) 
             bulk_fm(:,iisf) = bulk(:) 
             humrel_daily_fm(:,:,iisf) = humrel_daily(:,:)
             litterhum_daily_fm(:,iisf) = litterhum_daily(:)
             t2m_daily_fm(:,iisf) = t2m_daily(:)
             t2m_min_daily_fm(:,iisf) =t2m_min_daily(:)
             tsurf_daily_fm(:,iisf) = tsurf_daily(:)
             tsoil_daily_fm(:,:,iisf) =tsoil_daily(:,:)
             soilhum_daily_fm(:,:,iisf) =soilhum_daily(:,:)
             precip_fm(:,iisf) = precip_daily(:)
             gpp_daily_fm(:,:,iisf) =gpp_daily_x(:,:)
             veget_fm(:,:,iisf) = veget(:,:)
             veget_max_fm(:,:,iisf) =veget_max(:,:)
             lai_fm(:,:,iisf) =lai(:,:)
             drainage_fm(:,:,iisf) = drainage_daily(:,:) 
             runoff_fm  (:,:,iisf) = runoff_daily  (:,:) 
             max_eau_var_daily_fm(:,iisf) = max_eau_var_daily(:)
          ENDIF
          nf_cumul(isf(iisf)) = nf_cumul(isf(iisf))+1

          ! 5.8.3 Do we have to write the forcing states?
          IF (iisf == nsfm) THEN

             !! 5.8.3.1 Write these forcing states
             CALL forcing_write(forcing_id,1,nsfm)
             ! determine which forcing states must be read
             isf(1) = isf(nsfm)+1
             IF ( isf(1) > nsft ) isf(1) = 1
             DO iisf = 2, nsfm
                isf(iisf) = isf(iisf-1)+1
                IF (isf(iisf) > nsft)  isf(iisf) = 1
             ENDDO

             ! Read forcing variables - for debug use only
             ! CALL forcing_read(forcing_id,nsfm)
             iisf = 0

          ENDIF

       ENDIF

       !! 5.9 Compute daily CO2 flux diagnostics
       ! CO2 flux in @tex $gC m^{-2} s^{-1}$ @endtex (positive from atmosphere to land) is sum of:
       !   (1) co2 taken up by photosyntyhesis + (2) co2 taken up in the DGVM to establish saplings 
       ! - (3) plants maintenance respiration  - (4) plants growth respiration
       ! - (5) heterotrophic respiration from ground
       ! - (6) co2 emission from fire
       nep_daily(:,:)= gpp_daily(:,:)     + co2_to_bm_dgvm(:,:)  &
                     - resp_maint_d(:,:)  - resp_growth_d(:,:)   &
                     - resp_hetero_d(:,:) - co2_fire(:,:) 

       CALL xios_orchidee_send_field("nep",SUM(nep_daily*veget_cov_max,dim=2)/1e3/one_day)

       IF ( hist_id_stom_IPCC > 0 ) THEN
          vartmp(:) = SUM(nep_daily*veget_cov_max,dim=2)/1e3/one_day*contfrac
          CALL histwrite_p (hist_id_stom_IPCC, "nep", itime, &
               vartmp, kjpindex, hori_index)
       ENDIF

       ! Cumulate nep, harvest and land use change fluxes
       nep_monthly(:,:) = nep_monthly(:,:) + nep_daily(:,:)
       harvest_above_monthly(:) = harvest_above_monthly(:) + harvest_above(:)
       cflux_prod_monthly(:) = cflux_prod_monthly(:) + convflux(:) + & 
        & cflux_prod10(:) + cflux_prod100(:) + convflux_harvest(:) + & 
        & cflux_prod10_harvest(:) + cflux_prod100_harvest(:)
      

       IF(do_slow .AND. dsg_debug) THEN

            ! what is the change in total mass:
            mass_change(:,:,iphosphorus) =   SUM(p_input(:,:,:),DIM=1) - p_leaching_daily(:,:) - p_ssorption_daily(:,:)
            !DSGmass_dbg
            mass_change(:,:,iphosphorus) =   SUM(p_input(:,:,:),DIM=1) - p_leaching_daily(:,:) - p_ssorption_daily(:,:) + p_to_bm(:,:)
            !DSGmass_dbg
            ! DSG: his mass change lacks some fluxes:
            !mass_change(:,:,initrogen)   =   SUM(n_input(:,:,:),DIM=1) - n_leaching_daily(:,:) ! + ?
            ! DSG: this mass change lacks some fluxes:
            mass_change(:,:,icarbon)     =  nep_daily(:,:)

            ! Mass at the end of the timestep (accounting for mass which was removed from a pool but not added to another pool YET) 
            mass_after(:,:,:)            = SUM(som(:,:,:,:),DIM=2) +  SUM(SUM(litter(:,:,:,:,:),DIM=4),DIM=2) + SUM(biomass(:,:,:,:),DIM=3) 
            ! we need to correct for the litterfall which is already removed from biomass but not yet added to litter pools
            mass_after(:,:,:)            = mass_after(:,:,:) + SUM(turnover_daily(:,:,:,:),DIM=3) + SUM(bm_to_litter(:,:,:,:),DIM=3)
            ! add the mineal nutrient pools
            mass_after(:,:,iphosphorus)  = mass_after(:,:,iphosphorus) + SUM(soil_p_min(:,:,:),DIM=3)
            mass_after(:,:,initrogen)    = mass_after(:,:,initrogen)   + SUM(soil_n_min(:,:,:),DIM=3)

            ! we dont check C and N for now as I am not clear about the mass change (variables are lacking)
            mass_before_2(:,:,icarbon)   = zero
            mass_after   (:,:,icarbon)   = zero
            mass_change  (:,:,icarbon)   = zero
            mass_before_2(:,:,initrogen) = zero
            mass_after   (:,:,initrogen) = zero
            mass_change  (:,:,initrogen) = zero

            WRITE(6,*) 'location in code: '
            WRITE(6,*) 'stomate: end of stomate'
            CALL cons_mass(mass_before_2(:,:,:),             &  ! mass before
                  mass_after(:,:,:),                         &  ! mass after
                  mass_change(:,:,:),                        &  ! net of fluxes
                  'END of stomate') 
       ENDIF

       !! 5.10 Compute monthly CO2 fluxes 
       IF ( LastTsMonth ) THEN
          !! 5.10.1 Write history file for monthly fluxes
          CALL histwrite_p (hist_id_stomate, 'CO2FLUX', itime, &
               nep_monthly, kjpindex*nvm, horipft_index)
          CALL xios_orchidee_send_field('CO2FLUX',nep_monthly)
          
          ! Integrate nep_monthly over all grid-cells on local domain
          net_nep_monthly = zero
          DO ji=1,kjpindex
             DO j=2,nvm
                net_nep_monthly = net_nep_monthly + &
                     nep_monthly(ji,j)*resolution(ji,1)*resolution(ji,2)*contfrac(ji)*veget_cov_max(ji,j)
             ENDDO
          ENDDO
          ! Change unit from gC/m2 grid-cell into PgC/m2 grid-cell
          net_nep_monthly = net_nep_monthly*1e-15

     
          !! 5.10.2 Cumulative fluxes of land use cover change, harvest and net biosphere production
          ! Parallel processing, gather the information from different processors. first argument is the
          ! local variable, the second argument is the global variable. bcast send it to all processors.
          net_cflux_prod_monthly_sum = &
              &  SUM(cflux_prod_monthly(:)*resolution(:,1)*resolution(:,2)*contfrac(:))*1e-15
          CALL reduce_sum(net_cflux_prod_monthly_sum,net_cflux_prod_monthly_tot)
          CALL bcast(net_cflux_prod_monthly_tot)
          net_harvest_above_monthly_sum = &
             &   SUM(harvest_above_monthly(:)*resolution(:,1)*resolution(:,2)*contfrac(:))*1e-15
          CALL reduce_sum(net_harvest_above_monthly_sum,net_harvest_above_monthly_tot)
          CALL bcast(net_harvest_above_monthly_tot)
          CALL reduce_sum(net_nep_monthly,net_nep_monthly_sum)
          CALL bcast(net_nep_monthly_sum)
          net_biosp_prod_monthly_tot =  net_cflux_prod_monthly_tot + net_harvest_above_monthly_tot - net_nep_monthly_sum
          
          WRITE(numout,9010) 'GLOBAL net_cflux_prod_monthly    (Peta gC/month)  = ',net_cflux_prod_monthly_tot
          WRITE(numout,9010) 'GLOBAL net_harvest_above_monthly (Peta gC/month)  = ',net_harvest_above_monthly_tot
          WRITE(numout,9010) 'GLOBAL net_nep_monthly           (Peta gC/month)  = ',net_nep_monthly_sum
          WRITE(numout,9010) 'GLOBAL net_biosp_prod_monthly    (Peta gC/month)  = ',net_biosp_prod_monthly_tot

9010  FORMAT(A52,F17.14) 

          ! Reset Monthly values
          nep_monthly(:,:) = zero
          harvest_above_monthly(:) = zero
          cflux_prod_monthly(:)    = zero

       ENDIF ! Monthly processes - at the end of the month
       
       IF (spinup_analytic) THEN
          nbp_accu(:) = nbp_accu(:) + (SUM(nep_daily(:,:) * veget_max(:,:),dim=2) - (convflux(:) + cflux_prod10(:) + &
                    cflux_prod100(:)) - (convflux_harvest(:) + cflux_prod10_harvest(:) + &
                    cflux_prod100_harvest(:))  - harvest_above(:))/1e3 
       ENDIF

       !! 5.11 Reset daily variables
       humrel_daily(:,:)       = zero
       max_eau_var_daily(:)    = zero
       litterhum_daily(:)      = zero
       t2m_daily(:)            = zero
       t2m_min_daily(:)        = large_value
       tsurf_daily(:)          = zero
       tsoil_daily(:,:)        = zero
       soilhum_daily(:,:)      = zero
       precip_daily(:)         = zero
       gpp_daily(:,:)          = zero
       resp_maint_part(:,:,:)  = zero
       resp_hetero_d           = zero
       drainage_daily(:,:)     = zero 
       runoff_daily  (:,:)     = zero 
       n_uptake_daily(:,:,:)   = zero 
       p_uptake_daily(:,:)     = zero 
       p_ssorption_daily(:,:)  = zero 
       p_leaching_daily(:,:)   = zero 
       n_mineralisation_d(:,:) = zero 
       p_mineralisation_d(:,:) = zero 

       BNF_clover_daily(:,:)  = zero

       IF (printlev >= 3) THEN
          WRITE(numout,*) 'stomate_main: daily processes done'
       ENDIF

    ENDIF  ! Daily processes - at the end of the day

    
  !! 6. Outputs from Stomate

    ! co2_flux receives a value from STOMATE only if STOMATE is activated.
    ! Otherwise, the calling hydrological module must do this itself.

    !! 6.1 Respiration and fluxes
    resp_maint(:,:) = resp_maint_radia(:,:)*veget_cov_max(:,:)
    resp_maint(:,ibare_sechiba) = zero
    resp_growth(:,:)= resp_growth_d(:,:)*veget_cov_max(:,:)*dt_sechiba/one_day
    resp_hetero(:,:) = resp_hetero_radia(:,:)*veget_cov_max(:,:)
    
    !! 6.2 Derived CO2 fluxes
    ! CO2 flux in gC m^{-2} s^{-1} (positive towards the atmosphere) is sum of:
    ! (1) heterotrophic respiration from ground + (2) maintenance respiration 
    ! from the plants + (3) growth respiration from the plants + (4) co2 
    ! emissions from fire - (5) co2 taken up in the DGVM to establish 
    ! saplings - (6) co2 taken up by photosyntyhesis
    co2_flux(:,:) = resp_hetero(:,:) + resp_maint(:,:) + resp_growth(:,:) &
         & + (co2_fire(:,:)-co2_to_bm_dgvm(:,:))*veget_cov_max(:,:)/one_day &
         & - gpp(:,:)
    
    temp_growth(:)=t2m_month(:)-tp_00 



   
  !! 7. Analytical spinup

    IF (spinup_analytic) THEN

       tau_CN_longterm = MIN(tau_CN_longterm + dt_sechiba/one_day,spinup_period*one_year)

       !! 7.1. Update V and U at sechiba time step
       DO m = 2,nvm
          DO j = 1,kjpindex 
             ! V <- A * V
             MatrixV(j,m,:,:) = MATMUL(matrixA(j,m,:,:),MatrixV(j,m,:,:))
             ! U <- A*U + B
             VectorU(j,m,:)   = MATMUL(matrixA(j,m,:,:),VectorU(j,m,:)) + vectorB(j,m,:)
          ENDDO ! loop pixels
       ENDDO ! loop PFTS

       !! 7.2. What happened at the end of the year ?
       IF (LastTsYear) THEN

          !
          ! 7.2.1 Increase the years counter every LastTsYear which is the last sechiba time step of each year
          !
          global_years = global_years + 1 

          !
          ! 7.2.3 Is global_years is a multiple of the period time ?
          !

          !
          ! 3.2.1 When global_years is a multiple of the spinup_period, we calculate :
          !       1) the mean nbp flux over the period. This value is restarted
          !       2) we solve the matrix system by Gauss Jordan method
          !       3) We test if a point is at equilibrium : if yes, we mark the point (ok_equilibrium array)
          !       4) Then we reset the matrix 
          !       5) We erase the carbon_stock calculated by ORCHIDEE by the one found by the method
          IF( MOD(global_years, spinup_period) == 0 ) THEN
             WRITE(numout,*) 'Spinup analytic : Calculate if system is in equlibrium. global_years=',global_years
             ! The number total of days during the forcing period is given by :
             !    spinup_period*365 (we consider only the noleap calendar)
             nbp_flux(:) = nbp_accu(:) / ( spinup_period * 365.)
             ! Reset the values
             nbp_accu(:) = zero

             carbon_stock(:,ibare_sechiba,:) = zero
             ! Prepare the matrix for the resolution
             ! Add a temporary matrix W which contains I-MatrixV
             ! we should take the opposite of matrixV and add the identitiy : we solve (I-MatrixV)*C = VectorU
             MatrixW(:,:,:,:) = moins_un * MatrixV(:,:,:,:)
             DO jv = 1,nbpools
                MatrixW(:,:,jv,jv) =  MatrixW(:,:,jv,jv) + un
             ENDDO
             carbon_stock(:,:,:) = VectorU(:,:,:)

             !
             !  Solve the linear system
             !
             DO m = 2,nvm
                DO j = 1,kjpindex
                   ! the solution will be stored in VectorU : so it should be restarted before
                   ! loop over npts and nvm, so we solved npts*(nvm-1) (7,7) linear systems
                   CALL gauss_jordan_method(nbpools,MatrixW(j,m,:,:),carbon_stock(j,m,:))
                ENDDO ! loop pixels
             ENDDO ! loop PFTS

             ! Reset temporary matrixW
             MatrixW(:,:,:,:) = zero 

             previous_stock(:,:,:) = current_stock(:,:,:)
             current_stock(:,:,:)  = carbon_stock(:,:,:)  

             ! The relative error is calculated over the passive carbon pool (sum over the pfts) over the pixel.
             CALL error_L1_passive(kjpindex,nvm, nbpools, current_stock, previous_stock, veget_max, &
                  &                eps_carbon, carbon_eq)   

             !! ok_equilibrium is saved,
             WHERE( carbon_eq(:) .AND. .NOT.(ok_equilibrium(:)) )
                ok_equilibrium(:) = .TRUE.  
             ENDWHERE

             ! Reset matrixV for the pixel to the identity matrix and vectorU to zero
             MatrixV(:,:,:,:) = zero
             VectorU(:,:,:) = zero
             DO jv = 1,nbpools
                MatrixV(:,:,jv,jv) = un
             END DO
             IF (printlev >= 2) WRITE(numout,*) 'Reset for matrixV and VectorU done'    

             !! Write the values found in the standard outputs of ORCHIDEE
             litter(:,istructural,:,iabove,icarbon) = carbon_stock(:,:,istructural_above)
             litter(:,istructural,:,ibelow,icarbon) = carbon_stock(:,:,istructural_below)
             litter(:,imetabolic,:,iabove,icarbon)  = carbon_stock(:,:,imetabolic_above)
             litter(:,imetabolic,:,ibelow,icarbon)  = carbon_stock(:,:,imetabolic_below)
             litter(:,iwoody,:,iabove,icarbon)      = carbon_stock(:,:,iwoody_above)
             litter(:,iwoody,:,ibelow,icarbon)      = carbon_stock(:,:,iwoody_below)
             som(:,iactive,:,icarbon)   = carbon_stock(:,:,iactive_pool)
             som(:,isurface,:,icarbon)  = carbon_stock(:,:,isurface_pool)
             som(:,islow,:,icarbon)     = carbon_stock(:,:,islow_pool)
             som(:,ipassive,:,icarbon)  = carbon_stock(:,:,ipassive_pool) 

             !! Initialize the nitrogen & phosphorus stocks using the longterm CNP ratios and
             !! the new carbon stocks
             litter(:,:,:,:,initrogen)   = zero
             litter(:,:,:,:,iphosphorus) = zero
             som(:,:,:,initrogen)        = zero
             som(:,:,:,iphosphorus)      = zero
             
             WHERE( (CN_som_litter_longterm(:,:,istructural_above) .GT. zero) .AND. & 
                    (CP_som_litter_longterm(:,:,istructural_above) .GT. zero) )
                litter(:,istructural,:,iabove,initrogen) = carbon_stock(:,:,istructural_above) &
                     / CN_som_litter_longterm(:,:,istructural_above)

                litter(:,istructural,:,iabove,iphosphorus) = carbon_stock(:,:,istructural_above) & 
                     / CP_som_litter_longterm(:,:,istructural_above)
             ELSEWHERE
                litter(:,istructural,:,iabove,icarbon)     = zero
                litter(:,istructural,:,iabove,initrogen)   = zero
                litter(:,istructural,:,iabove,iphosphorus) = zero
             ENDWHERE

             WHERE( (CN_som_litter_longterm(:,:,istructural_below) .GT. zero) .AND. &
                    (CP_som_litter_longterm(:,:,istructural_below) .GT. zero) )
                litter(:,istructural,:,ibelow,initrogen) = carbon_stock(:,:,istructural_below) &
                     / CN_som_litter_longterm(:,:,istructural_below)
                litter(:,istructural,:,ibelow,iphosphorus) = carbon_stock(:,:,istructural_below) &
                     / CP_som_litter_longterm(:,:,istructural_below)
             ELSEWHERE
                litter(:,istructural,:,ibelow,icarbon)     = zero
                litter(:,istructural,:,ibelow,initrogen)   = zero
                litter(:,istructural,:,ibelow,iphosphorus) = zero
             ENDWHERE

             WHERE( (CN_som_litter_longterm(:,:,imetabolic_above).GT. zero) .AND. &
                    (CP_som_litter_longterm(:,:,imetabolic_above).GT. zero) )
                litter(:,imetabolic,:,iabove,initrogen)  = carbon_stock(:,:,imetabolic_above)  &
                     / CN_som_litter_longterm(:,:,imetabolic_above)
                litter(:,imetabolic,:,iabove,iphosphorus)  = carbon_stock(:,:,imetabolic_above)  &
                     / CP_som_litter_longterm(:,:,imetabolic_above)
             ELSEWHERE
                litter(:,imetabolic,:,iabove,icarbon)      = zero
                litter(:,imetabolic,:,iabove,initrogen)    = zero
                litter(:,imetabolic,:,iabove,iphosphorus)  = zero
             ENDWHERE

             WHERE( (CN_som_litter_longterm(:,:,imetabolic_below) .GT. zero ) .AND. &      
                    (CP_som_litter_longterm(:,:,imetabolic_below) .GT. zero ))
                litter(:,imetabolic,:,ibelow,initrogen)  = carbon_stock(:,:,imetabolic_below)  &
                     / CN_som_litter_longterm(:,:,imetabolic_below)
                litter(:,imetabolic,:,ibelow,iphosphorus)  = carbon_stock(:,:,imetabolic_below)  &
                     / CP_som_litter_longterm(:,:,imetabolic_below)
             ELSEWHERE
                litter(:,imetabolic,:,ibelow,icarbon)      = zero
                litter(:,imetabolic,:,ibelow,initrogen)    = zero
                litter(:,imetabolic,:,ibelow,iphosphorus)  = zero
             ENDWHERE

             WHERE( (CN_som_litter_longterm(:,:,iwoody_above) .GT. zero) .AND. &
                    (CP_som_litter_longterm(:,:,iwoody_above) .GT. zero) )
                litter(:,iwoody,:,iabove,initrogen)      = carbon_stock(:,:,iwoody_above)      &
                     / CN_som_litter_longterm(:,:,iwoody_above)    
                litter(:,iwoody,:,iabove,iphosphorus)      = carbon_stock(:,:,iwoody_above)      &
                     / CP_som_litter_longterm(:,:,iwoody_above)    
             ELSEWHERE
                litter(:,iwoody,:,iabove,icarbon)      = zero
                litter(:,iwoody,:,iabove,initrogen)    = zero
                litter(:,iwoody,:,iabove,iphosphorus)  = zero
             ENDWHERE

             WHERE( (CN_som_litter_longterm(:,:,iwoody_below) .GT. zero) .AND. &
                    (CP_som_litter_longterm(:,:,iwoody_below) .GT. zero) )
                litter(:,iwoody,:,ibelow,initrogen)      = carbon_stock(:,:,iwoody_below)      &
                     / CN_som_litter_longterm(:,:,iwoody_below)    
                litter(:,iwoody,:,ibelow,iphosphorus)      = carbon_stock(:,:,iwoody_below)      &
                     / CP_som_litter_longterm(:,:,iwoody_below)    
             ELSEWHERE
                litter(:,iwoody,:,ibelow,icarbon)      = zero
                litter(:,iwoody,:,ibelow,initrogen)    = zero
                litter(:,iwoody,:,ibelow,iphosphorus)  = zero
             ENDWHERE

             WHERE( (CN_som_litter_longterm(:,:,iactive_pool) .GT. zero) .AND. &
                    (CP_som_litter_longterm(:,:,iactive_pool) .GT. zero) )
                som(:,iactive,:,initrogen)               = carbon_stock(:,:,iactive_pool)      &
                     / CN_som_litter_longterm(:,:,iactive_pool)    
                som(:,iactive,:,iphosphorus)               = carbon_stock(:,:,iactive_pool)      &
                     / CP_som_litter_longterm(:,:,iactive_pool)    
             ELSEWHERE
                som(:,iactive,:,icarbon)               = zero
                som(:,iactive,:,initrogen)             = zero
                som(:,iactive,:,iphosphorus)           = zero
             ENDWHERE

             WHERE((CN_som_litter_longterm(:,:,isurface_pool) .GT. zero) .AND. &
                   (CP_som_litter_longterm(:,:,isurface_pool) .GT. zero) )
                som(:,isurface,:,initrogen)              = carbon_stock(:,:,isurface_pool)     &
                     / CN_som_litter_longterm(:,:,isurface_pool)   
                som(:,isurface,:,iphosphorus)              = carbon_stock(:,:,isurface_pool)     &
                     / CP_som_litter_longterm(:,:,isurface_pool)   
             ELSEWHERE
                som(:,isurface,:,icarbon)              = zero
                som(:,isurface,:,initrogen)            = zero
                som(:,isurface,:,iphosphorus)          = zero
             ENDWHERE

             WHERE((CN_som_litter_longterm(:,:,islow_pool) .GT. zero) .AND. &
                   (CP_som_litter_longterm(:,:,islow_pool) .GT. zero) )
                som(:,islow,:,initrogen)                 = carbon_stock(:,:,islow_pool)        &
                     / CN_som_litter_longterm(:,:,islow_pool)      
                som(:,islow,:,iphosphorus)                 = carbon_stock(:,:,islow_pool)        &
                     / CP_som_litter_longterm(:,:,islow_pool)      
             ELSEWHERE
                som(:,islow,:,icarbon)                 = zero
                som(:,islow,:,initrogen)               = zero
                som(:,islow,:,iphosphorus)             = zero
             ENDWHERE

             WHERE((CN_som_litter_longterm(:,:,ipassive_pool) .GT. zero).AND. &
                   (CP_som_litter_longterm(:,:,ipassive_pool) .GT. zero) )
                som(:,ipassive,:,initrogen)              = carbon_stock(:,:,ipassive_pool)     &
                     / CN_som_litter_longterm(:,:,ipassive_pool)     
                som(:,ipassive,:,iphosphorus)              = carbon_stock(:,:,ipassive_pool)     &
                     / CP_som_litter_longterm(:,:,ipassive_pool)     
             ELSEWHERE
                som(:,ipassive,:,icarbon)              = zero
                som(:,ipassive,:,initrogen)            = zero
                som(:,ipassive,:,iphosphorus)          = zero
             ENDWHERE

             ! DSG: we set the timer back to start 
             tau_CN_longterm = dt_sechiba/one_day

             ! Final step, test if all points at the local domain are at equilibrium
             ! The simulation can be stopped when all local domains have reached the equilibrium
             IF (printlev >=1) THEN
                IF (ALL(ok_equilibrium)) THEN
                   WRITE(numout,*) 'Spinup analytic : Equilibrium for carbon pools is reached for current local domain'
                ELSE
                   WRITE(numout,*) 'Spinup analytic : Equilibrium for carbon pools is not yet reached for current local domain'
                END IF
             END IF
             
          ENDIF ! ( MOD(global_years,spinup_period) == 0)



!       ELSE  ! useless

       ENDIF ! (LastTsYear)

    ENDIF !(spinup_analytic)
   
    IF (printlev >= 4) WRITE(numout,*) 'Leaving stomate_main'

  END SUBROUTINE stomate_main

!! ================================================================================================================================
!! SUBROUTINE 	: stomate_finalize
!!
!>\BRIEF        Write variables to restart file
!!
!! DESCRIPTION  : Write variables to restart file
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S): 
!!
!! REFERENCES	: 
!!
!! \n
!_ ================================================================================================================================

  SUBROUTINE stomate_finalize (kjit, kjpindex, index, clay, assim_param, silt, bulk) 
    
    IMPLICIT NONE
    
    !! 0. Variable and parameter declaration
    !! 0.1 Input variables
    INTEGER,INTENT(in)                       :: kjit              !! Time step number (unitless)
    INTEGER,INTENT(in)                       :: kjpindex          !! Domain size - terrestrial pixels only (unitless)
    INTEGER,DIMENSION(kjpindex),INTENT(in)   :: index             !! Indices of the terrestrial pixels only (unitless)
    REAL,DIMENSION(kjpindex),INTENT(in)      :: clay              !! Clay fraction of soil (0-1, unitless)
    REAL,DIMENSION(kjpindex),INTENT(in)      :: silt              !! Silt fraction of soil (0-1, unitless) 
    REAL,DIMENSION(kjpindex),INTENT(in)      :: bulk              !! Bulk density (kg/m**3) 
    REAL,DIMENSION(kjpindex,nvm,npco2),INTENT(in) :: assim_param    !! leaf nutrient derive vcmax,jmax and leafN,nue
                                                                            !! photosynthesis  
    !! 0.4 Local variables
    REAL                                   :: dt_days_read             !! STOMATE time step read in restart file (days)
    INTEGER                                :: l,k,ji, jv, i, j, m      !! indices    
    REAL,PARAMETER                         :: max_dt_days = 5.         !! Maximum STOMATE time step (days)
    REAL                                   :: hist_days                !! Writing frequency for history file (days)
    REAL,DIMENSION(0:nslm)                 :: z_soil                   !! Variable to store depth of the different soil layers (m)
    REAL,DIMENSION(kjpindex)               :: cvegtot                  !! Total "vegetation" cover (unitless)
    REAL,DIMENSION(kjpindex)               :: precip                   !! Total liquid and solid precipitation  
                                                                              !! @tex $(??mm dt_stomate^{-1})$ @endtex 
    REAL,DIMENSION(kjpindex,nvm)           :: gpp_d                    !! Gross primary productivity per ground area 
                                                                              !! @tex $(??gC m^{-2} dt_stomate^{-1})$ @endtex  
    REAL,DIMENSION(kjpindex,nvm)           :: gpp_daily_x              !! "Daily" gpp for teststomate  
                                                                              !! @tex $(??gC m^{-2} dt_stomate^{-1})$ @endtex 
    REAL,DIMENSION(kjpindex,nvm)           :: resp_hetero_litter       !! Litter heterotrophic respiration per ground area 
                                                                              !! @tex $(gC m^{-2} day^{-1})$ @endtex  
                                                                              !! ??Same variable is also used to 
                                                                              !! store heterotrophic respiration per ground area 
                                                                              !! over ::dt_sechiba?? 
    REAL,DIMENSION(kjpindex,nvm)           :: resp_hetero_soil         !! soil heterotrophic respiration  
                                                                              !! @tex $(gC m^{-2} day^{-1})$ @endtex
    REAL,DIMENSION(kjpindex,nvm)           :: veget_cov                !! Fractional coverage: actually share of the pixel 
                                                                              !! covered by a PFT (fraction of ground area), 
                                                                              !! taking into account LAI ??(= grid scale fpc)?? 
!DSG: I think this variable is not needed:
    REAL,DIMENSION(kjpindex,nvm,2)           :: vcmax                    !! Maximum rate of carboxylation
                                                                              !! @tex $(\mumol m^{-2} s^{-1})$ @endtex
!DSG
    REAL,DIMENSION(kjpindex,nlevs)         :: control_moist_inst       !! Moisture control of heterotrophic respiration 
                                                                              !! (0-1, unitless) 
    REAL,DIMENSION(kjpindex,nlevs)         :: control_temp_inst        !! Temperature control of heterotrophic 
                                                                              !! respiration, above and below (0-1, unitless) 
    REAL,DIMENSION(kjpindex,ncarb,nvm)     :: som_input_inst    !! Quantity of carbon going into carbon pools from 
                                                                              !! litter decomposition 
                                                                              !! @tex $(gC m^{-2} day^{-1})$ @endtex 
    
    INTEGER                                :: ier                      !! Check errors in netcdf call (unitless)
    REAL                                   :: sf_time                  !! Intermediate variable to calculate current time 
                                                                              !! step 
    INTEGER                                :: max_totsize              !! Memory management - maximum memory size (Mb)
    INTEGER                                :: totsize_1step            !! Memory management - memory required to store one 
                                                                              !! time step on one processor (Mb) 
    INTEGER                                :: totsize_tmp              !! Memory management - memory required to store one 
                                                                              !! time step on all processors(Mb) 
    REAL                                   :: xn                       !! How many times have we treated in this forcing 
    REAL, DIMENSION(kjpindex)              :: vartmp                   !! Temporary variable
    INTEGER                                :: vid                      !! Variable identifer of netCDF (unitless)
    INTEGER                                :: nneigh                   !! Number of neighbouring pixels
    INTEGER                                :: direct                   !! ??
    INTEGER,DIMENSION(ndm)                 :: d_id                     !! ??
    REAL,DIMENSION(nbp_glo)                :: clay_g                   !! Clay fraction of soil (0-1, unitless), parallel 
                                                                              !! computing 
    REAL,DIMENSION(nbp_glo)                :: silt_g                   !! Silt fraction of soil (0-1, unitless), parallel  
                                                                              !! computing 
    REAL,DIMENSION(nbp_glo)                :: bulk_g                   !! Bulk density (kg/m**3), parallel                                                                                                                             !! computing 
    REAL,ALLOCATABLE,DIMENSION(:,:,:,:)    :: carbon_input_g       !! Quantity of carbon going into carbon pools from 
                                                                              !! litter decomposition  
                                                                              !! @tex $(gC m^{-2} dt_sechiba^{-1})$ @endtex, parallel 
                                                                              !! computing 
    REAL,ALLOCATABLE,DIMENSION(:,:,:,:)    :: nitrogen_input_g     !! Quantity of nitrogen going into nitrogen pools from  
                                                                          !! litter decomposition   
                                                                          !! @tex $(gC m^{-2} dtradia^{-1})$ @endtex, parallel  
    REAL,ALLOCATABLE,DIMENSION(:,:,:)      :: control_moist_g          !! Moisture control of heterotrophic respiration 
                                                                              !! (0-1, unitless), parallel computing 
    REAL,ALLOCATABLE,DIMENSION(:,:,:)      :: control_temp_g           !! Temperature control of heterotrophic respiration 
                                                                              !! (0-1, unitless), parallel computing 
    REAL,ALLOCATABLE,DIMENSION(:,:)        :: npp_equil_g              !! Equilibrium NPP written to forcesoil 
                                                                              !! @tex $(gC m^{-2} year^{-1})$ @endtex, parallel 
                                                                              !! computing 

    REAL                                   :: net_cflux_prod_monthly_sum    !! AR5 output?? gC m2 month-1 (one variable for 
                                                                                   !! reduce_sum and one for bcast??), parallel 
                                                                                   !! computing 
    REAL                                   :: net_cflux_prod_monthly_tot    !! AR5 output?? gC m2 month-1 (one variable for 
                                                                                   !! reduce_sum and one for bcast??), parallel 
                                                                                   !! computing 
    REAL                                   :: net_harvest_above_monthly_sum !! AR5 output?? gC m2 month-1 (one variable for 
                                                                                   !! reduce_sum and one for bcast??), parallel 
                                                                                   !! computing 
    REAL                                   :: net_harvest_above_monthly_tot !! AR5 output?? gC m2 month-1 (one variable for 
                                                                                   !! reduce_sum and one for bcast??), parallel 
                                                                                   !! computing 
    REAL                                   :: net_biosp_prod_monthly_sum    !! AR5 output?? gC m2 month-1 (one variable for 
                                                                                   !! reduce_sum and one for bcast??), parallel 
                                                                                   !! computing 
    REAL                                   :: net_biosp_prod_monthly_tot    !! AR5 output?? gC m2 month-1 (one variable for 
                                                                                   !! reduce_sum and one for bcast??), parallel 
                                                                                   !! computing 
    REAL, DIMENSION(kjpindex,nvm,nbpools)  :: carbon_stock                  !! Array containing the carbon stock for each pool
                                                                                   !! used by ORCHIDEE

!_ ================================================================================================================================
    
    !! 1. Write restart file for stomate
    IF (printlev>=3) WRITE (numout,*) 'Write restart file for STOMATE'
       
    CALL writerestart &
         (kjpindex, index, &
         dt_days, days_since_beg, &
         ind, adapted, regenerate, &
         humrel_daily, max_eau_var_daily, gdd_init_date, litterhum_daily, &
         t2m_daily, t2m_min_daily, tsurf_daily, tsoil_daily,   &
         soilhum_daily, &
         drainage_longterm, runoff_longterm,                     &
         precip_daily, &
         gpp_daily, npp_daily, turnover_daily, &
         humrel_month, humrel_week, moiavail_growingseason, &
         t2m_longterm, tau_longterm, t2m_month, t2m_week, &
         tsoil_month, soilhum_month, fireindex, firelitter, &
         maxhumrel_lastyear, maxhumrel_thisyear, &
         minhumrel_lastyear, minhumrel_thisyear, &
         maxgppweek_lastyear, maxgppweek_thisyear, &
         gdd0_lastyear, gdd0_thisyear, &
         precip_lastyear, precip_thisyear, &
         gdd_m5_dormance, gdd_from_growthinit, gdd_midwinter, ncd_dormance, ngd_minus5, &
         PFTpresent, npp_longterm, lm_lastyearmax, lm_thisyearmax, &
         maxfpc_lastyear, maxfpc_thisyear, &
         turnover_longterm, gpp_week, npp_week, biomass, resp_maint_part, &
         leaf_age, leaf_frac, &
         senescence, when_growthinit, age, &
         resp_hetero_d, resp_maint_d, resp_growth_d, co2_fire, co2_to_bm_dgvm, &
         n_to_bm, p_to_bm, veget_lastlight, everywhere, need_adjacent, &
         RIP_time, &
         time_hum_min, hum_min_dormance, &
         litter, dead_leaves, &
         som, lignin_struc, lignin_wood,turnover_time,&
         prod10,prod100,flux10, flux100, &
         convflux, cflux_prod10, cflux_prod100, &
         prod10_harvest,prod100_harvest,flux10_harvest, flux100_harvest, &
         convflux_harvest, cflux_prod10_harvest, cflux_prod100_harvest, &
         woodharvestpft, bm_to_litter, carb_mass_total, &
         Tseason, Tseason_length, Tseason_tmp, &
         Tmin_spring_time, begin_leaves, onset_date, &
         global_years, ok_equilibrium, nbp_accu, nbp_flux, &
         MatrixV, VectorU, previous_stock, current_stock, assim_param, CN_som_litter_longterm, tau_CN_longterm, KF, k_latosa_adapt, &
         CP_som_litter_longterm, np_leaf_avg_season, soil_p_min,&
         f_Pdissolved, EW_grain,&
         rue_longterm,cn_leaf_avg_season,nstress_season, pstress_season, &
         lai_target_longterm, &
         soil_n_min,p_O2,bact, &
         record_reserve_grass, &
         wshtotsum, sr_ugb, sla_calc, nb_ani, grazed_frac, &
         import_yield, t2m_14, litter_not_avail, nb_grazingdays, &
         after_snow, after_wet, wet1day, wet2day, GRM_devstage, &
         days_senescence)
    
    !! 2. Write file with variables that force general processes in stomate
    IF (ok_co2 .AND. allow_forcing_write ) THEN
       IF ( TRIM(forcing_name) /= 'NONE' ) THEN  
          CALL forcing_write(forcing_id,1,iisf)
          ! Close forcing file
          IF (is_root_prc) ier = NF90_CLOSE (forcing_id)
          forcing_id=-1
       END IF
    END IF
    
    !! 3. Collect variables that force the soil processes in stomate
    IF (TRIM(Cforcing_name) /= 'NONE' ) THEN 
       
       !! Collet variables 
       IF (printlev >= 1) WRITE(numout,*) 'stomate: writing the forcing file for carbon spinup'
       DO iatt = 1, nparan*nbyear
          IF ( nforce(iatt) > 0 ) THEN
             carbon_input(:,:,:,iatt) = &
                  & carbon_input(:,:,:,iatt)/REAL(nforce(iatt),r_std)
             nitrogen_input(:,:,:,iatt) = &
                  & nitrogen_input(:,:,:,iatt)/REAL(nforce(iatt),r_std)
             control_moist(:,:,iatt) = &
                  & control_moist(:,:,iatt)/REAL(nforce(iatt),r_std)
             control_temp(:,:,iatt) = &
                  & control_temp(:,:,iatt)/REAL(nforce(iatt),r_std)
             npp_equil(:,iatt) = &
                  & npp_equil(:,iatt)/REAL(nforce(iatt),r_std)
             drainage(:,:,iatt) = & 
                  &  drainage(:,:,iatt)/REAL(nforce(iatt),r_std) 
             runoff(:,:,iatt) = & 
                  &  runoff(:,:,iatt)/REAL(nforce(iatt),r_std) 
          ELSE
             IF (printlev >= 1) THEN
                WRITE(numout,*) 'We have no soil carbon forcing data for this time step:', iatt
                WRITE(numout,*) ' -> we set them to zero'
             END IF
             carbon_input(:,:,:,iatt) = zero
             nitrogen_input(:,:,:,iatt) = zero
             control_moist(:,:,iatt) = zero
             control_temp(:,:,iatt) = zero
             npp_equil(:,iatt) = zero
          ENDIF
       ENDDO
       
       ! Allocate memory for parallel computing
       IF (is_root_prc) THEN
          ALLOCATE(carbon_input_g(nbp_glo,ncarb,nvm,nparan*nbyear))
          ALLOCATE(nitrogen_input_g(nbp_glo,ncarb,nvm,nparan*nbyear))
          ALLOCATE(control_moist_g(nbp_glo,nlevs,nparan*nbyear))
          ALLOCATE(control_temp_g(nbp_glo,nlevs,nparan*nbyear))
          ALLOCATE(npp_equil_g(nbp_glo,nparan*nbyear))
       ENDIF
       
       ! Gather distributed variables

       CALL gather(clay,clay_g)
       CALL gather(silt,silt_g) 
       CALL gather(bulk,bulk_g) 
       CALL gather(carbon_input,carbon_input_g)
       CALL gather(nitrogen_input,nitrogen_input_g)
       CALL gather(control_moist,control_moist_g)
       CALL gather(control_temp,control_temp_g)
       CALL gather(npp_equil,npp_equil_g)
       
       !! Create netcdf
       ! Create, define and populate a netcdf file containing the forcing data.
       ! For the root processor only (parallel computing). NF90_ are functions
       ! from and external library.  
       IF (is_root_prc) THEN
          IF (printlev>=2) WRITE (numout,*) 'Create Cforcing file : ',TRIM(Cforcing_name)
          ! Create new netCDF dataset
          ier = NF90_CREATE (TRIM(Cforcing_name),NF90_64BIT_OFFSET ,Cforcing_id)
          IF (ier /= NF90_NOERR) THEN
             WRITE (numout,*) 'Error in creating Cforcing file : ',TRIM(Cforcing_name)
             CALL ipslerr_p (3,'stomate_finalize', &
                  &        'PROBLEM creating Cforcing file', &
                  &        NF90_STRERROR(ier),'')
          END IF
          
          ! Add variable attribute
          ! Note ::nbp_glo is the number of global continental points
          ier = NF90_PUT_ATT (Cforcing_id,NF90_GLOBAL, &
               &                        'kjpindex',REAL(nbp_glo,r_std))
          ier = NF90_PUT_ATT (Cforcing_id,NF90_GLOBAL, &
               &                        'nparan',REAL(nparan,r_std))
          ier = NF90_PUT_ATT (Cforcing_id,NF90_GLOBAL, &
               &                        'nbyear',REAL(nbyear,r_std))
          
          ! Add new dimension
          ier = NF90_DEF_DIM (Cforcing_id,'points',nbp_glo,d_id(1))
          ier = NF90_DEF_DIM (Cforcing_id,'carbtype',ncarb,d_id(2))
          ier = NF90_DEF_DIM (Cforcing_id,'vegtype',nvm,d_id(3))
          ier = NF90_DEF_DIM (Cforcing_id,'level',nlevs,d_id(4))
          ier = NF90_DEF_DIM (Cforcing_id,'time_step',NF90_UNLIMITED,d_id(5))
          
          ! Add new variable
          ier = NF90_DEF_VAR (Cforcing_id,'points',    r_typ,d_id(1),vid)
          ier = NF90_DEF_VAR (Cforcing_id,'carbtype',  r_typ,d_id(2),vid)
          ier = NF90_DEF_VAR (Cforcing_id,'vegtype',   r_typ,d_id(3),vid)
          ier = NF90_DEF_VAR (Cforcing_id,'level',     r_typ,d_id(4),vid)
          ier = NF90_DEF_VAR (Cforcing_id,'time_step', r_typ,d_id(5),vid)
          ier = NF90_DEF_VAR (Cforcing_id,'index',     r_typ,d_id(1),vid)
          ier = NF90_DEF_VAR (Cforcing_id,'clay',      r_typ,d_id(1),vid)
          ier = NF90_DEF_VAR (Cforcing_id,'silt',      r_typ,d_id(1),vid) 
          ier = NF90_DEF_VAR (Cforcing_id,'bulk',      r_typ,d_id(1),vid) 
          ier = NF90_DEF_VAR (Cforcing_id,'carbon_input',r_typ, &
               &                        (/ d_id(1),d_id(2),d_id(3),d_id(5) /),vid)
          ier = NF90_DEF_VAR (Cforcing_id,'nitrogen_input',r_typ, &
               &                        (/ d_id(1),d_id(2),d_id(3),d_id(5) /),vid)
          ier = NF90_DEF_VAR (Cforcing_id,'control_moist',r_typ, &
               &                        (/ d_id(1),d_id(4),d_id(5) /),vid)
          ier = NF90_DEF_VAR (Cforcing_id,'control_temp',r_typ, &
               &                        (/ d_id(1),d_id(4),d_id(5) /),vid)
          ier = NF90_DEF_VAR (Cforcing_id,'npp_equil',r_typ, &
               &                        (/ d_id(1),d_id(5) /),vid)

          ier = NF90_DEF_VAR (Cforcing_id,'drainage',r_typ, & 
               &                        (/ d_id(1),d_id(3),d_id(6) /),vid) 
          ier = NF90_DEF_VAR (Cforcing_id,'runoff',r_typ, & 
               &                        (/ d_id(1),d_id(3),d_id(6) /),vid) 
          ier = NF90_ENDDEF (Cforcing_id)
          
          ! Given the name of a varaible, nf90_inq_varid finds the variable 
          ! ID (::vid). Put data value(s) into variable ::vid 
          ier = NF90_INQ_VARID (Cforcing_id,'points',vid)
          ier = NF90_PUT_VAR (Cforcing_id,vid, &
               &                          (/(REAL(i,r_std),i=1,nbp_glo)/))
          ier = NF90_INQ_VARID (Cforcing_id,'carbtype',vid)
          ier = NF90_PUT_VAR (Cforcing_id,vid, &
               &                        (/(REAL(i,r_std),i=1,ncarb)/))
          ier = NF90_INQ_VARID (Cforcing_id,'vegtype',vid)
          ier = NF90_PUT_VAR (Cforcing_id,vid, &
               &                            (/(REAL(i,r_std),i=1,nvm)/))
          ier = NF90_INQ_VARID (Cforcing_id,'level',vid)
          ier = NF90_PUT_VAR (Cforcing_id,vid, &
               &                          (/(REAL(i,r_std),i=1,nlevs)/))
          ier = NF90_INQ_VARID (Cforcing_id,'time_step',vid)
          ier = NF90_PUT_VAR (Cforcing_id,vid, &
               &                          (/(REAL(i,r_std),i=1,nparan*nbyear)/))
          ier = NF90_INQ_VARID (Cforcing_id,'index',vid)
          ier = NF90_PUT_VAR (Cforcing_id,vid, REAL(index_g,r_std) )
          ier = NF90_INQ_VARID (Cforcing_id,'clay',vid)
          ier = NF90_PUT_VAR (Cforcing_id,vid, clay_g )
          ier = NF90_INQ_VARID (Cforcing_id,'silt',vid) 
          ier = NF90_PUT_VAR (Cforcing_id,vid, silt_g ) 
          ier = NF90_INQ_VARID (Cforcing_id,'bulk',vid) 
          ier = NF90_PUT_VAR (Cforcing_id,vid, bulk_g ) 
          ier = NF90_INQ_VARID (Cforcing_id,'carbon_input',vid)
          ier = NF90_PUT_VAR (Cforcing_id,vid, carbon_input_g )
          ier = NF90_INQ_VARID (Cforcing_id,'nitrogen_input',vid)
          ier = NF90_PUT_VAR (Cforcing_id,vid, nitrogen_input_g )
          ier = NF90_INQ_VARID (Cforcing_id,'control_moist',vid)
          ier = NF90_PUT_VAR (Cforcing_id,vid, control_moist_g )
          ier = NF90_INQ_VARID (Cforcing_id,'control_temp',vid)
          ier = NF90_PUT_VAR (Cforcing_id,vid, control_temp_g )
          ier = NF90_INQ_VARID (Cforcing_id,'npp_equil',vid)
          ier = NF90_PUT_VAR (Cforcing_id,vid, npp_equil_g )
          
          ! Close netCDF
          ier = NF90_CLOSE (Cforcing_id)
          IF (ier /= NF90_NOERR) THEN
             CALL ipslerr_p (3,'stomate_finalize', &
                  &        'PROBLEM in closing Cforcing file', &
                  &        NF90_STRERROR(ier),'')
          END IF
          
          Cforcing_id = -1
       ENDIF

       ! Clear memory
       IF (is_root_prc) THEN
          DEALLOCATE(carbon_input_g)
          DEALLOCATE(nitrogen_input_g)
          DEALLOCATE(control_moist_g)
          DEALLOCATE(control_temp_g)
          DEALLOCATE(npp_equil_g)
       ENDIF
       
    ENDIF
  
  END SUBROUTINE stomate_finalize


!! ================================================================================================================================
!! SUBROUTINE 	: stomate_init
!!
!>\BRIEF        The routine is called only at the first simulation. At that 
!! time settings and flags are read and checked for internal consistency and 
!! memory is allocated for the variables in stomate.
!!
!! DESCRIPTION  : The routine reads the 
!! following flags from the run definition file:
!! -ipd (index of grid point for online diagnostics)\n
!! -ok_herbivores (flag to activate herbivores)\n
!! -treat_expansion (flag to activate PFT expansion across a pixel\n
!! -harvest_agri (flag to harvest aboveground biomass from agricultural PFTs)\n
!! \n
!! Check for inconsistent setting between the following flags:
!! -ok_stomate\n
!! -ok_dgvm\n
!! -ok_co2\n
!! \n
!! Memory is allocated for all the variables of stomate and new indexing tables 
!! are build. New indexing tables are needed because a single pixel can conatin 
!! several PFTs. The new indexing tables have separate indices for the different 
!! PFTs. Similar index tables are build for land use cover change.\n
!! \n
!! Several global variables and land cover change variables are initialized to 
!! zero.\n
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S): Strictly speaking the subroutine has no output 
!! variables. However, the routine allocates memory and builds new indexing 
!! variables for later use.\n 
!!
!! REFERENCE(S)	: None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE stomate_init &
       &  (kjpij, kjpindex, index, lalo, &
       &   rest_id_stom, hist_id_stom, hist_id_stom_IPCC)

  !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER,INTENT(in)                    :: kjpij             !! Total size of the un-compressed grid, including 
                                                                      !! oceans (unitless) 
    INTEGER,INTENT(in)                    :: kjpindex          !! Domain size - number of terrestrial pixels 
                                                                      !! (unitless) 
    INTEGER,INTENT(in)                    :: rest_id_stom      !! STOMATE's _Restart_ file identifier
    INTEGER,INTENT(in)                    :: hist_id_stom      !! STOMATE's _history_ file identifier
    INTEGER,INTENT(in)                    :: hist_id_stom_IPCC !! STOMATE's IPCC _history_ file identifier 
    INTEGER,DIMENSION(kjpindex),INTENT(in):: index             !! Indices of the terrestrial pixels on the global 
                                                                      !! map 
    REAL,DIMENSION(kjpindex,2),INTENT(in) :: lalo              !! Geogr. coordinates (latitude,longitude) (degrees)
   
    !! 0.2 Output variables

    !! 0.3 Modified variables

    !! 0.4 Local variables

    LOGICAL                                      :: l_error           !! Check errors in netcdf call
    INTEGER                               :: ier               !! Check errors in netcdf call
    INTEGER                               :: ji,j,ipd,l        !! Indices
!_ ================================================================================================================================
    
  !! 1. Online diagnostics

    IF ( kjpindex > 0 ) THEN
       !Config  Key  = STOMATE_DIAGPT
       !Config  Desc = Index of grid point for online diagnostics
       !Config If    = OK_STOMATE
       !Config  Def  = 1
       !Config  Help = This is the index of the grid point which
       !               will be used for online diagnostics.
       !Config Units = [-]
       ! By default ::ipd is set to 1
       ipd = 1
       ! Get ::ipd from run definition file
       CALL getin_p('STOMATE_DIAGPT',ipd)
       ipd = MIN( ipd, kjpindex )
       IF ( printlev >=3 ) THEN
          WRITE(numout,*) 'Stomate: '
          WRITE(numout,*) '  Index of grid point for online diagnostics: ',ipd
          WRITE(numout,*) '  Lon, lat:',lalo(ipd,2),lalo(ipd,1)
          WRITE(numout,*) '  Index of this point on GCM grid: ',index(ipd)
       END IF
    ENDIF
    
  !! 2. Check consistency of flags

    IF ( ( .NOT. ok_stomate ) .AND. ok_dgvm ) THEN
       WRITE(numout,*) 'Cannot do dynamical vegetation without STOMATE.'
       WRITE(numout,*) 'Inconsistency between ::ok_stomate and ::ok_dgvm'
       WRITE(numout,*) 'Stop: fatal error'
       STOP
    ENDIF

    IF ((.NOT.ok_co2).AND.ok_stomate) THEN
       WRITE(numout,*) 'Cannot call STOMATE without GPP.'
       WRITE(numout,*) 'Inconsistency between ::ok_stomate and ::ok_co2'
       WRITE(numout,*) 'Stop: fatal error'
       STOP
    ENDIF

  !! 3. Communicate settings
    
    IF (printlev >=2) THEN
       WRITE(numout,*) 'stomate first call - overview of the activated flags:'
       WRITE(numout,*) '  Photosynthesis: ', ok_co2
       WRITE(numout,*) '  STOMATE: ', ok_stomate
       WRITE(numout,*) '  LPJ: ', ok_dgvm
    END IF
  !! 4. Allocate memory for STOMATE's variables

    l_error = .FALSE.

    ALLOCATE(veget_cov_max(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for veget_cov_max. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(ind(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for ind. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(adapted(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for adapted. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(regenerate(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for regenerate. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(humrel_daily(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for humrel_daily. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(max_eau_var_daily(kjpindex),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for max_eau_var_daily. We stop. We need kjpindex words',kjpindex
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(litterhum_daily(kjpindex),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for litterhum_daily. We stop. We need kjpindex words',kjpindex
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(t2m_daily(kjpindex),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for t2m_daily. We stop. We need kjpindex words',kjpindex
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(t2m_min_daily(kjpindex),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for t2m_min_daily. We stop. We need kjpindex words',kjpindex
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(tsurf_daily(kjpindex),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for tsurf_daily. We stop. We need kjpindex words',kjpindex
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(tsoil_daily(kjpindex,nslm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for tsoil_daily. We stop. We need kjpindex*nslm words',kjpindex,nslm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(soilhum_daily(kjpindex,nslm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for soilhum_daily. We stop. We need kjpindex*nslm words',kjpindex,nslm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(precip_daily(kjpindex),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for precip_daily. We stop. We need kjpindex words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(gpp_daily(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for gpp_daily. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(npp_daily(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for npp_daily. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(turnover_daily(kjpindex,nvm,nparts,nelements),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for turnover_daily. We stop. We need kjpindex*nvm*nparts*nelements words', &
       &   kjpindex,nvm,nparts,nelements
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(turnover_littercalc(kjpindex,nvm,nparts,nelements),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for turnover_littercalc. We stop. We need kjpindex*nvm*nparts*nelements words', & 
        &  kjpindex,nvm,nparts,nelements
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(humrel_month(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for humrel_month. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(humrel_week(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for humrel_week. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(moiavail_growingseason(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for moiavail_growingseason. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(t2m_longterm(kjpindex),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for t2m_longterm. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(t2m_month(kjpindex),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for t2m_month. We stop. We need kjpindex words',kjpindex
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(Tseason(kjpindex),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for Tseason. We stop. We need kjpindex words',kjpindex
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(Tseason_length(kjpindex),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for Tseason_length. We stop. We need kjpindex words',kjpindex
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(Tseason_tmp(kjpindex),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for Tseason_tmp. We stop. We need kjpindex words',kjpindex
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(Tmin_spring_time(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for Tmin_spring_time. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(onset_date(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for onset_date. We stop. We need kjpindex*nvm*nparts words',kjpindex,nvm,2
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(t2m_week(kjpindex),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for t2m_week. We stop. We need kjpindex words',kjpindex
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(tsoil_month(kjpindex,nslm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for tsoil_month. We stop. We need kjpindex*nslm words',kjpindex,nslm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(soilhum_month(kjpindex,nslm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for soilhum_month. We stop. We need kjpindex*nslm words',kjpindex,nslm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(fireindex(kjpindex,nvm),stat=ier) 
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for fireindex. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(firelitter(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for firelitter. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(maxhumrel_lastyear(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for maxhumrel_lastyear. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(maxhumrel_thisyear(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for maxhumrel_thisyear. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(minhumrel_lastyear(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for minhumrel_lastyear. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(minhumrel_thisyear(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for minhumrel_thisyear. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(maxgppweek_lastyear(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for maxgppweek_lastyear. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(maxgppweek_thisyear(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for maxgppweek_thisyear. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(gdd0_lastyear(kjpindex),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for gdd0_lastyear. We stop. We need kjpindex words',kjpindex
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(gdd0_thisyear(kjpindex),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for gdd0_thisyear. We stop. We need kjpindex words',kjpindex
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(gdd_init_date(kjpindex,2),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for gdd_init_date. We stop. We need kjpindex*2 words',kjpindex,2
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(gdd_from_growthinit(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for gdd_from_growthinit. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(precip_lastyear(kjpindex),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for precip_lastyear. We stop. We need kjpindex*nvm words',kjpindex
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(precip_thisyear(kjpindex),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for precip_thisyear. We stop. We need kjpindex words',kjpindex
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(gdd_m5_dormance(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for gdd_m5_dormance. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(gdd_midwinter(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for gdd_midwinter. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(ncd_dormance(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for ncd_dormance. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(ngd_minus5(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for ngd_minus5. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(PFTpresent(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for PFTpresent. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(npp_longterm(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for npp_longterm. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(lm_lastyearmax(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for lm_lastyearmax. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(lm_thisyearmax(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for lm_thisyearmax. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(maxfpc_lastyear(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for maxfpc_lastyear. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(maxfpc_thisyear(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for maxfpc_thisyear. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(turnover_longterm(kjpindex,nvm,nparts,nelements),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for turnover_longterm. We stop. We need kjpindex*nvm*nparts*nelements words', & 
       &    kjpindex,nvm,nparts,nelements
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(gpp_week(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for gpp_week. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(npp_week(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for npp_week. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF


    ALLOCATE(biomass(kjpindex,nvm,nparts,nelements),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for biomass. We stop. We need kjpindex*nvm*nparts*nelements words', &
       &    kjpindex,nvm,nparts,nelements
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(senescence(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for senescence. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(begin_leaves(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for begin_leaves. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(when_growthinit(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for when_growthinit. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(age(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for age. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(resp_hetero_d(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for resp_hetero_d. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(resp_hetero_radia(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for resp_hetero_radia. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(resp_maint_d(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for resp_maint_d. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(resp_growth_d(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for resp_growth_d. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(co2_fire(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for co2_fire. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(co2_to_bm_dgvm(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for co2_to_bm_dgvm. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(n_to_bm(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for n_to_bm. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(p_to_bm(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for p_to_bm. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(veget_lastlight(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for veget_lastlight. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(everywhere(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for everywhere. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(need_adjacent(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for need_adjacent. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(leaf_age(kjpindex,nvm,nparts,nleafages),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for leaf_age. We stop. We need kjpindex*nvm*nleafages words', & 
       &      kjpindex,nvm,nparts,nleafages
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(leaf_frac(kjpindex,nvm,nparts,nleafages),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for leaf_frac. We stop. We need kjpindex*nvm*nleafages words', & 
       &      kjpindex,nvm,nparts,nleafages
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(RIP_time(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for RIP_time. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(time_hum_min(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for time_hum_min. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(hum_min_dormance(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for hum_min_dormance. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF


    ALLOCATE(litter(kjpindex,nlitt,nvm,nlevs,nelements),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for litter. We stop. We need kjpindex*nlitt*nvm*nlevs*nelements words', & 
       &    kjpindex,nlitt,nvm,nlevs,nelements
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(dead_leaves(kjpindex,nvm,nlitt),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for dead_leaves. We stop. We need kjpindex*nvm*nlitt words', & 
       &   kjpindex,nvm,nlitt
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(som(kjpindex,ncarb,nvm,nelements),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for som. We stop. We need kjpindex*ncarb*nvm*nelements words',kjpindex,ncarb,nvm,nelements
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(lignin_struc(kjpindex,nvm,nlevs),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for lignin_struc. We stop. We need kjpindex*nvm*nlevs words',kjpindex,nvm,nlevs
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(lignin_wood(kjpindex,nvm,nlevs),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for lignin_wood. We stop. We need kjpindex*nvm*nlevs words',kjpindex,nvm,nlevs
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(turnover_time(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for turnover_time. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(nep_daily(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for nep_daily. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(nep_monthly(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for nep_monthly. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE (cflux_prod_monthly(kjpindex), stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for cflux_prod_monthly. We stop. We need kjpindex words',kjpindex
       STOP 'stomate_init'
    ENDIF
 
    ALLOCATE (harvest_above_monthly(kjpindex), stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for harvest_above_monthly. We stop. We need kjpindex words',kjpindex
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(bm_to_litter(kjpindex,nvm,nparts,nelements),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for bm_to_litter. We stop. We need kjpindex*nvm*nparts*nelements words', & 
       &    kjpindex,nvm,nparts,nelements
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(bm_to_littercalc(kjpindex,nvm,nparts,nelements),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for bm_to_littercalc. We stop. We need kjpindex*nvm*nparts*nelements words', &
       &   kjpindex,nvm,nparts,nelements
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(herbivores(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for herbivores. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(hori_index(kjpindex),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for hori_index. We stop. We need kjpindex words',kjpindex
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(horipft_index(kjpindex*nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for horipft_index. We stop. We need kjpindex*nvm words',kjpindex*nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(resp_maint_part_radia(kjpindex,nvm,nparts),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for resp_maint_part_radia. We stop. We need kjpindex*nvm*nparts words', &
       &  kjpindex,nvm,nparts
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(resp_maint_radia(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for resp_maint_radia. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(resp_maint_part(kjpindex,nvm,nparts),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for resp_maint_part. We stop. We need kjpindex*nvm*nparts words', &
       &    kjpindex,nvm,nparts
       STOP 'stomate_init'
    ENDIF
    resp_maint_part(:,:,:) = zero

    ALLOCATE (horip10_index(kjpindex*10), stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for horip10_index. We stop. We need kjpindex*10 words',kjpindex,10
       STOP 'stomate_init'
    ENDIF

    ALLOCATE (horip100_index(kjpindex*100), stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for horip100_index. We stop. We need kjpindex*100 words',kjpindex,100
       STOP 'stomate_init'
    ENDIF

    ALLOCATE (horip11_index(kjpindex*11), stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for horip11_index. We stop. We need kjpindex*11 words',kjpindex,11
       STOP 'stomate_init'
    ENDIF

    ALLOCATE (horip101_index(kjpindex*101), stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for horip101_index. We stop. We need kjpindex*101 words',kjpindex,101
       STOP 'stomate_init'
    ENDIF

    ALLOCATE (prod10(kjpindex,0:10), stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for prod10. We stop. We need kjpindex*11 words',kjpindex,11
       STOP 'stomate_init'
    ENDIF

    ALLOCATE (prod100(kjpindex,0:100), stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for prod100. We stop. We need kjpindex*101 words',kjpindex,101
       STOP 'stomate_init'
    ENDIF

    ALLOCATE (flux10(kjpindex,10), stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for flux10. We stop. We need kjpindex*10 words',kjpindex,10
       STOP 'stomate_init'
    ENDIF

    ALLOCATE (flux100(kjpindex,100), stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for flux100. We stop. We need kjpindex*100 words',kjpindex,100
       STOP 'stomate_init'
    ENDIF

    ALLOCATE (convflux(kjpindex), stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for convflux. We stop. We need kjpindex words',kjpindex
       STOP 'stomate_init'
    ENDIF

    ALLOCATE (cflux_prod10(kjpindex), stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for cflux_prod10. We stop. We need kjpindex words',kjpindex
       STOP 'stomate_init'
    ENDIF

    ALLOCATE (cflux_prod100(kjpindex), stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for cflux_prod100. We stop. We need kjpindex words',kjpindex
       STOP 'stomate_init'
    ENDIF

    ALLOCATE (prod10_harvest(kjpindex,0:10), stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for prod10_harvest. We stop. We need kjpindex*11 words',kjpindex,11
       STOP 'stomate_init'
    ENDIF

    ALLOCATE (prod100_harvest(kjpindex,0:100), stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for prod100_harvest. We stop. We need kjpindex*101 words',kjpindex,101
       STOP 'stomate_init'
    ENDIF

    ALLOCATE (flux10_harvest(kjpindex,10), stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for flux10_harvest. We stop. We need kjpindex*10 words',kjpindex,10
       STOP 'stomate_init'
    ENDIF

    ALLOCATE (flux100_harvest(kjpindex,100), stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for flux100_harvest. We stop. We need kjpindex*100 words',kjpindex,100
       STOP 'stomate_init'
    ENDIF

    ALLOCATE (convflux_harvest(kjpindex), stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for convflux_harvest. We stop. We need kjpindex words',kjpindex
       STOP 'stomate_init'
    ENDIF

    ALLOCATE (cflux_prod10_harvest(kjpindex), stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for cflux_prod10_harvest. We stop. We need kjpindex words',kjpindex
       STOP 'stomate_init'
    ENDIF

    ALLOCATE (cflux_prod100_harvest(kjpindex), stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for cflux_prod100_harvest. We stop. We need kjpindex words',kjpindex
       STOP 'stomate_init'
    ENDIF

    ALLOCATE (woodharvestpft(kjpindex,nvm), stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for woodharvestpft. We stop. We need kjpindex*nvm words',kjpindex*nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE (harvest_above(kjpindex), stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for harvest_above. We stop. We need kjpindex words',kjpindex
       STOP 'stomate_init'
    ENDIF

    ALLOCATE (harvest_bioenergy(kjpindex,nvm), stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for harvest_bioenergy. We stop. We need kjpindex words',kjpindex*nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE (carb_mass_total(kjpindex), stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for carb_mass_total. We stop. We need kjpindex words',kjpindex
       STOP 'stomate_init'
    ENDIF

    ALLOCATE (som_input_daily(kjpindex,ncarb,nvm,nelements), stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for som_input_daily. We stop. We need kjpindex*ncarb*nvm*nelements words', & 
       &    kjpindex,ncarb,nvm,nelements
       STOP 'stomate_init'
    ENDIF

    ALLOCATE (control_temp_daily(kjpindex,nlevs), stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for control_temp_daily. We stop. We need kjpindex*nlevs words',kjpindex,nlevs
       STOP 'stomate_init'
    ENDIF

    ALLOCATE (control_moist_daily(kjpindex,nlevs), stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for control_moist_daily. We stop. We need kjpindex*nlevs words',kjpindex,nlevs
       STOP 'stomate_init'
    ENDIF

    ALLOCATE (fpc_max(kjpindex,nvm), stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for fpc_max. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(cn_leaf_avg_season(kjpindex,nvm),stat=ier) 
    l_error = l_error .OR. (ier /= 0) 
    IF (l_error) THEN 
       WRITE(numout,*) 'Memory allocation error for cn_leaf_avg_season. We stop. We need kjpindex*nvm words',kjpindex,nvm 
       STOP 'stomate_init' 
    ENDIF 
 
    ALLOCATE(nstress_season(kjpindex,nvm),stat=ier) 
    l_error = l_error .OR. (ier /= 0) 
    IF (l_error) THEN 
       WRITE(numout,*) 'Memory allocation error for nstress_season. We stop. We need kjpindex*nvm words',kjpindex,nvm 
       STOP 'stomate_init' 
    ENDIF 

    ALLOCATE(pstress_season(kjpindex,nvm),stat=ier) 
    l_error = l_error .OR. (ier /= 0) 
    IF (l_error) THEN 
       WRITE(numout,*) 'Memory allocation error for pstress_season. We stop. We need kjpindex*nvm words',kjpindex,nvm 
       STOP 'stomate_init' 
    ENDIF 

    ALLOCATE(lai_target_longterm(kjpindex,nvm),stat=ier) 
    l_error = l_error .OR. (ier /= 0) 
    IF (l_error) THEN 
       WRITE(numout,*) 'Memory allocation error for lai_target_longterm. We stop. We need kjpindex*nvm words',kjpindex,nvm 
       STOP 'stomate_init' 
    ENDIF 

    ALLOCATE(lai_target(kjpindex,nvm),stat=ier) 
    l_error = l_error .OR. (ier /= 0) 
    IF (l_error) THEN 
       WRITE(numout,*) 'Memory allocation error for lai_target. We stop. We need kjpindex*nvm words',kjpindex,nvm 
       STOP 'stomate_init' 
    ENDIF 
 
    ALLOCATE(soil_n_min(kjpindex,nvm,nnspec),stat=ier) 
    l_error = l_error .OR. (ier /= 0) 
    IF (l_error) THEN 
       WRITE(numout,*) 'Memory allocation error for soil_n_min. We stop. We need kjpindex*nvm words',kjpindex,nvm,nnspec 
       STOP 'stomate_init' 
    ENDIF 

    ALLOCATE(p_O2(kjpindex,nvm),stat=ier) 
    l_error = l_error .OR. (ier /= 0) 
    IF (l_error) THEN 
       WRITE(numout,*) 'Memory allocation error for p_O2. We stop. We need kjpindex*nvm words',kjpindex,nvm 
       STOP 'stomate_init' 
    ENDIF 

    ALLOCATE(bact(kjpindex,nvm),stat=ier) 
    l_error = l_error .OR. (ier /= 0) 
    IF (l_error) THEN 
       WRITE(numout,*) 'Memory allocation error for bact. We stop. We need kjpindex*nvm words',kjpindex,nvm 
       STOP 'stomate_init' 
    ENDIF 

    ALLOCATE(ok_equilibrium(kjpindex),stat=ier)
    l_error = l_error .OR. (ier /= 0) 
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for ok_equilibrium. We stop. We need kjpindex words',kjpindex
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(drainage_daily(kjpindex,nvm),stat=ier) 
    l_error = l_error .OR. (ier /= 0) 
    IF (l_error) THEN 
       WRITE(numout,*) ' Memory allocation error for drainage_daily. We stop. We need kjpindex words = ',kjpindex,nvm 
       STOP 'drainage_daily' 
    ENDIF 

    ALLOCATE(drainage_longterm(kjpindex),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for drainage_longterm. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'drainage_longterm'
    ENDIF

    ALLOCATE(runoff_daily(kjpindex,nvm),stat=ier) 
    l_error = l_error .OR. (ier /= 0) 
    IF (l_error) THEN 
       WRITE(numout,*) ' Memory allocation error for runoff_daily. We stop. We need kjpindex words = ',kjpindex,nvm 
       STOP 'runoff_daily' 
    ENDIF 

    ALLOCATE(runoff_longterm(kjpindex),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for runoff_longterm. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'runoff_longterm'
    ENDIF
 
    ALLOCATE (n_uptake_daily(kjpindex,nvm,nionspec), stat=ier) 
    l_error = l_error .OR. (ier.NE.0) 
    IF (l_error) THEN 
       WRITE(numout,*) ' Memory allocation error for n_uptake_daily. We stop. We need kjpindex words = ',kjpindex*nvm*nionspec 
       STOP 'n_uptake_daily' 
    ENDIF 
 
    ALLOCATE (n_mineralisation_d(kjpindex,nvm), stat=ier) 
    l_error = l_error .OR. (ier.NE.0) 
    IF (l_error) THEN 
       WRITE(numout,*) ' Memory allocation error for n_mineralisation_d. We stop. We need kjpindex words = ',kjpindex*nvm 
       STOP 'n_mineralisation_d' 
    ENDIF 

    ALLOCATE (p_uptake_daily(kjpindex,nvm), stat=ier) 
    l_error = l_error .OR. (ier.NE.0) 
    IF (l_error) THEN 
       WRITE(numout,*) ' Memory allocation error for p_uptake_daily. We stop. We need kjpindex words = ',kjpindex*nvm
       STOP 'p_uptake_daily' 
    ENDIF 

    ALLOCATE (p_leaching_daily(kjpindex,nvm), stat=ier) 
    l_error = l_error .OR. (ier.NE.0) 
    IF (l_error) THEN 
       WRITE(numout,*) ' Memory allocation error for p_leaching_daily. We stop. We need kjpindex words = ',kjpindex*nvm
       STOP 'p_leaching_daily' 
    ENDIF 

    ALLOCATE (p_ssorption_daily(kjpindex,nvm), stat=ier) 
    l_error = l_error .OR. (ier.NE.0) 
    IF (l_error) THEN 
       WRITE(numout,*) ' Memory allocation error for p_ssorption_daily. We stop. We need kjpindex words = ',kjpindex*nvm
       STOP 'p_ssorption_daily' 
    ENDIF 
 
    ALLOCATE (p_mineralisation_d(kjpindex,nvm), stat=ier) 
    l_error = l_error .OR. (ier.NE.0) 
    IF (l_error) THEN 
       WRITE(numout,*) ' Memory allocation error for p_mineralisation_d. We stop. We need kjpindex words = ',kjpindex*nvm 
       STOP 'p_mineralisation_d' 
    ENDIF 

    ALLOCATE (p_release_daily(kjpindex), stat=ier)
    l_error = l_error .OR. (ier.NE.0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for p_release_daily. We stop. We need kjpindex words = ',kjpindex
       STOP 'p_release_daily'
    ENDIF
    p_release_daily = zero

    ALLOCATE (EW_grain(kjpindex,nvm,nminerals,2), stat=ier)
    l_error = l_error .OR. (ier.NE.0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for EW_grain. We stop. We need kjpindex words = ',kjpindex*nvm*nminerals*2
       STOP 'EW_grain'
    ENDIF
    EW_grain(:,:,:,1) = zero
    EW_grain(:,:,:,2) = undef

    ALLOCATE (EW_flux(kjpindex,nvm,nECfluxes), stat=ier)
    l_error = l_error .OR. (ier.NE.0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for EW_flux. We stop. We need kjpindex words = ',kjpindex*nvm*nECfluxes
       STOP 'EW_flux'
    ENDIF
    EW_flux = zero

    ALLOCATE(carbon_eq(kjpindex),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for carbon_eq. We stop. We need kjpindex words',kjpindex
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(nbp_accu(kjpindex),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for nbp_accu. We stop. We need kjpindex words',kjpindex
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(nbp_flux(kjpindex),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for nbp_flux. We stop. We need kjpindex words',kjpindex
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(matrixA(kjpindex,nvm,nbpools,nbpools),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for matrixA. We stop. We need kjpindex*nvm*nbpools*nbpools words',  & 
       &     kjpindex, nvm, nbpools, nbpools
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(vectorB(kjpindex,nvm,nbpools),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for vectorB. We stop. We need kjpindex*nvm*nbpools words',  & 
       &     kjpindex, nvm, nbpools
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(VectorU(kjpindex,nvm,nbpools),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for VectorU. We stop. We need kjpindex*nvm*nbpools words',  & 
       &     kjpindex, nvm, nbpools
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(MatrixV(kjpindex,nvm,nbpools,nbpools),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for MatrixV. We stop. We need kjpindex*nvm*nbpools*nbpools words',  & 
       &     kjpindex, nvm, nbpools, nbpools
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(MatrixW(kjpindex,nvm,nbpools,nbpools),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for MatrixW. We stop. We need kjpindex*nvm*nbpools*nbpools words',  & 
       &     kjpindex, nvm, nbpools, nbpools
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(previous_stock(kjpindex,nvm,nbpools),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for previous_stock. We stop. We need kjpindex*nvm*nbpools words',  & 
       &     kjpindex, nvm, nbpools
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(current_stock(kjpindex,nvm,nbpools),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for current_stock. We stop. We need kjpindex*nvm*nbpools words',  & 
       &     kjpindex, nvm, nbpools
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(CN_som_litter_longterm(kjpindex,nvm,nbpools),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for CN_som_litter_longterm. We stop. We need kjpindex*nvm*nbpools words',  & 
       &     kjpindex, nvm, nbpools
       STOP 'stomate_init'
    ENDIF
    
    ALLOCATE(record_reserve_grass(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for record_reserve_grass. We stop. We need kjpindex*nvm words = ',kjpindex*nvm
       CALL ipslerr_p (3,'stomate_init', 'Memory allocation issue','','')
    ENDIF
    record_reserve_grass(:,:) = zero ! Is there a better place in the code for this?

    ALLOCATE(litter_avail(kjpindex,nlitt,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for litter_avail. We stop. We need kjpindex*nlitt*nvm words',  &
       &  kjpindex,nlitt,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(litter_not_avail(kjpindex,nlitt,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for litter_not_avail. We stop. We need kjpindex*nlitt*nvm words',  &
       &  kjpindex,nlitt,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(litter_avail_frac(kjpindex,nlitt,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for litter_avail_frac. We stop. We need kjpindex*nlitt*nvm words',  &
       &  kjpindex,nlitt,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(sla_calc(kjpindex,nvm),stat=ier)
   l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for sla_calc. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

   ALLOCATE(wshtotsum(kjpindex,nvm),stat=ier)
   l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for wshtotsum. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

   ALLOCATE(sr_ugb(kjpindex,nvm),stat=ier)
   l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for sr_ugb. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

   ALLOCATE(compt_ugb(kjpindex,nvm),stat=ier)
   l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for compt_ugb. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

   ALLOCATE(nb_ani(kjpindex,nvm),stat=ier)
   l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for nb_ani. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

   ALLOCATE(grazed_frac(kjpindex,nvm),stat=ier)
   l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for grazed_frac. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

   ALLOCATE(import_yield(kjpindex,nvm),stat=ier)
   l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for import_yield. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

   ALLOCATE(t2m_14(kjpindex),stat=ier)
   l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for t2m_14. We stop. We need kjpindex*nvm words',kjpindex
       STOP 'stomate_init'
    ENDIF

   ALLOCATE (sla_age1(kjpindex,nvm), stat=ier)
    l_error = l_error .OR. (ier.NE.0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for sla_age1. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

   ALLOCATE (when_growthinit_cut(kjpindex,nvm), stat=ier)
    l_error = l_error .OR. (ier.NE.0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for when_growthinit_cut. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF
   ALLOCATE(nb_grazingdays(kjpindex,nvm),stat=ier)
   l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for nb_grazingdays. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

   ALLOCATE(snowfall_daily(kjpindex),stat=ier)
   l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for snowfall_daily. We stop. We need kjpindex words',kjpindex
       STOP 'stomate_init'
    ENDIF

   ALLOCATE(snowmass_daily(kjpindex),stat=ier)
   l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for snowmass_daily. We stop. We need kjpindex words',kjpindex
       STOP 'stomate_init'
    ENDIF
! top 5 layer grassland soil moisture for grazing
    ALLOCATE(tmc_topgrass_daily(kjpindex),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for tmc_topgrass_daily. We stop. We need kjpindex words',kjpindex
       STOP 'stomate_init'
    ENDIF
    ALLOCATE(after_snow(kjpindex),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for after_snow. We stop. We need kjpindex words',kjpindex
       STOP 'stomate_init'
    ENDIF
    ALLOCATE(after_wet(kjpindex),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for after_wet. We stop. We need kjpindex words',kjpindex
       STOP 'stomate_init'
    ENDIF
    ALLOCATE(wet1day(kjpindex),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for wet1day. We stop. We need kjpindex words',kjpindex
       STOP 'stomate_init'
    ENDIF
    ALLOCATE(wet2day(kjpindex),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for wet2day. We stop. We need kjpindex words',kjpindex
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(grm_nfert(kjpindex,nvm,ninput),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for grm_nfert. We stop. We need kjpindex*nvm*ninput words',kjpindex,nvm,ninput
       STOP 'stomate_init'
    ENDIF
    ALLOCATE(grm_pfert(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for grm_pfert. We stop. We need kjpindex*nvm*ninput words',kjpindex,nvm,ninput
       STOP 'stomate_init'
    ENDIF
    ALLOCATE(BNF_clover(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for BNF_clover. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF
    ALLOCATE(BNF_clover_daily(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for BNF_clover_daily. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF
    ALLOCATE(GRM_devstage(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for GRM_devstage. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF
    ALLOCATE(p_input_dt(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for p_input_dt. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF
    ALLOCATE(agri_nfert(kjpindex,nvm,ninput),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for agri_nfert. We stop. We need kjpindex*nvm*ninput words',kjpindex,nvm,ninput
       STOP 'stomate_init'
    ENDIF
    ALLOCATE(agri_pfert(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for agri_pfert. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF
    ALLOCATE(rest_Ngrazing(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for rest_Ngrazing. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF
    ALLOCATE (grazed_above(kjpindex), stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for grazed_above. We stop. We need kjpindex words',kjpindex
       STOP 'stomate_init'
    ENDIF
    ALLOCATE(days_senescence(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for days_senescence. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF
    ALLOCATE(fclover(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for fclover. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(KF(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for KF. We stop. We need nvm words = ',kjpindex*nvm
       CALL ipslerr_p (3,'stomate_init', 'Memory allocation issue','','')
    ENDIF
    KF(:,:) = zero ! Is there a better place in the code for this?

    ALLOCATE(k_latosa_adapt(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for k_latosa_adapt. We stop. We need nvm words = ',kjpindex*nvm
       CALL ipslerr_p (3,'stomate_init', 'Memory allocation issue','','')
    ENDIF

    ALLOCATE (rue_longterm(kjpindex,nvm), stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for rue_longterm. We stop. We need kjpindex*nlevs words',kjpindex,nlevs
       CALL ipslerr_p (3,'stomate_init', 'Memory allocation issue','','')
    ENDIF
    rue_longterm(:,:) = un

! Phosphorus variables
    ALLOCATE(np_leaf_avg_season(kjpindex,nvm),stat=ier) 
    l_error = l_error .OR. (ier /= 0) 
    IF (l_error) THEN 
       WRITE(numout,*) 'Memory allocation error for np_leaf_avg_season. We stop. We need kjpindex*nvm words',kjpindex,nvm 
       STOP 'stomate_init' 
    ENDIF 

    ALLOCATE(CP_som_litter_longterm(kjpindex,nvm,nbpools),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for CP_som_litter_longterm. We stop. We need kjpindex*nvm*nbpools words',  & 
       &     kjpindex, nvm, nbpools
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(soil_p_min(kjpindex,nvm,npspec),stat=ier) 
    l_error = l_error .OR. (ier /= 0) 
    IF (l_error) THEN 
       WRITE(numout,*) 'Memory allocation error for soil_p_min. We stop. We need kjpindex*nvm words',kjpindex,nvm,npspec 
       STOP 'stomate_init' 
    ENDIF 

    ALLOCATE(f_Pdissolved(kjpindex,nvm),stat=ier) 
    l_error = l_error .OR. (ier /= 0) 
    IF (l_error) THEN 
       WRITE(numout,*) 'Memory allocation error for f_Pdissolved. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init' 
    ENDIF

  !! 5. File definitions

    ! Store history and restart files in common variables
    hist_id_stomate = hist_id_stom
    hist_id_stomate_IPCC = hist_id_stom_IPCC
    rest_id_stomate = rest_id_stom
    
    ! In STOMATE reduced grids are used containing only terrestrial pixels.
    ! Build a new indexing table for the vegetation fields separating 
    ! between the different PFTs. Note that ::index has dimension (kjpindex) 
    ! wheras ::indexpft has dimension (kjpindex*nvm). 

    hori_index(:) = index(:)

    DO j = 1, nvm
       DO ji = 1, kjpindex
          horipft_index((j-1)*kjpindex+ji) = index(ji)+(j-1)*kjpij + offset_omp - offset_mpi
       ENDDO
    ENDDO

    ! Similar index tables are build for the land cover change variables
    DO j = 1, 10
       DO ji = 1, kjpindex
          horip10_index((j-1)*kjpindex+ji) = index(ji)+(j-1)*kjpij + offset_omp - offset_mpi
       ENDDO
    ENDDO

    DO j = 1, 100
       DO ji = 1, kjpindex
          horip100_index((j-1)*kjpindex+ji) = index(ji)+(j-1)*kjpij + offset_omp - offset_mpi
       ENDDO
    ENDDO

    DO j = 1, 11
       DO ji = 1, kjpindex
          horip11_index((j-1)*kjpindex+ji) = index(ji)+(j-1)*kjpij + offset_omp - offset_mpi
       ENDDO
    ENDDO

    DO j = 1, 101
       DO ji = 1, kjpindex
          horip101_index((j-1)*kjpindex+ji) = index(ji)+(j-1)*kjpij + offset_omp - offset_mpi
       ENDDO
    ENDDO

  !! 6. Initialization of global and land cover change variables. 

    ! All variables are cumulative variables. bm_to_litter is not and is therefore
    ! excluded
    ! bm_to_litter(:,:,:) = zero
    turnover_daily(:,:,:,:) = zero
    resp_hetero_d(:,:) = zero
    nep_daily(:,:) = zero
    nep_monthly(:,:) = zero
    cflux_prod_monthly(:) = zero
    harvest_above_monthly(:) = zero
    control_moist_daily(:,:) = zero
    control_temp_daily(:,:) = zero
    som_input_daily(:,:,:,:) = zero
    drainage_daily(:,:) = zero 
    runoff_daily  (:,:) = zero 
    n_uptake_daily(:,:,:)=zero 
    p_uptake_daily(:,:)=zero 
    p_leaching_daily(:,:)=zero 
    p_ssorption_daily(:,:)=zero 
    n_mineralisation_d(:,:)=zero 
    p_mineralisation_d(:,:)=zero 

    ! Land cover change variables
    prod10(:,:)  = zero
    prod100(:,:) = zero
    flux10(:,:)  = zero
    flux100(:,:) = zero
    convflux(:)  = zero
    cflux_prod10(:) = zero
    cflux_prod100(:) = zero
    prod10_harvest(:,:)  = zero
    prod100_harvest(:,:) = zero
    flux10_harvest(:,:)  = zero
    flux100_harvest(:,:) = zero
    convflux_harvest(:)  = zero
    cflux_prod10_harvest(:) = zero
    cflux_prod100_harvest(:) = zero
    woodharvestpft(:,:) = zero
    fpc_max(:,:)=zero

   DO j = 1,nvm
     sla_calc(:,j) = sla(j)
     sla_age1(:,j) = sla_max(j)
   ENDDO
   wshtotsum(:,:)=zero
   sr_ugb(:,:)=zero
   compt_ugb(:,:)=zero
   nb_ani(:,:)=zero
   grazed_frac(:,:)=zero
   import_yield(:,:)=zero
   litter_not_avail(:,:,:) = 0.0
   litter_avail(:,:,:) = 0.0
   litter_avail_frac(:,:,:) = 1.0
   when_growthinit_cut(:,:)=20.0
   nb_grazingdays(:,:) = zero
   snowfall_daily(:) = zero
   snowmass_daily(:) = zero
   tmc_topgrass_daily(:) = zero
   after_snow(:) = zero
   after_wet(:) = zero
   wet1day(:) = 6
   wet2day(:) = 6
   grm_nfert(:,:,:) = zero
   grm_pfert(:,:) = zero
   BNF_clover(:,:) = zero
   BNF_clover_daily(:,:) = zero
   GRM_devstage(:,:) = 2.0
   p_input_dt(:,:) = zero
   agri_nfert(:,:,:) = zero
   agri_pfert(:,:) = zero
   rest_Ngrazing(:,:) = zero
   grazed_above(:) = zero
   days_senescence(:,:) = zero
   fclover(:,:) = zero
   ! It is now only assigned as universal value for site simulation
   ! it should be assigned a map for regional/global application
   DO j = 1, nvm
     IF (natural(j) .AND. .NOT. is_tree(j)) THEN
       fclover(:,j)=fix_legume_frac
     ENDIF
   ENDDO

   nstress_season(:,:) = zero 
   pstress_season(:,:) = zero 
   soil_n_min(:,:,:) = zero 
   soil_p_min(:,:,:) = zero 
   f_Pdissolved(:,:) = un
   lai_target(:,:) = zero
   lai_target_longterm(:,:) = zero

  END SUBROUTINE stomate_init


!! ================================================================================================================================
!! SUBROUTINE 	: stomate_clear
!!
!>\BRIEF        Deallocate memory of the stomate variables.
!!
!! DESCRIPTION  : None
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S): None
!!
!! REFERENCES	: None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
  
  SUBROUTINE stomate_clear

  !! 1. Deallocate all dynamics variables

    IF (ALLOCATED(veget_cov_max)) DEALLOCATE(veget_cov_max)
    IF (ALLOCATED(ind)) DEALLOCATE(ind)
    IF (ALLOCATED(adapted)) DEALLOCATE(adapted)
    IF (ALLOCATED(regenerate)) DEALLOCATE(regenerate)
    IF (ALLOCATED(humrel_daily)) DEALLOCATE(humrel_daily)
    IF (ALLOCATED(max_eau_var_daily)) DEALLOCATE(max_eau_var_daily)
    IF (ALLOCATED(gdd_init_date)) DEALLOCATE(gdd_init_date)
    IF (ALLOCATED(litterhum_daily)) DEALLOCATE(litterhum_daily)
    IF (ALLOCATED(t2m_daily))  DEALLOCATE(t2m_daily)
    IF (ALLOCATED(t2m_min_daily))  DEALLOCATE(t2m_min_daily)
    IF (ALLOCATED(tsurf_daily))  DEALLOCATE(tsurf_daily)
    IF (ALLOCATED(tsoil_daily)) DEALLOCATE(tsoil_daily)
    IF (ALLOCATED(soilhum_daily)) DEALLOCATE(soilhum_daily)
    IF (ALLOCATED(precip_daily)) DEALLOCATE(precip_daily)
    IF (ALLOCATED(gpp_daily)) DEALLOCATE(gpp_daily)
    IF (ALLOCATED(npp_daily)) DEALLOCATE(npp_daily)
    IF (ALLOCATED(turnover_daily)) DEALLOCATE(turnover_daily)
    IF (ALLOCATED(turnover_littercalc)) DEALLOCATE(turnover_littercalc)
    IF (ALLOCATED(humrel_month)) DEALLOCATE(humrel_month)
    IF (ALLOCATED(humrel_week)) DEALLOCATE(humrel_week)
    IF (ALLOCATED(moiavail_growingseason)) DEALLOCATE(moiavail_growingseason)
    IF (ALLOCATED(t2m_longterm)) DEALLOCATE(t2m_longterm)
    IF (ALLOCATED(t2m_month)) DEALLOCATE(t2m_month)
    IF (ALLOCATED(Tseason)) DEALLOCATE(Tseason)
    IF (ALLOCATED(Tseason_length)) DEALLOCATE(Tseason_length)
    IF (ALLOCATED(Tseason_tmp)) DEALLOCATE(Tseason_tmp)
    IF (ALLOCATED(Tmin_spring_time)) DEALLOCATE(Tmin_spring_time)
    IF (ALLOCATED(onset_date)) DEALLOCATE(onset_date)
    IF (ALLOCATED(begin_leaves)) DEALLOCATE(begin_leaves)
    IF (ALLOCATED(t2m_week)) DEALLOCATE(t2m_week)
    IF (ALLOCATED(tsoil_month)) DEALLOCATE(tsoil_month)
    IF (ALLOCATED(soilhum_month)) DEALLOCATE(soilhum_month)
    IF (ALLOCATED(fireindex)) DEALLOCATE(fireindex)
    IF (ALLOCATED(firelitter)) DEALLOCATE(firelitter)
    IF (ALLOCATED(maxhumrel_lastyear)) DEALLOCATE(maxhumrel_lastyear)
    IF (ALLOCATED(maxhumrel_thisyear)) DEALLOCATE(maxhumrel_thisyear)
    IF (ALLOCATED(minhumrel_lastyear)) DEALLOCATE(minhumrel_lastyear)
    IF (ALLOCATED(minhumrel_thisyear)) DEALLOCATE(minhumrel_thisyear)
    IF (ALLOCATED(maxgppweek_lastyear)) DEALLOCATE(maxgppweek_lastyear)
    IF (ALLOCATED(maxgppweek_thisyear)) DEALLOCATE(maxgppweek_thisyear)
    IF (ALLOCATED(gdd0_lastyear)) DEALLOCATE(gdd0_lastyear)
    IF (ALLOCATED(gdd0_thisyear)) DEALLOCATE(gdd0_thisyear)
    IF (ALLOCATED(precip_lastyear)) DEALLOCATE(precip_lastyear)
    IF (ALLOCATED(precip_thisyear)) DEALLOCATE(precip_thisyear)
    IF (ALLOCATED(gdd_m5_dormance)) DEALLOCATE(gdd_m5_dormance)
    IF (ALLOCATED(gdd_from_growthinit)) DEALLOCATE(gdd_from_growthinit)
    IF (ALLOCATED(gdd_midwinter)) DEALLOCATE(gdd_midwinter)
    IF (ALLOCATED(ncd_dormance)) DEALLOCATE(ncd_dormance)
    IF (ALLOCATED(ngd_minus5))  DEALLOCATE(ngd_minus5)
    IF (ALLOCATED(PFTpresent)) DEALLOCATE(PFTpresent)
    IF (ALLOCATED(npp_longterm)) DEALLOCATE(npp_longterm)
    IF (ALLOCATED(lm_lastyearmax)) DEALLOCATE(lm_lastyearmax)
    IF (ALLOCATED(lm_thisyearmax)) DEALLOCATE(lm_thisyearmax)
    IF (ALLOCATED(maxfpc_lastyear)) DEALLOCATE(maxfpc_lastyear)
    IF (ALLOCATED(maxfpc_thisyear)) DEALLOCATE(maxfpc_thisyear)
    IF (ALLOCATED(turnover_longterm)) DEALLOCATE(turnover_longterm)
    IF (ALLOCATED(gpp_week)) DEALLOCATE(gpp_week)
    IF (ALLOCATED(npp_week)) DEALLOCATE(npp_week)
    IF (ALLOCATED(biomass)) DEALLOCATE(biomass)
    IF (ALLOCATED(EW_grain)) DEALLOCATE(EW_grain)
    IF (ALLOCATED(senescence)) DEALLOCATE(senescence)
    IF (ALLOCATED(when_growthinit)) DEALLOCATE(when_growthinit)
    IF (ALLOCATED(age))  DEALLOCATE(age)
    IF (ALLOCATED(resp_hetero_d)) DEALLOCATE(resp_hetero_d)
    IF (ALLOCATED(resp_hetero_radia)) DEALLOCATE(resp_hetero_radia)
    IF (ALLOCATED(resp_maint_d)) DEALLOCATE(resp_maint_d)
    IF (ALLOCATED(resp_growth_d)) DEALLOCATE(resp_growth_d)
    IF (ALLOCATED(co2_fire)) DEALLOCATE(co2_fire)
    IF (ALLOCATED(co2_to_bm_dgvm)) DEALLOCATE(co2_to_bm_dgvm)
    IF (ALLOCATED(n_to_bm)) DEALLOCATE(n_to_bm)
    IF (ALLOCATED(p_to_bm)) DEALLOCATE(p_to_bm)
    IF (ALLOCATED(veget_lastlight)) DEALLOCATE(veget_lastlight)
    IF (ALLOCATED(everywhere)) DEALLOCATE(everywhere)
    IF (ALLOCATED(need_adjacent)) DEALLOCATE(need_adjacent)
    IF (ALLOCATED(leaf_age)) DEALLOCATE(leaf_age)
    IF (ALLOCATED(leaf_frac)) DEALLOCATE(leaf_frac)
    IF (ALLOCATED(RIP_time)) DEALLOCATE(RIP_time)
    IF (ALLOCATED(time_hum_min)) DEALLOCATE(time_hum_min)
    IF (ALLOCATED(hum_min_dormance)) DEALLOCATE(hum_min_dormance)
    IF (ALLOCATED(litter)) DEALLOCATE(litter)
    IF (ALLOCATED(dead_leaves)) DEALLOCATE(dead_leaves)
    IF (ALLOCATED(som)) DEALLOCATE(som)
    IF (ALLOCATED(lignin_struc)) DEALLOCATE(lignin_struc)
    IF (ALLOCATED(lignin_wood)) DEALLOCATE(lignin_wood)
    IF (ALLOCATED(turnover_time)) DEALLOCATE(turnover_time)
    IF (ALLOCATED(nep_daily)) DEALLOCATE(nep_daily)
    IF (ALLOCATED(nep_monthly)) DEALLOCATE(nep_monthly)
    IF (ALLOCATED(harvest_above_monthly)) DEALLOCATE (harvest_above_monthly)
    IF (ALLOCATED(cflux_prod_monthly)) DEALLOCATE (cflux_prod_monthly)
    IF (ALLOCATED(bm_to_litter)) DEALLOCATE(bm_to_litter)
    IF (ALLOCATED(bm_to_littercalc)) DEALLOCATE(bm_to_littercalc)
    IF (ALLOCATED(herbivores)) DEALLOCATE(herbivores)
    IF (ALLOCATED(resp_maint_part_radia)) DEALLOCATE(resp_maint_part_radia)
    IF (ALLOCATED(resp_maint_radia)) DEALLOCATE(resp_maint_radia)
    IF (ALLOCATED(resp_maint_part)) DEALLOCATE(resp_maint_part)
    IF (ALLOCATED(hori_index)) DEALLOCATE(hori_index)
    IF (ALLOCATED(horipft_index)) DEALLOCATE(horipft_index)
    IF (ALLOCATED(clay_fm)) DEALLOCATE(clay_fm)
    IF (ALLOCATED(silt_fm)) DEALLOCATE(silt_fm) 
    IF (ALLOCATED(bulk_fm)) DEALLOCATE(bulk_fm) 
    IF (ALLOCATED(humrel_daily_fm)) DEALLOCATE(humrel_daily_fm)
    IF (ALLOCATED(litterhum_daily_fm))  DEALLOCATE(litterhum_daily_fm)
    IF (ALLOCATED(t2m_daily_fm))  DEALLOCATE(t2m_daily_fm)
    IF (ALLOCATED(t2m_min_daily_fm))  DEALLOCATE(t2m_min_daily_fm)
    IF (ALLOCATED(tsurf_daily_fm)) DEALLOCATE(tsurf_daily_fm)
    IF (ALLOCATED(tsoil_daily_fm)) DEALLOCATE(tsoil_daily_fm)
    IF (ALLOCATED(soilhum_daily_fm))  DEALLOCATE(soilhum_daily_fm)
    IF (ALLOCATED(precip_fm)) DEALLOCATE(precip_fm)
    IF (ALLOCATED(gpp_daily_fm))  DEALLOCATE(gpp_daily_fm)
    IF (ALLOCATED(veget_fm)) DEALLOCATE(veget_fm)
    IF (ALLOCATED(veget_max_fm)) DEALLOCATE(veget_max_fm)
    IF (ALLOCATED(lai_fm))  DEALLOCATE(lai_fm)
    IF (ALLOCATED(drainage_fm)) DEALLOCATE(drainage_fm) 
    IF (ALLOCATED(runoff_fm)) DEALLOCATE(runoff_fm) 
    IF (ALLOCATED(max_eau_var_daily_fm)) DEALLOCATE(max_eau_var_daily_fm)
    !
    IF (ALLOCATED(ok_equilibrium)) DEALLOCATE(ok_equilibrium)
    IF (ALLOCATED(carbon_eq)) DEALLOCATE(carbon_eq)
    IF (ALLOCATED(matrixA)) DEALLOCATE(matrixA)
    IF (ALLOCATED(vectorB)) DEALLOCATE(vectorB)
    IF (ALLOCATED(MatrixV)) DEALLOCATE(MatrixV)
    IF (ALLOCATED(VectorU)) DEALLOCATE(VectorU)
    IF (ALLOCATED(MatrixW)) DEALLOCATE(MatrixW)
    IF (ALLOCATED(previous_stock)) DEALLOCATE(previous_stock)
    IF (ALLOCATED(current_stock)) DEALLOCATE(current_stock) 
    IF (ALLOCATED(CN_som_litter_longterm)) DEALLOCATE(CN_som_litter_longterm) 
    IF (ALLOCATED(KF)) DEALLOCATE (KF)
    IF (ALLOCATED(k_latosa_adapt)) DEALLOCATE (k_latosa_adapt)
    IF (ALLOCATED(rue_longterm)) DEALLOCATE (rue_longterm)
    IF (ALLOCATED(nbp_accu)) DEALLOCATE(nbp_accu)
    IF (ALLOCATED(nbp_flux)) DEALLOCATE(nbp_flux)

    IF (ALLOCATED(clay_fm_g)) DEALLOCATE(clay_fm_g)
    IF (ALLOCATED(silt_fm_g)) DEALLOCATE(silt_fm_g) 
    IF (ALLOCATED(bulk_fm_g)) DEALLOCATE(bulk_fm_g) 
    IF (ALLOCATED(humrel_daily_fm_g)) DEALLOCATE(humrel_daily_fm_g)
    IF (ALLOCATED(litterhum_daily_fm_g))  DEALLOCATE(litterhum_daily_fm_g)
    IF (ALLOCATED(t2m_daily_fm_g))  DEALLOCATE(t2m_daily_fm_g)
    IF (ALLOCATED(t2m_min_daily_fm_g))  DEALLOCATE(t2m_min_daily_fm_g)
    IF (ALLOCATED(tsurf_daily_fm_g)) DEALLOCATE(tsurf_daily_fm_g)
    IF (ALLOCATED(tsoil_daily_fm_g)) DEALLOCATE(tsoil_daily_fm_g)
    IF (ALLOCATED(soilhum_daily_fm_g))  DEALLOCATE(soilhum_daily_fm_g)
    IF (ALLOCATED(precip_fm_g)) DEALLOCATE(precip_fm_g)
    IF (ALLOCATED(gpp_daily_fm_g))  DEALLOCATE(gpp_daily_fm_g)
    IF (ALLOCATED(veget_fm_g)) DEALLOCATE(veget_fm_g)
    IF (ALLOCATED(veget_max_fm_g)) DEALLOCATE(veget_max_fm_g)
    IF (ALLOCATED(lai_fm_g))  DEALLOCATE(lai_fm_g)
    IF (ALLOCATED(drainage_fm_g))  DEALLOCATE(drainage_fm_g)    
    IF (ALLOCATED(runoff_fm_g))  DEALLOCATE(drainage_fm_g)    
    IF (ALLOCATED(max_eau_var_daily_fm_g)) DEALLOCATE(max_eau_var_daily_fm_g)
 
    IF (ALLOCATED(isf)) DEALLOCATE(isf)
    IF (ALLOCATED(nf_written)) DEALLOCATE(nf_written)
    IF (ALLOCATED(nf_cumul)) DEALLOCATE(nf_cumul)
    IF (ALLOCATED(nforce)) DEALLOCATE(nforce)
    IF (ALLOCATED(control_moist)) DEALLOCATE(control_moist)
    IF (ALLOCATED(control_temp)) DEALLOCATE(control_temp)
    IF (ALLOCATED(carbon_input)) DEALLOCATE(carbon_input)
    IF (ALLOCATED(nitrogen_input)) DEALLOCATE(nitrogen_input)
    IF ( ALLOCATED (horip10_index)) DEALLOCATE (horip10_index)
    IF ( ALLOCATED (horip100_index)) DEALLOCATE (horip100_index)
    IF ( ALLOCATED (horip11_index)) DEALLOCATE (horip11_index)
    IF ( ALLOCATED (horip101_index)) DEALLOCATE (horip101_index)
    IF ( ALLOCATED (prod10)) DEALLOCATE (prod10)
    IF ( ALLOCATED (prod100)) DEALLOCATE (prod100)
    IF ( ALLOCATED (flux10)) DEALLOCATE (flux10)
    IF ( ALLOCATED (flux100)) DEALLOCATE (flux100)
    IF ( ALLOCATED (convflux)) DEALLOCATE (convflux)
    IF ( ALLOCATED (cflux_prod10)) DEALLOCATE (cflux_prod10)
    IF ( ALLOCATED (cflux_prod100)) DEALLOCATE (cflux_prod100)
    IF ( ALLOCATED (prod10_harvest)) DEALLOCATE (prod10_harvest)
    IF ( ALLOCATED (prod100_harvest)) DEALLOCATE (prod100_harvest)
    IF ( ALLOCATED (flux10_harvest)) DEALLOCATE (flux10_harvest)
    IF ( ALLOCATED (flux100_harvest)) DEALLOCATE (flux100_harvest)
    IF ( ALLOCATED (convflux_harvest)) DEALLOCATE (convflux_harvest)
    IF ( ALLOCATED (cflux_prod10_harvest)) DEALLOCATE (cflux_prod10_harvest)
    IF ( ALLOCATED (cflux_prod100_harvest)) DEALLOCATE (cflux_prod100_harvest)
    IF ( ALLOCATED (woodharvestpft)) DEALLOCATE (woodharvestpft)
    IF ( ALLOCATED (harvest_above)) DEALLOCATE (harvest_above)
    IF ( ALLOCATED (harvest_bioenergy)) DEALLOCATE (harvest_bioenergy)
    IF ( ALLOCATED (som_input_daily)) DEALLOCATE (som_input_daily)
    IF ( ALLOCATED (control_temp_daily)) DEALLOCATE (control_temp_daily)
    IF ( ALLOCATED (control_moist_daily)) DEALLOCATE (control_moist_daily)


    IF ( ALLOCATED (drainage_daily)) DEALLOCATE(drainage_daily) 
    IF ( ALLOCATED (drainage_longterm)) DEALLOCATE(drainage_longterm)
    IF ( ALLOCATED (runoff_daily)) DEALLOCATE(runoff_daily) 
    IF ( ALLOCATED (runoff_longterm)) DEALLOCATE(runoff_longterm)
    IF ( ALLOCATED (n_uptake_daily)) DEALLOCATE(n_uptake_daily) 
    IF ( ALLOCATED (n_mineralisation_d)) DEALLOCATE(n_mineralisation_d) 
    IF ( ALLOCATED (cn_leaf_avg_season)) DEALLOCATE (cn_leaf_avg_season) 
    IF ( ALLOCATED (nstress_season)) DEALLOCATE (nstress_season) 
    IF ( ALLOCATED (soil_n_min)) DEALLOCATE (soil_n_min) 
    IF ( ALLOCATED (p_O2)) DEALLOCATE (p_O2) 
    IF ( ALLOCATED (bact)) DEALLOCATE (bact) 
    IF ( ALLOCATED (fpc_max)) DEALLOCATE (fpc_max)

    IF (ALLOCATED(record_reserve_grass)) DEALLOCATE (record_reserve_grass)
    IF (ALLOCATED(litter_avail)) DEALLOCATE(litter_avail)
    IF (ALLOCATED(litter_not_avail)) DEALLOCATE(litter_not_avail)
    IF (ALLOCATED(litter_avail_frac)) DEALLOCATE(litter_avail_frac)
    IF (ALLOCATED(sla_calc)) DEALLOCATE(sla_calc)
    IF (ALLOCATED(wshtotsum)) DEALLOCATE(wshtotsum)
    IF (ALLOCATED(sr_ugb)) DEALLOCATE(sr_ugb)
    IF (ALLOCATED(compt_ugb)) DEALLOCATE(compt_ugb)
    IF (ALLOCATED(nb_ani)) DEALLOCATE(nb_ani)
    IF (ALLOCATED(grazed_frac)) DEALLOCATE(grazed_frac)
    IF (ALLOCATED(import_yield)) DEALLOCATE(import_yield)
    IF (ALLOCATED(t2m_14)) DEALLOCATE(t2m_14)
    IF (ALLOCATED(sla_age1)) DEALLOCATE(sla_age1)
    IF (ALLOCATED(when_growthinit_cut)) DEALLOCATE(when_growthinit_cut)
    IF (ALLOCATED(nb_grazingdays)) DEALLOCATE(nb_grazingdays)
    IF (ALLOCATED(snowfall_daily)) DEALLOCATE(snowfall_daily)
    IF (ALLOCATED(snowmass_daily)) DEALLOCATE(snowmass_daily)
    IF (ALLOCATED(tmc_topgrass_daily)) DEALLOCATE(tmc_topgrass_daily)
    IF (ALLOCATED(after_snow)) DEALLOCATE(after_snow)
    IF (ALLOCATED(after_wet)) DEALLOCATE(after_wet)
    IF (ALLOCATED(wet1day)) DEALLOCATE(wet1day)
    IF (ALLOCATED(wet2day)) DEALLOCATE(wet2day)
    IF (ALLOCATED(grm_nfert)) DEALLOCATE(grm_nfert)
    IF (ALLOCATED(grm_pfert)) DEALLOCATE(grm_pfert)
    IF (ALLOCATED(BNF_clover)) DEALLOCATE(BNF_clover)
    IF (ALLOCATED(BNF_clover_daily)) DEALLOCATE(BNF_clover_daily)
    IF (ALLOCATED(GRM_devstage)) DEALLOCATE(GRM_devstage)
    IF (ALLOCATED(p_input_dt)) DEALLOCATE(p_input_dt)
    IF (ALLOCATED(agri_nfert)) DEALLOCATE(agri_nfert)
    IF (ALLOCATED(agri_pfert)) DEALLOCATE(agri_pfert)
    IF (ALLOCATED(rest_Ngrazing)) DEALLOCATE(rest_Ngrazing)
    IF (ALLOCATED(grazed_above)) DEALLOCATE (grazed_above)
    IF (ALLOCATED(days_senescence)) DEALLOCATE(days_senescence)
    IF (ALLOCATED(fclover)) DEALLOCATE(fclover)

    IF ( ALLOCATED (p_uptake_daily)) DEALLOCATE(p_uptake_daily)
    IF ( ALLOCATED (p_leaching_daily)) DEALLOCATE(p_leaching_daily)
    IF ( ALLOCATED (p_ssorption_daily)) DEALLOCATE(p_ssorption_daily)
    IF ( ALLOCATED (p_mineralisation_d)) DEALLOCATE(p_mineralisation_d)
    IF ( ALLOCATED (soil_p_min)) DEALLOCATE (soil_p_min)
    IF ( ALLOCATED (f_Pdissolved)) DEALLOCATE (f_Pdissolved)
    IF ( ALLOCATED (pstress_season)) DEALLOCATE (pstress_season) 

    IF ( ALLOCATED (lai_target_longterm)) DEALLOCATE (lai_target_longterm) 
    IF ( ALLOCATED (lai_target)) DEALLOCATE (lai_target)
 !! 2. reset l_first

    l_first_stomate=.TRUE.

 !! 3. call to clear functions

    CALL season_clear
    CALL stomatelpj_clear
    CALL littercalc_clear
    CALL vmax_clear
 
  END SUBROUTINE stomate_clear


!! ================================================================================================================================
!! SUBROUTINE 	: stomate_var_init
!!
!>\BRIEF        Initialize variables of stomate with a none-zero initial value.
!! Subroutine is called only if ::ok_stomate = .TRUE. STOMATE diagnoses some 
!! variables for SECHIBA : assim_param, deadleaf_cover, etc. These variables can 
!! be recalculated from STOMATE's prognostic variables. Note that height is
!! saved in SECHIBA.
!!
!! DESCRIPTION  : None
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S): leaf age (::leaf_age) and fraction of leaves in leaf 
!! age class (::leaf_frac). The maximum water on vegetation available for 
!! interception, fraction of soil covered by dead leaves
!! (::deadleaf_cover) and assimilation parameters (:: assim_param).
!!
!! REFERENCE(S)	: None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
  
  SUBROUTINE stomate_var_init &
       &  (kjpindex, veget_cov_max, leaf_age, leaf_frac, &
       &   dead_leaves, &
       &   veget, lai, deadleaf_cover, assim_param)


  !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER,INTENT(in)                             :: kjpindex        !! Domain size - terrestrial pixels only
    REAL,DIMENSION(kjpindex,nvm),INTENT(in)        :: veget           !! Fraction of pixel covered by PFT. Fraction 
                                                                             !! accounts for none-biological land covers 
                                                                             !! (unitless) 
    REAL,DIMENSION(kjpindex,nvm),INTENT(in)        :: veget_cov_max   !! Fractional coverage: maximum share of the pixel 
                                                                             !! covered by a PFT (unitless) 
    REAL,DIMENSION(kjpindex,nvm,nlitt),INTENT(in)  :: dead_leaves     !! Metabolic and structural fraction of dead leaves 
                                                                             !! per ground area 
                                                                             !! @tex $(gC m^{-2})$ @endtex 
    REAL,DIMENSION(kjpindex,nvm),INTENT(in)        :: lai             !! Leaf area index 
                                                                             !! @tex $(m^2 m{-2})$ @endtex 
    
    !! 0.2 Modified variables
    REAL,DIMENSION(kjpindex,nvm,npco2),INTENT(inout) :: assim_param   !! leaf nutrient derive vcmax,jmax and leafN,nue
                                                                             !! photosynthesis  
    REAL,DIMENSION(kjpindex,nvm,nparts,nleafages),INTENT(inout) :: leaf_age  !! Age of different leaf classes per PFT (days)
    REAL,DIMENSION(kjpindex,nvm,nparts,nleafages),INTENT(inout) :: leaf_frac !! Fraction of leaves in leaf age class per PFT 
                                                                             !! (unitless; 1) 
    
    !! 0.3 Output variables

    REAL,DIMENSION(kjpindex), INTENT (out)         :: deadleaf_cover  !! Fraction of soil covered by dead leaves 
                                                                             !! (unitless) 


    ! 0.4 Local variables
   
    REAL,PARAMETER                                 :: dt_0 = zero     !! Dummy time step, must be zero
    REAL,DIMENSION(kjpindex,nvm,2)                   :: vcmax         !! Dummy vcmax 
                                                                             !! @tex $(\mu mol m^{-2} s^{-1})$ @endtex
                                                                             !! third dimension: actual Vcmax=1 ; potential Vcmax=2
    REAL,DIMENSION(kjpindex,nvm,2)                   :: jmax          !! Dummy jmax 
                                                                             !! @tex $(\mu mol m^{-2} s^{-1})$ @endtex
                                                                             !! third dimension: actual Jmax=1 ; potential Jmax=2

    REAL,DIMENSION(kjpindex,nvm)                   :: nue             !! Nitrogen use Efficiency with impact of leaf age (umol CO2 (gN)-1 s-1) 
                                                                             !! @tex $(\mu mol m^{-2} s^{-1})$ @endtex 
    REAL,DIMENSION(kjpindex,nvm,nparts,nleafages)         :: leaf_age_tmp    !! Temporary variable
    REAL,DIMENSION(kjpindex,nvm,nparts,nleafages)         :: leaf_frac_tmp   !! Temporary variable
                                                                             !! (unitless; 1)     
    INTEGER                                        :: j               !! Index (untiless)
    
!_ ================================================================================================================================   


    !! 1.1 vcmax (stomate_vmax.f90)
    CALL vmax (kjpindex, dt_0, biomass, leaf_age, leaf_frac, vcmax, jmax, nue )

    ! Calculate assim_param if it was not found in the restart file
    IF (ALL(assim_param(:,:,:)==val_exp)) THEN
       ! Use temporary leaf_age_tmp and leaf_frac_tmp to preserve the input variables from being modified by the subroutine vmax.
       leaf_age_tmp(:,:,:,:)=leaf_age(:,:,:,:)
       leaf_frac_tmp(:,:,:,:)=leaf_frac(:,:,:,:)

       !! 1.2 Calculate a temporary vcmax (stomate_vmax.f90)
       CALL vmax (kjpindex, dt_0, biomass, leaf_age_tmp, leaf_frac_tmp, vcmax, jmax, nue )

    END IF
    
    !! Fill array which stores all parameter related to Assimilation to be
    !passed to sechiba
    assim_param(:,:,ivcmax) = vcmax(:,:,1)
    assim_param(:,:,ijmax)  = jmax(:,:,1)
    assim_param(:,:,inue)   = nue(:,:) 
    assim_param(:,:,ileafN) = biomass(:,:,ileaf,initrogen)
!DSGpot    
    assim_param(:,:,ivcmax_pot) = vcmax(:,:,2)
    assim_param(:,:,ijmax_pot)  = jmax(:,:,2)
!DSGpot    

    !! 2. Dead leaf cover (stomate_litter.f90)
    CALL deadleaf (kjpindex, veget_cov_max, dead_leaves, deadleaf_cover, &
                   sla_calc)
    
  END SUBROUTINE stomate_var_init


!! ================================================================================================================================
!! INTERFACE 	: stomate_accu
!!
!>\BRIEF        Accumulate a variable for the time period specified by 
!! dt_sechiba or calculate the mean value over the period of dt_stomate
!! 
!! DESCRIPTION : Accumulate a variable for the time period specified by 
!! dt_sechiba or calculate the mean value over the period of dt_stomate.
!! stomate_accu interface can be used for variables having 1, 2 or 3 dimensions.
!! The corresponding subruoutine stomate_accu_r1d, stomate_accu_r2d or
!! stomate_accu_r3d will be selected through the interface depending on the number of dimensions.
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S): accumulated or mean variable ::field_out:: 
!!
!! REFERENCE(S)	: None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
  
  SUBROUTINE stomate_accu_r1d (ldmean, field_in, field_out)
    
  !! 0. Variable and parameter declaration

    !! 0.1 Input variables
    LOGICAL,INTENT(in)                     :: ldmean    !! Flag to calculate the mean over
    REAL,DIMENSION(:),INTENT(in)    :: field_in  !! Field that needs to be accumulated
    
    !! 0.2 Modified variables
    REAL,DIMENSION(:),INTENT(inout) :: field_out !! Accumulated or mean field

!_ ================================================================================================================================

  !! 1. Accumulate field

    field_out(:) = field_out(:)+field_in(:)*dt_sechiba
   
  !! 2. Mean fields

    IF (ldmean) THEN
       field_out(:) = field_out(:)/dt_stomate
    ENDIF

  END SUBROUTINE stomate_accu_r1d

  SUBROUTINE stomate_accu_r2d (ldmean, field_in, field_out)
    
  !! 0. Variable and parameter declaration

    !! 0.1 Input variables
    LOGICAL,INTENT(in)                       :: ldmean    !! Flag to calculate the mean over
    REAL,DIMENSION(:,:),INTENT(in)    :: field_in  !! Field that needs to be accumulated
    
    !! 0.2 Modified variables
    REAL,DIMENSION(:,:),INTENT(inout) :: field_out !! Accumulated or mean field

!_ ================================================================================================================================

  !! 1. Accumulate field

    field_out(:,:) = field_out(:,:)+field_in(:,:)*dt_sechiba
   
  !! 2. Mean fields

    IF (ldmean) THEN
       field_out(:,:) = field_out(:,:)/dt_stomate
    ENDIF

  END SUBROUTINE stomate_accu_r2d

  SUBROUTINE stomate_accu_r3d (ldmean, field_in, field_out)
    
  !! 0. Variable and parameter declaration

    !! 0.1 Input variables
    LOGICAL,INTENT(in)                         :: ldmean    !! Flag to calculate the mean over
    REAL,DIMENSION(:,:,:),INTENT(in)    :: field_in  !! Field that needs to be accumulated
    
    !! 0.2 Modified variables
    REAL,DIMENSION(:,:,:),INTENT(inout) :: field_out !! Accumulated or mean field

!_ ================================================================================================================================

  !! 1. Accumulate field

    field_out(:,:,:) = field_out(:,:,:)+field_in(:,:,:)*dt_sechiba
   
  !! 2. Mean fields

    IF (ldmean) THEN
       field_out(:,:,:) = field_out(:,:,:)/dt_stomate
    ENDIF

  END SUBROUTINE stomate_accu_r3d

  SUBROUTINE stomate_accu_r4d (ldmean, field_in, field_out)
    
  !! 0. Variable and parameter declaration

    !! 0.1 Input variables
    LOGICAL,INTENT(in)                         :: ldmean    !! Flag to calculate the mean over
    REAL,DIMENSION(:,:,:,:),INTENT(in)    :: field_in  !! Field that needs to be accumulated
    
    !! 0.2 Modified variables
    REAL,DIMENSION(:,:,:,:),INTENT(inout) :: field_out !! Accumulated or mean field

!_ ================================================================================================================================

  !! 1. Accumulate field

    field_out(:,:,:,:) = field_out(:,:,:,:)+field_in(:,:,:,:)*dt_sechiba
   
  !! 2. Mean fields

    IF (ldmean) THEN
       field_out(:,:,:,:) = field_out(:,:,:,:)/dt_stomate
    ENDIF

  END SUBROUTINE stomate_accu_r4d

!! ================================================================================================================================
!! SUBROUTINE 	: init_forcing
!!
!>\BRIEF        Allocate memory for the variables containing the forcing data.
!! The maximum size of the allocated memory is specified in run definition file
!! (::max_totsize) and needs to be a compromise between charging the memory and 
!! accessing disks to get the forcing data.
!!
!! DESCRIPTION : None
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S): Strictly speaking the subroutine has no output 
!! variables. However, the routine allocates memory for later use. 
!!
!! REFERENCE(S)	: None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
  
  SUBROUTINE init_forcing (kjpindex,nsfm,nsft_loc)
    
  !! 0. Variable and parameter declaration

    !! 0.1 Input variables
    INTEGER,INTENT(in) :: kjpindex !! Domain size - terrestrial pixels only (unitless)
    INTEGER,INTENT(in) :: nsfm     !! Number of time steps that can be stored in memory (unitless)
    INTEGER,INTENT(in) :: nsft_loc !! Number of time steps in a year (unitless)

   !! 0.2 Output variables

   !! 0.3 Modified variables

   !! 0.4 Local variables

    LOGICAL                   :: l_error  !! Check errors in netcdf call
    INTEGER            :: ier      !! Check errors in netcdf call
!_ ================================================================================================================================
    
  !! 1. Allocate memory

    ! Note ::nvm is number of PFTs and ::nslm is number of soil layers
    l_error = .FALSE.
    ALLOCATE(clay_fm(kjpindex,nsfm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Problem with memory allocation: forcing variables clay_fm ',kjpindex,nsfm
       STOP 'init_forcing'
    ENDIF

    ALLOCATE(silt_fm(kjpindex,nsfm),stat=ier) 
    l_error = l_error .OR. (ier /= 0) 
    IF (l_error) THEN 
       WRITE(numout,*) 'Problem with memory allocation: forcing variables silt_fm ',kjpindex,nsfm 
       STOP 'init_forcing' 
    ENDIF 
    ALLOCATE(bulk_fm(kjpindex,nsfm),stat=ier) 
    l_error = l_error .OR. (ier /= 0) 
    IF (l_error) THEN 
       WRITE(numout,*) 'Problem with memory allocation: forcing variables bulk_fm ',kjpindex,nsfm 
       STOP 'init_forcing' 
    ENDIF 

    ALLOCATE(humrel_daily_fm(kjpindex,nvm,nsfm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Problem with memory allocation: forcing variables humrel_daily_fm ',kjpindex,nvm,nsfm
       STOP 'init_forcing'
    ENDIF
    ALLOCATE(litterhum_daily_fm(kjpindex,nsfm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Problem with memory allocation: forcing variables litterhum_daily_fm ',kjpindex,nsfm
       STOP 'init_forcing'
    ENDIF
    ALLOCATE(t2m_daily_fm(kjpindex,nsfm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Problem with memory allocation: forcing variables t2m_daily_fm ',kjpindex,nsfm
       STOP 'init_forcing'
    ENDIF
    ALLOCATE(t2m_min_daily_fm(kjpindex,nsfm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Problem with memory allocation: forcing variables t2m_min_daily_fm ',kjpindex,nsfm
       STOP 'init_forcing'
    ENDIF
    ALLOCATE(tsurf_daily_fm(kjpindex,nsfm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Problem with memory allocation: forcing variables tsurf_daily_fm ',kjpindex,nsfm
       STOP 'init_forcing'
    ENDIF
    ALLOCATE(tsoil_daily_fm(kjpindex,nslm,nsfm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Problem with memory allocation: forcing variables tsoil_daily_fm ',kjpindex,nslm,nsfm
       STOP 'init_forcing'
    ENDIF
    ALLOCATE(soilhum_daily_fm(kjpindex,nslm,nsfm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Problem with memory allocation: forcing variables soilhum_daily_fm ',kjpindex,nslm,nsfm
       STOP 'init_forcing'
    ENDIF
    ALLOCATE(precip_fm(kjpindex,nsfm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Problem with memory allocation: forcing variables precip_fm ',kjpindex,nsfm
       STOP 'init_forcing'
    ENDIF
    ALLOCATE(gpp_daily_fm(kjpindex,nvm,nsfm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Problem with memory allocation: forcing variables gpp_daily_fm ',kjpindex,nvm,nsfm
       STOP 'init_forcing'
    ENDIF
    ALLOCATE(veget_fm(kjpindex,nvm,nsfm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Problem with memory allocation: forcing variables veget_fm ',kjpindex,nvm,nsfm
       STOP 'init_forcing'
    ENDIF
    ALLOCATE(veget_max_fm(kjpindex,nvm,nsfm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Problem with memory allocation: forcing variables veget_max_fm ',kjpindex,nvm,nsfm
       STOP 'init_forcing'
    ENDIF
    ALLOCATE(lai_fm(kjpindex,nvm,nsfm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Problem with memory allocation: forcing variables lai_fm ',kjpindex,nvm,nsfm
       STOP 'init_forcing'
    ENDIF
    ALLOCATE(drainage_fm(kjpindex,nvm,nsfm),stat=ier) 
    l_error = l_error .OR. (ier /= 0) 
    IF (l_error) THEN 
       WRITE(numout,*) 'Problem with memory allocation: forcing variables drainage_fm ',kjpindex,nvm,nsfm 
       STOP 'init_forcing' 
    ENDIF 
    ALLOCATE(runoff_fm(kjpindex,nvm,nsfm),stat=ier) 
    l_error = l_error .OR. (ier /= 0) 
    IF (l_error) THEN 
       WRITE(numout,*) 'Problem with memory allocation: forcing variables runoff_fm ',kjpindex,nvm,nsfm 
       STOP 'init_forcing' 
    ENDIF 
    ALLOCATE(max_eau_var_daily_fm(kjpindex,nsfm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Problem with memory allocation: forcing variables humrel_daily_fm ',kjpindex,nsfm
       STOP 'init_forcing'
    ENDIF
    ALLOCATE(isf(nsfm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Problem with memory allocation: forcing variables isf ',nsfm
       STOP 'init_forcing'
    ENDIF
    ALLOCATE(nf_written(nsft_loc),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Problem with memory allocation: forcing variables nf_written ',nsft_loc
       STOP 'init_forcing'
    ENDIF
    ALLOCATE(nf_cumul(nsft_loc),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Problem with memory allocation: forcing variables nf_cumul ',nsft_loc
       STOP 'init_forcing'
    ENDIF
    
  !! 2. Allocate memory for the root processor only (parallel computing)

    ! Where, ::nbp_glo is the number of global continental points
    IF (is_root_prc) THEN
       ALLOCATE(clay_fm_g(nbp_glo,nsfm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) 'Problem with memory allocation: forcing variables clay_fm_g ',nbp_glo,nsfm
          STOP 'init_forcing'
       ENDIF
       ALLOCATE(silt_fm_g(nbp_glo,nsfm),stat=ier) 
       l_error = l_error .OR. (ier /= 0) 
       IF (l_error) THEN 
          WRITE(numout,*) 'Problem with memory allocation: forcing variables silt_fm_g ',nbp_glo,nsfm 
          STOP 'init_forcing' 
       ENDIF
       ALLOCATE(bulk_fm_g(nbp_glo,nsfm),stat=ier) 
       l_error = l_error .OR. (ier /= 0) 
       IF (l_error) THEN 
          WRITE(numout,*) 'Problem with memory allocation: forcing variables bulk_fm_g ',nbp_glo,nsfm 
          STOP 'init_forcing' 
       ENDIF
       ALLOCATE(humrel_daily_fm_g(nbp_glo,nvm,nsfm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) 'Problem with memory allocation: forcing variables humrel_daily_fm_g ',nbp_glo,nvm,nsfm
          STOP 'init_forcing'
       ENDIF
       ALLOCATE(litterhum_daily_fm_g(nbp_glo,nsfm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) 'Problem with memory allocation: forcing variables litterhum_daily_fm_g ',nbp_glo,nsfm
          STOP 'init_forcing'
       ENDIF
       ALLOCATE(t2m_daily_fm_g(nbp_glo,nsfm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) 'Problem with memory allocation: forcing variables t2m_daily_fm_g ',nbp_glo,nsfm
          STOP 'init_forcing'
       ENDIF
       ALLOCATE(t2m_min_daily_fm_g(nbp_glo,nsfm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) 'Problem with memory allocation: forcing variables t2m_min_daily_fm_g ',nbp_glo,nsfm
          STOP 'init_forcing'
       ENDIF
       ALLOCATE(tsurf_daily_fm_g(nbp_glo,nsfm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) 'Problem with memory allocation: forcing variables tsurf_daily_fm_g ',nbp_glo,nsfm
          STOP 'init_forcing'
       ENDIF
       ALLOCATE(tsoil_daily_fm_g(nbp_glo,nslm,nsfm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) 'Problem with memory allocation: forcing variables tsoil_daily_fm_g ',nbp_glo,nslm,nsfm
          STOP 'init_forcing'
       ENDIF
       ALLOCATE(soilhum_daily_fm_g(nbp_glo,nslm,nsfm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) 'Problem with memory allocation: forcing variables soilhum_daily_fm_g ',nbp_glo,nslm,nsfm
          STOP 'init_forcing'
       ENDIF
       ALLOCATE(precip_fm_g(nbp_glo,nsfm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) 'Problem with memory allocation: forcing variables precip_fm_g ',nbp_glo,nsfm
          STOP 'init_forcing'
       ENDIF
       ALLOCATE(gpp_daily_fm_g(nbp_glo,nvm,nsfm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) 'Problem with memory allocation: forcing variables gpp_daily_fm_g ',nbp_glo,nvm,nsfm
          STOP 'init_forcing'
       ENDIF
       ALLOCATE(veget_fm_g(nbp_glo,nvm,nsfm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) 'Problem with memory allocation: forcing variables veget_fm_g ',nbp_glo,nvm,nsfm
          STOP 'init_forcing'
       ENDIF
       ALLOCATE(veget_max_fm_g(nbp_glo,nvm,nsfm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) 'Problem with memory allocation: forcing variables veget_max_fm_g ',nbp_glo,nvm,nsfm
          STOP 'init_forcing'
       ENDIF
       ALLOCATE(lai_fm_g(nbp_glo,nvm,nsfm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) 'Problem with memory allocation: forcing variables lai_fm_g ',nbp_glo,nvm,nsfm
          STOP 'init_forcing'
       ENDIF

       ALLOCATE(drainage_fm_g(nbp_glo,nvm,nsfm),stat=ier) 
       l_error = l_error .OR. (ier /= 0) 
       IF (l_error) THEN 
          WRITE(numout,*) 'Problem with memory allocation: forcing variables drainage_fm_g ',nbp_glo,nvm,nsfm 
          STOP 'init_forcing' 
       ENDIF
       ALLOCATE(runoff_fm_g(nbp_glo,nvm,nsfm),stat=ier) 
       l_error = l_error .OR. (ier /= 0) 
       IF (l_error) THEN 
          WRITE(numout,*) 'Problem with memory allocation: forcing variables runoff_fm_g ',nbp_glo,nvm,nsfm 
          STOP 'init_forcing' 
       ENDIF
       ALLOCATE(max_eau_var_daily_fm_g(nbp_glo,nsfm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) 'Problem with memory allocation: forcing variables max_eau_var_daily_fm_g ',nbp_glo,nsfm
          STOP 'init_forcing'
       ENDIF
    ELSE
       ! Allocate memory for co-processors
       ALLOCATE(clay_fm_g(0,nsfm),stat=ier)
       ALLOCATE(silt_fm_g(0,nsfm),stat=ier) 
       ALLOCATE(bulk_fm_g(0,nsfm),stat=ier) 
       ALLOCATE(humrel_daily_fm_g(0,nvm,nsfm),stat=ier)
       ALLOCATE(litterhum_daily_fm_g(0,nsfm),stat=ier)
       ALLOCATE(t2m_daily_fm_g(0,nsfm),stat=ier)
       ALLOCATE(t2m_min_daily_fm_g(0,nsfm),stat=ier)
       ALLOCATE(tsurf_daily_fm_g(0,nsfm),stat=ier)
       ALLOCATE(tsoil_daily_fm_g(0,nslm,nsfm),stat=ier)
       ALLOCATE(soilhum_daily_fm_g(0,nslm,nsfm),stat=ier)
       ALLOCATE(precip_fm_g(0,nsfm),stat=ier)
       ALLOCATE(gpp_daily_fm_g(0,nvm,nsfm),stat=ier)
       ALLOCATE(veget_fm_g(0,nvm,nsfm),stat=ier)
       ALLOCATE(veget_max_fm_g(0,nvm,nsfm),stat=ier)
       ALLOCATE(lai_fm_g(0,nvm,nsfm),stat=ier)
       ALLOCATE(drainage_fm_g(0,nvm,nsfm),stat=ier) 
       ALLOCATE(runoff_fm_g  (0,nvm,nsfm),stat=ier) 
       ALLOCATE(max_eau_var_daily_fm_g(0,nsfm),stat=ier)
    ENDIF ! is_root_proc
    
    IF (l_error) THEN
       WRITE(numout,*) 'Problem with memory allocation: forcing variables'
       STOP 'init_forcing'
    ENDIF

  !! 3. Initilaize variables

    CALL forcing_zero
    
  END SUBROUTINE init_forcing


!! ================================================================================================================================
!! SUBROUTINE 	: forcing_zero
!!
!>\BRIEF        Initialize variables containing the forcing data; variables are 
!! set to zero.
!!
!! DESCRIPTION  : None
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S): None
!!
!! REFERENCES	: None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
  
  SUBROUTINE forcing_zero
    
    clay_fm(:,:) = zero
    silt_fm(:,:) = zero 
    bulk_fm(:,:) = zero 
    humrel_daily_fm(:,:,:) = zero
    litterhum_daily_fm(:,:) = zero
    t2m_daily_fm(:,:) = zero
    t2m_min_daily_fm(:,:) = zero
    tsurf_daily_fm(:,:) = zero
    tsoil_daily_fm(:,:,:) = zero
    soilhum_daily_fm(:,:,:) = zero
    precip_fm(:,:) = zero
    gpp_daily_fm(:,:,:) = zero
    veget_fm(:,:,:) = zero
    veget_max_fm(:,:,:) = zero
    lai_fm(:,:,:) = zero
    drainage_fm(:,:,:) = zero     
    runoff_fm  (:,:,:) = zero     
    max_eau_var_daily_fm(:,:) = zero
  END SUBROUTINE forcing_zero


!! ================================================================================================================================
!! SUBROUTINE 	: forcing_write
!!
!>\BRIEF        Appends data values to a netCDF file containing the forcing 
!! variables of the general processes in stomate.
!!
!! DESCRIPTION  : None
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S): netCDF file
!!
!! REFERENCES	: None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
  
  SUBROUTINE forcing_write(forcing_id,ibeg,iend)
    
  !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER,INTENT(in)      :: forcing_id  !! File identifer of forcing file, assigned when netcdf is created
    INTEGER,INTENT(in)      :: ibeg, iend  !! First and last time step to be written

    !! 0.2 Output variables

    !! 0.3 Modified variables

    !! 0.4 Local variables

    INTEGER                 :: ii          !! Index of isf where isf is the number of time steps that can be 
                                                  !! stored in memory 
    INTEGER                 :: iblocks     !! Index of block that is written
    INTEGER                 :: nblocks     !! Number of blocks that needs to be written
    INTEGER                 :: ier         !! Check errors in netcdf call
    INTEGER,DIMENSION(0:2)  :: ifirst      !! First block in memory - changes with iblocks
    INTEGER,DIMENSION(0:2)  :: ilast       !! Last block in memory - changes with iblocks
    INTEGER,PARAMETER       :: ndm = 10    !! Maximum number of dimensions
    INTEGER,DIMENSION(ndm)  :: start       !! First block to write
    INTEGER                 :: ndim        !! Dimensions of forcing to be added to the netCDF
    INTEGER,DIMENSION(ndm)  :: count_force !! Number of elements in each dimension  
    INTEGER                 :: vid         !! Variable identifer of netCDF
!_ ================================================================================================================================
    
  !! 1. Determine number of blocks of forcing variables that are stored in memory

    nblocks = 0
    ifirst(:) = 1
    ilast(:) = 1
    DO ii = ibeg, iend
       IF (     (nblocks /= 0) &
            &      .AND.(isf(ii) == isf(ilast(nblocks))+1)) THEN
          ! Last block found
          ilast(nblocks) = ii
       ELSE
          ! First block found
          nblocks = nblocks+1
          IF (nblocks > 2)  STOP 'Problem in forcing_write'
          ifirst(nblocks) = ii
          ilast(nblocks) = ii
       ENDIF
    ENDDO

  !! 2. Gather distributed variables (parallel computing)

    CALL gather(clay_fm,clay_fm_g)
    CALL gather(silt_fm,silt_fm_g) 
    CALL gather(bulk_fm,bulk_fm_g) 
    CALL gather(humrel_daily_fm,humrel_daily_fm_g)
    CALL gather(litterhum_daily_fm,litterhum_daily_fm_g)
    CALL gather(t2m_daily_fm,t2m_daily_fm_g)
    CALL gather(t2m_min_daily_fm,t2m_min_daily_fm_g)
    CALL gather(tsurf_daily_fm,tsurf_daily_fm_g)
    CALL gather(tsoil_daily_fm,tsoil_daily_fm_g)
    CALL gather(soilhum_daily_fm,soilhum_daily_fm_g)
    CALL gather(precip_fm,precip_fm_g)
    CALL gather(gpp_daily_fm,gpp_daily_fm_g)
    CALL gather(veget_fm,veget_fm_g)
    CALL gather(veget_max_fm,veget_max_fm_g)
    CALL gather(lai_fm,lai_fm_g)
    CALL gather(drainage_fm,drainage_fm_g)  
    CALL gather(runoff_fm,runoff_fm_g)  
    CALL gather(max_eau_var_daily_fm,max_eau_var_daily_fm_g)
 !! 3. Append data to netCDF file
   
    IF (is_root_prc) THEN
       ! The netCDF file has been created earlier in this module, a file ID is available 
       ! and variables and dimensions have already been defined
       DO iblocks = 1, nblocks
          IF (ifirst(iblocks) /= ilast(iblocks)) THEN
             ndim = 2
             start(1:ndim) = 1; start(ndim) = isf(ifirst(iblocks));
             count_force(1:ndim) = SHAPE(clay_fm_g)
             count_force(ndim) = isf(ilast(iblocks))-isf(ifirst(iblocks))+1
             ier = NF90_INQ_VARID (forcing_id,'clay',vid)
             ier = NF90_PUT_VAR (forcing_id,vid, &
                  &              clay_fm_g(:,ifirst(iblocks):ilast(iblocks)), &
                  & start=start(1:ndim), count=count_force(1:ndim))
             ndim = 2 
             start(1:ndim) = 1; start(ndim) = isf(ifirst(iblocks)); 
             count_force(1:ndim) = SHAPE(silt_fm_g) 
             count_force(ndim) = isf(ilast(iblocks))-isf(ifirst(iblocks))+1 
             ier = NF90_INQ_VARID (forcing_id,'silt',vid) 
             ier = NF90_PUT_VAR (forcing_id,vid, & 
                  &              silt_fm_g(:,ifirst(iblocks):ilast(iblocks)), & 
                  & start=start(1:ndim), count=count_force(1:ndim)) 
             ndim = 2 
             start(1:ndim) = 1; start(ndim) = isf(ifirst(iblocks)); 
             count_force(1:ndim) = SHAPE(bulk_fm_g) 
             count_force(ndim) = isf(ilast(iblocks))-isf(ifirst(iblocks))+1 
             ier = NF90_INQ_VARID (forcing_id,'bulk',vid) 
             ier = NF90_PUT_VAR (forcing_id,vid, & 
                  &              bulk_fm_g(:,ifirst(iblocks):ilast(iblocks)), & 
                  & start=start(1:ndim), count=count_force(1:ndim)) 

             ndim = 3;
             start(1:ndim) = 1; start(ndim) = isf(ifirst(iblocks));
             count_force(1:ndim) = SHAPE(humrel_daily_fm_g)
             count_force(ndim) = isf(ilast(iblocks))-isf(ifirst(iblocks))+1
             ier = NF90_INQ_VARID (forcing_id,'humrel',vid)
             ier = NF90_PUT_VAR (forcing_id, vid, &
                  &            humrel_daily_fm_g(:,:,ifirst(iblocks):ilast(iblocks)), &
                  &            start=start(1:ndim), count=count_force(1:ndim))
             ndim = 2;
             start(1:ndim) = 1; start(ndim) = isf(ifirst(iblocks));
             count_force(1:ndim) = SHAPE(litterhum_daily_fm_g)
             count_force(ndim) = isf(ilast(iblocks))-isf(ifirst(iblocks))+1
             ier = NF90_INQ_VARID (forcing_id,'litterhum',vid)
             ier = NF90_PUT_VAR (forcing_id, vid, &
                  &            litterhum_daily_fm_g(:,ifirst(iblocks):ilast(iblocks)), &
                  & start=start(1:ndim), count=count_force(1:ndim))
             ndim = 2;
             start(1:ndim) = 1; start(ndim) = isf(ifirst(iblocks));
             count_force(1:ndim) = SHAPE(t2m_daily_fm_g)
             count_force(ndim) = isf(ilast(iblocks))-isf(ifirst(iblocks))+1
             ier = NF90_INQ_VARID (forcing_id,'t2m',vid)
             ier = NF90_PUT_VAR (forcing_id, vid, &
                  &            t2m_daily_fm_g(:,ifirst(iblocks):ilast(iblocks)), &
                  & start=start(1:ndim), count=count_force(1:ndim))
             ndim = 2;
             start(1:ndim) = 1; start(ndim) = isf(ifirst(iblocks));
             count_force(1:ndim) = SHAPE(t2m_min_daily_fm_g)
             count_force(ndim) = isf(ilast(iblocks))-isf(ifirst(iblocks))+1
             ier = NF90_INQ_VARID (forcing_id,'t2m_min',vid)
             ier = NF90_PUT_VAR (forcing_id, vid, &
                  &            t2m_min_daily_fm_g(:,ifirst(iblocks):ilast(iblocks)), &
                  & start=start(1:ndim), count=count_force(1:ndim))
             ndim = 2;
             start(1:ndim) = 1; start(ndim) = isf(ifirst(iblocks));
             count_force(1:ndim) = SHAPE(tsurf_daily_fm_g)
             count_force(ndim) = isf(ilast(iblocks))-isf(ifirst(iblocks))+1
             ier = NF90_INQ_VARID (forcing_id,'tsurf',vid)
             ier = NF90_PUT_VAR (forcing_id, vid, &
                  &            tsurf_daily_fm_g(:,ifirst(iblocks):ilast(iblocks)), &
                  & start=start(1:ndim), count=count_force(1:ndim))
             ndim = 3;
             start(1:ndim) = 1; start(ndim) = isf(ifirst(iblocks));
             count_force(1:ndim) = SHAPE(tsoil_daily_fm_g)
             count_force(ndim) = isf(ilast(iblocks))-isf(ifirst(iblocks))+1
             ier = NF90_INQ_VARID (forcing_id,'tsoil',vid)
             ier = NF90_PUT_VAR (forcing_id, vid, &
                  &            tsoil_daily_fm_g(:,:,ifirst(iblocks):ilast(iblocks)), &
                  & start=start(1:ndim), count=count_force(1:ndim))
             ndim = 3;
             start(1:ndim) = 1; start(ndim) = isf(ifirst(iblocks));
             count_force(1:ndim) = SHAPE(soilhum_daily_fm_g)
             count_force(ndim) = isf(ilast(iblocks))-isf(ifirst(iblocks))+1
             ier = NF90_INQ_VARID (forcing_id,'soilhum',vid)
             ier = NF90_PUT_VAR (forcing_id, vid, &
                  &            soilhum_daily_fm_g(:,:,ifirst(iblocks):ilast(iblocks)), &
                  & start=start(1:ndim), count=count_force(1:ndim))
             ndim = 2;
             start(1:ndim) = 1; start(ndim) = isf(ifirst(iblocks));
             count_force(1:ndim) = SHAPE(precip_fm_g)
             count_force(ndim) = isf(ilast(iblocks))-isf(ifirst(iblocks))+1
             ier = NF90_INQ_VARID (forcing_id,'precip',vid)
             ier = NF90_PUT_VAR (forcing_id, vid, &
                  &            precip_fm_g(:,ifirst(iblocks):ilast(iblocks)), &
                  & start=start(1:ndim), count=count_force(1:ndim))
             ndim = 3;
             start(1:ndim) = 1; start(ndim) = isf(ifirst(iblocks));
             count_force(1:ndim) = SHAPE(gpp_daily_fm_g)
             count_force(ndim) = isf(ilast(iblocks))-isf(ifirst(iblocks))+1
             ier = NF90_INQ_VARID (forcing_id,'gpp',vid)
             ier = NF90_PUT_VAR (forcing_id, vid, &
                  &            gpp_daily_fm_g(:,:,ifirst(iblocks):ilast(iblocks)), &
                  &            start=start(1:ndim), count=count_force(1:ndim))
             ndim = 3;
             start(1:ndim) = 1; start(ndim) = isf(ifirst(iblocks));
             count_force(1:ndim) = SHAPE(veget_fm_g)
             count_force(ndim) = isf(ilast(iblocks))-isf(ifirst(iblocks))+1
             ier = NF90_INQ_VARID (forcing_id,'veget',vid)
             ier = NF90_PUT_VAR (forcing_id, vid, &
                  &            veget_fm_g(:,:,ifirst(iblocks):ilast(iblocks)), &
                  &            start=start(1:ndim), count=count_force(1:ndim))
             ndim = 3;
             start(1:ndim) = 1; start(ndim) = isf(ifirst(iblocks));
             count_force(1:ndim) = SHAPE(veget_max_fm_g)
             count_force(ndim) = isf(ilast(iblocks))-isf(ifirst(iblocks))+1
             ier = NF90_INQ_VARID (forcing_id,'veget_max',vid)
             ier = NF90_PUT_VAR (forcing_id, vid, &
                  &            veget_max_fm_g(:,:,ifirst(iblocks):ilast(iblocks)), &
                  &            start=start(1:ndim), count=count_force(1:ndim))
             ndim = 3;
             start(1:ndim) = 1; start(ndim) = isf(ifirst(iblocks));
             count_force(1:ndim) = SHAPE(lai_fm_g)
             count_force(ndim) = isf(ilast(iblocks))-isf(ifirst(iblocks))+1
             ier = NF90_INQ_VARID (forcing_id,'lai',vid)
             ier = NF90_PUT_VAR (forcing_id, vid, &
                  &            lai_fm_g(:,:,ifirst(iblocks):ilast(iblocks)), &
                  &            start=start(1:ndim), count=count_force(1:ndim))

             ndim = 2; 
             start(1:ndim) = 1; start(ndim) = isf(ifirst(iblocks)); 
             count_force(1:ndim) = SHAPE(drainage_fm_g) 
             count_force(ndim) = isf(ilast(iblocks))-isf(ifirst(iblocks))+1 
             ier = NF90_INQ_VARID (forcing_id,'drainage',vid) 
             ier = NF90_PUT_VAR (forcing_id, vid, & 
                  &            drainage_fm_g(:,:,ifirst(iblocks):ilast(iblocks)), & 
                  &            start=start(1:ndim), count=count_force(1:ndim)) 

             ndim = 2; 
             start(1:ndim) = 1; start(ndim) = isf(ifirst(iblocks)); 
             count_force(1:ndim) = SHAPE(runoff_fm_g) 
             count_force(ndim) = isf(ilast(iblocks))-isf(ifirst(iblocks))+1 
             ier = NF90_INQ_VARID (forcing_id,'runoff',vid) 
             ier = NF90_PUT_VAR (forcing_id, vid, & 
                  &            runoff_fm_g(:,:,ifirst(iblocks):ilast(iblocks)), & 
                  &            start=start(1:ndim), count=count_force(1:ndim)) 

             ndim = 2;
             start(1:ndim) = 1; start(ndim) = isf(ifirst(iblocks));
             count_force(1:ndim) = SHAPE(max_eau_var_daily_fm_g)
             count_force(ndim) = isf(ilast(iblocks))-isf(ifirst(iblocks))+1
             ier = NF90_INQ_VARID (forcing_id,'humrel',vid)
             ier = NF90_PUT_VAR (forcing_id, vid, &
                  &            max_eau_var_daily_fm_g(:,ifirst(iblocks):ilast(iblocks)), &
                  &            start=start(1:ndim), count=count_force(1:ndim))


          ENDIF
       ENDDO
    ENDIF
    
  !! 4. Adjust flag of forcing file
    nf_written(isf(:)) = .TRUE.

  END SUBROUTINE forcing_write

  
!! ================================================================================================================================
!! SUBROUTINE 	: stomate_forcing_read
!!
!>\BRIEF        Read forcing file.
!!
!! DESCRIPTION  : None
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S): None 
!!
!! REFERENCES	: None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
  
  SUBROUTINE stomate_forcing_read(forcing_id,nsfm)
   
  !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER,INTENT(in)  :: forcing_id           !! File identifer of forcing file, assigned when netcdf is created
    INTEGER,INTENT(in)  :: nsfm                 !! Number of time steps stored in memory        
    
    !! 0.2 Output variables

    !! 0.3 Modified variables

    !! 0.4 Local variables

    INTEGER                 :: ii                !! Index of isf where isf is the number of time steps that can be stored in 
                                                        !! memory 
    INTEGER                 :: iblocks           !! Index of block that is written
    INTEGER                 :: nblocks           !! Number of blocks that needs to be written
    INTEGER                 :: ier               !! Check error of netcdf call
    INTEGER,DIMENSION(0:2)  :: ifirst            !! First block in memory - changes with iblocks
    INTEGER,DIMENSION(0:2)  :: ilast             !! Last block in memory - changes with iblocks
    INTEGER,PARAMETER       :: ndm = 10          !! Maximum number of dimensions
    INTEGER,DIMENSION(ndm)  :: start             !! First block to write
    INTEGER                 :: ndim              !! Dimensions of forcing to be added to the netCDF
    INTEGER,DIMENSION(ndm)  :: count_force       !! Number of elements in each dimension
    INTEGER                 :: vid               !! Variable identifer of netCDF
    LOGICAL                        :: a_er=.FALSE.      !! Error catching from netcdf file
!_ ================================================================================================================================

    IF (printlev >= 4) WRITE(numout,*) "stomate_forcing_read "
    
  !! 1. Set to zero if the corresponding forcing state

    ! has not yet been written into the file  
    DO ii = 1, nsfm
       IF (.NOT.nf_written(isf(ii))) THEN
          clay_fm(:,ii) = zero
          silt_fm(:,iisf) = zero 
          bulk_fm(:,iisf) = zero 
          humrel_daily_fm(:,:,ii) = zero
          litterhum_daily_fm(:,ii) = zero
          t2m_daily_fm(:,ii) = zero
          t2m_min_daily_fm(:,ii) = zero
          tsurf_daily_fm(:,ii) = zero
          tsoil_daily_fm(:,:,ii) = zero
          soilhum_daily_fm(:,:,ii) = zero
          precip_fm(:,ii) = zero
          gpp_daily_fm(:,:,ii) = zero
          veget_fm(:,:,ii) = zero
          veget_max_fm(:,:,ii) = zero
          lai_fm(:,:,ii) = zero
          drainage_fm(:,:,iisf) = zero 
          runoff_fm  (:,:,iisf) = zero 
          max_eau_var_daily_fm(:,ii) = zero
       ENDIF
    ENDDO
    
  !! 2. determine blocks of forcing states that are contiguous in memory

    nblocks = 0
    ifirst(:) = 1
    ilast(:) = 1
    
    DO ii = 1, nsfm
       IF (nf_written(isf(ii))) THEN
          IF (     (nblocks /= 0) &
               &        .AND.(isf(ii) == isf(ilast(nblocks))+1)) THEN

             ! element is contiguous with last element found
             ilast(nblocks) = ii
          ELSE

             ! found first element of new block
             nblocks = nblocks+1
             IF (nblocks > 2)  STOP 'Problem in stomate_forcing_read'
             
             ifirst(nblocks) = ii
             ilast(nblocks) = ii
          ENDIF
       ENDIF
    ENDDO
    IF (printlev >= 4) WRITE(numout,*) "stomate_forcing_read nblocks, ifirst, ilast",nblocks, ifirst, ilast
    
  !! 3. Read variable values

    IF (is_root_prc) THEN
       DO iblocks = 1, nblocks
          IF (printlev >= 4) WRITE(numout,*) "stomate_forcing_read iblocks, ifirst(iblocks), ilast(iblocks)",iblocks, &
               ifirst(iblocks), ilast(iblocks)
          IF (ifirst(iblocks) /= ilast(iblocks)) THEN
             a_er=.FALSE.
             ndim = 2;
             start(1:ndim) = 1; start(ndim) = isf(ifirst(iblocks));
             count_force(1:ndim) = SHAPE(clay_fm_g)
             count_force(ndim) = isf(ilast(iblocks))-isf(ifirst(iblocks))+1
             ier = NF90_INQ_VARID (forcing_id,'clay',vid)
             a_er = a_er.OR.(ier /= 0)
             ier = NF90_GET_VAR (forcing_id, vid, &
                  &            clay_fm_g(:,ifirst(iblocks):ilast(iblocks)), &
                  &            start=start(1:ndim), count=count_force(1:ndim))
             a_er = a_er.OR.(ier /= 0)


             ndim = 2; 
             start(1:ndim) = 1; start(ndim) = isf(ifirst(iblocks)); 
             count_force(1:ndim) = SHAPE(silt_fm_g) 
             count_force(ndim) = isf(ilast(iblocks))-isf(ifirst(iblocks))+1 
             ier = NF90_INQ_VARID (forcing_id,'silt',vid) 
             a_er = a_er.OR.(ier /= 0) 
             ier = NF90_GET_VAR (forcing_id, vid, & 
                  &            silt_fm_g(:,ifirst(iblocks):ilast(iblocks)), & 
                  &            start=start(1:ndim), count=count_force(1:ndim)) 
             a_er = a_er.OR.(ier /= 0) 
             
             ndim = 2; 
             start(1:ndim) = 1; start(ndim) = isf(ifirst(iblocks)); 
             count_force(1:ndim) = SHAPE(bulk_fm_g) 
             count_force(ndim) = isf(ilast(iblocks))-isf(ifirst(iblocks))+1 
             ier = NF90_INQ_VARID (forcing_id,'bulk',vid) 
             a_er = a_er.OR.(ier /= 0) 
             ier = NF90_GET_VAR (forcing_id, vid, & 
                  &            bulk_fm_g(:,ifirst(iblocks):ilast(iblocks)), & 
                  &            start=start(1:ndim), count=count_force(1:ndim)) 
             a_er = a_er.OR.(ier /= 0) 
             
             ndim = 3;
             start(1:ndim) = 1; start(ndim) = isf(ifirst(iblocks));
             count_force(1:ndim) = SHAPE(humrel_daily_fm_g)
             count_force(ndim) = isf(ilast(iblocks))-isf(ifirst(iblocks))+1
             ier = NF90_INQ_VARID (forcing_id,'humrel',vid)
             a_er = a_er.OR.(ier /= 0)
             ier = NF90_GET_VAR (forcing_id, vid, &
                  &            humrel_daily_fm_g(:,:,ifirst(iblocks):ilast(iblocks)), &
                  &            start=start(1:ndim), count=count_force(1:ndim))
             a_er = a_er.OR.(ier /= 0)

             ndim = 2;
             start(1:ndim) = 1; start(ndim) = isf(ifirst(iblocks));
             count_force(1:ndim) = SHAPE(litterhum_daily_fm_g)
             count_force(ndim) = isf(ilast(iblocks))-isf(ifirst(iblocks))+1
             ier = NF90_INQ_VARID (forcing_id,'litterhum',vid)
             a_er = a_er.OR.(ier /= 0)
             ier = NF90_GET_VAR (forcing_id, vid, &
                  &              litterhum_daily_fm_g(:,ifirst(iblocks):ilast(iblocks)), &
                  &            start=start(1:ndim), count=count_force(1:ndim))
             a_er = a_er.OR.(ier /= 0)

             ndim = 2;
             start(1:ndim) = 1; start(ndim) = isf(ifirst(iblocks));
             count_force(1:ndim) = SHAPE(t2m_daily_fm_g)
             count_force(ndim) = isf(ilast(iblocks))-isf(ifirst(iblocks))+1
             ier = NF90_INQ_VARID (forcing_id,'t2m',vid)
             a_er = a_er.OR.(ier /= 0)
             ier = NF90_GET_VAR (forcing_id, vid, &
                  &              t2m_daily_fm_g(:,ifirst(iblocks):ilast(iblocks)), &
                  &            start=start(1:ndim), count=count_force(1:ndim))
             a_er = a_er.OR.(ier /= 0)

             ndim = 2;
             start(1:ndim) = 1; start(ndim) = isf(ifirst(iblocks));
             count_force(1:ndim) = SHAPE(t2m_min_daily_fm_g)
             count_force(ndim) = isf(ilast(iblocks))-isf(ifirst(iblocks))+1
             ier = NF90_INQ_VARID (forcing_id,'t2m_min',vid)
             a_er = a_er.OR.(ier /= 0)
             ier = NF90_GET_VAR (forcing_id, vid, &
                  &              t2m_min_daily_fm_g(:,ifirst(iblocks):ilast(iblocks)), &
                  &            start=start(1:ndim), count=count_force(1:ndim))
             a_er = a_er.OR.(ier /= 0)

             ndim = 2;
             start(1:ndim) = 1; start(ndim) = isf(ifirst(iblocks));
             count_force(1:ndim) = SHAPE(tsurf_daily_fm_g)
             count_force(ndim) = isf(ilast(iblocks))-isf(ifirst(iblocks))+1
             ier = NF90_INQ_VARID (forcing_id,'tsurf',vid)
             a_er = a_er.OR.(ier /= 0)
             ier = NF90_GET_VAR (forcing_id, vid, &
                  &              tsurf_daily_fm_g(:,ifirst(iblocks):ilast(iblocks)), &
                  &            start=start(1:ndim), count=count_force(1:ndim))
             a_er = a_er.OR.(ier /= 0)

             ndim = 3;
             start(1:ndim) = 1; start(ndim) = isf(ifirst(iblocks));
             count_force(1:ndim) = SHAPE(tsoil_daily_fm_g)
             count_force(ndim) = isf(ilast(iblocks))-isf(ifirst(iblocks))+1
             ier = NF90_INQ_VARID (forcing_id,'tsoil',vid)
             a_er = a_er.OR.(ier /= 0)
             ier = NF90_GET_VAR (forcing_id, vid, &
                  &              tsoil_daily_fm_g(:,:,ifirst(iblocks):ilast(iblocks)), &
                  &            start=start(1:ndim), count=count_force(1:ndim))
             a_er = a_er.OR.(ier /= 0)

             ndim = 3;
             start(1:ndim) = 1; start(ndim) = isf(ifirst(iblocks));
             count_force(1:ndim) = SHAPE(soilhum_daily_fm_g)
             count_force(ndim) = isf(ilast(iblocks))-isf(ifirst(iblocks))+1
             ier = NF90_INQ_VARID (forcing_id,'soilhum',vid)
             a_er = a_er.OR.(ier /= 0)
             ier = NF90_GET_VAR (forcing_id, vid, &
                  &              soilhum_daily_fm_g(:,:,ifirst(iblocks):ilast(iblocks)), &
                  &            start=start(1:ndim), count=count_force(1:ndim))
             a_er = a_er.OR.(ier /= 0)

             ndim = 2;
             start(1:ndim) = 1; start(ndim) = isf(ifirst(iblocks));
             count_force(1:ndim) = SHAPE(precip_fm_g)
             count_force(ndim) = isf(ilast(iblocks))-isf(ifirst(iblocks))+1
             ier = NF90_INQ_VARID (forcing_id,'precip',vid)
             a_er = a_er.OR.(ier /= 0)
             ier = NF90_GET_VAR (forcing_id, vid, &
                  &              precip_fm_g(:,ifirst(iblocks):ilast(iblocks)), &
                  &            start=start(1:ndim), count=count_force(1:ndim))
             a_er = a_er.OR.(ier /= 0)

             ndim = 3;
             start(1:ndim) = 1; start(ndim) = isf(ifirst(iblocks));
             count_force(1:ndim) = SHAPE(gpp_daily_fm_g)
             count_force(ndim) = isf(ilast(iblocks))-isf(ifirst(iblocks))+1
             ier = NF90_INQ_VARID (forcing_id,'gpp',vid)
             a_er = a_er.OR.(ier /= 0)
             ier = NF90_GET_VAR (forcing_id, vid, &
                  &            gpp_daily_fm_g(:,:,ifirst(iblocks):ilast(iblocks)), &
                  &            start=start(1:ndim), count=count_force(1:ndim))
             a_er = a_er.OR.(ier /= 0)

             ndim = 3;
             start(1:ndim) = 1; start(ndim) = isf(ifirst(iblocks));
             count_force(1:ndim) = SHAPE(veget_fm_g)
             count_force(ndim) = isf(ilast(iblocks))-isf(ifirst(iblocks))+1
             ier = NF90_INQ_VARID (forcing_id,'veget',vid)
             a_er = a_er.OR.(ier /= 0)
             ier = NF90_GET_VAR (forcing_id, vid, &
                  &            veget_fm_g(:,:,ifirst(iblocks):ilast(iblocks)), &
                  &            start=start(1:ndim), count=count_force(1:ndim))
             a_er = a_er.OR.(ier /= 0)

             ndim = 3;
             start(1:ndim) = 1; start(ndim) = isf(ifirst(iblocks));
             count_force(1:ndim) = SHAPE(veget_max_fm_g)
             count_force(ndim) = isf(ilast(iblocks))-isf(ifirst(iblocks))+1
             ier = NF90_INQ_VARID (forcing_id,'veget_max',vid)
             a_er = a_er.OR.(ier /= 0)
             ier = NF90_GET_VAR (forcing_id, vid, &
                  &            veget_max_fm_g(:,:,ifirst(iblocks):ilast(iblocks)), &
                  &            start=start(1:ndim), count=count_force(1:ndim))
             a_er = a_er.OR.(ier /= 0)

             ndim = 3;
             start(1:ndim) = 1; start(ndim) = isf(ifirst(iblocks));
             count_force(1:ndim) = SHAPE(lai_fm_g)
             count_force(ndim) = isf(ilast(iblocks))-isf(ifirst(iblocks))+1
             ier = NF90_INQ_VARID (forcing_id,'lai',vid)
             a_er = a_er.OR.(ier /= 0)
             ier = NF90_GET_VAR (forcing_id, vid, &
                  &            lai_fm_g(:,:,ifirst(iblocks):ilast(iblocks)), &
                  &            start=start(1:ndim), count=count_force(1:ndim))
             a_er = a_er.OR.(ier /= 0)


             ndim = 2; 
             start(1:ndim) = 1; start(ndim) = isf(ifirst(iblocks)); 
             count_force(1:ndim) = SHAPE(drainage_fm_g) 
             count_force(ndim) = isf(ilast(iblocks))-isf(ifirst(iblocks))+1 
             ier = NF90_INQ_VARID (forcing_id,'drainage',vid) 
             a_er = a_er.OR.(ier /= 0) 
             ier = NF90_GET_VAR (forcing_id, vid, & 
                  &            drainage_fm_g(:,:,ifirst(iblocks):ilast(iblocks)), & 
                  &            start=start(1:ndim), count=count_force(1:ndim)) 
             a_er = a_er.OR.(ier /= 0) 

             ndim = 2; 
             start(1:ndim) = 1; start(ndim) = isf(ifirst(iblocks)); 
             count_force(1:ndim) = SHAPE(runoff_fm_g) 
             count_force(ndim) = isf(ilast(iblocks))-isf(ifirst(iblocks))+1 
             ier = NF90_INQ_VARID (forcing_id,'runoff',vid) 
             a_er = a_er.OR.(ier /= 0) 
             ier = NF90_GET_VAR (forcing_id, vid, & 
                  &            runoff_fm_g(:,:,ifirst(iblocks):ilast(iblocks)), & 
                  &            start=start(1:ndim), count=count_force(1:ndim)) 
             a_er = a_er.OR.(ier /= 0) 

             ndim = 2;
             start(1:ndim) = 1; start(ndim) = isf(ifirst(iblocks));
             count_force(1:ndim) = SHAPE(max_eau_var_daily_fm_g)
             count_force(ndim) = isf(ilast(iblocks))-isf(ifirst(iblocks))+1
             ier = NF90_INQ_VARID (forcing_id,'max_eau_var',vid)
             a_er = a_er.OR.(ier /= 0)
             ier = NF90_GET_VAR (forcing_id, vid, &
                  &            max_eau_var_daily_fm_g(:,ifirst(iblocks):ilast(iblocks)), &
                  &            start=start(1:ndim), count=count_force(1:ndim))
             a_er = a_er.OR.(ier /= 0)

             IF (a_er) THEN
                CALL ipslerr_p (3,'stomate_forcing_read', &
                     &        'PROBLEM when read forcing file', &
                     &        '','')
             ENDIF

          ENDIF ! (ifirst(iblocks) /= ilast(iblocks))
       ENDDO ! iblocks
    ENDIF ! is_root_prc

  !! 4. Distribute the variable over several processors

    CALL scatter(clay_fm_g,clay_fm)
    CALL scatter(silt_fm_g,silt_fm) 
    CALL scatter(bulk_fm_g,bulk_fm) 
    CALL scatter(humrel_daily_fm_g,humrel_daily_fm)
    CALL scatter(litterhum_daily_fm_g,litterhum_daily_fm)
    CALL scatter(t2m_daily_fm_g,t2m_daily_fm)
    CALL scatter(t2m_min_daily_fm_g,t2m_min_daily_fm)
    CALL scatter(tsurf_daily_fm_g,tsurf_daily_fm)
    CALL scatter(tsoil_daily_fm_g,tsoil_daily_fm)
    CALL scatter(soilhum_daily_fm_g,soilhum_daily_fm)
    CALL scatter(precip_fm_g,precip_fm)
    CALL scatter(gpp_daily_fm_g,gpp_daily_fm)
    CALL scatter(veget_fm_g,veget_fm)
    CALL scatter(veget_max_fm_g,veget_max_fm)
    CALL scatter(lai_fm_g,lai_fm)
    CALL scatter(drainage_fm_g,drainage_fm) 
    CALL scatter(runoff_fm_g,runoff_fm) 
    CALL scatter(max_eau_var_daily_fm_g,max_eau_var_daily_fm)
  END SUBROUTINE stomate_forcing_read


!! ================================================================================================================================
!! SUBROUTINE 	: setlai
!!
!>\BRIEF        Routine to force the lai in STOMATE. The code in this routine
!! simply CALCULATES lai and is therefore not functional. The routine should be 
!! rewritten if one wants to force lai.
!!
!! DESCRIPTION  : None
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S): ::lai
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART : None
!! \n
!_ ================================================================================================================================
  
  SUBROUTINE setlai(npts,lai)

  !! 0 Variable and parameter declaration 
  
    !! 0.1 Input variables

    INTEGER,INTENT(in)                    :: npts !! Domain size - number of pixels (unitless)
    
    !! 0.2 Output variables

    REAL,DIMENSION(npts,nvm),INTENT(out)  :: lai  !! PFT leaf area index @tex $(m^{2} m^{-2})$ @endtex

    !! 0.3 Modified variables

    !! 0.4 Local variables

    INTEGER                               :: j    !! index (unitless)
!_ ================================================================================================================================
    
    !! 1. Set lai for bare soil to zero

    lai(:,ibare_sechiba) = zero

    !! 2. Multiply foliage biomass by sla to calculate lai for all PFTs and pixels

    DO j=2,nvm
       lai(:,j) = biomass(:,j,ileaf,icarbon)*sla(j)
    ENDDO
    
  END SUBROUTINE setlai

END MODULE stomate
