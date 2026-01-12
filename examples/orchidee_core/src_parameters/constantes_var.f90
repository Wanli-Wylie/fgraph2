! =================================================================================================================================
! MODULE       : constantes_var
!
! CONTACT      : orchidee-help _at_ listes.ipsl.fr
!
! LICENCE      : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF        constantes_var module contains most constantes like pi, Earth radius, etc...
!!              and all externalized parameters except pft-dependent constants.
!!
!!\n DESCRIPTION: This module contains most constantes and the externalized parameters of ORCHIDEE which 
!!                are not pft-dependent.\n
!!                In this module, you can set the flag diag_qsat in order to detect the pixel where the
!!                temperature is out of range (look qsatcalc and dev_qsatcalc in qsat_moisture.f90).\n
!!                The Earth radius is approximated by the Equatorial radius.The Earth's equatorial radius a,
!!                or semi-major axis, is the distance from its center to the equator and equals 6,378.1370 km.
!!                The equatorial radius is often used to compare Earth with other planets.\n
!!                The meridional mean is well approximated by the semicubic mean of the two axe yielding 
!!                6367.4491 km or less accurately by the quadratic mean of the two axes about 6,367.454 km
!!                or even just the mean of the two axes about 6,367.445 km.\n
!!                This module is already USE in module constantes. Therefor no need to USE it seperatly except
!!                if the subroutines in module constantes are not needed.\n
!!                
!! RECENT CHANGE(S):
!!
!! REFERENCE(S)	: 
!! - Louis, Jean-Francois (1979), A parametric model of vertical eddy fluxes in the atmosphere. 
!! Boundary Layer Meteorology, 187-202.\n
!!
!! SVN          :
!! $HeadURL: $
!! $Date: 2021-05-06 11:05:55 +0200 (四, 2021-05-06) $
!! $Revision: 7177 $
!! \n
!_ ================================================================================================================================

MODULE constantes_var

  USE defprec

  IMPLICIT NONE
!-

                         !-----------------------!
                         !  ORCHIDEE CONSTANTS   !
                         !-----------------------!

  !
  ! FLAGS 
  !
 LOGICAL, SAVE :: NC_COMPRESSION_ENABLE !! activate netcdf output compression 
!$OMP THREADPRIVATE(NC_COMPRESSION_ENABLE) 

  LOGICAL :: river_routing      !! activate river routing
!$OMP THREADPRIVATE(river_routing)
  LOGICAL :: hydrol_cwrr        !! activate 11 layers hydrolgy model
!$OMP THREADPRIVATE(hydrol_cwrr)
  LOGICAL, SAVE :: ok_nudge_mc  !! Activate nudging of soil moisture 
!$OMP THREADPRIVATE(ok_nudge_mc)
  LOGICAL, SAVE :: ok_nudge_snow!! Activate nudging of snow variables
!$OMP THREADPRIVATE(ok_nudge_snow)
  LOGICAL, SAVE :: nudge_interpol_with_xios  !! Activate reading and interpolation with XIOS for nudging fields
!$OMP THREADPRIVATE(nudge_interpol_with_xios)
  LOGICAL :: do_floodplains     !! activate flood plains
!$OMP THREADPRIVATE(do_floodplains)
  LOGICAL :: do_irrigation      !! activate computation of irrigation flux
!$OMP THREADPRIVATE(do_irrigation)
  LOGICAL :: ok_sechiba         !! activate physic of the model
!$OMP THREADPRIVATE(ok_sechiba)
  LOGICAL :: ok_co2             !! activate photosynthesis
!$OMP THREADPRIVATE(ok_co2)
  LOGICAL :: ok_stomate         !! activate carbon cycle
!$OMP THREADPRIVATE(ok_stomate)
  LOGICAL :: ok_ncycle          !! activate nitrogen cycle (true/false)
!$OMP THREADPRIVATE(ok_ncycle)
  LOGICAL :: ok_pcycle          !! activate phosphorus cycle (true/false)
!$OMP THREADPRIVATE(ok_pcycle)
  LOGICAL :: impose_cn          !! impose the CN ratio of leaves (true/false)
!$OMP THREADPRIVATE(impose_cn)
  LOGICAL :: impose_np          !! impose the NP ratio of leaves (true/false)
!$OMP THREADPRIVATE(impose_np)
  LOGICAL :: impose_cn_restart  !! impose the CN ratio of leaves from restart (true/false)
!$OMP THREADPRIVATE(impose_cn_restart)
  LOGICAL :: impose_np_restart  !! impose the NP ratio of leaves from restart (true/false)
!$OMP THREADPRIVATE(impose_np_restart)
  LOGICAL :: ok_dgvm            !! activate dynamic vegetation
!$OMP THREADPRIVATE(ok_dgvm)
  LOGICAL :: do_wood_harvest    !! activate wood harvest
!$OMP THREADPRIVATE(do_wood_harvest)
  LOGICAL :: ok_pheno           !! activate the calculation of lai using stomate rather than a prescription
!$OMP THREADPRIVATE(ok_pheno)
  LOGICAL :: ok_bvoc            !! activate biogenic volatile organic coumpounds
!$OMP THREADPRIVATE(ok_bvoc)
  LOGICAL :: ok_leafage         !! activate leafage
!$OMP THREADPRIVATE(ok_leafage)
  LOGICAL :: ok_radcanopy       !! use canopy radiative transfer model
!$OMP THREADPRIVATE(ok_radcanopy)
  LOGICAL :: ok_multilayer      !! use canopy radiative transfer model with multi-layers
!$OMP THREADPRIVATE(ok_multilayer)
  LOGICAL :: ok_pulse_NOx       !! calculate NOx emissions with pulse
!$OMP THREADPRIVATE(ok_pulse_NOx)
  LOGICAL :: ok_bbgfertil_NOx   !! calculate NOx emissions with bbg fertilizing effect
!$OMP THREADPRIVATE(ok_bbgfertil_NOx)
  LOGICAL :: ok_cropsfertil_NOx !! calculate NOx emissions with fertilizers use
!$OMP THREADPRIVATE(ok_cropsfertil_NOx)

  LOGICAL :: mass_conservation  !! crash the model at any violation of mass conservation; at the moment the phenology and prescribe are not tested as they designed to violate mass
  LOGICAL :: dsg_debug          !! test mass conservation, negative biomass, and stoichiometry


  LOGICAL :: ok_co2bvoc_poss    !! CO2 inhibition on isoprene activated following Possell et al. (2005) model
!$OMP THREADPRIVATE(ok_co2bvoc_poss)
  LOGICAL :: ok_co2bvoc_wilk    !! CO2 inhibition on isoprene activated following Wilkinson et al. (2006) model
!$OMP THREADPRIVATE(ok_co2bvoc_wilk)
  
  LOGICAL, SAVE :: OFF_LINE_MODE = .FALSE.  !! ORCHIDEE detects if it is coupled with a GCM or 
                                            !! just use with one driver in OFF-LINE. (true/false)
!$OMP THREADPRIVATE(OFF_LINE_MODE)  
  LOGICAL, SAVE :: impose_param = .TRUE.    !! Flag impos_param : read all the parameters in the run.def file
!$OMP THREADPRIVATE(impose_param)
  CHARACTER(LEN=80), SAVE     :: restname_in       = 'NONE'                 !! Input Restart files name for Sechiba component  
!$OMP THREADPRIVATE(restname_in)
  CHARACTER(LEN=80), SAVE     :: restname_out      = 'sechiba_rest_out.nc'  !! Output Restart files name for Sechiba component
!$OMP THREADPRIVATE(restname_out)
  CHARACTER(LEN=80), SAVE     :: stom_restname_in  = 'NONE'                 !! Input Restart files name for Stomate component
!$OMP THREADPRIVATE(stom_restname_in)
  CHARACTER(LEN=80), SAVE     :: stom_restname_out = 'stomate_rest_out.nc'  !! Output Restart files name for Stomate component
!$OMP THREADPRIVATE(stom_restname_out)
  INTEGER, SAVE :: printlev=2       !! Standard level for text output [0, 1, 2, 3]
!$OMP THREADPRIVATE(printlev)

  !
  ! TIME
  !
  INTEGER, PARAMETER  :: spring_days_max = 40  !! Maximum number of days during which we watch for possible spring frost damage
  !
  ! SPECIAL VALUES 
  !
  INTEGER, PARAMETER :: undef_int = 999999999     !! undef integer for integer arrays (unitless)
  !-
  REAL, SAVE :: val_exp = 999999.                 !! Specific value if no restart value  (unitless)
!$OMP THREADPRIVATE(val_exp)
  REAL, PARAMETER :: undef = -9999.               !! Special value for stomate (unitless)
  
  REAL, PARAMETER :: min_sechiba = 1.E-8_r_std    !! Epsilon to detect a near zero floating point (unitless)
  REAL, PARAMETER :: undef_sechiba = 1.E+20_r_std !! The undef value used in SECHIBA (unitless)
  
  REAL, PARAMETER :: min_stomate = 1.E-8_r_std    !! Epsilon to detect a near zero floating point (unitless)
  REAL, PARAMETER :: large_value = 1.E33_r_std    !! some large value (for stomate) (unitless)

  REAL, SAVE :: alpha_nudge_mc                    !! Nudging constant for soil moisture 
!$OMP THREADPRIVATE(alpha_nudge_mc)
  REAL, SAVE :: alpha_nudge_snow                  !! Nudging constant for snow variables
!$OMP THREADPRIVATE(alpha_nudge_snow)

  !
  !  DIMENSIONING AND INDICES PARAMETERS  
  !
  INTEGER, PARAMETER :: ibare_sechiba = 1 !! Index for bare soil in Sechiba (unitless)
  INTEGER, PARAMETER :: ivis = 1          !! index for albedo in visible range (unitless)
  INTEGER, PARAMETER :: inir = 2          !! index for albeod i near-infrared range (unitless) 
  INTEGER, PARAMETER :: nnobio = 1        !! Number of other surface types: land ice (lakes,cities, ...) (unitless)
  INTEGER, PARAMETER :: iice = 1          !! Index for land ice (see nnobio) (unitless)
  !-
  !! Soil
  INTEGER, PARAMETER :: classnb = 9       !! Levels of soil colour classification (unitless)
  !-
  INTEGER, PARAMETER :: nleafages = 4     !! leaf age discretisation ( 1 = no discretisation )(unitless)
  !-
  !! litter fractions: indices (unitless)
  INTEGER, PARAMETER :: ileaf = 1         !! Index for leaf compartment (unitless)
  INTEGER, PARAMETER :: isapabove = 2     !! Index for sapwood above compartment (unitless)
  INTEGER, PARAMETER :: isapbelow = 3     !! Index for sapwood below compartment (unitless)
  INTEGER, PARAMETER :: iheartabove = 4   !! Index for heartwood above compartment (unitless)
  INTEGER, PARAMETER :: iheartbelow = 5   !! Index for heartwood below compartment (unitless)
  INTEGER, PARAMETER :: iroot = 6         !! Index for roots compartment (unitless)
  INTEGER, PARAMETER :: ifruit = 7        !! Index for fruits compartment (unitless)
  INTEGER, PARAMETER :: icarbres = 8      !! Index for reserve compartment (unitless)
  INTEGER, PARAMETER :: ilabile = 9       !! Index for reserve compartment (unitless) 
  INTEGER, PARAMETER :: nparts = 9        !! Number of biomass compartments (unitless)
  !-
  !! indices for assimilation parameters 
  INTEGER, PARAMETER :: ivcmax = 1        !! Index for vcmax (assimilation parameters) (unitless)
  INTEGER, PARAMETER :: inue = 2          !! Index for nue (assimilationbn parameters) (unitless)
  INTEGER, PARAMETER :: ileafN = 3        !! Index for leaf N (assimilationbn parameters) (unitless)
  INTEGER, PARAMETER :: ijmax = 4         !! Index for jmax (assimilationbn parameters) (unitless)
!DSGpot
  INTEGER, PARAMETER :: ijmax_pot  = 5    !! Index for potential jmax (assimilationbn parameters) (unitless)
  INTEGER, PARAMETER :: ivcmax_pot = 6    !! Index for potential vcmax (assimilationbn parameters) (unitless)
  INTEGER, PARAMETER :: npco2      = 6    !! Number of assimilation parameters (unitless)
!DSGpot



  !-
  !! trees and litter: indices for the parts of heart-
  !! and sapwood above and below the ground 
  INTEGER, PARAMETER :: iabove = 1       !! Index for above part (unitless)
  INTEGER, PARAMETER :: ibelow = 2       !! Index for below part (unitless)
  INTEGER, PARAMETER :: nlevs = 2        !! Number of levels for trees and litter (unitless)
  !-
  !! litter: indices for metabolic and structural part
  INTEGER, PARAMETER :: imetabolic = 1   !! Index for metabolic litter (unitless)
  INTEGER, PARAMETER :: istructural = 2  !! Index for structural litter (unitless)
  INTEGER, PARAMETER :: iwoody = 3       !! Index for woody litter (unitless)
  INTEGER, PARAMETER :: nlitt = 3        !! Number of levels for litter compartments (unitless)
  !-
  !! carbon pools: indices
  INTEGER, PARAMETER :: iactive = 1      !! Index for active carbon pool (unitless)
  INTEGER, PARAMETER :: islow = 2        !! Index for slow carbon pool (unitless)
  INTEGER, PARAMETER :: ipassive = 3     !! Index for passive carbon pool (unitless)
  INTEGER, PARAMETER :: isurface = 4     !! Index for passive carbon pool (unitless)
  INTEGER, PARAMETER :: ncarb = 4        !! Number of soil carbon pools (unitless)
  !-
  !! For isotopes and nitrogen
  INTEGER, PARAMETER :: nelements = 3    !! Number of isotopes considered
  INTEGER, PARAMETER :: icarbon = 1      !! Index for carbon 
  INTEGER, PARAMETER :: initrogen = 2    !! Index for nitrogen 
  INTEGER, PARAMETER :: iphosphorus = 3  !! Index for phosphorus

  !! N-cycle : indices
  INTEGER, PARAMETER :: iammonium = 1    !! Index for Ammonium 
  INTEGER, PARAMETER :: initrate  = 2    !! Index for Nitrate
  INTEGER, PARAMETER :: inox      = 3    !! Index for NOX
  INTEGER, PARAMETER :: initrous  = 4    !! Index for N2O
  INTEGER, PARAMETER :: idinitro  = 5    !! Index for N2
  INTEGER, PARAMETER :: nionspec  = 2    !! Number of ionics form considered (ammonium, nitrate)
  INTEGER, PARAMETER :: nnspec    = 5    !! Number of N-species considered

  INTEGER, PARAMETER :: iatm_ammo = 1    !! Index for N input from Ammonium N atmospheric deposition
  INTEGER, PARAMETER :: iatm_nitr = 2    !! Index for N input from Nitrate N atmospheric deposition
  INTEGER, PARAMETER :: ibnf      = 3    !! Index for N input from BNF
  INTEGER, PARAMETER :: ninput    = 3    !! Number of N-input considered  

  INTEGER, PARAMETER :: i_nh4_to_no3 = 1 !! Index for NO3-production reaction 
  INTEGER, PARAMETER :: i_nh4_to_no  = 2 !! Index for NO -production reaction 
  INTEGER, PARAMETER :: i_nh4_to_n2o = 3 !! Index for N2O-production reaction 

  INTEGER, PARAMETER :: i_no3_to_n2 = 1 !! Index for NO3    consumption reaction
  INTEGER, PARAMETER :: i_no3_to_nox = 2 !! Index for NO/Nox consumption reaction
  INTEGER, PARAMETER :: i_no3_to_n2o  = 3 !! Index for N2O    consumption reaction

!  INTEGER, PARAMETER :: i_no3_to_nox = 1 !! Index for NO3    consumption reaction
!  INTEGER, PARAMETER :: i_nox_to_n2o = 2 !! Index for NO/Nox consumption reaction
!  INTEGER, PARAMETER :: i_n2o_to_n2  = 3 !! Index for N2O    consumption reaction


  !! P-cycle :  indices
  
  INTEGER, PARAMETER :: ipdissolved   = 1    !!  Index for P in solution
  INTEGER, PARAMETER :: ipsorbed      = 2    !!  Index for "sorbed" P
  ! the following fractions are sinks of P on time scale < 10,000 yr 
  !-> no need to compute their dynamics for now (Goll et al. (2012) 
  ! INTEGER, PARAMETER :: issorbed  = 3    !!  Index for strongly sorbed P
  ! INTEGER, PARAMETER :: ioccluded = 4    !!  Index for occluded P          

  INTEGER, PARAMETER :: npspec    = 2    !! Number of P-species considered

  INTEGER, PARAMETER :: iatm_p  = 1    !! Index for P input from  atmospheric deposition
  INTEGER, PARAMETER :: iweat_p = 2    !! Index for P input from chemical weathering
  INTEGER, PARAMETER :: pinput  = 2    !! Number of P-input considered  
  
  !REAL, PARAMETER    :: frac_labile_dissolved = 0.0001 !! fraction of labile P dissolved in soil solution []
  !                                                            !! value from Goll et al (2012)
  REAL, SAVE         :: tau_sorbed            = 54750.  !! turnover time of  P adsorbed to soil particles [days]
                                                               !! We choose 150yrs - estimates range from 25yr for Hawaii (Goll et al 2017) to 150yr (Wang et al 2010)
  REAL, PARAMETER    :: dissolved_init        = .01     !! initial value of P in soil solution when not read in from restart [gP/m2]
  REAL, PARAMETER    :: labile_init           = 2.      !! initial value of labile when not read in from restart [gP/m2]; value from Yang et al (2013); Figure3 plus 50% more to account for labile organic
                                                               


  !
  !! Indices used for analytical spin-up
  INTEGER, PARAMETER :: nbpools = 10              !! Total number of carbon pools (unitless)
  INTEGER, PARAMETER :: istructural_above = 1    !! Index for structural litter above (unitless)
  INTEGER, PARAMETER :: istructural_below = 2    !! Index for structural litter below (unitless)
  INTEGER, PARAMETER :: imetabolic_above = 3     !! Index for metabolic litter above (unitless)
  INTEGER, PARAMETER :: imetabolic_below = 4     !! Index for metabolic litter below (unitless)
  INTEGER, PARAMETER :: iwoody_above = 5         !! Index for woody litter above (unitless)
  INTEGER, PARAMETER :: iwoody_below = 6         !! Index for woody litter below (unitless)
  INTEGER, PARAMETER :: iactive_pool = 7         !! Index for active carbon pool (unitless)
  INTEGER, PARAMETER :: islow_pool   = 8         !! Index for slow carbon pool (unitless)
  INTEGER, PARAMETER :: ipassive_pool = 9        !! Index for passive carbon pool (unitless)
  INTEGER, PARAMETER :: isurface_pool = 10       !! Index for surface carbon pool (unitless)

  ! Enhanced Weathering
  LOGICAL, SAVE :: enhanced_weat = .FALSE.              !! Do we consider  enhanced weathering
!$OMP THREADPRIVATE(enhanced_weat)

  LOGICAL, SAVE :: EW_force_stock = .FALSE.               !! Do we force a fixed mineral stock size (instead of application rate)
!$OMP THREADPRIVATE(EW_force_stock)


  INTEGER, PARAMETER :: nminerals   = 2                 !! Number of different minerals considered (for now: basalt, forsterite)
  INTEGER, PARAMETER :: iforsterite = 1          !! Index for forsterite
  INTEGER, PARAMETER :: ibasalt     = 2              !! Index for basalt

  INTEGER, PARAMETER :: nECfluxes = 2                 !! Number of different weathering fluxes 
  INTEGER, PARAMETER :: ico2cons  = 1                 !! Index for CO2 consumption by EC
  INTEGER, PARAMETER :: iPrelease = 2                 !! Index for P release from EC
  !
  ! NUMERICAL AND PHYSICS CONSTANTS
  !
  !

  !-
  ! 1. Mathematical and numerical constants
  !-
  REAL, PARAMETER :: pi = 3.141592653589793238   !! pi souce : http://mathworld.wolfram.com/Pi.html (unitless)
  REAL, PARAMETER :: euler = 2.71828182845904523 !! e source : http://mathworld.wolfram.com/e.html (unitless)
  REAL, PARAMETER :: zero = 0._r_std             !! Numerical constant set to 0 (unitless)
  REAL, PARAMETER :: undemi = 0.5_r_std          !! Numerical constant set to 1/2 (unitless)
  REAL, PARAMETER :: un = 1._r_std               !! Numerical constant set to 1 (unitless)
  REAL, PARAMETER :: moins_un = -1._r_std        !! Numerical constant set to -1 (unitless)
  REAL, PARAMETER :: deux = 2._r_std             !! Numerical constant set to 2 (unitless)
  REAL, PARAMETER :: trois = 3._r_std            !! Numerical constant set to 3 (unitless)
  REAL, PARAMETER :: quatre = 4._r_std           !! Numerical constant set to 4 (unitless)
  REAL, PARAMETER :: cinq = 5._r_std             !![DISPENSABLE] Numerical constant set to 5 (unitless)
  REAL, PARAMETER :: six = 6._r_std              !![DISPENSABLE] Numerical constant set to 6 (unitless)
  REAL, PARAMETER :: huit = 8._r_std             !! Numerical constant set to 8 (unitless)
  REAL, PARAMETER :: mille = 1000._r_std         !! Numerical constant set to 1000 (unitless)

  !-
  ! 2 . Physics
  !-
  REAL, PARAMETER :: R_Earth = 6378000.              !! radius of the Earth : Earth radius ~= Equatorial radius (m)
  REAL, PARAMETER :: mincos  = 0.0001                !! Minimum cosine value used for interpolation (unitless) 
  REAL, PARAMETER :: pb_std = 1013.                  !! standard pressure (hPa)
  REAL, PARAMETER :: ZeroCelsius = 273.15            !! 0 degre Celsius in degre Kelvin (K)
  REAL, PARAMETER :: tp_00 = 273.15                  !! 0 degre Celsius in degre Kelvin (K)
  REAL, PARAMETER :: chalsu0 = 2.8345E06             !! Latent heat of sublimation (J.kg^{-1})
  REAL, PARAMETER :: chalev0 = 2.5008E06             !! Latent heat of evaporation (J.kg^{-1}) 
  REAL, PARAMETER :: chalfu0 = chalsu0-chalev0       !! Latent heat of fusion (J.kg^{-1}) 
  REAL, PARAMETER :: c_stefan = 5.6697E-8            !! Stefan-Boltzman constant (W.m^{-2}.K^{-4})
  REAL, PARAMETER :: cp_air = 1004.675               !! Specific heat of dry air (J.kg^{-1}.K^{-1}) 
  REAL, PARAMETER :: cte_molr = 287.05               !! Specific constant of dry air (kg.mol^{-1}) 
  REAL, PARAMETER :: kappa = cte_molr/cp_air         !! Kappa : ratio between specific constant and specific heat 
                                                            !! of dry air (unitless)
  REAL, PARAMETER :: msmlr_air = 28.964E-03          !! Molecular weight of dry air (kg.mol^{-1})
  REAL, PARAMETER :: msmlr_h2o = 18.02E-03           !! Molecular weight of water vapor (kg.mol^{-1}) 
  REAL, PARAMETER :: cp_h2o = &                      !! Specific heat of water vapor (J.kg^{-1}.K^{-1}) 
       & cp_air*(quatre*msmlr_air)/( 3.5_r_std*msmlr_h2o) 
  REAL, PARAMETER :: cte_molr_h2o = cte_molr/quatre  !! Specific constant of water vapor (J.kg^{-1}.K^{-1}) 
  REAL, PARAMETER :: retv = msmlr_air/msmlr_h2o-un   !! Ratio between molecular weight of dry air and water 
                                                            !! vapor minus 1(unitless)  
  REAL, PARAMETER :: rvtmp2 = cp_h2o/cp_air-un       !! Ratio between specific heat of water vapor and dry air
                                                            !! minus 1 (unitless)
  REAL, PARAMETER :: cepdu2 = (0.1_r_std)**2         !! Squared wind shear (m^2.s^{-2}) 
  REAL, PARAMETER :: ct_karman = 0.41_r_std          !! Van Karmann Constant (unitless)
  REAL, PARAMETER :: cte_grav = 9.80665_r_std        !! Acceleration of the gravity (m.s^{-2})
  REAL, PARAMETER :: pa_par_hpa = 100._r_std         !! Transform pascal into hectopascal (unitless)
  REAL, PARAMETER :: RR = 8.314                      !! Ideal gas constant (J.mol^{-1}.K^{-1})
  REAL, PARAMETER :: Sct = 1370.                     !! Solar constant (W.m^{-2}) 


  INTEGER, SAVE :: testpft = 6
  !-
  ! 3. Climatic constants
  !-
  !! Constantes of the Louis scheme 
  REAL, SAVE :: cb = 5._r_std              !! Constant of the Louis scheme (unitless);
                                                  !! reference to Louis (1979)
!$OMP THREADPRIVATE(cb)
  REAL, SAVE :: cc = 5._r_std              !! Constant of the Louis scheme (unitless);
                                                  !! reference to Louis (1979)
!$OMP THREADPRIVATE(cc)
  REAL, SAVE :: cd = 5._r_std              !! Constant of the Louis scheme (unitless);
                                                  !! reference to Louis (1979)
!$OMP THREADPRIVATE(cd)
  REAL, SAVE :: rayt_cste = 125.           !! Constant in the computation of surface resistance (W.m^{-2})
!$OMP THREADPRIVATE(rayt_cste)
  REAL, SAVE :: defc_plus = 23.E-3         !! Constant in the computation of surface resistance (K.W^{-1})
!$OMP THREADPRIVATE(defc_plus)
  REAL, SAVE :: defc_mult = 1.5            !! Constant in the computation of surface resistance (K.W^{-1})
!$OMP THREADPRIVATE(defc_mult)

  !-
  ! 4. Soil thermodynamics constants
  !-
  ! Look at constantes_soil.f90


  !
  ! OPTIONAL PARTS OF THE MODEL
  !
  LOGICAL,PARAMETER :: diag_qsat = .TRUE.         !! One of the most frequent problems is a temperature out of range
                                                  !! we provide here a way to catch that in the calling procedure. 
                                                  !! (from Jan Polcher)(true/false) 
  LOGICAL, SAVE     :: almaoutput =.FALSE.        !! Selects the type of output for the model.(true/false)
                                                  !! Value is read from run.def in intersurf_history
!$OMP THREADPRIVATE(almaoutput)

  !
  ! DIVERSE
  !
  CHARACTER(LEN=100), SAVE :: stomate_forcing_name='NONE'  !! NV080800 Name of STOMATE forcing file (unitless)
                                                           ! Compatibility with Nicolas Viovy driver.
!$OMP THREADPRIVATE(stomate_forcing_name)
  CHARACTER(LEN=100), SAVE :: stomate_Cforcing_name='NONE' !! NV080800 Name of soil forcing file (unitless)
                                                           ! Compatibility with Nicolas Viovy driver.
!$OMP THREADPRIVATE(stomate_Cforcing_name)
  INTEGER, SAVE :: forcing_id                 !! Index of the forcing file (unitless)
!$OMP THREADPRIVATE(forcing_id)
  LOGICAL, SAVE :: allow_forcing_write=.TRUE.        !! Allow writing of stomate_forcing file. 
                                                     !! This variable will be set to false for teststomate. 



                         !------------------------!
                         !  SECHIBA PARAMETERS    !
                         !------------------------!
 

  !
  ! GLOBAL PARAMETERS   
  !
  REAL, SAVE :: min_wind = 0.1      !! The minimum wind (m.s^{-1})
!$OMP THREADPRIVATE(min_wind)
  REAL, SAVE :: snowcri = 1.5       !! Sets the amount above which only sublimation occures (kg.m^{-2})
!$OMP THREADPRIVATE(snowcri)


  !
  ! FLAGS ACTIVATING SUB-MODELS
  !
  LOGICAL, SAVE :: treat_expansion = .FALSE.   !! Do we treat PFT expansion across a grid point after introduction? (true/false)
!$OMP THREADPRIVATE(treat_expansion)
  LOGICAL, SAVE :: ok_herbivores = .FALSE.     !! flag to activate herbivores (true/false)
!$OMP THREADPRIVATE(ok_herbivores)
  LOGICAL, SAVE :: harvest_agri = .TRUE.       !! flag to harvest aboveground biomass from agricultural PFTs)(true/false)
!$OMP THREADPRIVATE(harvest_agri)
  LOGICAL, SAVE :: harvest_be = .TRUE.       !! flag to harvest aboveground biomass from bioenergy grass PFTs)(true/false)
!$OMP THREADPRIVATE(harvest_be)
  LOGICAL, SAVE :: lpj_gap_const_mort          !! constant moratlity (true/false). Default value depend on OK_DGVM.
!$OMP THREADPRIVATE(lpj_gap_const_mort)
  LOGICAL, SAVE :: disable_fire = .FALSE.      !! flag that disable fire (true/false)
!$OMP THREADPRIVATE(disable_fire)
  LOGICAL, SAVE :: spinup_analytic = .FALSE.   !! Flag to activate analytical resolution for spinup (true/false)
!$OMP THREADPRIVATE(spinup_analytic)
  LOGICAL, SAVE :: spinup_cnp = .FALSE.        !! Flag to add nutrient in spinup simulations if immobilization exceeds nutrient supply
!$OMP THREADPRIVATE(spinup_cnp)

  LOGICAL, SAVE :: ok_explicitsnow             !! Flag to activate explicit snow scheme instead of default snow scheme
!$OMP THREADPRIVATE(ok_explicitsnow)

  !
  ! CONFIGURATION VEGETATION
  !
  LOGICAL, SAVE :: agriculture = .TRUE.    !! allow agricultural PFTs (true/false)
!$OMP THREADPRIVATE(agriculture)
  LOGICAL, SAVE :: impveg = .FALSE.        !! Impose vegetation ? (true/false)
!$OMP THREADPRIVATE(impveg)
  LOGICAL, SAVE :: impsoilt = .FALSE.      !! Impose soil texture? (true/false)
!$OMP THREADPRIVATE(impsoilt)
  LOGICAL, SAVE :: impsoils = .FALSE.      !! Impose soil order ? (true/false)
!$OMP THREADPRIVATE(impsoils)
  LOGICAL, SAVE :: impose_Ninput = .FALSE. !! Impose nutrient input values ? (true/false)
!$OMP THREADPRIVATE(impose_Ninput)
!  LOGICAL, SAVE :: impose_Pinput = .FALSE. !! Impose P input values ? (true/false)
!!$OMP THREADPRIVATE(impose_Pinput)
  LOGICAL, SAVE :: read_pweat = .FALSE. !! Impose P weathering input from netcdf map ? (true/false)
!$OMP THREADPRIVATE(read_pweat)
  LOGICAL, SAVE :: impose_Nmap = .FALSE. !! Impose N/P input from netcdf map ? (true/false)
!$OMP THREADPRIVATE(impose_Nmap)
  LOGICAL, SAVE :: read_bnf = .FALSE. !! Impose N/P input from netcdf map ? (true/false)
!$OMP THREADPRIVATE(read_bnf)
  LOGICAL, SAVE :: do_now_stomate_lcchange = .FALSE.  !! Time to call lcchange in stomate_lpj
!$OMP THREADPRIVATE(do_now_stomate_lcchange)
  LOGICAL, SAVE :: do_now_stomate_woodharvest = .FALSE.  !! Time to call woodharvest in stomate_lpj
!$OMP THREADPRIVATE(do_now_stomate_woodharvest)
  LOGICAL, SAVE :: done_stomate_lcchange = .FALSE.    !! If true, call lcchange in stomate_lpj has just been done. 
!$OMP THREADPRIVATE(done_stomate_lcchange)
  LOGICAL, SAVE :: read_lai = .FALSE.      !! Flag to read a map of LAI if STOMATE is not activated (true/false)
!$OMP THREADPRIVATE(read_lai)
  LOGICAL, SAVE :: map_pft_format = .TRUE. !! Read a land use vegetation map on PFT format (true/false)
!$OMP THREADPRIVATE(map_pft_format)
  LOGICAL, SAVE :: veget_reinit = .TRUE.   !! To change LAND USE file in a run. (true/false)
!$OMP THREADPRIVATE(veget_reinit)
  LOGICAL, SAVE :: ninput_reinit = .TRUE.  !! To change N INPUT file in a run. (true/false)
!$OMP THREADPRIVATE(ninput_reinit)
  LOGICAL, SAVE :: allow_agri_fert = .FALSE.  !! To read N P fertilization file in a run. (true/false)
!$OMP THREADPRIVATE(allow_agri_fert)
  !
  ! PARAMETERS USED BY BOTH HYDROLOGY MODELS
  !
  REAL, SAVE :: max_snow_age = 50._r_std !! Maximum period of snow aging (days)
!$OMP THREADPRIVATE(max_snow_age)
  REAL, SAVE :: snow_trans = 0.2_r_std   !! Transformation time constant for snow (m), reduced from the value 0.3 (04/07/2016)
!$OMP THREADPRIVATE(snow_trans)
  REAL, SAVE :: sneige                   !! Lower limit of snow amount (kg.m^{-2})
!$OMP THREADPRIVATE(sneige)
  REAL, SAVE :: maxmass_snow = 3000.     !! The maximum mass of snow (kg.m^{-2})
!$OMP THREADPRIVATE(maxmass_snow)

  !! Heat capacity
  REAL, PARAMETER :: capa_ice = 2.228*1.E3       !! Heat capacity of ice (J/kg/K)
  REAL, SAVE      :: so_capa_ice                 !! Heat capacity of saturated frozen soil (J/K/m3)
!$OMP THREADPRIVATE(so_capa_ice)
  REAL, PARAMETER :: rho_water = 1000.           !! Density of water (kg/m3)
  REAL, PARAMETER :: rho_ice = 920.              !! Density of ice (kg/m3)

  !! Thermal conductivities
  REAL, PARAMETER :: cond_water = 0.6            !! Thermal conductivity of liquid water (W/m/K)
  REAL, PARAMETER :: cond_ice = 2.2              !! Thermal conductivity of ice (W/m/K)
  REAL, PARAMETER :: cond_solid = 2.32           !! Thermal conductivity of mineral soil particles (W/m/K)

  !! Time constant of long-term soil humidity (s) 
  REAL, PARAMETER :: lhf = 0.3336*1.E6           !! Latent heat of fusion (J/kg)

  INTEGER, PARAMETER :: nsnow=3                  !! Number of levels in the snow for explicit snow scheme   
  REAL, PARAMETER    :: XMD    = 28.9644E-3 
  REAL, PARAMETER    :: XBOLTZ      = 1.380658E-23 
  REAL, PARAMETER    :: XAVOGADRO   = 6.0221367E+23 
  REAL, PARAMETER    :: XRD    = XAVOGADRO * XBOLTZ / XMD 
  REAL, PARAMETER    :: XCPD   = 7.* XRD /2. 
  REAL, PARAMETER    :: phigeoth = 0.057 ! 0. DKtest 
  REAL, PARAMETER    :: thick_min_snow = .01 

  !! The maximum snow density and water holding characterisicts 
  REAL, SAVE         :: xrhosmax = 750.  ! (kg m-3) 
  REAL, SAVE         :: xwsnowholdmax1   = 0.03  ! (-) 
  REAL, SAVE         :: xwsnowholdmax2   = 0.10  ! (-) 
  REAL, SAVE         :: xsnowrhohold     = 200.0 ! (kg/m3) 
  REAL, SAVE         :: xrhosmin = 50. 
  REAL, PARAMETER    :: xci = 2.106e+3 
  REAL, PARAMETER    :: xrv = 6.0221367e+23 * 1.380658e-23 /18.0153e-3 

  !! ISBA-ES Critical snow depth at which snow grid thicknesses constant 
  REAL, PARAMETER    :: xsnowcritd = 0.03  ! (m) 

  !! The threshold of snow depth used for preventing numerical problem in thermal calculations
  REAL, PARAMETER    :: snowcritd_thermal = 0.01  ! (m)  
  
  !! ISBA-ES CROCUS (Pahaut 1976): snowfall density coefficients: 
  REAL, PARAMETER       :: snowfall_a_sn = 109.0  !! (kg/m3) 
  REAL, PARAMETER       :: snowfall_b_sn =   6.0  !! (kg/m3/K) 
  REAL, PARAMETER       :: snowfall_c_sn =  26.0  !! [kg/(m7/2 s1/2)] 

  REAL, PARAMETER       :: dgrain_new_max=  2.0e-4!! (m) : Maximum grain size of new snowfall 
  
  !! Used in explicitsnow to prevent numerical problems as snow becomes vanishingly thin. 
  REAL, PARAMETER                :: psnowdzmin = .0001   ! m 
  REAL, PARAMETER                :: xsnowdmin = .000001  ! m 

  REAL, PARAMETER                :: ph2o = 1000.         !! Water density [kg/m3] 
  
  ! ISBA-ES Thermal conductivity coefficients from Anderson (1976): 
  ! see Boone, Meteo-France/CNRM Note de Centre No. 70 (2002) 
  REAL, SAVE                     :: ZSNOWTHRMCOND1 = 0.02    ! [W/m/K] 
  REAL, SAVE                     :: ZSNOWTHRMCOND2 = 2.5E-6  ! [W m5/(kg2 K)] 
  
  ! ISBA-ES Thermal conductivity: Implicit vapor diffn effects 
  ! (sig only for new snow OR high altitudes) 
  ! from Sun et al. (1999): based on data from Jordan (1991) 
  ! see Boone, Meteo-France/CNRM Note de Centre No. 70 (2002) 
  ! 
  REAL, SAVE                       :: ZSNOWTHRMCOND_AVAP  = -0.06023 ! (W/m/K) 
  REAL, SAVE                       :: ZSNOWTHRMCOND_BVAP  = -2.5425  ! (W/m) 
  REAL, SAVE                       :: ZSNOWTHRMCOND_CVAP  = -289.99  ! (K) 
  
  REAL,SAVE :: xansmax = 0.85      !! Maxmimum snow albedo
  REAL,SAVE :: xansmin = 0.50      !! Miniumum snow albedo
  REAL,SAVE :: xans_todry = 0.008  !! Albedo decay rate for dry snow
  REAL,SAVE :: xans_t = 0.240      !! Albedo decay rate for wet snow

  ! ISBA-ES Thermal conductivity coefficients from Anderson (1976):
  ! see Boone, Meteo-France/CNRM Note de Centre No. 70 (2002)
  REAL, PARAMETER                  :: XP00 = 1.E5

  ! ISBA-ES Thermal conductivity: Implicit vapor diffn effects
  ! (sig only for new snow OR high altitudes)
  ! from Sun et al. (1999): based on data from Jordan (1991)
  ! see Boone, Meteo-France/CNRM Note de Centre No. 70 (2002)
  !
  REAL, SAVE          :: ZSNOWCMPCT_RHOD  = 150.0        !! (kg/m3)
  REAL, SAVE          :: ZSNOWCMPCT_ACM   = 2.8e-6       !! (1/s)
  REAL, SAVE          :: ZSNOWCMPCT_BCM   = 0.04         !! (1/K) 
  REAL, SAVE          :: ZSNOWCMPCT_CCM   = 460.         !! (m3/kg)
  REAL, SAVE          :: ZSNOWCMPCT_V0    = 3.7e7        !! (Pa/s) 
  REAL, SAVE          :: ZSNOWCMPCT_VT    = 0.081        !! (1/K)
  REAL, SAVE          :: ZSNOWCMPCT_VR    = 0.018        !! (m3/kg)

  !
  ! BVOC : Biogenic activity  for each age class
  !
  REAL, SAVE, DIMENSION(nleafages) :: iso_activity = (/0.5, 1.5, 1.5, 0.5/)     !! Biogenic activity for each 
                                                                                       !! age class : isoprene (unitless)
!$OMP THREADPRIVATE(iso_activity)
  REAL, SAVE, DIMENSION(nleafages) :: methanol_activity = (/1., 1., 0.5, 0.5/)  !! Biogenic activity for each
                                                                                       !! age class : methanol (unnitless)
!$OMP THREADPRIVATE(methanol_activity)

  !
  ! condveg.f90
  !

  ! 1. Scalar

  ! 1.1 Flags used inside the module

  LOGICAL, SAVE :: alb_bare_model = .FALSE. !! Switch for choosing values of bare soil 
                                            !! albedo (see header of subroutine)
                                            !! (true/false)
!$OMP THREADPRIVATE(alb_bare_model)
  LOGICAL, SAVE :: alb_bg_modis = .FALSE.   !! Switch for choosing values of bare soil 
                                            !! albedo read from file
                                            !! (true/false)
!$OMP THREADPRIVATE(alb_bg_modis)
  LOGICAL, SAVE :: impaze = .FALSE.         !! Switch for choosing surface parameters
                                            !! (see header of subroutine).  
                                            !! (true/false)
!$OMP THREADPRIVATE(impaze)
  LOGICAL, SAVE :: rough_dyn = .FALSE.      !! Chooses between two methods to calculate the 
                                            !! the roughness height : static or dynamic (varying with LAI)
                                            !! (true/false)
!$OMP THREADPRIVATE(rough_dyn)

  LOGICAL, SAVE :: new_watstress = .FALSE.
!$OMP THREADPRIVATE(new_watstress)

  REAL, SAVE :: alpha_watstress = 1.
!$OMP THREADPRIVATE(alpha_watstress)

  ! 1.2 Others 


  REAL, SAVE :: height_displacement = 0.66        !! Factor to calculate the zero-plane displacement
                                                         !! height from vegetation height (m)
!$OMP THREADPRIVATE(height_displacement)
  REAL, SAVE :: z0_bare = 0.01                    !! bare soil roughness length (m)
!$OMP THREADPRIVATE(z0_bare)
  REAL, SAVE :: z0_ice = 0.001                    !! ice roughness length (m)
!$OMP THREADPRIVATE(z0_ice)
  REAL, SAVE :: tcst_snowa = 3.0                 !! Time constant of the albedo decay of snow (days), increased from the value 5.0 (04/07/2016)
!$OMP THREADPRIVATE(tcst_snowa)
  REAL, SAVE :: snowcri_alb = 10.                 !! Critical value for computation of snow albedo (cm)
!$OMP THREADPRIVATE(snowcri_alb)
  REAL, SAVE :: fixed_snow_albedo = undef_sechiba !! To choose a fixed snow albedo value (unitless)
!$OMP THREADPRIVATE(fixed_snow_albedo)
  REAL, SAVE :: z0_scal = 0.15                    !! Surface roughness height imposed (m)
!$OMP THREADPRIVATE(z0_scal)
  REAL, SAVE :: roughheight_scal = zero           !! Effective roughness Height depending on zero-plane 
                                                         !! displacement height (m) (imposed)
!$OMP THREADPRIVATE(roughheight_scal)
  REAL, SAVE :: emis_scal = 1.0                   !! Surface emissivity imposed (unitless)
!$OMP THREADPRIVATE(emis_scal)

  REAL, SAVE :: c1 = 0.32                         !! Constant used in the formulation of the ratio of 
!$OMP THREADPRIVATE(c1)                                  !! friction velocity to the wind speed at the canopy top
                                                         !! see Ershadi et al. (2015) for more info
  REAL, SAVE :: c2 = 0.264                        !! Constant used in the formulation of the ratio of 
!$OMP THREADPRIVATE(c2)                                  !! friction velocity to the wind speed at the canopy top
                                                         !! see Ershadi et al. (2015) for more info
  REAL, SAVE :: c3 = 15.1                         !! Constant used in the formulation of the ratio of 
!$OMP THREADPRIVATE(c3)                                  !! friction velocity to the wind speed at the canopy top
                                                         !! see Ershadi et al. (2015) for more info
  REAL, SAVE :: Cdrag_foliage = 0.2               !! Drag coefficient of the foliage
!$OMP THREADPRIVATE(Cdrag_foliage)                       !! See Ershadi et al. (2015) and Su et. al (2001) for more info
  REAL, SAVE :: Ct = 0.01                         !! Heat transfer coefficient of the leaf
!$OMP THREADPRIVATE(Ct)                                  !! See Ershadi et al. (2015) and Su et. al (2001) for more info
  REAL, SAVE :: Prandtl = 0.71                    !! Prandtl number used in the calculation of Ct_star
!$OMP THREADPRIVATE(Prandtl)                             !! See Su et. al (2001) for more info



  ! 2. Arrays

  REAL, SAVE, DIMENSION(2) :: alb_deadleaf = (/ .12, .35/)    !! albedo of dead leaves, VIS+NIR (unitless)
!$OMP THREADPRIVATE(alb_deadleaf)
  REAL, SAVE, DIMENSION(2) :: alb_ice = (/ .60, .20/)         !! albedo of ice, VIS+NIR (unitless)
!$OMP THREADPRIVATE(alb_ice)
  REAL, SAVE, DIMENSION(2) :: albedo_scal = (/ 0.25, 0.25 /)  !! Albedo values for visible and near-infrared 
                                                                     !! used imposed (unitless) 
!$OMP THREADPRIVATE(albedo_scal)
  REAL , SAVE, DIMENSION(classnb) :: vis_dry = (/0.24,&
       &0.22, 0.20, 0.18, 0.16, 0.14, 0.12, 0.10, 0.27/)  !! Soil albedo values to soil colour classification:
                                                          !! dry soil albedo values in visible range
!$OMP THREADPRIVATE(vis_dry)
  REAL, SAVE, DIMENSION(classnb) :: nir_dry = (/0.48,&
       &0.44, 0.40, 0.36, 0.32, 0.28, 0.24, 0.20, 0.55/)  !! Soil albedo values to soil colour classification:
                                                          !! dry soil albedo values in near-infrared range 
!$OMP THREADPRIVATE(nir_dry)
  REAL, SAVE, DIMENSION(classnb) :: vis_wet = (/0.12,&
       &0.11, 0.10, 0.09, 0.08, 0.07, 0.06, 0.05, 0.15/)  !! Soil albedo values to soil colour classification:
                                                          !! wet soil albedo values in visible range 
!$OMP THREADPRIVATE(vis_wet)
  REAL, SAVE, DIMENSION(classnb) :: nir_wet = (/0.24,&
       &0.22, 0.20, 0.18, 0.16, 0.14, 0.12, 0.10, 0.31/)  !! Soil albedo values to soil colour classification:
                                                          !! wet soil albedo values in near-infrared range
!$OMP THREADPRIVATE(nir_wet)
  REAL, SAVE, DIMENSION(classnb) :: albsoil_vis = (/ &
       &0.18, 0.16, 0.16, 0.15, 0.12, 0.105, 0.09, 0.075, 0.25/)   !! Soil albedo values to soil colour classification:
                                                                   !! Averaged of wet and dry soil albedo values
                                                                   !! in visible and near-infrared range
!$OMP THREADPRIVATE(albsoil_vis) 
  REAL, SAVE, DIMENSION(classnb) :: albsoil_nir = (/ &
       &0.36, 0.34, 0.34, 0.33, 0.30, 0.25, 0.20, 0.15, 0.45/)  !! Soil albedo values to soil colour classification:
                                                                !! Averaged of wet and dry soil albedo values
                                                                !! in visible and near-infrared range
!$OMP THREADPRIVATE(albsoil_nir)

  !
  ! diffuco.f90
  !

  ! 0. Constants

  REAL, PARAMETER :: Tetens_1 = 0.622         !! Ratio between molecular weight of water vapor and molecular weight  
                                                     !! of dry air (unitless)
  REAL, PARAMETER :: Tetens_2 = 0.378         !!
  REAL, PARAMETER :: ratio_H2O_to_CO2 = 1.6   !! Ratio of water vapor diffusivity to the CO2 diffusivity (unitless)
  REAL, PARAMETER :: mol_to_m_1 = 0.0244      !!
  REAL, PARAMETER :: RG_to_PAR = 0.5          !!
  REAL, PARAMETER :: W_to_mol = 4.6          !! W_to_mmol * RG_to_PAR = 2.3

  ! 1. Scalar

  INTEGER, SAVE :: nlai = 20             !! Number of LAI levels (unitless)
!$OMP THREADPRIVATE(nlai)
  LOGICAL, SAVE :: ldq_cdrag_from_gcm = .FALSE. !! Set to .TRUE. if you want q_cdrag coming from GCM
!$OMP THREADPRIVATE(ldq_cdrag_from_gcm)
  REAL, SAVE :: laimax = 12.             !! Maximal LAI used for splitting LAI into N layers (m^2.m^{-2})
!$OMP THREADPRIVATE(laimax)
  LOGICAL, SAVE :: downregulation_co2 = .FALSE.            !! Set to .TRUE. if you want CO2 downregulation.
!$OMP THREADPRIVATE(downregulation_co2)
  REAL, SAVE :: downregulation_co2_baselevel = 280. !! CO2 base level (ppm)
!$OMP THREADPRIVATE(downregulation_co2_baselevel)

  REAL, SAVE :: gb_ref = 1./25.                     !! Leaf bulk boundary layer resistance (s m-1)

  ! 3. Coefficients of equations

  REAL, SAVE :: lai_level_depth = 0.15  !!
!$OMP THREADPRIVATE(lai_level_depth)
!
  REAL, SAVE, DIMENSION(6) :: dew_veg_poly_coeff = &            !! coefficients of the 5 degree polynomomial used
  & (/ 0.887773, 0.205673, 0.110112, 0.014843,  0.000824,  0.000017 /) !! in the equation of coeff_dew_veg
!$OMP THREADPRIVATE(dew_veg_poly_coeff)
!
  REAL, SAVE               :: Oi=210000.    !! Intercellular oxygen partial pressure (ubar)
!$OMP THREADPRIVATE(Oi)
  !
  ! slowproc.f90 
  !

  ! 1. Scalar

  INTEGER, SAVE :: veget_year_orig = 0        !!  first year for landuse (number)
!$OMP THREADPRIVATE(veget_year_orig)
  INTEGER, SAVE :: ninput_year_orig = 0       !!  first year for N inputs (number)
!$OMP THREADPRIVATE(ninput_year_orig)
  REAL, SAVE :: clayfraction_default = 0.2    !! Default value for clay fraction (0-1, unitless)
!$OMP THREADPRIVATE(clayfraction_default)
  REAL, SAVE :: siltfraction_default = 0.5    !! Default value for silt fraction (0-1, unitless)
!$OMP THREADPRIVATE(siltfraction_default)
  REAL, SAVE :: bulk_default = 1000           !! Default value for bulk density of soil (kg/m3) DSG: value and units are consistent
!$OMP THREADPRIVATE(bulk_default)
  REAL, SAVE :: ph_default = 5.5              !! Default value for pH of soil (-)
!$OMP THREADPRIVATE(ph_default)
  REAL, SAVE , DIMENSION(3) :: p_input_default = (/ 0.0 , 0.0 ,0.0 /)         !! DEfault value for P input [g/m/yr]
!$OMP THREADPRIVATE(p_input_default)
  REAL, SAVE :: min_vegfrac = 0.001           !! Minimal fraction of mesh a vegetation type can occupy (0-1, unitless)
!$OMP THREADPRIVATE(min_vegfrac)
  REAL, SAVE :: frac_nobio_fixed_test_1 = 0.0 !! Value for frac_nobio for tests in 0-dim simulations (0-1, unitless)
!$OMP THREADPRIVATE(frac_nobio_fixed_test_1)
  
  REAL, SAVE :: stempdiag_bid = 280.          !! only needed for an initial LAI if there is no restart file
!$OMP THREADPRIVATE(stempdiag_bid)


                           !-----------------------------!
                           !  STOMATE AND LPJ PARAMETERS !
                           !-----------------------------!


  !
  ! lpj_constraints.f90
  !
  
  ! 1. Scalar

  REAL, SAVE  :: too_long = 5.      !! longest sustainable time without 
                                           !! regeneration (vernalization) (years)
!$OMP THREADPRIVATE(too_long)


  !
  ! lpj_establish.f90
  !

  ! 1. Scalar

  REAL, SAVE :: estab_max_tree = 0.12   !! Maximum tree establishment rate (ind/m2/dt_stomate)
!$OMP THREADPRIVATE(estab_max_tree)
  REAL, SAVE :: estab_max_grass = 0.12  !! Maximum grass establishment rate (ind/m2/dt_stomate)
!$OMP THREADPRIVATE(estab_max_grass)
  
  ! 3. Coefficients of equations

  REAL, SAVE :: establish_scal_fact = 5.  !!
!$OMP THREADPRIVATE(establish_scal_fact)
  REAL, SAVE :: max_tree_coverage = 0.98  !! (0-1, unitless)
!$OMP THREADPRIVATE(max_tree_coverage)
  REAL, SAVE :: ind_0_estab = 0.2         !! = ind_0 * 10.
!$OMP THREADPRIVATE(ind_0_estab)


  !
  ! lpj_fire.f90
  !

  ! 1. Scalar

  REAL, SAVE :: tau_fire = 30.           !! Time scale for memory of the fire index (days).
!$OMP THREADPRIVATE(tau_fire)
  REAL, SAVE :: litter_crit = 200.       !! Critical litter quantity for fire
                                                !! below which iginitions extinguish 
                                                !! @tex $(gC m^{-2})$ @endtex
!$OMP THREADPRIVATE(litter_crit)
  REAL, SAVE :: fire_resist_lignin = 0.5 !!
!$OMP THREADPRIVATE(fire_resist_lignin)
  ! 2. Arrays

  REAL, SAVE, DIMENSION(nparts) :: co2frac = &    !! The fraction of the different biomass 
       & (/ .95, .95, 0., 0.3, 0., 0., .95, .95, .95 /)       !! compartments emitted to the atmosphere 
!$OMP THREADPRIVATE(co2frac)                                                         !! when burned (unitless, 0-1)  

  ! 3. Coefficients of equations

  REAL, SAVE, DIMENSION(3) :: bcfrac_coeff = (/ .3,  1.3,  88.2 /)         !! (unitless)
!$OMP THREADPRIVATE(bcfrac_coeff)
  REAL, SAVE, DIMENSION(4) :: firefrac_coeff = (/ 0.45, 0.8, 0.6, 0.13 /)  !! (unitless)
!$OMP THREADPRIVATE(firefrac_coeff)

  !
  ! lpj_gap.f90
  !

  ! 1. Scalar

  REAL, SAVE :: ref_greff = 0.035         !! Asymptotic maximum mortality rate
                                                 !! @tex $(year^{-1})$ @endtex
!$OMP THREADPRIVATE(ref_greff)

  !               
  ! lpj_light.f90 
  !              

  ! 1. Scalar
  
  LOGICAL, SAVE :: annual_increase = .TRUE. !! for diagnosis of fpc increase, compare today's fpc to last year's maximum (T) or
                                            !! to fpc of last time step (F)? (true/false)
!$OMP THREADPRIVATE(annual_increase)
  REAL, SAVE :: min_cover = 0.05     !! For trees, minimum fraction of crown area occupied
                                            !! (due to its branches etc.) (0-1, unitless)
                                            !! This means that only a small fraction of its crown area
                                            !! can be invaded by other trees.
!$OMP THREADPRIVATE(min_cover)
  !
  ! lpj_pftinout.f90 
  !

  ! 1. Scalar

  REAL, SAVE :: min_avail = 0.01         !! minimum availability
!$OMP THREADPRIVATE(min_avail)
  REAL, SAVE :: ind_0 = 0.02             !! initial density of individuals
!$OMP THREADPRIVATE(ind_0)
  ! 3. Coefficients of equations
  
  REAL, SAVE :: RIP_time_min = 1.25      !! test whether the PFT has been eliminated lately (years)
!$OMP THREADPRIVATE(RIP_time_min)
  REAL, SAVE :: npp_longterm_init = 10.  !! Initialisation value for npp_longterm (gC.m^{-2}.year^{-1})
!$OMP THREADPRIVATE(npp_longterm_init)
  REAL, SAVE :: everywhere_init = 0.05   !!
!$OMP THREADPRIVATE(everywhere_init)




  !
  ! stomate_data.f90 
  !

  ! 1. Scalar 

  ! 1.2 climatic parameters 

  REAL, SAVE :: precip_crit = 100.        !! minimum precip, in (mm/year)
!$OMP THREADPRIVATE(precip_crit)
  REAL, SAVE :: gdd_crit_estab = 150.     !! minimum gdd for establishment of saplings
!$OMP THREADPRIVATE(gdd_crit_estab)
  REAL, SAVE :: fpc_crit = 0.95           !! critical fpc, needed for light competition and establishment (0-1, unitless)
!$OMP THREADPRIVATE(fpc_crit)

  ! 1.3 sapling characteristics

  REAL, SAVE :: alpha_grass = 0.5         !! alpha coefficient for grasses (unitless)
!$OMP THREADPRIVATE(alpha_grass)
  REAL, SAVE :: alpha_tree = 1.           !! alpha coefficient for trees (unitless)
!$OMP THREADPRIVATE(alpha_tree)
  REAL, SAVE :: struct_to_leaves = 0.05  !! Fraction of structural carbon in grass and crops as a share of the leaf
                                                !! carbon pool. Only used for grasses and crops (thus NOT for trees)
                                                !! (unitless)
!$OMP THREADPRIVATE(struct_to_leaves)

  REAL, SAVE :: labile_to_total = 0.01   !! Fraction of the labile pool in trees, grasses and crops as a share of the
                                                !! total carbon pool (accounting for the N-content of the different tissues).
                                                !! (unitless)
!$OMP THREADPRIVATE(labile_to_total)



  ! 1.4  time scales for phenology and other processes (in days)
  REAL, SAVE :: tau_hum_month = 20.        !! (days)       
!$OMP THREADPRIVATE(tau_hum_month)
  REAL, SAVE :: tau_hum_week = 7.          !! (days)  
!$OMP THREADPRIVATE(tau_hum_week)
  REAL, SAVE :: tau_t2m_month = 20.        !! (days)      
!$OMP THREADPRIVATE(tau_t2m_month)
  REAL, SAVE :: tau_t2m_week = 7.          !! (days)  
!$OMP THREADPRIVATE(tau_t2m_week)
  REAL, SAVE :: tau_tsoil_month = 20.      !! (days)     
!$OMP THREADPRIVATE(tau_tsoil_month)
  REAL, SAVE :: tau_soilhum_month = 20.    !! (days)     
!$OMP THREADPRIVATE(tau_soilhum_month)
  REAL, SAVE :: tau_gpp_week = 7.          !! (days)  
!$OMP THREADPRIVATE(tau_gpp_week)
  REAL, SAVE :: tau_gdd = 40.              !! (days)  
!$OMP THREADPRIVATE(tau_gdd)
  REAL, SAVE :: tau_ngd = 50.              !! (days)  
!$OMP THREADPRIVATE(tau_ngd)
  REAL, SAVE :: coeff_tau_longterm = 3.    !! (unitless)
!$OMP THREADPRIVATE(coeff_tau_longterm)
  REAL, SAVE :: tau_longterm_max           !! (days)  
!$OMP THREADPRIVATE(tau_longterm_max)

  ! 3. Coefficients of equations

  REAL, SAVE :: bm_sapl_carbres = 5.             !!
!$OMP THREADPRIVATE(bm_sapl_carbres)
  REAL, SAVE :: bm_sapl_labile = 5.             !!
!$OMP THREADPRIVATE(bm_sapl_carbres)
  REAL, SAVE :: bm_sapl_sapabove = 0.5           !!
!$OMP THREADPRIVATE(bm_sapl_sapabove)
  REAL, SAVE :: bm_sapl_heartabove = 2.          !!
!$OMP THREADPRIVATE(bm_sapl_heartabove)
  REAL, SAVE :: bm_sapl_heartbelow = 2.          !!
!$OMP THREADPRIVATE(bm_sapl_heartbelow)
  REAL, SAVE :: init_sapl_mass_leaf_nat = 0.1    !!
!$OMP THREADPRIVATE(init_sapl_mass_leaf_nat)
  REAL, SAVE :: init_sapl_mass_leaf_agri = 1.    !!
!$OMP THREADPRIVATE(init_sapl_mass_leaf_agri)
  REAL, SAVE :: init_sapl_mass_carbres = 5.      !!
!$OMP THREADPRIVATE(init_sapl_mass_carbres)
  REAL, SAVE :: init_sapl_mass_labile = 5.      !!
!$OMP THREADPRIVATE(init_sapl_mass_carbres)
  REAL, SAVE :: init_sapl_mass_root = 0.1        !!
!$OMP THREADPRIVATE(init_sapl_mass_root)
  REAL, SAVE :: init_sapl_mass_fruit = 0.3       !!  
!$OMP THREADPRIVATE(init_sapl_mass_fruit)
  REAL, SAVE :: cn_sapl_init = 0.5               !!
!$OMP THREADPRIVATE(cn_sapl_init)
  REAL, SAVE :: migrate_tree = 10.*1.E3          !!
!$OMP THREADPRIVATE(migrate_tree)
  REAL, SAVE :: migrate_grass = 10.*1.E3         !!
!$OMP THREADPRIVATE(migrate_grass)
  REAL, SAVE :: lai_initmin_tree = 0.3           !!
!$OMP THREADPRIVATE(lai_initmin_tree)
  REAL, SAVE :: lai_initmin_grass = 0.1          !!
!$OMP THREADPRIVATE(lai_initmin_grass)
  REAL, SAVE, DIMENSION(2) :: dia_coeff = (/ 4., 0.5 /)            !!
!$OMP THREADPRIVATE(dia_coeff)
  REAL, SAVE, DIMENSION(2) :: maxdia_coeff =(/ 100., 0.01/)        !!
!$OMP THREADPRIVATE(maxdia_coeff)
  REAL, SAVE, DIMENSION(4) :: bm_sapl_leaf = (/ 4., 4., 0.8, 5./)  !!
!$OMP THREADPRIVATE(bm_sapl_leaf)


  !
  ! stomate_litter.f90 
  !

  ! 0. Constants

  REAL, PARAMETER :: Q10 = 10.               !!

  ! 1. Scalar

  REAL, SAVE :: z_decomp = 0.2               !!  Maximum depth for soil decomposer's activity (m)
!$OMP THREADPRIVATE(z_decomp)

  ! 2. Arrays

  REAL, SAVE :: frac_soil_struct_aa = 0.55   !! corresponding to frac_soil(istructural,iactive,iabove) 
!$OMP THREADPRIVATE(frac_soil_struct_aa)
  REAL, SAVE :: frac_soil_struct_sua = 0.4    !! corresponding to frac_soil(istructural,isurface,iabove) 
!$OMP THREADPRIVATE(frac_soil_struct_sua)
  REAL, SAVE :: frac_soil_struct_ab = 0.45   !! corresponding to frac_soil(istructural,iactive,ibelow)
!$OMP THREADPRIVATE(frac_soil_struct_ab)
  REAL, SAVE :: frac_soil_struct_sa = 0.7    !! corresponding to frac_soil(istructural,islow,iabove)
!$OMP THREADPRIVATE(frac_soil_struct_sa)
  REAL, SAVE :: frac_soil_struct_sb = 0.7    !! corresponding to frac_soil(istructural,islow,ibelow)
!$OMP THREADPRIVATE(frac_soil_struct_sb)
  REAL, SAVE :: frac_soil_metab_aa = 0.45    !! corresponding to frac_soil(imetabolic,iactive,iabove)
!$OMP THREADPRIVATE(frac_soil_metab_aa)
  REAL, SAVE :: frac_soil_metab_sua = 0.4    !! corresponding to frac_soil(imetabolic,iactive,iabove)
!$OMP THREADPRIVATE(frac_soil_metab_sua)
  REAL, SAVE :: frac_soil_metab_ab = 0.45    !! corresponding to frac_soil(imetabolic,iactive,ibelow)
!$OMP THREADPRIVATE(frac_soil_metab_ab)
  REAL, SAVE, DIMENSION(nparts) :: CN_fix = &!! C/N ratio of each plant pool (0-100, unitless)
       & (/ 40., 40., 40., 40., 40., 40., 40., 40., 40./) 
!$OMP THREADPRIVATE(CN_fix)

  ! 3. Coefficients of equations

  REAL, SAVE :: metabolic_ref_frac = 0.85    !! used by litter and soilcarbon (0-1, unitless)
!$OMP THREADPRIVATE(metabolic_ref_frac)
  REAL, SAVE :: metabolic_LN_ratio = 0.018   !! (0-1, unitless)   
!$OMP THREADPRIVATE(metabolic_LN_ratio)
  REAL, SAVE :: tau_metabolic = 0.066        !!
!$OMP THREADPRIVATE(tau_metabolic)
  REAL, SAVE :: tau_struct = 0.245           !!
!$OMP THREADPRIVATE(tau_struct)
  ! Turnover rate (yr-1) - From Parton et al., 1993
  REAL, SAVE :: turn_metabolic = 15           !!
!$OMP THREADPRIVATE(turn_metabolic)
  REAL, SAVE :: turn_struct = 4                !!
!$OMP THREADPRIVATE(turn_struct)
  REAL, SAVE :: turn_woody = 1.33              !! from DOFOCO
!$OMP THREADPRIVATE(turn_woody)
  REAL, SAVE :: soil_Q10 = 0.69              !!= ln 2
!$OMP THREADPRIVATE(soil_Q10)
  REAL, SAVE :: tsoil_ref = 30.              !!
!$OMP THREADPRIVATE(tsoil_ref)
  REAL, SAVE :: litter_struct_coef = 3.      !! 
!$OMP THREADPRIVATE(litter_struct_coef)
  REAL, SAVE, DIMENSION(3) :: moist_coeff = (/ 1.1,  2.4,  0.29 /) !!
!$OMP THREADPRIVATE(moist_coeff)
  REAL, SAVE :: moistcont_min = 0.25  !! minimum soil wetness to limit the heterotrophic respiration
!$OMP THREADPRIVATE(moistcont_min)


  !
  ! stomate_lpj.f90
  !

  ! 1. Scalar

  REAL, SAVE :: frac_turnover_daily = 0.55  !! (0-1, unitless)
!$OMP THREADPRIVATE(frac_turnover_daily)


  !
  ! stomate_npp.f90 
  !

  ! 1. Scalar

  REAL, SAVE :: tax_max = 0.8 !! Maximum fraction of allocatable biomass used 
                                     !! for maintenance respiration (0-1, unitless)
!$OMP THREADPRIVATE(tax_max)


  !
  ! stomate_phenology.f90
  !

  ! 1. Scalar

  LOGICAL, SAVE :: always_init = .FALSE.           !! take carbon from atmosphere if carbohydrate reserve too small? (true/false)
!$OMP THREADPRIVATE(always_init)
  REAL, SAVE :: min_growthinit_time = 300.  !! minimum time since last beginning of a growing season (days)
!$OMP THREADPRIVATE(min_growthinit_time)
  REAL, SAVE :: moiavail_always_tree = 0.8  !! moisture monthly availability above which moisture tendency doesn't matter
                                                   !!  - for trees (0-1, unitless)
!$OMP THREADPRIVATE(moiavail_always_tree)
  REAL, SAVE :: moiavail_always_grass = 0.4 !! moisture monthly availability above which moisture tendency doesn't matter
                                                   !! - for grass (0-1, unitless)
!$OMP THREADPRIVATE(moiavail_always_grass)
  REAL, SAVE :: t_always                    !! monthly temp. above which temp. tendency doesn't matter
!$OMP THREADPRIVATE(t_always)
  REAL, SAVE :: t_always_add = 10.          !! monthly temp. above which temp. tendency doesn't matter (C)
!$OMP THREADPRIVATE(t_always_add)

  ! 3. Coefficients of equations
  
  REAL, SAVE :: gddncd_ref = 603.           !!
!$OMP THREADPRIVATE(gddncd_ref)
  REAL, SAVE :: gddncd_curve = 0.0091       !!
!$OMP THREADPRIVATE(gddncd_curve)
  REAL, SAVE :: gddncd_offset = 64.         !!
!$OMP THREADPRIVATE(gddncd_offset)


  !
  ! stomate_prescribe.f90
  !

  ! 3. Coefficients of equations

  REAL, SAVE :: bm_sapl_rescale = 40.       !!
!$OMP THREADPRIVATE(bm_sapl_rescale)


  !
  ! stomate_resp.f90
  !

  ! 3. Coefficients of equations

  REAL, SAVE :: maint_resp_min_vmax = 0.3   !!
!$OMP THREADPRIVATE(maint_resp_min_vmax)
  REAL, SAVE :: maint_resp_coeff = 1.4      !!
!$OMP THREADPRIVATE(maint_resp_coeff)


  !
  ! stomate_som_dynamics.f90 (in stomate_soilcarbon.f90)   
  !

  ! 2. Arrays 

  ! 2.1 frac_carb_coefficients

  REAL, SAVE :: frac_carb_ap = 0.004  !! from active pool: depends on clay content  (0-1, unitless)
                                             !! corresponding to frac_carb(:,iactive,ipassive)
!$OMP THREADPRIVATE(frac_carb_ap)
  REAL, SAVE :: frac_carb_sa = 0.42   !! from slow pool (0-1, unitless)
                                             !! corresponding to frac_carb(:,islow,iactive)
!$OMP THREADPRIVATE(frac_carb_sa)
  REAL, SAVE :: frac_carb_sp = 0.03   !! from slow pool (0-1, unitless) 
                                             !! corresponding to frac_carb(:,islow,ipassive)
!$OMP THREADPRIVATE(frac_carb_sp)
  REAL, SAVE :: frac_carb_pa = 0.45   !! from passive pool (0-1, unitless)
                                             !! corresponding to frac_carb(:,ipassive,iactive)
!$OMP THREADPRIVATE(frac_carb_pa)
  REAL, SAVE :: frac_carb_ps = 0.0    !! from passive pool (0-1, unitless)
                                             !! corresponding to frac_carb(:,ipassive,islow)
!$OMP THREADPRIVATE(frac_carb_ps)

  ! 2.1 Fixed fraction from one pool to another (or to CO2 emission)

  REAL, SAVE :: active_to_pass_ref_frac = 0.003  !! from active pool: depends on clay content  (0-1, unitless)
                                                        !! corresponding to frac_carb(:,iactive,ipassive)
  REAL, SAVE :: surf_to_slow_ref_frac = 0.4      !! from surface pool
                                                        !! corresponding to frac_carb(:,isurf,islow)
  REAL, SAVE :: active_to_CO2_ref_frac  = 0.85   !! from active pool: depends on clay content  (0-1, unitless)
                                                        !! corresponding to frac_resp(:,iactive)
!$OMP THREADPRIVATE(active_to_CO2_ref_frac)
  REAL, SAVE :: slow_to_pass_ref_frac   = 0.003  !! from slow pool: depends on clay content  (0-1, unitless) 
                                                        !! corresponding to frac_carb(:,islow,ipassive)
!$OMP THREADPRIVATE(slow_to_pass_ref_frac)
  REAL, SAVE :: slow_to_CO2_ref_frac    = 0.55   !! from slow pool (0-1, unitless) 
                                                        !! corresponding to frac_resp(:,islow)
!$OMP THREADPRIVATE(slow_to_CO2_ref_frac)
  REAL, SAVE :: pass_to_active_ref_frac = 0.45   !! from passive pool (0-1, unitless)
                                                        !! corresponding to frac_carb(:,ipassive,iactive)
!$OMP THREADPRIVATE(pass_to_active_ref_frac)
  REAL, SAVE :: pass_to_slow_ref_frac   = 0.0    !! from passive pool (0-1, unitless)
                                                        !! corresponding to frac_carb(:,ipassive,islow)
!$OMP THREADPRIVATE(pass_to_slow_ref_frac)

  ! 3. Define Variable fraction from one pool to another (function of silt and clay fraction)
  REAL, SAVE :: active_to_pass_clay_frac     = 0.032
!$OMP THREADPRIVATE(active_to_pass_clay_frac)
  !! residence times in carbon pools (days)
  REAL, SAVE :: carbon_tau_iactive = 0.149   !! residence times in active pool (days)
!$OMP THREADPRIVATE(carbon_tau_iactive)
  REAL, SAVE :: carbon_tau_islow = 5.48      !! residence times in slow pool (days)
!$OMP THREADPRIVATE(carbon_tau_islow)
  REAL, SAVE :: carbon_tau_ipassive = 241.   !! residence times in passive pool (days) !DSG: years?
!$OMP THREADPRIVATE(carbon_tau_ipassive)
  REAL, SAVE, DIMENSION(3) :: flux_tot_coeff = (/ 1.2, 1.4, .75/)
!$OMP THREADPRIVATE(flux_tot_coeff)
  !! residence times in carbon pools (days)

  REAL, SAVE :: active_to_CO2_clay_silt_frac = 0.68
!$OMP THREADPRIVATE(active_to_pass_clay_frac)
  REAL, SAVE :: slow_to_pass_clay_frac   = -0.009 
!$OMP THREADPRIVATE(slow_to_pass_clay_frac)

  ! C to N target ratios of differnt pools
  REAL, SAVE ::  CN_target_iactive_ref  = 15. !! CN target ratio of active pool for soil min N = 0
  REAL, SAVE ::  CN_target_islow_ref    = 20. !! CN target ratio of slow pool for soil min N = 0
  REAL, SAVE ::  CN_target_ipassive_ref = 10. !! CN target ratio of passive pool for soil min N = 0
  REAL, SAVE ::  CN_target_isurface_ref = 20. !! CN target ratio of surface pool for litter nitrogen content = 0

  REAL, SAVE ::  CN_target_iactive_Nmin  = -6.  !! CN target ratio change per mineral N unit (g m-2) for active pool 
  REAL, SAVE ::  CN_target_islow_Nmin    = -4.  !! CN target ratio change per mineral N unit (g m-2) for slow pool 
  REAL, SAVE ::  CN_target_ipassive_Nmin = -1.5 !! CN target ratio change per mineral N unit (g m-2) for passive pool 
  REAL, SAVE ::  CN_target_isurface_pnc  = -5.  !! CN target ratio change per plant nitrogen content unit (%) for surface pool


 !NHsink
  !REAL, SAVE ::  NP_target_iactive_ref  = 6.  !! NP target ratio of active pool from CENTURY model 
  REAL, SAVE ::  NP_target_iactive_ref  = 14.  !! NP target ratio of active pool from CENTURY model 
 !NHsink
  REAL, SAVE ::  NP_target_islow_ref    = 14. !! NP target ratio of slow pool 
 !NHsink
  !REAL, SAVE ::  NP_target_ipassive_ref = 6.  !! NP target ratio of passive pool 
  REAL, SAVE ::  NP_target_ipassive_ref = 14.  !! NP target ratio of passive pool 
 !NHsink
  REAL, SAVE ::  NP_target_isurface_ref = 14. !! NP target ratio of surface pool 


  !! Turnover in SOM pools (year-1)
  REAL, SAVE :: som_turn_isurface = 6.0           !! turnover of surface pool (year-1)
!$OMP THREADPRIVATE(som_turn_isurface)
  REAL, SAVE :: som_turn_iactive  = 7.3           !! turnover of active pool (year-1)
!$OMP THREADPRIVATE(som_turn_iactive)
  REAL, SAVE :: som_turn_islow    = 0.2           !! turnover of slow pool (year-1)
!$OMP THREADPRIVATE(som_turn_islow)
!DSG  REAL, SAVE :: som_turn_ipassive = 0.0045        !! turnover of passive pool (year-1)
  REAL, SAVE :: som_turn_ipassive = 0.002         !! turnover of passive pool (year-1)

!$OMP THREADPRIVATE(som_turn_ipassive)


  REAL, SAVE :: som_turn_iactive_clay_frac = 0.75 !! clay-dependant parameter impacting on turnover rate of active pool 
                                                         !! Tm parameter of Parton et al. 1993 (-)
!$OMP THREADPRIVATE(som_turn_iactive_clay_frac)

  !
  ! stomate_turnover.f90
  !

  ! 3. Coefficients of equations

  REAL, SAVE :: new_turnover_time_ref = 20. !!(days)
!$OMP THREADPRIVATE(new_turnover_time_ref)
  REAL, SAVE :: leaf_age_crit_tref = 20.    !! (C)
!$OMP THREADPRIVATE(leaf_age_crit_tref)
  REAL, SAVE, DIMENSION(3) :: leaf_age_crit_coeff = (/ 1.5, 0.75, 10./) !! (unitless)
!$OMP THREADPRIVATE(leaf_age_crit_coeff)


  !
  ! stomate_vmax.f90
  !
 
  ! 1. Scalar

  REAL, SAVE :: vmax_offset = 0.3        !! minimum leaf efficiency (unitless)
!$OMP THREADPRIVATE(vmax_offset)
  REAL, SAVE :: leafage_firstmax = 0.03  !! relative leaf age at which efficiency
                                                !! reaches 1 (unitless)
!$OMP THREADPRIVATE(leafage_firstmax)
  REAL, SAVE :: leafage_lastmax = 0.5    !! relative leaf age at which efficiency
                                                !! falls below 1 (unitless)
!$OMP THREADPRIVATE(leafage_lastmax)
  REAL, SAVE :: leafage_old = 1.         !! relative leaf age at which efficiency
                                                !! reaches its minimum (vmax_offset) 
                                                !! (unitless)
!$OMP THREADPRIVATE(leafage_old)


  !
  ! nitrogen_dynamics (in stomate_soilcarbon.f90) 
  !

  ! 0. Constants
  REAL, PARAMETER :: D_air = 1.73664     !! Oxygen diffusion rate in the air = 0.07236 m2/h
                                               !! from Table 2 of Li et al, 2000
                                               !! (m**2/day)

  REAL, PARAMETER :: C_molar_mass = 12  !! Carbon Molar mass (gC mol-1) 

  REAL, PARAMETER :: Pa_to_hPa    = 0.01      !! Conversion factor from Pa to hPa (-)
  REAL, PARAMETER :: V_O2         = 0.209476  !! Volumetric fraction of O2 in air (-)

  REAL, PARAMETER :: pk_NH4 = 9.25      !! The negative logarithm of the acid dissociation constant K_NH4     
                                               !! See Table 4 of Li et al. 1992 and Appendix A of Zhang et al. 2002    

!  REAL, SAVE                                      :: BNF_scal = 1.6   ! maximum BNF rate  see biologicalN2fixation for justification
!NHsink
  REAL, SAVE                                      :: BNF_scal = 3.2   ! maximum BNF rate  see biologicalN2fixation for justification
!NHsink

  REAL, SAVE                                      :: BNF_coef = -0.003  ! empirical factor

  !
  !  stomate_phosphorus.f90) 
  !
!DSG: the following turnover times are set to the ones of "normal" mineralization (=C decomposition)
!     exception of ipassive: there is evidence that the fraction of SOM which is bound to
!     minerals is rich in nutrients (= it cannot be accessed by biota ) (Tipping et al., 2016)
  !! Biochemical mineralization: 
  REAL, SAVE :: bcm_turn_isurface = 12.0           !! BCM turnover of P in surface pool (year-1)
!$OMP THREADPRIVATE(bcm_turn_isurface)
  REAL, SAVE :: bcm_turn_iactive  = 14.6           !! BCM turnover of P in active pool (year-1)
!$OMP THREADPRIVATE(bcm_turn_iactive)
  REAL, SAVE :: bcm_turn_islow    = 0.6           !! BCM turnover of P in slow pool (year-1)
!$OMP THREADPRIVATE(bcm_turn_islow)
  REAL, SAVE :: bcm_turn_ipassive = 0.006          !! BCM turnover of P in passive pool (year-1)
                                                         !! It's set to zero according to Tipping et al (201?)  
                                                         !! and to speed up spinup of global model
                                                         !
!$OMP THREADPRIVATE(bcm_turn_ipassive)

  REAL, SAVE :: sorb_tune   = 1.0                 !! tuning factor for sorbed fraction of labile P (dirty)
                                     
  ! 1. Scalar

  ! Coefficients for defining maximum porosity
  ! From Saxton, K.E., Rawls, W.J., Romberger, J.S., Papendick, R.I., 1986
  ! Estimationg generalized soil-water characteristics from texture. 
  ! Soil Sci. Soc. Am. J. 50, 1031-1036
  ! Cited in Table 5 (page 444) of
  ! Y. Pachepsky, W.J. Rawls
  ! Development of Pedotransfer Functions in Soil Hydrology
  ! Elsevier, 23 nov. 2004 - 542 pages
  ! http://books.google.fr/books?id=ar_lPXaJ8QkC&printsec=frontcover&hl=fr#v=onepage&q&f=false
  REAL, SAVE :: h_saxton = 0.332          !! h coefficient 
!$OMP THREADPRIVATE(h_saxton)
  REAL, SAVE :: j_saxton = -7.251*1e-4    !! j coefficient 
!$OMP THREADPRIVATE(j_saxton)
  REAL, SAVE :: k_saxton = 0.1276         !! k coefficient
!$OMP THREADPRIVATE(k_saxton)

  ! Values of the power used in the equation defining the diffusion of oxygen in soil
  ! from Table 2 of Li et al, 2000
  REAL, SAVE :: diffusionO2_power_1 = 3.33 !! (unitless)
!$OMP THREADPRIVATE(diffusionO2_power_1)
  REAL, SAVE :: diffusionO2_power_2 = 2.0  !! (unitless) 
!$OMP THREADPRIVATE(diffusionO2_power_2)

  ! Temperature-related Factors impacting on Oxygen diffusion rate
  ! From eq. 2 of Table 2 (Li et al, 2000)
  REAL, SAVE ::   F_nofrost = 1.2          !! (unitless)
!$OMP THREADPRIVATE(F_nofrost)
  REAL, SAVE ::   F_frost   = 0.8          !! (unitless)
!$OMP THREADPRIVATE(F_frost)

  ! Coefficients used in the calculation of Volumetric fraction of anaerobic microsites 
  ! a and b constants are not specified in Li et al., 2000
  ! S. Zaehle used a=0.85 and b=1 without mention to any publication
  REAL, SAVE ::   a_anvf    = 0.85 !! (-)
!$OMP THREADPRIVATE(a_anvf)
  REAL, SAVE ::   b_anvf    = 1.   !! (-)
!$OMP THREADPRIVATE(b_anvf)

  ! Coefficients used in the calculation of the Fraction of adsorbed NH4+
  ! Li et al. 1992, JGR, Table 4
  REAL, SAVE ::   a_FixNH4 = 0.41  !! (-)
!$OMP THREADPRIVATE(a_FixNH4)
  REAL, SAVE ::   b_FixNH4 = -0.47 !! (-)
!$OMP THREADPRIVATE(b_FixNH4)
  REAL, SAVE ::   clay_max = 0.63  !! (-)
!$OMP THREADPRIVATE(clay_max)

  ! Coefficients used in the calculation of the Response of Nitrification
  ! to soil moisture
  ! Zhang et al. 2002, Ecological Modelling, appendix A, page 101
  REAL, SAVE ::   fw_0 =  -0.0243  !! (-)
!$OMP THREADPRIVATE(fw_0)
  REAL, SAVE ::   fw_1 =   0.9975  !! (-)
!$OMP THREADPRIVATE(fw_1)
  REAL, SAVE ::   fw_2 =  -5.5368  !! (-)
!$OMP THREADPRIVATE(fw_2)
  REAL, SAVE ::   fw_3 =  17.651   !! (-)
!$OMP THREADPRIVATE(fw_3)
  REAL, SAVE ::   fw_4 = -12.904   !! (-)
!$OMP THREADPRIVATE(fw_4)

  ! Coefficients used in the calculation of the Response of Nitrification
  ! to Temperature
  ! Zhang et al. 2002, Ecological Modelling, appendix A, page 101
  REAL, SAVE ::   ft_nit_0 =  -0.0233 !! (-)
!$OMP THREADPRIVATE(ft_nit_0)
  REAL, SAVE ::   ft_nit_1 =   0.3094 !! (-)
!$OMP THREADPRIVATE(ft_nit_1)
  REAL, SAVE ::   ft_nit_2 =  -0.2234 !! (-)
!$OMP THREADPRIVATE(ft_nit_2)
  REAL, SAVE ::   ft_nit_3 =   0.1566 !! (-)
!$OMP THREADPRIVATE(ft_nit_3)
  REAL, SAVE ::   ft_nit_4 =  -0.0272 !! (-)
!$OMP THREADPRIVATE(ft_nit_4)

  ! Coefficients used in the calculation of the Response of Nitrification
  ! to pH
  ! Zhang et al. 2002, Ecological Modelling, appendix A, page 101
  REAL, SAVE ::   fph_0 = -1.2314  !! (-)
!$OMP THREADPRIVATE(fph_0)
  REAL, SAVE ::   fph_1 = 0.7347   !! (-)
!$OMP THREADPRIVATE(fph_1)
  REAL, SAVE ::   fph_2 = -0.0604  !! (-)
!$OMP THREADPRIVATE(fph_2)

  ! Coefficients used in the calculation of the response of NO2 or NO 
  ! production during nitrificationof to Temperature
  ! Zhang et al. 2002, Ecological Modelling, appendix A, page 102
  REAL, SAVE ::   ftv_0 = 2.72   !! (-)
!$OMP THREADPRIVATE(ftv_0)
  REAL, SAVE ::   ftv_1 = 34.6   !! (-)
!$OMP THREADPRIVATE(ftv_1)
  REAL, SAVE ::   ftv_2 = 9615.  !! (-) 
!$OMP THREADPRIVATE(ftv_2)

 !DSGminN REAL, SAVE ::   k_nitrif = 0.2         !! Nitrification rate at 20 ◦C and field capacity (day-1)
 !DSGminN                                               !! from Schmid et al., 2001
 REAL, SAVE ::   k_nitrif = 1.2         !! Nitrification rate at 20 ◦C and field capacity (day-1)
                                               !! from OCN
 !DSGminN 
!$OMP THREADPRIVATE(k_nitrif)
 REAL, SAVE ::   Rn2oN = 0.0008
!$OMP THREADPRIVATE(Rn2on)
 REAL, SAVE ::   RnoN = 0.02
!$OMP THREADPRIVATE(Rnon)
 REAL, SAVE ::   scal_anvf = 0.8
 REAL, SAVE ::   scal_ph = 0.0

  REAL, SAVE ::   n2o_nitrif_p = 0.0006  !! Reference n2o production per N-NO3 produced g N-N2O  (g N-NO3)-1
                                                !! From Zhang et al., 2002 - Appendix A p. 102
!$OMP THREADPRIVATE(n2o_nitrif_p)
  REAL, SAVE ::   no_nitrif_p = 0.0025   !! Reference NO production per N-NO3 produced g N-NO  (g N-NO3)-1
                                                !! From Zhang et al., 2002 - Appendix A p. 102
!$OMP THREADPRIVATE(no_nitrif_p)

  ! NO production from chemodenitrification
  ! based on Kesik et al., 2005, Biogeosciences
  ! Coefficients used in the calculation of the Response to Temperature
  REAL, SAVE ::   chemo_t0  = -31494. !! (-)
!$OMP THREADPRIVATE(chemo_t0)
  ! Coefficients use in the calculation of the Response to pH
  REAL, SAVE ::   chemo_ph0 = -1.62   !! (-)
!$OMP THREADPRIVATE(chemo_ph0)
  ! Coefficients used in the calculation of NO production from chemodenitrification
  REAL, SAVE ::   chemo_0   = 30.     !! (-)
!$OMP THREADPRIVATE(chemo_0)
  REAL, SAVE ::   chemo_1   = 16565.  !! (-)
!$OMP THREADPRIVATE(chemo_1)

  ! Denitrification processes
  ! Li et al, 2000, JGR Table 4 eq 1, 2 and 4
  !
  ! Coefficients used in the Temperature response of 
  ! relative growth rate of total denitrifiers - Eq. 2 Table 4 of Li et al., 2000
  REAL, SAVE ::   ft_denit_0 = 2.     !! (-)
!$OMP THREADPRIVATE(ft_denit_0)
  REAL, SAVE ::   ft_denit_1 = 22.5   !! (-)
!$OMP THREADPRIVATE(ft_denit_1)
  REAL, SAVE ::   ft_denit_2 = 10.    !! (-)
!$OMP THREADPRIVATE(ft_denit_2)
  !
  ! Coefficients used in the pH response of 
  ! relative growth rate of total denitrifiers - Eq. 2 Table 4 of Li et al., 2000
  REAL, SAVE ::   fph_no3_0  = 4.25    !! (-) 
!$OMP THREADPRIVATE(fph_no3_0)
  REAL, SAVE ::   fph_no3_1  = 0.5     !! (-)
!$OMP THREADPRIVATE(fph_no3_1)
  REAL, SAVE ::   fph_no_0  = 5.25     !! (-)
!$OMP THREADPRIVATE(fph_no_0)
  REAL, SAVE ::   fph_no_1  = 1.       !! (-)
!$OMP THREADPRIVATE(fph_no_1)
  REAL, SAVE ::   fph_n2o_0  = 6.25    !! (-)
!$OMP THREADPRIVATE(fph_n2o_0)
  REAL, SAVE ::   fph_n2o_1  = 1.5     !! (-)
!$OMP THREADPRIVATE(fph_n2o_1)

  REAL, SAVE ::   Kn = 0.083           !! Half Saturation of N oxydes (kgN/m3)
                                              !! Table 4 of Li et al., 2000
!$OMP THREADPRIVATE(Kn)

  ! Maximum Relative growth rate of Nox denitrifiers
  ! Eq.1 Table 4 Li et al., 2000 
  REAL, SAVE ::   mu_no3_max = 0.67   !! (hour-1)
!$OMP THREADPRIVATE(mu_no3_max)
  REAL, SAVE ::   mu_no_max  = 0.34   !! (hour-1)
!$OMP THREADPRIVATE(mu_no_max)
  REAL, SAVE ::   mu_n2o_max = 0.34   !! (hour-1)
!$OMP THREADPRIVATE(mu_n2o_max)

  ! Maximum growth yield of NOx denitrifiers on N oxydes
  ! Table 4 Li et al., 2000
  REAL, SAVE ::   Y_no3 = 0.401 !! (kgC / kgN)
!$OMP THREADPRIVATE(Y_no3)
  REAL, SAVE ::   Y_no  = 0.428 !! (kgC / kgN)
!$OMP THREADPRIVATE(Y_no)
  REAL, SAVE ::   Y_n2o = 0.151 !! (kgC / kgN)
!$OMP THREADPRIVATE(Y_n2o)

  ! Maintenance coefficient on N oxyde
  ! Table 4 Li et al., 2000
  REAL, SAVE ::   M_no3 = 0.09   !! (kgN / kgC / hour)
!$OMP THREADPRIVATE(M_no3)
  REAL, SAVE ::   M_no  = 0.035  !! (kgN / kgC / hour)
!$OMP THREADPRIVATE(M_no)
  REAL, SAVE ::   M_n2o = 0.079  !! (kgN / kgC / hour)
!$OMP THREADPRIVATE(M_n2o)

        
  REAL, SAVE ::   Maint_c = 0.0076    !! Maintenance coefficient of carbon (kgC/kgC/h)
                                        !! Table 4 Li et al., 2000
!$OMP THREADPRIVATE(Maint_c)
  REAL, SAVE ::   Yc = 0.503     !! Maximum growth yield on soluble carbon (kgC/kgC)
                                        !! Table 4 Li et al., 2000
!$OMP THREADPRIVATE(Yc)

  !! Coefficients used in the eq. defining the response of N-emission to clay fraction (-)
  !! from  Table 4, Li et al. 2000
  REAL, SAVE ::   F_clay_0 = 0.13   
!$OMP THREADPRIVATE(F_clay_0) 
  REAL, SAVE ::   F_clay_1 = -0.079
!$OMP THREADPRIVATE(F_clay_1)


  REAL, SAVE ::   ratio_nh4_fert = 0.875  !! Proportion of ammonium in the fertilizers (ammo-nitrate) 
                                                 !! = 7./8. (-)
!$OMP THREADPRIVATE(ratio_nh4_fert)

  ! 2. Arrays
  REAL, SAVE, DIMENSION(2)  :: vmax_n_uptake = (/ 5.4 , 5.4 /) !! Vmax of nitrogen uptake by plants 
                                                                      !! for Ammonium (ind.1) and Nitrate (ind.2)
                                                                      !! (in umol (g DryWeight_root)-1 h-1)
                                                                      !! from Zaehle & Friend (2010) "calibrated"
 !REAL, SAVE, DIMENSION(2)  :: vmax_n_uptake = (/ 3. , 3. /)  !! original Kronzucker et al. (1995,1996)
                                                                      
!$OMP THREADPRIVATE(vmax_n_uptake)

  REAL, SAVE, DIMENSION(2)  :: K_N_min = (/ 98., 98. /)    !! [NH4+] (resp. [NO3-]) for which the Nuptake 
                                                                  !! equals vmax/2.   (umol per litter)
                                                                  !! from Zaehle & Friend (2010) "calibrated"

 ! REAL, SAVE, DIMENSION(2)  :: K_N_min = (/ 30., 30. /)    !! original values from Kronzucker, 1995,1996
!$OMP THREADPRIVATE(K_N_min)

  REAL, SAVE, DIMENSION(2)  :: low_K_N_min = (/ 0.0002, 0.0002 /) !! Rate of N uptake not associated with 
                                                                         !! Michaelis- Menten Kinetics for Ammonium 
                                                                         !! (ind.1) and Nitrate (ind.2)
                                                                         !! from Kronzucker, 1995 ((umol)-1)
!$OMP THREADPRIVATE(low_K_N_min)


  !! Other N-related parameters
  REAL, SAVE                                  :: Dmax = 0.25      !! Parameter te be clarified (what it is, units, ...)
                                                                         !! used in stomate_growth_fun_all

  REAL, SAVE :: reserve_time_tree = 30.     !! Maximum number of days during which
                                                   !! carbohydrate reserve may be used for 
                                                   !! trees (days)
!$OMP THREADPRIVATE(reserve_time_tree)
  REAL, SAVE :: reserve_time_grass = 20.    !! Maximum number of days during which
                                                   !! carbohydrate reserve may be used for 
                                                   !! grasses (days)
!$OMP THREADPRIVATE(reserve_time_grass)

  !
  ! stomate_phosphorus.f90
  !

  ! 1. Scalar
  REAL, SAVE                :: vmax_P_uptake =  1.39              !! Vmax of phosphorus uptake by plants 
                                                                         !! by high-affinty/activity system 
                                                                         !! (in umol (g DryWeight_root)-1 h-1)
                                                                         !!  HAS:  from Bouma et al. (2001) (average of both species) 3.96 
                                  
!$OMP THREADPRIVATE(vmax_P_uptake)

  REAL, SAVE                :: K_P_min = 3.                       !! dissolved P concentraion for which the Puptake 
                                                                         !! equals vmax/2.   (umol per litter)
                                                                         !! by low-affinty/activity system (ind.1) 
                                                                         !! and high-affinty/activity system (ind.2)
                                                                         !! from Schachtman et al. (1998), Page 448: 
                                                                         !! averages of range given for
                                                                         !! high affinity system (3-7) 
!$OMP THREADPRIVATE(K_P_min)

  REAL, SAVE                :: low_K_P_min = 0.01                 !! linear increase in uptake for high P concentrations (chosen to fit the uptake
                                                                         !! at high concentrations from Zhang et al. (2009)

                                                                        
!$OMP THREADPRIVATE(low_K_P_min)


  !
  ! stomate_season.f90 
  !

  ! 1. Scalar

  REAL, SAVE :: gppfrac_dormance = 0.2  !! report maximal GPP/GGP_max for dormance (0-1, unitless)
!$OMP THREADPRIVATE(gppfrac_dormance)
  REAL, SAVE :: tau_climatology = 20.   !! tau for "climatologic variables (years)
!$OMP THREADPRIVATE(tau_climatology)
  REAL, SAVE :: hvc1 = 0.019            !! parameters for herbivore activity (unitless)
!$OMP THREADPRIVATE(hvc1)
  REAL, SAVE :: hvc2 = 1.38             !! parameters for herbivore activity (unitless)
!$OMP THREADPRIVATE(hvc2)
  REAL, SAVE :: leaf_frac_hvc = 0.33    !! leaf fraction (0-1, unitless)
!$OMP THREADPRIVATE(leaf_frac_hvc)
  REAL, SAVE :: tlong_ref_max = 303.1   !! maximum reference long term temperature (K)
!$OMP THREADPRIVATE(tlong_ref_max)
  REAL, SAVE :: tlong_ref_min = 253.1   !! minimum reference long term temperature (K)
!$OMP THREADPRIVATE(tlong_ref_min)

  ! 3. Coefficients of equations

  REAL, SAVE :: ncd_max_year = 3.
!$OMP THREADPRIVATE(ncd_max_year)
  REAL, SAVE :: gdd_threshold = 5.
!$OMP THREADPRIVATE(gdd_threshold)
  REAL, SAVE :: green_age_ever = 2.
!$OMP THREADPRIVATE(green_age_ever)
  REAL, SAVE :: green_age_dec = 0.5
!$OMP THREADPRIVATE(green_age_dec)

  INTEGER, SAVE :: ncirc = 1                  !! Number of circumference classes used to calculate C allocation
!$OMP THREADPRIVATE(ncirc)
  REAL, PARAMETER :: kilo_to_unit = 1.0E03   !! Convert Kilo to unit


LOGICAL, SAVE :: lbypass_cc = .FALSE.                      !! Set to true for a temporary patch of a known bug, though the underlying
!$OMP THREADPRIVATE(lbypass_cc)
!DSG LOGICAL, SAVE :: ld_fake_height=.FALSE.  ! a flag to turn on the statements ;DSG:  I wouldnt have guess
!DSG !$OMP THREADPRIVATE(ld_fake_height)
REAL, SAVE :: sync_threshold = 0.0001     !! The threshold above which a warning is generated when the
!$OMP THREADPRIVATE(sync_threshold)
LOGICAL,PARAMETER :: ld_biomass=.FALSE.   ! a flag to turn on debug statements
INTEGER, SAVE        :: test_pft = 4                              !! Number of PFT for which detailed output 
!$OMP THREADPRIVATE(test_pft)

INTEGER, SAVE        :: test_grid = 1                                 !! Number of the grid square for which detailed output 
!$OMP THREADPRIVATE(test_grid)

LOGICAL,PARAMETER                       :: ld_stop=.FALSE.      ! a flag to turn on some stop statements.
LOGICAL,SAVE                            :: ld_alloc=.FALSE.     ! a flag to turn on debug statements
LOGICAL,PARAMETER                       :: ld_warn=.FALSE.      ! a flag to turn on various warnings
LOGICAL,PARAMETER                       :: jr_nextstep = .FALSE.   ! set this to TRUE to activate the 
LOGICAL,PARAMETER                       :: ld_massbal=.FALSE.   ! a flag to turn on debug statements
INTEGER, PARAMETER :: ipoolchange = 5           !! change in pool size i.e. change in biomass

INTEGER, PARAMETER :: ilat2in = 4               !! incoming lateral flux i.e. N deposition for the land
INTEGER, PARAMETER :: ilat2out = 3              !! outgoing lateral flux i.e. DOC leaching for the litter routine
INTEGER, PARAMETER :: iatm2land = 1             !! atmosphere to land fluxes such as GPP and co2_2_bm
INTEGER, PARAMETER :: iland2atm = 2             !! land to atmosphere fluxes such as Rh, Ra and product decomposition
INTEGER, PARAMETER :: nmbcomp = 5               !! The total nomber of components in our mass balance check


REAL, SAVE :: max_delta_KF = 0.1          !! Maximum change in KF from one time step to another (m)
                                                   !! This is a bit arbitrary. 
!$OMP THREADPRIVATE(max_delta_KF)

 REAL, SAVE :: maint_from_gpp = 0.8        !! Some carbon needs to remain to support the growth, hence, 
                                                   !! respiration will be limited. In this case resp_maint 
                                                   !! (gC m-2 dt-1) should not be more than 80% (::maint_from_gpp) 
                                                   !! of the GPP (gC m-2 s-1)
!$OMP THREADPRIVATE(maint_from_gpp)
  REAL, PARAMETER :: m2_to_ha = 10000.       !! Conversion from m2 to hectares
  REAL, PARAMETER :: ha_to_m2 = 0.0001       !! Conversion from hectares (forestry) to m2 (rest of the code)


 ! JC MOD045 wstress_lowgrass for grasses only
 REAL, SAVE :: wstress_lowgrass = 0.1       !! minimum water stress factor for grass allocation
 REAL, SAVE :: sstress_lowgrass = 0.1       !! minimum total stress factor for grass allocation
 ! JC MOD050 root resp suppression
 REAL, SAVE :: coeff_suppress_resp = 0.1    !! default coeff makes no suppression
 REAL, SAVE :: days_senesc_crit = 10.0      !! days start to suppress maint_resp for dormancy
 REAL, SAVE :: GRM_RtoL_turn = 0.0          !! relative grass root turnover compared to leaf turnover
                                                   !! during senescence (default 1 = annual grasses)
 
 !gmjc constantes for GRM module
 LOGICAL, SAVE      ::  GRM_allow_BNF = .FALSE.      !! allowed legume BNF?
 LOGICAL, SAVE      ::  GRM_BNF_newmethod = .FALSE.  !! use new method (Lazzarotto et al. 2009) for BNF?
 LOGICAL, SAVE      ::  GRM_enable_grazing = .FALSE. !! activate GRM module?
 LOGICAL, SAVE      ::  GRM_allow_DEVSTAGE = .FALSE. !! allow grass allocation follows exactly devstage as PaSim? (should not activate)
 REAL, SAVE  ::  fix_legume_frac = zero       !! fixed legume fraction for BNF
 REAL, SAVE  ::  reserve_time_cut = 20.
 REAL, SAVE  ::  lai_happy_cut = 0.25
 REAL, SAVE  ::  tau_leafinit_cut = 10
 REAL, SAVE  ::  tau_t2m_14 = 14.
 REAL, SAVE  ::  N_effect = 0.65         !! this is the maximum additive effect of nitrogen fertilization on Vcmax

END MODULE constantes_var
