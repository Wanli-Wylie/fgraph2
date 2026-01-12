! =================================================================================================================================
! MODULE 	: constantes_soil_var
!
! CONTACT       : orchidee-help _at_ listes.ipsl.fr
!
! LICENCE       : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF         "constantes_soil_var" module contains the parameters related to soil and hydrology.
!!
!!\n DESCRIPTION : The non saturated hydraulic properties are defined from the  
!!                 formulations of van Genuchten (1980) and Mualem (1976), combined as  
!!                 explained in d'Orgeval (2006). \n
!!                 The related parameters for main soil textures (coarse, medium and fine if "fao", 
!!                 12 USDA testures if "usda") come from Carsel and Parrish (1988).
!!
!! RECENT CHANGE(S): Sonke Zaehle changed hcrit_litter value according to Shilong Piao
!!                   from 0.03 to 0.08, 080806
!!                   AD: mcw and mcf depend now on soil texture, based on Van Genuchten equations 
!!                   and classical matric potential values, and pcent is adapted
!!
!! REFERENCE(S)	:
!!- Roger A.Pielke, (2002), Mesoscale meteorological modeling, Academic Press Inc. 
!!- Polcher, J., Laval, K., Dümenil, L., Lean, J., et Rowntree, P. R. (1996).
!! Comparing three land surface schemes used in general circulation models. Journal of Hydrology, 180(1-4), 373--394.
!!- Ducharne, A., Laval, K., et Polcher, J. (1998). Sensitivity of the hydrological cycle
!! to the parametrization of soil hydrology in a GCM. Climate Dynamics, 14, 307--327. 
!!- Rosnay, P. de et Polcher, J. (1999). Modelling root water uptake in a complex land surface
!! scheme coupled to a GCM. Hydrol. Earth Syst. Sci., 2(2/3), 239--255.
!!- d'Orgeval, T. et Polcher, J. (2008). Impacts of precipitation events and land-use changes
!! on West African river discharges during the years 1951--2000. Climate Dynamics, 31(2), 249--262. 
!!- Carsel, R. and Parrish, R.: Developing joint probability distributions of soil water
!! retention characteristics, Water Resour. Res.,24, 755–769, 1988.
!!- Mualem Y (1976) A new model for predicting the hydraulic conductivity  
!! of unsaturated porous media. Water Resources Research 12(3):513-522
!!- Van Genuchten M (1980) A closed-form equation for predicting the  
!! hydraulic conductivity of unsaturated soils. Soil Sci Soc Am J, 44(5):892-898
!!
!! SVN          :
!! $HeadURL: $
!! $Date: $
!! $Revision: $
!! \n
!_ ================================================================================================================================

MODULE constantes_soil_var

  USE defprec
  !USE constantes
  USE constantes_var
  USE vertical_soil_var

  IMPLICIT NONE

  LOGICAL, SAVE             :: check_cwrr          !! To check the water balance in hydrol (true/false)
!$OMP THREADPRIVATE(check_cwrr)
  LOGICAL, SAVE             :: check_cwrr2         !! Calculate diagnostics to check the water balance in hydrol (true/false)
!$OMP THREADPRIVATE(check_cwrr2)
  LOGICAL, SAVE             :: check_waterbal      !! The check the water balance (true/false)
!$OMP THREADPRIVATE(check_waterbal)


  !! Number of soil classes

  INTEGER, PARAMETER :: ntext=3                  !! Number of soil textures (Silt, Sand, Clay)
  INTEGER, PARAMETER :: nstm=3                   !! Number of soil tiles (unitless)
  CHARACTER(LEN=30)         :: soil_classif             !! Type of classification used for the map of soil types.
                                                        !! It must be consistent with soil file given by 
                                                        !! SOILCLASS_FILE parameter.
!$OMP THREADPRIVATE(soil_classif)
  INTEGER, PARAMETER :: nscm_fao=3               !! For FAO Classification (unitless)
  INTEGER, PARAMETER :: nscm_usda=12             !! For USDA Classification (unitless)
  INTEGER, SAVE      :: nscm=nscm_fao            !! Default value for nscm
!$OMP THREADPRIVATE(nscm)

  INTEGER, SAVE      :: p_nscm=nscm_usda         !! Default value for p_nscm
!$OMP THREADPRIVATE(p_nscm)

  !! Number of lithological classes
  INTEGER, PARAMETER :: nlcm_GLiM=16              !! For GLiM Classification (unitless)
  INTEGER, SAVE      :: nlcm=nlcm_GLiM            !! Default value for nlcm
!$OMP THREADPRIVATE(nlcm)

  !! Parameters for soil thermodynamics

  REAL, SAVE :: so_capa_dry = 1.80e+6            !! Dry soil Heat capacity of soils 
                                                        !! @tex $(J.m^{-3}.K^{-1})$ @endtex 
!$OMP THREADPRIVATE(so_capa_dry)
  REAL, SAVE :: so_cond_dry = 0.40               !! Dry soil Thermal Conductivity of soils
                                                        !! @tex $(W.m^{-2}.K^{-1})$ @endtex
!$OMP THREADPRIVATE(so_cond_dry)
  REAL, SAVE :: so_capa_wet = 3.03e+6            !! Wet soil Heat capacity of soils 
                                                        !! @tex $(J.m^{-3}.K^{-1})$ @endtex
!$OMP THREADPRIVATE(so_capa_wet)
  REAL, SAVE :: so_cond_wet = 1.89               !! Wet soil Thermal Conductivity of soils 
                                                        !! @tex $(W.m^{-2}.K^{-1})$ @endtex 
!$OMP THREADPRIVATE(so_cond_wet)
  REAL, SAVE :: sn_cond = 0.3                    !! Thermal Conductivity of snow 
                                                        !! @tex $(W.m^{-2}.K^{-1})$ @endtex  
!$OMP THREADPRIVATE(sn_cond)
  REAL, SAVE :: sn_dens = 330.0                  !! Snow density for the soil thermodynamics
                                                        !! (kg/m3)
!$OMP THREADPRIVATE(sn_dens)
  REAL, SAVE :: sn_capa                          !! Heat capacity for snow 
                                                        !! @tex $(J.m^{-3}.K^{-1})$ @endtex
!$OMP THREADPRIVATE(sn_capa)
  REAL, SAVE :: water_capa = 4.18e+6             !! Water heat capacity 
                                                        !! @tex $(J.m^{-3}.K^{-1})$ @endtex
!$OMP THREADPRIVATE(water_capa)
  REAL, SAVE :: brk_capa = 2.0e+6                !! Heat capacity of generic rock
                                                        !! @tex $(J.m^{-3}.K^{-1})$ @endtex
!$OMP THREADPRIVATE(brk_capa)
  REAL, SAVE :: brk_cond = 3.0                   !! Thermal conductivity of saturated granitic rock
                                                        !! @tex $(W.m^{-1}.K^{-1})$ @endtex
!$OMP THREADPRIVATE(brk_cond)


  !! Specific parameters for the Choisnel hydrology

  REAL, SAVE :: min_drain = 0.001                !! Diffusion constant for the slow regime
                                                        !! (This is for the diffusion between reservoirs)
                                                        !! @tex $(kg.m^{-2}.dt^{-1})$ @endtex
!$OMP THREADPRIVATE(min_drain)
  REAL, SAVE :: max_drain = 0.1                  !! Diffusion constant for the fast regime 
                                                        !! @tex $(kg.m^{-2}.dt^{-1})$ @endtex
!$OMP THREADPRIVATE(max_drain)
  REAL, SAVE :: exp_drain = 1.5                  !! The exponential in the diffusion law (unitless)
!$OMP THREADPRIVATE(exp_drain)
  REAL, SAVE :: qsintcst = 0.1                   !! Transforms leaf area index into size of interception reservoir
                                                        !! (unitless)
!$OMP THREADPRIVATE(qsintcst)
  REAL, SAVE :: mx_eau_nobio = 150.              !! Volumetric available soil water capacity in nobio fractions
                                                        !! @tex $(kg.m^{-3} of soil)$ @endtex
!$OMP THREADPRIVATE(mx_eau_nobio)
  REAL, SAVE :: rsol_cste = 33.E3                !! Constant in the computation of resistance for bare soil evaporation
                                                        !! @tex $(s.m^{-2})$ @endtex
!$OMP THREADPRIVATE(rsol_cste)
  REAL, SAVE :: hcrit_litter=0.08                !! Scaling depth for litter humidity (m)
!$OMP THREADPRIVATE(hcrit_litter)


  !! Parameters specific for the CWRR hydrology.

  !!  1. Parameters for FAO Classification

  !! Parameters for soil type distribution

  REAL,DIMENSION(nscm_fao),SAVE :: soilclass_default_fao = &   !! Default soil texture distribution for fao :
 & (/ 0.28, 0.52, 0.20 /)                                             !! in the following order : COARSE, MEDIUM, FINE (unitless)
!$OMP THREADPRIVATE(soilclass_default_fao)

  REAL,PARAMETER,DIMENSION(nscm_fao) :: nvan_fao = &            !! Van Genuchten coefficient n (unitless)
 & (/ 1.89      , 1.56      , 1.31       /)                             !  RK: 1/n=1-m

  REAL,PARAMETER,DIMENSION(nscm_fao) :: avan_fao = &            !! Van Genuchten coefficient a 
  & (/ 0.0075      , 0.0036      , 0.0019       /)                     !!  @tex $(mm^{-1})$ @endtex

  REAL,PARAMETER,DIMENSION(nscm_fao) :: mcr_fao = &             !! Residual volumetric water content 
 & (/ 0.065      , 0.078      , 0.095       /)                         !!  @tex $(m^{3} m^{-3})$ @endtex

  REAL,PARAMETER,DIMENSION(nscm_fao) :: mcs_fao = &             !! Saturated volumetric water content 
 & (/ 0.41      , 0.43      , 0.41       /)                            !!  @tex $(m^{3} m^{-3})$ @endtex

  REAL,PARAMETER,DIMENSION(nscm_fao) :: ks_fao = &              !! Hydraulic conductivity at saturation 
 & (/ 1060.8      , 249.6      , 62.4       /)                         !!  @tex $(mm d^{-1})$ @endtex

! The max available water content is smaller when mcw and mcf depend on texture,
! so we increase pcent to a classical value of 80%
  REAL,PARAMETER,DIMENSION(nscm_fao) :: pcent_fao = &           !! Fraction of saturated volumetric soil moisture 
 & (/ 0.8      , 0.8      , 0.8       /)                               !! above which transpir is max (0-1, unitless)

  REAL,PARAMETER,DIMENSION(nscm_fao) :: free_drain_max_fao = &  !! Max=default value of the permeability coeff  
 & (/ 1.0      , 1.0      , 1.0       /)                               !! at the bottom of the soil (0-1, unitless)

!! We use the VG relationships to derive mcw and mcf depending on soil texture
!! assuming that the matric potential for wilting point and field capacity is
!! -150m (permanent WP) and -3.3m respectively
!! (-1m for FC for the three sandy soils following Richards, L.A. and Weaver, L.R. (1944)
!! Note that mcw GE mcr
  REAL,PARAMETER,DIMENSION(nscm_fao) :: mcf_fao = &             !! Volumetric water content at field capacity 
 & (/ 0.1218      , 0.1654      , 0.2697       /)                      !!  @tex $(m^{3} m^{-3})$ @endtex

  REAL,PARAMETER,DIMENSION(nscm_fao) :: mcw_fao = &             !! Volumetric water content at wilting point 
 & (/ 0.0657      ,  0.0884      , 0.1496      /)                      !!  @tex $(m^{3} m^{-3})$ @endtex

  REAL,PARAMETER,DIMENSION(nscm_fao) :: mc_awet_fao = &         !! Vol. wat. cont. above which albedo is cst 
 & (/ 0.25      , 0.25      , 0.25       /)                            !!  @tex $(m^{3} m^{-3})$ @endtex

  REAL,PARAMETER,DIMENSION(nscm_fao) :: mc_adry_fao = &         !! Vol. wat. cont. below which albedo is cst
 & (/ 0.1      , 0.1      , 0.1       /)                               !!  @tex $(m^{3} m^{-3})$ @endtex

  REAL,PARAMETER,DIMENSION(nscm_fao) :: SMCMAX_fao = &          !! porosity
 & (/ 0.41      , 0.43      , 0.41       /)                            !! & (/ 0.434      , 0.439      , 0.465       /) !!noah lsm

  REAL,PARAMETER,DIMENSION(nscm_fao) :: QZ_fao = &              !! QUARTZ CONTENT (SOIL TYPE DEPENDENT)
 & (/ 0.60      , 0.40      , 0.35       /)                            !! Peters et al [1998]

  REAL,PARAMETER,DIMENSION(nscm_fao) :: so_capa_dry_ns_fao = &  !! Dry soil Heat capacity of soils,J.m^{-3}.K^{-1}
 & (/ 1.34e+6      , 1.21e+6      , 1.23e+6       /)                   !! Pielke [2002, 2013]

  !!  2. Parameters for USDA Classification

  !! Parameters for soil type distribution :
  !! Sand, Loamy Sand, Sandy Loam, Silt Loam, Silt, Loam, Sandy Clay Loam, Silty Clay Loam, Clay Loam, Sandy Clay, Silty Clay, Clay

  REAL,DIMENSION(nscm_usda),SAVE :: soilclass_default_usda = &    !! Default soil texture distribution in the above order :
 & (/ 0.28, 0.52, 0.20, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)   !! Thus different from "FAO"'s COARSE, MEDIUM, FINE
                                                                         !! which have indices 3,6,9 in the 12-texture vector
  !$OMP THREADPRIVATE(soilclass_default_usda)

  REAL,PARAMETER,DIMENSION(nscm_usda) :: nvan_usda = &            !! Van Genuchten coefficient n (unitless)
 & (/ 2.68      , 2.28      , 1.89      , 1.41      , &                   !  RK: 1/n=1-m
 &    1.37      , 1.56      , 1.48      , 1.23      , &
 &    1.31      , 1.23      , 1.09      , 1.09       /)

  REAL,PARAMETER,DIMENSION(nscm_usda) :: avan_usda = &            !! Van Genuchten coefficient a 
 & (/ 0.0145      , 0.0124      , 0.0075      , 0.0020      , &          !!  @tex $(mm^{-1})$ @endtex
 &    0.0016      , 0.0036      , 0.0059      , 0.0010      , &
 &    0.0019      , 0.0027      , 0.0005      , 0.0008       /)

  REAL,PARAMETER,DIMENSION(nscm_usda) :: mcr_usda = &             !! Residual volumetric water content 
 & (/ 0.045      , 0.057      , 0.065      , 0.067      , &              !!  @tex $(m^{3} m^{-3})$ @endtex
 &    0.034      , 0.078      , 0.100      , 0.089      , &
 &    0.095      , 0.100      , 0.070      , 0.068       /)

  REAL,PARAMETER,DIMENSION(nscm_usda) :: mcs_usda = &             !! Saturated volumetric water content 
 & (/ 0.43      , 0.41      , 0.41      , 0.45      , &                  !!  @tex $(m^{3} m^{-3})$ @endtex
 &    0.46      , 0.43      , 0.39      , 0.43      , &
 &    0.41      , 0.38      , 0.36      , 0.38       /)

  REAL,PARAMETER,DIMENSION(nscm_usda) :: ks_usda = &              !! Hydraulic conductivity at saturation
 & (/ 7128.0      , 3501.6      , 1060.8      , 108.0      , &           !!  @tex $(mm d^{-1})$ @endtex
 &    60.0      , 249.6      , 314.4      , 16.8      , &
 &    62.4      , 28.8      , 4.8      , 48.0       /)

  REAL,PARAMETER,DIMENSION(nscm_usda) :: pcent_usda = &           !! Fraction of saturated volumetric soil moisture
 & (/ 0.8      , 0.8      , 0.8      , 0.8      , &                      !! above which transpir is max (0-1, unitless)
 &    0.8      , 0.8      , 0.8      , 0.8      , &
 &    0.8      , 0.8      , 0.8      , 0.8       /)

  REAL,PARAMETER,DIMENSION(nscm_usda) :: free_drain_max_usda = &  !! Max=default value of the permeability coeff 
 & (/ 1.0      , 1.0      , 1.0      , 1.0      , &                      !! at the bottom of the soil (0-1, unitless)
 &    1.0      , 1.0      , 1.0      , 1.0      , &
 &    1.0      , 1.0      , 1.0      , 1.0       /)

  REAL,PARAMETER,DIMENSION(nscm_usda) :: mcf_usda = &             !! Volumetric water content at field capacity
 & (/ 0.0493      , 0.0710      , 0.1218      , 0.2402      , &          !!  @tex $(m^{3} m^{-3})$ @endtex
      0.2582      , 0.1654      , 0.1695      , 0.3383      , &
      0.2697      , 0.2672      , 0.3370      , 0.3469       /)
  
  REAL,PARAMETER,DIMENSION(nscm_usda) :: mcw_usda = &             !! Volumetric water content at wilting point
 & (/ 0.0450      , 0.0570      , 0.0657      , 0.1039      , &          !!  @tex $(m^{3} m^{-3})$ @endtex
      0.0901      , 0.0884      , 0.1112      , 0.1967      , &
      0.1496      , 0.1704      , 0.2665      , 0.2707       /)

  REAL,PARAMETER,DIMENSION(nscm_usda) :: mc_awet_usda = &         !! Vol. wat. cont. above which albedo is cst
 & (/ 0.25      , 0.25      , 0.25      , 0.25      , &                  !!  @tex $(m^{3} m^{-3})$ @endtex
 &    0.25      , 0.25      , 0.25      , 0.25      , &
 &    0.25      , 0.25      , 0.25      , 0.25       /)

  REAL,PARAMETER,DIMENSION(nscm_usda) :: mc_adry_usda = &         !! Vol. wat. cont. below which albedo is cst
 & (/ 0.1      , 0.1      , 0.1      , 0.1      , &                      !!  @tex $(m^{3} m^{-3})$ @endtex
 &    0.1      , 0.1      , 0.1      , 0.1      , &
 &    0.1      , 0.1      , 0.1      , 0.1       /)


  !!  2.1 Parameters for USDA Classification (this time using soil order classification and not silt/clay)

  !! Parameters for soil orders (n=12):
  ! alphabetical sequence: 1=Alfisols, 2=Andisols, 3=Aridisols, 4=Entisols, 
  !                        5=Gelisols, 6=Histosols, 7=Inceptisols, 8=Mollisols, 
  !                        9=Oxisols, 10=Spodosols, 11=Ultisols, 12=Vertisols

  ! DSG Attention: the nscm_usda here correspond to the 12 USDA soil order not to
  ! 12 texture classe (like for the hydrological parameters). This inconsistency could
  ! lead to soil order related parameters which are not consistent with each
  ! other. The handling of the hydroligical (& other) parameters should also use USDA soil
  ! order and derive in a second step texture information to avoid problems.

! PARAMETERS FOR FREUNDLICH P SORPTION: from Kavkic et al. (in prep.) based on
! Achat et al (2016) data compilation incl croplands.
! Marko distinguishs three classes: Molisols, Oxisols, and rest. - This could be
! revised later 
  REAL,PARAMETER,DIMENSION(nscm_usda) :: K_freund_usda   = &                          !!   Freundlich Isotherm coefficient [L(water) kg-1(soil)]
 & (/ 72.3452646248      , 72.3452646248      , 72.3452646248      , 72.3452646248      , &  !!   Values from Kvacik et al (in prep)
 &    72.3452646248      , 72.3452646248      , 72.3452646248      , 185.458103612      , &
 &    348.899878516      , 72.3452646248      , 72.3452646248      , 72.3452646248       /)

  REAL,PARAMETER,DIMENSION(nscm_usda) :: n_freund_usda = &        !! Freundlich Isotherm exponent  [ ]
 & (/ 1.590280      ,1.590280      , 1.590280      , 1.590280      , &   !! Values from Kvacik et al (in prep)
 &    1.590280      ,1.590280      , 1.590280      , 1.590280      , &  
 &    1.590280      ,1.590280      , 1.590280      , 1.590280       /)
! PARAMETERS FOR FREUNDLICH P SORPTION - END

  REAL,PARAMETER,DIMENSION(nscm_usda) :: tau_sorb_usda = &      !!  turnover constant for P losses to stronlgy sorption 
!tauSorb_XtrTopics_slow:
! & (/ 36500.,36500., 36500., 36500., &                                   !!  [days]
! &    36500.,36500., 36500., 36500., &                                   !!  calibrated/tuned !! 
! &    9125. ,36500., 36500., 36500. /)
! REF:
! & (/ 18250.,18250., 18250., 18250., &                                   !!  [days]
! &    18250.,18250., 18250., 18250., &                                   !!  calibrated/tuned !! 
! &    9125. ,18250., 18250., 18250. /)
! Tausorb_New
 & (/ 36500.,36500., 36500., 36500., &                                   !!  [days]
 &    36500.,36500., 36500., 36500., &                                   !!  calibrated/tuned !! 
 &    18250. ,36500., 36500., 36500. /)


  !!  3.Parameters of GLiM lithology ( n=16)
  !! DSG: Currently there are only parameters for the GLiM classification; 
  !!      nonetheless, I follow the soil parameters to ensure new lithological 
  !!      classifications could be addded in the future.

! alphabetical sequence of lithological classes: 
!    ! 1  - evaporites
!    ! 2  - ice & glaciers
!    ! 3  - metamorphics
!    ! 4  - no data
!    ! 5  - acid plutonic rocks
!    ! 6  - basic plutonic rocks
!    ! 7  - intermediat plutonic rocks
!    ! 8  - pyroclastics
!    ! 9  - carbonate sedimentary rocks 
!    ! 10 - mixed sedimentary rocks 
!    ! 11 - siliciclastic sedimentary rocks
!    ! 12 - unconsolidated sediments
!    ! 13 - acid volcanic rocks
!    ! 14 - basic volcanic rocks
!    ! 15 - intermediate volcanic rocks
!    ! 16 - water bodies

  REAL,PARAMETER,DIMENSION(nlcm_GLiM) :: Ea_GLiM = &           !! Activation energy for Arrhenius term [J/mol]
 & (/ undef     , undef     , undef     , 6.E4      ,     &           !! in the weathering model by Hartmann et al. (2014)
 &    6.E4      , 4.E4      , 6.E4      , 4.E4      ,     &           !! ;parameter values from Goll et al. (2014)
 &    undef     , 6.E4      , 6.E4      , 6.E4      ,     & 
 &    6.E4      , 5.E4      , 5.E4      , undef      /)



  REAL,PARAMETER,DIMENSION(nlcm_GLiM) :: bS_emp_GLiM =      &  !! empirical factor to relate P release by silicate weathering to runoff
 & (/ undef        , undef        , undef        ,  1.25E-5      , &  !! ;parameter values from Goll et al. (2014)
 &    1.77E-5      , 1.42E-5      , 1.77E-5      , 16.16E-5      , &  !! [g(P)/m(H2O]? <- this needs to be confirmed from global maps; I got it w/o units 
 &    undef        , 1.51E-5      , 1.76E-5      ,  1.74E-5      , & 
 &    0.79E-5      , 8.52E-5      , 8.52E-5      , undef      /)


  REAL,PARAMETER,DIMENSION(nlcm_GLiM) :: bC_emp_GLiM =   &     !! empirical factor to relate P release by carbonate weathering to runoff
 & (/ 3.96E-5      , undef        , undef        ,  undef,      &     !! ;parameter values from Goll et al. (2014)
 &    undef        , undef        , undef        ,  undef,      &
 &    7.26E-5      , 1.72E-5      , undef        ,  undef,      &
 &    undef        , undef        , undef        ,  undef/)


  REAL,PARAMETER,DIMENSION(nscm_usda) :: SMCMAX_usda = &          !! porosity
 & (/ 0.43      , 0.41      , 0.41      , 0.45      , &
 &    0.46      , 0.43      , 0.39      , 0.43      , &
 &    0.41      , 0.38      , 0.36      , 0.38       /)
 
  REAL,PARAMETER,DIMENSION(nscm_usda) :: QZ_usda = &              !! QUARTZ CONTENT (SOIL TYPE DEPENDENT)
 & (/ 0.92      , 0.82      , 0.60      , 0.25      , &
 &    0.10      , 0.40      , 0.60      , 0.10      , &
 &    0.35      , 0.52      , 0.10      , 0.25       /)                  !! Peters et al [1998]

  REAL,PARAMETER,DIMENSION(nscm_usda) :: so_capa_dry_ns_usda = &  !! Dry soil Heat capacity of soils,J.m^{-3}.K^{-1}
 & (/ 1.47e+6      , 1.41e+6      , 1.34e+6      , 1.27e+6      , &
 &    1.21e+6      , 1.21e+6      , 1.18e+6      , 1.32e+6      , &
 &    1.23e+6      , 1.18e+6      , 1.15e+6      , 1.09e+6       /)      !! Pielke [2002, 2013]
  
  !! Parameters for the numerical scheme used by CWRR

  INTEGER, PARAMETER :: imin = 1                                 !! Start for CWRR linearisation (unitless)
  INTEGER, PARAMETER :: nbint = 50                               !! Number of interval for CWRR linearisation (unitless)
  INTEGER, PARAMETER :: imax = nbint+1                           !! Number of points for CWRR linearisation (unitless)
  REAL, PARAMETER    :: w_time = 1.0                             !! Time weighting for CWRR numerical integration (unitless)


  !! Variables related to soil freezing, in thermosoil : 
  LOGICAL, SAVE        :: ok_Ecorr                    !! Flag for energy conservation correction
  LOGICAL, SAVE        :: ok_freeze_thermix           !! Flag to activate thermal part of the soil freezing scheme
  LOGICAL, SAVE        :: read_reftemp                !! Flag to initialize soil temperature using climatological temperature
  REAL, SAVE    :: poros                       !! Soil porosity (from USDA classification, mean value)(-)
  REAL, SAVE    :: fr_dT                       !! Freezing window (K)

  !! Variables related to soil freezing, in hydrol : 
  LOGICAL, SAVE        :: ok_freeze_cwrr              !! CWRR freezing scheme by I. Gouttevin
  LOGICAL, SAVE        :: ok_thermodynamical_freezing !! Calculate frozen fraction thermodynamically


! Enhanced weathering:
  REAL, PARAMETER, DIMENSION(nminerals) :: EW_targetStock =  &    ! Mineral stock size [g(rock)/m2]
 & (/ 0.0      , 3000.0      /)
  REAL, PARAMETER, DIMENSION(nminerals) :: EW_targetSize  =  &    ! Mineral grain size  [mu m]
 & (/ 20.      , 20.0      /)
! Enhanced weathering:


END MODULE constantes_soil_var
