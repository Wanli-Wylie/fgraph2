! =================================================================================================================================
! MODULE       : constantes_mtc
!
! CONTACT      : orchidee-help _at_ listes.ipsl.fr
!
! LICENCE      : IPSL (2011)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF         This module contains the standard values of the parameters for the 13 metaclasses of vegetation used by ORCHIDEE.
!!
!!\n DESCRIPTION: None
!!
!! RECENT CHANGE(S): 
!!
!! REFERENCE(S)	: 
!! - Kuppel, S. (2012): Doctoral Thesis, Assimilation de mesures de flux turbulents d'eau et de carbone dans un modèle de la biosphère 
!! continentale 
!! - Kuppel, S., Peylin, P., Chevallier, F., Bacour, C., Maignan, F., and Richardson, A. D. (2012). Constraining a global ecosystem
!! model with multi-site eddy-covariance data, Biogeosciences, 9, 3757-3776, DOI 10.5194/bg-9-3757-2012.
!! - Wohlfahrt, G., M. Bahn, E. Haubner, I. Horak, W. Michaeler, K.Rottmar, U. Tappeiner, and A. Cemusca, 1999: Inter-specific 
!! variation of the biochemical limitation to photosynthesis and related leaf traits of 30 species from mountain grassland 
!! ecosystems under different land use. Plant Cell Environ., 22, 12811296.
!! - Malhi, Y., Doughty, C., and Galbraith, D. (2011). The allocation of ecosystem net primary productivity in tropical forests, 
!! Philosophical Transactions of the Royal Society B-Biological Sciences, 366, 3225-3245, DOI 10.1098/rstb.2011.0062.
!! - Earles, J. M., Yeh, S., and Skog, K. E. (2012). Timing of carbon emissions from global forest clearance, Nature Climate Change, 2, 
!! 682-685, Doi 10.1038/Nclimate1535.
!! - Piao, S. L., Luyssaert, S., Ciais, P., Janssens, I. A., Chen, A. P., Cao, C., Fang, J. Y., Friedlingstein, P., Luo, Y. Q., and 
!! Wang, S. P. (2010). Forest annual carbon cost: A global-scale analysis of autotrophic respiration, Ecology, 91, 652-661, 
!! Doi 10.1890/08-2176.1.
!! - Verbeeck, H., Peylin, P., Bacour, C., Bonal, D., Steppe, K., and Ciais, P. (2011). Seasonal patterns of co2 fluxes in amazon 
!! forests: Fusion of eddy covariance data and the orchidee model, Journal of Geophysical Research-Biogeosciences, 116, 
!! Artn G02018, Doi 10.1029/2010jg001544.
!! - MacBean, N., Maignan, F., Peylin, P., Bacour, C., Breon, F. M., & Ciais, P. (2015). Using satellite data to improve the leaf 
!! phenology of a global terrestrial biosphere model. Biogeosciences, 12(23), 7185-7208.
!!
!! SVN          :
!! $HeadURL: $
!! $Date: 2021-05-06 11:05:55 +0200 (四, 2021-05-06) $
!! $Revision: 7177 $
!! \n
!_ ================================================================================================================================

MODULE constantes_mtc

  USE defprec
  USE constantes

  IMPLICIT NONE

  !
  ! METACLASSES CHARACTERISTICS
  !

  INTEGER, PARAMETER :: nvmc = 13                         !! Number of MTCS fixed in the code (unitless)

  CHARACTER(len=34), PARAMETER, DIMENSION(nvmc) :: MTC_name = &  !! description of the MTC (unitless)
  & (/ 'bare ground                       ', &          !  1
  &    'tropical  broad-leaved evergreen  ', &          !  2
  &    'tropical  broad-leaved raingreen  ', &          !  3
  &    'temperate needleleaf   evergreen  ', &          !  4
  &    'temperate broad-leaved evergreen  ', &          !  5
  &    'temperate broad-leaved summergreen', &          !  6
  &    'boreal    needleleaf   evergreen  ', &          !  7
  &    'boreal    broad-leaved summergreen', &          !  8
  &    'boreal    needleleaf   summergreen', &          !  9
  &    '          C3           grass      ', &          ! 10
  &    '          C4           grass      ', &          ! 11
  &    '          C3           agriculture', &          ! 12
  &    '          C4           agriculture'  /)         ! 13


  !
  ! VEGETATION STRUCTURE
  !
  INTEGER,PARAMETER, DIMENSION(nvmc) :: leaf_tab_mtc  =  &                 !! leaf type (1-4, unitless)
  & (/  4,   1,   1,   2,   1,   1,   2,   &                                      !! 1=broad leaved tree, 2=needle leaved tree
  &     1,   2,   3,   3,   3,   3   /)                                           !! 3=grass 4=bare ground

  CHARACTER(len=6), PARAMETER, DIMENSION(nvmc) :: pheno_model_mtc  =  &  !! which phenology model is used? (tabulated) 
  & (/  'none  ',   'none  ',   'moi   ',   'none  ',   'none  ',  &
  &     'ncdgdd',   'none  ',   'ncdgdd',   'ngd   ',   'moigdd',  &
  &     'moigdd',   'moigdd',   'moigdd'  /) 

  LOGICAL, PARAMETER, DIMENSION(nvmc) :: is_tropical_mtc  =  &                       !! Is PFT tropical ? (true/false)
  & (/ .FALSE.,   .TRUE.,    .TRUE.,    .FALSE.,   .FALSE.,   .FALSE.,   .FALSE.,  &
  &    .FALSE.,   .FALSE.,   .FALSE.,   .FALSE.,   .FALSE.,   .FALSE. /)    

!gmjc
  ! Is the vegetation a grassland that can be managed ?
  LOGICAL,PARAMETER, DIMENSION(nvmc) :: is_grassland_manag_mtc   =  &
  & (/ .FALSE., .FALSE.,  .FALSE., .FALSE., .FALSE.,  &
  &    .FALSE.,  .FALSE.,  .FALSE., .FALSE., .FALSE., &
  &    .FALSE., .FALSE., .FALSE. /)
  ! Is the vegetation a grassland that can be cut ?
  LOGICAL,PARAMETER, DIMENSION(nvmc) :: is_grassland_cut_mtc   =  &
  & (/ .FALSE., .FALSE.,  .FALSE., .FALSE., .FALSE.,  &
  &    .FALSE.,  .FALSE.,  .FALSE., .FALSE., .FALSE., &
  &    .FALSE., .FALSE., .FALSE. /)
  ! Is the vegetation a grassland that can be grazed ?
  LOGICAL,PARAMETER, DIMENSION(nvmc) :: is_grassland_grazed_mtc   =  &
  & (/ .FALSE., .FALSE.,  .FALSE., .FALSE., .FALSE.,  &
  &    .FALSE.,  .FALSE.,  .FALSE., .FALSE., .FALSE., &
  &    .FALSE., .FALSE., .FALSE. /)
  ! Management Intensity reading in MANAGEMENT_MAP
  INTEGER, PARAMETER, DIMENSION(nvmc) :: management_intensity_mtc  =  &
  & (/ 0,   0,   0,   0,   0,   0,   0,  &
  &    0,   0,   0,   0,   0,   0 /)
  ! Start year of management reading in MANAGEMENT_MAP & GRAZING_MAP
  INTEGER, PARAMETER, DIMENSION(nvmc) :: management_start_mtc  =  &
  & (/ 0,   0,   0,   0,   0,   0,   0,  &
  &    0,   0,   0,   0,   0,   0 /)
  ! Start year of N deposition reading in DEPOSITION_MAP
  INTEGER, PARAMETER, DIMENSION(nvmc) :: deposition_start_mtc  =  &
  & (/ 0,   0,   0,   0,   0,   0,   0,  &
  &    0,   0,   0,   0,   0,   0 /)
  ! Number of year that management should be read
  INTEGER, PARAMETER, DIMENSION(nvmc) :: nb_year_management_mtc  =  &
  & (/ 1,   1,   1,   1,   1,   1,   1,  &
  &    1,   1,   1,   1,   1,   1 /)
  ! maximum specific leaf area (m**2/gC)
   REAL, PARAMETER, DIMENSION(nvmc)     ::    sla_max_mtc = &  !m2/gC
  & (/ 1.5E-2, 1.53E-2, 2.6E-2, 9.26E-3, 2.0E-2, 2.6E-2, 9.26E-3, &
  &    2.6E-2,  1.9E-2, 4.8E-2,  4.5E-2, 4.8E-2, 4.5E-2 /)
  ! miniimum specific leaf area (m**2/gC)
   REAL, PARAMETER, DIMENSION(nvmc)     ::    sla_min_mtc  = &
  &(/ 1.5E-2, 1.53E-2, 2.6E-2, 9.26E-3,   2E-2, 2.6E-2, 9.26E-3, &
  &   2.6E-2, 1.9E-2,  3.84E-2, 3.6E-2,  3.84E-2, 3.6E-2 /)
!end gmjc
  CHARACTER(LEN=5), PARAMETER, DIMENSION(nvmc) :: type_of_lai_mtc  =  &  !! Type of behaviour of the LAI evolution algorithm
  & (/ 'inter', 'inter', 'inter', 'inter', 'inter',  &                   !! for each vegetation type. (unitless)
  &    'inter', 'inter', 'inter', 'inter', 'inter',  &                   !! Value of type_of_lai : mean or interp
  &    'inter', 'inter', 'inter' /)

  LOGICAL, PARAMETER, DIMENSION(nvmc) :: natural_mtc =  &                         !! natural?  (true/false)
  & (/ .TRUE.,   .TRUE.,   .TRUE.,   .TRUE.,   .TRUE.,    .TRUE.,   .TRUE.,  &
  &    .TRUE.,   .TRUE.,   .TRUE.,   .TRUE.,   .FALSE.,   .FALSE.  /)
! dgvmjc
  LOGICAL, PARAMETER, DIMENSION(nvmc) :: pasture_mtc =  & !! pasture? (true/false)
  & (/ .FALSE.,   .FALSE.,   .FALSE.,   .FALSE.,   .FALSE.,    .FALSE., .FALSE.,  &
  &    .FALSE.,   .FALSE.,   .FALSE.,   .FALSE.,   .FALSE.,   .FALSE. /)
  LOGICAL, PARAMETER, DIMENSION(nvmc) :: rangeland_mtc =  & !! rangeland? (true/false)
  & (/ .FALSE.,   .FALSE.,   .FALSE.,   .FALSE.,   .FALSE.,    .FALSE., .FALSE.,  &
  &    .FALSE.,   .FALSE.,   .FALSE.,   .FALSE.,   .FALSE.,   .FALSE. /)
  LOGICAL, PARAMETER, DIMENSION(nvmc) :: bioenergy_mtc =  & !! bioenergy crop? (true/false)
  & (/ .FALSE.,   .FALSE.,   .FALSE.,   .FALSE.,   .FALSE.,    .FALSE., .FALSE., &
  &    .FALSE.,   .FALSE.,   .FALSE.,   .FALSE.,   .FALSE.,   .FALSE. /)
  LOGICAL, PARAMETER, DIMENSION(nvmc) :: is_agri_fert_mtc =  & !! fertilized? (true/false)
  & (/ .FALSE.,   .FALSE.,   .FALSE.,   .FALSE.,   .FALSE.,    .FALSE., .FALSE., &
  &    .FALSE.,   .FALSE.,   .FALSE.,   .FALSE.,   .FALSE.,   .FALSE. /)
!end dgvmjc

  REAL, PARAMETER, DIMENSION(nvmc) :: veget_ori_fixed_mtc  =  &  !! Value for veget_ori for tests in
  & (/ 0.2,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  &                !! 0-dim simulations (0-1, unitless)
  &    0.0,   0.0,   0.8,   0.0,   0.0,   0.0  /)

  REAL, PARAMETER, DIMENSION(nvmc) :: llaimax_mtc  =  &          !! laimax for maximum
  & (/ 0.0,   8.0,   8.0,   4.0,   4.5,   4.5,   4.0,  &                !! See also type of lai interpolation
  &    4.5,   4.0,   2.0,   2.0,   2.0,   2.0  /)                       !! @tex $(m^2.m^{-2})$ @endtex

  REAL, PARAMETER, DIMENSION(nvmc) :: llaimin_mtc  = &           !! laimin for minimum lai
  & (/ 0.0,   8.0,   0.0,   4.0,   4.5,   0.0,   4.0,  &                !! See also type of lai interpolation (m^2.m^{-2})
  &    0.0,   0.0,   0.0,   0.0,   0.0,   0.0  /)                       !! @tex $(m^2.m^{-2})$ @endtex

  REAL, PARAMETER, DIMENSION(nvmc) :: height_presc_mtc  =  &     !! prescribed height of vegetation (m)
  & (/  0.0,   30.0,   30.0,   20.0,   20.0,   20.0,   15.0,  &         !! Value for height_presc : one for each vegetation type
  &    15.0,   15.0,    0.5,    0.6,    1.0,    1.0  /)

  REAL, PARAMETER, DIMENSION(nvmc) :: z0_over_height_mtc = &         !! Factor to calculate roughness height from 
  & (/  0.0, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625,  &         !! vegetation height (unitless)   
  &  0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625  /)

  REAL, PARAMETER, DIMENSION(nvmc) :: ratio_z0m_z0h_mtc = &      !! Ratio between z0m and z0h values (roughness height for momentum and for heat)
  & (/  1.0,    1.0,    1.0,    1.0,    1.0,    1.0,    1.0,  &         
  &     1.0,    1.0,    1.0,    1.0,    1.0,    1.0  /)


  REAL, PARAMETER, DIMENSION(nvmc) :: rveg_mtc  =  &             !! Potentiometer to set vegetation resistance (unitless)
  & (/ 1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,  &                !! Nathalie on March 28th, 2006 - from Fred Hourdin,
  &    1.0,   1.0,   1.0,   1.0,   1.0,   1.0   /)

  REAL, PARAMETER, DIMENSION(nvmc) :: sla_mtc  =  &                       !! specif leaf area @tex $(m^2.gC^{-1})$ @endtex
  & (/ 1.5E-2,   1.53E-2,   2.6E-2,   9.26E-3,   2.6E-2,   2.6E-2,   9.26E-3,  &
  &    2.6E-2,    1.9E-2,   4.1E-2,    4.0E-2,   4.1E-2,   4.0E-2  /) 

  REAL, PARAMETER, DIMENSION(nvmc) :: availability_fact_mtc  =  &     !! calculate mortality in lpj_gap
  & (/ undef,   0.14,  0.14,   0.10,   0.10,   0.10,   0.05,  &
  &     0.05,   0.05,  undef,  undef,  undef,  undef  /)

  !
  ! EVAPOTRANSPIRATION (sechiba)
  !
  REAL, PARAMETER, DIMENSION(nvmc) :: rstruct_const_mtc  =  &  !! Structural resistance.
  & (/ 0.0,   25.0,   25.0,   25.0,   25.0,   25.0,   25.0,  &        !! @tex $(s.m^{-1})$ @endtex
  &   25.0,   25.0,    2.5,    2.0,    2.0,    2.0   /)

  REAL, PARAMETER, DIMENSION(nvmc) :: kzero_mtc  =  &                  !! A vegetation dependent constant used in the 
  & (/    0.0,   12.E-5,   12.E-5,   12.E-5,   12.E-5,   25.E-5,   12.E-5,  & !! calculation  of the surface resistance. 
  &    25.E-5,   25.E-5,   30.E-5,   30.E-5,   30.E-5,   30.E-5  /)           !! @tex $(kg.m^2.s^{-1})$ @endtex


  !
  ! WATER (sechiba)
  !
  REAL, PARAMETER, DIMENSION(nvmc) :: wmax_veg_mtc  =  &        !! Volumetric available soil water capacity in each PFT
  & (/ 150.0,   150.0,   150.0,   150.0,   150.0,   150.0,   150.0,  & !! @tex $(kg.m^{-3} of soil)$ @endtex
  &    150.0,   150.0,   150.0,   150.0,   150.0,   150.0  /)         
                                                                      

  REAL, PARAMETER, DIMENSION(nvmc) :: humcste_ref4m  =  &       !! Root profile description for the different 
  & (/ 5.0,   0.4,   0.4,   1.0,   0.8,   0.8,   1.0,  &               !! vegetations types. @tex $(m^{-1})$ @endtex
  &    1.0,   0.8,   4.0,   1.0,   4.0,   1.0  /)                      !! These are the factor in the exponential which gets       
                                                                       !! the root density as a function of depth
                                                                       !! Values for zmaxh = 4.0  
  
  REAL, PARAMETER, DIMENSION(nvmc) :: humcste_ref2m  =  &       !! Root profile description for the different
  & (/ 5.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,  &               !! vegetations types.  @tex $(m^{-1})$ @endtex
  &    1.0,   4.0,   0.6,   0.6,   4.0,   4.0  /)                      !! These are the factor in the exponential which gets       
                                                                       !! the root density as a function of depth
                                                                       !! Values for zmaxh = 2.0

  REAL, PARAMETER, DIMENSION(nvmc) :: throughfall_by_mtc  =  &  !! Fraction of rain intercepted by the canopy
  & (/ 30.0,   30.0,   30.0,   30.0,   30.0,   30.0,   30.0,  &        !! (0-100, unitless)
  &    30.0,   30.0,   30.0,   30.0,   30.0,   30.0  /)


  !
  ! ALBEDO (sechiba)
  !
  REAL, PARAMETER, DIMENSION(nvmc) :: snowa_aged_vis_mtc  =  &  !! Minimum snow albedo value for each vegetation type
  & (/ 0.50,    0.0,    0.0,   0.15,   0.14,   0.14,   0.15,  &        !! after aging (dirty old snow) (unitless), visible albedo
  &    0.14,   0.22,   0.35,   0.35,   0.35,   0.35  /)                !! Source : Values are from the Thesis of S. Chalita (1992), optimized on 04/07/2016

  REAL, PARAMETER, DIMENSION(nvmc) :: snowa_aged_nir_mtc  =  &  !! Minimum snow albedo value for each vegetation type
  & (/ 0.35,    0.0,    0.0,   0.14,   0.14,   0.14,   0.14,  &        !! after aging (dirty old snow) (unitless), near infrared albedo
  &    0.14,   0.14,   0.18,   0.18,   0.18,   0.18  /)                !! Source : Values are from the Thesis of S. Chalita (1992)

  REAL, PARAMETER, DIMENSION(nvmc) :: snowa_dec_vis_mtc  =  &   !! Decay rate of snow albedo value for each vegetation type
  & (/ 0.45,    0.0,    0.0,   0.10,   0.06,   0.11,   0.10,  &        !! as it will be used in condveg_snow (unitless), visible albedo
  &    0.11,   0.18,   0.60,   0.60,   0.60,   0.60  /)                !! Source : Values are from the Thesis of S. Chalita (1992), optimized on 04/07/2016

  REAL, PARAMETER, DIMENSION(nvmc) :: snowa_dec_nir_mtc  =  &   !! Decay rate of snow albedo value for each vegetation type
  & (/ 0.45,    0.0,    0.0,   0.06,   0.06,   0.11,   0.06,  &        !! as it will be used in condveg_snow (unitless), near infrared albedo
  &    0.11,   0.11,   0.52,   0.52,   0.52,   0.52  /)                !! Source : Values are from the Thesis of S. Chalita (1992)

  REAL, PARAMETER, DIMENSION(nvmc) :: alb_leaf_vis_mtc  =  &    !! leaf albedo of vegetation type, visible albedo, optimized on 04/07/2016
  & (/ 0.00,   0.0397, 0.0474, 0.0386, 0.0484, 0.0411, 0.041,  &       !! (unitless)
  &    0.0541, 0.0435, 0.0524, 0.0508, 0.0509, 0.0606  /)

  REAL, PARAMETER, DIMENSION(nvmc) :: alb_leaf_nir_mtc  =  &    !! leaf albedo of vegetation type, near infrared albedo, optimized on 04/07/2016
  & (/ 0.00,   0.227,  0.214,  0.193,  0.208,  0.244,  0.177,  &       !! (unitless)
  &    0.218,  0.213,  0.252,  0.265,  0.272,  0.244  /)

  !
  ! SOIL - VEGETATION
  !
  INTEGER, PARAMETER, DIMENSION(nvmc) :: pref_soil_veg_mtc  =  &       !! The soil tile number for each vegetation
  & (/ 1,   2,   2,   2,   2,   2,   2,  &                                    
  &    2,   2,   3,   3,   3,   3  /)                                         


  !
  ! PHOTOSYNTHESIS
  !
  !-
  ! 1 .CO2
  !-
  LOGICAL, PARAMETER, DIMENSION(nvmc) :: is_c4_mtc  =  &                            !! flag for C4 vegetation types (true/false)
  & (/ .FALSE.,  .FALSE.,   .FALSE.,   .FALSE.,   .FALSE.,   .FALSE.,   .FALSE.,  &
  &    .FALSE.,  .FALSE.,   .FALSE.,   .TRUE.,    .FALSE.,   .TRUE.  /)

  REAL, PARAMETER, DIMENSION(nvmc) :: vcmax_fix_mtc  =  &     !! values used for vcmax when STOMATE is not
  & (/  0.0,   40.0,   50.0,   30.0,   35.0,   40.0,   30.0,  &      !! activated @tex $(\mu mol.m^{-2}.s^{-1})$ @endtex
  &    40.0,   35.0,   60.0,   60.0,   70.0,   70.0  /)

! For C4 plant we define a very small downregulation effect as C4 plant are
! currently saturate with respect to CO2 impact on vcmax
  REAL, PARAMETER, DIMENSION(nvmc) :: downregulation_co2_coeff_mtc  =  &  !! coefficient for CO2 downregulation 
  & (/  0.0,   0.38,   0.38,   0.28,   0.28,   0.28,   0.22,  &
  &     0.22,  0.22,   0.26,   0.03,   0.26,   0.03 /)

  REAL, PARAMETER, DIMENSION(nvmc) :: E_KmC_mtc  = &            !! Energy of activation for KmC (J mol-1)
  & (/undef,  79430.,  79430.,  79430.,  79430.,  79430.,  79430.,  &  !! See Medlyn et al. (2002) 
  &  79430.,  79430.,  79430.,  79430.,  79430.,  79430.  /)           !! from Bernacchi al. (2001)

  REAL, PARAMETER, DIMENSION(nvmc) :: E_KmO_mtc  = &            !! Energy of activation for KmO (J mol-1)
  & (/undef,  36380.,  36380.,  36380.,  36380.,  36380.,  36380.,  &  !! See Medlyn et al. (2002) 
  &  36380.,  36380.,  36380.,  36380.,  36380.,  36380.  /)           !! from Bernacchi al. (2001)

REAL, PARAMETER, DIMENSION(nvmc) :: E_Sco_mtc  = &            !! Energy of activation for Sco (J mol-1)
  & (/undef, -24460., -24460., -24460., -24460., -24460., -24460.,  &  !! See Table 2 of Yin et al. (2009)
  & -24460., -24460., -24460., -24460., -24460., -24460.  /)           !! Value for C4 plants is not mentioned - We use C3 for all plants


  REAL, PARAMETER, DIMENSION(nvmc) :: E_gamma_star_mtc  = &     !! Energy of activation for gamma_star (J mol-1)
  & (/undef,  37830.,  37830.,  37830.,  37830.,  37830.,  37830.,  &  !! See Medlyn et al. (2002) from Bernacchi al. (2001) 
  &  37830.,  37830.,  37830.,  37830.,  37830.,  37830.  /)           !! for C3 plants - We use the same values for C4 plants

  REAL, PARAMETER, DIMENSION(nvmc) :: E_Vcmax_mtc  = &          !! Energy of activation for Vcmax (J mol-1)
  & (/undef,  71513.,  71513.,  71513.,  71513.,  71513.,  71513.,  &  !! See Table 2 of Yin et al. (2009) for C4 plants
  &  71513.,  71513.,  71513.,  67300.,  71513.,  67300.  /)           !! and Kattge & Knorr (2007) for C3 plants (table 3)

  REAL, PARAMETER, DIMENSION(nvmc) :: E_Jmax_mtc  = &            !! Energy of activation for Jmax (J mol-1)
  & (/undef,  49884.,  49884.,  49884.,  49884.,  49884.,  49884.,  &   !! See Table 2 of Yin et al. (2009) for C4 plants
  &  49884.,  49884.,  49884.,  77900.,  49884.,  77900.  /)            !! and Kattge & Knorr (2007) for C3 plants (table 3)

  REAL, PARAMETER, DIMENSION(nvmc) :: aSV_mtc     = &            !! a coefficient of the linear regression (a+bT) defining the Entropy term for Vcmax (J K-1 mol-1)
  & (/undef,  668.39,  668.39,  668.39,  668.39,  668.39,  668.39,  &   !! See Table 3 of Kattge & Knorr (2007)
  &  668.39,  668.39,  668.39,  641.64,  668.39,  641.64  /)            !! For C4 plants, we assume that there is no
                                                                        !! acclimation and that at for a temperature of 25°C, aSV is the same for both C4 and C3 plants (no strong jusitification - need further parametrization)

  REAL, PARAMETER, DIMENSION(nvmc) :: bSV_mtc     = &            !! b coefficient of the linear regression (a+bT) defining the Entropy term for Vcmax (J K-1 mol-1 °C-1)
  & (/undef,   -1.07,   -1.07,   -1.07,   -1.07,   -1.07,   -1.07,  &   !! See Table 3 of Kattge & Knorr (2007)
  &   -1.07,   -1.07,   -1.07,      0.,   -1.07,      0.  /)            !! We assume No acclimation term for C4 plants

  REAL, PARAMETER, DIMENSION(nvmc) :: tphoto_min_mtc  =  &  !! minimum photosynthesis temperature (deg C) 
  & (/  undef,   -4.0,    -4.0,   -4.0,   -4.0,   -4.0,   -4.0,  & 
  &      -4.0,   -4.0,    -4.0,   -4.0,   -4.0,   -4.0  /)

  REAL, PARAMETER, DIMENSION(nvmc) :: tphoto_max_mtc  =  &  !! maximum photosynthesis temperature (deg C) 
  & (/  undef,   55.0,    55.0,   55.0,   55.0,   55.0,   55.0,  & 
  &      55.0,   55.0,    55.0,   55.0,   55.0,   55.0  /)

  REAL, PARAMETER, DIMENSION(nvmc) :: aSJ_mtc     = &            !! a coefficient of the linear regression (a+bT) defining the Entropy term for Jmax (J K-1 mol-1)
  & (/undef,  659.70,  659.70,  659.70,  659.70,  659.70,  659.70,  &   !! See Table 3 of Kattge & Knorr (2007)
  &  659.70,  659.70,  659.70,    630.,  659.70,    630.  /)            !! and Table 2 of Yin et al. (2009) for C4 plants

  REAL, PARAMETER, DIMENSION(nvmc) :: bSJ_mtc     = &            !! b coefficient of the linear regression (a+bT) defining the Entropy term for Jmax (J K-1 mol-1 °C-1)
  & (/undef,   -0.75,   -0.75,   -0.75,   -0.75,   -0.75,   -0.75,  &   !! See Table 3 of Kattge & Knorr (2007)
  &   -0.75,   -0.75,   -0.75,      0.,   -0.75,      0.  /)            !! We assume no acclimation term for C4 plants

  REAL, PARAMETER, DIMENSION(nvmc) :: D_Vcmax_mtc  = &           !! Energy of deactivation for Vcmax (J mol-1)
  & (/undef, 200000., 200000., 200000., 200000., 200000., 200000.,  &   !! Medlyn et al. (2002) also uses 200000. for C3 plants (same value than D_Jmax)
  & 200000., 200000., 200000., 192000., 200000., 192000.  /)            !! 'Consequently', we use the value of D_Jmax for C4 plants

  REAL, PARAMETER, DIMENSION(nvmc) :: D_Jmax_mtc  = &            !! Energy of deactivation for Jmax (J mol-1)
  & (/undef, 200000., 200000., 200000., 200000., 200000., 200000.,  &   !! See Table 2 of Yin et al. (2009)
  & 200000., 200000., 200000., 192000., 200000., 192000.  /)            !! Medlyn et al. (2002) also uses 200000. for C3 plants

  REAL, PARAMETER, DIMENSION(nvmc) :: E_gm_mtc  = &              !! Energy of activation for gm (J mol-1) 
  & (/undef,  49600.,  49600.,  49600.,  49600.,  49600.,  49600.,  &   !! See Table 2 of Yin et al. (2009) 
  &  49600.,  49600.,  49600.,   undef,  49600.,   undef  /)            
	 	 
  REAL, PARAMETER, DIMENSION(nvmc) :: S_gm_mtc  = &              !! Entropy term for gm (J K-1 mol-1) 
  & (/undef,   1400.,   1400.,   1400.,   1400.,   1400.,   1400.,  &   !! See Table 2 of Yin et al. (2009) 
  &   1400.,   1400.,   1400.,   undef,   1400.,   undef  /) 
	 	 
  REAL, PARAMETER, DIMENSION(nvmc) :: D_gm_mtc  = &              !! Energy of deactivation for gm (J mol-1) 
  & (/undef, 437400., 437400., 437400., 437400., 437400., 437400.,  &   !! See Table 2 of Yin et al. (2009) 
  & 437400., 437400., 437400.,   undef, 437400.,   undef  /)            

  REAL, PARAMETER, DIMENSION(nvmc) :: E_Rd_mtc  = &              !! Energy of activation for Rd (J mol-1)
  & (/undef,  46390.,  46390.,  46390.,  46390.,  46390.,  46390.,  &   !! See Table 2 of Yin et al. (2009)
  &  46390.,  46390.,  46390.,  46390.,  46390.,  46390.  /)           

  REAL, PARAMETER, DIMENSION(nvmc) :: Vcmax25_mtc  =  &          !! Maximum rate of Rubisco activity-limited carboxylation at 25°C
  & (/ undef,   50.0,    65.0,    45.0,   45.0,   55.0,   45.0,  &      !! @tex $(\mu mol.m^{-2}.s^{-1})$ @endtex
  &     45.0,   35.0,    70.0,    70.0,   70.0,   70.0  /)

  REAL, PARAMETER, DIMENSION(nvmc) :: arJV_mtc    = &            !! a coefficient of the linear regression (a+bT) defining the Jmax25/Vcmax25 ratio (mu mol e- (mu mol CO2)-1)
  & (/undef,    2.00,    2.00,    2.04,    2.22,    2.18,    2.28,  &   !! comes from ORCHIDEE-CAN (no reference given)
  &    2.18,    2.00,    2.00,   2.00,    2.00,   2.00  /)              !! 
                                                                        !! 

  REAL, PARAMETER, DIMENSION(nvmc) :: brJV_mtc    = &            !! b coefficient of the linear regression (a+bT) defining the Jmax25/Vcmax25 ratio ((mu mol e- (mu mol CO2)-1) (°C)-1)
  & (/undef,  -0.035,  -0.035,  -0.035,  -0.035,  -0.035,  -0.035,  &   !! See Table 3 of Kattge & Knorr (2007)
  &  -0.035,  -0.035,  -0.035,      0.,  -0.035,      0.  /)            !! We assume No acclimation term for C4 plants

  REAL, PARAMETER, DIMENSION(nvmc) :: KmC25_mtc  = &             !! Michaelis–Menten constant of Rubisco for CO2 at 25°C (ubar)
  & (/undef,   404.9,   404.9,   404.9,   404.9,  404.9,   404.9,  &    !! See Table 2 of Yin et al. (2009) for C4
  &   404.9,   404.9,   404.9,    650.,   404.9,   650.  /)             !! and Medlyn et al (2002) for C3

  REAL, PARAMETER, DIMENSION(nvmc) :: KmO25_mtc  = &             !! Michaelis–Menten constant of Rubisco for O2 at 25°C (ubar)
  & (/undef, 278400., 278400., 278400., 278400., 278400., 278400.,  &   !! See Table 2 of Yin et al. (2009) for C4 plants and Medlyn et al. (2002) for C3
  & 278400., 278400., 278400., 450000., 278400., 450000.  /)           

REAL, PARAMETER, DIMENSION(nvmc) :: Sco25_mtc  = &             !! Relative CO2 /O2 specificity factor for Rubisco at 25Â°C (bar bar-1)
  & (/undef,   2800.,   2800.,   2800.,   2800.,   2800.,   2800.,  &   !! See Table 2 of Yin et al. (2009)
  &   2800.,   2800.,   2800.,   2590.,   2800.,   2590.  /)           

  REAL, PARAMETER, DIMENSION(nvmc) :: gm25_mtc  = &              !! Mesophyll diffusion conductance at 25Â°C (mol m-2 s-1 bar-1) 
  & (/undef,     0.4,     0.4,     0.4,     0.4,    0.4,      0.4,  &   !! See legend of Figure 6 of Yin et al. (2009) 
  &     0.4,     0.4,     0.4,   undef,     0.4,  undef  /)             !! and review by Flexas et al. (2008) - gm is not used for C4 plants 

  REAL, PARAMETER, DIMENSION(nvmc) :: gamma_star25_mtc  = &      !! Ci-based CO2 compensation point in the absence of Rd at 25°C (ubar)
  & (/undef,   42.75,   42.75,   42.75,   42.75,   42.75,   42.75,  &   !! See Medlyn et al. (2002) for C3 plants - For C4 plants, we use the same value (probably uncorrect)
  &   42.75,   42.75,   42.75,   42.75,   42.75,   42.75  /)    

  REAL, PARAMETER, DIMENSION(nvmc) :: a1_mtc  = &                !! Empirical factor involved in the calculation of fvpd (-)
  & (/undef,    0.85,    0.85,    0.85,    0.85,    0.85,  0.85,  &     !! See Table 2 of Yin et al. (2009)
  &    0.85,    0.85,    0.85,    0.72,    0.85,    0.72  /)           

  REAL, PARAMETER, DIMENSION(nvmc) :: b1_mtc  = &                !! Empirical factor involved in the calculation of fvpd (-)
  & (/undef,    0.14,    0.14,    0.14,    0.14,    0.14,  0.14,  &     !! See Table 2 of Yin et al. (2009)
  &    0.14,    0.14,    0.14,    0.20,    0.14,    0.20  /)           

  REAL, PARAMETER, DIMENSION(nvmc) :: g0_mtc  = &                !! Residual stomatal conductance when irradiance approaches zero (mol CO2 m−2 s−1 bar−1)
  & (/undef, 0.00625, 0.00625, 0.00625, 0.00625, 0.00625, 0.00625,  &   !! Value from ORCHIDEE - No other reference.
  & 0.00625, 0.00625, 0.00625, 0.01875, 0.00625, 0.01875  /)            !! modofy to account for the conversion for conductance to H2O to CO2 

  REAL, PARAMETER, DIMENSION(nvmc) :: h_protons_mtc  = &         !! Number of protons required to produce one ATP (mol mol-1)
  & (/undef,      4.,      4.,      4.,      4.,      4.,    4.,  &     !! See Table 2 of Yin et al. (2009) - h parameter
  &      4.,      4.,      4.,      4.,      4.,      4.  /)           

  REAL, PARAMETER, DIMENSION(nvmc) :: fpsir_mtc = &              !! Fraction of PSII e− transport rate 
  & (/undef,   undef,   undef,   undef,   undef,  undef,  undef,  &     !! partitioned to the C4 cycle (-)
  &   undef,   undef,   undef,     0.4,   undef,    0.4  /)             !! See Table 2 of Yin et al. (2009) - x parameter        
 
  REAL, PARAMETER, DIMENSION(nvmc) :: fQ_mtc = &                 !! Fraction of electrons at reduced plastoquinone 
  & (/undef,   undef,   undef,   undef,   undef,  undef,  undef,  &     !! that follow the Q-cycle (-) - Values for C3 platns are not used
  &   undef,   undef,   undef,      1.,   undef,     1.  /)             !! See Table 2 of Yin et al. (2009)          

  REAL, PARAMETER, DIMENSION(nvmc) :: fpseudo_mtc = &            !! Fraction of electrons at PSI that follow 
  & (/undef,   undef,   undef,   undef,   undef,  undef,  undef,  &     !! pseudocyclic transport (-) - Values for C3 platns are not used
  &   undef,   undef,   undef,     0.1,   undef,    0.1  /)             !! See Table 2 of Yin et al. (2009)    

  REAL, PARAMETER, DIMENSION(nvmc) :: kp_mtc = &                 !! Initial carboxylation efficiency of the PEP carboxylase (mol m−2 s−1 bar−1) 
  & (/undef,   undef,   undef,   undef,   undef,  undef,  undef,  &     !! See Table 2 of Yin et al. (2009)
  &   undef,   undef,   undef,     0.7,   undef,    0.7  /)                 

  REAL, PARAMETER, DIMENSION(nvmc) :: alpha_mtc = &              !! Fraction of PSII activity in the bundle sheath (-)
  & (/undef,   undef,   undef,   undef,   undef,  undef,  undef,  &     !! See legend of Figure 6 of Yin et al. (2009)
  &   undef,   undef,   undef,     0.1,   undef,    0.1  /)                 

  REAL, PARAMETER, DIMENSION(nvmc) :: gbs_mtc = &                !! Bundle-sheath conductance (mol m−2 s−1 bar−1)
  & (/undef,   undef,   undef,   undef,   undef,  undef,  undef,  &     !! See legend of Figure 6 of Yin et al. (2009)
  &   undef,   undef,   undef,   0.003,   undef,  0.003  /)    

  REAL, PARAMETER, DIMENSION(nvmc) :: theta_mtc = &              !! Convexity factor for response of J to irradiance (-)
  & (/undef,     0.7,     0.7,     0.7,     0.7,    0.7,    0.7,  &     !! See Table 2 of Yin et al. (2009) 
  &     0.7,     0.7,     0.7,     0.7,     0.7,    0.7  /)

  REAL, PARAMETER, DIMENSION(nvmc) :: alpha_LL_mtc = &           !! Conversion efficiency of absorbed light into J at strictly limiting light (mol e− (mol photon)−1)
  & (/undef,     0.53,     0.53,     0.3,     0.3,    0.3,    0.3,  &   !! See comment from Yin et al. (2009) after eq. 4
  &     0.3,     0.3,     0.3,     0.3,     0.3,    0.3  /)            !! alpha value from Medlyn et al. (2002)   
                                                                        !! 0.093 mol CO2 fixed per mol absorbed photons
                                                                        !! times 4 mol e- per mol CO2 produced

  REAL, PARAMETER, DIMENSION(nvmc) :: stress_vcmax_mtc = &       !! Water stress on vcmax 
  & (/    1.,     1.,     1.,       1.,      1.,     1.,      1., &     !! IF 1 then stress is disabled on Vcmax
  &      1.,     1.,     1.,       0.,      1.,     0.  /)              !! NVrecommneds to keep stres for C4 plants

  REAL, PARAMETER, DIMENSION(nvmc) :: stress_gs_mtc = &          !! Water stress on gs
  & (/    0.,     0.,     0.,       0.,      0.,     0.,      0., &     !! IF 1 then stress is disabled on stomatal conductance
  &      0.,     0.,     0.,       1.,      0.,     1.  /)

  REAL, PARAMETER, DIMENSION(nvmc) :: stress_gm_mtc = &          !! Water stress on gm
  & (/    0.,     0.,     0.,       0.,      0.,     0.,      0., &     !! IF 1 then stress is disabled on mesophyl conductance
  &      0.,     0.,     0.,       1.,      0.,     1.  /)
    
  !-
  ! 2 .Stomate
  !-
  REAL, PARAMETER, DIMENSION(nvmc) :: ext_coeff_mtc  =  &     !! extinction coefficient of the Monsi&Saeki
  & (/ 0.5,   0.5,   0.5,   0.5,   0.5,   0.5,   0.5,  &             !! relationship (1953) ((m2[ground]) (m-2[leaf])) 
  &    0.5,   0.5,   0.5,   0.5,   0.5,   0.5  /)

  REAL, PARAMETER, DIMENSION(nvmc) :: nue_opt_mtc = &           !! Nitrogen use efficiency of Vcmax 
  & (/ undef,   22.,    22.,    20.,    33.,    33.,    20., &         !! ((mumol[CO2] s-1) (gN[leaf])-1); 
  &      33.,   22.,    45.,    45.,    62.75,    62.75     /)             !! based on the work of Kattge et al. (2009, GCB)

  REAL, PARAMETER, DIMENSION(nvmc) :: ext_coeff_vegetfrac_mtc  =  &     !! extinction coefficient used for defining the fraction
  & (/ 1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,  &                       !!  of bare soil (unitless)
  &    1.0,   1.0,   1.0,   1.0,   1.0,   1.0  /)

!DSG: back from 0.18 to the OCN value of 0.11 as Vcmax is most of the time too high for
!realistic leaf C:N concentrations; now its better
  REAL, PARAMETER, DIMENSION(nvmc) :: ext_coeff_N_mtc  =  &    !! extinction coefficient of the leaf N content profile within the canopy
  & (/ 0.11,  0.11,  0.11,  0.11,  0.11,  0.11,  0.11,  &             !! ((m2[ground]) (m-2[leaf]))
  &    0.11,  0.11,  0.11,  0.11,  0.11,  0.11  /)                    !! based on Dewar et al. (2012, value of 0.18), on Carswell et al. (2000, value of 0.11 used in OCN) 
!DSG
 
  !
  ! ALLOCATION (stomate)
  ! 
  REAL, PARAMETER, DIMENSION(nvmc) :: R0_mtc = &              !! Default root allocation (0-1, unitless)
  & (/ undef,   0.30,   0.30,   0.30,   0.30,  0.30,    0.30, &
  &     0.30,   0.30,   0.30,   0.30,   0.30,  0.30 /)                   

  REAL, PARAMETER, DIMENSION(nvmc) :: S0_mtc = &              !! Default sapwood allocation (0-1, unitless)
  & (/ undef,   0.25,   0.25,   0.30,   0.30,  0.30,    0.30, &
  &     0.30,   0.30,   0.30,   0.30,   0.30,  0.30 /)                   

  !
  ! RESPIRATION (stomate)
  !
  REAL, PARAMETER, DIMENSION(nvmc) :: frac_growthresp_mtc  =  &  !! fraction of GPP which is lost as growth respiration
  & (/  undef,   0.28,   0.28,   0.28,   0.28,   0.28,   0.28,  &
  &      0.28,   0.28,   0.28,   0.28,   0.28,   0.28  /)

!=====================================
!=====================================
!DSG: these parameter seem not be in use (anymore): START
  REAL, PARAMETER, DIMENSION(nvmc) :: maint_resp_slope_c_mtc  =  &  !! slope of maintenance respiration coefficient (1/K),
  & (/  undef,   0.20,   0.20,   0.16,   0.16,   0.16,   0.16,  &          !! constant c of aT^2+bT+c, tabulated
  &      0.16,   0.16,   0.16,   0.12,   0.16,   0.12  /)

  REAL, PARAMETER, DIMENSION(nvmc) :: maint_resp_slope_b_mtc  =  &  !! slope of maintenance respiration coefficient (1/K),
  & (/  undef,   0.0,        0.0,   0.0,        0.0,   0.0,   0.0,  &      !! constant b of aT^2+bT+c, tabulated
  &       0.0,   0.0,   -0.00133,   0.0,   -0.00133,   0.0  /)

  REAL, PARAMETER, DIMENSION(nvmc) :: maint_resp_slope_a_mtc  =  &  !! slope of maintenance respiration coefficient (1/K),
  & (/  undef,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  &                !! constant a of aT^2+bT+c, tabulated
  &       0.0,   0.0,   0.0,   0.0,   0.0,   0.0  /)

  REAL, PARAMETER, DIMENSION(nvmc) :: cm_zero_leaf_mtc  =   &                  !! maintenance respiration coefficient
  & (/   undef,   2.35E-3,   2.62E-3,   1.01E-3,   2.35E-3,   2.62E-3,   1.01E-3,  &  !! at 0 deg C,for leaves, tabulated, 
  &    2.62E-3,   2.05E-3,   2.62E-3,   2.62E-3,   2.62E-3,   2.62E-3  /)             !! @tex $(gC.gC^{-1}.day^{-1})$ @endtex

  REAL, PARAMETER, DIMENSION(nvmc) :: cm_zero_sapabove_mtc =  &                !! maintenance respiration coefficient 
  & (/   undef,   1.19E-4,   1.19E-4,   1.19E-4,   1.19E-4,   1.19E-4,   1.19E-4,  &  !! at 0 deg C, for sapwood above,
  &    1.19E-4,   1.19E-4,   1.19E-4,   1.19E-4,   1.19E-4,   1.19E-4  /)             !! tabulated, @tex $(gC.gC^{-1}.day^{-1})$ @endtex

  REAL, PARAMETER, DIMENSION(nvmc) :: cm_zero_sapbelow_mtc  =  &               !! maintenance respiration coefficient
  & (/   undef,   1.19E-4,   1.19E-4,   1.19E-4,   1.19E-4,   1.19E-4,   1.19E-4,  &  !! at 0 deg C, for sapwood below, 
  &    1.19E-4,   1.19E-4,   1.19E-4,   1.19E-4,   1.19E-4,   1.19E-4  /)             !! tabulated, @tex $(gC.gC^{-1}.day^{-1})$ @endtex 

  REAL, PARAMETER, DIMENSION(nvmc) :: cm_zero_heartabove_mtc  =  &             !! maintenance respiration coefficient
  & (/  undef,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  &                           !! at 0 deg C, for heartwood above,
  &       0.0,   0.0,   0.0,   0.0,   0.0,   0.0  /)                                  !! tabulated, @tex $(gC.gC^{-1}.day^{-1})$ @endtex 

  REAL, PARAMETER, DIMENSION(nvmc) :: cm_zero_heartbelow_mtc  =  &             !! maintenance respiration coefficient
  & (/  undef,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  &                           !! at 0 deg C, for heartwood below, 
  &       0.0,   0.0,   0.0,   0.0,   0.0,   0.0  /)                                  !! tabulated, @tex $(gC.gC^{-1}.day^{-1})$ @endtex 

  REAL, PARAMETER, DIMENSION(nvmc) :: cm_zero_root_mtc  =  &                   !! maintenance respiration coefficient
  & (/   undef,   1.67E-3,   1.67E-3,   1.67E-3,   1.67E-3,   1.67E-3,   1.67E-3,  &  !! at 0 deg C, for roots, tabulated,
  &    1.67E-3,   1.67E-3,   1.67E-3,   1.67E-3,   1.67E-3,   1.67E-3  /)             !! @tex $(gC.gC^{-1}.day^{-1})$ @endtex 

  REAL, PARAMETER, DIMENSION(nvmc) :: cm_zero_fruit_mtc  =  &                  !! maintenance respiration coefficient
  & (/   undef,   1.19E-4,   1.19E-4,   1.19E-4,   1.19E-4,   1.19E-4,   1.19E-4,  &  !! at 0 deg C, for fruits, tabulated,
  &    1.19E-4,   1.19E-4,   1.19E-4,   1.19E-4,   1.19E-4,   1.19E-4  /)             !!  @tex $(gC.gC^{-1}.day^{-1})$ @endtex

  REAL, PARAMETER, DIMENSION(nvmc) :: cm_zero_carbres_mtc  =  &                !! maintenance respiration coefficient
  & (/   undef,   1.19E-4,   1.19E-4,   1.19E-4,   1.19E-4,   1.19E-4,   1.19E-4,  &  !! at 0 deg C, for carbohydrate reserve,
  &    1.19E-4,   1.19E-4,   1.19E-4,   1.19E-4,   1.19E-4,   1.19E-4  /)             !! tabulated, @tex $(gC.gC^{-1}.day^{-1})$ @endtex

!DSG: these are not used:   END
!=====================================
!=====================================
  REAL, PARAMETER, DIMENSION(nvmc) :: coeff_maint_init_mtc  =   &              !! maintenance respiration coefficient
  & (/   undef,   2.50E-2,   2.50E-2,   3.84E-2,   3.84E-2,   3.84E-2,  13.80E-2,  &  !! at 10 deg C - adjusted accounting for T acclimation 
  &   13.80E-2,  13.80E-2,   4.90E-2,   13.80E-2,   3.84E-2,   3.84E-2  /)            !! trop<temp<boreal; grass C4 > C3

  REAL, PARAMETER, DIMENSION(nvmc) :: tref_maint_resp_mtc  =   &              !! maintenance respiration Temperature coefficient (deg C)
  & (/   undef,     10.02,     10.02,     10.02,     10.02,     10.02,     10.02,  &  
  &      10.02,     10.02,     10.02,     10.02,     10.02,     10.02  /)             
  REAL, PARAMETER, DIMENSION(nvmc) :: tmin_maint_resp_mtc  =   &              !! maintenance respiration Temperature coefficient (deg C)
  & (/   undef,    -46.02,    -46.02,    -46.02,    -46.02,    -46.02,    -46.02,  &  
  &     -46.02,    -46.02,    -46.02,    -46.02,    -46.02,    -46.02  /)             
  REAL, PARAMETER, DIMENSION(nvmc) :: e0_maint_resp_mtc  =   &              !! maintenance respiration Temperature coefficient (unitless)
  & (/   undef,   308.56,     308.56,   308.56,     308.56,   308.56,     308.56,  &  
  &     308.56,   308.56,     308.56,   308.56,     308.56,   308.56  /)             

  !
  ! SOM decomposition (stomate)
  !
  REAL, PARAMETER, DIMENSION(nvmc) :: LC_leaf_mtc = &        !! Lignin/C ratio of leaf pool (unitless)
  & (/  undef,   0.18,   0.18,   0.24,   0.18,   0.18,   0.24,  &   !! based on CN from White et al. (2000)       
  &      0.18,   0.24,   0.09,   0.09,   0.09,   0.09  /)

  REAL, PARAMETER, DIMENSION(nvmc) :: LC_sapabove_mtc = &    !! Lignin/C ratio of sapabove pool (unitless)
  & (/  undef,   0.23,   0.23,   0.29,   0.23,   0.23,   0.29,  &   !! based on CN from White et al. (2000)       
  &      0.23,   0.29,   0.09,   0.09,   0.09,   0.09  /)

  REAL, PARAMETER, DIMENSION(nvmc) :: LC_sapbelow_mtc = &    !! Lignin/C ratio of sapbelow pool (unitless)
  & (/  undef,   0.23,   0.23,   0.29,   0.23,   0.23,   0.29,  &   !! based on CN from White et al. (2000)       
  &      0.23,   0.29,   0.09,   0.09,   0.09,   0.09  /)

  REAL, PARAMETER, DIMENSION(nvmc) :: LC_heartabove_mtc = &  !! Lignin/C ratio of heartabove pool (unitless)
  & (/  undef,   0.23,   0.23,   0.29,   0.23,   0.23,   0.29,  &   !! based on CN from White et al. (2000)       
  &      0.23,   0.29,   0.09,   0.09,   0.09,   0.09  /)

  REAL, PARAMETER, DIMENSION(nvmc) :: LC_heartbelow_mtc = &  !! Lignin/C ratio of heartbelow pool (unitless)
  & (/  undef,   0.23,   0.23,   0.29,   0.23,   0.23,   0.29,  &   !! based on CN from White et al. (2000)       
  &      0.23,   0.29,   0.09,   0.09,   0.09,   0.09  /)

  REAL, PARAMETER, DIMENSION(nvmc) :: LC_fruit_mtc = &       !! Lignin/C ratio of fruit pool (unitless)
  & (/  undef,   0.09,   0.09,   0.09,   0.09,   0.09,   0.09,  &   !! based on CN from White et al. (2000)       
  &      0.09,   0.09,   0.09,   0.09,   0.09,   0.09  /)

  REAL, PARAMETER, DIMENSION(nvmc) :: LC_root_mtc = &        !! Lignin/C ratio of root pool (unitless)
  & (/  undef,   0.22,   0.22,   0.22,   0.22,   0.22,   0.22,  &   !! based on CN from White et al. (2000)       
  &      0.22,   0.22,   0.22,   0.22,   0.22,   0.22  /)

  REAL, PARAMETER, DIMENSION(nvmc) :: LC_carbres_mtc = &     !! Lignin/C ratio of carbres pool (unitless)
  & (/  undef,   0.18,   0.18,   0.24,   0.18,   0.18,   0.24,  &   !! based on CN from White et al. (2000)       
  &      0.18,   0.24,   0.09,   0.09,   0.09,   0.09  /)

  REAL, PARAMETER, DIMENSION(nvmc) :: LC_labile_mtc = &      !! Lignin/C ratio of labile pool (unitless)
  & (/  undef,   0.18,   0.18,   0.24,   0.18,   0.18,   0.24,  &   !! based on CN from White et al. (2000)       
  &      0.18,   0.24,   0.09,   0.09,   0.09,   0.09  /)

  REAL, PARAMETER, DIMENSION(nvmc) :: decomp_factor_mtc  =  &  !! Multpliactive factor modifying the standard decomposition factor for each SOM pool
  & (/  undef,     1.,     1.,     1.,     1.,     1.,     1.,  &          
  &        1.,     1.,     1.,     1.,    1.2,    1.4  /)
  !
  ! STAND STRUCTURE
  !
  REAL, PARAMETER, DIMENSION(nvmc) :: tree_ff_mtc = &                          !! Tree form factor to reduce
  &(/ undef,  0.7,    0.7,   0.7,   0.7,   0.7,   0.75, &                             !! the volume of a cylinder
  &     0.7,  0.7,  undef, undef, undef, undef /)                                     !! to the volume of the real tree shape

  REAL, PARAMETER, DIMENSION(nvmc) :: pipe_density_mtc = &                     !! Wood density @tex $(gC.m^{-3})$ @endtex
  &(/   undef,   3.e5,   3.e5,   2.083e5,   4.8e5,   2.38e5,   1.95e5,  &               !!
  &    2.38e5,  2.4875e5,    2.e5,   2.e5,   2.e5,   2.e5  /)                               !!
                                                                                      !!

  REAL, PARAMETER, DIMENSION(nvmc) :: pipe_tune1_mtc = &                       !! cn_area = pipe_tune1*...
  &(/ undef,  200.,  200.,  83.,  113.,  91.,  55., &                              !!    stem diameter**pipe_tune_exp_coeff
  &    77.,  83., undef, undef, undef, undef /) 
       
  REAL, PARAMETER, DIMENSION(nvmc) :: pipe_tune2_mtc = &                       !! height=pipe_tune2 * diameter**pipe_tune3
  &(/ undef,   40.,   40.,   32.,   14.,   35.,   43., &
  &     41.,   30., undef, undef, undef, undef /) 

  REAL, PARAMETER, DIMENSION(nvmc) :: pipe_tune3_mtc = &                       !! height=pipe_tune2 * diameter**pipe_tune3
  &(/ undef,   0.65,   0.65,   0.57,   0.33,   0.66,  0.58, &
  &     0.66,   0.52, undef, undef, undef, undef /)   
      
  REAL, PARAMETER, DIMENSION(nvmc) :: pipe_tune4_mtc = &                       !! CHECK - needed for stem diameter
  &(/ undef,   0.3,   0.3,   0.3,   0.3,   0.3,   0.3, &
        0.3,   0.3, undef, undef, undef, undef /)

  REAL, PARAMETER, DIMENSION(nvmc) :: pipe_k1_mtc = &                          !! CHECK 
  &(/ undef,  8.e3,  8.e3,  8.e3,  8.e3,  8.e3,  8.e3, &
  &    8.e3,  8.e3, undef, undef, undef, undef /) 

  REAL, PARAMETER, DIMENSION(nvmc) :: pipe_tune_exp_coeff_mtc = &              !! cn_area = pipe_tune1*...  
  &(/ undef,   1.6,   1.6,   1.33,   .85,   1.16,   1.3, &                              !!    stem diameter**pipe_tune_exp_coeff
  &     1.21,   1.4, undef, undef, undef, undef /)  
      
  REAL, PARAMETER, DIMENSION(nvmc) :: mass_ratio_heart_sap_mtc = &             !! mass ratio (heartwood+sapwood)/heartwood 
  &(/ undef,    3.,    3.,    3.,    3.,    3.,    3., &
  &      3.,    3., undef, undef, undef, undef /)

  REAL, PARAMETER, DIMENSION(nvmc) :: canopy_cover_mtc = &                     !! Prescribed canopy cover (1-gap fraction) 
  & (/ undef,    0.9,   0.9,   0.7,   0.7,   0.7,   0.6, &                            !! of a canopy (unitless)
  &      0.5,    0.5,   0.9,   0.9,   0.9,   0.9 /)  

  INTEGER, PARAMETER, DIMENSION(nvmc) :: nmaxtrees_mtc = &                     !! Initial number of trees per ha. This parameter is 
  & (/  -9999,   2000,   2000,   7333,   8000,  15000, 15000,  &                      !! used at .firstcall. and after clearcuts
  &     15000,  15000,  10000,  10000,  10000,  10000 /)                              !! the value is used by the allometric allocation
                                                                                      !! and forestry subroutines.
                                                                                      !! Values from DOFOCO run.def

  REAL, PARAMETER, DIMENSION(nvmc) :: height_init_min_mtc = &                  !! The minimum height (m) of a tree sapling when a forest 
  &(/ undef,    2.,    2.,    2.,    2.,    3.,    2., &                              !! stand is established. Owing to the allometric 
  &      2.,    4.,   0.3,   0.3,   0.3,   0.3 /)                                     !! relationship this setting determines all
                                                                                      !! biomass components of a newly establised stand

  REAL, PARAMETER, DIMENSION(nvmc) :: height_init_max_mtc = &                  !! The maximum height (m) of a tree sapling when a forest 
  &(/ undef,   15.,   15.,    3.,    3.,    4.,    3., &                              !! stand is established. 
  &      3.,    5.,   0.3,   0.3,   0.4,   0.6 /)                                     
                                                                                      

  !
  ! FIRE (stomate)
  !
  REAL,PARAMETER, DIMENSION(nvmc) :: flam_mtc  =  &         !! flamability: critical fraction of water 
  & (/  undef,   0.15,   0.25,   0.25,   0.25,   0.25,   0.25,  &  !! holding capacity (0-1, unitless)
  &      0.25,   0.25,   0.25,   0.25,   0.35,   0.35  /)

  REAL,PARAMETER, DIMENSION(nvmc) :: resist_mtc  =  &       !! fire resistance (0-1, unitless)
  & (/ undef,   0.95,   0.90,   0.90,   0.90,   0.90,   0.90,  &
  &    0.90,    0.90,    0.0,    0.0,    0.0,    0.0 /) 

  REAL,PARAMETER, DIMENSION(nvmc) :: frac_co_fire_mtc  =  & !! CO fraction in fire C emission 
  & (/ undef,   0.08105155,   0.08105155,   0.07705976,   0.07705976, 0.07705976,   0.11705334,  & 
  &    0.11705334,    0.11705334,    0.05529691,    0.05529691, 0.05529691,    0.05529691 /) 

  REAL,PARAMETER, DIMENSION(nvmc) :: frac_ch4_fire_mtc  =  & !! CH4 fraction in fire C emission 
  & (/ undef,   0.007732579,   0.007732579,   0.005148993,  0.005148993, 0.005148993,   0.009613121,  & 
  &    0.009613121,    0.009613121,   0.002979889,   0.002979889, 0.009087086,   0.009087086 /) 
      
  REAL,PARAMETER, DIMENSION(nvmc) :: frac_Nemit_fire_mtc  =  & !! N:C ratio in fire emission  
  & (/ undef,   0.006652489,   0.006652489,   0.004542138, 0.004542138,   0.004542138,   0.008614297,  & 
  &    0.008614297,    0.008614297,   0.005588204,   0.005588204, 0.007446571,   0.007446571 /) 

  REAL,PARAMETER, DIMENSION(nvmc) :: frac_Pemit_fire_mtc  =  & !! P:C ratio in fire emission 
  & (/ undef,   3.689501e-06,   3.689501e-06,   0.0,   0.0,   0.0, 0.0,  & 
  &    0.0,    0.0,   3.007377e-06,   3.007377e-06,   0.0,   0.0 /) 

  REAL,PARAMETER, DIMENSION(nvmc) :: frac_nh3_fire_mtc  =  &  !! NH3 emission fraction in fire emission 
  & (/ undef,   0.3348122,   0.3348122,   0.3111859,   0.3111859, 0.3111859,   0.5592236,  & 
  &    0.5592236,    0.5592236,   0.1569450,   0.1569450, 0.4996005,   0.4996005 /) 

  REAL,PARAMETER, DIMENSION(nvmc) :: frac_nox_fire_mtc  =  &  !!  NOx emission fraction in fire emission 
  & (/ undef,   0.3637622,   0.3637622,   0.4030598,   0.4030598, 0.4030598,   0.1048544,  & 
  &    0.1048544,    0.1048544,   0.6670163,   0.6670163, 0.4057432,   0.4057432 /) 

  REAL,PARAMETER, DIMENSION(nvmc) :: frac_n2o_fire_mtc  =  &  !! N2O emission fraction in fire emission 
  & (/ undef,   0.03890504,   0.03890504,   0.04580225,   0.04580225, 0.04580225,   0.06513684,  & 
  &    0.06513684,    0.06513684,   0.04664450,   0.04664450, 0.01779055,   0.01779055 /) 

  !
  ! FLUX - LUC
  !
  REAL, PARAMETER, DIMENSION(nvmc) :: coeff_lcchange_1_mtc  =  &   !! Coeff of biomass export for the year
  & (/  undef,   0.897,   0.897,   0.597,   0.597,   0.597,   0.597,  &   !! (0-1, unitless)
  &     0.597,   0.597,   0.597,   0.597,   0.597,   0.597  /)

  REAL, PARAMETER, DIMENSION(nvmc) :: coeff_lcchange_10_mtc  =  &  !! Coeff of biomass export for the decade 
  & (/  undef,   0.103,   0.103,   0.299,   0.299,   0.299,   0.299,  &   !! (0-1, unitless)
  &     0.299,   0.299,   0.299,   0.403,   0.299,   0.403  /) 

  REAL, PARAMETER, DIMENSION(nvmc) :: coeff_lcchange_100_mtc  =  & !! Coeff of biomass export for the century
  & (/  undef,     0.0,     0.0,   0.104,   0.104,   0.104,   0.104,  &   !! (0-1, unitless)
  &     0.104,   0.104,   0.104,     0.0,   0.104,     0.0  /)


  !
  ! PHENOLOGY
  !
  ! The latest modifications regarding leafagecrit, senescence_temp_c, leaffall, hum_min_time and nosenescence_hum are inspired by
  ! MacBean et al. (2015), following the optimization of phenology parameters using MODIS NDVI (FM/PP).
  !-
  ! 1. Stomate
  !-
  REAL, PARAMETER, DIMENSION(nvmc) :: lai_max_to_happy_mtc  =  &  !! threshold of LAI below which plant uses carbohydrate reserves
  & (/  undef,   0.5,   0.5,   0.5,   0.5,   0.5,   0.5,  &
  &       0.5,   0.5,   0.5,   0.5,   0.5,   0.5  /)

  REAL, PARAMETER, DIMENSION (nvmc) :: lai_max_mtc  =  &          !! maximum LAI, PFT-specific 
  & (/ undef,   7.0,   7.0,   5.0,   5.0,   5.0,   4.5,  &               !! @tex $(m^2.m^{-2})$ @endtex
  &      4.5,   3.0,   2.5,   2.5,   5.0,   5.0  /)

  INTEGER, PARAMETER, DIMENSION(nvmc) :: pheno_type_mtc  =  &     !! type of phenology (0-4, unitless)
  & (/  0,   1,   3,   1,   1,   2,   1,  &                              !! 0=bare ground 1=evergreen,  2=summergreen, 
  &     2,   2,   4,   4,   2,   3  /)                                   !! 3=raingreen,  4=perennial
  REAL, PARAMETER, DIMENSION (nvmc) :: frac_bioenergy_harvest_mtc  =  & !!for bioenergy PFTs, fraction of biomass to harvest, PFT-specific
  & (/ undef,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  &                   !!@tex $(m^2.m^{-2})$ @endtex
  &      0.0,   0.0,   0.0,   0.0,   0.0, 0.0  /)
  !-
  ! 2. Leaf Onset
  !-
  REAL, PARAMETER, DIMENSION(nvmc) :: pheno_gdd_crit_c_mtc  =  &    !! critical gdd, tabulated (C),
  & (/  undef,   undef,   undef,   undef,   undef,   undef,   undef,  &    !! constant c of aT^2+bT+c
  &     undef,   undef,   320.0,   400.0,   300.0,   400.0  /)

  REAL, PARAMETER, DIMENSION(nvmc) :: pheno_gdd_crit_b_mtc  =  &    !! critical gdd, tabulated (C),
  & (/  undef,   undef,   undef,   undef,   undef,   undef,   undef,  &    !! constant b of aT^2+bT+c
  &     undef,   undef,    6.25,     0.0,    6.25,     0.0  /)

  REAL, PARAMETER, DIMENSION(nvmc) :: pheno_gdd_crit_a_mtc  =  &    !! critical gdd, tabulated (C),
  & (/  undef,   undef,     undef,   undef,   undef,   undef,   undef,  &  !! constant a of aT^2+bT+c
  &     undef,   undef,   0.03125,     0.0,  0.0315,   0.0  /)

  REAL, PARAMETER, DIMENSION(nvmc) :: pheno_moigdd_t_crit_mtc  = &  !! temperature threshold for C4 grass(C)
  & (/  undef,   undef,     undef,   undef,   undef,   undef,   undef,  &  
  &     undef,   undef,     undef,    22.0,   undef,   undef  /)

  REAL, PARAMETER, DIMENSION(nvmc) :: ngd_crit_mtc  =  &            !! critical ngd, tabulated. 
  & (/  undef,   undef,   undef,   undef,   undef,   undef,   undef,  &    !! Threshold -5 degrees (days)
  &     undef,    17.0,   undef,   undef,   undef,   undef  /)

  REAL, PARAMETER, DIMENSION(nvmc) :: ncdgdd_temp_mtc  =  &         !! critical temperature for the ncd vs. gdd 
  & (/  undef,   undef,   undef,   undef,   undef,     5.0,   undef,  &    !! function in phenology (C)
  &       0.0,   undef,   undef,   undef,   undef,   undef  /)

  REAL, PARAMETER, DIMENSION(nvmc) :: hum_frac_mtc  =  &            !! critical humidity (relative to min/max) 
  & (/  undef,   undef,   0.5,   undef,   undef,   undef,   undef, &       !! for phenology (unitless)
  &     undef,   undef,   0.5,     0.5,     0.5,     0.5  /)

  REAL, PARAMETER, DIMENSION(nvmc) :: hum_min_time_mtc  =  &        !! minimum time elapsed since
  & (/  undef,   undef,   50.0,   undef,   undef,   undef,   undef,  &     !! moisture minimum (days)
  &     undef,   undef,   36.0,    35.0,    75.0,    75.0  /) 

!  REAL, PARAMETER, DIMENSION(nvmc) :: tau_root_mtc  =  &            !! roots longevity (days)
!  & (/  undef,  365.,   365.,    365.,    365.,   365.,   365.,  &
!  &      365.,  365.,   365.,    365.,    365.,   365.  /)/7               !!DSG: there is a "/7" here! The values I use in run.def from SL are factor!of 7 larger?
  REAL, PARAMETER, DIMENSION(nvmc) :: tau_root_mtc  =  &            !! roots longevity (days)
  & (/  undef,  730.,   730.,    730.,    730.,   730.,   730.,  &
  &      730.,  730.,   70.,    70.,    70.,   70.  /)
!DSG: recalibration
!  REAL, PARAMETER, DIMENSION(nvmc) :: tau_sap_mtc  =  &             !! time (days)  
!  & (/  undef,   730.0,   730.0,   730.0,   730.0,   730.0,   730.0,  &
!  &     730.0,   730.0,   undef,   undef,   undef,   undef  /)
  REAL, PARAMETER, DIMENSION(nvmc) :: tau_sap_mtc  =  &             !! time (days)  
  & (/  undef,  20000.0,  20000.0,   15148.0,   11680.0,   8760.0,   11680.0,  &
  &     11680.0,   13870.0,   70.,   280.,   280.,   280.  /)
!DSG: recalibration

  REAL, PARAMETER, DIMENSION(nvmc) :: tau_leafinit_mtc  =  &  !! time to attain the initial foliage using the carbohydrate reserve
  & (/  undef,   10.,   10.,   10.,   10.,   10.,   10.,  &
  &       10.,   10.,   10.,   10.,   10.,   10.  /)

  REAL, PARAMETER, DIMENSION(nvmc) :: tau_fruit_mtc  =  &           !! fruit lifetime (days)
  & (/  undef,  90.0,    90.0,    90.0,    90.0,   90.0,   90.0,  &
  &      90.0,  90.0,   undef,   undef,   undef,   undef  /)

  REAL, PARAMETER, DIMENSION(nvmc) :: ecureuil_mtc  =  &            !! fraction of primary leaf and root allocation
  & (/  undef,   0.0,   1.0,   0.0,   0.0,   1.0,   0.0,  &                !! put into reserve (0-1, unitless)
  &       1.0,   1.0,   1.0,   1.0,   1.0,   1.0  /)

  REAL, PARAMETER, DIMENSION(nvmc) :: alloc_min_mtc  =  &           !! NEW - allocation above/below = f(age) 
  & (/  undef,   0.2,     0.2,     0.2,     0.2,    0.2,   0.35,  &        !! O-CAN
  &       0.7,   0.2,     0.2,     0.2,     0.2,    0.2   /)

  REAL, PARAMETER, DIMENSION(nvmc) :: alloc_max_mtc  =  &           !! NEW - allocation above/below = f(age) 
  & (/  undef,   0.8,     0.8,     0.8,     0.8,    0.8,   0.8,  &         !! O-CAN
  &       0.8,   0.8,     0.8,     0.8,     0.8,    0.8  /)

  REAL, PARAMETER, DIMENSION(nvmc) :: demi_alloc_mtc  =  &          !! NEW - allocation above/below = f(age) 
  & (/  undef,   5.0,     5.0,     5.0,     5.0,    5.0,   5.0,  &         !! - 30/01/04 NV/JO/PF
  &       5.0,   5.0,   undef,   undef,   undef,   undef  /)

  REAL, PARAMETER, DIMENSION(nvmc) :: leaflife_mtc  =  &            !! leaf longevity, tabulated (??units??)
  & (/  undef,   0.5,   2.0,   0.33,   1.0,   2.0,   0.33,  &
  &       2.0,   2.0,   2.0,   2.0,    2.0,   2.0  /)
  !-
  ! 3. Senescence
  !-
  REAL, PARAMETER, DIMENSION(nvmc) :: leaffall_mtc  =  &             !! length of death of leaves, tabulated (days)
  & (/  undef,   undef,   5.0,   undef,   undef,   30.0,   undef,  &
  &       5.0,    5.0,   10.0,    10.0,    10.0,   10.0  /)

  REAL, PARAMETER, DIMENSION(nvmc) :: leafagecrit_mtc  =  &          !! critical leaf age, tabulated (days)
  & (/  undef,   730.0,   180.0,   910.0,   730.0,   160.0,   910.0,  &
  &     220.0,   120.0,    80.0,   80.0,   120.0,   120.0  /)

  REAL, PARAMETER, DIMENSION(nvmc) :: rootagecrit_mtc  =  &          !! critical leaf age, tabulated (days)
  & (/  undef,   730.0,   180.0,   910.0,   730.0,   160.0,   910.0,  &
  &     220.0,   120.0,   160.0,   160.0,   120.0,   120.0  /)

  CHARACTER(LEN=6), PARAMETER, DIMENSION(nvmc) :: senescence_type_mtc  =  & !! type of senescence, tabulated (unitless)
  & (/  'none  ',  'none  ',   'dry   ',  'none  ',  'none  ',  &
  &     'cold  ',  'none  ',   'cold  ',  'cold  ',  'mixed ',  &
  &     'mixed ',  'crop  ',   'crop  '            /)

  REAL, PARAMETER, DIMENSION(nvmc) :: senescence_hum_mtc  =  &       !! critical relative moisture availability
  & (/  undef,   undef,   0.3,   undef,   undef,   undef,   undef,  &       !! for senescence (0-1, unitless)
  &     undef,   undef,   0.2,     0.2,     0.3,     0.2  /)

  REAL, PARAMETER, DIMENSION(nvmc) :: nosenescence_hum_mtc  =  &     !! relative moisture availability above which 
  & (/  undef,   undef,   0.8,   undef,   undef,   undef,   undef,  &       !! there is no humidity-related senescence
  &     undef,   undef,   0.6,     0.3,     0.3,     0.3  /)                !! (0-1, unitless)

  REAL, PARAMETER, DIMENSION(nvmc) :: max_turnover_time_mtc  =  &    !! maximum turnover time for grasses (days)
  & (/  undef,   undef,   undef,   undef,   undef,   undef,   undef,  &
  &     undef,   undef,    80.0,    80.0,    80.0,    80.0  /)

  REAL, PARAMETER, DIMENSION(nvmc) :: min_turnover_time_mtc  =  &    !! minimum turnover time for grasses (days)
  & (/  undef,   undef,   undef,   undef,   undef,   undef,   undef,  &
  &     undef,   undef,    10.0,    10.0,    10.0,    10.0  /)
 
  REAL, PARAMETER, DIMENSION(nvmc) :: recycle_leaf_mtc = &          !! Fraction of N leaf that is recycled when leaves are senescent
  & (/  undef,  0.561,   0.612,   0.627,   0.561,  0.612,  0.627,  &       !! Vergutz et al (2012) Fig.4
  &     0.612,  0.627,   0.746,   0.746,   0.746,  0.746  /)

  REAL, PARAMETER, DIMENSION(nvmc) :: recycle_root_mtc = &        !! Fraction of N leaf that is recycled when roots are senescent
  & (/  undef,  0.275, 0.275,  0.275,  0.275,  0.275,  0.275,   &        !! Kunkle et al (2009) , Freschet et al. 2010 
  &     0.275,  0.275, 0.275,  0.275,  0.275,  0.275  /)                 !! both report 27-28% for temperate andsubartic species

  REAL, PARAMETER, DIMENSION(nvmc) :: min_leaf_age_for_senescence_mtc  =  &  !! minimum leaf age to allow 
  & (/  undef,   undef,   90.0,   undef,   undef,   90.0,   undef,  &               !! senescence g (days)
  &      60.0,    60.0,   30.0,    30.0,   120.0,  120.0  /)

  REAL, PARAMETER, DIMENSION(nvmc) :: senescence_temp_c_mtc  =  &    !! critical temperature for senescence (C)
  & (/  undef,   undef,    undef,   undef,   undef,   16.0,   undef,  &     !! constant c of aT^2+bT+c, tabulated
  &      14.0,    10.0,      5.0,     5.0,     5.0,    10.0  /)             !! (unitless)

  REAL, PARAMETER, DIMENSION(nvmc) :: senescence_temp_b_mtc  =  &    !! critical temperature for senescence (C), 
  & (/  undef,   undef,   undef,   undef,   undef,   0.0,   undef,  &       !! constant b of aT^2+bT+c, tabulated
  &       0.0,     0.0,     0.1,     0.0,     0.0,   0.0  /)                !! (unitless)

  REAL, PARAMETER, DIMENSION(nvmc) :: senescence_temp_a_mtc  =  &    !! critical temperature for senescence (C), 
  & (/  undef,   undef,     undef,   undef,   undef,   0.0,   undef,  &     !! constant a of aT^2+bT+c, tabulated
  &       0.0,     0.0,   0.00375,     0.0,     0.0,   0.0  /)              !! (unitless)

  REAL, PARAMETER, DIMENSION(nvmc) :: gdd_senescence_mtc  =  &       !! minimum gdd to allow senescence of crops (days)
  & (/  undef,   undef,    undef,   undef,     undef,    undef,    undef,  &
  &     undef,   undef,    undef,   undef,     1800.,    2700.  /)

  !-
  ! 4. N cycle
  !-
  REAL, PARAMETER, DIMENSION(nvmc) :: cn_leaf_min_mtc  = &                !! minimum CN ratio of leaves 
  & (/  undef,     12.5,      12.5,     28.,       16.,      16.,      28., &    !! (gC/gN)
  &       16.,     28.,      16.,     22.,       16.,      22.   /)              !! 

  REAL, PARAMETER, DIMENSION(nvmc) :: cn_leaf_max_mtc  = &            !! maximum CN ratio of leaves 
 & (/   undef,     45.,      45.,     75.,       45.,      45.,      75.,  & !! (gC/gN) 
 &        45.,     75.,      45.,     60.,       45.,      60.   /)

  REAL, PARAMETER, DIMENSION(nvmc) :: max_soil_n_bnf_mtc   = &        !! Value of total N (NH4+NO3) 
 & (/     0.0,     1.5,      1.5,     1.5,       1.5,      1.5,      1.5,  & 
 &        1.5,     1.5,       2.,      2.,        2.,       2.   /)          !! above which we stop adding N via BNF
                                                                             !! (gN/m**2)
  REAL, PARAMETER, DIMENSION(nvmc) :: fixed_fsoiln_mtc   = &              !! fixed scalar to vary soil C:N range
 & (/     1.0,     .75,      .75,     .5,       .5,      .5,      0.0,  & !!  1=max n content, 0 min n content
         0.0,     0.0,      1.0,     1.0,       1.0,      1.0   /)          !! 

  !-
  ! 5. P cycle
  !- 
  REAL, PARAMETER, DIMENSION(nvmc) :: np_leaf_min_mtc  =           &  !! minimum NP ratio of leaves 
  & (/  undef,     12.,      12.,     12.,      12.,     12.,    12.,   &          !! (gN/gP); 
  &    12.,    12.,     12.,    12.,      12.,     12.   /)                        !! 

  REAL, PARAMETER, DIMENSION(nvmc) :: np_leaf_max_mtc  =                      & !! maximum NP ratio of leaves 
 & (/   undef,     30.00,      30.00,     30.00,       30.00,      30.00,    30.00,  & !! (gN/gP) 
 &      30.00,     30.00,      30.00,     30.00,       30.00,      30.00   /)          !! 

 REAL, PARAMETER, DIMENSION(nvmc) :: p_recycle_leaf_mtc  = &             !! Leaf P resoprtion efficiency
 & (/   undef,     .65,      .65,     .65,       .65,      .65,      .65,  &    !! Vergutz et al (2012) ( intially we used PFT values from Fig.4)
 &        .65,     .65,      .821,     .821,       .821,      .821   /)         !! This lead to sever overestimation of conifer GPP, thus we now take the !global average for tree PFTs

 REAL, PARAMETER, DIMENSION(nvmc) :: p_recycle_root_mtc  = &                  !! Freschet et al 2010
 & (/   undef,  .57,   .57,  .57,   .57,  .57,   .57,  &         !! 
 &        .57,  .57,   .57,  .57,   .57,  .57   /)


  !
  ! DGVM
  !
!DSG
  REAL, PARAMETER, DIMENSION(nvmc) :: residence_time_mtc  =  &    !! residence time of trees (years)
  & (/  undef,   50.0,   50.0,   50.0,   50.0,  50.0,   50.0,  &         !!
  &      50.0,   50.0,    0.0,    0.0,    0.0,    0.0  /)                !!
!DSG: recalibration

  REAL, PARAMETER, DIMENSION(nvmc) :: tmin_crit_mtc  =  &
  & (/  undef,     0.0,     0.0,   -30.0,   -14.0,   -30.0,   -45.0,  &  !! critical tmin, tabulated (C)
  &     -45.0,   -60.0,   undef,   undef,   undef,   undef  /)

  REAL, PARAMETER, DIMENSION(nvmc) :: tcm_crit_mtc  =  &
  & (/  undef,   undef,   undef,     5.0,    15.5,    15.5,   -8.0,  &   !! critical tcm, tabulated (C)
  &      -8.0,    -8.0,   undef,   undef,   undef,   undef  /)



  !
  ! Biogenic Volatile Organic Compounds
  !
  REAL, PARAMETER, DIMENSION(nvmc) :: em_factor_isoprene_mtc = &     !! Isoprene emission factor 
  & (/  0.,    24.,   24.,    8.,   16.,   45.,   8.,  &                    !!
  &    18.,    0.5,   12.,   18.,    5.,    5.  /)

  REAL, PARAMETER, DIMENSION(nvmc) :: em_factor_monoterpene_mtc = &  !! Monoterpene emission factor
  & (/   0.,   2.0,    2.0,   1.8,    1.4,    1.6,    1.8,  &               !! @tex $(\mu gC.g^{-1}.h^{-1})$ @endtex
  &    1.4,    1.8,    0.8,   0.8,    0.22,     0.22  /)

  REAL, PARAMETER :: LDF_mono_mtc = 0.6                                  !! monoterpenes fraction dependancy to light
  REAL, PARAMETER :: LDF_sesq_mtc = 0.5                                  !! sesquiterpenes fraction dependancy to light
  REAL, PARAMETER :: LDF_meth_mtc = 0.8                                  !! methanol fraction dependancy to light
  REAL, PARAMETER :: LDF_acet_mtc = 0.2                                  !! acetone fraction dependancy to light

  REAL, PARAMETER, DIMENSION(nvmc) :: em_factor_apinene_mtc = &      !! Alfa pinene emission factor percentage
  & (/   0.,   0.395,   0.395,   0.354,   0.463,   0.326,   0.354, &        !! ATTENTION: for each PFT they are PERCENTAGE of monoterpene EF
  &   0.316,   0.662,   0.231,   0.200,   0.277,   0.277 /)

  REAL, PARAMETER, DIMENSION(nvmc) :: em_factor_bpinene_mtc = &      !! Beta pinene emission factor  percentage      
  & (/   0.,   0.110,   0.110,   0.146,   0.122,   0.087,   0.146, &        !! ATTENTION: for each PFT they are PERCENTAGE of monoterpene EF
  &   0.063,   0.150,   0.123,   0.080,   0.154,   0.154  /)

  REAL, PARAMETER, DIMENSION(nvmc) :: em_factor_limonene_mtc = &     !! Limonene emission factor percentage
  & (/   0.,   0.092,   0.092,   0.083,   0.122,   0.061,   0.083, &        !! ATTENTION: for each PFT they are PERCENTAGE of monoterpene EF
  &   0.071,   0.037,   0.146,   0.280,   0.092,   0.092  /)

  REAL, PARAMETER, DIMENSION(nvmc) :: em_factor_myrcene_mtc = &      !! Myrcene emission factor percentage
  & (/   0.,   0.073,   0.073,   0.050,   0.054,   0.028,   0.050, &        !! ATTENTION: for each PFT they are PERCENTAGE of monoterpene EF
  &   0.019,   0.025,   0.062,   0.057,   0.046,   0.046  /)

  REAL, PARAMETER, DIMENSION(nvmc) :: em_factor_sabinene_mtc = &     !! Sabinene emission factor percentage
  & (/   0.,   0.073,   0.073,   0.050,   0.083,   0.304,   0.050, &        !! ATTENTION: for each PFT they are PERCENTAGE of monoterpene EF
  &   0.263,   0.030,   0.065,   0.050,   0.062,   0.062  /)

  REAL, PARAMETER, DIMENSION(nvmc) :: em_factor_camphene_mtc = &     !! Camphene emission factor percentage
  & (/   0.,   0.055,   0.055,   0.042,   0.049,   0.004,   0.042, &        !! ATTENTION: for each PFT they are PERCENTAGE of monoterpene EF
  &   0.005,   0.023,   0.054,   0.053,   0.031,   0.031  /)

  REAL, PARAMETER, DIMENSION(nvmc) :: em_factor_3carene_mtc = &      !! 3-carene emission factor percentage
  & (/   0.,   0.048,   0.048,   0.175,   0.010,   0.024,   0.175, &        !! ATTENTION: for each PFT they are PERCENTAGE of monoterpene EF
  &   0.013,   0.042,   0.065,   0.057,   0.200,   0.200  /)

  REAL, PARAMETER, DIMENSION(nvmc) :: em_factor_tbocimene_mtc = &    !! T-beta-ocimene emission factor percentage
  & (/   0.,   0.092,   0.092,   0.054,   0.044,   0.113,   0.054, &        !! ATTENTION: for each PFT they are PERCENTAGE of monoterpene EF
  &   0.105,   0.028,   0.138,   0.120,   0.031,   0.031  /)

  REAL, PARAMETER, DIMENSION(nvmc) :: em_factor_othermonot_mtc = &   !! Other monoterpenes emission factor percentage
  & (/   0.,   0.062,   0.062,   0.046,   0.054,   0.052,   0.046, &        !! ATTENTION: for each PFT they are PERCENTAGE of monoterpene EF
  &   0.144,   0.003,   0.115,   0.103,   0.108,   0.108  /)

  REAL, PARAMETER, DIMENSION(nvmc) :: em_factor_sesquiterp_mtc = &   !! Sesquiterpene emission factor
  & (/   0.,  0.45,   0.45,   0.13,   0.30,   0.36,   0.15, &               !! @tex $(\mu gC.g^{-1}.h^{-1})$ @endtex
  &    0.30,  0.25,   0.60,   0.60,   0.08,   0.08  /)

  REAL, PARAMETER :: beta_mono_mtc = 0.10                            !! Monoterpenes temperature dependency coefficient 
  REAL, PARAMETER :: beta_sesq_mtc = 0.17                            !! Sesquiterpenes temperature dependency coefficient 
  REAL, PARAMETER :: beta_meth_mtc = 0.08                            !! Methanol temperature dependency coefficient 
  REAL, PARAMETER :: beta_acet_mtc = 0.10                            !! Acetone temperature dependency coefficient 
  REAL, PARAMETER :: beta_oxyVOC_mtc = 0.13                          !! Other oxygenated BVOC temperature dependency coefficient 


  REAL, PARAMETER, DIMENSION(nvmc) :: em_factor_ORVOC_mtc = &        !! ORVOC emissions factor 
  &  (/  0.,    1.5,    1.5,    1.5,    1.5,   1.5,    1.5,  &              !! @tex $(\mu gC.g^{-1}.h^{-1})$ @endtex
  &     1.5,    1.5,    1.5,    1.5,    1.5,   1.5  /) 

  REAL, PARAMETER, DIMENSION(nvmc) :: em_factor_OVOC_mtc = &         !! OVOC emissions factor 
  &  (/  0.,    1.5,    1.5,    1.5,    1.5,   1.5,    1.5,  &              !! @tex $(\mu gC.g^{-1}.h^{-1})$ @endtex
  &     1.5,    1.5,    1.5,    1.5,    1.5,   1.5  /)
  
  REAL, PARAMETER, DIMENSION(nvmc) :: em_factor_MBO_mtc = &          !! MBO emissions factor
  & (/     0., 2.e-5, 2.e-5,   1.4, 2.e-5, 2.e-5, 0.14,  &                  !! @tex $(\mu gC.g^{-1}.h^{-1})$ @endtex
  &     2.e-5, 2.e-5, 2.e-5, 2.e-5, 2.e-5, 2.e-5  /)  
  
  REAL, PARAMETER, DIMENSION(nvmc) :: em_factor_methanol_mtc = &     !! Methanol emissions factor 
  & (/  0.,    0.8,   0.8,   1.8,   0.9,   1.9,   1.8,  &                   !! @tex $(\mu gC.g^{-1}.h^{-1})$ @endtex
  &    1.8,    1.8,   0.7,   0.9,    2.,     2.  /)  
  
  REAL, PARAMETER, DIMENSION(nvmc) :: em_factor_acetone_mtc = &      !! Acetone emissions factor
  & (/  0.,   0.25,   0.25,   0.30,   0.20,   0.33,   0.30,  &              !! @tex $(\mu gC.g^{-1}.h^{-1})$ @endtex
  &   0.25,   0.25,   0.20,   0.20,   0.08,   0.08  /)
  
  REAL, PARAMETER, DIMENSION(nvmc) :: em_factor_acetal_mtc = &       !! Acetaldehyde emissions factor 
  & (/  0.,   0.2,    0.2,     0.2,   0.2,   0.25,   0.25,   0.16,   &      !! @tex $(\mu gC.g^{-1}.h^{-1})$ @endtex
  &   0.16,   0.12,   0.12,   0.035,   0.020  /)  
  
  REAL, PARAMETER, DIMENSION(nvmc) :: em_factor_formal_mtc = &       !! Formaldehyde emissions factor
  & (/  0.,   0.04,   0.04,  0.08,    0.04,    0.04,  0.04,  &              !! @tex $(\mu gC.g^{-1}.h^{-1})$ @endtex
  &   0.04,   0.04,  0.025, 0.025,   0.013,   0.013  /)  

  REAL, PARAMETER, DIMENSION(nvmc) :: em_factor_acetic_mtc = &       !! Acetic Acid emissions factor
  & (/   0.,   0.025,   0.025,   0.025,   0.022,   0.08,   0.025,   &      !! @tex $(\mu gC.g^{-1}.h^{-1})$ @endtex
  &   0.022,   0.013,   0.012,   0.012,   0.008,   0.008  /)  

  REAL, PARAMETER, DIMENSION(nvmc) :: em_factor_formic_mtc = &       !! Formic Acid emissions factor
  & (/  0.,  0.015,  0.015,   0.02,    0.02,   0.025,  0.025,  &            !! @tex $(\mu gC.g^{-1}.h^{-1})$ @endtex
  &  0.015,  0.015,  0.010,  0.010,   0.008,   0.008  /)  

  REAL,PARAMETER, DIMENSION(nvmc) :: em_factor_no_wet_mtc = &        !! NOx emissions factor soil emissions and exponential
  & (/  0.,   2.6,   0.06,   0.03,   0.03,   0.03,   0.03,  &               !! dependancy factor for wet soils
  &  0.03,   0.03,   0.36,   0.36,   0.36,   0.36  /)                       !! @tex $(ngN.m^{-2}.s^{-1})$ @endtex

  REAL,PARAMETER, DIMENSION(nvmc) :: em_factor_no_dry_mtc = &        !! NOx emissions factor soil emissions and exponential
  & (/  0.,   8.60,   0.40,   0.22,   0.22,   0.22,   0.22,  &              !! dependancy factor for dry soils
  &   0.22,   0.22,   2.65,   2.65,   2.65,   2.65  /)                      !! @tex $(ngN.m^{-2}.s^{-1})$ @endtex 

  REAL, PARAMETER, DIMENSION(nvmc) :: Larch_mtc = &                  !! Larcher 1991 SAI/LAI ratio (unitless)
  & (/   0.,   0.015,   0.015,   0.003,   0.005,   0.005,   0.003,  &
  &   0.005,   0.003,   0.005,   0.005,   0.008,   0.008  /)  

! DSG: recalibration. Because of the unrealistic CUE, biomass turnover, and
! biomass stocks, I had to recalibrate the C cycle by trail and error using some
! rough range of CUE, Biomass and biomass turnover from Bloom et al
! (2016),Thurner et al (2015).

! One issue here is that some of the PFT tend to be at one of the boundaries for
! latosa which makes the calibration not straight forward.  Other PFTs go the whole spread 
! and thereby the biomass stocks can vary from super low to way too high for a
! given LAI. The original calibration as documented in Naudts et al (2015) uses no biomass
! (turnover) constraint for the calibration of LATOSA which should be done in a
! future recalibration. I can get a perfect LAImax using different sets of
! parameter values resulting biomass stocks and biomass turnover varying an order of magnitude.

! This certainly needs to be improved, but for now it gives a more or less acceptable status of C cycle
! This is updated from JC global simulations
  REAL, PARAMETER, DIMENSION(nvmc) :: k_latosa_max_mtc = &  !! Maximum leaf-to-sapwood area ratio 
  & (/ undef,  4000.,  4000.,  1600.,  1800.,  2100.,  1500.,  &   !! mixture of values from NV/SL and result of recalibration
  &    3100.,  2500.,  2.0,  1.0,  2.5,  3.0 /)                     

  REAL, PARAMETER, DIMENSION(nvmc) :: k_latosa_min_mtc = &  !! Minimum leaf-to-sapwood area ratio as defined in McDowell et al
  & (/ undef,  1475.,  1475.,  417.,  675.,  1600.,  400.,  &     !! 2002, Oecologia 
  &    2400.,  1500.,  1.5,  0.5,  1.75,  1.5  /)                  !! 

! DSG: recalibration by JC

  REAL, PARAMETER, DIMENSION(nvmc) :: k_root_mtc = &        !! Fine root specific conductivity. Values 
  & (/ undef,   10.,   10.,   3.3,   25.6,   2.51,   2.38,    &    !! from  ORCHIDEE-CAN
  &       2.51,    3.,   1.,  1.,  1.,  1.   /)*1.e-7         
                                                                       
  REAL, PARAMETER, DIMENSION(nvmc) :: k_sap_mtc = &         !! Maximal sapwood specific conductivity. Values compiled in T. Hickler 
  & (/-9999., 1.E-3, 1.E-3, 5.34E-4, 1.08E-4, 3.E-3, 6.25E-4,&     !! et al. 2006. @tex $(m^{2} s^{-1} MPa^{-1})$ @endtex 
  &    3.E-3, 5.82E-4, 3.E-4, 3.E-4, 3.E-4, 3.E-4    /)            !! Values from DOFOCO run.def

  REAL, PARAMETER, DIMENSION(nvmc) :: lai_to_height_mtc = &                    !! Convertion from lai to height for grasses
  &(/ undef, undef, undef, undef, undef, undef, undef, &                              !! and cropland. Convert lai because that way a dynamic
  &   undef, undef,   0.1,   0.2,   0.1,   0.1 /)                                     !! sla is accounted for  

  REAL, PARAMETER, DIMENSION(nvmc) :: deleuze_a_mtc = &                          !! intercept of the intra-tree competition within a stand
  & (/ undef,  0.23,  0.23,  0.23,  0.23,  0.23,  0.23, &                               !! based on the competion rule of Deleuze and Dhote 2004
  &     0.23,  0.23, undef, undef, undef, undef /)                                      !! Used when n_circ > 6

  REAL, PARAMETER, DIMENSION(nvmc) :: deleuze_b_mtc = &                          !! slope of the intra-tree competition within a stand
  & (/ undef,  0.58,  0.58,  0.58,  0.58,  0.58,  0.58, &                               !! based on the competion rule of Deleuze and Dhote 2004
  &     0.58,  0.58, undef, undef, undef, undef /)                                      !! Used when n_circ > 6

  REAL, PARAMETER, DIMENSION(nvmc) :: deleuze_p_all_mtc = &                      !! Percentile of the circumferences that receives photosynthates
  & (/ undef,  0.50,  0.50,  0.99,  0.99,  0.99,  0.99, &                               !! based on the competion rule of Deleuze and Dhote 2004
  &     0.99,  0.99, undef, undef, undef, undef /)                                      !! Used when n_circ > 6 for FM1, FM2 and FM4

  REAL, PARAMETER, DIMENSION(nvmc) :: m_dv_mtc = &                               !! Parameter in the Deleuze & Dhote allocation 
  & (/ undef,  1.05,  1.05,  1.05,  1.05,  1.05,  1.05, &                               !! rule that relaxes the cut-off imposed by
  &     1.05,  1.05,    0.,    0.,    0.,    0. /)                                      !! ::sigma. Owing to m_relax trees still grow
                                                                                        !! a little when their ::circ is below ::sigma
  
  REAL, PARAMETER, DIMENSION(nvmc) :: fruit_alloc_mtc = &   !! Fraction of biomass allocated to fruit production (0-1)
  & (/ undef,   0.1,    0.1,    0.1,     0.1,     0.05,    0.1, &   !! currently only parameterized for forest PFTs 
  &      0.05,   0.05,  0.1,    0.1,     0.1,     0.1 /)  

  REAL, PARAMETER, DIMENSION(nvmc) :: labile_reserve_mtc = &                   !! The lab_fac is divided by this value to obtain 
  &(/  0.,  60.,  30.,  60.,  60.,  30.,  60., &                                      !! a new parameter. This new parameter is a fraction 
  &    30.,  30.,  30.,  30.,  30.,  30.  /)                                          !! that is multiplied with the plant biomass to obatin 
                                                                                      !! the optimal size of the labile pool. The dependency
                                                                                      !! on lab_fac is a nice feature but the whole 
                                                                                      !! parameterization is arbitrary

  REAL, PARAMETER, DIMENSION(nvmc) :: evergreen_reserve_mtc = &                !! Fraction of sapwood mass stored in the reserve pool of evergreen 
  &(/ undef, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, &                                    !! trees (unitless, 0-1)
  &    0.05,  0.05, 0.05, 0.05, 0.05, 0.05 /)

  REAL, PARAMETER, DIMENSION(nvmc) :: senescense_reserve_mtc = &               !! Fraction of sapwood mass stored in the reserve pool of deciduous 
  &(/ undef, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, &                                    !! trees during senescense(unitless, 0-1)
  &    0.15,  0.15, 0.15, 0.15, 0.15, 0.15 /)

  REAL, PARAMETER, DIMENSION(nvmc) :: deciduous_reserve_mtc = &                !! Fraction of sapwood mass stored in the reserve pool of deciduous 
  &(/ undef, 0.12, 0.12, 0.12, 0.12, 0.12, 0.12, &                                    !! trees during the growing season (unitless, 0-1) 
  &    0.12,  0.12, 0.12, 0.12, 0.12, 0.12 /)
 
  REAL, PARAMETER, DIMENSION(nvmc) :: fcn_root_mtc = &      !! N/C of "root" for allocation relative to leaf N/C  according 
  & (/ undef,   .86,    .86,    .86,    .86,    .86,   .86,   &    !! to stich et al 2003
  &      .86,   .86,    .86,    .86,    .86,    .86   /)

   REAL, PARAMETER, DIMENSION(nvmc) :: fcn_wood_mtc = &      !! C/N of "wood" for allocation relative to leaf C/N  according 
   & (/ undef,    .108,   .110,   .134,   .128,   .110,  .134,  &   !! to minimum C/N of leaves(cn_leaf_min) 
   &     .110,    .077,     1.,     1.,     1.,     1.   /)         !! values are from Wang et al (2010)

  REAL, PARAMETER, DIMENSION(nvmc) :: branch_ratio_mtc =  &  !! Ratio of branches to total woody biomass (unitless)
  & (/  0.0,   0.38,   0.38,   0.25,   0.38,   0.38,   0.25,  &
  &    0.38,   0.25,    0.0,    0.0,    0.0,    0.0 /)

  REAL, PARAMETER, DIMENSION(nvmc) :: cn_leaf_init_mtc = &       !! in case of impose_np, we assume this leaf C:N ratios
  & (/  undef,  18.,  21.,  40.,  30.,  25.,  35.0, &                   !! the closer these values are to the dynamical computed C:N ratio the faster the
  &       30.,  30.,  35.,  35.,  35.,  35.  /)                         !! spinp, need to be ajusted when we have a global run

  REAL, PARAMETER, DIMENSION(nvmc) :: np_leaf_init_mtc =      &  !! in case of impose_np, we assume this leaf N:P ratios
  & (/  undef, 20.,  20.,  20.,  20.,  20.,  20.,  &                    !! the closer these values are to the dynamical computed N:P ratio the faster the
  &     20.,  20.,  20.,  20.,  20.,  20.   /)                          !! spinp, need to be ajusted when we have a global run

  REAL, PARAMETER, DIMENSION(nvmc) :: fnp_wood_mtc =       &     !! N/P of "wood" for allocation relative to np_leaf_min=5 (rigid_woodNP=TRUE) according 
  & (/ undef,    1.1,   1.1,   1.1,   1.1,   1.1,  1.1,           &     !! Sardans (2015) reports a value of 11.2 for N:P wood
  &       1.1,    1.1,   1.1,   1.1,  1.1,   1.1   /)                   !! 
  
  REAL, PARAMETER, DIMENSION(nvmc) :: fnp_root_mtc = &      !! N/P of "root" for allocation relative to actual leaf
  & (/ undef,    1.,   1.,   1.,   1.,   1.,  1.,           &      !! due to lack of data, we assume it to be same as leaves
  &       1.,    1.,   1.,   1.,   1.,   1.   /)                   !! 

              
END MODULE constantes_mtc
