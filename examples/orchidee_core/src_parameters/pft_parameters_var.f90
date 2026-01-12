! =================================================================================================================================
! MODULE       : pft_parameters_var
!
! CONTACT      : orchidee-help _at_ listes.ipsl.fr
!
! LICENCE      : IPSL (2011)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF        This module contains the variables in function of plant funtional type (pft).
!!
!!\n DESCRIPTION: This module contains the declarations for the externalized variables in function of the 
!!                plant foncional type(pft). \n
!!                The module is already USE in module pft_parameters. Therefor no need to USE it seperatly except
!!                if the subroutines in module pft_parameters are not needed.\n
!!
!! RECENT CHANGE(S): None
!!
!! REFERENCE(S)	: None
!!
!! SVN          :
!! $HeadURL: $
!! $Date: 2021-05-06 11:05:55 +0200 (四, 2021-05-06) $
!! $Revision: 7177 $
!! \n
!_ ================================================================================================================================

MODULE pft_parameters_var

  USE defprec
  
  IMPLICIT NONE


  !
  ! PFT GLOBAL
  !
  INTEGER, SAVE :: nvm = 13                               !! Number of vegetation types (2-N, unitless)
!$OMP THREADPRIVATE(nvm)

  INTEGER, ALLOCATABLE, SAVE, DIMENSION(:) :: pft_to_mtc  !! Table of conversion : we associate one pft to one metaclass 
                                                                 !! (1-13, unitless)
!$OMP THREADPRIVATE(pft_to_mtc)

  CHARACTER(LEN=34), ALLOCATABLE, SAVE, DIMENSION(:) :: PFT_name !! Description of the PFT (unitless)
!$OMP THREADPRIVATE(PFT_name)

  LOGICAL, SAVE   :: l_first_pft_parameters = .TRUE.             !! To keep first call trace of the module (true/false)
!$OMP THREADPRIVATE(l_first_pft_parameters)

  !
  ! VEGETATION STRUCTURE
  !
  INTEGER, ALLOCATABLE, SAVE, DIMENSION(:) :: leaf_tab       !! leaf type (1-4, unitless)
                                                                    !! 1=broad leaved tree, 2=needle leaved tree, 
                                                                    !! 3=grass 4=bare ground 
!$OMP THREADPRIVATE(leaf_tab)

  CHARACTER(len=6), ALLOCATABLE, SAVE, DIMENSION(:) :: pheno_model  !! which phenology model is used? (tabulated) (unitless)
!$OMP THREADPRIVATE(pheno_model)

  LOGICAL, ALLOCATABLE, SAVE, DIMENSION(:) :: is_tree               !! Is the vegetation type a tree ? (true/false)
!$OMP THREADPRIVATE(is_tree)

  LOGICAL, ALLOCATABLE, SAVE, DIMENSION(:) :: is_deciduous          !! Is PFT deciduous ? (true/false)
!$OMP THREADPRIVATE(is_deciduous)

  LOGICAL, ALLOCATABLE, SAVE, DIMENSION(:) :: is_evergreen          !! Is PFT evegreen ? (true/false)
!$OMP THREADPRIVATE(is_evergreen)

  LOGICAL, ALLOCATABLE, SAVE, DIMENSION(:) :: is_needleleaf         !! Is PFT needleleaf ? (true/false)
!$OMP THREADPRIVATE(is_needleleaf)
 
  LOGICAL, ALLOCATABLE, SAVE, DIMENSION(:) :: is_tropical           !! Is PFT tropical ? (true/false)
!$OMP THREADPRIVATE(is_tropical)

  LOGICAL, ALLOCATABLE, SAVE, DIMENSION(:) :: natural               !! natural? (true/false)
!$OMP THREADPRIVATE(natural)
! dgvmjc
  LOGICAL, ALLOCATABLE, SAVE, DIMENSION(:) :: pasture               !! pasture? (true/false)
!$OMP THREADPRIVATE(pasture)
  LOGICAL, ALLOCATABLE, SAVE, DIMENSION(:) :: rangeland               !! rangeland? (true/false)
!$OMP THREADPRIVATE(rangeland)
  LOGICAL, ALLOCATABLE, SAVE, DIMENSION(:) :: bioenergy             !! bioenergy plantation? (true/false)
!$OMP THREADPRIVATE(bioenergy)
  LOGICAL, ALLOCATABLE, SAVE, DIMENSION(:) :: is_agri_fert          !! is agricultural land fertilized? (true/false)
! end dgvmjc
! firejc
  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: frac_co_fire       !! CO fraction in fire C emission
  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: frac_ch4_fire       !! CH4 fraction in fire C emission
  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: frac_Nemit_fire     !! N:C ratio during fire emission
  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: frac_Pemit_fire     !! P:C ratio during fire emission
  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: frac_nh3_fire       !! NH3 fraction in fire N emission
  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: frac_nox_fire       !! NOx fraction in fire N emission
  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: frac_n2o_fire       !! N2O fraction in fire N emission
! end firejc

  CHARACTER(len=5), ALLOCATABLE, SAVE, DIMENSION(:) :: type_of_lai  !! Type of behaviour of the LAI evolution algorithm 
                                                                    !! for each vegetation type.
                                                                    !! Value of type_of_lai, one for each vegetation type :
                                                                    !! mean or interp
!$OMP THREADPRIVATE(type_of_lai)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: veget_ori_fixed_test_1 !! Value for veget_ori for tests in 0-dim simulations 
                                                                         !! (0-1, unitless)
!$OMP THREADPRIVATE(veget_ori_fixed_test_1)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: llaimax                !! laimax for maximum lai see also type of lai 
                                                                         !! interpolation
                                                                         !! @tex $(m^2.m^{-2})$ @endtex 
!$OMP THREADPRIVATE(llaimax)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: llaimin                !! laimin for minimum lai see also type of lai 
                                                                         !! interpolation 
                                                                         !! @tex $(m^2.m^{-2})$ @endtex
!$OMP THREADPRIVATE(llaimin)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: height_presc           !! prescribed height of vegetation.(m)
                                                                         !! Value for height_presc : one for each vegetation type
!$OMP THREADPRIVATE(height_presc)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: z0_over_height        !! Factor to calculate roughness height from 
                                                                        !! vegetation height (unitless)   
!$OMP THREADPRIVATE(z0_over_height)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: ratio_z0m_z0h         !! Ratio between z0m and z0h
!$OMP THREADPRIVATE(ratio_z0m_z0h)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) ::  rveg_pft              !! Potentiometer to set vegetation resistance (unitless)
                                                                         !! Nathalie on March 28th, 2006 - from Fred Hourdin,
!$OMP THREADPRIVATE(rveg_pft)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: sla                    !! specif leaf area @tex $(m^2.gC^{-1})$ @endtex 
!$OMP THREADPRIVATE(sla)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: availability_fact      !! calculate dynamic mortality in lpj_gap
!$OMP THREADPRIVATE(availability_fact)

  !
  ! EVAPOTRANSPIRATION (sechiba)
  !
  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: rstruct_const          !! Structural resistance. 
                                                                         !! Value for rstruct_const : one for each vegetation type
                                                                         !! @tex $(s.m^{-1})$ @endtex
!$OMP THREADPRIVATE(rstruct_const)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: kzero                  !! A vegetation dependent constant used in the calculation
                                                                         !! of the surface resistance. 
                                                                         !! Value for kzero one for each vegetation type
                                                                         !! @tex $(kg.m^2.s^{-1})$ @endtex
!$OMP THREADPRIVATE(kzero)

  !
  ! WATER (sechiba)
  !
  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: wmax_veg  !! Volumetric available soil water capacity in each PFT
                                                            !! @tex $(kg.m^{-3} of soil)$ @endtex
!$OMP THREADPRIVATE(wmax_veg)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: humcste   !! Root profile description for the different vegetation types.
                                                            !! These are the factor in the exponential which gets
                                                            !! the root density as a function of depth 
                                                            !! @tex $(m^{-1})$ @endtex
!$OMP THREADPRIVATE(humcste)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: throughfall_by_pft !! Fraction of rain intercepted by the canopy (0-100, unitless)
!$OMP THREADPRIVATE(throughfall_by_pft)

  !
  ! ALBEDO (sechiba)
  !
  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: snowa_aged_vis !! Minimum snow albedo value for each vegetation type
                                                                 !! after aging (dirty old snow), visible albedo (unitless)
                                                                 !! Source : Values are from the Thesis of S. Chalita (1992)
!$OMP THREADPRIVATE(snowa_aged_vis)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: snowa_aged_nir !! Minimum snow albedo value for each vegetation type
                                                                 !! after aging (dirty old snow), near infrared albedo (unitless)
                                                                 !! Source : Values are from the Thesis of S. Chalita (1992)
!$OMP THREADPRIVATE(snowa_aged_nir)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: snowa_dec_vis  !! Decay rate of snow albedo value for each vegetation type
                                                                 !! as it will be used in condveg_snow, visible albedo (unitless)
                                                                 !! Source : Values are from the Thesis of S. Chalita (1992)
!$OMP THREADPRIVATE(snowa_dec_vis)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: snowa_dec_nir  !! Decay rate of snow albedo value for each vegetation type
                                                                 !! as it will be used in condveg_snow, near infrared albedo (unitless)
                                                                 !! Source : Values are from the Thesis of S. Chalita (1992)
!$OMP THREADPRIVATE(snowa_dec_nir)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: alb_leaf_vis  !! leaf albedo of vegetation type, visible albedo (unitless)
!$OMP THREADPRIVATE(alb_leaf_vis)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: alb_leaf_nir  !! leaf albedo of vegetation type, near infrared albedo (unitless)
!$OMP THREADPRIVATE(alb_leaf_nir)

  !
  ! SOIL - VEGETATION
  !
  INTEGER, ALLOCATABLE, SAVE, DIMENSION(:) :: pref_soil_veg      !! Table which contains the correlation between the soil
                                                                        !! types and vegetation type. Two modes exist :
                                                                        !! 1) pref_soil_veg = 0 then we have an equidistribution
                                                                        !!    of vegetation on soil types
                                                                        !! 2) Else for each pft the prefered soil type is given :
                                                                        !!    1=sand, 2=loan, 3=clay
                                                                        !! This variable is initialized in slowproc.(1-3, unitless)
!$OMP THREADPRIVATE(pref_soil_veg)

  !
  ! PHOTOSYNTHESIS
  !
  !-
  ! 1. CO2
  !-
  LOGICAL, ALLOCATABLE, SAVE, DIMENSION(:) :: is_c4             !! flag for C4 vegetation types (true/false)
!$OMP THREADPRIVATE(is_c4)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: vcmax_fix     !! values used for vcmax when STOMATE is not activated
                                                                !! @tex $(\mu mol.m^{-2}.s^{-1})$ @endtex
!$OMP THREADPRIVATE(vcmax_fix)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: downregulation_co2_coeff !! Coefficient for CO2 downregulation (unitless)
!$OMP THREADPRIVATE(downregulation_co2_coeff)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: E_KmC         !! Energy of activation for KmC (J mol-1)
!$OMP THREADPRIVATE(E_KmC)                                                               
  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: E_KmO         !! Energy of activation for KmO (J mol-1)
!$OMP THREADPRIVATE(E_KmO)          
REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: E_Sco           !! Energy of activation for Sco (J mol-1)
!$OMP THREADPRIVATE(E_Sco)            
  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: E_gamma_star  !! Energy of activation for gamma_star (J mol-1)
!$OMP THREADPRIVATE(E_gamma_star)    
  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: E_Vcmax       !! Energy of activation for Vcmax (J mol-1) 
!$OMP THREADPRIVATE(E_Vcmax)                                                              
  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: E_Jmax        !! Energy of activation for Jmax (J mol-1)
!$OMP THREADPRIVATE(E_Jmax) 
  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: aSV           !! a coefficient of the linear regression (a+bT) defining the Entropy term for Vcmax (J K-1 mol-1)
!$OMP THREADPRIVATE(aSV)    
  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: bSV           !! b coefficient of the linear regression (a+bT) defining the Entropy term for Vcmax (J K-1 mol-1 °C-1)
!$OMP THREADPRIVATE(bSV)
  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: tphoto_min   !! minimum photosynthesis temperature (deg C)
!$OMP THREADPRIVATE(tphoto_min) 
  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: tphoto_max   !! maximum photosynthesis temperature (deg C)
!$OMP THREADPRIVATE(tphoto_max) 
  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: aSJ           !! a coefficient of the linear regression (a+bT) defining the Entropy term for Jmax (J K-1 mol-1)
!$OMP THREADPRIVATE(aSJ)    
  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: bSJ           !! b coefficient of the linear regression (a+bT) defining the Entropy term for Jmax (J K-1 mol-1 °C-1)
!$OMP THREADPRIVATE(bSJ)    
  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: D_Vcmax       !! Energy of deactivation for Vcmax (J mol-1)
!$OMP THREADPRIVATE(D_Vcmax)                     
  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: D_Jmax        !! Energy of deactivation for Jmax (J mol-1)
!$OMP THREADPRIVATE(D_Jmax)                           

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: E_gm          !! Energy of activation for gm (J mol-1) 
!$OMP THREADPRIVATE(E_gm)                                       
  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: S_gm          !! Entropy term for gm (J K-1 mol-1) 
!$OMP THREADPRIVATE(S_gm)                                       
  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: D_gm          !! Energy of deactivation for gm (J mol-1) 
!$OMP THREADPRIVATE(D_gm)                                       
         
  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: E_Rd          !! Energy of activation for Rd (J mol-1)
!$OMP THREADPRIVATE(E_Rd)                                      
  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: Vcmax25       !! Maximum rate of Rubisco activity-limited carboxylation at 25°C
                                                                !! @tex $(\mu mol.m^{-2}.s^{-1})$ @endtex
!$OMP THREADPRIVATE(Vcmax25)
  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: arJV          !! a coefficient of the linear regression (a+bT) defining the Jmax25/Vcmax25 ratio (mu mol e- (mu mol CO2)-1)
!$OMP THREADPRIVATE(arJV)
  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: brJV          !! b coefficient of the linear regression (a+bT) defining the Jmax25/Vcmax25 ratio (mu mol e- (mu mol CO2)-1)
!$OMP THREADPRIVATE(brJV)
  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: KmC25         !! Michaelis–Menten constant of Rubisco for CO2 at 25°C (ubar)
!$OMP THREADPRIVATE(KmC25)                                     
  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: KmO25         !! Michaelis–Menten constant of Rubisco for O2 at 25°C (ubar)
!$OMP THREADPRIVATE(KmO25)                
REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: Sco25           !! Relative CO2 /O2 specificity factor for Rubisco at 25Â°C (bar bar-1)
!$OMP THREADPRIVATE(Sco25)      
  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: gamma_star25  !! Ci-based CO2 compensation point in the absence of Rd at 25°C (ubar)
!$OMP THREADPRIVATE(gamma_star25)       
  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: gm25         !! Mesophyll diffusion conductance at 25°C (mol m−2 s−1 bar−1)
!$OMP THREADPRIVATE(gm25)      
  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: a1            !! Empirical factor involved in the calculation of fvpd (-)
!$OMP THREADPRIVATE(a1)                                        
  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: b1            !! Empirical factor involved in the calculation of fvpd (-)
!$OMP THREADPRIVATE(b1)                                        
  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: g0            !! Residual stomatal conductance when irradiance approaches zero (mol m−2 s−1 bar−1)
!$OMP THREADPRIVATE(g0)                                        
  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: h_protons     !! Number of protons required to produce one ATP (mol mol-1)
!$OMP THREADPRIVATE(h_protons)                                 
  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: fpsir         !! Fraction of PSII e− transport rate partitioned to the C4 cycle (-)
!$OMP THREADPRIVATE(fpsir)                                         
  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: fQ            !! Fraction of electrons at reduced plastoquinone that follow the Q-cycle (-) - Values for C3 platns are not used
!$OMP THREADPRIVATE(fQ)                                        
  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: fpseudo       !! Fraction of electrons at PSI that follow  pseudocyclic transport (-) - Values for C3 platns are not used
!$OMP THREADPRIVATE(fpseudo)                                   
  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: kp            !! Initial carboxylation efficiency of the PEP carboxylase (mol m−2 s−1 bar−1) 
!$OMP THREADPRIVATE(kp)                                        
  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: alpha         !! Fraction of PSII activity in the bundle sheath (-)
!$OMP THREADPRIVATE(alpha)                                     
  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: gbs           !! Bundle-sheath conductance (mol m−2 s−1 bar−1)
!$OMP THREADPRIVATE(gbs)                                       
  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: theta         !! Convexity factor for response of J to irradiance (-)
!$OMP THREADPRIVATE(theta)                                     
  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: alpha_LL      !! Conversion efficiency of absorbed light into J at strictly limiting light (mol e− (mol photon)−1)
!$OMP THREADPRIVATE(alpha_LL)
  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: stress_vcmax  !! Stress on vcmax
!$OMP THREADPRIVATE(stress_vcmax)
  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: stress_gs     !! Stress on vcmax
!$OMP THREADPRIVATE(stress_gs)
  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: stress_gm     !! Stress on vcmax
!$OMP THREADPRIVATE(stress_gm)

  !-
  ! 2. Stomate
  !-
  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: ext_coeff     !! extinction coefficient of the Monsi&Saeki relationship (1953)
                                                                !! (unitless)
!$OMP THREADPRIVATE(ext_coeff)
  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: nue_opt       !! Nitrogen use efficiency of Vcmax 
                                                                !! ((mumol[CO2] s-1) (gN[leaf])-1)
                                                                !! based on the work of Kattge et al. (2009, GCB)
!$OMP THREADPRIVATE(nue_opt)
  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: ext_coeff_vegetfrac     !! extinction coefficient used for the calculation of the 
                                                                !! bare soil fraction (unitless)
!$OMP THREADPRIVATE(ext_coeff_vegetfrac)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: ext_coeff_N   !! extinction coefficient of the leaf N content profile within the canopy
                                                                !! ((m2[ground]) (m-2[leaf]))
                                                                !! based on Dewar et al. (2012, value of 0.18), on Carswell et al. (2000, value of 0.11 used in OCN) 
!$OMP THREADPRIVATE(ext_coeff_N)



  !
  ! ALLOCATION (stomate)
  !
  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: R0            !! Default root allocation (0-1, unitless)
!$OMP THREADPRIVATE(R0)
  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: S0            !! Default sapwood allocation (0-1, unitless)
!$OMP THREADPRIVATE(S0)
  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: L0            !! Default leaf allocation (0-1, unitless)
!$OMP THREADPRIVATE(L0)


  !
  ! RESPIRATION (stomate)
  !
  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: frac_growthresp  !! fraction of GPP which is lost as growth respiration

!$OMP THREADPRIVATE(frac_growthresp)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:) :: maint_resp_slope  !! slope of maintenance respiration coefficient 
                                                                      !! (1/K, 1/K^2, 1/K^3), used in the code
!$OMP THREADPRIVATE(maint_resp_slope)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: maint_resp_slope_c  !! slope of maintenance respiration coefficient (1/K),
                                                                      !! constant c of aT^2+bT+c , tabulated
!$OMP THREADPRIVATE(maint_resp_slope_c)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: maint_resp_slope_b  !! slope of maintenance respiration coefficient (1/K), 
                                                                      !! constant b of aT^2+bT+c , tabulated
!$OMP THREADPRIVATE(maint_resp_slope_b)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: maint_resp_slope_a  !! slope of maintenance respiration coefficient (1/K), 
                                                                      !! constant a of aT^2+bT+c , tabulated
!$OMP THREADPRIVATE(maint_resp_slope_a)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:) :: coeff_maint_zero  !! maintenance respiration coefficient at 0 deg C, 
                                                                      !! @tex $(gC.gC^{-1}.day^{-1})$ @endtex
!$OMP THREADPRIVATE(coeff_maint_zero)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: cm_zero_leaf        !! maintenance respiration coefficient at 0 deg C,
                                                                      !! for leaves, tabulated 
                                                                      !! @tex $(gC.gC^{-1}.day^{-1})$ @endtex
!$OMP THREADPRIVATE(cm_zero_leaf)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: cm_zero_sapabove    !! maintenance respiration coefficient at 0 deg C,
                                                                      !! for sapwood above, tabulated
                                                                      !! @tex $(gC.gC^{-1}.day^{-1})$ @endtex
!$OMP THREADPRIVATE(cm_zero_sapabove)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: cm_zero_sapbelow    !! maintenance respiration coefficient at 0 deg C,
                                                                      !! for sapwood below, tabulated 
                                                                      !! @tex $(gC.gC^{-1}.day^{-1})$ @endtex
!$OMP THREADPRIVATE(cm_zero_sapbelow)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: cm_zero_heartabove  !! maintenance respiration coefficient at 0 deg C
                                                                      !! for heartwood above, tabulated 
                                                                      !! @tex $(gC.gC^{-1}.day^{-1})$ @endtex
!$OMP THREADPRIVATE(cm_zero_heartabove)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: cm_zero_heartbelow  !! maintenance respiration coefficient at 0 deg C,
                                                                      !! for heartwood below, tabulated 
                                                                      !! @tex $(gC.gC^{-1}.day^{-1})$ @endtex
!$OMP THREADPRIVATE(cm_zero_heartbelow)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: cm_zero_root        !! maintenance respiration coefficient at 0 deg C,
                                                                      !! for roots, tabulated
                                                                      !! @tex $(gC.gC^{-1}.day^{-1})$ @endtex
!$OMP THREADPRIVATE(cm_zero_root)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: cm_zero_fruit       !! maintenance respiration coefficient  at 0 deg C,
                                                                      !! for fruits, tabulated
                                                                      !! @tex $(gC.gC^{-1}.day^{-1})$ @endtex
!$OMP THREADPRIVATE(cm_zero_fruit)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: cm_zero_carbres     !! maintenance respiration coefficient at 0 deg C,
                                                                      !! for carbohydrate reserve, tabulated
                                                                      !! @tex $(gC.gC^{-1}.day^{-1})$ @endtex
!$OMP THREADPRIVATE(cm_zero_carbres)
  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: coeff_maint_init    !! maintenance respiration coefficient at 10 deg C,
                                                                      !! @tex $(gC.gN^{-1}.day^{-1})$ @endtex
!$OMP THREADPRIVATE(coeff_maint_init)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: tref_maint_resp     !! maintenance respiration Temperature coefficient,
                                                                      !! @tex $(degC)$ @endtex
!$OMP THREADPRIVATE(tref_maint_resp)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: tmin_maint_resp     !! maintenance respiration Temperature coefficient,
                                                                      !! @tex $(degC)$ @endtex
!$OMP THREADPRIVATE(tmin_maint_resp)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: e0_maint_resp     !! maintenance respiration Temperature coefficient,
                                                                      !! @tex $(unitless)$ @endtex
!$OMP THREADPRIVATE(e0_maint_resp)

 
  !
  ! FIRE (stomate)
  !
  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: flam              !! flamability : critical fraction of water holding 
                                                                    !! capacity (0-1, unitless)
!$OMP THREADPRIVATE(flam)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: resist            !! fire resistance (0-1, unitless)
!$OMP THREADPRIVATE(resist)
  !
  ! grassland management
  !
  !gmjc
  ! Is the vegetation type a managed grassland ?
  LOGICAL,ALLOCATABLE, SAVE, DIMENSION (:) :: is_grassland_manag
  ! Is the vegetation type a cut grassland for management adaptation ?
  LOGICAL,ALLOCATABLE, SAVE, DIMENSION (:) :: is_grassland_cut
  ! Is the vegetation type a grazed grassland for management adaptation ?
  LOGICAL,ALLOCATABLE, SAVE, DIMENSION (:) :: is_grassland_grazed
  ! Management Intensity reading in MANAGEMENT_MAP
  INTEGER, ALLOCATABLE, SAVE, DIMENSION(:) :: management_intensity
  ! Start year of management reading in MANAGEMENT_MAP & GRAZING_MAP
  INTEGER, ALLOCATABLE, SAVE, DIMENSION(:) :: management_start
  ! Start year of N deposition reading in DEPOSITION_MAP
  INTEGER, ALLOCATABLE, SAVE, DIMENSION(:) :: deposition_start
  ! Number of year that management should be read
  INTEGER, ALLOCATABLE, SAVE, DIMENSION(:) :: nb_year_management
  ! maximum specific leaf area (m**2/gC)
  REAL,  ALLOCATABLE, SAVE, DIMENSION(:) :: sla_max
  ! minimum specific leaf area (m**2/gC)
  REAL, ALLOCATABLE, SAVE, DIMENSION(:)  :: sla_min
  !end gmjc

  !
  ! FLUX - LUC (Land Use Change)
  !
  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: coeff_lcchange_1   !! Coeff of biomass export for the year (unitless)
!$OMP THREADPRIVATE(coeff_lcchange_1)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: coeff_lcchange_10  !! Coeff of biomass export for the decade (unitless)
!$OMP THREADPRIVATE(coeff_lcchange_10)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: coeff_lcchange_100 !! Coeff of biomass export for the century (unitless)
!$OMP THREADPRIVATE(coeff_lcchange_100)
 
 
  !
  ! PHENOLOGY
  !
  !-
  ! 1. Stomate
  !-
  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: lai_max_to_happy  !! threshold of LAI below which plant uses carbohydrate reserves
!$OMP THREADPRIVATE(lai_max_to_happy)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: lai_max           !! maximum LAI, PFT-specific @tex $(m^2.m^{-2})$ @endtex
!$OMP THREADPRIVATE(lai_max)

  INTEGER, ALLOCATABLE, SAVE, DIMENSION(:) :: pheno_type     !! type of phenology (0-4, unitless)
                                                                    !! 0=bare ground 1=evergreen,  2=summergreen, 
                                                                    !! 3=raingreen,  4=perennial
                                                                    !! For the moment, the bare ground phenotype is not managed, 
                                                                    !! so it is considered as "evergreen"
!$OMP THREADPRIVATE(pheno_type)
  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: frac_bioenergy_harvest !! = bioenergy crop harvest fraction, PFT-specific @tex @endtex
!$OMP THREADPRIVATE(frac_bioenergy_harvest)

  !-
  ! 2. Leaf Onset
  !-
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:) :: pheno_gdd_crit   !! critical gdd,tabulated (C), used in the code
!$OMP THREADPRIVATE(pheno_gdd_crit)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: pheno_gdd_crit_c   !! critical gdd,tabulated (C), 
                                                                     !! constant c of aT^2+bT+c (unitless)
!$OMP THREADPRIVATE(pheno_gdd_crit_c)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: pheno_gdd_crit_b   !! critical gdd,tabulated (C), 
                                                                     !! constant b of aT^2+bT+c (unitless)
!$OMP THREADPRIVATE(pheno_gdd_crit_b)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: pheno_gdd_crit_a   !! critical gdd,tabulated (C), 
                                                                     !! constant a of aT^2+bT+c (unitless)
!$OMP THREADPRIVATE(pheno_gdd_crit_a)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: pheno_moigdd_t_crit!! Monthly avearage temperature treashold used for C4 grass (C)
!$OMP THREADPRIVATE(pheno_moigdd_t_crit)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: ngd_crit           !! critical ngd,tabulated. Threshold -5 degrees (days)
!$OMP THREADPRIVATE(ngd_crit)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: ncdgdd_temp        !! critical temperature for the ncd vs. gdd function
                                                                     !! in phenology (C)
!$OMP THREADPRIVATE(ncdgdd_temp)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: hum_frac           !! critical humidity (relative to min/max) for phenology
                                                                     !! (0-1, unitless)
!$OMP THREADPRIVATE(hum_frac)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: hum_min_time       !! minimum time elapsed since moisture minimum (days)
!$OMP THREADPRIVATE(hum_min_time)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: tau_sap            !! sapwood -> heartwood conversion time (days)
!$OMP THREADPRIVATE(tau_sap)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: tau_fruit          !! fruit lifetime (days)
!$OMP THREADPRIVATE(tau_fruit)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: tau_root           !! root turnover (1/days)
!$OMP THREADPRIVATE(tau_root)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: tau_leafinit  !! time to attain the initial foliage using the carbohydrate reserve
!$OMP THREADPRIVATE(tau_leafinit)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: ecureuil           !! fraction of primary leaf and root allocation put
                                                                     !! into reserve (0-1, unitless)
!$OMP THREADPRIVATE(ecureuil)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: alloc_min          !! NEW - allocation above/below = f(age) - 30/01/04 NV/JO/PF
!$OMP THREADPRIVATE(alloc_min)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: alloc_max          !! NEW - allocation above/below = f(age) - 30/01/04 NV/JO/PF
!$OMP THREADPRIVATE(alloc_max)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: demi_alloc         !! NEW - allocation above/below = f(age) - 30/01/04 NV/JO/PF
!$OMP THREADPRIVATE(demi_alloc)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: leaflife_tab       !! leaf longevity, tabulated (??units??)
!$OMP THREADPRIVATE(leaflife_tab)

  !-
  ! 3. Senescence
  !-
  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: leaffall              !! length of death of leaves,tabulated (days)
!$OMP THREADPRIVATE(leaffall)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: leafagecrit           !! critical leaf age,tabulated (days)
!$OMP THREADPRIVATE(leafagecrit)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: rootagecrit           !! critical leaf age,tabulated (days)
!$OMP THREADPRIVATE(rootagecrit)

  CHARACTER(len=6), ALLOCATABLE, SAVE, DIMENSION(:) :: senescence_type  !! type of senescence,tabulated (unitless)
                                                                        !! List of avaible types of senescence :
                                                                        !! 'cold  ', 'dry   ', 'mixed ', 'none  '
!$OMP THREADPRIVATE(senescence_type)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: senescence_hum        !! critical relative moisture availability for senescence
                                                                        !! (0-1, unitless)
!$OMP THREADPRIVATE(senescence_hum)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: nosenescence_hum      !! relative moisture availability above which there is
                                                                        !! no humidity-related senescence (0-1, unitless)
!$OMP THREADPRIVATE(nosenescence_hum)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: max_turnover_time     !! maximum turnover time for grasses (days)
!$OMP THREADPRIVATE(max_turnover_time)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: min_turnover_time     !! minimum turnover time for grasses (days)
!$OMP THREADPRIVATE(min_turnover_time)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: recycle_leaf          !! Fraction of N leaf that is recycled when leaves are senescent
!$OMP THREADPRIVATE(recycle_leaf)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: recycle_root          !! Fraction of N root that is recycled when leaves are senescent
!$OMP THREADPRIVATE(recycle_root)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: p_recycle_leaf          !! Fraction of P leaf that is recycled when leaves are senescent
!$OMP THREADPRIVATE(p_recycle_leaf)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: p_recycle_root          !! Fraction of P root that is recycled when leaves are senescent
!$OMP THREADPRIVATE(p_recycle_root)


  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: min_leaf_age_for_senescence  !! minimum leaf age to allow senescence g (days)
!$OMP THREADPRIVATE(min_leaf_age_for_senescence)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:) :: senescence_temp     !! critical temperature for senescence (C),
                                                                        !! used in the code
!$OMP THREADPRIVATE(senescence_temp)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: senescence_temp_c     !! critical temperature for senescence (C), 
                                                                        !! constant c of aT^2+bT+c , tabulated (unitless)
!$OMP THREADPRIVATE(senescence_temp_c)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: senescence_temp_b     !! critical temperature for senescence (C), 
                                                                        !! constant b of aT^2+bT+c , tabulated (unitless)
!$OMP THREADPRIVATE(senescence_temp_b)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: senescence_temp_a     !! critical temperature for senescence (C),
                                                                        !! constant a of aT^2+bT+c , tabulated (unitless)
!$OMP THREADPRIVATE(senescence_temp_a)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: gdd_senescence        !! minimum gdd to allow senescence of crops (days)
!$OMP THREADPRIVATE(gdd_senescence)

  !-
  ! 4. N cycle
  !-
  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: cn_leaf_min           !! minimum CN ratio of leaves (gC/gN)
!$OMP THREADPRIVATE(cn_leaf_min)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: cn_leaf_max           !! maximum CN ratio of leaves (gC/gN) 
!$OMP THREADPRIVATE(cn_leaf_max)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: max_soil_n_bnf        !! Value of total N (NH4+NO3) 
                                                                        !! above which we stop adding N via BNF
                                                                        !! (gN/m**2)
!$OMP THREADPRIVATE(max_soil_n_bnf)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: fixed_fsoiln        !! 
!$OMP THREADPRIVATE(fixed_fsoiln)

  !-
  ! 5. P cycle
  !-
  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: np_leaf_min           !! minimum NP ratio of leaves (gN/gP)
!$OMP THREADPRIVATE(np_leaf_min)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: np_leaf_max           !! maximum NP ratio of leaves (gN/gP) 
!$OMP THREADPRIVATE(np_leaf_max)

  !
  ! DGVM
  !

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: residence_time        !! residence time of trees (y) 
!$OMP THREADPRIVATE(residence_time)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: tmin_crit             !! critical tmin, tabulated (C)
!$OMP THREADPRIVATE(tmin_crit)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: tcm_crit              !! critical tcm, tabulated (C)
!$OMP THREADPRIVATE(tcm_crit)

  !
  ! Biogenic Volatile Organic Compounds
  !

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: em_factor_isoprene       !! Isoprene emission factor
                                                                           !! @tex $(\mu gC.g^{-1}.h^{-1})$ @endtex
!$OMP THREADPRIVATE(em_factor_isoprene)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: em_factor_monoterpene    !! Monoterpene emission factor
                                                                           !! @tex $(\mu gC.g^{-1}.h^{-1})$ @endtex
!$OMP THREADPRIVATE(em_factor_monoterpene)

  REAL, SAVE :: LDF_mono                                            !! monoterpenes fraction dependancy to light
!$OMP THREADPRIVATE(LDF_mono) 
  REAL, SAVE :: LDF_sesq                                            !! sesquiterpenes fraction dependancy to light
!$OMP THREADPRIVATE(LDF_sesq)
  REAL, SAVE :: LDF_meth                                            !! methanol fraction dependancy to light
!$OMP THREADPRIVATE(LDF_meth)
  REAL, SAVE :: LDF_acet                                            !! acetone fraction dependancy to light
!$OMP THREADPRIVATE(LDF_acet)
  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: em_factor_apinene        !! Alfa pinene emission factor
                                                                           !! @tex $(\mu gC.g^{-1}.h^{-1})$ @endtex
!$OMP THREADPRIVATE(em_factor_apinene)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: em_factor_bpinene        !! Beta pinene emission factor
                                                                           !! @tex $(\mu gC.g^{-1}.h^{-1})$ @endtex
!$OMP THREADPRIVATE(em_factor_bpinene)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: em_factor_limonene       !! Limonene emission factor
                                                                           !! @tex $(\mu gC.g^{-1}.h^{-1})$ @endtex
!$OMP THREADPRIVATE(em_factor_limonene)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: em_factor_myrcene        !! Myrcene emission factor
                                                                           !! @tex $(\mu gC.g^{-1}.h^{-1})$ @endtex
!$OMP THREADPRIVATE(em_factor_myrcene)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: em_factor_sabinene       !! Sabinene emission factor
                                                                           !! @tex $(\mu gC.g^{-1}.h^{-1})$ @endtex
!$OMP THREADPRIVATE(em_factor_sabinene)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: em_factor_camphene       !! Camphene emission factor
                                                                           !! @tex $(\mu gC.g^{-1}.h^{-1})$ @endtex
!$OMP THREADPRIVATE(em_factor_camphene)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: em_factor_3carene        !! 3-carene emission factor
                                                                           !! @tex $(\mu gC.g^{-1}.h^{-1})$ @endtex
!$OMP THREADPRIVATE(em_factor_3carene)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: em_factor_tbocimene      !! T-beta-ocimene emission factor
                                                                           !! @tex $(\mu gC.g^{-1}.h^{-1})$ @endtex
!$OMP THREADPRIVATE(em_factor_tbocimene)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: em_factor_othermonot     !! Other monoterpenes emission factor
                                                                           !! @tex $(\mu gC.g^{-1}.h^{-1})$ @endtex
!$OMP THREADPRIVATE(em_factor_othermonot)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: em_factor_sesquiterp     !! Sesquiterpene emission factor
                                                                           !! @tex $(\mu gC.g^{-1}.h^{-1})$ @endtex
!$OMP THREADPRIVATE(em_factor_sesquiterp)

  REAL, SAVE :: beta_mono                                           !! Monoterpenes temperature dependency coefficient 
!$OMP THREADPRIVATE(beta_mono)
  REAL, SAVE :: beta_sesq                                           !! Sesquiterpenes temperature dependency coefficient
!$OMP THREADPRIVATE(beta_sesq)
  REAL, SAVE :: beta_meth                                           !! Methanol temperature dependency coefficient 
!$OMP THREADPRIVATE(beta_meth)
  REAL, SAVE :: beta_acet                                           !! Acetone temperature dependency coefficient
!$OMP THREADPRIVATE(beta_acet)
  REAL, SAVE :: beta_oxyVOC                                         !! Other oxygenated BVOC temperature dependency coefficient
!$OMP THREADPRIVATE(beta_oxyVOC)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: em_factor_ORVOC          !! ORVOC emissions factor
                                                                           !! @tex $(\mu gC.g^{-1}.h^{-1})$ @endtex
!$OMP THREADPRIVATE(em_factor_ORVOC)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: em_factor_OVOC           !! OVOC emissions factor
                                                                           !! @tex $(\mu gC.g^{-1}.h^{-1})$ @endtex
!$OMP THREADPRIVATE(em_factor_OVOC)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: em_factor_MBO            !! MBO emissions factor 
                                                                           !! @tex $(\mu gC.g^{-1}.h^{-1})$ @endtex
!$OMP THREADPRIVATE(em_factor_MBO)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: em_factor_methanol       !! Methanol emissions factor
                                                                           !! @tex $(\mu gC.g^{-1}.h^{-1})$ @endtex
!$OMP THREADPRIVATE(em_factor_methanol)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: em_factor_acetone        !! Acetone emissions factor 
                                                                           !! @tex $(\mu gC.g^{-1}.h^{-1})$ @endtex
!$OMP THREADPRIVATE(em_factor_acetone)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: em_factor_acetal         !! Acetaldehyde emissions factor
                                                                           !! @tex $(\mu gC.g^{-1}.h^{-1})$ @endtex
!$OMP THREADPRIVATE(em_factor_acetal)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: em_factor_formal         !! Formaldehyde emissions factor
                                                                           !! @tex $(\mu gC.g^{-1}.h^{-1})$ @endtex
!$OMP THREADPRIVATE(em_factor_formal)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: em_factor_acetic         !! Acetic Acid emissions factor 
                                                                           !! @tex $(\mu gC.g^{-1}.h^{-1})$ @endtex
!$OMP THREADPRIVATE(em_factor_acetic)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: em_factor_formic         !! Formic Acid emissions factor 
                                                                           !! @tex $(\mu gC.g^{-1}.h^{-1})$ @endtex
!$OMP THREADPRIVATE(em_factor_formic)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: em_factor_no_wet         !! NOx emissions factor soil emissions and 
                                                                           !! exponential dependancy factor for wet soils
                                                                           !! @tex $(ngN.m^{-2}.s^{-1})$ @endtex
!$OMP THREADPRIVATE(em_factor_no_wet)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: em_factor_no_dry         !! NOx emissions factor soil emissions and
                                                                           !! exponential dependancy factor for dry soils
                                                                           !! @tex $(ngN.m^{-2}.s^{-1})$ @endtex
!$OMP THREADPRIVATE(em_factor_no_dry)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: Larch                    !! Larcher 1991 SAI/LAI ratio (unitless)
!$OMP THREADPRIVATE(Larch)

  !
  ! INTERNAL PARAMETERS USED IN STOMATE_DATA
  !

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: lai_initmin   !! Initial lai for trees/grass 
                                                                !! @tex $(m^2.m^{-2})$ @endtex
!$OMP THREADPRIVATE(lai_initmin)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: bm_sapl   !! sapling biomass @tex $(gC.ind^{-1})$ @endtex
!$OMP THREADPRIVATE(bm_sapl)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: migrate       !! migration speed @tex $(m.year^{-1})$ @endtex
!$OMP THREADPRIVATE(migrate)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: maxdia        !! maximum stem diameter from which on crown area no longer 
                                                                !! increases (m)
!$OMP THREADPRIVATE(maxdia)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: cn_sapl       !! crown of tree when sapling  @tex $(m^2$)$ @endtex
!$OMP THREADPRIVATE(cn_sapl)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: leaf_timecst  !! time constant for leaf age discretisation (days)
!$OMP THREADPRIVATE(leaf_timecst)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: k_latosa_max       !! Maximum leaf-to-sapwood area ratio (unitless)
!$OMP THREADPRIVATE(k_latosa_max)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: k_latosa_min       !! Minimum leaf-to-sapwood area ratio (unitless)
!$OMP THREADPRIVATE(k_latosa_min)

 REAL, ALLOCATABLE, SAVE, DIMENSION(:,:) :: LC                 !! Lignin/C ratio of the different biomass pools and PFTs (unitless)
                                                                      !! based on CN from White et al. (2000) 
!$OMP THREADPRIVATE(LC)

 REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: LC_leaf              !! Lignin/C ratio of leaf pool (unitless)
                                                                      !! based on CN from White et al. (2000) 
!$OMP THREADPRIVATE(LC_leaf)

 REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: LC_sapabove          !! Lignin/C ratio of sapabove pool (unitless)
                                                                      !! based on CN from White et al. (2000) 
!$OMP THREADPRIVATE(LC_sapabove)

 REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: LC_sapbelow          !! Lignin/C ratio of sapbelow pool (unitless)
                                                                      !! based on CN from White et al. (2000) 
!$OMP THREADPRIVATE(LC_sapbelow)

 REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: LC_heartabove        !! Lignin/C ratio of heartabove pool (unitless)
                                                                      !! based on CN from White et al. (2000) 
!$OMP THREADPRIVATE(LC_heartabove)

 REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: LC_heartbelow        !! Lignin/C ratio of heartbelow pool (unitless)
                                                                      !! based on CN from White et al. (2000) 
!$OMP THREADPRIVATE(LC_heartbelow)

 REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: LC_fruit             !! Lignin/C ratio of fruit pool (unitless)
                                                                      !! based on CN from White et al. (2000) 
!$OMP THREADPRIVATE(LC_fruit)

 REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: LC_root              !! Lignin/C ratio of root pool (unitless)
                                                                      !! based on CN from White et al. (2000) 
!$OMP THREADPRIVATE(LC_root)

 REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: LC_carbres           !! Lignin/C ratio of carbres pool (unitless)
                                                                      !! based on CN from White et al. (2000) 
!$OMP THREADPRIVATE(LC_carbres)

 REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: LC_labile            !! Lignin/C ratio of labile pool (unitless)
                                                                      !! based on CN from White et al. (2000) 
!$OMP THREADPRIVATE(LC_labile)

 REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: decomp_factor        !! Multpliactive factor modifying 
                                                                      !! the standard decomposition factor for each SOM pool
!$OMP THREADPRIVATE(decomp_factor)



  !
  ! STAND STRUCTURE (stomate) 
  !
  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: pipe_density        !! Wood density in @tex $(gC.m^{-3})$ @endtex
!$OMP THREADPRIVATE(pipe_density)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: pipe_tune1          !! crown area = pipe_tune1*stem diameter**pipe_tune_exp_coeff
!$OMP THREADPRIVATE(pipe_tune1)
  
  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: pipe_tune2          !! height=pipe_tune2 * diameter**pipe_tune3
!$OMP THREADPRIVATE(pipe_tune2)      
  
  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: pipe_tune3          !! height=pipe_tune2 * diameter**pipe_tune3
!$OMP THREADPRIVATE(pipe_tune3)      
  
  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: pipe_tune4          !! ???needed for stem diameter
!$OMP THREADPRIVATE(pipe_tune4)      
  
  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: pipe_k1             !! ???
!$OMP THREADPRIVATE(pipe_k1)       
  
  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: pipe_tune_exp_coeff !! crown area = pipe_tune1*stem diameter**pipe_tune_exp_coeff
!$OMP THREADPRIVATE(pipe_tune_exp_coeff)      
  

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: tree_ff             !! Volume reduction factor from cylinder to real tree shape (inc.branches)
!$OMP THREADPRIVATE(tree_ff)
  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: mass_ratio_heart_sap!! mass ratio (heartwood+sapwood)/heartwood 
!$OMP THREADPRIVATE(mass_ratio_heart_sap)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: canopy_cover        !! Canopy cover - current values are guesses for testing
                                                                      !! could tune this variable to match MODIS albedo
!$OMP THREADPRIVATE(canopy_cover)
   
  INTEGER, ALLOCATABLE, SAVE, DIMENSION(:) :: nmaxtrees      !! Intial number of seedlings per hectare. Used 
                                                                    !! in prescribe to initialize the model and after 
                                                                    !! every clearcut
!$OMP THREADPRIVATE(nmaxtrees)
  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: height_init_min   !! The minimum height of a tree sapling when a forest 
                                                                    !! stand is established. Owing to the allometric 
                                                                    !! relationship this setting determines all
                                                                    !! biomass components of a newly establised stand
                                                                    !! @tex $(m)$ @endtex
!$OMP THREADPRIVATE(height_init_min)
  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: height_init_max   !! The maximum height of a tree sapling when a forest 
                                                                    !! stand is established. 
                                                                    !! @tex $(m)$ @endtex
!$OMP THREADPRIVATE(height_init_max)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: k_root              !! Fine root specific conductivity 
                                                                      !! @tex $(m^{3} kg^{-1} s^{-1} MPa^{-1})$ @endtex 
!$OMP THREADPRIVATE(k_root)
 
  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: k_sap               !! Sapwood specific conductivity
                                                                      !! @tex $(m^{3} kg^{-1} s^{-1} MPa^{-1})$ @endtex 
!$OMP THREADPRIVATE(k_sap)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: lai_to_height       !! Covert lai into vegetation height for grasses and crops
!$OMP THREADPRIVATE(lai_to_height)    

  REAL, ALLOCATABLE, SAVE, DIMENSION(:)    :: deleuze_a      !! intercept of the intra-tree competition within a stand
                                                                    !! based on the competion rule of Deleuze and Dhote 2004
                                                                    !! Used when n_circ > 6
!$OMP THREADPRIVATE(deleuze_a)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:)    :: deleuze_b      !! slope of the intra-tree competition within a stand
                                                                    !! based on the competion rule of Deleuze and Dhote 2004
                                                                    !! Used when n_circ > 6
!$OMP THREADPRIVATE(deleuze_b)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:)    :: deleuze_p_all  !! Percentile of the circumferences that receives photosynthates
                                                                    !! based on the competion rule of Deleuze and Dhote 2004
                                                                    !! Used when n_circ > 6 for FM1, FM2 and FM4 
!$OMP THREADPRIVATE(deleuze_p_all)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: m_dv                !! Parameter in the Deleuze & Dhote allocation rule that
                                                                      !! relaxes the cut-off imposed by ::sigma. Owing to m_relax 
                                                                      !! trees still grow a little when their ::circ is below 
                                                                      !! ::sigma
!$OMP THREADPRIVATE(m_dv) 

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: fruit_alloc         !! Fraction of biomass allocated to fruit production (0-1)

!$OMP THREADPRIVATE(fruit_alloc) 

 REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: labile_reserve       !! The size of the labile pool as a fraction of the
                                                                      !! weekly gpp (-). For example, 3 indicates that the 
                                                                      !! is 3 times the weekly gpp. 
!$OMP THREADPRIVATE(labile_reserve)

 REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: deciduous_reserve    !! Fraction of sapwood mass stored in the reserve pool of deciduous 
                                                                      !! trees during the growing season (unitless, 0-1) 

!$OMP THREADPRIVATE(deciduous_reserve)

 REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: evergreen_reserve    !! Fraction of sapwood mass stored in the reserve pool of evergreen 
                                                                      !! trees (unitless, 0-1) 

!$OMP THREADPRIVATE(evergreen_reserve)

 REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: senescense_reserve   !! Fraction of sapwood mass stored in the reserve pool of deciduous 
                                                                      !! trees during senescense(unitless, 0-1)

!$OMP THREADPRIVATE(senescense_reserve)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: fcn_wood            !! CN of wood for allocation, relative to leaf CN according
                                                                      !! to stich et al 2003
!$OMP THREADPRIVATE(fcn_wood)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: fcn_root            !! CN of roots for allocation, relative to leaf CN according
                                                                      !! to stich et al 2003
!$OMP THREADPRIVATE(fcn_root)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: fnp_wood            !! NP of wood for allocation, relative to leaf NP
                                                                      !! McGroddy et al. (2004) and Sardans et al. (2015) 
                                                                      
!$OMP THREADPRIVATE(fnp_wood)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: fnp_root            !! CN of roots for allocation, relative to leaf NP
                                                                      !! McGroddy et al. (2004) and Sardans et al. (2015) 
!$OMP THREADPRIVATE(fnp_root)




  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: branch_ratio      !! branches/total aboveground biomass ratio
                                                                    !! (cf carbofor for CITEPA inventory, these 
                                                                    !! Guerric, Lim 2004, Peischl 2007, 
                                                                    !! Schulp 2008: 15-30% slash after harvest,
                                                                    !! Zaehle 2007: 30% slash after harvest)
!$OMP THREADPRIVATE(branch_ratio)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: cn_leaf_init      !! CN of foliage for allocation, according to stich et al 2003
!$OMP THREADPRIVATE(cn_leaf_init)

  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: np_leaf_init      !! CN of foliage for allocation, 
                                                                    !! according to McGroddy et al. (2004) and Sardans et al. (2015)
!$OMP THREADPRIVATE(np_leaf_init)

END MODULE pft_parameters_var

