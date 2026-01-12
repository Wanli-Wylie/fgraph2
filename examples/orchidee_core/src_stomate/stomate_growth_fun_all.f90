!=================================================================================================================================
!
! CONTACT      : orchidee-help _at_ ipsl.jussieu.fr
!
! LICENCE      : IPSL (2006)
!                This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF         Plant growth and C-allocation among the biomass components (leaves, wood, roots, fruit, reserves, labile) 
!!               is calculated making use of functional allocation which combines the pipe model and allometric relationships 
!!               proposed by Sitch et al 2003 and adjusted by Zaehle et al 2010.
!!
!!\n DESCRIPTION: This module calculates three processes: (1) daily maintenance respiration based on the half-hourly
!!               respiration calculated in stomate_resp.f90, (2) the absolute allocation to the different biomass
!!               components based on functional allocation approach and (3) the allocatable biomass as the residual
!!               of GPP-Ra. Multiplication of the allocation fractions and allocatable biomass given the changes in
!!               biomass pools.
!!
!! RECENT CHANGE(S): Until 1.9.6 only one allocation scheme was available (now contained in stomate_grwoth_res_lim.f90). 
!!               This module consists of an alternative formalization of plant growth.
!!                             
!! REFERENCE(S)	: - Sitch, S., Smith, B., Prentice, I.C., Arneth, A., Bondeau, A., Cramer, W.,
!!               Kaplan, J.O., Levis, S., Lucht, W., Sykes, M.T., Thonicke, K., Venevsky, S. (2003), Evaluation of 
!!               ecosystem dynamics, plant geography and terrestrial carbon cycling in the LPJ Dynamic Global Vegetation 
!!               Model, Global Change Biology, 9, 161-185.\n 
!!               - Zaehle, S. and Friend, A.D. (2010), Carbon and nitrogen cycle dynamics in the O-CN land surface model: 1. 
!!               Model description, site-scale evaluation, and sensitivity to parameter estimates, Global Biogeochemical 
!!               Cycles, 24, GB1005.\n
!!               - Magnani F., Mencuccini M. & Grace J. 2000. Age-related decline in stand productivity: the role of 
!!               structural acclimation under hydraulic constraints Plant, Cell and Environment 23, 251–263.\n
!!               - Bloom A.J., Chapin F.S. & Mooney H.A. (1985) Resource limitation in plants. An economic analogy.  
!!               Annual Review Ecology Systematics 16, 363–392.\n
!!               - Case K.E. & Fair R.C. (1989) Principles of Economics. Prentice Hall, London.\n
!!               - McDowell, N., Barnard, H., Bond, B.J., Hinckley, T., Hubbard, R.M., Ishii, H., Köstner, B., 
!!               Magnani, F. Marshall, J.D., Meinzer, F.C., Phillips, N., Ryan, M.G., Whitehead D. 2002. The 
!!               relationship between tree height and leaf area: sapwood area ratio. Oecologia, 132:12–20.\n
!!               - Novick, K., Oren, R., Stoy, P., Juang, F.-Y., Siqueira, M., Katul, G. 2009. The relationship between
!!               reference canopy conductance and simplified hydraulic architecture. Advances in water resources 32,
!!               809-819.     
!!
!! SVN          :
!! $HeadURL$
!! $Date$
!! $Revision$
!! \n
!_ ===============================================================================================================================
MODULE stomate_growth_fun_all

  ! Modules used:
  USE ioipsl_para
  USE xios_orchidee
  USE pft_parameters
  USE stomate_data
  USE constantes
  USE constantes_soil
  USE stomate_phosphorus
  USE function_library,    ONLY: wood_to_qmdia, wood_to_qmheight, wood_to_ba_eff, biomass_to_lai,&
                                 calculate_c0_alloc, wood_to_volume, wood_to_ba

  IMPLICIT NONE

  ! Private & public routines

  PRIVATE
  PUBLIC growth_fun_all_clear, growth_fun_all

 ! Variables shared by all subroutines in this module

  LOGICAL, SAVE                                             :: firstcall = .TRUE.  !! Is this the first call? (true/false)
 
!!$  !+++TEMP+++
  INTEGER, SAVE                                             :: istep = 0

CONTAINS



!! ================================================================================================================================
!! SUBROUTINE   : growth_fun_all_clear
!!
!>\BRIEF          Set the flag ::firstcall to .TRUE. and as such activate section 
!! 1.1 of the subroutine alloc (see below).\n
!!
!_ ================================================================================================================================

  SUBROUTINE growth_fun_all_clear
    firstcall = .TRUE.
  END SUBROUTINE growth_fun_all_clear



!! ================================================================================================================================
!! SUBROUTINE 	: growth_fun_all
!!
!>\BRIEF          Allocate net primary production (= photosynthesis
!!                minus autothrophic respiration) to: labile carbon pool carbon reserves, aboveground sapwood,
!!                belowground sapwood, root, fruits and leaves following the pipe model and allometric constraints.
!!
!! DESCRIPTION  : Total maintenance respiration for the whole plant is calculated by summing maintenance 
!!                respiration of the different plant compartments. Maintenance respiration is subtracted 
!!                from whole-plant allocatable biomass (up to a maximum fraction of the total allocatable biomass). 
!!                Growth respiration is then calculated as a prescribed (0.75) fraction of the allocatable
!!                biomass. Subsequently NPP is calculated by substracting total autotrophic  respiration from GPP 
!!                i.e. NPP = GPP - maintenance resp - growth resp.
!!
!!                The pipe model assumes that one unit of leaf mass requires a proportional amount of sapwood to 
!!                transport water from the roots to the leaves. Also a proportional fraction of roots is needed to 
!!                take up the water from the soil. The proportional amounts between leaves, sapwood and roots are 
!!                given by allocation factors. These allocation factors are PFT specific and depend on a parameter 
!!                quantifying the leaf to sapwood area (::k_latosa_target), the specific laeaf area (::sla), wood 
!!                density (::pipe_density) and a scaling parameter between leaf and root mass.
!! 
!!                Lai is optimised for mean annual radiation use efficiency and the C cost for producing the 
!!                canopy. The cost-benefit ratio is optimised when the marginal gain / marginal cost = 1 lai target 
!!                is used to calculate whether the reserves are used. This approach allows plants to get out of 
!!                senescence and to start developping a canopy in early spring. 
!!                  
!!                As soon as a canopy has emerged, C (b_inc_tot) becomes available at the stand level through 
!!                photosynthesis and, C is allocated at the tree level (b_inc) following both the pipe model and 
!!                allometric constraints. Mass conservation requires:
!!                (1) Cs_inc + Cr_inc + Cl_inc = b_inc
!!                (2) sum(b_inc) = b_inc_tot
!!
!!                Wood allocation depends on tree basal area following the rule of Deleuze & Dhote
!!                delta_ba = gammas*(circ - m*sigmas + sqrt((m*sigmas + circ).^2 - (4*sigmas*circ)))/2
!!                (3) <=> delta_ba = circ_class_dba*gammas
!!                Where circ_class_dba = (circ - m*sigmas + sqrt((m*sigmas + circ).^2 - (4*sigmas*circ)))/2
!! 
!!                Allometric relationships
!!                height = pipe_tune2*(dia.^pipe_tune3)
!!                Re-write this relationship as a function of ba 
!!                (4) height = pipe_tune2 * (4/pi*ba)^(pipe_tune3/2)
!!                (5a) Cl/Cs = KF/height for trees
!!                (5b) Cs = Cl / KF
!!                (6) Cl = Cr * LF
!!
!!                Use a linear approximation to avoid iterations. Given that allocation is calculated daily, a 
!!                local lineair assumption is fair. Eq (4) can thus be rewritten as:
!!                s = step/(pipe_tune2*(4/pi*(ba+step)).^(pipe_tune3/2)-pipe_tune2*(4/pi*ba).^(pipe_tune3/2))
!!                Where step is a small but realistic (for the time step) change in ba
!!                (7)  <=> delta_height = delta_ba/s
!!
!!                Calculate Cs_inc from allometric relationships
!!                Cs_inc = tree_ff*pipe_density*(ba+delta_ba)*(height+delta_height) - Cs - Ch         
!!                Cs_inc = tree_ff*pipe_density*(ba+delta_ba)*(height+delta_ba/s) - Cs - Ch
!!                (8)  <=> Cs_inc = tree_ff*pipe_density*(ba+a*gammas)*(height+(a/s*gammas)) - Cs - Ch
!!
!!                Rewrite (5) as 
!!                Cl_inc = KF*(Cs_inc+Cs)/(height+delta_height) - Cl
!!                Substitute (7) in (4) and solve for Cl_inc
!!                Cl_inc = KF*(tree_ff*pipe_density*(ba+circ_class_dba*gammas)*(height+(circ_class_dba/s*gammas)) - Ch)/ & 
!!                   (height+(circ_class_dba/s*gammas)) - Cl  
!!                (9)  <=> Cl_inc = KF*tree_ff*pipe_density*(ba+circ_class_dba*gammas) - &
!!                            (KF*Ch)/(height+(circ_class_dba/s*gammas)) - Cl
!!
!!                Rewrite (6) as 
!!                Cr_inc = (Cl_inc+Cl)/LF - Cr
!!                Substitute (9) in (6)
!!                (10)  <=> Cr_inc = KF/LF*tree_ff*pipe_density*(ba+circ_class_dba*gammas) - &
!!                            (KF*Ch/LF)/(height+(circ_class_dba/s*gammas)) - Cr
!!
!!                Depending on the specific case that needs to be solved equations (1) takes one of the following forms:
!!                (a) b_inc = Cl_inc + Cr_inc + Cs_inc, (b) b_inc = Cl_inc + Cr_inc, (c) b_inc = Cl_inc + Cs_inc or 
!!                (d) b_inc = Cr_inc + Cs_inc. One of these alternative forms of eq. 1 are then combined with 
!!                eqs 8, 9 and 10 and solved for gammas. The details for the solution of these four cases are given in the
!!                code. Once gammas is know, eqs 6 - 10 are used to calculate the allocation to leaves (Cl_inc), 
!!                roots (Cr_inc) and sapwood (Cs_inc).
!!
!!                Because of turnover, biomass pools are not all the time in balance according to rules prescribed
!!                by the pipe model. To test whether biomass pools are balanced, the target biomasses are calculated 
!!                and balance is restored whenever needed up to the level that the biomass pools for leaves, sapwood
!!                and roots are balanced according to the pipe model. Once the balance is restored C is allocated to
!!                fruits, leaves, sapwood and roots by making use of the pipe model (below this called ordinary 
!!                allocation).
!!
!!                Although strictly speaking allocation factors are not necessary in this scheme (Cl_inc could simply 
!!                be added to biomass(:,:,ileaf,icarbon), Cr_inc to biomass(:,:,iroot,icarbon), etc.), they are 
!!                nevertheless calculated because using allocation factors facilitates comparison to the resource 
!!                limited allocation scheme (stomate_growth_res_lim.f90) and it comes in handy for model-data comparison.
!!
!!                Effective basal area, height and circumferences are use in the allocation scheme because their 
!!                calculations make use of the total (above and belowground) biomass. In forestry the same measures
!!                exist (and they are also calculated in ORCHIDEE) but only account for the aboverground biomass.
!! 
!!
!! RECENT CHANGE(S): - The code by Sonke Zaehle made use of ::Cl_target that was derived from ::lai_target which in turn 
!!                was a function of ::rue_longterm. Cl_target was then used as a threshold value to decide whether there 
!!                was only phenological growth (just leaves and roots) or whether there was full allometric growth to the 
!!                leaves, roots and sapwood. This approach was inconsistent with the pipe model because full allometric 
!!                growth can only occur if all three biomass pools are in balance. ::lai_target is no longer used as a 
!!                criterion to switch between phenological and full allometric growth. Its use is now restricted to trigger
!!                the use of reserves in spring. 
!!
!! RECENT CHANGE(S): - Reserve pools are included again and the several modification to the carbon dynamics were made 
!!                     DSG
!!
!! MAIN OUTPUT VARIABLE(S): ::npp and :: biomass. Seven different biomass compartments (leaves, roots, above and 
!!                belowground wood, carbohydrate reserves, labile and fruits).
!!
!! REFERENCE(S)	:- Sitch, S., Smith, B., Prentice, I.C., Arneth, A., Bondeau, A., Cramer, W.,
!!                Kaplan, J.O., Levis, S., Lucht, W., Sykes, M.T., Thonicke, K., Venevsky, S. (2003), Evaluation of 
!!                ecosystem dynamics, plant geography and terrestrial carbon cycling in the LPJ Dynamic Global Vegetation 
!!                Model, Global Change Biology, 9, 161-185.\n 
!!                - Zaehle, S. and Friend, A.D. (2010), Carbon and nitrogen cycle dynamics in the O-CN land surface model: 1. 
!!                Model description, site-scale evaluation, and sensitivity to parameter estimates, Global Biogeochemical 
!!                Cycles, 24, GB1005.\n
!!                - Magnani F., Mencuccini M. & Grace J. 2000. Age-related decline in stand productivity: the role of 
!!                structural acclimation under hydraulic constraints Plant, Cell and Environment 23, 251–263.
!!                - Bloom A.J., Chapin F.S. & Mooney H.A. (1985) Resource limitation in plants. An economic analogy. Annual 
!!                Review Ecology Systematics 16, 363–392.
!!                - Case K.E. & Fair R.C. (1989) Principles of Economics. Prentice Hall, London.
!!                - McDowell, N., Barnard, H., Bond, B.J., Hinckley, T., Hubbard, R.M., Ishii, H., Köstner, B., 
!!                Magnani, F. Marshall, J.D., Meinzer, F.C., Phillips, N., Ryan, M.G., Whitehead D. 2002. The 
!!                relationship between tree height and leaf area: sapwood area ratio. Oecologia, 132:12–20
!!                - Novick, K., Oren, R., Stoy, P., Juang, F.-Y., Siqueira, M., Katul, G. 2009. The relationship between
!!                reference canopy conductance and simplified hydraulic architecture. Advances in water resources 32,
!!                809-819.   
!!
!! +++++++++++++++++++++
!! MAKE A NEW FLOW CHART
!! +++++++++++++++++++++
!! FLOWCHART    : 
!!
!_ ================================================================================================================================

  SUBROUTINE growth_fun_all (npts, dt, veget_max, veget, PFTpresent, &
       senescence, when_growthinit, t2m, &
       nstress_season, pstress_season, moiavail_growingseason, &
       gpp_daily, resp_maint_part, resp_maint, &
!DSGlabfac
       gpp_week, &
!DSGlabfac
       resp_growth, npp, biomass, age, &
       leaf_age, leaf_frac, use_reserve, &
       ind, rue_longterm, KF,  k_latosa_adapt, &
       cn_leaf_avg_season, n_uptake_daily, N_support, &
       np_leaf_avg_season, p_uptake_daily, P_support, &
       lai_target, lai_target_longterm, &
       record_reserve_grass, &       
       sla_calc, sla_age1, BNF_clover_daily, &
       GRM_devstage, when_growthinit_cut,days_senescence)

 !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER, INTENT(in)                        :: npts                   !! Domain size - number of grid cells 
                                                                                !! (unitless)
    REAL, INTENT(in)                           :: dt                     !! Time step of the simulations for stomate
                                                                                !! (days)
    REAL, DIMENSION(:), INTENT(in)             :: t2m                    !! Temperature at 2 meter (K)
    REAL, DIMENSION(:,:), INTENT(in)           :: ind                    !! Number of individuals at the stand level
                                                                                !! @tex $(m^{-2})$ @endtex     
    REAL, DIMENSION(:,:), INTENT(in)           :: veget_max              !! PFT "Maximal" coverage fraction of a PFT 
                                                                                !! (= ind*cn_ind) 
                                                                                !! @tex $(m^2 m^{-2})$ @endtex
    REAL, DIMENSION(:,:), INTENT(in)           :: veget                  !! Fraction of vegetation type including 
                                                                                !! non-biological fraction (unitless)   
    REAL, DIMENSION(:,:), INTENT(in)           :: when_growthinit        !! Days since beginning of growing season 
                                                                                !! (days)
    REAL, DIMENSION(:,:), INTENT(in)           :: rue_longterm           !! Longterm "radiation use efficicency"
                                                                                !! calculated as the ratio of GPP over 
                                                                                !! the fraction of radiation absorbed
                                                                                !! by the canopy 
                                                                                !! @tex $(gC.m^{-2}day^{-1})$ @endtex
    REAL, DIMENSION(:,:,:), INTENT(in)         :: resp_maint_part        !! Maintenance respiration of different  
                                                                                !! plant parts
                                                                                !! @tex $(gC.m^{-2}dt^{-1})$ @endtex
    LOGICAL, DIMENSION(:,:), INTENT(in)               :: senescence             !! Is the PFT senescent?  - only for 
                                                                                !! deciduous trees (true/false)
    LOGICAL, DIMENSION(:,:), INTENT(in)               :: PFTpresent             !! PFT exists (true/false)    
    REAL, DIMENSION(npts,nvm), INTENT(in)      :: nstress_season         !! N-related seasonal stress (used for allocation)    
    REAL, DIMENSION(npts,nvm), INTENT(in)      :: pstress_season         !! P-related seasonal stress (used for allocation)    
    REAL, DIMENSION(npts,nvm), INTENT(in)      :: moiavail_growingseason !! mean growing season moisture availability 
    REAL, DIMENSION(npts,nvm), INTENT(in)     :: lai_target_longterm      !! Target LAI used to calcuate the maximum reserve pool
   
    !! 0.2 Output variables

    REAL, DIMENSION(:,:), INTENT(out)          :: resp_maint             !! PFT maintenance respiration 
                                                                                !! @tex $(gC.m^{-2}dt^{-1})$ @endtex    
    REAL, DIMENSION(:,:), INTENT(out)          :: resp_growth            !! PFT growth respiration 
                                                                                !! @tex $(gC.m^{-2}dt^{-1})$ @endtex
    REAL, DIMENSION(:,:), INTENT(out)          :: npp                    !! PFT net primary productivity 
                                                                                !! @tex $(gC.m^{-2}dt^{-1})$ @endtex

    REAL, DIMENSION(npts,nvm), INTENT(out)     :: lai_target             !! Target LAI @tex $(m^{2}m^{-2})$ @endtex


    !! 0.3 Modified variables

    REAL, DIMENSION(:,:), INTENT(inout)        :: gpp_daily              !! PFT gross primary productivity 
                                                                                !! @tex $(gC.m^{-2}dt^{-1})$ @endtex
    REAL, DIMENSION(:,:), INTENT(in)           :: gpp_week               !! Weekly PFT gross primary productivity 
                                                                                !! @tex $(gC.m^{-2}dt^{-1})$ @endtex
    REAL, DIMENSION(:,:), INTENT(inout)        :: use_reserve            !! Flag to use the reserves to support
                                                                                !! phenological growth (0 or 1)
    REAL, DIMENSION(:,:), INTENT(inout)        :: age                    !! PFT age (days)      
    REAL, DIMENSION(:,:,:,:), INTENT(inout)    :: biomass                !! PFT total biomass 
                                                                                !! @tex $(gC m^{-2})$ @endtex
!JC separate tissue age
!    REAL, DIMENSION(:,:,:), INTENT(inout)      :: leaf_age               !! PFT age of different leaf classes 
!                                                                                !! (days)
!    REAL, DIMENSION(:,:,:), INTENT(inout)      :: leaf_frac              !! PFT fraction of leaves in leaf age class
!                                                                                !! (0-1, unitless)
    REAL, DIMENSION(npts,nvm,nparts,nleafages), INTENT(inout)      :: leaf_age
    REAL, DIMENSION(npts,nvm,nparts,nleafages), INTENT(inout)      :: leaf_frac
!End JC separate tissue age
    REAL, DIMENSION(:,:), INTENT(inout)        :: KF                     !! Scaling factor to convert sapwood mass
                                                                                !! into leaf mass (m)
    REAL, DIMENSION(:,:), INTENT(in)           :: k_latosa_adapt         !! Leaf to sapwood area adapted for long 
                                                                                !! term water stress (m)
    REAL, DIMENSION(:,:), INTENT(in)           :: cn_leaf_avg_season     !! Average leaf nitrogen concentration (C:N) of the growing season  
                                                                                !! (gC/gN) 
    REAL, DIMENSION(npts,nvm,nionspec), INTENT(in) :: n_uptake_daily     !! Uptake of soil N by plants  
    REAL, DIMENSION(npts,nvm), INTENT(inout)       :: N_support          !! Nitrogen which is added to the ecosystem to support vegetation growth (only used when IMPOSE_CN .EQ. true) (gC/m2/day)

    REAL, DIMENSION(:,:), INTENT(in)               :: np_leaf_avg_season !! Average leaf nitrogen concentration (N:P) of the growing season  
                                                                                !! (gN/gP) 
    REAL, DIMENSION(npts,nvm), INTENT(in)          :: p_uptake_daily     !! Uptake of soil N by plants  
    REAL, DIMENSION(npts,nvm), INTENT(inout)       :: P_support          !! Phosphorus which is added to the ecosystem to support vegetation growth (only used when IMPOSE_NP .EQ. true) (gP/m2/day)
    REAL, DIMENSION(npts,nvm), INTENT(inout)       :: record_reserve_grass  !! record maximum grass reserve pool                                                                              
    ! N fertilization limitation factor for grassland Vcmax and SLA
    REAL, DIMENSION(npts,nvm), INTENT(inout)       :: sla_calc           !! leaf age-related SLA
    REAL, DIMENSION(npts,nvm), INTENT(inout)       :: sla_age1           !! maximum SLA for youngest leaves
    REAL, DIMENSION(npts,nvm), INTENT(in)          :: BNF_clover_daily   !! daily clover BNF in grass-legume mixtures
    REAL, DIMENSION(npts,nvm), INTENT(in)          :: GRM_devstage       !! development stage of grasses
    REAL, DIMENSION(npts,nvm), INTENT(in)          :: when_growthinit_cut!! growing days after grass harvest
    REAL, DIMENSION(npts,nvm), INTENT(in)          :: days_senescence    !! day counter after senescence is detected
    
    !! 0.4 Local variables
    CHARACTER(30)                                     :: var_name               !! To store variable names for I/O
    REAL, DIMENSION(npts,nvm)                  :: c0_alloc               !! Root to sapwood tradeoff parameter
    LOGICAL                                           :: grow_wood=.TRUE.       !! Flag to grow wood
    INTEGER                                    :: ipts,j,k,l,m,icirc,imed!! Indeces(unitless)
    INTEGER                                    :: ipar, iele, imbc       !! Indeces(unitless)
    INTEGER                                    :: ilev                   !! Indeces(unitless)
    REAL                                       :: frac                   !! No idea??
    REAL                                       :: a,b,c                  !! Temporary variables to solve a
                                                                                !! quadratic equation (unitless)
    ! Stand level
    REAL                                       :: gtemp                  !! Turnover coefficient of labile C pool
                                                                                !! (0-1??, units??)
    REAL                                       :: reserve_pool           !! Intentional size of the reserve pool
                                                                                !! @tex $(gC.m^{-2})$ @endtex
    REAL                                       :: labile_pool            !! Intentional size of the labile apool
                                                                                !! @tex $(gC.m^{-2})$ @endtex
    REAL                                       :: reserve_scal           !! Protection of the reserve against 
                                                                                !! overuse (unitless)
    REAL                                       :: use_lab                !! Availability of labile biomass 
                                                                                !! @tex $(gC.m^{-2})$ @endtex
    REAL                                       :: use_res                !! Availability of resource biomass 
                                                                                !! @tex $(gC.m^{-2})$ @endtex
    REAL                                       :: use_max                !! Maximum use of labile and resource pool 
                                                                                !! @tex $(gC.m^{-2})$ @endtex

    REAL, DIMENSION(nvm)                       :: root_reserve           !! fraction of max root biomass which are covered by carbon reserve
    REAL                                       :: leaf_meanage           !! Mean age of the leaves (days?)
    REAL                                       :: reserve_time           !! Maximum number of days during which 
                                                                                !! carbohydrate reserve may be used (days)
    REAL                                       :: b_inc_tot              !! Carbon that needs to allocated in the
                                                                                !! fixed number of trees (gC)
    REAL                                       :: b_inc_temp             !! Temporary b_inc at the stand-level
                                                                                !! @tex $(gC.plant^{-1})$ @endtex
    REAL, DIMENSION(npts,nvm)                  :: scal                   !! Scaling factor between average 
                                                                                !! individual and individual plant
                                                                                !! @tex $(plant.m^{-2})$ @endtex
    REAL                                       :: total_inc              !! Total biomass increase
                                                                                !! @tex $(gC.plant^{-1})$ @endtex
    REAL                                       :: KF_old                 !! Scaling factor to convert sapwood mass
                                                                                !! into leaf mass (m) at the previous
                                                                                !! time step
    REAL, DIMENSION(nvm)                       :: lai_happy              !! Lai threshold below which carbohydrate 
                                                                                !! reserve may be used 
                                                                                !! @tex $(m^2 m^{-2})$ @endtex
    REAL, DIMENSION(nvm)                       :: deleuze_p              !! Percentile of trees that will receive 
                                                                                !! photosynthates. The proxy for intra stand
                                                                                !! competition. Depends on the management
                                                                                !! strategy when ncirc < 6
    REAL, DIMENSION(npts)                      :: tl                     !! Long term annual mean temperature (C)
    REAL, DIMENSION(npts)                      :: bm_add                 !! Biomass increase
                                                                                !! @tex $(gC.m^{-2})$ @endtex
    REAL, DIMENSION(npts)                      :: bm_new                 !! New biomass @tex $(gC.m^{-2})$ @endtex
    REAL                                       :: alloc_sap_above        !! Fraction allocated to sapwood above 
                                                                                !! ground
    REAL, DIMENSION(npts,nvm)                  :: residual               !! Copy of b_inc_tot after all C has been
                                                                                !! allocated @tex $(gC.m^{-2})$ @endtex
                                                                                !! if all went well the value should be zero

    REAL, DIMENSION(npts,nvm)                  :: ltor                   !! Leaf to root ratio (unitless)   
    REAL, DIMENSION(npts,nvm)                  :: lstress_fac            !! Light stress factor, based on total
                                                                                !! transmitted light (unitless, 0-1)
    REAL, DIMENSION(npts,nvm)                  :: k_latosa               !! Target leaf to sapwood area ratio
    REAL, DIMENSION(npts,nvm)                  :: LF                     !! Scaling factor to convert sapwood mass
                                                                                !! into root mass (unitless)
    REAL, DIMENSION(npts,nvm)                  :: bm_alloc_tot           !! Allocatable biomass for the whole plant
!JC separate tissue age
! add nparts to lm_old and leaf_mass_young
    REAL, DIMENSION(npts,nvm,nparts)                  :: lm_old
    REAL, DIMENSION(npts,nvm,nparts)                  :: leaf_mass_young
!End JC separate tissue age
    REAL, DIMENSION(npts,nvm)                  :: lai                    !! PFT leaf area index 
                                                                                !! @tex $(m^2 m^{-2})$ @endtex
    REAL, DIMENSION(npts,nvm)                  :: qm_dia                 !! Quadratic mean diameter of the stand (m)
    REAL, DIMENSION(npts,nvm)                  :: qm_height              !! Height of a tree with the quadratic mean
                                                                                !! diameter (m)
    REAL, DIMENSION(npts,nvm)                  :: ba                     !! Basal area. variable for histwrite_p only (m2)
    REAL, DIMENSION(npts,nvm)                  :: wood_volume            !! wood_volume (m3 m-2)
    REAL, DIMENSION(npts,nvm,nparts)           :: f_alloc                !! PFT fraction of NPP that is allocated to
                                                                                !! the different components (0-1, unitless)
    REAL, DIMENSION(npts,nvm,nparts,nelements) :: bm_alloc               !! PFT biomass increase, i.e. NPP per plant  
                                                                                !! part @tex $(gC.m^{-2}dt^{-1})$ @endtex

    REAL, DIMENSION(npts,nvm,2)                :: bm_res_lab             !! to be able to keep track of changes in labile and reserve C pools
                                                                                !! I introduce this variable instead of keeping track of 
                                                                                !! all the deficit,residual calculations
    ! Tree level	
    REAL, DIMENSION(ncirc)                     :: step                   !! Temporary variables to solve a
                                                                                !! quadratic equation (unitless)
    REAL, DIMENSION(ncirc)                     :: s                      !! tree-level linear relationship between
                                                                                !! basal area and height. This variable is
                                                                                !! used to linearize the allocation scheme
    REAL, DIMENSION(ncirc)                     :: Cs_inc_est             !! Initial value estimate for Cs_inc. The 
                                                                                !! value is used to linearize the ba~height
                                                                                !! relationship 
                                                                                !! @tex $(gC.plant^{-1})$ @endtex
    REAL, DIMENSION(ncirc)                     :: Cl                     !! Individual plant, leaf compartment
                                                                                !! @tex $(gC.plant^{-1})$ @endtex
    REAL, DIMENSION(ncirc)                     :: Cr                     !! Individual plant, root compartment
                                                                                !! @tex $(gC.plant^{-1})$ @endtex
    REAL, DIMENSION(ncirc)                     :: Cs                     !! Individual plant, sapwood compartment
                                                                                !! @tex $(gC.plant^{-1})$ @endtex
    REAL, DIMENSION(ncirc)                     :: Ch                     !! Individual plant, heartwood compartment
                                                                                !! @tex $(gC.plant^{-1})$ @endtex
    REAL, DIMENSION(ncirc)                     :: Cl_inc                 !! Individual plant increase in leaf 
                                                                                !! compartment
                                                                                !! @tex $(gC.plant^{-1})$ @endtex 
    REAL, DIMENSION(ncirc)                     :: Cr_inc                 !! Individual plant increase in root
                                                                                !! compartment
                                                                                !! @tex $(gC.plant^{-1})$ @endtex 
    REAL, DIMENSION(ncirc)                     :: Cs_inc                 !! Individual plant increase in sapwood
                                                                                !! compartment
                                                                                !! @tex $(gC.plant^{-1})$ @endtex 
    REAL, DIMENSION(ncirc)                     :: Cf_inc                 !! Individual plant increase in fruit
                                                                                !! compartment
                                                                                !! @tex $(gC.plant^{-1})$ @endtex
    REAL, DIMENSION(ncirc)                     :: Cl_incp                !! Phenology related individual plant 
                                                                                !! increase in leaf compartment
                                                                                !! @tex $(gC.plant^{-1})$ @endtex
    REAL, DIMENSION(ncirc)                     :: Cr_incp                !! Phenology related individual plant 
                                                                                !! increase in leaf compartment
                                                                                !! @tex $(gC.plant^{-1})$ @endtex
    REAL, DIMENSION(ncirc)                     :: Cs_incp                !! Phenology related individual plant 
                                                                                !! increase in sapwood compartment
                                                                                !! @tex $(gC.plant^{-1})$ @endtex
    REAL, DIMENSION(ncirc)                     :: Cl_target              !! Individual plant maximum leaf mass given 
                                                                                !! its current sapwood mass
                                                                                !! @tex $(gC.plant^{-1})$ @endtex
    REAL, DIMENSION(ncirc)                     :: Cr_target              !! Individual plant maximum root mass given 
                                                                                !! its current sapwood mass
                                                                                !! @tex $(gC.plant^{-1})$ @endtex
    REAL, DIMENSION(ncirc)                     :: Cs_target              !! Individual plant maximum sapwood mass 
                                                                                !! given its current leaf mass or root mass
                                                                                !! @tex $(gC.plant^{-1})$ @endtex     
    REAL, DIMENSION(ncirc)                     :: delta_ba               !! Change in basal area for a unit
                                                                                !! investment into sapwood mass (m)
    REAL, DIMENSION(ncirc)                     :: delta_height           !! Change in height for a unit
                                                                                !! investment into sapwood mass (m)
    REAL, DIMENSION(ncirc)                     :: circ_class_ba          !! Basal area (forestry definition) of the model 
                                                                                !! tree in each circ class 
                                                                                !! @tex $(m^{2} m^{-2})$ @endtex 
    REAL, DIMENSION(ncirc)                     :: circ_class_ba_eff      !! Effective basal area of the model tree in each
                                                                                !! circ class @tex $(m^{2} m^{-2})$ @endtex 
    REAL, DIMENSION(ncirc)                     :: circ_class_dba         !! Share of an individual tree in delta_ba
                                                                                !! thus, circ_class_dba*gammas = delta_ba
    REAL, DIMENSION(ncirc)                     :: circ_class_height_eff  !! Effective tree height calculated from allometric 
                                                                                !! relationships (m)
    REAL, DIMENSION(ncirc)                     :: circ_class_circ_eff    !! Effective circumference of individual trees (m)
    REAL                                       :: woody_biomass          !! Woody biomass. Temporary variable to 
                                                                                !! calculate wood volume (gC m-2)
    REAL                                       :: temp_share             !! Temporary variable to store the share
                                                                                !! of biomass of each circumference class
                                                                                !! to the total biomass            
    REAL                                       :: temp_class_biomass     !! Biomass across parts for a single circ
                                                                                !! class @tex $(gC m^{-2})$ @endtex
    REAL                                       :: temp_total_biomass     !! Biomass across parts and circ classes
                                                                                !! @tex $(gC m^{-2})$ @endtex
    REAL, DIMENSION(npts,nvm,ncirc)            :: store_delta_ba         !! Store delta_ba in this variable before writing
                                                                                !! to the output file (m). Adding this variable
                                                                                !! was faster than changing the dimensions
                                                                                !! of delta_ba which would have been the same
    REAL, DIMENSION(npts,nvm,ncirc)            :: store_circ_class_ba    !! Store circ_class_ba in this variable before
                                                                                !! writing to the output file (m). Adding this
                                                                                !! variable was faster than changing the 
                                                                                !! dimensions of circ_class_ba_ba which would 
                                                                                !! have been the same
    REAL, DIMENSION(npts,nvm,nmbcomp,nelements):: check_intern           !! Contains the components of the internal
                                                                                !! mass balance chech for this routine
                                                                                !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    REAL, DIMENSION(npts,nvm,nelements)        :: closure_intern         !! Check closure of internal mass balance
                                                                                !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    REAL, DIMENSION(npts,nvm,nelements)        :: pool_start, pool_end   !! Start and end pool of this routine 
                                                                                !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    REAL                                       :: median_circ            !! Median circumference (m)
    REAL                                       :: deficit                !! Carbon that needs to be respired in 
                                                                                !! excess of todays gpp  
                                                                                !! @tex $(gC.m^{-2}dt^{-1})$ @endtex
    REAL                                       :: excess                 !! Carbon that needs to be re-allocated
                                                                                !! after the needs of the reserve and 
                                                                                !! labile pool are satisfied  
                                                                                !! @tex $(gC.m^{-2}dt^{-1})$ @endtex
    REAL                                       :: shortage               !! Shortage in the reserves that needs to 
                                                                                !! be re-allocated after to minimise the
                                                                                !! tension between required and available
                                                                                !! reserves
                                                                                !! @tex $(gC.m^{-2}dt^{-1})$ @endtex
    
    INTEGER                                           :: i,tempi
!   (temp variables for impose intraseasonal LAI dynamic)   
    INTEGER                                           :: month_id               !! index of month 
    REAL                                       :: ratio_move             !! tmperal variable to move the allocatable carbon
                                                                                !! from leaf to sapwood
    REAL                                       :: impose_lai             !! get value of impose maximun LAI
    REAL, DIMENSION(13)                        :: lai_scale              !! monthly lai scaling facter    
    REAL                                       :: daily_lai              !! Daily LAI value interpolated by impose lai & lai_scale 
    CHARACTER(len=256)                                :: temp_text              !! dummy text variable exchange

    REAL, DIMENSION(npts,nvm,ncirc)                  :: circ_class_n           !! Number of individuals in the single circ class
    REAL, DIMENSION(npts,nvm,ncirc,nparts,nelements) :: circ_class_biomass     !! Biomass components of the model tree  
                                                                                !! within the single circumference class
                                                                                !! @tex $(g C ind^{-1})$ @endtex
    REAL, DIMENSION(npts,nvm)                       :: sigma                  !! Threshold for indivudal tree growth in 
                                                                                !! the equation of Deleuze & Dhote (2004)(m).
                                                                                !! Trees whose circumference is smaller than 
                                                                                !! sigma don't grow much
    REAL, DIMENSION(npts,nvm)                      :: gammas                 !! Slope for individual tree growth in the 
                                                                                !! equation of Deleuze & Dhote (2004) (m)
    REAL, DIMENSION(npts,nvm)                      :: lab_fac                !! Activity of labile pool factor (units??)

    ! Nitrogen cycle 
    REAL, DIMENSION(npts,nvm)                  :: n_alloc_tot            !! nitrogen growth (gN/m2/dt) 
    REAL , DIMENSION(npts,nvm)                 :: cn_leaf                !! nitrogen concentration in leaves (gC/gN) 
    REAL , DIMENSION(npts,nvm)                 :: transloc_N             !! Transloc variables Nitrogen
    REAL, DIMENSION(npts,nelements)            :: alloc_c,alloc_d,alloc_e                  !! allocation coefficients of nitrogen to leaves, roots and wood 
    REAL, DIMENSION(npts)                      :: sum_sap,sum_oth        !! carbon growth of wood and root+fruits (gC/m2) 
    REAL                                       :: costf_N                !! nitrogen cost of a unit carbon growth given current C partitioning and nitrogen concentration 
    REAL                                       :: deltacn,deltacnmax     !! (maximum) change in leaf nitrogen concentration 
    REAL                                       :: n_avail                !! nitrogen available to growth (dummy) 
    REAL                                       :: bm_supply_n            !! carbon growth sustainable by n_avail, considering costf 
    REAL, DIMENSION(npts,nvm)                  :: fcn_wood_act           !! Actual wood CN to leaf CN ratio; needed for tree PFTs which have a rigid wood CNP concentration

    ! Phosphorus cycle
    REAL, DIMENSION(npts,nvm)                  :: p_alloc_tot            !! phosphorus growth (g(P)/m2/dt)
    REAL , DIMENSION(npts,nvm)                 :: np_leaf                !! phosphorus concentration in leaves (gC/g(P))
    REAL , DIMENSION(npts,nvm)                 :: transloc_P             !! Transloc variables Phosphorus
    REAL                                       :: costf_P                !! phosphorus cost of a unit carbon growth given current C partitioning and phosphorus concentration
    REAL                                       :: deltanp,deltanpmax     !! (maximum) change in leaf phosphorus concentration
    REAL                                       :: p_avail                !! phosphorus available to growth (dummy)
    REAL                                       :: bm_supply_p            !! carbon growth sustainable by p_avail, considering costf_P
    REAL, DIMENSION(npts,nvm)                  :: fnp_wood_act           !! Actual wood NP to leaf NP ratio; needed for tree PFTs which have a rigid wood CNP concentration

    ! Carbon, Nitrogen & Phosphorus cycle
    REAL, DIMENSION(npts,nvm,nelements)        :: alloc_tot              !! carbon & nitrogen & phosphorus growth (g(P)/m2/dt)

    ! Functional allocation scheme
    ! water stress factor, based on water potential gradient in plant relative to max
    REAL, DIMENSION(npts,nvm)                  :: wstress_fac
    REAL, DIMENSION(npts,nvm)                  :: sstress_fac            !! soilstress_fac replaces to account also for nutrient stresss
    REAL, DIMENSION(npts,nvm,nparts)           :: tmp_resp_maint_part

    CHARACTER(LEN=2), DIMENSION(nelements)            :: element_str            !! string suffix indicating element 
    REAL                                       :: frac_growthresp_dyn    !! Fraction of gpp used for growth respiration (-)
                                                                                !! considers the special case at the leaf onset : frac_growthresp_dyn=0
                                                                                !! DSG: WRONG! it is fraction of allocated carbon not GPP
    REAL                                       :: lai_temp,veget_temp                                              
    REAL,PARAMETER                             :: lab_fac_fixed=0.95     !! the labile fraction of the labile pool is kept constant on 95%

!gmjc
    ! weighted leaf age (leaf_age * fraction of the age class)
    REAL, DIMENSION(npts,nvm)                      :: leaf_age_w
    REAL, DIMENSION(npts,nvm)                      :: sla_age2
    REAL, DIMENSION(npts,nvm)                      :: sla_age3
    REAL, DIMENSION(npts,nvm)                      :: sla_age4
    ! SLA max and SLA min affected by N fertilization
    REAL, DIMENSION(npts,nvm)                      :: sla_max_Nfert
    REAL, DIMENSION(npts,nvm)                      :: sla_min_Nfert
    REAL :: wear2wsh 
    REAL :: lamrep
    REAL :: fdev_lam
    REAL :: fdev_ear
    REAL :: fdev_st
    REAL :: fdev_above
    REAL :: fdev_below
!end gmjc
    REAL, DIMENSION(nparts) :: npartsvar     !! Biomass components of the model tree  
    REAL, DIMENSION(ncirc) :: ncircvar      !! 
      
!_ ================================================================================================================================

    IF (printlev.GE.3) THEN
       WRITE(numout,*) 'Entering functional allocation growth' 
       WRITE(numout,*) 'biomass(test_grid,test_pft,:,icarbon)   :',biomass(test_grid,test_pft,:,icarbon)
       WRITE(numout,*) 'biomass(test_grid,test_pft,:,initrogen) :',biomass(test_grid,test_pft,:,initrogen)
       WRITE(numout,*) 'biomass(test_grid,test_pft,:,iphosphorus) :',biomass(test_grid,test_pft,:,iphosphorus)
       WRITE(numout,*) 'n_uptake_daily(test_grid,test_pft,iammonium)   : ',  n_uptake_daily(test_grid,test_pft,iammonium)
       WRITE(numout,*) 'n_uptake_daily(test_grid,test_pft,initrate)    : ',  n_uptake_daily(test_grid,test_pft,initrate)
       WRITE(numout,*) 'p_uptake_daily(test_grid,test_pft)             : ',  p_uptake_daily(test_grid,test_pft)
    ENDIF

!! 1. Initialize

    !! 1.1 First call only
    IF ( firstcall ) THEN
      
       firstcall = .FALSE.

    ENDIF ! first call

    !! 1.2 Initialize variables at every call
    bm_alloc(:,:,:,:)          = zero 
    bm_alloc_tot(:,:)          = zero 
    n_alloc_tot(:,:)           = zero 
    p_alloc_tot(:,:)           = zero 
    qm_height(:,:)             = zero
    delta_ba                   = zero
    lai_target(:,:)            = zero
    resp_maint(:,:)            = zero
    resp_growth(:,:)           = zero
    lstress_fac(:,:)           = zero
    sigma(:,:)                 = zero
    gammas(:,:)                = zero
    k_latosa(:,:)              = zero
    store_circ_class_ba(:,:,:) = zero
    store_delta_ba(:,:,:)      = zero
    ltor(:,:)                  = val_exp
    npp(:,:)                   = zero
    residual(:,:)              = zero
    lab_fac(:,1)               = zero
    tmp_resp_maint_part(:,:,:) = zero

    ! remember the size of labile and reserve carbon pool to diagnose the
    ! bm_allo to them
    bm_res_lab(:,:,1)       = biomass(:,:,icarbres,icarbon)
    bm_res_lab(:,:,2)       = biomass(:,:,ilabile ,icarbon)

    root_reserve(:)           = 1.0          !! fraction of max root biomass which are covered by carbon reserve

    !! 1.2.2 Initialize check for mass balance closure
    !  The mass balance is calculated at the end of this routine
    !  in section 8
    pool_start = zero
    DO ipar = 1,nparts
       DO iele = 1,nelements
          !  Initial biomass pool
          pool_start(:,:,iele) = pool_start(:,:,iele) + &
               (biomass(:,:,ipar,iele) * veget_max(:,:))
       ENDDO
    ENDDO

    circ_class_n(:,:,1)=ind(:,:)

    DO ipar= 1,nparts
       DO iele = 1,nelements
          WHERE (circ_class_n(:,:,1).GT.zero) 
             circ_class_biomass(:,:,1,ipar,iele) = biomass(:,:,ipar,iele)/circ_class_n(:,:,1)
          ELSEWHERE   
             circ_class_biomass(:,:,1,ipar,iele) = zero
          ENDWHERE
       ENDDO
    ENDDO

    !! 1.2.3 Calculate LAI threshold below which carbohydrate reserve is used.
    !  Lai_max is a PFT-dependent parameter specified in stomate_constants.f90
    ! +++CHECK+++
    ! Can we make this a function of Cs or rue_longterm? this double prescribed value does not make
    ! to much sense to me. It is not really dynamic.
    lai_happy(:) = lai_max(:) * lai_max_to_happy(:)
    ! +++++++++++
  
    ! some initisalisations
    ! nitrogen concentration in leaves as CN
    WHERE(biomass(:,:,ileaf,initrogen).GT.min_stomate) 
       cn_leaf(:,:)=biomass(:,:,ileaf,icarbon)/biomass(:,:,ileaf,initrogen)
    ELSEWHERE
       cn_leaf(:,:)=cn_leaf_avg_season(:,:)
    ENDWHERE
    

    ! phosphorus concentration in leaves as NP
    WHERE(biomass(:,:,ileaf,iphosphorus) .GT. min_stomate) 
       np_leaf(:,:)=biomass(:,:,ileaf,initrogen)/biomass(:,:,ileaf,iphosphorus)
    ELSEWHERE
       np_leaf(:,:)=np_leaf_avg_season(:,:)
    ENDWHERE

    ! for the situation of impose_cn and impose_cn_restart
    ! cn_leaf should be forced to cn_leaf_avg_season
    IF (impose_cn .AND. impose_cn_restart) THEN
      cn_leaf(:,:)=cn_leaf_avg_season(:,:)
    ENDIF
    ! for the situation of impose_np and impose_np_restart
    ! np_leaf should be forced to np_leaf_avg_season
    IF (impose_np .AND. impose_np_restart) THEN
      np_leaf(:,:)=np_leaf_avg_season(:,:)
    ENDIF

    ! Calculate the change in fcn_wood due to changes in leaf CN for woody PFTs
    DO j=2,nvm
      IF ( is_tree(j) ) THEN
          fcn_wood_act(:,j) = fcn_wood(j) * cn_leaf(:,j)/cn_leaf_min(j)
          fnp_wood_act(:,j) = fnp_wood(j) * np_leaf(:,j)/np_leaf_min(j)
      ELSE
          ! for grass woody CNP varies with leaf
          fcn_wood_act(:,j) = fcn_wood(j) * un
          fnp_wood_act(:,j) = fnp_wood(j) * un
      ENDIF
    ENDDO

    !! Add mineral nutrients taken up from the soil to the plant labile nutrient pools
    biomass(:,:,ilabile,initrogen) = biomass(:,:,ilabile,initrogen) & 
         + n_uptake_daily(:,:,iammonium)+n_uptake_daily(:,:,initrate) &
! JCADD 20170427 BNF for grassland
         + BNF_clover_daily(:,:)
! END JCADD

    biomass(:,:,ilabile,iphosphorus) = biomass(:,:,ilabile,iphosphorus) & 
         + p_uptake_daily(:,:)


    !DSG debug: ===================================================
    IF (dsg_debug) THEN
       WRITE(6,*) 'location in code: '
       WRITE(6,*) 'stomate_growth_fun: after nutrient uptake'
       IF (printlev.GE.3) THEN
          WRITE(numout,*) 'location in code: '
          WRITE(numout,*) 'stomate_growth_fun: after nutrient uptake'
          WRITE(numout,*) 'biomass           (test_grid,test_pft ,ilabile,initrogen)',       biomass(test_grid,test_pft  ,ilabile,initrogen)
          WRITE(numout,*) 'biomass           (test_grid,test_pft ,ilabile,iphosphorus)',     biomass(test_grid,test_pft  ,ilabile,iphosphorus)
       ENDIF
       CALL check_mass(npts,biomass(:,:,:,:),'stomate_growth_fun: after nutrient uptake')
    ENDIF


 !! 2. Use carbohydrate reserve to support growth

    ! Save old leaf mass, biomass got last updated in stomate_phenology.f90
!JC separate tissue age
!    lm_old(:,:) = biomass(:,:,ileaf,icarbon)
    lm_old(:,:,:) = biomass(:,:,:,icarbon)
!End JC separate tissue age

    ! lai for bare soil is by definition zero
    lai(:,ibare_sechiba) = zero

    DO j = 2, nvm ! Loop over # PFTs 

       !! 2.1 Calculate demand for carbohydrate reserve to support leaf and root growth.
       !  Maximum time (days) since start of the growing season during which carbohydrate 
       !  may be used
       IF ( is_tree(j) ) THEN

          reserve_time = reserve_time_tree    

       ELSE

          reserve_time = reserve_time_grass

       ENDIF

       ! Growth is only supported by the use of carbohydrate reserves if the following 
       ! conditions are  statisfied:\n
       ! - PFT is not senescent;\n
       ! - LAI must be low (i.e. below ::lai_happy) and\n
       ! - Day of year of the simulation is in the beginning of the growing season.
   
       DO ipts = 1,npts

          ! Calculate lai
          lai(ipts,j) = biomass_to_lai(biomass(ipts,j,ileaf,icarbon),j)

          ! We might need the c0_alloc factor, so let's calculate it.
          c0_alloc(ipts,j) = calculate_c0_alloc(ipts, j, tau_root(j), &
               tau_sap(j))

       ENDDO

!JC MOD021 use reserve for regrowth after harvest
! in case some times the grasses just die after harvest
       WHERE (( ( biomass(:,j,ileaf,icarbon) .GT. min_stomate ) .AND. & 
            ( .NOT. senescence(:,j) ) .AND. &
            ( lai(:,j) .LT. lai_happy(j) ) .AND. &
            ( when_growthinit(:,j) .LT. reserve_time ) ) .OR. &
            ( when_growthinit_cut(:,j) .LT. 2.0))
          
          ! Tell the labile and resource pool to use its reserve
          use_reserve(:,j) = 1.0

       ENDWHERE
 
    ENDDO ! loop over # PFTs
 
 !! 3. Initialize allocation

    DO j = 2, nvm ! Loop over # PFTs 
 
       !!  3.1 Calculate scaling factors, temperature sensitivity, target lai to decide on
       !!  reserve use, labile fraction, labile biomass and total allocatable biomass 
       
       ! Convert long term temperature from K to C
       tl(:) = t2m(:) - ZeroCelsius

       DO ipts = 1, npts

          IF (veget_max(ipts,j) .LE. min_stomate) THEN

              ! this vegetation type is not present, so no reason to do the 
              ! calculation. CYCLE will take us out of the innermost DO loop
               CYCLE

          ENDIF

          !! 3.1 Water stress
          !  The waterstress factor varies between 0.1 and 1 and is calculated
          !  from ::moiavail_growingseason. The latter is only used in the allometric 
          !  allocation and its time integral is determined by tau_sap for trees
          !  (see constantes_mtc.f90 for tau_sap and see pft_constantes.f90 for 
          !  the definition of tau_hum_growingseason). The time integral for 
          !  grasses and crops is a prescribed constant (see constantes.f90). For
          !  trees KF (and indirecrtly LF) and for grasses LF are multiplied
          !  by wstress. Because the calculated values are too low for its purpose
          !  Sonke Zhaele multiply it by two in the N-branch (see stomate_season.f90). 
          !  This approach maintains the physiological basis of KF while combining it 
          !  with a simple multiplicative factor for water stress. Clearly after 
          !  multiplication with 2, wstress is closer to 1 and will thus result in a 
          !  KF values closer to the physiologically expected KF. We did not see the 
          !  need to multiply by 2 because the way we now calculate ::moiavail_growingseason
          !  is less volatile than before. Before it ranged between 0 and 1, now the 
          !  range is more like 0.4 to 0.9.

!DSGlie_detector 
!         that is not true: moiavail_growingseason is still multiplies by a
!         factor of two in stomate_season. So if you want larger effect of water
!         stress on ltor, go to stomate_season and remove the two
!DSGlie_detector 


          ! Veget is now calculated from Pgap to be fully consistent within the model. Hence
          ! dividing by veget_max gives a value between 0 and 1 that denotes the amount of 
          ! light reaching the forest floor.
           IF (veget_max(ipts,j) .GT. min_stomate) THEN
             lstress_fac(ipts,j) = un - (veget(ipts,j) / veget_max(ipts,j))
          ELSE
             lstress_fac(ipts,j) = zero
          ENDIF

          !! 3.2 Initialize scaling factors
          ! Stand level scaling factors
          LF(ipts,j) = 1._r_std

          ! Tree level scaling factors
          ltor(ipts,j) = 1._r_std
          circ_class_height_eff(:) = 1._r_std


          !! 3.3 Calculate structural characteristics

          ! Target lai is calculated at the stand level for the tree height of a 
          ! virtual tree with the mean basal area or the so called quadratic mean diameter
          npartsvar = circ_class_biomass(ipts,j,1,:,icarbon)
          ncircvar = circ_class_n(ipts,j,:)

          qm_dia(ipts,j) = wood_to_qmdia(npartsvar , ncircvar, j)
          qm_height(ipts,j) = wood_to_qmheight(npartsvar, ncircvar, j)

          IF ( is_tree(j) ) THEN
            wstress_fac(ipts,j)=MAX(moiavail_growingseason(ipts,j),0.1)
            sstress_fac(ipts,j)=MAX(MIN(moiavail_growingseason(ipts,j),  &
                                            pstress_season(ipts,j),         &
                                            nstress_season(ipts,j) )        &
                                       ,0.1)
          ELSE
            ! JC MOD044 grasses should be more tolerant for water stress
            ! Thus a threshold is introduced to avoid too high water stress
            wstress_fac(ipts,j)=MAX(moiavail_growingseason(ipts,j),wstress_lowgrass)
            sstress_fac(ipts,j)=MAX(MIN(wstress_fac(ipts,j),  &
                                            pstress_season(ipts,j),         &
                                            nstress_season(ipts,j) )        &
                                       ,sstress_lowgrass)
          ENDIF

          IF (dsg_debug) THEN
                      IF (j.EQ.test_pft .AND. ipts==test_grid) THEN
                        WRITE(numout,*) 'wstress_fac(ipts,j)',wstress_fac(ipts,j)
                        WRITE(numout,*) 'sstress_fac(ipts,j)',sstress_fac(ipts,j)
                      ENDIF 
          ENDIF

          !! 3.4 Calculate allocation factors for trees and grasses 
          IF ( SUM(biomass(ipts,j,:,icarbon)) .GT. min_stomate ) THEN

             ! Trees
             IF (is_tree(j)) THEN

                ! Note that KF may already be calculated in stomate_prescribe.f90 (if called)
                ! it is recalculated because the biomass pools for grasses and crops
                ! may have been changed in stomate_phenology.f90. Trees were added to this
                ! calculation just to be consistent.
                
                ! To be fully consistent with the hydraulic limitations and pipe theory,
                ! k_latosa_zero should be calculated from equation (18) in Magnani et al.
                ! To do so, total hydraulic resistance and tree height need to known. This
                ! poses a problem as the resistance depends on the leaf area and the leaf 
                ! area on the resistance. There is no independent equation and equations 12
                ! and 18 depend on each other and substitution would be circular. Hence 
                ! prescribed k_latosa_zero values were obtained from observational records 
                ! and are given in mtc_parameters.f90.

                ! The relationship between height and k_latosa as reported in McDowell 
                ! et al 2002 and Novick et al 2009 is implemented to adjust k_latosa for 
                ! the height of the stand.  The slope of the relationship is calculated in 
                ! stomate_data.f90 This did NOT result in a realistic model behavior.
                !!$ k_latosa(ipts,j) = wstress_fac(ipts,j) * & 
                !!$     (k_latosa_max(j) - latosa_height(j) * qm_height(ipts,j))

                ! Alternatively, k_latosa is also reported to be a function of diameter 
                ! (i.e. stand thinning, Simonin et al 2006, Tree Physiology, 26:493-503).
                ! Here the relationship with thinning was interpreted as a realtionship with
                ! light stress.

                k_latosa(ipts,j) = (k_latosa_adapt(ipts,j) + (lstress_fac(ipts,j) * &
                                   (k_latosa_max(j)-k_latosa_min(j))))*wstress_fac(ipts,j)
                ! +++++++++++
               
                ! Also k_latosa has been reported to be a function of CO2 concentration 
                ! (Atwell et al. 2003, Tree Physiology, 23:13-21 and Pakati et al. 2000, 
                ! Global Change Biology, 6:889-897). This effect is not accounted for in 
                ! the current code

                ! Scaling factor to convert sapwood mass into leaf mass (KF)
                ! derived from
                ! LA_ind = k1 * SA_ind, k1=latosa (pipe-model)
                ! <=> Cl * vm/ind * sla = k1 * Cs * vm/ind / wooddens / tree_ff / height_new
                ! <=> Cl = Cs * k1 / wooddens / tree_ff/ height_new /sla
                ! <=> Cl = Cs * KF / height_new, where KF = k1 / (wooddens * sla * tree_ff)
                ! (1) Cl = Cs * KF / height_new
                KF_old = KF(ipts,j) 
                KF(ipts,j) = k_latosa(ipts,j) / (sla_calc(ipts,j) * pipe_density(j) * tree_ff(j))
                
                ! KF of the previous time step was stored in ::KF_old to check its absolute 
                ! change. If this absolute change is too big the whole allocation will crash
                ! because it will calculate negative increments which are compensated by 
                ! positive increments that exceed the available carbon for allocation. This 
                ! would suggest that for examples the plant destructs leaves and uses the 
                ! available carbon to produce more roots. This is would repesent an unwanted 
                ! outcome. Large changes from time step to another makes its difficult for 
                ! the scheme to ever reach allometric balance. This balance is needed for the
                ! allocation scheme to allow 'ordinary allocation', which in turn is needed
                ! to make use of the allocation rule of Dhote and Deleuze. It needs to be 
                ! avoided that the code spends too much time in phenological growth and the 
                ! if-then statements that help to restore allometric balance. For this reason
                ! the absolute change in KF from one time step to another are truncated.
                IF (KF_old - KF(ipts,j) .GT. max_delta_KF ) THEN
        
                   IF(ld_warn)THEN
                      WRITE(numout,*) 'WARNING 2: KF was truncated'
                      WRITE(numout,*) 'WARNING 2: PFT, ipts: ',j,ipts
                      WRITE(numout,'(A,3F20.10)') 'WARNING 2: KF_old, KF(ipts,j), max_delta_KF: ',&
                           KF_old, KF(ipts,j), max_delta_KF
                   ENDIF

                   ! Add maximum absolute change
                   KF(ipts,j) = KF_old - max_delta_KF
                   
                   IF(ld_warn)THEN
                      WRITE(numout,'(A,3F20.10)') 'WARNING 2: Reset, KF_old, KF(ipts,j): ',&
                           KF_old, KF(ipts,j)
                   ENDIF
                ELSEIF (KF_old - KF(ipts,j) .LT. -max_delta_KF) THEN
                   
                   IF(ld_warn)THEN
                      WRITE(numout,*) 'WARNING 3: KF was truncated'
                      WRITE(numout,*) 'WARNING 3: PFT, ipts: ',j,ipts
                      WRITE(numout,'(A,3F20.10)') 'WARNING 3: KF_old, KF(ipts,j), max_delta_KF: ',&
                           KF_old, KF(ipts,j), -max_delta_KF
                   ENDIF

                   ! Remove maximum absolute change
                   KF(ipts,j) = KF_old + max_delta_KF
                   
                   IF(ld_warn)THEN
                      WRITE(numout,'(A,3F20.10)') 'WARNING 3: Reset, KF_old, KF(ipts,j): ',&
                        KF_old, KF(ipts,j)
                   ENDIF
                ELSE
                   ! The change in KF is acceptable no action required
                ENDIF
    
                ! Scaling factor to convert sapwood mass into root mass  (LF) 
                ! derived from 
                ! Cs = c0 * height * Cr (Magnani 2000)
                ! Cr = Cs / c0 / height_new
                ! scaling parameter between leaf and root mass, derived from
                ! Cr = Cs / c0 / height_new
                ! let Cs = Cl / KF * height_new
                ! <=> Cr = ( Cl * height_new / KF ) / ( c0 * height_new )
                ! <=> Cl = Cr * KF * c0
                ! <=> Cl = Cr * LF, where LF = KF * c0
                ! (2) Cl = Cr * LF
                ! +++CHECK+++
                ! How do we want to account for waterstress? wstress is accounted for in c0_alloc
                ! DSGlie_detector: is that true? I dont think so
                LF(ipts,j) = c0_alloc(ipts,j) * KF(ipts,j) *sstress_fac(ipts,j)  !&


                ! Calculate non-nitrogen stressed leaf to root ratio to calculate the
                ! allocation to the reserves. Should be multiplied by a nitrogen stress
                ! have a look in OCN. This code should be considered as a placeholder
                
                ! DSG: it MUST be multiplied with nstress; if not, then plants under N
                ! stress increase their CN ratio and can thereby store more
                ! carbon in biomass/higher LAI -> achieve higher growth rates than under
                ! optimal N conditions; which is nonesense
                ltor(ipts,j) = c0_alloc(ipts,j) * KF(ipts,j)

                !---TEMP---
                IF (j.EQ.test_pft .AND. ld_alloc .AND. ipts == test_grid) THEN
                   WRITE(numout,*) 'KF, LF, ', KF(ipts,j), LF(ipts,j)
                   WRITE(numout,*) 'c0_alloc, ', c0_alloc(ipts,j)
                   WRITE(numout,*) 'tau_root, tau_sap, ', tau_root(j), tau_sap(j)
                   WRITE(numout,*) 'k_root, k_sap, ', k_root(j), k_sap(j)
                   WRITE(numout,*) 'ltor, ', ltor(ipts,j)
                   
                ENDIF
                !----------

!JC GRASS KF LF ltor start
             ! Grasses and crops
             ELSEIF (.NOT. is_tree(j)) THEN
                
                !+++CHECK+++
                ! Similar to ::k_latosa for trees we defined it for grasses. Note that for trees
                ! the definition is supported by some observations. For grasses we didn't look very
                ! hard to check the literature. Someone interested in grasses should invest some
                ! time in this issue and replace this code by a parameter that can be derived
                ! from observations. In the end it was decided to use the same variable name for 
                ! grasses, crops and trees as that allowed us to optimize this parameter.
! JC MOD042 for grasses, there will be no water stress to k_latosa, grasses do
! not actually need sapwood for water transportation
!                k_latosa(ipts,j) = (k_latosa_adapt(ipts,j) + (lstress_fac(ipts,j) * &
!                    (k_latosa_max(j)-k_latosa_min(j))))*wstress_fac(ipts,j)
                k_latosa(ipts,j) = (k_latosa_adapt(ipts,j) + (lstress_fac(ipts,j) * &
                    (k_latosa_max(j)-k_latosa_min(j))))

                IF (GRM_devstage(ipts,j) .LE. 1.0) THEN
                  k_latosa(ipts,j) = (k_latosa_adapt(ipts,j) + &
                      (k_latosa_max(j)-k_latosa_min(j)))!*wstress_fac(ipts,j)
                ELSE IF (GRM_devstage(ipts,j) .LT. 2.0) THEN
                  k_latosa(ipts,j) = (k_latosa_adapt(ipts,j) + &
                      (2.0-GRM_devstage(ipts,j)) * &
                      (k_latosa_max(j)-k_latosa_min(j)))!*wstress_fac(ipts,j)
                ELSE
                  k_latosa(ipts,j) = k_latosa_adapt(ipts,j)!*wstress_fac(ipts,j)
                ENDIF
                ! The mass of the structural carbon relates to the mass of the leaves through
                ! a prescribed parameter ::k_latosa
                ! waterstress on KF is accounted for in k_latosa
                   KF(ipts,j) = k_latosa(ipts,j)
                !++++++++++

!JC MOD039 No stem before GRM_devstage < 1.0
                IF (.NOT. natural(j)) THEN
                  IF (GRM_devstage(ipts,j) .LE. 1.0) THEN
                    k_latosa(ipts,j) = 100.
                    KF(ipts,j) = k_latosa(ipts,j)
                  ENDIF
                ENDIF               
                ! Stressed root allocation, wstress is accounted for in c0_alloc  
                ! DSGlie_detector: is that true? Waterstress is in KF
!DSGsstress     The leaf to root ratio (LF) should be lower when nutrient stressed
!JC MOD028 since KF decrease from 2 to 0.5, the LF follows KF is not realistic,
! this will cause too much root, and prevent regrowth.
! so I fixed the LF to a more stable rate, which is more close to PaSim alloc
! the leaf root ratio is generally more stable accross the year.
! BUT not sure if need double count wstress_fac
!                LF(ipts,j) = c0_alloc(ipts,j) * KF(ipts,j)*sstress_fac(ipts,j)
                LF(ipts,j) = c0_alloc(ipts,j) * k_latosa_max(j) * sstress_fac(ipts,j)

!DSGwtf?! : why do you need LF and ltor? and why do you do the exactly same calculation again
!again?
                ! Calculate nutrient stressed leaf to root ratio to calculate the allocation to the reserves           
                ltor(ipts,j) = LF(ipts,j)
!DSGwtf?!
                !---TEMP---
                IF (j.EQ.test_pft .AND. ld_alloc .AND. ipts == test_grid) THEN
                   WRITE(numout,*) 'KF, LF, ', KF(ipts,j), LF(ipts,j)
                   WRITE(numout,*) 'k_latosa_adapt',k_latosa_adapt(ipts,j)
                   WRITE(numout,*) 'c0_alloc, ', c0_alloc(ipts,j)
                   WRITE(numout,*) 'lstress_fac, ',lstress_fac(ipts,j)
                   WRITE(numout,*) 'sstress_fac, ',sstress_fac(ipts,j)
                   WRITE(numout,*) 'tau_root, tau_sap, ', tau_root(j), tau_sap(j)
                   WRITE(numout,*) 'k_root, k_sap, ', k_root(j), k_sap(j)
                   WRITE(numout,*) 'ltor, ', ltor(ipts,j)
                   
                ENDIF
                !----------
                
             ENDIF

          ENDIF

          !+++CHECK+++
          !! 3.5 Calculate optimal LAI
          !  The calculation of the optimal LAI was copied and adjusted from O-CN. In O-CN it
          !  was also used in the allocation but that seems to be inconsistent with the allometric
          !  rules that are implemented. Say that the actual LAI is below the optimal LAI then
          !  the O-CN approach will keep pumping carbon to grow the optimal LAI. If we would apply
          !  the same method it means that during this phase the rule of Deleuze and Dhote would
          !  not be used. For that reason we dropped the use of LAI_optimal and replaced it by
          !  an allometric-based Cl_target value. Initially, lai_target was still calculated as 
          !  described below and used in the calculation of the reserves.
          !  Further testing showed that for some parameter sets lai_target was over 8 whereas the
          !  realized lai was close to 4. This leaves us with a frustrated plant that will invest a
          !  lot in its reserves but can never use them because it is constrained by the allometric
          !  rules. To grow an LAI of 8 it would need to have a crazy sapwoodmass.
          !  At a more fundamental level it is clear why the plant's LAI should not exceed lai_target
          !  because then it costs more to produce and maintain the leaf than that the new leaf can
          !  produce but there is no reason why the plant should try to reach lai_target. For these
          !  reasons it was decided to abandon this approach to lai_target and simply replace 
          !  lai_target by Cl_target * sla

          !! 3.5.1 Scaling factor
          !  Scaling factor to convert variables to the individual plant
          !  Different approach between the DGVM and statitic approach
          IF (ok_dgvm) THEN

             ! The DGVM does currently NOT work with the new allocation, consider this as
             ! placeholder. The original code had two different transformations to 
             ! calculate the scalars. Both could be used but the units will differ.
             ! When fixing the DGVM check which quantities need to be multiplied by scal 
             ! scal = ind(ipts,j) * cn_ind(ipts,j) / veget_max(ipts,j)
             scal(ipts,j) = veget_max(ipts,j) / ind(ipts,j)  
          ELSE

             ! circ_class_biomass contain the data at the tree level
             ! no conversion required
             scal(ipts,j) = 1.

          ENDIF

          !! 3.5.2 Calculate lai_target based on the allometric rules
          IF ( SUM(biomass(ipts,j,:,icarbon)) .GT. min_stomate ) THEN

          IF ( is_tree(j)) THEN

             ! Basal area at the tree level (m2 tree-1)
             circ_class_ba_eff(:) = wood_to_ba_eff(circ_class_biomass(ipts,j,1,:,icarbon),j)

             ! Current biomass pools per tree (gC tree^-1) 
             ! We will have different trees so this has to be calculated from the diameter relationships            
             Cs(:) = ( circ_class_biomass(ipts,j,:,isapabove,icarbon) + &
                  circ_class_biomass(ipts,j,:,isapbelow,icarbon) ) * scal(ipts,j)
             Cr(:) = circ_class_biomass(ipts,j,:,iroot,icarbon) * scal(ipts,j)
             Cl(:) = circ_class_biomass(ipts,j,:,ileaf,icarbon) * scal(ipts,j)
             Ch(:) = ( circ_class_biomass(ipts,j,:,iheartabove,icarbon) + &
                  circ_class_biomass(ipts,j,:,iheartbelow,icarbon) ) * scal(ipts,j)

             DO l = 1,ncirc 

                !  Calculate tree height
                circ_class_height_eff(l) = pipe_tune2(j)*(4/pi*circ_class_ba_eff(l))**(pipe_tune3(j)/2)

                !  Do the biomass pools respect the pipe model?
                !  Do the current leaf, sapwood and root components respect the allometric 
                !  constraints? Due to plant phenology it is possible that we have too much 
                !  sapwood compared to the leaf and root mass (i.e. in early spring). 
                !  Calculate the optimal root and leaf mass, given the current wood mass 
                !  by using the basic allometric relationships. Calculate the optimal sapwood
                !  mass as a function of the current leaf and root mass.
                Cl_target(l) = MAX( KF(ipts,j) * Cs(l) / circ_class_height_eff(l), &
                     Cr(l) * LF(ipts,j) , Cl(l))
                Cs_target(l) = MAX( Cl(l) / KF(ipts,j) * circ_class_height_eff(l), &
                     Cr(l) * LF(ipts,j) / KF(ipts,j) * circ_class_height_eff(l) , Cs(l))

                ! Check dimensions of the trees
                ! If Cs = Cs_target then ba and height are correct, else calculate the correct dimensions
                IF ( Cs_target(l) - Cs(l) .GT. min_stomate ) THEN

                   ! If Cs = Cs_target then dia and height are correct. However, if Cl = Cl_target
                   ! or Cr = Cr_target then dia and height need to be re-estimated. Cs_target should
                   ! satify the relationship Cl/Cs = KF/height where height is a function of Cs_target
                   ! <=> (KF*Cs_target)/(pipe_tune2*(Cs_target+Ch)/pi
                   ! /4)**(pipe_tune3/(2+pipe_tune3)) = Cl_target
                   ! Search Cs needed to sustain the max of Cl or Cr.
                   !  Search max of Cl and Cr first 
                   Cl_target(l) = MAX(Cl(l), Cr(l)*LF(ipts,j))

                   !---TEMP---
                   IF (j.EQ.test_pft .AND. ld_alloc .AND. ipts==test_grid) THEN
                      WRITE(numout,*) 'circ_class_height_eff, ', circ_class_height_eff(l)
                      WRITE(numout,*) 'KF, LF, ', KF(ipts,j), LF(ipts,j)
                      WRITE(numout,*) 'Cl_target,Cl, ', Cl_target(l), Cl(l)
                      WRITE(numout,*) 'Cs, ', Cs(l)
                      WRITE(numout,*) 'Cr, ', Cr(l)
                      WRITE(numout,*) 'Ch, ', Ch(l)
                      WRITE(numout,*) 'Cs_target, ', Cs_target(l)
                      WRITE(numout,*) 'PRELIM. Cl_target, ', Cl_target(l)
                      WRITE(numout,*) 'PRELIM. Cr_target, ', Cr_target(l)       
                   ENDIF

                   !----------
                   !DSG: dont call the random number generator:!
!                   Cs_target(l) =  newX(KF(ipts,j), Ch(l),&
!                        & pipe_tune2(j), pipe_tune3(j), Cl_target(l),&
!                        & tree_ff(j)*pipe_density(j)*pi/4*pipe_tune2(j), Cs(l),&
!                        & 2*Cs(l), 2, j, ipts)

                   ! Recalculate height and ba from the correct
                   !  Cs_target
                   circ_class_height_eff(l) = Cs_target(l)*KF(ipts,j)&
                        &/Cl_target(l)
                   circ_class_ba_eff(l) = pi/4*(circ_class_height_eff(l)&
                        &/pipe_tune2(j))**(2/pipe_tune3(j))
                   Cl_target(l) = KF(ipts,j) * Cs_target(l) /&
                        & circ_class_height_eff(l)
                   Cr_target(l) = Cl_target(l) / LF(ipts,j)

                   !---TEMP---
                   IF (j.EQ.test_pft .AND. ld_alloc .AND. ipts==test_grid) THEN
                      WRITE(numout,*) 'New Cl_target, ', Cl_target(l)
                      WRITE(numout,*) 'New Cr_target, ', Cr_target(l)       
                   ENDIF
                   !---------- 

                ENDIF

             ENDDO

!!! ========================================= !!!             
!DSG: lai_target can get negative or super large e19.
             ! Calculate lai_target
             lai_target(ipts,j) = SUM(Cl_target(:)*circ_class_n(ipts,j,:)) * sla_calc(ipts,j)
!!! ========================================= !!!             
!!! ========================================= !!!             

             !---TEMP---
             IF (j.EQ.test_pft .AND. ld_alloc .AND. ipts==test_grid) THEN
                WRITE(numout,*) 'lai_target, ',lai_target(ipts,j)
                WRITE(numout,*) 'circ_class_n',circ_class_n(ipts,j,:)
                WRITE(numout,*) 'Cl_target',Cl_target(:)
                WRITE(numout,*) 'KF, LF, ', KF(ipts,j), LF(ipts,j)
             ENDIF
             !----------

          ! Grasses and croplands
          ELSEIF ( .NOT. is_tree(j)) THEN
             
             ! Current biomass pools per grass/crop (gC ind^-1)
             ! Cs has too many dimensions for grass/crops. To have a consistent notation the same variables
             ! are used as for trees but the dimension of Cs, Cl and Cr i.e. ::ncirc should be ignored            
             Cs(1) = biomass(ipts,j,isapabove,icarbon) * scal(ipts,j)
             Cr(1) = biomass(ipts,j,iroot,icarbon) * scal(ipts,j)
             Cl(1) = biomass(ipts,j,ileaf,icarbon) * scal(ipts,j)
             Ch(1) = zero
   
             ! Do the biomass pools respect the pipe model?
             ! Do the current leaf, sapwood and root components respect the allometric 
             ! constraints? Calculate the optimal root and leaf mass, given the current wood mass 
             ! by using the basic allometric relationships. Calculate the optimal sapwood
             ! mass as a function of the current leaf and root mass.
             Cl_target(1) = MAX( Cs(1) * KF(ipts,j) , Cr(1) * LF(ipts,j), Cl(1) )
             Cs_target(1) = MAX( Cl_target(1) / KF(ipts,j), Cr(1) * LF(ipts,j) / KF(ipts,j), Cs(1) ) 
             Cr_target(1) = MAX( Cl_target(1) / LF(ipts,j), Cs_target(1) * KF(ipts,j) / LF(ipts,j), Cr(1) )

             ! Calculate lai_target
! JC MOD sla_calc
             lai_target(ipts,j) = Cl_target(1) * circ_class_n(ipts,j,1) * sla_calc(ipts,j)

             !---TEMP---
             IF (j.EQ.test_pft .AND. ld_alloc .AND. ipts==test_grid) THEN
                WRITE(numout,*) 'lai_target, ',lai_target(ipts,j)
                WRITE(numout,*) 'circ_class_n',circ_class_n(ipts,j,1)
                WRITE(numout,*) 'Cl_target',Cl_target(1)
                WRITE(numout,*) 'KF, LF, ', KF(ipts,j), LF(ipts,j)
             ENDIF
             !----------

          ENDIF

          ENDIF ! biomass > min_stomate
          

          !! 3.6 Calculate mean leaf age
          leaf_meanage = zero
          DO m = 1,nleafages
          
             leaf_meanage = leaf_meanage + leaf_age(ipts,j,ileaf,m) * leaf_frac(ipts,j,ileaf,m)
          
          ENDDO
             
          !! 3.7 Calculate labile fraction 
          !  Use constant labile fraction (to initiate on-set of leaves in spring)
          !  The dynamic lab_fac in ORCHIDEE-CAN decouples NPP from GPP to a
          !  degree which is not realistic anymore.
            
          ! DSG: lab_fac of (1,.95,) don't work with the allometric calculation as
          ! they are not precise enough or there are bugs. which lead to
          ! negative biomass due to overspending.
          lab_fac= lab_fac_fixed
        
          !! 3.8 Calculate allocatable carbon
          !  Total allocatable biomass during this time step determined from GPP.
          !  GPP was calculated as CO2 assimilation in enerbil.f90
          !  Under some exceptional conditions :gpp could be negative when
          !  the dark respiration exceeds the photosynthesis. When this happens
          !  the dark respiration is paid for by the labile and carbres pools

          deficit = zero

          IF ( (biomass(ipts,j,ilabile,icarbon) + gpp_daily(ipts,j) * dt) .LT. zero ) THEN

             deficit = (biomass(ipts,j,ilabile,icarbon) + gpp_daily(ipts,j) * dt)

           ! The deficit is less than the carbon reserve
             IF (-deficit .LE. biomass(ipts,j,icarbres,icarbon)) THEN

              ! Pay the deficit from the reserve pool
                biomass(ipts,j,icarbres,icarbon) = &
                     biomass(ipts,j,icarbres,icarbon) + deficit
                biomass(ipts,j,ilabile,icarbon) = &
                     biomass(ipts,j,ilabile,icarbon) - deficit

             ELSE

                ! Not enough carbon to pay the deficit, the individual 
                ! is going to die at the end of this day
                biomass(ipts,j,ilabile,icarbon) = &
                     biomass(ipts,j,ilabile,icarbon) + biomass(ipts,j,icarbres,icarbon) 
                biomass(ipts,j,icarbres,icarbon) = zero

                ! Truncate the dark respiration to the available carbon.  Now we
                ! should use up all the reserves.  If the plant has no leaves, it
                ! will die quickly after this.
                gpp_daily(ipts,j) = - biomass(ipts,j,ilabile,icarbon)/dt 

             ENDIF

          ENDIF
       
          ! Labile carbon pool after possible correction for dark respiration
          biomass(ipts,j,ilabile,icarbon) = biomass(ipts,j,ilabile,icarbon) + &
               gpp_daily(ipts,j) * dt

          IF (printlev.GE.3) THEN
             IF(ipts == test_grid .AND. j == test_pft)THEN
                WRITE(numout,*) 'Adding gpp to labile pool'
                WRITE(numout,*) 'deficit(already taken)',deficit !DSG: here the deficit is  already paied
                WRITE(numout,*) 'biomass(test_grid,test_pft,:,icarbon)',biomass(ipts,j,:,icarbon)
                WRITE(numout,*) 'gpp_daily(test_grid,test_pft)',gpp_daily(ipts,j)
                WRITE(numout,*) 'SUM(biomass(test_grid,test_pft,:,icarbon))',SUM(biomass(ipts,j,:,icarbon))
             ENDIF
          ENDIF

          !! 3.9 Calculate activity of labile carbon pool  
          !  Similar realtionship as that used for the temperature response of 
          !  maintenance respiration.  The parameters in the equation were calibrated
          !  to give a fraction of 0.1 of GPP at reference temperature tl (i.e. 10°C)
!DSGliedetector: you mean gtemp=1 for tl=10C, so the fraction will be 1*lab_fac
          !  Note that the temperature response has a lower slope than for respiration
          !  to avoid too large turnover rates at high temperature.
          IF (tl(ipts) .GT. -2.) THEN

             gtemp = EXP(308.56/4.*(1.0/56.02-1.0/(tl(ipts)+46.02)))

          ELSE 

             gtemp = zero            

          ENDIF
           
          ! ==========================================
          !DSGtriggerXYR938
          !  If there is a plant, and we are either at the very start or in the growing season
          !  not during senescences, calculate labile pool use for growth
          IF (ind(ipts,j) .GT. min_stomate .AND. & ! if there is a plant ...
               .NOT.senescence(ipts,j) &           ! ... which is not senescence ...
               .AND. ( biomass(ipts,j,ileaf,icarbon) .GT. min_stomate .OR. & !... and it has leaves or ...
               use_reserve(ipts,j) .GT. min_stomate ) &                      ! ... is allowed to grow 
               ) THEN  

             IF ((use_reserve(ipts,j) .GT. min_stomate) .AND. &
                 (biomass(ipts,j,icarbres,icarbon) .GT. min_stomate)) THEN
             !DSGtriggerXYR938
             ! ==========================================

                biomass(ipts,j,ilabile,icarbon) = biomass(ipts,j,ilabile,icarbon) + 0.10 * &
                                                   biomass(ipts,j,icarbres,icarbon)

                biomass(ipts,j,icarbres,icarbon) = biomass(ipts,j,icarbres,icarbon) * 0.90
             ENDIF

                ! The labile pool is filled. Re-calculate turnover of the labile pool. 
                ! Only if the labile pool is very small turnover will exceed 0.75 
                ! and the pool will thus be almost entirely emptied
             IF (biomass(ipts,j,ilabile,icarbon) .GT. min_stomate) THEN
                 gtemp = MAX(MIN( gtemp * lab_fac(ipts,j), un ), zero)
             ELSE
                 gtemp = zero
             ENDIF


          ENDIF

          IF (printlev.GE.3) THEN
             IF(ipts == test_grid .AND. j == test_pft)THEN
                WRITE(numout,*) 'After moving carbon between carbres and labile pool'
                WRITE(numout,*) 'SUM(biomass(test_grid,test_pft,:,icarbon))',SUM(biomass(ipts,j,:,icarbon))
                WRITE(numout,*) 'biomass(test_grid,test_pft,:,icarbon)',biomass(ipts,j,:,icarbon)
             ENDIF
          ENDIF


          !! 3.10 Calculate allocatable part of the labile pool
          !  If there is a plant, and we are 
          !  either at the very start or in the growing season not during 
          !  senescences and after senescence only if use_reserve > 0., calculate
          !  labile pool use for growth.
          IF (ind(ipts,j) .GT. min_stomate .AND. &
               .NOT.senescence(ipts,j)&
               .AND.( biomass(ipts,j,ileaf,icarbon) .GT. min_stomate .OR. &
               use_reserve(ipts,j) .GT. min_stomate ) &
             ) THEN

             ! Use carbon from the labile pool to allocate. The allometric (or 
             ! functional) allocation scheme transfers gpp to the labile pool 
             ! (see above) and then uses the labile pool (gpp + labile(t-1)) to sustain 
             ! growth. The fraction of the labile pool that can be used is a 
             ! function of the temperature and phenology. bm_alloc_tot in 
             ! gC m-2 dt-1
             bm_alloc_tot(ipts,j) = gtemp * biomass(ipts,j,ilabile,icarbon)

          ! The conditions do not support growth
          ELSE

             bm_alloc_tot(ipts,j) = zero
             
          ENDIF

          IF (printlev.GE.3) THEN
             IF(ipts == test_grid .AND. j == test_pft)THEN
                WRITE(numout,*) 'allocatable carbon:'
                WRITE(numout,*) 'bm_alloc_tot(test_grid,test_pft)',bm_alloc_tot(ipts,j)
             ENDIF
          ENDIF

          !! 3.10 Maintenance respiration
!JC MOD050 suppress resp_maint of grasses when detecting dormancy
          tmp_resp_maint_part(ipts,j,:) = resp_maint_part(ipts,j,:)
          IF (.NOT. is_tree(j) .AND. natural(j)) THEN
            ! suppression starts when 1) no gpp
            ! 2) marginal leaves (i.e., senescence almost completed)
            ! 3) senescence have been carried out for more than days_senesc_crit days 
            ! (due to drought or low T)  
            IF (gpp_week(ipts,j) .LT. SUM(tmp_resp_maint_part(ipts,j,:)) .AND. &
                (biomass(ipts,j,ileaf,icarbon) .LE. (lai_initmin(j) / 2.)/sla_calc(ipts,j) .OR. &
                days_senescence(ipts,j) .GT. days_senesc_crit)) THEN

              ! just suppress all maintence respiration for dormant period
                tmp_resp_maint_part(ipts,j,iroot) = tmp_resp_maint_part(ipts,j,iroot) * &
                   coeff_suppress_resp
                tmp_resp_maint_part(ipts,j,ileaf) = tmp_resp_maint_part(ipts,j,ileaf) * &
                   coeff_suppress_resp
                tmp_resp_maint_part(ipts,j,isapabove) = tmp_resp_maint_part(ipts,j,isapabove) * &
                   coeff_suppress_resp
                tmp_resp_maint_part(ipts,j,ifruit) = tmp_resp_maint_part(ipts,j,ifruit) * &
                   coeff_suppress_resp

            ENDIF
          ENDIF

          !  First, total maintenance respiration for the whole plant is calculated by 
          !  summing maintenance respiration of the different plant compartments.
          !  This simply recalculates the maintenance respiration from stomate_resp.f90   
          !  Maintenance respiration of the different plant parts is calculated in 
          !  stomate_resp.f90 as a function of the plant's temperature, the long term 
          !  temperature and plant coefficients:
          !  The unit of ::resp_maint is gC m-2 dt-1
          resp_maint(ipts,j) = resp_maint(ipts,j) + SUM(tmp_resp_maint_part(ipts,j,:))
           
          ! Following the calculation of hourly maintenance respiration, verify that 
          ! the PFT has not been killed after calcul of resp_maint_part in stomate.
          ! Can this generaly calculated ::resp_maint be use under the given
          ! conditions? Surpress the respiration for deciduous
          !  PFTs as long as they haven't carried leaves at least once. When 
          !  starting from scratch there is no budburst in the first year because 
          !  the longterm phenological parameters are not initialized yet. If
          !  not surpressed respiration consumes all the reserves before the PFT
          !  can start growing. The code would establish a new PFT but it was
          !  decided to surpress this respiration because it has no physiological
          !  bases.
          IF (ind(ipts,j) .GT. min_stomate .AND. &
               rue_longterm(ipts,j) .NE. un) THEN

             !+++CHECK+++
             ! Can the calculated maintenance respiration be used ? Or 
             ! does it has to be adjusted for special cases. Maintenance 
             ! respiration should be positive. In case it is very low, use 20%
             ! (::maint_from_labile) of the active labile carbon pool (gC m-2 dt-1)
             ! resp_maint(ipts,j) = MAX(zero, MAX(maint_from_labile * gtemp * 
             ! biomass(ipts,j,ilabile,icarbon), resp_maint(ipts,j)))
         
             ! Calculate resp_maint for the labile pool as well, no need to have the
             ! above threshold. Make sure resp_maint is not zero
             resp_maint(ipts,j) = MAX(zero, resp_maint(ipts,j))
             !+++++++++++

             ! Phenological growth makes use of the reserves. Some carbon needs to remain
             ! to support the growth, hence, respiration will be limited. In this case 
             ! resp_maint ((gC m-2 dt-1) should not be more than 80% (::maint_from_gpp) 
             ! of the GPP (gC m-2 s-1) 

          ELSE

             ! No plants, no respiration
             resp_maint(ipts,j) = zero

          ENDIF

          ! The calculation of ::resp_maint is solely based on the demand i.e.
          ! given the biomass and the condition of the plant, how much should be
          ! respired. It is not sure that this demand can be satisfied i.e. the 
          ! calculated maintenance respiration may exceed the available carbon

          deficit = zero

          IF ( bm_alloc_tot(ipts,j) - resp_maint(ipts,j) .LT. zero ) THEN

             deficit = bm_alloc_tot(ipts,j) - resp_maint(ipts,j)

             ! The deficit is less than the carbon reserve
             IF (-deficit .LE. biomass(ipts,j,icarbres,icarbon)) THEN

                ! Pay the deficit from the reserve pool
                biomass(ipts,j,icarbres,icarbon) = &
                     biomass(ipts,j,icarbres,icarbon) + deficit

                ! add the defict to bm_alloc_tot 
                bm_alloc_tot(ipts,j) = bm_alloc_tot(ipts,j) - deficit

             ELSE
                ! Not enough carbon to pay the deficit, the individual 
                ! is going to die at the end of this day

                ! truncate deficit to what is available in reserve
                deficit = - biomass(ipts,j,icarbres,icarbon)
                ! deplete the reserve pool
                biomass(ipts,j,icarbres,icarbon) = zero
                ! add the deficit to bm_alloc_tot
                bm_alloc_tot(ipts,j) = bm_alloc_tot(ipts,j) - deficit
                ! truncate the maintenance respiration to the available carbon
                resp_maint(ipts,j) = bm_alloc_tot(ipts,j)
             ENDIF


          ENDIF

          ! Final ::resp_maint is known
          bm_alloc_tot(ipts,j) = bm_alloc_tot(ipts,j) - resp_maint(ipts,j)

          IF(ld_alloc .AND. ipts == test_grid .AND. j == test_pft)THEN
             WRITE(numout,*) "resp_maint ", resp_maint(ipts,j)
          ENDIF

          biomass(ipts,j,ilabile,icarbon) = MAX(biomass(ipts,j,ilabile,icarbon) - &
                                             (resp_maint(ipts,j) + deficit),zero)

          !  When starting from scratch there is no budburst in the first year because 
          !  the longterm phenological parameters are not initialized yet. If
          !  not surpressed respiration consumes all the reserves before the PFT
          !  can start growing. The code would establish a new PFT but it was
          !  decided to surpress this respiration because it has no physiological
          !  bases.
          IF (ind(ipts,j) .GT. min_stomate .AND. &
               rue_longterm(ipts,j) .NE. un) THEN
             
             frac_growthresp_dyn = frac_growthresp(j)
          ELSE
             frac_growthresp_dyn = zero
          ENDIF

          !! 3.11 Growth respiration
          !  Calculate total growth respiration and update allocatable carbon
          !  Growth respiration is a tax on productivity, not actual allocation
          !  Total growth respiration has be calculated before the allocation 
          !  takes place because the allocation itself is not linear. After 
          !  the allocation has been calculated, growth respiration can be 
          !  calculated for each biomass component separatly. The unit of
          !  resp_growth is gC m-2 dt-1. Surpress the respiration for deciduous
          !  PFTs as long as they haven't carried leaves at least once.
          resp_growth(ipts,j) = frac_growthresp_dyn * MAX(zero, gpp_daily(ipts,j)*dt)


          ! The calculation of ::resp_growth is solely based on the demand
          ! i.e.
          ! given the biomass and the condition of the plant, how much should be
          ! respired. It is not sure that this demand can be satisfied i.e. the
          ! calculated maintenance respiration may exceed the available carbon

          deficit = zero

          IF ( bm_alloc_tot(ipts,j) - resp_growth(ipts,j) .LT. zero ) THEN
    
             deficit = bm_alloc_tot(ipts,j) - resp_growth(ipts,j)
    
             ! The deficit is less than the carbon reserve
             IF (-deficit .LE. biomass(ipts,j,icarbres,icarbon)) THEN
    
                ! Pay the deficit from the reserve pool
                biomass(ipts,j,icarbres,icarbon) = &
                     biomass(ipts,j,icarbres,icarbon) + deficit
                bm_alloc_tot(ipts,j) = bm_alloc_tot(ipts,j) - deficit
    
             ELSE
    
                ! truncate deficit to what is available in reserve
                deficit = - biomass(ipts,j,icarbres,icarbon)
                ! deplete the reserve pool
                biomass(ipts,j,icarbres,icarbon) = zero
                ! add the deficit to bm_alloc_tot
                bm_alloc_tot(ipts,j) = bm_alloc_tot(ipts,j) - deficit
                ! truncate the maintenance respiration to the available carbon
                resp_growth(ipts,j) = bm_alloc_tot(ipts,j)
    
             ENDIF
    
          ENDIF

          ! as the code below is not yet cleaned from frac_growthresp_dyn, we set frac_growthresp_dyn to zero
          frac_growthresp_dyn = zero

          ! Final ::resp_maint is known
          bm_alloc_tot(ipts,j) = bm_alloc_tot(ipts,j) - resp_growth(ipts,j)     

          biomass(ipts,j,ilabile,icarbon) = MAX(biomass(ipts,j,ilabile,icarbon) - &
                                             (resp_growth(ipts,j) + deficit),zero)

          IF (printlev.GE.3) THEN
             IF(ipts == test_grid .AND. j == test_pft)THEN
                WRITE(numout,*) 'resp_maint(test_grid,test_pft)',resp_maint(test_grid,test_pft)
                WRITE(numout,*) 'resp_maint deduced from bm_alloc_tot:'
                WRITE(numout,*) 'bm_alloc_tot(test_grid,test_pft)',bm_alloc_tot(ipts,j)
                WRITE(numout,*) 'deficit',deficit
                WRITE(numout,*) 'resp_maint and deficit deduced from labile:'
                WRITE(numout,*) 'biomass(test_grid,test_pft,ilabile,icarbon)',biomass(ipts,j,ilabile,icarbon)
             ENDIF
          ENDIF

          !! 3.12 Distribute stand level ilabile and icarbres at the tree level
          !  The labile and carbres pools are calculated at the stand level but
          !  are then redistributed at the tree level. This has the advantage 
          !  that biomass and circ_class_biomass have the same dimensions for 
          !  nparts which comes in handy when phenology and mortality are
          !  calculated.
          IF ( is_tree(j) ) THEN

             ! Reset to zero to enable a loop over nparts
             circ_class_biomass(ipts,j,:,ilabile,:) = zero
             circ_class_biomass(ipts,j,:,icarbres,:) = zero

             ! Distribute labile and reserve pools over the circumference classes 
             DO m = 1,nelements

                ! Total biomass across parts and circumference classes
                temp_total_biomass = zero

                DO l = 1,ncirc 

                   DO k = 1,nparts

                      temp_total_biomass = temp_total_biomass + &
                           circ_class_biomass(ipts,j,l,k,icarbon) * circ_class_n(ipts,j,l)

                   ENDDO

                ENDDO

                ! Total biomass across parts but for a specific circumference class
                DO l = 1,ncirc

                   temp_class_biomass = zero

                   DO k = 1,nparts

                      temp_class_biomass = temp_class_biomass + &
                           circ_class_biomass(ipts,j,l,k,icarbon) * circ_class_n(ipts,j,l)

                   ENDDO
    
                   IF (temp_total_biomass .NE. zero) THEN

                      ! Share of this circumference class to the total biomass                      
                      temp_share = temp_class_biomass / temp_total_biomass

                      ! Allocation of ilabile at the tree level (gC tree-1)
                      circ_class_biomass(ipts,j,l,ilabile,m) = temp_share * &
                           biomass(ipts,j,ilabile,m) / circ_class_n(ipts,j,l)

                      ! Allocation of icarbres at the tree level (gC tree-1)
                      circ_class_biomass(ipts,j,l,icarbres,m) = temp_share * &
                           biomass(ipts,j,icarbres,m) / circ_class_n(ipts,j,l)
                      IF((m .EQ. initrogen).OR.(m .EQ. iphosphorus)) THEN
                         circ_class_biomass(ipts,j,l,ileaf,m) = temp_share * &
                           biomass(ipts,j,ileaf,m) / circ_class_n(ipts,j,l)
                         circ_class_biomass(ipts,j,l,isapabove,m) = temp_share * &
                           biomass(ipts,j,isapabove,m) / circ_class_n(ipts,j,l)
                         circ_class_biomass(ipts,j,l,isapbelow,m) = temp_share * &
                           biomass(ipts,j,isapbelow,m) / circ_class_n(ipts,j,l)
                         circ_class_biomass(ipts,j,l,iheartabove,m) = temp_share * &
                           biomass(ipts,j,iheartabove,m) / circ_class_n(ipts,j,l)
                         circ_class_biomass(ipts,j,l,iheartbelow,m) = temp_share * &
                           biomass(ipts,j,iheartbelow,m) / circ_class_n(ipts,j,l)
                         circ_class_biomass(ipts,j,l,iroot,m) = temp_share * &
                           biomass(ipts,j,iroot,m) / circ_class_n(ipts,j,l)
                      ENDIF


                   ELSE

                      circ_class_biomass(ipts,j,l,ilabile,m) = zero
                      circ_class_biomass(ipts,j,l,icarbres,m) = zero
                      IF((m .EQ. initrogen).OR.(m .EQ. iphosphorus)) THEN
                         circ_class_biomass(ipts,j,l,ileaf,m) = zero
                         circ_class_biomass(ipts,j,l,isapabove,m) = zero
                         circ_class_biomass(ipts,j,l,isapbelow,m) = zero
                         circ_class_biomass(ipts,j,l,iheartabove,m) = zero
                         circ_class_biomass(ipts,j,l,iheartbelow,m) = zero
                         circ_class_biomass(ipts,j,l,iroot,m) = zero
                      ENDIF

                   ENDIF

                ENDDO ! ncirc

             ENDDO  ! nelements

          ! Grasses and crops
          ELSE

             DO m = 1,nelements

                ! synchronize biomass and circ_class_biomass
                IF (ind(ipts,j) .GT. zero) THEN

                   circ_class_biomass(ipts,j,1,:,m) = biomass(ipts,j,:,m) / ind(ipts,j)

                ELSE

                   circ_class_biomass(ipts,j,1,:,m) = zero

                ENDIF

             ENDDO

          ENDIF ! is_tree

       ENDDO ! pnts

 !! 5. Allometric allocation

       DO ipts = 1, npts

          !!  5.1 Initialize allocated biomass pools
          f_alloc(ipts,j,:) = zero
          Cl_inc(:)         = zero
          Cs_inc(:)         = zero
          Cr_inc(:)         = zero
          Cf_inc(:)         = zero
          Cl_incp(:)        = zero
          Cs_incp(:)        = zero
          Cr_incp(:)        = zero 
          Cs_inc_est(:)     = zero
          Cl_target(:)      = zero
          Cr_target(:)      = zero
          Cs_target(:)      = zero
 
          !! 5.2 Calculate allocated biomass pools for trees

!JC tree allocated biomass
          !! 5.2.1 Stand to tree allocation rule of Deleuze & Dhote
          IF ( is_tree(j) .AND. bm_alloc_tot(ipts,j) .GT. min_stomate ) THEN

             !  Basal area at the tree level (m2 tree-1)
             circ_class_ba_eff(:) = wood_to_ba_eff(circ_class_biomass(ipts,j,1,:,icarbon),j)
             circ_class_circ_eff(:) = 2 * pi * SQRT(circ_class_ba_eff(:)/pi)

             ! According to equation (-) in Bellasen et al 2010. 
             ! ln(sigmas) = a_sig * ln(circ_med) + b_sig
             ! sigmas = exp(a_sig*log(median(circ_med))+b_sig);
             ! However, in the code (sapiens_forestry.f90) a different expression was used
             ! sigmas = 0.023+0.58*prctile(circ_med,0.05);
             ! Any of these implementations could work but seem to be more suited for 
             ! continues or nearly continuous diameter distributions, say n_circ > 10
             ! For a small number of diameter classes sigma depends on a prescribed
             ! circumference percentile.
             IF (ncirc .GE. 6) THEN          
               
                ! Calculate the median circumference
                DO l = 1,ncirc
                   
                      IF (SUM(circ_class_n(ipts,j,1:l)) .GE. 0.5 * ind(ipts,j)) THEN
                         
                         median_circ  = circ_class_circ_eff(l) - 5 * min_stomate
                      
                         EXIT
                      
                      ENDIF
                   
                ENDDO

                sigma(ipts,j) = deleuze_a(j) + deleuze_b(j) * median_circ
                
             ELSE
                deleuze_p(j) = deleuze_p_all(j)
                ! Search for the X percentile, where X is given by ::deleuze_p
                ! Substract a very small number (5*min_stomate) just to be sure that 
                ! the circ_class will be corectly accounted for in GE or LE statements
                DO l = 1,ncirc
                   
                      IF (SUM(circ_class_n(ipts,j,1:l)) .GE. deleuze_p(j) * ind(ipts,j)) THEN
                         
                         sigma(ipts,j) = circ_class_circ_eff(l) - 5 * min_stomate
                      
                         EXIT
                      
                      ENDIF
                   
                ENDDO
                
             ENDIF

          
             !! 5.2 Calculate allocated biomass pools for trees
             !  Only possible if there is biomass to allocate
             !  Use sigma and m_dv to calculate a single coefficient that can be
             !  used in the subsequent allocation scheme.
             circ_class_dba(:) = (circ_class_circ_eff(:) - m_dv(j)*sigma(ipts,j) + &
                  SQRT((m_dv(j)*sigma(ipts,j) + circ_class_circ_eff(:))**2 - &
                  (4*sigma(ipts,j)*circ_class_circ_eff(:)))) / 2

             !! 5.2.2 Scaling factor to convert variables to the individual plant
             !  Allocation is on an individual basis. Stand-level variables need to convert to a 
             !  single individual. Different approach between the DGVM and statitic approach
             IF (ok_dgvm) THEN

                ! The DGVM does currently NOT work with the new allocation, consider this as
                ! placeholder. The original code had two different transformations to 
                ! calculate the scalars. Both could be used but the units will differ.
                ! When fixing the DGVM check which quantities need to be multiplied by scal 
                ! scal = ind(ipts,j) * cn_ind(ipts,j) / veget_max(ipts,j)
                scal(ipts,j) = veget_max(ipts,j) / ind(ipts,j)  
             
             ELSE
             
                ! circ_class_biomass contain the data at the tree level
                ! no conversion required
                scal(ipts,j) = 1.
             
             ENDIF
                      

             !! 5.2.3 Current biomass pools per tree (gC tree^-1) 
             ! We will have different trees so this has to be calculated from the 
             ! diameter relationships            
             Cs(:) = ( circ_class_biomass(ipts,j,:,isapabove,icarbon) + &
                  circ_class_biomass(ipts,j,:,isapbelow,icarbon) ) * scal(ipts,j)
             Cr(:) = circ_class_biomass(ipts,j,:,iroot,icarbon) * scal(ipts,j)
             Cl(:) = circ_class_biomass(ipts,j,:,ileaf,icarbon) * scal(ipts,j)
             Ch(:) = ( circ_class_biomass(ipts,j,:,iheartabove,icarbon) + &
                  circ_class_biomass(ipts,j,:,iheartbelow,icarbon) ) * scal(ipts,j)
 
             ! Total amount of carbon that needs to ba allocated (::bm_alloc_tot). 
             ! bm_alloc_tot is in gC m-2 day-1. At 1 m2 there are ::ind number of 
             ! trees. We calculate the allocation for ::ncirc trees. Hence b_inc_tot 
             ! needs to be scaled in the allocation routines. For all cases were 
             ! allocation takes place for a single circumference class, scaling 
             ! could be done before the allocation. In the ordinary allocation 
             ! allocation takes place to all circumference classes at the same time. 
             ! Hence scaling takes place in that step for consistency we scale during 
             ! allocation. Note that b_inc (the carbon allocated to an individual 
             ! circumference class cannot be estimates at this point.
             b_inc_tot = bm_alloc_tot(ipts,j)
             

             !! 5.2.4 C-allocation for trees
             !  The mass conservation equations are detailed in the header of this subroutine.
             !  The scheme assumes a functional relationships between leaves, sapwood and 
             !  roots. When carbon is added to the leaf biomass pool, an increase in the root
             !  biomass is to be expected to sustain water transport from the roots to the 
             !  leaves. Also sapwood is needed to sustain this water transport and to support
             !  the leaves.
             DO l = 1,ncirc 

                !! 5.2.4.1 Calculate tree height
                circ_class_height_eff(l) = pipe_tune2(j)* & 
                     (4/pi*circ_class_ba_eff(l))**(pipe_tune3(j)/2)


                !! 5.2.4.2 Do the biomass pools respect the pipe model?
                !  Do the current leaf, sapwood and root components respect the allometric 
                !  constraints? Due to plant phenology it is possible that we have too much 
                !  sapwood compared to the leaf and root mass (i.e. in early spring). 
                !  Calculate the optimal root and leaf mass, given the current wood mass 
                !  by using the basic allometric relationships. Calculate the optimal sapwood
                !  mass as a function of the current leaf and root mass.
                Cl_target(l) = MAX( KF(ipts,j) * Cs(l) / circ_class_height_eff(l), &
                     Cr(l) * LF(ipts,j) , Cl(l))
                Cr_target(l) = MAX( Cl_target(l) / LF(ipts,j), &
                     Cs(l) * KF(ipts,j) / LF(ipts,j) / circ_class_height_eff(l) , Cr(l))
                Cs_target(l) = MAX( Cl(l) / KF(ipts,j) * circ_class_height_eff(l), &
                     Cr(l) * LF(ipts,j) / KF(ipts,j) * circ_class_height_eff(l) , Cs(l))

                !---TEMP---
                IF (j.EQ.test_pft .AND. ld_alloc .AND. ipts==test_grid) THEN
                   WRITE(numout,*) 'bm_alloc_tot, ', bm_alloc_tot(ipts,j)
                   WRITE(numout,*) 'Does the tree needs reshaping? Class: ',l
                   WRITE(numout,*) 'circ_class_height_eff, ', circ_class_height_eff(l)
                   WRITE(numout,*) 'KF, LF, ', KF(ipts,j), LF(ipts,j)
                   WRITE(numout,*) 'Cl_target-Cl, ', Cl_target(l)-Cl(l), Cl_target(l), Cl(l)
                   WRITE(numout,*) 'Cs_target-Cs, ', Cs_target(l)-Cs(l), Cs_target(l), Cs(l)
                   WRITE(numout,*) 'Cr_target-Cr, ', Cr_target(l)-Cr(l), Cr_target(l), Cr(l)
                ENDIF
                !----------

                !! 5.2.4.2 Check dimensions of the trees
                ! If Cs = Cs_target then ba and height are correct, else calculate 
                ! the correct dimensions
                IF ( Cs_target(l) - Cs(l) .GT. min_stomate ) THEN

                   ! If Cs = Cs_target then dia and height are correct. However, 
                   ! if Cl = Cl_target or Cr = Cr_target then dia and height need 
                   ! to be re-estimated. Cs_target should satify the relationship 
                   ! Cl/Cs = KF/height where height is a function of Cs_target
                   ! <=> (KF*Cs_target)/(pipe_tune2*(Cs_target+Ch)/pi/4)**&
                   ! (pipe_tune3/(2+pipe_tune3)) = Cl_target. Search Cs needed to 
                   ! sustain the max of Cl or Cr. Search max of Cl and Cr first 
                   Cl_target(l) = MAX(Cl(l), Cr(l)*LF(ipts,j))

                   ! Recalculate height and ba from the correct Cs_target
                   circ_class_height_eff(l) = Cs_target(l)*KF(ipts,j)/Cl_target(l)
                   circ_class_ba_eff(l) = pi/4*(circ_class_height_eff(l)/ & 
                        pipe_tune2(j))**(2/pipe_tune3(j))
                   Cl_target(l) = KF(ipts,j) * Cs_target(l) / circ_class_height_eff(l)
                   Cr_target(l) = Cl_target(l) / LF(ipts,j)

                ENDIF

                !---TEMP---
                IF (j.EQ.test_pft .AND. ld_alloc .AND. ipts==test_grid) THEN
                   WRITE(numout,*) 'height_fin, ba_fin, ', circ_class_height_eff(:), &
                        circ_class_ba_eff(:)
                   WRITE(numout,*) 'Cl_target, Cs_target, Cr_target, ', Cl_target(:), &
                        Cs_target(:), Cr_target(:)
                   WRITE(numout,*) 'New target values'
                   WRITE(numout,*) 'Cl_target-Cl, ', Cl_target(l)-Cl(l), Cl_target(l), Cl(l)
                   WRITE(numout,*) 'Cs_target-Cs, ', Cs_target(l)-Cs(l), Cs_target(l), Cs(l)
                   WRITE(numout,*) 'Cr_target-Cr, ', Cr_target(l)-Cr(l), Cr_target(l), Cr(l)
                ENDIF
                !-----------

             ENDDO

             ! The step estimate is used to linearalize the diameter vs height relationship.
             ! Use a prior to distribute b_inc_tot over the individual trees. The share of
             ! the total sapwood mass is used as a prior. Subsequently, estimate the change in 
             ! diameter by assuming all the available C for allocation will be used in Cs.  
             ! Hence, this represents the maximum possible diameter increase. It was not tested
             ! whether this is the best prior but it seems to work OK although it often results
             ! in very small (1e-8) negative values, with even more rare 1e-6 negative values. 
             ! A C-balance closure check could reveal 
             ! whether this is a real issue and requires to change the prior or not.
             ! Calculate the linear slope (::s) of the relationship between ba and h as 
             ! (1) s = (ba2-ba)/(height2-height). 
             ! The goal is to approximate the ba2 that
             ! is predicted through the non-linear ordinary allocation approach, as this will 
             ! keep the trees in allometric balance. In the next time step, allometric 
             ! balance is recalculated and can be corrected through the so-called phenological 
             ! growth; hence, small deviations resulting from the linearization will not 
             ! accumulate with time.
             ! Note that ba2 = ba + delta_ba and that height and ba are related as
             ! (2) height = k2*(4*ba/pi)**(k3/2)
             ! At this stage the only information we have is that there is b_inc_tot (gC m-2) 
             ! available for allocation. There are two obvious approximations both making use
             ! of the same assumption, i.e. that for the initial estimate of delta_ba height is
             ! constant. The first approximation is crude and assumes that all the available C
             ! is used in Cs_inc (thus Cs_inc = b_inc_tot / ind ). The second approximation,
             ! implemented here, makes use of the allometric rules and thus accounts for the
             ! knowledge that allocating one unit the sapwood comes with a cost in leaves and 
             ! roots thus:
             ! b_inc_temp = Cs_inc+Cl_inc+Cr_inc
             ! (3) <=> b_inc_temp ~= (Cs_inc_est+Cs) + KF*(Cs_inc_est+Cs)/H + ...
             !    KF/LF*(Cs_inc_est+Cs)/H - Cs - Cl - Cr
             ! b_inc_temp is the amount of carbon that can be allocated to each diameter class.
             ! However, only the total amount i.e. b_inc_tot is known. Total allocatable carbon
             ! is distributed over the different diameter classes proportional to their share
             ! of the total wood biomass. Divide by circ_class_n to get the correct units 
             ! (gC tree-1)
             ! (4)  b_inc_temp ~= b_inc_tot / circ_class_n * (circ_class_n * ba**(1+k3)) / ...
             !    sum(circ_class_n * ba**(1+k3))
             ! By substituting (4) in (3) an expression is obtained to approximate the carbon
             ! that will be allocated to sapwood growth per diameter class ::Cs_inc_est. This
             ! estimate is then used to calculate delta_ba (called ::step) as
             ! step = (Cs+Ch+Cs_inc_set)/(tree_ff*pipe_density*height) - ba where height is 
             ! calculated from (2) after replacing ba by ba+delta_ba
             
             !+++CHECK+++
             !Alternative decribed in the documentation - most complete
             Cs_inc_est(:) = ( b_inc_tot / circ_class_n(ipts,j,:) * &
                  (circ_class_n(ipts,j,:) * circ_class_ba_eff(:)**(un+pipe_tune3(j))) / &
                  (SUM(circ_class_n(ipts,j,:) * circ_class_ba_eff(:)**(un+pipe_tune3(j)))) + &
                  Cs(:) + Cl(:) + Cr(:)) * circ_class_height_eff(:) / &
                  (circ_class_height_eff(:) + KF(ipts,j) + KF(ipts,j)/LF(ipts,j)) - Cs(:)

             !Keep it simple
             !Cs_inc_est(:) = ( b_inc_tot / circ_class_n(ipts,j,:) * &
             !     (circ_class_n(ipts,j,:) * circ_class_ba_eff(:)**(un+pipe_tune3(j))) / &
             !     (SUM(circ_class_n(ipts,j,:) * circ_class_ba_eff(:)**(un+pipe_tune3(j)))))
             !+++++++++++

             step(:) = ((Ch(:)+Cs(:)+Cs_inc_est(:)) / (tree_ff(j)*pipe_density(j)* &
                  circ_class_height_eff(:))) - circ_class_ba_eff(:)

             !++++ CHECK ++++++
             ! It can happen that step is equal to zero sometimes.  I'm not sure why, but 
             ! there was a case where it was nonzero for circ classes 1 and 3, and zero 
             ! for 2.  This causes s to be zero and provokes a divide by zero error later 
             ! on.  What if we make it not zero?  This might cause a small mass balance
             ! error for this timestep, but I would rather have that than getting an 
             ! infinite biomass, which is what happened in the other case.  These limits 
             ! are arbitrary and adjusted by hand.  If the output file doesn't show this 
             ! warning very often, I think we're okay, since the amount of carbon is really
             ! small.
             DO l=1,ncirc
                IF(step(l) .LT. min_stomate*0.01 .AND. step(l) .GT. zero)THEN
                   step(l)=min_stomate*0.02
                   IF (ld_alloc) THEN
                      WRITE(numout,*) 'WARNING: Might cause mass balance problems in fun_all, position 1'
                      WRITE(numout,*) 'WARNING: ips,j ',ipts,j
                   END IF
                ELSEIF(step(l) .GT. -min_stomate*0.01 .AND. step(l) .LT. zero)THEN
                   step(l)=-min_stomate*0.02
                   IF (ld_alloc) THEN
                      WRITE(numout,*) 'WARNING: Might cause mass balance problems in fun_all, position 2'
                      WRITE(numout,*) 'WARNING: ips,j ',ipts,j
                   END IF
                ENDIF
             ENDDO
             !+++++++++++++++++
             s(:) = step(:)/(pipe_tune2(j)*(4.0_r_std/pi*(circ_class_ba_eff(:)+step(:)))**&
                  (pipe_tune3(j)/deux) - &
                  pipe_tune2(j)*(4.0_r_std/pi*circ_class_ba_eff(:))**(pipe_tune3(j)/deux))

             !! 5.2.4.3 Phenological growth
             !  Phenological growth and reshaping of the tree in line with the pipe model. 
             !  Turnover removes C from the different plant components but at a 
             !  component-specific rate, as such the allometric constraints are distorted 
             !  at every time step and should be restored before ordinary growth can 
             !  take place
             l = ncirc
             DO WHILE ( (l .GT. zero) .AND. (b_inc_tot .GT. min_stomate) )

                !! 5.2.4.3.1 The available wood can sustain the available leaves and roots
                !  Calculate whether the wood is in allometric balance. The target values 
                !  should always be larger than the current pools so the use of ABS is 
                !  redundant but was used to be on the safe side (here and in the rest 
                !  of the module) as it could help to find logical flaws.
                IF ( ABS(Cs_target(l) - Cs(l)) .LT. min_stomate ) THEN

                   ! Use the difference between the target and the actual to
                   ! ensure mass balance closure because l times a values 
                   ! smaller than min_stomate can still add up to a value 
                   ! exceeding min_stomate.
                   Cs_incp(l) = MAX(zero, Cs_target(l) - Cs(l))

                   ! Enough leaves and wood, only grow roots
                   IF ( ABS(Cl_target(l) - Cl(l))  .LT. min_stomate ) THEN

                      ! Allocate at the tree level to restore allometric balance
                      ! Some carbon may have been used for Cs_incp and Cl_incp
                      ! adjust the total allocatable carbon
                      Cl_incp(l) = MAX(zero, Cl_target(l) - Cl(l))
                      Cr_incp(l) = MAX( MIN(b_inc_tot / circ_class_n(ipts,j,l) - &
                           Cs_incp(l) - Cl_incp(l), Cr_target(l) - Cr(l)), zero )

                      ! Write debug comments to output file
                      IF ((j.EQ.test_pft .AND. ld_alloc) .OR. ld_warn) THEN
                         CALL comment(npts, Cl_target, Cl, Cs_target, Cs, Cr_target, Cr, &
                              delta_ba, ipts, j, l, b_inc_tot, Cl_incp, Cs_incp, Cr_incp, &
                              KF, LF, Cl_inc, Cs_inc, Cr_inc, Cf_inc, grow_wood, &
                              circ_class_n, ind, 1)
                      ENDIF

                   ! Sufficient wood and roots, allocate C to leaves
                   ELSEIF ( ABS(Cr_target(l) - Cr(l)) .LT. min_stomate ) THEN

                      ! Allocate at the tree level to restore allometric balance
                      ! Some carbon may have been used for Cs_incp and Cr_incp
                      ! adjust the total allocatable carbon
                      Cr_incp(l) = MAX(zero, Cr_target(l) - Cr(l))
                      Cl_incp(l) = MAX( MIN(b_inc_tot / circ_class_n(ipts,j,l) - &
                           Cs_incp(l) - Cr_incp(l), Cl_target(l) - Cl(l)), zero )

                      ! Write debug comments to output file
                      IF ((j.EQ.test_pft .AND. ld_alloc) .OR. ld_warn) THEN
                         CALL comment(npts, Cl_target, Cl, Cs_target, Cs, Cr_target, Cr, &
                              delta_ba, ipts, j, l, b_inc_tot, Cl_incp, Cs_incp, Cr_incp, &
                              KF, LF, Cl_inc, Cs_inc, Cr_inc, Cf_inc, &
                              grow_wood, circ_class_n, ind, 2)
                      ENDIF
                      
                   ! Both leaves and roots are needed to restore the allometric relationships
                   ELSEIF ( ABS(Cl_target(l) - Cl(l)) .GT. min_stomate .AND. &
                        ABS(Cr_target(l) - Cr(l)) .GT. min_stomate ) THEN                 

                      ! Allocate at the tree level to restore allometric balance
                      !  The equations can be rearanged and written as
                      !  (i) b_inc = Cl_inc + Cr_inc
                      !  (ii) Cr_inc = (Cl_inc+Cl)/LF - Cr
                      !  Substitue (ii) in (i) and solve for Cl_inc
                      !  <=> Cl_inc = (LF*(b_inc+Cr)-Cl)/(1+LF)
                      Cl_incp(l) = MIN( ((LF(ipts,j) * ((b_inc_tot/circ_class_n(ipts,j,l) - &
                           Cs_incp(l)) + Cr(l))) - Cl(l)) / & 
                           (1 + LF(ipts,j)), Cl_target(l) - Cl(l) )
                      Cr_incp(l) = MIN ( ((Cl_incp(l) + Cl(l)) / LF(ipts,j)) - Cr(l), &
                           Cr_target(l) - Cr(l))

                      ! The imbalance between Cr and Cl can be so big that (Cl+Cl_inc)/LF 
                      ! is still less then the available root carbon (observed!). This would 
                      ! result in a negative Cr_incp
                      IF ( Cr_incp(l) .LT. zero ) THEN

                         Cl_incp(l) = MIN( b_inc_tot/circ_class_n(ipts,j,l) - Cs_incp(l), &
                              Cl_target(l) - Cl(l) )
                         Cr_incp(l) = (b_inc_tot/circ_class_n(ipts,j,l)) - Cs_incp(l) - &
                              Cl_incp(l)

                      ELSEIF (Cl_incp(l) .LT. zero) THEN

                         Cr_incp(l) = MIN( b_inc_tot/circ_class_n(ipts,j,l) - Cs_incp(l), &
                              Cr_target(l) - Cr(l) )
                         Cl_incp(l) = (b_inc_tot/circ_class_n(ipts,j,l)) - &
                              Cs_incp(l) - Cr_incp(l)

                      ENDIF                          

                      ! Write debug comments to output file
                      IF ((j.EQ.test_pft .AND. ld_alloc) .OR. ld_warn) THEN
                         CALL comment(npts, Cl_target, Cl, Cs_target, Cs, Cr_target, Cr, &
                              delta_ba, ipts, j, l, b_inc_tot, Cl_incp, Cs_incp, Cr_incp, &
                              KF, LF, Cl_inc, Cs_inc, Cr_inc, Cf_inc, &
                              grow_wood, circ_class_n, ind, 3)
                      ENDIF      

                   ELSE

                      WRITE(numout,*) 'Exc 1-3: unexpected exception' 
                      IF(ld_stop)THEN
                         CALL ipslerr_p (3,'growth_fun_all',&
                              'Exc 1-3: unexpected exception','','')
                      ENDIF

                   ENDIF

                !! 5.2.4.3.2 Enough leaves to sustain the wood and roots
                ELSEIF ( ABS(Cl_target(l) - Cl(l)) .LT. min_stomate ) THEN

                   ! Use the difference between the target and the actual to
                   ! ensure mass balance closure because l times a values 
                   ! smaller than min_stomate can still add up to a value 
                   ! exceeding min_stomate.
                   Cl_incp(l) = MAX(zero, Cl_target(l) - Cl(l))

                   ! Enough leaves and wood, only grow roots
                   ! This duplicates Exc 1 and these lines should never be called 
                   IF ( ABS(Cs_target(l) - Cs(l)) .LT. min_stomate ) THEN

                      ! Allocate at the tree level to restore allometric balance
                      ! Some carbon may have been used for Cs_incp and Cl_incp
                      ! adjust the total allocatable carbon
                      Cs_incp(l) = MAX(zero, ABS(Cs_target(l) - Cs(l)))
                      Cr_incp(l) = MAX( MIN(b_inc_tot/circ_class_n(ipts,j,l) - &
                           Cl_incp(l) - Cs_incp(l), Cr_target(l) - Cr(l)), zero )

                      ! Write debug comments to output file
                      IF ((j.EQ.test_pft .AND. ld_alloc) .OR. ld_warn) THEN
                         CALL comment(npts, Cl_target, Cl, Cs_target, Cs, Cr_target, Cr, &
                              delta_ba, ipts, j, l, b_inc_tot, Cl_incp, Cs_incp, Cr_incp, &
                              KF, LF, Cl_inc, Cs_inc, Cr_inc, Cf_inc, grow_wood, &
                              circ_class_n, ind, 4)
                      ENDIF

                   ! Enough leaves and roots. Need to grow sapwood to support the available 
                   ! canopy and roots
                   ELSEIF ( ABS(Cr_target(l) - Cr(l)) .LT. min_stomate ) THEN

                      ! In truth, there might be a little root carbon to allocate here,
                      ! since min_stomate is not equal to zero.  If there is
                      ! enough of this small carbon in every circ class, and there
                      ! are enough circ classes, ordinary allocation will be skipped
                      ! below and we might try to force allocation, which is silly
                      ! if the different in the root masses is around 1e-8. This
                      ! means we will allocate a tiny amount to the roots to make
                      ! sure they are exactly in balance.                  
                      Cr_incp(l) = MAX(zero, ABS(Cr_target(l) - Cr(l)))
                      Cs_incp(l) = MAX( MIN(b_inc_tot/circ_class_n(ipts,j,l) - &
                           Cl_incp(l) - Cr_incp(l), Cs_target(l) - Cs(l)), zero )

                      ! Write debug comments to output file
                      IF ((j.EQ.test_pft .AND. ld_alloc) .OR. ld_warn) THEN
                         CALL comment(npts, Cl_target, Cl, Cs_target, Cs, Cr_target, Cr, &
                              delta_ba, ipts, j, l, b_inc_tot, Cl_incp, Cs_incp, Cr_incp, &
                              KF, LF, Cl_inc, Cs_inc, Cr_inc, Cf_inc, &
                              grow_wood, circ_class_n, ind, 5)
                      ENDIF                     

                   ! Need both wood and roots to restore the allometric relationships
                   ELSEIF ( ABS(Cs_target(l) - Cs(l) ) .GT. min_stomate .AND. &
                        ABS(Cr_target(l) - Cr(l)) .GT. min_stomate ) THEN

                      ! circ_class_ba_eff and circ_class_height_eff are already calculated
                      ! for a tree in balance. It would be rather complicated to follow
                      ! the allometric rules for wood allocation (implying changes in height 
                      ! and basal area) because the tree is not in balance yet. First try 
                      ! if we can simply satisfy the allocation needs
                      IF (Cs_target(l) - Cs(l) + Cr_target(l) - Cr(l) .LE. &
                           b_inc_tot/circ_class_n(ipts,j,l) - Cl_incp(l)) THEN
                         
                         Cr_incp(l) = Cr_target(l) - Cr(l)
                         Cs_incp(l) = Cs_target(l) - Cs(l)

                      ! Try to satisfy the need for roots
                      ELSEIF (Cr_target(l) - Cr(l) .LE. b_inc_tot/circ_class_n(ipts,j,l) - &
                           Cl_incp(l)) THEN

                         Cr_incp(l) = Cr_target(l) - Cr(l)
                         Cs_incp(l) = b_inc_tot/circ_class_n(ipts,j,l) - &
                              Cl_incp(l) - Cr_incp(l)
                         
                      ! There is not enough use whatever is available
                      ELSE
                         
                         Cr_incp(l) = b_inc_tot/circ_class_n(ipts,j,l) - Cl_incp(l)
                         Cs_incp(l) = zero
                         
                      ENDIF

                      ! Write debug comments to output file
                      IF ((j.EQ.test_pft .AND. ld_alloc) .OR. ld_warn) THEN
                         CALL comment(npts, Cl_target, Cl, Cs_target, Cs, Cr_target, Cr, &
                              delta_ba, ipts, j, l, b_inc_tot, Cl_incp, Cs_incp, Cr_incp, &
                              KF, LF, Cl_inc, Cs_inc, Cr_inc, Cf_inc, grow_wood, &
                              circ_class_n, ind, 6)
                      ENDIF   

                   ELSE

                      WRITE(numout,*) 'Exc 4-6: unexpected exception'
                      IF(ld_stop)THEN
                         CALL ipslerr_p (3,'growth_fun_all',&
                              'Exc 4-6: unexpected exception','','')
                      ENDIF
                      
                   ENDIF


                !! 5.2.4.3.3 Enough roots to sustain the wood and leaves
                ELSEIF ( ABS(Cr_target(l) - Cr(l)) .LT. min_stomate ) THEN

                   ! Use the difference between the target and the actual to
                   ! ensure mass balance closure because l times a values 
                   ! smaller than min_stomate can still add up to a value 
                   ! exceeding min_stomate.
                   Cr_incp(l) = MAX(zero, Cr_target(l) - Cr(l))

                   ! Enough roots and wood, only grow leaves
                   ! This duplicates Exc 2 and these lines should thus never be called 
                   IF ( ABS(Cs_target(l) - Cs(l)) .LT. min_stomate ) THEN

                      ! Allocate at the tree level to restore allometric balance
                      Cs_incp(l) = MAX(zero, Cs_target(l) - Cs(l))
                      Cl_incp(l) = MAX( MIN(b_inc_tot/circ_class_n(ipts,j,l) - &
                           Cs_incp(l) - Cr_incp(l), &
                           Cl_target(l) - Cl(l)), zero )

                      ! Write debug comments to output file
                      IF ((j.EQ.test_pft .AND. ld_alloc) .OR. ld_warn) THEN
                         CALL comment(npts, Cl_target, Cl, Cs_target, Cs, Cr_target, Cr, &
                              delta_ba, ipts, j, l, b_inc_tot, Cl_incp, Cs_incp, Cr_incp, &
                              KF, LF, Cl_inc, Cs_inc, Cr_inc, Cf_inc, grow_wood, &
                              circ_class_n, ind, 7)
                      ENDIF

                   ! Enough leaves and roots. Need to grow sapwood to support the 
                   ! available canopy and roots. Duplicates Exc. 4 and these lines 
                   ! should thus never be called 
                   ELSEIF ( ABS(Cl_target(l) - Cl(l)) .LT. min_stomate ) THEN

                      ! Allocate at the tree level to restore allometric balance
                      Cl_incp(l) = MAX(zero, Cl_target(l) - Cl(l))
                      Cs_incp(l) = MAX( MIN(b_inc_tot/circ_class_n(ipts,j,l) - &
                           Cr_incp(l) - Cl_incp(l), Cs_target(l) - Cs(l) ), zero )

                      ! Write debug comments to output file
                      IF ((j.EQ.test_pft .AND. ld_alloc) .OR. ld_warn) THEN
                         CALL comment(npts, Cl_target, Cl, Cs_target, Cs, Cr_target, Cr, &
                              delta_ba, ipts, j, l, b_inc_tot, Cl_incp, Cs_incp, Cr_incp, &
                              KF, LF, Cl_inc, Cs_inc, Cr_inc, Cf_inc, grow_wood, &
                              circ_class_n, ind, 8)
                      ENDIF

                   ! Need both wood and leaves to restore the allometric relationships
                   ELSEIF ( ABS(Cs_target(l) - Cs(l)) .GT. min_stomate .AND. &
                        ABS(Cl_target(l) - Cl(l)) .GT. min_stomate ) THEN

                      ! circ_class_ba_eff and circ_class_height_eff are already calculated
                      ! for a tree in balance. It would be rather complicated to follow
                      ! the allometric rules for wood allocation (implying changes in height 
                      ! and basal area) because the tree is not in balance.First try if we 
                      ! can simply satisfy the allocation needs
                      IF (Cs_target(l) - Cs(l) + Cl_target(l) - Cl(l) .LE. &
                           b_inc_tot/circ_class_n(ipts,j,l) - Cr_incp(l)) THEN

                         Cl_incp(l) = Cl_target(l) - Cl(l)
                         Cs_incp(l) = Cs_target(l) - Cs(l)

                      ! Try to satisfy the need for leaves
                      ELSEIF (Cl_target(l) - Cl(l) .LE. b_inc_tot/circ_class_n(ipts,j,l) - &
                           Cr_incp(l)) THEN

                         Cl_incp(l) = Cl_target(l) - Cl(l)
                         Cs_incp(l) = b_inc_tot/circ_class_n(ipts,j,l) - &
                              Cr_incp(l) - Cl_incp(l)

                      ! There is not enough use whatever is available
                      ELSE

                         Cl_incp(l) = b_inc_tot/circ_class_n(ipts,j,l) - Cr_incp(l)
                         Cs_incp(l) = zero

                      ENDIF

                      ! Write debug comments to output file
                      IF ((j.EQ.test_pft .AND. ld_alloc) .OR. ld_warn) THEN
                         CALL comment(npts, Cl_target, Cl, Cs_target, Cs, Cr_target, Cr, &
                              delta_ba, ipts, j, l, b_inc_tot, Cl_incp, Cs_incp, Cr_incp, &
                              KF, LF, Cl_inc, Cs_inc, Cr_inc, Cf_inc, grow_wood, &
                              circ_class_n, ind, 9)
                      ENDIF

                   ELSE

                      WRITE(numout,*) 'Exc 7-9: unexpected exception'
                      IF(ld_stop)THEN
                         CALL ipslerr_p (3,'growth_fun_all',&
                              'Exc 7-9: unexpected exception','','')
                      ENDIF

                   ENDIF

                ! Either Cl_target, Cs_target or Cr_target should be zero
                ELSE

                   ! Something possibly important was overlooked
                   WRITE(numout,*) 'WARNING 4: logical flaw in the phenological allocation, PFT, class: ', j, l
                   WRITE(numout,*) 'WARNING 4: PFT, ipts: ',j,ipts
                   WRITE(numout,*) 'Cs - Cs_target', Cs(l), Cs_target(l)
                   WRITE(numout,*) 'Cl - Cl_target', Cl(l), Cl_target(l)
                   WRITE(numout,*) 'Cr - Cr_target', Cr(l), Cr_target(l)
                   IF(ld_stop)THEN
                       CALL ipslerr_p (3,'growth_fun_all',&
                            'WARNING 4: logical flaw in the phenological allocation','','')
                   ENDIF

                ENDIF

                !! 5.2.4.4 Wrap-up phenological allocation
                IF ( Cl_incp(l) .GE. zero .OR. Cr_incp(l) .GE. zero .OR. &
                     Cs_incp(l) .GE. zero) THEN

                   ! Fake allocation for less messy equations in next case, 
                   ! incp needs to be added to inc at the end
                   Cl(l) = Cl(l) + Cl_incp(l)
                   Cr(l) = Cr(l) + Cr_incp(l)
                   Cs(l) = Cs(l) + Cs_incp(l)
                   b_inc_tot = b_inc_tot - (circ_class_n(ipts,j,l) * &
                        (Cl_incp(l) + Cr_incp(l) + Cs_incp(l)))             

                   ! Something is wrong with the calculations
                   IF (b_inc_tot .LT. -min_stomate) THEN

                      WRITE(numout,*) 'WARNING 5: numerical problem, overspending in phenological allocation'
                      WRITE(numout,*) 'WARNING 5: PFT, ipts: ',j,ipts
                      CALL ipslerr_p (3,'growth_fun_all',&
                           'WARNING 5: numerical problem, overspending in phenological allocation','','')
                   ENDIF

                ELSE

                   ! The code was written such that the increment pools should be 
                   ! greater than or equal to zero. If this is not the case, something 
                   ! fundamental is wrong with the if-then constructs under §5.2.4.3
                   WRITE(numout,*) 'WARNING 6: PFT, ipts: ',j,ipts
                   CALL ipslerr_p (3,'growth_fun_all',&
                        'WARNING 6: numerical problem, one of the increment pools is less than zero','','')
                ENDIF

                ! Set counter for next circumference class
                l = l-1

             ENDDO ! DO WHILE l.GE.1 .AND. b_inc_tot .GT. min_stomate


             !! 5.2.5 Calculate the expected size of the reserve pool 
             !  use the minimum of either (1) 2% of the total sapwood biomass or 
             !  (2) the amount of carbon needed to develop the optimal LAI and "30% of" (DSG)the roots
             !  This reserve pool estimate is only used to decide whether wood should be
             !  grown or not. When really dealing with the reserves the reserve pool is 
             !  recalculated. See further below §7.1. 
           
             ! DSG: this kamikaze behaviour of tress doesn't work well for TeDBF PFT;
             !      at the start of the growing season these trees rely too much
             !      on C from photosynthesis to built leave & roots.  This is not in line with
             !      observations (Hartmann & Trumbore, 2016) and a bit risky for
             !      the average behaviour of an ecosystem.
             !      The minimum  size of the reserve pool should be sufficient to built up
             !      a second set of leaves at the start of the growing season
             !      (see Hartmann & Trumbore, 2016).
             !      For a start, I increase the reserves for roots to 100% (OCN).

             ! DSG: Deciduousness is connected to a higher defoliation
             ! resistancy (Piper and Fajardo (2014)); so they should have higher
             ! storage capacity as evergreen ones (as reflected in the variable
             ! for the upper limit of reserves "deciduous_reserve / evergreen_reserve",
             ! but it is not reflected in the dynamic part of the equation in
             ! the same extent (only after DSGlabXX 100% vs 30% of root mass); I
             ! will increase that further by allowing them 150% of total root
             ! and leaf mass to be stored.

             !DSGlabXX reserve_pool = MIN( 0.02 * ( biomass(ipts,j,isapabove,icarbon) + &
             ! the increase from 2-12% didn't help much as the other limit is
             ! usually hit with a reserve / sapwood ratio of ~5-8%.
             ! we need to allow also to save N for more of the roots 
             reserve_pool = MIN( deciduous_reserve(j) * ( biomass(ipts,j,isapabove,icarbon) + &
                  biomass(ipts,j,isapbelow,icarbon)), &
                  lai_target_longterm(ipts,j)/sla(j)*(1.+root_reserve(j)/ltor(ipts,j))) 
             grow_wood = .TRUE.

             ! If the carbohydrate pool is too small, don't grow wood
             IF ( (pheno_type(j) .NE. 1) .AND. &
                  (biomass(ipts,j,icarbres,icarbon) .LE. reserve_pool) ) THEN

                grow_wood = .FALSE.

             ENDIF

             ! Write debug comments to output file
             IF ((j.EQ.test_pft .AND. ld_alloc) .OR. ld_warn) THEN
                CALL comment(npts, Cl_target, Cl, Cs_target, Cs, Cr_target, Cr, &
                     delta_ba, ipts, j, l, b_inc_tot, Cl_incp, Cs_incp, Cr_incp, &
                     KF, LF, Cl_inc, Cs_inc, Cr_inc, Cf_inc, grow_wood, &
                     circ_class_n, ind, 10)
             ENDIF


             !! 5.2.6 Ordinary growth
             !  Allometric relationship between components is respected, sustain 
             !  ordinary growth and allocate
             !  biomass to leaves, wood, roots and fruits.
             IF ( (SUM( ABS(Cl_target(:) - Cl(:)) ) .LE. min_stomate) .AND. &
                  (SUM( ABS(Cs_target(:) - Cs(:)) ) .LE. min_stomate) .AND. &
                  (SUM( ABS(Cr_target(:) - Cr(:)) ) .LE. min_stomate) .AND. &
                  (grow_wood) .AND. (b_inc_tot .GT. min_stomate) ) THEN

                ! Allocate fraction of carbon to fruit production (at the tree level)
                Cf_inc(:) = b_inc_tot / SUM(circ_class_n(ipts,j,:)) * fruit_alloc(j)

                ! Residual carbon is allocated to the other components (b_inc_tot is
                ! at the stand level)
                b_inc_tot = b_inc_tot * (un-fruit_alloc(j))

                ! Substitute (7), (8) and (9) in (1) 
                ! b_inc = tree_ff*pipe_density*(ba+circ_class_dba*gammas)*...
                ! (height+(circ_class_dba/s*gammas)) - Cs - Ch + ...
                !    KF*tree_ff*pipe_density*(ba+circ_class_dba*gammas) - ... 
                !    (KF*Ch)/(height+(circ_class_dba/s*gammas)) - Cl + ...
                !    KF/LF*tree_ff*pipe_density*(ba+circ_class_dba*gammas) - ...
                !    (KF*Ch/LF)/(height+(circ_class_dba/s*gammas)) - Cr
                !
                ! b_inc+Cs+Ch+Cl+Cr = tree_ff*pipe_density*(ba+circ_class_dba*gammas)*...
                !    (height+(circ_class_dba/s*gammas))  + ...
                !    KF*tree_ff*pipe_density*(ba+circ_class_dba*gammas) - ...
                !    (KF*Ch)/(height+(circ_class_dba/s*gammas)) + ...
                !    KF/LF*tree_ff*pipe_density*(ba+circ_class_dba*gammas) - ...
                !    (KF*Ch/LF)/(height+(circ_class_dba/s*gammas))
                ! <=> b_inc+Cs+Ch+Cl+Cr = circ_class_dba^2/s*tree_ff*...
                !    pipe_density*gammas^2 + circ_class_dba/s*ba*tree_ff*...
                !    pipe_density*gammas + ...
                !    circ_class_dba*height*tree_ff*pipe_density*gammas + ...
                !    bcirc_class_dba*height*tree_ff*pipe_density - ...
                !    (Ch*KF*s)/(circ_class_dba*gammas+height*s) + ...
                !    circ_class_dba*KF*tree_ff*pipe_density*gammas + ...
                !    ba*KF*tree_ff*pipe_density - ...
                !    (Ch*KF*s)/(LF*(circ_class_dba*gammas+height*s)) + ...
                !    circ_class_dba*KF/LF*tree_ff*pipe_density*gammas + ...
                !    ba*KF/LF*tree_ff*pipe_density
                ! (10) b_inc+Cs+Ch+Cl+Cr = (circ_class_dba^2/s*tree_ff*...
                !    pipe_density)*gammas^2 + ...
                !    (circ_class_dba/s*ba*tree_ff*pipe_density + ...
                !    circ_class_dba*height*tree_ff*pipe_density + ...
                !    circ_class_dba*KF*tree_ff*pipe_density + ...
                !    circ_class_dba*KF/LF*tree_ff*pipe_density)*gammas - ...
                !    (Ch*KF*s)(1+1/LF)/(circ_class_dba*gammas+height*s) + ...
                !    bcirc_class_dba*height*tree_ff*pipe_density + ...
                !    ba*KF*tree_ff*pipe_density + ba*KF/LF*tree_ff*pipe_density
                !
                ! Note that b_inc is not known, only b_inc_tot (= sum(b_inc) is known. 
                ! The above equations are for individual trees, at the stand level we 
                ! have to take the sum over the individuals which is 
                ! equivalant to substituting (10) in (2) 
                ! (11) sum(b_inc) + sum(Cs+Ch+Cl+Cr) = ...
                !    sum(circ_class_dba^2/s*tree_ff*pipe_density) * gammas^2 + ...
                !    sum(circ_class_dba/s*ba*tree_ff*pipe_density + ...
                !    circ_class_dba*height*tree_ff*pipe_density + ...
                !    circ_class_dba*KF*tree_ff*pipe_density + ...
                !    circ_class_dba*KF/LF*tree_ff*pipe_density) * gammas - ...
                !    sum[(Ch*KF*s)(1+1/LF)/(circ_class_dba*gammas+height*s)] + ...
                !    sum(bcirc_class_dba*height*tree_ff*pipe_density + ...
                !    ba*KF*tree_ff*pipe_density + ba*KF/LF*tree_ff*pipe_density)
                !
                ! The term sum[(Ch*KF*s)(1+1/LF)/(circ_class_dba*gammas+height*s)] 
                ! can be approximated by a series expansion
                ! (12) sum((Ch*KF*s)(1+1/LF)/(height*s) + ...
                !    sum((Ch*KF*s)(1+1/LF)*circ_class_dba/(height*s)^2)*gammas + ...
                !    sum((Ch*KF*s)(1+1/LF)*circ_class_dba^2/(height*s)^3)*gammas^2
                !
                ! Substitute (12) in (11)
                ! sum(b_inc) + sum(Cs+Ch+Cl+Cr) = ...
                !    sum(circ_class_dba^2/s*tree_ff*pipe_density - ...
                !    (Ch*KF*s)*(1+1/LF)*circ_class_dba^2/(height*s)^3) * gammas^2 + ...
                !    sum(circ_class_dba/s*ba*tree_ff*pipe_density + ...
                !    circ_class_dba*height*tree_ff*pipe_density + ...
                !    circ_class_dba*KF*tree_ff*pipe_density + ...
                !    circ_class_dba*KF/LF*tree_ff*pipe_density + ...
                !    (Ch*KF*s)*(1+1/LF)*circ_class_dba/(height*s)^2) * gammas + ...
                !    sum(bcirc_class_dba*height*tree_ff*pipe_density + ...
                !    ba*KF*tree_ff*pipe_density + ba*KF/LF*tree_ff*pipe_density - ...
                !    (Ch*KF*s)*(1+1/LF)/(height*s))
                !
                ! Solve this quadratic equation for gammas.
                a = SUM( circ_class_n(ipts,j,:) * &
                     (circ_class_dba(:)**2/s(:)*tree_ff(j)*pipe_density(j) - &
                     (Ch(:)*KF(ipts,j)*s(:))*(1+1/LF(ipts,j))*&
                     (circ_class_dba(:)**2/(circ_class_height_eff(:)*s(:))**3)) )
                b = SUM( circ_class_n(ipts,j,:) * &
                     (circ_class_dba(:)/s(:)*circ_class_ba_eff(:)*tree_ff(j)*pipe_density(j) + &
                     circ_class_dba(:)*circ_class_height_eff(:)*tree_ff(j)*pipe_density(j) + &
                     circ_class_dba(:)*KF(ipts,j)*tree_ff(j)*pipe_density(j) + &
                     circ_class_dba(:)*KF(ipts,j)/LF(ipts,j)*tree_ff(j)*pipe_density(j) + &
                     (Ch(:)*KF(ipts,j)*s(:))*(1+1/LF(ipts,j))*circ_class_dba(:)/&
                     (circ_class_height_eff(:)*s(:))**2) )
                c = SUM( circ_class_n(ipts,j,:) * &
                     (circ_class_ba_eff(:)*circ_class_height_eff(:)*&
                     tree_ff(j)*pipe_density(j) + &
                     circ_class_ba_eff(:)*KF(ipts,j)*tree_ff(j)*pipe_density(j) + &
                     circ_class_ba_eff(:)*KF(ipts,j)/LF(ipts,j)*tree_ff(j)*pipe_density(j) - &
                     (Ch(:)*KF(ipts,j)*s(:))*(1+1/LF(ipts,j))/&
                     (circ_class_height_eff(:)*s(:)) - &
                     (Cs(:) + Ch(:) + Cl(:) + Cr(:))) ) - b_inc_tot

                ! Solve the quadratic equation a*gammas2 + b*gammas + c = 0, for gammas.
                gammas(ipts,j) = (-b + sqrt(b**2-4*a*c)) / (2*a)  

                !++++++ TEMP ++++++
!!$                IF(test_pft == j .AND. test_grid ==ipts)THEN
!!$                   WRITE(numout,*) 'Testing for the slope'
!!$                   DO i=1,100
!!$                      tempi=1
!!$                      delta_ba(tempi) = circ_class_dba(tempi) * gammas * REAL(i)/50.0
!!$                      delta_height(tempi) = delta_ba(tempi)/s(tempi)
!!$                      Cs_inc(tempi) = tree_ff(j)*pipe_density(j)*(circ_class_ba_eff(tempi) + delta_ba(tempi))*(circ_class_height_eff(tempi) + &
!!$                           delta_height(tempi)) - Cs(tempi) - Ch(tempi)
!!$                      Cl_inc(tempi) = KF(ipts,j)*tree_ff(j)*pipe_density(j)*(circ_class_ba_eff(tempi)+delta_ba(tempi)) - &
!!$                           (KF(ipts,j)*Ch(tempi))/(circ_class_height_eff(tempi)+delta_height(tempi)) - Cl(tempi)
!!$                      Cr_inc(tempi) = KF(ipts,j)/LF(ipts,j)*tree_ff(j)*pipe_density(j)*(circ_class_ba_eff(tempi)+delta_ba(tempi)) - &
!!$                           (KF(ipts,j)*Ch(tempi)/LF(ipts,j))/(circ_class_height_eff(tempi)+delta_height(tempi)) - Cr(tempi)    
!!$                      WRITE(numout,'(10F20.10)') delta_height(tempi),delta_ba(tempi),Cs_inc(tempi),Cl_inc(tempi),Cr_inc(tempi)
!!$                   ENDDO
!!$                   WRITE(numout,*) 'End testing for the slope'
!!$                END IF
                !+++++++++++++++++
                
                ! The solution for gammas is then used to calculate delta_ba (eq. 3), 
                ! delta_height (eq. 6), Cs_inc (eq. 7), Cl_inc (eq. 8) and Cr_inc (eq. 9). 
                ! See comment on the calculation of delta_height and its implications on 
                ! numerical consistency at the similar statement in §5.2.4.3.1
                delta_ba(:) = circ_class_dba(:) * gammas(ipts,j)
                store_delta_ba(ipts,j,:) = delta_ba(:)
                delta_height(:) = delta_ba(:)/s(:)              
                Cs_inc(:) = tree_ff(j)*pipe_density(j)*(circ_class_ba_eff(:) + &
                     delta_ba(:))*(circ_class_height_eff(:) + &
                     delta_height(:)) - Cs(:) - Ch(:)
                Cl_inc(:) = KF(ipts,j)*tree_ff(j)*pipe_density(j)*&
                     (circ_class_ba_eff(:)+delta_ba(:)) - &
                     (KF(ipts,j)*Ch(:))/(circ_class_height_eff(:)+delta_height(:)) - Cl(:)
                Cr_inc = KF(ipts,j)/LF(ipts,j)*tree_ff(j)*pipe_density(j)*&
                     (circ_class_ba_eff(:)+delta_ba(:)) - &
                     (KF(ipts,j)*Ch(:)/LF(ipts,j))/(circ_class_height_eff(:)+&
                     delta_height(:)) - Cr(:)

                ! After thousands of simulation years we had a single pixel where 
                ! one of the three circ_class got a negative growth. The cause is not
                ! is not entirely clear but could be related to the fact that KF
                ! changes from one day to another in combination with a low b_inc.
                ! If this happens, we don't allocate and simply leave the carbon
                ! in the labile pool. We will try again with more carbon the next day.
                ! This case is dealt with later in the code - see warning 10
                

                ! Write debug comments to output file
                IF ((j.EQ.test_pft .AND. ld_alloc) .OR. ld_warn) THEN
                   WRITE(numout,*) 'a, b, c, gammas, ', a, b, c, gammas(ipts,j)
                   WRITE(numout,*) 'delta_height, ', delta_height(:)
                   CALL comment(npts, Cl_target, Cl, Cs_target, Cs, Cr_target, Cr, &
                        delta_ba, ipts, j, l, b_inc_tot, Cl_incp, Cs_incp, Cr_incp, &
                        KF, LF, Cl_inc, Cs_inc, Cr_inc, Cf_inc, grow_wood, &
                        circ_class_n, ind, 11)
                ENDIF

                ! Wrap-up ordinary growth  
                ! Calculate C that was not allocated, note that Cf_inc was already substracted
                b_inc_tot = b_inc_tot - &
                     SUM(circ_class_n(ipts,j,:)*(Cl_inc(:) + Cr_inc(:) + Cs_inc(:)))

                !---TEMP---
                IF (j.EQ.test_pft .AND. ld_alloc.AND. ipts==test_grid) THEN
                   WRITE(numout,*) 'wrap-up ordinary allocation, left b_in_tot, ', b_inc_tot 
                ENDIF
                !----------


             !! 5.2.7 Don't grow wood, use C to fill labile pool
             ELSEIF ( (.NOT. grow_wood) .AND. (b_inc_tot .GT. min_stomate) ) THEN

                ! Calculate the C that needs to be distributed to the 
                ! labile pool. The fraction is proportional to the ratio 
                ! between the total allocatable biomass and the unallocated 
                ! biomass per tree (b_inc now contains the unallocated 
                ! biomass). At the end of the allocation scheme bm_alloc_tot 
                ! is substracted from the labile biomass pool to update the 
                ! biomass pool (biomass(:,:,ilabile) = biomass(:,:,ilabile) - 
                ! bm_alloc_tot(:,:)). At that point, the scheme puts the 
                ! unallocated b_inc into the labile pool. What we 
                ! want is that the unallocated fraction is removed from 
                ! ::bm_alloc_tot such that only the allocated C is removed 
                ! from the labile pool. DSGdebug_01: 

                bm_alloc_tot(ipts,j) = bm_alloc_tot(ipts,j) - b_inc_tot

                ! Wrap-up ordinary growth  
                ! Calculate C that was not allocated (b_inc_tot), the 
                ! equation should read b_inc_tot = b_inc_tot - b_inc_tot 
                ! note that Cf_inc was already substracted
                b_inc_tot = zero 

                !---TEMP---
                IF (j.EQ.test_pft .AND. ld_alloc.AND. ipts==test_grid) THEN
                   WRITE(numout,*) 'No wood growth, move remaining C to labile pool'
                   WRITE(numout,*) 'bm_alloc_tot_new, ',bm_alloc_tot(ipts,j)
                   WRITE(numout,*) 'wrap-up ordinary allocation, left b_inc_tot, ', b_inc_tot 
                ENDIF
                !---------- 


             !! 5.2.8 Error - the allocation scheme is overspending 
             ELSEIF (b_inc_tot .LT. min_stomate) THEN

                IF (b_inc_tot .LT. -min_stomate) THEN

                   ! Something is wrong with the calculations
                   WRITE(numout,*) 'WARNING 7: numerical problem overspending in ordinary allocation'
                   WRITE(numout,*) 'WARNING 7: PFT, ipts: ',j,ipts
                   WRITE(numout,*) 'WARNING 7: b_inc_tot', b_inc_tot
                   IF(ld_stop)THEN
                      CALL ipslerr_p (3,'growth_fun_all',&
                           'WARNING 7: numerical problem overspending in ordinary allocation','','')
                   ENDIF

                ELSE

                   IF (j.EQ.test_pft .AND. ld_alloc.AND. ipts==test_grid) THEN

                      ! Succesful allocation
                      WRITE(numout,*) 'Successful allocation'

                   ENDIF

                ENDIF

                ! Althought the biomass components respect the allometric relationships, there
                ! is no carbon left to allocate                      
                b_inc_tot = zero
                Cl_inc(:) = zero
                Cs_inc(:) = zero
                Cr_inc(:) = zero
                Cf_inc(:) = zero

             ENDIF ! Ordinary allocation

             !! 5.2.9 Forced allocation
             !  Although this should not happen, in case the functional allocation did not 
             !  consume all the allocatable carbon, the remaining C is left for the next day, 
             !  and some of the biomass is used to produce fruits (tuned). The numerical 
             !  precision of the allocation scheme (i.e. the linearisation) is similar to 
             !  min_stomate (i.e. 10-8) resulting in 'false' warnings. In the latter case 
             !  forced allocation is applied but only for very small amounts of carbon 
             !  i.e. between 10-5 and 10-8.  
             IF ( b_inc_tot .GT. min_stomate) THEN

                WRITE(numout,*) 'WARNING 8: b_inc_tot greater than min_stomate force allocation'
                WRITE(numout,*) 'WARNING 8: PFT, ipts: ',j,ipts
                WRITE(numout,*) 'WARNING 8: b_inc_tot, ', b_inc_tot
                IF(ld_stop)THEN
                   CALL ipslerr_p (3,'growth_fun_all',&
                        'WARNING 8: b_inc_tot greater than min_stomate force allocation','','')
                ENDIF

                !+++CHECK+++
                ! We should not end-up here. We need some code to break the conditions 
                ! that made us end-up here. The current code will do this job.
!!$                ! Calculate fraction that will be allocated to fruit. The fraction is proportional to the 
!!$                ! ratio between the total allocatable biomass and the unallocated biomass per tree
!!$                frac = fruit_alloc(j) * MIN(1., bm_alloc_tot(ipts,j) / b_inc_tot)  
!!$                Cf_inc(:) = Cf_inc(:) + b_inc_tot * frac
!!$                b_inc_tot = b_inc_tot * (1 - frac)
!!$
!!$                ! Calculate the C that needs to be distributed to the labile pool. The fraction is proportional 
!!$                ! to the ratio between the total allocatable biomass and the unallocated biomass per tree (b_inc 
!!$                ! now contains the unallocated biomass). At the end of the allocation scheme bm_alloc_tot is 
!!$                ! substracted from the labile biomass pool to update the biomass pool (biomass(:,:,ilabile) = 
!!$                ! biomass(:,:,ilabile,icarbon) - bm_alloc_tot(:,:)). At that point, the scheme puts the 
!!$                ! unallocated b_inc into the labile pool.
!!$                bm_alloc_tot(ipts,j) = bm_alloc_tot(ipts,j) * &
!!$                     ( 1. - (1.-frac) * b_inc_tot / bm_alloc_tot(ipts,j) )
                !+++++++++++

             ELSEIF ( (b_inc_tot .LT. min_stomate) .AND. (b_inc_tot .GE. -min_stomate) ) THEN

                ! Successful allocation
                IF (j.EQ.test_pft .AND. ld_alloc.AND. ipts==test_grid) THEN
                   WRITE(numout,*) 'Successful allocation'
                ENDIF

             ELSE

                ! Something possibly important was overlooked
                IF ( (b_inc_tot .LT. 100*min_stomate) .AND. (b_inc_tot .GE. -100*min_stomate) ) THEN
                   IF (j.EQ.test_pft .AND. ld_alloc.AND. ipts==test_grid) THEN
                      WRITE(numout,*) 'Marginally successful allocation - precision better than 10-6'
                      WRITE(numout,*) 'PFT, b_inc_tot', j, b_inc_tot
                   ENDIF
                ELSE
                   WRITE(numout,*) 'WARNING 9: Logical flaw unexpected result in the ordinary allocation'
                   WRITE(numout,*) 'WARNING 9: b_inc_tot, ',b_inc_tot
                   WRITE(numout,*) 'WARNING 9: PFT, ipts: ',j,ipts
                  IF(ld_stop)THEN
                     CALL ipslerr_p (3,'growth_fun_all',&
                          'WARNING 9: Logical flaw unexpected result in the ordinary allocation','','')
                  ENDIF
                ENDIF

             ENDIF

             ! The second problem we need to catch is when one of the increment pools is 
             ! negative. This is an undesired outcome (see comment where ::KF_old is 
             ! calculated in this routine. In that case we write a warning, set all increment
             ! pools to zero and try it again at the next time step. A likely cause of this 
             ! problem is a too large change in KF from one time step to another. Try decreasing
             ! the acceptable value for an absolute increase in KF.
             IF (MINVAL(Cs_inc(:)) .LT. zero .OR. MINVAL(Cr_inc(:)) .LT. zero .OR. &
                  MINVAL(Cs_inc(:)) .LT. zero) THEN

                ! Do not allocate - save the carbon for the next time step
                WRITE(numout,*) 'WARNING 10: numerical problem, one of the increment pools is less than zero'
                WRITE(numout,*) 'WARNING 10: PFT, ipts: ',j,ipts
                WRITE(numout,*) 'WARNING 10: Cl_inc(:): ',Cl_inc(:)
                WRITE(numout,*) 'WARNING 10: Cr_inc(:): ',Cr_inc(:)
                WRITE(numout,*) 'WARNING 10: Cs_inc(:): ',Cs_inc(:)
                WRITE(numout,*) 'WARNING 10: PFT, ipts: ',j,ipts
                WRITE(numout,*) 'WARNING 10: We will undo the allocation'
                WRITE(numout,*) ' and save the carbon for the next day'

                Cl_inc(:) = zero
                Cr_inc(:) = zero
                Cs_inc(:) = zero

             ENDIF

             !! 5.2.10 Wrap-up phenological and ordinary allocation
             Cl_inc(:) = Cl_inc(:) + Cl_incp(:)
             Cr_inc(:) = Cr_inc(:) + Cr_incp(:)
             Cs_inc(:) = Cs_inc(:) + Cs_incp(:)
             residual(ipts,j) = b_inc_tot

             !---TEST---
             IF (j.EQ.test_pft .AND. ld_alloc.AND. ipts==test_grid) THEN
                WRITE(numout,*) 'Final allocation', ipts, j
                WRITE(numout,*) 'Cl, Cs, Cr', Cl(:), Cs(:), Cr(:) 
                WRITE(numout,*) 'Cl_incp, Cs_incp, Cr_incp, ', Cl_incp(:), Cs_incp(:), Cr_incp(:)
                WRITE(numout,*) 'Cl_inc, Cs_ins, Cr_inc, Cf_inc, ', Cl_inc(:), Cs_inc(:), Cr_inc(:), Cf_inc(:)
                WRITE(numout,*) 'unallocated/residual, ', b_inc_tot
                WRITE(numout,*) 'Old ba, delta_ba, new ba, ', circ_class_ba_eff(:), delta_ba(:), circ_class_ba_eff(:)+delta_ba(:)
                DO l=1,ncirc
                   WRITE(numout,*) 'Circ_class_biomass, ',circ_class_biomass(ipts,j,l,:,icarbon)
                ENDDO
             ENDIF
             !----------


!============================================================================================================
!! DSG: This is residual is totally avoidable if one would not have several variables
!for the same thing :(
! DSG: the following bit of code appears twice (grep WARNING 23)
!

             !! 5.2.11 Account for the residual
             !  The residual is usually around ::min_stomate but we deal 
             !  with it anyway to make sure the mass balance is closed
             !  and as a way to detect errors. Move the unallocated carbon
             !  back into the labile pool


             deficit = zero 

             IF (ABS((biomass(ipts,j,ilabile,icarbon) -   bm_alloc_tot(ipts,j))        &   
                     + residual(ipts,j)) .GT. min_stomate) THEN

                deficit = (biomass(ipts,j,ilabile,icarbon) - bm_alloc_tot(ipts,j) ) &
                            + residual(ipts,j)

                ! The deficit is less than the carbon reserve
                IF (-deficit .LE. biomass(ipts,j,icarbres,icarbon)) THEN

                   ! Pay the deficit from the reserve pool
                   biomass(ipts,j,icarbres,icarbon) = &
                        biomass(ipts,j,icarbres,icarbon) + deficit
                   biomass(ipts,j,ilabile,icarbon)  = &
                        biomass(ipts,j,ilabile,icarbon) - deficit

                ! If not, try to reduce it from growth respiration:
                ELSEIF ( ABS(deficit) .LE. (resp_growth(ipts,j) + resp_maint(ipts,j))) THEN

                   biomass(ipts,j,ilabile,icarbon)  = &
                        biomass(ipts,j,ilabile,icarbon) - deficit
! arbitrary removal from respiration components
                   ! remove it from maint_respiration
                   resp_maint(ipts,j) = resp_maint(ipts,j) + deficit
                   ! if deficit is more than maint resp remove the rest from
                   ! growth resp
                   IF (resp_maint(ipts,j) .LT. zero) THEN
                       resp_growth(ipts,j) = resp_growth(ipts,j) + resp_maint(ipts,j)
                       resp_maint(ipts,j)  = zero
                   ENDIF
! arbitrary removal from respiration components
                ELSE

                   ! Not enough carbon to pay the deficit
                   ! There is likely a bigger problem somewhere in
                   ! this routine
!DSGtempFIX
                   !  DSG: There isn't a bigger problem in the routine, but the
                   !  routine itself is the problem. 
                   !  This error can be trigger when tuning the model on-the-fly.
                   !  Thus, we won't stop the model but prevent any growth this time
                   !  step. Allocation routine can try next timestep again to
                   !  get things right. 
                   bm_alloc_tot(ipts,j) = zero
!DSGtempFIX

!!DSG: nice, why didn't you fix it then?!?!
!                   WRITE(numout,*) 'WARNING 11: PFT, ipts: ',j,ipts
!                   CALL ipslerr_p (3,'growth_fun_all',&
!                        'WARNING 11: numerical problem overspending ',&
!                        'when trying to account for unallocatable C ','')
!!DSG: nice, why didn't you fix it then?!?!

                ENDIF

             ENDIF

!============================================================================================================

             !! 5.2.12 Standardise allocation factors
             !  Strictly speaking the allocation factors do not need to be 
             !  calculated because the functional allocation scheme allocates 
             !  absolute amounts of carbon. Hence, Cl_inc could simply be 
             !  added to biomass(:,:,ileaf,icarbon), Cr_inc to 
             !  biomass(:,:,iroot,icarbon), etc. However, using allocation 
             !  factors bears some elegance in respect to distributing the 
             !  growth respiration if this would be required. Further it 
             !  facilitates comparison to the resource limited allocation 
             !  scheme (stomate_growth_res_lim.f90) and it comes in handy 
             !  for model-data comparison. This allocation takes place at 
             !  the tree level - note that ::biomass is the only prognostic 
             !  variable from the tree-based allocation

             !  Allocation   
             Cl_inc(:) = MAX(zero, circ_class_n(ipts,j,:) * Cl_inc(:))
             Cr_inc(:) = MAX(zero, circ_class_n(ipts,j,:) * Cr_inc(:))
             Cs_inc(:) = MAX(zero, circ_class_n(ipts,j,:) * Cs_inc(:))
             Cf_inc(:) = MAX(zero, circ_class_n(ipts,j,:) * Cf_inc(:))

             ! Total_inc is based on the updated Cl_inc, Cr_inc, Cs_inc and Cf_inc. Therefore, do not multiply
             ! circ_class_n(ipts,j,:) again
             total_inc = SUM(Cf_inc(:) + Cl_inc(:) + Cs_inc(:) + Cr_inc(:))

             ! Relative allocation
             IF ( total_inc .GT. min_stomate ) THEN

                Cl_inc(:) = Cl_inc(:) / total_inc
                Cs_inc(:) = Cs_inc(:) / total_inc
                Cr_inc(:) = Cr_inc(:) / total_inc
                Cf_inc(:) = Cf_inc(:) / total_inc

             ELSE

                bm_alloc_tot(ipts,j) = zero
                Cl_inc(:) = zero
                Cs_inc(:) = zero
                Cr_inc(:) = zero
                Cf_inc(:) = zero

             ENDIF


             !! 5.2.13 Convert allocation to allocation facors
             !  Convert allocation of individuals to ORCHIDEE's allocation 
             !  factors - see comment for 5.2.5. Aboveground sapwood 
             !  allocation is age dependent in trees. ::alloc_min and 
             !  ::alloc_max must range between 0 and 1. 
             alloc_sap_above = alloc_min(j) + ( alloc_max(j) - alloc_min(j) ) * &
                  ( 1. - EXP( -age(ipts,j) / demi_alloc(j) ) )

             ! Leaf, wood, root and fruit allocation
             f_alloc(ipts,j,ileaf)     = SUM(Cl_inc(:))
             f_alloc(ipts,j,isapabove) = SUM(Cs_inc(:)*alloc_sap_above)
             f_alloc(ipts,j,isapbelow) = SUM(Cs_inc(:)*(1.-alloc_sap_above))
             f_alloc(ipts,j,iroot)     = SUM(Cr_inc(:))
             f_alloc(ipts,j,ifruit)    = SUM(Cf_inc(:))
             
             ! Absolute allocation at the tree level and for an individual tree (gC tree-1)
             ! The labile and reserve pools are not allocated at the tree level. However,
             ! stand level ilabile and icarbres biomass will be redistributed at the tree
             ! level later in this subroutine. This is done after the relative allocation
             ! beacuse now ::alloc_sap_above is known
             circ_class_biomass(ipts,j,:,ileaf,icarbon) = &
                  circ_class_biomass(ipts,j,:,ileaf,icarbon) + &
                  ( Cl_inc(:) * total_inc / circ_class_n(ipts,j,:) )
             circ_class_biomass(ipts,j,:,isapabove,icarbon) = &
                  circ_class_biomass(ipts,j,:,isapabove,icarbon) + &
                  ( Cs_inc(:) * alloc_sap_above * total_inc / &
                  circ_class_n(ipts,j,:) )
             circ_class_biomass(ipts,j,:,isapbelow,icarbon) = &
                  circ_class_biomass(ipts,j,:,isapbelow,icarbon) + &
                  ( Cs_inc(:) * (un - alloc_sap_above) * &
                  total_inc / circ_class_n(ipts,j,:) )
             circ_class_biomass(ipts,j,:,iroot,icarbon) = &
                  circ_class_biomass(ipts,j,:,iroot,icarbon) + &
                  ( Cr_inc(:) * total_inc / circ_class_n(ipts,j,:) )
             circ_class_biomass(ipts,j,:,ifruit,icarbon) = &
                  circ_class_biomass(ipts,j,:,ifruit,icarbon) + &
                  ( Cf_inc(:) * total_inc / circ_class_n(ipts,j,:) )

             circ_class_biomass(ipts,j,:,ileaf,initrogen)     = circ_class_biomass(ipts,j,:,ileaf,initrogen) + &
                  ( Cl_inc(:) * total_inc / circ_class_n(ipts,j,:) )/cn_leaf(ipts,j)
             circ_class_biomass(ipts,j,:,isapabove,initrogen) = circ_class_biomass(ipts,j,:,isapabove,initrogen) + &
                  ( Cs_inc(:) * alloc_sap_above * total_inc / circ_class_n(ipts,j,:) )/cn_leaf(ipts,j)*fcn_wood_act(ipts,j)
             circ_class_biomass(ipts,j,:,isapbelow,initrogen) = circ_class_biomass(ipts,j,:,isapbelow,initrogen) + &
                  ( Cs_inc(:) * (un - alloc_sap_above) * total_inc / circ_class_n(ipts,j,:) )/cn_leaf(ipts,j)*fcn_wood_act(ipts,j)
             circ_class_biomass(ipts,j,:,iroot,initrogen)     = circ_class_biomass(ipts,j,:,iroot,initrogen) + &
                  ( Cr_inc(:) * total_inc / circ_class_n(ipts,j,:) )/cn_leaf(ipts,j)*fcn_root(j)
             circ_class_biomass(ipts,j,:,ifruit,initrogen)    = circ_class_biomass(ipts,j,:,ifruit,initrogen) + &
                  ( Cf_inc(:) * total_inc / circ_class_n(ipts,j,:) )/cn_leaf(ipts,j)*fcn_root(j)


             circ_class_biomass(ipts,j,:,ileaf,iphosphorus)     = circ_class_biomass(ipts,j,:,ileaf,iphosphorus)  + &
                  ( Cl_inc(:) * total_inc / circ_class_n(ipts,j,:) )/(cn_leaf(ipts,j)*np_leaf(ipts,j))
             circ_class_biomass(ipts,j,:,isapabove,iphosphorus) = circ_class_biomass(ipts,j,:,isapabove,iphosphorus) + &
                  ( Cs_inc(:) * alloc_sap_above * total_inc / circ_class_n(ipts,j,:) )/                                &
                  (cn_leaf(ipts,j)*np_leaf(ipts,j))*(fcn_wood_act(ipts,j)*fnp_wood_act(ipts,j))
             circ_class_biomass(ipts,j,:,isapbelow,iphosphorus) = circ_class_biomass(ipts,j,:,isapbelow,iphosphorus) + &
                  ( Cs_inc(:) * (un - alloc_sap_above) * total_inc / circ_class_n(ipts,j,:) )/                         &
                  (cn_leaf(ipts,j)*np_leaf(ipts,j))*(fcn_wood_act(ipts,j)*fnp_wood_act(ipts,j))
             circ_class_biomass(ipts,j,:,iroot,iphosphorus)     = circ_class_biomass(ipts,j,:,iroot,iphosphorus) + &
                  ( Cr_inc(:) * total_inc / circ_class_n(ipts,j,:) )/                       &  
                  (cn_leaf(ipts,j)*np_leaf(ipts,j))*(fcn_root(j)*fnp_root(j))
             circ_class_biomass(ipts,j,:,ifruit,iphosphorus)    = circ_class_biomass(ipts,j,:,ifruit,iphosphorus) + &
                  ( Cf_inc(:) * total_inc / circ_class_n(ipts,j,:) )/                       &
                  (cn_leaf(ipts,j)*np_leaf(ipts,j))*(fcn_root(j)*fnp_root(j))


!             !+++TEMP+++
!             IF (ld_alloc) THEN
!                tempi = zero
!                DO icirc = 1,ncirc
!                   IF ( Cl_inc(icirc) .LT. zero) THEN
!                      WRITE(numout,*) 'Cl_inc, ', j, Cl_inc(icirc)
!                      tempi = un
!                   ENDIF
!                   IF ( Cs_inc(icirc) * alloc_sap_above .LT. zero) THEN
!                      WRITE(numout,*) 'Cs_inc aboveground, ', j, &
!                           Cs_inc(icirc) * alloc_sap_above
!                      tempi = un
!                   ENDIF
!                   IF ( Cs_inc(icirc) * (un - alloc_sap_above) .LT. zero) THEN
!                      WRITE(numout,*) 'Cs_inc aboveground, ', j, &
!                           Cs_inc(icirc) * (un-alloc_sap_above)
!                      tempi = un
!                   ENDIF
!                   IF ( Cr_inc(icirc) .LT. zero) THEN
!                      WRITE(numout,*) 'Cr_inc, ', j, Cr_inc(icirc)
!                      tempi = un
!                   ENDIF
!                   IF ( Cf_inc(icirc) .LT. zero) THEN
!                      WRITE(numout,*) 'Cf_inc, ', j, Cf_inc(icirc)
!                      tempi = un
!                   ENDIF
!                   IF ( total_inc .LT. zero) THEN
!                      WRITE(numout,*) 'total_inc, ', j, total_inc
!                      tempi = un
!                   ENDIF
!                   IF ( circ_class_n(ipts,j,icirc) .LT. zero) THEN
!                      WRITE(numout,*) 'circ_class_n, ', j, circ_class_n(ipts,j,icirc)
!                      tempi = un
!                   ENDIF
!                ENDDO
!                IF (tempi == un) CALL ipslerr_p (3,'growth_fun_all',&
!                     'WARNING 11bis: the solution has negative values',&
!                     'None of these variables should be negative','')   
!             ENDIF
!             !++++++++++

          ELSEIF (is_tree(j)) THEN
       
             IF(ld_alloc) WRITE(numout,*) 'there is no tree biomass to allocate, PFT, ', j 

          ENDIF ! Is there biomass to allocate (§5.2 - far far up)

!!!!! JC GRASS START
          !! 5.3 Calculate allocated biomass pools for grasses and crops
          !  Only possible if there is biomass to allocate
          IF ( .NOT. is_tree(j) .AND. bm_alloc_tot(ipts,j) .GT. min_stomate ) THEN

             !! 5.3.1 Scaling factor to convert variables to the individual plant
             !  Allocation is on an individual basis (gC ind-1). Stand-level variables 
             !  need to convert to a single individual. The absence of sapwood makes 
             !  this irrelevant because the allocation reduces to a linear function 
             !  (contrary to the non-linearity of tree allocation). For the
             !  beauty of consistency, the transformations will be implemented.
             !  Different approach between the DGVM and statitic approach
             IF (ok_dgvm) THEN

                ! The DGVM does NOT work with the functional allocation. Consider
                ! this code as a placeholder. The original code had two different 
                ! transformations to calculate the scalars. Both could be used but 
                ! the units will differ. For consistency only one was retained 
                ! scal = ind(ipts,j) * cn_ind(ipts,j) / veget_max(ipts,j)
                scal(ipts,j) = veget_max(ipts,j) / ind(ipts,j)

             ELSE

                ! By dividing the actual biomass by the number of individuals
                ! the biomass of an individual is obtained. Note that a grass/crop 
                ! individual was defined as 1m-2 of vegetation 
                scal(ipts,j) = 1./ ind(ipts,j)

             ENDIF


             !! 5.3.2 Current biomass pools per grass/crop (gC ind^-1)
             !  Cs has too many dimensions for grass/crops. To have a consistent notation the same variables
             !  are used as for trees but the dimension of Cs, Cl and Cr i.e. ::ncirc should be ignored            
             Cs(:) = biomass(ipts,j,isapabove,icarbon) * scal(ipts,j)
             Cr(:) = biomass(ipts,j,iroot,icarbon) * scal(ipts,j)
             Cl(:) = biomass(ipts,j,ileaf,icarbon) * scal(ipts,j)
             Ch(:) = zero

             ! Total amount of carbon that needs to be allocated (::bm_alloc_tot). bm_alloc_tot is
             ! in gC m-2 day-1. At 1 m2 there are ::ind number of grasses.
             b_inc_tot = bm_alloc_tot(ipts,j)

             !! 5.3.3 C-allocation for crops and grasses
             !  The mass conservation equations are detailed in the header of this subroutine.
             !  The scheme assumes a functional relationships between leaves and roots for grasses and crops. 
             !  When carbon is added to the leaf biomass pool, an increase in the root biomass is to be 
             !  expected to sustain water transport from the roots to the leaves.

             !! 5.3.3.1 Do the biomass pools respect the pipe model?
             !  Do the current leaf, sapwood and root components respect the allometric 
             !  constraints? Calculate the optimal root and leaf mass, given the current wood mass 
             !  by using the basic allometric relationships. Calculate the optimal sapwood
             !  mass as a function of the current leaf and root mass.
             Cl_target(1) = MAX( Cs(1) * KF(ipts,j) , Cr(1) * LF(ipts,j), Cl(1) )
             Cs_target(1) = MAX( Cl_target(1) / KF(ipts,j), Cr(1) * LF(ipts,j) / KF(ipts,j), Cs(1) ) 
             Cr_target(1) = MAX( Cl_target(1) / LF(ipts,j), Cs_target(1) * KF(ipts,j) / LF(ipts,j), Cr(1) )
             
             ! Write debug comments to output file
             IF (j .EQ. test_pft .AND. ld_alloc.AND. ipts==test_grid) THEN
                WRITE(numout,*) 'bm_alloc_tot, ',bm_alloc_tot(ipts,j)
                WRITE(numout,*) 'Does the grass/crop needs reshaping?'
                WRITE(numout,*) 'KF, LF, ', KF(ipts,j), LF(ipts,j)
                WRITE(numout,*) 'qm_height, ', (Cl_target(1) * sla(j) * lai_to_height(j))
                WRITE(numout,*) 'Cl_target-Cl, ', Cl_target(1)-Cl(1), Cl_target(1), Cl(1)
                WRITE(numout,*) 'Cs_target-Cs, ', Cs_target(1)-Cs(1), Cs_target(1), Cs(1)
                WRITE(numout,*) 'Cr_target-Cr, ', Cr_target(1)-Cr(1), Cr_target(1), Cr(1)
             ENDIF


             !! 5.3.3.2 Phenological growth
             !  Phenological growth and reshaping of the grass/crop in line with the pipe model. Turnover removes 
             !  C from the different plant components but at a component-specific rate, as such the allometric 
             !  constraints are distorted at every time step and should be restored before ordinary growth can 
             !  take place

             !! 5.3.3.2.1 The available structural C can sustain the available leaves and roots
             !  Calculate whether the structural c is in allometric balance. The target values should 
             !  always be larger than the current pools so the use of ABS is redundant but was used to 
             !  be on the safe side (here and in the rest of the module) as it could help to find 
             !  logical flaws.        
             IF ( ABS(Cs_target(1) - Cs(1)) .LT. min_stomate ) THEN

                Cs_incp(1) = MAX(zero, Cs_target(1) - Cs(1))

                ! Enough leaves and structural biomass, only grow roots
                IF ( ABS(Cl_target(1) - Cl(1))  .LT. min_stomate ) THEN

                   ! Allocate at the tree level to restore allometric balance
                   Cl_incp(1) = MAX(zero, Cl_target(1) - Cl(1))
                   Cr_incp(1) = MAX( MIN(b_inc_tot / ind(ipts,j) - Cs_incp(1) - Cl_incp(1), &
                        Cr_target(1) - Cr(1)), zero )

                   ! Write debug comments to output file
                   IF ((j.EQ.test_pft .AND. ld_alloc) .OR. ld_warn) THEN
                      CALL comment(npts, Cl_target, Cl, Cs_target, Cs, Cr_target, Cr, delta_ba, ipts, j, l, &
                           b_inc_tot, Cl_incp, Cs_incp, Cr_incp, KF, LF, Cl_inc, Cs_inc, Cr_inc, Cf_inc, &
                           grow_wood, circ_class_n, ind, 12)
                   ENDIF

                ! Sufficient structural C and roots, allocate C to leaves
                ELSEIF ( ABS(Cr_target(1) - Cr(1)) .LT. min_stomate ) THEN

                   ! Allocate at the tree level to restore allometric balance 
                   Cr_incp(1) = MAX(zero, Cr_target(1) - Cr(1))
                   Cl_incp(1) = MAX( MIN(b_inc_tot / ind(ipts,j) - Cs_incp(1) - Cr_incp(1), &
                        Cl_target(1) - Cl(1)), zero )

                   ! Update vegetation height
                   qm_height = (Cl(1) + Cl_incp(1)) * sla_calc(ipts,j) * lai_to_height(j)
               
                   ! Write debug comments to output file
                   IF ((j.EQ.test_pft .AND. ld_alloc) .OR. ld_warn) THEN
                      CALL comment(npts, Cl_target, Cl, Cs_target, Cs, Cr_target, Cr, delta_ba, ipts, j, l, &
                           b_inc_tot, Cl_incp, Cs_incp, Cr_incp, KF, LF, Cl_inc, Cs_inc, Cr_inc, Cf_inc, &
                           grow_wood, circ_class_n, ind, 13)
                   ENDIF

                ! Both leaves and roots are needed to restore the allometric relationships
                ELSEIF ( ABS(Cl_target(1) - Cl(1)) .GT. min_stomate .AND. &
                     ABS(Cr_target(1) - Cr(1)) .GT. min_stomate ) THEN                 

                   ! Allocate at the tree level to restore allometric balance
                   !  The equations can be rearanged and written as
                   !  (i) b_inc = Cl_inc + Cr_inc
                   !  (ii) Cr_inc = (Cl_inc+Cl)/LF - Cr
                   !  Substitue (ii) in (i) and solve for Cl_inc
                   !  <=> Cl_inc = (LF*(b_inc+Cr)-Cl)/(1+LF)
                   Cl_incp(1) = MIN( ((LF(ipts,j) * ((b_inc_tot/ind(ipts,j)) - Cs_incp(1) + Cr(1))) - Cl(1)) / & 
                        (1 + LF(ipts,j)), Cl_target(1) - Cl(1) )
                   Cr_incp(1) = MIN ( ((Cl_incp(1) + Cl(1)) / LF(ipts,j)) - Cr(1), &
                        Cr_target(1) - Cr(1))

                   ! The imbalance between Cr and Cl can be so big that (Cl+Cl_inc)/LF is still less
                   ! then the available root carbon (observed!). This would result in a negative Cr_incp
                   IF ( Cr_incp(1) .LT. zero ) THEN

                      Cl_incp(1) = MIN( b_inc_tot/ind(ipts,j) - Cs_incp(1), Cl_target(1) - Cl(1) )
                      Cr_incp(1) = b_inc_tot/ind(ipts,j) - Cs_incp(1) - Cl_incp(1)

                   ELSEIF (Cl_incp(1) .LT. zero) THEN

                      Cr_incp(1) = MIN( b_inc_tot/ind(ipts,j) - Cs_incp(1), Cr_target(1) - Cr(1) )
                      Cl_incp(1) = (b_inc_tot/ind(ipts,j)) - Cs_incp(1) - Cr_incp(1)

                   ENDIF

                   ! Update vegetation height
                   qm_height = (Cl(1) + Cl_incp(1)) * sla_calc(ipts,j) * lai_to_height(j)

                   ! Write debug comments to output file
                   IF ((j.EQ.test_pft .AND. ld_alloc) .OR. ld_warn) THEN
                      CALL comment(npts, Cl_target, Cl, Cs_target, Cs, Cr_target, Cr, delta_ba, ipts, j, l, &
                           b_inc_tot, Cl_incp, Cs_incp, Cr_incp, KF, LF, Cl_inc, Cs_inc, Cr_inc, Cf_inc, &
                           grow_wood, circ_class_n, ind, 14)
                   ENDIF    

                ELSE

                   WRITE(numout,*) 'WARNING 12: Exc 1-3 unexpected exception'
                   WRITE(numout,*) 'WARNING 12: PFT, ipts: ',j,ipts
                   IF(ld_stop)THEN
                      CALL ipslerr_p (3,'growth_fun_all',&
                          'WARNING 12: Exc 1-3 unexpected exception','','') 
                   ENDIF

                ENDIF


             !! 5.3.3.3.2 Enough leaves to sustain the structural C and roots
             ELSEIF ( ABS(Cl_target(1) - Cl(1)) .LT. min_stomate ) THEN
                
                Cl_incp(1) = MAX(zero, Cl_target(1) - Cl(1))

                ! Enough leaves and structural C, only grow roots
                ! This duplicates Exc 1 and these lines should never be called 
                IF ( ABS(Cs_target(1) - Cs(1)) .LT. min_stomate ) THEN

                   ! Allocate at the tree level to restore allometric balance
                   Cs_incp(1) = MAX(zero, Cs_target(1) - Cs(1))
                   Cr_incp(1) = MAX( MIN(b_inc_tot/ind(ipts,j) - Cl_incp(1) - Cs_incp(1), &
                        Cr_target(1) - Cr(1)), zero )

                   ! Write debug comments to output file 
                   IF ((j.EQ.test_pft .AND. ld_alloc) .OR. ld_warn) THEN
                      CALL comment(npts, Cl_target, Cl, Cs_target, Cs, Cr_target, Cr, delta_ba, ipts, j, l, &
                           b_inc_tot, Cl_incp, Cs_incp, Cr_incp, KF, LF, Cl_inc, Cs_inc, Cr_inc, Cf_inc, &
                           grow_wood, circ_class_n, ind, 15)
                   ENDIF 

                ! Enough leaves and roots. Need to grow structural C to support the available canopy and roots
                ELSEIF ( ABS(Cr_target(1) - Cr(1)) .LT. min_stomate ) THEN

                   Cr_incp(1) = MAX(zero, Cr_target(1) - Cr(1))
                   Cs_incp(1) = MAX( MIN(b_inc_tot/ind(ipts,j) - Cr_incp(1) - Cl_incp(1), &
                        Cs_target(1) - Cs(1)), zero )

                   ! Write debug comments to output file 
                   IF ((j.EQ.test_pft .AND. ld_alloc) .OR. ld_warn) THEN
                      CALL comment(npts, Cl_target, Cl, Cs_target, Cs, Cr_target, Cr, delta_ba, ipts, j, l, &
                           b_inc_tot, Cl_incp, Cs_incp, Cr_incp, KF, LF, Cl_inc, Cs_inc, Cr_inc, Cf_inc, &
                           grow_wood, circ_class_n, ind, 16)
                   ENDIF                  

                ! Need both structural C and roots to restore the allometric relationships
                ELSEIF ( ABS(Cs_target(1) - Cs(1) ) .GT. min_stomate .AND. &
                     ABS(Cr_target(1) - Cr(1)) .GT. min_stomate ) THEN

                   !  First try if we can simply satisfy the allocation needs
                   IF (Cs_target(1) - Cs(1) + Cr_target(1) - Cr(1) .LE. &
                           b_inc_tot/ind(ipts,j) - Cl_incp(1)) THEN
                         
                         Cr_incp(1) = Cr_target(1) - Cr(1)
                         Cs_incp(1) = Cs_target(1) - Cs(1)

                      ! Try to satisfy the need for the roots
                      ELSEIF (Cr_target(1) - Cr(1) .LE. b_inc_tot/ind(ipts,j) - Cl_incp(1)) THEN

                         Cr_incp(1) = Cr_target(1) - Cr(1)
                         Cs_incp(1) = b_inc_tot/ind(ipts,j) - Cl_incp(1) - Cr_incp(1)
                         

                      ! There is not enough use whatever is available
                      ELSE
                         
                         Cr_incp(1) = b_inc_tot/ind(ipts,j) - Cl_incp(1)
                         Cs_incp(1) = zero
                         
                      ENDIF

                   ! Write debug comments to output file 
                   IF ((j.EQ.test_pft .AND. ld_alloc) .OR. ld_warn) THEN
                      CALL comment(npts, Cl_target, Cl, Cs_target, Cs, Cr_target, Cr, delta_ba, ipts, j, l, &
                           b_inc_tot, Cl_incp, Cs_incp, Cr_incp, KF, LF, Cl_inc, Cs_inc, Cr_inc, Cf_inc, &
                           grow_wood, circ_class_n, ind, 17)
                   ENDIF 

                ELSE

                   WRITE(numout,*) 'WARNING 13: Exc 4-6 unexpected exception'
                   WRITE(numout,*) 'WARNING 13: PFT, ipts: ',j,ipts
                   IF(ld_stop)THEN
                      CALL ipslerr_p (3,'growth_fun_all',&
                           'WARNING 13: Exc 4-6 unexpected exception','','')
                   ENDIF

                ENDIF


             !! 5.3.3.3.3 Enough roots to sustain the wood and leaves
             ELSEIF ( ABS(Cr_target(1) - Cr(1)) .LT. min_stomate ) THEN

                Cr_incp(1) = MAX(zero, Cr_target(1) - Cr(1)) 

                ! Enough roots and wood, only grow leaves
                ! This duplicates Exc 2 and these lines should thus never be called 
                IF ( ABS(Cs_target(1) - Cs(1)) .LT. min_stomate ) THEN

                   ! Allocate at the tree level to restore allometric balance
                   Cs_incp(1) = MAX(zero, Cs_target(1) - Cs(1)) 
                   Cl_incp(1) = MAX( MIN(b_inc_tot/ind(ipts,j) - Cr_incp(1) - Cs_incp(1), &
                        Cl_target(1) - Cl(1)), zero )

                   ! Update vegetation height
                   qm_height = (Cl(1) + Cl_incp(1)) * sla_calc(ipts,j) * lai_to_height(j)

                   ! Write debug comments to output file 
                   IF ((j.EQ.test_pft .AND. ld_alloc) .OR. ld_warn) THEN
                      CALL comment(npts, Cl_target, Cl, Cs_target, Cs, Cr_target, Cr, delta_ba, ipts, j, l, &
                           b_inc_tot, Cl_incp, Cs_incp, Cr_incp, KF, LF, Cl_inc, Cs_inc, Cr_inc, Cf_inc, &
                           grow_wood, circ_class_n, ind, 18)
                   ENDIF 

                ! Enough leaves and roots. Need to grow sapwood to support the available canopy and roots
                ! Duplicates Exc. 4 and these lines should thus never be called 
                ELSEIF ( ABS(Cl_target(1) - Cl(1)) .LT. min_stomate ) THEN

                   ! Allocate at the tree level to restore allometric balance
                   Cl_incp(1) = MAX(zero, Cl_target(1) - Cl(1)) 
                   Cs_incp(1) = MAX( MIN(b_inc_tot/ind(ipts,j) - Cr_incp(1) - Cl_incp(1), &
                        Cs_target(1) - Cs(1) ), zero )

                   ! Write debug comments to output file 
                   IF ((j.EQ.test_pft .AND. ld_alloc) .OR. ld_warn) THEN
                      CALL comment(npts, Cl_target, Cl, Cs_target, Cs, Cr_target, Cr, delta_ba, ipts, j, l, &
                           b_inc_tot, Cl_incp, Cs_incp, Cr_incp, KF, LF, Cl_inc, Cs_inc, Cr_inc, Cf_inc, &
                           grow_wood, circ_class_n, ind, 19)
                   ENDIF                       

                ! Need both wood and leaves to restore the allometric relationships
                ELSEIF ( ABS(Cs_target(1) - Cs(1)) .GT. min_stomate .AND. &
                     ABS(Cl_target(1) - Cl(1)) .GT. min_stomate ) THEN

                   ! circ_class_ba_eff and circ_class_height_eff are already calculated
                   ! for a tree in balance. It would be rather complicated to follow
                   ! the allometric rules for wood allocation (implying changes in height 
                   ! and basal area) because the tree is not in balance.First try if we 
                   ! can simply satisfy the allocation needs
                   IF (Cs_target(1) - Cs(1) + Cl_target(1) - Cl(1) .LE. &
                        b_inc_tot/ind(ipts,j) - Cr_incp(1)) THEN
                      
                      Cl_incp(1) = Cl_target(1) - Cl(1)
                      Cs_incp(1) = Cs_target(1) - Cs(1)
                      
                   ! Try to satisfy the need for leaves
                   ELSEIF (Cl_target(1) - Cl(1) .LE. b_inc_tot/ind(ipts,j) - Cr_incp(1)) THEN

                      Cl_incp(1) = Cl_target(1) - Cl(1)
                      Cs_incp(1) = b_inc_tot/ind(ipts,j) - Cr_incp(1) - Cl_incp(1)

                   ! There is not enough use whatever is available
                   ELSE

                      Cl_incp(1) = b_inc_tot/ind(ipts,j) - Cr_incp(1)
                      Cs_incp(1) = zero
                      
                   ENDIF
                      
                   ! Calculate the height of the expanded canopy
                   qm_height(ipts,j) = (Cl(1) + Cl_inc(1)) * sla_calc(ipts,j) * lai_to_height(j)

                   ! Write debug comments to output file 
                   IF ((j.EQ.test_pft .AND. ld_alloc) .OR. ld_warn) THEN
                      CALL comment(npts, Cl_target, Cl, Cs_target, Cs, Cr_target, Cr, delta_ba, ipts, j, l, &
                           b_inc_tot, Cl_incp, Cs_incp, Cr_incp, KF, LF, Cl_inc, Cs_inc, Cr_inc, Cf_inc, &
                           grow_wood, circ_class_n, ind, 20)
                   ENDIF

                ELSE

                   WRITE(numout,*) 'WARNING 14: Exc 7-9 unexpected exception'
                   WRITE(numout,*) 'WARNING 14: PFT, ipts: ',j, ipts
                   IF(ld_stop)THEN
                      CALL ipslerr_p (3,'growth_fun_all',&
                           'WARNING 14: Exc 7-9 unexpected exception','','')
                   ENDIF

                ENDIF

             ! Either Cl_target, Cs_target or Cr_target should be zero
             ELSE

                ! Something possibly important was overlooked
                WRITE(numout,*) 'WARNING 15: Logical flaw in phenological allocation '
                WRITE(numout,*) 'WARNING 15: PFT, ipts: ',j, ipts
                WRITE(numout,*) 'Cs - Cs_target', Cs(1), Cs_target(1)
                WRITE(numout,*) 'Cl - Cl_target', Cl(1), Cl_target(1)
                WRITE(numout,*) 'Cr - Cr_target', Cr(1), Cr_target(1)
                IF(ld_stop)THEN
                   CALL ipslerr_p (3,'growth_fun_all',&
                        'WARNING 15: Logical flaw in phenological allocation','','')
                ENDIF

             ENDIF


             !! 5.3.4 Wrap-up phenological allocation
             IF ( Cl_incp(1) .GE. zero .OR. Cr_incp(1) .GE. zero .OR. Cs_incp(1) .GE. zero) THEN

                ! Fake allocation for less messy equations in next 
                ! case, incp needs to be added to inc at the end
                Cl(1) = Cl(1) + Cl_incp(1)
                Cr(1) = Cr(1) + Cr_incp(1)
                Cs(1) = Cs(1) + Cs_incp(1)
                b_inc_tot = b_inc_tot - (ind(ipts,j) * (Cl_incp(1) + Cr_incp(1) + Cs_incp(1)))

             ELSE

                ! The code was written such that the increment pools should be greater than or equal
                ! to zero. If this is not the case, something fundamental is wrong with the if-then 
                ! constructs under §5.3.3.2
                WRITE(numout,*) 'WARNING 16: numerical problem, one of the increment pools is less than zero'
                WRITE(numout,*) 'WARNING 16: Cl_incp(1), Cr_incp(1), Cs_incp(1), j, ipts',&
                     Cl_incp(1), Cr_incp(1), Cs_incp(1), j, ipts
                IF(ld_stop)THEN
                   CALL ipslerr_p (3,'growth_fun_all',&
                        'WARNING 16: numerical problem, one of the increment pools is less than zero','','')
                ENDIF

             ENDIF

             ! Height depends on Cl, so update height when Cl gets updated
             qm_height(ipts,j) = Cl(1) * sla_calc(ipts,j) * lai_to_height(j) 

             ! Something is wrong with the calculations
             IF (b_inc_tot .LT. -min_stomate) THEN

                WRITE(numout,*) 'WARNING 17: numerical problem overspending in the phenological allocation'
                WRITE(numout,*) 'WARNING 17: b_inc_tot, j, ipts',b_inc_tot, j, ipts 
                WRITE(numout,*) 'WARNING 17: Cl_incp, Cr_incp, Cs_incp, ', Cl_incp(1), Cr_incp(1), Cs_incp(1)
                IF(ld_stop)THEN
                    CALL ipslerr_p (3,'growth_fun_all',&
                         'WARNING 17: numerical problem overspending in the phenological allocation','','')
                ENDIF

             ENDIF

    
             !! 5.3.5 Calculate the expected size of the reserve pool 
             !  use the minimum of either (1) 20% of the total sapwood biomass or 
             !  (2) the amount of carbon needed to develop the optimal LAI and the roots
             !  This reserve pool estimate is only used to decide whether wood should be
             !  grown or not. When really dealing with the reserves the reserve pool is 
             !  recalculated. See further below §7.1.


             reserve_pool = MIN( 0.50 * ( biomass(ipts,j,isapabove,icarbon)+ &
             biomass(ipts,j,iroot,icarbon) + &
             biomass(ipts,j,isapbelow,icarbon)), 1.5*lai_target_longterm(ipts,j)/sla(j)*&  !   had no effect on storage!
             (1.+root_reserve(j)/ltor(ipts,j)) )

             grow_wood = .TRUE.
!!JC Debug
!             IF ((j.EQ.test_pft .AND. ld_alloc .AND. j==test_grid)) THEN
!               WRITE(numout,*) 'biomass(icarbres): ',biomass(ipts,j,icarbres,icarbon)
!               WRITE(numout,*) 'reserve_pool: ',reserve_pool
!               WRITE(numout,*) 'lai_target: ',lai_target(ipts,j)
!               WRITE(numout,*) 'grow_wood: ',grow_wood
!             ENDIF
!!End JC Debug
             ! If the carbohydrate pool is too small, don't grow structural C and 
             ! thus skip ordinary allocation
!JC Debug MOD020 
! here, when target Cl is reached, and reserve_pool is not full
! labile biomass will fill the reserves rather than grow tissue
! However, at the beginning of growing season, there should be no biomass
! been allocated to reserve
! So either assumed no allocation to reserves before reserve_time_grass
! or before devstage = 1 or 0.75 
             IF ( (pheno_type(j) .NE. 1) .AND. &
                  (biomass(ipts,j,icarbres,icarbon) .LT. reserve_pool) &!) THEN
                  .AND. (when_growthinit(ipts,j) .GT. reserve_time_grass) &
                  .AND. (GRM_devstage(ipts,j) .GT. 1.0) ) THEN
!!JC Debug MOD017
!! assuming no ordinary growth limitation for grassea
!                IF (is_tree(j)) THEN
!!JC Debug If icarbres is not full, we should put labile carbon to reserves
!rather than allocate to tissue to avoid deplete icarbres when lab_fac = 0.95
                  grow_wood = .FALSE.
!                ENDIF
!!End JC Debug
             ENDIF
! JC MOD040 crops always use their reserves
             IF (.NOT. natural(j) ) THEN
               grow_wood = .TRUE.
             ENDIF

             ! Write debug comments to output file 
             IF ((j.EQ.test_pft .AND. ld_alloc) .OR. ld_warn) THEN
                CALL comment(npts, Cl_target, Cl, Cs_target, Cs, Cr_target, Cr, &
                     delta_ba, ipts, j, l, b_inc_tot, Cl_incp, Cs_incp, Cr_incp, &
                     KF, LF, Cl_inc, Cs_inc, Cr_inc, Cf_inc, grow_wood, &
                     circ_class_n, ind, 21)
             ENDIF

!JC MOD026 fruit growth started from devear
             IF ((GRM_devstage(ipts,j) .GT. 0.77) .AND. &
                 (b_inc_tot/ind(ipts,j) .GT. min_stomate)) THEN
                  Cf_inc(:) = b_inc_tot * fruit_alloc(j)
                  b_inc_tot = b_inc_tot * (1-fruit_alloc(j))
             ENDIF


             !! 5.3.6 Ordinary growth
             !  Allometric relationship between components is respected, sustain 
             !  ordinary growth and allocate biomass to leaves, wood, roots and fruits.
             IF ( (ABS(Cl_target(1) - Cl(1) ) .LE. min_stomate) .AND. &
                  (ABS(Cs_target(1) - Cs(1) ) .LE. min_stomate) .AND. &
                  (ABS(Cr_target(1) - Cr(1) ) .LE. min_stomate) .AND. &
                  (grow_wood) .AND. (b_inc_tot/ind(ipts,j) .GT. min_stomate) ) THEN 

                ! Allocate fraction of carbon to fruit production (at the tree level)
!JC Debug MOD008
! there should not be allocation to fruit/seed
! at the beginning of growing season 
! thus ordinary allocation to fruit is avoid 
!                  Cf_inc(:) = b_inc_tot * fruit_alloc(j)
!!                IF (when_growthinit(ipts,j) .GT. reserve_time) THEN
!!                  Cf_inc(:) = b_inc_tot * fruit_alloc(j)
!!                ! Residual carbon is allocated to the other components (b_inc_tot is 
!!                ! at the stand level)
!!                  b_inc_tot = b_inc_tot * (1-fruit_alloc(j))
!!                ENDIF
!End JC Debug
                ! Following allometric allocation
                ! (i) b_inc = Cl_inc + Cr_inc + Cs_inc
                ! (ii) Cr_inc = (Cl + Cl_inc)/LF - Cr
                ! (iii) Cs_inc = (Cl + Cl_inc) / KF - Cs
                ! Substitue (ii) and (iii) in (i) and solve for Cl_inc 
                ! <=> b_inc = Cl_inc + ( Cl_inc + Cl ) / KF - Cs + ( Cl_inc + Cl ) / LF - Cr
	        ! <=> b_inc = Cl_inc * ( 1.+ 1/KF + 1./LF ) + Cl/LF - Cs - Cr
	        ! <=> Cl_inc = ( b_inc - Cl/LF + Cs + Cr ) / ( 1.+ 1/KF + 1./LF )
                Cl_inc(1) = MAX( (b_inc_tot/ind(ipts,j) - Cl(1)/LF(ipts,j) - &
                     Cl(1)/KF(ipts,j) + Cs(1) + Cr(1)) / &
                     (1. + 1./KF(ipts,j) + 1./LF(ipts,j)), zero)
               
                IF (Cl_inc(1) .LE. zero) THEN

                   Cr_inc(:) = zero
                   Cs_inc(:) = zero

                ELSE

                   ! Calculate the height of the expanded canopy
                   qm_height(ipts,j) = (Cl(1) + Cl_inc(1)) * sla_calc(ipts,j) * lai_to_height(j)

                   ! Use the solution for Cl_inc to calculate Cr_inc and Cs_inc according to (ii) and (iii)
                   Cr_inc(1) = (Cl(1) + Cl_inc(1)) / LF(ipts,j) - Cr(1)                
                   Cs_inc(1) = (Cl(1)+Cl_inc(1)) / KF(ipts,j) - Cs(1)
               
                ENDIF

                ! Write debug comments to output file 
                IF ((j.EQ.test_pft .AND. ld_alloc) .OR. ld_warn) THEN
                   CALL comment(npts, Cl_target, Cl, Cs_target, Cs, Cr_target, Cr, delta_ba, ipts, j, l, &
                        b_inc_tot, Cl_incp, Cs_incp, Cr_incp, KF, LF, Cl_inc, Cs_inc, Cr_inc, Cf_inc, &
                        grow_wood, circ_class_n, ind, 22)
                ENDIF
                
                ! Wrap-up ordinary growth  
                ! Calculate C that was not allocated, note that Cf_inc was already substracted
                b_inc_tot = b_inc_tot - (ind(ipts,j) * (Cl_inc(1) + Cr_inc(1) + Cs_inc(1)))

                !---TEMP---
                IF (j.EQ.test_pft .AND. ld_alloc) THEN
                   WRITE(numout,*) 'wrap-up ordinary allocation, left b_in_tot, ', b_inc_tot 
                ENDIF
                !----------


             !! 5.3.7 Don't grow wood, use C to fill labile pool
             ELSEIF ( (.NOT. grow_wood) .AND. (b_inc_tot .GT. min_stomate) ) THEN

                ! Calculate the C that needs to be distributed to the 
                ! labile pool. The fraction is proportional to the ratio 
                ! between the total allocatable biomass and the unallocated 
                ! biomass per tree (b_inc now contains the unallocated 
                ! biomass). At the end of the allocation scheme bm_alloc_tot 
                ! is substracted from the labile biomass pool to update the 
                ! biomass pool (biomass(:,:,ilabile) = biomass(:,:,ilabile) - 
                ! bm_alloc_tot(:,:)). At that point, the scheme puts the 
                ! unallocated b_inc into the labile pool. What we 
                ! want is that the unallocated fraction is removed from 
                ! ::bm_alloc_tot such that only the allocated C is removed 
                ! from the labile pool. b_inc_tot will be moved back into
                ! the labile pool in 5.2.11

                bm_alloc_tot(ipts,j) = bm_alloc_tot(ipts,j) - b_inc_tot
                                
                ! Wrap-up ordinary growth  
                ! Calculate C that was not allocated (b_inc_tot), the 
                ! equation should read b_inc_tot = b_inc_tot - b_inc_tot 
                ! note that Cf_inc was already substracted
                b_inc_tot = zero 

                !---TEMP---
                !IF (j.EQ.test_pft .AND. ld_alloc) THEN
                !   WRITE(numout,*) 'No wood growth, move remaining C to labile pool'
                !   WRITE(numout,*) 'bm_alloc_tot_new, ',bm_alloc_tot(ipts,j)
                !   WRITE(numout,*) 'wrap-up ordinary allocation, left b_inc_tot, ', b_inc_tot 
                !ENDIF
                !---------- 

             !! 5.3.8 Error - the allocation scheme is overspending 
             ELSEIF (b_inc_tot .LE. min_stomate) THEN  
      
                IF (b_inc_tot .LT. -min_stomate) THEN

                   ! Something is wrong with the calculations
                   WRITE(numout,*) 'WARNING 18: numerical problem overspending in ordinary allocation'
                   WRITE(numout,*) 'WARNING 18: PFT, ipts, b_inc_tot: ', j, ipts,b_inc_tot
                   IF(ld_stop)THEN
                      CALL ipslerr_p (3,'growth_fun_all',&
                           'WARNING 18: numerical problem overspending in ordinary allocation','','')
                   ENDIF

                ELSE

                   IF (j .EQ. test_pft .AND. ld_alloc) THEN

                      ! Succesful allocation
                      WRITE(numout,*) 'Successful allocation'

                   ENDIF

                ENDIF
              
                ! Althought the biomass components respect the allometric relationships, there
                ! is no carbon left to allocate                      
                b_inc_tot = zero
                Cl_inc(1) = zero
                Cs_inc(1) = zero
                Cr_inc(1) = zero
                Cf_inc(1) = zero

             ELSE

                WRITE(numout,*) 'WARNING 19: Logical flaw unexpected result in ordinary allocation'
                WRITE(numout,*) 'WARNING 19: PFT, ipts: ', j, ipts
                WRITE(numout,*) 'WARNING 19: ',ABS(Cl_target(1) - Cl(1) ) , Cl(1)
                WRITE(numout,*) 'WARNING 19: ',ABS(Cs_target(1) - Cs(1) )  , Cs(1)
                WRITE(numout,*) 'WARNING 19: ',ABS(Cr_target(1) - Cr(1) )  , Cr(1)
                WRITE(numout,*) 'WARNING 19: ',grow_wood
                WRITE(numout,*) 'WARNING 19: ',b_inc_tot,ind(ipts,j),b_inc_tot/ind(ipts,j)
                IF(ld_stop)THEN
                   CALL ipslerr_p (3,'growth_fun_all',&
                        'WARNING 19: Logical flaw unexpected result in ordinary allocation','','')
                ENDIF
                
             ENDIF ! Ordinary allocation


             !! 5.3.9 Forced allocation
             !  Although this should not happen, in case the functional allocation did not consume 
             !  all the allocatable carbon, the remaining C is left for the next day, and some of 
             !  the biomass is used to produce fruits (tuned)
             IF ( b_inc_tot .GT. min_stomate ) THEN 

                WRITE(numout,*) 'WARNING 20: unexpected outcome force allocation'
                WRITE(numout,*) 'WARNING 20: grow_wood, b_inc_tot: ', grow_wood, b_inc_tot
                WRITE(numout,*) 'WARNING 20: PFT, ipts: ',j,ipts
                IF(ld_stop)THEN
                   CALL ipslerr_p (3,'growth_fun_all',&
                        'WARNING 20: unexpected outcome force allocation','','')
                ENDIF

                !+++CHECK+++
!!$                ! Calculate fraction that will be allocated to fruit. The fraction is proportional to the 
!!$                ! ratio between the total allocatable biomass and the unallocated biomass per tree
!!$                frac = 0.1 * MIN(1., bm_alloc_tot(ipts,j) / b_inc_tot)
!!$                Cf_inc(:) = Cf_inc(:) + b_inc_tot * frac
!!$                b_inc_tot = b_inc_tot * (1 - frac)
!!$
!!$                ! Calculate the C that needs to be distributed to the labile pool. The fraction is proportional 
!!$                ! to the ratio between the total allocatable biomass and the unallocated biomass per tree (b_inc 
!!$                ! now contains the unallocated biomass). At the end of the allocation scheme bm_alloc_tot is 
!!$                ! substracted from the labile biomass pool to update the biomass pool (biomass(:,:,ilabile) = 
!!$                ! biomass(:,:,ilabile,icarbon) - bm_alloc_tot(:,:)). At that point, the scheme puts the 
!!$                ! unallocated b_inc into the labile pool. 
!!$                bm_alloc_tot(ipts,j) = bm_alloc_tot(ipts,j) * ( 1. - (1.-frac) * b_inc_tot / bm_alloc_tot(ipts,j) )
                !++++++++++

             ELSEIF ( (b_inc_tot .LT. min_stomate) .AND. (b_inc_tot .GE. -min_stomate) ) THEN

                ! Successful allocation
                !---TEMP---
                IF (j.EQ.test_pft .AND. ld_alloc) THEN
                   WRITE(numout,*) 'Successful allocation'
                ENDIF
                !----------

             ELSE

                ! Something possibly important was overlooked
                IF ( (b_inc_tot .LT. 100*min_stomate) .AND. (b_inc_tot .GE. -100*min_stomate) ) THEN
                   IF (j.EQ.test_pft .AND. ld_alloc) THEN
                      WRITE(numout,*) 'Marginally successful allocation - precision is better than 10-6', j
                   ENDIF
                ELSE
                   WRITE(numout,*) 'WARNING 21: Logical flaw unexpected result in ordinary allocation'
                   WRITE(numout,*) 'WARNING 21: b_inc_tot', b_inc_tot
                   WRITE(numout,*) 'WARNING 21: PFT, ipts: ',j,ipts
                   CALL ipslerr_p (3,'growth_fun_all',&
                        'WARNING 21: Logical flaw unexpected result in ordinary allocation','','')
                ENDIF

             ENDIF

             ! The second problem we need to catch is when one of the increment pools is 
             ! negative. This is an undesired outcome (see comment where ::KF_old is 
             ! calculated in this routine. In that case we write a warning, set all increment
             ! pools to zero and try it again at the next time step. A likely cause of this 
             ! problem is a too large change in KF from one time step to another. Try decreasing
             ! the acceptable value for an absolute increase in KF.
             IF (Cs_inc(1) .LT. zero .OR. Cr_inc(1) .LT. zero .OR. Cs_inc(1) .LT. zero) THEN
             
                ! Do not allocate - save the carbon for the next time step
                Cl_inc(1) = zero
                Cr_inc(1) = zero
                Cs_inc(1) = zero
                WRITE(numout,*) 'WARNING 22: numerical problem, one of the increment pools is less than zero'
                WRITE(numout,*) 'WARNING 22: PFT, ipts: ',j,ipts
               
             ENDIF


             !! 5.3.10 Wrap-up phenological and ordinary allocation
             Cl_inc(1) = Cl_inc(1) + Cl_incp(1)
             Cr_inc(1) = Cr_inc(1) + Cr_incp(1)
             Cs_inc(1) = Cs_inc(1) + Cs_incp(1)
             residual(ipts,j) = b_inc_tot

             !---TEMP---
             IF (j.EQ.test_pft .AND. ld_alloc) THEN
                WRITE(numout,*) 'Final allocation', j 
                WRITE(numout,*) 'Cl, Cs, Cr', Cl(1), Cs(1), Cr(1) 
                WRITE(numout,*) 'Cl_incp, Cs_incp, Cr_incp, ', Cl_incp(1), Cs_incp(1), Cr_incp(1)
                WRITE(numout,*) 'Cl_inc, Cs_ins, Cr_inc, Cf_inc, ', Cl_inc(1), Cs_inc(1), Cr_inc(1), Cf_inc(1)
                WRITE(numout,*) 'unallocated/residual, ', b_inc_tot
             ENDIF
             !----------
            
!===========================================================================
!! This is residual is totally avoidable if one would not have several variables
!for the same thing :(
! DSG: the following bit of code appears twice (grep WARNING 11)
             !! 5.3.11 Account for the residual
             !  The residual is usually around ::min_stomate but we deal 
             !  with it anyway to make sure the mass balance is closed
             !  and as a way to detect errors. Move the unallocated carbon
             !  back into the labile poolA

             deficit = zero 

             IF (ABS((biomass(ipts,j,ilabile,icarbon) -   bm_alloc_tot(ipts,j))        &
                     + residual(ipts,j)) .GT. min_stomate) THEN

                deficit = (biomass(ipts,j,ilabile,icarbon) - bm_alloc_tot(ipts,j) ) &
                            + residual(ipts,j)

                ! The deficit is less than the carbon reserve
                IF (-deficit .LE. biomass(ipts,j,icarbres,icarbon)) THEN

                   ! Pay the deficit from the reserve pool
                   biomass(ipts,j,icarbres,icarbon) = &
                        biomass(ipts,j,icarbres,icarbon) + deficit
                   biomass(ipts,j,ilabile,icarbon)  = &
                        biomass(ipts,j,ilabile,icarbon) - deficit
                
                ! If not, try to reduce it from growth respiration:
                ELSEIF ( ABS(deficit) .LE. (resp_growth(ipts,j) + resp_maint(ipts,j))) THEN

                   biomass(ipts,j,ilabile,icarbon)  = &
                        biomass(ipts,j,ilabile,icarbon) - deficit
! arbitrary removal from respiration components
                   ! remove it from maint_respiration
                   resp_maint(ipts,j) = resp_maint(ipts,j) + deficit
                   ! if deficit is more than maint resp remove the rest from
                   ! growth resp
                   IF (resp_maint(ipts,j) .LT. zero) THEN
                       resp_growth(ipts,j) = resp_growth(ipts,j) + resp_maint(ipts,j)
                       resp_maint(ipts,j)  = zero
                   ENDIF
! arbitrary removal from respiration components

                ELSE
                   ! Not enough carbon to pay the deficit
                   ! There is likely a bigger problem somewhere in
                   ! this routine.
                   WRITE(numout,*) 'WARNING 23: PFT, ipts: ',j,ipts
!DSGtempFIX
                   !  DSG: There isn't a bigger problem in the routine, but the
                   !  routine itself is the problem. 
                   !  This error can be trigger when tuning the model on-the-fly.
                   !  Thus, we won't stop the model but prevent any growth this time
                   !  step. Allocation routine can try next timestep again to
                   !  get things right:
                   bm_alloc_tot(ipts,j) = zero
  !                 CALL ipslerr_p (3,'growth_fun_all',&
  !                      'WARNING 23: numerical problem overspending ',&
  !                      'when trying to account for unallocatable C ','')

!DSGtempFIX
                ENDIF
                
             ENDIF
!===========================================================================
                

             !! 5.3.12 Standardise allocation factors
             !  Strictly speaking the allocation factors do not need to be calculated because the functional
             !  allocation scheme allocates absolute amounts of carbon. Hence, Cl_inc could simply be added to
             !  biomass(:,:,ileaf,icarbon), Cr_inc to biomass(:,:,iroot,icarbon), etc. However, using allocation 
             !  factors bears some elegance in respect to distributing the growth respiration if this would be 
             !  required. Further it facilitates comparison to the resource limited allocation scheme 
             !  (stomate_growth_res_lim.f90) and it comes in handy for model-data comparison. This allocation
             !  takes place at the tree level - note that ::biomass is the only prognostic variable from the tree-based
             !  allocation
                          
             !  Allocation   
             Cl_inc(1) = MAX(zero, ind(ipts,j) * Cl_inc(1))
             Cr_inc(1) = MAX(zero, ind(ipts,j) * Cr_inc(1))
             Cs_inc(1) = MAX(zero, ind(ipts,j) * Cs_inc(1))
             Cf_inc(1) = MAX(zero, ind(ipts,j) * Cf_inc(1))
             
             ! Total_inc is based on the updated Cl_inc, Cr_inc, Cs_inc and Cf_inc. Therefore, do not multiply
             ! ind(ipts,j) again
             total_inc = (Cf_inc(1) + Cl_inc(1) + Cs_inc(1) + Cr_inc(1))
             
             ! Relative allocation
             IF ( total_inc .GT. min_stomate ) THEN

                Cl_inc(1) = Cl_inc(1) / total_inc
                Cs_inc(1) = Cs_inc(1) / total_inc
                Cr_inc(1) = Cr_inc(1) / total_inc
                Cf_inc(1) = Cf_inc(1) / total_inc

             ELSE

                bm_alloc_tot(ipts,j) = zero
                Cl_inc(1) = zero
                Cs_inc(1) = zero
                Cr_inc(1) = zero
                Cf_inc(1) = zero

             ENDIF

 
             !! 5.3.13 Convert allocation to allocation facors
             !  Convert allocation of individuals to ORCHIDEE's allocation factors - see comment for 5.2.5
             !  Aboveground sapwood allocation is age dependent in trees, but there is only aboveground
             !  allocation in grasses 
             alloc_sap_above = un

             ! Leaf, wood, root and fruit allocation
             f_alloc(ipts,j,ileaf)     = Cl_inc(1)
             f_alloc(ipts,j,isapabove) = Cs_inc(1)*alloc_sap_above
             f_alloc(ipts,j,isapbelow) = Cs_inc(1)*(1.-alloc_sap_above)
             f_alloc(ipts,j,iroot)     = Cr_inc(1)
             f_alloc(ipts,j,ifruit)    = Cf_inc(1)

             ! JC ADD for a different allocation scheme used by grasses only
             ! NOTE: it should not be used
             IF (GRM_allow_DEVSTAGE) THEN
               IF ( (biomass(ipts,j,ifruit,icarbon)+biomass(ipts,j,isapabove,icarbon)+ &
                    biomass(ipts,j,ileaf,icarbon)) .GT. min_stomate) THEN
                 wear2wsh = biomass(ipts,j,ifruit,icarbon) / &
                    (biomass(ipts,j,ifruit,icarbon)+biomass(ipts,j,isapabove,icarbon)+ &
                    biomass(ipts,j,ileaf,icarbon))
               ELSE
                 wear2wsh = 0.0
               ENDIF
               IF (GRM_devstage(ipts,j) .GT. 0.0 .AND. GRM_devstage(ipts,j) .LE. 0.77) THEN 
                 lamrep = 1.0 - (1.0 - 0.25) * GRM_devstage(ipts,j)/0.77
                 fdev_lam = 0.3 + (0.8 - 0.3) * (0.77 - GRM_devstage(ipts,j)) / 0.77
                 fdev_ear = 0.0
                 fdev_st = 1.0 - fdev_lam - fdev_ear
               ELSE IF (GRM_devstage(ipts,j) .GT. 0.77 .AND. GRM_devstage(ipts,j) .LE. 1.0) THEN
                 lamrep = 0.25 + (1.0 - 0.25) * (GRM_devstage(ipts,j) - 0.77) / (1.0 - 0.77)
                 fdev_lam = 0.3
                 fdev_ear = MAX(0.0,MIN(0.1, 0.6 * (GRM_devstage(ipts,j) - 0.77) / (1.0 - 0.77) * &
                                (0.1 - wear2wsh)/ 0.1) )
                 fdev_st = 1.0 - fdev_lam - fdev_ear
               ELSE IF (GRM_devstage(ipts,j) .GT. 1.0 .AND. GRM_devstage(ipts,j) .LE. 2.0) THEN
                 lamrep = 1.0
                 fdev_lam = 0.3
                 fdev_ear = MAX(0.0,MIN(0.1, 0.6 * (GRM_devstage(ipts,j) - 0.77) / (1.0 - 0.77) * &
                                (0.1 - wear2wsh)/ 0.1) )
                 fdev_st = 1.0 - fdev_lam - fdev_ear
               ELSE
                 lamrep = 1.0
                 fdev_lam = 0.7
                 fdev_ear = 0.0
                 fdev_st = 1.0 - fdev_lam - fdev_ear
               ENDIF
               fdev_above = 0.6 + 1.0 * (1.0 - lamrep) * 0.4
               fdev_below = 0.4 - 1.0 * (1.0 - lamrep) * 0.4
             f_alloc(ipts,j,ileaf)     = fdev_above * fdev_lam
             f_alloc(ipts,j,isapabove) = fdev_above * fdev_st
             f_alloc(ipts,j,isapbelow) = 0.0
             f_alloc(ipts,j,iroot)     = fdev_below
             f_alloc(ipts,j,ifruit)    = fdev_above * fdev_ear               

             ENDIF

          ELSEIF (.NOT. is_tree(j)) THEN 
       
             IF(ld_alloc) WRITE(numout,*) 'there is no non-tree biomass to allocate, PFT, ', j
             f_alloc(ipts,j,ileaf)     = zero
             f_alloc(ipts,j,isapabove) = zero
             f_alloc(ipts,j,isapbelow) = zero
             f_alloc(ipts,j,iroot)     = zero
             f_alloc(ipts,j,ifruit)    = zero
  
          ENDIF ! .NOT. is_tree(j) and there is biomass to allocate (§5.3 - far far up)

       ENDDO ! npts

       ! correct the biomass to be allocated by the residual.
       bm_alloc_tot(:,j)= bm_alloc_tot(:,j) - residual(:,j)

    IF (dsg_debug) THEN
       IF ( j == test_pft) THEN
            WRITE(numout,*) 'bm_alloc_tot is harmonized with allocation constraints:'
            WRITE(numout,*) 'bm_alloc_tot(test_grid,j)',bm_alloc_tot(test_grid,j)
            WRITE(numout,*) 'residual(test_grid,j)',residual(test_grid,j)
            WRITE(numout,*) 'biomass(test_grid,test_pft,ilabile,icarbon)',biomass(test_grid,test_pft,ilabile,icarbon)
            WRITE(numout,*) 'biomass(test_grid,test_pft,icarbres,icarbon)',biomass(test_grid,test_pft,icarbres,icarbon)
       ENDIF
     ENDIF

! JC Nutrient start
       ! X nutrient limitation
       ! X.1 diagnose how much of the growth can be sustained by the available nutrients

       DO ipts = 1 , npts 

          ! allow easier debugging
          p_avail    = zero
          n_avail    = zero
          costf_N    = zero
          costf_P    = zero
          deltacn    = zero
          deltacnmax = zero
          deltanp    = zero
          deltanpmax = zero
          bm_supply_n = zero
          bm_supply_p = zero

          !DSG useless deltacn=1.0

          ! diagnose nutrient dynamics only if there is growth

          IF ( bm_alloc_tot(ipts,j).GT.min_stomate ) THEN
             ! nitrogen cost, given required N for a given allocatable biomass C, and an
             ! intended leaf CN as N = C * costf_N / C:N
             costf_N = f_alloc(ipts,j,ileaf) +  &
                       fcn_wood_act(ipts,j) * (f_alloc(ipts,j,isapabove) + f_alloc(ipts,j,isapbelow)) + &
                       fcn_root(j) * (f_alloc(ipts,j,iroot)     + f_alloc(ipts,j,ifruit) )

             ! fraction of labile N allocatable for growth
             !no growth respiration calculated here!
             n_avail = & !MIN( &
                       MAX(biomass(ipts,j,ilabile,initrogen)*0.9-min_stomate,zero)!,&

             ! carbon growth which can be sustained by the given nitrogen availability and keeping current leaf C:N stoichiometry
             bm_supply_n = n_avail / costf_N / (un-frac_growthresp_dyn) * cn_leaf(ipts,j) 

             ! elasticity of leaf nitrogen concentration
             !DSG: This if elsif statements are clearly not the nicest way but 
             !     I have to make sure the stoichiometric limits are not violated
             !     by machine precision issues

             ! Elasticity terms tend to give deltanpmax < 1 for cn = cn_max;
             ! thus stoichiometry can run out of range if min_stomate is use
             ! here:
             ! instead we use 0.5
             !IF (ABS(cn_leaf_max(j)-cn_leaf(ipts,j)).LT.min_stomate) THEN

             IF (ABS(cn_leaf_max(j)-cn_leaf(ipts,j)).LT..5) THEN 
                 ! if we at the upper bound
                 deltacnmax = un
             ELSEIF (ABS(cn_leaf_min(j)-cn_leaf(ipts,j)).LT.min_stomate) THEN
                 ! if we at the lower bound
                 deltacnmax = zero
             ELSE
               !DSG: this is different from OCN: 
               deltacnmax=un - exp(-((1.6 * MIN((un/cn_leaf(ipts,j))-(un/cn_leaf_min(j)),zero) / &
                     ( (un/cn_leaf_max(j)) - (un/cn_leaf_min(j)) ) )**4.1))
             ENDIF
             !DSG:end

             !DSG: below is the original elasticity equation;
             ! I do not know the rationale of the change; the term original term
             ! was made up by Soenke  - so it is in general fine to replace it
             ! deltacnmax=exp(-(un/(cn_leaf(ipts,j)*0.5*(un/cn_leaf_max(j)+un/cn_leaf_min(j))))**8)

             ! to exactly the same for phosphorus
             ! phosphorus cost, given required P for a given allocatable biomass C, and an
             ! intended leaf C:P as P = C * costf_P / (C:P) ; C:P = C:N * N:P
             costf_P = f_alloc(ipts,j,ileaf) +                                                            &
                      fcn_wood_act(ipts,j)*fnp_wood_act(ipts,j) * (f_alloc(ipts,j,isapabove) + f_alloc(ipts,j,isapbelow)) + &
                      fcn_root(j)*fnp_root(j) * (f_alloc(ipts,j,iroot)     + f_alloc(ipts,j,ifruit))

             ! fraction of labile P allocatable for growth  (10%); assumption taken from N cycle; "expert choice"
             p_avail = MAX(biomass(ipts,j,ilabile,iphosphorus)*0.9-min_stomate,zero) 
                  
             ! carbon growth which can be sustained by the given phosphorus availability and keeping current leaf C:N:P stoichiometry
             bm_supply_p = p_avail / costf_P / (un-frac_growthresp_dyn) * cn_leaf(ipts,j) * np_leaf(ipts,j)

             ! Elasticity terms tend to give deltanpmax < 1 for np = np_max;
             ! thus stoichiometry can run out of range if min_stomate is use
             ! here:
             ! instead we use 0.2
            ! IF (ABS(np_leaf_max(j)-np_leaf(ipts,j)).LT. min_stomate) THEN 
             IF (ABS(np_leaf_max(j)-np_leaf(ipts,j)).LT. .5) THEN 
                 ! if we at the upper bound
                 deltanpmax = un
           !  ELSEIF (ABS(np_leaf_min(j)-np_leaf(ipts,j)).LT. min_stomate) THEN
             ELSEIF (ABS(np_leaf_min(j)-np_leaf(ipts,j)).LT. .5) THEN
                 ! if we at the lower bound
                 deltanpmax = zero
             ELSE
                ! elasticity of leaf phosphorus concentration (adopted from OCN)
                ! DSG: here is the new formulation of NV


                ! DSG: new term(s) for xtra-wide N:P ratios 5-60:
                ! deltanpmax=un - exp(-((1.3 * MIN((un/np_leaf(ipts,j))-(un/np_leaf_min(j)),zero) / &
                !      ( (un/np_leaf_max(j)) - (un/np_leaf_min(j)) ) )**7.1))
                ! deltanpmax=un - exp(-((1.2 * MIN((un/np_leaf(ipts,j))-(un/np_leaf_min(j)),zero) / &
                !      ( (un/np_leaf_max(j)) - (un/np_leaf_min(j)) ) )**10.1))


                ! more narrow range 12:30 or 12:35
                deltanpmax=un - exp(-((1.6 * MIN((un/np_leaf(ipts,j))-(un/np_leaf_min(j)),zero) / &
                     ( (un/np_leaf_max(j)) - (un/np_leaf_min(j)) ) )**4.1))
                ! more narrow range less felixibility
                ! deltanpmax=un - exp(-((1.6 * MIN((un/np_leaf(ipts,j))-(un/np_leaf_min(j)),zero) / &
                !      ( (un/np_leaf_max(j)) - (un/np_leaf_min(j)) ) )**3.8))

                ! 


                ! DSG: here is the old formulation of SZ:
           !     deltanpmax=exp(-(un/(np_leaf(ipts,j)*0.5*(un/np_leaf_max(j)+un/np_leaf_min(j))))**8)
             ENDIF



             ! we now know the maximum growth sustainable by available nutrients

             ! if we impose a CN or NP ratio we add nutrients if they are not in
             ! sufficient supply
             IF ((impose_cn).AND.( bm_alloc_tot(ipts,j) .GT. bm_supply_n )) THEN

                IF ((printlev>=3).AND.( ipts == test_grid).AND.( j == test_pft)) THEN
                   WRITE(numout,*)  'not enough N, we add N to keep imposed CN concentration'
                ENDIF

                N_support(ipts,j) = N_support(ipts,j)+(bm_alloc_tot(ipts,j)-bm_supply_n)&
                     *costf_N*(un-frac_growthresp_dyn) / cn_leaf(ipts,j)/0.9 + min_stomate
                biomass(ipts,j,ilabile,initrogen)=biomass(ipts,j,ilabile,initrogen) + &
                     (bm_alloc_tot(ipts,j)-bm_supply_n)&
                     *costf_N*(un-frac_growthresp_dyn) / cn_leaf(ipts,j)/0.9 + min_stomate
                n_avail = & !MIN( &
                     MAX(biomass(ipts,j,ilabile,initrogen)*0.9,0.0)!,&
                IF((test_grid == ipts).AND.(test_pft==j)) WRITE(numout,*) "bm_alloc_tot(ipts,j):",bm_alloc_tot(ipts,j)
                IF((test_grid == ipts).AND.(test_pft==j)) WRITE(numout,*) "bm_supply_n:",bm_supply_n
                ! carbon growth possible given nitrogen availability and current nitrogen concentration
                bm_supply_n = bm_alloc_tot(ipts,j)
                
                IF((test_grid == ipts).AND.(test_pft==j)) WRITE(numout,*) "bm_supply_n:",bm_supply_n

            ENDIF

            ! same for phosphorus
            IF ((impose_np).AND.( bm_alloc_tot(ipts,j) .GT. bm_supply_p )) THEN

                IF ((printlev>=3).AND.( ipts == test_grid).AND.( j == test_pft)) THEN
                   WRITE(numout,*)  'not enough P, we add P to keep imposed NP concentration'
                ENDIF

                P_support(ipts,j) = P_support(ipts,j)+(bm_alloc_tot(ipts,j)-bm_supply_p)&
                     *costf_P*(un-frac_growthresp_dyn) / (cn_leaf(ipts,j)*np_leaf(ipts,j))/0.9 & ! needs to be cofirmed
                     + min_stomate 

                biomass(ipts,j,ilabile,iphosphorus)=biomass(ipts,j,ilabile,iphosphorus) + &
                     (bm_alloc_tot(ipts,j)-bm_supply_p)&
                     *costf_P*(un-frac_growthresp_dyn) / (cn_leaf(ipts,j)*np_leaf(ipts,j))/0.9 + min_stomate ! needs to be confirmed
                p_avail = & !MIN( &
                     MAX(biomass(ipts,j,ilabile,iphosphorus)*0.9,zero)

                IF((test_grid == ipts).AND.(test_pft==j)) WRITE(numout,*) "bm_alloc_tot(ipts,j):",bm_alloc_tot(ipts,j)
                ! carbon growth possible given phosphorus availability and current phosphorus concentration
                bm_supply_p = bm_alloc_tot(ipts,j)
                IF((test_grid == ipts).AND.(test_pft==j)) WRITE(numout,*) "bm_supply_p:",bm_supply_p
            ENDIF

            IF ( bm_alloc_tot(ipts,j) .GT. MIN(bm_supply_n, bm_supply_p)) THEN  

                 IF ((printlev>=3).AND.( ipts == test_grid).AND.( j == test_pft)) THEN
                    WRITE(numout,*)  'not enough nutrients'
                 ENDIF
   
                 ! case of NOT enough nitrogen and/or phosphorus to sustain intended growth, 
                 ! reduce carbon allocation to meet nitrogen/phosphorus availability
                 ! taking account of the maximal change of leaf stoichiometry
    
                 ! 1. computed the increase in C to nutrient ratio (here C:N) to reduce nutrient demand per allocated C
                 !  - delta of nitrogen concentrations in response to nitrogen & phosphorus deficit
                 !  - delta of phosphorus concentrations in response to nitrogen & phosphorus deficit
                 !  - cannot change leaf C:N:P in bud burst period
    
    
                 ! delta of nitrogen concetrations in response to nitrogen deficit

                 IF ((printlev>=3).AND.(ipts == test_grid).AND.(j==test_pft)) THEN
                     WRITE(numout,*)  'we try to increase the C:N ratio'
                 ENDIF

                 deltacnmax = Dmax * (un-deltacnmax) 
                 ! DSG: take the minimum of Nsupply/Ndemand and Psupply/Pdemand
                 ! to derive the "wished" change in CN ratio
                 deltacn    = MIN(n_avail /  ( bm_alloc_tot(ipts,j) * (un-frac_growthresp_dyn)       &
                                               * costf_N * un/cn_leaf(ipts,j) )                   , &
                                  p_avail /  ( bm_alloc_tot(ipts,j) * (un-frac_growthresp_dyn)       &
                                               * costf_P * un/(cn_leaf(ipts,j)*np_leaf(ipts,j) )))
                 ! increase CN ratio:                              
                 deltacn    = MIN(MAX(deltacn,un-deltacnmax),un)

                 IF(is_tree(j)) THEN
                    ! update fcn_wood_act 
                    ! as deltacn < 1 here(!), we need to divide to increase fcn_wood:
                    fcn_wood_act(ipts,j) = fcn_wood_act(ipts,j) /deltacn
                    ! update costf_N
                    costf_N = f_alloc(ipts,j,ileaf) +                                &
                           fcn_wood_act(ipts,j) *                                    &
                           (f_alloc(ipts,j,isapabove) + f_alloc(ipts,j,isapbelow)) + &
                           fcn_root(j) * (f_alloc(ipts,j,iroot) + f_alloc(ipts,j,ifruit) )

                    !  update costf_P
                     costf_P = f_alloc(ipts,j,ileaf) +                              &
                          fcn_wood_act(ipts,j)*fnp_wood_act(ipts,j) *       &
                          (f_alloc(ipts,j,isapabove) + f_alloc(ipts,j,isapbelow)) + &
                          fcn_root(j)*fnp_root(j) * (f_alloc(ipts,j,iroot)     + f_alloc(ipts,j,ifruit))
                 ENDIF

                 ! 1.1 nitrogen demand given by nitrogen available or lowered leaf C:N ratio
                 n_alloc_tot(ipts,j) =  MIN( n_avail , & 
                      bm_alloc_tot(ipts,j) * (un-frac_growthresp_dyn) * costf_N * un/cn_leaf(ipts,j)*deltacn )
    
                 ! 1.1 biomass supported by nitrogen availability and lowered leaf C:N ratio:
                 bm_supply_n = MIN( bm_alloc_tot(ipts,j) , &
                      n_alloc_tot(ipts,j) / costf_N / (un-frac_growthresp_dyn) / (un/cn_leaf(ipts,j) * deltacn) ) 
    
                 ! 1.2 phosphorus demand given by phosphorus available or lowered C:N ratio (keeping current leaf N:P ratio) (EQ CHECKED!)
                 p_alloc_tot(ipts,j) =  MIN( p_avail , & 
                      bm_alloc_tot(ipts,j) * (un-frac_growthresp_dyn) * costf_P * un/cn_leaf(ipts,j) * deltacn * un/np_leaf(ipts,j))

                 ! 1.2 biomass supported by phosphorus availability and lowered C:N ratio (keeping current leaf N:P ratio): (EQ CHECKED!)
                 bm_supply_p = MIN( bm_alloc_tot(ipts,j) , &
                            p_alloc_tot(ipts,j) / costf_P / (1.-frac_growthresp_dyn) / (un/cn_leaf(ipts,j)*deltacn *un/np_leaf(ipts,j)) ) 

                 IF((printlev>=3).AND.(ipts==test_grid).AND.(j==test_pft)) THEN
                     WRITE(numout,*) 'deltacn',deltacn
                     WRITE(numout,*) 'deltacnmax',deltacnmax
                     WRITE(numout,*) 'n_alloc_tot(test_grid,test_pft)',n_alloc_tot(ipts,j) 
                     WRITE(numout,*) 'p_alloc_tot(test_grid,test_pft)',p_alloc_tot(ipts,j) 
                     WRITE(numout,*) 'bm_supply_n(test_grid,test_pft)',bm_supply_n
                     WRITE(numout,*) 'bm_supply_p(test_grid,test_pft)',bm_supply_p
                     WRITE(numout,*) 'cn_leaf(test_grid,test_pft)',cn_leaf(ipts,j)
                     WRITE(numout,*) 'np_leaf(test_grid,test_pft)',np_leaf(ipts,j)
                     WRITE(numout,*) 'costf_N',costf_N
                     WRITE(numout,*) 'costf_P',costf_P
                     WRITE(numout,*) 'n_avail', n_avail
                     WRITE(numout,*) 'p_avail', p_avail
                 ENDIF

                 IF (impose_np) THEN
                 ! in this case no variation of NP ratio is allowed
                    deltanp=un 
    
                    ! set growth to minimum which can be supported: 
                    bm_alloc_tot(ipts,j) = MIN(bm_supply_n, bm_supply_p) 
    
                    ! revised nutrient allocation 
                    p_alloc_tot(ipts,j) =  MIN( p_avail , & 
                            bm_alloc_tot(ipts,j) * (un-frac_growthresp_dyn) * costf_P * un/cn_leaf(ipts,j)*deltacn*un/np_leaf(ipts,j)*deltanp)

                    n_alloc_tot(ipts,j) =  MIN( n_avail , & 
                      bm_alloc_tot(ipts,j) * (un-frac_growthresp_dyn) * costf_N * un/cn_leaf(ipts,j)*deltacn)
    
                 ELSE 
                    ! adjustment of N:P stoichiometry allowed (calcuation of changes in the N:P ratio follow the approach of C:N stoichiometric changes)
                        
                    ! First, check which element is more scarce to decide if the
                    ! N:P ratio should be increased or decreased

                    IF ( bm_supply_p .GT. bm_supply_n) THEN  ! N limits more than P: try to decrease N:P ratio

                       IF ((printlev>=4).AND.( ipts == test_grid).AND.( j == test_pft)) THEN
                         WRITE (numout,*) ' we try to decrease N:P ratio (not enough nutrient)' ! which means allocate more P than the current NP ratio dictates
                         WRITE (numout,*) ' BEFORE: p_alloc_tot(ipts,j)',p_alloc_tot(ipts,j)
                         WRITE (numout,*) ' n_alloc_tot(ipts,j)',n_alloc_tot(ipts,j)
                         WRITE (numout,*) 'p_avail',p_avail
                         WRITE (numout,*) 'costf_N',costf_N 
                         WRITE (numout,*) 'costf_P',costf_P
                         WRITE (numout,*) 'bm_supply_p, bm_supply_n',bm_supply_p,bm_supply_n
                       ENDIF

                       ! now, try to decrease N:P ratio (accounting for the an allowed maximum change)
                       ! maximum change
                       deltanpmax = Dmax * deltanpmax 

                       ! use biomass supported by P vs N to scale the change in leaf NP ratio
                       IF ( bm_supply_n .GT. zero ) THEN
                          deltanp    = MIN(bm_supply_p/bm_supply_n,un+deltanpmax)
                       ELSE
                          deltanp    = un+deltanpmax
                       ENDIF

                       ! set growth to the lower - the N supported - growth
                       ! (assuming that the f_allocs are piecewise linear
                       ! with bm_alloc_tot, which is first-order correct: 
                       bm_alloc_tot(ipts,j) = MIN( bm_alloc_tot(ipts,j) , &
                                n_alloc_tot(ipts,j) / (un-frac_growthresp_dyn) / costf_N *cn_leaf(ipts,j)/deltacn)

                       IF(is_tree(j)) THEN
                          ! update fnp_wood_act 
                          ! as deltanp > 1 here(!), we need to divide to decrease fnp_wood_act:
                          fnp_wood_act(ipts,j) = fnp_wood_act(ipts,j)/deltanp

                          costf_P = f_alloc(ipts,j,ileaf) +                                  &
                                   fcn_wood_act(ipts,j)*fnp_wood_act(ipts,j) *               &
                                   (f_alloc(ipts,j,isapabove) + f_alloc(ipts,j,isapbelow)) + &
                                   fcn_root(j)*fnp_root(j) * (f_alloc(ipts,j,iroot)   + f_alloc(ipts,j,ifruit))
                       ENDIF

                       ! revised P allocation (accounting for surplus of available phosphorus to be used to decrease N:P ratio) 
                       p_alloc_tot(ipts,j) =  MIN( p_avail , &
                            bm_alloc_tot(ipts,j) * (un-frac_growthresp_dyn) * costf_P * un/cn_leaf(ipts,j)*deltacn*un/np_leaf(ipts,j)*deltanp)


                       IF ((printlev>=4).AND.( ipts == test_grid).AND.( j == test_pft)) THEN
                         WRITE (numout,*) 'delta_np ,delta_np_max',deltanp,  un+deltanpmax
                         WRITE (numout,*) 'AFTER: p_alloc_tot(ipts,j)',p_alloc_tot(ipts,j)
                         WRITE (numout,*) 'AFTER: costf_P',costf_P
                         WRITE (numout,*) 'AFTER: bm_alloc_tot',  bm_alloc_tot(ipts,j)
                       ENDIF
    
                    ELSEIF ( bm_supply_p .LT. bm_supply_n) THEN ! P limits more than N: try to increase N:P ratio

                       IF ((printlev>=4).AND.( ipts == test_grid).AND.( j == test_pft)) THEN
                         WRITE (numout,*) ' we try to increase the N:P ratio (not enough nutrient)' ! which means allocate more N than the current CN ratio dictates
                         WRITE (numout,*) ' BEFORE: p_alloc_tot(ipts,j)',p_alloc_tot(ipts,j)
                         WRITE (numout,*) ' n_alloc_tot(ipts,j)',n_alloc_tot(ipts,j)
                         WRITE (numout,*) 'p_avail',p_avail
                         WRITE (numout,*) 'costf_N',costf_N 
                         WRITE (numout,*) 'costf_P',costf_P
                         WRITE (numout,*) 'bm_supply_p, bm_supply_n',bm_supply_p,bm_supply_n
                       ENDIF

                       ! now, try to increase N:P ratio (accounting for the an allowed maximum change)
                       ! maximum change

                       !DSG_2lowGPP
                       ! as leaf N is more important than leaf P for PS, the N:P
                       ! ratio should be allowed to increase faster than the C:N
                       !deltanpmax = Dmax * (un-deltanpmax) 
                       deltanpmax  = 2. * Dmax * (un-deltanpmax) 
                       !DSG_2lowGPP

                       ! use biomass supported by P vs N to scale change in NP ratio
                       IF ( bm_supply_p .GT. zero ) THEN
                          deltanp    = MIN(bm_supply_n/bm_supply_p,un+deltanpmax)
                       ELSE
                          deltanp    = un+deltanpmax
                       ENDIF

                       IF(is_tree(j)) THEN
                          ! update fnp_wood_act 
                          ! as deltanp > 1 here(!), we need to multiply to increase fnp_wood_act:
                          fnp_wood_act(ipts,j) = fnp_wood_act(ipts,j)*deltanp
                          
                          !here calculate the costf_P for P costs per C
                          ! ( as done when computing changes in CN ratio) 
                          costf_P = f_alloc(ipts,j,ileaf) +                                  &
                                   fcn_wood_act(ipts,j)*fnp_wood_act(ipts,j) *               &
                                   (f_alloc(ipts,j,isapabove) + f_alloc(ipts,j,isapbelow)) + &
                                   fcn_root(j)*fnp_root(j) * (f_alloc(ipts,j,iroot)   + f_alloc(ipts,j,ifruit))
                       ENDIF
    
                       ! set growth to the lower - the P supported growth
                       ! (assuming that the f_allocs are piecewise linear
                       ! with bm_alloc_tot, which is first-order correct: 
                       bm_alloc_tot(ipts,j) = MIN( bm_alloc_tot(ipts,j) , &
                           p_alloc_tot(ipts,j) / (un-frac_growthresp_dyn) / costf_P * cn_leaf(ipts,j) / deltacn  * np_leaf(ipts,j) *deltanp)
    
                       ! revised N allocation (accounting for surplus of available nitrogen to be used to increase N:P ratio) 
                       n_alloc_tot(ipts,j) =  MIN( n_avail , &
                                               p_alloc_tot(ipts,j)  *costf_N /(un/np_leaf(ipts,j))*deltanp /costf_P)

                       IF ((printlev>=4).AND.( ipts == test_grid).AND.( j == test_pft)) THEN
                         WRITE (numout,*) 'delta_np ,delta_np_max',deltanp,  un+deltanpmax
                         WRITE (numout,*) 'fnp_wood_act(ipts,j)' ,  fnp_wood_act(ipts,j)
                         WRITE (numout,*) 'AFTER: p_alloc_tot(ipts,j)',p_alloc_tot(ipts,j)
                         WRITE (numout,*) 'AFTER: costf_P',costf_P
                         WRITE (numout,*) 'AFTER: bm_alloc_tot',  bm_alloc_tot(ipts,j)
                       ENDIF

                    ELSE !! there are exactly the same (most likely zero) ...
                       ! ... set allocated C to C supported by nutrients
                       bm_alloc_tot(ipts,j) = MIN( bm_alloc_tot(ipts,j) ,bm_supply_n, bm_supply_p)

                       deltanp=un ! just for consistency/debugging
                    ENDIF !( bm_supply_p .GT. bm_supply_n) 

                 ENDIF ! (impose_np)

            ELSE ! DO the following only when neither N nor P is limiting:

                 ! sufficient nitrogen and phosphorus, 
                 ! increase of leaf nutrient concentration dependent on 
                 ! distance to maximal leaf nitrogen concentration
                 ! nitrogen constrained such that nitrogen concentration can only increase 1% per day

                IF (impose_cn) THEN
                   deltacn=un 
                   deltanp=un 
                ELSE

                   IF ((printlev>=3).AND.( ipts == test_grid).AND.( j == test_pft)) THEN
                     WRITE (numout,*) ' we try to decrease the C:N ratio (sufficient nutrients)'
                   ENDIF

                   !DSG_2lowGPP
                   ! as leaf N is more important than leaf P for PS, the N:P
                   ! ratio should be allowed to increase faster than the C:N
                   !deltanpmax = Dmax * deltanpmax
                   deltacnmax = 2. * Dmax * deltacnmax
                   !DSG_2lowGPP

                   ! Use the minimum of Nsupply/Ndemand and Psupply/Pdemand
                   deltacn    = MIN(n_avail /                                                                &
                        ( bm_alloc_tot(ipts,j) * (un-frac_growthresp_dyn) * costf_N * un/cn_leaf(ipts,j) ),   & 
                                    p_avail /                                                                &
                        ( bm_alloc_tot(ipts,j) * (un-frac_growthresp_dyn) * costf_P * un/(cn_leaf(ipts,j)*np_leaf(ipts,j)) ))

                   deltacn    = MIN(deltacn,un+deltacnmax)

                ENDIF

                IF(is_tree(j)) THEN
                   ! update fcn_wood_act 
                   ! as deltacn > 1 here(!), we need to divide to decrease fcn_wood:
                   fcn_wood_act(ipts,j) = fcn_wood_act(ipts,j) /deltacn

                   ! update costf_N
                   costf_N = f_alloc(ipts,j,ileaf) +                                &
                          fcn_wood_act(ipts,j) *                                    &
                          (f_alloc(ipts,j,isapabove) + f_alloc(ipts,j,isapbelow)) + &
                          fcn_root(j) * (f_alloc(ipts,j,iroot) + f_alloc(ipts,j,ifruit) )

                   ! update costf_P
                    costf_P = f_alloc(ipts,j,ileaf) +                              &
                         fcn_wood_act(ipts,j)*fnp_wood_act(ipts,j) *       & 
                         (f_alloc(ipts,j,isapabove) + f_alloc(ipts,j,isapbelow)) + &
                         fcn_root(j)*fnp_root(j) * (f_alloc(ipts,j,iroot)     + f_alloc(ipts,j,ifruit))
                ENDIF

                ! N allocation 
                n_alloc_tot(ipts,j) =  MIN( n_avail , & 
                     bm_alloc_tot(ipts,j) * (un-frac_growthresp_dyn) * costf_N &
                     * un/cn_leaf(ipts,j)*deltacn )

                ! P allocation 
                p_alloc_tot(ipts,j) =  MIN( p_avail , & 
                     bm_alloc_tot(ipts,j) * (un-frac_growthresp_dyn) * costf_P   &
                     * un/cn_leaf(ipts,j)*deltacn   * un/np_leaf(ipts,j) ) 

                IF ((printlev>=4).AND.( ipts == test_grid).AND.( j == test_pft)) THEN
                  WRITE (numout,*) 'delta_cn ,delta_cn_max',deltacn,  un+deltacnmax
                  WRITE (numout,*) 'fcn_wood_act(ipts,j)' ,  fcn_wood_act(ipts,j)
                  WRITE (numout,*) 'AFTER: p_alloc_tot(ipts,j)',p_alloc_tot(ipts,j)
                  WRITE (numout,*) 'AFTER: n_alloc_tot(ipts,j)',n_alloc_tot(ipts,j)
                  WRITE (numout,*) 'AFTER: costf_N',costf_N
                  WRITE (numout,*) 'AFTER: costf_P',costf_P
                ENDIF

                IF(.NOT.impose_np) THEN
                   ! try to correct deviation of leaf_np from an minimum NP (np_leaf_prescribed) 
                   ! here is room for improvement: it would be better to use
                   ! an optimal leaf NP in respect to photosynthesis; as the
                   ! leaf P - photosynthesis dependence is yet not developed, we 
                   ! here let the plant adjust the leaf NP to an average 

!DSGfkn3k_dbg
! The threshold doesn't work in case of very small p_alloc_tot values:
!                   IF ((np_leaf(ipts,j) .GT. np_leaf_min(j)) .AND.   & 
!                   ! also check if p_alloc << p_avail to avoid inconsistencies between p_alloc and costf_P in case new p_alloc >> p_avail
!                       ((p_alloc_tot(ipts,j)-p_avail) .GT. (p_alloc_tot(ipts,j)*2.))) THEN
! Thus for now, we accept small deviations in leaf N:P from inconsistencies between p_alloc and costf_P in case new p_alloc >> p_avail
                   IF (np_leaf(ipts,j) .GT. np_leaf_min(j)) THEN
!DSGfkn3k_dbg

                      IF (( ipts == test_grid).AND.(j == test_pft)) THEN
                         WRITE (numout,*) ' we try to decrease the N:P ratio (sufficient nutrients):'
                      ENDIF 

                      ! use the distance between actual and target NP ratio to
                      ! to scale correction
                      deltanp    = np_leaf(ipts,j)/np_leaf_min(j) ! gets a value larger than 1
                      deltanpmax = Dmax * deltanpmax
                      deltanp    = MIN(deltanp,un+deltanpmax)

                      IF(is_tree(j)) THEN
                         ! update fnp_wood_act 
                         ! as deltanp > 1 here(!), we need to multiply to increase fnp_wood_act:
                         fnp_wood_act(ipts,j) = fnp_wood_act(ipts,j)/deltanp

                         costf_P = f_alloc(ipts,j,ileaf) +                                  &
                                  fcn_wood_act(ipts,j)*fnp_wood_act(ipts,j) *               &
                                  (f_alloc(ipts,j,isapabove) + f_alloc(ipts,j,isapbelow)) + &
                                  fcn_root(j)*fnp_root(j) * (f_alloc(ipts,j,iroot)   + f_alloc(ipts,j,ifruit))
                      ENDIF

                          p_alloc_tot(ipts,j) =  MIN( p_avail , & 
                               bm_alloc_tot(ipts,j) * (un-frac_growthresp_dyn) * costf_P * un/cn_leaf(ipts,j)*deltacn*un/np_leaf(ipts,j)*deltanp)

                      IF ((printlev>=4).AND.( ipts == test_grid).AND.( j == test_pft)) THEN
                        WRITE (numout,*) 'delta_np ,delta_np_max',deltanp,  un+deltanpmax
                        WRITE (numout,*) 'fnp_wood_act(ipts,j)' ,  fnp_wood_act(ipts,j)
                        WRITE (numout,*) 'AFTER: p_alloc_tot(ipts,j)',p_alloc_tot(ipts,j)
                        WRITE (numout,*) 'p_avail(ipts,j)',p_avail
                        WRITE (numout,*) 'AFTER: costf_P',costf_P
                      ENDIF
                   ENDIF
                ENDIF

             ENDIF!( bm_alloc_tot(ipts,j) .GT. MIN(bm_supply_n, bm_supply_p)) 

          ENDIF! ( costf_N.GT.min_stomate ) 

       ENDDO

       IF((printlev>=3).AND.(j == test_pft))THEN
          WRITE(numout,*) 'bm_alloc_tot harmonized with nutrient constraints:'
          WRITE(numout,*) 'bm_alloc_tot(test_grid,test_pft)',bm_alloc_tot(test_grid,j)
          WRITE(numout,*) 'growth_resp(test_grid,test_pft), ', resp_growth(test_grid,j) 
          WRITE(numout,*) 'labile after growth resp is deduced'
          WRITE(numout,*) 'biomass(test_grid,test_pft,ilabile,icarbon)',biomass(test_grid,j,ilabile,icarbon)
       ENDIF

       !! 5.4 Allocate allocatable biomass to different plant compartments
       !  The amount of allocatable biomass to each compartment is a fraction ::f_alloc of the total 
       !  allocatable biomass - see comment for 5.2.6
       DO k = 1, nparts

          bm_alloc(:,j,k,icarbon) = f_alloc(:,j,k) * bm_alloc_tot(:,j) 

       ENDDO

       !---TEMP---
       IF ((printlev>=3).AND.(j .EQ. test_pft)) THEN
                WRITE(numout,*) 'bm_alloc_tot(ipts,j), ', j, bm_alloc_tot(test_grid,j)
                WRITE(numout,*) 'f_alloc(ipts,j,:), ', f_alloc(test_grid,j,:)
                WRITE(numout,*) 'residual(ipts,j), ', residual(test_grid,j)
                WRITE(numout,*) 'bm_alloc(test_grid,j,:,icarbon)',bm_alloc(test_grid,j,:,icarbon)
       ENDIF
       !---------- 

    ENDDO ! # End Loop over # of PFTs   

    ! we need to zero the array for PFT 1, since it has not been calculated but it
    ! is used in implict loops below
    bm_alloc_tot(:,1)       = zero
    bm_alloc(:,1,:,icarbon) = zero

    ! 
    ! 4. calculate nitrogen & phosphorus fluxes associated with biomass growth 
    !

    ! set the allocated C,N & P in a single variable; to minimize code in the
    ! next stepts:
    alloc_tot(:,:,icarbon)     = bm_alloc_tot(:,:)
    alloc_tot(:,:,initrogen)   = n_alloc_tot(:,:)
    alloc_tot(:,:,iphosphorus) = p_alloc_tot(:,:)

    ! this is the case of dynamic allocation of nitrogen taken up from the soil
    ! The principles of this allocation are
    ! 1) nitrogen allocated to the plant tissue is dependent on the labile nitrogen/carbon ratio
    !    i.e. carbon_alloc*(N/C)_labile
    ! 2) the proportions of C/N ratios between different compartments are prescribed as in Hybrid 3
    ! 3) since different to Hybrid, grasses have sapwood (...) the proportions have be adjusted for
    !    grasses, since their tillers are not lignified...

    DO j = 2, nvm
       DO ipts = 1, npts

          ! only allocate nitrogen when there is construction of new biomass
          IF(bm_alloc_tot(ipts,j).GT.min_stomate) THEN

!             IF (j ==test_pft) WRITE(6,*) 'CASE we grow something'

             alloc_c(ipts,:) = zero
             alloc_d(ipts,:) = zero
             alloc_e(ipts,:) = zero

             ! pool sapwood and roots+fruits into each on pool
             sum_sap(ipts)= bm_alloc(ipts,j,isapabove,icarbon)+bm_alloc(ipts,j,isapbelow,icarbon)
             sum_oth(ipts)= bm_alloc(ipts,j,iroot,icarbon)+bm_alloc(ipts,j,ifruit,icarbon)

             IF((sum_sap(ipts)+sum_oth(ipts)) .GT. min_stomate) THEN
                ! in case there is new allocation to leaves
                IF(bm_alloc(ipts,j,ileaf,icarbon).GT.min_stomate) THEN

!                   IF ((j ==test_pft).AND.(ipts==test_grid)) WRITE(6,*) 'CASE we grow leaves'

                ! nitrogen
                   alloc_c(ipts,initrogen) = un/(un+(fcn_wood_act(ipts,j)*sum_sap(ipts)+fcn_root(j)*sum_oth(ipts))/bm_alloc(ipts,j,ileaf,icarbon))
                   alloc_d(ipts,initrogen) = fcn_wood_act(ipts,j)*sum_sap(ipts)/bm_alloc(ipts,j,ileaf,icarbon)*alloc_c(ipts,initrogen)
                   alloc_e(ipts,initrogen) = fcn_root(j)*sum_oth(ipts)/bm_alloc(ipts,j,ileaf,icarbon)*alloc_c(ipts,initrogen)

                ! phosphorus
                   alloc_c(ipts,iphosphorus) =   &
                               un/(un+ ( fcn_wood_act(ipts,j)*fnp_wood_act(ipts,j)*sum_sap(ipts)+  &
                                         fcn_root(j)*fnp_root(j)*sum_oth(ipts) ) &
                                       / bm_alloc(ipts,j,ileaf,icarbon))

                   alloc_d(ipts,iphosphorus) = fcn_wood_act(ipts,j)*fnp_wood_act(ipts,j)* sum_sap(ipts)/       &
                                   bm_alloc(ipts,j,ileaf,icarbon) * alloc_c(ipts,iphosphorus)
                   alloc_e(ipts,iphosphorus) = fcn_root(j)*fnp_root(j)* sum_oth(ipts)/       &
                                   bm_alloc(ipts,j,ileaf,icarbon) * alloc_c(ipts,iphosphorus)
                   ! otherwise add no nutrients to leaves and alot to whatever is constructed
                ELSE
!                   IF (j ==test_pft) WRITE(6,*) 'CASE we dont grow leaves'
                   bm_alloc(ipts,j,ilabile,icarbon) = bm_alloc(ipts,j,ilabile,icarbon) + bm_alloc(ipts,j,ileaf,icarbon)
                   bm_alloc(ipts,j,ileaf  ,icarbon) = zero 

                   alloc_c(ipts,initrogen) = zero

                   alloc_d(ipts,initrogen) = sum_sap(ipts)/(sum_sap(ipts)+fcn_root(j)/fcn_wood_act(ipts,j)*sum_oth(ipts))
                   alloc_e(ipts,initrogen) = un-alloc_d(ipts,initrogen)

                   alloc_c(ipts,iphosphorus) = zero
                   alloc_d(ipts,iphosphorus) = sum_sap(ipts)/(sum_sap(ipts) + &
                                     (fcn_root(j)*fnp_root(j))/(fcn_wood_act(ipts,j)*fnp_wood_act(ipts,j))*sum_oth(ipts))
                   !
                   alloc_e(ipts,iphosphorus) = un-alloc_d(ipts,iphosphorus)
                ENDIF
                ! case of only allocation to leaves!
             ELSEIF(bm_alloc(ipts,j,ileaf,icarbon) .GT. min_stomate)THEN

!                IF(j==test_pft) WRITE(6,*) 'CASE XX2'
                alloc_c(ipts,initrogen)   = un 
                alloc_c(ipts,iphosphorus) = un 
             ELSE
                alloc_tot(ipts,j,:)      = zero
                bm_alloc(ipts,j,:,:)     = zero
             ENDIF


             !calculate allocation of nitrogen and phosphorus
             DO m = initrogen,iphosphorus

                ! a) to leaves:
                bm_alloc(ipts,j,ileaf,m) = alloc_c(ipts,m)*alloc_tot(ipts,j,m)

                ! b) to sapwood:
                IF(sum_sap(ipts).GT.min_stomate) THEN
                    ! aboveground
                    bm_alloc(ipts,j,isapabove,m) = alloc_d(ipts,m) * & 
                         bm_alloc(ipts,j,isapabove,icarbon)/sum_sap(ipts)*alloc_tot(ipts,j,m)
                    ! belowground
                    bm_alloc(ipts,j,isapbelow,m) = alloc_d(ipts,m) * &
                         bm_alloc(ipts,j,isapbelow,icarbon)/sum_sap(ipts)*alloc_tot(ipts,j,m)
                ELSE
                    ! if we dont allocate nutrients we cannot allocate carbon neither
                    bm_alloc(ipts,j,ilabile,icarbon)   = bm_alloc(ipts,j,ilabile,icarbon)   + &
                                                         bm_alloc(ipts,j,isapabove,icarbon) + &
                                                         bm_alloc(ipts,j,isapbelow,icarbon)
                    bm_alloc(ipts,j,isapabove,icarbon) = zero
                    bm_alloc(ipts,j,isapbelow,icarbon) = zero
                ENDIF

                ! c) to fruits and roots:
                IF(sum_oth(ipts).GT.min_stomate) THEN
                   ! roots
                   bm_alloc(ipts,j,iroot,m)     = alloc_e(ipts,m) * & 
                        bm_alloc(ipts,j,iroot,icarbon)    /sum_oth(ipts)*alloc_tot(ipts,j,m)
                   ! fruits
                   bm_alloc(ipts,j,ifruit,m)    = alloc_e(ipts,m) * & 
                        bm_alloc(ipts,j,ifruit,icarbon)   /sum_oth(ipts)*alloc_tot(ipts,j,m)
                ELSE
                    ! if we dont allocate nutrients we cannot allocate carbon neither
                    bm_alloc(ipts,j,ilabile,icarbon)   = bm_alloc(ipts,j,ilabile,icarbon) + &
                                                         bm_alloc(ipts,j,iroot,icarbon)   + &
                                                         bm_alloc(ipts,j,ifruit,icarbon)
                    bm_alloc(ipts,j,iroot,icarbon)   = zero
                    bm_alloc(ipts,j,ifruit,icarbon)  = zero
                ENDIF

                !necessary because bm_alloc_tot can be positive in deciduous, but all C is put into the reserve
                !in this case, all fractions are set to zero, thus no nitrogen is allocated
                !alternatively formulation is minus (c+d+e)*n_alloc_tot
                alloc_tot(ipts,j,m)=(alloc_c(ipts,m)+alloc_d(ipts,m)+alloc_e(ipts,m))*alloc_tot(ipts,j,m)
             ENDDO

          ELSE ! no biomass growth; thus set allocated nutrients to zero

             IF((printlev>=3).AND.(test_pft==j) .AND. ipts==test_grid)THEN
                 WRITE(numout,*)'CASE: we grow nothing'
             ENDIF

             alloc_tot(ipts,j,:)      = zero
             bm_alloc(ipts,j,:,:)     = zero
          ENDIF
       ENDDO !npts

       IF((printlev>=3).AND.(test_pft==j).AND.(ipts==test_grid))THEN
         WRITE(numout,*) 'alloc_c(test_grid,initrogen)=',alloc_c(test_grid,initrogen) 
         WRITE(numout,*) 'alloc_d(test_grid,initrogen)=',alloc_d(test_grid,initrogen) 
         WRITE(numout,*) 'alloc_e(test_grid,initrogen)=',alloc_e(test_grid,initrogen) 
         WRITE(numout,*) 'n_alloc_tot(test_grid,j)=',n_alloc_tot(test_grid,j)

       ENDIF

    ENDDO !nvm 

    IF(dsg_debug) THEN
       !DSG debug: ===================================================
       WRITE(6,*) 'location in code: '
       WRITE(6,*) 'stomate_growth_fun: after nutrient allocation'
       IF (printlev>=4) THEN
          WRITE(numout,*) 'location in code: '
          WRITE(numout,*) 'stomate_growth_fun: after nutrient allocation'
          WRITE(numout,*) 'biomass(test_grid,test_pft,:,icarbon)     : ',  biomass(test_grid,test_pft,:,icarbon)
          WRITE(numout,*) 'biomass(test_grid,test_pft,:,initrogen)   : ',  biomass(test_grid,test_pft,:,initrogen)
          WRITE(numout,*) 'biomass(test_grid,test_pft,:,iphosphorus) : ',  biomass(test_grid,test_pft,:,iphosphorus)
          WRITE(numout,*) 'SUM of bm_alloc must be alloc_tot'
       ENDIF
       CALL check_mass(npts,biomass(:,:,:,:),'stomate_growth_fun: after nutrient allocation')
    ENDIF

    !! 6. Update the biomass with newly allocated biomass after respiration
    ! Finally, we can update 
    ! a) the donor biomass pool

!==   TEMPORAY   ===============================================================================
!DSGdebug: this is really nasty, but I must do it.  Because of calculating everything
! a couple of times separately there is an issue of overspending of labile
! carbon (WARNING 11 and 23); we allow this if its less than min_stomate by violating mass
! conservation. The following MAX function should be removed as soon as the
! allometric calculation are done in a reasonable way.
    biomass(:,:,ilabile,:)   = MAX(biomass(:,:,ilabile,:)   - alloc_tot(:,:,:),zero)
!==   TEMPORAY  ================================================================================

    ! we update the variable bm_alloc_tot as it is use later on in the code; it
    ! would be best to kick this variable out and just use alloc_tot
    bm_alloc_tot(:,:)= alloc_tot(:,:,icarbon)

    ! b) the receiving biomass pools 
    biomass(:,:,:,:)     = biomass(:,:,:,:)     + bm_alloc(:,:,:,:)
    

    IF(dsg_debug)THEN
       !DSG debug: ===================================================
       WRITE(numout,*) 'we added bm_alloc to biomass pools'
       WRITE(numout,*) 'biomass(test_grid,test_pft,:,icarbon)       :', biomass(test_grid,test_pft,:,icarbon)
       WRITE(numout,*) 'biomass(test_grid,test_pft,:,initrogen)     :', biomass(test_grid,test_pft,:,initrogen)
       WRITE(numout,*) 'biomass(test_grid,test_pft,:,iphosphorus)   :', biomass(test_grid,test_pft,:,iphosphorus)
       WRITE(numout,*) 'bm_alloc(test_grid,test_pft,:,icarbon)      :', bm_alloc(test_grid,test_pft,:,icarbon)
       WRITE(numout,*) 'bm_alloc(test_grid,test_pft,:,initrogen)    :', bm_alloc(test_grid,test_pft,:,initrogen)
       WRITE(numout,*) 'bm_alloc(test_grid,test_pft,:,iphosphorus)  :', bm_alloc(test_grid,test_pft,:,iphosphorus)
       WRITE(numout,*) 'fcn_wood_act(test_grid,test_pft)            :', fcn_wood_act(test_grid,test_pft)
       WRITE(numout,*) 'fnp_wood_act(test_grid,test_pft)            :', fnp_wood_act(test_grid,test_pft)

       CALL check_mass(npts,biomass(:,:,:,:),'stomate_growth_fun: we added bm_alloc to biomass pools')
    ENDIF

    !=========================================================== 
    !  DSG: the following part shifts N from leaves & roots to the labile pool and vice versa depending on
    !  the differences in stoichiometry between new biomass and existing  biomass:
    !  High N availability:  in case of incoming leaf biomass with a
    !  higher NC ratio compared to existing leaf biomass, the N content of
    !  existing leaves is tried to be increased to decrease the difference in
    !  stoichiometry between existing leaves and new leaves by 5%.  
    !  DSG: the relaxation is done to allow plants to adjust their leave stoichiometry
    !  during the vegetation in case nutrient stress changes

!DSGfkn4k_dbg
!DSG: this doesn't work (and it is not needed)
!    IF( Dmax.GT.zero )THEN
!       DO j=2,nvm
!          ! where there are leaf and leaf growth
!          WHERE((biomass(:,j,ileaf,icarbon).GT.min_stomate).AND.(bm_alloc(:,j,ileaf,initrogen).GT.min_stomate))
!             ! calculate the amount of N needed to relax the difference between
!             ! new and existing leaf C:N ratio by 5%
!             transloc_N(:,j) = biomass(:,j,ileaf,icarbon) * 0.05 * & 
!                  (bm_alloc(:,j,ileaf,initrogen)/bm_alloc(:,j,ileaf,icarbon) - & 
!                   un/cn_leaf(:,j))
!             transloc_N(:,j) = MAX(biomass(:,j,ileaf,icarbon)/cn_leaf_max(j)-biomass(:,j,ileaf,initrogen) &
!                                   ,transloc_N(:,j))
!             transloc_N(:,j) = MIN(biomass(:,j,ileaf,icarbon)/cn_leaf_min(j)-biomass(:,j,ileaf,initrogen) &
!                                   ,transloc_N(:,j))
!
!             ! the size of labile pool is the limit to the N being shifted - in case N is added to leaves (transloc_N > zero)
!             transloc_N(:,j) = MIN(biomass(:,j,ilabile,initrogen)*0.7,transloc_N(:,j))
!
!             ! check if there is enough labile P to ensure the N:P of leaves doesn't change; 
!             ! P is a higher turnover in cells than N, so there is no issue with
!             ! letting it passively follow N changes
!             transloc_N(:,j) = MIN(biomass(:,j,ilabile,iphosphorus)*0.7*np_leaf(:,j), transloc_N(:,j))
!
!             ! add/remove the nitrogen to/from the leaf pool ...
!             biomass(:,j,ileaf,initrogen)   = biomass(:,j,ileaf,initrogen)   + transloc_N(:,j) 
!             ! remove/add it from/to the labile pool 
!             biomass(:,j,ilabile,initrogen) = biomass(:,j,ilabile,initrogen) - transloc_N(:,j)
!             ! update the variable containing nitrogen allocated to leaves
!             bm_alloc(:,j,ileaf,initrogen)  = bm_alloc(:,j,ileaf,initrogen)  + transloc_N(:,j) 
!
!
!             ! take care about P in a two step way 
!             ! first, shift P assuming no change in N:P ratio
!             ! second, shift P assuming N:P can adjust, too (the coding approach
!             ! with "first" and "second" might be not the most elegant one, but
!             ! it is more easy to follow as the calculations are similiar)
!
!             ! first, shift P assuming no change in N:P ratio
!             ! calculate the amount of to be shifted
!             transloc_P(:,j) = transloc_N(:,j)/np_leaf(:,j)
!             ! add/remove the phosphorus to/from the leaf pool ...
!             biomass(:,j,ileaf,iphosphorus)   = biomass(:,j,ileaf,iphosphorus)   + transloc_P(:,j) 
!             ! remove/add it from/to the labile pool 
!             biomass(:,j,ilabile,iphosphorus) = biomass(:,j,ilabile,iphosphorus) - transloc_P(:,j)
!             ! update the variable containing phosphorus allocated to leaves
!             bm_alloc(:,j,ileaf,iphosphorus)  = bm_alloc(:,j,ileaf,iphosphorus)  + transloc_P(:,j) 
!
!             ! second, shift P assuming N:P can adjust, too (the coding approach
!             ! calculate the amount of P needed to relax the difference between
!             ! new and existing leaf N:P ratio by 5% (P has turnover is a higher
!             ! than N)
!             transloc_P(:,j) = biomass(:,j,ileaf,initrogen) * 0.05 * & 
!                  (bm_alloc(:,j,ileaf,iphosphorus)/bm_alloc(:,j,ileaf,initrogen) - & 
!                   un/np_leaf(:,j))
!
!             !prevent N:P > NP_MAX & prevent overuse of labile P
!             transloc_P(:,j) = MIN(MAX(biomass(:,j,ileaf,initrogen)/np_leaf_max(j)-biomass(:,j,ileaf,iphosphorus) &
!                                   ,transloc_P(:,j)),biomass(:,j,ilabile,iphosphorus)*0.7)
!             !prevent N:P < NP_MIN
!             transloc_P(:,j) = MIN(biomass(:,j,ileaf,initrogen)/np_leaf_min(j)-biomass(:,j,ileaf,iphosphorus) &
!                                   ,transloc_P(:,j))
!
!             ! add/remove the phosphorus to/from the leaf pool ...
!             biomass(:,j,ileaf,iphosphorus)   = biomass(:,j,ileaf,iphosphorus)   + transloc_P(:,j) 
!             ! remove/add it from/to the labile pool 
!             biomass(:,j,ilabile,iphosphorus) = biomass(:,j,ilabile,iphosphorus) - transloc_P(:,j)
!             ! update the variable containing phosphorus allocated to leaves
!             bm_alloc(:,j,ileaf,iphosphorus)  = bm_alloc(:,j,ileaf,iphosphorus)  + transloc_P(:,j) 
!
!          ENDWHERE
!          
!          ! we do exactly the same with the nutrient concentration of roots
!          WHERE((biomass(:,j,iroot,icarbon).GT.min_stomate).AND.(bm_alloc(:,j,iroot,initrogen).GT.min_stomate))
!             transloc_N(:,j) = biomass(:,j,iroot,icarbon) * 0.05 * & 
!                  (fcn_root(j)/cn_leaf(:,j)-biomass(:,j,iroot,initrogen)/biomass(:,j,iroot,icarbon))
!             transloc_N(:,j) = MIN(biomass(:,j,iroot,icarbon)/cn_leaf_min(j)*fcn_root(j) &
!                                   - biomass(:,j,iroot,initrogen)                        &
!                                    ,transloc_N(:,j))
!             transloc_N(:,j) = MAX(biomass(:,j,iroot,icarbon)/cn_leaf_max(j)*fcn_root(j) &
!                                   - biomass(:,j,iroot,initrogen)                        &
!                                     ,transloc_N(:,j))
!
!             ! the size of labile pool is the limit to the N being shifted - in case N is added to roots (transloc_N > zero)
!             transloc_N(:,j) = MIN(biomass(:,j,ilabile,initrogen)*0.7,transloc_N(:,j))
!
!             ! check if there is enough labile P to ensure the N:P doesn't change
!             transloc_N(:,j) = MIN(biomass(:,j,ilabile,iphosphorus)*0.7/fnp_root(j)*np_leaf(:,j), transloc_N(:,j))
!             transloc_P(:,j) = transloc_N(:,j)/np_leaf(:,j)*fnp_root(j)
!             ! add/remove the phosphorus to/from the leaf pool ...
!             biomass(:,j,iroot,iphosphorus)   = biomass(:,j,iroot,iphosphorus)   + transloc_P(:,j) 
!             ! remove/add it from/to the labile pool 
!             biomass(:,j,ilabile,iphosphorus) = biomass(:,j,ilabile,iphosphorus) - transloc_P(:,j)
!             ! update the variable containing phosphorus allocated to leaves
!             bm_alloc(:,j,iroot,iphosphorus)  = bm_alloc(:,j,iroot,iphosphorus)  + transloc_P(:,j) 
!
!             ! add the nitrogen to the root biomass
!             biomass(:,j,iroot,initrogen)  = biomass(:,j,iroot,initrogen) +transloc_N(:,j) 
!             ! ... and remove it from the labile nitrogen
!             biomass(:,j,ilabile,initrogen)= biomass(:,j,ilabile,initrogen)-transloc_N(:,j)
!             ! update the variable containing nitrogen allocated to roots
!             bm_alloc(:,j,iroot,initrogen) = bm_alloc(:,j,iroot,initrogen)+transloc_N(:,j)
!
!             ! relax the N:P by add/removing P (P has a higher turnover than N)
!             transloc_P(:,j) = biomass(:,j,iroot,initrogen) * 0.05 * & 
!                  (bm_alloc(:,j,iroot,iphosphorus)/bm_alloc(:,j,iroot,initrogen) - & 
!                   un/np_leaf(:,j))
!
!             !prevent N:P > NP_MAX & do not overuse labile P
!             transloc_P(:,j) = MIN(MAX(biomass(:,j,iroot,initrogen)/np_leaf_max(j)*fnp_root(j)-biomass(:,j,iroot,iphosphorus) &
!                                   ,transloc_P(:,j)),biomass(:,j,ilabile,iphosphorus)*0.7)
!
!             !prevent N:P < NP_MIN
!             transloc_P(:,j) = MIN(biomass(:,j,iroot,initrogen)/np_leaf_min(j)*fnp_root(j)-biomass(:,j,iroot,iphosphorus) &
!                                   ,transloc_P(:,j))
!
!
!             ! add/remove the phosphorus to/from the root pool ...
!             biomass(:,j,iroot,iphosphorus)   = biomass(:,j,iroot,iphosphorus)   + transloc_P(:,j) 
!             ! remove/add it from/to the labile pool 
!             biomass(:,j,ilabile,iphosphorus) = biomass(:,j,ilabile,iphosphorus) - transloc_P(:,j)
!             ! update the variable containing phosphorus allocated to leaves
!             bm_alloc(:,j,iroot,iphosphorus)  = bm_alloc(:,j,iroot,iphosphorus)  + transloc_P(:,j) 
!          ENDWHERE
!          ! We should do the same for fruits to make sure fruit CNP is the same
!          ! as root (needed for allocation)
!
!       ENDDO
!    ENDIF
!DSGfkn4k_dbg

    !DSG_debug_X6 : 
    IF(dsg_debug) THEN
      !DSG debug: ===================================================
       WRITE(6,*) 'location in code: '
       WRITE(6,*) 'stomate_growth_fun: after shifting'
       IF (printlev>=4) THEN
         WRITE(numout,*) 'location in code: '
         WRITE(numout,*) 'stomate_growth_fun: after shifting'
         WRITE(numout,*) 'biomass(test_grid,test_pft,:,icarbon)       : ',  biomass(test_grid   ,test_pft,:,icarbon)
         WRITE(numout,*) 'biomass(test_grid,test_pft,:,initrogen)     : ',  biomass(test_grid   ,test_pft,:,initrogen)
         WRITE(numout,*) 'biomass(test_grid,test_pft,:,iphosphorus)   : ',  biomass(test_grid   ,test_pft,:,iphosphorus)
         WRITE(numout,*) 'transloc_N(test_grid,test_pft)              : ',  transloc_N(test_grid,test_pft)
       ENDIF
       !WRITE(numout,*) 'fcn_root(j)*cn_leaf(:,j)',fcn_root(test_pft)*un/cn_leaf(test_grid,test_pft)
       CALL check_mass(npts,biomass(:,:,:,:),'stomate_growth_fun: after shift nutrients')
    ENDIF

! JC tree grass start labile vs. reserve pool
 !! 7. Use or fill reserve pools depending on relative size of the labile and reserve C pool

    ! +++ CHECK +++
    ! Externalize all the hard coded values i.e. 0.3

    ! Calculate the labile pool for all plants and also the reserve pool for trees      
    DO j = 2,nvm

       DO ipts = 1,npts
         
          ! set them to val_exp to faciliate debugging output
          labile_pool   = val_exp 
          reserve_pool  = val_exp 
          excess        = val_exp 
          use_max       = val_exp 
          use_lab       = val_exp 

         IF (printlev>=4 .AND. ipts==test_grid .AND. j==test_pft) THEN
             WRITE(numout,*) 'location in code: '
             WRITE(numout,*) 'stomate: growth_fun: 6'
             WRITE(numout,*) 'ipts=',ipts
             WRITE(numout,*) 'pft=',j
             WRITE(numout,*) 'biomass(ipts,j,:,icarbon)    : ', biomass(ipts,j,:,icarbon)
         ENDIF


          IF (veget_max(ipts,j) .LE. min_stomate .OR. &
              SUM(biomass(ipts,j,:,icarbon)) .LE. min_stomate) THEN

             ! this vegetation type is not present, so no reason to do the 
             ! calculation. CYCLE will take us out of the innermost DO loop
             CYCLE

          ENDIF

          ! There is vegetation present and has started growing. The second and third condition
          ! required to make the PFT survive the first year during which the long term climate
          ! variables are initialized for the phenology. If these conditions are not added, the 
          ! reserves are respired well before growth ever starts
          IF ( veget_max(ipts,j) .GT. min_stomate .AND. &
               rue_longterm(ipts,j) .GE. zero .AND. &                  
               rue_longterm(ipts,j) .NE. un) THEN
             
             !! 7.1 Calculate the optimal size of the pools   
             !  The size of the labile pool is proportional to the assumed activity of living tissues 
             !  and its relative nitrogen content). The numerical value of ::lab_fac is already a
             !  tuning variable, the division by 10 (default for ::labile_reserve) stresses the importance 
             !  of this variable to scale other processes.
             !+++CHECK+++
             ! There is an inconsistency in the calculation - most pools are in gN but leaves is in gC
             ! The correction is proposed, that implies that the parameter labile_reserve will need to 
             ! be tuned
!!$          VERSION WITH CONSISTENT UNITS
!!$             labile_pool = lab_fac(ipts,j)/labile_reserve * &
!!$                  ( biomass(ipts,j,ileaf,icarbon) / cn_leaf_prescribed(j) + & 
!!$                  fcn_root(j) * ( biomass(ipts,j,iroot,icarbon) + biomass(ipts,j,ifruit,icarbon) ) + & 
!!$                  fcn_wood(j) * ( biomass(ipts,j,isapabove,icarbon) + biomass(ipts,j,isapbelow,icarbon) + &
!!$                  biomass(ipts,j,icarbres,icarbon) ) )
!!$          ORIGINAL VERSION WITH INCONSISTENT UNITS
!!$             labile_pool = lab_fac(ipts,j)/labile_reserve * ( biomass(ipts,j,ileaf,icarbon)  + & 
!!$                  fcn_root(j) * ( biomass(ipts,j,iroot,icarbon) + biomass(ipts,j,ifruit,icarbon) ) + & 
!!$                  fcn_wood(j) * ( biomass(ipts,j,isapabove,icarbon) + biomass(ipts,j,isapbelow,icarbon) + &
!!$                  biomass(ipts,j,icarbres,icarbon) ) )
!!$             labile_pool = MAX ( labile_pool, gpp_to_labile(j) * gpp_week(ipts,j) )
             
             ! We had an endless series of problems which were often difficult to
             ! understand but which always seemed to be related to a sudden drop
             ! in biomass(ilabile). This drop was often the result of a sudden 
             ! change in labile_pool. 
!DSG:        these drops likely occured due to the way the if cases were triggered for
!            the transfer from labile to reserve, vice versa, and to excess
!            respiraiton.
!DSG: 
             !Given that there is not much science behind
             ! this approach it seems a good idea to remove this max statement to
             ! avoid sudden changes. Rather than using the actual biomass we propose
             ! to use the target biomass. This assumes that the tree would like to 
             ! fill its labile pool to be optimal when it would be in allometric
             ! balance.
             IF (is_tree(j)) THEN
             
                ! We will make use of the REAL sapwood, heartwood and effective height
                ! and then calculate the target leaves and roots. This approach gives
                ! us a target for a labile_pool of a tree in allometric balance. 
                ! Basal area at the tree level (m2 tree-1)
                circ_class_ba_eff(:) = wood_to_ba_eff(circ_class_biomass(ipts,j,1,:,icarbon),j)

                ! Current biomass pools per tree (gC tree^-1) 
                ! We will have different trees so this has to be calculated from the 
                ! diameter relationships            
                Cs(:) = ( circ_class_biomass(ipts,j,:,isapabove,icarbon) + &
                     circ_class_biomass(ipts,j,:,isapbelow,icarbon) ) * scal(ipts,j)

                DO l = 1,ncirc 

                   !  Calculate tree height
                   circ_class_height_eff(l) = pipe_tune2(j)*(4/pi*circ_class_ba_eff(l))**(pipe_tune3(j)/2)
                   
                   ! Use the pipe model to calculate the target leaf and root
                   ! biomasses
                   Cl_target(l) = KF(ipts,j) * Cs(l) / circ_class_height_eff(l)
                   Cr_target(l) = Cl_target(l) / LF(ipts,j)
                   
                ENDDO

             ! grasses and crops
             ELSEIF ( .NOT. is_tree(j)) THEN
             
                Cs(:) = zero
                Cl_target(:) = zero 
                Cr_target(:) = zero

                ! Current biomass pools per grass/crop (gC ind^-1)
                ! Cs has too many dimensions for grass/crops. To have a consistent notation the same variables
                ! are used as for trees but the dimension of Cs, Cl and Cr i.e. ::ncirc should be ignored            
                Cs(1) = biomass(ipts,j,isapabove,icarbon) * scal(ipts,j)
   
                ! Use the pipe model to calculate the target leaf and root
                ! biomasses
                Cl_target(1) = Cs(1) * KF(ipts,j)
                Cr_target(1) = Cl_target(1) / LF(ipts,j)

             ENDIF !is_tree


             IF (j .EQ. test_pft .AND. ld_alloc .AND. ipts == test_grid) THEN
                WRITE(numout,*) 'labile_reserve(j)=',labile_reserve(j)
                WRITE(numout,*) 'Cl_target(:)=',Cl_target(:)
                WRITE(numout,*) 'Cr_target(:)=',Cr_target(:)
                WRITE(numout,*) 'cn_leaf_avg_season(ipts,j)=',cn_leaf_avg_season(ipts,j)
                WRITE(numout,*) 'labile_pool_1',(SUM(Cl_target(:)) + &
                        SUM(Cr_target(:))* fcn_root(j))  &
                        /cn_leaf_avg_season(ipts,j)
                WRITE(numout,*) 'gpp_week(ipts,j)',gpp_week(ipts,j)
                WRITE(numout,*) 'lab_fac, labile pool, ', lab_fac(ipts,j),labile_pool
                WRITE(numout,*) 'circ_class_biomass(ipts,j,:,isapbelow,icarbon)=',circ_class_biomass(ipts,j,:,isapbelow,icarbon)
                WRITE(numout,*) 'circ_class_biomass(ipts,j,:,isapabove,icarbon)=',circ_class_biomass(ipts,j,:,isapabove,icarbon)
                WRITE(numout,*) 'scal=',scal(ipts,j)
                WRITE(numout,*) 'ind(ipts,j)=',ind(ipts,j)
                WRITE(numout,*) 'labile_reserve(j)=',labile_reserve(j)
                WRITE(numout,*) 'Cl_target(:)=',Cl_target(:)
                WRITE(numout,*) 'Cr_target(:)=',Cr_target(:)
                WRITE(numout,*) 'Cs(:)=',Cs(:)
                WRITE(numout,*) 'fcn_root(j)=',fcn_root(j)
                WRITE(numout,*) 'fcn_wood_act(j)=',fcn_wood_act(ipts,j)
                WRITE(numout,*) 'cn_leaf(ipts,j)=',cn_leaf(ipts,j)
                WRITE(numout,*) 'lab_fac, labile pool, ', lab_fac(ipts,j), labile_pool
             ENDIF

                 labile_pool = lab_fac(ipts,j) / un    * & !MAX(lab_fac(ipts,j) / un    *       &
                       ( SUM(Cl_target(:)) + SUM(Cr_target(:))* fcn_root(j))  &
                       !the GPP threshold might be a problem (10* seems to high)
                        /cn_leaf_avg_season(ipts,j)  !, 10. * gpp_week(ipts,j))

             ! The max size of reserve pool is proportional to the size of the storage organ (the sapwood)
             ! and a the leaf functional trait of the PFT (::phene_type_tab). The reserve pool is 
             ! constrained by the mass needed to replace foliage and roots. This constraint prevents the
             ! scheme from putting too much reserves in big trees (which have a lot of sapwood compared to
             ! small trees). Exessive storage would hamper tree growth and would make mortality less likely.
             IF(is_tree(j)) THEN

                IF (pheno_type(j).EQ.1) THEN 

                   ! Evergreen trees are not very conservative with respect to C-storage. Therefore, only maximum of 5% 
                   ! of their sapwood mass is stored in their reserve pool.
                   ! DSGlabXX: I agree. no change fro evergreens
                   reserve_pool = MIN(evergreen_reserve(j) * ( biomass(ipts,j,isapabove,icarbon) + &
                       biomass(ipts,j,isapbelow,icarbon)), 0.25*lai_target_longterm(ipts,j)/sla(j)*(1.+root_reserve(j)/ltor(ipts,j)))
                ELSE
                   ! Deciduous trees are more conservative and 12% of their sapwood mass is stored in the 
                   ! reserve pool. The scheme avoids that during the growing season too much reserve are 
                   ! accumulated (which would hamper growth), therefore, the reduced rate of 12% is used 
                   ! until scenecence.
                   IF (bm_alloc_tot(ipts,j) .GT. min_stomate) THEN
                      reserve_pool = MIN(deciduous_reserve(j) * ( biomass(ipts,j,isapabove,icarbon) + &
                           biomass(ipts,j,isapbelow,icarbon)), lai_target_longterm(ipts,j)/sla(j)*(1.+root_reserve(j)/ltor(ipts,j)))

                  ELSE  
                      ! If the plant is scenecent, allow for a higher reserve mass. Plants can then use the 
                      ! excess labile C, that is no longer used for growth and would be respired otherwise,
                      ! to regrow leaves after the dormant period.
                      reserve_pool = MIN(senescense_reserve(j) * ( biomass(ipts,j,isapabove,icarbon) + &
                           biomass(ipts,j,isapbelow,icarbon)), lai_target_longterm(ipts,j)/sla(j)*(1.+root_reserve(j)/ltor(ipts,j)))
                   ENDIF ! Scenecent

                ENDIF ! Phenology type

!JC GRASS reserve_pool
             ! Grasses and Crops
             ELSE

                !+++CHECK+++
                ! The min criterion results in the reserves being zero because isapabove goes to zero 
                ! when the reserves are most needed
!JC Debug MOD005
                   ! Grasses conservertive: 50% of their
                   ! root and sapwood mass is stored in the
                   ! reserve pool. 
                   ! usually the reserve_pool is limited by root+sap biomass
                   ! during growing season
                   ! here we reset reserve_pool reserve_time_grass days
                   ! after regrowth

                 IF (natural(j)) THEN
                 reserve_pool = MIN(0.5 * ( biomass(ipts,j,iroot,icarbon) + biomass(ipts,j,isapabove,icarbon) + & 
!DSGlaitarget: although grass is often annual, I use the season lai_target to
!avoid the drops in NPP, the grass community itself has a memory via seeds and
!we are prescribing community behaviour and not of a single plant

                     biomass(ipts,j,isapbelow,icarbon)), 1.5*lai_target_longterm(ipts,j)/sla(j)*(1.+root_reserve(j)/ltor(ipts,j)))
                 ELSE
                   reserve_pool = MIN(0.5 * ( biomass(ipts,j,iroot,icarbon) + &
                                  biomass(ipts,j,isapabove,icarbon) + &
                                  biomass(ipts,j,isapbelow,icarbon)), &
                                  MAX(2.0,lai_target_longterm(ipts,j))/sla(j)*(1.+root_reserve(j)/ltor(ipts,j)))

                 ENDIF
!JC Debug MOD005
! opt1 add a conservative variable recording reserve_pool?
! opt2 conserve root and sapabove biomass during winter to mimic the perennual
! grassland?
               ! reset record_reserve_grass at the reserve_time_grass days
               ! after grass leaf-onset (BUT how about regrowth after cut?)
               IF (when_growthinit(ipts,j) .EQ. reserve_time_grass) THEN
                  record_reserve_grass(ipts,j) = reserve_pool
               ENDIF

               IF (record_reserve_grass(ipts,j) .LT. reserve_pool) THEN
                   record_reserve_grass(ipts,j) = reserve_pool             
               ELSE
                  ! do nothing if record is the max
                  record_reserve_grass(ipts,j) = record_reserve_grass(ipts,j)
               ENDIF
               ! normal situation: choose maximum of record_reserve_grass and
               ! reserve_pool
                 ! force reserve_pool to be the maximum of record_reserve_grass
                 ! and reserve_pool just calculated
                 ! which means we do not change reserve_pool anymore when
                 ! senecence start
                 reserve_pool = MAX(record_reserve_grass(ipts,j),reserve_pool)

                 !+++TEMP+++ 
                 IF (j .EQ. test_pft .AND. ld_alloc .AND. ipts == test_grid) THEN
                    WRITE(numout,*) 'bm_alloc_tot, ',bm_alloc_tot(ipts,j)
!                    WRITE(numout,*) 'lai_target, ',lai_target(ipts,j)
                    WRITE(numout,*) '50% r+s, ',0.5 * (biomass(ipts,j,iroot,icarbon) +biomass(ipts,j,isapabove,icarbon) +biomass(ipts,j,isapbelow,icarbon))
                    WRITE(numout,*) 'target2,', lai_target_longterm(ipts,j)/sla(j)*(1.+root_reserve(j)/ltor(ipts,j))
                    WRITE(numout,*) 'reserve pool, 50%', reserve_pool
                 ENDIF
                 !++++++++++

                 IF (j .EQ. test_pft .AND. ld_alloc .AND. ipts == test_grid &
                    .AND. record_reserve_grass(ipts,j) .GT. reserve_pool) THEN
                    WRITE(numout,*) 'reserve pool, forcedrecord', reserve_pool
                 ENDIF
!End JC Debug MOD005
             ENDIF ! tree and grasses/crop

             !! 7.2 Move carbon between the reserve pools
             !  Fill the reserve pools up to their optimal level or until the min/max limits are reached
             !  The original approcah in OCN resulted in instabilities and sometimes oscilations. For
             !  this reason a more simple and straightforward transfer between the pools has been 
             !  implemented.

             !! 7.2.1 Burn excess reserves
             !  The actual reserve and/or labile pool exceed the required pools (as calculated in 7.1)
             !  The excessive reserve pools respires C which needs to be accounted for in the 
             !  growth respiration. Because of this line of code, npp cannot be calculated at the
             !  start of the of this subroutine (i.e. between section 4 and 5). Last, correct the
             !  labile and reserve pool for this respiration flux. Note that instead of respiring
             !  this carbon it could be used for leaching, feeding mycorrhizae, producing DOCs,
             !  producing VOCs or any other component of NPP that does not end up in the biomass.
             IF ( (biomass(ipts,j,icarbres,icarbon) .GE. reserve_pool) .AND. &
                  (biomass(ipts,j,ilabile,icarbon) .GE. labile_pool) ) THEN
                          
                ! reserves are full
                IF ( biomass(ipts,j,icarbres,icarbon) .GT. reserve_pool ) THEN

                   excess = biomass(ipts,j,icarbres,icarbon) - reserve_pool

                   resp_growth(ipts,j) = resp_growth(ipts,j) + 0.1 * excess
                   biomass(ipts,j,icarbres,icarbon) = biomass(ipts,j,icarbres,icarbon) - &
                        0.1 * excess

                ENDIF

                ! DSG: labile is full, but we don't burn from labile, as we made sure
                ! surplus labile goes to reserve fast


             !! 7.2.2 Enough reserves, not enough labile
             ELSEIF ( (biomass(ipts,j,icarbres,icarbon) .GE. reserve_pool) .AND. &
                  (biomass(ipts,j,ilabile,icarbon) .LT. labile_pool) ) THEN

                ! We need to move carbon from the reserve pool to the labile pool. We will move
                ! this gradually. Calculate the maximum flow of the carbon reserve pool to labile 
                IF (is_tree(j)) THEN

                   ! During the dormant season for evergreens or the growing season for both functional leaf
                   ! traits. 
                   IF ( (biomass(ipts,j,ileaf,icarbon) .GT. min_stomate .AND. & ! if there are leaves ...
                        f_alloc(ipts,j,isapabove) .LE. min_stomate) .OR.      & ! ... and we dont grow anything,

!DSGlabfac: the following trigger will be always YES; 
!           this shouldn't be an issue if we make sure there is no reserve
!           coming into labile (which might trigger excess respiration); thus
!           check DSGtriggerXYR938
!DSGlabfac                
                        (lab_fac(ipts,j) .GT. 0.1) ) THEN                       ! ... or in growing season
                      
                      ! Don't move more than an arbitrary 5% from the carbohydrate reserve pool
                      use_max = biomass(ipts,j,icarbres,icarbon) * 0.05
                   ! Dormant season 
                   ELSE 

                      ! Don't move any carbon between the reserve and the labile pool
                      use_max = zero

                   ENDIF ! Growing or dormant season
! JC grass move reserve to labile
                ! Grasses
                ELSE 

                   ! During the growing season
                   IF (biomass(ipts,j,ileaf,icarbon).GT.min_stomate) THEN 

                      ! Don't move more than an arbitrary 5% from the carbohydrate reserve pool
!DSG: increase this 5% fraction? It would take about 1.5 months to get reserve to 1/10 of its intial size with 5%
!     this could be a bit long for grass 
                      use_max = biomass(ipts,j,icarbres,icarbon) * 0.05

                   ! Dormant season
                   ELSE 

                      ! Don't move any carbon between the reserve and the labile pool
                      use_max = zero

                   ENDIF ! Growing or dormant season

                ENDIF ! Trees or grasses

                ! Calculate the required flow of the carbon reserve pool to labile  
                ! Propose to use what can be moved from the reserve pool (::use_max) or 
                ! the amount required to fill the pool.
                use_res = MAX(zero, MIN(use_max, labile_pool-biomass(ipts,j,ilabile,icarbon)))
                   
                ! Update labile pool and reserve
                bm_alloc(ipts,j,icarbres,icarbon) =  use_res
                biomass(ipts,j,ilabile,icarbon) = biomass(ipts,j,ilabile,icarbon) + &
                     bm_alloc(ipts,j,icarbres,icarbon)
                biomass(ipts,j,icarbres,icarbon) = biomass(ipts,j,icarbres,icarbon) - &
                     bm_alloc(ipts,j,icarbres,icarbon)
                
             ! 7.2.3 Enough labile, not enough reserves
             ELSEIF ( (biomass(ipts,j,icarbres,icarbon) .LT. reserve_pool) .AND. &
                  (biomass(ipts,j,ilabile,icarbon) .GE. labile_pool) ) THEN
                
                ! The labile carbon is more mobile than the reserve pool but
                ! it is also more important in the allocation scheme because 
                ! it generates growth respiration and it is used to calculate
                ! bm_alloc. Therefore, the mobility of the labile pool was
                ! restricted to an arbitrary 15% 
                !DSG: one needs to make sure to get ridd of labile or (case1) excess
                !     respiration of labile pool can be huge messing up NP if
                !     activate excess respiration from labile pool / (case2) if not
                !     activate the labile pools gets way too big under nutrient
                !     stress
                
                ! Propose to use what can be moved from the reserve labile pool(::use_max) or 
                ! the amount required to fill the pool.
!DSG: we need to get ridd of carbon from the labile pool if it is above its
!optimal size:
!     a) because labile carbon is osmotic active thus potentially
!     harmful -> plants store it away immediately
!     b) we must avoid that the labile P pool keeps growing under nutrient
!     stress. if reserve pools hit their limits, we have huge excess respiration
!    -> thus:
                use_lab =MAX(zero,0.25*(biomass(ipts,j,ilabile,icarbon)-labile_pool))
                
                ! Update labile pool and reserve
                bm_alloc(ipts,j,icarbres,icarbon) =  use_lab
                biomass(ipts,j,ilabile,icarbon) = biomass(ipts,j,ilabile,icarbon) - &
                     bm_alloc(ipts,j,icarbres,icarbon)
                biomass(ipts,j,icarbres,icarbon) = biomass(ipts,j,icarbres,icarbon) + &
                     bm_alloc(ipts,j,icarbres,icarbon)
             

             !! 7.2.3 We don't have enough carbon in both reserve pools. We will 
             !  redistribute what we have to minimise the tension between available and 
             !  required
             ELSEIF ( (biomass(ipts,j,icarbres,icarbon) .LT. reserve_pool) .AND. &
                  (biomass(ipts,j,ilabile,icarbon) .LT. labile_pool) ) THEN 
             
              !DSG: not needed


             !!  7.4 Unexpected condition
             ELSE
                
                IF(ld_warn .OR. ld_alloc) THEN
                   WRITE(numout,*) 'An unexpected condition occured for the reserve pools'
                   WRITE(numout,*) 'required: reserve and labile pool, ',reserve_pool, labile_pool
                   WRITE(numout,*) 'available: reserve and labile, ', biomass(ipts,j,icarbres,icarbon), &
                        biomass(ipts,j,ilabile,icarbon)
                   CALL ipslerr_p (3,'growth_fun_all',&
                        'An unexpected condition occured for the reserve pools','','')
                ENDIF

             ENDIF
               
          ELSEIF ( veget_max(ipts,j) .GT. min_stomate .AND. &
               rue_longterm(ipts,j) .EQ. un) THEN

             ! There hasn't been any photosynthesis yet. This happens when a new vegetation
             ! is prescribed and the longterm phenology variables are not initialized yet. 
             ! These conditions happen when the model is started from scratch (no restart files).
             ! Because the plants are very small, they contain little reserves. We increased the 
             ! amount of reserves by a factor ::tune_r_in_sapling where r stands for reserves.
             ! However, this amount gets simply respired before it is needed because the 
             ! reserve_pool is calculated as a function of the sapwood biomass which is very
             ! low because the plants are really small. Here we skip recalculating the 
             ! reserve_pool until the day we start using it.

          ELSE

             ! No reason to be here
             WRITE(numout,*) 'Error: unexpected condition for the reserve pools, pft, ',j
             WRITE(numout,*) 'veget_max, rue_longterm, ', veget_max(ipts,j), rue_longterm(ipts,j)

          ENDIF ! rue_longterm

          !! 7.8 Calculate NPP
          !  Calculate the NPP @tex $(gC.m^{-2}dt^{-1})$ @endtex as the difference between GPP and the two 
          !  components of autotrophic respiration (maintenance and growth respiration). GPP, R_maint and R_growth 
          !  are prognostic variables, NPP is calculated as the residuals and is thus a diagnostic variable. 
          !  Note that NPP is not used in the allocation scheme, instead bm_alloc_tot is allocated. The 
          !  physiological difference between both is that bm_alloc_tot does no longer contain the reserves and
          !  labile pools and is only the carbon that needs to go into the biomass pools. NPP contains the reserves
          !  and labile carbon. Note that GPP is in gC m-2 s-1 whereas the respiration terms were calculated in
          !  gC m-2 dt-1

          npp(ipts,j) = gpp_daily(ipts,j) - resp_growth(ipts,j)/dt - resp_maint(ipts,j)/dt


         ! overwrite the allocation fraction to labile and carboresrves with
         ! pool changes to account for the reshuffeling of C between both pools
          bm_alloc(ipts,j,icarbres,icarbon) = biomass(ipts,j,icarbres,icarbon) - bm_res_lab(ipts,j,1)
          bm_alloc(ipts,j,ilabile,icarbon)  = biomass(ipts,j,ilabile,icarbon)  - bm_res_lab(ipts,j,2)

          ! as this could introdcuce mass conservation issue we check here that
          ! bm_alloc is consistent with npp
!DSGdebug: because of the min_stomate threshold being used all around the code
!here, the accuracy is a bit less than min_stomate; here we use +25%
!          IF ((ABS(npp(ipts,j) - SUM(bm_alloc(ipts,j,:,icarbon))).GT.min_stomate)) THEN
          IF ((ABS(npp(ipts,j) - SUM(bm_alloc(ipts,j,:,icarbon))).GT. (1.25*min_stomate))) THEN
!DSG
            WRITE(6,*) "npp is not equal the allocated carbon"
            WRITE(6,*) "grid=",ipts
            WRITE(6,*) "pft=",j
            WRITE(6,*) "npp(ipts,j)", npp(ipts,j)
            WRITE(6,*) "SUM(bm_alloc(ipts,j,:,icarbon))",SUM(bm_alloc(ipts,j,:,icarbon))
            WRITE(6,*) "ABS delta",ABS(npp(ipts,j) - SUM(bm_alloc(ipts,j,:,icarbon)))
            CALL ipslerr_p ( 3, 'npp not equal allocated carbon ','model STOP', &
                               &         '','')
          ENDIF

          ! 5.3  use or fill reserve pools depending on relative size of the labile and reserve N pool

          IF(veget_max(ipts,j) .GT. min_stomate) THEN
             ! nitrogen
             costf_N = f_alloc(ipts,j,ileaf) + fcn_wood_act(ipts,j) * (f_alloc(ipts,j,isapabove)+f_alloc(ipts,j,isapbelow)) + &
                  fcn_root(j) * ( f_alloc(ipts,j,iroot) + f_alloc(ipts,j,ifruit) )

             ! phosphorus
             costf_P = f_alloc(ipts,j,ileaf)                                             & 
                        + (fcn_wood_act(ipts,j)*fnp_wood_act(ipts,j))                    & 
                           * (f_alloc(ipts,j,isapabove)+f_alloc(ipts,j,isapbelow))      & 
                        + (fcn_root(j)*fnp_root(j))                                      &
                            * ( f_alloc(ipts,j,iroot) + f_alloc(ipts,j,ifruit) ) 

             IF (costf_N.EQ.zero) costf_N=un
             IF (costf_P.EQ.zero) costf_P=un

             ! 5.3.1 Nitrogen
             
             ! target pool size of labile nitrogen is derived from the size of labile carbon using the
             ! current leaf CN; DSG: it's better to use a target labile N to be able to allocate leaves of a minimal CN 

             labile_pool = biomass(ipts,j,ilabile,icarbon)/cn_leaf_min(j) * costf_N

             ! excess or deficit of nitrogen in the labile pool
             ! use_lab is negative if labile pool is less than optimal
             use_lab=biomass(ipts,j,ilabile,initrogen)-labile_pool

             ! Target reserve_pool size
             !   is assumed to be the nitrogen needed to alloacte 12% of the total sapwood biomass according to the CN of new growth;
             !   DSG: this is dangerous, plants should store reserve N according
             !   to actual reserve C, not the optimal reserve C stock in case
             !   they are C limited. DSGlab03a fixes that
             
             !   in case plants are N limited they should not store any reserve,
             !   but should give priority to the labile pool. To ensure this, we
             !   use leaf_cn_min for the optimal labile pool but cn_leaf for
             !   optimal reserve (DSGlab03b)
 



                reserve_pool = biomass(ipts,j,icarbres,icarbon)/cn_leaf(ipts,j) * & 
                                (1.+fcn_root(j)/ltor(ipts,j))/(1.+root_reserve(j)/ltor(ipts,j))
             
             ! Calculate the maximum flow of the nitrogen reserve pool to labile (use_res)
             ! the fuller the reserve is the more we give, the fuller the labile
             ! is the less it wants 

             IF(reserve_pool.GT.min_stomate)THEN ! if we have reserves... 
!JC Debug MOD013
! According to the explanation following the calculation of reserve_scal
! the correct function should be
                reserve_scal = MAX(MIN(biomass(ipts,j,icarbres,initrogen) / reserve_pool,1.), 0.1)
!                reserve_scal = MIN(MAX(1. - biomass(ipts,j,icarbres,initrogen) / reserve_pool,1.),0.1) ! calculate a scaling factor 
                                                                                                       ! which is 1 if reserve > optimal
                                                                                                       ! and between <1 and 0.1 depending on the deviation of the reserve from optimal; therbey the reserve pool can be depleted my maximal 10% if very low

!End JC Debug
                IF(labile_pool.GT.min_stomate)THEN 
                   use_res=biomass(ipts,j,icarbres,initrogen) * reserve_scal * &
                        MAX(1.-biomass(ipts,j,ilabile,initrogen) / labile_pool,0.0) ! scale the flow with the deviation from optimal size of labile N
                ELSE
                   use_res=zero
                ENDIF
             ELSE
                use_res=biomass(ipts,j,icarbres,initrogen) ! just use it up if it the reserves are marginal
             ENDIF

             IF(use_lab.GT.0.0)THEN ! if the labile pool is larger than optimal ...
                ! ... we assign N to fill the reserve to optimal size if needed
                use_max=MAX(-use_lab, biomass(ipts,j,icarbres,initrogen)-reserve_pool)

             ! if the labile pool is smaller than optimal ...
             ELSEIF(use_res.GT.-use_lab)THEN  ! ... (a) and the reserve pool can give as much as needed to fill labile N to opitmal
                use_max=-use_lab              !        ... then give it all to the labile pool
             ELSE                             ! ... (b) and the reserve pool cannot give as much as needed to fill the labile N to optimal size
                use_max=use_res               !        ... then give as much as we as the reserve pool allows
             ENDIF

             biomass(ipts,j,icarbres,initrogen) = biomass(ipts,j,icarbres,initrogen)-use_max
             biomass(ipts,j,ilabile,initrogen)  = biomass(ipts,j,ilabile,initrogen)+use_max
             bm_alloc(ipts,j,icarbres,initrogen)= bm_alloc(ipts,j,icarbres,initrogen)-use_max

             !IF((ipts == test_grid).AND.(j == test_pft)) THEN
             !  WRITE(numout,*) 'biomass(test_grid,test_pft,:,initrogen)    : ',  biomass(test_grid,test_pft,:,initrogen)
             !ENDIF

             ! we do the same for P as for N; so please check the comments above
             ! to understand what is done here
             ! 5.3.2 Phosphorus:

             labile_pool = biomass(ipts,j,ilabile,icarbon)/  &
                                       (cn_leaf_min(j) * &
                                        np_leaf(ipts,j) * costf_P)

             ! excess or deficit of phosphorus in the labile pool
             use_lab=biomass(ipts,j,ilabile,iphosphorus) - labile_pool

             reserve_pool = biomass(ipts,j,icarbres,icarbon)/cn_leaf(ipts,j) * & 
                               (un+fcn_root(j)/ltor(ipts,j))/(un+(root_reserve(j))/ltor(ipts,j)) / &
                                np_leaf(ipts,j) *                                & 
                              (un+fnp_root(j)/ltor(ipts,j))/(un+(root_reserve(j))/ltor(ipts,j))

!JC Debug
                 IF (j .EQ. test_pft .AND. ld_alloc .AND. ipts == test_grid) THEN
                    WRITE(numout,*) 'costf_N,costf_P',costf_N,costf_P
                    WRITE(numout,*) 'iphosphorus labile_pool,reserve_pool',labile_pool,reserve_pool
                    WRITE(numout,*) 'iphosphorus ilabile icarbres',biomass(ipts,j,ilabile,iphosphorus),biomass(ipts,j,icarbres,iphosphorus)
                    WRITE(numout,*) 'np_leaf_min,np_leaf',np_leaf_min(j),np_leaf(ipts,j)
                    WRITE(numout,*) 'iphosphorus use_lab',use_lab
                 ENDIF

!End JC Debug


             IF(reserve_pool.GT.min_stomate)THEN
!JC Debug MOD013
! According to the explanation following the calculation of reserve_scal
! the correct function should be
                reserve_scal = MAX(MIN(biomass(ipts,j,icarbres,iphosphorus) / reserve_pool,1.), 0.1)
!                reserve_scal = MIN(MAX(1. - biomass(ipts,j,icarbres,iphosphorus) / reserve_pool,1.),0.1)
!End JC Debug
                IF(labile_pool.GT.min_stomate)THEN
                   use_res = biomass(ipts,j,icarbres,iphosphorus) * reserve_scal * &
                        MAX(un-biomass(ipts,j,ilabile,iphosphorus) / labile_pool,zero)
                ELSE
                   use_res=zero
                ENDIF
             ELSE
                use_res=biomass(ipts,j,icarbres,iphosphorus)
             ENDIF

             IF(use_lab.GT.0.0)THEN
                use_max=MAX(-use_lab,biomass(ipts,j,icarbres,iphosphorus)-reserve_pool)
             ELSEIF(use_res.GT.-use_lab)THEN
                use_max=-use_lab
             ELSE
                use_max=use_res
             ENDIF

             biomass(ipts,j,icarbres,iphosphorus) = biomass(ipts,j,icarbres,iphosphorus) - use_max
             biomass(ipts,j,ilabile,iphosphorus)  = biomass(ipts,j,ilabile,iphosphorus)   + use_max
             bm_alloc(ipts,j,icarbres,iphosphorus)= bm_alloc(ipts,j,icarbres,iphosphorus) - use_max

             IF(printlev>=4)THEN
                 IF((ipts==test_grid).AND.(j==test_pft)) THEN
                    WRITE(numout,*) 'CHECK values in growth_fun_all'
                    WRITE(numout,*) 'biomass(test_grid,test_pft,icarb,iphosphorus): ',  biomass(test_grid,test_pft,icarbres,iphosphorus)
                    WRITE(numout,*) 'biomass(test_grid,test_pft,icarb,initrogen): ',  biomass(test_grid,test_pft,icarbres,initrogen)
                    WRITE(numout,*) 'biomass(test_grid,test_pft,ilabile,iphosphorus): ',  biomass(test_grid,test_pft,ilabile,iphosphorus)
                    WRITE(numout,*) 'use_max ', use_max
                    WRITE(numout,*) 'use_res ', use_res
                    WRITE(numout,*) 'use_lab ', use_lab
                    WRITE(numout,*) 'reserve_pool ', reserve_pool
                    WRITE(numout,*) 'labile_pool ', labile_pool
                    WRITE(numout,*) 'costf_P', costf_P
                 ENDIF
              ENDIF

!DSGlabXY
                 ! We need to kept the reserve nitrogen & phosphorus pool filled
                 ! according to the filling status of carbon reserve to allow
                 ! decidiuous PFTs to grow when we initialize a impose_cn=n with the
                 ! restarts from a impose_cn=y simulation; if we dont ensure
                 ! that there is enough reserve N the PFT will get extinct, 
                 ! so ...
!JC Debug MOD015
! The grass need the same reserve pool to make the grasses growth
! when restart from impose_cn=y
! Though it may be less important as that for deciduous forest
!                IF(is_tree(j) ) THEN
                 IF(natural(j)) THEN
                   IF (pheno_type(j).NE.1) THEN 
                     IF (impose_cn) THEN
                        ! take what is needed to keep nitrogen reserves optimal
                        N_support(ipts,j) = N_support(ipts,j)+ &
                                           MAX((biomass(ipts,j,icarbres,icarbon)/cn_leaf(ipts,j) * &  
                                           (1.+fcn_root(j)/ltor(ipts,j))/(1.+root_reserve(j)/ltor(ipts,j))) &
                                            - biomass(ipts,j,icarbres,initrogen),zero)
                        ! fill the pool
                        biomass(ipts,j,icarbres,initrogen)=MAX(biomass(ipts,j,icarbres,initrogen), &
                                                           biomass(ipts,j,icarbres,icarbon)/cn_leaf(ipts,j) * &  
                                                           (1.+fcn_root(j)/ltor(ipts,j))/(1.+root_reserve(j)/ltor(ipts,j)))

                     ELSEIF (impose_np) THEN
                        P_support(ipts,j) = P_support(ipts,j)+ &
                                            MAX((biomass(ipts,j,icarbres,icarbon)/cn_leaf(ipts,j) * & 
                                                (un+fcn_root(j)/ltor(ipts,j))/(un+(root_reserve(j))/ltor(ipts,j)) / &
                                                 np_leaf(ipts,j) *                                & 
                                                (un+fnp_root(j)/ltor(ipts,j))/(un+(root_reserve(j))/ltor(ipts,j))) &
                                                - biomass(ipts,j,icarbres,iphosphorus),zero)

                        ! fill the pool
                        biomass(ipts,j,icarbres,iphosphorus)=MAX(biomass(ipts,j,icarbres,iphosphorus), &
                                                             biomass(ipts,j,icarbres,icarbon)/cn_leaf(ipts,j) * & 
                                                             (un+fcn_root(j)/ltor(ipts,j))/(un+(root_reserve(j))/ltor(ipts,j)) / &
                                                              np_leaf(ipts,j) *                                & 
                                                             (un+fnp_root(j)/ltor(ipts,j))/(un+(root_reserve(j))/ltor(ipts,j)))
                     ENDIF 
                   ENDIF
                 ENDIF

!DSGlabXY


          ENDIF

!DSGremove
       IF (printlev>=4 .AND. ipts==test_grid .AND. j==test_pft) THEN
         WRITE(numout,*) 'location in code: '
         WRITE(numout,*) 'stomate_growth_fun: before shifting'
         WRITE(numout,*) 'biomass(test_grid,test_pft,:,icarbon)       : ',  biomass(test_grid   ,test_pft,:,icarbon)
       ENDIF
!DSGremove


          !! 7.9 Distribute stand level ilabile and icarbres at the tree level
          !  The labile and carbres pools are calculated at the stand level but are then redistributed at the
          !  tree level. This has the advantage that biomass and circ_class_biomass have the same dimensions
          !  for nparts which comes in handy when phenology and mortality are calculated.
                
          IF (is_tree(j)) THEN
 
             ! Initialize to enable a loop over nparts
             circ_class_biomass(ipts,j,:,ilabile,:) = zero
             circ_class_biomass(ipts,j,:,icarbres,:) = zero
 
             ! Distribute labile and reserve pools over the circumference classes 
             DO m = 1,nelements

                ! Total biomass across parts and circumference classes
                temp_total_biomass = zero

                DO l = 1,ncirc 
                   
                   DO k = 1,nparts
                      
                      temp_total_biomass = temp_total_biomass + circ_class_biomass(ipts,j,l,k,icarbon) * circ_class_n(ipts,j,l)
                                
                   ENDDO

                ENDDO
              
                ! Total biomass across parts but for a specific circumference class
                DO l = 1,ncirc

                   temp_class_biomass = zero
  
                   DO k = 1,nparts

                      temp_class_biomass = temp_class_biomass + circ_class_biomass(ipts,j,l,k,icarbon) * circ_class_n(ipts,j,l)
               
                   ENDDO

                   IF(  temp_total_biomass .NE. zero) THEN

                      ! Share of this circumference class to the total biomass
                      temp_share = temp_class_biomass / temp_total_biomass
 
                      ! Allocation of ilabile at the tree level (gC tree-1)
                      circ_class_biomass(ipts,j,l,ilabile,m) = temp_share * &
                           biomass(ipts,j,ilabile,m) / circ_class_n(ipts,j,l)
                      
                      ! Allocation of icarbres at the tree level (gC tree-1)
                      circ_class_biomass(ipts,j,l,icarbres,m) = temp_share * &
                           biomass(ipts,j,icarbres,m) / circ_class_n(ipts,j,l)

                      ! we need to synchronize all the circ_class pools if it is
                      ! carbon or nutrients:
                      circ_class_biomass(ipts,j,l,ileaf,m) = temp_share * &
                        biomass(ipts,j,ileaf,m) / circ_class_n(ipts,j,l)
                      circ_class_biomass(ipts,j,l,isapabove,m) = temp_share * &
                        biomass(ipts,j,isapabove,m) / circ_class_n(ipts,j,l)
                      circ_class_biomass(ipts,j,l,isapbelow,m) = temp_share * &
                        biomass(ipts,j,isapbelow,m) / circ_class_n(ipts,j,l)
                      circ_class_biomass(ipts,j,l,iheartabove,m) = temp_share * &
                        biomass(ipts,j,iheartabove,m) / circ_class_n(ipts,j,l)
                      circ_class_biomass(ipts,j,l,iheartbelow,m) = temp_share * &
                        biomass(ipts,j,iheartbelow,m) / circ_class_n(ipts,j,l)
                      circ_class_biomass(ipts,j,l,iroot,m) = temp_share * &
                        biomass(ipts,j,iroot,m) / circ_class_n(ipts,j,l)
                      circ_class_biomass(ipts,j,l,ifruit,m) = temp_share * &
                        biomass(ipts,j,ifruit,m) / circ_class_n(ipts,j,l)
                   ELSE

                      circ_class_biomass(ipts,j,l,ilabile,m)  = zero
                      circ_class_biomass(ipts,j,l,icarbres,m) = zero

                      IF((m .EQ. initrogen).OR.(m .EQ. iphosphorus)) THEN
                          !DSG: mayeb this IF statement is wrong liek the one
                          !before?!    
                         circ_class_biomass(ipts,j,l,ileaf,m)       = zero
                         circ_class_biomass(ipts,j,l,isapabove,m)   = zero
                         circ_class_biomass(ipts,j,l,isapbelow,m)   = zero
                         circ_class_biomass(ipts,j,l,iheartabove,m) = zero
                         circ_class_biomass(ipts,j,l,iheartbelow,m) = zero
                         circ_class_biomass(ipts,j,l,iroot,m)       = zero
                         circ_class_biomass(ipts,j,l,ifruit,m)      = zero
                      ENDIF

                   ENDIF

                ENDDO ! ncirc

             ENDDO  ! nelements

          ! Grasses and crops
          ELSE

             DO m = 1,nelements

                ! synchronize biomass and circ_class_biomass
                IF (ind(ipts,j) .GT. zero) THEN

                   circ_class_biomass(ipts,j,1,:,m) = biomass(ipts,j,:,m) / ind(ipts,j)

                ELSE

                   circ_class_biomass(ipts,j,1,:,m) = zero

                ENDIF

             ENDDO  
          
          ENDIF ! is_tree
  
       ENDDO ! pnts

    ENDDO ! PFTs

    IF(dsg_debug) THEN
       WRITE(numout,*) 'AFTER distributing stand level ilabile'
       IF (printlev>=4) THEN
         WRITE(numout,*) 'location in code: '
         WRITE(numout,*) 'stomate_growth_fun: after shifting'
         WRITE(numout,*) 'biomass(test_grid,test_pft,:,icarbon)       : ',  biomass(test_grid   ,test_pft,:,icarbon)
       ENDIF
       CALL check_mass(npts,biomass(:,:,:,:),'stomte_growth_fun: after shift standlevel')
    ENDIF

 !! 8. Check mass balance closure
    
    ! Calculate pools at the end of the routine
    pool_end = zero
    DO ipar = 1,nparts
       DO iele = 1,nelements
          pool_end(:,:,iele) = pool_end(:,:,iele) + &
               (biomass(:,:,ipar,iele) * veget_max(:,:))
       ENDDO
    ENDDO

    ! Calculate components of the mass balance
    check_intern(:,:,iatm2land,icarbon) = gpp_daily(:,:) * dt * veget_max(:,:)
    check_intern(:,:,iland2atm,icarbon) = -un * (resp_maint(:,:) + resp_growth(:,:)) * &
         veget_max(:,:)
    check_intern(:,:,ilat2out,icarbon) = -un * zero
    check_intern(:,:,ilat2in,icarbon) = un * zero
    check_intern(:,:,ipoolchange,icarbon) = -un * (pool_end(:,:,icarbon) - &
         pool_start(:,:,icarbon))
    closure_intern = zero
    DO imbc = 1,nmbcomp
       closure_intern(:,:,icarbon) = closure_intern(:,:,icarbon) + &
            check_intern(:,:,imbc,icarbon)
    ENDDO

    ! Write conclusion
    DO ipts=1,npts
       DO j=1,nvm
          IF(ABS(closure_intern(ipts,j,icarbon)) .LE. min_stomate)THEN
             IF (ld_massbal) WRITE(numout,*) 'Mass balance closure in stomate_growth_fun_all.f90'
          ELSE
             WRITE(numout,*) 'Error: mass balance is not closed in stomate_growth_fun_all.f90'
             WRITE(numout,*) '   ipts,j; ', ipts,j
             WRITE(numout,*) '   Difference is, ', closure_intern(ipts,j,icarbon)
             WRITE(numout,*) '   pool_end,pool_start: ', pool_end(ipts,j,icarbon), pool_start(ipts,j,icarbon)
             WRITE(numout,*) '   gpp_daily, veget_max: ', gpp_daily(ipts,j),veget_max(ipts,j)
             WRITE(numout,*) '   resp_maint,resp_growth: ', resp_maint(ipts,j),resp_growth(ipts,j)
             IF(ld_stop)THEN
                CALL ipslerr_p (3,'growth_fun_all', 'Mass balance error.','','')
             ENDIF
          ENDIF
       ENDDO
    ENDDO

!    !+++HACK++++
!    ! JR 080315 temporary hardwire for testing PFTs 4
!    ! comment out this sections for tree growth profiles 
!    IF (ld_fake_height) THEN
!       !Do nothing
!    ELSE 
!       !Go to James' Hardwire 
!       IF (jr_nextstep) THEN
!         ! this was a simple means to hardwire a canopy profile without having to use a spin-up file
!         ! for testing. It is now commented out to avoid compilation errors for the INPUT variables
!         ! ind and circ_class_n
!
!         ! ind(1,4) = 0.6d0     
!         ! biomass(1,4,1,1) = 10919.1492214105     
!         ! biomass(1,4,2,1) = 119092.584535974     
!         ! biomass(1,4,3,1) = 119092.584535974     
!         ! biomass(1,4,4,1) = 23818.5169071949   
!         ! biomass(1,4,5,1) = 23818.5169071949     
!         ! biomass(1,4,6,1) = 3717.22198731141     
!         ! biomass(1,4,7,1) = 0.000000000000000E+000
!         ! biomass(1,4,8,1) = 0.000000000000000E+000
!         ! biomass(1,4,9,1) = 348.380698397579     
!         ! circ_class_n(1,4,1) = 0.570200000000000     
!         ! circ_class_n(1,4,2) = 2.840000000000000E-002
!         ! circ_class_n(1,4,3) = 1.400000000000000E-003
!
!         ! circ_class_biomass(1,4,1,1,1) =  18198.5820356842     
!         ! circ_class_biomass(1,4,1,2,1) =  198487.640893291     
!         ! circ_class_biomass(1,4,1,3,1) =  198487.640893291     
!         ! circ_class_biomass(1,4,1,4,1) =  39697.5281786581     
!         ! circ_class_biomass(1,4,1,5,1) =  39697.5281786581     
!         ! circ_class_biomass(1,4,1,6,1) =  6195.36997885235     
!         ! circ_class_biomass(1,4,1,7,1) = 0.000000000000000E+000
!         ! circ_class_biomass(1,4,1,8,1) =  0.000000000000000E+000
!         ! circ_class_biomass(1,4,1,9,1) =  580.634497329298     
!         ! circ_class_biomass(1,4,2,1,1) =  18198.5820356842     
!         ! circ_class_biomass(1,4,2,2,1) =  198487.640893291     
!         ! circ_class_biomass(1,4,2,3,1) =  198487.640893291     
!         ! circ_class_biomass(1,4,2,4,1) =  39697.5281786581     
!         ! circ_class_biomass(1,4,2,5,1) =  39697.5281786581     
!         ! circ_class_biomass(1,4,2,6,1) =  6195.36997885235     
!         ! circ_class_biomass(1,4,2,7,1) =  0.000000000000000E+000
!         ! circ_class_biomass(1,4,2,8,1) =   0.000000000000000E+000
!         ! circ_class_biomass(1,4,2,9,1) =  580.634497329298     
!         ! circ_class_biomass(1,4,3,1,1) =  18198.5820356842     
!         ! circ_class_biomass(1,4,3,2,1) =  198487.640893291     
!         ! circ_class_biomass(1,4,3,3,1) =  198487.640893291     
!         ! circ_class_biomass(1,4,3,4,1) =  39697.5281786581     
!         ! circ_class_biomass(1,4,3,5,1) =  39697.5281786581     
!         ! circ_class_biomass(1,4,3,6,1) =  6195.36997885235     
!         ! circ_class_biomass(1,4,3,7,1) =  0.000000000000000E+000
!         ! circ_class_biomass(1,4,3,8,1) =  0.000000000000000E+000
!         ! circ_class_biomass(1,4,3,9,1) =  580.634497329298     
!       END IF ! (jr_nextstep)
!    END IF !(ld_fake_height)
    !----------


 !! 9. Update leaf age

    !  Leaf age is needed to calculate the turnover and vmax in the stomate_turnover.f90 and stomate_vmax.f90 routines. 
    !  Leaf biomass is distributed according to its age into several "age classes" with age class=1 representing the
    !  youngest class, and consisting of the most newly allocated leaf biomass. 
    
    !! 9.1 Update quantity and age of the leaf biomass in the youngest class
    !  The new amount of leaf biomass in the youngest age class (leaf_mass_young) is the sum of :
    !  - the leaf biomass that was already in the youngest age class (leaf_frac(:,j,1) * lm_old(:,j)) with the 
    !  leaf age given in leaf_age(:,j,1) 
    !  - and the new biomass allocated to leaves (bm_alloc(:,j,ileaf,icarbon)) with a leaf age of zero.
!JC separate tissue age
! add nparts young mass
!    leaf_mass_young(:,:) = leaf_frac(:,:,1) * lm_old(:,:) + bm_alloc(:,:,ileaf,icarbon)
    leaf_mass_young(:,:,:) = leaf_frac(:,:,:,1) * lm_old(:,:,:) + bm_alloc(:,:,:,icarbon)
!End JC separate tissue age
    ! The age of the updated youngest age class is the average of the ages of its 2 components: bm_alloc(leaf) of age
    ! '0', and leaf_frac*lm_old(=leaf_mass_young-bm_alloc) of age 'leaf_age(:,:,1)' 
    DO ipts=1,npts

       DO j=1,nvm

          ! IF(veget_max(ipts,j) == zero)THEN
          !     ! this vegetation type is not present, so no reason to do the 
          !     ! calculation
          !     CYCLE
          ! ENDIF
!JC separate tissue age
!          IF( (bm_alloc(ipts,j,ileaf,icarbon) .GT. min_stomate ) .AND. &
!               ( leaf_mass_young(ipts,j) .GT. min_stomate ) )THEN
!             
!
!             leaf_age(ipts,j,1) = MAX ( zero, leaf_age(ipts,j,1) * &
!                  ( leaf_mass_young(ipts,j) - bm_alloc(ipts,j,ileaf,icarbon) ) / &
!                  & leaf_mass_young(ipts,j) )
!             
!          ENDIF
! calculate for all tissues
          DO ipar = 1,nparts
            IF( (bm_alloc(ipts,j,ipar,icarbon) .GT. min_stomate ) .AND. &
                 ( leaf_mass_young(ipts,j,ipar) .GT. min_stomate ) )THEN
  
  
               leaf_age(ipts,j,ipar,1) = MAX ( zero, leaf_age(ipts,j,ipar,1) * &
                    ( leaf_mass_young(ipts,j,ipar) - bm_alloc(ipts,j,ipar,icarbon) ) / &
                    & leaf_mass_young(ipts,j,ipar) )
  
            ENDIF
          ENDDO
!End JC separate tissue age          

          !+++TEMP+++
!!$          IF(j == test_pft .AND. ipts == test_grid)THEN
!!$             WRITE(numout,*) 'VCMAX: leaf_age growth: ',leaf_age(ipts,j,1),&
!!$                  bm_alloc(ipts,j,ileaf,icarbon),leaf_mass_young(ipts,j),&
!!$                  biomass(ipts,j,ileaf,icarbon)
!!$          ENDIF
          !++++++++++

       ENDDO

    ENDDO
          
    !! 8.2 Decrease reduction of photosynthesis 
    !  Decrease reduction of photosynthesis from new (undamaged) foliage
!!$    WHERE(biomass(:,:,ileaf,icarbon).GT.min_stomate)
!!$    
!!$      t_photo_stress(:,:) =  (t_photo_stress(:,:) * lm_old(:,:) + & 
!!$            bm_alloc(:,:,ileaf,icarbon))/biomass(:,:,ileaf,icarbon)
!!$   
!!$    ENDWHERE


    !! 9.3 Update leaf age
    !  Update fractions of leaf biomass in each age class (fraction in youngest class increases)

    !! 9.3.1 Update age of youngest leaves
    !  For age class 1 (youngest class), because we have added biomass to the youngest class, we need to update
    !  the fraction of total leaf biomass that belongs to the youngest age class : updated mass in class divided
    !  by new total leaf mass
!JC separate tissue age
!    WHERE ( biomass(:,:,ileaf,icarbon) .GT. min_stomate )
!
!          leaf_frac(:,:,1) = leaf_mass_young(:,:) / biomass(:,:,ileaf,icarbon)
!
!    ENDWHERE
! update fraction of youngest tissues
    WHERE ( biomass(:,:,:,icarbon) .GT. min_stomate )

          leaf_frac(:,:,:,1) = leaf_mass_young(:,:,:) / biomass(:,:,:,icarbon)

    ENDWHERE
!End JC separate tissue age

    !! 9.3.2 Update age of other age classes
    !  Because the total leaf biomass has changed, we need to update the fraction of leaves in each age class:
    !  mass in leaf age class (from previous fraction of leaves in this class and previous total leaf biomass) 
    !  divided by new total mass
!JC separate tissue age
!    DO m = 2, nleafages ! Loop over # leaf age classes
!
!       WHERE ( biomass(:,:,ileaf,icarbon) .GT. min_stomate )
!
!          leaf_frac(:,:,m) = leaf_frac(:,:,m) * lm_old(:,:) / biomass(:,:,ileaf,icarbon)
!
!       ENDWHERE
!
!    ENDDO       ! Loop over # leaf age classes
! update fraction of other age class tissues
    DO m = 2, nleafages ! Loop over # leaf age classes

       WHERE ( biomass(:,:,:,icarbon) .GT. min_stomate )

          leaf_frac(:,:,:,m) = leaf_frac(:,:,:,m) * lm_old(:,:,:) / biomass(:,:,:,icarbon)

       ENDWHERE

    ENDDO       ! Loop over # leaf age classes
!End JC separate tissue age

!gmjc varied sla for managed grassland
    leaf_age_w = 0.0
    DO j = 2,nvm
      IF (is_grassland_manag(j)) THEN
         sla_max_Nfert(:,j)=sla_max(j)
         sla_min_Nfert(:,j)=sla_min(j)

      ELSE
         sla_max_Nfert(:,j)=sla_max(j)
         sla_min_Nfert(:,j)=sla_min(j)
      ENDIF

      WHERE ( ( bm_alloc(:,j,ileaf,icarbon) .GT. 0.0 ) .AND. &
!JC separate tissue age
!             ( leaf_mass_young(:,j) .GT. 0.0 ) )
             ( leaf_mass_young(:,j,ileaf) .GT. 0.0 ) )

       sla_age1(:,j) = (sla_age1(:,j) * &
!JC separate tissue age
!                       (leaf_mass_young(:,j)-bm_alloc(:,j,ileaf,icarbon)) + &
!                       sla_max(j) * bm_alloc(:,j,ileaf,icarbon)) / leaf_mass_young(:,j)
                       (leaf_mass_young(:,j,ileaf)-bm_alloc(:,j,ileaf,icarbon)) + &
                       sla_max(j) * bm_alloc(:,j,ileaf,icarbon)) / leaf_mass_young(:,j,ileaf)
       sla_age2(:,j) = sla_max(j)*0.9
       sla_age3(:,j) = sla_max(j)*0.85
       sla_age4(:,j) = sla_max(j)*0.8
!JC separate tissue age
!       sla_calc(:,j) = sla_age1(:,j) * leaf_frac(:,j,1) + &
!                       sla_age2(:,j) * leaf_frac(:,j,2) +  &
!                       sla_age3(:,j) * leaf_frac(:,j,3) + &
!                       sla_age4(:,j) * leaf_frac(:,j,4)
       sla_calc(:,j) = sla_age1(:,j) * leaf_frac(:,j,ileaf,1) + &
                       sla_age2(:,j) * leaf_frac(:,j,ileaf,2) +  &
                       sla_age3(:,j) * leaf_frac(:,j,ileaf,3) + &
                       sla_age4(:,j) * leaf_frac(:,j,ileaf,4)
      ENDWHERE

      leaf_age_w(:,j) = 0.0
      DO m = 1, nleafages
!JC separate tissue age
!        leaf_age_w(:,j) = leaf_age_w(:,j)+ leaf_age(:,j,m)*leaf_frac(:,j,m)
        leaf_age_w(:,j) = leaf_age_w(:,j)+ leaf_age(:,j,ileaf,m)*leaf_frac(:,j,ileaf,m)
      END DO

      ! sla_calc can not be greater than sla_max or less than sla_min, and sla
      ! will be at maximum when age< 10
      WHERE (sla_calc(:,j) .GT. sla_max(j))
        sla_calc(:,j) = sla_max(j)
      ELSE WHERE (sla_calc(:,j) .LT. sla_min(j))
        sla_calc(:,j) = sla_min(j)
      ENDWHERE
    END DO
!end gmjc

    IF(dsg_debug) THEN
       !DSG debug: ===================================================
       WRITE(6,*) 'location in code: '
       WRITE(6,*) 'stomate: growth_fun: 10'
    !   WRITE(6,*) 'biomass(test_grid,test_pft,icarbres,icarbon)    : ',  biomass(test_grid,test_pft,icarbres,icarbon)
       CALL check_mass(npts,biomass(:,:,:,:), 'stomate: growth_fun: 10')
    ENDIF


 !! 10. Update whole-plant age 
    
    !! 10.1 PFT age
    !  At every time step, increase age of the biomass that was already present at previous time step. 
    !  Age is expressed in years, and the time step 'dt' in days so age increase is: dt divided by number 
    !  of days in a year.
    WHERE ( PFTpresent(:,:) )

       age(:,:) = age(:,:) + dt/one_year

    ELSEWHERE

       age(:,:) = zero

    ENDWHERE


    !! 10.2 Age of grasses and crops
    !  For grasses and crops, biomass with age 0 has been added to the whole plant with age 'age'. New biomass is the sum of 
    !  the current total biomass in all plant parts (bm_new), bm_new(:) = SUM( biomass(:,j,:), DIM=2 ). The biomass that has 
    !  just been added is the sum of the allocatable biomass of all plant parts (bm_add), its age is zero. bm_add(:) = 
    !  SUM( bm_alloc(:,j,:,icarbon), DIM=2 ). Before allocation, the plant biomass is bm_new-bm_add, its age is "age(:,j)". 
    !  The age of the new biomass is the average of the ages of previous and added biomass.
    !  For trees, age is treated in "establish" if vegetation is dynamic, and in turnover routines if it is static (in this 
    !  case, only the age of the heartwood is accounted for).
    DO j = 2,nvm

       IF ( .NOT. is_tree(j) ) THEN

          bm_new(:) = biomass(:,j,ileaf,icarbon) + biomass(:,j,isapabove,icarbon) + &
               biomass(:,j,iroot,icarbon) + biomass(:,j,ifruit,icarbon)
          bm_add(:) = bm_alloc(:,j,ileaf,icarbon) + bm_alloc(:,j,isapabove,icarbon) + &
               bm_alloc(:,j,iroot,icarbon) + bm_alloc(:,j,ifruit,icarbon)

          WHERE ( ( bm_new(:) .GT. min_stomate ) .AND. ( bm_add(:) .GT. min_stomate ) )
             
             age(:,j) = age(:,j) * ( bm_new(:) - bm_add(:) ) / bm_new(:)
          
          ENDWHERE

       ENDIF ! is .NOT. tree

    ENDDO  ! Loop over #PFTs 

!!$    ! +++HACK+++
!!$    !  10.3  This is only for the model validation 
!!$    !        LAI is imposed with "IMPOSE_LAI" (the maximun value within a year)     
!!$    !        
!!$    !        Reset the Biomass in leaf, labile and carbonate pools
!!$    !        In order to impose the LAI value & profile we need to recalculate the biomass
!!$    !        based on impose lai and default SLA(specific of leaf area index)
!!$    !        A monthy LAI scaling factor LAI_SCALE was also introduced for descrption the dynamic of LAI  
!!$    !        >>>  This will arise a mass balance issue in stomate_lpj routine 
!!$    !        >>>  Make sure you change the ERROR level of "check_biomass_sync" in funtion_library
!!$    !        >>>  to "2" to avoid model stop  
!!$    !        So, the model is suggested not to run over than one year when the flage turned on.
!DSG     istep=istep+1
!DSG     IF (ld_fake_height) THEN
!DSG        WRITE(numout,'(A,I6,I6)') '!=== CALL FUNTIONAL ALLOCATION STEP & RESET BIOMASS ====',test_pft, istep
!DSG 
!DSG !!$       !SET ind number of the individual in a square meter (trees/m2)
!DSG !!$       !CALL getin_p('STAND_DENSITY',ind(test_grid,test_pft))
!DSG !!$       !Recalculate the stand density 
!DSG !!$       !DO icirc=1,ncirc
!DSG !!$       !   circ_class_n(test_grid,test_pft,icirc)=ind(test_grid,test_pft)*circ_class_dist(icirc)
!DSG !!$       !   WRITE(numout,*) 'Stand density(trees/m2):',ind(test_grid,test_pft)
!DSG !!$       !   WRITE(numout,*) 'circ_class_n:', circ_class_n(test_grid,test_pft,icirc)
!DSG !!$       !ENDDO   
!DSG !!$
!DSG        !GET impose LAI & LAI_SCALE
!DSG        CALL getin_p('IMPOSE_LAI',impose_lai) 
!DSG        DO j=1, 13
!DSG           WRITE(temp_text,'(A11,I5.5)') 'LAI_SCALE__',j
!DSG           CALL getin_p(trim(temp_text),lai_scale(j))
!DSG           WRITE(numout,*) trim(temp_text), lai_scale(j)
!DSG        ENDDO
!DSG        !START a simple linear interpolation based on impose lai and lai scale 
!DSG        !Here, we use a simple 30 days for cycling of one month
!DSG        month_id = INT(istep/30.) + 1
!DSG        IF (month_id .GT. 12) THEN
!DSG           ! only for final 5 days set to a constant as impose_lai*lai_scale(13)
!DSG           daily_lai = (impose_lai)*lai_scale(month_id) 
!DSG        ELSE 
!DSG           daily_lai = ( (impose_lai)*lai_scale(month_id) + &
!DSG                ((lai_scale(month_id+1)-lai_scale(month_id))*impose_lai/30) * (MOD(istep,30)+1) )  
!DSG        ENDIF
!DSG        IF (daily_lai .LT. 0.) daily_lai=0.
!DSG        WRITE(numout,*) 'MONTH_ID:',month_id ,'MONTH_DAY:', (MOD(istep,30)+1)
!DSG        WRITE(numout,*) '!=== Dai=ly LAI ====:', daily_lai
!DSG        WRITE(numout,*) 'BIOMASS_IN_LEAF:', daily_lai/sla(test_pft)
!DSG        !Covert the daily lai to biomass in leaf, labile ans carbres pools
!DSG        biomass(test_grid,test_pft,ileaf,icarbon)=   daily_lai/sla(test_pft)
!DSG        !For these two carbon pools we can simply set them as a coinstant value
!DSG        !labile and carbre pools are for sustaining the photothesis during bad weather conditions.  
!DSG        biomass(test_grid,test_pft,ilabile,icarbon)=   1000.
!DSG        biomass(test_grid,test_pft,icarbres,icarbon)=  1000.
!DSG        !Calculate the biomass in each circ_class again.  
!DSG        DO icirc=1,ncirc
!DSG           circ_class_biomass(test_grid,test_pft,icirc,:,icarbon)= &
!DSG                (biomass(test_grid,test_pft,:,icarbon)/float(ncirc))/circ_class_n(test_grid,test_pft,icirc)
!DSG           WRITE(numout,*) 'circ_class_n:',circ_class_n(test_grid,test_pft,icirc)
!DSG           WRITE(numout,*) 'leaf_biomass:',biomass(test_grid,test_pft,ileaf,icarbon)
!DSG           WRITE(numout,*) 'circ_biomass:',circ_class_biomass(test_grid,test_pft,icirc,ileaf,icarbon)
!DSG        ENDDO
!DSG     ENDIF  !ld_fake_height
!DSG !!$    !++++++++ 



 !! 11. Write history files

    !---TEMP---
!!$    DO ipts = 1, npts
!!$       height_out(ipts,1) = zero
!!$       DO j = 2, nvm
!!$          height_out(ipts,j) = SUM(circ_class_height_eff(:))/ncirc
!!$       ENDDO
!!$    ENDDO
    !---------

!!$    !+++++++++ TEMP ++++++++++
!!$    ! Just for testing.  Set the labile and reserve pools to zero to see if it dies.
!!$    istep=istep+1
!!$    IF(istep == 600)THEN
!!$       WRITE(numout,'(A,I6,I6)') '!********** KILLING PFT ',test_pft,istep
!!$       biomass(test_grid,test_pft,ileaf,icarbon)=zero
!!$       biomass(test_grid,test_pft,ilabile,icarbon)=zero
!!$       biomass(test_grid,test_pft,icarbres,icarbon)=zero
!!$       circ_class_biomass(test_grid,test_pft,:,ileaf,icarbon)=zero
!!$       circ_class_biomass(test_grid,test_pft,:,ilabile,icarbon)=zero
!!$       circ_class_biomass(test_grid,test_pft,:,icarbres,icarbon)=zero
!!$    ENDIF
!!$    !+++++++++++++++++++++++++

    ! Save in history file the variables describing the biomass allocated to the plant parts

!DSGbalance
!    WRITE(6,*) 'DSG01: BM_ALLOC_LEAF:', bm_alloc(test_grid,test_pft,ileaf,:)
!    WRITE(6,*) 'DSG01: BIOMASS(LEAF):', biomass(test_grid,test_pft,ileaf,:)

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
       CALL histwrite_p (hist_id_stomate, 'BM_ALLOC_LEAF'//TRIM(element_str(l)), itime, &
            bm_alloc(:,:,ileaf,l), npts*nvm, horipft_index)
       CALL histwrite_p (hist_id_stomate, 'BM_ALLOC_SAP_AB'//TRIM(element_str(l)), itime, &
            bm_alloc(:,:,isapabove,l), npts*nvm, horipft_index)
       CALL histwrite_p (hist_id_stomate, 'BM_ALLOC_SAP_BE'//TRIM(element_str(l)), itime, &
            bm_alloc(:,:,isapbelow,l), npts*nvm, horipft_index)

       ! this is given out in stomate_turnover:
 !DSG: no NPP goes to heart      CALL histwrite_p (hist_id_stomate, 'BM_ALLOC_HEART_AB'//TRIM(element_str(l)), itime, &
 !DSG: no NPP goes to heart           bm_alloc(:,:,iheartabove,l), npts*nvm, horipft_index)
 !DSG: no NPP goes to heart      CALL histwrite_p (hist_id_stomate, 'BM_ALLOC_HEART_BE'//TRIM(element_str(l)), itime, &
 !DSG: no NPP goes to heart           bm_alloc(:,:,iheartbelow,l), npts*nvm, horipft_index)

       CALL histwrite_p (hist_id_stomate, 'BM_ALLOC_ROOT'//TRIM(element_str(l)), itime, &
            bm_alloc(:,:,iroot,l), npts*nvm, horipft_index)
       CALL histwrite_p (hist_id_stomate, 'BM_ALLOC_FRUIT'//TRIM(element_str(l)), itime, &
            bm_alloc(:,:,ifruit,l), npts*nvm, horipft_index)
      !DSG added: 
       CALL histwrite_p (hist_id_stomate, 'BM_ALLOC_RES'//TRIM(element_str(l)), itime, &
            bm_alloc(:,:,icarbres,l), npts*nvm, horipft_index)
       CALL histwrite_p (hist_id_stomate, 'BM_ALLOC_LAB'//TRIM(element_str(l)), itime, &
            bm_alloc(:,:,ilabile,l), npts*nvm, horipft_index)
      !DSG added: 

       CALL xios_orchidee_send_field('BM_ALLOC_LEAF'//TRIM(element_str(l)), bm_alloc(:,:,ileaf,l))
       CALL xios_orchidee_send_field('BM_ALLOC_SAP_AB'//TRIM(element_str(l)), bm_alloc(:,:,isapabove,l))
       CALL xios_orchidee_send_field('BM_ALLOC_SAP_BE'//TRIM(element_str(l)), bm_alloc(:,:,isapbelow,l))
       ! this is given out in stomate_turnover:
  !     CALL xios_orchidee_send_field('BM_ALLOC_HEART_AB'//TRIM(element_str(l)), bm_alloc(:,:,iheartabove,l))
  !     CALL xios_orchidee_send_field('BM_ALLOC_HEART_BE'//TRIM(element_str(l)), bm_alloc(:,:,iheartbelow,l))
       CALL xios_orchidee_send_field('BM_ALLOC_ROOT'//TRIM(element_str(l)), bm_alloc(:,:,iroot,l))
       CALL xios_orchidee_send_field('BM_ALLOC_FRUIT'//TRIM(element_str(l)), bm_alloc(:,:,ifruit,l))
       CALL xios_orchidee_send_field('BM_ALLOC_RES'//TRIM(element_str(l)), bm_alloc(:,:,icarbres,l))
       CALL xios_orchidee_send_field('BM_ALLOC_LAB'//TRIM(element_str(l)), bm_alloc(:,:,ilabile,l))

    ENDDO

    CALL histwrite_p (hist_id_stomate, 'RUE_LONGTERM', itime, &
         rue_longterm(:,:), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'KF', itime, &
         k_latosa(:,:), npts*nvm, horipft_index)
!gmjc
    CALL histwrite_p(hist_id_stomate, 'SLA_CALC', itime, &
                    sla_calc(:,:), npts*nvm, horipft_index)
    CALL xios_orchidee_send_field('SLA_CALC',sla_calc(:,:))
!end gmjc
    CALL xios_orchidee_send_field('RUE_LONGTERM', rue_longterm(:,:))
    CALL xios_orchidee_send_field('KF', k_latosa(:,:))

    IF (printlev.GE.3) THEN
       WRITE(numout,*) 'Leaving functional allocation growth'
       WRITE(numout,*) 'leaf_biomass C:',biomass(test_grid,test_pft,ileaf,icarbon)
       WRITE(numout,*) 'leaf_biomass N:',biomass(test_grid,test_pft,ileaf,initrogen)
       WRITE(numout,*) 'reserve_biomass C:',biomass(test_grid,test_pft,icarbres,icarbon)
       WRITE(numout,*) 'reserve_biomass N:',biomass(test_grid,test_pft,icarbres,initrogen)
    ENDIF     
    IF (dsg_debug) THEN
       !DSG debug: ===================================================
       WRITE(6,*) 'location in code: '
       WRITE(6,*) 'stomate: end growth_fun'
       !WRITE(6,*) 'biomass(test_grid,test_pft,icarbres,icarbon)    : ',  biomass(test_grid,test_pft,icarbres,icarbon)
       CALL check_mass(npts,biomass(:,:,:,:), 'stomate: end growth_fun')
    ENDIF
END SUBROUTINE growth_fun_all



!! ================================================================================================================================
!! FUNCTION	: func_derfunc
!!
!>\BRIEF        Calculate value for a function and its derivative
!!
!!
!! DESCRIPTION  : the routine describes the function and its derivative. Both function and derivative are used
!!              by the optimisation scheme. Hence, this function is part of the optimisation scheme and is only 
!!              called by the optimisation 
!!
!! RECENT CHANGE(S): 
!!
!! MAIN OUTPUT VARIABLE(S): f, df
!!
!! REFERENCE(S)	: Numerical recipies in Fortran 77
!! 
!! FLOWCHART : 
!! \n
!_ ================================================================================================================================
 
 SUBROUTINE func_derfunc(x, n, o, p, q, r, t, eq_num, f, df)

!! 0. Variable and parameter declaration

    !! 0.1 Input variables
    REAL, INTENT(in)                :: x           !! x value for which the function f(x) will be evaluated
    REAL, INTENT(in)                :: n,o,p,q,r,t !! Coefficients of the equation. Not all equations use all coefficients
    INTEGER, INTENT(in)             :: eq_num      !! Function i.e. f(x), g(x), ...

    !! 0.2 Output variables
    REAL, INTENT(out)               :: f           !! Value y for f(x)
    REAL, INTENT(out)               :: df          !! Value y for derivative[f(x)]   

    !! 0.3 Modified variables

    !! 0.4 Local variables
!_ ================================================================================================================================

!! 1. Calculate f(x) and df(x)

    IF (eq_num .EQ. 1) THEN

       !f = n*x**4 + o*x**3 + p*x**2 + q*x + r  
       !df = 4*n*x**3 + 3*o*x**2 + 2*p*x + q
    
    ELSEIF (eq_num .EQ. 2) THEN
    
       f = ( (n*x)/(p*((x+o)/t)**(q/(2+q))) ) - r
       df = ( n*(o*(q+2)+2*x)*((o+x)/t)**(-q/(q+2)) ) / ( p*(q+2)*(o+x) )
    
    ENDIF

 END SUBROUTINE func_derfunc


!! ================================================================================================================================
!! FUNCTION	: iterative_solver
!!
!>\BRIEF        find best fitting x for f(x)
!!
!!
!! DESCRIPTION  : The function makes use of an iterative approach to optimise the value for X. The solver
!!              splits the search region in two but there is an additional check to ensure that bounds are not
!!              exceeded.   
!!
!! RECENT CHANGE(S): 
!!
!! MAIN OUTPUT VARIABLE(S): x
!!
!! REFERENCE(S)	: Numerical recipies in Fortran 77
!! 
!! FLOWCHART : 
!! \n
!_ ================================================================================================================================

  FUNCTION newX(n, o, p, q, r, t, x1, x2, eq_num, j, ipts)

!! 0. Variable and parameter declaration

    !! 0.1 Input variables
    REAL, INTENT(in)        :: n,o,p,q,r,t      !! Coefficients of the equation. Not all 
                                                       !! equations use all coefficients
    REAL, INTENT(in)        :: x1               !! Lower boundary off search range
    REAL, INTENT(in)        :: x2               !! Upper boundary off search range
    INTEGER, INTENT(in)     :: eq_num           !! Function for which an iterative solution is 
                                                       !! searched
    INTEGER, INTENT(in)     :: j                !! Number of PFT
    INTEGER, INTENT(in)     :: ipts             !! Number of grdi square...for debugging
    
    !! 0.2 Output variables
   
    !! 0.3 Modified variables

    !! 0.4 Local variables
    INTEGER, PARAMETER      :: maxit = 20       !! Maximum number of iterations
    INTEGER, PARAMETER      :: max_attempt = 5  !! Maximum number of iterations
    INTEGER                 :: i, attempt       !! Index
    REAL                    :: newX             !! New estimate for X
    REAL                    :: fl, fh, f        !! Value of the function for the lower bound (x1), 
                                                       !! upper bound (x2) and the new value (newX)
    REAL                    :: xh, xl           !! Checked lower and upper bounds
    REAL                    :: df               !! Value of the derivative of the function for newX
    REAL                    :: dx, dxold        !! Slope of improvement
    REAL                    :: temp             !! Dummy variable for value swaps
    REAL                    :: low, high        !! temporary variables for x1 and x2 to avoid
                                                       !! intent in/out conflicts with Cs
    LOGICAL                        :: found_range      !! Flag indicating whether the range in which
                                                       !! a solution exists was identified.
    
    
!_ ================================================================================================================================   

!! 1. Find solution for X
 
    ! Not sure whether our initial range is large enough. We will
    ! start with a narrow range so we are more likely to fine the
    ! solution witin ::maxit iterations. If there is no solution
    ! in the initial range we will expande the range and try again

    ! Initilaze flags and counters
    attempt = 2
    found_range = .FALSE.
    low = x1
    high = x2
    

    ! Calculate y for the upper and lower bound
    DO WHILE (.NOT. found_range .AND. attempt .LT. max_attempt)
 
       CALL func_derfunc(low, n, o, p, q, r, t, eq_num, fl, df)
       CALL func_derfunc(high, n, o, p, q, r, t, eq_num, fh, df)
 
       IF ((fl .GT. 0.0 .AND. fh .GT. 0.0) .OR. &
            (fl .LT. 0.0 .AND. fh .LT. 0.0)) THEN
       
          IF (attempt .GT. max_attempt) THEN

             ! If the sign of y does not changes between the upper 
             ! and lower bound there no solution with the specified range
             WRITE(numout,*) 'Iterative procedure - tried really hard but' 
             WRITE(numout,*) 'no solution exists within the specified range'
             WRITE(numout,*) 'PFT, grid square: ',j,ipts
             CALL ipslerr_p (3,'growth_fun_all','newX',&
                  'Iterative procedure - tried really hard but failed','')

          ELSE
 
             ! Update counter
             attempt = attempt + 1
             
             ! Use previous upper boundary as the lower
             ! boundary for the next range search. Increase
             ! the upper boundary
             temp = high
             high = x1 * attempt
             low = temp
             
             ! Enlarge the search range
!!$             WRITE(numout,*) 'Iterative procedure - enlarge the search range'
!!$             WRITE(numout,*) 'New range: ', x1, x2
!!$             WRITE(numout,*) 'PFT, grid square, range: ',j,ipts,attempt

          ENDIF
       
       ELSE

          found_range = .TRUE.
          
       ENDIF
      
    ENDDO

    ! Only when we found a range we will search for the solution 
    IF (found_range) THEN
       
       ! If the sign of y changes between the upper and lower bound there is a solution
       IF ( ABS(fl) .LT. min_stomate ) THEN          

          ! The lower bound is the solution - most likely the lower bound is too high
          newX = x1
          RETURN

       ELSEIF ( ABS(fh) .LT. min_stomate ) THEN

          ! The upper bound is the solution - most likely the upper bound is too low
          newX = x2
          RETURN

       ELSEIF (fl .LT. 0.0) THEN
          
          ! Accept the lower and upper bounds as specified
          xl = x1
          xh = x2
       ELSE

          ! Lower and upper bounds were swapped, correct their ranking 
          xh = x1
          xl = x2
       ENDIF

       ! Estimate the initial newX value 	
       newX = 0.5 * (x1+x2)
       dxold = ABS(x2-x1)
       dx = dxold

    ENDIF
   
    ! Calculate y=f(x) and df(x) for initial guess of newX
    CALL func_derfunc(newX, n, o, p, q, r, t, eq_num, f, df)
    
    ! Evaluate for the maximum number of iterations  
    DO  i = 1,maxit

       IF ( ((newX-xh)*df-f)*((newX-xl)*df-f) .GT. 0.0 .OR. ABS(deux*f) > ABS(dxold*df) ) THEN
             
          ! Bisection
          dxold = dx
          dx = 0.5 * (xh-xl)
          newX = xl+dx
          IF (xl .EQ. newX) RETURN

       ELSE
             
          ! Newton
          dxold = dx
          dx = f/df
          temp = newX
          newX = newX-dx
          IF (temp .EQ. newX) RETURN

       ENDIF

       ! Precision reached
       IF ( ABS(dx) .LT. min_stomate) RETURN
          
       ! Precision was not reached calculate f(x) and df(x) for newX
       CALL func_derfunc(newX, n, o, p, q, r, t, eq_num, f, df)
          
       ! Narrow down the range
       IF (f .LT. 0.0) then
          xl = newX
       ELSE
          xh = newX
       ENDIF

    ENDDO ! maximum number of iterations
       
    !---TEMP---
    IF (j.EQ. test_pft) THEN 
       WRITE(numout,*) 'Iterative procedure: exceeded maximum iterations'
    ENDIF
    !----------

  END FUNCTION newX


!! ================================================================================================================================
!! SUBROUTINE	: comments
!!
!>\BRIEF        Contains all comments to check the code
!!
!!
!! DESCRIPTION  : contains all comments to check the code. By setting pft_test to 0, this routine is not called 
!!
!! RECENT CHANGE(S): 
!!
!! MAIN OUTPUT VARIABLE(S): none
!!
!! REFERENCE(S)	: none
!! 
!! FLOWCHART : 
!! \n
!_ ================================================================================================================================

  SUBROUTINE comment(npts, Cl_target, Cl, Cs_target, & 
       Cs, Cr_target, Cr, delta_ba, &
       ipts, j, l, b_inc_tot, & 
       Cl_incp, Cs_incp, Cr_incp, KF, LF, &
       Cl_inc, Cs_inc, Cr_inc, Cf_inc, &
       grow_wood, circ_class_n, ind, n_comment)

    !! 0. Variable and parameter declaration

    !! 0.1 Input variables
    INTEGER, INTENT(in)                         :: npts		                    !! Defined in stomate_growth_fun_all
    REAL, DIMENSION(:), INTENT(in)              :: Cl_target, Cs_target, Cr_target   !! Defined in stomate_growth_fun_all
    REAL, DIMENSION(:), INTENT(in)              :: Cl_incp, Cs_incp, Cr_incp         !! Defined in stomate_growth_fun_all
    REAL, DIMENSION(:), INTENT(in)              :: Cl_inc, Cs_inc, Cr_inc, Cf_inc    !! Defined in stomate_growth_fun_all
    REAL, DIMENSION(:,:,:), INTENT(in)          :: circ_class_n                      !! Defined in stomate_growth_fun_all
    REAL, DIMENSION(:), INTENT(in)              :: Cl, Cs, Cr                        !! Defined in stomate_growth_fun_all
    REAL, DIMENSION(:), INTENT(in)              :: delta_ba                          !! Defined in stomate_growth_fun_all
    REAL, DIMENSION(:,:), INTENT(in)            :: KF, LF, ind                       !! Defined in stomate_growth_fun_all
    REAL, INTENT(in)                            :: b_inc_tot                         !! Defined in stomate_growth_fun_all
    INTEGER, INTENT(in)                         :: ipts, j, l                        !! Defined in stomate_growth_fun_all
    LOGICAL, INTENT(in)                                :: grow_wood                         !! Defined in stomate_growth_fun_all

    !! 0.2 Output variables

    !! 0.3 Modified variables

    !! 0.4 Local variables
    INTEGER                                     :: n_comment                         !! Comment number  
    !_ ================================================================================================================================

    SELECT CASE (n_comment)
    CASE (1)
       ! Enough leaves and wood, grow roots
       WRITE(numout,*) 'Exc 1: Cl_incp(=0), Cs_incp (=0), Cr_incp (<>0), unallocated, class, '
       WRITE(numout,*) Cl_incp(l), Cs_incp(l), Cr_incp(l), &
            b_inc_tot - (circ_class_n(ipts,j,l)*(Cl_incp(l)+Cs_incp(l)+Cr_incp(l)) ), l
       IF (b_inc_tot - (circ_class_n(ipts,j,l)*(Cl_incp(l)+Cs_incp(l)+Cr_incp(l)) ) .LT. -min_stomate) THEN
          WRITE(numout,*) 'Exc 1.1: unallocated less then 0: overspending, ', b_inc_tot - &
               (circ_class_n(ipts,j,l)*(Cl_incp(l)+Cs_incp(l)+Cr_incp(l)))
       ELSE
          IF ( (b_inc_tot - ( circ_class_n(ipts,j,l)*(Cl_incp(l)+Cs_incp(l)+Cr_incp(l)) ) .GE. zero) .AND. &
               ((circ_class_n(ipts,j,l) * ABS(Cr_target(l)-Cr(l)-Cr_incp(l))) .LE. min_stomate) ) THEN
             WRITE(numout,*) 'Exc 1.2: unallocated <>= 0 but tree is in good shape: successful allocation'
          ELSEIF ( (b_inc_tot - ( circ_class_n(ipts,j,l)*(Cl_incp(l)+Cs_incp(l)+Cr_incp(l)) ) &
               .LE. min_stomate) .AND. &
               (circ_class_n(ipts,j,l) * ABS(Cr_target(l)-Cr(l)-Cr_incp(l)) .GT. min_stomate) ) THEN
             WRITE(numout,*) 'Exc 1.3: unallocated = 0, but the tree needs more reshaping: successful allocation'
          ELSE
             WRITE(numout,*) 'WARNING 24: Exc 1.4 unexpected result'
             WRITE(numout,*) 'WARNING 24: PFT, ipts: ',j,ipts
          ENDIF
       ENDIF

    CASE (2)
       ! Enough wood and roots, grow leaves 
       WRITE(numout,*) 'Exc 2: Cl_incp(<>0), Cs_incp (=0), Cr_incp (=0), unallocated, class, '
       WRITE(numout,*) Cl_incp(l), Cs_incp(l), Cr_incp(l), &
            b_inc_tot - (circ_class_n(ipts,j,l)*(Cl_incp(l)+Cs_incp(l)+Cr_incp(l)) ), l
       IF (b_inc_tot - (circ_class_n(ipts,j,l)*(Cl_incp(l)+Cs_incp(l)+Cr_incp(l)) ) .LT. -min_stomate) THEN
          WRITE(numout,*) 'Exc 2.1: unallocated less then 0: overspending, ', b_inc_tot - &
               (circ_class_n(ipts,j,l)*(Cl_incp(l)+Cs_incp(l)+Cr_incp(l)))
       ELSE
          IF ( (b_inc_tot - ( circ_class_n(ipts,j,l)*(Cl_incp(l)+Cs_incp(l)+Cr_incp(l)) ) .GE. zero) .AND. &
               ((circ_class_n(ipts,j,l) * ABS(Cl_target(l)-Cl(l)-Cl_incp(l))) .LE. min_stomate) ) THEN
             WRITE(numout,*) 'Exc 2.2: unallocated <>= 0 but tree is in good shape: successful allocation'
          ELSEIF ( (b_inc_tot - ( circ_class_n(ipts,j,l)*(Cl_incp(l)+Cs_incp(l)+Cr_incp(l)) ) &
               .LE. min_stomate) .AND. &
               ((circ_class_n(ipts,j,l) * ABS(Cl_target(l)-Cl(l)-Cl_incp(l))) .GT. min_stomate) ) THEN
             WRITE(numout,*) 'Exc 2.3: unallocated = 0, but the tree needs more reshaping: successful allocation'
          ELSE
             WRITE(numout,*) 'WARNING 25: Exc 2.4 unexpected result'
             WRITE(numout,*) 'WARNING 25: PFT, ipts: ',j,ipts
          ENDIF
       ENDIF


    CASE (3)

       ! Enough wood, grow leaves and roots
       WRITE(numout,*) 'Exc 3: Cl_incp(<>0), Cs_incp(=0), Cr_incp(<>0), unallocated, class, '
       WRITE(numout,*) Cl_incp(l), Cs_incp(l), Cr_incp(l), b_inc_tot - & 
            (circ_class_n(ipts,j,l) * (Cl_incp(l)+Cs_incp(l)+Cr_incp(l))), l
       IF (b_inc_tot - circ_class_n(ipts,j,l) * (Cl_incp(l) + Cs_incp(l) + Cr_incp(l))  &
            .LT. -min_stomate) THEN
          WRITE(numout,*) 'Exc 3.1: unallocated less then 0: overspending, ', b_inc_tot - &
               (circ_class_n(ipts,j,l) * (Cl_incp(l)+Cs_incp(l)+Cr_incp(l)) )
       ELSE
          IF ( (b_inc_tot - (circ_class_n(ipts,j,l)*(Cl_incp(l)+Cs_incp(l)+Cr_incp(l))) &
               .GE. min_stomate) .AND. &
               ((circ_class_n(ipts,j,l) * ABS(Cl_target(l)-Cl(l)-Cl_incp(l)) ) .LE. min_stomate) .AND. &
               ((circ_class_n(ipts,j,l) * ABS(Cr_target(l)-Cr(l)-Cr_incp(l)) ) .LE. min_stomate) ) THEN
             WRITE(numout,*) 'Exc 3.2: unallocated <>= 0 but tree is in good shape: successful allocation'
          ELSEIF ( (b_inc_tot - (circ_class_n(ipts,j,l)*(Cl_incp(l)+Cs_incp(l)+Cr_incp(l)) ) &
               .LE. min_stomate) .AND. &
               ((circ_class_n(ipts,j,l) * ABS(Cl_target(l)-Cl(l)-Cl_incp(l)) ) .GT. min_stomate) .AND. &
               ((circ_class_n(ipts,j,l) * ABS(Cr_target(l)-Cr(l)-Cr_incp(l)) ) .GT. min_stomate) ) THEN
             WRITE(numout,*) 'Exc 3.3: unallocated = 0, but the tree needs more reshaping: successful allocation'
          ELSE
             WRITE(numout,*) 'WARNING 26: Exc 3.4 unexpected result'
             WRITE(numout,*) 'WARNING 26: PFT, ipts: ',j,ipts
          ENDIF
       ENDIF

    CASE(4)
       ! Enough leaves and wood, grow roots
       WRITE(numout,*) 'Exc 4: Cl_incp(=0), Cs_incp (=0), Cr_incp (<>0), unallocated, class, '
       WRITE(numout,*) Cl_incp(l), Cs_incp(l), Cr_incp(l), &
            b_inc_tot - (circ_class_n(ipts,j,l)*(Cl_incp(l)+Cs_incp(l)+Cr_incp(l)) ), l
       IF (b_inc_tot - (circ_class_n(ipts,j,l)*(Cl_incp(l)+Cs_incp(l)+Cr_incp(l)) ) .LT. -min_stomate) THEN
          WRITE(numout,*) 'Exc 4.1: unallocated less then 0: overspending, ', b_inc_tot - &
               (circ_class_n(ipts,j,l)*(Cl_incp(l)+Cs_incp(l)+Cr_incp(l)))
       ELSE
          IF ( (b_inc_tot - ( circ_class_n(ipts,j,l)*(Cl_incp(l)+Cs_incp(l)+Cr_incp(l)) ) .GE. zero) .AND. &
               ((circ_class_n(ipts,j,l) * ABS(Cr_target(l)-Cr(l)-Cr_incp(l))) .LE. min_stomate) ) THEN
             WRITE(numout,*) 'Exc 4.2: unallocated <>= 0 but tree is in good shape: successful allocation'
          ELSEIF ( (b_inc_tot - ( circ_class_n(ipts,j,l)*(Cl_incp(l)+Cs_incp(l)+Cr_incp(l)) ) &
               .LE. min_stomate) .AND. &
               ((circ_class_n(ipts,j,l) * ABS(Cr_target(l)-Cr(l)-Cr_incp(l))) .GT. min_stomate) ) THEN
             WRITE(numout,*) 'Exc 4.3: unallocated = 0, but the tree needs more reshaping: successful allocation'
          ELSE
             WRITE(numout,*) 'WARNING 27: Exc 4.4 unexpected result'
             WRITE(numout,*) 'WARNING 27: PFT, ipts: ',j,ipts
          ENDIF
       ENDIF

    CASE(5)
       ! Enough leaves and roots, grow wood
       WRITE(numout,*) 'Exc 5: Cl_incp(=0), Cs_incp (<>0), Cr_incp (=0), unallocated, class, '
       WRITE(numout,*) Cl_incp(l), Cs_incp(l), Cr_incp(l), &
            b_inc_tot - (circ_class_n(ipts,j,l)*(Cl_incp(l)+Cs_incp(l)+Cr_incp(l)) ), l
       IF (b_inc_tot - (circ_class_n(ipts,j,l)*(Cl_incp(l)+Cs_incp(l)+Cr_incp(l)) ) .LT. -min_stomate) THEN
          WRITE(numout,*) 'Exc 5.1: unallocated less then 0: overspending, ', b_inc_tot - &
               (circ_class_n(ipts,j,l)*(Cl_incp(l)+Cs_incp(l)+Cr_incp(l)))
       ELSE
          IF ( (b_inc_tot - ( circ_class_n(ipts,j,l)*(Cl_incp(l)+Cs_incp(l)+Cr_incp(l)) ) .GE. zero) .AND. &
               ((circ_class_n(ipts,j,l) * ABS(Cs_target(l)-Cs(l)-Cs_incp(l))) .LE. min_stomate) ) THEN
             WRITE(numout,*) 'Exc 5.2: unallocated <>= 0 but tree is in good shape: successful allocation'
          ELSEIF ( (b_inc_tot - ( circ_class_n(ipts,j,l)*(Cl_incp(l)+Cs_incp(l)+Cr_incp(l)) ) &
               .LE. min_stomate) .AND. &
               ((circ_class_n(ipts,j,l) * ABS(Cs_target(l)-Cs(l)-Cs_incp(l))) .GT. min_stomate) ) THEN
             WRITE(numout,*) 'Exc 5.3: unallocated = 0, but the tree needs more reshaping: successful allocation'
          ELSE
             WRITE(numout,*) 'WARNING 28: Exc 5.4 unexpected result'
             WRITE(numout,*) 'WARNING 28: PFT, ipts: ',j,ipts
          ENDIF
       ENDIF

    CASE(6)
       ! Enough leaves, grow wood and roots
       WRITE(numout,*) 'Exc 6: Cl_incp(=0), Cs_incp(<>0), Cr_incp(<>0), unallocated'
       WRITE(numout,*) Cl_incp(l), Cs_incp(l), Cr_incp(l), &
            b_inc_tot - (circ_class_n(ipts,j,l) * (Cl_incp(l)+Cs_incp(l)+Cr_incp(l)))
       IF (b_inc_tot - (circ_class_n(ipts,j,l) * (Cl_incp(l)+Cs_incp(l)+Cr_incp(l))) .LT. -min_stomate) THEN
          WRITE(numout,*) 'Exc 6.1: unallocated less then 0: overspending, ', &
               b_inc_tot - (circ_class_n(ipts,j,l) * (Cl_incp(l)+Cs_incp(l)+Cr_incp(l)))
       ELSE
          IF ( (b_inc_tot - (circ_class_n(ipts,j,l)*(Cl_incp(l)+Cs_incp(l)+Cr_incp(l))) .GE. zero) .AND. &
               ((circ_class_n(ipts,j,l) * ABS(Cs_target(l)-Cs(l)-Cs_incp(l))) .LE. min_stomate) .AND. &
               ((circ_class_n(ipts,j,l) * ABS(Cr_target(l)-Cr(l)-Cr_incp(l))) .LE. min_stomate) ) THEN
             WRITE(numout,*) 'Exc 6.2: unallocated <>= 0 but tree is in good shape: successful allocation'
          ELSEIF ( (b_inc_tot - (circ_class_n(ipts,j,l)*(Cl_incp(l)+Cs_incp(l)+Cr_incp(l))) .LE. min_stomate) .AND. &
               ((circ_class_n(ipts,j,l) * ABS(Cs_target(l)-Cs(l)-Cs_incp(l))) .GT. min_stomate) .OR. &
               ((circ_class_n(ipts,j,l) * ABS(Cr_target(l)-Cr(l)-Cr_incp(l))) .GT. min_stomate) ) THEN
             WRITE(numout,*) &
                  'Exc 6.3: unallocated = 0, but the tree needs more reshaping: successful allocation'
          ELSE
             WRITE(numout,*) 'WARNING 29: Exc 6.4 unexpected result'
             WRITE(numout,*) 'WARNING 29: PFT, ipts: ',j,ipts
          ENDIF
       ENDIF

    CASE(7)
       ! Enough leaves and wood, grow roots
       WRITE(numout,*) 'Exc 7: Cl_incp(=0), Cs_incp (=0), Cr_incp (<>0), unallocated, class, '
       WRITE(numout,*) Cl_incp(l), Cs_incp(l), Cr_incp(l), &
            b_inc_tot - (circ_class_n(ipts,j,l)*(Cl_incp(l)+Cs_incp(l)+Cr_incp(l)) ), l
       IF (b_inc_tot - (circ_class_n(ipts,j,l)*(Cl_incp(l)+Cs_incp(l)+Cr_incp(l)) ) .LT. -min_stomate) THEN
          WRITE(numout,*) 'Exc 7.1: unallocated less then 0: overspending, ', b_inc_tot - &
               (circ_class_n(ipts,j,l)*(Cl_incp(l)+Cs_incp(l)+Cr_incp(l)))
       ELSE
          IF ( (b_inc_tot - ( circ_class_n(ipts,j,l)*(Cl_incp(l)+Cs_incp(l)+Cr_incp(l)) ) .GE. zero) .AND. &
               ((circ_class_n(ipts,j,l) * ABS(Cl_target(l)-Cl(l)-Cl_incp(l))) .LE. min_stomate) ) THEN
             WRITE(numout,*) 'Exc 7.2: unallocated <>= 0 but tree is in good shape: successful allocation'
          ELSEIF ( (b_inc_tot - ( circ_class_n(ipts,j,l)*(Cl_incp(l)+Cs_incp(l)+Cr_incp(l)) ) &
               .LE. min_stomate) .AND. &
               ((circ_class_n(ipts,j,l) * ABS(Cl_target(l)-Cl(l)-Cl_incp(l))) .GT. min_stomate) ) THEN
             WRITE(numout,*) 'Exc 7.3: unallocated = 0, but the tree needs more reshaping: successful allocation'
          ELSE
             WRITE(numout,*) 'WARNING 30: Exc 7.4 unexpected result'
             WRITE(numout,*) 'WARNING 30: PFT, ipts: ',j,ipts
          ENDIF
       ENDIF

    CASE(8)
       ! Enough leaves and roots, grow wood
       WRITE(numout,*) 'Exc 8: Cl_incp(=0), Cs_incp (<>0), Cr_incp (=0), unallocated, class, '
       WRITE(numout,*) Cl_incp(l), Cs_incp(l), Cr_incp(l), &
            b_inc_tot - (circ_class_n(ipts,j,l)*(Cl_incp(l)+Cs_incp(l)+Cr_incp(l)) ), l
       IF (b_inc_tot - (circ_class_n(ipts,j,l)*(Cl_incp(l)+Cs_incp(l)+Cr_incp(l)) ) .LT. -min_stomate) THEN
          WRITE(numout,*) 'Exc 8.1: unallocated less then 0: overspending, ', b_inc_tot - &
               (circ_class_n(ipts,j,l)*(Cl_incp(l)+Cs_incp(l)+Cr_incp(l)))
       ELSE
          IF ( (b_inc_tot - ( circ_class_n(ipts,j,l)*(Cl_incp(l)+Cs_incp(l)+Cr_incp(l)) ) .GE. zero) .AND. &
               ((circ_class_n(ipts,j,l) * ABS(Cs_target(l)-Cs(l)-Cs_incp(l))) .LE. min_stomate) ) THEN
             WRITE(numout,*) 'Exc 8.2: unallocated <>= 0 but tree is in good shape: successful allocation'
          ELSEIF ( (b_inc_tot - ( circ_class_n(ipts,j,l)*(Cl_incp(l)+Cs_incp(l)+Cr_incp(l)) ) &
               .LE. min_stomate) .AND. &
               ((circ_class_n(ipts,j,l) * ABS(Cs_target(l)-Cs(l)-Cs_incp(l))) .GT. min_stomate) ) THEN
             WRITE(numout,*) 'Exc 8.3: unallocated = 0, but the tree needs more reshaping: successful allocation'
          ELSE
             WRITE(numout,*) 'WARNING 31: Exc 8.4 unexpected result'
             WRITE(numout,*) 'WARNING 31: PFT, ipts: ',j,ipts
          ENDIF
       ENDIF

    CASE(9)
       ! Enough roots, grow leaves and wood
       WRITE(numout,*) 'Exc 9: delta_ba, Cl_incp(<>0), Cs_incp(<>0), Cr_incp(=0), unallocated, class, '
       WRITE(numout,*) delta_ba(:), Cl_incp(l), Cs_incp(l), Cr_incp(l), &
            b_inc_tot - (circ_class_n(ipts,j,l) * (Cl_incp(l)+Cs_incp(l)+Cr_incp(l))), l
       IF (b_inc_tot - (circ_class_n(ipts,j,l) * (Cl_incp(l)+Cs_incp(l)+Cr_incp(l))) .LT. -min_stomate) THEN
          WRITE(numout,*) 'Exc 9.1: unallocated less then 0: overspending, ', &
               b_inc_tot - (circ_class_n(ipts,j,l) * (Cl_incp(l)+Cs_incp(l)+Cr_incp(l)))
       ELSE
          IF ( (b_inc_tot - (circ_class_n(ipts,j,l) * (Cl_incp(l)+Cs_incp(l)+Cr_incp(l))) .GE. zero) .AND. &
               ((circ_class_n(ipts,j,l) * ABS(Cl_target(l)-Cl(l)-Cl_incp(l))) .LE. min_stomate) .AND. &
               ((circ_class_n(ipts,j,l) * ABS(Cs_target(l)-Cs(l)-Cs_incp(l))) .LE. min_stomate) ) THEN
             WRITE(numout,*) 'Exc 9.2: unallocated <>= 0 but tree is in good shape: successful allocation'
          ELSEIF ( (b_inc_tot - (circ_class_n(ipts,j,l) * (Cl_incp(l)+Cs_incp(l)+Cr_incp(l))) &
               .LE. min_stomate) .AND. &
               ((circ_class_n(ipts,j,l) * ABS(Cl_target(l)-Cl(l)-Cl_incp(l))) .GT. min_stomate) .OR. &
               ((circ_class_n(ipts,j,l) * ABS(Cs_target(l)-Cs(l)-Cs_incp(l))) .GT. min_stomate) ) THEN
             WRITE(numout,*) 'Exc 9.3: unallocated = 0, but the tree needs more reshaping: successful allocation'
          ELSE
             WRITE(numout,*) 'WARNING 32: Exc 9.4 unexpected result'
             WRITE(numout,*) 'WARNING 32: PFT, ipts: ',j,ipts
          ENDIF
       ENDIF

    CASE(10)
       ! Ready for ordinary allocation
       WRITE(numout,*) 'Ready for ordinary allocation?'
       WRITE(numout,*) 'KF, LF, ', KF(ipts,j), LF(ipts,j)
       WRITE(numout,*) 'b_inc_tot, ', b_inc_tot
       WRITE(numout,*) 'Cl, Cs, Cr', Cl(:), Cs(:), Cr(:)
       WRITE(numout,*) 'Cl_target-Cl, ', Cl_target(:)-Cl(:)
       WRITE(numout,*) 'Cs_target-Cs, ', Cs_target(:)-Cs(:)
       WRITE(numout,*) 'Cr_target-Cr, ', Cr_target(:)-Cr(:)
       IF (b_inc_tot .GT. min_stomate) THEN
          IF (SUM(ABS(Cl_target(:)-Cl(:))) .LE. min_stomate) THEN
             IF (SUM(ABS(Cs_target(:)-Cs(:))) .LE. min_stomate) THEN
                IF (SUM(ABS(Cr_target(:)-Cr(:))) .LE. min_stomate) THEN
                   IF (grow_wood) THEN
                      WRITE(numout,*) 'should result in exc 10.1 or 10.2'
                   ELSE
                      WRITE(numout,*) 'No wood growth.  Not a problem!  Just an observation.'
                   ENDIF
                ELSE
                   WRITE(numout,*) 'WARNING 34: problem with Cr_target'
                   WRITE(numout,*) 'WARNING 34: PFT, ipts: ',j,ipts
                ENDIF
             ELSE
                WRITE(numout,*) 'WARNING 35: problem with Cs_target'
                WRITE(numout,*) 'WARNING 35: PFT, ipts: ',j,ipts
             ENDIF
          ELSE
             WRITE(numout,*) 'WARNING 36: problem with Cl_target'
             WRITE(numout,*) 'WARNING 36: PFT, ipts: ',j,ipts
          ENDIF
       ELSEIF(b_inc_tot .LT. -min_stomate) THEN 
          WRITE(numout,*) 'WARNING 37: problem with b_inc_tot'
          WRITE(numout,*) 'WARNING 37: PFT, ipts: ',j,ipts
       ELSE
          WRITE(numout,*) 'no unallocated fraction'
       ENDIF

    CASE(11)
       ! Ordinary allocation
       WRITE(numout,*) 'delta_ba, ', delta_ba
       IF ( (SUM(Cl_inc(:)) .GE. zero) .AND. (SUM(Cs_inc(:)) .GE. zero) .AND. &
            (SUM(Cr_inc(:)) .GE. zero) .AND. &
            ( b_inc_tot - SUM(circ_class_n(ipts,j,:) * (Cl_inc(:)+Cs_inc(:)+Cr_inc(:))) .GT. -1*min_stomate) .AND. &
            ( b_inc_tot - SUM(circ_class_n(ipts,j,:) * (Cl_inc(:)+Cs_inc(:)+Cr_inc(:))) .LT. min_stomate ) ) THEN
          WRITE(numout,*) 'Exc 10.1: Ordinary allocation was succesful'
          WRITE(numout,*) 'Cl_inc, Cs_inc, Cr_inc, unallocated', Cl_inc(:), Cs_inc(:), Cr_inc(:), & 
               b_inc_tot - SUM(circ_class_n(ipts,j,:) * (Cl_inc(:)+Cs_inc(:)+Cr_inc(:)))
       ELSE
          WRITE(numout,*) 'WARNING 38: Exc 10.2 problem with ordinary allocation'
          WRITE(numout,*) 'WARNING 38: PFT, ipts: ',j,ipts
          WRITE(numout,*) 'Cl_inc, Cs_inc, Cr_inc, unallocated', Cl_inc(:), Cs_inc(:), Cr_inc(:), & 
               b_inc_tot - SUM(circ_class_n(ipts,j,:) * (Cl_inc(:)+Cs_inc(:)+Cr_inc(:)))
       ENDIF

    CASE(12)
       ! Enough leaves and structure, grow roots
       WRITE(numout,*) 'Exc 1: Cl_incp(=0), Cs_incp (=0), Cr_incp (<>0), unallocated, '
       WRITE(numout,*) Cl_incp(1), Cs_incp(1), Cr_incp(1), &
            b_inc_tot - (ind(ipts,j)*(Cl_incp(1)+Cs_incp(1)+Cr_incp(1)) )
       IF (b_inc_tot - (ind(ipts,j)*(Cl_incp(1)+Cs_incp(1)+Cr_incp(1)) ) .LT. -min_stomate) THEN
          WRITE(numout,*) 'Exc 1.1: unallocated less then 0: overspending, ', b_inc_tot - &
               (ind(ipts,j)*(Cl_incp(1)+Cs_incp(1)+Cr_incp(1)))
       ELSE
          IF ( (b_inc_tot - ( ind(ipts,j)*(Cl_incp(1)+Cs_incp(1)+Cr_incp(1)) ) .GE. zero) .AND. &
               ((ind(ipts,j) * ABS(Cr_target(1)-Cr(1)-Cr_incp(1))) .LE. min_stomate) ) THEN
             WRITE(numout,*) 'Exc 1.2: unallocated <>= 0 but tree is in good shape: successful allocation'
          ELSEIF ( (b_inc_tot - ( ind(ipts,j)*(Cl_incp(1)+Cs_incp(1)+Cr_incp(1)) ) &
               .LE. min_stomate) .AND. &
               (ind(ipts,j) * ABS(Cr_target(1)-Cr(1)-Cr_incp(1)) .GT. min_stomate) ) THEN
             WRITE(numout,*) 'Exc 1.3: unallocated = 0, but the tree needs more reshaping: successful allocation'
          ELSE
             WRITE(numout,*) 'WARNING 39: Exc 1.4 unexpected result'
             WRITE(numout,*) 'WARNING 39: PFT, ipts: ',j,ipts
          ENDIF
       ENDIF

    CASE(13)
       ! Enough structural C and roots, grow leaves
       WRITE(numout,*) 'Exc 2: Cl_incp(<>0), Cs_incp (=0), Cr_incp (=0), unallocated, '
       WRITE(numout,*) Cl_incp(1), Cs_incp(1), Cr_incp(1), &
            b_inc_tot - (ind(ipts,j)*(Cl_incp(1)+Cs_incp(1)+Cr_incp(1)) )
       IF (b_inc_tot - (ind(ipts,j)*(Cl_incp(1)+Cs_incp(1)+Cr_incp(1)) ) .LT. -min_stomate) THEN
          WRITE(numout,*) 'Exc 2.1: unallocated less then 0: overspending, ', b_inc_tot - &
               (ind(ipts,j)*(Cl_incp(1)+Cs_incp(1)+Cr_incp(1)))
       ELSE
          IF ( (b_inc_tot - ( ind(ipts,j)*(Cl_incp(1)+Cs_incp(1)+Cr_incp(1)) ) .GE. zero) .AND. &
               ((ind(ipts,j) * ABS(Cl_target(1)-Cl(1)-Cl_incp(1))) .LE. min_stomate) ) THEN
             WRITE(numout,*) 'Exc 2.2: unallocated <>= 0 but tree is in good shape: successful allocation'
          ELSEIF ( (b_inc_tot - ( ind(ipts,j)*(Cl_incp(1)+Cs_incp(1)+Cr_incp(1)) ) &
               .LE. min_stomate) .AND. &
               ((ind(ipts,j) * ABS(Cl_target(1)-Cl(1)-Cl_incp(1))) .GT. min_stomate) ) THEN
             WRITE(numout,*) 'Exc 2.3: unallocated = 0, but the tree needs more reshaping: successful allocation'
          ELSE
             WRITE(numout,*) 'WARNING 40: Exc l.4 unexpected result'
             WRITE(numout,*) 'WARNING 40: PFT, ipts: ',j,ipts
          ENDIF
       ENDIF

    CASE(14)
       ! Enough structural C and root, grow leaves
       WRITE(numout,*) 'Exc 3: Cl_incp(<>0), Cs_incp(=0), Cr_incp(<>0), unallocated, '
       WRITE(numout,*) Cl_incp(1), Cs_incp(1), Cr_incp(1), b_inc_tot - & 
            (ind(ipts,j) * (Cl_incp(1)+Cs_incp(1)+Cr_incp(1)))
       IF (b_inc_tot - ind(ipts,j) * (Cl_incp(1) + Cs_incp(1) + Cr_incp(1))  &
            .LT. -min_stomate) THEN
          WRITE(numout,*) 'Exc 3.1: unallocated less then 0: overspending, ', b_inc_tot - &
               (ind(ipts,j) * (Cl_incp(1)+Cs_incp(1)+Cr_incp(1)) )
       ELSE
          IF ( (b_inc_tot - (ind(ipts,j)*(Cl_incp(1)+Cs_incp(1)+Cr_incp(1))) &
               .GE. min_stomate) .AND. &
               ((ind(ipts,j) * ABS(Cl_target(1)-Cl(1)-Cl_incp(1)) ) .LE. min_stomate) .AND. &
               ((ind(ipts,j) * ABS(Cr_target(1)-Cr(1)-Cr_incp(1)) ) .LE. min_stomate) ) THEN
             WRITE(numout,*) 'Exc 3.2: unallocated <>= 0 but tree is in good shape: successful allocation'
          ELSEIF ( (b_inc_tot - (ind(ipts,j)*(Cl_incp(1)+Cs_incp(1)+Cr_incp(1)) ) &
               .LE. min_stomate) .AND. &
               ((ind(ipts,j) * ABS(Cl_target(1)-Cl(1)-Cl_incp(1)) ) .GT. min_stomate) .AND. &
               ((ind(ipts,j) * ABS(Cr_target(1)-Cr(1)-Cr_incp(1)) ) .GT. min_stomate) ) THEN
             WRITE(numout,*) 'Exc 3.3: unallocated = 0, but the tree needs more reshaping: successful allocation'
          ELSE
             WRITE(numout,*) 'WARNING 41: Exc 3.4 unexpected result'
             WRITE(numout,*) 'WARNING 41: PFT, ipts: ',j,ipts
          ENDIF
       ENDIF

    CASE(15)
       ! Enough leaves and structural C, grow roots
       WRITE(numout,*) 'Exc 4: Cl_incp(=0), Cs_incp (=0), Cr_incp (<>0), unallocated, '
       WRITE(numout,*) Cl_incp(1), Cs_incp(1), Cr_incp(1), &
            b_inc_tot - (ind(ipts,j)*(Cl_incp(1)+Cs_incp(1)+Cr_incp(1)) )
       IF (b_inc_tot - (ind(ipts,j)*(Cl_incp(1)+Cs_incp(1)+Cr_incp(1)) ) .LT. -min_stomate) THEN
          WRITE(numout,*) 'Exc 4.1: unallocated less then 0: overspending, ', b_inc_tot - &
               (ind(ipts,j)*(Cl_incp(1)+Cs_incp(1)+Cr_incp(1)))
       ELSE
          IF ( (b_inc_tot - ( ind(ipts,j)*(Cl_incp(1)+Cs_incp(1)+Cr_incp(1)) ) .GE. zero) .AND. &
               ((ind(ipts,j) * ABS(Cr_target(1)-Cr(1)-Cr_incp(1))) .LE. min_stomate) ) THEN
             WRITE(numout,*) 'Exc 4.2: unallocated <>= 0 but tree is in good shape: successful allocation'
          ELSEIF ( (b_inc_tot - ( ind(ipts,j)*(Cl_incp(1)+Cs_incp(1)+Cr_incp(1)) ) &
               .LE. min_stomate) .AND. &
               ((ind(ipts,j) * ABS(Cr_target(1)-Cr(1)-Cr_incp(1))) .GT. min_stomate) ) THEN
             WRITE(numout,*) 'Exc 4.3: unallocated = 0, but the tree needs more reshaping: successful allocation'
          ELSE
             WRITE(numout,*) 'WARNING 42: Exc 4.4 unexpected result'
             WRITE(numout,*) 'WARNING 42: PFT, ipts: ',j,ipts
          ENDIF
       ENDIF

    CASE(16)
       ! Enough leaves and roots, grow structural C            
       WRITE(numout,*) 'Exc 5: Cl_incp(=0), Cs_incp (<>0), Cr_incp (=0), unallocated, '
       WRITE(numout,*) Cl_incp(1), Cs_incp(1), Cr_incp(1), &
            b_inc_tot - (ind(ipts,j)*(Cl_incp(1)+Cs_incp(1)+Cr_incp(1)) )
       IF (b_inc_tot - (ind(ipts,j)*(Cl_incp(1)+Cs_incp(1)+Cr_incp(1)) ) .LT. -min_stomate) THEN
          WRITE(numout,*) 'Exc 5.1: unallocated less then 0: overspending, ', b_inc_tot - &
               (ind(ipts,j)*(Cl_incp(1)+Cs_incp(1)+Cr_incp(1)))
       ELSE
          IF ( (b_inc_tot - ( ind(ipts,j)*(Cl_incp(1)+Cs_incp(1)+Cr_incp(1)) ) .GE. zero) .AND. &
               ((ind(ipts,j) * ABS(Cs_target(1)-Cs(1)-Cs_incp(1))) .LE. min_stomate) ) THEN
             WRITE(numout,*) 'Exc 5.2: unallocated <>= 0 but tree is in good shape: successful allocation'
          ELSEIF ( (b_inc_tot - ( ind(ipts,j)*(Cl_incp(1)+Cs_incp(1)+Cr_incp(1)) ) &
               .LE. min_stomate) .AND. &
               ((ind(ipts,j) * ABS(Cs_target(1)-Cs(1)-Cs_incp(1))) .GT. min_stomate) ) THEN
             WRITE(numout,*) 'Exc 5.3: unallocated = 0, but the tree needs more reshaping: successful allocation'
          ELSE
             WRITE(numout,*) 'WARNING 43: Exc 5.4 unexpected result'
             WRITE(numout,*) 'WARNING 43: PFT, ipts: ',j,ipts
          ENDIF
       ENDIF

    CASE(17)
       ! Enough leaves, grow structural C and roots
       WRITE(numout,*) 'Exc 6: Cl_incp(=0), Cs_incp(<>0), Cr_incp(<>0), unallocated'
       WRITE(numout,*) Cl_incp(1), Cs_incp(1), Cr_incp(1), &
            b_inc_tot - (ind(ipts,j) * (Cl_incp(1)+Cs_incp(1)+Cr_incp(1)))
       IF (b_inc_tot - (ind(ipts,j) * (Cl_incp(1)+Cs_incp(1)+Cr_incp(1))) .LT. -min_stomate) THEN
          WRITE(numout,*) 'Exc 6.1: unallocated less then 0: overspending, ', &
               b_inc_tot - (ind(ipts,j) * (Cl_incp(1)+Cs_incp(1)+Cr_incp(1)))
       ELSE
          IF ( (b_inc_tot - ind(ipts,j)*(Cl_incp(1)+Cs_incp(1)+Cr_incp(1)) .GE. zero) .AND. &
               ((ind(ipts,j) * ABS(Cs_target(1)-Cs(1)-Cs_incp(1))) .LE. min_stomate) .AND. &
               ((ind(ipts,j) * ABS(Cr_target(1)-Cr(1)-Cr_incp(1))) .LE. min_stomate) ) THEN
             WRITE(numout,*) 'Exc 6.2: unallocated <>= 0 but tree is in good shape: successful allocation'
          ELSEIF ( (b_inc_tot - ind(ipts,j)*(Cl_incp(1)+Cs_incp(1)+Cr_incp(1)) .LE. min_stomate) .AND. &
               ((ind(ipts,j) * ABS(Cs_target(1)-Cs(1)-Cs_incp(1))) .GT. min_stomate) .AND. &
               ((ind(ipts,j) * ABS(Cr_target(1)-Cr(1)-Cr_incp(1))) .GT. min_stomate) ) THEN
             WRITE(numout,*) &
                  'Exc 6.3: unallocated = 0, but the tree needs more reshaping: successful allocation'
          ELSE
             WRITE(numout,*) 'WARNING 44: Exc 6.4 unexpected result'
             WRITE(numout,*) 'WARNING 44: PFT, ipts: ',j,ipts
          ENDIF
       ENDIF

    CASE(18)
       ! Enough leaves and structural C, grow roots
       WRITE(numout,*) 'Exc 7: Cl_incp(=0), Cs_incp (=0), Cr_incp (<>0), unallocated, '
       WRITE(numout,*) Cl_incp(1), Cs_incp(1), Cr_incp(1), &
            b_inc_tot - (ind(ipts,j)*(Cl_incp(1)+Cs_incp(1)+Cr_incp(1)) )
       IF (b_inc_tot - (ind(ipts,j)*(Cl_incp(1)+Cs_incp(1)+Cr_incp(1)) ) .LT. -min_stomate) THEN
          WRITE(numout,*) 'Exc 7.1: unallocated less then 0: overspending, ', b_inc_tot - &
               (ind(ipts,j)*(Cl_incp(1)+Cs_incp(1)+Cr_incp(1)))
       ELSE
          IF ( (b_inc_tot - ( ind(ipts,j)*(Cl_incp(1)+Cs_incp(1)+Cr_incp(1)) ) .GE. zero) .AND. &
               ((ind(ipts,j) * ABS(Cl_target(1)-Cl(1)-Cl_incp(1))) .LE. min_stomate) ) THEN
             WRITE(numout,*) 'Exc 7.2: unallocated <>= 0 but tree is in good shape: successful allocation'
          ELSEIF ( (b_inc_tot - ( ind(ipts,j)*(Cl_incp(1)+Cs_incp(1)+Cr_incp(1)) ) &
               .LE. min_stomate) .AND. &
               ((ind(ipts,j) * ABS(Cl_target(1)-Cl(1)-Cl_incp(1))) .GT. min_stomate) ) THEN
             WRITE(numout,*) 'Exc 7.3: unallocated = 0, but the tree needs more reshaping: successful allocation'
          ELSE
             WRITE(numout,*) 'WARNING 45: Exc 7.4 unexpected result'
             WRITE(numout,*) 'WARNING 45: PFT, ipts: ',j,ipts
          ENDIF
       ENDIF

    CASE(19)
       ! Enough leaves and roots, grow structural C
       WRITE(numout,*) 'Exc 8: Cl_incp(=0), Cs_incp (<>0), Cr_incp (=0), unallocated, '
       WRITE(numout,*) Cl_incp(1), Cs_incp(1), Cr_incp(1), &
            b_inc_tot - (ind(ipts,j)*(Cl_incp(1)+Cs_incp(1)+Cr_incp(1)) )
       IF (b_inc_tot - (ind(ipts,j)*(Cl_incp(1)+Cs_incp(1)+Cr_incp(1)) ) .LT. -min_stomate) THEN
          WRITE(numout,*) 'Exc 8.1: unallocated less then 0: overspending, ', b_inc_tot - &
               (ind(ipts,j)*(Cl_incp(1)+Cs_incp(1)+Cr_incp(1)))
       ELSE
          IF ( (b_inc_tot - ( ind(ipts,j)*(Cl_incp(1)+Cs_incp(1)+Cr_incp(1)) ) .GE. zero) .AND. &
               ((ind(ipts,j) * ABS(Cs_target(1)-Cs(1)-Cs_incp(1))) .LE. min_stomate) ) THEN
             WRITE(numout,*) 'Exc 8.2: unallocated <>= 0 but tree is in good shape: successful allocation'
          ELSEIF ( (b_inc_tot - ( ind(ipts,j)*(Cl_incp(1)+Cs_incp(1)+Cr_incp(1)) ) &
               .LE. min_stomate) .AND. &
               ((ind(ipts,j) * ABS(Cs_target(1)-Cs(1)-Cs_incp(1))) .GT. min_stomate) ) THEN
             WRITE(numout,*) 'Exc 8.3: unallocated = 0, but the tree needs more reshaping: successful allocation'
          ELSE
             WRITE(numout,*) 'WARNING 46: Exc 8.4 unexpected result'
             WRITE(numout,*) 'WARNING 46: PFT, ipts: ',j,ipts
          ENDIF
       ENDIF

    CASE(20)
       ! Enough roots, grow structural C and leaves
       WRITE(numout,*) 'Exc 9: Cl_incp(<>0), Cs_incp(<>0), Cr_incp(=0), unallocated, '
       WRITE(numout,*) Cl_incp(1), Cs_incp(1), Cr_incp(1), &
            b_inc_tot - (ind(ipts,j) * (Cl_incp(1)+Cs_incp(1)+Cr_incp(1)))
       WRITE(numout,*) 'term 1', b_inc_tot - (ind(ipts,j) * (Cl_incp(1)+Cs_incp(1)+Cr_incp(1)))
       WRITE(numout,*) 'term 2', (ind(ipts,j) * ABS(Cl_target(1)-Cl(1)-Cl_incp(1)))
       WRITE(numout,*) 'term 3', (ind(ipts,j) * ABS(Cs_target(1)-Cs(1)-Cs_incp(1)))
       IF (b_inc_tot - (ind(ipts,j) * (Cl_incp(1)+Cs_incp(1)+Cr_incp(1))) .LT. -min_stomate) THEN
          WRITE(numout,*) 'Exc 9.1: unallocated less then 0: overspending, ', &
               b_inc_tot - (ind(ipts,j) * (Cl_incp(1)+Cs_incp(1)+Cr_incp(1)))
       ELSE
          IF ( (b_inc_tot - (ind(ipts,j) * (Cl_incp(1)+Cs_incp(1)+Cr_incp(1))) .GE. zero) .AND. &
               ((ind(ipts,j) * ABS(Cl_target(1)-Cl(1)-Cl_incp(1))) .LE. min_stomate) .AND. &
               ((ind(ipts,j) * ABS(Cs_target(1)-Cs(1)-Cs_incp(1))) .LE. min_stomate) ) THEN
             WRITE(numout,*) 'Exc 9.2: unallocated <>= 0 but tree is in good shape: successful allocation'
          ELSEIF ( (b_inc_tot - (ind(ipts,j) * (Cl_incp(1)+Cs_incp(1)+Cr_incp(1))) .LE. min_stomate) .AND. &
               (((ind(ipts,j) * ABS(Cl_target(1)-Cl(1)-Cl_incp(1))) .GT. min_stomate) .OR. &
               ((ind(ipts,j) * ABS(Cs_target(1)-Cs(1)-Cs_incp(1))) .GT. min_stomate) ) ) THEN
             WRITE(numout,*) 'Exc 9.3: unallocated = 0, but the tree needs more reshaping: successful allocation'
          ELSE
             WRITE(numout,*) 'WARNING 47: Exc 9.4 unexpected result'
             WRITE(numout,*) 'WARNING 47: PFT, ipts: ',j,ipts
          ENDIF
       ENDIF

    CASE(21)
       ! Ready for ordinary allocation
       WRITE(numout,*) 'Ready for ordinary allocation?'
       WRITE(numout,*) 'KF, LF, ', KF(ipts,j), LF(ipts,j)
       WRITE(numout,*) 'b_inc_tot, ', b_inc_tot
       WRITE(numout,*) 'Cl, Cs, Cr', Cl(1), Cs(1), Cr(1)
       WRITE(numout,*) 'Cl_target-Cl, ', Cl_target(1)-Cl(1)
       WRITE(numout,*) 'Cs_target-Cs, ', Cs_target(1)-Cs(1)
       WRITE(numout,*) 'Cr_target-Cr, ', Cr_target(1)-Cr(1)
       IF (b_inc_tot .GT. min_stomate) THEN
          IF (ABS(Cl_target(1)-Cl(1)) .LE. min_stomate) THEN
             IF (ABS(Cs_target(1)-Cs(1)) .LE. min_stomate) THEN
                IF (ABS(Cr_target(1)-Cr(1)) .LE. min_stomate) THEN
                   IF (b_inc_tot .GT. min_stomate) THEN
                      IF (grow_wood) THEN
                         WRITE(numout,*) 'should result in exc 10.1 or 10.2'
                      ELSE
                         WRITE(numout,*) 'WARNING 48: no wood growth'
                         WRITE(numout,*) 'WARNING 48: PFT, ipts: ',j,ipts
                      ENDIF
                   ENDIF
                ELSE
                   WRITE(numout,*) 'WARNING 49: problem with Cr_target'
                   WRITE(numout,*) 'WARNING 49: PFT, ipts: ',j,ipts
                ENDIF
             ELSE
                WRITE(numout,*) 'WARNING 50: problem with Cs_target'
                WRITE(numout,*) 'WARNING 50: PFT, ipts: ',j,ipts
             ENDIF
          ELSE
             WRITE(numout,*) 'WARNING 51: problem with Cl_target'
             WRITE(numout,*) 'WARNING 51: PFT, ipts: ',j,ipts
          ENDIF
       ELSEIF(b_inc_tot .LT. -min_stomate) THEN 
          WRITE(numout,*) 'WARNING 52: problem with b_inc_tot'
          WRITE(numout,*) 'WARNING 52: PFT, ipts: ',j,ipts
       ELSE
          WRITE(numout,*) 'no unallocated fraction'
       ENDIF

    CASE(22)
       ! Ordinary allocation
       IF ( ((Cl_inc(1)) .GE. zero) .AND. ((Cs_inc(1)) .GE. zero) .AND. &
            ((Cr_inc(1)) .GE. zero) .AND. &
            ( b_inc_tot - (ind(ipts,j) * (Cl_inc(1)+Cs_inc(1)+Cr_inc(1))) .GT. -1*min_stomate) .AND. &
            ( b_inc_tot - (ind(ipts,j) * (Cl_inc(1)+Cs_inc(1)+Cr_inc(1))) .LT. min_stomate ) ) THEN
          WRITE(numout,*) 'Exc 10.1: Ordinary allocation was succesful'
          WRITE(numout,*) 'Cl_inc, Cs_inc, Cr_inc, unallocated', Cl_inc(1), Cs_inc(1), Cr_inc(1), & 
               b_inc_tot - (ind(ipts,j) * (Cl_inc(1)+Cs_inc(1)+Cr_inc(1)))
       ELSE
          WRITE(numout,*) 'WARNING 53: Exc 10.2 problem with ordinary allocation'
          WRITE(numout,*) 'WARNING 53: PFT, ipts: ',j,ipts
          WRITE(numout,*) 'Cl_inc, Cs_inc, Cr_inc, unallocated', Cl_inc(1), Cs_inc(1), Cr_inc(1), & 
               b_inc_tot - (ind(ipts,j) * (Cl_inc(1)+Cs_inc(1)+Cr_inc(1)))
       ENDIF

    END SELECT

  END SUBROUTINE comment

END MODULE stomate_growth_fun_all
