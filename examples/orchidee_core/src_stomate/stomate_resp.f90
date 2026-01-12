! =================================================================================================================================
! MODULE           : stomate_resp
!
! CONTACT          : orchidee-help _at_ listes.ipsl.fr
!
! LICENCE          : IPSL (2006)
!                  This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF           Calculates maintenance respiration for different plant components
!!
!!\n DESCRIPTION   : None
!!
!! RECENT CHANGE(S): None
!!
!! REFERENCE(S)	   :
!!- McCree KJ. An equation for the respiration of white clover plants grown under controlled conditions. 
!! In: Setlik I, editor. Prediction and measurement of photosynthetic productivity. Wageningen, The Netherlands: Pudoc; 1970. p. 221-229.
!! - Krinner G, Viovy N, de Noblet-Ducoudre N, Ogee J, Polcher J, Friedlingstein P,
!! Ciais P, Sitch S, Prentice I C (2005) A dynamic global vegetation model for studies
!! of the coupled atmosphere-biosphere system. Global Biogeochemical Cycles, 19, GB1015,
!! doi: 10.1029/2003GB002199.\n
!! Ruimy A., Dedieu G., Saugier B. (1996), TURC: A diagnostic model
!! of continental gross primary productivity and net primary productivity,
!! Global Biogeochemical Cycles, 10, 269-285.\n

!! SVN :
!! $HeadURL: svn://forge.ipsl.jussieu.fr/orchidee/branches/ORCHIDEE-CN-P/ORCHIDEE/src_stomate/stomate_resp.f90 $
!! $Date: 2017-11-29 11:07:35 +0100 (三, 2017-11-29) $
!! $Revision: 4792 $
!! \n
!_ ================================================================================================================================
 
MODULE stomate_resp

  ! modules used:
  USE stomate_data
  USE pft_parameters
  USE constantes  
  USE constantes_soil 

  IMPLICIT NONE

  ! private & public routines
  PRIVATE
  PUBLIC maint_respiration,maint_respiration_clear

  LOGICAL, SAVE                                              :: firstcall_resp = .TRUE.                 !! first call
!$OMP THREADPRIVATE(firstcall_resp)

CONTAINS


!! ================================================================================================================================
!! SUBROUTINE 	: maint_respiration_clear
!!
!>\BRIEF        : Set the flag ::firstcall_resp to .TRUE. and as such activate section 
!!                1.1 of the subroutine maint_respiration (see below).
!_ ================================================================================================================================

  SUBROUTINE maint_respiration_clear
    firstcall_resp=.TRUE.
  END SUBROUTINE maint_respiration_clear


!! ================================================================================================================================
!! SUBROUTINE 	: maint_respiration
!!
!>\BRIEF         Calculate PFT maintenance respiration of each living plant part by 
!! multiplying the biomass of plant part by maintenance respiration coefficient which
!! depends on long term mean annual temperature. PFT maintenance respiration is carbon flux 
!! with the units @tex $(gC.m^{-2}dt_sechiba^{-1})$ @endtex, and the convention is from plants to the 
!! atmosphere.
!!
!! DESCRIPTION : The maintenance respiration of each plant part for each PFT is the biomass of the plant 
!! part multiplied by maintenance respiration coefficient. The biomass allocation to different 
!! plant parts is done in routine stomate_alloc.f90. The maintenance respiration coefficient is 
!! calculated in this routine.\n
!!
!! The maintenance respiration coefficient is the fraction of biomass that is lost during 
!! each time step, which increases linearly with temperature (2-meter air temperature for aboveground plant
!! tissues; root-zone temperature for below-ground tissues). Air temperature is an input forcing variable. 
!! Root-zone temperature is a convolution of root and soil temperature profiles and also calculated 
!! in this routine.\n
!!
!! The calculation of maintenance respiration coefficient (fraction of biomass respired) depends linearly
!! on temperature:
!! - the relevant temperature for different plant parts (air temperature or root-zone temperature)\n
!! - intercept: prescribed maintenance respiration coefficients at 0 Degree Celsius for 
!!   different plant parts for each PFT in routine stomate_constants.f90\n
!! - slope: calculated with a quadratic polynomial with the multi-annual mean air temperature 
!! (the constants are in routine stomate_constants.f90) as follows\n 
!!    \latexonly
!!      \input{resp3.tex} 
!!    \endlatexonly
!!   Where, maint_resp_slope1, maint_resp_slope2, maint_resp_slope3 are constant in stomate_constants.f90.
!!   Then coeff_maint is calculated as follows:\n
!!    \latexonly
!!      \input{resp4.tex} 
!!    \endlatexonly  
!! If the calculation result is negative, maintenance respiration coefficient will take the value 0.
!! Therefore the maintenance respiration will also be 0.\n
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): PFT maintenance respiration of different plant parts (::resp_maint_part_radia)
!!
!! REFERENCE(S)	:
!! McCree KJ. An equation for the respiration of white clover plants grown under controlled conditions. In: 
!! Setlik I, editor. Prediction and measurement of photosynthetic productivity. Wageningen, 
!! The Netherlands: Pudoc; 1970. p. 221-229.
!! Krinner G, Viovy N, de Noblet-Ducoudre N, Ogee J, Polcher J, Friedlingstein P,
!! Ciais P, Sitch S, Prentice I C (2005) A dynamic global vegetation model for studies
!! of the coupled atmosphere-biosphere system. Global Biogeochemical Cycles, 19, GB1015,
!! doi: 10.1029/2003GB002199.\n
!! Ruimy A., Dedieu G., Saugier B. (1996), TURC: A diagnostic model
!! of continental gross primary productivity and net primary productivity,
!! Global Biogeochemical Cycles, 10, 269-285.\n
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE maint_respiration ( npts, t2m,stempdiag,&
       rprof,biomass,resp_maint_part_radia)

!! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER, INTENT(in)                         :: npts      !! Domain size - number of grid cells (unitless)
    REAL, DIMENSION(npts), INTENT(in)           :: t2m       !! 2 meter air temperature - forcing variable (K)
    REAL, DIMENSION(npts,nslm), INTENT (in)     :: stempdiag !! Soil temperature of each soil layer (K)
    REAL, DIMENSION(npts,nvm), INTENT(in)       :: rprof     !! PFT root depth as calculated in stomate.f90 from parameter 
                                                                    !! humcste which is root profile for different PFTs 
                                                                    !! in slowproc.f90 (m)
    REAL, DIMENSION(npts,nvm,nparts,nelements),INTENT(in) :: biomass   !! PFT total biomass calculated in stomate_alloc.f90 
                                                                    !! @tex $(gC.m^{-2})$ @endtex

    !! 0.2 Output variables

    REAL, DIMENSION(npts,nvm,nparts), INTENT(out) :: resp_maint_part_radia !! PFT maintenance respiration of different plant
                                                                                  !! parts @tex $(gC.m^{-2}dt_sechiba^{-1} )$ @endtex

    !! 0.3 Modified variables
 
    !! 0.4 Local variables

    REAL, SAVE, ALLOCATABLE, DIMENSION(:)    :: z_soil       !! Variable to store depth of the different soil layers (m)
!$OMP THREADPRIVATE(z_soil)
    REAL, DIMENSION(npts,nvm)        :: t_root               !! PFT root temperature (convolution of root and soil 
                                                                    !! temperature profiles) (K)
    REAL, DIMENSION(npts,nvm,nparts) :: coeff_maint          !! PFT maintenance respiration coefficients of different 
                                                                    !! plant compartments at 0 deg C 
                                                                    !! @tex $(g.g^{-1}dt_sechiba^{-1})$ @endtex
    REAL, DIMENSION(npts)            :: rpc                  !! Scaling factor for integrating vertical soil 
                                                                    !! profiles (unitless)
    REAL, DIMENSION(npts,nparts)     :: t_maint_radia        !! Temperature which is pertinent for maintenance respiration, 
                                                                    !! which is air/root temperature for above/below-ground 
                                                                    !! compartments (K) 
    INTEGER                          :: i,j,k,l,m            !! Indeces (unitless)
    REAL, DIMENSION(npts,nvm)        :: gtemp                !! Temperature response of respiration in the
                                                                    !! Lloyd-Taylor Model (-)
    REAL, DIMENSION(npts)            :: cn                   !! CN ratio of a biomass pool ((gC)(gN)-1) 
    REAL                             :: limit_cn             !! Calculate limiting C/N ratio ((gC)(gN)-1)
    INTEGER                          :: ier                  !! Error handling 

!EXTERNALIZE: it is also used in stomate.f90
    REAL, PARAMETER                  :: SNC=0.004            !! Structural nitrogen  [gN g-1C] based on C:N of dead wood (White et al., 2000)
!EXTERNALIZE
                                                                    !! assuming carbon dry matter ratio of 0.5  
                                                                   

!_ ================================================================================================================================
    
    
    IF (printlev>=3) WRITE(numout,*) 'Entering respiration'
    
 !! 1. Initializations

    ! initialize output variable
    resp_maint_part_radia(:,:,:) = zero 
    
    IF ( firstcall_resp ) THEN

       !! 1.1. Soil levels (first call only)
       !       Set the depth of the different soil layers (number of layers: nslm) 
       !       previously calculated as variable diaglev in routines sechiba.f90 and slowproc.f90
       ALLOCATE(z_soil(0:nslm), stat=ier)
       IF ( ier /= 0 ) CALL ipslerr_p(3,'maint_respiration','Pb in allocate of z_soil','','')
       z_soil(0) = zero
       z_soil(1:nslm) = diaglev(1:nslm)

       firstcall_resp = .FALSE.
    ENDIF

    
    
    !! 1.2. Calculate root temperature
    !       Calculate root temperature as the convolution of root and soil temperature profiles
    DO j = 2,nvm ! Loop over # PFTs

       !! 1.2.1 Calculate rpc
       !  - rpc is an integration constant to make the integral over the root profile is equal 'one', 
       !    calculated as follows:\n
       !  \latexonly
       !    \input{resp1.tex} 
       !  \endlatexonly
       rpc(:) = un / ( un - EXP( -z_soil(nslm) / rprof(:,j) ) )

       !! 1.2.2 Calculate root temperature
       !        - Integrate root profile temperature (K) over soil layers (number of layers = nslm)
       !          with rpc and soil temperature (K) of each soil layer as follows:\n
       !        \latexonly
       !          \input{resp2.tex} 
       !        \endlatexonly
       !        Where, stempdiag is diagnostic temperature profile of soil (K)\n
       t_root(:,j) = zero

       DO l = 1, nslm ! Loop over # soil layers

          t_root(:,j) = &
               t_root(:,j) + stempdiag(:,l) * rpc(:) * &
               ( EXP( -z_soil(l-1)/rprof(:,j) ) - EXP( -z_soil(l)/rprof(:,j) ) )

       ENDDO ! Loop over # soil layers

    ENDDO ! Loop over # PFTs

 !! 2. Define maintenance respiration coefficients

    DO j = 2,nvm ! Loop over # PFTs

       !! 2.1 Temperature for maintenanace respiration
       !      Temperature which is used to calculate maintenance respiration for different plant compartments
       !      (above- and belowground)\n
       !      - for aboveground parts, we use 2-meter air temperature, t2m\n
       !      - for belowground parts, we use root temperature calculated in section 1.2 of this subroutine\n
       
       ! 2.1.1 Aboveground biomass
       t_maint_radia(:,ileaf) = t2m(:)
       t_maint_radia(:,isapabove) = t2m(:)
       t_maint_radia(:,ifruit) = t2m(:)

       ! 2.1.2 Belowground biomass
       t_maint_radia(:,isapbelow) = t_root(:,j)
       t_maint_radia(:,iroot) = t_root(:,j)

       !! 2.1.3 Heartwood biomass
       !        Heartwood does does not respire (coeff_maint_zero is set to zero)

       t_maint_radia(:,iheartbelow) = t_root(:,j)
       t_maint_radia(:,iheartabove) = t2m(:)

       t_maint_radia(:,ilabile) = t2m(:)
       !! 2.1.4 Reserve biomass
       !        Use aboveground temperature for trees and belowground temeperature for grasses
       IF ( is_tree(j) ) THEN
          t_maint_radia(:,icarbres) = t2m(:)
       ELSE
          t_maint_radia(:,icarbres) = t_root(:,j)
       ENDIF

       DO k = 1, nparts ! Loop over # plant parts

          ! Comment taken from DOFOCO: The original code refers to the Sitch et al 2003 as the source of the equations and parameters
          ! for modelling maintenance respiration. The equations below are consistent with the paper. The
          ! parameter setting for coeff_maint is in the range of 0.066 to 0.011 as reported in the paper but
          ! exact values are not given. Although, the principle of a climate correction for coeff_maint is
          ! mentioned in Sitch et al 2003, the reduction factors themselves were not given. As it appears
          ! now this block of code pretends much more knowledge then we actually have. Rather than using a
          ! baseline coeff_maint that is later corrected for the climate region, the parameter values for
          ! coeff_maint could be simply prescribed and made pft-specific.
          ! Further down in the code C/N ratios are used to constrain respiration. The C/N ratios were reset
          ! to the values presented in Sitch et al 2003 but still seem on the low side. If pft-specific
          ! values are to be used, changes in respiration could be compensated for by changing
          ! coeff_maint_init. Given Vicca et al 2012 (Ecology Letters) an NPP/GPP ratio of 0.5 is 'universal'
          ! for forests given a sufficient nutrient supply and strictly defining NPP as solely its biomass
          ! components (thus excluding VOC, exudation and subsidies to myccorrhizae as is the case in ORCHIDEE).
          ! Unless observation based values are available for coeff_maint_init, these values could be adjusted
          ! within the range of 0.066 to 0.011 to obtain an NPP/GPP of 0.5 in the absence of nutrient
          ! limitations.
           
          !DSG: completely agree, this piece of code is totally overparametrized
          !     and leads to higher CUE under nutrient stress thant under
          !     optimal nutrient conditions which contradicts observation
          
          ! LPJ respiration factors based on Sitch et al. 2003 - first part of the calculation
          IF ( k.EQ.iheartabove .OR. k.EQ.iheartbelow .OR. k.EQ.icarbres .OR. &
!DSGrevise_Ra
! PROBLEM:maintenance respiration is about a factor of 2 too high. 
! The seperation between maintenance and growth respiration is arbitrary (see Cannell et al., 2000), based on Ryan et al. 1991.
! In ORCHIDEE growth respiration is a fixed fraction of GPP not growth, which
! contradicts Ryan et al. 1991. Using growth would make no sense with nutrients
! as nutrient stress decreases growth AND decreases CUE.
! This whole respiration thing makes not much sense to me, but
! to fix this issue properly I lack the time (as this affects many
! parts of the model). 

! So what I did is 
! 1)to exclude sapwood N in the calculation of
! maintenance respiration. The coeff_maint which relates tissue N to respiration is anyway made up and there is a functional
! relationship between leaf and sapwood mass. There is further no evidence that
! sapwood N is good predictor of respiration while root and leaf N is. As growth
! respiration is not even derived from GPP this totally fine.
! According to Ryan et al 1991 maintenace respiration is ~80% of total Ra for
! grass and forests, which is in contradiction to the stuff cited in Ali et al.
! 2016 where they claim its 50:50 w/o giving much evidence. We with this fix somewhere in
! between 50 and 80% so I guess we are fine.

               (k.EQ.isapabove).OR.(k.EQ.isapbelow).OR.(k.EQ.ilabile)) THEN ! DSG:nitrogen labile should  per definition not respire 
             coeff_maint(:,j,k) = zero
          ELSE
             ! Use a  PFT-specific value - Values from OCN are used
             coeff_maint(:,j,k) = coeff_maint_init(j)* dt_sechiba/one_day
                     
          ENDIF

       ENDDO

    
 !! 3. Calculate maintenance respiration
 !! DSG: the use of tissue N to calculate maintenance respiration is a bit
 !       problematic, as respiration will go down when tissue N concentration declines
 !       due to nutrient stress. This contradicts observations which show an
 !       increase in respiration under nutrient stress. The whole autotrophic
 !       respiration formulation in ORCHIDEE should be revised at some point.

       DO k = 1, nparts ! Loop over # plant parts 
           
          ! LPJ respiration factors based on Sitch et al. 2003 - second part of the calculation
          ! Temperature response, LLoyd and Taylor, 1994. E0 = 308.56 comes from the paper of
          ! Lloyd and Taylor but was fitted for soil respiration which only partly consists of
          ! authotrophic (root) respiration.   
          WHERE(t_maint_radia(:,k)-ZeroCelsius-tmin_maint_resp(j) .GT.min_stomate)
             gtemp(:,j) = EXP(e0_maint_resp(j)*(1.0/(tref_maint_resp(j)-tmin_maint_resp(j))-1.0/(t_maint_radia(:,k)-ZeroCelsius-tmin_maint_resp(j))))
          ! No gtemp below -46.01 degrees Celsius
          ELSEWHERE
             gtemp(:,j) = 0.0
          ENDWHERE

          ! Following Ali et al. (2016) we account for structural N which does
          ! not contribute to respiration
          resp_maint_part_radia(:,j,k) = coeff_maint(:,j,k) * gtemp(:,j) * & 
                                         (biomass(:,j,k,initrogen) - SNC*biomass(:,j,k,icarbon))

       ENDDO ! Loop over # plant parts

    ENDDO

    IF (printlev>=4) WRITE(numout,*) 'Leaving respiration'

  END SUBROUTINE maint_respiration

END MODULE stomate_resp
