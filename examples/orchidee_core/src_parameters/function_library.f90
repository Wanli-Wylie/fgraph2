! =================================================================================================================================
! MODULE       : function_library
!
! CONTACT      : orchidee-help _at_ ipsl.jussieu.fr
!
! LICENCE      : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF       Collection of functions that are used throughout the ORCHIDEE code
!!
!!\n DESCRIPTION: Collection of modules to : (1) convert one variable into another i.e. basal area
!! to diameter, diamter to tree height, diameter to crown area, etc. (2) ...
!!
!! RECENT CHANGE(S): None
!!
!! REFERENCE(S)	:
!!
!! SVN          :
!! $HeadURL: svn://forge.ipsl.jussieu.fr/orchidee/branches/ORCHIDEE-DOFOCO/ORCHIDEE/src_stomate/stomate_prescribe.f90 $
!! $Date: 2013-01-04 16:50:56 +0100 (Fri, 04 Jan 2013) $
!! $Revision: 1126 $
!! \n
!_ ================================================================================================================================

MODULE function_library

  ! modules used:

!$  USE ioipsl_para
  USE pft_parameters
  USE constantes

  IMPLICIT NONE

  ! private & public routines

  PRIVATE
  PUBLIC calculate_c0_alloc, wood_to_ba_eff, wood_to_ba, &
       wood_to_qmheight, &
       wood_to_qmdia,  &
       wood_to_volume, &
       biomass_to_lai,  &
       check_biomass_sync

  CONTAINS

  
!! ================================================================================================================================
!! FUNCTION     : calculate_c0_alloc
!!
!>\BRIEF        Calculate the baseline root vs sapwood allocation
!!
!! DESCRIPTION : Calculates the baseline root vs sapwood allocation based on the 
!! parameters of the pipe model (hydraulic conductivities) and the
!! turnover of the different components              
!! 
!! RECENT CHANGE(S): None
!!
!! RETURN VALUE : ::c0_alloc (m)
!!
!! REFERENCE(S)	: 
!!
!! FLOWCHART	: None
!! \n
!_ ================================================================================================================================
   
  FUNCTION calculate_c0_alloc(pts, pft, tau_eff_root, tau_eff_sap)

 !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER                             :: pts               !! Pixel number (-)
    INTEGER                             :: pft               !! PFT number (-)
    REAL                                :: tau_eff_root      !! Effective longivety for leaves (days)
    REAL                                :: tau_eff_sap       !! Effective longivety for leaves (days)
    
    !! 0.2 Output variables
                
    REAL                                :: calculate_c0_alloc  !! quadratic mean height (m) 

    !! 0.3 Modified variables

    !! 0.4 Local variables
    REAL                                :: sapwood_density
    REAL                                :: qm_dia            !! quadratic mean diameter (m)

!_ ================================================================================================================================
 
    !! 1. Calculate c0_alloc
    IF ( is_tree(pft) ) THEN

       sapwood_density = deux * pipe_density(pft) / kilo_to_unit
       calculate_c0_alloc = sqrt(k_root(pft)/k_sap(pft)*tau_eff_sap/tau_eff_root*sapwood_density)

    ! Grasses and croplands   
    ELSE

       !+++CHECK+++
       ! Simply copied the same formulation as for trees but note
       ! that the sapwood in trees vs grasses and crops has a very
       ! meaning. In grasses and crops is structural carbon to ensure
       ! that the allocation works. In trees it really is the sapwood
       sapwood_density = deux * pipe_density(pft) / kilo_to_unit
       calculate_c0_alloc = sqrt(k_root(pft)/k_sap(pft)*tau_eff_sap/tau_eff_root*sapwood_density)
       !+++++++++++

    ENDIF ! is_tree(j)

 END FUNCTION calculate_c0_alloc





!! ================================================================================================================================
!! FUNCTION     : wood_to_ba_eff
!!
!>\BRIEF        Calculate effective basal area from woody biomass making use of allometric relationships
!!
!! DESCRIPTION :  Calculate basal area of an individual tree from the woody biomass of that tree making 
!! use of allometric relationships. Effective basal area accounts for both above and below ground carbon
!! and is the basis for the application of the rule of Deleuze and Dhote.
!! (i) woodmass = tree_ff * pipe_density*ba*height
!! (ii) height = pipe_tune2 * sqrt(4/pi*ba) ** pipe_tune_3  
!! 
!! RECENT CHANGE(S): None
!!
!! RETURN VALUE : effective basal area (m2 ind-1)
!!
!! REFERENCE(S)	: 
!!
!! FLOWCHART	: None
!! \n
!_ ================================================================================================================================
 
 FUNCTION wood_to_ba_eff(biomass_temp, pft)

 !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER                                          :: pft                  !! PFT number (-)
    REAL, DIMENSION(:)                             :: biomass_temp         !! Biomass of an individual tree within a circ 
                                                                                    !! class @tex $(m^{2} ind^{-1})$ @endtex 

    !! 0.2 Output variables
                
    REAL, DIMENSION(ncirc)                           :: wood_to_ba_eff       !! Effective basal area of an individual tree within a circ
                                                                                    !! class @tex $(m^{2} ind^{-1})$ @endtex

    !! 0.3 Modified variables

    !! 0.4 Local variables
    INTEGER                                          :: l                    !! Index
    REAL                                             :: woodmass_ind         !! Woodmass of an individual tree
                                                                                    !! @tex $(gC ind^{-1})$ @endtex

!_ ================================================================================================================================
 
 !! 1. Calculate basal area from woodmass

    IF ( is_tree(pft) ) THEN
       
       DO l = 1,ncirc
         
          ! Woodmass of an individual tree
          woodmass_ind = biomass_temp(isapabove) + biomass_temp(isapbelow) + &
             biomass_temp(iheartabove) + biomass_temp(iheartbelow)

          ! Basal area of that individual (m2 ind-1)  
          wood_to_ba_eff(l) = (pi/4*(woodmass_ind/(tree_ff(pft)*pipe_density(pft)*pipe_tune2(pft))) &
                  **(2./pipe_tune3(pft)))**(pipe_tune3(pft)/(pipe_tune3(pft)+2))
                 
       ENDDO

    ELSE

       WRITE(numout,*) 'pft ',pft
       CALL ipslerr_p (3,'wood_to_ba_eff', &
            'wood_to_ba_eff is not defined for this PFT.', &
            'See the output file for more details.','')

    ENDIF

 END FUNCTION wood_to_ba_eff



!! ================================================================================================================================
!! FUNCTION     : wood_to_ba
!!
!>\BRIEF        Calculate basal area from woody biomass making use of allometric relationships
!!
!! DESCRIPTION : Calculate basal area of an individual tree from the woody biomass of that tree making 
!! use of allometric relationships given below. Here basal area is defined in line with its classical 
!! forestry meaning.
!! (i) woodmass = tree_ff * pipe_density*ba*height
!! (ii) height = pipe_tune2 * sqrt(4/pi*ba) ** pipe_tune_3  
!! 
!! RECENT CHANGE(S): None
!!
!! RETURN VALUE : basal area (m2 ind-1)
!!
!! REFERENCE(S)	: 
!!
!! FLOWCHART	: None
!! \n
!_ ================================================================================================================================
 
 FUNCTION wood_to_ba(biomass_temp, pft)

 !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER                                          :: pft               !! PFT number (-)
    REAL, DIMENSION(:)                             :: biomass_temp      !! Biomass of an individual tree within a circ 
                                                                                 !! class @tex $(m^{2} ind^{-1})$ @endtex 

    !! 0.2 Output variables
                
    REAL                                            :: wood_to_ba        !! Basal area of an individual tree within a circ
                                                                                 !! class @tex $(m^{2} ind^{-1})$ @endtex

    !! 0.3 Modified variables

    !! 0.4 Local variables
    REAL                                             :: woodmass_ind      !! Woodmass of an individual tree
                                                                                 !! @tex $(gC ind^{-1})$ @endtex

!_ ================================================================================================================================
 
 !! 1. Calculate basal area from woodmass

    IF ( is_tree(pft) ) THEN
       
         
          ! Woodmass of an individual tree
          woodmass_ind = biomass_temp(iheartabove) + biomass_temp(isapabove)


          ! Basal area of that individual (m2 ind-1)  
          wood_to_ba = (pi/4*(woodmass_ind/(tree_ff(pft)*pipe_density(pft)*pipe_tune2(pft))) &
                  **(2./pipe_tune3(pft)))**(pipe_tune3(pft)/(pipe_tune3(pft)+2))
                 

    ELSE

       WRITE(numout,*) 'pft ',pft
       CALL ipslerr_p (3,'wood_to_ba', &
            'wood_to_ba is not defined for this PFT.', &
            'See the output file for more details.','')

    ENDIF

 END FUNCTION wood_to_ba





!! ================================================================================================================================
!! FUNCTION     : wood_to_qmheight
!!
!>\BRIEF        Calculate the quadratic mean height from the biomass
!!
!! DESCRIPTION : Calculates the quadratic mean height from the biomass
!! 
!! RECENT CHANGE(S): None
!!
!! RETURN VALUE : ::qm_height (m)
!!
!! REFERENCE(S)	: 
!!
!! FLOWCHART	: None
!! \n
!_ ================================================================================================================================
   
  FUNCTION wood_to_qmheight(biomass_temp, ind, pft)

 !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER                             :: pft               !! PFT number (-)
    REAL, DIMENSION(nparts)       :: biomass_temp      !! Biomass of the leaves @tex $(gC m^{-2})$ @endtex
    REAL, DIMENSION(ncirc)              :: ind               !! Number of individuals @tex $(m^{-2})$ @endtex


    !! 0.2 Output variables
                
    REAL                                :: wood_to_qmheight  !! quadratic mean height (m) 

    !! 0.3 Modified variables

    !! 0.4 Local variables
    REAL, DIMENSION(ncirc)              :: circ_class_ba     !! basal area for each circ_class @tex $(m^{2})$ @endtex
    REAL                                :: qm_dia            !! quadratic mean diameter (m)

!_ ================================================================================================================================
 
    !! 1. Calculate qm_height from the biomass
    IF ( is_tree(pft) ) THEN

       ! Basal area at the tree level (m2 tree-1)
       circ_class_ba(:) = wood_to_ba(biomass_temp(:),pft) 
       
       IF (SUM(ind(:)) .NE. zero) THEN

          qm_dia = SQRT( 4/pi*SUM(circ_class_ba(:)*ind(:))/SUM(ind(:)) )
          
       ELSE
          
          qm_dia = zero

       ENDIF
       
       wood_to_qmheight = pipe_tune2(pft)*(qm_dia**pipe_tune3(pft))
             

    ! Grasses and croplands   
    ELSE

       ! Calculate height as a function of the leaf and structural biomass. Use structural
       ! biomass to make sure that the grasslands have a roughness length during the winter
       ! If the biomass increases, vegetation height will increase as well. Divide by 
       ! ind(ipts,j) to obtain the height of the individual. biomass(ileaf) is in gC m-2 
       ! whereas qm is the height of the individual.
       IF (SUM(ind(:)) .NE. zero) THEN

          wood_to_qmheight = (biomass_temp(ileaf) + biomass_temp(isapabove)) / &
               SUM(ind(:)) * sla(pft) * lai_to_height(pft)

       ELSE

          wood_to_qmheight = zero

       ENDIF

    ENDIF ! is_tree(j)

 END FUNCTION wood_to_qmheight






!! ================================================================================================================================
!! FUNCTION     : wood_to_qmdia
!!
!>\BRIEF        Calculate the quadratic mean diameter from the biomass
!!
!! DESCRIPTION : Calculates the quadratic mean diameter from the aboveground biomss
!! 
!! RECENT CHANGE(S): None
!!
!! RETURN VALUE : ::qm_dia (m)
!!
!! REFERENCE(S)	: 
!!
!! FLOWCHART	: None
!! \n
!_ ================================================================================================================================
   
  FUNCTION wood_to_qmdia(biomass_temp, ind, pft)

 !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER                             :: pft               !! PFT number (-)
    REAL, DIMENSION(nparts)       :: biomass_temp      !! Biomass of the leaves @tex $(gC m^{-2})$ @endtex 
    REAL, DIMENSION(ncirc)              :: ind               !! Number of individuals @tex $(m^{-2})$ @endtex

    !! 0.2 Output variables
                
    REAL                                :: wood_to_qmdia     !! quadratic mean diameter (m) 

    !! 0.3 Modified variables

    !! 0.4 Local variables
    REAL, DIMENSION(ncirc)              :: circ_class_ba     !! basal area for each circ_class @tex $(m^{2})$ @endtex

!_ ================================================================================================================================
 
    !! 1. Calculate qm_dia from the biomass
    IF ( is_tree(pft) ) THEN

       ! Basal area at the tree level (m2 tree-1)
       circ_class_ba(:) = wood_to_ba(biomass_temp(:),pft) 
       
       IF (SUM(ind(:)) .NE. zero) THEN

          wood_to_qmdia = SQRT( 4/pi*SUM(circ_class_ba(:)*ind(:))/SUM(ind(:)) )
          
       ELSE
          
          wood_to_qmdia = zero

       ENDIF


    ! Grasses and croplands   
    ELSE

       wood_to_qmdia = zero   
       
    ENDIF ! is_tree(pft)

 END FUNCTION wood_to_qmdia


!! ================================================================================================================================
!! FUNCTION     : wood_to_volume
!!
!>\BRIEF        This allometric function computes volume as a function of 
!! biomass at stand scale. Volume \f$(m^3 m^{-2}) = f(biomass (gC m^{-2}))\f$
!!
!! DESCRIPTION : None
!!
!! RECENT CHANGE(S): None
!! 
!! RETURN VALUE : biomass_to_volume
!!
!! REFERENCE(S)	: See above, module description.
!!
!! FLOWCHART	: None
!! \n
!_ ================================================================================================================================

  FUNCTION wood_to_volume(biomass,pft,branch_ratio,inc_branches)

 !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    REAL, DIMENSION(:)   :: biomass              !! Stand biomass @tex $(gC m^{-2})$ @endtex  
    REAL                 :: branch_ratio         !! Branch ratio of sap and heartwood biomass
                                                        !! unitless
    INTEGER              :: pft                  !! Plant functional type (unitless)
    INTEGER              :: inc_branches         !! Include the branches in the volume calculation?
                                                        !! 0: exclude the branches from the volume calculation 
                                                        !! (thus correct the biomass for the branch ratio)
                                                        !! 1: include the branches in the volume calculation
                                                        !! (thus use all aboveground biomass)
                                    
    

    !! 0.2 Output variables

    REAL                 :: wood_to_volume       !! The volume of wood per square meter 
                                                        !!  @tex $(m^3 m^{-2})$ @endtex

    !! 0.3 Modified variables

    !! 0.4 Local variables
    
    REAL                 :: woody_biomass        !! Woody biomass at the stand level 
                                                        !! @tex $(gC m^{-2})$ @endtex

!_ ================================================================================================================================

 !! 1. Volume to biomass

    ! Woody biomass used in the calculation
    IF (inc_branches .EQ. 0) THEN

       woody_biomass=(biomass(isapabove)+biomass(iheartabove))*(un - branch_ratio)

    ELSEIF (inc_branches .EQ. 1) THEN

       woody_biomass=(biomass(isapabove)+biomass(iheartabove))

    ELSE

    ENDIF

    ! Wood volume expressed in m**3 / m**2
    wood_to_volume = woody_biomass/(pipe_density(pft))

  END FUNCTION wood_to_volume



!! ================================================================================================================================
!! FUNCTION     : biomass_to_lai
!!
!>\BRIEF        Calculate the LAI based on the leaf biomass
!!
!! DESCRIPTION : Calculates the LAI of a PFT/grid square based on the leaf biomass
!! 
!! RECENT CHANGE(S): None
!!
!! RETURN VALUE : ::LAI [m**2 m**{-2}]
!!
!! REFERENCE(S)	: 
!!
!! FLOWCHART	: None
!! \n
!_ ================================================================================================================================
   
  FUNCTION biomass_to_lai(leaf_biomass, pft)

 !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER                                          :: pft               !! PFT number (-)
    REAL                                             :: leaf_biomass      !! Biomass of the leaves
                                                                                 !! @tex $(gC m^{-2})$ @endtex 


    !! 0.2 Output variables
                
    REAL                                             :: biomass_to_lai    !! Leaf area index 
                                                                                 !! @tex $(m^{2} m^{-2})$ @endtex 

    !! 0.3 Modified variables

    !! 0.4 Local variables
    REAL                                             :: impose_lai         !! LAI read from run.def
!_ ================================================================================================================================
 
    !! 1. Calculate the LAI from the leaf biomass

    biomass_to_lai = leaf_biomass * sla(pft)
     
!!$    !+++++++++ TEMP ++++++++++
!!$    ! This is a perfect place to hack the code to make it run with
!!$    ! constant lai
!!$    WRITE(numout,*) 'WARNING ERROR: Using fake lai values for testing!'
!!$    biomass_to_lai=3.79052
!!$    !+++++++++++++++++++++++++

            !+++++++ TEMP ++++++++++
            ! This code is only used evaluation of the performance of the multi-layer energy budget. 
            ! To reduce the complexity of the tests we want to impose the LAI and its vertical distribution. 
            ! The solution is not very elegant but it works. 
            !  IF (ld_fake_height) THEN
            ! In order to imposed lai, we read the TOTAL_LAI from run.def
            !  CALL getin_p('TOTAL_LAI', impose_lai)
            ! This part of code reset the sla vale to match which alow modeled LAI equal to TOTAL LAI. 
            ! Althought this is ugly way to match the modeled LAI and impose LAI. 
            ! You probably need to go to your ORCHIDEE out put file to find out the suitable SLA value 
            ! and reset it agin in the run.def.
            ! So, we impose LAI & structure for a quick testing the performance of multilayer energy budget
            ! without changing the leaf_biomass.    
            !     IF ( leaf_biomass .GT. 0.0) THEN
            !        sla(pft)=impose_lai/leaf_biomass
            !          WRITE(numout,'(A,F20.8)') 'USE A FAKE SLA BASED ON imposed LAI/LEAFMASS=', sla(pft)
            !     ENDIF 
            !        biomass_to_lai=leaf_biomass*sla(pft) 
            !  ENDIF 
            !++++++++++++++++++++++++

     END FUNCTION biomass_to_lai





!! ================================================================================================================================
!! SUBROUTINE  : check_biomass_sync
!!
!>\BRIEF       
!!
!! DESCRIPTION : 
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S):
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
  SUBROUTINE check_biomass_sync ( check_point, npts, biomass, &
       circ_class_biomass, circ_class_n , ind, &
       lsync, bm_sync)

 !! 0. Variable and parameter description

    !! 0.1 Input variables
    INTEGER, INTENT(in)                                 :: npts                 !! Domain size (unitless)
    REAL, DIMENSION(:,:,:,:,:), INTENT(in)              :: circ_class_biomass   !! Biomass of the componets of the model  
                                                                                       !! tree within a circumference
                                                                                       !! class @tex $(gC ind^{-1})$ @endtex  
    REAL, DIMENSION(:,:,:), INTENT(in)                  :: circ_class_n         !! Number of individuals in each circ class
                                                                                       !! @tex $(m^{-2})$ @endtex
    REAL, DIMENSION(:,:,:,:), INTENT(in)                :: biomass              !! Stand level biomass 
                                                                                       !! @tex $(gC m^{-2})$ @endtex
    CHARACTER(*),INTENT(in)                                    :: check_point          !! A flag to indicate at which
                                                                                       !! point in the code we're doing
                                                                                       !! this check
    REAL, DIMENSION(:,:), INTENT(in)                    :: ind                  !! Density of individuals 
                                                                                       !! @tex $(m^{-2})$ @endtex 

    !! 0.2 Output variables
    LOGICAL,INTENT(out)                                        :: lsync
    REAL, DIMENSION(:,:,:), INTENT(out)                 :: bm_sync              !! The difference betweeen the
                                                                                       !! biomass in the circ_classes and
                                                                                       !! the total biomass 
                                                                                       !! @tex $(gC m^{-2})$ @endtex
    !! 0.3 Modified variables

    !! 0.4 Local variables
    INTEGER                                                    :: iele,ipts,ivm,ipar,icir
    REAL                                                :: total_circ_class_biomass
    REAL,DIMENSION(ncirc)                               :: tree_size
    LOGICAL                                                    :: lnegative
    
!_ ================================================================================================================================

    lsync=.TRUE.
    lnegative=.FALSE.

    bm_sync(:,:,:)=zero

    !++++++ TEMP ++++++
    ! We gain 5-10% speed by skipping this routine
       
    !++++++++++++

    ! Check to see if the biomass is not equal to the total biomass
    ! in circ_class_biomass anywhere.
    DO ipts=1,npts

       DO ivm=1,nvm

          ! Only woody PFTs have circumference classes therefore
          ! only woody PFTs need to be syncronized
          IF(.NOT. lbypass_cc)THEN
             IF(is_tree(ivm)) THEN
                tree_size(:)=zero
                DO icir=1,ncirc
                   tree_size(icir)=SUM(circ_class_biomass(ipts,ivm,icir,:,1))
                ENDDO
                DO icir=2,ncirc
                   IF(tree_size(icir) .LT. tree_size(icir-1)-min_stomate)THEN
                      WRITE(numout,*) 'ERROR: stopping in sync'
                      WRITE(numout,*) check_point
                      WRITE(numout,*) 'ipts,ivm: ',ipts,ivm
                      WRITE(numout,*) 'tree_size(icir), tree_size(icir-1), ',&
                           tree_size(icir), tree_size(icir-1), tree_size(icir) - tree_size(icir-1)  
                      WRITE(numout,*) 'icir, tree_size: ',icir, tree_size(:)
                      !+++ TEMP +++
                      !This would not STOP the ORCHIDEE beacause the mass balance is due to imposed LAI   
                     ! IF(ld_fake_height)  THEN 
                     !      CALL ipslerr_p (2,'check_biomass_sync', &
                     !      'The size of the trees in the circ class are not monotonically increasing!',&
                     !      'Look in the output file for more details.',&
                     !      '')
                     ! ELSE
                           CALL ipslerr_p (3,'check_biomass_sync', &
                           'The size of the trees in the circ class are not monotonically increasing!',&
                           'Look in the output file for more details.',&
                           '')
                     ! ENDIF
                      !++++++++++++    
                   ENDIF
                ENDDO
             ENDIF
          ENDIF
          
          DO iele=1,icarbon
                
             DO ipar=1,nparts
                
                total_circ_class_biomass=zero
                DO icir=1,ncirc
                   
                   total_circ_class_biomass=total_circ_class_biomass+&
                        circ_class_biomass(ipts,ivm,icir,ipar,iele)*circ_class_n(ipts,ivm,icir)

                   ! Check as well to see if our biomass is ever negative.
                   ! It really should not be.
                   IF(circ_class_biomass(ipts,ivm,icir,ipar,iele) .LT. -min_stomate)THEN

                      lnegative=.TRUE.
                      WRITE(numout,*) '!***********************************'
                      WRITE(numout,*) 'Error: Negative biomass component!'
                      WRITE(numout,*) 'Check point: ',TRIM(check_point)
                      WRITE(numout,*) 'circ_class_biomass(ipts,ivm,icir,ipar,iele) ',&
                           circ_class_biomass(ipts,ivm,icir,ipar,iele)
                      WRITE(numout,'(A,5I5)') 'ipts,ivm,icir,ipar,iele',ipts,ivm,icir,ipar,iele
                      WRITE(numout,*) '!***********************************'

                   ENDIF
                ENDDO

                IF(ABS(biomass(ipts,ivm,ipar,iele) -  &
                     total_circ_class_biomass) .GT. sync_threshold)THEN

                   WRITE(numout,*) '!***********************************'
                   WRITE(numout,*) 'Biomass and circ_class_biomass are not equal!'
                   WRITE(numout,*) 'Check point: ',TRIM(check_point)
                   WRITE(numout,100) 'biomass(ipts,ivm,ipar,iele) ',&
                        biomass(ipts,ivm,ipar,iele)
                   WRITE(numout,100) 'total_circ_class_biomass ',&
                        total_circ_class_biomass
                   WRITE(numout,100) 'Difference: ',&
                        ABS(biomass(ipts,ivm,ipar,iele) - total_circ_class_biomass)
                   WRITE(numout,*) 'ipts,ivm,ipar,iele',ipts,ivm,ipar,iele
                   WRITE(numout,*) '!***********************************'
100                FORMAT(A,E20.10)
                   lsync=.FALSE.

                ENDIF

             ENDDO

             ! we are not going to save the biomass for every component right now,
             ! just the total
             bm_sync(ipts,ivm,iele)=zero

             DO ipar=1,nparts
                   

                DO icir=1,ncirc

                   bm_sync(ipts,ivm,iele)=bm_sync(ipts,ivm,iele)+&
                        circ_class_biomass(ipts,ivm,icir,ipar,iele)*circ_class_n(ipts,ivm,icir)
                ENDDO ! ncirc

             ENDDO ! nparts

             bm_sync(ipts,ivm,iele)=ABS(bm_sync(ipts,ivm,iele)-&
                  SUM(biomass(ipts,ivm,:,iele)))

          ENDDO ! nelements



       ENDDO ! loop over PFTs

    ENDDO ! loop over points
 
    !---TEMP---
    IF(ld_biomass)THEN
       WRITE(numout,*) 'Check point: ',TRIM(check_point)
       WRITE(numout,*) 'test_pft, test_grid: ',test_pft,test_grid
       WRITE(numout,*) 'biomass (ileaf), ', biomass(test_grid,test_pft,ileaf,icarbon)
       WRITE(numout,*) 'biomass (iwood), ', biomass(test_grid,test_pft,isapabove,icarbon) + &
           biomass(test_grid,test_pft,isapbelow,icarbon) + biomass(test_grid,test_pft,iheartabove,icarbon) + &
           biomass(test_grid,test_pft,iheartbelow,icarbon)
       WRITE(numout,*) 'biomass (iroot), ', biomass(test_grid,test_pft,iroot,icarbon)
       WRITE(numout,'(A,20F14.6)') 'biomassHHH, ',biomass(test_grid,test_pft,:,icarbon)
       DO icir=1,ncirc
          WRITE(numout,'(A,I1,20F14.6)') 'ccbiomass',icir,circ_class_biomass(test_grid,test_pft,icir,:,icarbon)
       ENDDO
       WRITE(numout,*) 'circ_class_biomass, ',&
            SUM (SUM(circ_class_biomass(test_grid,test_pft,:,:,icarbon),2) * &
            circ_class_n(test_grid,test_pft,:))
       WRITE(numout,*) 'circ_class_n, ', SUM(circ_class_n(test_grid,test_pft,:))
       WRITE(numout,*) 'circ_class_n(:), ', circ_class_n(test_grid,test_pft,:)
       WRITE(numout,*) 'ind, ', ind(test_grid,test_pft)
    ENDIF

!!$    !----------

    IF(.NOT. lsync) THEN
       WRITE(numout,*) 'ERROR: stopping in sync #2'
       WRITE(numout,*) 'Stopping'
       CALL ipslerr_p (3,'check_biomass_sync', &
            'circ_class_biomass*circ_class_n is not equal to the total biomass',&
            'Look in the output file for more details.',&
            '')

    ENDIF
    IF(lnegative) THEN
       WRITE(numout,*) 'ERROR: negative biomass'
       WRITE(numout,*) 'Stopping'
       CALL ipslerr_p (3,'check_biomass_sync', &
            'One of the biomass pools is negative!',&
            'Look in the output file for more details.',&
            '')
    ENDIF

  END SUBROUTINE check_biomass_sync

END MODULE function_library

