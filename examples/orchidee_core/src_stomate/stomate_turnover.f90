! =================================================================================================================================
! MODULE       : stomate_turnover.f90
!
! CONTACT      : orchidee-help _at_ listes.ipsl.fr
!
! LICENCE      : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF        This module manages the end of the growing season and calculates herbivory and turnover of leaves, fruits, fine roots.
!! 
!!\n DESCRIPTION: This subroutine calculates leaf senescence due to climatic conditions or as a 
!! function of leaf age and new LAI, and subsequent turnover of the different plant biomass compartments (sections 1 to 6), 
!! herbivory (section 7), fruit turnover for trees (section 8) and sapwood conversion (section 9). 
!!
!! RECENT CHANGE(S): None
!!
!! SVN          :
!! $HeadURL: svn://forge.ipsl.jussieu.fr/orchidee/branches/ORCHIDEE-CN-P/ORCHIDEE/src_stomate/stomate_turnover.f90 $
!! $Date: 2020-08-11 09:58:42 +0200 (二, 2020-08-11) $
!! $Revision: 6859 $
!! \n
!_ ================================================================================================================================

MODULE stomate_turnover

  ! modules used:
  USE xios_orchidee
  USE ioipsl_para
  USE stomate_data
  USE constantes
  USE pft_parameters

  IMPLICIT NONE

  ! private & public routines

  PRIVATE
  PUBLIC turn, turn_clear

  LOGICAL, SAVE                          :: firstcall_turnover = .TRUE.           !! first call (true/false)
!$OMP THREADPRIVATE(firstcall_turnover)

CONTAINS


!! ================================================================================================================================
!! SUBROUTINE   : turn_clear
!!
!>\BRIEF        Set flag ::firstcall_turnover to .TRUE., and therefore activate section 1 
!!              of subroutine turn which writes a message to the output.
!!                
!_ ================================================================================================================================

  SUBROUTINE turn_clear
    firstcall_turnover=.TRUE.
  END SUBROUTINE turn_clear


!! ================================================================================================================================
!! SUBROUTINE    : turn
!!
!>\BRIEF         Calculate turnover of leaves, roots, fruits and sapwood due to aging or climatic 
!!               induced senescence. Calculate herbivory.
!!
!! DESCRIPTION : This subroutine determines the turnover of leaves and fine roots (and stems for grasses)
!! and simulates following processes:
!! 1. Mean leaf age is calculated from leaf ages of separate leaf age classes. Should actually 
!!    be recalculated at the end of the routine, but it does not change too fast. The mean leaf 
!!    age is calculated using the following equation:
!!    \latexonly
!!    \input{turnover_lma_update_eqn1.tex}
!!    \endlatexonly
!!    \n 
!! 2. Meteorological senescence: the detection of the end of the growing season and shedding 
!!    of leaves, fruits and fine roots due to unfavourable meteorological conditions.
!!    The model distinguishes three different types of "climatic" leaf senescence, that do not 
!!    change the age structure: sensitivity to cold temperatures, to lack of water, or both. 
!!    If meteorological conditions are fulfilled, a flag ::senescence is set to TRUE. Note
!!    that evergreen species do not experience climatic senescence.
!!    Climatic senescence is triggered by sensitivity to cold temperatures where the critical 
!!    temperature for senescence is calculated using the following equation:
!!    \latexonly
!!    \input{turnover_temp_crit_eqn2.tex}
!!    \endlatexonly
!!    \n
!!    Climatic senescence is triggered by sensitivity to lack of water availability where the 
!!    moisture availability critical level is calculated using the following equation:
!!    \latexonly
!!    \input{turnover_moist_crit_eqn3.tex}
!!    \endlatexonly
!!    \n
!!    Climatic senescence is triggered by sensitivity to temperature or to lack of water where
!!    critical temperature and moisture availability are calculated as above.\n
!!    Trees in climatic senescence lose their fine roots at the same rate as they lose their leaves. 
!!    The rate of biomass loss of both fine roots and leaves is presribed through the equation:
!!    \latexonly
!!    \input{turnover_clim_senes_biomass_eqn4.tex}
!!    \endlatexonly
!!    \n
!!    with ::leaffall(j) a PFT-dependent time constant which is given in 
!!    ::stomate_constants. In grasses, leaf senescence is extended to the whole plant 
!!    (all carbon pools) except to its carbohydrate reserve.    
!! 3. Senescence due to aging: the loss of leaves, fruits and  biomass due to aging
!!    At a certain age, leaves fall off, even if the climate would allow a green plant
!!    all year round. Even if the meteorological conditions are favorable for leaf maintenance,
!!    plants, and in particular, evergreen trees, have to renew their leaves simply because the 
!!    old leaves become inefficient. Roots, fruits (and stems for grasses) follow leaves. 
!!    The ??senescence?? rate varies with leaf age. Note that plant is not declared senescent 
!!    in this case (wchich is important for allocation: if the plant loses leaves because of 
!!    their age, it can renew them). The leaf turnover rate due to aging of leaves is calculated
!!    using the following equation:
!!    \latexonly
!!    \input{turnover_age_senes_biomass_eqn5.tex}
!!    \endlatexonly
!!    \n
!!    Drop all leaves if there is a very low leaf mass during senescence. After this, the biomass 
!!    of different carbon pools both for trees and grasses is set to zero and the mean leaf age 
!!    is reset to zero. Finally, the leaf fraction and leaf age of the different leaf age classes 
!!    is set to zero. For deciduous trees: next to leaves, also fruits and fine roots are dropped.
!!    For grasses: all aboveground carbon pools, except the carbohydrate reserves are affected:
!! 4. Update the leaf biomass, leaf age class fraction and the LAI
!!    Older leaves will fall more frequently than younger leaves and therefore the leaf age 
!!    distribution needs to be recalculated after turnover. The fraction of biomass in each 
!!    leaf class is updated using the following equation:
!!    \latexonly
!!    \input{turnover_update_LeafAgeDistribution_eqn6.tex}
!!    \endlatexonly
!!    \n
!! 5. Simulate herbivory activity and update leaf and fruits biomass. Herbivore activity 
!!    affects the biomass of leaves and fruits as well as stalks (only for grasses).
!!    However, herbivores do not modify leaf age structure.
!! 6. Calculates fruit turnover for trees. Trees simply lose their fruits with a time 
!!    constant ::tau_fruit(j), that is set to 90 days for all PFTs in ::stomate_constants 
!! 7. Convert sapwood to heartwood for trees and update heart and softwood above and 
!!    belowground biomass. Sapwood biomass is converted into heartwood biomass 
!!    with a time constant tau ::tau_sap(j) of 1 year. Note that this biomass conversion 
!!    is not added to "turnover" as the biomass is not lost. For the updated heartwood, 
!!    the sum of new heartwood above and new heartwood below after converting sapwood to 
!!    heartwood, is saved as ::hw_new(:). Creation of new heartwood decreases the age of 
!!    the plant ??carbon?? with a factor that is determined by: old heartwood ::hw_old(:) 
!!    divided by the new heartwood ::hw_new(:)
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLES: ::Biomass of leaves, fruits, fine roots and sapwood above (latter for grasses only), 
!!                        ::Update LAI, ::Update leaf age distribution with new leaf age class fraction 
!!
!! REFERENCE(S) : 
!! - Krinner, G., N. Viovy, N. de Noblet-Ducoudre, J. Ogee, J. Polcher, P. 
!! Friedlingstein, P. Ciais, S. Sitch and I.C. Prentice (2005), A dynamic global
!! vegetation model for studies of the coupled atmosphere-biosphere system, Global
!! Biogeochemical Cycles, 19, doi:10.1029/2003GB002199. 
!! - McNaughton, S. J., M. Oesterheld, D. A. Frank and K. J. Williams (1989), 
!! Ecosystem-level patterns of primary productivity and herbivory in terrestrial habitats, 
!! Nature, 341, 142-144, 1989. 
!! - Sitch, S., C. Huntingford, N. Gedney, P. E. Levy, M. Lomas, S. L. Piao, , Betts, R., Ciais, P., Cox, P., 
!! Friedlingstein, P., Jones, C. D., Prentice, I. C. and F. I. Woodward : Evaluation of the terrestrial carbon  
!! cycle, future plant geography and climate-carbon cycle feedbacks using 5 dynamic global vegetation 
!! models (dgvms), Global Change Biology, 14(9), 2015–2039, 2008. 
!!
!! FLOWCHART    : 
!! \latexonly
!! \includegraphics[scale=0.5]{turnover_flowchart_1.png}
!! \includegraphics[scale=0.5]{turnover_flowchart_2.png}
!! \endlatexonly
!! \n
!_ ================================================================================================================================

  SUBROUTINE turn (npts, dt, PFTpresent, &
       herbivores, &
       maxmoiavail_lastyear, minmoiavail_lastyear, &
       moiavail_week, moiavail_month, t2m_longterm, t2m_month, t2m_week, veget_cov_max, &
       gdd_from_growthinit, leaf_age, leaf_frac, age, lai, biomass, &
       turnover, recycling, senescence,turnover_time, &
       sla_calc,GRM_enable_grazing,when_growthinit,days_senescence)

    !! 0. Variable and parameter declaration

    !! 0.1 Input variables 

    INTEGER, INTENT(in)                                 :: npts                 !! Domain size - number of grid cells 
                                                                                       !! (unitless) 
    REAL, INTENT(in)                                    :: dt                   !! time step (dt_days)
    LOGICAL, DIMENSION(npts,nvm), INTENT(in)                   :: PFTpresent           !! PFT exists (true/false)
    REAL, DIMENSION(npts,nvm), INTENT(in)               :: herbivores           !! time constant of probability of a leaf to 
                                                                                       !! be eaten by a herbivore (days) 
    REAL, DIMENSION(npts,nvm), INTENT(in)               :: maxmoiavail_lastyear !! last year's maximum moisture availability 
                                                                                       !! (0-1, unitless)
    REAL, DIMENSION(npts,nvm), INTENT(in)               :: minmoiavail_lastyear !! last year's minimum moisture availability 
                                                                                       !! (0-1, unitless)
    REAL, DIMENSION(npts,nvm), INTENT(in)               :: moiavail_week        !! "weekly" moisture availability 
                                                                                       !! (0-1, unitless)
    REAL, DIMENSION(npts,nvm), INTENT(in)               :: moiavail_month       !! "monthly" moisture availability 
                                                                                       !! (0-1, unitless)
    REAL, DIMENSION(npts), INTENT(in)                   :: t2m_longterm         !! "long term" 2 meter reference 
                                                                                       !! temperatures (K) 
    REAL, DIMENSION(npts), INTENT(in)                   :: t2m_month            !! "monthly" 2-meter temperatures (K)
    REAL, DIMENSION(npts), INTENT(in)                   :: t2m_week             !! "weekly" 2 meter temperatures (K)
    REAL, DIMENSION(npts,nvm), INTENT(in)               :: veget_cov_max        !! "maximal" coverage fraction of a PFT (LAI 
                                                                                       !! -> infinity) on ground (unitless) 
    REAL, DIMENSION(npts,nvm), INTENT(in)               :: gdd_from_growthinit  !! gdd senescence for crop
    REAL, DIMENSION(npts,nvm), INTENT(in)               :: sla_calc             !! Leaf age-related SLA
    LOGICAL, INTENT(in)                                        :: GRM_enable_grazing   !! whether activate GRM
    REAL, DIMENSION(npts,nvm), INTENT(in)               :: when_growthinit      !! days after growth init

    !! 0.2 Output variables

    REAL, DIMENSION(npts,nvm,nparts,nelements), INTENT(out) :: turnover         !! Turnover @tex ($gC m^{-2}$) @endtex (corrected for recycled nutrients)
    REAL, DIMENSION(npts,nvm,nparts,nelements), INTENT(out) :: recycling        !! Recycled nutrients 

    LOGICAL, DIMENSION(npts,nvm), INTENT(out)                      :: senescence       !! is the plant senescent? (true/false) 
                                                                                       !! (interesting only for deciduous trees: 
                                                                                       !! carbohydrate reserve) 
    !! 0.3 Modified variables

    REAL, DIMENSION(npts,nvm,nparts,nleafages), INTENT(inout) :: leaf_age       !! age of the leaves (days)
    REAL, DIMENSION(npts,nvm,nparts,nleafages), INTENT(inout) :: leaf_frac      !! fraction of leaves in leaf age class 
                                                                                       !! (0-1, unitless)
    REAL, DIMENSION(npts,nvm), INTENT(inout)                  :: age            !! age (years)
    REAL, DIMENSION(npts,nvm), INTENT(inout)                  :: lai            !! leaf area index @tex ($m^2 m^{-2}$) 
                                                                                       !! @endtex 
    REAL, DIMENSION(npts,nvm,nparts,nelements), INTENT(inout) :: biomass        !! biomass @tex ($gC m^{-2}$) @endtex
    REAL, DIMENSION(npts,nvm), INTENT(inout)                  :: turnover_time  !! turnover_time of grasses (days)
    REAL, DIMENSION(npts,nvm), INTENT(inout)                  :: days_senescence!! accumulated days of senescence

    !! 0.4 Local  variables

    REAL, DIMENSION(npts,nvm,nparts)                    :: leaf_meanage         !! mean age of the leaves (days)
    REAL, DIMENSION(npts,nvm,nelements,2)                :: bm_sap_2_heartw      !! Heartwood growth (g m-2)
    CHARACTER(LEN=2), DIMENSION(nelements)            :: element_str            !! string suffix indicating element 
    REAL, DIMENSION(npts)                               :: dturnover            !! Intermediate variable for turnover ??
                                                                                       !! @tex ($gC m^{-2}$) @endtex 
    REAL, DIMENSION(npts)                               :: moiavail_crit        !! critical moisture availability, function 
                                                                                       !! of last year's moisture availability 
                                                                                       !! (0-1, unitless)
    REAL, DIMENSION(npts)                               :: tl                   !! long term annual mean temperature, (C)
    REAL, DIMENSION(npts)                               :: t_crit               !! critical senescence temperature, function 
                                                                                       !! of long term annual temperature (K) 
    LOGICAL, DIMENSION(npts)                                   :: shed_rest            !! shed the remaining leaves? (true/false)
    LOGICAL, DIMENSION(npts)                                   :: shed_minimum         !! shed the remaining biomass? (true/false)
    REAL, DIMENSION(npts)                               :: tot_bio              !! total biomass for detecting marginal biomass
    REAL, DIMENSION(npts)                               :: sapconv              !! Sapwood conversion @tex ($gC m^{-2}$) 
                                                                                       !! @endtex 
    REAL, DIMENSION(npts)                               :: hw_old               !! old heartwood mass @tex ($gC m^{-2}$) 
                                                                                       !! @endtex 
    REAL, DIMENSION(npts)                               :: hw_new               !! new heartwood mass @tex ($gC m^{-2}$) 
                                                                                       !! @endtex 
    REAL, DIMENSION(npts,nparts)                        :: lm_old               !! old leaf mass @tex ($gC m^{-2}$) @endtex
    REAL, DIMENSION(npts,nparts,nleafages)              :: delta_lm             !! leaf mass change for each age class @tex 
                                                                                       !! ($gC m^{-2}$) @endtex 
    REAL, DIMENSION(npts,nparts)                        :: turnover_rate        !! turnover rate (unitless) 
    REAL, DIMENSION(npts,nvm)                           :: leaf_age_crit        !! critical leaf age (days)
    REAL, DIMENSION(npts,nvm)                           :: root_age_crit        !! critical root age (days)
    REAL, DIMENSION(npts,nvm)                           :: new_turnover_time    !! instantaneous turnover time (days)
    INTEGER                                             :: j,m,k,ipar           !! Index (unitless)

!_ ================================================================================================================================

    IF (printlev>=3) WRITE(numout,*) 'Entering turnover'

    !! 1. first call - output messages

    IF ( firstcall_turnover ) THEN

       IF (printlev >=2 ) THEN
          WRITE(numout,*) 'turnover:'
          WRITE(numout,*) ' > minimum mean leaf age for senescence (days) (::min_leaf_age_for_senescence) : '&
               ,min_leaf_age_for_senescence
       END IF
       firstcall_turnover = .FALSE.

    ENDIF

    !! 2. Initializations 

    !! 2.1 set output to zero
    turnover(:,:,:,:)      = zero
    recycling(:,:,:,:)     = zero
    new_turnover_time(:,:) = zero
    senescence(:,:)        = .FALSE.

    !! 2.2 Recalculate mean leaf age
    !      Mean leaf age is recalculated from leaf ages of separate leaf age classes. 
    !      Should actually be recalculated at the 
    !      end of this routine, but it does not change too fast.
    !      The mean leaf age is calculated using the following equation:
    !      \latexonly
    !      \input{turnover_lma_update_eqn1.tex}
    !      \endlatexonly
    !      \n
    leaf_meanage(:,:,:) = zero

    DO m = 1, nleafages
    !WRITE(numout,*) 'invalid at', m
       leaf_meanage(:,:,:) = leaf_meanage(:,:,:) + leaf_age(:,:,:,m) * leaf_frac(:,:,:,m)
    ENDDO

    !! 3. Climatic senescence

    !     Three different types of "climatic" leaf senescence,
    !     that do not change the age structure. 
    DO j = 2,nvm ! Loop over # PFTs

       !! 3.1 Determine if there is climatic senescence. 
       !      The climatic senescence can be of three types:
       !      sensitivity to cold temperatures, to lack of water, or both. If meteorological conditions are 
       !      fulfilled, a flag senescence is set to TRUE.
       !      Evergreen species do not experience climatic senescence.

       SELECT CASE ( senescence_type(j) )


       CASE ('crop' )!for crop senescence is based on a GDD criterium as in crop models
          WHERE ( ( biomass(:,j,ileaf,icarbon) .GT. zero ) .AND. &
               (( leaf_meanage(:,j,ileaf) .GT. min_leaf_age_for_senescence(j) ) .AND. &
               ( gdd_from_growthinit(:,j) .GT.  gdd_senescence(j)) .OR. &
               !JC MOD041 avoid no senescence for crops, force to senescence after 280 days
               ! in some cold climate GDD and leaf age crit will nevel achieved
               (when_growthinit(:,j) .GT. min_growthinit_time-20)))
             senescence(:,j) = .TRUE.
          ENDWHERE

       CASE ( 'cold' )

          !! 3.1.1 Summergreen species: Climatic senescence is triggered by sensitivity to cold temperatures
          !        Climatic senescence is triggered by sensitivity to cold temperatures as follows: 
          !        If biomass is large enough (i.e. when it is greater than zero), 
          !        AND (i.e. when leaf mean age is above a certain PFT-dependent treshold ::min_leaf_age_for_senescence,
          !        which is given in ::stomate_constants),      
          !        AND the monthly temperature is low enough (i.e. when monthly temperature ::t2m_month(:) is below a critical 
          !        temperature ::t_crit(:), which is calculated in this module),
          !        AND the temperature tendency is negative (i.e. when weekly temperatures ::t2m_week(:) are lower than monthly 
          !        temperatures ::t2m_month(:))
          !        If these conditions are met, senescence is set to TRUE.
          !
          !        The critical temperature for senescence is calculated using the following equation:
          !        \latexonly
          !        \input{turnover_temp_crit_eqn2.tex}
          !        \endlatexonly
          !        \n
          !
          ! Critical temperature for senescence may depend on long term annual mean temperature
          tl(:) = t2m_longterm(:) - ZeroCelsius
          t_crit(:) = ZeroCelsius + senescence_temp(j,1) + &
               tl(:) * senescence_temp(j,2) + &
               tl(:)*tl(:) * senescence_temp(j,3)

          WHERE ( ( biomass(:,j,ileaf,icarbon) .GT. zero ) .AND. &
               ( leaf_meanage(:,j,ileaf) .GT. min_leaf_age_for_senescence(j) ) .AND. &
               ( t2m_month(:) .LT. t_crit(:) ) .AND. ( t2m_week(:) .LT. t2m_month(:) ) )

             senescence(:,j) = .TRUE.

          ENDWHERE

       CASE ( 'dry' )

          !! 3.1.2 Raingreen species: Climatic senescence is triggered by sensitivity to lack of water availability 
          !        Climatic senescence is triggered by sensitivity to lack of water availability as follows:  
          !        If biomass is large enough (i.e. when it is greater than zero), 
          !        AND (i.e. when leaf mean age is above a certain PFT-dependent treshold ::min_leaf_age_for_senescence,
          !        which is given in ::stomate_constants),      
          !        AND the moisture availability drops below a critical level (i.e. when weekly moisture availability 
          !        ::moiavail_week(:,j) is below a critical moisture availability ::moiavail_crit(:),
          !        which is calculated in this module), 
          !        If these conditions are met, senescence is set to TRUE.
          !
          !        The moisture availability critical level is calculated using the following equation:
          !        \latexonly
          !        \input{turnover_moist_crit_eqn3.tex}
          !        \endlatexonly
          !        \n
          moiavail_crit(:) = &
               MIN( MAX( minmoiavail_lastyear(:,j) + hum_frac(j) * &
               ( maxmoiavail_lastyear(:,j) - minmoiavail_lastyear(:,j) ), &
               senescence_hum(j) ), &
               nosenescence_hum(j) )

          WHERE ( ( biomass(:,j,ileaf,icarbon) .GT. zero ) .AND. &
               ( leaf_meanage(:,j,ileaf) .GT. min_leaf_age_for_senescence(j) ) .AND. &
               ( moiavail_week(:,j) .LT. moiavail_crit(:) ) )

             senescence(:,j) = .TRUE.

          ENDWHERE

       CASE ( 'mixed' )

          !! 3.1.3 Mixed criterion: Climatic senescence is triggered by sensitivity to temperature or to lack of water  
          !        Climatic senescence is triggered by sensitivity to temperature or to lack of water availability as follows:
          !        If biomass is large enough (i.e. when it is greater than zero), 
          !        AND (i.e. when leaf mean age is above a certain PFT-dependent treshold ::min_leaf_age_for_senescence,
          !        which is given in ::stomate_constants),      
          !        AND the moisture availability drops below a critical level (i.e. when weekly moisture availability 
          !        ::moiavail_week(:,j) is below a critical moisture availability ::moiavail_crit(:), calculated in this module), 
          !        OR 
          !        the monthly temperature is low enough (i.e. when monthly temperature ::t2m_month(:) is below a critical 
          !        temperature ::t_crit(:), calculated in this module),
          !        AND the temperature tendency is negative (i.e. when weekly temperatures ::t2m_week(:) are lower than 
          !        monthly temperatures ::t2m_month(:)).
          !        If these conditions are met, senescence is set to TRUE.
          moiavail_crit(:) = &
               MIN( MAX( minmoiavail_lastyear(:,j) + hum_frac(j) * &
               (maxmoiavail_lastyear(:,j) - minmoiavail_lastyear(:,j) ), &
               senescence_hum(j) ), &
               nosenescence_hum(j) )

          tl(:) = t2m_longterm(:) - ZeroCelsius
          t_crit(:) = ZeroCelsius + senescence_temp(j,1) + &
               tl(:) * senescence_temp(j,2) + &
               tl(:)*tl(:) * senescence_temp(j,3)

          IF ( is_tree(j) ) THEN
             ! critical temperature for senescence may depend on long term annual mean temperature
             WHERE ( ( biomass(:,j,ileaf,icarbon) .GT. zero ) .AND. &
                  ( leaf_meanage(:,j,ileaf) .GT. min_leaf_age_for_senescence(j) ) .AND. &
                  ( ( moiavail_week(:,j) .LT. moiavail_crit(:) ) .OR. &
                  ( ( t2m_month(:) .LT. t_crit(:) ) .AND. ( t2m_week(:) .LT. t2m_month(:) ) ) ) )
                senescence(:,j) = .TRUE.
             ENDWHERE
          ELSE

            turnover_time(:,j) = max_turnover_time(j)

            WHERE ( ( biomass(:,j,ileaf,icarbon) .GT. zero ) .AND. &
                 ( leaf_meanage(:,j,ileaf) .GT. min_leaf_age_for_senescence(j) ) .AND. &
                 ( ( moiavail_week(:,j) .LT. moiavail_crit(:) )))
                turnover_time(:,j) = max_turnover_time(j) * &
                     (1.-   (1.- (moiavail_week(:,j)/  moiavail_crit(:)))**2)
                senescence(:,j) = .TRUE.
            ENDWHERE

            WHERE ( turnover_time(:,j) .LT. min_turnover_time(j) ) 
                turnover_time(:,j) = min_turnover_time(j) 
            ENDWHERE  

            WHERE ((( biomass(:,j,ileaf,icarbon) .GT. zero ) .AND. &
                ( leaf_meanage(:,j,ileaf) .GT. min_leaf_age_for_senescence(j) ) .AND. &
                ((t2m_month(:) .LT. t_crit(:)) .AND. (lai(:,j) .GT. lai_max(j)/4.) .OR. &
                (t2m_month(:) .LT. ZeroCelsius)) .AND. ( t2m_week(:) .LT. t2m_month(:) )))
               turnover_time(:,j)= leaffall(j)
               senescence(:,j) = .TRUE. 
            ENDWHERE
            
         ENDIF


       !! Evergreen species do not experience climatic senescence
       CASE ( 'none' )

          
       !! In case no climatic senescence type is recognized.
       CASE default

          WRITE(numout,*) '  turnover: don''t know how to treat this PFT.'
          WRITE(numout,*) '  number (::j) : ',j
          WRITE(numout,*) '  senescence type (::senescence_type(j)) : ',senescence_type(j)

          CALL ipslerr_p(3,"turn","Dont know how to treat this PFT.","","")

       END SELECT

       !JC MOD050 suppress resp_maint of grasses when detecting dormancy
       ! calculate accumulated days where grasses senescence happened
       ! due to low soil moisture or low temperature
       WHERE(senescence(:,j) .EQV. .TRUE.)
         days_senescence(:,j) = days_senescence(:,j) + un
       ELSEWHERE
         days_senescence(:,j) = days_senescence(:,j) - un
       ENDWHERE
       WHERE(days_senescence(:,j) .LT. min_stomate)
         days_senescence(:,j) = zero
       ENDWHERE

       !! 3.2 Drop leaves and roots, plus stems and fruits for grasses
       DO m=1,nelements

          IF ( is_tree(j) ) THEN

             !! 3.2.1 Trees in climatic senescence lose their fine roots at the same rate as they lose their leaves. 
             !        The rate of biomass loss of both fine roots and leaves is presribed through the equation:
             !        \latexonly
             !        \input{turnover_clim_senes_biomass_eqn4.tex}
             !        \endlatexonly
             !        \n
             !         with ::leaffall(j) a PFT-dependent time constant which is given in ::stomate_constants),
             WHERE ( senescence(:,j) )
             turnover(:,j,ileaf,m)   = biomass(:,j,ileaf,m)   * dt / leaffall(j)
             turnover(:,j,iroot,m)   = biomass(:,j,iroot,m)   * dt / leaffall(j)
             ENDWHERE

          ELSE

             !! 3.2.2 In grasses, leaf senescence is extended to the whole plant 
             !        In grasses, leaf senescence is extended to the whole plant (all carbon pools) except to its
             !        carbohydrate reserve.      
             !JC separate tissue age
             ! for crops, I donot change the turnover by leaffall
             IF (senescence_type(j) .EQ. 'crop') THEN
                ! 3.2.2.1 crops with 'crop' phenological model
                WHERE ( senescence(:,j) )
                turnover(:,j,ileaf    ,m) = biomass(:,j,ileaf    ,m) * dt / leaffall(j)
                turnover(:,j,iroot    ,m) = biomass(:,j,iroot    ,m) * dt / leaffall(j)
                turnover(:,j,isapabove,m) = biomass(:,j,isapabove,m) * dt / leaffall(j)
                turnover(:,j,ifruit   ,m) = biomass(:,j,ifruit   ,m) * dt / leaffall(j)
                !JC MOD029 add ilabile senescence, otherwise, due to no burn of ilabile,
                !labile will increase abnormally when reserve becomes full
                turnover(:,j,ilabile  ,m) = biomass(:,j,ilabile  ,m) * dt / leaffall(j)    
                !JC MOD036 add reserve senescence for new crop phenology strategy
                turnover(:,j,icarbres  ,m) = biomass(:,j,icarbres  ,m) * dt / leaffall(j)
                ENDWHERE
             ELSE
                ! for grasses, I will not apply root turnover eventually
                ! BUT now I keep it to only test the new implementation of different tissue age
                ! NOW, I remove the root turnover during leaf fall
                ! 3.2.2.2 grass or crops based on 'mixed' phenological model
                WHERE (turnover_time(:,j) .LT. max_turnover_time(j)) 
                   turnover(:,j,ileaf    ,m) = biomass(:,j,ileaf    ,m) * dt / turnover_time(:,j)
                   turnover(:,j,isapabove,m) = biomass(:,j,isapabove,m) * dt / turnover_time(:,j)
                   ! JC MOD044 reverse the root turnover during leaffall
                   ! that's because the possible large root biomass causing maint_resp during
                   ! non-growing season that depleting reserves in dry regions (long non-growing
                   ! season)
                   ! JC 09Nov2017 remove root fall again for testing overwinter/overdrought survival
                   ! default GRM_RtoL_turn is 0.0
                   turnover(:,j,iroot    ,m) = biomass(:,j,iroot    ,m) * dt / turnover_time(:,j) * GRM_RtoL_turn 

                   turnover(:,j,ifruit   ,m) = biomass(:,j,ifruit   ,m) * dt / turnover_time(:,j)
                ENDWHERE
             ENDIF
          ENDIF      ! tree/grass
       ENDDO ! nelements

       ! update biomass pools
       biomass(:,j,ileaf    ,:) = biomass(:,j,ileaf    ,:) - turnover(:,j,ileaf    ,:)
       biomass(:,j,isapabove,:) = biomass(:,j,isapabove,:) - turnover(:,j,isapabove,:)
       biomass(:,j,iroot    ,:) = biomass(:,j,iroot    ,:) - turnover(:,j,iroot    ,:)
       biomass(:,j,ifruit   ,:) = biomass(:,j,ifruit   ,:) - turnover(:,j,ifruit   ,:)
       !JC MOD029 add ilabile senescence
       biomass(:,j,ilabile  ,:) = biomass(:,j,ilabile  ,:) - turnover(:,j,ilabile  ,:)
       !JC MOD036 add reserve senescence for new crop phenology strategy
       IF (senescence_type(j) .EQ. 'crop') THEN
         biomass(:,j,icarbres  ,:) = biomass(:,j,icarbres  ,:) - turnover(:,j,icarbres  ,:)
       ENDIF
       !JC MOD030 In reality, there should be no nutrient recycle for crops or annual
       !grasses
       ! currently, to be sure the crops can be growth, I keep the nutrient recycling
       ! for leaf and root, but not for ilabile
       ! but subject to be changed after discussion.
       
       !JC MOD036 no recycling for crops
       IF (senescence_type(j) .NE. 'crop') THEN
       ! compute the recycled fraction of nitrogen
       recycling(:,j,ileaf,initrogen) = turnover(:,j,ileaf,initrogen) * recycle_leaf(j) 
       recycling(:,j,iroot,initrogen) = turnover(:,j,iroot,initrogen) * recycle_root(j) 
       ! correct turnover for the recycled nitrogen
       turnover(:,j,ileaf,initrogen) = turnover(:,j,ileaf,initrogen) - recycling(:,j,ileaf,initrogen)
       turnover(:,j,iroot,initrogen) = turnover(:,j,iroot,initrogen) - recycling(:,j,iroot,initrogen)
       ! add the recycled nutrients to the labile pool
       biomass(:,j,ilabile,initrogen) = biomass(:,j,ilabile,initrogen)    &
                                          + recycling(:,j,ileaf,initrogen) &
                                          + recycling(:,j,iroot,initrogen)

       ! compute the recycled fraction of phosphorus
       recycling(:,j,ileaf,iphosphorus) = turnover(:,j,ileaf,iphosphorus) * p_recycle_leaf(j) 
       recycling(:,j,iroot,iphosphorus) = turnover(:,j,iroot,iphosphorus) * p_recycle_root(j) 
       ! correct turnover for the recycled phosphorus
       turnover(:,j,ileaf,iphosphorus) = turnover(:,j,ileaf,iphosphorus) - recycling(:,j,ileaf,iphosphorus)
       turnover(:,j,iroot,iphosphorus) = turnover(:,j,iroot,iphosphorus) - recycling(:,j,iroot,iphosphorus)
       ! add the recycled nutrients to the labile pool 
       biomass(:,j,ilabile,iphosphorus) = biomass(:,j,ilabile,iphosphorus)  &
                                          + recycling(:,j,ileaf,iphosphorus) &
                                          + recycling(:,j,iroot,iphosphorus)
       ENDIF

       ! JCADD 20170501 CNP allow stem nutrient resorption,
       ! stem has the same resorption rate as leaf
       IF (.NOT. is_tree(j) .AND. natural(j)) THEN       
         recycling(:,j,isapabove,initrogen) = turnover(:,j,isapabove,initrogen) * recycle_leaf(j)
         turnover(:,j,isapabove,initrogen) = turnover(:,j,isapabove,initrogen) - recycling(:,j,isapabove,initrogen)
         biomass(:,j,ilabile,initrogen) = biomass(:,j,ilabile,initrogen)    &
                                            + recycling(:,j,isapabove,initrogen)
  
         recycling(:,j,isapabove,iphosphorus) = turnover(:,j,isapabove,iphosphorus) * p_recycle_leaf(j)
         turnover(:,j,isapabove,iphosphorus) = turnover(:,j,isapabove,iphosphorus) - recycling(:,j,isapabove,iphosphorus)
         biomass(:,j,ilabile,iphosphorus) = biomass(:,j,ilabile,iphosphorus)  &
                                            + recycling(:,j,isapabove,iphosphorus)
       ENDIF ! grasses

    ENDDO        ! loop over PFTs


    !! 4. Leaf fall
    !     At a certain age, leaves fall off, even if the climate would allow a green plant
    !     all year round. Even if the meteorological conditions are favorable for leaf maintenance,
    !     plants, and in particular, evergreen trees, have to renew their leaves simply because the 
    !     old leaves become inefficient.   
    !     Roots, fruits (and stems) follow leaves. The decay rate varies with leaf age.
    !     Note that plant is not declared senescent in this case (wchich is important for allocation:
    !     if the plant loses leaves because of their age, it can renew them).
    !
    !     The leaf turnover rate due to aging of leaves is calculated using the following equation:
    !     \latexonly
    !     \input{turnover_age_senes_biomass_eqn5.tex}
    !     \endlatexonly
    !     \n
    DO j = 2,nvm ! Loop over # PFTs

       !! save old leaf mass
       !JC separate tissue age
       ! save all mass
       lm_old(:,:) = biomass(:,j,:,icarbon)

       !! initialize leaf mass change in age class
       delta_lm(:,:,:) = zero

       IF ( is_tree(j) .OR. (.NOT. natural(j)) ) THEN

          !! 4.1 Trees: leaves, roots, fruits roots and fruits follow leaves.

          !! 4.1.1 Critical age: prescribed for trees
          leaf_age_crit(:,j) = leafagecrit(j)
          root_age_crit(:,j) = rootagecrit(j)
        
       ELSE

          !! 4.2 Grasses: leaves, roots, fruits, sap follow leaves.

          !! 4.2.1 Critical leaf age depends on long-term temperature
          !        Critical leaf age depends on long-term temperature
          !        generally, lower turnover in cooler climates.
          leaf_age_crit(:,j) = &
               MIN( leafagecrit(j) * leaf_age_crit_coeff(1) , &
               MAX( leafagecrit(j) * leaf_age_crit_coeff(2) , &
               leafagecrit(j) - leaf_age_crit_coeff(3) * &
               ( t2m_longterm(:)-ZeroCelsius - leaf_age_crit_tref ) ) )

          root_age_crit(:,j) = &
               MIN( rootagecrit(j) * leaf_age_crit_coeff(1) , &
               MAX( rootagecrit(j) * leaf_age_crit_coeff(2) , &
               rootagecrit(j) - leaf_age_crit_coeff(3) * &
               ( t2m_longterm(:)-ZeroCelsius - leaf_age_crit_tref ) ) )

       END IF

       ! 4.2.2 Loop over leaf age classes
       DO m = 1, nleafages

          turnover_rate(:,:) = zero
          ! calculate turnover_rate for each grid each part
          DO ipar = 1,nparts
            IF ((ipar .EQ. ileaf) .OR. (ipar .EQ. isapabove)) THEN
              WHERE ( leaf_age(:,j,ipar,m) .GT. leaf_age_crit(:,j)/2. )

                turnover_rate(:,ipar) =  &
                    MIN( 0.99_r_std, dt / ( leaf_age_crit(:,j) * &
                    ( leaf_age_crit(:,j) / leaf_age(:,j,ipar,m) )**quatre ) )

              ENDWHERE
            ! root life span of grasses is much longer than leaf
            ! for trees, it is different
            ! now rootagecrit can be read to control
            ELSE IF (ipar .EQ. iroot) THEN
              WHERE ( leaf_age(:,j,ipar,m) .GT. root_age_crit(:,j)/2. )

                turnover_rate(:,ipar) =  &
                    MIN( 0.99_r_std, dt / ( root_age_crit(:,j) * &
                    ( root_age_crit(:,j) / leaf_age(:,j,ipar,m) )**quatre ) )

              ENDWHERE
            ELSE 
              turnover_rate(:,ipar) = zero
            ENDIF
          ENDDO
          !! for C pools
          ! leaf
          dturnover(:)          = biomass(:,j,ileaf,icarbon) * leaf_frac(:,j,ileaf,m) * turnover_rate(:,ileaf)
          turnover(:,j,ileaf,icarbon) = turnover(:,j,ileaf,icarbon) + dturnover(:)
          biomass (:,j,ileaf,icarbon) = biomass (:,j,ileaf,icarbon) - dturnover(:)
          ! save leaf mass change
          delta_lm(:,ileaf,m) = - dturnover(:)
          ! root
          dturnover(:) = biomass(:,j,iroot,icarbon) * leaf_frac(:,j,iroot,m) * turnover_rate(:,iroot)
          turnover(:,j,iroot,icarbon) = turnover(:,j,iroot,icarbon) + dturnover(:)
          biomass(:,j,iroot,icarbon) = biomass(:,j,iroot,icarbon) - dturnover(:)
          ! save leaf mass change
          delta_lm(:,iroot,m) = - dturnover(:)
          ! fruit
          dturnover(:) = biomass(:,j,ifruit,icarbon) * leaf_frac(:,j,ileaf,m) * turnover_rate(:,ileaf)
          turnover(:,j,ifruit,icarbon) = turnover(:,j,ifruit,icarbon) + dturnover(:)
          biomass(:,j,ifruit,icarbon) = biomass(:,j,ifruit,icarbon) - dturnover(:)
          ! save leaf mass change
          delta_lm(:,ifruit,m) = - dturnover(:)
          ! for grasses/crops, sapabove can also turnover
          IF (.NOT. is_tree(j)) THEN
            dturnover(:) = biomass(:,j,isapabove,icarbon) * leaf_frac(:,j,isapabove,m) * turnover_rate(:,isapabove)
            turnover(:,j,isapabove,icarbon) = turnover(:,j,isapabove,icarbon) + dturnover(:)
            biomass(:,j,isapabove,icarbon) = biomass(:,j,isapabove,icarbon) - dturnover(:)
          ! save leaf mass change
          delta_lm(:,isapabove,m) = - dturnover(:)
          ENDIF

          !! for N pools there is resorption
          ! leaf
          dturnover(:) = biomass(:,j,ileaf,initrogen) * leaf_frac(:,j,ileaf,m) * turnover_rate(:,ileaf)
          biomass(:,j,ilabile,initrogen) = biomass(:,j,ilabile,initrogen) + recycle_leaf(j) * dturnover(:)
          biomass(:,j,ileaf,initrogen)   = biomass(:,j,ileaf,initrogen)   - dturnover(:)
          !! add these fluxes to the recycling and turnover fluxes:
          recycling(:,j,ileaf,initrogen) = recycling(:,j,ileaf,initrogen) + recycle_leaf(j) * dturnover(:)
          turnover(:,j,ileaf,initrogen)  = turnover(:,j,ileaf,initrogen)  + ( 1.0 - recycle_leaf(j) ) * dturnover(:)
          ! root
          dturnover(:) = biomass(:,j,iroot,initrogen) * leaf_frac(:,j,iroot,m) * turnover_rate(:,iroot)
          biomass(:,j,ilabile,initrogen) = biomass(:,j,ilabile,initrogen) + recycle_root(j) * dturnover(:)
          biomass(:,j,iroot,initrogen)   = biomass(:,j,iroot,initrogen)   - dturnover(:)
          !! add these fluxes to the recycling and turnover fluxes:
          recycling(:,j,iroot,initrogen) = recycling(:,j,iroot,initrogen) + recycle_root(j) * dturnover(:)
          turnover(:,j,iroot,initrogen)  = turnover(:,j,iroot,initrogen)  + ( 1.0 - recycle_root(j) ) * dturnover(:)
          ! for grasses, there should be no fruit turnover until scenescene
          ! fruit (no recycling)
          dturnover(:) = biomass(:,j,ifruit,initrogen) * leaf_frac(:,j,ileaf,m) * turnover_rate(:,ileaf)
          turnover(:,j,ifruit,initrogen) = turnover(:,j,ifruit,initrogen) + dturnover(:)
          biomass(:,j,ifruit,initrogen) = biomass(:,j,ifruit,initrogen) - dturnover(:)
          ! for grasses/crops, sapabove can also turnover
          IF (.NOT. is_tree(j)) THEN
            dturnover(:) = biomass(:,j,isapabove,initrogen) * leaf_frac(:,j,isapabove,m) * turnover_rate(:,isapabove)
            biomass(:,j,ilabile,initrogen) = biomass(:,j,ilabile,initrogen) + recycle_leaf(j) * dturnover(:)
            biomass(:,j,isapabove,initrogen)   = biomass(:,j,isapabove,initrogen)   - dturnover(:)
            !! add these fluxes to the recycling and turnover fluxes:
            recycling(:,j,isapabove,initrogen) = recycling(:,j,isapabove,initrogen) + recycle_leaf(j) * dturnover(:)
            turnover(:,j,isapabove,initrogen)  = turnover(:,j,isapabove,initrogen)  + ( 1.0 - recycle_leaf(j) ) * dturnover(:)
          ENDIF

          !! for P pools there is resorption
          ! leaf
          dturnover(:) = biomass(:,j,ileaf,iphosphorus) * leaf_frac(:,j,ileaf,m) * turnover_rate(:,ileaf)
          biomass(:,j,ilabile,iphosphorus) = biomass(:,j,ilabile,iphosphorus) + recycle_leaf(j) * dturnover(:)
          biomass(:,j,ileaf,iphosphorus)   = biomass(:,j,ileaf,iphosphorus)   - dturnover(:)
          !! add these fluxes to the recycling and turnover fluxes:
          recycling(:,j,ileaf,iphosphorus) = recycling(:,j,ileaf,iphosphorus) + recycle_leaf(j) * dturnover(:)
          turnover(:,j,ileaf,iphosphorus)  = turnover(:,j,ileaf,iphosphorus)  + ( 1.0 - recycle_leaf(j) ) * dturnover(:)
          ! root
          dturnover(:) = biomass(:,j,iroot,iphosphorus) * leaf_frac(:,j,iroot,m) * turnover_rate(:,iroot)
          biomass(:,j,ilabile,iphosphorus) = biomass(:,j,ilabile,iphosphorus) + recycle_root(j) * dturnover(:)
          biomass(:,j,iroot,iphosphorus)   = biomass(:,j,iroot,iphosphorus)   - dturnover(:)
          !! add these fluxes to the recycling and turnover fluxes:
          recycling(:,j,iroot,iphosphorus) = recycling(:,j,iroot,iphosphorus) + recycle_root(j) * dturnover(:)
          turnover(:,j,iroot,iphosphorus)  = turnover(:,j,iroot,iphosphorus)  + ( 1.0 - recycle_root(j) ) * dturnover(:)
! for grasses, there should be no fruit turnover until scenescene
          ! fruit (no recycling)
          dturnover(:) = biomass(:,j,ifruit,iphosphorus) * leaf_frac(:,j,ileaf,m) * turnover_rate(:,ileaf)
          turnover(:,j,ifruit,iphosphorus) = turnover(:,j,ifruit,iphosphorus) + dturnover(:)
          biomass(:,j,ifruit,iphosphorus) = biomass(:,j,ifruit,iphosphorus) - dturnover(:)
          ! for grasses/crops, sapabove can also turnover
          IF (.NOT. is_tree(j)) THEN
            dturnover(:) = biomass(:,j,isapabove,iphosphorus) * leaf_frac(:,j,isapabove,m) * turnover_rate(:,isapabove)
            biomass(:,j,ilabile,iphosphorus) = biomass(:,j,ilabile,iphosphorus) + recycle_leaf(j) * dturnover(:)
            biomass(:,j,isapabove,iphosphorus)   = biomass(:,j,isapabove,iphosphorus)   - dturnover(:)
            !! add these fluxes to the recycling and turnover fluxes:
            recycling(:,j,isapabove,iphosphorus) = recycling(:,j,isapabove,iphosphorus) + recycle_leaf(j) * dturnover(:)
            turnover(:,j,isapabove,iphosphorus)  = turnover(:,j,isapabove,iphosphorus)  + ( 1.0 - recycle_leaf(j) ) * dturnover(:)
          ENDIF
          
       ENDDO

       !! 4.3 Recalculate the fraction of leaf biomass in each leaf age class.
       !      Older leaves will fall more fast than younger leaves and therefore 
       !      the leaf age distribution needs to be recalculated after turnover. 
       !      The fraction of biomass in each leaf class is updated using the following equation:
       !      \latexonly
       !      \input{turnover_update_LeafAgeDistribution_eqn6.tex}
       !      \endlatexonly
       !      \n
       !
       !      new fraction = new leaf mass of that fraction / new total leaf mass
       !                   = (old fraction*old total leaf mass ::lm_old(:) + biomass change of that fraction ::delta_lm(:,m)  ) /
       !                     new total leaf mass ::biomass(:,j,ileaf
       DO m = 1, nleafages
          
          WHERE ( biomass(:,j,:,icarbon) .GT. min_sechiba )
             leaf_frac(:,j,:,m) = ( leaf_frac(:,j,:,m)*lm_old(:,:) + delta_lm(:,:,m) ) / biomass(:,j,:,icarbon)
          ELSEWHERE
             leaf_frac(:,j,:,m) = zero
          ENDWHERE
        
       ENDDO
       
    ENDDO         ! loop over PFTs

    !! 5. New (provisional) LAI 
    !     ::lai(:,j) is determined from the leaf biomass ::biomass(:,j,ileaf,icarbon) and the 
    !     specific leaf surface :: sla(j) (m^2 gC^{-1})
    !     The leaf area index is updated using the following equation:
    !     \latexonly
    !     \input{turnover_update_LAI_eqn7.tex}
    !     \endlatexonly
    !     \n

    !    lai(:,ibare_sechiba) = zero
    !    DO j = 2, nvm ! Loop over # PFTs
    !        lai(:,j) = biomass(:,j,ileaf,icarbon) * sla(j)
    !    ENDDO

    !! 6. Definitely drop all leaves if there is a very low leaf mass during senescence.

    !     Both for deciduous trees and grasses same conditions are checked:
    !     If biomass is large enough (i.e. when it is greater than zero), 
    !     AND when senescence is set to true
    !     AND the leaf biomass drops below a critical minimum biomass level (i.e. when it is lower than half
    !     the minimum initial LAI ::lai_initmin(j) divided by the specific leaf area ::sla(j),
    !     ::lai_initmin(j) is set to 0.3 in stomate_data.f90 and sla is a constant that is set to 0.015366 m2/gC), 
    !     If these conditions are met, the flag ::shed_rest(:) is set to TRUE.
    !
    !     After this, the biomass of different carbon pools both for trees and grasses is set to zero
    !     and the mean leaf age is reset to zero.
    !     Finally, the leaf fraction and leaf age of the different leaf age classes is set to zero.
    DO j = 2,nvm ! Loop over # PFTs

       shed_rest(:) = .FALSE.

       !! 6.1 For deciduous trees: next to leaves, also fruits and fine roots are dropped 
       !      For deciduous trees: next to leaves, also fruits and fine roots are dropped: fruit ::biomass(:,j,ifruit) 
       !      and fine root ::biomass(:,j,iroot) carbon pools are set to zero.
       IF ( is_tree(j) .AND. ( senescence_type(j) .NE. 'none' ) ) THEN

          ! check whether we shed the remaining leaves
          WHERE ( ( biomass(:,j,ileaf,icarbon) .GT. zero ) .AND. senescence(:,j) .AND. &
               ( biomass(:,j,ileaf,icarbon) .LT. (lai_initmin(j) / 2.)/sla_calc(:,j) ) )

             shed_rest(:) = .TRUE.
            
             ! first leaves and root for which nutrients are partly recycled:
             turnover(:,j,ileaf,icarbon)  = turnover(:,j,ileaf,icarbon) + biomass(:,j,ileaf,icarbon)
             turnover(:,j,iroot,icarbon)  = turnover(:,j,iroot,icarbon) + biomass(:,j,iroot,icarbon)

             ! for N pools
             biomass(:,j,ilabile,initrogen) = biomass(:,j,ilabile,initrogen) &
                  + recycle_leaf(j) * biomass(:,j,ileaf,initrogen) &
                  + recycle_root(j) * biomass(:,j,iroot,initrogen)             
             ! track turnover vs recycling fluxes:
             turnover(:,j,ileaf,initrogen) = turnover(:,j,ileaf,initrogen) &
                  + (un-recycle_leaf(j)) * biomass(:,j,ileaf,initrogen)
             turnover(:,j,iroot,initrogen) = turnover(:,j,iroot,initrogen) &
                  + (un-recycle_root(j)) * biomass(:,j,iroot,initrogen)

             recycling(:,j,ileaf,initrogen) = recycling(:,j,ileaf,initrogen) &
                  + recycle_leaf(j) * biomass(:,j,ileaf,initrogen)
             recycling(:,j,iroot,initrogen) = recycling(:,j,iroot,initrogen) &
                  + recycle_root(j) * biomass(:,j,iroot,initrogen)

             ! for P pools
             biomass(:,j,ilabile,iphosphorus) = biomass(:,j,ilabile,iphosphorus) &
                  + p_recycle_leaf(j) * biomass(:,j,ileaf,iphosphorus) &
                  + p_recycle_root(j) * biomass(:,j,iroot,iphosphorus)             
             ! track turnover vs recycling fluxes:
             turnover(:,j,ileaf,iphosphorus) = turnover(:,j,ileaf,iphosphorus) &
                  + (un-p_recycle_leaf(j)) * biomass(:,j,ileaf,iphosphorus)
             turnover(:,j,iroot,iphosphorus) = turnover(:,j,iroot,iphosphorus) &
                  + (un-p_recycle_root(j)) * biomass(:,j,iroot,iphosphorus)

             recycling(:,j,ileaf,iphosphorus) = recycling(:,j,ileaf,iphosphorus) &
                  + p_recycle_leaf(j) * biomass(:,j,ileaf,iphosphorus)
             recycling(:,j,iroot,iphosphorus) = recycling(:,j,iroot,iphosphorus) &
                  + p_recycle_root(j) * biomass(:,j,iroot,iphosphorus)
             
             ! fruits dont recycle 
             turnover(:,j,ifruit,icarbon)     = turnover(:,j,ifruit,icarbon)     + biomass(:,j,ifruit,icarbon)
             turnover(:,j,ifruit,initrogen)   = turnover(:,j,ifruit,initrogen)   + biomass(:,j,ifruit,initrogen)
             turnover(:,j,ifruit,iphosphorus) = turnover(:,j,ifruit,iphosphorus) + biomass(:,j,ifruit,iphosphorus)
             ! update leaf, root, fruit biomass
             biomass(:,j,ileaf, icarbon)  = zero
             biomass(:,j,iroot, icarbon)  = zero
             biomass(:,j,ifruit,icarbon)  = zero

             biomass(:,j,ileaf, initrogen)  = zero
             biomass(:,j,iroot, initrogen)  = zero
             biomass(:,j,ifruit,initrogen)  = zero

             biomass(:,j,ileaf, iphosphorus)  = zero
             biomass(:,j,iroot, iphosphorus)  = zero
             biomass(:,j,ifruit,iphosphorus)  = zero
             ! reset leaf age
             leaf_meanage(:,j,ileaf) = zero
             leaf_meanage(:,j,iroot) = zero
             leaf_meanage(:,j,ifruit) = zero
             lai(:,j) = zero
          ENDWHERE

       ENDIF

       !! 6.2 For grasses: all aboveground carbon pools, except the carbohydrate reserves are affected: 
       !      For grasses: all aboveground carbon pools, except the carbohydrate reserves are affected: 
       !      fruit ::biomass(:,j,ifruit,icarbon), fine root ::biomass(:,j,iroot,icarbon) and sapwood above 
       !      ::biomass(:,j,isapabove,icarbon) carbon pools are set to zero. 
       ! JC MOD023 only shed root for crops but not for grasses
       IF ( .NOT. is_tree(j) .AND. .NOT. natural(j)) THEN

          ! Shed the remaining leaves if LAI very low.
          WHERE ( ( biomass(:,j,ileaf,icarbon) .GT. zero ) .AND. senescence(:,j) .AND. &
               ! lai_initmin for grass is 0.1
               (  biomass(:,j,ileaf,icarbon) .LT. (lai_initmin(j) / 2.)/sla_calc(:,j) ))

             shed_rest(:) = .TRUE.
             !JC MOD036 new crop strategy, shed all tissues and reserves
             ! no recycling:
             turnover(:,j,ileaf,icarbon)         = turnover(:,j,ileaf,icarbon)        + biomass(:,j,ileaf,icarbon)
             turnover(:,j,iroot,icarbon)         = turnover(:,j,iroot,icarbon)        + biomass(:,j,iroot,icarbon)
             turnover(:,j,ileaf,initrogen)       = turnover(:,j,ileaf,initrogen)      + biomass(:,j,ileaf,initrogen)
             turnover(:,j,iroot,initrogen)       = turnover(:,j,iroot,initrogen)      + biomass(:,j,iroot,initrogen)
             turnover(:,j,ileaf,iphosphorus)     = turnover(:,j,ileaf,iphosphorus)    + biomass(:,j,ileaf,iphosphorus)
             turnover(:,j,iroot,iphosphorus)     = turnover(:,j,iroot,iphosphorus)    + biomass(:,j,iroot,iphosphorus)

             turnover(:,j,isapabove,icarbon)     = turnover(:,j,isapabove,icarbon)    + biomass(:,j,isapabove,icarbon)
             turnover(:,j,ifruit   ,icarbon)     = turnover(:,j,ifruit   ,icarbon)    + biomass(:,j,ifruit   ,icarbon)
             turnover(:,j,isapabove,initrogen)   = turnover(:,j,isapabove,initrogen)  + biomass(:,j,isapabove,initrogen)
             turnover(:,j,ifruit   ,initrogen)   = turnover(:,j,ifruit   ,initrogen)  + biomass(:,j,ifruit   ,initrogen)
             turnover(:,j,isapabove,iphosphorus) = turnover(:,j,isapabove,iphosphorus)+ biomass(:,j,isapabove,iphosphorus)
             turnover(:,j,ifruit   ,iphosphorus) = turnover(:,j,ifruit   ,iphosphorus)+ biomass(:,j,ifruit   ,iphosphorus)

             turnover(:,j,ilabile,icarbon)       = turnover(:,j,ilabile,icarbon)      + biomass(:,j,ilabile,icarbon)
             turnover(:,j,icarbres,icarbon)      = turnover(:,j,icarbres,icarbon)     + biomass(:,j,icarbres,icarbon)
             turnover(:,j,ilabile,initrogen)     = turnover(:,j,ilabile,initrogen)    + biomass(:,j,ilabile,initrogen)
             turnover(:,j,icarbres,initrogen)    = turnover(:,j,icarbres,initrogen)   + biomass(:,j,icarbres,initrogen)
             turnover(:,j,ilabile,iphosphorus)   = turnover(:,j,ilabile,iphosphorus)  + biomass(:,j,ilabile,iphosphorus)
             turnover(:,j,icarbres,iphosphorus)  = turnover(:,j,icarbres,iphosphorus) + biomass(:,j,icarbres,iphosphorus)
             ! update tissues which were shed
             biomass(:,j,ileaf    ,icarbon)  = zero
             biomass(:,j,isapabove,icarbon)  = zero
             biomass(:,j,iroot    ,icarbon)  = zero
             biomass(:,j,ifruit   ,icarbon)  = zero
             biomass(:,j,ilabile  ,icarbon)  = zero
             biomass(:,j,icarbres,icarbon)  = zero

             biomass(:,j,ileaf    ,initrogen)  = zero
             biomass(:,j,isapabove,initrogen)  = zero
             biomass(:,j,iroot    ,initrogen)  = zero
             biomass(:,j,ifruit   ,initrogen)  = zero
             biomass(:,j,ilabile  ,initrogen)  = zero
             biomass(:,j,icarbres,initrogen)  = zero

             biomass(:,j,ileaf    ,iphosphorus)  = zero
             biomass(:,j,isapabove,iphosphorus)  = zero
             biomass(:,j,iroot    ,iphosphorus)  = zero
             biomass(:,j,ifruit   ,iphosphorus)  = zero
             biomass(:,j,ilabile  ,iphosphorus)  = zero
             biomass(:,j,icarbres,iphosphorus)  = zero
             ! reset leaf age
             leaf_meanage(:,j,ileaf) = zero
             leaf_meanage(:,j,isapabove) = zero
             leaf_meanage(:,j,iroot) = zero
             leaf_meanage(:,j,ifruit) = zero
             lai(:,j) = zero
          ENDWHERE

       ELSE IF ( .NOT. is_tree(j) .AND. natural(j)) THEN
          ! Shed the remaining leaves if LAI very low.
          WHERE ( ( biomass(:,j,ileaf,icarbon) .GT. zero ) .AND. senescence(:,j) .AND. &
               ! lai_initmin for grass is 0.1
               (  biomass(:,j,ileaf,icarbon) .LT. (lai_initmin(j) / 2.)/sla_calc(:,j) ))

             shed_rest(:) = .TRUE.

             ! first take about tissues from which nutrients are recycled:

             turnover(:,j,ileaf,icarbon)     = turnover(:,j,ileaf,icarbon)     + biomass(:,j,ileaf,icarbon)
             ! nitrogen
             turnover(:,j,ileaf,initrogen) = turnover(:,j,ileaf,initrogen) &
                  + (un-recycle_leaf(j)) * biomass(:,j,ileaf,initrogen)
             recycling(:,j,ileaf,initrogen) = recycling(:,j,ileaf,initrogen) &
                  + recycle_leaf(j) * biomass(:,j,ileaf,initrogen)
             biomass(:,j,ilabile,initrogen) = biomass(:,j,ilabile,initrogen) &
                  + recycle_leaf(j) * biomass(:,j,ileaf,initrogen) 
             ! phosphorus
             turnover(:,j,ileaf,iphosphorus) = turnover(:,j,ileaf,iphosphorus) &
                  + (un-p_recycle_leaf(j)) * biomass(:,j,ileaf,iphosphorus)
             recycling(:,j,ileaf,iphosphorus) = recycling(:,j,ileaf,iphosphorus) &
                  + p_recycle_leaf(j) * biomass(:,j,ileaf,iphosphorus)
             biomass(:,j,ilabile,iphosphorus) = biomass(:,j,ilabile,iphosphorus) &
                  + p_recycle_leaf(j) * biomass(:,j,ileaf,iphosphorus) 

             ! no recycling:
             turnover(:,j,isapabove,icarbon)     = turnover(:,j,isapabove,icarbon) + biomass(:,j,isapabove,icarbon)
             turnover(:,j,ifruit   ,icarbon)     = turnover(:,j,ifruit   ,icarbon) + biomass(:,j,ifruit   ,icarbon)
             turnover(:,j,isapabove,initrogen)   = turnover(:,j,isapabove,initrogen) + biomass(:,j,isapabove,initrogen)
             turnover(:,j,ifruit   ,initrogen)   = turnover(:,j,ifruit   ,initrogen) + biomass(:,j,ifruit   ,initrogen)
             turnover(:,j,isapabove,iphosphorus) = turnover(:,j,isapabove,iphosphorus) + biomass(:,j,isapabove,iphosphorus)
             turnover(:,j,ifruit   ,iphosphorus) = turnover(:,j,ifruit   ,iphosphorus) + biomass(:,j,ifruit   ,iphosphorus)

             ! update tissues which were shed
             biomass(:,j,ileaf    ,icarbon)  = zero
             biomass(:,j,isapabove,icarbon)  = zero
             biomass(:,j,ifruit   ,icarbon)  = zero

             biomass(:,j,ileaf    ,initrogen)  = zero
             biomass(:,j,isapabove,initrogen)  = zero
             biomass(:,j,ifruit   ,initrogen)  = zero

             biomass(:,j,ileaf    ,iphosphorus)  = zero
             biomass(:,j,isapabove,iphosphorus)  = zero
             biomass(:,j,ifruit   ,iphosphorus)  = zero

             ! reset leaf age
             leaf_meanage(:,j,ileaf) = zero
             leaf_meanage(:,j,isapabove) = zero
             leaf_meanage(:,j,ifruit) = zero

          ENDWHERE

       ENDIF

       ! JC debug for detecting marginal biomass
       ! we should shed all tissues because it will nevel grow again
       ! and it will violate phenology leaf onset, because possiblly no nutrient at all
       shed_minimum(:) = .FALSE.
       tot_bio(:) = zero
       tot_bio(:) = SUM(biomass(:,j,:,icarbon),DIM=2)
       DO ipar = 1, nparts
         DO m = 1, nelements
           WHERE (tot_bio(:) .LE. min_stomate)
             shed_minimum(:) = .TRUE.
             turnover(:,j,ipar,m) = turnover(:,j,ipar,m) + biomass(:,j,ipar,m)
             biomass(:,j,ipar,m) = zero
           ENDWHERE
         ENDDO
       ENDDO
 
       !IF (ANY(shed_minimum(:))) THEN
       !  WRITE(numout,*) 'we shed marginal biomass (sum(biomass) < min_stomate)'
       !  WRITE(numout,*) 'in : ', COUNT(shed_minimum)
       !  WRITE(numout,*) 'of (npts): ', npts
       !  WRITE(numout,*) 'for PFT: ', j
       !ENDIF
       !! 6.3 Reset the leaf age structure: the leaf fraction and leaf age of the different leaf age classes is set to zero.
      
       DO m = 1, nleafages
         DO ipar = 1,nparts
           IF (ipar .NE. iroot) THEN
             WHERE ( shed_rest(:) .OR. shed_minimum(:) )
   
                leaf_age(:,j,ipar,m) = zero
                leaf_frac(:,j,ipar,m) = zero
   
             ENDWHERE
           ENDIF
         ENDDO
       ENDDO

    ENDDO          ! loop over PFTs

    !! 7. Herbivore activity: elephants, cows, gazelles but no lions.
 
    !     Herbivore activity affects the biomass of leaves and fruits as well 
    !     as stalks (only for grasses). Herbivore activity does not modify leaf 
    !     age structure. Herbivores ::herbivores(:,j) is the time constant of 
    !     probability of a leaf to be eaten by a herbivore, and is calculated in 
    !     ::stomate_season. following Mc Naughton et al. [1989].

    IF ( ok_herbivores ) THEN

       ! If the herbivore activity is allowed (if ::ok_herbivores is true, which is set in run.def), 
       ! remove the amount of biomass consumed by herbivory from the leaf biomass ::biomass(:,j,ileaf,icarbon) and 
       ! the fruit biomass ::biomass(:,j,ifruit,icarbon).
       ! The daily amount consumed equals the biomass multiplied by 1 day divided by the time constant ::herbivores(:,j).
       DO j = 2,nvm ! Loop over # PFTs

          IF ( is_tree(j) ) THEN

             !! For trees: only the leaves and fruit carbon pools are affected

             WHERE (biomass(:,j,ileaf,icarbon) .GT. zero)
                ! added by shilong
                WHERE (herbivores(:,j).GT. min_sechiba)
                   !leaf carbon
                   dturnover(:) = biomass(:,j,ileaf,icarbon) * dt / herbivores(:,j)
                   turnover(:,j,ileaf,icarbon) = turnover(:,j,ileaf,icarbon) + dturnover(:)
                   biomass(:,j,ileaf,icarbon) = biomass(:,j,ileaf,icarbon) - dturnover(:)
                   !leaf nitrogen
                   dturnover(:) = biomass(:,j,ileaf,initrogen) * dt / herbivores(:,j) 
                   turnover(:,j,ileaf,initrogen) = turnover(:,j,ileaf,initrogen) + dturnover(:)
                   biomass(:,j,ileaf,initrogen) = biomass(:,j,ileaf,initrogen) - dturnover(:) 
                   !leaf phosphorus
                   dturnover(:) = biomass(:,j,ileaf,iphosphorus) * dt / herbivores(:,j) 
                   turnover(:,j,ileaf,iphosphorus) = turnover(:,j,ileaf,iphosphorus) + dturnover(:)
                   biomass(:,j,ileaf,iphosphorus) = biomass(:,j,ileaf,iphosphorus) - dturnover(:) 
                   !fruit carbon
                   dturnover(:) = biomass(:,j,ifruit,icarbon) * dt / herbivores(:,j)
                   turnover(:,j,ifruit,icarbon) = turnover(:,j,ifruit,icarbon) + dturnover(:)
                   biomass(:,j,ifruit,icarbon) = biomass(:,j,ifruit,icarbon) - dturnover(:)
                   !fruit nitrogen
                   dturnover(:) = biomass(:,j,ifruit,initrogen) * dt / herbivores(:,j)
                   turnover(:,j,ifruit,initrogen) = turnover(:,j,ifruit,initrogen) + dturnover(:) 
                   biomass(:,j,ifruit,initrogen) = biomass(:,j,ifruit,initrogen) - dturnover(:) 
                   !fruit phosphorus
                   dturnover(:) = biomass(:,j,ifruit,iphosphorus) * dt / herbivores(:,j)
                   turnover(:,j,ifruit,iphosphorus) = turnover(:,j,ifruit,iphosphorus) + dturnover(:) 
                   biomass(:,j,ifruit,iphosphorus) = biomass(:,j,ifruit,iphosphorus) - dturnover(:) 

                   ! DSGpiss: the nutrients consumed by herbivores should mostly
                   ! go to readily available nutrient pools (ammonia, nitrate,
                   ! dissolved P)
                   ! status: not done yet

                ENDWHERE
             ENDWHERE

          ELSE

             ! For grasses: all aboveground carbon pools are affected: leaves, fruits and sapwood above
             WHERE ( biomass(:,j,ileaf,icarbon) .GT. zero )
                ! added by shilong
                WHERE (herbivores(:,j) .GT. min_sechiba)

                   dturnover(:) = biomass(:,j,ileaf,icarbon) * dt / herbivores(:,j)
                   turnover(:,j,ileaf,icarbon) = turnover(:,j,ileaf,icarbon) + dturnover(:)
                   biomass(:,j,ileaf,icarbon) = biomass(:,j,ileaf,icarbon) - dturnover(:)

                   dturnover(:) = biomass(:,j,ileaf,initrogen) * dt / herbivores(:,j) 
                   turnover(:,j,ileaf,initrogen) = turnover(:,j,ileaf,initrogen) + dturnover(:)
                   biomass(:,j,ileaf,initrogen) = biomass(:,j,ileaf,initrogen) - dturnover(:) 

                   dturnover(:) = biomass(:,j,isapabove,icarbon) * dt / herbivores(:,j)
                   turnover(:,j,isapabove,icarbon) = turnover(:,j,isapabove,icarbon) + dturnover(:)
                   biomass(:,j,isapabove,icarbon) = biomass(:,j,isapabove,icarbon) - dturnover(:)

                   dturnover(:) = biomass(:,j,isapabove,initrogen)  * dt / herbivores(:,j)
                   turnover(:,j,isapabove,initrogen) = turnover(:,j,isapabove,initrogen) + dturnover(:) 
                   biomass(:,j,isapabove,initrogen) = biomass(:,j,isapabove,initrogen) - dturnover(:) 

                   dturnover(:) = biomass(:,j,ifruit,icarbon) * dt / herbivores(:,j)
                   turnover(:,j,ifruit,icarbon) = turnover(:,j,ifruit,icarbon) + dturnover(:)
                   biomass(:,j,ifruit,icarbon) = biomass(:,j,ifruit,icarbon) - dturnover(:)

                   dturnover(:) = biomass(:,j,ifruit,initrogen) * dt / herbivores(:,j)
                   turnover(:,j,ifruit,initrogen) = turnover(:,j,ifruit,initrogen) + dturnover(:) 
                   biomass(:,j,ifruit,initrogen) = biomass(:,j,ifruit,initrogen) - dturnover(:) 

                   ! same fore P
                   dturnover(:)                        = biomass(:,j,ileaf,iphosphorus) * dt / herbivores(:,j) 
                   turnover(:,j,ileaf,iphosphorus)     = turnover(:,j,ileaf,iphosphorus) + dturnover(:)
                   biomass(:,j,ileaf,iphosphorus)      = biomass(:,j,ileaf,iphosphorus) - dturnover(:) 
                   dturnover(:)                        = biomass(:,j,isapabove,iphosphorus)  * dt / herbivores(:,j)
                   turnover(:,j,isapabove,iphosphorus) = turnover(:,j,isapabove,iphosphorus) + dturnover(:) 
                   biomass(:,j,isapabove,iphosphorus)  = biomass(:,j,isapabove,iphosphorus) - dturnover(:) 
                   dturnover(:)                        = biomass(:,j,ifruit,iphosphorus) * dt / herbivores(:,j)
                   turnover(:,j,ifruit,iphosphorus)    = turnover(:,j,ifruit,iphosphorus) + dturnover(:) 
                   biomass(:,j,ifruit,iphosphorus)     = biomass(:,j,ifruit,iphosphorus) - dturnover(:) 

                   ! DSGpiss: the nutrients consumed by herbivores should mostly
                   ! go to readily available nutrient pools (ammonia, nitrate,
                   ! dissolved P)
                   ! status: not done yet
                   
                ENDWHERE

             ENDWHERE

          ENDIF  ! tree/grass?

       ENDDO    ! loop over PFTs

    ENDIF ! end herbivores

    !! 8. Fruit turnover for trees

    !     Fruit turnover for trees: trees simply lose their fruits with a time constant ::tau_fruit(j), 
    !     that is set to 90 days for all PFTs in ::stomate_constants

    DO k = 1,nelements 
       DO j = 2,nvm ! Loop over # PFTs
          IF ( is_tree(j) ) THEN

             dturnover(:) = biomass(:,j,ifruit,k) * dt / tau_fruit(j)
             turnover(:,j,ifruit,k) = turnover(:,j,ifruit,k) + dturnover(:)
             biomass(:,j,ifruit,k) = biomass(:,j,ifruit,k) - dturnover(:)
             
          ENDIF
       ENDDO       ! loop over PFTs
    END DO

    !! 9 Conversion of sapwood to heartwood both for aboveground and belowground sapwood and heartwood.

    !   Following LPJ (Sitch et al., 2003), sapwood biomass is converted into heartwood biomass 
    !   with a time constant tau ::tau_sap(j) of 1 year. DSG: not anymore 1 yr
    !   Note that this biomass conversion is not added to "turnover" as the biomass is not lost!

    ! init:
    bm_sap_2_heartw =  zero 

    DO j = 2,nvm ! Loop over # PFTs

       IF ( is_tree(j) ) THEN

          !! For the recalculation of age in 9.2 (in case the vegetation is not dynamic ie. ::ok_dgvm is FALSE), 
          !! the heartwood above and below is stored in ::hw_old(:).
          IF ( .NOT. ok_dgvm ) THEN
             hw_old(:) = biomass(:,j,iheartabove,icarbon) + biomass(:,j,iheartbelow,icarbon)
          ENDIF

          !! 9.1 Calculate the rate of sapwood to heartwood conversion 
          !      Calculate the rate of sapwood to heartwood conversion with the time constant ::tau_sap(j) 
          !      and update aboveground and belowground sapwood ::biomass(:,j,isapabove) and ::biomass(:,j,isapbelow)
          !      and heartwood ::biomass(:,j,iheartabove) and ::biomass(:,j,iheartbelow).

          DO k = 1,nelements

             ! Above the ground
             sapconv(:) = biomass(:,j,isapabove,k) * dt / tau_sap(j)
             biomass(:,j,isapabove,k) = biomass(:,j,isapabove,k) - sapconv(:)
             biomass(:,j,iheartabove,k) =  biomass(:,j,iheartabove,k) + sapconv(:)
             ! for output:
             bm_sap_2_heartw(:,j,k,1) = sapconv(:)
             
             ! Below the ground
             sapconv(:) = biomass(:,j,isapbelow,k) * dt / tau_sap(j)
             biomass(:,j,isapbelow,k) = biomass(:,j,isapbelow,k) - sapconv(:)
             biomass(:,j,iheartbelow,k) =  biomass(:,j,iheartbelow,k) + sapconv(:)
             ! for output:
             bm_sap_2_heartw(:,j,k,2) = sapconv(:)

          END DO

          !! 9.2 If the vegetation is not dynamic, the age of the plant is decreased. 
          !      The updated heartwood, the sum of new heartwood above and new heartwood below after 
          !      converting sapwood to heartwood, is saved as ::hw_new(:) .
          !      Creation of new heartwood decreases the age of the plant with a factor that is determined by: 
          !      old heartwood ::hw_old(:) divided by the new heartwood ::hw_new(:)
          IF ( .NOT. ok_dgvm ) THEN

             hw_new(:) = biomass(:,j,iheartabove,icarbon) + biomass(:,j,iheartbelow,icarbon)

             WHERE ( hw_new(:) .GT. min_sechiba )

                age(:,j) = age(:,j) * hw_old(:)/hw_new(:)

             ENDWHERE

          ENDIF

       ENDIF

    ENDDO       ! loop over PFTs

    DO k=1,nelements
       IF     (k == icarbon) THEN
          element_str(k) = '_c'
       ELSEIF (k == initrogen) THEN
          element_str(k) = '_n'
       ELSEIF (k == iphosphorus) THEN
          element_str(k) = '_p'
       ELSE
          STOP 'Define element_str'
       ENDIF

       CALL xios_orchidee_send_field('BM_ALLOC_HEART_AB'//TRIM(element_str(k)), bm_sap_2_heartw(:,:,k,1))
       CALL xios_orchidee_send_field('BM_ALLOC_HEART_BE'//TRIM(element_str(k)), bm_sap_2_heartw(:,:,k,2))
    ENDDO

    CALL xios_orchidee_send_field("HERBIVORES",herbivores)
    CALL xios_orchidee_send_field("LEAF_AGE",leaf_meanage(:,:,ileaf))
    CALL xios_orchidee_send_field("SAP_AB_AGE",leaf_meanage(:,:,isapabove))
    CALL xios_orchidee_send_field("ROOT_AGE",leaf_meanage(:,:,iroot))    

    ! Write mean leaf age and time constant of probability of a leaf to be eaten by a herbivore 
    ! to the stomate output file.
    CALL histwrite_p (hist_id_stomate, 'LEAF_AGE', itime, &
         leaf_meanage(:,:,ileaf), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'SAP_AB_AGE', itime, &
         leaf_meanage(:,:,isapabove), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'ROOT_AGE', itime, &
         leaf_meanage(:,:,iroot), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'HERBIVORES', itime, &
         herbivores, npts*nvm, horipft_index)

    IF (printlev>=4) WRITE(numout,*) 'Leaving turnover'

  END SUBROUTINE turn

END MODULE stomate_turnover
