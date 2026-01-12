! =================================================================================================================================
! MODULE       : grassland_cutting
!
! CONTACT      : orchidee-help _at_ ipsl.jussieu.fr
!
! LICENCE      : IPSL (2006)
! This software is governed by the CeCILL licence see
! ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF       This module execute grassland harvest process, including
!! calculate grass biomass after harvest, carbon/nitrogen loss during 
!! harvest
!!
!!\n DESCRIPTION : None
!!
!! RECENT CHANGE(S) : None
!!
!! REFERENCE(S) : None
!!
!! \n
!_
!================================================================================================================================

! Excution of cutting
MODULE grassland_cutting

  ! modules used:

  USE grassland_fonctions
  USE grassland_constantes
  USE stomate_data
  USE constantes
  USE pft_parameters
  USE ioipsl
  ! USE Balances 

  IMPLICIT NONE
  REAL, SAVE                 :: mem_tjulian
  REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: wshtotsumprev

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!! NEW CUTTING SUBROUTINE for CNP version
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Excute the same process for grids where flag_cutting=1
  SUBROUTINE cutting_spa_cnp(&
     npts              , &
     tjulian           , &
     flag_cutting      , &
     wshtotcutinit     , &
     biomass           , &
     loss              , &
     losscnp           , &
     grassharvestcnp   , &
     devstage          , &
     tasum             , &
     tgrowth           , &
     gmean             , &
     regcount          , &
     tlossstart        , &
     wshtotsum         , &
     wshtotsumprev     , &
     tcut              , &
     tcut_modif, &
     flam_left)

    !
    ! 0 declarations
    !

    ! 0.1 input

    INTEGER                              , INTENT(in)    :: npts
    INTEGER                              , INTENT(in)  :: tjulian
    INTEGER, DIMENSION(npts,nvm)         , INTENT(in)    :: flag_cutting
    REAL, DIMENSION(npts,nvm,ncut)       , INTENT(in)    :: wshtotcutinit
    REAL, DIMENSION(npts,nvm,nparts,nelements)     , INTENT(inout)   :: biomass
    REAL, DIMENSION(npts,nvm)            , INTENT(inout)   :: loss
    REAL, DIMENSION(npts,nvm,nparts,nelements)  , INTENT(inout)   :: losscnp
    REAL, DIMENSION(npts,nvm,nelements)  , INTENT(inout)   :: grassharvestcnp
    REAL, DIMENSION(npts,nvm)            , INTENT(inout) :: devstage
    ! stade of developpment of plant (-)
    REAL, DIMENSION(npts,nvm)            , INTENT(inout)   :: tasum
    REAL, DIMENSION(npts,nvm)            , INTENT(inout)   :: tgrowth
    REAL, DIMENSION(npts,nvm,ngmean)     , INTENT(inout)   :: gmean
    INTEGER, DIMENSION(npts,nvm)         , INTENT(inout) :: regcount
    ! number of cut realized(-)
    REAL, DIMENSION(npts,nvm)            , INTENT(inout)   :: tlossstart
    REAL, DIMENSION(npts,nvm)     , intent(inout) :: wshtotsum
    ! yield = total (substrate + structural) shoot dry matter (kg DM m-2)
    REAL, DIMENSION(npts,nvm)     , intent(inout) :: wshtotsumprev
    REAL, DIMENSION(npts,nvm,ncut), INTENT(inout)   :: tcut
    ! date of cutting
    REAL, DIMENSION(npts,nvm,ncut), INTENT(inout)   :: tcut_modif
    REAL, INTENT(in)    :: flam_left

    ! local
    REAL, PARAMETER       :: yieldloss    = 0.05
    INTEGER                     :: h, i, j
    REAL, DIMENSION(npts,nvm) :: wlam          ! leaf mass (kg m-2)
    REAL, DIMENSION(npts,nvm) :: wst           ! stem mass (kg m-2)
    REAL, DIMENSION(npts,nvm) :: wear          ! ear mass (kg m-2)
    REAL, DIMENSION(npts,nvm) :: wroot          ! ear mass (kg m-2)
    REAL, DIMENSION(npts,nvm) :: wshtotreg  ! yield of regrowth
    REAL, DIMENSION(npts,nvm) :: prev_wlam          ! saved leaf mass (kg m-2)
    REAL, DIMENSION(npts,nvm) :: prev_wst           ! saved stem mass (kg m-2)
    REAL, DIMENSION(npts,nvm) :: prev_wear          ! saved ear mass (kg m-2)
    REAL, DIMENSION(npts,nvm) :: prev_wroot          ! saved ear mass (kg m-2)
    REAL, DIMENSION(npts,nvm) :: frac_wlam          ! left fraction leaf mass (kg m-2)
    REAL, DIMENSION(npts,nvm) :: frac_wst           ! left fraction stem mass (kg m-2)
    REAL, DIMENSION(npts,nvm) :: frac_wear          ! left fraction ear mass (kg m-2)
    REAL, DIMENSION(npts,nvm) :: frac_Hlabile
    REAL, DIMENSION(npts,nvm) :: wshcutinit

    DO j=2,nvm
!      WHERE(flag_cutting(:,j) .EQ. 1)
        !transform biomass to dry matter kgDM/m2
        wlam(:,j) = ( biomass(:,j,ileaf,icarbon)/(1000*CtoDM) )
        wst(:,j)  = ( biomass(:,j,isapabove,icarbon)/(1000*CtoDM) )
        wear(:,j) = ( biomass(:,j,ifruit,icarbon)/(1000*CtoDM) )
        wroot(:,j) = ( biomass(:,j,iroot,icarbon)/(1000*CtoDM) )
        prev_wlam(:,j) = wlam(:,j)
        prev_wst(:,j) = wst(:,j)
        prev_wear(:,j) = wear(:,j)
        prev_wroot(:,j) = wroot(:,j)
!      END WHERE

      ! wshcutinit is the residual shoot dry matter (kgDM/m2) after cut
      ! it depends on the available shoot (leaf + stem) and wshtotcutinit
      ! wshtotcutinit is read from input file ~0.1 - 0.15 kgDM/m2
      DO i=1,npts
        IF (flag_cutting(i,j) .EQ. 1) THEN
          WRITE(numout,*) 'execute grassland cut',tjulian, i, j
          wshcutinit(i,j) = MIN (wlam(i,j)+wst(i,j), wshtotcutinit(i,j,regcount(i,j) + 1))
          IF ((wlam(i,j)+wst(i,j)) .GT. min_stomate) THEN
            ! start to cut
            ! all ear was cut
            wear(i,j) = zero
            ! leafs was cut in priority in practice due to its location (higher)
            ! Idealy, there are 10% of leaf and 90% of stem in biomass after cut
            IF (wst(i,j) .GE. ((un-flam_left)*wshcutinit(i,j)) .AND. &
                wlam(i,j) .GE. (flam_left*wshcutinit(i,j)) ) THEN
              wst(i,j) = (un-flam_left)*wshcutinit(i,j)
              wlam(i,j) = flam_left*wshcutinit(i,j) 
            ! while wst can not higher than 0.9*wshcutinit(i,j)
            ! more leaf need to be left 
            ELSE IF (wst(i,j) .LT. ((un-flam_left)*wshcutinit(i,j)) ) THEN
              wst(i,j) = wst(i,j)
              wlam(i,j) = wshcutinit(i,j) - wst(i,j)
            ! while wst >= 0.9*wshcutinit(i,j) 
            ! but wlam can not higher than 0.1*wshcutinit(i,j)
            ! more wst left
            ELSE ! (wst(i,j) .GE. (0.9*wshcutinit(i,j)) .AND.
                 ! wlam(i,j) .LT. (0.1*wshcutinit(i,j)) ) THEN
              wlam(i,j) = wlam(i,j)
              wst(i,j) = wshcutinit(i,j) - wlam(i,j)
            END IF

            IF (wlam(i,j) .LT. -min_stomate .OR. &
               wst(i,j) .LT. -min_stomate .OR. &
               wear(i,j) .LT. -min_stomate) THEN
              WRITE(numout,*)  'ERROR cutting negative residue',i,j, &
                wlam(i,j),wst(i,j),wear(i,j)
            END IF

            ! calculate labile pool that been harvested
            ! according to the fraction of total harvested biomass in total
            ! biomass
            frac_Hlabile(i,j) = (wlam(i,j)+wst(i,j)+wear(i,j)+wroot(i,j))/ &
                 (prev_wlam(i,j)+prev_wst(i,j)+prev_wear(i,j)+prev_wroot(i,j))
            IF (frac_Hlabile(i,j) .GT. un+min_stomate .OR. &
                frac_Hlabile(i,j) .LT. -min_stomate) THEN
              WRITE(numout,*)  'ERROR cutting negative Hlabile',i,j, &
                frac_Hlabile(i,j)
            END IF
            losscnp(i,j,ilabile,:) = biomass(i,j,ilabile,:)*(un-frac_Hlabile(i,j))*yieldloss
            grassharvestcnp(i,j,:) = grassharvestcnp(i,j,:) + &
                           biomass(i,j,ilabile,:)*(un-frac_Hlabile(i,j))* (un-yieldloss)
            biomass(i,j,ilabile,:) = biomass(i,j,ilabile,:) * frac_Hlabile(i,j)
                      
  
 
            ! calculate loss losscnp and biomass left
            ! dry matter loss
            loss(i,j) = MAX(zero,( (prev_wlam(i,j)+prev_wst(i,j)+prev_wear(i,j)) - &
                                  (wlam(i,j)+wst(i,j)+wear(i,j)) )*yieldloss + &
                                   losscnp(i,j,ilabile,icarbon)/(1000*CtoDM) )            

            ! same CNP fraction has been cut
            frac_wear(i,j) = zero
            losscnp(i,j,ifruit,:) = biomass(i,j,ifruit,:)*yieldloss
            grassharvestcnp(i,j,:) = grassharvestcnp(i,j,:) + &
                           biomass(i,j,ifruit,:)* (un-yieldloss)
            biomass(i,j,ifruit,:) = zero

            IF (prev_wlam(i,j) .GT. min_stomate) THEN
              frac_wlam(i,j) = wlam(i,j)/prev_wlam(i,j)
              losscnp(i,j,ileaf,:) = biomass(i,j,ileaf,:) * (un - frac_wlam(i,j))*yieldloss
              grassharvestcnp(i,j,:) = grassharvestcnp(i,j,:) + &
                             ( biomass(i,j,ileaf,:) * (un - frac_wlam(i,j)) * &
                               (un-yieldloss) )
              biomass(i,j,ileaf,:) = biomass(i,j,ileaf,:) * frac_wlam(i,j)
            END IF
            IF (prev_wst(i,j) .GT. min_stomate) THEN
              frac_wst(i,j) = wst(i,j)/prev_wst(i,j)
              losscnp(i,j,isapabove,:) = biomass(i,j,isapabove,:) * (un - frac_wst(i,j))*yieldloss
              grassharvestcnp(i,j,:) = grassharvestcnp(i,j,:) + &
                             ( biomass(i,j,isapabove,:) * (un - frac_wst(i,j)) * &
                               (un-yieldloss) )
              biomass(i,j,isapabove,:) = biomass(i,j,isapabove,:) * frac_wst(i,j)
            END IF
            ! yield of regrowth
            wshtotreg(i,j) = MAX(zero,loss(i,j)/yieldloss*(un-yieldloss))
            wshtotsum(i,j) = wshtotsum(i,j) + wshtotreg(i,j)
            wshtotsumprev(i,j) = wshtotsum(i,j)

          END IF ! enough wlam+wst

          ! reset status variables
          IF ((f_autogestion .EQ. 1) .OR. (f_postauto .NE. 0) .OR. &
           & (f_autogestion .EQ. 3) .OR. (f_autogestion .EQ. 4)) THEN
            tcut(i,j,regcount(i,j))       = tjulian
            tcut_modif(i,j,regcount(i,j)) = tjulian
          END IF

          tlossstart(i,j)  = tjulian
         
          IF ((devstage(i,j) .LT. devsecond) .AND. (regcount(i,j) .EQ. 1)) THEN
            devstage(i,j) = 0.0
          ELSE
            ! do change devstage if not set to 0?
! JC MOD024
            devstage(i,j) = 2.0
!            devstage(i,j) = devstage(i,j)
          END IF

          ! reinitialize accumulated temperature (>5 degree)
          tasum(i,j)  = 0.0

          ! time of growthing = 0
          tgrowth(i,j)  = 0.0          

          ! annulation gmean
          gmean(i,j,:) = 0.0

          ! increase count number of cut
          regcount(i,j)  = regcount(i,j)  + 1

        END IF ! flag_cutting
      END DO ! npts

    END DO ! nvm
  END SUBROUTINE cutting_spa_cnp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!! CUTTING
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Excute the same process for grids where flag_cutting=1
  SUBROUTINE cutting_spa(&
     npts              , &
     tjulian           , &
     flag_cutting      , &
     wshtotcutinit     , &
     lcutinit          , &
     wsh               , &
     wshtot            , &
     wr                , &
     c                 , &
     n                 , &
     napo              , &
     nsym              , &
     fn                , &
     t                 , &
     nel               , &
     biomass           , &
     devstage          , &
     regcount          , &
     gmean             , &
     wc_frac           , &
     wnapo             , &
     wnsym             , &
     wgn               , &
     tasum             , &
     tgrowth           , &
     loss              , &
     lossc             , &
     lossn             , &
     tlossstart        , &
     lai               , &
     tcut              , &
     tcut_modif        , &
     wshtotsum         , &
     controle_azote_sum)         

    ! 
    ! 0 declarations
    !

    ! 0.1 input

    INTEGER                              , INTENT(in)    :: npts
    INTEGER                              , INTENT(in)  :: tjulian
    INTEGER, DIMENSION(npts,nvm)         , INTENT(in)    :: flag_cutting
    REAL, DIMENSION(npts,nvm,ncut)       , INTENT(in)    :: wshtotcutinit
    REAL, DIMENSION(npts,nvm,ncut)       , INTENT(in)    :: lcutinit
    REAL, DIMENSION(npts,nvm)            , INTENT(in)    :: wsh  
    ! total dry mass of shoots
    REAL, DIMENSION(npts,nvm)            , INTENT(in)    :: wshtot   
    REAL, DIMENSION(npts,nvm)            , INTENT(in)    :: wr     
    REAL, DIMENSION(npts,nvm)            , INTENT(in)    :: c  
    ! substrat C concentration(kg C/kg)
    REAL, DIMENSION(npts,nvm)            , INTENT(in)    :: n   
    ! substrat N concentration (kg N/kg)
    REAL, DIMENSION(npts,nvm)            , INTENT(in)    :: napo 
    ! n concentration of apoplast (kg N m-2)
    REAL, DIMENSION(npts,nvm)            , INTENT(in)    :: nsym 
    ! n concentration of symplast (kg N m-2)
    REAL, DIMENSION(npts,nvm)            , INTENT(in)    :: fn  
    ! structral N concentration kgN/kg
    INTEGER                              , INTENT(in)    :: t   
    ! time (d)
    REAL, DIMENSION(npts,nvm)            , INTENT(in)    :: nel  
    ! net lactation energy (mj/kg)
    REAL, DIMENSION(npts,nvm,nparts,nelements)     , INTENT(inout)   :: biomass  
    REAL, DIMENSION(npts,nvm)            , INTENT(inout) :: devstage
    ! stade of developpment of plant (-)
    INTEGER, DIMENSION(npts,nvm)         , INTENT(inout) :: regcount
    ! number of cut realized(-)
    REAL, DIMENSION(npts,nvm,ngmean)     , INTENT(out)   :: gmean
    REAL, DIMENSION(npts,nvm)            , INTENT(out)   :: wc_frac
    REAL, DIMENSION(npts,nvm)            , INTENT(out)   :: wnapo
    REAL, DIMENSION(npts,nvm)            , INTENT(out)   :: wnsym
    REAL, DIMENSION(npts,nvm)            , INTENT(out)   :: wgn
    REAL, DIMENSION(npts,nvm)            , INTENT(out)   :: tasum
    REAL, DIMENSION(npts,nvm)            , INTENT(out)   :: tgrowth
    REAL, DIMENSION(npts,nvm)            , INTENT(out)   :: loss
    REAL, DIMENSION(npts,nvm)            , INTENT(out)   :: lossc
    REAL, DIMENSION(npts,nvm)            , INTENT(out)   :: lossn
    REAL, DIMENSION(npts,nvm)            , INTENT(out)   :: tlossstart
    REAL, DIMENSION(npts,nvm)            , INTENT(inout)   :: lai 
    ! leaf area index of an individual plant
    REAL, DIMENSION(npts,nvm,ncut), INTENT(out)   :: tcut     
    ! date of cutting
    REAL, DIMENSION(npts,nvm,ncut), INTENT(out)   :: tcut_modif
    REAL, DIMENSION(npts,nvm)     , intent(inout) :: wshtotsum  
    ! yield = total (substrate + structural) shoot dry matter (kg DM m-2)
    REAL, DIMENSION(npts,nvm), INTENT(inout)      :: controle_azote_sum

    ! local
    REAL, PARAMETER       :: yieldloss    = 0.05
    REAL, DIMENSION(npts,nvm) :: proportion_wsh, proportion_llam
    INTEGER                     :: h, i, j
    REAL, DIMENSION(npts,nvm) :: wlam          ! leaf mass (kg m-2)
    REAL, DIMENSION(npts,nvm) :: wst           ! stem mass (kg m-2)
    REAL, DIMENSION(npts,nvm) :: wear          ! ear mass (kg m-2)
    REAL, DIMENSION(npts,nvm) :: wshtotreg  ! yield of regrowth
    REAL, DIMENSION(npts,nvm) :: w_postecut
    REAL, DIMENSION(npts,nvm) :: wshcutinit
    DO j=2,nvm
      WHERE(flag_cutting(:,j) .EQ. 1) 
        !transform biomass to dry matter kgDM/m2
        wlam(:,j) = ( biomass(:,j,ileaf,icarbon)/(1000*CtoDM) ) / &
                 & (1.0 + (mc /12.0)*c(:,j) + (mn /14.0)*n(:,j) )
        wst(:,j)  = ( biomass(:,j,isapabove,icarbon)/(1000*CtoDM) ) / &
                 & (1.0 + (mc /12.0)*c(:,j) + (mn /14.0)*n(:,j) )
        wear(:,j) = ( biomass(:,j,ifruit,icarbon)/(1000*CtoDM) ) / &
                 & (1.0 + (mc /12.0)*c(:,j) + (mn /14.0)*n(:,j) )
      END WHERE
!070904 AIG definitions
!
! # wshcutinit: residual structural shoot dry mass
! the less, it equals simulated structural shoot dry mass in input file (when it is less
! than the value in input file)
! in such a case, wshcutinit=wsh=wlam+wstem and proportion_wsh=1 
! the more, it equals residual structural shoot dry mass
! # wshtotcutinit: residual total shoot dry mass (in input file)
!   in the CALL of the subroutine cutting_spa in module 'plantes', wshtotcutinit_Pl 
!   replaces wshtotcutinit which is an argument of the subroutine cutting_spa
!   wshtotcutinit_Pl >=0.12
! # proportion_wsh: proportion of the residual in comparison with the structural shoot 
! dry mass simulated
! # proportion_llam : proportion of the residual in comparison with the structural leaf
! aera simulated

      DO i=1,npts
        IF (flag_cutting(i,j) .EQ. 1) THEN
          wshcutinit(i,j) = MIN (wsh(i,j), wshtotcutinit(i,j,regcount(i,j) + 1) / &
                           & (1.0 + (mc/12.0)*c(i,j) + (mn/14.0)*n(i,j)))
        END IF
      END DO
      WHERE (flag_cutting(:,j) .EQ. 1) 
        WHERE (ABS(wlam(:,j) + wst(:,j)) .GT. 10e-15) 
          proportion_wsh(:,j)  = wshcutinit(:,j)/(wlam(:,j) + wst(:,j))
        ELSEWHERE
        !JCmodif 070904 AIG 
        ! I think there is one error here
        ! If ABS(wlam + wst)< minimum_vital, there is not enough shoot biomass for the cut
        ! to occur, so that proportion_wsh =1 and any dry mass remain the same value(see after)
            !proportion_wsh(:)  = 0.0
          proportion_wsh(:,j)  = 1.0
        !ENDJCmodif 070904 AIG end
        END WHERE
      END WHERE

    ! 070725 AIG confirm
    !-----------------------------
    ! calculation of new biomasses, LAI and substrate concentrations for each component
    ! of different ages after cutting
    !-----------------------------

      WHERE (flag_cutting(:,j) .EQ. 1) 

        wlam(:,j)   = proportion_wsh(:,j) * wlam(:,j)
        wst(:,j)    = proportion_wsh(:,j) * wst(:,j)
        wear(:,j)   = 0.0
        wc_frac(:,j)     = c(:,j)    * (wr(:,j) + wshcutinit(:,j))
        wnapo(:,j)  = napo(:,j) * (wr(:,j) + wshcutinit(:,j))
        wnsym(:,j)  = nsym(:,j) * (wr(:,j) + wshcutinit(:,j))
        wgn(:,j)    = fn(:,j)   * (wr(:,j) + wshcutinit(:,j))
        ! assumed the composition (leaf/stem) in rest of the grass
        w_postecut(:,j) = wlam(:,j) + wst(:,j)
        wlam(:,j) = 0.1 * w_postecut(:,j)
        wst(:,j) = 0.9 * w_postecut(:,j)

        WHERE((devstage(:,j) .LT. devsecond) .AND. (regcount(:,j) .EQ. 1))
            devstage(:,j) = 0.0
        ELSEWHERE 
!            devstage(:,j) = 2.0
            devstage(:,j) = devstage(:,j)
        END WHERE

        ! 070725 AIG confirm
        !-----------------------------
        ! reinitialization of air sum temperature and counter for regrowth after last cutting
        !-----------------------------

        ! reinitialize accumulated temperature (>5 degree)
        tasum(:,j)  = 0.0
        
        ! time of growthing = 0
        tgrowth(:,j)  = 0.0
      END WHERE

       ! 070725 AIG confirm
        !-----------------------------
        ! calculations of 
        ! - biomass and substrate losses
        ! - total yield(without in or_pa) and LAI
        ! after last cutting   
        !-----------------------------

      DO i=1,npts

        IF (flag_cutting(i,j) .EQ. 1) THEN
          !ernteverluste, ms 1999
          !constant  yieldloss  = 0.05
          loss(i,j)        = MAX (0.0, wshtot(i,j) - wshtotcutinit(i,j,regcount(i,j)+1)) * yieldloss 
          ! lossc(i)       = (c(i) + fcsh) * (wsh(i) - wshcutinit(i)) * yieldloss 
          lossc(i,j)       = loss(i,j)*CtoDM!*8
          lossn(i,j)       = loss(i,j)* (n(i,j) + fn(i,j))
          ! lossn(i)       = (n(i) + fn(i))   * (wsh(i) - wshcutinit(i)) * yieldloss 
          tlossstart(i,j)  = t 

          ! yield of regrowth
          wshtotreg(i,j) = MAX(0.0 ,wshtot(i,j) - wshtotcutinit(i,j,regcount(i,j)+1)) * (1 - yieldloss)
        END IF
      END DO

      ! RUN WITH GRASS AUTOGESTION (SATURANT OR NONLIMITANT)
      ! it is important to set it for auto cut
      IF ((f_autogestion .EQ. 1) .OR. (f_postauto .NE. 0) .OR. &
        & (f_autogestion .EQ. 3) .OR. (f_autogestion .EQ. 4)) THEN
        DO i=1,npts
          IF  (flag_cutting(i,j) .EQ. 1) THEN   
            tcut(i,j,regcount(i,j))       = tjulian
            tcut_modif(i,j,regcount(i,j)) = tjulian
          END IF
        END DO 
      END IF

      DO i=1,npts
        IF (flag_cutting(i,j) .EQ. 1) THEN
          ! annulation gmean
          gmean(i,j,:) = 0.0

          ! increase count number of cut
          regcount(i,j)  = regcount(i,j)  + 1 
        END IF
      END DO

        ! 070725 AIG confirm
        !-----------------------------
        ! calculations of yield, total LAI, total carbon and nitrogen and nel 
        ! after regrowth
        !-----------------------------

      WHERE(flag_cutting(:,j) .EQ. 1) 
        wshtotsum(:,j)     = wshtotsum(:,j) + wshtotreg(:,j)
        biomass(:,j,ileaf,icarbon)     = (wlam(:,j) * 1000*CtoDM) * (1.0 + (mc /12.0)*c(:,j) + &
                                        & (mn /14.0)*n(:,j) )
        biomass(:,j,isapabove,icarbon) = (wst(:,j)  * 1000*CtoDM) * (1.0 + (mc /12.0)*c(:,j) + &
                                        & (mn /14.0)*n(:,j) )
        biomass(:,j,ifruit,icarbon)    = (wear(:,j) * 1000*CtoDM) * (1.0 + (mc /12.0)*c(:,j) + &
                                        & (mn /14.0)*n(:,j) )

      END WHERE
      IF (f_autogestion .LT. 2) THEN
        controle_azote_sum(:,j)=controle_azote_sum(:,j)+wshtotsum(:,j)-wshtotsumprev(:,j)
      ENDIF

        wshtotsumprev(:,j)=wshtotsum(:,j)
    END DO ! nvm

  END SUBROUTINE cutting_spa

END MODULE grassland_cutting
