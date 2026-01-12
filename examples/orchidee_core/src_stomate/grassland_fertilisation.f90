! =================================================================================================================================
! MODULE       : grassland_fertilisation
!
! CONTACT      : orchidee-help _at_ ipsl.jussieu.fr
!
! LICENCE      : IPSL (2006)
! This software is governed by the CeCILL licence see
! ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF       This module execute grassland fertilization process,
!! calculate carbon/nitrogen in fertilizer (organic/mineral),
!! spatialize the fertilizer into corresponding grid
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
MODULE grassland_fertilisation 

  USE grassland_fonctions
  USE grassland_constantes
  USE pft_parameters
  USE constantes
  USE ioipsl

  IMPLICIT NONE
  REAL, PARAMETER :: c2Nsolidmanure   = 15.0
  REAL, PARAMETER :: c2Nslurry        = 10.0
!  ! flag for verify that without twice fertilisation at same time
!  LOGICAL     , ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: tfert_verif 
!  LOGICAL     , ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: tfert_verif2  ! idem
!  LOGICAL     , ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: tfert_verif3  ! idem


CONTAINS

  ! two function of fertilization
  ! one is for each time
  ! the other one is for strategy spatialization of fertilization in management


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!  FERTILISATION a chaque pas de temps
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE fertilisation_pas_temps(&
     npts                           , &
     fertcount                      , &
     dt                             , &
     tjulian                        , &
     deltat                         , &
     tfert                          , &
     Nliquidmanure                  , &
     Nslurry                        , &
     Nsolidmanure                   , &
     fcOrganicFertmetabolic         , &
     fcOrganicFertstruct            , &
     fnOrganicFerturine             , &
     fnOrganicFertmetabolic         , &
     fnOrganicFertstruct            , &
     fpOrganicFerturine             , &
     fpOrganicFertmetabolic         , &
     fpOrganicFertstruct            , &
     c2nratiostruct                 , &
     Fert_on)

    INTEGER                          , INTENT(in)     :: npts               
     ! Number of spatial points
    INTEGER, DIMENSION(npts,nvm)     , INTENT(inout)  :: fertcount           
    ! counter for fertilizer application (-)
    REAL                             , INTENT(in)     :: dt
    INTEGER                          , INTENT(in)     :: tjulian              
    ! Julian day (-)
    REAL                             , INTENT(in)     :: deltat
    REAL, DIMENSION(npts,nvm, nfert) , INTENT(in)     :: tfert               
    ! zeitpunkt der duengung h (1,..,nfert) (d)
    REAL, DIMENSION(npts,nvm, nfert) , INTENT(in)     :: Nliquidmanure        
    ! Nitrogen in liquid manure (kg N m-2) (lisier)
    REAL, DIMENSION(npts,nvm, nfert) , INTENT(in)     :: Nslurry              
    ! Nitrogen in slurry (kg N m-2) (boues)
    REAL, DIMENSION(npts,nvm, nfert) , INTENT(in)     :: Nsolidmanure         
    ! Nitrogen in solid manure (kg N m-2)
    REAL, DIMENSION(npts,nvm)        , INTENT(out)    :: fcOrganicFertmetabolic 
    ! Carbon flux from slurry and manure to metabolic SOM pool (kg c m-2 d-1)
    REAL, DIMENSION(npts,nvm)        , INTENT(out)    :: fcOrganicFertstruct  
    ! Carbon flux from slurry and manure to structural SOM pool (kg c m-2 d-1)
    REAL, DIMENSION(npts,nvm)        , INTENT(out)    :: fnOrganicFerturine   
    ! Nitrogen flux from Organic ferilization to urine N pool (kg N m-2 d-1)
    REAL, DIMENSION(npts,nvm)        , INTENT(out)    :: fnOrganicFertmetabolic
    ! Nitrogen flux from Organic ferilization to metabolic soil N pool (kg N m-2 d-1)
    REAL, DIMENSION(npts,nvm)        , INTENT(out)    :: fnOrganicFertstruct
    ! Nitrogen flux from Organic ferilization to struct soil N pool (kg N m-2 d-1)
    REAL, DIMENSION(npts,nvm)        , INTENT(out)    :: fpOrganicFerturine
    ! Nitrogen flux from Organic ferilization to urine N pool (kg N m-2 d-1)
    REAL, DIMENSION(npts,nvm)        , INTENT(out)    :: fpOrganicFertmetabolic
    ! Nitrogen flux from Organic ferilization to metabolic soil N pool (kg N m-2
    ! d-1)
    REAL, DIMENSION(npts,nvm)        , INTENT(out)    :: fpOrganicFertstruct
    ! Nitrogen flux from Organic ferilization to struct soil N pool (kg N m-2
    ! d-1)
    REAL, DIMENSION(npts,nvm)        , INTENT(in)     :: c2nratiostruct

    REAL, DIMENSION(npts,nvm)        , INTENT(inout)     :: Fert_on
    ! variables locales
    INTEGER :: i,j
    REAL, DIMENSION(npts,nvm) :: fnfert_resultat1   
    REAL, DIMENSION(npts,nvm) :: fnfert_resultat2  
    REAL, DIMENSION(npts,nvm) :: fnfert_resultat3  

    fnfert_resultat1 = 0.0
    fnfert_resultat2 = 0.0
    fnfert_resultat3 = 0.0

    !fertilization
    !application rates
          fcOrganicFertmetabolic(:,:)  = 0.0
          fcOrganicFertstruct(:,:)     = 0.0
          fnOrganicFerturine(:,:)      = 0.0
          fnOrganicFertmetabolic(:,:)  = 0.0
          fnOrganicFertstruct(:,:)     = 0.0
          fpOrganicFerturine(:,:)      = 0.0
          fpOrganicFertmetabolic(:,:)  = 0.0
          fpOrganicFertstruct(:,:)     = 0.0
    ! engrais de ferme ! valeu estimated based on fertilisation 
    ! in chambers ms 1999
    DO j = 2, nvm
      DO i= 1, npts
        IF (fertcount(i,j) .GT. 0) THEN
          CALL fnfert(fnfert_resultat1(i,j), tjulian, deltat, &
             tfert(i,j,fertcount(i,j)), tapplvg, Nliquidmanure(i,j,fertcount(i,j)))
          CALL fnfert(fnfert_resultat2(i,j), tjulian, deltat, &
             tfert(i,j,fertcount(i,j)), tapplka, Nslurry(i,j,fertcount(i,j)))
          CALL fnfert(fnfert_resultat3(i,j), tjulian, deltat, &
             tfert(i,j,fertcount(i,j)), tapplmist, Nsolidmanure(i,j,fertcount(i,j)))

         
          fcOrganicFertmetabolic(i,j)  =  &
             fvgcmetabolic * (1.0 - fvgnurine) * c2Nslurry * fnfert_resultat1(i,j) + &
             fkacmetabolic * (1.0 - fkanurine) * c2Nslurry * fnfert_resultat2(i,j) + &
             fmistcmetabolic * c2Nsolidmanure * fnfert_resultat3(i,j)

          fcOrganicFertstruct(i,j)  =  &
             (1.0 - fvgcmetabolic) * (1.0 - fvgnurine) * c2Nslurry * fnfert_resultat1(i,j) + &
             (1.0 - fkacmetabolic) * (1.0 - fkanurine) * c2Nslurry * fnfert_resultat2(i,j) + &
             (1.0 - fmistcmetabolic) * c2Nsolidmanure *fnfert_resultat3(i,j)

          fnOrganicFerturine(i,j)  =  fvgnurine*fnfert_resultat1(i,j) + &
                                      fkanurine*fnfert_resultat2(i,j)
          
          ! in case negative fnOrganicFertXXX appears
          fnOrganicFertstruct(i,j)  = fcOrganicFertstruct(i,j)/c2nratiostruct(i,j)
    
          fnOrganicFertmetabolic(i,j)  = &
             (1.0 - fvgnurine)*fnfert_resultat1(i,j) + &
             (1.0 - fkanurine)*fnfert_resultat2(i,j) + &
             fnfert_resultat3(i,j) - fcOrganicFertstruct(i,j)/c2nratiostruct(i,j)  !ms 1999
          IF (fnOrganicFertmetabolic(i,j) .LT. -1e-10) THEN
            WRITE(numout,*) 'ERROR negative fnOrganicFertmetabolic',fnOrganicFertmetabolic(i,j)
            STOP 'ERROR negative fnOrganicFertmetabolic'
          END IF

          ! Phosphorus excreted via urine by ruminants is usually negligible
          ! because a possible surplus of digested P is recirculated to the
          ! rumen and from there excreted in faeces instead of in urine.
          ! http://www.fao.org/WAIRDOCS/LEAD/X6113E/x6113e06.htm
          fpOrganicFerturine(i,j)      = 0.0
          ! 10:3 of N:P ratio is assumed 
          ! according to US states data
          ! https://www.epa.gov/nutrient-policy-data/estimated-animal-agriculture-nitrogen-and-phosphorus-manure
          fpOrganicFertmetabolic(i,j)  = 0.3 * fnOrganicFertmetabolic(i,j)
          fpOrganicFertstruct(i,j)     = 0.3 * fnOrganicFertstruct(i,j)

          Fert_on(i,j) = (fnOrganicFerturine(i,j) + fnOrganicFertmetabolic(i,j) + &
                         fnOrganicFertstruct(i,j))
 
         ELSE
          fcOrganicFertmetabolic(i,j)  = 0.0
          fcOrganicFertstruct(i,j)     = 0.0
          fnOrganicFerturine(i,j)      = 0.0
          fnOrganicFertmetabolic(i,j)  = 0.0
          fnOrganicFertstruct(i,j)     = 0.0
          fpOrganicFerturine(i,j)      = 0.0
          fpOrganicFertmetabolic(i,j)  = 0.0
          fpOrganicFertstruct(i,j)     = 0.0
          Fert_on(i,j) = 0.0
         ENDIF
      END DO ! i npts
    END DO ! j nvm

  END SUBROUTINE fertilisation_pas_temps


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!! FNFERT 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE Fnfert(&
!     npts           , &
     Fnfertform     , &
     tjulianform    , &
     deltatform     , &
     tfertform      , &
     tapplform      , &
     Nfertform)

!    INTEGER (i_std)                   , INTENT(in)  :: npts
    INTEGER                 , INTENT(in)  :: tjulianform
    REAL                 , INTENT(in)  :: deltatform
    REAL                 , INTENT(in)  :: tapplform
!    REAL, DIMENSION(npts), INTENT(in)  :: tfertform
!    REAL, DIMENSION(npts), INTENT(in)  :: Nfertform
!    REAL, DIMENSION(npts), INTENT(out) :: Fnfertform
    REAL, INTENT(in)  :: tfertform
    REAL, INTENT(in)  :: Nfertform
    REAL, INTENT(out) :: Fnfertform

    ! variables locales :
!    INTEGER :: i


    IF (deltatform .EQ. 0.0) THEN
        Fnfertform = 0.0
    ELSE
      IF (tfertform-0.5 .LE. REAL(tjulianform) .AND. &
             tfertform+0.9 .GE. REAL(tjulianform)) THEN
        Fnfertform = Nfertform
      ELSE
        Fnfertform = 0.0
      END IF
! JC comment
! the original statement is for multiple application within one day for PaSim
!      WHERE ((tjulianform + deltatform) .LE. tfertform(:)) 
!
!        Fnfertform(:) = 0.0
!
!      ELSEWHERE ((tjulianform .LE. tfertform(:)) .AND. ((tjulianform + deltatform) .GT. tfertform(:))) 
!
!        Fnfertform(:) = &
!          Nfertform(:)*(MIN((tjulianform + deltatform - tfertform(:)), tapplform)) / &
!          deltatform/tapplform
!
!      ELSEWHERE (tjulianform .LT. (tfertform(:) + tapplform)) 
!
!        Fnfertform(:) = &
!          Nfertform(:)*(MIN((tfertform(:) + tapplform - tjulianform),deltatform)) / &
!          deltatform/tapplform
!
!      ELSEWHERE
!
!        Fnfertform(:) = 0.0
!      END WHERE
    ENDIF

  END SUBROUTINE Fnfert

! JC comment 13Feb2017
! do not need this subroutine anymore
! the processes has been added in the grassland management
!!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  !!!!!! SPATIALISATION DE LA FERTILISATION GEREE PAR MANAGEMENT
!!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!  SUBROUTINE Fertilisation_spa(&
!!     npts                  , &
!!     flag_fertilisation    , &
!!     fertcount_start       , &
!!     tjulian               , &
!!     tfert                 , &
!!     nfertnittotprevyear   , &
!!     nfertammtotprevyear   , &
!!     nfertnit              , &
!!     nfertamm              , &
!!     fertcount             , &
!!     nfertammtot           , &
!!     nfertnittot           , &
!!     nfertammtotyear       , &
!!     nfertnittotyear       , &
!!     controle_azote_sum    , &
!!     controle_azote_sum_mem)
!!
!!    INTEGER (i_std)                         , INTENT(in)    :: npts
!!    INTEGER   , DIMENSION(npts,nvm)      , INTENT(in)    :: flag_fertilisation
!!    REAL, DIMENSION(npts,nvm)      , INTENT(inout)    :: nfertnittotprevyear
!!    REAL, DIMENSION(npts,nvm)      , INTENT(inout)    :: nfertammtotprevyear
!!    REAL, DIMENSION(npts,nvm,nfert), INTENT(in)    :: nfertnit
!!    REAL, DIMENSION(npts,nvm,nfert), INTENT(in)    :: nfertamm
!!    INTEGER, DIMENSION(npts,nvm)      , INTENT(inout) :: fertcount
!!    REAL, DIMENSION(npts,nvm)      , INTENT(inout) :: nfertammtot
!!    REAL, DIMENSION(npts,nvm)      , INTENT(inout) :: nfertnittot
!!    REAL, DIMENSION(npts,nvm)      , INTENT(out)   :: nfertammtotyear
!!    REAL, DIMENSION(npts,nvm)      , INTENT(out)   :: nfertnittotyear
!!    INTEGER, DIMENSION(npts,nvm)      , INTENT(in)    :: fertcount_start
!!    INTEGER                       , INTENT(in)    :: tjulian       ! Julian day (-)
!!    REAL, DIMENSION(npts,nvm,nfert), INTENT(in)    :: tfert
!!    REAL, DIMENSION(npts,nvm)      , INTENT(inout) :: controle_azote_sum
!!    REAL, DIMENSION(npts,nvm)      , INTENT(out)   :: controle_azote_sum_mem
!!
!!    INTEGER :: i,j
!!
!!
!!    IF (blabla_pasim) PRINT *, 'PASIM main grassland : call fertilisation_spa'
!!
!!    DO j = 2, nvm
!!
!!      DO i=1,npts
!!
!!        IF (flag_fertilisation(i,j) .EQ. 1) THEN 
!!          !counter for fertilizer application
!!          fertcount(i,j)  = fertcount(i,j)  + 1
!!
!!          !mineral fertilization
!!          nfertammtot(i,j)  = nfertammtot(i,j)  + nfertamm(i,j,fertcount(i,j))
!!          nfertnittot(i,j)  = nfertnittot(i,j)  + nfertnit(i,j,fertcount(i,j))
!!          
!!          nfertammtotyear(i,j)  = nfertammtot(i,j)  - nfertammtotprevyear(i,j) 
!!          nfertnittotyear(i,j)  = nfertnittot(i,j)  - nfertnittotprevyear(i,j) 
!!
!!        END IF
!!      END DO
!!    END DO
!!    !*****RUN NONLIMITANT
!!    IF (f_nonlimitant .EQ. 1.) THEN
!!      DO j = 2, nvm
!!        DO i=1,npts
!!          IF (flag_fertilisation(i,j) .EQ. 1) THEN
!!
!!            controle_azote_sum_mem(i,j) = controle_azote_sum(i,j)
!!
!!            IF ((tjulian .GE. tfert(i,j,fertcount_start(i,j)))  .AND. &
!!               (tjulian .LE. tfert(i,j,fertcount_start(i,j))+0.9) .AND. &
!!               (.NOT. (tfert_verif3(i,j,fertcount_start(i,j))) )) THEN
!!
!!              tfert_verif3(i,j,fertcount_start(i,j))= .TRUE.
!!              controle_azote_sum(i,j) = 0.
!!
!!            ENDIF
!!          END IF
!!        END DO
!!      END DO
!!    ENDIF
!!
!!  END SUBROUTINE Fertilisation_spa

END MODULE grassland_fertilisation
