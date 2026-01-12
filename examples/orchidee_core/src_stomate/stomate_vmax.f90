! =================================================================================================================================
! MODULE 	: stomate_vmax
!
! CONTACT	: orchidee-help _at_ listes.ipsl.fr
!
! LICENCE      	: IPSL (2006). This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF        calculates the leaf efficiency.
!!	
!!\n DESCRIPTION: None
!!
!! RECENT CHANGE(S): None
!!
!! SVN		:
!! $HeadURL: svn://forge.ipsl.jussieu.fr/orchidee/branches/ORCHIDEE-CN-P/ORCHIDEE/src_stomate/stomate_vmax.f90 $ 
!! $Date: 2019-08-27 16:52:54 +0200 (二, 2019-08-27) $
!! $Revision: 6176 $
!! \n
!_ =================================================================================================================================

MODULE stomate_vmax

  ! modules used:

  USE ioipsl_para
  USE stomate_data
  USE constantes
  USE pft_parameters

  IMPLICIT NONE

  ! private & public routines

  PRIVATE
  PUBLIC vmax, vmax_clear

  ! first call
  LOGICAL, SAVE                                              :: firstcall_vmax = .TRUE.
!$OMP THREADPRIVATE(firstcall_vmax)

CONTAINS

!! ================================================================================================================================
!! SUBROUTINE	: vmax_clear
!!
!>\BRIEF	  Flag setting 
!!
!!\n DESCRIPTION: This subroutine sets flags ::firstcall_vmax, to .TRUE., and therefore activates   
!!		  section 1.1 of the ::vmax subroutine which writes messages to the output. \n
!!		  This subroutine is called at the end of the subroutine ::stomate_clear, in the 
!!		  module ::stomate.
!!
!! RECENT CHANGE(S):None
!!
!! MAIN OUTPUT VARIABLE(S): ::firstcall_vmax
!!
!! REFERENCE(S)  : None 
!!
!! FLOWCHART     : None
!! \n		  
!_ =================================================================================================================================

  SUBROUTINE vmax_clear
    firstcall_vmax=.TRUE.
  END SUBROUTINE vmax_clear



!! ================================================================================================================================
!! SUBROUTINE    : vmax
!!
!>\BRIEF         This subroutine computes vcmax photosynthesis parameters 
!! given optimal vcmax parameter values and a leaf age-related efficiency.
!!
!! DESCRIPTION (functional, design, flags): 
!! Leaf age classes are introduced to take into account the fact that photosynthetic activity depends on leaf age
!! (Ishida et al., 1999). There are \f$nleafages\f$ classes (constant defined in stomate_constants.f90).
!! This subroutine first calculates the new age of each leaf age-class based on fraction of leaf 
!! that goes from one to another class.                                              
!! Then calculation of the new fraction of leaf in each class is performed.      
!! Last, leaf efficiency is calculated for each PFT and for each leaf age class.
!! vcmax is defined as vcmax25 and vjmax_opt weighted by a mean leaf
!! efficiency. vcmax25 is PFT-dependent constants defined in constants_mtc.f90.
!! DSG: vcmax25 is not used anymore but vcmax_opt instead which is a function of leaf nutrient content
!!
!! This routine is called once at the beginning by stomate_var_init and then at each stomate time step by stomateLpj.
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): vcmax
!!
!! REFERENCE(S)	: 
!! - Ishida, A., A. Uemura, N. Koike, Y. Matsumoto, and A. Lai Hoe (1999),
!! Interactive effects of leaf age and self-shading on leaf structure, photosynthetic
!! capacity and chlorophyll fluorescence in the rain forest tree,
!! dryobalanops aromatica, Tree Physiol., 19, 741-747
!!
!! FLOWCHART    : None
!!
!! REVISION(S)	: None
!! \n
!_ ================================================================================================================================

  SUBROUTINE vmax (npts, dt, &
       biomass,              &
       leaf_age, leaf_frac,  &
       vcmax, jmax, nue)

    !
    !! 0. Variable and parameter declaration
    !

    !
    !! 0.1 Input variables
    !
    INTEGER, INTENT(in)                                     :: npts               !! Domain size (unitless)
    REAL, INTENT(in)                                        :: dt                 !! time step of stomate (days)

    REAL, DIMENSION(npts,nvm,nparts,nelements), INTENT(in)  :: biomass            !! 

    !
    !! 0.2 Output variables 
    !
    REAL, DIMENSION(npts,nvm,2), INTENT(out)                  :: vcmax            !! Maximum rate of carboxylation 
                                                                                         !! @tex ($\mu mol m^{-2} s^{-1}$) @endtex
                                                                                         !! third dimension: actual Vmax=1 ; potential Vcmax=2
    REAL, DIMENSION(npts,nvm,2), INTENT(out)                  :: jmax               !! Maximum rate of electron transport 
                                                                                         !! @tex ($\mu mol m^{-2} s^{-1}$) @endtex
                                                                                         !! third dimension: actual Jmax=1 ; potential Jmax=2
    REAL,DIMENSION (npts,nvm), INTENT(out)                  :: nue                !! Nitrogen use Efficiency with impact of leaf age (umol CO2 (gN)-1 s-1) 
    !
    !! 0.3 Modified variables
    !
    REAL, DIMENSION(npts,nvm,nparts,nleafages), INTENT(inout) :: leaf_age         !! plant tissue age for each age class
    REAL, DIMENSION(npts,nvm,nparts,nleafages), INTENT(inout) :: leaf_frac        !! plant tissue fraction for each age class

    !
    !! 0.4 Local variables
    !

    REAL, DIMENSION(npts,nvm,2)                             :: leaf_nutrient      !! leaf N(=1) and P(=2)  [10-3 g(N) g-1(DW)]

    REAL, DIMENSION(npts,nvm,2)                               :: vcmax_opt          !! vcmax_opt derived from leaf nutrient content @tex ($\mu mol m^{-2} s^{-1}$) @endtex
                                                                                         !! third dimension: actual Jmax=1 ; potential Jmax=2
    REAL, DIMENSION(npts,nvm,2)                               :: jmax_opt           !! jmax_opt derived from leaf nutrient content @tex ($\mu mol m^{-2} s^{-1}$) @endtex
                                                                                         !! third dimension: actual Jmax=1 ; potential Jmax=2

    REAL, DIMENSION(npts)                                   :: leaf_efficiency    !! leaf efficiency (vcmax/vcmax25)
                                                                                         !! (unitless)
    REAL, DIMENSION(npts,nvm,nparts,nleafages)              :: d_leaf_frac        !! delta leaf fraction 
    REAL, DIMENSION(npts,nparts,nleafages)                  :: leaf_age_new       !! updated leaf age
    REAL, DIMENSION(npts,nparts)                            :: sumfrac            !! temporary variable for total leaf fraction

    REAL, DIMENSION(npts)                                   :: rel_age            !! relative leaf age (age/critical age)
                                                                                         !! (unitless)
    INTEGER                                                 :: j,m                !! indices (unitless)
    INTEGER                                                 :: ipar               !! indices (unitless)


!_ ================================================================================================================================

    IF (printlev>=3) WRITE(numout,*) 'Entering vmax'

    !
    !! 1 Initialization
    !

    !
    !! 1.1 first call: info about flags and parameters.
    !

    IF ( firstcall_vmax ) THEN
       
       IF (printlev >= 2) THEN
          WRITE(numout,*) 'vmax:'

          WRITE(numout,*) '   > offset (minimum vcmax/vmax_opt):' , vmax_offset
          WRITE(numout,*) '   > relative leaf age at which vmax reaches vcmax_opt:', leafage_firstmax 
          WRITE(numout,*) '   > relative leaf age at which vmax falls below vcmax_opt:', leafage_lastmax
          WRITE(numout,*) '   > relative leaf age at which vmax reaches its minimum:', leafage_old
       END IF
       firstcall_vmax = .FALSE.

    ENDIF

    !
    !! 1.2 initialize output
    !

    vcmax(:,:,:)     = zero
    vcmax_opt(:,:,:) = zero
    jmax(:,:,:)      = zero
    jmax_opt(:,:,:)  = zero
    nue(:,:)       = zero 

    !
    !! 1.3 Derive vcmax_opt  and jmax_opt from leaf nutrient content following Ellsworth et al. (in prep)
    ! 
    DO j=2,nvm  ! loop over PFT
         WHERE ( biomass(:,j,ileaf,icarbon) .GT. min_stomate ) 
            ! MINIMUM NUTRIENT CONCENTRATION
            !! 1.3.0 Convert g(N)/m2(ground) -> 10-3g(N)/g(DW) assuming 50%C in DW 
!max leaf N:            leaf_nutrient(:,j,1) = un/cn_leaf_min(j)/2*1.E3
!act leaf N: 
            leaf_nutrient(:,j,1) = (biomass(:,j,ileaf,initrogen)   / (2.*biomass(:,j,ileaf,icarbon)))*1.E3
!act leaf N END

            leaf_nutrient(:,j,2) =  un/(cn_leaf_min(j)*np_leaf_min(j))/2*1.E3
            !! 1.3.1 Derive potential Vcmax and Jmax from maximum leaf nutrient
! old relationship
            vcmax_opt(:,j,2) =  EXP(4.44904 + (0.3472 *LOG(leaf_nutrient(:,j,2))) + (0.49078*LOG(leaf_nutrient(:,j,1))))
            jmax_opt (:,j,2) =  EXP(5.49435 + (0.37345*LOG(leaf_nutrient(:,j,2))) + (0.41435*LOG(leaf_nutrient(:,j,1))))
! new relationship
!            vcmax_opt(:,j,2) =  EXP(4.810 + (0.363 *LOG(leaf_nutrient(:,j,2))) + (0.402*LOG(leaf_nutrient(:,j,1))))
!            jmax_opt (:,j,2) =  EXP(5.703 + (0.398*LOG(leaf_nutrient(:,j,2))) + (0.340*LOG(leaf_nutrient(:,j,1))))

            ! ACTUAL NUTRIENT CONCENTRATION
            !! 1.3.2 Convert g(N)/m2(ground) -> 10-3g(N)/g(DW) assuming 50%C in DW 
            leaf_nutrient(:,j,1) = (biomass(:,j,ileaf,initrogen)   / (2.*biomass(:,j,ileaf,icarbon)))*1.E3
            leaf_nutrient(:,j,2) = (biomass(:,j,ileaf,iphosphorus) / (2.*biomass(:,j,ileaf,icarbon)))*1.E3
            !! 1.3.3 Derive Vcmax and Jmax from leaf nutrients
            vcmax_opt(:,j,1) =  EXP(4.44904 + (0.3472 *LOG(leaf_nutrient(:,j,2))) + (0.49078*LOG(leaf_nutrient(:,j,1))))
            jmax_opt (:,j,1) =  EXP(5.49435 + (0.37345*LOG(leaf_nutrient(:,j,2))) + (0.41435*LOG(leaf_nutrient(:,j,1))))
! new relationship
!            vcmax_opt(:,j,1) =  EXP(4.810 + (0.363 *LOG(leaf_nutrient(:,j,2))) + (0.402*LOG(leaf_nutrient(:,j,1))))
!            jmax_opt (:,j,1) =  EXP(5.703 + (0.398*LOG(leaf_nutrient(:,j,2))) + (0.340*LOG(leaf_nutrient(:,j,1))))
         ELSEWHERE
            vcmax_opt(:,j,1)     = zero
            jmax_opt (:,j,1)     = zero
            vcmax_opt(:,j,2)     = zero
            jmax_opt (:,j,2)     = zero

            leaf_nutrient(:,j,1) = zero
            leaf_nutrient(:,j,2) = zero
         ENDWHERE
    ENDDO

    !! 1.3.3 Convert Vcmax and Jmax to ORCHIDEE units 
    ! [nmol CO2 (g drymass)-1 s-1] -> [mu mol CO2 m2-1 s-1]
    vcmax_opt(:,:,1) = vcmax_opt(:,:,1) *2/SPREAD(sla(:),DIM=1,NCOPIES=npts) /1.E3
    jmax_opt (:,:,1) = jmax_opt (:,:,1) *2/SPREAD(sla(:),DIM=1,NCOPIES=npts) /1.E3
    vcmax_opt(:,:,2) = vcmax_opt(:,:,2) *2/SPREAD(sla(:),DIM=1,NCOPIES=npts) /1.E3
    jmax_opt (:,:,2) = jmax_opt (:,:,2) *2/SPREAD(sla(:),DIM=1,NCOPIES=npts) /1.E3

    ! As needle leaves have a well-documented lower slope between leaf N and Vcmax than other
    ! broadleaf (Kattge et al.2009; we must correc the Vcmax and Jmax derived
    ! from David's relationship for broadleaves for needle leaves; we use the
    ! ratio of NUE from Kattges paper which is 0.66 (derived from the PFT
    ! averages shown in his paper table 3); This is pretty rough approximation
    ! and needs to be replaced by something more sophisticated as soon as
    ! available:
    DO j = 2,nvm ! Loop over # PFTs
       IF (is_needleleaf(j)) THEN
          vcmax_opt(:,j,:) = vcmax_opt(:,j,:) * 0.66
          jmax_opt (:,j,:) = jmax_opt (:,j,:) * 0.66
      ENDIF
    ENDDO
    
    !
    !! 2 leaf age: general increase and turnover between age classes.
    !

    !
    !! 2.1 increase leaf age
    !
    !! The age of the leaves in each leaf-age-class increases by 1 time step.
    DO m = 1, nleafages ! Loop over # leaf age classes
       DO j = 2,nvm ! Loop over # PFTs
          WHERE ( leaf_frac(:,j,:,m) .GT. min_stomate )

             leaf_age(:,j,:,m) = leaf_age(:,j,:,m) + dt

          ENDWHERE
       ENDDO    ! Loop over # PFTs

    ENDDO   ! Loop over # leaf age classes

    !
    !! 2.2 turnover between leaf age classes
    !     d_leaf_frac(:,:,m) = what leaves m-1 and goes into m
    !

    DO j = 2,nvm   ! Loop over # PFTs

       !! 2.2.1 fluxes

       !! nothing goes into first age class
       d_leaf_frac(:,j,:,1) = zero

       !! for others age classes (what goes from m-1 to m)
       DO m = 2, nleafages 
          !! leaf_timecst is defined in stomate_constants.f90 as the quotient of 
          !! the critical leaf age per the number of age classes.
          !! The critical leaf age is a PFT-dependent constant defined in stomate_constants.f90, 
          !! that represents the leaf life span.
          !! This time constant (leaf_timecst) determines the turnover between the nleafages different leaf age classes
          !! (see section [118] in Krinner et al. (2005)).
          ! JC comments
          ! now I assumed that all tissues have the same life span as leaf
          ! BUT it is not the truth, and subject to change in the future!!!
          d_leaf_frac(:,j,:,m) = leaf_frac(:,j,:,m-1) * dt/leaf_timecst(j)
       ENDDO

       !! 2.2.2 new leaf age in class
       !!       new age = ( old age * (old fraction - fractional loss) + fractional increase * age of the source class ) / new fraction
       !!       The leaf age of the youngest class (m=1) is updated into stomate_alloc          
       leaf_age_new(:,:,:) = zero

       DO m = 2, nleafages-1       ! Loop over age classes
          !! For all age classes except first and last 
          WHERE ( d_leaf_frac(:,j,:,m) .GT. min_stomate .AND. &
            ( leaf_frac(:,j,:,m) + d_leaf_frac(:,j,:,m)- d_leaf_frac(:,j,:,m+1) ) .GT. min_stomate)

             leaf_age_new(:,:,m) = ( ( (leaf_frac(:,j,:,m)- d_leaf_frac(:,j,:,m+1)) * leaf_age(:,j,:,m) )  + &
                  ( d_leaf_frac(:,j,:,m) * leaf_age(:,j,:,m-1) ) ) / &
                  ( leaf_frac(:,j,:,m) + d_leaf_frac(:,j,:,m)- d_leaf_frac(:,j,:,m+1) )

          ENDWHERE
       ENDDO       ! Loop over age classes

       !! For last age class, there is no leaf fraction leaving the class. 

       WHERE ( d_leaf_frac(:,j,:,nleafages) .GT. min_stomate  .AND. &
         ( leaf_frac(:,j,:,nleafages) + d_leaf_frac(:,j,:,nleafages) ) .GT. min_stomate )

          leaf_age_new(:,:,nleafages) = ( ( leaf_frac(:,j,:,nleafages) * leaf_age(:,j,:,nleafages) )  + &
               ( d_leaf_frac(:,j,:,nleafages) * leaf_age(:,j,:,nleafages-1) ) ) / &
               ( leaf_frac(:,j,:,nleafages) + d_leaf_frac(:,j,:,nleafages) )

       ENDWHERE

       DO m = 2, nleafages       ! Loop over age classes

          WHERE ( d_leaf_frac(:,j,:,m) .GT. min_stomate )

             leaf_age(:,j,:,m) = leaf_age_new(:,:,m)

          ENDWHERE

       ENDDO       ! Loop over age classes

       !! 2.2.3 calculate new fraction

       DO m = 2, nleafages       ! Loop over age classes

          ! where the change comes from
          leaf_frac(:,j,:,m-1) = leaf_frac(:,j,:,m-1) - d_leaf_frac(:,j,:,m)

          ! where it goes to
          leaf_frac(:,j,:,m) = leaf_frac(:,j,:,m) + d_leaf_frac(:,j,:,m)

       ENDDO       ! Loop over age classes

       !! 2.2.4 renormalize fractions in order to prevent accumulation 
       !       of numerical errors

       ! correct small negative values

       DO m = 1, nleafages
          leaf_frac(:,j,:,m) = MAX( zero, leaf_frac(:,j,:,m) )
       ENDDO

       ! total of fractions, should be very close to one where there is leaf mass

       sumfrac(:,:) = zero

       DO m = 1, nleafages       ! Loop over age classes

          sumfrac(:,:) = sumfrac(:,:) + leaf_frac(:,j,:,m)

       ENDDO       ! Loop over age classes

       ! normalize

       DO m = 1, nleafages       ! Loop over age classes

          WHERE ( sumfrac(:,:) .GT. min_stomate )

             leaf_frac(:,j,:,m) = leaf_frac(:,j,:,m) / sumfrac(:,:) 

          ELSEWHERE

             leaf_frac(:,j,:,m) = zero

          ENDWHERE

       ENDDO       ! Loop over age classes
       !JC add to avoid stop in debug mode
       ! due to the invalid float by fruit age
       WHERE (leaf_frac(:,j,:,:) .LE. min_stomate)
         leaf_age(:,j,:,:) = zero
       ENDWHERE
       ! fruit age-related turnover follows the ileaf
       WHERE (leaf_age(:,j,ifruit,:) .LE. zero .OR. &
              leaf_age(:,j,ifruit,:) .GT. 730) 
         leaf_age(:,j,ifruit,:) = zero
       ENDWHERE

    ENDDO         ! Loop over PFTs



    !==============================
    !DSG: age effect disabled START
    ! WARNING: the following age related decline in Vcmax 
    !          is reverted at the end of the routine.
    !          I keep this code for now, as it might be 
    !          useful later.
    !==============================
    ! TO MAKE IT WORK AGAIN ADJUST THE NUMBER OF DIMENSION FOR VCMAX, JMAX

    !
    !! 3 calculate vmax as a function of the age
    !

!    DO j = 2,nvm
!
!       vcmax(:,j) = zero
!       jmax(:,j)  = zero
!       nue(:,j)   = zero 
!
!       ! sum up over the different age classes
!       ! In case of OK_DGVM, evergreen needleleaf tress get a special treatment:
!       IF (ok_dgvm .AND. pheno_type(j)==1 .AND. leaf_tab(j)==2) THEN
!          ! pheno_typ=evergreen and leaf_tab=needleleaf
!          !vcmax(:,j) = Vcmax25(j)
!          vcmax(:,j) = vcmax_opt(:,j)
!          jmax(:,j)  = jmax_opt(:,j)
!
!       ELSE
!       ! all other cases:
!
!          DO m = 1, nleafages       ! Loop over age classes
!
!             !
!             !! 3.1 efficiency in each of the age classes
!             !!     it varies from vmax_offset to 1 
!             !!     linearly increases from vmax_offset to 1 for 0 < rel_age < leafage_firstmax
!             !!     is 1 when leafage_firstmax < rel_age < leafage_lastmax
!             !!     linearly decreases from 1 to vmax_offset for leafage_lastmax < rel_age < leafage_firstmax
!             !!     vmax_offset for rel_age >= leafage_old
!             !!     (Ishida et al., 1999)
!             rel_age(:) = leaf_age(:,j,ileaf,m) / leafagecrit(j)
!
!             leaf_efficiency(:) = MAX( vmax_offset, MIN( un, &
!                  vmax_offset + (un - vmax_offset) * rel_age(:) / leafage_firstmax, &
!                  un - (un - vmax_offset) * ( rel_age(:) - leafage_lastmax ) / &
!                  ( leafage_old - leafage_lastmax ) ) )
!
!             !
!             !! 3.2 add to mean vmax
!             !             
!             vcmax(:,j) = vcmax(:,j) + vcmax_opt(:,j) * leaf_efficiency(:) * leaf_frac(:,j,ileaf,m)
!             jmax(:,j)  = jmax(:,j)  + jmax_opt(:,j)  * leaf_efficiency(:) * leaf_frac(:,j,ileaf,m)
!             nue(:,j)   = nue(:,j)   + nue_opt(j)     * leaf_efficiency(:) * leaf_frac(:,j,ileaf,m)
!
!          ENDDO     ! loop over age classes
!       ENDIF
!
!    ENDDO       ! loop over PFTs

    !DSG: we bypass the leaf age related downregulation of Vcmax and Jmax
    vcmax(:,:,:) = vcmax_opt(:,:,:)
    jmax(:,:,:)  = jmax_opt(:,:,:)
    DO j=2,nvm  ! loop over PFT
         nue(:,j) = nue_opt(j)
    END DO
    


    !============================
    !DSG: age effect disabled END
    !============================

    IF (printlev>=4) WRITE(numout,*) 'Leaving vmax'

  END SUBROUTINE vmax

END MODULE stomate_vmax
