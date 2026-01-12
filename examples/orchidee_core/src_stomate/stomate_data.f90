! =================================================================================================================================
! MODULE 	: stomate_data
!
! CONTACT      : orchidee-help _at_ listes.ipsl.fr
!
! LICENCE      : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF         "stomate_data" module defines the values about the PFT parameters. It will print
!! the values of the parameters for STOMATE in the standard outputs. 
!!
!!\n DESCRIPTION: None 
!!
!! RECENT CHANGE(S): Sonke Zaehle: Reich et al, 1992 find no statistically significant differences 
!!                  between broadleaved and coniferous forests, specifically, the assumption that grasses grow 
!!                  needles is not justified. Replacing the function with the one based on Reich et al. 1997. 
!!                  Given that sla=100cm2/gDW at 9 months, sla is:
!!                  sla=exp(5.615-0.46*ln(leaflon in months))
!!
!! REFERENCE(S)	: None
!!
!! SVN          :
!! $HeadURL: svn://forge.ipsl.jussieu.fr/orchidee/branches/ORCHIDEE-CN-P/ORCHIDEE/src_stomate/stomate_data.f90 $
!! $Date: 2017-10-26 14:32:36 +0200 (四, 2017-10-26) $
!! $Revision: 4717 $
!! \n
!_ ================================================================================================================================

MODULE stomate_data

  ! modules used:

  USE constantes
  USE time, ONLY : one_day, dt_sechiba, one_year
  USE pft_parameters
  USE defprec
  

  IMPLICIT NONE

  INTEGER,ALLOCATABLE,SAVE,DIMENSION(:) :: hori_index     !! Move to Horizontal indices
!$OMP THREADPRIVATE(hori_index)

  INTEGER,ALLOCATABLE,SAVE,DIMENSION(:) :: horipft_index  !! Horizontal + PFT indices
!$OMP THREADPRIVATE(horipft_index)

  ! Land cover change

  INTEGER, ALLOCATABLE, SAVE, DIMENSION(:) :: horip10_index   !! Horizontal + P10 indices
!$OMP THREADPRIVATE(horip10_index)
  INTEGER, ALLOCATABLE, SAVE, DIMENSION(:) :: horip100_index  !! Horizontal + P100 indice
!$OMP THREADPRIVATE(horip100_index)
  INTEGER, ALLOCATABLE, SAVE, DIMENSION(:) :: horip11_index   !! Horizontal + P11 indices
!$OMP THREADPRIVATE(horip11_index)
  INTEGER, ALLOCATABLE, SAVE, DIMENSION(:) :: horip101_index  !! Horizontal + P101 indices
!$OMP THREADPRIVATE(horip101_index)

  INTEGER,SAVE :: itime                 !! time step
!$OMP THREADPRIVATE(itime)
  INTEGER,SAVE :: hist_id_stomate       !! STOMATE history file ID
!$OMP THREADPRIVATE(hist_id_stomate)
  INTEGER,SAVE :: hist_id_stomate_IPCC  !! STOMATE history file ID for IPCC output
!$OMP THREADPRIVATE(hist_id_stomate_IPCC)
  INTEGER,SAVE :: rest_id_stomate       !! STOMATE restart file ID
!$OMP THREADPRIVATE(rest_id_stomate)

  REAL,PARAMETER :: adapted_crit = 1. - ( 1. / euler ) !! critical value for being adapted (1-1/e) (unitless)
  REAL,PARAMETER :: regenerate_crit = 1. / euler       !! critical value for being regenerative (1/e) (unitless)


  ! private & public routines

  PUBLIC data

CONTAINS

!! ================================================================================================================================
!! SUBROUTINE 	: data
!!
!>\BRIEF         This routine defines the values of the PFT parameters. It will print the values of the parameters for STOMATE
!!               in the standard outputs of ORCHIDEE. 
!!
!! DESCRIPTION : This routine defines PFT parameters. It initializes the pheno_crit structure by tabulated parameters.\n
!!               Some initializations are done for parameters. The SLA is calculated according *to* Reich et al (1992).\n
!!               Another formulation by Reich et al(1997) could be used for the computation of the SLA.
!!               The geographical coordinates might be used for defining some additional parameters
!!               (e.g. frequency of anthropogenic fires, irrigation of agricultural surfaces, etc.). \n
!!               For the moment, this possibility is not used. \n
!!               The specifc leaf area (SLA) is calculated according Reich et al, 1992 by :
!!               \latexonly
!!               \input{stomate_data_SLA.tex}
!!               \endlatexonly
!!               The sapling (young) biomass for trees and for each compartment of biomass is calculated by :
!!               \latexonly
!!               \input{stomate_data_sapl_tree.tex}
!!               \endlatexonly
!!               The sapling biomass for grasses and for each compartment of biomass is calculated by :
!!               \latexonly
!!               \input{stomate_data_sapl_grass.tex}
!!               \endlatexonly
!!               The critical stem diameter is given by the following formula :
!!               \latexonly
!!               \input{stomate_data_stem_diameter.tex}
!!               \endlatexonly
!!
!! RECENT CHANGE(S): Sonke Zaehle: Reich et al, 1992 find no statistically significant differences 
!!                  between broadleaved and coniferous forests, specifically, the assumption that grasses grow 
!!                  needles is not justified. Replacing the function with the one based on Reich et al. 1997. 
!!                  Given that sla=100cm2/gDW at 9 months, sla is:
!!                  sla=exp(5.615-0.46*ln(leaflon in months)) 
!!                   \latexonly
!!                   \input{stomate_data_SLA_Reich_97.tex}
!!                   \endlatexonly
!!
!! MAIN OUTPUT VARIABLE(S): 
!!
!! REFERENCE(S) :
!! - Reich PB, Walters MB, Ellsworth DS, (1992), Leaf life-span in relation to leaf, plant and 
!! stand characteristics among diverse ecosystems. Ecological Monographs, Vol 62, pp 365-392.
!! - Reich PB, Walters MB, Ellsworth DS (1997) From tropics to tundra: global convergence in plant 
!!  functioning. Proc Natl Acad Sci USA, 94:13730 13734
!!
!! FLOWCHART    :
!! \n
!_ ================================================================================================================================

  SUBROUTINE data (npts, lalo)


    !! 0. Variables and parameter declaration


    !! 0.1 Input variables

    INTEGER, INTENT(in)                   :: npts    !! [DISPENSABLE] Domain size (unitless)
    REAL,DIMENSION (npts,2), INTENT (in)  :: lalo    !! [DISPENSABLE] Geographical coordinates (latitude,longitude)

    !! 0.4 Local variables

    INTEGER                               :: i,j     !! Index (unitless)
    REAL                                  :: alpha   !! alpha's : (unitless)
    REAL                                  :: dia     !! stem diameter (m)
    REAL                                  :: csa_sap !! Crown specific area sapling @tex $(m^2.ind^{-1})$ @endtex
    REAL                                  :: cn_leaf         !! C to N ratio of Leaf pool (gC per gN) 
    REAL                                  :: cn_wood         !! C to N ratio of Woody pools (gC per gN) 
    REAL                                  :: cn_root         !! C to N ratio of Root pool (gC per gN) 
    REAL                                  :: np_leaf         !! N to P ratio of Leaf pool (gN per gP)
    REAL                                  :: np_wood         !! N to P ratio of Woody pools (gN per gP)
    REAL                                  :: np_root         !! N to P ratio of Root pool (gN per gP)
!_ ================================================================================================================================

    !- pheno_gdd_crit
    pheno_gdd_crit(:,1) = pheno_gdd_crit_c(:)
    pheno_gdd_crit(:,2) = pheno_gdd_crit_b(:)         
    pheno_gdd_crit(:,3) = pheno_gdd_crit_a(:) 
    !
    !- senescence_temp
    senescence_temp(:,1) = senescence_temp_c(:)
    senescence_temp(:,2) = senescence_temp_b(:)
    senescence_temp(:,3) = senescence_temp_a(:)

    ! 
    !-LC 
    LC(:,ileaf) = LC_leaf(:) 
    LC(:,isapabove) = LC_sapabove(:) 
    LC(:,isapbelow) = LC_sapbelow(:) 
    LC(:,iheartabove) = LC_heartabove(:) 
    LC(:,iheartbelow) = LC_heartbelow(:) 
    LC(:,iroot) = LC_root(:) 
    LC(:,ifruit) = LC_fruit(:) 
    LC(:,icarbres) = LC_carbres(:) 
    LC(:,ilabile) = LC_labile(:) 

    IF ( printlev >= 2 ) WRITE(numout,*) 'data: PFT characteristics'

    DO j = 2,nvm ! Loop over # PFTS 

       IF ( printlev >= 2 ) WRITE(numout,'(a,i3,a,a)') '    > PFT#',j,': ', PFT_name(j)

       !
       ! 1 tree? (true/false)
       !
       IF ( printlev >= 2 ) WRITE(numout,*) '       tree: (::is_tree) ', is_tree(j)

       !
       ! 2 flamability (0-1, unitless)
       !

       IF ( printlev >= 2 ) WRITE(numout,*) '       litter flamability (::flam) :', flam(j)

       !
       ! 3 fire resistance (unitless)
       !

       IF ( printlev >= 2 ) WRITE(numout,*) '       fire resistance (::resist) :', resist(j)

       !
       ! 4 specific leaf area per mass carbon = 2 * sla / dry mass (m^2.gC^{-1})
       !

       ! S. Zaehle: Reich et al, 1992 find no statistically significant differences between broadleaved and coniferous
       ! forests, specifically, the assumption that grasses grow needles is not justified. Replacing the function
       ! with the one based on Reich et al. 1997. Given that sla=100cm2/gDW at 9 months, sla is:
       ! sla=exp(5.615-0.46*ln(leaflon in months))

       ! Oct 2010 : sla values are prescribed by values given by N.Viovy 

       ! includes conversion from 
       !!       sla(j) = 2. * 1e-4 * EXP(5.615 - 0.46 * log(12./leaflife_tab(j)))
       !!\latexonly
       !!\input{stomate_data_SLA.tex}
       !!\endlatexonly
!       IF ( leaf_tab(j) .EQ. 2 ) THEN
!
!          ! needle leaved tree
!          sla(j) = 2. * ( 10. ** ( 2.29 - 0.4 * LOG10(12./leaflife_tab(j)) ) ) *1e-4
!
!       ELSE
!
!          ! broad leaved tree or grass (Reich et al 1992)
!          sla(j) = 2. * ( 10. ** ( 2.41 - 0.38 * LOG10(12./leaflife_tab(j)) ) ) *1e-4
!
!       ENDIF

!!!$      IF ( leaf_tab(j) .EQ. 1 ) THEN
!!!$
!!!$        ! broad leaved tree
!!!$
!!!$        sla(j) = 2. * ( 10. ** ( 2.41 - 0.38 * LOG10(12./leaflife_tab(j)) ) ) *1e-4
!!!$
!!!$      ELSE
!!!$
!!!$        ! needle leaved or grass (Reich et al 1992)
!!!$
!!!$        sla(j) = 2. * ( 10. ** ( 2.29 - 0.4 * LOG10(12./leaflife_tab(j)) ) ) *1e-4
!!!$
!!!$      ENDIF
!!!$
!!!$      IF ( ( leaf_tab(j) .EQ. 2 ) .AND. ( pheno_type_tab(j) .EQ. 2 ) ) THEN
!!!$
!!!$        ! summergreen needle leaf
!!!$
!!!$        sla(j) = 1.25 * sla(j)
!!!$
!!!$      ENDIF

       IF ( printlev >= 2 ) WRITE(numout,*) '       specific leaf area (m**2/gC) (::sla):', sla(j), 12./leaflife_tab(j)

       !
       ! 5 sapling characteristics
       !

       IF ( is_tree(j) ) THEN

          !> 5.1 trees

          !!\latexonly
          !!\input{stomate_data_sapl_tree.tex}
          !!\endlatexonly

          alpha = alpha_tree

          bm_sapl(j,ileaf,icarbon) = &
               &     ((bm_sapl_leaf(1)*pipe_tune1(j)*(mass_ratio_heart_sap(j) *bm_sapl_leaf(2)*sla(j)/(pi*pipe_k1(j))) & 
               &     **bm_sapl_leaf(3))/sla(j))**bm_sapl_leaf(4)

          IF ( pheno_type(j) .NE. 1 ) THEN
             ! not evergreen
             bm_sapl(j,icarbres,icarbon) = bm_sapl_carbres * bm_sapl(j,ileaf,icarbon)
             bm_sapl(j,ilabile,icarbon) = bm_sapl_labile * bm_sapl(j,ileaf,icarbon)
          ELSE
             bm_sapl(j,icarbres,icarbon) = zero
             bm_sapl(j,ilabile,icarbon) = zero
          ENDIF ! (pheno_type_tab(j) .NE. 1 )

          csa_sap = bm_sapl(j,ileaf,icarbon) / ( pipe_k1(j) / sla(j) )

          dia = (mass_ratio_heart_sap(j) * csa_sap * dia_coeff(1) / pi ) ** dia_coeff(2)

          bm_sapl(j,isapabove,icarbon) = &
               bm_sapl_sapabove * pipe_density(j) * csa_sap * pipe_tune2(j) * dia ** pipe_tune3(j)
          bm_sapl(j,isapbelow,icarbon) = bm_sapl(j,isapabove,icarbon)

          bm_sapl(j,iheartabove,icarbon) = bm_sapl_heartabove * bm_sapl(j,isapabove,icarbon)
          bm_sapl(j,iheartbelow,icarbon) = bm_sapl_heartbelow * bm_sapl(j,isapbelow,icarbon)

       ELSE

          !> 5.2 grasses

          !!\latexonly
          !!\input{stomate_data_sapl_grass.tex}
          !!\endlatexonly

          alpha = alpha_grass

          IF ( natural(j) ) THEN
             bm_sapl(j,ileaf,icarbon) = init_sapl_mass_leaf_nat / sla(j)
          ELSE
             bm_sapl(j,ileaf,icarbon) = init_sapl_mass_leaf_agri / sla(j)
          ENDIF

          bm_sapl(j,icarbres,icarbon) = init_sapl_mass_carbres *bm_sapl(j,ileaf,icarbon)
          bm_sapl(j,ilabile,icarbon) = init_sapl_mass_labile *bm_sapl(j,ileaf,icarbon)

          bm_sapl(j,isapabove,icarbon) = zero
          bm_sapl(j,isapbelow,icarbon) = zero

          bm_sapl(j,iheartabove,icarbon) = zero
          bm_sapl(j,iheartbelow,icarbon) = zero

       ENDIF !( is_tree(j) )

       bm_sapl(j,iroot,icarbon) = init_sapl_mass_root * (1./alpha) * bm_sapl(j,ileaf,icarbon)

       bm_sapl(j,ifruit,icarbon) = init_sapl_mass_fruit  * bm_sapl(j,ileaf,icarbon)


       cn_leaf=cn_leaf_init(j)
       cn_wood=cn_leaf/fcn_wood(j)
       cn_root=cn_leaf/fcn_root(j)

       np_leaf=np_leaf_init(j)
       np_wood=np_leaf/fnp_wood(j)
       np_root=np_leaf/fnp_root(j)
       
       DO i=1,nparts
          IF (i.EQ.ileaf) THEN
             bm_sapl(j,i,initrogen) = bm_sapl(j,i,icarbon) / cn_leaf
             bm_sapl(j,i,iphosphorus) = bm_sapl(j,i,initrogen) / np_leaf
          ELSE IF (i.LT.iroot) THEN
             bm_sapl(j,i,initrogen) = bm_sapl(j,i,icarbon) / cn_wood
             bm_sapl(j,i,iphosphorus) = bm_sapl(j,i,initrogen) / np_wood
          ELSE
             bm_sapl(j,i,initrogen) = bm_sapl(j,i,icarbon) / cn_root
             bm_sapl(j,i,iphosphorus) = bm_sapl(j,i,initrogen) / np_root
          ENDIF
       ENDDO

       IF ( printlev >= 2 ) THEN
          WRITE(numout,*) '       sapling biomass (gC):'
          WRITE(numout,*) '         leaves: (::bm_sapl(j,ileaf,icarbon))',bm_sapl(j,ileaf,icarbon)
          WRITE(numout,*) '         sap above ground: (::bm_sapl(j,ispabove,icarbon)):',bm_sapl(j,isapabove,icarbon)
          WRITE(numout,*) '         sap below ground: (::bm_sapl(j,isapbelow,icarbon))',bm_sapl(j,isapbelow,icarbon)
          WRITE(numout,*) '         heartwood above ground: (::bm_sapl(j,iheartabove,icarbon))',bm_sapl(j,iheartabove,icarbon)
          WRITE(numout,*) '         heartwood below ground: (::bm_sapl(j,iheartbelow,icarbon))',bm_sapl(j,iheartbelow,icarbon)
          WRITE(numout,*) '         roots: (::bm_sapl(j,iroot,icarbon))',bm_sapl(j,iroot,icarbon)
          WRITE(numout,*) '         fruits: (::bm_sapl(j,ifruit,icarbon))',bm_sapl(j,ifruit,icarbon)
          WRITE(numout,*) '         carbohydrate reserve: (::bm_sapl(j,icarbres,icarbon))',bm_sapl(j,icarbres,icarbon)
          WRITE(numout,*) '         labile reserve: (::bm_sapl(j,ilabile,icarbon))',bm_sapl(j,ilabile,icarbon)
       ENDIF

       !
       ! 6 migration speed (m/year)
       !

       IF ( is_tree(j) ) THEN

          migrate(j) = migrate_tree

       ELSE

          ! can be any value as grasses are, per *definition*, everywhere (big leaf).
          migrate(j) = migrate_grass

       ENDIF !( is_tree(j) )

       IF ( printlev >= 2 ) WRITE(numout,*) '       migration speed (m/year): (::migrate(j))', migrate(j)

       !
       ! 7 critical stem diameter: beyond this diameter, the crown area no longer
       !     increases (m)
       !

       IF ( is_tree(j) ) THEN

          !!\latexonly
          !!\input{stomate_data_stem_diameter.tex}
          !!\endlatexonly

          maxdia(j) = ( ( pipe_tune4(j) / ((pipe_tune2(j)*pipe_tune3(j))/(maxdia_coeff(1)**pipe_tune3(j))) ) &
               ** ( un / ( pipe_tune3(j) - un ) ) ) * maxdia_coeff(2)
          cn_sapl(j) = cn_sapl_init !crown of individual tree, first year

       ELSE

          maxdia(j) = undef
          cn_sapl(j)=1

       ENDIF !( is_tree(j) )

       IF ( printlev >= 2 ) WRITE(numout,*) '       critical stem diameter (m): (::maxdia(j))', maxdia(j)

       !
       ! 8 Coldest tolerable temperature (K)
       !

       IF ( ABS( tmin_crit(j) - undef ) .GT. min_stomate ) THEN
          tmin_crit(j) = tmin_crit(j) + ZeroCelsius
       ELSE
          tmin_crit(j) = undef
       ENDIF 

       IF ( printlev >= 2 ) &
            WRITE(numout,*) '       coldest tolerable temperature (K): (::tmin_crit(j))', tmin_crit(j)

       !
       ! 9 Maximum temperature of the coldest month: need to be below this temperature
       !      for a certain time to regrow leaves next spring *(vernalization)* (K)
       !

       IF ( ABS ( tcm_crit(j) - undef ) .GT. min_stomate ) THEN
          tcm_crit(j) = tcm_crit(j) + ZeroCelsius
       ELSE
          tcm_crit(j) = undef
       ENDIF

       IF ( printlev >= 2 ) &
            WRITE(numout,*) '       vernalization temperature (K): (::tcm_crit(j))', tcm_crit(j)

       !
       ! 10 critical values for phenology
       !

       ! 10.1 model used

       IF ( printlev >= 2 ) &
            WRITE(numout,*) '       phenology model used: (::pheno_model(j)) ',pheno_model(j)

       ! 10.2 growing degree days. What kind of gdd is meant (i.e. threshold 0 or -5 deg C
       !        or whatever), depends on how this is used in stomate_phenology.


       IF ( ( printlev >= 2 ) .AND. ( ALL(pheno_gdd_crit(j,:) .NE. undef) ) ) THEN
          WRITE(numout,*) '         critical GDD is a function of long term T (C): (::gdd)'
          WRITE(numout,*) '          ',pheno_gdd_crit(j,1), &
               ' + T *',pheno_gdd_crit(j,2), &
               ' + T^2 *',pheno_gdd_crit(j,3)
       ENDIF

       ! consistency check

       IF ( ( ( pheno_model(j) .EQ. 'moigdd' ) .OR. &
            ( pheno_model(j) .EQ. 'humgdd' )       ) .AND. &
            ( ANY(pheno_gdd_crit(j,:) .EQ. undef) )                      ) THEN
          CALL ipslerr_p(3,'stomate_data','problem with phenology parameters, critical GDD. (::pheno_model)','','')
       ENDIF

       ! 10.3 number of growing days

       IF ( ( printlev >= 2 ) .AND. ( ngd_crit(j) .NE. undef ) ) &
            WRITE(numout,*) '         critical NGD: (::ngd_crit(j))', ngd_crit(j)

       ! 10.4 critical temperature for ncd vs. gdd function in phenology (C)

       IF ( ( printlev >= 2 ) .AND. ( ncdgdd_temp(j) .NE. undef ) ) &
            WRITE(numout,*) '         critical temperature for NCD vs. GDD (C): (::ncdgdd_temp(j))', &
            ncdgdd_temp(j)

       ! 10.5 humidity fractions (0-1, unitless)

       IF ( ( printlev >= 2 ) .AND. ( hum_frac(j) .NE. undef ) ) &
            WRITE(numout,*) '         critical humidity fraction: (::hum_frac(j))', &
            &  hum_frac(j)


       ! 10.6 minimum time elapsed since moisture minimum (days)

       IF ( ( printlev >= 2 ) .AND. ( hum_min_time(j) .NE. undef ) ) &
            WRITE(numout,*) '         time to wait after moisture min (d): (::hum_min_time(j))', &
        &    hum_min_time(j)

       !
       ! 11 critical values for senescence
       !

       ! 11.1 type of senescence

       IF ( printlev >= 2 ) &
            WRITE(numout,*) '       type of senescence: (::senescence_type(j))',senescence_type(j)

       ! 11.2 critical temperature for senescence (C)

       IF ( ( printlev >= 2 ) .AND. ( ALL(senescence_temp(j,:) .NE. undef) ) ) THEN
          WRITE(numout,*) '         critical temperature for senescence (C) is'
          WRITE(numout,*) '          a function of long term T (C): (::senescence_temp)'
          WRITE(numout,*) '          ',senescence_temp(j,1), &
               ' + T *',senescence_temp(j,2), &
               ' + T^2 *',senescence_temp(j,3)
       ENDIF

       ! consistency check

       IF ( ( ( senescence_type(j) .EQ. 'cold' ) .OR. &
            ( senescence_type(j) .EQ. 'mixed' )      ) .AND. &
            ( ANY(senescence_temp(j,:) .EQ. undef ) )           ) THEN
          CALL ipslerr_p(3,'stomate_data','Problem with senescence parameters, temperature. (::senescence_type)','','')
       ENDIF

       ! 11.3 critical relative moisture availability for senescence

       IF ( ( printlev >= 2 ) .AND. ( senescence_hum(j) .NE. undef ) ) THEN
          WRITE(numout,*)  '  max. critical relative moisture availability for' 
          WRITE(numout,*)  '  senescence: (::senescence_hum(j))',  &
               & senescence_hum(j)
       END IF

       ! consistency check

       IF ( ( ( senescence_type(j) .EQ. 'dry' ) .OR. &
            ( senescence_type(j) .EQ. 'mixed' )     ) .AND. &
            ( senescence_hum(j) .EQ. undef )                   ) THEN
          CALL ipslerr_p(3,'stomate_data','Problem with senescence parameters, humidity.(::senescence_type)','','')
       ENDIF

       ! 14.3 relative moisture availability above which there is no moisture-related
       !      senescence (0-1, unitless)

       IF ( ( printlev >= 2 ) .AND. ( nosenescence_hum(j) .NE. undef ) ) THEN
          WRITE(numout,*) '         relative moisture availability above which there is' 
          WRITE(numout,*) '             no moisture-related senescence: (::nosenescence_hum(j))', nosenescence_hum(j)
       END IF
       ! consistency check

       IF ( ( ( senescence_type(j) .EQ. 'dry' ) .OR. &
            ( senescence_type(j) .EQ. 'mixed' )     ) .AND. &
            ( nosenescence_hum(j) .EQ. undef )                   ) THEN
          CALL ipslerr_p(3,'stomate_data','Problem with senescence parameters, humidity. (::senescence_type)','','')
       ENDIF

       !
       ! 12 sapwood -> heartwood conversion time (days)
       !

       IF ( printlev >= 2 ) &
            WRITE(numout,*) '       sapwood -> heartwood conversion time (d): (::tau_sap(j))', tau_sap(j)

       !
       ! 13 fruit lifetime (days)
       !

       IF ( printlev >= 2 ) WRITE(numout,*) '       fruit lifetime (d): (::tau_fruit(j))', tau_fruit(j)

       !
       ! 14 length of leaf death (days)
       !      For evergreen trees, this variable determines the lifetime of the leaves.
       !      Note that it is different from the value given in leaflife_tab.
       !

       IF ( printlev >= 2 ) &
            WRITE(numout,*) '       length of leaf death (d): (::leaffall(j))', leaffall(j)

       !
       ! 15 maximum lifetime of leaves (days)
       !

       IF ( ( printlev >= 2 ) .AND. ( leafagecrit(j) .NE. undef ) ) &
            WRITE(numout,*) '       critical leaf age (d): (::leafagecrit(j))', leafagecrit(j)

       IF ( ( printlev >= 1 ) .AND. ( rootagecrit(j) .NE. undef ) ) &
            WRITE(numout,*) '       critical leaf age (d): (::rootagecrit(j))', rootagecrit(j)

       !
       ! 16 time constant for leaf age discretisation (days)
       !

       leaf_timecst(j) = leafagecrit(j) / REAL( nleafages,r_std )

       IF ( printlev >= 2 ) &
            WRITE(numout,*) '       time constant for leaf age discretisation (d): (::leaf_timecst(j))', &
            leaf_timecst(j)

       !
       ! 17 minimum lai, initial (m^2.m^{-2})
       !

       IF ( is_tree(j) ) THEN
          lai_initmin(j) = lai_initmin_tree
       ELSE
          lai_initmin(j) = lai_initmin_grass
       ENDIF !( is_tree(j) )

       IF ( printlev >= 2 ) &
            WRITE(numout,*) '       initial LAI: (::lai_initmin(j))', lai_initmin(j)

       !
       ! 19 maximum LAI (m^2.m^{-2})
       !

       IF ( printlev >= 2 ) &
            WRITE(numout,*) '       critical LAI above which no leaf allocation: (::lai_max(j))', lai_max(j)

       !
       ! 20 fraction of primary leaf and root allocation put into reserve (0-1, unitless)
       !

       IF ( printlev >= 2 ) &
            WRITE(numout,*) '       reserve allocation factor: (::ecureuil(j))', ecureuil(j)


       !
       ! 23 natural ?
       !

       IF ( printlev >= 2 ) &
            WRITE(numout,*) '       Natural: (::natural(j))', natural(j)

       !
       ! 24 Vcmax et Vjmax (umol.m^{-2}.s^{-1}) 
       !

       IF ( printlev >= 2 ) &
            WRITE(numout,*) '       Maximum rate of carboxylation: (::Vcmax_25(j))', vcmax25(j)
       !
       ! 25 constants for photosynthesis temperatures
       !

       IF ( printlev >= 2 ) THEN


          !
          ! 26 Properties
          !

          WRITE(numout,*) '       C4 photosynthesis: (::is_c4(j))', is_c4(j)
          WRITE(numout,*) '       Depth constant for root profile (m): (::1./humcste(j))', 1./humcste(j)

       ENDIF

       !
       ! 27 extinction coefficient of the Monsi and Saeki (1953) relationship 
       !
       IF ( printlev >= 2 ) THEN
          WRITE(numout,*) '       extinction coefficient: (::ext_coeff(j))', ext_coeff(j)
       ENDIF

       !
       ! 30 fraction of allocatable biomass which is lost as growth respiration (0-1, unitless)
       !
       IF ( printlev >= 2 ) &
            WRITE(numout,*) '       growth respiration fraction: (::frac_growthresp(j))', frac_growthresp(j)

    ENDDO ! Loop over # PFTS 

    !
    ! 29 time scales for phenology and other processes (in days)
    !

    tau_longterm_max = coeff_tau_longterm * one_year

    IF ( printlev >= 2 ) THEN

       WRITE(numout,*) '   > time scale for ''monthly'' moisture availability (d): (::tau_hum_month)', &
            tau_hum_month
       WRITE(numout,*) '   > time scale for ''weekly'' moisture availability (d): (::tau_hum_week)', &
           tau_hum_week
       WRITE(numout,*) '   > time scale for ''monthly'' 2 meter temperature (d): (::tau_t2m_month)', &
            tau_t2m_month
       WRITE(numout,*) '   > time scale for ''weekly'' 2 meter temperature (d): (::tau_t2m_week)', &
            tau_t2m_week
       WRITE(numout,*) '   > time scale for ''weekly'' GPP (d): (::tau_gpp_week)', &
            tau_gpp_week
       WRITE(numout,*) '   > time scale for ''monthly'' soil temperature (d): (::tau_tsoil_month)', &
            tau_tsoil_month
       WRITE(numout,*) '   > time scale for ''monthly'' soil humidity (d): (::tau_soilhum_month)', &
            tau_soilhum_month
       WRITE(numout,*) '   > time scale for vigour calculations (y): (::tau_longterm_max / one_year)', &
            tau_longterm_max / one_year

    ENDIF

    IF (printlev >= 4) WRITE(numout,*) 'Leaving stomate_data'

  END SUBROUTINE data

END MODULE stomate_data
