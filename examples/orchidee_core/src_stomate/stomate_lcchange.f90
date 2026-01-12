! =================================================================================================================================
! MODULE       : stomate_lcchange
!
! CONTACT      : orchidee-help _at_ listes.ipsl.fr
!
! LICENCE      : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF       Impact of land cover change on carbon stocks
!!
!!\n DESCRIPTION: None
!!
!! RECENT CHANGE(S): None
!!
!! REFERENCE(S)	: None
!!
!! SVN          :
!! $HeadURL: svn://forge.ipsl.jussieu.fr/orchidee/branches/ORCHIDEE-CN-P/ORCHIDEE/src_stomate/stomate_lcchange.f90 $
!! $Date: 2021-04-29 14:43:02 +0200 (四, 2021-04-29) $
!! $Revision: 7172 $
!! \n
!_ ================================================================================================================================


MODULE stomate_lcchange

  ! modules used:
  
  USE ioipsl_para
  USE stomate_data
  USE pft_parameters
  USE constantes
  
  IMPLICIT NONE
  
  PRIVATE
  PUBLIC lcchange_main
  
CONTAINS


!! ================================================================================================================================
!! SUBROUTINE   : lcchange_main
!!
!>\BRIEF        Impact of land cover change on carbon stocks
!!
!! DESCRIPTION  : This subroutine is always activate if VEGET_UPDATE>0Y in the configuration file, which means that the 
!! vegetation map is updated regulary. lcchange_main is called from stomateLpj the first time step after the vegetation 
!! map has been changed. 
!! The impact of land cover change on carbon stocks is computed in this subroutine. The land cover change is written
!! by the difference of current and previous "maximal" coverage fraction of a PFT. 
!! On the basis of this difference, the amount of 'new establishment'/'biomass export',
!! and increase/decrease of each component, are estimated.\n
!!
!! Main structure of lpj_establish.f90 is:
!! 1. Initialization
!! 2. Calculation of changes in carbon stocks and biomass by land cover change
!! 3. Update 10 year- and 100 year-turnover pool contents
!! 4. History
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S) : ::prod10, ::prod100, ::flux10, ::flux100,
!!   :: cflux_prod10 and :: cflux_prod100 
!!
!! REFERENCES   : None
!!
!! FLOWCHART    : 
!! \latexonly 
!!     \includegraphics[scale=0.5]{lcchange.png}
!! \endlatexonly
!! \n
!_ ================================================================================================================================

  
  SUBROUTINE lcchange_main ( npts, dt_days, veget_cov_max_old, veget_cov_max_new, &
       biomass, ind, age, PFTpresent, senescence, when_growthinit, everywhere, &        
       co2_to_bm, n_to_bm, p_to_bm, &
       bm_to_litter, turnover_daily, bm_sapl, cn_ind,flux10,flux100, &
       prod10,prod100, convflux, cflux_prod10,cflux_prod100,&
       nflux_prod_total, pflux_prod_total, leaf_frac,&
       flux10_harvest,flux100_harvest, prod10_harvest,prod100_harvest,&
       convflux_harvest,cflux_prod10_harvest,cflux_prod100_harvest,&
       npp_longterm, lm_lastyearmax, litter, &
       litter_avail, litter_not_avail, &
       som, soil_n_min, soil_p_min, KF, k_latosa_adapt, rue_longterm, &
       lignin_struc, lignin_wood, &
       harvestwood,harvestpft)  ! do_wood_harvest 

    IMPLICIT NONE
    
  !! 0. Variable and parameter declaration 
    
    !! 0.1 Input variables
    
    INTEGER, INTENT(in)                                       :: npts             !! Domain size - number of pixels (unitless)
    REAL, INTENT(in)                                   :: dt_days          !! Time step of vegetation dynamics for stomate
                                                                                  !! (days)
    REAL, DIMENSION(nvm, nparts,nelements), INTENT(in) :: bm_sapl          !! biomass of sapling 
                                                                                  !! @tex ($gC individual^{-1}$) @endtex
    REAL, DIMENSION(npts,nvm), INTENT(in)              :: veget_cov_max_old!! Current "maximal" coverage fraction of a PFT (LAI
                                                                                  !! -> infinity) on ground
    REAL, DIMENSION(npts,nvm), INTENT(in)              :: veget_cov_max_new!! New "maximal" coverage fraction of a PFT (LAI ->
                                                                                  !! infinity) on ground (unitless) 
    REAL, DIMENSION(npts), INTENT(in)                  :: harvestwood      !! harvested wood (gC/m2/year)
 
    !! 0.2 Output variables

    REAL, DIMENSION(npts), INTENT(out)                 :: convflux         !! release during first year following land cover
                                                                                  !! change
    REAL, DIMENSION(npts), INTENT(out)                 :: cflux_prod10     !! total annual release from the 10 year-turnover
                                                                                  !! pool @tex ($gC m^{-2}$) @endtex
    REAL, DIMENSION(npts), INTENT(out)                 :: cflux_prod100    !! total annual release from the 100 year-
                                                                                  !! turnover pool @tex ($gC m^{-2}$) @endtex
    REAL, DIMENSION(npts), INTENT(out)                 :: nflux_prod_total !! release of N associated to land cover change  @tex ($gN m^{-2}$) @endtex 
    REAL, DIMENSION(npts), INTENT(out)                 :: pflux_prod_total !! release of P associated to land cover change  @tex ($gPm^{-2}$) @endtex
    REAL, DIMENSION(npts,nvm,nparts,nelements), INTENT(inout):: turnover_daily   !! Turnover rates 
 
    REAL, DIMENSION(npts), INTENT(out)                 :: convflux_harvest      !! flux release during first year following wood harvest           
    REAL, DIMENSION(npts), INTENT(out)                 :: cflux_prod10_harvest  !! total annual release from the 10 year-turnover
                                                                                       !! pool @tex ($gC m^{-2}$) @endtex
    REAL, DIMENSION(npts), INTENT(out)                 :: cflux_prod100_harvest !! total annual release from the 100 year-
                                                                                       !! turnover pool @tex ($gC m^{-2}$) @endtex
    REAL, DIMENSION(npts,nvm), INTENT(out)             :: harvestpft            !! wood harvested for each PFT
                                                                                       !! @tex ($gC m^{-2} day^{-1}$) @endtex

    !! 0.3 Modified variables   
    
    REAL, DIMENSION(npts,nvm,nparts,nelements), INTENT(inout):: biomass    !! biomass @tex ($gC m^{-2}$) @endtex
    REAL, DIMENSION(npts,nvm), INTENT(inout)           :: ind              !! Number of individuals @tex ($m^{-2}$) @endtex
    REAL, DIMENSION(npts,nvm), INTENT(inout)           :: age              !! mean age (years)
    LOGICAL, DIMENSION(npts,nvm), INTENT(inout)               :: senescence       !! plant senescent (only for deciduous trees) Set
                                                                                  !! to .FALSE. if PFT is introduced or killed
    LOGICAL, DIMENSION(npts,nvm), INTENT(inout)               :: PFTpresent       !! Is pft there (unitless)
    REAL, DIMENSION(npts,nvm), INTENT(inout)           :: everywhere       !! is the PFT everywhere in the grid box or very 
                                                                                  !! localized (unitless)
    REAL, DIMENSION(npts,nvm), INTENT(inout)           :: when_growthinit  !! how many days ago was the beginning of the 
                                                                                  !! growing season (days)
    REAL, DIMENSION(npts,nvm), INTENT(inout)           :: co2_to_bm        !! biomass uptaken 
                                                                                  !! @tex ($gC m^{-2} day^{-1}$) @endtex
    REAL, DIMENSION(npts,nvm), INTENT(inout)           :: n_to_bm          !! N taken to create the seedlings 
    REAL, DIMENSION(npts,nvm), INTENT(inout)           :: p_to_bm          !! P taken to create the seedlings
    REAL, DIMENSION(npts,nvm,nparts,nelements), INTENT(inout) :: bm_to_litter !! conversion of biomass to litter 
                                                                                  !! @tex ($gC m^{-2} day^{-1}$) @endtex
    REAL, DIMENSION(npts,nvm), INTENT(inout)           :: cn_ind           !! crown area of individuals 
                                                                                  !! @tex ($m^{2}$) @endtex
    REAL, DIMENSION(npts,0:10), INTENT(inout)          :: prod10           !! products remaining in the 10 year-turnover
                                                                                  !! pool after the annual release for each 
                                                                                  !! compartment (10 + 1 : input from year of land
                                                                                  !! cover change)
    REAL, DIMENSION(npts,0:100), INTENT(inout)         :: prod100          !! products remaining in the 100 year-turnover
                                                                                  !! pool after the annual release for each 
                                                                                  !! compartment (100 + 1 : input from year of land
                                                                                  !! cover change)
    REAL, DIMENSION(npts,10), INTENT(inout)            :: flux10           !! annual release from the 10/100 year-turnover 
                                                                                  !! pool compartments
    REAL, DIMENSION(npts,100), INTENT(inout)           :: flux100          !! annual release from the 10/100 year-turnover
                                                                                  !! pool compartments
    REAL, DIMENSION(npts,nvm,nparts,nleafages), INTENT(inout) :: leaf_frac        !! fraction of leaves in leaf age class 
                                                                                  !! (unitless)
    REAL, DIMENSION(npts,nvm), INTENT(inout)           :: lm_lastyearmax   !! last year's maximum leaf mass for each PFT 
                                                                                  !! @tex ($gC m^{-2}$) @endtex
    REAL, DIMENSION(npts,nvm), INTENT(inout)           :: npp_longterm     !! "long term" net primary productivity 
                                                                                  !! @tex ($gC m^{-2} year^{-1}$) @endtex
    REAL,DIMENSION(npts,nlitt,nvm,nlevs,nelements), INTENT(inout):: litter !! metabolic and structural litter, above and 
                                                                                  !! below ground @tex ($gC m^{-2}$) @endtex
    REAL,DIMENSION(npts,ncarb,nvm,nelements), INTENT(inout)      :: som    !! SOM pool: active, slow, or passive  
  !! @tex ($g(C or N) m^{-2}$) @endtex 
    REAL,DIMENSION(npts,nvm,nnspec), INTENT(inout)     :: soil_n_min       !!  soil mineral nitrogen pool 
    REAL,DIMENSION(npts,nvm,npspec), INTENT(inout)     :: soil_p_min       !!  soil mineral phosphorus pool
    REAL, DIMENSION(:,:), INTENT(inout)                :: KF               !! Scaling factor to convert sapwood mass into leaf 
                                                                                  !! mass (m)
    REAL, DIMENSION(:,:), INTENT(inout)                :: k_latosa_adapt   !! Leaf to sapwood area adapted for long 
                                                                                  !! term water stress (m)
    REAL, DIMENSION(:,:), INTENT(inout)                :: rue_longterm     !! Longterm radiation use efficiency 
                                                                                  !! (??units??) 
    REAL, DIMENSION(npts,nvm,nlevs), INTENT(inout)     :: lignin_struc     !! ratio Lignine/Carbon in structural litter,
                                                                                  !! above and below ground
    REAL, DIMENSION(npts,nvm,nlevs), INTENT(inout)     :: lignin_wood      !! ratio Lignine/Carbon in woody litter,
                                                                                  !! above and below ground
    REAL, DIMENSION(npts,nlitt,nvm), INTENT(inout)     :: litter_avail     !! litter available for grazing     
    REAL, DIMENSION(npts,nlitt,nvm) , INTENT(inout)    :: litter_not_avail !! litter not available for grazing (e.g., manure)

    ! products remaining in the 10/100 year-turnover pool after the annual release for each compartment
    ! (10 or 100 + 1 : input from year of land cover change)
    REAL, DIMENSION(npts,0:10), INTENT(inout)                   :: prod10_harvest
    REAL, DIMENSION(npts,0:100), INTENT(inout)                  :: prod100_harvest

    ! annual release from the 10/100 year-turnover pool compartments
    REAL, DIMENSION(npts,10), INTENT(inout)                     :: flux10_harvest
    REAL, DIMENSION(npts,100), INTENT(inout)                    :: flux100_harvest
                                                                                 !! @tex ($gC m^{-2}$) @endtex

    !! 0.4 Local variables

    INTEGER                                            :: i, j, k, l, m    !! indices (unitless)
    REAL,DIMENSION(npts,nelements)                     :: bm_new           !! biomass increase @tex ($gC m^{-2}$) @endtex
    REAL,DIMENSION(npts,nparts,nelements)              :: biomass_loss     !! biomass loss @tex ($gC m^{-2}$) @endtex
    REAL                                               :: above            !! aboveground biomass @tex ($gC m^{-2}$) @endtex
    REAL,DIMENSION(npts)                               :: dilu_KF         !! KF dilution 
    REAL,DIMENSION(npts)                               :: dilu_k_latosa_adapt        !! k_latosa_adapt dilution 
    REAL,DIMENSION(npts)                               :: dilu_rue_longterm        !! rue_longterm dilution 
    REAL,DIMENSION(npts,nlitt,nlevs,nelements)         :: dilu_lit         !! Litter dilution @tex ($gC m^{-2}$) @endtex
    REAL,DIMENSION(npts,ncarb,nelements)               :: dilu_som         !! SOM dilution @tex ($g(C or N) m^{-2}$) @endtex 
    REAL,DIMENSION(npts,nnspec)                        :: dilu_sin         !! Soil Inorganic Nitrogen dilution @tex ($gN m^{-2}$) @endtex     
  !! Soil Inorganic Nitrogen dilution @tex ($gN m^{-2}$) @endtex
    REAL,DIMENSION(npts,npspec)                        :: dilu_sin_p
!! Soil Inorganic Phosphorus dilution @tex ($gN m^{-2}$) @endtex
    REAL,DIMENSION(npts,nlevs)                         :: dilu_lf_struc    !! fraction of structural litter that is lignin
                                                                                  !! (0-1,unitless)
    REAL,DIMENSION(npts,nlevs)                         :: dilu_lf_wood     !! fraction of woody litter that is lignin
                                                                                  !! (0-1,unitless)
    REAL,DIMENSION(nvm)                                :: delta_veg        !! changes in "maximal" coverage fraction of PFT 
    REAL                                               :: delta_veg_sum    !! sum of delta_veg
    REAL,DIMENSION(npts,nvm)                           :: delta_ind        !! change in number of individuals  
    REAL                                               :: harvest          !! wood harvest renomalized by forest fraction
    REAL                                               :: relharvest       !! relative harvest fraction of each pft 
    REAL                                               :: relbiomass       !! relative biomass of each PFT 
    REAL                                               :: vegettree        !! vegfrac of trees
    REAL                                               :: abovemean        !! mean aboveground biomass @tex ($gC m^{-2}$) @endtex

!_ ================================================================================================================================

    IF (printlev>=3) WRITE(numout,*) 'Entering lcchange_main'
    
  !! 1. initialization
    
    prod10(:,0)         = zero
    prod100(:,0)        = zero   
    above               = zero
    convflux(:)         = zero
    cflux_prod10(:)     = zero
    cflux_prod100(:)    = zero
    delta_ind(:,:)      = zero
    delta_veg(:)        = zero
    nflux_prod_total(:) = zero 
    pflux_prod_total(:) = zero
 

    !! 2. calculation of changes in carbon stocks and biomass by wood harvest \n
    !! 2.1 initialization of carbon stocks\n
    prod10_harvest(:,0)           = zero
    prod100_harvest(:,0)          = zero   
    convflux_harvest(:)           = zero
    cflux_prod10_harvest(:)       = zero
    cflux_prod100_harvest(:)      = zero
    harvestpft(:,:)               = zero

  !! DSG: Wood harvest is now in own routine: wood_harvest
    
  !! 3. calculation of changes in carbon stocks and biomass by land cover change\n
    
    DO i = 1, npts ! Loop over # pixels - domain size
       
       !! 3.1 initialization of carbon stocks\n
       delta_veg(:) = veget_cov_max_new(i,:)-veget_cov_max_old(i,:)
       delta_veg_sum = SUM(delta_veg,MASK=delta_veg.LT.0.)
      
       dilu_lit(i,:,:,:) = zero
       dilu_som(i,:,:) = zero
       biomass_loss(i,:,:) = zero
       dilu_sin(i,:) = zero 
       dilu_sin_p(i,:) = zero
       dilu_lf_struc(i,:) = zero
       dilu_lf_wood(i,:) = zero
 
       dilu_KF(:)=zero
       dilu_k_latosa_adapt(:)=zero
       dilu_rue_longterm(:)=zero
       !! 3.2 if vegetation coverage decreases, compute dilution of litter, soil carbon, and biomass.\n
       DO j=2, nvm
          IF ( delta_veg(j) < -min_stomate ) THEN 
             dilu_lit(i,:,:,:) = dilu_lit(i,:,:,:) + delta_veg(j)*litter(i,:,j,:,:) / delta_veg_sum
             dilu_som(i,:,:) =  dilu_som(i,:,:) + delta_veg(j) * som(i,:,j,:) / delta_veg_sum 
             dilu_sin(i,:)=  dilu_sin(i,:) + delta_veg(j) * soil_n_min(i,j,:) / delta_veg_sum 
             dilu_sin_p(i,:)=  dilu_sin_p(i,:) + delta_veg(j) * soil_p_min(i,j,:) / delta_veg_sum
             dilu_lf_struc(i,:) = dilu_lf_struc(i,:) + &
                  delta_veg(j) * lignin_struc(i,j,:)* litter(i,istructural,j,:,icarbon) / delta_veg_sum
             dilu_lf_wood(i,:) = dilu_lf_wood(i,:) + &
                  delta_veg(j) * lignin_wood(i,j,:)*litter(i,iwoody,j,:,icarbon) / delta_veg_sum
             biomass_loss(i,:,:) = biomass_loss(i,:,:) + biomass(i,j,:,:)*delta_veg(j) / delta_veg_sum

             ! JC Debug I don't know why ther should be dilution of KF k_latosa_adapt and
             ! rue_longterm, they are PFT specific and should not be diluted
    !DSG         dilu_KF(i) = dilu_KF(i) + delta_veg(j) * KF(i,j) / delta_veg_sum
    !DSG         dilu_k_latosa_adapt(i) = dilu_k_latosa_adapt(i) + delta_veg(j) * k_latosa_adapt(i,j) / delta_veg_sum
    !DSG         dilu_rue_longterm(i) = dilu_rue_longterm(i) + delta_veg(j) * rue_longterm(i,j) / delta_veg_sum
          ENDIF
       ENDDO
       
       !! 3.3 
       DO j=2, nvm ! Loop over # PFTs

          !! 3.3.1 The case that vegetation coverage of PFTj increases
          IF ( delta_veg(j) > min_stomate) THEN

             !! 3.3.1.1 Initial setting of new establishment
             IF (veget_cov_max_old(i,j) .LT. min_stomate) THEN 

                ! JC Debug for new PFT, we don't establish here using bm_sapl
                ! from stomate_data. Instead, we will let stomate_prescribe
                ! to initialize vegetation biomass
                ! so it is necessary to make ind = 0 
                IF (is_tree(j)) THEN

                   ! cn_sapl(j)=0.5; stomate_data.f90
                   cn_ind(i,j) = cn_sapl(j) 
                ELSE
                   cn_ind(i,j) = un
                ENDIF
                ind(i,j)= delta_veg(j) / cn_ind(i,j)
                PFTpresent(i,j) = .TRUE.
                everywhere(i,j) = 1.
                senescence(i,j) = .FALSE.
                age(i,j) = zero

                ! large_value = 1.E33_r_std
                when_growthinit(i,j) = large_value 
                leaf_frac(i,j,:,:) = zero
                leaf_frac(i,j,:,1) = 1.0
                npp_longterm(i,j) = npp_longterm_init
                lm_lastyearmax(i,j) = bm_sapl(j,ileaf,icarbon) * ind(i,j)
                ind(i,j)=zero 
                ! JC Debug If the PFT has just been established, 
                ! initialize KF, k_latosa_adapt and rue_longterm
                KF(i,j) = k_latosa_min(j)
                k_latosa_adapt(i,j) = k_latosa_min(j)
                ! rue_longterm should be un as initialized (no gpp)
                rue_longterm(i,j) = un 

             ENDIF
             IF ( cn_ind(i,j) > min_stomate ) THEN
                delta_ind(i,j) = delta_veg(j) / cn_ind(i,j) 
             ENDIF
             
             !! 3.3.1.2 Update of biomass in each each carbon stock component 
             !!         Update of biomass in each each carbon stock component (leaf, sapabove, sapbelow,
             !>         heartabove, heartbelow, root, fruit, and carbres)\n
             DO k = 1, nparts ! loop over # carbon stock components, nparts = 8; stomate_constant.f90 
                DO l = 1,nelements ! loop over # elements

                   bm_new(i,l) = delta_ind(i,j) * bm_sapl(j,k,l) 
                   IF (veget_cov_max_old(i,j) .GT. min_stomate) THEN

                      ! in the case that bm_new is overestimated compared with biomass?
                      IF ((bm_new(i,l)/delta_veg(j)) > biomass(i,j,k,l)) THEN
                         bm_new(i,l) = biomass(i,j,k,l)*delta_veg(j)
                      ENDIF
                   ENDIF

                   IF (veget_cov_max_old(i,j) .GE. min_stomate) THEN 
                   biomass(i,j,k,l) = ( biomass(i,j,k,l) * veget_cov_max_old(i,j) + bm_new(i,l) ) / veget_cov_max_new(i,j)
                   co2_to_bm(i,j) = co2_to_bm(i,j) + (bm_new(i,icarbon)* dt_days) / (one_year * veget_cov_max_new(i,j))
                   n_to_bm(i,j) = n_to_bm(i,j) + (bm_new(i,initrogen)* dt_days) / (one_year * veget_cov_max_new(i,j))
                   p_to_bm(i,j) = p_to_bm(i,j) + (bm_new(i,iphosphorus)* dt_days) / (one_year * veget_cov_max_new(i,j))
                   ELSE
                   ! newly introduced PFT, it will be prescribed in prescribe
                   biomass(i,j,k,l) = zero

                   ENDIF ! (veget_cov_max_old(i,j) .GE. min_stomate)
                END DO ! loop over # elements
             ENDDO ! loop over # carbon stock components

             !! 3.3.1.3 Calculation of dilution in litter, soil carbon, and  input of litter
             !!        In this 'IF statement', dilu_* is zero. Formulas for litter and soil carbon
             !!         could be shortend?? Are the following formulas correct?

             ! JC Debug I don't know why ther should be dilution of KF k_latosa_adapt and
             ! rue_longterm, they are PFT specific and should not be diluted
             !!KF(i,j) = ( KF(i,j) * veget_cov_max_old(i,j) + &
             !!     dilu_KF(i) *  delta_veg(j)) / veget_cov_max_new(i,j)
             !!
             !!k_latosa_adapt(i,j) = ( k_latosa_adapt(i,j) * veget_cov_max_old(i,j) + &
             !!     dilu_k_latosa_adapt(i) *  delta_veg(j)) / veget_cov_max_new(i,j)
             !!
             !!rue_longterm(i,j) = ( rue_longterm(i,j) * veget_cov_max_old(i,j) + &
             !!     dilu_rue_longterm(i) *  delta_veg(j)) / veget_cov_max_new(i,j)

             ! Lignin fraction of structural litter
             lignin_struc(i,j,:)=(lignin_struc(i,j,:) * veget_cov_max_old(i,j)* litter(i,istructural,j,:,icarbon) + & 
                  dilu_lf_struc(i,:) * delta_veg(j)) / veget_cov_max_new(i,j) 

             ! Lignin fraction of woody litter
             lignin_wood(i,j,:)=(lignin_wood(i,j,:) * veget_cov_max_old(i,j)* litter(i,iwoody,j,:,icarbon) + & 
                  dilu_lf_wood(i,:) * delta_veg(j)) / veget_cov_max_new(i,j)

             ! Litter
             litter(i,:,j,:,:)=(litter(i,:,j,:,:) * veget_cov_max_old(i,j) + &
                  dilu_lit(i,:,:,:) * delta_veg(j)) / veget_cov_max_new(i,j)
                !gmjc available and not available litter for grazing
                ! only not available litter increase/decrease, available litter
                ! will not
                ! change, due to tree litter can not be eaten
               IF (is_grassland_manag(j) .AND. is_grassland_grazed(j)) THEN
                 litter_avail(i,:,j) = litter_avail(i,:,j) * veget_cov_max_old(i,j) / veget_cov_max_new(i,j)
                 litter_not_avail(i,:,j) = litter(i,:,j,iabove,icarbon) - litter_avail(i,:,j)
               ENDIF
                !end gmjc
            
             WHERE ( litter(i,istructural,j,:,icarbon) > min_stomate )
                lignin_struc(i,j,:) = lignin_struc(i,j,:)/litter(i,istructural,j,:,icarbon)
             ELSEWHERE
                lignin_struc(i,j,:) = LC_leaf(j)
             ENDWHERE


             WHERE ( litter(i,iwoody,j,:,icarbon) > min_stomate )
                lignin_wood(i,j,:) = lignin_wood(i,j,:)/litter(i,iwoody,j,:,icarbon)
             ELSEWHERE
                lignin_wood(i,j,:) = LC_heartabove(j)
             ENDWHERE

             ! Soil Organic Matter 
             som(i,:,j,:)=(som(i,:,j,:) * veget_cov_max_old(i,j) + dilu_som(i,:,:) * delta_veg(j)) / veget_cov_max_new(i,j) 
 
             ! Soil inorganic nitrogen  
             soil_n_min(i,j,:)=(soil_n_min(i,j,:) * veget_cov_max_old(i,j) + &  
                  dilu_sin(i,:) * delta_veg(j)) / veget_cov_max_new(i,j) 
             ! Soil inorganic phosphorus
             soil_p_min(i,j,:)=(soil_p_min(i,j,:) * veget_cov_max_old(i,j) + &
                  dilu_sin_p(i,:) * delta_veg(j)) / veget_cov_max_new(i,j)

             DO l = 1,nelements

                ! Litter input
                bm_to_litter(i,j,isapbelow,l) = (bm_to_litter(i,j,isapbelow,l) * veget_cov_max_old(i,j) + &
                     biomass_loss(i,isapbelow,l) * delta_veg(j)) / veget_cov_max_new(i,j)
                bm_to_litter(i,j,iheartbelow,l) = (bm_to_litter(i,j,iheartbelow,l) * veget_cov_max_old(i,j) + &
                     biomass_loss(i,iheartbelow,l) * delta_veg(j)) / veget_cov_max_new(i,j)
                bm_to_litter(i,j,iroot,l) = (bm_to_litter(i,j,iroot,l) * veget_cov_max_old(i,j) + &
                     biomass_loss(i,iroot,l) * delta_veg(j)) / veget_cov_max_new(i,j)
                bm_to_litter(i,j,ifruit,l) = (bm_to_litter(i,j,ifruit,l) * veget_cov_max_old(i,j) + &
                     biomass_loss(i,ifruit,l) * delta_veg(j)) / veget_cov_max_new(i,j)
                bm_to_litter(i,j,icarbres,l) = (bm_to_litter(i,j,icarbres,l) * veget_cov_max_old(i,j) + &
                     biomass_loss(i,icarbres,l) * delta_veg(j)) / veget_cov_max_new(i,j)
                bm_to_litter(i,j,ilabile,l) = (bm_to_litter(i,j,ilabile,l) * veget_cov_max_old(i,j) + &
                     biomass_loss(i,ilabile,l) * delta_veg(j)) / veget_cov_max_new(i,j)
                bm_to_litter(i,j,ileaf,l) = (bm_to_litter(i,j,ileaf,l) * veget_cov_max_old(i,j) + &
                     biomass_loss(i,ileaf,l) * delta_veg(j)) / veget_cov_max_new(i,j)

             END DO

             age(i,j)=age(i,j)*veget_cov_max_old(i,j)/veget_cov_max_new(i,j)
             
          !! 3.3.2 The case that vegetation coverage of PFTj is no change or decreases
          ELSE 
 
             !! 3.3.2.1 Biomass export
             ! coeff_lcchange_*:  Coeff of biomass export for the year, decade, and century
             above = biomass(i,j,isapabove,icarbon) + biomass(i,j,iheartabove,icarbon)
             convflux(i)  = convflux(i)  - ( coeff_lcchange_1(j) * above * delta_veg(j) ) 
             prod10(i,0)  = prod10(i,0)  - ( coeff_lcchange_10(j) * above * delta_veg(j) )
             prod100(i,0) = prod100(i,0) - ( coeff_lcchange_100(j) * above * delta_veg(j) )
             above = biomass(i,j,isapabove,initrogen) + biomass(i,j,iheartabove,initrogen) 
             nflux_prod_total(i) = nflux_prod_total(i) - above* delta_veg(j) 
             above = biomass(i,j,isapabove,iphosphorus) + biomass(i,j,iheartabove,iphosphorus)
             pflux_prod_total(i) = pflux_prod_total(i) - above* delta_veg(j)

             !! 3.3.2.2 Total reduction
             !! If the vegetation is to small, it has been set to 0.
             IF ( veget_cov_max_new(i,j) .LT. min_stomate ) THEN 
                
                ind(i,j) = zero
                biomass(i,j,:,:) = zero
                PFTpresent(i,j) = .FALSE.
                senescence(i,j) = .FALSE.
                age(i,j) = zero
                when_growthinit(i,j) = undef
                everywhere(i,j) = zero
                litter(i,:,j,:,:) = zero
                bm_to_litter(i,j,:,:) = zero
                turnover_daily(i,j,:,:) = zero
                som(i,:,j,:) = zero
                soil_n_min(i,j,:) = zero
                soil_p_min(i,j,:) = zero               
             ENDIF
 
          ENDIF ! End if PFT's coverage reduction
          
       ENDDO ! Loop over # PFTs
       
       !! 3.4 update 10 year-turnover pool content following flux emission
       !!     (linear decay (10%) of the initial carbon input)
       DO  l = 0, 8
          m = 10 - l
          cflux_prod10(i) =  cflux_prod10(i) + flux10(i,m)
          prod10(i,m)     =  prod10(i,m-1)   - flux10(i,m-1)
          flux10(i,m)     =  flux10(i,m-1)
          
          IF (prod10(i,m) .LT. 1.0) prod10(i,m) = zero
          ! Similar treatment for wood harvest
          IF (do_wood_harvest) THEN
             cflux_prod10_harvest(i) =  cflux_prod10_harvest(i) + flux10_harvest(i,m)
             prod10_harvest(i,m)   =  prod10_harvest(i,m-1)   - flux10_harvest(i,m-1)
             flux10_harvest(i,m)   =  flux10(i,m-1)
             IF (prod10_harvest(i,m) .LT. 1.0) prod10_harvest(i,m) = zero
          ENDIF

       ENDDO
       
       cflux_prod10(i) = cflux_prod10(i) + flux10(i,1) 
       flux10(i,1)     = 0.1 * prod10(i,0)
       prod10(i,1)     = prod10(i,0)
       ! Similar treatment for wood harvest as for LU
       IF (do_wood_harvest) THEN
          cflux_prod10_harvest(i) = cflux_prod10_harvest(i) + flux10_harvest(i,1) 
          flux10_harvest(i,1)     = 0.1 * prod10_harvest(i,0)
          prod10_harvest(i,1)     = prod10_harvest(i,0)
       ENDIF
       
       !! 3.5 update 100 year-turnover pool content following flux emission\n
       DO   l = 0, 98
          m = 100 - l
          cflux_prod100(i)  =  cflux_prod100(i) + flux100(i,m)
          prod100(i,m)      =  prod100(i,m-1)   - flux100(i,m-1)
          flux100(i,m)      =  flux100(i,m-1)
          
          IF (prod100(i,m).LT.1.0) prod100(i,m) = zero
          ! Similar treatment for wood harvest as for LU
          IF (do_wood_harvest) THEN
             cflux_prod100_harvest(i)  =  cflux_prod100_harvest(i) + flux100_harvest(i,m)
             prod100_harvest(i,m)      =  prod100_harvest(i,m-1)   - flux100_harvest(i,m-1)
             flux100_harvest(i,m)      =  flux100_harvest(i,m-1)
             IF (prod100_harvest(i,m).LT.1.0) prod100_harvest(i,m) = zero
          ENDIF
       ENDDO
       
       cflux_prod100(i)  = cflux_prod100(i) + flux100(i,1) 
       flux100(i,1)      = 0.01 * prod100(i,0)
       prod100(i,1)      = prod100(i,0)
       prod10(i,0)        = zero
       prod100(i,0)       = zero 

       ! Similar treatment for wood harvest as for LU
       IF (do_wood_harvest) THEN
          cflux_prod100_harvest(i)  = cflux_prod100_harvest(i) + flux100_harvest(i,1) 
          flux100_harvest(i,1)      = 0.01 * prod100_harvest(i,0)
          prod100_harvest(i,1)      = prod100_harvest(i,0)
          prod10_harvest(i,0)        = 0.0
          prod100_harvest(i,0)       = 0.0 
       ENDIF
       
    ENDDO ! Loop over # pixels - domain size
    
  !! 4. history
    convflux        = convflux/one_year*dt_days
    cflux_prod10    = cflux_prod10/one_year*dt_days
    cflux_prod100   = cflux_prod100/one_year*dt_days

    ! Similar treatment for wood harvest as for LU
    IF (do_wood_harvest) THEN
       convflux_harvest        = convflux_harvest/one_year*dt_days
       cflux_prod10_harvest    = cflux_prod10_harvest/one_year*dt_days
       cflux_prod100_harvest   = cflux_prod100_harvest/one_year*dt_days
    ENDIF
   
    IF (printlev>=4) WRITE(numout,*) 'Leaving lcchange_main'
    
  END SUBROUTINE lcchange_main
  
END MODULE stomate_lcchange
