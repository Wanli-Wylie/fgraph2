! =================================================================================================================================
! MODULE       : lpj_cover
!
! CONTACT      : orchidee-help _at_ listes.ipsl.fr
!
! LICENCE      : IPSL (2006)
!                This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF        Recalculate vegetation cover and LAI
!!
!!\n DESCRIPTION : None
!!
!! RECENT CHANGE(S) : None
!!
!! REFERENCE(S) : 
!!        Sitch, S., B. Smith, et al. (2003), Evaluation of ecosystem dynamics,
!!        plant geography and terrestrial carbon cycling in the LPJ dynamic 
!!        global vegetation model, Global Change Biology, 9, 161-185.\n
!!        Smith, B., I. C. Prentice, et al. (2001), Representation of vegetation
!!        dynamics in the modelling of terrestrial ecosystems: comparing two
!!        contrasting approaches within European climate space,
!!        Global Ecology and Biogeography, 10, 621-637.\n
!!
!! SVN :
!! $HeadURL: svn://forge.ipsl.jussieu.fr/orchidee/branches/ORCHIDEE-CN-P/ORCHIDEE/src_stomate/lpj_cover.f90 $
!! $Date: 2017-11-29 10:57:22 +0100 (三, 2017-11-29) $
!! $Revision: 4790 $
!! \n
!_ ================================================================================================================================

MODULE lpj_cover

  ! modules used:

  USE ioipsl_para
  USE xios_orchidee
  USE stomate_data
  USE pft_parameters

  IMPLICIT NONE

  ! private & public routines

  PRIVATE
  PUBLIC cover

CONTAINS

!! ================================================================================================================================
!! SUBROUTINE     : lpj_cover
!!
!>\BRIEF          Recalculate vegetation cover and LAI
!!
!!\n DESCRIPTION : Veget_cov_max is first renewed here according to newly calculated foliage biomass in this calculation step 
!! Then, litter, soil carbon, and biomass are also recalcuted with taking into account the changes in Veget_cov_max (i.e. delta_veg)
!! Grid-scale fpc (foliage projected coverage) is calculated to obtain the shadede ground area by leaf's light capture
!! Finally, grid-scale fpc is adjusted not to exceed 1.0
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S) : ::lai (leaf area index, @tex $(m^2 m^{-2})$ @endtex), 
!! :: veget (fractional vegetation cover, unitless)
!!
!! REFERENCE(S)   : None
!! 
!! FLOWCHART :
!! \latexonly 
!!     \includegraphics[scale=0.5]{lpj_cover_flowchart.png}
!! \endlatexonly
!! \n
!_ ================================================================================================================================

  SUBROUTINE cover (npts, cn_ind, ind, biomass, &
       veget_cov_max, veget_cov_max_old, lai, litter, &
!gmjc
       litter_avail, litter_not_avail, &
!end gmjc
       som, turnover_daily, bm_to_litter, &
       co2_to_bm, co2_fire, resp_hetero, resp_maint, resp_growth, gpp_daily,&
       lignin_struc, lignin_wood, soil_n_min, soil_p_min)

!! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER, INTENT(in)                                  :: npts             !! Domain size (unitless)  
    REAL, DIMENSION(npts,nvm), INTENT(in)                :: cn_ind           !! Crown area 
                                                                                    !! @tex $(m^2)$ @endtex per individual
    REAL, DIMENSION(npts,nvm), INTENT(in)                :: ind              !! Number of individuals 
                                                                                    !! @tex $(m^{-2})$ @endtex
    REAL, DIMENSION(npts,nvm), INTENT(in)                :: veget_cov_max_old!! "Maximal" coverage fraction of a PFT (LAI-> 
                                                                                    !! infinity) on ground at beginning of time 

    !! 0.2 Output variables

    !! 0.3 Modified variables

    REAL, DIMENSION(npts,nvm), INTENT(inout)             :: lai                 !! Leaf area index OF AN INDIVIDUAL PLANT 
                                                                                       !! @tex $(m^2 m^{-2})$ @endtex
    REAL, DIMENSION(npts,nlitt,nvm,nlevs,nelements), INTENT(inout) :: litter    !! Metabolic and structural litter, above and 
                                                                                       !! below ground @tex $(gC m^{-2})$ @endtex
!gmjc
    REAL, DIMENSION(npts,nlitt,nvm), INTENT(inout):: litter_avail
    REAL, DIMENSION(npts,nlitt,nvm) , INTENT(inout):: litter_not_avail
!end gmjc
    REAL, DIMENSION(npts,ncarb,nvm,nelements), INTENT(inout) :: som            !! Carbon pool: active, slow, or passive @tex $(gC m^{-2})$ @endtex
    REAL, DIMENSION(npts,nvm,nparts,nelements), INTENT(inout) :: biomass        !! Biomass @tex $(gC m^{-2})$ @endtex
    REAL, DIMENSION(npts,nvm), INTENT(inout)                  :: veget_cov_max  !! "Maximal" coverage fraction of a PFT (LAI->
                                                                                       !! infinity) on ground (unitless)
    REAL, DIMENSION(npts,nvm,nparts,nelements), INTENT(inout) :: turnover_daily !! Turnover rates (gC m^{-2} day^{-1})
    REAL, DIMENSION(npts,nvm,nparts,nelements), INTENT(inout) :: bm_to_litter   !! Conversion of biomass to litter 
                                                                                       !! @tex $(gC m^{-2} day^{-1})$ @endtex
    REAL, DIMENSION(npts,nvm), INTENT(inout)                  :: co2_to_bm      !! biomass up take for establishment           
    REAL, DIMENSION(npts,nvm), INTENT(inout)                  :: co2_fire       !! Carbon emitted to the atmosphere by fire(living
                                                                                       !! and dead biomass)
    REAL, DIMENSION(npts,nvm), INTENT(inout)                  :: resp_hetero    !! Heterotrophic respiration
    REAL, DIMENSION(npts,nvm), INTENT(inout)                  :: resp_maint     !! Maintenance respiration 
    REAL, DIMENSION(npts,nvm), INTENT(inout)                  :: resp_growth    !! Growth respiration
    REAL, DIMENSION(npts,nvm), INTENT(inout)                  :: gpp_daily      !! Daily gross primary productivity
    REAL, DIMENSION(npts,nvm,nlevs), INTENT(inout)            :: lignin_struc   !! ratio Lignine/Carbon in structural litter,
                                                                                       !! above and below ground
    REAL, DIMENSION(npts,nvm,nlevs), INTENT(inout)            :: lignin_wood    !! ratio Lignine/Carbon in woody litter,
                                                                                       !! above and below ground
    REAL,DIMENSION(npts,nvm,nnspec), INTENT(inout)            :: soil_n_min     !! mineral nitrogen in the soil (gN/m**2)  
    REAL,DIMENSION(npts,nvm,npspec), INTENT(inout)            :: soil_p_min     !! mineral phosphorus in the soil (gP/m**2)
    !! 0.4 Local variables

    INTEGER                                              :: i,j,k,m             !! Index (unitless)
    REAL, DIMENSION(npts,nlitt,nlevs,nelements)          :: dilu_lit            !! Litter dilution @tex $(gC m^{-2})$ @endtex
    REAL, DIMENSION(npts,ncarb,nelements)                :: dilu_som            !! Soil Organic Matter ndilution 
                                                                                       !! @tex $(gC(or N) m^{-2})$ @endtex
    REAL, DIMENSION(npts,nnspec)                         :: dilu_sin            !! Soil Inorganic Nitrogen dilution 
                                                                                       !! @tex $(gN m^{-2})$ @endtex 
    REAL,DIMENSION(npts,npspec)                        :: dilu_sin_p
!! Soil Inorganic Phosphorus dilution @tex ($gN m^{-2}$) @endtex
    REAL, DIMENSION(npts,nlevs)                          :: dilu_lf_struc       !! fraction of structural litter that is lignin
                                                                                       !! (0-1,unitless)
    REAL, DIMENSION(npts,nlevs)                          :: dilu_lf_wood        !! fraction of woody litter that is lignin
                                                                                       !! (0-1,unitless)
    REAL, DIMENSION(npts,nparts,nelements)               :: dilu_bio            !! Dilution for carbon variables
    REAL, DIMENSION(npts)                                :: dilu_TCarbon        !! Dilution for carbon variables
    REAL, DIMENSION(npts,nparts,nelements)               :: dilu_turnover_daily !! Dilution for carbon variables
    REAL, DIMENSION(npts,nparts,nelements)               :: dilu_bm_to_litter   !! Dilution for carbon variables
    REAL, DIMENSION(npts)                                :: dilu_co2flux_new    !! Dilution for carbon variables
    REAL, DIMENSION(npts)                                :: dilu_gpp_daily      !! Dilution for carbon variables
    REAL, DIMENSION(npts)                                :: dilu_resp_growth    !! Dilution for carbon variables
    REAL, DIMENSION(npts)                                :: dilu_resp_maint     !! Dilution for carbon variables
    REAL, DIMENSION(npts)                                :: dilu_resp_hetero    !! Dilution for carbon variables
    REAL, DIMENSION(npts)                                :: dilu_co2_to_bm      !! Dilution for carbon variables
    REAL, DIMENSION(npts)                                :: dilu_co2_fire       !! Dilution for carbon variables
    REAL, DIMENSION(npts,nvm)                            :: TCarbon             !! Total carbon
    REAL, DIMENSION(npts,nvm)                            :: co2flux_new         !! NBP after re-calculation in order to conserve carbon 
    REAL, DIMENSION(npts,nvm)                            :: co2flux_old         !! NBP before re-calculation 
    REAL, DIMENSION(nvm)                                 :: delta_veg           !! Conversion factors (unitless)
    REAL, DIMENSION(nvm)                                 :: reduct              !! Conversion factors (unitless)
    REAL                                                 :: delta_veg_sum       !! Conversion factors (unitless)
    REAL                                                 :: diff                !! Conversion factors (unitless)
    REAL                                                 :: sr                  !! Conversion factors (unitless)
    REAL, DIMENSION(npts)                                :: frac_nat            !! Conversion factors (unitless)
    REAL, DIMENSION(npts)                                :: sum_vegettree       !! Conversion factors (unitless)
    REAL, DIMENSION(npts)                                :: sum_vegetgrass      !! Conversion factors (unitless) 
    REAL, DIMENSION(npts)                                :: sum_veget_natveg    !! Conversion factors (unitless)
    REAL, DIMENSION(npts)                                :: vartmp              !! Temporary variable used to add history

!_ ================================================================================================================================

 !! 1. If the vegetation is dynamic, calculate new maximum vegetation cover for natural plants
  
    IF ( ok_dgvm ) THEN

       !! 1.1  Calculate initial values of vegetation cover
       frac_nat(:) = un
       sum_veget_natveg(:) = zero
       veget_cov_max(:,ibare_sechiba) = un

       DO j = 2,nvm ! loop over PFTs

          IF ( natural(j) ) THEN
	     
             ! Summation of individual tree crown area to get total foliar projected coverage
             veget_cov_max(:,j) = ind(:,j) * cn_ind(:,j)
             sum_veget_natveg(:) = sum_veget_natveg(:) + veget_cov_max(:,j)

          ELSE
             
             !fraction occupied by agriculture needs to be substracted for the DGVM
             !this is used below to constrain veget for natural vegetation, see below
             frac_nat(:) = frac_nat(:) - veget_cov_max(:,j)

          ENDIF

       ENDDO ! loop over PFTs

       DO i = 1, npts ! loop over grid points
  
          ! Recalculation of vegetation projected coverage when ::frac_nat was below ::sum_veget_natveg
          ! It means that non-natural vegetation will recover ::veget_cov_max as natural vegetation
          IF (sum_veget_natveg(i) .GT. frac_nat(i) .AND. frac_nat(i) .GT. min_stomate) THEN

             DO j = 2,nvm ! loop over PFTs
                IF( natural(j) ) THEN
                   veget_cov_max(i,j) =  veget_cov_max(i,j) * frac_nat(i) / sum_veget_natveg(i)
                ENDIF
             ENDDO ! loop over PFTs

          ENDIF
       ENDDO ! loop over grid points

       ! Renew veget_cov_max of bare soil as 0 to difference of veget_cov_max (ibare_sechiba) 
       ! to current veget_cov_max
       DO j = 2,nvm ! loop over PFTs
          veget_cov_max(:,ibare_sechiba) = veget_cov_max(:,ibare_sechiba) - veget_cov_max(:,j)
       ENDDO ! loop over PFTs
       veget_cov_max(:,ibare_sechiba) = MAX( veget_cov_max(:,ibare_sechiba), zero )

       !! 1.2 Calculate carbon fluxes between PFTs to maintain mass balance
       !! Assure carbon closure when veget_cov_max changes(delta_veg): if veget_cov_max of some PFTs decrease, we use "dilu" to 
       !! record the corresponding lost in carbon (biomass, litter, soil carbon, gpp, respiration etc.) for 
       !! these PFTs, and re-allocate "dilu" to those PFTs with increasing veget_cov_max.
       DO i = 1, npts ! loop over grid points
          
          ! Calculate the change in veget_cov_max between previous time step and current time step
          delta_veg(:) = veget_cov_max(i,:)-veget_cov_max_old(i,:)
          delta_veg_sum = SUM(delta_veg,MASK=delta_veg.LT.zero)

          dilu_lit(i,:,:,:) = zero
          dilu_som(i,:,:) = zero
          dilu_sin(i,:) = zero
          dilu_sin_p(i,:) = zero
          dilu_lf_struc(i,:) = zero
          dilu_lf_wood(i,:) = zero
          dilu_bio(i,:,:) = zero
          dilu_TCarbon(i)=zero
          dilu_turnover_daily(i,:,:)=zero
          dilu_bm_to_litter(i,:,:)=zero
          dilu_co2flux_new(i)=zero
          dilu_gpp_daily(i)=zero
          dilu_resp_growth(i)=zero
          dilu_resp_maint(i)=zero
          dilu_resp_hetero(i)=zero
          dilu_co2_to_bm(i)=zero
          dilu_co2_fire(i)=zero

          ! Calculate TCarbon: total carbon including biomass, litter and soil carbon, as well as "today's" turnover and 
          ! bm_to_litter due to mortality, because today's turnover and bm_to_litter are not yet added into "litter" until tomorrow. 
          DO j=1, nvm
                TCarbon(i,j)=SUM(biomass(i,j,:,icarbon))+SUM(som(i,:,j,icarbon))+SUM(litter(i,:,j,:,icarbon))+&
                     SUM(turnover_daily(i,j,:,icarbon))+SUM(bm_to_litter(i,j,:,icarbon))
                co2flux_old(i,j)=resp_maint(i,j)+resp_growth(i,j)+resp_hetero(i,j)+co2_fire(i,j)-co2_to_bm(i,j)-gpp_daily(i,j)
                co2flux_new(i,j)=resp_maint(i,j)+resp_growth(i,j)+resp_hetero(i,j)+co2_fire(i,j)-co2_to_bm(i,j)-gpp_daily(i,j)
          ENDDO

          DO j=1, nvm ! loop over PFTs
             IF ( delta_veg(j) < -min_stomate ) THEN 
                dilu_lit(i,:,:,:) =  dilu_lit(i,:,:,:) + delta_veg(j) * litter(i,:,j,:,:) / delta_veg_sum
                dilu_som(i,:,:) =  dilu_som(i,:,:) + delta_veg(j) * som(i,:,j,:) / delta_veg_sum
                dilu_sin(i,:) =  dilu_sin(i,:) + delta_veg(j) * soil_n_min(i,j,:) / delta_veg_sum 
                dilu_sin_p(i,:)=  dilu_sin_p(i,:) + delta_veg(j) * soil_p_min(i,j,:) / delta_veg_sum
                dilu_lf_struc(i,:) = dilu_lf_struc(i,:) + & 
                     delta_veg(j) * lignin_struc(i,j,:) * litter(i,istructural,j,:,icarbon)/ delta_veg_sum
                dilu_lf_wood(i,:) = dilu_lf_wood(i,:) + &
                     delta_veg(j) * lignin_wood(i,j,:)*litter(i,iwoody,j,:,icarbon)  / delta_veg_sum
                dilu_TCarbon(i)= dilu_TCarbon(i) + delta_veg(j) * TCarbon(i,j) / delta_veg_sum
                dilu_turnover_daily(i,:,:)=dilu_turnover_daily(i,:,:)+delta_veg(j)*turnover_daily(i,j,:,:)/delta_veg_sum
                dilu_bm_to_litter(i,:,:)=dilu_bm_to_litter(i,:,:)+delta_veg(j)*bm_to_litter(i,j,:,:)/delta_veg_sum
                dilu_co2flux_new(i)=dilu_co2flux_new(i)+delta_veg(j)*co2flux_old(i,j)/delta_veg_sum
                dilu_gpp_daily(i)=dilu_gpp_daily(i)+delta_veg(j)*gpp_daily(i,j)/delta_veg_sum
                dilu_resp_growth(i)=dilu_resp_growth(i)+delta_veg(j)*resp_growth(i,j)/delta_veg_sum
                dilu_resp_maint(i)=dilu_resp_maint(i)+delta_veg(j)*resp_maint(i,j)/delta_veg_sum
                dilu_resp_hetero(i)=dilu_resp_hetero(i)+delta_veg(j)*resp_hetero(i,j)/delta_veg_sum
                dilu_co2_to_bm(i)=dilu_co2_to_bm(i)+delta_veg(j)*co2_to_bm(i,j)/delta_veg_sum
                dilu_co2_fire(i)=dilu_co2_fire(i)+delta_veg(j)*co2_fire(i,j)/delta_veg_sum
             ENDIF
          ENDDO ! loop over PFTs

          DO j=1, nvm ! loop over PFTs
             IF ( delta_veg(j) > min_stomate) THEN

                ! Dilution of reservoirs
                ! Recalculate the litter and soil carbon with taking into accout the change in 
                ! veget_cov_max (delta_veg)
                ! Litter
                ! Lignin fraction of structural litter
                lignin_struc(i,j,:)=(lignin_struc(i,j,:) * veget_cov_max_old(i,j)* litter(i,istructural,j,:,icarbon) + & 
                     dilu_lf_struc(i,:) * delta_veg(j)) / veget_cov_max(i,j) 

                ! Lignin fraction of woody litter
                lignin_wood(i,j,:)=(lignin_wood(i,j,:) * veget_cov_max_old(i,j)* litter(i,iwoody,j,:,icarbon) + & 
                     dilu_lf_wood(i,:) * delta_veg(j)) / veget_cov_max(i,j)

                litter(i,:,j,:,:)=(litter(i,:,j,:,:) * veget_cov_max_old(i,j) + dilu_lit(i,:,:,:) * delta_veg(j)) / veget_cov_max(i,j)

                WHERE ( litter(i,istructural,j,:,icarbon) > min_stomate )
                   lignin_struc(i,j,:) = lignin_struc(i,j,:)/litter(i,istructural,j,:,icarbon)
                ELSEWHERE
                   lignin_struc(i,j,:) = LC_leaf(j)
                ENDWHERE

                WHERE ( litter(i,iwoody,j,:,icarbon) > min_stomate )
                   lignin_wood(i,j,:) = lignin_struc(i,j,:)/litter(i,iwoody,j,:,icarbon)
                ELSEWHERE
                   lignin_wood(i,j,:) = LC_heartabove(j)
                ENDWHERE
                !gmjc available and not available litter for grazing
                ! only not available litter change, available litter will not
                ! change, because tree litter can not be eaten
               IF (is_grassland_manag(j) .AND. is_grassland_grazed(j)) THEN
                 litter_avail(i,:,j) = litter_avail(i,:,j) * veget_cov_max_old(i,j) / veget_cov_max(i,j)
                 litter_not_avail(i,:,j) = litter(i,:,j,iabove,icarbon) - litter_avail(i,:,j)
               ENDIF
                !end gmjc
                ! Soil carbon
                som(i,:,j,:)=(som(i,:,j,:) * veget_cov_max_old(i,j) + dilu_som(i,:,:) * delta_veg(j)) / veget_cov_max(i,j)
                ! Soil inorganic nitrogen  
                soil_n_min(i,j,:)=(soil_n_min(i,j,:) * veget_cov_max_old(i,j) + dilu_sin(i,:) &  
                     * delta_veg(j)) / veget_cov_max(i,j) 
                ! Soil inorganic phosphorus
                soil_p_min(i,j,:)=(soil_p_min(i,j,:) * veget_cov_max_old(i,j) + &
                     dilu_sin_p(i,:) * delta_veg(j)) / veget_cov_max(i,j)
                TCarbon(i,j)=(TCarbon(i,j) * veget_cov_max_old(i,j) + dilu_TCarbon(i) * delta_veg(j)) / veget_cov_max(i,j)
                turnover_daily(i,j,:,:)=(turnover_daily(i,j,:,:)*veget_cov_max_old(i,j)+&
                     dilu_turnover_daily(i,:,:)*delta_veg(j))/veget_cov_max(i,j)
                bm_to_litter(i,j,:,:)=(bm_to_litter(i,j,:,:)*veget_cov_max_old(i,j)+&
                     dilu_bm_to_litter(i,:,:)*delta_veg(j))/veget_cov_max(i,j)
                co2flux_new(i,j)=(co2flux_old(i,j)*veget_cov_max_old(i,j)+dilu_co2flux_new(i)*delta_veg(j))/veget_cov_max(i,j)
                gpp_daily(i,j)=(gpp_daily(i,j)*veget_cov_max_old(i,j)+dilu_gpp_daily(i)*delta_veg(j))/veget_cov_max(i,j)
                resp_growth(i,j)=(resp_growth(i,j)*veget_cov_max_old(i,j)+dilu_resp_growth(i)*delta_veg(j))/veget_cov_max(i,j)
                resp_maint(i,j)=(resp_maint(i,j)*veget_cov_max_old(i,j)+dilu_resp_maint(i)*delta_veg(j))/veget_cov_max(i,j)
                resp_hetero(i,j)=(resp_hetero(i,j)*veget_cov_max_old(i,j)+dilu_resp_hetero(i)*delta_veg(j))/veget_cov_max(i,j)
                co2_to_bm(i,j)=(co2_to_bm(i,j)*veget_cov_max_old(i,j)+dilu_co2_to_bm(i)*delta_veg(j))/veget_cov_max(i,j)
                co2_fire(i,j)=(co2_fire(i,j)*veget_cov_max_old(i,j)+dilu_co2_fire(i)*delta_veg(j))/veget_cov_max(i,j)
             ENDIF

             IF(veget_cov_max(i,j).GT.min_stomate) THEN

                ! Correct biomass densities to conserve mass
                ! since it's defined on veget_cov_max
                biomass(i,j,:,:) = biomass(i,j,:,:) * veget_cov_max_old(i,j) / veget_cov_max(i,j)

             ENDIF

          ENDDO ! loop over PFTs
       ENDDO ! loop over grid points
    ELSE !DSGdebug: the following variables are in the IO also when ok_dgvm-FALSE; so set them
       DO i = 1, npts ! loop over grid points
          Do j=1, nvm ! loop over PFTs
                TCarbon(i,j)=SUM(biomass(i,j,:,icarbon))+SUM(som(i,:,j,icarbon))+SUM(litter(i,:,j,:,icarbon))+&
                     SUM(turnover_daily(i,j,:,icarbon))+SUM(bm_to_litter(i,j,:,icarbon))
                co2flux_old(i,j)=resp_maint(i,j)+resp_growth(i,j)+resp_hetero(i,j)+co2_fire(i,j)-co2_to_bm(i,j)-gpp_daily(i,j)
                co2flux_new(i,j)=resp_maint(i,j)+resp_growth(i,j)+resp_hetero(i,j)+co2_fire(i,j)-co2_to_bm(i,j)-gpp_daily(i,j)
          ENDDO
       ENDDO
       !DSGdebug
    ENDIF

    vartmp(:)=SUM(gpp_daily*veget_cov_max,dim=2)
    CALL histwrite_p (hist_id_stomate, "tGPP", itime, vartmp, npts, hori_index)
    CALL xios_orchidee_send_field("tGPP",vartmp(:))

    vartmp(:)=SUM(resp_growth*veget_cov_max,dim=2)
    CALL histwrite_p (hist_id_stomate, "tRESP_GROWTH", itime, vartmp, npts, hori_index)
    CALL xios_orchidee_send_field("tRESP_GROWTH",vartmp(:))

    vartmp(:)=SUM(resp_maint*veget_cov_max,dim=2)
    CALL histwrite_p (hist_id_stomate, "tRESP_MAINT", itime, vartmp, npts, hori_index)
    CALL xios_orchidee_send_field("tRESP_MAINT",vartmp(:))

    vartmp(:)=SUM(resp_hetero*veget_cov_max,dim=2)
    CALL histwrite_p (hist_id_stomate, "tRESP_HETERO", itime, vartmp, npts, hori_index)
    CALL xios_orchidee_send_field("tRESP_HETERO",vartmp(:))

    vartmp(:)=SUM(co2_to_bm*veget_cov_max,dim=2)
    CALL histwrite_p (hist_id_stomate, "tCO2_TAKEN", itime, vartmp, npts, hori_index)
    CALL xios_orchidee_send_field("tCO2_TAKEN",vartmp(:))

    vartmp(:)=SUM(co2_fire*veget_cov_max,dim=2)
    CALL histwrite_p (hist_id_stomate, "tCO2_FIRE", itime, vartmp, npts, hori_index)
    CALL xios_orchidee_send_field("tCO2_FIRE",vartmp(:))

    vartmp(:)=SUM(co2flux_new*veget_cov_max,dim=2)
    CALL histwrite_p (hist_id_stomate, "tCO2FLUX", itime, vartmp, npts, hori_index)
    CALL xios_orchidee_send_field("tCO2FLUX",vartmp(:))

    vartmp(:)=SUM(co2flux_old*veget_cov_max_old,dim=2)
    CALL histwrite_p (hist_id_stomate, "tCO2FLUX_OLD", itime, vartmp, npts, hori_index)
    CALL xios_orchidee_send_field("tCO2FLUX_OLD",vartmp(:))

    vartmp(:)=SUM(TCarbon*veget_cov_max,dim=2)
    CALL histwrite_p (hist_id_stomate, "tCARBON", itime, vartmp, npts, hori_index)
    CALL xios_orchidee_send_field("tCARBON",vartmp(:))

    vartmp(:)=SUM(SUM(biomass(:,:,:,icarbon),dim=3)*veget_cov_max,dim=2)
    CALL histwrite_p (hist_id_stomate, "tBIOMASS", itime, vartmp, npts, hori_index)
    CALL xios_orchidee_send_field("tBIOMASS",vartmp(:))

    vartmp(:)=SUM(SUM(SUM(litter(:,:,:,:,icarbon),dim=4),dim=2)*veget_cov_max,dim=2)
    CALL histwrite_p (hist_id_stomate, "tLITTER", itime, vartmp, npts, hori_index)
    CALL xios_orchidee_send_field("tLITTER",vartmp(:))

    vartmp(:)=SUM(SUM(som(:,:,:,icarbon),dim=2)*veget_cov_max,dim=2)
    CALL histwrite_p (hist_id_stomate, "tSOILC", itime, vartmp, npts, hori_index)
    CALL xios_orchidee_send_field("tSOILC",vartmp(:))

  END SUBROUTINE cover

END MODULE lpj_cover
