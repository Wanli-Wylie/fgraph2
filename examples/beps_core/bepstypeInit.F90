! This module is used to initialize bepstype variables
! Editted by J.Wang
! Date: 10May2017

module bepstypeInit
   use shr_kind_mod,only: r8=>shr_kind_r8
   use beps_par,only:npoints,PFT
   use bepstype
   implicit none
   save

   public  :: Initbepstype
#ifdef COUP_CSM
   private :: initlnd2atm
#endif
   private :: initatm2lnd
   private :: InitSurf
   private :: InitOuput
contains

   subroutine Initbepstype()
      implicit none

      allocate(v2last(npoints,0:40,PFT))
      v2last  = 0.

      call initatm2lnd()
#ifdef COUP_CSM
      call initlnd2atm()
#endif

      call InitSurf()

      call InitSoilstat()

      call InitOuput()
      return
   end subroutine


   subroutine initatm2lnd()
      implicit none
      type(forc),pointer ::p

      p=>clim
      allocate(p%Temp(npoints))
      allocate(p%Tempmx(npoints))
      allocate(p%Tempmn(npoints))
      allocate(p%Wind(npoints))

#ifdef COUP_CSM
      allocate(p%Zref(npoints))
      allocate(p%Rain(npoints))
      allocate(p%Snow(npoints))
      allocate(p%Swndr(npoints))
      allocate(p%Swvdr(npoints))
      allocate(p%Swndf(npoints))
      allocate(p%Swvdf(npoints))
      allocate(p%Swdr(npoints))
      allocate(p%Swdf(npoints))
      allocate(p%Lwdn(npoints))
      allocate(p%shum(npoints))
      allocate(p%pres(npoints))
#else
      allocate(p%Srad(npoints))
      allocate(p%Rh(npoints))
      allocate(p%Rain(npoints))
      allocate(p%Snow(npoints))
      allocate(p%Swdr(npoints))
      allocate(p%Swdf(npoints))
#endif

      p%Temp(:)      = 0.
      p%Tempmx(:)      = 0.
      p%Tempmn(:)      = 0.
      p%Wind(:)      = 0.
#ifdef COUP_CSM
      p%Zref(:)     = 0.
      p%Rain(:)     = 0.
      p%Snow(:)     = 0.
      p%Swndr(:)    = 0.
      p%Swvdr(:)    = 0.
      p%Swndf(:)    = 0.
      p%Swvdf(:)    = 0.
      p%Swdr(:)     = 0.
      p%Swdf(:)     = 0.
      p%Lwdn(:)     = 0.
      p%shum(:)     = 0.
      p%pres(:)     = 0.
#else
      p%Srad(:)     = 0.
      p%Rh(:)       = 0.
      p%Rain(:)     = 0.
      p%Snow(:)     = 0.
      p%Swdr(:)     = 0.
      p%Swdf(:)     = 0.
#endif

      return
   end subroutine

#ifdef COUP_CSM

   subroutine initlnd2atm()
      implicit none

      type(CPL),pointer :: p
      p=>lnd2atm

      allocate(p%rofliq(npoints))
      allocate(p%rofice(npoints))
      allocate(p%t_Rad(npoints))
      allocate(p%tref(npoints))
      allocate(p%qref(npoints))
      allocate(p%avsdr(npoints))
      allocate(p%anidr(npoints))
      allocate(p%avsdf(npoints))
      allocate(p%anidf(npoints))
      allocate(p%snowh(npoints))
      allocate(p%u10(npoints))
      allocate(p%ddvel(npoints))
      allocate(p%fv(npoints))
      allocate(p%ram1(npoints))
      allocate(p%soilw(npoints))
      allocate(p%taux(npoints))
      allocate(p%tauy(npoints))
      allocate(p%LH(npoints))
      allocate(p%SH(npoints))
      allocate(p%lwup(npoints))
      allocate(p%evap(npoints))
      allocate(p%swnet(npoints))
      allocate(p%fco2(npoints))
      allocate(p%flxdst1(npoints))
      allocate(p%flxdst2(npoints))
      allocate(p%flxdst3(npoints))
      allocate(p%flxdst4(npoints))
      allocate(p%flxvoc(npoints))

      p%rofliq(:)     = 0.
      p%rofice(:)     = 0.
      p%t_Rad(:)      = 0.
      p%tref(:)       = 0.
      p%qref(:)       = 0.
      p%avsdr(:)      = 0.
      p%anidr(:)      = 0.
      p%avsdf(:)      = 0.
      p%anidf(:)      = 0.
      p%snowh(:)      = 0.
      p%u10(:)        = 0.
      p%ddvel(:)      = 0.
      p%fv(:)         = 0.
      p%ram1(:)       = 0.
      p%soilw(:)      = 0.
      p%taux(:)       = 0.
      p%tauy(:)       = 0.
      p%LH(:)         = 0.
      p%SH(:)         = 0.
      p%lwup(:)       = 0.
      p%evap(:)       = 0.
      p%swnet(:)      = 0.
      p%fco2(:)       = 0.
      p%flxdst1(:)    = 0.
      p%flxdst2(:)    = 0.
      p%flxdst3(:)    = 0.
      p%flxdst4(:)    = 0.
      p%flxvoc(:)     = 0.

   end subroutine

#endif

   subroutine InitSurf()
      implicit none
      type(surf),pointer::p
      p=>bound

      allocate(p%lcno(npoints,PFT))
      allocate(p%stext(npoints))
      allocate(p%PCT_PFT(npoints,PFT))
      allocate(p%clumping(npoints))
      allocate(p%longitude(npoints))
      allocate(p%latitude(npoints))

      allocate(p%sdp(npoints))
      allocate(p%st(npoints))
      allocate(p%sw(npoints))

! 2025/02/16
      allocate(p%HeightC(npoints))

      allocate(p%laiyr(npoints,PFT))
      allocate(p%nppyr(npoints,PFT))
      allocate(p%ccd(npoints,PFT))
      allocate(p%cfmd(npoints,PFT))
      allocate(p%cfsd(npoints,PFT))
      allocate(p%cm(npoints,PFT))
      allocate(p%cp(npoints,PFT))
      allocate(p%cs(npoints,PFT))
      allocate(p%csm(npoints,PFT))
      allocate(p%csmd(npoints,PFT))
      allocate(p%cssd(npoints,PFT))
      allocate(p%lai(npoints,PFT))
!allocate(p%Vcmax(npoints,PFT))
!*********************** crop ********************
      allocate(p%tt_veg(npoints,PFT))
      allocate(p%tt_rep(npoints,PFT))
      allocate(p%phot_type(npoints,PFT))
      allocate(p%emer_doy(npoints,PFT))
      allocate(p%har_doy(npoints,PFT))
!*********************** crop ********************
! selected parameters
      allocate(p%p_Vcmax(npoints,PFT))
      allocate(p%p_b_h2o(npoints,PFT))
      allocate(p%p_m_h2o(npoints,PFT))
      allocate(p%p_Ksat_scalar(npoints))
      allocate(p%p_b_scalar(npoints))
      allocate(p%p_fei_min(npoints,PFT))
      allocate(p%p_Lr(npoints,PFT))
      allocate(p%p_r_r(npoints,PFT))
      allocate(p%p_ppslh(npoints,PFT))


      p%lcno      = 0.
      p%stext     = 0.
      p%PCT_PFT   = 0.
      p%clumping  = 0.
      p%longitude = 0.
      p%latitude  = 0.

      p%sdp       = 0.
      p%st        = 0.
      p%sw        = 0.

      p%laiyr     = 0.
      p%nppyr     = 0.
      p%ccd       = 0.
      p%cfmd      = 0.
      p%cfsd      = 0.
      p%cm        = 0.
      p%cp        = 0.
      p%cs        = 0.
      p%csm       = 0.
      p%csmd      = 0.
      p%cssd      = 0.
!2025/02/16
      p%HeightC   = 0.

      p%lai       = 0.
!p%Vcmax     = 0.
!*********************** crop ********************
      p%tt_veg(:,:)      = 0.
      p%tt_rep(:,:)     = 0.
      p%phot_type(:,:)   = 0.
      p%emer_doy(:,:)  = 0.
      p%har_doy(:,:) = 0.
!*********************** crop ********************
! selected parameters
      p%p_Vcmax(:,:) = 0.
      p%p_b_h2o(:,:) = 0.
      p%p_m_h2o(:,:) = 0.
      p%p_Ksat_scalar(:) = 0.
      p%p_b_scalar(:) = 0.
      p%p_fei_min(:,:) = 0.
      p%p_Lr(:,:) = 0.
      p%p_r_r(:,:) = 0.
      p%p_ppslh(:,:) = 0.



   end subroutine

   subroutine InitSoilstat()
      implicit none
      type(soils),pointer  :: p
      p => soilstat

      allocate(p%n_layer(npoints))
      allocate(p%Zp(npoints,PFT))
      allocate(p%Zsp(npoints,PFT))
      allocate(p%r_rain_g(npoints,PFT))
      allocate(p%r_drainage(npoints,PFT))
      allocate(p%r_root_decay(npoints,PFT))
      allocate(p%psi_min(npoints,PFT))
      allocate(p%alpha(npoints,PFT))
      allocate(p%f_soilwater(npoints,PFT))
      allocate(p%d_soil(npoints,0:max_layers-1))
      allocate(p%f_root(npoints,0:max_layers-1,PFT))
      allocate(p%dt(npoints,0:max_layers-1,PFT))
      allocate(p%thermal_cond(npoints,0:max_layers-1,PFT))
      allocate(p%theta_vfc(npoints,0:max_layers-1,PFT))
      allocate(p%theta_vwp(npoints,0:max_layers-1,PFT))
      allocate(p%fei(npoints,0:max_layers-1,PFT))
      allocate(p%Ksat(npoints,0:max_layers-1,PFT))
      allocate(p%psi_sat(npoints,0:max_layers-1,PFT))
      allocate(p%b(npoints,0:max_layers-1,PFT))
      allocate(p%density_soil(npoints,0:max_layers-1))
      allocate(p%f_org(npoints,0:max_layers-1,PFT))
      allocate(p%ice_ratio(npoints,0:max_layers-1,PFT))
      allocate(p%thetam(npoints,0:max_layers-1,PFT))
      allocate(p%thetam_prev(npoints,0:max_layers-1,PFT))
      allocate(p%temp_soil_p(npoints,0:max_layers-1,PFT))
      allocate(p%temp_soil_c(npoints,0:max_layers-1,PFT))
      allocate(p%f_ice(npoints,0:max_layers-1,PFT))
      allocate(p%psim(npoints,0:max_layers-1,PFT))
      allocate(p%thetab(npoints,0:max_layers-1,PFT))
      allocate(p%psib(npoints,0:max_layers-1,PFT))
      allocate(p%r_waterflow(npoints,0:max_layers-1,PFT))
      allocate(p%km(npoints,0:max_layers-1,PFT))
      allocate(p%kb(npoints,0:max_layers-1,PFT))
      allocate(p%KK(npoints,0:max_layers-1,PFT))
      allocate(p%Cs(npoints,0:max_layers-1,PFT))
      allocate(p%lambda(npoints,0:max_layers-1,PFT))
      allocate(p%Ett(npoints,0:max_layers-1,PFT))
      allocate(p%G(npoints,0:max_layers-1,PFT))
! 2025/02/17
      allocate(p%f_feileaf(npoints,PFT))
      allocate(p%Sp(npoints,PFT))
      allocate(p%psim_prev(npoints,0:max_layers-1,PFT))

      p%f_feileaf(:,:)       = 0.
      p%Sp(:,:)              = 0.
      p%psim_prev(:,:,:)     = 0.
!...................

      p%n_layer(:)           = 0
      p%Zp(:,:)              = 0.
      p%Zsp(:,:)             = 0.
      p%r_rain_g(:,:)        = 0.
      p%r_drainage(:,:)      = 0.
      p%r_root_decay(:,:)    = 0.
      p%psi_min(:,:)         = 0.
      p%alpha(:,:)           = 0.
      p%f_soilwater(:,:)     = 0.

      p%d_soil(:,:)          = 0.
      p%f_root(:,:,:)        = 0.
      p%dt(:,:,:)            = 0.
      p%thermal_cond(:,:,:)  = 0.
      p%theta_vfc(:,:,:)     = 0.
      p%theta_vwp(:,:,:)     = 0.
      p%fei(:,:,:)           = 0.
      p%Ksat(:,:,:)          = 0.
      p%psi_sat(:,:,:)       = 0.
      p%b(:,:,:)             = 0.
      p%density_soil(:,:)    = 0.
      p%f_org(:,:,:)         = 0.
      p%ice_ratio(:,:,:)     = 0.
      p%thetam(:,:,:)        = 0.
      p%thetam_prev(:,:,:)   = 0.
      p%temp_soil_p(:,:,:)   = 0.
      p%temp_soil_c(:,:,:)   = 0.
      p%f_ice(:,:,:)         = 0.
      p%psim(:,:,:)          = 0.
      p%thetab(:,:,:)        = 0.
      p%psib(:,:,:)           = 0.
      p%r_waterflow(:,:,:)   = 0.
      p%km(:,:,:)            = 0.
      p%kb(:,:,:)            = 0.
      p%KK(:,:,:)            = 0.
      p%Cs(:,:,:)            = 0.
      p%lambda(:,:,:)        = 0.
      p%Ett(:,:,:)           = 0.
      p%G(:,:,:)             = 0.

   end subroutine

   subroutine InitOuput()
      implicit none
      type(res),pointer::p
      p=>output

      allocate(p%GPPpft(npoints,PFT))
      allocate(p%SIFpft(npoints,PFT))
      allocate(p%SIFpft_sat(npoints,PFT))
      allocate(p%NPPpft(npoints,PFT))
      allocate(p%NEPpft(npoints,PFT))
      allocate(p%SHpft(npoints,PFT))
      allocate(p%LHpft(npoints,PFT))
      allocate(p%Transpft(npoints,PFT))
      allocate(p%Evappft(npoints,PFT))
      allocate(p%Net_Radpft(npoints,PFT))
      allocate(p%GPP(npoints))
      allocate(p%SIF(npoints))
      allocate(p%SIF_sat(npoints))
      allocate(p%NPP(npoints))
      allocate(p%NEP(npoints))
      allocate(p%LAIpft(npoints,PFT))
      allocate(p%LAI(npoints))
      allocate(p%SH(npoints))
      allocate(p%LH(npoints))
      allocate(p%Trans(npoints))
      allocate(p%Evap(npoints))
      allocate(p%Net_Rad(npoints))
      allocate(p%Thetampft(npoints,PFT))
      allocate(p%Thetam(npoints))
      allocate(p%fAPARpft(npoints,PFT))
      allocate(p%fAPAR(npoints))
      allocate(p%VODpft(npoints,PFT))
      allocate(p%VOD(npoints))
      allocate(p%COS_fluxpft(npoints,PFT))
      allocate(p%COS_flux(npoints))
!********************************** lm
      allocate(p%Harvcpft(npoints,PFT))
      allocate(p%Harvc(npoints))
!2025/02/18
      allocate(p%fei_leafpft(npoints,PFT))
      allocate(p%fei_leaf(npoints))
      allocate(p%ETapft(npoints,PFT))
      allocate(p%ETa(npoints))
      allocate(p%LHapft(npoints,PFT))
      allocate(p%LHa(npoints))
      allocate(p%Quptpft(npoints,PFT))
      allocate(p%Qupt(npoints))
      allocate(p%Thetam2_pft(npoints,PFT))
      allocate(p%Thetam_layer2(npoints))
      allocate(p%Thetam3_pft(npoints,PFT))
      allocate(p%Thetam_layer3(npoints))
      allocate(p%Thetam4_pft(npoints,PFT))
      allocate(p%Thetam_layer4(npoints))
      allocate(p%Thetam5_pft(npoints,PFT))
      allocate(p%Thetam_layer5(npoints))
! 2025/02/18
      p%fei_leafpft(:,:) = 0.
      p%fei_leaf(:) = 0.
      p%ETapft(:,:) = 0.
      p%ETa(:) = 0.
      p%LHapft(:,:) = 0.
      p%LHa(:) = 0.
      p%Quptpft(:,:) = 0.
      p%Qupt(:) = 0.
      p%Thetam2_pft(:,:) = 0.
      p%Thetam_layer2(:) = 0.
      p%Thetam3_pft(:,:) = 0.
      p%Thetam_layer3(:) = 0.
      p%Thetam4_pft(:,:) = 0.
      p%Thetam_layer4(:) = 0.
      p%Thetam5_pft(:,:) = 0.
      p%Thetam_layer5(:) = 0.

      p%GPPpft(:,:)  = 0.
      p%SIFpft(:,:)  = 0.
      p%SIFpft(:,:)  = 0.
      p%NPPpft(:,:)  = 0.
      p%NEPpft(:,:)  = 0.
      p%SHpft(:,:)   = 0.
      p%LHpft(:,:)   = 0.
      p%Transpft(:,:)= 0.
      p%Evappft(:,:) = 0.
      p%Net_Radpft(:,:) = 0.
      p%LAIpft(:,:)  = 0.
      p%Thetampft(:,:)  = 0.
      p%fAPARpft(:,:)  = 0.
      p%VODpft(:,:)  = 0.
      p%COS_fluxpft(:,:)  = 0.
!****************************** lm
      p%Harvcpft(:,:)   = 0.

      p%GPP(:)       = 0.
      p%SIF(:)       = 0.
      p%SIF_sat(:)   = 0.
      p%NPP(:)       = 0.
      p%NEP(:)       = 0.
      p%LAI(:)       = 0.
      p%SH(:)        = 0.
      p%LH(:)        = 0.
      p%Trans(:)     = 0.
      p%Evap(:)      = 0.
      p%Net_Rad(:)   = 0.
      p%Thetam(:)   = 0.
      p%fAPAR(:)    = 0.
      p%COS_flux(:) = 0.
      p%VOD(:)      = 0.
!****************************** lm
      p%Harvc(:)    = 0.
   end subroutine

end module
