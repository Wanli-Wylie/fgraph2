module mid_results
use shr_kind_mod,only: r8=>shr_kind_r8
implicit none

type :: results
   real(r8)  :: gpp_o_sunlit
   real(r8)  :: gpp_u_sunlit
   real(r8)  :: gpp_o_shaded
   real(r8)  :: gpp_u_shaded
   real(r8)  :: plant_resp
   real(r8)  :: npp_o
   real(r8)  :: npp_u
   real(r8)  :: GPP
   real(r8)  :: SIF
   real(r8)  :: NPP
   real(r8)  :: NEP
   real(r8)  :: soil_resp
   real(r8)  :: Net_Rad
   real(r8)  :: SH
   real(r8)  :: LH
   real(r8)  :: Trans
   real(r8)  :: Evap
   real(r8)  :: thetam_surf
   real(r8)  :: COS_flux
   real(r8)  :: lai
   real(r8)  :: lai_old
   real(r8)  :: lai_new
   real(r8)  :: COS_plant
   real(r8)  :: COS_grnd
   real(r8)  :: fAPAR
    !******************************* lm
   real(r8)  :: Harvc
   real(r8)  :: thetam2
   real(r8)  :: thetam3
   real(r8)  :: thetam4
   real(r8)  :: thetam5
   real(r8)  :: ETa  ! trans from SAPC
   real(r8)  :: LHa  !
   real(r8)  :: VOD
   real(r8)  :: fei_leaf
   real(r8)  :: Qupt ! water uptake from SPAC
   real(r8)  :: res_H
   real(r8)  :: res_A
end type results

type :: climatedata
   real(r8) :: temp      ! temperatre  (oC)
   real(r8)  :: Srad     !solar radiation
   real(r8)  :: LR       !downward longwave radiation
   real(r8)  :: rainfall !liquid water rainfall
   real(r8)  :: snow     !snow
   real(r8)  :: S_dff    !diffuse solar radiation
   real(r8)  :: S_dir    !direct solar radiation
   real(r8)  :: rh
   real(r8)  :: wind
   real(r8)  :: tempmx
   real(r8)  :: tempmn
end type climatedata


!! These types are adapted from DB.h, but many variables are not used @J.Wang
!! Importantly, these three types are only used by photosynthesis module.
type :: meteorology
   real(r8)  :: ustar                 !friction velocity, m s-1
   real(r8)  :: ustarnew              !updated friction velocity with new H, m s-1
   real(r8)  :: rhova_g               !absolute humidity, g m-3
   real(r8)  :: rhova_kg              !absolute humidity, kg m-3
   real(r8)  :: sensible_heat_flux    !sensible heat flux, W M-2
   real(r8)  :: H_old                 !old sensible heat flux, W m-2
   real(r8)  :: air_density           !air density, kg m-3
   real(r8)  :: T_Kelvin              !absolute air temperature, K
   real(r8)  :: press_kpa             !station pressure, kPa
   real(r8)  :: press_bars            !station pressure, bars
   real(r8)  :: press_Pa              !pressure, Pa
   real(r8)  :: pstat273              !gas constant computations
   real(r8)  :: air_density_mole      !air density, mole m-3
   real(r8)  :: relative_humidity     !relative humidity, ea/es(T)
   real(r8)  :: vpd                   !vapor pressure deficit
   real(r8)  :: ir_in                 !infrared flux density
end type meteorology

type  :: factors
   real(r8) :: latent                 !latent heat of vaporization, J kg-1
   real(r8) :: latent18               !latent heat of vaporization times molecular mass of vapor, 18 g mol-1
   real(r8) :: heatcoef               !factor for sensible heat flux density
   real(r8) :: a_filt                 !filter coefficients
   real(r8) :: b_filt                 !filter coefficients
   real(r8) :: co2                    !CO2 factor, ma/mc * rhoa (mole m-3)
end type factors

type :: boundary_layer_resistances
   real(r8) :: vapor                  !resistance for water vapor, s/m
   real(r8) :: heat                   !resistance for heat, s/m
   real(r8) :: co2                    !resistance for CO2, s/m
end type boundary_layer_resistances

end module

