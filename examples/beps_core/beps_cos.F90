!!!Calculation of uptake of carbonyl sulfide by plants and soil.
!!!Equations presented here follow:
!
!         Berry J.A. et al., 2013: A coupled model of the global cycles of carbonyl sulfide
!                 and CO2: A possible new window on the carbon cycle. Journal of
!                 Geophysical Research, Biogeosciences, 118, doi:10.1002/jgrg.20068.
!
!
subroutine cos_calc(cosa,g_sw,g_bw,vmax0,soilp,mid_res)
    use shr_kind_mod, only: r8 =>shr_kind_r8
    use beps_soilMod
    implicit none

    !Input Variables
    type(soil), intent(in)           ::soilp
    type(results), intent(inout)     ::mid_res

    real(r8) :: vmax0, g_sw, g_bw

    !Local Variables
    real(r8) :: cosa   ! CAS COS concentration (mol COS/mol air)
    real(r8) :: gcosm
    real(r8) :: gtcos
    real(r8) :: cos_grnd1,cos_grnd2 ! local ground COS flux (mol/m2/sec)
    real(r8) :: freeze  ! term to restrict COS uptake by frozen ground
    real(r8) :: moist   ! term to restrict COS uptake by saturated ground
    real(r8) :: wfract    ! fraction of saturation in top 3 soil layers
    real(r8) :: F_opt, S_opt, F_g, S_g, a
    real(r8), parameter :: k_cos_soil = 1.2E-4

    !Misc Variables
    integer :: j
    real(r8):: soil_T, soil_S, dsoil, soil_ice
    real(r8):: cos_soil_abiotic, cos_soil_biotic
    intrinsic log

    soil_T    = 0.
    soil_S    = 0.
    soil_ice  = 0.
    dsoil     = 0.

!    cosa = 450. * 1.0e-12  ! unit: mol COS per mol air, needs to set for every month or year
    gcosm = dble(1.40e3) * vmax0 * 1.0e-6
    gtcos = dble(1.0) / ((dble(1.94) / g_sw) + (dble(1.56) / (g_bw/1.6)) &
            + (dble(1.0) / gcosm))
    mid_res%COS_plant = gtcos * cosa

    !...ground uptake of COS
    do j = 1,3
        soil_T = soil_T + soilp%temp_soil_c(j-1)*soilp%d_soil(j-1)
        soil_S = soil_S + soilp%thetam(j-1)*soilp%d_soil(j-1)
        soil_ice = soil_ice + soilp%ice_ratio(j-1)*soilp%d_soil(j-1)
        dsoil = dsoil + soilp%d_soil(j-1)
    end do

    soil_T = soil_T/dsoil
    soil_S = soil_S/dsoil * 100.
    soil_ice = soil_ice/dsoil

    cos_soil_abiotic = 0.437 *exp(0.0984 * soil_T)

    F_opt = -0.00986 * soil_T * soil_T + 0.197 * soil_T - 9.32
    S_opt = 0.28 * soil_T + 14.5
    F_g = -0.0119 * soil_T * soil_T + 0.110 * soil_T -1.18
    S_g  = 35.0
    a = log(F_opt/F_g) * (log(S_opt/S_g) + (S_g/S_opt - 1.))**(-1)

    cos_soil_biotic = F_opt * (soil_S/S_opt)**a * exp(-a * (soil_S/S_opt -  1.))

    cos_grnd1 = cos_soil_abiotic + cos_soil_biotic

    !...calculations limited to top 3 soil layers, but
    !...scaling contributions to total soil respiration
    if(soil_ice < 1.e-8 .and. soil_S < 1.e-8) then
        freeze = 1.0
        else
            freeze = (soil_ice * soil_S/100.)/(soil_ice * soil_S/100. + soil_S/100.)
    end if

    wfract = soil_S/100.
    moist = min(1.0, 1.42 - (wfract * 1.42))

    cos_grnd2 = (mid_res%NPP - mid_res%NEP)/12. * cosa * k_cos_soil * freeze * moist

    mid_res%COS_grnd = max(cos_grnd1, cos_grnd2)

end subroutine cos_calc


