! =================================================================================================================================
! MODULE       : thermosoilc
!
! CONTACT      : orchidee-help _at_ listes.ipsl.fr
!
! LICENCE      : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF        Calculates the soil temperatures by solving the heat
!! diffusion equation within the soil. This module is only used with Choisnel hydrology.
!!
!!\n DESCRIPTION : General important informations about the numerical scheme and
!!                 the soil vertical discretization:\n
!!               - the soil is divided into "ngrnd" (=7 by default) layers, reaching to as
!!                 deep as 5.5m down within the soil, with thiscknesses
!!                 following a geometric series of ration 2.\n
!!               - "jg" is usually used as the index going from 1 to ngrnd to describe the
!!                  layers, from top (jg=1) to bottom (jg=ngrnd)\n
!!               - the thermal numerical scheme is implicit finite differences.\n
!!                 -- When it is resolved in thermosoil_profile at the present timestep t, the
!!                 dependancy from the previous timestep (t-1) is hidden in the
!!                 integration coefficients cgrnd and dgrnd, which are therefore
!!                 calculated at the very end of thermosoil_main (call to
!!                 thermosoil_coef) for use in the next timestep.\n
!!                 -- At timestep t, the system becomes :\n 
!!
!!                              T(k+1)=cgrnd(k)+dgrnd(k)*T(k) \n
!!                                      -- EQ1 -- \n
!!
!!                 (the bottom boundary condition has been used to obtained this equation).\n
!!                 To solve it, the uppermost soil temperature T(1) is required.
!!                 It is obtained from the surface temperature Ts, which is
!!                 considered a linear extrapolation of T(1) and T(2)\n
!!
!!                           Ts=(1-lambda)*T(1) -lambda*T(2) \n 
!!                                      -- EQ2--\n
!!
!!                 -- caveat 1 : Ts is called 'temp_soil_new' in this routine,
!!                 don' t act.\n
!!                 -- caveat 2 : actually, the surface temperature at time t Ts
!!                 depends on the soil temperature at time t through the
!!                 ground heat flux. This is again implicitly solved, with Ts(t)
!!                 expressed as :\n
!!
!!                 soilcap*(Ts(t)-Ts(t-1))/dt=soilflux+otherfluxes(Ts(t))\n 
!!                                      -- EQ3 --\n
!!
!!                 and the dependency from the previous timestep is hidden in
!!                 soilcap and soilflux (apparent surface heat capacity and heat
!!                 flux respectively). Soilcap and soilflux are therefore
!!                 calculated at the previsou timestep, at the very end of thermosoil
!!                 (final call to thermosoil_coef) and stored to be used at the next time step.
!!                 At timestep t, EQ3 is solved for Ts in enerbil, and Ts
!!                 is used in thermosoil to get T(1) and solve EQ1.\n
!!
!! - lambda is the @tex $\mu$ @endtex of F. Hourdin' s PhD thesis, equation (A28); ie the
!! coefficient of the linear extrapolation of Ts (surface temperature) from T1 and T2 (ptn(jg=1) and ptn(jg=2)), so that:\n
!! Ts= (1+lambda)*T(1)-lambda*T(2) --EQ2-- \n
!! lambda = (zz_coeff(1))/((zz_coef(2)-zz_coef(1))) \n
!!
!! - cstgrnd is the attenuation depth of the diurnal temperature signal
!! (period : one_day) as a result of the heat conduction equation
!! with no coefficients :
!!\latexonly
!!\input{thermosoil_var_init0.tex}
!!\endlatexonly
!!  -- EQ4 --\n
!! This equation results from the change of variables :
!! z' =z*sqrt(Cp/K) where z' is the new depth (homogeneous
!! to sqrt(time) ), z the real depth (in m), Cp and K the soil heat
!! capacity and conductivity respectively.\n
!!
!! the attenuation depth of a diurnal thermal signal for EQ4 is therefore homogeneous to sqrt(time) and
!! equals : \n
!! cstgrnd = sqrt(oneday/Pi)
!!
!! - lskin is the attenuation depth of the diurnal temperature signal
!! (period : one_day) within the soil for the complete heat conduction equation
!! (ie : with coefficients)
!!\latexonly
!!\input{thermosoil_var_init00.tex}
!!\endlatexonly
!! -- EQ5 --  \n
!! it can be retrieved from cstgrnd using the change of variable  z' =z*sqrt(Cp/K):\n
!! lskin = sqrt(K/Cp)*cstgrnd =  sqrt(K/Cp)*sqrt(oneday//Pi)\n
!! 
!! In thermosoil, the ratio lskin/cstgrnd is frequently used as the
!! multiplicative factor to go from
!!'adimensional' depths (like z' ) to real depths (z). z' is not really
!! adimensional but is reffered to like this in the code.
!!
!!
!! RECENT CHANGE(S) : thermosoilc module is a copy of thermosoil before the vertical discretization was changed.
!!                    This module is used only for the Choisnel hydrology. 
!!
!! REFERENCE(S) : None
!!
!! SVN          :
!! $HeadURL: svn://forge.ipsl.jussieu.fr/orchidee/branches/ORCHIDEE-CN-P/ORCHIDEE/src_sechiba/thermosoilc.f90 $
!! $Date: 2017-10-26 14:32:36 +0200 (四, 2017-10-26) $
!! $Revision: 4717 $
!! \n
!_ ================================================================================================================================

MODULE thermosoilc

  ! modules used :
  USE ioipsl
  USE ioipsl_para
  USE xios_orchidee
  USE constantes
  USE time, ONLY : one_day, dt_sechiba
  USE constantes_soil
  USE sechiba_io_p
  USE grid


  IMPLICIT NONE

  !private and public routines :
  PRIVATE
  PUBLIC :: thermosoilc_main, thermosoilc_clear, thermosoilc_levels, thermosoilc_initialize, thermosoilc_finalize

  REAL, SAVE                               :: lambda, cstgrnd, lskin   !! See Module description
!$OMP THREADPRIVATE(lambda, cstgrnd, lskin)

  REAL, SAVE                               :: fz1, zalph               !! usefull constants for diverse use
!$OMP THREADPRIVATE(fz1, zalph)
  REAL, ALLOCATABLE, SAVE, DIMENSION (:,:) :: ptn                      !! vertically discretized 
                                                                              !! soil temperatures @tex ($K$) @endtex. 
!$OMP THREADPRIVATE(ptn)
  REAL, ALLOCATABLE, SAVE, DIMENSION (:)   :: zz                       !! depths of the soil thermal numerical nodes. 
                                                                              !! Caveats: they are not exactly the centers of the
                                                                              !! thermal layers, see the calculation in 
                                                                              !! ::thermosoilc_var_init  @tex ($m$) @endtex.
!$OMP THREADPRIVATE(zz)
  REAL, ALLOCATABLE, SAVE, DIMENSION (:)   :: zz_coef		      !! depths of the boundaries of the thermal layers,
                                                                              !! see the calculation in 
                                                                              !! thermosoilc_var_init  @tex ($m$) @endtex.
!$OMP THREADPRIVATE(zz_coef)
  REAL, ALLOCATABLE, SAVE, DIMENSION (:)   :: dz1                      !! numerical constant used in the thermal numerical
                                                                              !! scheme  @tex ($m^{-1}$) @endtex. ; it corresponds
                                                                              !! to the coefficient  @tex $d_k$ @endtex of equation
                                                                              !! (A.12) in F. Hourdin PhD thesis.
!$OMP THREADPRIVATE(dz1)
  REAL, ALLOCATABLE, SAVE, DIMENSION (:)   :: dz2                      !! thicknesses of the thermal layers  @tex ($m$)
                                                                              !! @endtex; typically: 
                                                                              !! dz2(jg)=zz_coef(jg+1)-zz_coef(jg); calculated once 
                                                                              !! and for all in thermosoilc_var_init
!$OMP THREADPRIVATE(dz2)
  REAL, ALLOCATABLE, SAVE, DIMENSION (:,:) :: cgrnd                    !! integration coefficient for the numerical scheme,
                                                                              !! see eq.1
!$OMP THREADPRIVATE(cgrnd)
  REAL, ALLOCATABLE, SAVE, DIMENSION (:,:) :: dgrnd                    !! integration coefficient for the numerical scheme,
                                                                              !! see eq.1
!$OMP THREADPRIVATE(dgrnd)

  REAL, ALLOCATABLE, SAVE, DIMENSION (:,:) :: pcapa                    !! volumetric vertically discretized soil heat 
                                                                              !! capacity  @tex ($J K^{-1} m^{-3}$) @endtex. 
                                                                              !! It depends on the soil
                                                                              !! moisture content (shum_ngrnd_perma) and is calculated at 
                                                                              !! each time step in thermosoilc_coef.
!$OMP THREADPRIVATE(pcapa)
  REAL, ALLOCATABLE, SAVE, DIMENSION (:,:) :: pkappa                   !! vertically discretized soil thermal conductivity 
                                                                              !!  @tex ($W K^{-1} m^{-1}$) @endtex. Same as pcapa.
!$OMP THREADPRIVATE(pkappa)
  REAL, ALLOCATABLE, SAVE, DIMENSION (:,:) :: pcapa_snow               !! volumetric vertically discretized snow heat 
                                                                              !! capacity @tex ($J K^{-1} m^{-3}$) @endtex.
!$OMP THREADPRIVATE(pcapa_snow)
  REAL, ALLOCATABLE, SAVE, DIMENSION (:,:) :: pkappa_snow              !! vertically discretized snow thermal conductivity 
                                                                              !! @tex ($W K^{-1} m^{-1}$) @endtex.
!$OMP THREADPRIVATE(pkappa_snow)

  REAL, ALLOCATABLE, SAVE, DIMENSION (:,:) :: pcapa_en                 !! heat capacity used for surfheat_incr and 
                                                                              !! coldcont_incr 
!$OMP THREADPRIVATE(pcapa_en)
  REAL, ALLOCATABLE, SAVE, DIMENSION (:,:) :: ptn_beg                  !! vertically discretized temperature at the 
                                                                              !! beginning of the time step  @tex ($K$) @endtex; 
                                                                              !! is used in 
                                                                              !! thermosoilc_energy for energy-related diagnostic of
                                                                              !! the routine.
!$OMP THREADPRIVATE(ptn_beg)
  REAL, ALLOCATABLE, SAVE, DIMENSION (:)   :: temp_sol_beg             !! Surface temperature at the beginning of the 
                                                                              !! timestep  @tex ($K$) @endtex
!$OMP THREADPRIVATE(temp_sol_beg)
  REAL, ALLOCATABLE, SAVE, DIMENSION (:)   :: surfheat_incr            !! Change in soil heat content during the timestep 
                                                                              !!  @tex ($J$) @endtex.
!$OMP THREADPRIVATE(surfheat_incr)
  REAL, ALLOCATABLE, SAVE, DIMENSION (:)   :: coldcont_incr            !! Change in snow heat content  @tex ($J$) @endtex.
!$OMP THREADPRIVATE(coldcont_incr)
  REAL, ALLOCATABLE, SAVE, DIMENSION (:,:) :: shum_ngrnd_perma         !! Saturation degree on the thermal axes (0-1, dimensionless)
!$OMP THREADPRIVATE(shum_ngrnd_perma)

  !  Variables related to soil freezing
  REAL, ALLOCATABLE, SAVE, DIMENSION (:,:) :: profil_froz              !! Frozen fraction of the soil on hydrological levels (-)
!$OMP THREADPRIVATE(profil_froz)
  REAL,ALLOCATABLE, SAVE, DIMENSION (:) :: e_soil_lat                  !! Accumulated latent heat for the whole soil (J)
!$OMP THREADPRIVATE(e_soil_lat)
  REAL,ALLOCATABLE, SAVE, DIMENSION (:,:) :: pcappa_supp               !! Additional heat capacity due to soil freezing for each soil layer (J/K)
!$OMP THREADPRIVATE(pcappa_supp)    


CONTAINS
  !!  =============================================================================================================================
  !! SUBROUTINE		 		    : thermosoilc_initialize
  !!
  !>\BRIEF			            Allocate module variables, read from restart file or initialize with default values
  !!
  !! DESCRIPTION			    : Allocate module variables, read from restart file or initialize with default values.
  !!                                          Call thermosoilc_var_init to calculate physical constants. 
  !!                                          Call thermosoilc_coef to calculate thermal soil properties.
  !!
  !! RECENT CHANGE(S)			    : None
  !!
  !! REFERENCE(S)			    : None
  !! 
  !! FLOWCHART                              : None
  !! \n
  !_ ==============================================================================================================================
  SUBROUTINE thermosoilc_initialize (kjit,          kjpindex,   rest_id,                   &
                                    temp_sol_new,  snow,       shumdiag_perma,            &
                                    soilcap,       soilflx,    stempdiag,                 &
                                    gtemp, &
                                    frac_snow_veg,frac_snow_nobio,totfrac_nobio, &
                                    snowdz, snowrho, snowtemp, lambda_snow, cgrnd_snow, dgrnd_snow, pb)

    !! 0. Variable and parameter declaration
    !! 0.1 Input variables
    INTEGER, INTENT(in)                            :: kjit             !! Time step number (unitless) 
    INTEGER, INTENT(in)                            :: kjpindex         !! Domain size (unitless)
    INTEGER,INTENT (in)                            :: rest_id          !! Restart file identifier (unitless)
    REAL,DIMENSION (kjpindex), INTENT (in)         :: temp_sol_new     !! Surface temperature at the present time-step,
    REAL,DIMENSION (kjpindex), INTENT (in)         :: snow             !! Snow mass (kg)
    REAL,DIMENSION (kjpindex,nslm), INTENT (in)    :: shumdiag_perma   !! Soil saturation degree on the diagnostic axis (0-1, unitless)
    REAL,DIMENSION (kjpindex), INTENT(in)          :: frac_snow_veg    !! Snow cover fraction on vegeted area
    REAL,DIMENSION (kjpindex,nnobio), INTENT(in)   :: frac_snow_nobio  !! Snow cover fraction on non-vegeted area
    REAL,DIMENSION (kjpindex),INTENT(in)           :: totfrac_nobio    !! Total fraction of continental ice+lakes+cities+...
                                                                              !! (unitless,0-1)
    REAL,DIMENSION (kjpindex,nsnow),INTENT (in)    :: snowdz           !! Snow depth
    REAL, DIMENSION (kjpindex,nsnow), INTENT(in)    :: snowrho         !! Snow density
    REAL, DIMENSION (kjpindex,nsnow), INTENT(in) :: snowtemp        !! Snow temperature
    REAL, DIMENSION (kjpindex), INTENT (in)         :: pb              !! Surface presure (hPa)


    !! 0.2 Output variables
    REAL,DIMENSION (kjpindex), INTENT (out)        :: soilcap          !! apparent surface heat capacity (J m-2 K-1)
    REAL,DIMENSION (kjpindex), INTENT (out)        :: soilflx          !! apparent soil heat flux (W m-2)
    REAL,DIMENSION (kjpindex,nslm), INTENT (out)   :: stempdiag        !! diagnostic temperature profile (K)
    REAL,DIMENSION (kjpindex),INTENT(out)          :: gtemp            !! First soil layer temperature

    !! 0.3 Modified variables
    REAL, DIMENSION (kjpindex), INTENT(inout)       :: lambda_snow     !! Coefficient of the linear extrapolation of surface temperature 
    REAL, DIMENSION (kjpindex,nsnow), INTENT (inout):: cgrnd_snow      !! Integration coefficient for snow numerical scheme
    REAL, DIMENSION (kjpindex,nsnow), INTENT (inout):: dgrnd_snow      !! Integration coefficient for snow numerical scheme

    !! 0.4 Local variables
    INTEGER                                        :: ier, i
    LOGICAL                                               :: calculate_coef   !! Local flag to initialize variables by call to thermosoilc_coef   
!_ ================================================================================================================================
    
    IF (printlev >= 3) WRITE (numout,*) 'Start thermosoilc_initialize '

    !! 1. Allocate soil temperatures variables
    ALLOCATE (ptn(kjpindex,ngrnd),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'thermosoilc_initialize', 'Error in allocation of ptn','','')

    ALLOCATE (zz(ngrnd),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'thermosoilc_initialize', 'Error in allocation of zz','','')

    ALLOCATE (zz_coef(ngrnd),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'thermosoilc_initialize', 'Error in allocation of zz_coef','','')

    ALLOCATE (dz1(ngrnd),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'thermosoilc_initialize', 'Error in allocation of dz1','','')

    ALLOCATE (dz2(ngrnd),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'thermosoilc_initialize', 'Error in allocation of dz2','','')

    ALLOCATE (cgrnd(kjpindex,ngrnd-1),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'thermosoilc_initialize', 'Error in allocation of cgrnd','','')

    ALLOCATE (dgrnd(kjpindex,ngrnd-1),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'thermosoilc_initialize', 'Error in allocation of dgrnd','','')

    ALLOCATE (pkappa_snow(kjpindex,nsnow),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'thermosoilc_initialize', 'Error in allocation of pkappa_snow','','')

    ALLOCATE (pcapa_snow(kjpindex,nsnow),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'thermosoilc_initialize', 'Error in allocation of pcapa_snow','','')

    ALLOCATE (pcapa(kjpindex,ngrnd),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'thermosoilc_initialize', 'Error in allocation of pcapa','','')

    ALLOCATE (pkappa(kjpindex,ngrnd),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'thermosoilc_initialize', 'Error in allocation of pkappa','','')

    ! Temporary fix: Initialize following variable because they are output to xios before the first calculation
    pcapa  = 0
    pkappa = 0
    pcapa_snow  = 0
    pkappa_snow = 0

    ALLOCATE (surfheat_incr(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'thermosoilc_initialize', 'Error in allocation of surfheat_incr','','')

    ALLOCATE (coldcont_incr(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'thermosoilc_initialize', 'Error in allocation of coldcont_incr','','')

    ALLOCATE (pcapa_en(kjpindex,ngrnd),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'thermosoilc_initialize', 'Error in allocation of pcapa_en','','')

    ALLOCATE (ptn_beg(kjpindex,ngrnd),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'thermosoilc_initialize', 'Error in allocation of ptn_beg','','')

    ALLOCATE (temp_sol_beg(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'thermosoilc_initialize', 'Error in allocation of temp_sol_beg','','')

    ALLOCATE (shum_ngrnd_perma(kjpindex,ngrnd),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'thermosoilc_initialize', 'Error in allocation of shum_ngrnd_perma','','')

    ALLOCATE (profil_froz(kjpindex,ngrnd),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'thermosoilc_initialize', 'Error in allocation of profil_froz','','')

    IF (ok_freeze_thermix) THEN
       ALLOCATE (pcappa_supp(kjpindex,ngrnd),stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'thermosoilc_initialize', 'Error in allocation of ok_freeze_termix','','')
    END IF
    IF (ok_Ecorr) THEN
       ALLOCATE (e_soil_lat(kjpindex),stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'thermosoilc_initialize', 'Error in allocation of e_soil_lat','','')
    END IF

    !! 2. Initialize variable from restart file or with default values 
    
    !! Reads restart files for soil temperatures only. If no restart file is
    !! found,  the initial soil temperature is by default set to 280K at all depths. The user
    !! can decide to initialize soil temperatures at an other value, in which case he should set the flag THERMOSOIL_TPRO
    !! to this specific value in the run.def.
    IF (printlev>=3) WRITE (numout,*) 'Read restart file for THERMOSOIL variables'

    CALL ioconf_setatt_p('UNITS', 'K')
    CALL ioconf_setatt_p('LONG_NAME','Soil Temperature profile')
    CALL restget_p (rest_id, 'ptn', nbp_glo, ngrnd, 1, kjit, .TRUE., ptn, "gather", nbp_glo, index_g)

    ! Initialize ptn if it was not found in restart file
    IF (ALL(ptn(:,:)==val_exp)) THEN 
       ! ptn was not found in restart file

       IF (read_reftemp) THEN
          ! Read variable ptn from file
          CALL thermosoilc_read_reftempfile(kjpindex,lalo,ptn)
       ELSE
          ! Initialize ptn with a constant value which can be set in run.def
          
          !Config Key   = THERMOSOIL_TPRO
          !Config Desc  = Initial soil temperature profile if not found in restart
          !Config Def   = 280.
          !Config If    = OK_SECHIBA
          !Config Help  = The initial value of the temperature profile in the soil if 
          !Config         its value is not found in the restart file. Here
          !Config         we only require one value as we will assume a constant 
          !Config         throughout the column.
          !Config Units = Kelvin [K]
          CALL setvar_p (ptn, val_exp,'THERMOSOIL_TPRO',280.)
       END IF
    END IF
    
    ! Initialize ptn_beg (variable needed in thermosoilc_coef before calucation in thermosoilc_energy)
    ptn_beg(:,:) = ptn(:,:)
    
    ! Initialize temp_sol_beg with values from previous time-step
    temp_sol_beg(:) = temp_sol_new(:) 
    
    ! Read e_soil_lat from restart file or initialize
    IF (ok_Ecorr) THEN
       CALL restget_p (rest_id, 'e_soil_lat', nbp_glo, 1, 1, kjit, .TRUE., &
            e_soil_lat, "gather", nbp_glo, index_g)
       CALL setvar_p (e_soil_lat, val_exp,'NO_KEYWORD',zero)
    END IF

    ! Read gtemp from restart file
    CALL restget_p (rest_id, 'gtemp', nbp_glo, 1, 1, kjit, .TRUE., &
         gtemp, "gather", nbp_glo, index_g)
    CALL setvar_p (gtemp, val_exp,'NO_KEYWORD',zero)
    
    ! Read variables calculated in thermosoilc_coef from restart file
    ! If the variables were not found in the restart file, the logical 
    ! calculate_coef will be true and thermosoilc_coef will be called further below.
    ! These variables need to be in the restart file to avoid a time shift that
    ! would be done using thermosoilc_coef at this stage.
    calculate_coef=.FALSE.
    CALL ioconf_setatt_p('UNITS', 'J m-2 K-1')
    CALL ioconf_setatt_p('LONG_NAME','Apparent surface heat capacity')
    CALL restget_p (rest_id, 'soilcap', nbp_glo, 1, 1, kjit, .TRUE., &
         soilcap, "gather", nbp_glo, index_g)
    IF (ALL(soilcap(:)==val_exp)) calculate_coef=.TRUE.

    CALL ioconf_setatt_p('UNITS', 'W m-2')
    CALL ioconf_setatt_p('LONG_NAME','Apparent soil heat flux')
    CALL restget_p (rest_id, 'soilflx', nbp_glo, 1, 1, kjit, .TRUE., &
         soilflx, "gather", nbp_glo, index_g)
    IF (ALL(soilflx(:)==val_exp)) calculate_coef=.TRUE.

    CALL ioconf_setatt_p('UNITS', '')
    CALL ioconf_setatt_p('LONG_NAME','Integration coefficient for the numerical scheme')
    CALL restget_p (rest_id, 'cgrnd', nbp_glo, ngrnd-1, 1, kjit, .TRUE., &
         cgrnd, "gather", nbp_glo, index_g)
    IF (ALL(cgrnd(:,:)==val_exp)) calculate_coef=.TRUE.

    CALL ioconf_setatt_p('UNITS', '')
    CALL ioconf_setatt_p('LONG_NAME','Integration coefficient for the numerical scheme')
    CALL restget_p (rest_id, 'dgrnd', nbp_glo, ngrnd-1, 1, kjit, .TRUE., &
         dgrnd, "gather", nbp_glo, index_g)
    IF (ALL(dgrnd(:,:)==val_exp)) calculate_coef=.TRUE.

    CALL ioconf_setatt_p('UNITS', '')
    CALL ioconf_setatt_p('LONG_NAME','Integration coefficient for the numerical scheme')
    CALL restget_p (rest_id, 'cgrnd_snow', nbp_glo, nsnow, 1, kjit, .TRUE., &
         cgrnd_snow, "gather", nbp_glo, index_g)
    IF (ALL(cgrnd_snow(:,:)==val_exp)) calculate_coef=.TRUE.

    CALL ioconf_setatt_p('UNITS', '')
    CALL ioconf_setatt_p('LONG_NAME','Integration coefficient for the numerical scheme')
    CALL restget_p (rest_id, 'dgrnd_snow', nbp_glo, nsnow, 1, kjit, .TRUE., &
         dgrnd_snow, "gather", nbp_glo, index_g)
    IF (ALL(dgrnd_snow(:,:)==val_exp)) calculate_coef=.TRUE.

    CALL ioconf_setatt_p('UNITS', '')
    CALL ioconf_setatt_p('LONG_NAME','Coefficient of the linear extrapolation of surface temperature')
    CALL restget_p (rest_id, 'lambda_snow', nbp_glo, 1, 1, kjit, .TRUE., &
         lambda_snow, "gather", nbp_glo, index_g)
    IF (ALL(lambda_snow(:)==val_exp)) calculate_coef=.TRUE.


    !! 2.2.Computes physical constants and arrays; initializes soil thermal properties; produces the first stempdiag
    !!  Computes some physical constants and arrays depending on the soil vertical discretization 
    !! (lskin, cstgrnd, zz, zz_coef, dz1, dz2); get the vertical humidity onto the thermal levels
    CALL thermosoilc_var_init (kjpindex, shumdiag_perma, snow, &
         zz, zz_coef, dz1, dz2, stempdiag)
    
    !! 2.3. Computes cgrd, dgrd, soilflx and soilcap coefficients from restart values or initialisation values.
    !!      This is done only if they were not found in restart file.
    IF (calculate_coef) THEN
       IF (printlev>=3) WRITE (numout,*) 'thermosoilc_coef will be called in the intialization phase'
       CALL thermosoilc_coef (&
            kjpindex, temp_sol_new, snow, &
            ptn,      soilcap, soilflx,      zz, &
            dz1,      dz2,           &
            cgrnd,    dgrnd, &
            frac_snow_veg,frac_snow_nobio,totfrac_nobio, &
            snowdz,snowrho,  snowtemp,     pb,         &
	    lambda_snow,    cgrnd_snow,    dgrnd_snow)
    END IF

  END SUBROUTINE thermosoilc_initialize


!! ================================================================================================================================
!! SUBROUTINE   : thermosoilc_main
!!
!>\BRIEF        Thermosoilc_main computes the soil thermal properties and dynamics, ie solves
!! the heat diffusion equation within the soil. The soil temperature profile is
!! then interpolated onto the diagnostic axis.
!!
!! DESCRIPTION : The resolution of the soil heat diffusion equation 
!! relies on a numerical finite-difference implicit scheme
!! fully described in the reference and in the header of the thermosoilc module.
!! - The dependency of the previous timestep hidden in the 
!! integration coefficients cgrnd and dgrnd (EQ1), calculated in thermosoilc_coef, and 
!! called at the end of the routine to prepare for the next timestep.
!! - The effective computation of the new soil temperatures is performed in thermosoilc_profile. 
!!
!! - thermosoilc_coef calculates the coefficients for the numerical scheme for the very first iteration of thermosoilc;
!! after that, thermosoilc_coef is called only at the end of the module to calculate the coefficients for the next timestep.
!! - thermosoilc_profile solves the numerical scheme.\n
!!
!! - Flags : one unique flag : THERMOSOIL_TPRO (to be set to the desired initial soil in-depth temperature in K; by default 280K)
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S): vertically discretized soil temperatures ptn, soil
!! thermal properties (pcapa, pkappa), apparent surface heat capacity (soilcap)
!! and heat flux (soilflux) to be used in enerbil at the next timestep to solve
!! the surface energy balance.
!!
!! REFERENCE(S) : 
!! - Hourdin, F. (1992). Study and numerical simulation of the general circulation of planetary atmospheres,
!!  Ph.D. thesis, Paris VII University. Remark: the part of F. Hourdin' s PhD thesis relative to the thermal
!!  integration scheme has been scanned and is provided along with the documentation, with name : 
!!  Hourdin_1992_PhD_thermal_scheme.pdf
!!
!! FLOWCHART    : 
!! \latexonly
!! \includegraphics[scale = 1]{thermosoil_flowchart.png}
!! \endlatexonly
!! 
!! \n
!_ ================================================================================================================================

  SUBROUTINE thermosoilc_main (kjit, kjpindex, &
       index, indexgrnd, &
       indexnslm, &
       temp_sol_new, snow, soilcap, soilflx, &
       shumdiag_perma, stempdiag, ptnlev1, rest_id, hist_id, hist2_id, &
       snowdz,snowrho,snowtemp,gtemp,pb,&
       frac_snow_veg,frac_snow_nobio,totfrac_nobio,temp_sol_add, &
       lambda_snow, cgrnd_snow, dgrnd_snow)

    !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER, INTENT(in)                            :: kjit             !! Time step number (unitless) 
    INTEGER, INTENT(in)                            :: kjpindex         !! Domain size (unitless)
    INTEGER,INTENT (in)                            :: rest_id,hist_id  !! Restart_ file and history file identifier 
                                                                              !! (unitless)
    INTEGER,INTENT (in)                            :: hist2_id         !! history file 2 identifier (unitless)
    INTEGER,DIMENSION (kjpindex), INTENT (in)      :: index            !! Indeces of the points on the map (unitless)
    INTEGER,DIMENSION (kjpindex*ngrnd), INTENT (in):: indexgrnd        !! Indeces of the points on the 3D map (vertical 
                                                                              !! dimension towards the ground) (unitless)
    REAL,DIMENSION (kjpindex), INTENT (inout)      :: temp_sol_new     !! Surface temperature at the present time-step,
                                                                              !! temp_sol_new is only modified for the case ok_explicitsnow
                                                                              !! Ts @tex ($K$) @endtex
    REAL,DIMENSION (kjpindex), INTENT (in)         :: snow             !! Snow mass @tex ($kg$) @endtex.
                                                                              !! Caveat: when there is snow on the
                                                                              !! ground, the snow is integrated into the soil for
                                                                              !! the calculation of the thermal dynamics. It means
                                                                              !! that the uppermost soil layers can completely or 
                                                                              !! partially consist in snow. In the second case, zx1
                                                                              !! and zx2 are the fraction of the soil layer 
                                                                              !! consisting in snow and 'normal' soil, respectively
                                                                              !! This is calculated in thermosoilc_coef.
    REAL,DIMENSION (kjpindex,nslm), INTENT (in)    :: shumdiag_perma   !! Soil saturation degree on the diagnostic axis (0-1, unitless)
    INTEGER,DIMENSION (kjpindex*nslm), INTENT (in) :: indexnslm        !! Indeces of the points on the 3D map
    REAL,DIMENSION (kjpindex,nsnow),INTENT (in)    :: snowdz           !! Snow depth
    REAL,DIMENSION (kjpindex,nsnow),INTENT (in)    :: snowrho          !! Snow density
    REAL, DIMENSION (kjpindex), INTENT (in)        :: pb               !! Surface presure (hPa)

    REAL,DIMENSION (kjpindex,nsnow),INTENT (inout) :: snowtemp         !! Snow temperature
    REAL,DIMENSION (kjpindex), INTENT(in)          :: frac_snow_veg    !! Snow cover fraction on vegeted area
    REAL,DIMENSION (kjpindex,nnobio), INTENT(in)   :: frac_snow_nobio  !! Snow cover fraction on non-vegeted area
    REAL,DIMENSION (kjpindex),INTENT(in)           :: totfrac_nobio    !! Total fraction of continental ice+lakes+cities+...
                                                                              !!(unitless,0-1)
    REAL,DIMENSION (kjpindex), INTENT (inout)      :: temp_sol_add     !! additional surface temperature due to the melt of first layer
                                                                              !!at the present time-step @tex ($K$) @endtex

    !! 0.2 Output variables

    REAL,DIMENSION (kjpindex), INTENT (out)        :: ptnlev1          !! 1st level soil temperature   
    REAL,DIMENSION (kjpindex),INTENT(out)          :: gtemp            !! First soil layer temperature


    !! 0.3 Modified variables

    REAL,DIMENSION (kjpindex), INTENT (inout)      :: soilcap          !! apparent surface heat capacity
                                                                              !! @tex ($J m^{-2} K^{-1}$) @endtex
    REAL,DIMENSION (kjpindex), INTENT (inout)      :: soilflx          !! apparent soil heat flux @tex ($W m^{-2}$) @endtex
                                                                              !! , positive 
                                                                              !! towards the soil, writen as Qg (ground heat flux) 
                                                                              !! in the history files, and computed at the end of 
                                                                              !! thermosoilc for the calculation of Ts in enerbil, 
                                                                              !! see EQ3.
    REAL,DIMENSION (kjpindex,nslm), INTENT (inout) :: stempdiag        !! diagnostic temperature profile @tex ($K$) @endtex
                                                                              !! , eg on the 
                                                                              !! diagnostic axis (levels:1:nslm). The soil 
                                                                              !! temperature is put on this diagnostic axis to be
                                                                              !! used by other modules (slowproc.f90; routing.f90;
                                                                              !! hydrol or hydrolc when a frozen soil 
                                                                              !! parametrization is used..)
    REAL, DIMENSION (kjpindex), INTENT(inout)       :: lambda_snow     !! Coefficient of the linear extrapolation of surface temperature 
    REAL, DIMENSION (kjpindex,nsnow), INTENT (inout):: cgrnd_snow      !! Integration coefficient for snow numerical scheme
    REAL, DIMENSION (kjpindex,nsnow), INTENT (inout):: dgrnd_snow      !! Integration coefficient for snow numerical scheme

    !! 0.4 Local variables

    INTEGER                                        :: jv,ji,ii

!_ ================================================================================================================================
    
  !! 3. Put the soil wetness diagnostic on the levels of the soil temperature

    !!?? this could logically be put just before the last call to
    !!thermosoil_coef, as the results are used there...
    CALL thermosoilc_humlev(kjpindex, shumdiag_perma, snow)

    
  !! 4. Effective computation of the soil temperatures profile, using the cgrd and dgrd coefficients from previsou tstep.
    
    CALL thermosoilc_profile (kjpindex, temp_sol_new, ptn,  &
         stempdiag,snowtemp, frac_snow_veg, frac_snow_nobio,totfrac_nobio, &
         cgrnd_snow, dgrnd_snow)

  !! 5. Call to thermosoilc_energy, still to be clarified..

    CALL thermosoilc_energy (kjpindex, temp_sol_new, soilcap)

  !! 6. Writing the history files according to the ALMA standards (or not..)

    !in only one file (hist2_id <=0) or in 2 different files (hist2_id >0).
    CALL xios_orchidee_send_field("ptn",ptn)
    CALL xios_orchidee_send_field("soilflx",soilflx)
    CALL xios_orchidee_send_field("surfheat_incr",surfheat_incr)
    CALL xios_orchidee_send_field("coldcont_incr",coldcont_incr)
    CALL xios_orchidee_send_field("pkappa",pkappa)
    CALL xios_orchidee_send_field("pkappa_snow",pkappa_snow)
    CALL xios_orchidee_send_field("pcapa",pcapa)
    CALL xios_orchidee_send_field("pcapa_snow",pcapa_snow)
    CALL xios_orchidee_send_field("snowtemp",snowtemp)

    IF ( .NOT. almaoutput ) THEN
      CALL histwrite_p(hist_id, 'ptn', kjit, ptn, kjpindex*ngrnd, indexgrnd)
      CALL histwrite_p(hist_id, 'Qg', kjit, soilflx, kjpindex, index)
    ELSE
      CALL histwrite_p(hist_id, 'SoilTemp', kjit, ptn, kjpindex*ngrnd, indexgrnd)
      CALL histwrite_p(hist_id, 'Qg', kjit, soilflx, kjpindex, index)
      CALL histwrite_p(hist_id, 'DelSurfHeat', kjit, surfheat_incr, kjpindex, index)
      CALL histwrite_p(hist_id, 'DelColdCont', kjit, coldcont_incr, kjpindex, index)
    ENDIF
    IF ( hist2_id > 0 ) THEN
       IF ( .NOT. almaoutput ) THEN
          CALL histwrite_p(hist2_id, 'ptn', kjit, ptn, kjpindex*ngrnd, indexgrnd)
       ELSE
          CALL histwrite_p(hist2_id, 'SoilTemp', kjit, ptn, kjpindex*ngrnd, indexgrnd)
          CALL histwrite_p(hist2_id, 'Qg', kjit, soilflx, kjpindex, index)
          CALL histwrite_p(hist2_id, 'DelSurfHeat', kjit, surfheat_incr, kjpindex, index)
          CALL histwrite_p(hist2_id, 'DelColdCont', kjit, coldcont_incr, kjpindex, index)
       ENDIF
    ENDIF
    
  !! 7. A last final call to thermosoilc_coef
 
    !! A last final call to thermosoilc_coef, which calculates the different
    !!coefficients (cgrnd, dgrnd, dz1, soilcap, soilflx) from this time step to be
    !!used at the next time step, either in the surface temperature calculation
    !!(soilcap, soilflx) or in the soil thermal numerical scheme.
    CALL thermosoilc_coef (&
         kjpindex, temp_sol_new, snow, &
         ptn,      soilcap, soilflx,      zz,   &
         dz1,      dz2,                     &
         cgrnd,   dgrnd, &
         frac_snow_veg,frac_snow_nobio,totfrac_nobio,&
         snowdz,snowrho,  snowtemp,     pb, &
	 lambda_snow,    cgrnd_snow,    dgrnd_snow)

    ! Save variables for explicit snow model
    gtemp(:) = ptn(:,1)

    !! Initialize output arguments to be used in sechiba
    ptnlev1(:) = ptn(:,1)


    !! Surface temperature is forced to zero celcius if its value is larger than melting point, only for explicit snow scheme
    IF  (ok_explicitsnow) THEN
       DO ji=1,kjpindex
          IF  (SUM(snowdz(ji,:)) .GT. 0.0) THEN
             IF (temp_sol_new(ji) .GE. tp_00) THEN
                temp_sol_new(ji) = tp_00
             ENDIF
          END IF
       END DO
    ENDIF


    IF (printlev>=3) WRITE (numout,*) ' thermosoilc_main done '

  END SUBROUTINE thermosoilc_main

  !!  =============================================================================================================================
  !! SUBROUTINE		 		    : thermosoilc_finalize
  !!
  !>\BRIEF                                    Write to restart file
  !!
  !! DESCRIPTION			    : This subroutine writes the module variables and variables calculated in thermosoilc
  !!                                          to restart file
  !! \n
  !_ ==============================================================================================================================
  SUBROUTINE thermosoilc_finalize (kjit,    kjpindex, rest_id,   gtemp,     &
                                  soilcap, soilflx, lambda_snow, cgrnd_snow, dgrnd_snow)

    !! 0. Variable and parameter declaration
    !! 0.1 Input variables
    INTEGER, INTENT(in)                            :: kjit             !! Time step number (unitless) 
    INTEGER, INTENT(in)                            :: kjpindex         !! Domain size (unitless)
    INTEGER,INTENT (in)                            :: rest_id          !! Restart file identifier(unitless)
    REAL,DIMENSION (kjpindex),INTENT(in)           :: gtemp            !! First soil layer temperature
    REAL,DIMENSION (kjpindex),INTENT(in)           :: soilcap
    REAL,DIMENSION (kjpindex),INTENT(in)           :: soilflx
    REAL, DIMENSION (kjpindex), INTENT(in)         :: lambda_snow      !! Coefficient of the linear extrapolation of surface temperature 
    REAL, DIMENSION (kjpindex,nsnow), INTENT (in)  :: cgrnd_snow       !! Integration coefficient for snow numerical scheme
    REAL, DIMENSION (kjpindex,nsnow), INTENT (in)  :: dgrnd_snow       !! Integration coefficient for snow numerical scheme

!_ ================================================================================================================================
    
    !! 1. Write variables to restart file to be used for the next simulation
    IF (printlev>=3) WRITE (numout,*) 'Write restart file with THERMOSOIL variables'
    
    CALL restput_p(rest_id, 'ptn', nbp_glo, ngrnd, 1, kjit, ptn, 'scatter', nbp_glo, index_g)
    
    IF (ok_Ecorr) THEN
       CALL restput_p(rest_id, 'e_soil_lat', nbp_glo, 1 , 1, kjit, e_soil_lat, 'scatter', nbp_glo, index_g)
    END IF
    
    CALL restput_p(rest_id, 'gtemp', nbp_glo, 1, 1, kjit, gtemp, 'scatter', nbp_glo, index_g)
    CALL restput_p(rest_id, 'soilcap', nbp_glo, 1, 1, kjit, soilcap, 'scatter', nbp_glo, index_g)
    CALL restput_p(rest_id, 'soilflx', nbp_glo, 1, 1, kjit, soilflx, 'scatter', nbp_glo, index_g)
    CALL restput_p(rest_id, 'cgrnd', nbp_glo, ngrnd-1, 1, kjit, cgrnd, 'scatter', nbp_glo, index_g)
    CALL restput_p(rest_id, 'dgrnd', nbp_glo, ngrnd-1, 1, kjit, dgrnd, 'scatter', nbp_glo, index_g)
    CALL restput_p(rest_id, 'cgrnd_snow', nbp_glo, nsnow, 1, kjit, cgrnd_snow, 'scatter', nbp_glo, index_g)
    CALL restput_p(rest_id, 'dgrnd_snow', nbp_glo, nsnow, 1, kjit, dgrnd_snow, 'scatter', nbp_glo, index_g)
    CALL restput_p(rest_id, 'lambda_snow', nbp_glo, 1, 1, kjit, lambda_snow, 'scatter', nbp_glo, index_g)

  END SUBROUTINE thermosoilc_finalize


!! ================================================================================================================================
!! SUBROUTINE   : thermosoilc_clear
!!
!>\BRIEF        Desallocates the allocated arrays.
!! The call of thermosoilc_clear originates from sechiba_clear but the calling sequence and 
!! its purpose require further investigation.
!!
!! DESCRIPTION  : None
!!
!! RECENT CHANGE(S) : None 
!!
!! MAIN OUTPUT VARIABLE(S): None
!!
!! REFERENCE(S) : None 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

 SUBROUTINE thermosoilc_clear()

        IF ( ALLOCATED (ptn)) DEALLOCATE (ptn)
        IF ( ALLOCATED (cgrnd)) DEALLOCATE (cgrnd) 
        IF ( ALLOCATED (dgrnd)) DEALLOCATE (dgrnd) 
        IF ( ALLOCATED (pcapa)) DEALLOCATE (pcapa)
        IF ( ALLOCATED (pkappa))  DEALLOCATE (pkappa)
        IF ( ALLOCATED (pcapa_en)) DEALLOCATE (pcapa_en)
        IF ( ALLOCATED (ptn_beg)) DEALLOCATE (ptn_beg)
        IF ( ALLOCATED (temp_sol_beg)) DEALLOCATE (temp_sol_beg)
        IF ( ALLOCATED (surfheat_incr)) DEALLOCATE (surfheat_incr)
        IF ( ALLOCATED (coldcont_incr)) DEALLOCATE (coldcont_incr)
        IF ( ALLOCATED (shum_ngrnd_perma)) DEALLOCATE (shum_ngrnd_perma)
        IF ( ALLOCATED (profil_froz)) DEALLOCATE (profil_froz)

  END SUBROUTINE thermosoilc_clear


!! ================================================================================================================================
!! FUNCTION     : fz
!!
!>\BRIEF        fz(rk), the function's result, is the "rk"th element of a geometric series 
!! with first element fz1 and ration zalph.
!!
!! DESCRIPTION  : This function is used to calculate the depths of the boudaries of the thermal layers (zz_coef) and 
!! of the numerical nodes (zz) of the thermal scheme. Formulae to get the adimensional depths are followings :
!!      zz(jg)      = fz(REAL(jg,r_std) - undemi); \n
!!      zz_coef(jg) = fz(REAL(jg,r_std))
!!
!! RECENT CHANGE(S) : None
!!
!! RETURN VALUE : fz(rk)
!!
!! REFERENCE(S) : None 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  FUNCTION fz(rk) RESULT (fz_result)

  !! 0. Variables and parameter declaration

    !! 0.1 Input variables

    REAL, INTENT(in)                        :: rk
    
    !! 0.2 Output variables

    REAL                                    :: fz_result
    
    !! 0.3 Modified variables

    !! 0.4 Local variables

!_ ================================================================================================================================

    fz_result = fz1 * (zalph ** rk - un) / (zalph - un)

  END FUNCTION fz


!! ================================================================================================================================
!! FUNCTION     : thermosoilc_levels
!!
!>\BRIEF          Depth of nodes for the thermal layers in meters. 
!!
!! DESCRIPTION  : Calculate and return the depth in meters of the nodes of the soil layers. This calculation is the same 
!!                as done in thermosoilc_var_init for zz. See thermosoilc_var_init for more details. 
!!
!! RECENT CHANGE(S) : None
!!
!! RETURN VALUE : Vector of soil depth for the nodes in meters
!!
!! REFERENCE(S) : None 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  FUNCTION thermosoilc_levels() RESULT (zz_out)
    
    !! 0.1 Return variable
    
    REAL, DIMENSION (ngrnd)  :: zz_out      !! Depth of soil layers in meters
    
    !! 0.2 Local variables
    INTEGER                  :: jg
    REAL                     :: so_capa
    REAL                     :: so_cond
    
!_ ================================================================================================================================

    !! 1. Define some parameters
    so_capa = (so_capa_dry + so_capa_wet)/deux
    so_cond = (so_cond_dry + so_cond_wet)/deux
    
    cstgrnd=SQRT(one_day / pi)
    lskin = SQRT(so_cond / so_capa * one_day / pi)

    !! Parameters needed by fz function
    fz1 = 0.3_r_std * cstgrnd
    zalph = deux

    !!  2. Get adimentional depth of the numerical nodes
    DO jg=1,ngrnd
       zz_out(jg) = fz(jg - undemi)
    ENDDO
    
    !! 3. Convert to meters
    DO jg=1,ngrnd
       zz_out(jg) = zz_out(jg) /  cstgrnd * lskin
    END DO

  END FUNCTION thermosoilc_levels


!! ================================================================================================================================
!! SUBROUTINE   : thermosoilc_var_init
!!
!>\BRIEF        Define and initializes the soil thermal parameters
!!		  
!! DESCRIPTION	: This routine\n
!! 1. Defines the parameters ruling the vertical grid of the thermal scheme (fz1, zalpha).\n
!! 2. Defines the scaling coefficients for adimensional depths (lskin, cstgrnd, see explanation in the 
!!    variables description of thermosoilc_main). \n
!! 3. Calculates the vertical discretization of the soil (zz, zz_coef, dz2) and the constants used
!!    in the numerical scheme and which depend only on the discretization (dz1, lambda).\n
!! 4. Initializes the soil thermal parameters (capacity, conductivity) based on initial soil moisture content.\n
!! 5. Produces a first temperature diagnostic based on temperature initialization.\n
!!
!! The scheme comprizes ngrnd=7 layers by default.
!! The layer' s boundaries depths (zz_coef) follow a geometric series of ratio zalph=2 and first term fz1.\n
!! zz_coef(jg)=fz1.(1-zalph^jg)/(1-zalph) \n
!! The layers' boudaries depths are first calculated 'adimensionally', ie with a
!! discretization adapted to EQ5. This discretization is chosen for its ability at
!! reproducing a thermal signal with periods ranging from days to centuries. (see
!! Hourdin, 1992). Typically, fz1 is chosen as : fz1=0.3*cstgrnd (with cstgrnd the
!! adimensional attenuation depth). \n
!! The factor lskin/cstgrnd is then used to go from adimensional depths to
!! depths in m.\n
!! zz(real)=lskin/cstgrnd*zz(adimensional)\n
!! Similarly, the depths of the numerical nodes are first calculated
!! adimensionally, then the conversion factor is applied.\n
!! the numerical nodes (zz) are not exactly the layers' centers : their depths are calculated as follows:\n
!! zz(jg)=fz1.(1-zalph^(jg-1/2))/(1-zalph)\n
!! The values of zz and zz_coef used in the default thermal discretization are in the following table.
!! \latexonly
!! \includegraphics{thermosoilc_var_init1.jpg}
!! \endlatexonly\n
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S) : None
!!
!! REFERENCE(S)	:
!! - Hourdin, F. (1992). Study and numerical simulation of the general circulation of 
!! planetary atmospheres, Ph.D. thesis, Paris VII University.
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE thermosoilc_var_init(kjpindex, shumdiag_perma, snow, &
                                 zz, zz_coef, dz1, dz2, stempdiag)

  !! 0. Variables and parameter declaration

    !! 0.1 Input variables

    INTEGER, INTENT(in)                               :: kjpindex          !! Domain size (unitless)
    REAL, DIMENSION (kjpindex,nslm), INTENT (in)      :: shumdiag_perma    !! Relative soil humidity on the diagnostic axis 
                                                                                  !! (unitless), [0,1]. (see description of the 
                                                                                  !! variables of thermosoilc_main for more 
                                                                                  !! explanations) 
    REAL, DIMENSION (kjpindex), INTENT (in)     	     :: snow              !! Snow quantity   
   
    !! 0.2 Output variables

    REAL, DIMENSION (ngrnd), INTENT(out)              :: zz                !! depths of the layers'numerical nodes 
                                                                                  !! @tex ($m$)@endtex 
    REAL, DIMENSION (ngrnd), INTENT(out)              :: zz_coef		  !! depths of the layers'boundaries 
                                                                                  !! @tex ($m$)@endtex 
    REAL, DIMENSION (ngrnd), INTENT(out)              :: dz1               !! numerical constant depending on the vertical
                                                                                  !! thermal grid only @tex  ($m^{-1}$) @endtex. 
                                                                                  !! (see description
                                                                                  !! of the variables of thermosoilc_main for more
                                                                                  !! explanations)
    REAL, DIMENSION (ngrnd), INTENT(out)              :: dz2               !! thicknesses of the soil thermal layers 
                                                                                  !! @tex ($m$) @endtex
    REAL, DIMENSION (kjpindex,nslm), INTENT (out)     :: stempdiag         !! Diagnostic temperature profile @tex ($K$)
                                                                                  !! @endtex

    !! 0.3 Modified variables

    !! 0.4 Local variables

    INTEGER                                           :: ier, ji, jg       !! Index (unitless)
    REAL                                              :: sum               !! Temporary variable
    REAL                                              :: so_capa           !! Average Thermal Conductivity of soils 
                                                                                  !! @tex $(W.m^{-2}.K^{-1})$ @endtex
    REAL                                              :: so_cond           !! Average Thermal Conductivity of soils 
                                                                                  !! @tex $(W.m^{-2}.K^{-1})$ @endtex

!_ ================================================================================================================================

  !! 1. Initialization of the parameters of the vertical discretization and of the attenuation depths
    
    !! so_capa and so_cond are temporary variables which contain the average values of soil conductivity 
    !! and soil conductivity and are only needed in thermosoilc_var_init to set the vertical layering. 
    so_capa = (so_capa_dry + so_capa_wet)/deux
    so_cond = (so_cond_dry + so_cond_wet)/deux

    cstgrnd=SQRT(one_day / pi)
    lskin = SQRT(so_cond / so_capa * one_day / pi)
    fz1 = 0.3_r_std * cstgrnd
    zalph = deux

    !! Calculate so_capa_ice
    so_capa_ice = so_capa_dry + poros*capa_ice*rho_ice
    WRITE(numout,*) 'Calculation of so_capa_ice=', so_capa_ice,' using poros=',poros,' and capa_ice=',capa_ice
    
  !! 2.  Computing the depth of the thermal levels (numerical nodes) and the layers boundaries
   
    !! Computing the depth of the thermal levels (numerical nodes) and 
    !! the layers boundariesusing the so-called
    !! adimentional variable z' = z/lskin*cstgrnd (with z in m)
    
    !! 2.1 adimensional thicknesses of the layers
    DO jg=1,ngrnd

    !!?? code simplification hopefully possible here with up-to-date compilers !
    !!! This needs to be solved soon. Either we allow CPP options in SECHIBA or the VPP
    !!! fixes its compiler 
    !!!#ifdef VPP5000
      dz2(jg) = fz(jg-undemi+undemi) - fz(jg-1-undemi+undemi)
    !!!#else
    !!!      dz2(jg) = fz(REAL(jg,r_std)) - fz(REAL(jg-1,r_std))
    !!!#endif
    ENDDO
    
    !! 2.2 adimentional depth of the numerical nodes and layers' boudaries
    DO jg=1,ngrnd
      zz(jg)      = fz(jg - undemi)
      zz_coef(jg) = fz(jg -undemi+undemi)
    ENDDO

    !! 2.3 Converting to meters
    DO jg=1,ngrnd
      zz(jg)      = zz(jg) /  cstgrnd * lskin
      zz_coef(jg) = zz_coef(jg) / cstgrnd * lskin 
      dz2(jg)     = dz2(jg) /  cstgrnd * lskin
    ENDDO

    !! 2.4 Computing some usefull constants for the numerical scheme
    DO jg=1,ngrnd-1
      dz1(jg)  = un / (zz(jg+1) - zz(jg))
    ENDDO
    lambda = zz(1) * dz1(1)

    !! 2.5 Get the wetness profile on the thermal vertical grid from the diagnostic axis
    CALL thermosoilc_humlev(kjpindex, shumdiag_perma, snow)

  !! 3. Diagnostics : consistency checks on the vertical grid.

    WRITE (numout,*) 'diagnostic des niveaux dans le sol' !!?? to be changed,
    WRITE (numout,*) 'niveaux intermediaires et pleins'
    sum = zero
    DO jg=1,ngrnd
      sum = sum + dz2(jg)
      WRITE (numout,*) zz(jg),sum
    ENDDO

    
  !! 4. Compute a first diagnostic temperature profile
    
    CALL thermosoilc_diaglev(kjpindex, stempdiag)

    IF (printlev>=3) WRITE (numout,*) ' thermosoilc_var_init done '

  END SUBROUTINE thermosoilc_var_init
  

!! ================================================================================================================================
!! SUBROUTINE   : thermosoilc_coef
!!
!>\BRIEF        Calculate soil thermal properties, integration coefficients, apparent heat flux,
!! surface heat capacity,  
!!
!! DESCRIPTION	: This routine computes : \n
!!		1. the soil thermal properties. \n 
!!		2. the integration coefficients of the thermal numerical scheme, cgrnd and dgrnd,
!!              which depend on the vertical grid and on soil properties, and are used at the next 
!!              timestep.\n
!!              3. the soil apparent heat flux and surface heat capacity soilflux
!!              and soilcap, used by enerbil to compute the surface temperature at the next
!!              timestep.\n
!!             -  The soil thermal properties depend on water content (shum_ngrnd_perma, shumdiag_perma) and on the presence 
!!              of snow : snow is integrated into the soil for the thermal calculations, ie if there 
!!              is snow on the ground, the first thermal layer(s) consist in snow, depending on the 
!!              snow-depth. If a layer consists out of snow and soil, wheighed soil properties are 
!!              calculated\n
!!             - The coefficients cgrnd and dgrnd are the integration
!!              coefficients for the thermal scheme \n
!!                              T(k+1)=cgrnd(k)+dgrnd(k)*T(k) \n
!!                                      -- EQ1 -- \n
!!              They correspond respectively to $\beta$ and $\alpha$ from F. Hourdin\'s thesis and 
!!              their expression can be found in this document (eq A19 and A20)
!!             - soilcap and soilflux are the apparent surface heat capacity and flux
!!               used in enerbil at the next timestep to solve the surface
!!               balance for Ts (EQ3); they correspond to $C_s$ and $F_s$ in F.
!!               Hourdin\'s PhD thesis and are expressed in eq. A30 and A31. \n
!!                 soilcap*(Ts(t)-Ts(t-1))/dt=soilflux+otherfluxes(Ts(t)) \n
!!                                      -- EQ3 --\n
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S): cgrnd, dgrnd, pcapa, pkappa, soilcap, soilflx
!!
!! REFERENCE(S) :
!! - Hourdin, F. (1992). Study and numerical simulation of the general circulation of planetary atmospheres,
!! Ph.D. thesis, Paris VII University. Remark: the part of F. Hourdin's PhD thesis relative to the thermal
!! integration scheme has been scanned and is provided along with the documentation, with name : 
!! Hourdin_1992_PhD_thermal_scheme.pdf
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE thermosoilc_coef (kjpindex, temp_sol_new, snow,            &
                              ptn,      soilcap, soilflx,      zz,     &
                              dz1,      dz2,                           &
                              cgrnd,    dgrnd,                         &
                              frac_snow_veg,frac_snow_nobio,totfrac_nobio, &
                              snowdz,snowrho,  snowtemp,     pb,         &
			      lambda_snow,    cgrnd_snow,    dgrnd_snow)

    !! 0. Variables and parameter declaration

    !! 0.1 Input variables

    INTEGER, INTENT(in)                             :: kjpindex     !! Domain size (unitless)
    REAL, DIMENSION (kjpindex), INTENT (in)         :: temp_sol_new !! soil surface temperature @tex ($K$) @endtex
    REAL, DIMENSION (kjpindex), INTENT (in)         :: snow         !! snow mass @tex ($Kg$) @endtex
    REAL, DIMENSION (ngrnd), INTENT(in)             :: zz           !! depths of the soil thermal numerical nodes 
                                                                           !! @tex ($m$) @endtex 
    REAL, DIMENSION (ngrnd), INTENT(in)             :: dz1          !! numerical constant depending on the vertical 
                                                                           !! thermal grid only @tex ($m^{-1}$) @endtex 
    REAL, DIMENSION (ngrnd), INTENT(in)             :: dz2          !! thicknesses of the soil thermal layers
                                                                           !! @tex ($m$) @endtex 
    REAL,DIMENSION (kjpindex), INTENT(in)           :: frac_snow_veg   !! Snow cover fraction on vegeted area
    REAL,DIMENSION (kjpindex,nnobio), INTENT(in)    :: frac_snow_nobio !! Snow cover fraction on non-vegeted area
    REAL,DIMENSION (kjpindex),INTENT(in)            :: totfrac_nobio   !! Total fraction of continental ice+lakes+cities+...
                                                                              !!(unitless,0-1)
    REAL, DIMENSION (kjpindex,nsnow), INTENT(in)    :: snowdz          !! Snow depth
    REAL, DIMENSION (kjpindex,nsnow), INTENT(in)    :: snowrho         !! Snow density
    REAL, DIMENSION (kjpindex,nsnow), INTENT(in)    :: snowtemp        !! Snow temperature
    REAL, DIMENSION (kjpindex), INTENT (in)         :: pb              !! Surface presure (hPa)
    
    !! 0.2 Output variables

    REAL, DIMENSION (kjpindex), INTENT (out)        :: soilcap      !! surface heat capacity
                                                                           !! @tex ($J m^{-2} K^{-1}$) @endtex
    REAL, DIMENSION (kjpindex), INTENT (out)        :: soilflx      !! surface heat flux @tex ($W m^{-2}$) @endtex,
                                                                           !! positive towards the 
                                                                           !! soil, writen as Qg (ground heat flux) in the history 
                                                                           !! files.
    REAL, DIMENSION (kjpindex,ngrnd-1), INTENT(out) :: cgrnd        !! matrix coefficient for the computation of soil 
                                                                           !! temperatures (beta in F. Hourdin thesis)
    REAL, DIMENSION (kjpindex,ngrnd-1), INTENT(out) :: dgrnd        !! matrix coefficient for the computation of soil 
                                                                           !! temperatures (alpha in F. Hourdin thesis)


    !! 0.3 Modified variable

    REAL, DIMENSION (kjpindex,ngrnd), INTENT (inout):: ptn          !! vertically discretized soil temperatures
                                                                           !! @tex ($K$) @endtex
    REAL, DIMENSION (kjpindex), INTENT(inout)       :: lambda_snow  !! Coefficient of the linear extrapolation of surface temperature 
    REAL, DIMENSION (kjpindex,nsnow), INTENT (inout):: cgrnd_snow   !! Integration coefficient for snow numerical scheme
    REAL, DIMENSION (kjpindex,nsnow), INTENT (inout):: dgrnd_snow   !! Integration coefficient for snow numerical scheme

    !! 0.4 Local variables

    INTEGER                                         :: ji, jg
    REAL, DIMENSION (kjpindex,ngrnd-1)              :: zdz1         !! numerical (buffer) constant 
                                                                           !! @tex ($W m^{-1} K^{-1}$) @endtex
    REAL, DIMENSION (kjpindex,ngrnd)                :: zdz2         !! numerical (buffer) constant  
                                                                           !! @tex ($W m^{-1} K^{-1}$) @endtex
    REAL, DIMENSION (kjpindex)                      :: z1           !! numerical constant @tex ($W m^{-1} K^{-1}$) @endtex
    REAL, DIMENSION (kjpindex)                      :: soilcap_nosnow !! surface heat capacity
                                                                             !! @tex ($J m^{-2} K^{-1}$) @endtex
    REAL, DIMENSION (kjpindex)                      :: soilflx_nosnow !! surface heat flux @tex ($W m^{-2}$) @endtex,
                                                                             !! positive towards the soil, written as Qg (ground heat flux in the history files).
    REAL, DIMENSION (kjpindex)                      :: snowcap        !! apparent snow heat capacity @tex ($J m^{-2} K^{-1}$)
    REAL, DIMENSION (kjpindex)                      :: snowflx        !! apparent snow-atmosphere heat flux @tex ($W m^{-2}$) @endtex
    REAL, DIMENSION (kjpindex,nsnow)                :: dz1_snow
    REAL, DIMENSION (kjpindex,nsnow)                :: ZSNOWDZM
    REAL, DIMENSION (kjpindex,nsnow)                :: dz2_snow
    REAL, DIMENSION (kjpindex,nsnow)                :: zdz1_snow
    REAL, DIMENSION (kjpindex,nsnow)                :: zdz2_snow
    REAL, DIMENSION (kjpindex)                      :: z1_snow

!_ ================================================================================================================================

  !! 1. Computation of the soil thermal properties
   
    ! Computation of the soil thermal properties; snow properties are also accounted for
    IF (ok_freeze_thermix) THEN
       CALL thermosoilc_getdiff( kjpindex, snow, ptn, snowrho, snowtemp, pb )
    ELSE IF (ok_explicitsnow) THEN
       ! Special case with explicit snow model without soil freezing
       CALL thermosoilc_getdiff_old_thermix_without_snow( kjpindex, snowrho, snowtemp, pb )
    ELSE
       ! Special case with old snow without soil freezing
       CALL thermosoilc_getdiff_old_thermix_with_snow( kjpindex, snow )
    ENDIF

    ! Energy conservation : Correction to make sure that the same latent heat is released and 
    ! consumed during freezing and thawing
    IF (ok_Ecorr) THEN
       CALL thermosoilc_readjust(kjpindex, ptn)
    ENDIF
    

    !! 2. Computation of the coefficients of the numerical integration scheme for the soil layers

    !! 2.1 Calculate numerical coefficients zdz1 and zdz2
    DO jg=1,ngrnd
      DO ji=1,kjpindex
        zdz2(ji,jg)=pcapa(ji,jg) * dz2(jg)/dt_sechiba
      ENDDO
    ENDDO
    
    DO jg=1,ngrnd-1
      DO ji=1,kjpindex
        zdz1(ji,jg) = dz1(jg) * pkappa(ji,jg)
      ENDDO
    ENDDO
    
    !! 2.2 Calculate coefficients cgrnd and dgrnd for soil
    DO ji = 1,kjpindex
      z1(ji) = zdz2(ji,ngrnd) + zdz1(ji,ngrnd-1)
      cgrnd(ji,ngrnd-1) = zdz2(ji,ngrnd) * ptn(ji,ngrnd) / z1(ji)
      dgrnd(ji,ngrnd-1) = zdz1(ji,ngrnd-1) / z1(ji)
    ENDDO

    DO jg = ngrnd-1,2,-1
      DO ji = 1,kjpindex
        z1(ji) = un / (zdz2(ji,jg) + zdz1(ji,jg-1) + zdz1(ji,jg) * (un - dgrnd(ji,jg)))
        cgrnd(ji,jg-1) = (ptn(ji,jg) * zdz2(ji,jg) + zdz1(ji,jg) * cgrnd(ji,jg)) * z1(ji)
        dgrnd(ji,jg-1) = zdz1(ji,jg-1) * z1(ji)
      ENDDO
    ENDDO


    !! 3. Computation of the coefficients of the numerical integration scheme for the snow layers

    !! 3.1 Calculate numerical coefficients zdz1_snow, zdz2_snow and lambda_snow
    DO ji = 1, kjpindex

       dz2_snow(ji,:) = 0

       IF ( ok_explicitsnow ) THEN

          ! Calculate internal values
          DO jg = 1, nsnow
             ZSNOWDZM(ji,jg) = MAX(snowdz(ji,jg),psnowdzmin)
          ENDDO
          dz2_snow(ji,:)=ZSNOWDZM(ji,:)

          DO jg = 1, nsnow-1
             dz1_snow(ji,jg)  = 2.0 / (dz2_snow(ji,jg+1)+dz2_snow(ji,jg))
          ENDDO

          lambda_snow(ji) = dz2_snow(ji,1)/2.0 * dz1_snow(ji,1)

          DO jg=1,nsnow
             zdz2_snow(ji,jg)=pcapa_snow(ji,jg) * dz2_snow(ji,jg)/dt_sechiba
          ENDDO

          DO jg=1,nsnow-1
             zdz1_snow(ji,jg) = dz1_snow(ji,jg) * pkappa_snow(ji,jg)
          ENDDO

          ! the bottom snow
          zdz1_snow(ji,nsnow) = pkappa_snow(ji,nsnow) / ( zz_coef(1) + dz2_snow(ji,nsnow)/2 )

       ELSE
          ! There is no snow, these values will never be used
          lambda_snow(ji) = lambda
       ENDIF

    ENDDO

    !! 3.2 Calculate coefficients cgrnd_snow and dgrnd_snow for snow
    DO ji = 1,kjpindex
       IF ( ok_explicitsnow ) THEN
          ! bottom level
          z1_snow(ji) = zdz2(ji,1)+(un-dgrnd(ji,1))*zdz1(ji,1)+zdz1_snow(ji,nsnow)
          cgrnd_snow(ji,nsnow) = (zdz2(ji,1) * ptn(ji,1) + zdz1(ji,1) * cgrnd(ji,1) ) / z1_snow(ji)
          dgrnd_snow(ji,nsnow) = zdz1_snow(ji,nsnow) / z1_snow(ji)

          ! next-to-bottom level
          z1_snow(ji) = zdz2_snow(ji,nsnow)+(un-dgrnd_snow(ji,nsnow))*zdz1_snow(ji,nsnow)+zdz1_snow(ji,nsnow-1)
          cgrnd_snow(ji,nsnow-1) = (zdz2_snow(ji,nsnow)*snowtemp(ji,nsnow)+&
               zdz1_snow(ji,nsnow)*cgrnd_snow(ji,nsnow))/z1_snow(ji)
          dgrnd_snow(ji,nsnow-1) = zdz1_snow(ji,nsnow-1) / z1_snow(ji)

          DO jg = nsnow-1,2,-1
             z1_snow(ji) = un / (zdz2_snow(ji,jg) + zdz1_snow(ji,jg-1) + zdz1_snow(ji,jg) * (un - dgrnd_snow(ji,jg)))
             cgrnd_snow(ji,jg-1) = (snowtemp(ji,jg) * zdz2_snow(ji,jg) + zdz1_snow(ji,jg) * cgrnd_snow(ji,jg)) * z1_snow(ji)
             dgrnd_snow(ji,jg-1) = zdz1_snow(ji,jg-1) * z1_snow(ji)
          ENDDO
       ELSE
          ! There is no snow, these values will never be used
          cgrnd_snow(ji,:) = cgrnd(ji,1)
          dgrnd_snow(ji,:) = dgrnd(ji,1)
       ENDIF
    ENDDO


    !! 4. Computation of the apparent ground heat flux 
    !! Computation of apparent snow-atmosphere flux  
    DO ji = 1,kjpindex
       IF ( ok_explicitsnow ) THEN
          snowflx(ji) = zdz1_snow(ji,1) * (cgrnd_snow(ji,1) + (dgrnd_snow(ji,1)-1.) * snowtemp(ji,1))
          snowcap(ji) = (zdz2_snow(ji,1) * dt_sechiba + dt_sechiba * (un - dgrnd_snow(ji,1)) * zdz1_snow(ji,1))
          z1_snow(ji) = lambda_snow(ji) * (un - dgrnd_snow(ji,1)) + un 
          snowcap(ji) = snowcap(ji) / z1_snow(ji)
          snowflx(ji) = snowflx(ji) + &
               & snowcap(ji) * (snowtemp(ji,1) * z1_snow(ji) - lambda_snow(ji) * cgrnd_snow(ji,1) - temp_sol_new(ji)) / dt_sechiba
       ELSE
          snowflx(ji) = zero
          snowcap(ji) = zero
       ENDIF
    ENDDO


    !! Computation of the apparent ground heat flux (> towards the soil) and
    !! apparent surface heat capacity, used at the next timestep by enerbil to
    !! compute the surface temperature.
    DO ji = 1,kjpindex
      soilflx_nosnow(ji) = zdz1(ji,1) * (cgrnd(ji,1) + (dgrnd(ji,1)-1.) * ptn(ji,1))
      soilcap_nosnow(ji) = (zdz2(ji,1) * dt_sechiba + dt_sechiba * (un - dgrnd(ji,1)) * zdz1(ji,1))
      z1(ji) = lambda * (un - dgrnd(ji,1)) + un
      soilcap_nosnow(ji) = soilcap_nosnow(ji) / z1(ji)
      soilflx_nosnow(ji) = soilflx_nosnow(ji) + &
         & soilcap_nosnow(ji) * (ptn(ji,1) * z1(ji) - lambda * cgrnd(ji,1) - temp_sol_new(ji)) / dt_sechiba
    ENDDO

    !! Add snow fraction
    IF ( ok_explicitsnow ) THEN
       ! Using an effective heat capacity and heat flux by a simple pondering of snow and soil fraction
       DO ji = 1, kjpindex
          soilcap(ji) = snowcap(ji)*frac_snow_veg(ji)*(1-totfrac_nobio(ji))+    & ! weights related to snow cover fraction on vegetation  
               soilcap_nosnow(ji)*SUM(frac_snow_nobio(ji,:))*totfrac_nobio(ji)+ & ! weights related to SCF on nobio
               soilcap_nosnow(ji)*(1-(frac_snow_veg(ji)*(1-totfrac_nobio(ji))+SUM(frac_snow_nobio(ji,:))*totfrac_nobio(ji))) ! weights related to non snow fraction
          soilflx(ji) = snowflx(ji)*frac_snow_veg(ji)*(1-totfrac_nobio(ji))+    & ! weights related to snow cover fraction on vegetation  
               soilflx_nosnow(ji)*SUM(frac_snow_nobio(ji,:))*totfrac_nobio(ji)+ & ! weights related to SCF on nobio
               soilflx_nosnow(ji)*(1-(frac_snow_veg(ji)*(1-totfrac_nobio(ji))+SUM(frac_snow_nobio(ji,:))*totfrac_nobio(ji))) ! weights related to non snow fraction
       ENDDO
    ELSE
       ! Do not consider snow fraction
       soilcap(:)=soilcap_nosnow(:)
       soilflx(:)=soilflx_nosnow(:)
    END IF

    IF (printlev>=3) WRITE (numout,*) ' thermosoilc_coef done '

  END SUBROUTINE thermosoilc_coef
 
 
!! ================================================================================================================================
!! SUBROUTINE   : thermosoilc_profile
!!
!>\BRIEF        In this routine solves the numerical soil thermal scheme, ie calculates the new soil temperature profile; 
!! This profile is then exported onto the diagnostic axis (call thermosoilc_diaglev)
!!
!! DESCRIPTION	: The calculation of the new soil temperature profile is based on
!! the cgrnd and dgrnd values from the previous timestep and the surface temperature Ts aka temp_sol_new. (see detailed
!! explanation in the header of the thermosoilc module or in the reference).\n
!!                              T(k+1)=cgrnd(k)+dgrnd(k)*T(k)\n
!!                                      -- EQ1 --\n
!!                           Ts=(1-lambda)*T(1) -lambda*T(2)\n 
!!                                      -- EQ2--\n
!!
!! RECENT CHANGE(S) : None
!! 
!! MAIN OUTPUT VARIABLE(S): ptn (soil temperature profile on the thermal axis), 
!!                          stempdiag (soil temperature profile on the diagnostic axis)
!!
!! REFERENCE(S) :
!! - Hourdin, F. (1992). Study and numerical simulation of the general circulation of planetary atmospheres,
!! Ph.D. thesis, Paris VII University. Remark: the part of F. Hourdin's PhD thesis relative to the thermal
!! integration scheme has been scanned and is provided along with the documentation, with name : 
!! Hourdin_1992_PhD_thermal_scheme.pdf
!!
!! FLOWCHART    : None 
!! \n 
!_ ================================================================================================================================

 SUBROUTINE thermosoilc_profile (kjpindex, temp_sol_new, ptn, stempdiag, &
                                 snowtemp, frac_snow_veg,frac_snow_nobio,totfrac_nobio, &
                                 cgrnd_snow,dgrnd_snow)

  !! 0. Variables and parameter declaration

    !! 0.1 Input variables

    INTEGER, INTENT(in)                               :: kjpindex       !! Domain size (unitless)
    REAL,DIMENSION (kjpindex), INTENT (in)            :: temp_sol_new   !! Surface temperature at the present time-step 
                                                                               !! @tex ($K$) @endtex
    REAL,DIMENSION (kjpindex), INTENT(in)             :: frac_snow_veg  !! Snow cover fraction on vegeted area
    REAL,DIMENSION (kjpindex,nnobio), INTENT(in)      :: frac_snow_nobio!! Snow cover fraction on non-vegeted area
    REAL,DIMENSION (kjpindex),INTENT(in)              :: totfrac_nobio  !! Total fraction of continental ice+lakes+cities+...(unitless,0-1)
    REAL,DIMENSION (kjpindex,nsnow), INTENT(in)       :: snowtemp       !! Snow temperature
    REAL,DIMENSION (kjpindex,nsnow), INTENT(in)       :: cgrnd_snow     !! Integration coefficient for snow numerical scheme
    REAL,DIMENSION (kjpindex,nsnow), INTENT(in)       :: dgrnd_snow     !! Integration coefficient for snow numerical scheme

   
    !! 0.2 Output variables
    REAL,DIMENSION (kjpindex,nslm), INTENT (out)      :: stempdiag      !! diagnostic temperature profile 
                                                                               !! @tex ($K$) @endtex
    REAL,DIMENSION (kjpindex,ngrnd), INTENT (out)     :: ptn            !! vertically discretized soil temperatures 
                                                                               !! @tex ($K$) @endtex

    !! 0.3 Modified variables


    !! 0.4 Local variables

    INTEGER                                           :: ji, jg
    REAL                                              :: temp_sol_eff
     
!_ ================================================================================================================================
    
  !! 1. Computes the soil temperatures ptn.

    !! 1.1. ptn(jg=1) using EQ1 and EQ2
    DO ji = 1,kjpindex

       IF ( ok_explicitsnow ) THEN
          ! Soil temperature calculation with explicit snow if there is snow on the ground
          ! using an effective surface temperature by a simple pondering 
          temp_sol_eff=snowtemp(ji,nsnow)*frac_snow_veg(ji)*(1-totfrac_nobio(ji))+ & ! weights related to snow cover fraction on vegetation  
               temp_sol_new(ji)*SUM(frac_snow_nobio(ji,:))*totfrac_nobio(ji)+ & ! weights related to SCF on nobio
               temp_sol_new(ji)*(1-(frac_snow_veg(ji)*(1-totfrac_nobio(ji))+SUM(frac_snow_nobio(ji,:))*totfrac_nobio(ji))) ! weights related to non snow fraction
          ! Soil temperature calculation with explicit snow if there is snow on the ground
          ptn(ji,1) = cgrnd_snow(ji,nsnow) + dgrnd_snow(ji,nsnow) * temp_sol_eff
       ELSE
          ! Standard soil temperature calculation
          ptn(ji,1) = (lambda * cgrnd(ji,1) + temp_sol_new(ji)) / (lambda *(un - dgrnd(ji,1)) + un)
       ENDIF
    ENDDO

    !! 1.2. ptn(jg=2:ngrnd) using EQ1.
    DO jg = 1,ngrnd-1
      DO ji = 1,kjpindex
        ptn(ji,jg+1) = cgrnd(ji,jg) + dgrnd(ji,jg) * ptn(ji,jg)
      ENDDO
    ENDDO

  !! 2. Put the soil temperatures onto the diagnostic axis 
  
    !! Put the soil temperatures onto the diagnostic axis for convenient
    !! use in other routines (stomate..)
    CALL thermosoilc_diaglev(kjpindex, stempdiag)

    IF (printlev>=3) WRITE (numout,*) ' thermosoilc_profile done '

  END SUBROUTINE thermosoilc_profile


!! ================================================================================================================================
!! SUBROUTINE   : thermosoilc_diaglev
!!
!>\BRIEF        Interpolation of the soil in-depth temperatures onto the diagnostic profile.
!!
!! DESCRIPTION  : This is a very easy linear interpolation, with intfact(jsl, jg) the fraction
!! the thermal layer jg comprised within the diagnostic layer jsl. The depths of
!! the diagnostic levels are diaglev(1:nslm), computed in slowproc.f90.
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S): stempdiag (soil temperature profile on the diagnostic axis)
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None 
!! \n 
!_ ================================================================================================================================

  SUBROUTINE thermosoilc_diaglev(kjpindex, stempdiag)
    
  !! 0. Variables and parameter declaration

    !! 0.1 Input variables
 
    INTEGER, INTENT(in)                          :: kjpindex       !! Domain size (unitless)
   
    !! 0.2 Output variables

    REAL,DIMENSION (kjpindex,nslm), INTENT (out) :: stempdiag      !! Diagnostoc soil temperature profile @tex ($K$) @endtex
    
    !! 0.3 Modified variables

    !! 0.4 Local variables

    INTEGER                                      :: ji, jsl, jg
    REAL                                         :: lev_diag, prev_diag, lev_prog, prev_prog
    REAL, SAVE, ALLOCATABLE, DIMENSION(:,:)      :: intfact
!$OMP THREADPRIVATE(intfact)
    LOGICAL, PARAMETER                                  :: check=.FALSE.
!_ ================================================================================================================================
    
  !! 1. Computes intfact(jsl, jg)

    !! Computes intfact(jsl, jg), the fraction
    !! the thermal layer jg comprised within the diagnostic layer jsl.
    IF ( .NOT. ALLOCATED(intfact)) THEN
        
        ALLOCATE(intfact(nslm, ngrnd))
        
        prev_diag = zero
        DO jsl = 1, nslm
          lev_diag = diaglev(jsl)
          prev_prog = zero
          DO jg = 1, ngrnd
             IF ( jg == ngrnd .AND. (prev_prog + dz2(jg)) < lev_diag ) THEN
                lev_prog = lev_diag
             ELSE
                lev_prog = prev_prog + dz2(jg)
             ENDIF
            intfact(jsl,jg) = MAX(MIN(lev_diag,lev_prog)-MAX(prev_diag, prev_prog), zero)/(lev_diag-prev_diag)
            prev_prog = lev_prog
          ENDDO
           prev_diag = lev_diag
        ENDDO
        
        IF ( check ) THEN
           WRITE(numout,*) 'thermosoilc_diagev -- thermosoilc_diaglev -- thermosoilc_diaglev --' 
           DO jsl = 1, nslm
              WRITE(numout,*) jsl, '-', intfact(jsl,1:ngrnd)
           ENDDO
           WRITE(numout,*) "SUM -- SUM -- SUM SUM -- SUM -- SUM"
           DO jsl = 1, nslm
              WRITE(numout,*) jsl, '-', SUM(intfact(jsl,1:ngrnd))
           ENDDO
           WRITE(numout,*) 'thermosoilc_diaglev -- thermosoilc_diaglev -- thermosoilc_diaglev --' 
        ENDIF
        
    ENDIF

 !! 2. does the interpolation

    stempdiag(:,:) = zero
    DO jg = 1, ngrnd
      DO jsl = 1, nslm
        DO ji = 1, kjpindex
          stempdiag(ji,jsl) = stempdiag(ji,jsl) + ptn(ji,jg)*intfact(jsl,jg)
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE thermosoilc_diaglev
  

!! ================================================================================================================================
!! SUBROUTINE   : thermosoilc_humlev
!!
!>\BRIEF        Interpolates the diagnostic soil humidity profile shumdiag_perma(nslm, diagnostic axis) onto 
!!              the thermal axis, which gives shum_ngrnd_perma(ngrnd, thermal axis).
!!
!! DESCRIPTION  : Same as in thermosoilc_diaglev : This is a very easy linear interpolation, with intfactw(jsl, jg) the fraction
!! the thermal layer jsl comprised within the diagnostic layer jg. 
!! 
!! The depths of the diagnostic levels are diaglev(1:nslm), computed in slowproc.f90.
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S): shum_ngrnd_perma (soil saturation degree on the thermal axis)
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None 
!! \n 
!_ ================================================================================================================================
  SUBROUTINE thermosoilc_humlev(kjpindex, shumdiag_perma, snow)
  
  !! 0. Variables and parameter declaration

    !! 0.1 Input variables
 
    INTEGER, INTENT(in)                            :: kjpindex    !! Domain size (unitless)
    REAL,DIMENSION (kjpindex,nslm), INTENT (in)    :: shumdiag_perma  
    REAL,DIMENSION (kjpindex), INTENT (in)         :: snow 
    
    !! 0.2 Output variables

    !! 0.3 Modified variables

    !! 0.4 Local variables
    INTEGER                                       :: ji, jsl, jg
    REAL                                          :: lev_diag, prev_diag, lev_prog, prev_prog
    REAL, DIMENSION(ngrnd,nslm)                   :: intfactw     !! fraction of each diagnostic layer (jsl) comprized within
                                                                         !! a given thermal layer (jg)(0-1, unitless) 
    REAL, DIMENSION(kjpindex)               :: snow_h
    LOGICAL, PARAMETER :: check=.FALSE.

!_ ================================================================================================================================
    
  !! 1. computes intfactw(jsl,jg), the fraction of each diagnostic layer (jg) comprized within a given thermal layer (jsl)
    IF ( check ) &
         WRITE(numout,*) 'thermosoilc_humlev --' 

    ! Snow height
    snow_h(:)=snow(:)/sn_dens
    !
    shum_ngrnd_perma(:,:) = zero
    DO ji=1,kjpindex
       prev_diag = zero
       DO jsl = 1, ngrnd
          lev_diag = prev_diag + dz2(jsl)
          prev_prog = snow_h(ji)
          DO jg = 1, nslm
             IF ( jg == nslm .AND. diaglev(jg)+snow_h(ji) < lev_diag ) THEN
                lev_prog = lev_diag+snow_h(ji)
             ELSE
                lev_prog = diaglev(jg)+snow_h(ji)
             ENDIF
             intfactw(jsl,jg) = MAX(MIN(lev_diag,lev_prog)-MAX(prev_diag, prev_prog), zero)/(lev_diag-prev_diag)
             prev_prog = lev_prog
          ENDDO
          prev_diag = lev_diag
       ENDDO

       DO jsl = 1, nslm
          DO jg = 1, ngrnd
             shum_ngrnd_perma(ji,jg) = shum_ngrnd_perma(ji,jg) + shumdiag_perma(ji,jsl)*intfactw(jg,jsl)
          ENDDO
       ENDDO
       !
       IF ( check ) THEN
          DO jg = 1, ngrnd
             WRITE(numout,*) ji,jg, '-', intfactw(jg,1:nslm),'-sum-', SUM(intfactw(jg,1:nslm))
          ENDDO
       ENDIF
    ENDDO

    IF ( check ) &
         WRITE(numout,*) 'thermosoilc_humlev --' 

  END SUBROUTINE thermosoilc_humlev


!! ================================================================================================================================
!! SUBROUTINE   : thermosoilc_energy
!!
!>\BRIEF         Energy check-up.
!!
!! DESCRIPTION  : I didn\'t comment this routine since at do not understand its use, please
!! ask initial designers (Jan ? Nathalie ?).
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S) : ??
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None 
!! \n 
!_ ================================================================================================================================

  SUBROUTINE thermosoilc_energy(kjpindex, temp_sol_new, soilcap)
  
   !! 0. Variables and parameter declaration

    !! 0.1 Input variables

    INTEGER, INTENT(in)                     :: kjpindex     !! Domain size (unitless)
    REAL,DIMENSION (kjpindex), INTENT (in)  :: temp_sol_new !! Surface temperature at the present time-step, Ts 
                                                                   !! @tex ($K$) @endtex
    REAL,DIMENSION (kjpindex), INTENT (in)  :: soilcap      !! Apparent surface heat capacity 
                                                                   !! @tex ($J m^{-2} K^{-1}$) @endtex, 
                                                                   !! see eq. A29 of F. Hourdin\'s PhD thesis.
    
    !! 0.2 Output variables

    !! 0.3 Modified variables
    
    !! 0.4 Local variables
    
    INTEGER                                 :: ji, jg
!_ ================================================================================================================================
   

     DO ji = 1, kjpindex
      surfheat_incr(ji) = zero
      coldcont_incr(ji) = zero
     ENDDO
    
    !  Sum up the energy content of all layers in the soil.
    DO ji = 1, kjpindex
   
       IF (pcapa_en(ji,1) .LE. sn_capa) THEN
          
          ! Verify the energy conservation in the surface layer
          coldcont_incr(ji) = soilcap(ji) * (temp_sol_new(ji) - temp_sol_beg(ji))
          surfheat_incr(ji) = zero
       ELSE
          
          ! Verify the energy conservation in the surface layer
          surfheat_incr(ji) = soilcap(ji) * (temp_sol_new(ji) - temp_sol_beg(ji))
          coldcont_incr(ji) = zero
       ENDIF
    ENDDO
    
    ptn_beg(:,:)      = ptn(:,:)
    temp_sol_beg(:)   = temp_sol_new(:)

  END SUBROUTINE thermosoilc_energy



!! ================================================================================================================================
!! SUBROUTINE   : thermosoilc_readjust
!!
!>\BRIEF        
!!
!! DESCRIPTION	: Energy conservation : Correction to make sure that the same latent heat is released and 
!!                consumed during freezing and thawing  
!!
!! RECENT CHANGE(S) : None
!! 
!! MAIN OUTPUT VARIABLE(S): ptn (soil temperature profile on the thermal axis), 
!!                          
!! REFERENCE(S) :
!!
!! FLOWCHART    : None 
!! \n 
!_ ================================================================================================================================

  SUBROUTINE thermosoilc_readjust(kjpindex, ptn)

   !! 0. Variables and parameter declaration

    !! 0.1 Input variables
    INTEGER, INTENT(in)                             :: kjpindex

    !! 0.2 Modified variables
    REAL,DIMENSION(kjpindex,ngrnd),INTENT(inout)    :: ptn

    !! 0.3 Local variables
    INTEGER  :: ji, jg
    REAL :: ptn_tmp

    DO jg=1, ngrnd
       DO ji=1, kjpindex
          ! All soil latent energy is put into e_soil_lat(ji)
          ! because the variable soil layers make it difficult to keep track of all
          ! layers in this version
          ! NOTE : pcapa has unit J/K/m3 and pcappa_supp has J/K
          e_soil_lat(ji)=e_soil_lat(ji)+pcappa_supp(ji,jg)*(ptn(ji,jg)-ptn_beg(ji,jg))
       END DO
   END DO

   DO ji=1, kjpindex
      IF (e_soil_lat(ji).GT.min_sechiba.AND.MINVAL(ptn(ji,:)).GT.ZeroCelsius+fr_dT/2.) THEN
         ! The soil is thawed: we spread the excess of energy over the uppermost 6 levels e.g. 2.7m
         ! Here we increase the temperatures
         DO jg=1,6
            ptn_tmp=ptn(ji,jg)
            
            ptn(ji,jg)=ptn(ji,jg)+MIN(e_soil_lat(ji)/pcapa(ji,jg)/zz_coef(6), 0.5)
            e_soil_lat(ji)=e_soil_lat(ji)-(ptn(ji,jg)-ptn_tmp)*pcapa(ji,jg)*dz2(jg)
         ENDDO
      ELSE IF (e_soil_lat(ji).LT.-min_sechiba.AND.MINVAL(ptn(ji,:)).GT.ZeroCelsius+fr_dT/2.) THEN
         ! The soil is thawed
         ! Here we decrease the temperatures
         DO jg=1,6
            ptn_tmp=ptn(ji,jg)
            ptn(ji,jg)=MAX(ZeroCelsius+fr_dT/2., ptn_tmp+e_soil_lat(ji)/pcapa(ji,jg)/zz_coef(6))
            e_soil_lat(ji)=e_soil_lat(ji)+(ptn_tmp-ptn(ji,jg))*pcapa(ji,jg)*dz2(jg)
         END DO
      END IF
   END DO

  END SUBROUTINE thermosoilc_readjust
   
!-------------------------------------------------------------------



!! ================================================================================================================================
!! SUBROUTINE   : thermosoilc_getdiff
!!
!>\BRIEF          Computes soil and snow heat capacity and conductivity    
!!
!! DESCRIPTION	: Computation of the soil and snow thermal properties; snow properties are also accounted for
!!
!! RECENT CHANGE(S) : None
!! 
!! MAIN OUTPUT VARIABLE(S): 
!!                          
!! REFERENCE(S) :
!!
!! FLOWCHART    : None 
!! \n 
!_ ================================================================================================================================

  SUBROUTINE thermosoilc_getdiff( kjpindex, snow, ptn, snowrho, snowtemp, pb )

   !! 0. Variables and parameter declaration

    !! 0.1 Input variables
    INTEGER,INTENT(in)				:: kjpindex
    REAL,DIMENSION(kjpindex),INTENT (in)	        :: snow       !! Snow mass
    REAL,DIMENSION(kjpindex,nsnow),INTENT(in)    :: snowrho    !! Snow density
    REAL,DIMENSION(kjpindex,nsnow),INTENT(in)    :: snowtemp   !! Snow temperature
    REAL,DIMENSION(kjpindex),INTENT(in)          :: pb         !! Surface pressure
    REAL,DIMENSION(kjpindex,ngrnd),INTENT(in)	:: ptn        !! Soil temperature profile

    !! 0.3 Local variables
    REAL						:: xx         !! Unfrozen fraction of the soil
    REAL, DIMENSION(kjpindex)             	:: snow_h
    REAL, DIMENSION(kjpindex,ngrnd) 		:: zx1, zx2    
    REAL			   			:: cap_iw     !! Heat capacity of ice/water mixture 
    REAL						:: csat       !! Thermal conductivity for saturated soil
    INTEGER						:: ji,jg


    DO ji = 1,kjpindex
       IF (.NOT. ok_explicitsnow) THEN      
          ! 1. Determine the fractions of snow and soil
          snow_h(ji) = snow(ji) / sn_dens
      
          !
          !  1.1. The first level
          !
          IF ( snow_h(ji) .GT. zz_coef(1) ) THEN
             ! the 1st level is in the snow => the 1st layer is entirely snow
             zx1(ji,1) = 1.
             zx2(ji,1) = 0.
          ELSE IF ( snow_h(ji) .GT. zero ) THEN      
             ! the 1st level is beyond the snow and the snow is present
             zx1(ji,1) = snow_h(ji) / zz_coef(1)
             zx2(ji,1) = ( zz_coef(1) - snow_h(ji)) / zz_coef(1)	
          ELSE
             ! there is no snow at all, quoi ;-)
             zx1(ji,1) = 0.
             zx2(ji,1) = 1.       
          ENDIF
      
          !
          !  1.2. The other levels
          !
          DO jg = 2, ngrnd
             IF ( snow_h(ji) .GT. zz_coef(jg) ) THEN
                ! the current level is in the snow => the current layer is entirely snow
                zx1(ji,jg) = 1.
                zx2(ji,jg) = 0.
             ELSE IF ( snow_h(ji) .GT. zz_coef(jg-1) ) THEN
                ! the current layer is partially snow
                zx1(ji,jg) = (snow_h(ji) - zz_coef(jg-1)) / (zz_coef(jg) - zz_coef(jg-1))
                zx2(ji,jg) = ( zz_coef(jg) - snow_h(ji)) / (zz_coef(jg) - zz_coef(jg-1))
             ELSE
                ! both levels are out of snow => the current layer is entirely soil	  
                zx1(ji,jg) = 0.
                zx2(ji,jg) = 1.       
             ENDIF
          ENDDO
       ELSE
          zx1(ji,:) = 0.
          zx2(ji,:) = 1.
       END IF

      DO jg = 1, ngrnd
         !
         ! 2. Calculate heat capacity with allowance for permafrost
         ! 2.1. soil heat capacity depending on temperature and humidity
	 
         IF (ptn(ji,jg) .LT. ZeroCelsius-fr_dT/2.) THEN
	    ! frozen soil
            pcapa(ji,jg) = so_capa_dry + shum_ngrnd_perma(ji,jg)*(so_capa_ice - so_capa_dry)!Isa : old version, proved to be correct
	    profil_froz(ji,jg) = 1.
 	    pcappa_supp(ji,jg)= 0.

     	 ELSEIF (ptn(ji,jg) .GT. ZeroCelsius+fr_dT/2.) THEN
	    ! unfrozen soil	 
     	    pcapa(ji,jg) =  so_capa_dry + shum_ngrnd_perma(ji,jg)*(so_capa_wet - so_capa_dry)!Isa : old version, proved to be correct
	    profil_froz(ji,jg) = 0.
 	    pcappa_supp(ji,jg)= 0.
     	 ELSE
     	   ! xx is the unfrozen fraction of soil water	   	   
     	   xx = (ptn(ji,jg)-(ZeroCelsius-fr_dT/2.)) / fr_dT
           profil_froz(ji,jg) = (1. - xx)

	   ! net heat capacity of the ice/water mixture
     	   cap_iw = xx * so_capa_wet + (1.-xx) * so_capa_ice
	   pcapa(ji,jg) = so_capa_dry + shum_ngrnd_perma(ji,jg)*(cap_iw-so_capa_dry) + shum_ngrnd_perma(ji,jg)*poros*lhf*rho_water/fr_dT
 	   pcappa_supp(ji,jg)= shum_ngrnd_perma(ji,jg)*poros*lhf*rho_water/fr_dT*zx2(ji,jg)*dz2(jg)

         ENDIF

         !
	 ! 2.2. Take into account the snow and soil fractions in the layer
	 !
         pcapa(ji,jg) = zx1(ji,jg) * sn_capa + zx2(ji,jg) * pcapa(ji,jg)

	 !
	 ! 2.3. Calculate the heat capacity for energy conservation check 
	 IF ( zx1(ji,jg).GT.0. ) THEN
            pcapa_en(ji,jg) = sn_capa
	 ELSE
            pcapa_en(ji,jg) = pcapa(ji,jg)
	 ENDIF
 
         !
         ! 3. Calculate the heat conductivity with allowance for permafrost (Farouki, 1981, Cold Reg. Sci. Technol.)
         !
	 ! 3.1. unfrozen fraction
         xx = (ptn(ji,jg)-(ZeroCelsius-fr_dT/2.)) / fr_dT * poros
         xx = MIN( poros, MAX( 0., xx ) )

         ! 3.2. saturated conductivity
      	 csat = cond_solid**(1.-poros) * cond_ice**(poros-xx) * cond_water**xx
     	
       	 ! 3.3. unsaturated conductivity
     	 pkappa(ji,jg) =(csat - so_cond_dry)*shum_ngrnd_perma(ji,jg) + so_cond_dry

         !
	 ! 3.4. Take into account the snow and soil fractions in the layer
         pkappa(ji,jg) = un / ( zx1(ji,jg) / sn_cond + zx2(ji,jg) / pkappa(ji,jg) )
      END DO            
    ENDDO   

    ! 4. Calculate snow heat capacity and conductivity
     DO ji = 1,kjpindex
         pcapa_snow(ji,:) = snowrho(ji,:) * xci
         pkappa_snow(ji,:) = (ZSNOWTHRMCOND1 + ZSNOWTHRMCOND2*snowrho(ji,:)*snowrho(ji,:)) +      &
               MAX(0.0,(ZSNOWTHRMCOND_AVAP+(ZSNOWTHRMCOND_BVAP/(snowtemp(ji,:)+ &
                ZSNOWTHRMCOND_CVAP)))*(XP00/(pb(ji)*100.)))
    END DO


   
   END SUBROUTINE thermosoilc_getdiff

!! ================================================================================================================================
!! SUBROUTINE   : thermosoilc_getdiff_old_thermix_with_snow
!!
!>\BRIEF          Computes soil heat capacity and conductivity    
!!
!! DESCRIPTION	: Computes soil heat capacity and conductivity
!!                Special case with old snow without soil freezing
!!
!! RECENT CHANGE(S) : None
!! 
!! MAIN OUTPUT VARIABLE(S):
!!                          
!! REFERENCE(S) :
!!
!! FLOWCHART    : None 
!! \n 
!_ ================================================================================================================================


   SUBROUTINE thermosoilc_getdiff_old_thermix_with_snow( kjpindex, snow )


   !! 0. Variables and parameter declaration

    !! 0.1 Input variables
    INTEGER, INTENT(in) :: kjpindex
    REAL,DIMENSION(kjpindex),INTENT (in) :: snow

    !! 0.2 Local variables
    INTEGER					:: ji,jg
    REAL                                 :: snow_h       !! snow_h is the snow height @tex ($m$) @endtex 
    REAL                                 :: zx1, zx2     !! zx1 and zx2 are the layer fraction consisting in snow and soil respectively.

     
    ! Computation of the soil thermal properties; snow properties are also accounted for
    DO ji = 1,kjpindex
      snow_h = snow(ji) / sn_dens

      ! First layer
      IF ( snow_h .GT. zz_coef(1) ) THEN
          pcapa(ji,1) = sn_capa
          pcapa_en(ji,1) = sn_capa
          pkappa(ji,1) = sn_cond
      ELSE IF ( snow_h .GT. zero ) THEN
          pcapa_en(ji,1) = sn_capa
          zx1 = snow_h / zz_coef(1)
          zx2 = ( zz_coef(1) - snow_h) / zz_coef(1)
          pcapa(ji,1) = zx1 * sn_capa + zx2 * so_capa_wet
          pkappa(ji,1) = un / ( zx1 / sn_cond + zx2 / so_cond_wet )
      ELSE
          pcapa(ji,1) = so_capa_dry + shum_ngrnd_perma(ji,1)*(so_capa_wet - so_capa_dry)
          pkappa(ji,1) = so_cond_dry + shum_ngrnd_perma(ji,1)*(so_cond_wet - so_cond_dry)
          pcapa_en(ji,1) = so_capa_dry + shum_ngrnd_perma(ji,1)*(so_capa_wet - so_capa_dry)
      ENDIF

      ! Mid layers
      DO jg = 2, ngrnd - 2
        IF ( snow_h .GT. zz_coef(jg) ) THEN
            pcapa(ji,jg) = sn_capa
            pkappa(ji,jg) = sn_cond
            pcapa_en(ji,jg) = sn_capa
        ELSE IF ( snow_h .GT. zz_coef(jg-1) ) THEN
            zx1 = (snow_h - zz_coef(jg-1)) / (zz_coef(jg) - zz_coef(jg-1))
            zx2 = ( zz_coef(jg) - snow_h) / (zz_coef(jg) - zz_coef(jg-1))
            pcapa(ji,jg) = zx1 * sn_capa + zx2 * so_capa_wet
            pkappa(ji,jg) = un / ( zx1 / sn_cond + zx2 / so_cond_wet )
            pcapa_en(ji,jg) = sn_capa
        ELSE
            pcapa(ji,jg) = so_capa_dry + shum_ngrnd_perma(ji,jg)*(so_capa_wet - so_capa_dry)
            pkappa(ji,jg) = so_cond_dry + shum_ngrnd_perma(ji,jg)*(so_cond_wet - so_cond_dry)
            pcapa_en(ji,jg) = so_capa_dry + shum_ngrnd_perma(ji,jg)*(so_capa_wet - so_capa_dry)
        ENDIF
      ENDDO

      ! Last two layers: These layers can not be filled with snow
      DO jg = ngrnd - 1, ngrnd
         pcapa(ji,jg) = so_capa_dry
         pkappa(ji,jg) = so_cond_dry
         pcapa_en(ji,jg) = so_capa_dry
      END DO
      
    ENDDO ! DO ji = 1,kjpindex


    END SUBROUTINE thermosoilc_getdiff_old_thermix_with_snow



!! ================================================================================================================================
!! SUBROUTINE   : thermosoilc_getdiff_old_thermix_without_snow
!!
!>\BRIEF          Computes soil and snow heat capacity and conductivity    
!!
!! DESCRIPTION	: Calculations of soil and snow thermal properties without effect of freezing. This subroutine is only called
!!                when explicitsnow is activated.
!!                
!!
!! RECENT CHANGE(S) : None
!! 
!! MAIN OUTPUT VARIABLE(S):
!!                          
!! REFERENCE(S) :
!!
!! FLOWCHART    : None 
!! \n 
!_ ================================================================================================================================

    SUBROUTINE thermosoilc_getdiff_old_thermix_without_snow( kjpindex,snowrho, snowtemp, pb )

   !! 0. Variables and parameter declaration

    !! 0.1 Input variables
      INTEGER, INTENT(in) :: kjpindex
      REAL, DIMENSION (kjpindex,nsnow), INTENT(in) :: snowrho  !! Snow density
      REAL, DIMENSION (kjpindex,nsnow), INTENT(in) :: snowtemp !! Snow temperature
      REAL,DIMENSION (kjpindex), INTENT (in)       :: pb       !! Surface presure (hPa)


    !! 0.1 Local variables
      INTEGER						  :: ji,jg

     ! soil
      DO jg = 1,ngrnd
         DO ji = 1,kjpindex
            pkappa(ji,jg) = so_cond_dry + shum_ngrnd_perma(ji,jg)*(so_cond_wet - so_cond_dry)
            pcapa(ji,jg) = so_capa_dry + shum_ngrnd_perma(ji,jg)*(so_capa_wet - so_capa_dry)
            pcapa_en(ji,jg) = so_capa_dry + shum_ngrnd_perma(ji,jg)*(so_capa_wet - so_capa_dry)
         ENDDO
      ENDDO

     ! snow 
     DO ji = 1,kjpindex
         pcapa_snow(ji,:) = snowrho(ji,:) * xci
         pkappa_snow(ji,:) = (ZSNOWTHRMCOND1 + ZSNOWTHRMCOND2*snowrho(ji,:)*snowrho(ji,:)) +      & 
              MAX(0.0,(ZSNOWTHRMCOND_AVAP+(ZSNOWTHRMCOND_BVAP/(snowtemp(ji,:)+ &
               ZSNOWTHRMCOND_CVAP)))*(XP00/(pb(ji)*100.)))
     END DO


   END SUBROUTINE thermosoilc_getdiff_old_thermix_without_snow


!! ================================================================================================================================
!! SUBROUTINE   : thermosoilc_read_reftempfile
!!
!>\BRIEF          
!!
!! DESCRIPTION	: Read file with longterm temperature
!!                
!!
!! RECENT CHANGE(S) : None
!! 
!! MAIN OUTPUT VARIABLE(S): reftemp : Reference temerature
!!                          
!! REFERENCE(S) :
!!
!! FLOWCHART    : None 
!! \n 
!_ ================================================================================================================================

 SUBROUTINE thermosoilc_read_reftempfile(kjpindex,lalo,reftemp)
    
    USE interpweight

    IMPLICIT NONE
    !! 0. Variables and parameter declaration

    !! 0.1 Input variables
    INTEGER, INTENT(in) :: kjpindex
    REAL, DIMENSION(kjpindex,2), INTENT(in) :: lalo
    
    !! 0.2 Output variables
    REAL, DIMENSION(kjpindex, ngrnd), INTENT(out) :: reftemp
    REAL, DIMENSION(kjpindex)                     :: areftemp         !! Availability of data for  the interpolation

    !! 0.3 Local variables
    INTEGER :: ib
    CHARACTER(LEN=80) :: filename
    REAL,DIMENSION(kjpindex)                      :: reftemp_file     !! Horizontal temperature field interpolated from file [C]
    REAL                                          :: vmin, vmax       !! min/max values to use for the
                                                                             !!   renormalization
    CHARACTER(LEN=80)                                    :: variablename     !! Variable to interpolate
                                                                             !!   the file
    CHARACTER(LEN=80)                                    :: lonname, latname !! lon, lat names in input file
    REAL, DIMENSION(:), ALLOCATABLE               :: variabletypevals !! Values for all the types of the variable
                                                                             !!   (variabletypevals(1) = -un, not used)
    CHARACTER(LEN=50)                                    :: fractype         !! method of calculation of fraction
                                                                             !!   'XYKindTime': Input values are kinds 
                                                                             !!     of something with a temporal 
                                                                             !!     evolution on the dx*dy matrix'
    LOGICAL                                              :: nonegative       !! whether negative values should be removed
    CHARACTER(LEN=50)                                    :: maskingtype      !! Type of masking
                                                                             !!   'nomask': no-mask is applied
                                                                             !!   'mbelow': take values below maskvals(1)
                                                                             !!   'mabove': take values above maskvals(1)
                                                                             !!   'msumrange': take values within 2 ranges;
                                                                             !!      maskvals(2) <= SUM(vals(k)) <= maskvals(1)
                                                                             !!      maskvals(1) < SUM(vals(k)) <= maskvals(3)
                                                                             !!       (normalized by maskvals(3))
                                                                             !!   'var': mask values are taken from a 
                                                                             !!     variable inside the file (>0)
    REAL, DIMENSION(3)                            :: maskvals         !! values to use to mask (according to 
                                                                             !!   `maskingtype') 
    CHARACTER(LEN=250)                                   :: namemaskvar      !! name of the variable to use to mask 
    REAL                                          :: reftemp_norefinf
    REAL                                          :: reftemp_default  !! Default value
       
    
    !Config Key   = SOIL_REFTEMP_FILE
    !Config Desc  = File with climatological soil temperature
    !Config If    = READ_REFTEMP
    !Config Def   = reftemp.nc
    !Config Help  = 
    !Config Units = [FILE]
    filename = 'reftemp.nc'
    CALL getin_p('REFTEMP_FILE',filename)

    variablename = 'temperature'

    IF (printlev >= 1) WRITE(numout,*) "thermosoilc_read_reftempfile: Read and interpolate file " &
         // TRIM(filename) //" for variable " //TRIM(variablename)

    ! For this case there are not types/categories. We have 'only' a continuos
    ! field
    ! Assigning values to vmin, vmax

    vmin = 0.
    vmax = 9999.

!   For this file we do not need neightbours!
    neighbours = 0

    !! Variables for interpweight
    ! Type of calculation of cell fractions
    fractype = 'default'
    ! Name of the longitude and latitude in the input file
    lonname = 'nav_lon'
    latname = 'nav_lat'
    ! Default value when no value is get from input file
    reftemp_default = 1.
    ! Reference value when no value is get from input file
    reftemp_norefinf = 1.
    ! Should negative values be set to zero from input file?
    nonegative = .FALSE.
    ! Type of mask to apply to the input data (see header for more details)
    maskingtype = 'nomask'
    ! Values to use for the masking (here not used)
    maskvals = (/ undef_sechiba, undef_sechiba, undef_sechiba /)
    ! Name of the variable with the values for the mask in the input file (only if maskkingtype='var') (here not used)
    namemaskvar = ''

    CALL interpweight_2Dcont(kjpindex, 0, 0, lalo, resolution, neighbours,                            &
      contfrac, filename, variablename, lonname, latname, vmin, vmax, nonegative, maskingtype,        &
      maskvals, namemaskvar, -1, fractype, reftemp_default, reftemp_norefinf,                         &
      reftemp_file, areftemp)
    IF (printlev >= 5) WRITE(numout,*)'  thermosoilc_read_reftempfile after interpweight2D_cont'

    ! Copy reftemp_file temperature to all ground levels and transform into Kelvin
    DO ib=1, kjpindex
      reftemp(ib, :) = reftemp_file(ib)+ZeroCelsius
    END DO

    ! Write diagnostics
    CALL xios_orchidee_send_field("areftemp",areftemp)
       
  END SUBROUTINE thermosoilc_read_reftempfile
  

END MODULE thermosoilc
