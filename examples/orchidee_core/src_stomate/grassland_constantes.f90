! =================================================================================================================================
! MODULE       : grassland_constantes
!
! CONTACT      : orchidee-help _at_ ipsl.jussieu.fr
!
! LICENCE      : IPSL (2006)
! This software is governed by the CeCILL licence see
! ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF       This module defined all constantes used in
!! grassland management module
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
MODULE grassland_constantes

  USE constantes

  LOGICAL, PARAMETER      :: blabla_pasim    = .FALSE. 
  LOGICAL, PARAMETER      :: DEBUG_anne      = .FALSE.
  LOGICAL, PARAMETER      :: BIG_DEBUG_anne  = .FALSE.
  INTEGER, PARAMETER      :: taille      = SELECTED_REAL_KIND(8)
  INTEGER, PARAMETER      :: taille_stdn = SELECTED_REAL_KIND(6,30)
  INTEGER, PARAMETER      :: nmaxprinc   = 8760
  INTEGER, PARAMETER      :: maxvar      = 100
  INTEGER, PARAMETER      :: NCUT        = 10
  INTEGER, PARAMETER      :: NSOILMAX    = 6
  INTEGER, PARAMETER      :: NFERT       = 10
  INTEGER, PARAMETER      :: NCANOPY     = 10
  INTEGER, PARAMETER      :: NGMEAN      = 10
  INTEGER, PARAMETER      :: NSTOCKING   = 10
  INTEGER, PARAMETER      :: NGAUSS      = 5
  INTEGER, PARAMETER      :: igrass      = 10
  LOGICAL, PARAMETER      :: Yearly_1snowzero          = .FALSE.
  REAL, PARAMETER :: minimum_vital             = 1e-30
  REAL, PARAMETER :: DEFAULTMaxt               = 0.02
  REAL, PARAMETER :: lambda                    = 2.49
  REAL, PARAMETER :: gamma                     = 0.065
  REAL, PARAMETER :: rho                       = 1.23
  REAL, PARAMETER :: alpha350                  = 0.145
  REAL, PARAMETER :: devsecond                 = 0.77
  REAL, PARAMETER :: sigma                     = 5.67e-8
  LOGICAL     , PARAMETER :: soilnh4pmobile            = .TRUE.
  REAL, PARAMETER :: nh3mol2nh3ug              = 17.0e6
  REAL, PARAMETER :: qm2mol                    = 44.6
  REAL, PARAMETER :: nn2omol2kg                = 0.028
  REAL, PARAMETER :: lf                        = 333750.0
  REAL, PARAMETER :: cp                        = 1.01e3
  REAL, PARAMETER :: rwt                       = 0
  REAL, PARAMETER :: dailysunshine             = 0
  REAL, PARAMETER :: pa2atm                    = 101325.0
  ! konvertierung von umol co2 nach kg c (kg c/umol co2)
  REAL, PARAMETER :: cumol2ckg                 = 12.0e-9 
  ! convertion de mol n en kg n 
  REAL, PARAMETER :: nmol2kg                   = 0.014 
  ! konvertierung von m nach mm (mm/m)  
  REAL, PARAMETER :: m2mm                      = 1000.0  
  REAL, PARAMETER :: highresoutput             = 0
  REAL, PARAMETER :: tfr                       = 273.15
  REAL, PARAMETER :: d0                        = 2.15e6
  REAL, PARAMETER :: nstructfix                = 0
  ! dichte von wasser @293.16k (kg h2o/m**3)
  REAL, PARAMETER :: rhow                      = 998.206  
  REAL, PARAMETER :: cw                        = 4.18e6
  REAL, PARAMETER :: kg2mg                     = 1.0e6
  ! anzahl sekunden pro tag (s/d)
  REAL, PARAMETER :: d2s                       = 86400   
  LOGICAL   , PARAMETER :: unanaerob                 = .TRUE.
  REAL, PARAMETER :: karman                    = 0.41
  REAL, PARAMETER :: dailyweatherdata          = 0
  REAL, PARAMETER :: rgas                      = 8.3143
  REAL, PARAMETER :: solar                     = 1370
  REAL, PARAMETER :: infinit                   = 1.0e20
  REAL, PARAMETER :: grav                      = 9.81
  ! ! konvertierung von l nach m**3 (m**3/l)
  REAL, PARAMETER :: l2m3                      = 1.0e-3 
  REAL, PARAMETER :: ug2kg                     = 1.0e-9
  REAL, PARAMETER :: cmol2kg                   = 0.012
  REAL, PARAMETER :: mw                        = 0.018
  REAL, PARAMETER :: mc                        = 28.5
  REAL, PARAMETER :: fcsh                      = 0.39
  REAL, PARAMETER :: ssta                      = 6.6
  REAL, PARAMETER :: fwnapo                    = 0.05
  REAL, PARAMETER :: mn                        = 62.0
  REAL, PARAMETER :: slam                      = 33.5   ! = 33.5
  REAL, PARAMETER :: fcr                       = 0.50
  REAL, PARAMETER :: tconstraintmax            = 10000.0
  REAL, PARAMETER :: nage                      = 4.0
  ! parameter for the calculation of enteric methane emission (kg CH4 (kg life weight)-1 d-1) 
  REAL, PARAMETER :: aCH4                      = 0.0002867   ! -0.0087  @equation constantes::aCH4
  ! parameter for the calculation of enteric methane emission (kg CH4 (kg life weight)-1 d-1)
  REAL, PARAMETER :: bCH4                      = 0.000045    ! 0.0541   @equation constantes::bCH4
  ! parameter for the calculation of enteric methane emission !!! 
  REAL, PARAMETER :: CH4toC                    = 0.75
  !        @equation constantes::CH4toC
  ! parameter for transform the biomass(gC/m2) to wsh (kgDM/m2) 
  ! for grass C ratio in Dry Matter usually 0.45 for trees usually 0.5
  REAL, PARAMETER :: CtoDM                    = 0.45  
  REAL, PARAMETER :: zeta                      = 10.0
  !facteur de conversion de MS a MO
  REAL, PARAMETER :: dm2om                     = 0.9    
  REAL, PARAMETER :: vcmaxadap                 = 1.0
  REAL, PARAMETER :: absorvl                   = 0.85
  ! fraction of net energy converted to metabolizable energy (-)
  REAL(Taille), PARAMETER :: k_CH4                     = 0.6   
  ! flag d'activation de l'effet des temperatures elevees sur l'ingestion    
  LOGICAL     , PARAMETER :: f_temperature_DMI         = .TRUE.
  ! TRUE : CH4 is calculated from N.Vuichard equation    
  LOGICAL     , PARAMETER :: f_CH4_methode             = .FALSE.   

  REAL, PARAMETER :: devear                    = 0.52
  REAL, PARAMETER :: nanimaltotmax             = 0.001
  REAL, PARAMETER :: tbase         = 278.15 !278.0
  REAL, PARAMETER :: trep          = 278.15 !278.0
  REAL, PARAMETER :: tasumrep      = 225.0

  ! variables for version 3.7 
  REAL, PARAMETER ::  l_min = 0.1
  REAL, PARAMETER ::  age_init_1 = 0
  REAL, PARAMETER ::  age_init_2 = 0
  REAL, PARAMETER ::  age_init_3 = 0
  REAL, PARAMETER ::  age_init_4 = 0

  REAL, PARAMETER :: fligninstructinit = 0.54
  REAL, PARAMETER :: fligninresidue    = 0.25

  REAL, PARAMETER :: devstagemin    = 0.6
  REAL, PARAMETER :: tgrowthmin     = 45 !30.0
  REAL, PARAMETER :: tgrowthmax     = 70.0
  REAL, PARAMETER :: misval         = -99.999
  REAL, PARAMETER :: gmeansloperel  = -0.002
  REAL, PARAMETER :: devstocking    = 0.001
  ! original tseasonendmin  = 250.0 is for France grassland
  ! JC change it to 0 for global application
  REAL, PARAMETER :: tseasonendmin  = -1.0

  ! constante for soil physics module
  REAL, PARAMETER    :: Temp_1csoil           = 2.4E6
  REAL, PARAMETER    :: Temp_1eta             = 0.6
  REAL, PARAMETER    :: sswmax                = 10000.0
  REAL, PARAMETER    :: tspring               = 100.0
  REAL, PARAMETER    :: tautomn               = 300.0

  ! constante for principal module
  REAL, PARAMETER :: psifc             = -612.
  REAL, PARAMETER :: psipwp            = -152905.

  ! constante for fertilisation module
  REAL, PARAMETER :: tapplmist        = 1.0
  REAL, PARAMETER :: fmistcmetabolic  = 0.7
  REAL, PARAMETER :: tapplvg          = 1.0
  REAL, PARAMETER :: asunshine        = 0.19
  REAL, PARAMETER :: fvgcmetabolic    = 0.8
  REAL, PARAMETER :: bsunshine        = 0.56
  REAL, PARAMETER :: tapplka          = 1.0
  REAL, PARAMETER :: c2nmist          = 15.0
  REAL, PARAMETER :: fvgnurine        = 0.61
  REAL, PARAMETER :: fkacmetabolic    = 0.9
  REAL, PARAMETER :: c2nkot           = 10.0
  REAL, PARAMETER :: fkanurine        = 0.9
  ! si DMI/DMIpot < intake_tolerance alors on sort les animaux ou on les complemente
  REAL, SAVE      :: intake_tolerance    = 0.80
  ! Maximum quantity of forage or concentrate to supplement animals when auto-supplementation (kg)
  REAL, SAVE      :: q_max_complement    = 6  

  INTEGER, SAVE :: f_saturant     = 0.
  INTEGER, SAVE :: f_nonlimitant  = 0.
  INTEGER, SAVE :: f_autogestion  = 0.     !080110 AIG fin
  INTEGER, SAVE :: f_complementation  = 0. !11 janvier 2009 AIG
  ! 0 : pas de fertilisation
  ! 1 : fertilisation geree par le modele pour un niveau de satisfaction des besoins en N
  INTEGER, SAVE :: f_fertilization = 0
  INTEGER, SAVE :: f_postauto = 0
  !071129 AIG
  ! constante for auto management modifications by nicolas vuichard
  !  INTEGER     , PARAMETER :: n_somitermin        = 2
  !  INTEGER     , PARAMETER :: n_somitermax        = 400
  !  REAL(Taille), PARAMETER :: nsatur_somerrormax  = 0.01
  ! Minimum threshold of Wshtot, below which animals are kept indoors (kg.m-2),
  ! either used when management is optimized by the model or for usual runs.
  REAL, PARAMETER :: min_grazing         = 0.08 
  !  INTEGER     , PARAMETER :: cte_nb_high_out     = 1
  !  INTEGER     , PARAMETER :: cte_nb_low_out      = 1
  !  INTEGER     , PARAMETER :: nb_last_year        = 1
  !20/03/2009 AIG & MG
  INTEGER, SAVE :: n_out          = 3 
  !Is the number of output for autogestion mangement run in import Yield file. 
  ! Dim1 : Fraction of grazed aera (F)   Dim2: Ratio F/(1-F)  Dim3: number of grazing days

END MODULE grassland_constantes
