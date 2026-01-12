! =================================================================================================================================
! MODULE       : vertical_soil_var
!
! CONTACT      : orchidee-help _at_ listes.ipsl.fr
!
! LICENCE      : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF        
!!
!!\n DESCRIPTION: 
!!                
!! RECENT CHANGE(S):
!!
!! REFERENCE(S)	: 
!!
!! SVN          :
!! $HeadURL: svn://forge.ipsl.jussieu.fr/orchidee/branches/ORCHIDEE-CN-P/ORCHIDEE/src_parameters/vertical_soil_var.f90 $
!! $Date: 2017-10-26 14:32:36 +0200 (四, 2017-10-26) $
!! $Revision: 4717 $
!! \n
!_ ================================================================================================================================

MODULE vertical_soil_var

  USE defprec

  IMPLICIT NONE
  PUBLIC

  !! Dimensioning parameters
  INTEGER, SAVE      :: ngrnd     !! Number of soil layer for thermo (unitless)
!$OMP THREADPRIVATE(ngrnd)
  INTEGER, SAVE      :: nslm      !! Number of levels in CWRR (unitless)
!$OMP THREADPRIVATE(nslm)
  REAL, SAVE         :: zmaxh     !! Maximum depth of soil reservoir in hydrol (m). Old name dpu_max or depth_Wmax
!$OMP THREADPRIVATE(zmaxh)
  REAL, SAVE         :: zmaxt     !! Maximum depth of the soil thermodynamics (m)
!$OMP THREADPRIVATE(zmaxt)

  !! Variables defining the vertical layering in soil moisture and temperature
  REAL, SAVE, ALLOCATABLE, DIMENSION(:) :: znt          !! Depth of nodes for thermal (m) 
!$OMP THREADPRIVATE(znt)
  REAL, SAVE, ALLOCATABLE, DIMENSION(:) :: znh          !! Depth of nodes for hydrology (m)
!$OMP THREADPRIVATE(znh)
  REAL, SAVE, ALLOCATABLE, DIMENSION(:) :: dnh          !! Distance between the current node and the one above for hydrology (m)
!$OMP THREADPRIVATE(dnh)
  REAL, SAVE, ALLOCATABLE, DIMENSION(:) :: dlh          !! Soil layer thickness for hydrology (m) 
!$OMP THREADPRIVATE(dlh)
  REAL, SAVE, ALLOCATABLE, DIMENSION(:) :: dlt          !! Soil layer thickness for thermal (m)
!$OMP THREADPRIVATE(dlt)
  REAL, SAVE, ALLOCATABLE, DIMENSION(:) :: zlh          !! Depth of lower layer-interface for hydrology (m)
!$OMP THREADPRIVATE(zlh)
  REAL, SAVE, ALLOCATABLE, DIMENSION(:) :: zlt          !! Depth of lower layer-interface for thermal (m)
!$OMP THREADPRIVATE(zlt)

  REAL,ALLOCATABLE, DIMENSION(:),SAVE :: diaglev        !! The lower limit of the layer on which soil moisture
                                                               !! (relative) and temperature are going to be diagnosed.
                                                               !! These variables are made for transfering the information
                                                               !! to the biogeophyical processes modelled in STOMATE. 
!$OMP THREADPRIVATE(diaglev)

END MODULE vertical_soil_var
