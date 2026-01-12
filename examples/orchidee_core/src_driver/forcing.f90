!  ==============================================================================================================================\n
!  MODULE 	: 
! 
!  CONTACT      : orchidee-help _at_ listes.ipsl.fr
!
!  LICENCE      : IPSL (2006)
!  This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF        This module contains subroutines to manage the forcing list
!!
!!
!!\n DESCRIPTION  : This module contains subroutines for reading the forcing file for the dim2_driver.
!!                  Following subroutines are public and called from dim2_driver :
!!                    - forcing_info : Open the forcing file and return information about the grid in the forcing file. 
!!                                     Prepare for a zoom if needed. 
!!                                     Initialization of parallelization related to the grid. 
!!                    - forcing_read : Return the forcing data for the current time step of the model. The forcing file will
!!                                     be read if it has not already been done for the current time-step in the forcing file.
!!                    - forcing_grid : Calculate the longitudes and latitudes of the model grid.
!! 
!! RECENT CHANGE(S): None 
!! 
!! REFERENCE(S) : None
!!   
!! SVN     :
!! $HeadURL: svn://forge.ipsl.jussieu.fr/orchidee/branches/ORCHIDEE-CN-P/ORCHIDEE/src_driver/readdim2.f90 $ 
!! $Date: 2019-01-16 11:11:16 +0100 (mer. 16 janv. 2019) $
!! $Revision: 5700 $
!! \n
!_ ================================================================================================================================

MODULE forcing 

   USE string

   USE ioipsl_para
   USE ioipsl_lists

   IMPLICIT NONE

   PRIVATE
   PUBLIC  :: T_FORCING_LIST

   TYPE T_FORCING_LIST
       CHARACTER(LEN=200), DIMENSION(:), ALLOCATABLE :: file_list
       LOGICAL :: use_list ! Forcings are red from file
       LOGICAL :: has_cycle ! start the list again once it's at the end of the list
       INTEGER :: selected_file ! Index to point the file_list filename

   CONTAINS
       PROCEDURE :: get
       PROCEDURE :: total
       PROCEDURE :: from_list
       PROCEDURE :: next
   END TYPE T_FORCING_LIST

   INTERFACE T_FORCING_LIST
     PROCEDURE :: init_t_forcing_list
   END INTERFACE 

CONTAINS

!! ==============================================================================================================================\n
!! SUBROUTINE 	: 
!!
!>\BRIEF        
!!
!!\n DESCRIPTION : 
!!
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): 
!!
!! REFERENCE(S)	: 
!!
!_ ================================================================================================================================
  FUNCTION init_t_forcing_list(forcing_file_flag, forcing_cycle) RESULT ( this )
    TYPE(t_forcing_list) :: this
    CHARACTER(LEN=*) :: forcing_file_flag     !! Name of the file to be opened
    LOGICAL, INTENT(in) :: forcing_cycle

    CHARACTER(LEN=:), ALLOCATABLE :: filename     !! Name of the file to be opened
!    CHARACTER(LEN=300), DIMENSION(:), ALLOCATABLE :: space_list     !! Name of the file to be opened
    LOGICAL :: is_filelist, is_ext_found
    INTEGER :: ier, nfiles, fidx, file_counter, elem
    TYPE(t_file_def_list) :: flist

    ! Init
    this%use_list = .FALSE.
    nfiles = 1
    this%selected_file = 1
    this%has_cycle = forcing_cycle 

    ! From run.def
    filename = TRIM(forcing_file_flag)
!    WRITE(*,*) "init_t_forcing_list:: filename=", forcing_file_flag
!    WRITE(*,*) "init_t_forcing_list:: filename=", filename

    ! Discover which type: simple filename or a forcing file list
    is_filelist = string_endswith(filename, FILE_LIST_EXT)
!    WRITE(*,*) "init_t_forcing_list:: is_filelist=", is_filelist

    IF (is_filelist) THEN
        flist = t_file_def_list(filename)
        this%use_list = .TRUE.
        nfiles = flist%get_num_files()

        ALLOCATE(this%file_list(nfiles), stat=ier)
        IF (ier /= 0) CALL ipslerr_p(3, 'init_t_forcing_list','Allocation memory error', &
                        'Error code:', 'TODO: to add later')

        ! Copy files into array
        DO fidx = 1, nfiles
          this%file_list(fidx) = flist%get_file(fidx)
        ENDDO

    ELSE
!        WRITE(*,*) "init_t_forcing_list:: filename=", filename
        ! how many files are defined in the key?
        nfiles = countsubstring(filename, ".nc")
!        WRITE(*,*) "init_t_forcing_list:: nfiles=", nfiles

        ALLOCATE(this%file_list(nfiles), stat=ier)
        IF (ier /= 0) CALL ipslerr_p(3, 'init_t_forcing_list','Allocation memory error', &
                        'Error code:', 'TODO: to add later')

        IF (nfiles .EQ. 1) THEN ! just one
            WRITE(*,*) "init_t_forcing_list:: single forcing detected"
            this%file_list(1)=filename
        ELSE IF (nfiles .LT. 1) THEN ! No one
            CALL ipslerr_p(3, 'init_t_forcing_list','No netcdf file found', &
                        'Found:', filename)
        ELSE
            CALL ipslerr_p(3, 'init_t_forcing_list','More than 1 netcdf file found', &
                        'Option not available', "Found: "//filename)

        ENDIF
!        ELSE ! More than 1
!            WRITE(*,*) "init_t_forcing_list:: multiple forcings detected"
!            this%use_list = .TRUE.
!
!            ! separate 
!            CALL string_split(filename, " ", space_list, ier)
!
!            ! filter only keep filenames
!            file_counter = 1
 !           DO elem=1, SIZE(space_list)
 !!             is_ext_found = countsubstring(space_list(elem), ".nc")
 !             IF (is_ext_found) THEN
  !!              this%file_list(file_counter) = TRIM(space_list(elem))
  !              file_counter = file_counter + 1
  !            ENDIF
  !          ENDDO

!            CALL string_split(filename, ".nc ", this%file_list, ier)
!        ENDIF
    ENDIF

!    WRITE(*,*) "init_t_forcing_list:: forcing files found=", this%file_list
      
  END FUNCTION init_t_forcing_list

!! ==============================================================================================================================\n
!! SUBROUTINE 	: 
!!
!>\BRIEF        
!!
!!\n DESCRIPTION : 
!!
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): 
!!
!! REFERENCE(S)	: 
!!
!_ ================================================================================================================================
  SUBROUTINE next(this) 
    CLASS(t_forcing_list) :: this

    this%selected_file = this%selected_file + 1
    IF (this%selected_file .GT. SIZE(this%file_list)) THEN
      IF (this%has_cycle) CALL ipslerr_p(3,'forcing::next','There are no more forcing files to read','','')
      this%selected_file = 1
    ENDIF
  END SUBROUTINE next

!! ==============================================================================================================================\n
!! SUBROUTINE 	: 
!!
!>\BRIEF        
!!
!!\n DESCRIPTION : 
!!
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): 
!!
!! REFERENCE(S)	: 
!!
!_ ================================================================================================================================
  FUNCTION get(this, idx) RESULT ( filename )
    CLASS(t_forcing_list) :: this
    INTEGER :: idx

    CHARACTER(LEN=:), ALLOCATABLE :: filename

    filename = this%file_list(idx)
  END FUNCTION get

!! ==============================================================================================================================\n
!! SUBROUTINE 	: 
!!
!>\BRIEF        
!!
!!\n DESCRIPTION : 
!!
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): 
!!
!! REFERENCE(S)	: 
!!
!_ ================================================================================================================================
  FUNCTION total(this) RESULT ( total_files )
    CLASS(t_forcing_list) :: this

    INTEGER :: total_files

    total_files = SIZE(this%file_list) 
  END FUNCTION total

!! ==============================================================================================================================\n
!! SUBROUTINE 	: 
!!
!>\BRIEF        
!!
!!\n DESCRIPTION : 
!!
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): 
!!
!! REFERENCE(S)	: 
!!
!_ ================================================================================================================================
  FUNCTION from_list(this) RESULT ( has_list )
    CLASS(t_forcing_list) :: this

    LOGICAL :: has_list

    has_list = this%use_list
  END FUNCTION from_list

END MODULE forcing
