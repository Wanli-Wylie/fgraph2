! ==============================================================================================================================
! MODULE   : ioipsl_lists
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
!!
!! RECENT CHANGE(S): None
!!
!! REFERENCES(S)    : None
!!
!! SVN              :
!! $HeadURL: svn://forge.ipsl.jussieu.fr/orchidee/branches/ORCHIDEE-CN-P/ORCHIDEE/src_parallel/ioipsl_para.f90 $ 
!! $Date: 2017-10-26 14:32:36 +0200 (Thu, 26 Oct 2017) $
!! $Revision: 4717 $
!! \n
!_ ================================================================================================================================

MODULE ioipsl_lists

  USE ioipsl
!  USE mod_orchidee_para_var
!  USE mod_orchidee_transfert_para
!-
  IMPLICIT NONE

  PUBLIC :: FILE_LIST_EXT, t_file_def_list

  CHARACTER(LEN=5), PARAMETER :: FILE_LIST_EXT = ".list" ! File list extension. E.g: something.list

  TYPE t_file_def_list
      CHARACTER(LEN=:), ALLOCATABLE :: filename
      INTEGER :: file_id
      
      CHARACTER(LEN=200), DIMENSION(:), ALLOCATABLE :: file_list
  CONTAINS
      PROCEDURE :: read_filelist
      PROCEDURE :: get_file
      PROCEDURE :: get_num_files
  END TYPE t_file_def_list

  INTERFACE t_file_def_list
    PROCEDURE :: init_file_def_list
  END INTERFACE t_file_def_list

  LOGICAL, PARAMETER :: is_debug = .FALSE.

CONTAINS

  !!  =============================================================================================================================
  !! SUBROUTINE:  init_file_def_list
  !!
  !>\BRIEF	    Type initialization by providing a filepath 
  !!
  !! DESCRIPTION: Initialize the type by providing a filepath. It will read the
  !!            file list and load the filenames.
  !!
  !! \n
  !_ ==============================================================================================================================
  FUNCTION init_file_def_list(f_name) RESULT ( this )
    TYPE(t_file_def_list) :: this
    CHARACTER(LEN=*), INTENT(in) :: f_name

    LOGICAL :: file_exists
    INTEGER :: nb_files, ier

    this%filename = TRIM(f_name)

    ! Does file exists? If NO, crash
    INQUIRE(FILE=this%filename, EXIST=file_exists)
    IF (.NOT. file_exists) THEN
        CALL ipslerr(3, 'init_file_def_lists','Error','File not found:',this%filename)
    ENDIF

    ! Discover how many files are in the list
    nb_files = this%get_num_files()
    ALLOCATE(this%file_list(nb_files), stat=ier)
    IF (ier/= 0) CALL ipslerr(3,'init_file_def_list','Allocation memory error','Error found:','TODO: add here error code')

    ! Load them
    CALL this%read_filelist()

  END FUNCTION init_file_def_list

  !!  =============================================================================================================================
  !! SUBROUTINE:  get_file
  !!
  !>\BRIEF	    Get a file name from an index
  !!
  !! DESCRIPTION: Provides a file name from the list
  !!
  !! \n
  !_ ==============================================================================================================================
  FUNCTION get_file(this, idx) RESULT (f_name)
    CLASS(t_file_def_list), INTENT(inout) :: this
    INTEGER, INTENT(in) :: idx

    CHARACTER(LEN=:), ALLOCATABLE :: f_name

    IF (idx .GT. LEN(this%file_list)) THEN
        CALL ipslerr(3,"get_file","Index not valid","The selected index is bigger than the existing number of files ","")
    ENDIF

    f_name = this%file_list(idx)

  END FUNCTION get_file

  !!  =============================================================================================================================
  !! SUBROUTINE:  read_filelist
  !!
  !>\BRIEF	    Read the list file and load data
  !!
  !! DESCRIPTION: Read the list file and load data
  !!
  !! \n
  !_ ==============================================================================================================================
  SUBROUTINE read_filelist(this) 
!-
    CLASS(t_file_def_list) :: this
    INTEGER :: n_files
!-
    INTEGER, PARAMETER :: fid = 22 ! file id
    INTEGER :: eof, io_err, nb_lastkey, len_str, ier, nb_files
    CHARACTER(LEN=:), ALLOCATABLE :: READ_str
!-
    nb_files = 0
    eof = 0
    nb_lastkey = 0
    n_files = 0

    IF (.NOT. ALLOCATED(this%file_list)) THEN
      CALL ipslerr(3,'read_filelist','this%file_list is not allocated','','')
    ENDIF
!-
    OPEN (UNIT=fid,FILE=this%filename,STATUS="OLD",IOSTAT=io_err)
    IF (io_err /= 0) THEN
      CALL ipslerr (3,'read_filelist', &
         &  'Could not open file:', this%filename,'')
      RETURN
    ENDIF
!-
    DO WHILE (eof /= 1)
!---
      CALL ioipsl_lists_skipafew (fid,READ_str,eof,nb_lastkey)
      len_str = LEN_TRIM(READ_str)
      IF (len_str > 0) THEN
          n_files = n_files + 1
          this%file_list(n_files) = TRIM(READ_str)
      ENDIF
      IF (is_debug) WRITE(*,*) "read_filelist:: len_str=", len_str, READ_str
!---
    ENDDO 

    CLOSE(fid)

    RETURN
!-
  END SUBROUTINE read_filelist

  !!  =============================================================================================================================
  !! SUBROUTINE: get_num_files
  !!
  !>\BRIEF	    Number of files to read
  !!
  !! DESCRIPTION: Number of files to read
  !!
  !! \n
  !_ ==============================================================================================================================
  FUNCTION get_num_files(this) RESULT (n_files)
!-
    CLASS(t_file_def_list) :: this
    INTEGER :: n_files
!-
    INTEGER, PARAMETER :: fid = 22 ! file id
    INTEGER :: eof, io_err, nb_lastkey, len_str
    CHARACTER(LEN=:), ALLOCATABLE :: READ_str
!-
    n_files = 0
    eof = 0
    nb_lastkey = 0
!-
    ! File is alread loaded, no need to read the whole file again
    IF ( ALLOCATED(this%file_list )) THEN
        n_files = SIZE(this%file_list)
        RETURN
    ENDIF
!-
    OPEN (UNIT=fid,FILE=this%filename,STATUS="OLD",IOSTAT=io_err)
    IF (io_err /= 0) THEN
      CALL ipslerr (3,'get_num_files', &
         &  'Could not open file:', this%filename,'')
      RETURN
    ENDIF
!-
    DO WHILE (eof /= 1)
!---
      CALL ioipsl_lists_skipafew (fid,READ_str,eof,nb_lastkey)
      len_str = LEN_TRIM(READ_str)
      IF (len_str > 0) n_files = n_files + 1
      IF (is_debug) WRITE(*,*) "get_num_files:: len_str=", len_str, READ_str
!---
    ENDDO 

    CLOSE(fid)

    RETURN
!-
  END FUNCTION get_num_files
!-
! FROM ioispl
!-
SUBROUTINE ioipsl_lists_readline(unitf, out_string, is_eof)
!---------------------------------------------------------------------
  USE ISO_FORTRAN_ENV,ONLY : IOSTAT_EOR,IOSTAT_END
!-
  IMPLICIT NONE
!-
  INTEGER, PARAMETER :: CHARLEN = 100     ! buffer size
  INTEGER, INTENT(in) :: unitf
  INTEGER, INTENT(out) :: is_eof
  CHARACTER(LEN=:),INTENT(out),ALLOCATABLE :: out_string
!-
  CHARACTER(LEN=CHARLEN)  :: dummy
!-
  CHARACTER(LEN=:), ALLOCATABLE :: buff1  ! buffer
  INTEGER :: ioerr                  ! error code
  INTEGER :: readlength                   ! number of chars read from file 
  LOGICAL :: is_eol, is_first_ite         ! end of line? 

  is_eof = 0
  is_eol = .FALSE.
  buff1 = ""

  DO WHILE (.NOT. is_eol)
!-
      dummy = ""
      READ (UNIT=unitf,FMT='(A)', ADVANCE='NO', SIZE=readlength,ERR=9998,END=7778,IOSTAT=ioerr) dummy
      IF ((ioerr==IOSTAT_EOR).OR.(ioerr==IOSTAT_END)) ioerr = 0
!-
      ! keep looping if line is commented
      dummy = TRIM(ADJUSTL(dummy))
!-
      ! is end of line?
      is_eol = (readlength .LT. CHARLEN)
!-
      ! merge with previous buffer
      buff1 = TRIM(buff1)//TRIM(dummy)
  ENDDO
!-
  out_string=TRIM(buff1)
!-
  RETURN
!-
9998 CONTINUE
  CALL ipslerr (3,'ioipsl_lists_readline','Error while reading file',' ',' ')
!-
7778 CONTINUE
  out_string = TRIM(dummy)
  is_eof = 1
    
END SUBROUTINE ioipsl_lists_readline
!-
! ioipsl_lists_skipafew: reads  
! FROM ioispl
!-
SUBROUTINE ioipsl_lists_skipafew (unit,out_string,is_eof,nb_lastkey)
!---------------------------------------------------------------------
  USE ISO_FORTRAN_ENV,ONLY : IOSTAT_EOR,IOSTAT_END
!-
  IMPLICIT NONE
!-
  INTEGER :: unit,is_eof,nb_lastkey
  CHARACTER(LEN=:),INTENT(out),ALLOCATABLE :: out_string
!-
  CHARACTER(LEN=1) :: first
  CHARACTER(LEN=:), ALLOCATABLE :: dummy
!---------------------------------------------------------------------
  first=COMMENT_TAG
  is_eof = 0
  dummy = ""
!-
! Loop until a non commented line is found
  DO WHILE (first == COMMENT_TAG .AND. is_eof == 0)
!-
     CALL ioipsl_lists_readline(unit, dummy, is_eof)
!-
!    Is first char a comment? #
     IF (LEN(dummy) > 0) THEN
         first=dummy(1:1)
         IF (first == COMMENT_TAG) THEN
           nb_lastkey = 0
         ENDIF
     ENDIF
!-
  ENDDO
!-
  CALL nocomment(dummy)
  out_string = TRIM(dummy)
!----------------------------
END SUBROUTINE ioipsl_lists_skipafew

END MODULE ioipsl_lists
