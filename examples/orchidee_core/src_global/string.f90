!  ==============================================================================================================================\n
!  MODULE 	: 
! 
!  CONTACT      : orchidee-help _at_ listes.ipsl.fr
!
!  LICENCE      : IPSL (2006)
!  This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF        
!!
!!
!!\n DESCRIPTION  : 
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

MODULE string 

   IMPLICIT NONE

   PRIVATE
   PUBLIC  :: string_endswith, countsubstring, string_split

CONTAINS

!! ==============================================================================================================================\n
!! SUBROUTINE 	: 
!!
!>\BRIEF        
!!
!!\n DESCRIPTION : string and with 'ending'
!!              true : in_str ends with ENDING
!!              false: in_str does not end with ENDING
!!
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): 
!!
!! REFERENCE(S)	: 
!!
!_ ================================================================================================================================
  FUNCTION string_endswith(in_str, ending) RESULT (is_ending)
    CHARACTER(LEN=*), INTENT(in) :: in_str
    CHARACTER(LEN=*), INTENT(in) :: ending

    CHARACTER(LEN=:), ALLOCATABLE :: clean_in_str
    INTEGER :: idx

    LOGICAL :: is_ending

    is_ending = .FALSE.

    clean_in_str =TRIM(in_str) 
    idx = INDEX(clean_in_str, ending)
    !WRITE(*,*) "string_endswith:: idx=", idx

    ! substring found
    IF (idx .GT. 0) THEN
        ! is it at the back?
        !WRITE(*,*) "string_endswith:: idx+len(ending)=", idx+LEN(ending)-1
        !WRITE(*,*) "string_endswith:: LEN(clean_in_str)=", LEN(clean_in_str)
        IF (idx + LEN(ending) - 1 .EQ. LEN(clean_in_str) ) THEN
          ! yes
          is_ending = .TRUE.
        ENDIF
    ENDIF ! (idx .EQ. 0) 

  END FUNCTION string_endswith


!! ==============================================================================================================================\n
!! SUBROUTINE 	: 
!!
!>\BRIEF        
!!
!!\n DESCRIPTION : string and with 'ending'
!!              https://www.rosettacode.org/wiki/Count_occurrences_of_a_substring#Fortran
!!
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): 
!!
!! REFERENCE(S)	: 
!!
!_ ================================================================================================================================
  FUNCTION countsubstring(s1, s2) RESULT(c)
    character(*), intent(in) :: s1, s2
    integer :: c, p, posn
           
    c = 0
    if(len(s2) == 0) return
    p = 1
    do 
      posn = index(s1(p:), s2)
!      WRITE(*,*) "countsubstring:: s1=", s1(p:)
!      WRITE(*,*) "countsubstring:: pos=", posn
      if(posn == 0) exit
      c = c + 1
      p = p + posn + len(s2)-1
!      WRITE(*,*) "countsubstring:: next p=", s1(p:)
    end do
  END FUNCTION countsubstring

!! ==============================================================================================================================\n
!! SUBROUTINE 	: 
!!
!>\BRIEF        
!!
!!\n DESCRIPTION : string and with 'ending'
!!              https://www.rosettacode.org/wiki/Count_occurrences_of_a_substring#Fortran
!!
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): 
!!
!! REFERENCE(S)	: 
!!
!_ ================================================================================================================================
  SUBROUTINE string_split(str, token, out_str, stat) 
    CHARACTER(*), INTENT(IN) :: str
    CHARACTER(*), INTENT(IN) :: token

    CHARACTER(LEN=300), ALLOCATABLE, DIMENSION(:), INTENT(out) :: out_str
    INTEGER, INTENT(out) :: stat
    INTEGER :: c, p, posn, num_tok, startpos, endpos, ier
    
    num_tok = countsubstring(str,token) + 1
!    WRITE(*,*) "string_split:: found parts=", num_tok

    stat = 0

    IF (.NOT. ALLOCATED(out_str) ) THEN
        ALLOCATE(out_str(num_tok), stat=ier)
        !CALL iplserr(3,"string_split",'Memory allocation error for out_str','Error code=',ier)
        stat = ier
    ELSE
        !IF (SIZE(out_str) .LT. num_tok) THEN
        !    CALL ipslerr(3,'string_split','The size of the allocated output array','is not big enough. found=',num_tok)
        !ENDIF
        stat = -2
    ENDIF

    startpos=1
    endpos=1
!    WRITE(*,*) "string_split:: str=", str

    c = 0
    if(len(token) == 0) return
    do 
!      WRITE(*,*) ""
!      WRITE(*,*) "string_split:: index=", str(endpos:)
      posn = index(str(endpos:), token)
!      WRITE(*,*) "string_split:: posn=", posn
      if(posn == 0) exit
      c = c + 1

!      WRITE(*,*) "string_split:: startpos=", startpos
!      WRITE(*,*) "string_split:: endpos=", startpos+posn-1
      out_str(c) = str(startpos:startpos+posn-1)
!      WRITE(*,*) "string_split:: out_str=", TRIM(out_str(c))

      startpos = startpos + posn 
      endpos = endpos + posn + len(token) - 1
    end do

!    WRITE(*,*) "string_split:: end, startpos=", startpos
!    WRITE(*,*) "string_split:: end, endpos=", len(str)
    ! last iteration
    out_str(c+1) = str(startpos:len(str))

!    WRITE(*,*) "string_split:: out_str=", out_str(c+1)
  END SUBROUTINE string_split

END MODULE string 
