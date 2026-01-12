!************************************************
! original code was written by W. Ju on July 2004
! in C program.
! This Fortran version is updated by Jun Wang on
!  5/12/2016
!************************************************

subroutine s_coszs(lat,lon,CosZs,hr,Hsolar1)
   use shr_kind_mod,only: r8=>shr_kind_r8
   use beps_con,only: PI
   use beps_time_manager
   implicit none

!!integer,intent(in)  :: day,hour    ! the hour in day and the day in year (hour begins from 0, ends at 23)
   real(r8),intent(in) :: lat,lon
   real(r8),intent(out):: CosZs       ! Solar Zenith
   real(r8),intent(out):: hr
   real(r8),intent(out):: Hsolar1

   real(r8)            :: lon1        ! change 0.5~359.5 into -179.5~179.5
   integer             :: day,hour
   integer             :: yr,mn,dd,tod
   character(len = 30)  :: calendar
   real(r8)  :: Delta,Lat_arc         !  Hsolar1
   real(r8)  :: doy      !!days of the year

   call get_curr_date(yr,mn,dd,tod)
   hour   = tod/3600

   day       = get_curr_calday()
   calendar  = get_calendar()

!! determine whether it is leap year
   if(trim(calendar) == trim(NO_LEAP_C)) then
      doy   = 365.
   else if(trim(calendar) == trim(GREGORIAN_C)) then
      if((mod(yr,4) ==0 .and. mod(yr,100) /= 0 ).or. mod(yr,400) ==0) then
         doy  = 366.
      else
         doy  = 365.
      end if
   end if


   Delta = 0.006918-0.399912*cos(day*2.0*PI/doy)+0.070257*sin(day*2.0*PI/doy) &
      -0.006758*cos(day*4.0*PI/doy)+0.000907*sin(day*4.0*PI/doy)
!delta is the declination angle of sun.

!! longitude 0.5~359.5  => -179.5~179.5

   if(lon > 180 .and. lon < 360) then
      lon1 = lon - 360.
   else
      lon1 = lon
   end if

   hr    = hour + lon1/15.0
   if(hr > 24) hr = 24 - hr
   if(hr < 0 ) hr = 24 + hr

   Lat_arc   = PI*lat/180.0
   Hsolar1   = (hr - 12.0)*2.0*PI/24.0  !local hour angle in arc.
   CosZs     = cos(Delta)*cos(Lat_arc)*cos(Hsolar1)+sin(Delta)*sin(Lat_arc)

end subroutine

