! 091110 by hongbo
! This is a Monte Carlo approach to do integration.
! This program is the simplest and most naive just for pratice
!-----------------------------------------------------------------
!module intergral_range
!  implicit none
!  real :: x,y,integral
!  real,parameter :: xi=0.0,xf=1.0,ymin=-1,ymax=2
! range for scattering points. In this problem, f(x)=2x, x in (0,1)
! we can safely choose ymin=-1 and ymax=2
!  real :: area
!  area=(xf-xi)*(ymax-ymin)
!
!end module
! HOW TO USE MODULE IN LINUX ?
!-------------------------------------------------------------------

program montecarlo
! use intergral_range
  implicit none

  real,external :: fx
  real,external :: random
  integer, parameter :: n=100000000
  integer :: s=0,i
  real :: x,y,integral
  real,parameter :: xi=0.0,xf=1.0,ymin=-1,ymax=2
! range for scattering points. In this problem, f(x)=2x, x in (0,1)
! we can safely choose ymin=-1 and ymax=2
  real :: area
  area=(xf-xi)*(ymax-ymin)

  
  call random_seed()
  do i=1,n,1
     x=random(xi,xf)
     y=random(ymin,ymax)
     if ( y <= fx(x) .and. y>=0) then 
! in this problem, f(x)>0 for all x. this maybe not true in another problem, 
        s=s+1
     end if
  end do
  
  integral=area*real(s)/real(n)
  write(*,*) integral

end program

!-------------------------------------------------------------------

function fx(x)
  implicit none

  real :: x
  real :: fx

!  fx=2.0*x
   fx=exp(x)

  return
end function

!------------------------------------------------------------------

function random(low,up)
  implicit none
  
  real :: random,low,up,lens,t
  
  lens=up-low
  call random_number(t)
  random=low+lens*t
  return

end function
