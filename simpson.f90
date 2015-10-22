  program simpson
  implicit none

  real,external :: fx
  real :: xmin,xmax,s,h
  integer :: i, seg=10


  xmin = 0.0
  xmax = 10.0
  s = 0.0
  h = (xmax-xmin)/seg

  do i=0,seg-1,1
     s=s+h*(fx(xmin+i*h) + 4*fx(xmin+i*h+h/2.0) + fx(xmin+i*h+h))/6.0
  end do

  write(*,*) s,log(11.0)

  end program

!--------------------------------

  function fx(x)
  implicit none

  real :: fx,x
 
!  fx = x**3
   fx = 1/(1+x)

  return
  end function
