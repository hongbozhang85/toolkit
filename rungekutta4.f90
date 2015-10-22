! 091109 by hongbo
! This program is runge kutta of order 4 for solving ODE with initial condition.
! NOT for system of ODEs
! This program is the simplest rk4 just for pratice

program  rk4
implicit none

  real, parameter :: h=0.01   ! h is step length, h = x_{i+1} - x_{i}
  integer,parameter :: n=100  ! n is the number of steps 
  integer :: step=1
  real :: energy(n)           ! energy(i) is y(x_{i}) 
  real :: x=0.0,y=0.0
  real :: slope1, slope2, slope3, slope4

  real, external :: slope

  open(unit=44, file='data')

  do step=1,n,1
    slope1=slope(x,y)
    slope2=slope(x+0.5*h,y+0.5*h*slope1)
    slope3=slope(x+0.5*h,y+0.5*h*slope2)
    slope4=slope(x+h,y+h*slope3)
    energy(step)=y
!    write(44, *) energy(step)
    write(44, *) energy(step),x
    y=y+(slope1+2.0*slope2+2.0*slope3+slope4)*h/6.0
    x=x+h  
  end do

end program


function slope(x,y)
implicit none
  real :: x,y
  real :: slope

  slope = 2.0*x               ! slope = y'

return
end function

