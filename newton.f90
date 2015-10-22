!091210 by hongbo
!a simple exercise of Newton-Raphson Method. Or Tangent Method to find root
!Picard Algorithm with g(x)=x-f(x)/f'(x)


program Newton
implicit none

real, external :: func,derifunc
real :: x,p, tolx=1.0E-4
integer,parameter :: time=1000
integer :: i=0

  write(*,*) "input initial data"
  read(*,*) x

  p=x-func(x)/derifunc(x)

  do while(abs(p-x)>tolx .and. i<time)
     x=p
     p=x-func(x)/derifunc(x)
     i=i+1
  end do

  if (i==time .or. i>time) then
     write(*,*) "Method Failed"
  end if

  write(*,*) "x=",p,"f(x)=",func(p),"iteration time = ",i


end program

!------------------------

function func(x)
implicit none

real :: func,x

  func = exp(x)+x

return
end function

!------------------------

function derifunc(x)
implicit none

real :: derifunc,x,eps=1.0E-3
real,external :: func
 
 derifunc = exp(x)+1

! derifunc = (func(x+eps)-func(x))/eps

return
end function
