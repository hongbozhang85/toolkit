!091210 by hongbo
!to find root by half cut method 区间半分法

program half_cut
implicit none

real,external :: func
real :: a=0,b=0,p, tolf=1.0E-4, tolx=1.0E-4 !tolerant
integer :: time=1000, i=0

!--read initial data---!

  do while (func(a)*func(b)>0)
    write(*,*) "f(a)=",func(a),"f(b)=",func(b),"input initial a & b, b>a"
    read(*,*) a,b
  end do

  p=(a+b)/2 

!---------------!
!---main code---!
!---------------!

  do while (abs(func(p))>tolf .and. b-a>tolx .and. i<time)
     if (func(p)>0) then
        b=p
     else
        a=p
     end if
     p=(a+b)/2
     i=i+1
!!!  write(*,*) "x=",p,"f(x)=",func(p),"iteration time = ",i
  end do

!---main code ends---!
  
  if (i==time) then
     write(*,*) "Method Failed"
  end if
  
  write(*,*) "x=",p,"f(x)=",func(p),"iteration time = ",i

end program

!-----------------

function func(x)
implicit none

real :: func,x

!  func=x-5
   func=exp(x)+x

return
end function
