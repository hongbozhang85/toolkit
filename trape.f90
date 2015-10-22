! 091111 by hongbo
! trape method to do integration

program trape
implicit none

  real,external :: fx
  real,parameter :: pi=3.1415,xmin=0,xmax=1 ! range of integral
  integer,parameter :: seg=100000  ! dived the range of integral into 100 pieces
  real :: datalist(seg+1)   ! 100 pieces--> 101 data
  real :: width  ! length of every piece
  real :: s=0,ans
  integer :: i

!---generate fx(x_i)---!
  width=(xmax-xmin)/seg
  do i=0,seg,1
     datalist(i+1)=fx(xmin+i*width)
  end do

!---trape integrate---!
  do i=0,seg,1
     s=s+datalist(i+1)
  end do

  ans=width*(s-datalist(1)/2.0-datalist(seg+1)/2.0)
  
  write(*,*) ans

end program

!----------------------------------------

function fx(x)
implicit none
  real :: fx,x
  
!  fx=sin(x)
   fx=exp(x)
!  fx=2*x

return
end function
