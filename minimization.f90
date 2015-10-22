! 091127 by hongbo
! a very simple toy code which will be used in estimating the best-fit
! gradient way
! step is very important, it cannot be large nor small, or the do while cycle will not stop

  program minimization
  implicit none

  real(kind=8),external :: f
  integer,parameter :: dimen = 2
  real(kind=8) :: beginx(dimen), gradient(dimen),stepvector(dimen,dimen)
  real(kind=8),parameter :: step=0.01,zero=1.0D-4
!step = 1.0,  zero = 1.0D-4 !step=0.01,zero=1.0D-3
  integer :: i,j,countcycle=0
!-------------------------!
!---input initial value---!
!-------------------------!

  beginx(1) = 5.0
  beginx(2) = 6.0

!!!write(*,*) "input the initial point"  
!!!read(*,*) beginx(:)

!-Unit Matrix-!

  do i=1,dimen,1
     do j=1,dimen,1
        if (i==j) then
           stepvector(i,j) = step
        else
           stepvector(i,j) = 0.0
        end if
     end do
  end do

!!! stepvector = step*reshape...

!---get gradient---!
	do i=1,dimen,1
	  gradient(i) = (f(beginx+stepvector(i,:))-f(beginx))/step
	end do

!---------------------!
!---gradient method---!
!---------------------!

!--------useful code begin----------!

  do while (dot_product(gradient,gradient)>zero)
	beginx=beginx-step*gradient
!!!  beginx=beginx-step*gradient/dot_product(gradient,gradient)
	do i=1,dimen,1
	  gradient(i) = (f(beginx+stepvector(i,:))-f(beginx))/step
	end do
        countcycle = countcycle + 1
  end do

!--------useful code end----------!
  
  write(*,*) "The minimal is", f(beginx), "at (", beginx,")"
  write(*,*) "cycle is", countcycle 

  end program

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  function f(x)
  implicit none

  integer,parameter :: dimen = 2
  real(kind=8) :: f,x(dimen)

!!!  f = x(1)**2 + x(2)**2
     f = x(1)**2 + 2*x(1) + x(2)**2

  return 
  end function
