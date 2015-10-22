! 100518 by hongbo Hongbo Zhang 19/05/2010 04:00:08 
! an example program to do cubic spline interpolation
! tri-diagonal method [google: cubic spline/tri-diagonal]
! one must know the 2nd derivative at the end points to do cubic spline in this program.
! C_S algorithm will also work if 1st derivative at end points are known, but not for this program.
! This program is only for natural boundary condition, or 2nd derivative are known.
! referrence http://www.vbforums.com/showthread.php?t=480806

  program main
  implicit none

!  character(80) :: datafile="spline_data.txt"
  interface
    subroutine cubic_spline(n,datalist,result,startpoint,endpoint)
      implicit none
      integer :: n
      real(kind=8) :: datalist(2,n)
      real(kind=8) :: result(4,n-1)
      real(kind=8), optional :: startpoint, endpoint
    end subroutine cubic_spline
  end interface

  integer,parameter :: n = 3
  real(kind=8) :: datalist(2,n)
  real(kind=8) :: result(4,n-1)
  data datalist /0.0,0.0,0.5,0.125,1.0,1.0/
  real(kind=8) :: startpoint=0.0, endpoint=6.0

!!!  write(*,*) datalist(1,2)

!  call cubic_spline(n,datalist,result)
   call cubic_spline(n,datalist,result,startpoint,endpoint)

  end program main

!============================================================!
!                  cubic spline subroutine                   !
!============================================================!

  subroutine cubic_spline(n,datalist,result,startpoint,endpoint)
  implicit none

! n        : the number of pairs in datalist
! datalist : a list of data (x_i,y_i) to be interpolated
!            datalist(:,i) = (x_i,y_i)
! result   : the result returned. the coefficients of n-1 cubic polynomials 
!            y_i= a_i x^3 + b_i x^2 + c_i x + d_i 
!            result(:,i)= (a_i,b_i,c_i,d_i)
! startpoint, endpoint (optional): 2nd derivative at x_0 and x_n. if not present, natrual boundary condition will be used.

  integer :: n
  real(kind=8) :: datalist(2,n)
  real(kind=8) :: result(4,n-1)
  real(kind=8), optional :: startpoint, endpoint
  real(kind=8) :: start2de,end2de
  real(kind=8) :: A(n,n)
  real(kind=8) :: B(n)
  real(kind=8) :: h(n-1)
  integer :: i
  real(kind=8) :: p(n), q(n) , y(n), y2(n)
! in Crout decomposition of tri-diagonal matrix, A=LU, p(n) is diagonal of L, and q(n) is upper diagonal of U
! y2(i) is the 2nd derivative at the ith point, y2(i)=f''(x_i).
! y(i) is an auxiliary varible
  
! Ax=6B where x is array of y''_i; A is tri-diagonal matrix
! A(i,j) means ith line (HANG), jth row (LIE)

  if (present(startpoint)) then
     start2de = startpoint
  else
     start2de = 0.0
  end if

  if (present(endpoint)) then
     end2de = endpoint
  else
     end2de = 0.0
  end if

  A = 0.0

  do i=1, n-1, 1
     h(i) = datalist(1,i+1)-datalist(1,i)
  end do
  
  B(1) = start2de/6.0
  B(n) = end2de/6.0

  do i=2, n-1, 1
     B(i) = (datalist(2,i+1)-datalist(2,i))/h(i)-(datalist(2,i)-datalist(2,i-1))/h(i-1)
  end do

  A(1,1)=1.0
  A(n,n)=1.0

  do i=2, n-1, 1
     A(i,i-1) = h(i-1)
     A(i,i) = 2*(h(i-1)+h(i))
     A(i,i+1) = h(i)
  end do

! Ax=6B => x=6*A^(-1)*B inverse matrix will cause instability or cost amount of memory
! because A is tri-diagonal, we can use some special method to solve this tri-diagonal linear system
! referrence 林成森1998 P79
!=================================================!
!              tri-diagonal  system               !
!=================================================!
  
  p(1)=A(1,1)        !1.0
  q(1)=A(1,2)/A(1,1) !0.0
  
  do i=2,n-1,1
     p(i) = A(i,i)-A(i,i-1)*q(i-1)
     q(i) = A(i,i+1)/p(i)
  end do

  p(n)=A(n,n)-A(n,n-1)*q(n-1)

  y(1)=6.0*B(1)/A(1,1)
  
  do i=2,n,1
     y(i) = (6.0*B(i)-A(i,i-1)*y(i-1))/p(i)
  end do

  y2(n)=y(n)

  do i=n-1,1,-1
     y2(i) = y(i)-q(i)*y2(i+1)
  end do

!!!  write(*,*) y2

!=================================================!
!               tri-diagonal   end                !
!=================================================!
  
  do i=1,n-1,1
     result(1,i) = (y2(i+1)-y2(i))/6.0/h(i)
     result(2,i) = y2(i)/2.0
     result(3,i) = (datalist(2,i+1)-datalist(2,i))/h(i)-y2(i+1)*h(i)/6.0-y2(i)*h(i)/3.0
     result(4,i) = datalist(2,i)
  end do
  
  write(*,*) result

!================================================!
!                   Draw Figure                  !
!================================================!
!  call system()

  end subroutine cubic_spline




  
