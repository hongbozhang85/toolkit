!091210 by hongbo
!a simple exercise

!cccccccccccccccccccccccccccccccccccccccccccccc!
!ccc         Metropolis algorithm           ccc!
!ccc     with symmetry proposal density     ccc!
!cccccccccccccccccccccccccccccccccccccccccccccc!


program Metropolis
implicit none

  real,external :: distribution,proposal
  integer,parameter :: chain_length=100000
  real :: chain(chain_length),y,ran
  character(80) :: chain_file = "chain.txt"
  integer :: i

 call random_seed()
!--input the initial point--!
 write(*,*) "input the initial point"
 read(*,*) chain(1)

 open(unit=44,file=chain_file)

 write(*,*) chain(1)

 write(44,*) 2 !chain(1)

 do i=2, chain_length, 1
    call random_walk(chain(i-1),y)
    if (distribution(y)>distribution(chain(i-1))) then
       chain(i)=y
    else
        call random_number(ran)
        if(ran < distribution(y)/distribution(chain(i-1))) then
           chain(i)=y
        else
           chain(i)=chain(i-1)
        end if
    end if
    write(44,*) chain(i)
 end do


 close(44)
end program

!----------------------------------

function distribution(var)
implicit none

real :: distribution, var

  distribution = exp(-var**2)

return
end function

!---------------------------------

function proposal(x,y)
implicit none

real :: proposal,x,y

  proposal = exp(-(x-y)**2)/sqrt(3.1415926)

return
end function

!--------------------------------
!-- solve equation by newton method
!-- /int_(-/inf)^(new) dx*proposal(x,origin) = t
!-- to get 'new'. where t is random number

subroutine random_walk(origin,new)
implicit none

  real :: origin, new,t
  real,external :: proposal,intproposal
  real :: x,p, tolx=1.0E-3
  integer,parameter :: time=10000
  integer :: i


  i=0
  call random_number(t)

  x=origin

  p=x-(intproposal(x,origin)-t)/proposal(x,origin)

  do while(abs(p-x)>tolx .and. i<time)
     x=p
     p=x-(intproposal(x,origin)-t)/proposal(x,origin)
     i=i+1
  end do

  if (i==time .or. i>time) then
     write(*,*) "Method Failed"
     write(*,*) x,"new=",p,"f(x)=",intproposal(x,origin)-t,"iteration time = ",i
  end if

  new=p



end subroutine

!-------------------------------
!--integrate probability---
!--/int_(-/inf)^(new) dx*proposal(x,origin)

function intproposal(x,origin)
implicit none

  real,external :: proposal
  real :: intproposal,x,origin
  real :: xmin,xmax ! range of integral
  integer,parameter :: seg=500  ! dived the range of integral into 100 pieces
  real :: datalist(seg+1)   ! 100 pieces--> 101 data
  real :: width  ! length of every piece
  real :: s
  integer :: i

  s=0
  xmin=origin-5
  xmax=x

!---generate proposal(x_i)---!
  width=(xmax-xmin)/seg
  do i=0,seg,1
     datalist(i+1)=proposal(xmin+i*width,origin)
  end do

!---trape integrate---!
  do i=0,seg,1
     s=s+datalist(i+1)
  end do

  intproposal=width*(s-datalist(1)/2.0-datalist(seg+1)/2.0)

end function
