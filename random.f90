! 091110 by hongbo
! random(a,b) generates random number between a and b
! PAY ATTENTION: make sure a<b or you will make mistake. 
! Of course, you can avoid mistake by slightly modifying this code
! 
! PAY ATTENTION need "call random_seed()" between use this function


function random(low,up)
  implicit none
  
  real :: random,low,up,lens,t
  
  lens=up-low
  call random_number(t)
  random=low+lens*t
  return

end function
