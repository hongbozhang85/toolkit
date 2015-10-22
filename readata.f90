! 091114 by hongbo
! a simple example to show how to read data

program readdata
implicit none

  type snadata
    real(kind=8) :: z,mu,sigma
  end type

  integer,parameter :: sn_num = 307
  type(snadata) :: sd(sn_num)
  character(80) :: filename = "sn_z_mu_dmu.txt"
  logical :: alive
  character(6) :: sn_name
  integer :: i

!---read sna data---!
  write(*,*) 'Reading Supernovea Data, Uion Set'

!---whether the data file exists
  inquire (file=filename, exist=alive)
    if (.not.alive) then
       write(*,*) filename,"doesn't exist"
       stop
    else 
       write(*,*) "data file exists"
    end if

  open (unit=44, file=filename)
  do i=1,sn_num,1
     read(44,*) sn_name, sd(i)
     write(*,*) sn_name, sd(i)
  end do
  close(44)

end
