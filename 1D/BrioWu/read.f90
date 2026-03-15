program main
implicit none
integer :: i
real(8) :: x, rho, vx, vy, vz, pre, Bx, By, Bz


 open(1,file="briowu_nonregsol.dat",action="read")
 open(2,file="briowu_nonregsol_new.dat",action="write")

 do 
    read(1,*,end=100) x, rho, pre, vx, vy, vz, Bx, By, Bz
    write(2,"(9e17.8)") x-0.5d0, rho, vx, vy, vz, pre, Bx, By, Bz
 enddo
 100 print*, "end reading"
 close(1)
 close(2)

end program main
