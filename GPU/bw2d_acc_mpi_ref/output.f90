
!=============================================================
! Output
! Description:
!   Write snapshot output of the solution (cell centers and primitive variables).
!   The routine decides whether to output based on an internal snapshot clock
!   (e.g., tsnap and dtsnap) and a logical flag argument.
!
! Inputs:
!   flag     If true, check output condition and write if needed
!   xv(:)    Cell-center coordinates
!   yv(:)    Cell-center coordinates
!   Q(:,:,:)   Primitive variables to be written
!
! Notes:
!   With variable dt, consider using a "while (time >= next_output_time)" style
!   to avoid missing outputs when the simulation time jumps over an output time.
!=============================================================
subroutine Output( time, flag, xv, yv, Q )
use params, only : IDN, IVX, IVY, IVZ, IPR, IBX, IBY, IBZ, IPS, NVAR, & 
                   nxtot, nytot, nx, ny, is, ie, js, je, gam, &
                   flag_binary, dirname, dtsnap, unitsnap
use mpimod
implicit none
real(8), intent(in) :: time! false --> output per dtsnap, true --> force to output
logical, intent(in) :: flag
real(8), intent(in) :: xv(nxtot), yv(nytot), Q(NVAR,nxtot,nytot)

integer::i,j
character(100)::filename
real(8), save :: tsnap = - dtsnap
integer, save :: nsnap = 0

real(8),dimension(is:ie) :: xout
real(8),dimension(js:je) :: yout
real(4),dimension(1:5,is:ie,js:je) :: Qouthyd
real(4),dimension(6:Nvar,is:ie,js:je) :: Qoutmhd


    if( .not.flag) then
        if( time + 1.0d-14.lt. tsnap+dtsnap) return
    endif
!$acc update host(Q)
    xout(is:ie) = xv(is:ie)
    yout(is:ie) = yv(is:ie)
    Qouthyd(1:5,is:ie,js:je) = real(Q(1:5,is:ie,js:je)) ! single precision
    Qoutmhd(6:Nvar,is:ie,js:je) = real(Q(6:Nvar,is:ie,js:je))! single precision
    write(filename,'((i2.2)(a1)(i5.5))') myid,"-",nsnap
    if( flag_binary ) then
        filename = trim(dirname)//"/snap"//trim(filename)//".bin"
        open(unitsnap,file=filename,form='unformatted',access="stream",action="write")
        write(unitsnap) time
        write(unitsnap) nx
        write(unitsnap) ny
        write(unitsnap) 5
        write(unitsnap) NVAR - 5
        write(unitsnap) xout
        write(unitsnap) yout
        write(unitsnap) Qouthyd
        write(unitsnap) Qoutmhd ! single precision
        close(unitsnap)
        print *, "output binary file:  ",filename,time
    else 
        filename = trim(dirname)//"/snap"//trim(filename)//".dat"
        open(unitsnap,file=filename,form='formatted',action="write")
        write(unitsnap,*) "# time = ",time
        write(unitsnap,*) "#nx, ny = ", nx, ny
          do j=js,je
          do i=is,ie
              write(unitsnap,*) xv(i), yv(j), Q(IDN,i,j), Q(IVX,i,j), Q(IVY,i,j), Q(IVZ,i,j), &
                                Q(IPR,i,j), Q(IBX,i,j), Q(IBY,i,j), Q(IBZ,i,j) , Q(IPS,i,j)

          enddo
          enddo
          close(unitsnap)
          print *, "output ascii file:  ",filename,time
       endif


    nsnap=nsnap+1
    tsnap=tsnap + dtsnap

return
end subroutine Output
!=============================================================
! makedirs
! Description:
!   Create an output directory if it does not exist.
!   This routine runs the OS command:
!     mkdir -p 'outdir'
!   so that:
!     - If the directory already exists: it succeeds and does nothing.
!     - If parent directories are missing: they are created.
!=============================================================
subroutine makedirs(outdir)
implicit none
integer :: istat
character(len=*), intent(in) :: outdir
character(len=1024) :: cmd = ""
    if (len_trim(outdir) == 0) then 
       write(*,*) "makedirs: outdir is empty" 
       stop 1 
    end if

    cmd = "mkdir -p '" // trim(outdir) // "'"
    call execute_command_line(trim(cmd),exitstat=istat)
    if( istat .ne. 0 ) then
        print*, "makedirs: command failed, status=", istat
        print*, "cmd: ", trim(cmd)
    endif

end subroutine makedirs
!=============================================================
! RealtimeAnalysis
! Description:
!   Perform on-the-fly diagnostics during the simulation loop.
!   This routine is intended to be called every chosen interval to 
!   monitor the run without generating heavy I/O.
!
! Inputs:
!   time     Current simulation time
!   step     Current step index
!   xv(:)    Cell-center coordinates
!   yv(:)    Cell-center coordinates
!   Q(:,:,:) Primitive variables Q=(rho, v, p)
!   U(:,:,:) Conservative variables U=(rho, mom, E) (optional but useful)
!=============================================================
subroutine RealtimeAnalysis(xv,yv,Q,phys_evo)
use params, only : IDN, IVX, IVY, IVZ, IPR, IBX, IBY, IBZ, &
                   IMX, IMY, IMZ, IEN, IPS, NVAR, nxtot, nytot, &
                   is, ie, js, je, gam, nevo
implicit none
real(8), intent(in)  :: xv(nxtot), yv(nytot), Q(NVAR,nxtot,nytot)
real(8), intent(out) :: phys_evo(nevo)
integer::i,j
real(8) :: tmp

!$acc update host(Q)      
      tmp = 0.0d0
      do j=js,je
      do i=is,ie
          tmp = tmp + Q(IDN,i,j)*Q(IPR,i,j)/(xv(i)+yv(j))
      enddo
      enddo
      phys_evo(1:nevo) = 0.0d0
      
return
end subroutine
