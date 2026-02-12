module params

! time
real(8),parameter:: timemax=0.2d0 ! simulation end time

! coordinate 
integer, parameter :: nx = 64*4       ! the number of grids in the simulation box
integer, parameter :: ngh = 2            ! the number of ghost cells
integer, parameter :: nxtot = nx+2*ngh+1 ! the total number of face-centered grids including ghost cells
integer, parameter :: is = ngh+1         ! the index of the leftmost grid
integer, parameter :: ie = nx+ngh     ! the index of the rightmost grid
real(8), parameter :: x1min = -0.5d0, x1max = 0.5d0

! indices of the primitive variables
integer, parameter :: IDN = 1
integer, parameter :: IMX = 2
integer, parameter :: IPR = 3
integer, parameter :: NVAR = 3

! indices of the conservative variables
integer, parameter :: IVX = 2
integer, parameter :: IEN = 3

! thermodynamics
real(8),parameter::gam=1.4d0 !! adiabatic index

! Riemann solver
integer, parameter :: flag_flux = 2 ! 1 (Lax), 2 (HLL)

! output 
character(20),parameter::dirname="lax" ! directory name

! snapshot
integer, parameter :: unitsnap = 17

! realtime analysis 
integer, parameter :: unitevo = 11
integer, parameter :: nevo = 1    ! the number of variables derived in the realtime analysis
end module

!===============================================================================
!
!===============================================================================
program main
use params, only: nxtot, NVAR, dirname, timemax, nevo, unitevo
implicit none

! time evolution
integer :: ntime = 0    ! counter of the timestep
real(8) :: time = 0.0d0  ! time 
real(8) :: dt   = 0.0d0  ! time width

! definition of arrays 
real(8),dimension(nxtot)::xf
real(8),dimension(nxtot)::xv
real(8),dimension(NVAR,nxtot) :: U ! conservative variables
real(8),dimension(NVAR,nxtot) :: Q ! primitive variables
real(8),dimension(NVAR,nxtot) :: F ! numerical flux

real(8) ::  phys_evo(nevo)    ! variables derived in the realtime analysis

real(8), external :: TimestepControl


     ! make the directory for output
      call makedirs(trim(dirname))

      write(6,*) "setup grids and initial condition"
      call GenerateGrid(xf, xv)
      call GenerateProblem(xv, Q)
      call Prim2Consv(Q, U)
      call Output( time, .TRUE., xv, Q )

      write(6,*) "Start the simulation"
      !open file to output the result of the realtime analysis
      open(unitevo,file=trim(dirname)//'/'//'ana.dat', action="write")
! main loop
      do 
         dt = TimestepControl(xf, Q)
         if( time + dt > timemax ) dt = timemax - time
         call BoundaryCondition(Q)
         call NumericalFlux(dt, xv, Q, F)
         call UpdateConsv( dt, xf, F, U )
         call Consv2Prim( U, Q )
         time=time+dt
         ntime=ntime+1
         call Output( time, .FALSE., xv, Q )

         if( mod(ntime,10) .eq. 0 ) then
            call RealtimeAnalysis(xv,Q,phys_evo)
            write(unitevo,*) time, phys_evo(1:nevo)
         endif

         if(time >= timemax) exit 
      enddo 
      call Output( time, .TRUE., xv, Q )
end program main
!=============================================================
! GenerateGrid
! Description:
!   Generate a 1D uniform grid for a finite-volume scheme.
!   This routine fills:
!     - xf(:): face (cell-boundary) coordinates
!     - xv(:): cell-center coordinates
!   The grid uses global parameters (x1min, x1max, nx, ngh, ...).
!
! Notes:
!   Be careful about array sizes when mixing cell-centered and face-centered
!   quantities. The number of faces is (number of cells + 1).
!=============================================================
subroutine GenerateGrid(xf, xv)
use params, only : x1min, x1max, nxtot, ngh, nx
implicit none
real(8), intent(out) :: xf(nxtot), xv(nxtot)
real(8) :: dx
integer::i

    dx=(x1max-x1min)/nx
    do i=1,nxtot
         xf(i) = dx*(i-(ngh+1))+x1min
    enddo
    do i=1,nxtot-1
         xv(i) = 0.5d0*(xf(i+1)+xf(i))
    enddo

return
end subroutine GenerateGrid
!=============================================================
! GenerateProblem
! Description:
!   Set initial conditions for a 1D Riemann problem (Sod shock tube).
!   The primitive variables Q(:,i) = (rho, v, p) are assigned based on xv(i):
!     - left state (x < 0):   rho = 1.0,   v = 0.0, p = 1.0
!     - right state (x >= 0): rho = 0.125, v = 0.0, p = 0.1
!   The routine typically initializes only the active zone (i=is:ie).
!   Ghost zones are filled later by BoundaryCondition().
!=============================================================
subroutine GenerateProblem(xv, Q)
use params, only: IDN, IVX, IPR, NVAR, is, ie, nxtot
implicit none
integer::i
real(8), intent(in ) :: xv(nxtot)
real(8), intent(out) :: Q(NVAR,nxtot)

      do i=is,ie
         if( xv(i) < 0.0d0 ) then 
             Q(IDN,i) = 1.0d0
             Q(IVX,i) = 0.0d0
             Q(IPR,i) = 1.0d0
        else 
             Q(IDN,i) = 0.125d0
             Q(IVX,i) = 0.0d0
             Q(IPR,i) = 0.1d0
         endif
      enddo

return
end subroutine GenerateProblem
!=============================================================
! BoundaryCondition
! Description:
!   Apply boundary conditions by filling ghost cells ( i < is and i > ie ) 
!   of the primitive array Q.
!=============================================================
subroutine BoundaryCondition(Q)
use params, only : IDN, IVX, IPR, NVAR, is, ie, ngh, nxtot
implicit none
real(8), intent(inout) :: Q(NVAR,nxtot)
integer::i

    do i=1,ngh 
         Q(IDN,is-i)  = Q(IDN,is-1+i)
         Q(IVX,is-i)  = Q(IVX,is-1+i)
         Q(IPR,is-i)  = Q(IPR,is-1+i)
    enddo

    do i=1,ngh
         Q(IDN,ie+i) = Q(IDN,ie-i+1)
         Q(IVX,ie+i) = Q(IVX,ie-i+1)
         Q(IPR,ie+i) = Q(IPR,ie-i+1)
    enddo

return
end subroutine BoundaryCondition
!=============================================================
! Prim2Consv
! Description:
!   Convert primitive variables Q = (rho, v, p) to conservative variables
!   U = (rho, mom, E) for the 1D Euler equations with an ideal-gas EOS.
!     - mom = rho * v
!     - E   = 0.5 * rho * v^2 + p / (gam - 1)
!   Operates on the active zone (i=is:ie).
!=============================================================
subroutine Prim2Consv(Q, U)
use params, only : IDN, IMX, IEN, IVX, IPR, NVAR, is, ie, nxtot, gam
implicit none
real(8), intent(in) :: Q(NVAR,nxtot)
real(8), intent(out) :: U(NVAR,nxtot)
integer::i

    do i=is,ie
        U(IDN,i) = Q(IDN,i)
        U(IMX,i) = Q(IDN,i)*Q(IVX,i)
        U(IEN,i)  = 0.5d0*Q(IDN,i)*Q(IVX,i)**2 + Q(IPR,i)/(gam - 1.0d0)
    enddo
      
return
end subroutine Prim2Consv
!=============================================================
! Consv2Prim
! Description:
!   Convert conservative variables U = (rho, mom, E) to primitive variables
!   Q = (rho, v, p) for the 1D Euler equations with an ideal-gas EOS.
!     - v = mom / rho
!     - p = (E - 0.5 * mom^2 / rho) * (gam - 1)
!   Operates on the active zone (i=is:ie).
!=============================================================
subroutine Consv2Prim( U, Q )
use params, only : IDN, IVX, IPR, IMX, IEN, NVAR, gam, is, ie, nxtot
implicit none
real(8), intent(in) :: U(NVAR,nxtot)
real(8), intent(out) :: Q(NVAR,nxtot)
integer::i

    do i=is,ie
        Q(IDN,i) = U(IDN,i)
        Q(IVX,i) = U(IMX,i)/U(IDN,i)
        Q(IPR,i) = ( U(IEN,i) - 0.5d0*U(IMX,i)**2/U(IDN,i) )*(gam-1.0d0)
    enddo

return
end subroutine Consv2Prim
!!=============================================================
! TimestepControl
! Description:
!   Compute a stable time step dt based on a CFL condition for the 1D Euler
!   equations.
!
! Inputs:
!   xf(:)    Face coordinates (used to compute cell widths dx)
!   Q(:,:)   Primitive variables Q=(rho, v, p) in the active zone
!
! Output:
!   dt       Time step satisfying
!              dt = CFL * min_i [ dx_i / (|v_i| + c_s,i) ]
!            where c_s = sqrt(gam * p / rho).
!=============================================================
Real(8) Function TimestepControl(xf, Q) 
use params, only : IDN, IVX, IPR, NVAR, nxtot, is, ie, gam 
implicit none
real(8), intent(in) :: xf(nxtot), Q(NVAR,nxtot)
real(8)::dtlocal
real(8)::dtmin
integer::i

    dtmin=1.0d90

    do i=is,ie
        dtlocal =(xf(i+1)-xf(i))/(abs(Q(IVX,i)) + dsqrt(gam*Q(IPR,i)/Q(IDN,i)))
         if(dtlocal .lt. dtmin) dtmin = dtlocal
    enddo

    TimestepControl = 0.3d0 * dtmin

return
end function TimestepControl
!=============================================================
! NumericalFlux
! Description:
!   Compute numerical fluxes at cell faces from cell-centered primitive states.
!   Steps:
!     1) Reconstruct left/right states at each interface (here: 1st-order,
!        piecewise-constant reconstruction).
!     2) Solve an approximate Riemann problem to obtain the interface flux.
!=============================================================
subroutine NumericalFlux( dt, xv, Q, F )
use params, only : nxtot, nvar, is, ie, flag_flux
implicit none
real(8), intent(in)  :: dt
real(8), intent(in)  :: xv(nxtot)
real(8), intent(in)  :: Q(NVAR,nxtot)
real(8), intent(out) :: F(NVAR,nxtot)
integer :: i
real(8),dimension(NVAR,nxtot):: Ql,Qr
real(8),dimension(NVAR):: flx

    do i=is-1,ie+1
        Ql(:,i+1) = Q(:,i) 
        Qr(:,i  ) = Q(:,i) 
    enddo

    do i=is,ie+1
       if( flag_flux == 1 ) then 
           call Lax((xv(i) - xv(i-1))/dt,Ql(:,i),Qr(:,i),flx)
       else if (flag_flux == 2) then
           call HLL(Ql(:,i),Qr(:,i),flx)
       endif
        F(:,i)  = flx(:)
    enddo

return
end subroutine Numericalflux
!=============================================================
! Lax
! Description:
!   Compute the numerical flux at a cell interface using the
!   Lax (Lax–Friedrichs / Rusanov-type) flux formula for the 1D Euler equations.
!
! Inputs:
!   a        Numerical dissipation speed (typically dx/dt, or a maximum wave speed)
!   Ql(:)    Left primitive state  (rho, v, p)
!   Qr(:)    Right primitive state (rho, v, p)
!
! Output:
!   flx(:)   Conservative flux (mass, momentum, energy) at the interface.
!
! Method:
!   1) Convert Ql and Qr to conservative variables Ul, Ur.
!   2) Compute physical fluxes Fl(Ul), Fr(Ur).
!   3) Return
!        flx = 0.5*(Fl + Fr) - 0.5*a*(Ur - Ul)
!
! Notes:
!   - a controls numerical viscosity. Using a = dx/dt corresponds to the classic
!     Lax–Friedrichs scheme; using a = max(|v|+cs) gives the local Lax–Friedrichs
!     (Rusanov) flux.
!=============================================================
subroutine Lax(dxdt,Ql,Qr,flx)
use params, only : IDN, IVX, IPR, IMX, IEN, NVAR, gam
implicit none
real(8),intent(in) :: dxdt
real(8),intent(in) :: Ql(NVAR), Qr(NVAR)
real(8),intent(out) :: flx(NVAR)
integer :: n
real(8):: Ul(NVAR), Ur(NVAR)
real(8):: Fl(NVAR), Fr(NVAR)

    ! conserved variables in the left and right states
    Ul(IDN) = Ql(IDN)
    Ur(IDN) = Qr(IDN)

    Ul(IMX) = Ql(IDN)*Ql(IVX)
    Ur(IMX) = Qr(IDN)*Qr(IVX)

    Ul(IEN) = 0.5d0*Ql(IDN)*Ql(IVX)**2 + Ql(IPR)/(gam - 1.0d0)
    Ur(IEN) = 0.5d0*Qr(IDN)*Qr(IVX)**2 + Qr(IPR)/(gam - 1.0d0)

    ! flux in the left and right states
    Fl(IDN) = Ul(IMX)
    Fr(IDN) = Ur(IMX)

    Fl(IMX) = Ql(IPR) + Ql(IDN)*Ql(IVX)**2 
    Fr(IMX) = Qr(IPR) + Qr(IDN)*Qr(IVX)**2 

    Fl(IEN) = ( gam*Ql(IPR)/(gam - 1.0d0) + 0.5d0*Ql(IDN)*Ql(IVX)**2)*Ql(IVX)
    Fr(IEN) = ( gam*Qr(IPR)/(gam - 1.0d0) + 0.5d0*Qr(IDN)*Qr(IVX)**2)*Qr(IVX)

    do n=1,NVAR 
        flx(n)  = 0.5d0*(Fl(n) + Fr(n)) - 0.5d0*dxdt*(Ur(n) - Ul(n))
    enddo


return
end subroutine Lax
!=============================================================
! HLL
! Description:
!   Compute the interface flux using the HLL approximate Riemann solver for
!   the 1D Euler equations.
!
! Inputs:
!   Ql(:)  Left primitive state  (rho, v, p)
!   Qr(:)  Right primitive state (rho, v, p)
!
! Output:
!   flx(:) Conservative flux (mass, momentum, energy)
!=============================================================
subroutine HLL(Ql,Qr,flx)
use params, only : IDN, IVX, IPR, IMX, IEN, NVAR, gam
implicit none
real(8),intent(in)  :: Ql(NVAR), Qr(NVAR)
real(8),intent(out) :: flx(NVAR)
integer :: n
real(8):: Ul(NVAR), Ur(NVAR)
real(8):: Fl(NVAR), Fr(NVAR)
real(8):: csl,csr
real(8):: sl, sr

    ! conserved variables in the left and right states
    Ul(IDN) = Ql(IDN)
    Ur(IDN) = Qr(IDN)

    Ul(IMX) = Ql(IDN)*Ql(IVX)
    Ur(IMX) = Qr(IDN)*Qr(IVX)

    Ul(IEN) = 0.5d0*Ql(IDN)*Ql(IVX)**2 + Ql(IPR)/(gam - 1.0d0)
    Ur(IEN) = 0.5d0*Qr(IDN)*Qr(IVX)**2 + Qr(IPR)/(gam - 1.0d0)

    ! flux in the left and right states
    Fl(IDN) = Ul(IMX)
    Fr(IDN) = Ur(IMX)

    Fl(IMX) = Ql(IPR) + Ql(IDN)*Ql(IVX)**2 
    Fr(IMX) = Qr(IPR) + Qr(IDN)*Qr(IVX)**2 


    Fl(IEN) = ( gam*Ql(IPR)/(gam - 1.0d0) + 0.5d0*Ql(IDN)*Ql(IVX)**2)*Ql(IVX)
    Fr(IEN) = ( gam*Qr(IPR)/(gam - 1.0d0) + 0.5d0*Qr(IDN)*Qr(IVX)**2)*Qr(IVX)

    csl = dsqrt(gam*Ql(IPR)/Ql(IDN))
    csr = dsqrt(gam*Qr(IPR)/Qr(IDN))

    sl = min(Ql(IVX),Qr(IVX)) - max(csl,csr)
    sr = max(Ql(IVX),Qr(IVX)) + max(csl,csr)

    if( sl > 0.0d0 ) then 
        do n=1,NVAR
            flx(n) = Fl(n)
        enddo
    else if (sr < 0.0d0 ) then
        do n=1,NVAR
             flx(n) = Fr(n)
        enddo
    else 
        do n=1,NVAR
            flx(n)  = (sr*Fl(n) - sl*Fr(n) + sl*sr*( Ur(n) - Ul(n) ))/(sr - sl)
        enddo
    endif

return
end subroutine HLL
!=============================================================
! UpdateConsv
! Description:
!   Update conservative variables U by one finite-volume time step:
!     U_i^{n+1} = U_i^{n} - dt * (F_{i+1/2} - F_{i-1/2}) / dx_i
!   where dx_i is computed from face coordinates xf.
!
! Inputs:
!   dt     Time step
!   xf(:)  Face coordinates
!   F(:,:) Interface fluxes (at faces)
!
! In/Out:
!   U(:,:) Conservative variables updated in-place (active zone i=is:ie).
!=============================================================
subroutine UpdateConsv( dt, xf, F, U )
use params, only : NVAR, is, ie, nxtot
implicit none
real(8), intent(in)  :: F(NVAR,nxtot), dt, xf(nxtot)
real(8), intent(inout) :: U(NVAR,nxtot)
integer::i,n

    do i=is,ie
        do n=1,NVAR
            U(n,i) = U(n,i) + dt*(- F(n,i+1) + F(n,i))/(xf(i+1)-xf(i)) 
        enddo
    enddo

return
end subroutine UpdateConsv
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
!   Q(:,:)   Primitive variables to be written
!
! Notes:
!   With variable dt, consider using a "while (time >= next_output_time)" style
!   to avoid missing outputs when the simulation time jumps over an output time.
!=============================================================
subroutine Output( time, flag, xv, Q )
use params, only: IDN, IVX, IPR, NVAR, nxtot, dirname, is, ie, unitsnap
implicit none
real(8), intent(in) :: time
logical, intent(in) :: flag
real(8), intent(in) :: xv(nxtot), Q(NVAR,nxtot)
real(8), parameter:: dtsnap=5.0d-3
integer::i
character(40) :: filename
real(8), save :: tsnap = - dtsnap
integer :: nsnap = 0

    if( .not.flag) then
          if( time + 1.0d-14.lt. tsnap+dtsnap) return
    endif

    write(filename,'(i5.5)') nsnap
    filename = trim(dirname)//"/snap"//trim(filename)//".dat"
    open(unitsnap,file=filename,form='formatted',action="write")
    write(unitsnap,"(a2,f6.4)") "# ",time
    do i=is,ie
         write(unitsnap,*) xv(i), Q(IDN,i), Q(IVX,i), Q(IPR,i)
    enddo
    close(unitsnap)

    write(6,*) "output:",nsnap,time

    nsnap=nsnap+1
    tsnap=tsnap + dtsnap

return
end subroutine Output
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
!   Q(:,:)   Primitive variables Q=(rho, v, p)
!   U(:,:)   Conservative variables U=(rho, mom, E) (optional but useful)
!=============================================================
subroutine RealtimeAnalysis(xv,Q,phys_evo)
use params, only : IDN, IVX, IPR, NVAR, is, ie, gam, nevo, nx, nxtot
implicit none
real(8), intent(in)  :: xv(nxtot), Q(NVAR,nxtot)
real(8), intent(out) :: phys_evo(nevo)
integer :: i
real(8) :: tmp

      tmp = 0.0d0
      do i=is,ie
           tmp = tmp + Q(IDN,i)*Q(IPR,i)*xv(i)
      enddo
      phys_evo(1:nevo) = tmp/dble(nx)
      
return
end subroutine
!!=============================================================
!! AnalysisAfterSimu
!! Description:
!!   Perform post-processing after the time integration has finished.
!!   This routine is intended to be called once at the end of the run to produce
!!   summary diagnostics and/or derived outputs.
!!
!! Inputs:
!!   dirname  Output directory name where snapshots/logs are stored
!!   xv(:)    Cell-center coordinates
!!   Q(:,:)   Final primitive variables Q=(rho, v, p)
!!   U(:,:)   Final conservative variables U=(rho, mom, E)
!!=============================================================
!!subroutine AnalysisAfterSimu(time,xf,xv,Q)
!!real(8), intent(in)  :: xf(:), xv(:), Q(:,:)
!!real(8), intent(in)  :: time
!!integer :: i
!!real(8) :: tmp
!!
!!      tmp = 0.0d0
!!      do i=is,ie
!!           tmp = tmp + 1.0d0
!!      enddo
!!      
!!return
!!end subroutine
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
    istat = system(trim(cmd))
    if( istat .ne. 0 ) then
        print*, "makedirs: command failed, status=", istat
        print*, "cmd: ", trim(cmd)
    endif

end subroutine makedirs

