module params
real(8),parameter:: timemax=0.1d0 ! simulation end time

integer,parameter::nx=256*1       ! the number of grids in the simulation box
integer,parameter::ngh=2            ! the number of ghost cells
integer,parameter::nxtot=nx+2*ngh+1 ! the total number of grids including ghost cells
integer,parameter::is=ngh+1         ! the index of the leftmost grid
integer,parameter::ie=nx+ngh     ! the index of the rightmost grid
real(8),parameter:: xmin=-0.5d0,xmax=0.5d0

real(8),parameter:: Bx=0.75

! indices of the primitive variables
integer, parameter :: IDN = 1
integer, parameter :: IMX = 2
integer, parameter :: IMY = 3
integer, parameter :: IMZ = 4
integer, parameter :: IPR = 5
integer, parameter :: IBY = 6
integer, parameter :: IBZ = 7
integer, parameter :: NVAR = 7

! indices of the conservative variables
integer, parameter :: IVX = 2
integer, parameter :: IVY = 3
integer, parameter :: IVZ = 4
integer, parameter :: IEN = 5

! adiabatic index
real(8),parameter::gam = 5.0d0/3.0d0 

! CFL number
real(8),parameter :: cfl_number = 0.1d0

! output 
character(20),parameter::dirname="lax" ! directory name

! snapshot
integer, parameter :: unitsnap = 17
real(8), parameter :: dtsnap=5.0d-3

! realtime analysis 
integer, parameter :: unitevo = 11
integer, parameter :: nevo = 1    ! the number of variables derived in the realtime analysis

end module
!===============================================================================
!
!===============================================================================
program main
use params, only: nxtot, NVAR, dirname, timemax, unitevo, nevo
implicit none

! time evolution
integer :: ntime = 0    ! counter of the timestep
real(8) :: time = 0.0d0  ! time 
real(8) :: dt = 0.0d0  ! time width


! definition of arrays 
real(8),dimension(nxtot)::xf,xv
real(8),dimension(NVAR,nxtot) :: U
real(8),dimension(NVAR,nxtot) :: Q
real(8),dimension(NVAR,nxtot) :: F

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
         call BoundaryCondition( Q)
         call NumericalFlux( dt, xv, Q, F )
         call UpdateConsv( dt, xf, F, U)
         call Consv2Prim( U, Q )
         time=time+dt
         print*,"time = ",time, "dt = ",dt
         ntime = ntime + 1
         call Output( time, .FALSE., xv, Q )

         if( mod(ntime,10) .eq. 0 ) then
            call RealtimeAnalysis(xv,Q,phys_evo)
            write(unitevo,*) time, phys_evo(1:nevo)
         endif

         if(time >= timemax) exit 
      enddo 
      call Output( time, .TRUE., xv, Q )
!      call AnalysisAfterSimu(time,xv,Q)

      write(6,*) "program has been finished"

end program main
!=============================================================
! GenerateGrid
! Description:
!   Generate a 1D uniform grid for a finite-volume scheme.
!   This routine fills:
!     - xf(:): face (cell-boundary) coordinates
!     - xv(:): cell-center coordinates
!   The grid uses global parameters (xmin, xmax, nx, ngh, ...).
!
! Notes:
!   Be careful about array sizes when mixing cell-centered and face-centered
!   quantities. The number of faces is (number of cells + 1).
!=============================================================
subroutine GenerateGrid(xf, xv)
use params, only : xmin, xmax, nxtot, ngh, nx
implicit none
real(8), intent(out) :: xf(nxtot), xv(nxtot)
real(8) :: dx,dy
integer::i

    dx=(xmax-xmin)/nx
    do i=1,nxtot 
        xf(i) = dx*(i-(ngh+1))+xmin
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
!   
!   The routine typically initializes only the active zone (i=is:ie).
!   Ghost zones are filled later by BoundaryCondition().
!=============================================================
subroutine GenerateProblem(xv, Q)
use params, only: IDN, IVX, IVY, IVZ, IPR, IBY, IBZ, Bx, NVAR, is, ie, nxtot
implicit none
integer::i
real(8), intent(in ) :: xv(nxtot)
real(8), intent(out) :: Q(NVAR,nxtot)

    do i=is,ie
        if( xv(i) < 0.0d0 ) then 
             Q(IDN,i) = 1.0d0
             Q(IVX,i) = 0.0d0
             Q(IVY,i) = 0.0d0
             Q(IVZ,i) = 0.0d0
             Q(IPR,i) = 1.0d0
             Q(IBY,i) = 1.0d0
             Q(IBZ,i) = 0.0d0
        else 
             Q(IDN,i) = 0.125d0
             Q(IVX,i) = 0.0d0
             Q(IVY,i) = 0.0d0
             Q(IVZ,i) = 0.0d0
             Q(IPR,i) = 0.1d0
             Q(IBY,i) = -1.0d0
             Q(IBZ,i) = 0.0d0
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
use params, only : IDN, IVX, IVY, IVZ, IPR, IBY, IBZ, NVAR, is, ie, ngh, nxtot
implicit none
real(8), intent(inout) :: Q(NVAR,nxtot)
integer::i,ihy

      do i=1,ngh
      do ihy=1,NVAR
          Q(ihy,is-i)  = Q(ihy,is-1+i)
      enddo
      enddo

      do i=1,ngh
      do ihy=1,NVAR
          Q(ihy,ie+i) = Q(ihy,ie-i+1)
      enddo
      enddo

return
end subroutine BoundaryCondition
!-------------------------------------------------------------------
!       Primitive variables ===> Conservative variables
!       Input  : Q
!       Output : U
!-------------------------------------------------------------------
subroutine Prim2Consv(Q, U)
use params, only : IDN, IVX, IVY, IVZ, IPR, IBY, IBZ, &
                   IMX, IMY, IMZ, IEN, NVAR, Bx, is, ie, nxtot, gam
implicit none
real(8), intent(in)  :: Q(NVAR,nxtot)
real(8), intent(out) :: U(NVAR,nxtot)
integer::i

      do i=is,ie
          U(IDN,i) = Q(IDN,i)
          U(IMX,i) = Q(IDN,i)*Q(IVX,i)
          U(IMY,i) = Q(IDN,i)*Q(IVY,i)
          U(IMZ,i) = Q(IDN,i)*Q(IVZ,i)
          U(IEN,i) = 0.5d0*Q(IDN,i)*( Q(IVX,i)**2 + Q(IVY,i)**2 + Q(IVZ,i)**2 ) &
                   + 0.5d0*( Bx**2 + Q(IBY,i)**2 + Q(IBZ,i)**2 ) &
                   + Q(IPR,i)/(gam - 1.0d0)
          U(IBY,i) = Q(IBY,i)
          U(IBZ,i) = Q(IBZ,i)
      enddo
      
return
end subroutine Prim2Consv
!-------------------------------------------------------------------
!       Conservative variables ===> Primitive variables
!       Input  : U
!       Output : Q
!-------------------------------------------------------------------
subroutine Consv2Prim( U, Q )
use params, only : IDN, IVX, IVY, IVZ, IPR, IBY, IBZ, &
                   IMX, IMY, IMZ, IEN, NVAR, Bx, is, ie, nxtot, gam
implicit none
real(8), intent(in) :: U(NVAR,nxtot)
real(8), intent(out) :: Q(NVAR,nxtot)
integer::i
real(8) :: inv_d;

    do i=is,ie
        Q(IDN,i) = U(IDN,i)
        inv_d = 1.0d0/U(IDN,i)
        Q(IVX,i) = U(IMX,i)*inv_d
        Q(IVY,i) = U(IMY,i)*inv_d
        Q(IVZ,i) = U(IMZ,i)*inv_d
        Q(IPR,i) = ( U(IEN,i) &
                    - 0.5d0*(U(IMX,i)**2 + U(IMY,i)**2 + U(IMZ,i)**2)*inv_d  &
                    - 0.5d0*(Bx**2 + U(IBY,i)**2 + U(IBZ,i)**2) )*(gam-1.0d0)
        Q(IBY,i) = U(IBY,i)
        Q(IBZ,i) = U(IBZ,i)
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
use params, only : IDN, IVX, IPR, IBY, IBZ, Bx, NVAR, nxtot, is, ie, gam, &
                   cfl_number
implicit none
real(8), intent(in) :: xf(nxtot), Q(NVAR,nxtot)
real(8)::dtl1
real(8)::dtl2
real(8)::dtl3
real(8)::dtlocal
real(8)::dtmin,cf
integer::i

    dtmin=1.0d90

    do i=is,ie
         cf = dsqrt( (gam*Q(IPR,i) + Bx**2 + Q(IBY,i)**2 + Q(IBZ,i)**2)/Q(IDN,i))
         dtlocal =(xf(i+1)-xf(i))/(abs(Q(IVX,i)) + cf)
         if(dtlocal .lt. dtmin) dtmin = dtlocal
    enddo

    TimestepControl = cfl_number*dtmin

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
subroutine NumericalFlux(dt, xv, Q, F)
use params, only : nxtot, NVAR, is, ie
implicit none
real(8), intent(in) :: dt
real(8), intent(in) :: xv(nxtot)
real(8), intent(in) :: Q(NVAR,nxtot)
real(8), intent(out) :: F(NVAR,nxtot)
real(8) :: Ql(NVAR,nxtot),Qr(NVAR,nxtot)
real(8) :: flx(NVAR)
integer :: i,ihy
real(8) :: dQ

      do i=is-1,ie+1
      do ihy=1,NVAR
         dQ = 0.0d0

         Ql(ihy,i+1) = Q(ihy,i) + 0.5d0*dQ
         Qr(ihy,i  ) = Q(ihy,i) - 0.5d0*dQ
      enddo
      enddo

      do i=is,ie+1
!         call Lax((xv(i) - xv(i-1))/dt,Ql(:,i),Qr(:,i),flx(:))
         call HLL(Ql(:,i),Qr(:,i),flx(:))
!         call HLLD(Ql(:,i),Qr(:,i),flx(:))
         do ihy=1,NVAR
            F(ihy,i)  = flx(ihy)
         enddo
      enddo


return
end subroutine Numericalflux
!=============================================================
! Lax
! Description:
!   Compute the numerical flux at a cell interface using the
!   Lax (Lax–Friedrichs / Rusanov-type) flux formula for the 1D MHD equations.
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
use params, only : IDN, IVX, IVY, IVZ, IPR, IBY, IBZ, &
                   IMX, IMY, IMZ, IEN, NVAR, Bx, gam
implicit none
real(8),intent(in)  :: Ql(NVAR), Qr(NVAR)
real(8),intent(in)  :: dxdt
real(8),intent(out) :: flx(NVAR)
integer :: i, n
real(8) :: Ul(NVAR), Ur(NVAR)
real(8) :: Fl(NVAR), Fr(NVAR)
real(8) :: ptotl, ptotr, pbl, pbr

    pbl = 0.5d0*(Bx**2 + Ql(IBY)**2 + Ql(IBZ)**2)
    pbr = 0.5d0*(Bx**2 + Qr(IBY)**2 + Qr(IBZ)**2)
    ptotl = Ql(IPR) + pbl
    ptotr = Qr(IPR) + pbr

    ! conserved variables in the left and right states
    Ul(IDN) = Ql(IDN)
    Ul(IMX) = Ql(IDN)*Ql(IVX)
    Ul(IMY) = Ql(IDN)*Ql(IVY)
    Ul(IMZ) = Ql(IDN)*Ql(IVZ)
    Ul(IEN) = 0.5d0*Ql(IDN)*( Ql(IVX)**2 + Ql(IVY)**2 + Ql(IVZ)**2) & 
                  + pbl + Ql(IPR)/(gam - 1.0d0)
    Ul(IBY) = Ql(IBY)
    Ul(IBZ) = Ql(IBZ)

    Ur(IDN) = Qr(IDN)
    Ur(IMX) = Qr(IDN)*Qr(IVX)
    Ur(IMY) = Qr(IDN)*Qr(IVY)
    Ur(IMZ) = Qr(IDN)*Qr(IVZ)
    Ur(IEN) = 0.5d0*Qr(IDN)*( Qr(IVX)**2 + Qr(IVY)**2 + Qr(IVZ)**2) & 
                  + pbr + Qr(IPR)/(gam - 1.0d0)
    Ur(IBY) = Qr(IBY)
    Ur(IBZ) = Qr(IBZ)

    ! flux in the left and right states
    Fl(IDN) = Ul(IMX)
    Fl(IMX) = Ul(IMX)*Ql(IVX) + ptotl - Bx**2
    Fl(IMY) = Ul(IMY)*Ql(IVX) - Bx*Ql(IBY)
    Fl(IMZ) = Ul(IMZ)*Ql(IVX) - Bx*Ql(IBZ)
    Fl(IEN) = ( Ul(IEN) + ptotl )*Ql(IVX) &
            - Bx*( Bx*Ql(IVX) + Ql(IBY)*Ql(IVY) + Ql(IBZ)*Ql(IVZ) )
    Fl(IBY) = Ql(IBY)*Ql(IVX) - Bx*Ql(IVY)
    Fl(IBZ) = Ql(IBZ)*Ql(IVX) - Bx*Ql(IVZ)

    Fr(IDN) = Ur(IMX)
    Fr(IMX) = Ur(IMX)*Qr(IVX) + ptotr - Bx**2
    Fr(IMY) = Ur(IMY)*Qr(IVX) - Bx*Qr(IBY)
    Fr(IMZ) = Ur(IMZ)*Qr(IVX) - Bx*Qr(IBZ)
    Fr(IEN) = ( Ur(IEN) + ptotr )*Qr(IVX) &
                  - Bx*( Bx*Qr(IVX) + Qr(IBY)*Qr(IVY) + Qr(IBZ)*Qr(IVZ) )
    Fr(IBY) = Qr(IBY)*Qr(IVX) - Bx*Qr(IVY)
    Fr(IBZ) = Qr(IBZ)*Qr(IVX) - Bx*Qr(IVZ)

    do n=1,NVAR 
        flx(n)  = 0.5d0*(Fl(n) + Fr(n)) - 0.5d0*dxdt*(Ur(n) - Ul(n))
    enddo


return
end subroutine Lax
!=============================================================
! HLL
! Description:
!   Compute the interface flux using the HLL approximate Riemann solver for
!   the 1D MHD equations.
!
! Inputs:
!   Ql(:)  Left primitive state  (rho, v, p)
!   Qr(:)  Right primitive state (rho, v, p)
!
! Output:
!   flx(:) Conservative flux (mass, momentum, energy)
!=============================================================
subroutine HLL(Ql,Qr,flx)
use params, only : IDN, IVX, IVY, IVZ, IPR, IBY, IBZ, &
                   IMX, IMY, IMZ, IEN, NVAR, Bx, gam
implicit none
real(8),intent(in) :: Ql(NVAR), Qr(NVAR)
real(8),intent(out):: flx(NVAR)
real(8):: Ul(NVAR), Ur(NVAR)
real(8):: Fl(NVAR), Fr(NVAR)



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
integer::i,ihy

      do i=is,ie
      do ihy=1,NVAR
         U(ihy,i) = U(ihy,i) + dt*(- F(ihy,i+1) + F(ihy,i))/(xf(i+1)-xf(i)) 
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
use params, only : nxtot, NVAR, IDN, IVX, IVY, IVZ, IPR, Bx, &
                   IBY, IBZ, dtsnap, unitsnap, is, ie, dirname
implicit none
real(8), intent(in) :: time
logical, intent(in) :: flag
real(8), intent(in) :: xv(nxtot), Q(NVAR,nxtot)
integer :: i
character(40) :: filename
real(8), save :: tsnap = - dtsnap
integer, save :: nsnap = 0

    if( .not.flag) then
          if( time + 1.0d-14.lt. tsnap+dtsnap) return
    endif

    write(filename,'(i5.5)') nsnap
    filename = trim(dirname)//"/snap"//trim(filename)//".dat"
    open(unitsnap,file=filename,form='formatted',action="write")
    write(unitsnap,"(a2,f6.4)") "# ",time
    do i=is,ie
          write(unitsnap,'(1p,9(es24.16e3,1x))') xv(i), Q(IDN,i), Q(IVX,i), Q(IVY,i), Q(IVZ,i), & 
                                               Q(IPR,i), Bx, Q(IBY,i), Q(IBZ,i)
    enddo
    close(unitsnap)

    print*, "output time=", time, trim(filename)

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
!   xv(:)    Cell-center coordinates
!   Q(:,:)   Primitive variables Q=(rho, vx, vy, vz, p, By, Bz)
!=============================================================
subroutine RealtimeAnalysis(xv,Q,phys_evo)
use params, only : IDN, IVX, IVY, IBZ, IPR, Bx, IBY, IBZ, NVAR, &
                    is, ie, gam, nx, nxtot, nevo
implicit none
real(8), intent(in)  :: xv(nxtot), Q(NVAR,nxtot)
real(8), intent(out) :: phys_evo(nevo)
integer :: i,j,k
real(8) :: tmp

      phys_evo(1) = 0.0d0
      
return
end subroutine
!=============================================================
! AnalysisAfterSimu
! Description:
!   Perform post-processing after the time integration has finished.
!   This routine is intended to be called once at the end of the run to produce
!   summary diagnostics and/or derived outputs.
!
! Inputs:
!   xv(:)    Cell-center coordinates
!   Q(:,:)   Primitive variables Q=(rho, vx, vy, vz, p, By, Bz)
!=============================================================
subroutine AnalysisAfterSimu(time,xv,Q)
use params, only : IDN, IVX, IVY, IBZ, IPR, Bx, IBY, IBZ, NVAR, &
                    is, ie, gam, nx, nxtot, nevo
implicit none
real(8), intent(in)  :: xv(nxtot), Q(NVAR,nxtot)
real(8), intent(in)  :: time
integer :: i
real(8) :: error

      error = 0.0d0
      do i=is,ie
           error = error + 1.0d0
      enddo

      open(1,file="error.dat",action="write",position="append")
      write(1,*) nx, error
      print*,"nx = ",nx, "error = ",error
      close(1)
      
return
end subroutine
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

