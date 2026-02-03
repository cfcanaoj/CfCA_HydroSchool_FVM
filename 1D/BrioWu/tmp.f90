program main
implicit none

! time evolution
integer :: ntime = 0    ! counter of the timestep
real(8) :: time = 0.0d0  ! time 
real(8) :: dt   = 0.0d0  ! time width
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
integer, parameter :: NFLX = 7

! indices of the conservative variables
integer, parameter :: IVX = 2
integer, parameter :: IVY = 3
integer, parameter :: IVZ = 4
integer, parameter :: IEN = 5

! adiabatic index
real(8),parameter::gam=5.0d0/3.0d0 


! definition of arrays 
real(8),dimension(nxtot)::xf,xv
real(8),dimension(nxtot,NVAR) :: U
real(8),dimension(nxtot,NVAR) :: Q
real(8),dimension(nxtot,NFLX) :: flux

! output 
character(20),parameter::dirname="lax" ! directory name

! snapshot
integer, parameter :: unitsnap = 17

! realtime analysis 
!integer, parameter :: unitevo = 11
!integer, parameter :: nevo = 3    ! the number of variables derived in the realtime analysis
!real(8) ::  phys_evo(nevo)    ! variables derived in the realtime analysis

     ! make the directory for output
     call makedirs(trim(dirname))

     call GenerateGrid(xf, xv)
     call GenerateProblem(xv, Q)
     call PrIMYConsv(Q, U)
     call Output( .TRUE., dirname, xv, Q )

      write(6,*) "Start the simulation"
      !open file to output the result of the realtime analysis
!      open(unitevo,file=trim(dirname)//'/'//'ana.dat', action="write")
! main loop
      mloop: do 
         dt = TimestepControl(xf, Q)
         if( time + dt > timemax ) dt = timemax - time
         call BoundaryCondition( Q)
         call NumericalFlux( Q, flux )
         call UpdateConsv( dt, xf, flux, U)
         call Consv2Prim( U, Q )
         time=time+dt
         print*,"time = ",time, "dt = ",dt
         ntime = ntime + 1
         call Output( .FALSE., dirname, xv, Q )

!         if( mod(ntime,10) .eq. 0 ) then
!            call RealtimeAnalysis(xf,xv,Q,phys_evo)
!            write(unitevo,*) time, phys_evo(1:nevo)
!         endif


         if(time >= timemax) exit mloop
      enddo mloop
      call Output( .TRUE., dirname, xv, Q )
!      call AnalysisAfterSimu(time,xf,xv,Q)

!      write(6,*) "program has been finished"
contains

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
implicit none
real(8), intent(out) :: xf(:), xv(:)
real(8) :: dx
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
!   Set initial conditions for a 1D Riemann problem (Brio-Wu shock tube).
!   The primitive variables Q(i,:) = (rho, v, p) are assigned based on xv(i):
!     - left  state (x < 0):   rho = 1.0,   v = (0,0,0), p = 1.0, B = (0.75,1,0)
!     - right state (x > 0):   rho = 0.125, v = (0,0,0), p = 0.1, B = (0.75,-1,0)
!   The routine typically initializes only the active zone (i=is:ie).
!   Ghost zones are filled later by BoundaryCondition().
!=============================================================
subroutine GenerateProblem(xv, Q)
implicit none
integer::i
real(8), intent(in ) :: xv(:)
real(8), intent(out) :: Q(:,:)

    do i=is,ie
        if( xv(i) < 0.0d0 ) then 
             Q(i,IDN) = 1.0d0
             Q(i,IVX) = 0.0d0
             Q(i,IVY) = 0.0d0
             Q(i,IVZ) = 0.0d0
             Q(i,IPR) = 1.0d0
             Q(i,IBY) = 1.0d0
             Q(i,IBZ) = 0.0d0
        else 
             Q(i,IDN) = 0.125d0
             Q(i,IVX) = 0.0d0
             Q(i,IVY) = 0.0d0
             Q(i,IVZ) = 0.0d0
             Q(i,IPR) = 0.1d0
             Q(i,IBY) = -1.0d0
             Q(i,IBZ) = 0.0d0
         endif
    enddo

      
!      call BoundaryCondition

return
end subroutine GenerateProblem
!=============================================================
! BoundaryCondition
! Description:
!   Apply boundary conditions by filling ghost cells ( i < is and i > ie ) 
!   of the primitive array Q.
!
!=============================================================
subroutine BoundaryCondition(Q)
implicit none
real(8), intent(inout) :: Q(:,:)
integer::i,ihy


 
      do ihy=1,NVAR
      do i=1,ngh
          Q(is-i,ihy)  = Q(is-1+i,ihy)
      enddo
      enddo

      do ihy=1,NVAR
      do i=1,ngh
          Q(ie+i,ihy) = Q(ie-i+1,ihy)
      enddo
      enddo

return
end subroutine BoundaryCondition
!=============================================================
! Prim2Consv
! Description:
!   Convert primitive variables Q = (rho, v, p) to conservative variables
!   U = (rho, mom, E) for the 1D Euler equations with an ideal-gas EOS.
!     - mom = rho * v
!     - E   = 0.5 * rho * v^2 + p / (gam - 1) + 0.5 * B^2
!   Operates on the active zone (i=is:ie).
!=============================================================
subroutine PrIMYConsv(Q, U)
implicit none
real(8), intent(in) :: Q(:,:)
real(8), intent(out) :: U(:,:)
integer::i

      do i=is,ie
          U(i,IDN) = Q(i,IDN)
          U(i,IMX) = Q(i,IDN)*Q(i,IVX)
          U(i,IMY) = Q(i,IDN)*Q(i,IVY)
          U(i,IMZ) = Q(i,IDN)*Q(i,IVZ)
          U(i,IEN) = 0.5d0*Q(i,IDN)*( Q(i,IVX)**2 + Q(i,IVY)**2 + Q(i,IVZ)**2 ) &
                   + 0.5d0*( Bx**2 + Q(i,IBY)**2 + Q(i,IBZ)**2 ) &
                   + Q(i,IPR)/(gam - 1.0d0)
          U(i,IBY) = Q(i,IBY)
          U(i,IBZ) = Q(i,IBZ)
      enddo
      
return
end subroutine PrIMYConsv
!=============================================================
! Consv2Prim
! Description:
!   Convert conservative variables U = (rho, mom, E) to primitive variables
!   Q = (rho, v, p) for the 1D Euler equations with an ideal-gas EOS.
!     - v = mom / rho
!     - p = (E - 0.5 * mom^2 / rho - 0.5 * B^2) * (gam - 1)
!   Operates on the active zone (i=is:ie).
!=============================================================
subroutine Consv2Prim( U, Q )
implicit none
real(8), intent(in) :: U(:,:)
real(8), intent(out) :: Q(:,:)
integer::i
real(8) :: inv_d;

    do i=is,ie
        Q(i,IDN) = U(i,IDN)
        inv_d = 1.0d0/U(i,IDN)
        Q(i,IVX) = U(i,IMX)*inv_d
        Q(i,IVY) = U(i,IMY)*inv_d
        Q(i,IVZ) = U(i,IMZ)*inv_d
        Q(i,IPR) = ( U(i,IEN) &
                    - 0.5d0*(U(i,IMX)**2 + U(i,IMY)**2 + U(i,IMZ)**2)*inv_d  &
                    - 0.5d0*(Bx**2 + U(i,IBY)**2 + U(i,IBZ)**2) )*(gam-1.0d0)
        Q(i,IBY) = U(i,IBY)
        Q(i,IBZ) = U(i,IBZ)
    enddo

return
end subroutine Consv2Prim
!=============================================================
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
implicit none
real(8), intent(in) :: xf(:), Q(:,:)
real(8)::dtlocal
real(8)::dtmin,cf
integer::i

    dtmin=1.0d90

    do i=is,ie
         cf = dsqrt( (gam*Q(i,IPR) + Bx**2 + Q(i,IBY)**2 + Q(i,IBZ)**2)/Q(i,IDN))
         dtlocal =(xf(i+1)-xf(i))/(abs(Q(i,IVX)) + cf)
         if(dtlocal .lt. dtmin) dtmin = dtlocal
    enddo

      TimestepControl = 0.1d0 * dtmin
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
subroutine NumericalFlux( Q, flux )
implicit none
real(8), intent(in) :: Q(:,:)
real(8), intent(out) :: flux(:,:)
real(8),dimension(nxtot,NFLX):: Ql,Qr
real(8),dimension(NFLX):: flx
integer::i,ihy
real(8) :: dQm, dQp, dQ

      do ihy=1,NVAR
      do i=is-1,ie+1
         dQp = Q(i+1,ihy) - Q(i  ,ihy)
         dQm = Q(i  ,ihy) - Q(i-1,ihy)

         if(dQp*dQm .gt. 0.0d0)then
            dQ = 2.0d0*dQp*dQm/(dQp+dQm)
         else
            dQ = 0.0d0
         endif

         Ql(i+1,ihy) = Q(i,ihy) + 0.5d0*dQ
         Qr(i  ,ihy) = Q(i,ihy) - 0.5d0*dQ
      enddo
      enddo

      do i=is,ie+1
         call Lax((xv(i) - xv(i-1))/dt,Ql(i,:),Qr(i,:),flx(:))
!         call HLL(Ql(i,:),Qr(i,:),flx)
!         call HLLD(Ql(i,:),Qr(i,:),flx)
         flux(i,:)  = flx(:)
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
!   Ql(:)    Left primitive state  (rho, v, p, B)
!   Qr(:)    Right primitive state (rho, v, p, B)
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
implicit none
real(8),intent(in)  :: Ql(:), Qr(:)
real(8),intent(in)  :: dxdt
real(8),intent(out) :: flx(:)
integer :: n
real(8):: Ul(NVAR), Ur(NVAR)
real(8):: Fl(NVAR), Fr(NVAR)
real(8):: pbl,pbr,ptotl,ptotr

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
implicit none
real(8),intent(in)::Ql(:), Qr(:)
real(8),intent(out) :: flx(:)
real(8):: Ul(NFLX), Ur(NFLX)
real(8):: Fl(NFLX), Fr(NFLX)
real(8):: cfl,cfr
real(8):: sl, sr
real(8):: pbl, pbr, ptotl, ptotr

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

!    cfl = dsqrt( (gam*Ql(IPR) + Ql(IBY)**2 + Ql(IBZ)**2 + Bx**2)/Ql(IDN))
!    cfr = dsqrt( (gam*Qr(IPR) + Qr(IBY)**2 + Qr(IBZ)**2 + Bx**2)/Qr(IDN))
    cfl = sqrt( 0.5d0*( 2.0d0*pbl + gam*Ql(IPR) &
                     + dsqrt( (2.0d0*pbl - gam*Ql(IPR))**2 &
                     + 4.0d0*gam*Ql(IPR)*( Ql(IBY)**2 + Ql(IBZ)**2 ) ) )/Ql(IDN) )
    cfr = sqrt( 0.5d0*( 2.0d0*pbr + gam*Qr(IPR) &
                     + dsqrt( (2.0d0*pbr - gam*Qr(IPR))**2 &
                     + 4.0d0*gam*Qr(IPR)*( Qr(IBY)**2 + Qr(IBZ)**2 ) ) )/Qr(IDN) )

   sl = min(Ql(IVX),Qr(IVX)) - max(cfl,cfr)
   sr = max(Ql(IVX),Qr(IVX)) + max(cfl,cfr)
!    sl = min(Ql(IVX) - cfl,Qr(IVX) - cfr)
!    sr = max(Ql(IVX) + cfl,Qr(IVX) + cfr)


    if( sl > 0.0d0 ) then
        flx(:)  = Fl(:) 
    else if (sr <= 0.0d0 ) then
        flx(:)  = Fr(:) 
    else 
        flx(:)  = (sr*Fl(:) - sl*Fr(:) + sl*sr*( Ur(:) - Ul(:) ))/(sr - sl)
    endif

return
end subroutine HLL
!=============================================================
! HLLD
! Description:
!   Compute the interface flux using the HLLD approximate Riemann solver for
!   the 1D Euler equations.
!
! Inputs:
!   Ql(:)  Left primitive state  (rho, v, p)
!   Qr(:)  Right primitive state (rho, v, p)
!
! Output:
!   flx(:) Conservative flux (mass, momentum, energy)
!=============================================================
subroutine HLLD(Ql,Qr,flx)
implicit none
real(8),intent(in)  ::Ql(:), Qr(:)
real(8),intent(out) :: flx(:)
real(8):: Ul(NFLX), Ur(NFLX)
real(8):: Ulst(NFLX), Urst(NFLX)
real(8):: Uldst(NFLX), Urdst(NFLX)
real(8):: Fl(NFLX), Fr(NFLX)
real(8):: cfl,cfr
real(8):: S0, S1, S2, S3, S4
real(8):: pbl, pbr, ptotl, ptotr
real(8) :: sqrtdl, sqrtdr, v_dot_B_stl, v_dot_B_str
real(8) :: Ulst_d_inv, Urst_d_inv, sum_sqrtd_inv, tmp
real(8) :: ptot_stl, ptot_str,ptot_st, Cl, Cr, Cml, Cmr, Cml_inv, Cmr_inv, bxsgn

    pbl = 0.5d0*(Bx**2 + Ql(IBY)**2 + Ql(IBZ)**2)
    pbr = 0.5d0*(Bx**2 + Qr(IBY)**2 + Qr(IBZ)**2)
    ptotl = Ql(IPR) + pbl
    ptotr = Qr(IPR) + pbr

    cfl = sqrt( 0.5d0*( 2.0d0*pbl + gam*Ql(IPR) &
                     + dsqrt( (2.0d0*pbl - gam*Ql(IPR))**2 &
                     + 4.0d0*gam*Ql(IPR)*( Ql(IBY)**2 + Ql(IBZ)**2 ) ) )/Ql(IDN) )
    cfr = sqrt( 0.5d0*( 2.0d0*pbr + gam*Qr(IPR) &
                     + dsqrt( (2.0d0*pbr - gam*Qr(IPR))**2 &
                     + 4.0d0*gam*Qr(IPR)*( Qr(IBY)**2 + Qr(IBZ)**2 ) ) )/Qr(IDN) )

    S0 = min( Ql(IVX) - cfl, Qr(IVX) - cfr)
    S4 = max( Ql(IVX) + cfl, Qr(IVX) + cfr)

          ! conserved variables in the left and right states
    Ul(IDN) = Ql(IDN)
    Ul(IVX) = Ql(IDN)*Ql(IVX)
    Ul(IVY) = Ql(IDN)*Ql(IVY)
    Ul(IVZ) = Ql(IDN)*Ql(IVZ)
    Ul(IEN) = 0.5d0*Ql(IDN)*( Ql(IVX)**2 + Ql(IVY)**2 + Ql(IVZ)**2) & 
            + pbl + Ql(IPR)/(gam - 1.0d0)
    Ul(IBY) = Ql(IBY)
    Ul(IBZ) = Ql(IBZ)

    Ur(IDN) = Qr(IDN)
    Ur(IVX) = Qr(IDN)*Qr(IVX)
    Ur(IVY) = Qr(IDN)*Qr(IVY)
    Ur(IVZ) = Qr(IDN)*Qr(IVZ)
    Ur(IEN) = 0.5d0*Qr(IDN)*( Qr(IVX)**2 + Qr(IVY)**2 + Qr(IVZ)**2) & 
            + pbr + Qr(IPR)/(gam - 1.0d0)
    Ur(IBY) = Qr(IBY)
    Ur(IBZ) = Qr(IBZ)

    !--- Step 3.  Compute L/R fluxes
    Fl(IDN) = Ul(IVX)
    Fl(IVX) = Ul(IVX)*Ql(IVX) + ptotl - Bx**2
    Fl(IVY) = Ul(IVY)*Ql(IVX) - Bx*Ql(IBY)
    Fl(IVZ) = Ul(IVZ)*Ql(IVX) - Bx*Ql(IBZ)
    Fl(IEN) = ( Ul(IEN) + ptotl - Bx**2 )*Ql(IVX) &
            - Bx*( Ql(IBY)*Ql(IVY) + Ql(IBZ)*Ql(IVZ) )
    Fl(IBY) = Ql(IBY)*Ql(IVX) - Bx*Ql(IVY)
    Fl(IBZ) = Ql(IBZ)*Ql(IVX) - Bx*Ql(IVZ)

    Fr(IDN) = Ur(IVX)
    Fr(IVX) = Ur(IVX)*Qr(IVX) + ptotr - Bx**2
    Fr(IVY) = Ur(IVY)*Qr(IVX) - Bx*Qr(IBY)
    Fr(IVZ) = Ur(IVZ)*Qr(IVX) - Bx*Qr(IBZ)
    Fr(IEN) = ( Ur(IEN) + ptotr - Bx**2 )*Qr(IVX) &
            - Bx*( Qr(IBY)*Qr(IVY) + Qr(IBZ)*Qr(IVZ) )
    Fr(IBY) = Qr(IBY)*Qr(IVX) - Bx*Qr(IVY)
    Fr(IBZ) = Qr(IBZ)*Qr(IVX) - Bx*Qr(IVZ)

    !--- Step 4.  Compute middle and Alfven wave speeds
    Cl = S0 - Ql(IVX)
    Cr = S4 - Qr(IVX)

    S2 = ( Cr*Ur(IVX) - Cl*Ul(IVX) + (ptotl - ptotr) ) &
           /( Cr*Ur(IDN) - Cl*Ul(IDN) )

    Cml = S0 - S2
    Cmr = S4 - S2
    Cml_inv = 1.0d0/Cml
    Cmr_inv = 1.0d0/Cmr

    Ulst(IDN) = Ul(IDN)*Cl*Cml_inv
    Urst(IDN) = Ur(IDN)*Cr*Cmr_inv
    Ulst_d_inv = 1.0d0/Ulst(IDN)
    Urst_d_inv = 1.0d0/Urst(IDN)
    sqrtdl = dsqrt(Ulst(IDN))
    sqrtdr = dsqrt(Urst(IDN))

!          if( sqrtdr .ne. sqrtdr ) then
!              print*, "sqrtdr",sqrtdr, Cr, Cmr
!             print*,"S", S0,S2,S4
!              stop
!          endif

    S1 = S2 - dabs(Bx)/sqrtdl
    S3 = S2 + dabs(Bx)/sqrtdr

    !--- Step 5.  Compute intermediate states
   ptot_stl = ptotl + Ul(IDN)*Cl*(S2 - Ql(IVX))
   ptot_str = ptotr + Ur(IDN)*Cr*(S2 - Qr(IVX))

   ptot_st = 0.5d0*(ptot_stl + ptot_str)

   Ulst(IVX) = Ulst(IDN)*S2
   if( dabs( Ul(IDN)*Cl*Cml-Bx**2) < 1.0d-8*ptot_st ) then
       Ulst(IVY) = Ulst(IDN)*Ql(IVY)
       Ulst(IVZ) = Ulst(IDN)*Ql(IVZ)

       Ulst(IBY) = Ul(IBY)
       Ulst(IBZ) = Ul(IBZ)
   else 
       tmp = Bx*( Cl - Cml )/(Ul(IDN)*Cl*Cml - Bx**2)
       Ulst(IVY) = Ulst(IDN)*( Ql(IVY) - Ul(IBY)*tmp )
       Ulst(IVZ) = Ulst(IDN)*( Ql(IVZ) - Ul(IBZ)*tmp )

       tmp = (Ul(IDN)*Cl**2 - Bx**2)/( Ul(IDN)*Cl*Cml - Bx**2)
       Ulst(IBY) = Ul(IBY)*tmp
       Ulst(IBZ) = Ul(IBZ)*tmp
   endif

   v_dot_B_stl = ( Ulst(IVX)*Bx + Ulst(IVY)*Ulst(IBY) + Ulst(IVZ)*Ulst(IBZ) )*Ulst_d_inv
   Ulst(IEN) = ( Cl*Ul(IEN) - ptotl*Ql(IVX) + ptot_st*S2 &
               + Bx*( Ql(IVX)*Bx + Ql(IVY)*Ul(IBY) + Ql(IVZ)*Ul(IBZ) - v_dot_B_stl) )*Cml_inv

   Urst(IVX) = Urst(IDN)*S2
   if( dabs( Ur(IDN)*Cr*Cmr-Bx**2) < 1.0d-8*ptot_st ) then
       Urst(IVY) = Urst(IDN)*Qr(IVY)
       Urst(IVZ) = Urst(IDN)*Qr(IVZ)

       Urst(IBY) = Ur(IBY)
       Urst(IBZ) = Ur(IBZ)
   else 
       tmp = Bx*( Cr - Cmr )/(Ur(IDN)*Cr*Cmr - Bx**2)
       Urst(IVY) = Urst(IDN)*( Qr(IVY) - Ur(IBY)*tmp )
       Urst(IVZ) = Urst(IDN)*( Qr(IVZ) - Ur(IBZ)*tmp )

       tmp = (Ur(IDN)*Cr**2 - Bx**2)/( Ur(IDN)*Cr*Cmr - Bx**2)
       Urst(IBY) = Ur(IBY)*tmp
       Urst(IBZ) = Ur(IBZ)*tmp
   endif

   v_dot_B_str = ( Urst(IVX)*Bx + Urst(IVY)*Urst(IBY) + Urst(IVZ)*Urst(IBZ) )*Urst_d_inv
   Urst(IEN) = ( Cr*Ur(IEN) - ptotr*Qr(IVX) + ptot_st*S2 &
               + Bx*( Qr(IVX)*Bx + Qr(IVY)*Ur(IBY) + Qr(IVZ)*Ur(IBZ) - v_dot_B_str) )*Cmr_inv
 
   if( 0.5d0*Bx**2 < 1.0d-8*ptot_st )then
       Uldst(:) = Ulst(:)
       Urdst(:) = Urst(:)
   else 
       sum_sqrtd_inv = 1.0d0/(sqrtdl + sqrtdr) 
       if (Bx > 0.0d0 ) then 
           bxsgn = 1.0d0
       else 
           bxsgn = -1.0d0
       endif

       Uldst(IDN) = Ulst(IDN)
       Urdst(IDN) = Urst(IDN)

       Uldst(IVX) = Ulst(IVX)
       Urdst(IVX) = Urst(IVX)

       tmp = sum_sqrtd_inv*(  sqrtdl*(Ulst(IVY)*Ulst_d_inv) + sqrtdr*(Urst(IVY)*Urst_d_inv) &
                            + bxsgn*(Urst(IBY) - Ulst(IBY)) )
       Uldst(IVY) = Uldst(IDN)*tmp
       Urdst(IVY) = Urdst(IDN)*tmp
!

       tmp = sum_sqrtd_inv*(  sqrtdl*(Ulst(IVZ)*Ulst_d_inv) + sqrtdr*(Urst(IVZ)*Urst_d_inv) &
                            + bxsgn*(Urst(IBZ) - Ulst(IBZ)) )
       Uldst(IVZ) = Uldst(IDN)*tmp
       Urdst(IVZ) = Urdst(IDN)*tmp

       tmp = sum_sqrtd_inv*(  sqrtdl*Urst(IBY) + sqrtdr*Ulst(IBY) &
                 + bxsgn*sqrtdl*sqrtdr*( (Urst(IVY)*Urst_d_inv) - (Ulst(IVY)*Ulst_d_inv) ) )
       Uldst(IBY) = tmp
       Urdst(IBY) = tmp

       tmp = sum_sqrtd_inv*(  sqrtdl*Urst(IBZ) + sqrtdr*Ulst(IBZ) &
                 + bxsgn*sqrtdl*sqrtdr*( (Urst(IVZ)*Urst_d_inv) - (Ulst(IVZ)*Ulst_d_inv) ) )
       Uldst(IBZ) = tmp
       Urdst(IBZ) = tmp
!
       tmp = S2*Bx + (Uldst(IVY)*Uldst(IBY) + Uldst(IVZ)*Uldst(IBZ))/Uldst(IDN)
       Uldst(IEN) = Ulst(IEN) - sqrtdl*bxsgn*(v_dot_B_stl - tmp)
       Urdst(IEN) = Urst(IEN) + sqrtdr*bxsgn*(v_dot_B_str - tmp)
   endif

    !--- Step 6.  Compute flux
    if( S0 >= 0.0d0 ) then
         flx(:) = Fl(:)
    else if (S4 <= 0.0d0 ) then
         flx(:) = Fr(:)
    else  if  (S1 >= 0.0d0) then
         flx(:) = Fl(:) + S0*(Ulst(:) - Ul(:))
     else if (S2 >= 0.0d0) then
         flx(:) = Fl(:) + S0*(Ulst(:) - Ul(:)) + S1*(Uldst(:) - Ulst(:))
     else if (S3 > 0.0d0 ) then
         flx(:) = Fr(:) + S4*(Urst(:) - Ur(:)) + S3*(Urdst(:) - Urst(:))
     else 
         flx(:) = Fr(:) + S4*(Urst(:) - Ur(:)) 
     endif

return
end subroutine HLLD
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
subroutine UpdateConsv( dt, xf, flux, U )
implicit none
real(8), intent(in)  :: flux(:,:), dt, xf(:)
real(8), intent(out) :: U(:,:)
integer::i,n

      do n=1,NVAR
      do i=is,ie
         U(i,n) = U(i,n) + dt*(- flux(i+1,n) + flux(i,n))/(xf(i+1)-xf(i)) 
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
!   dirname  Output directory name
!   xv(:)    Cell-center coordinates
!   Q(:,:)   Primitive variables to be written
!
! Notes:
!   With variable dt, consider using a "while (time >= next_output_time)" style
!   to avoid missing outputs when the simulation time jumps over an output time.
!=============================================================
subroutine Output( flag, dirname, xv, Q )
implicit none
logical,       intent(in) :: flag
character(20), intent(in) :: dirname 
real(8),       intent(in) :: xv(:), Q(:,:)
real(8), parameter:: dtsnap=5.0d-3
integer::i
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
          write(unitsnap,*) xv(i), Q(i,IDN), Q(i,IVX), Q(i,IVY), Q(i,IVZ), Q(i,IPR), Bx, &
          Q(i,IBY), Q(i,IBZ)
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
!   Q(:,:)   Primitive variables Q=(rho, v, p)
!   U(:,:)   Conservative variables U=(rho, mom, E) (optional but useful)
!=============================================================
!subroutine RealtimeAnalysis(xf,xv,Q,phys_evo)
!real(8), intent(in)  :: xf(:), xv(:), Q(:,:)
!real(8), intent(out) :: phys_evo(:)
!integer :: i
!real(8) :: tmp
!
!      tmp = 0.0d0
!      do i=is,ie
!           tmp = tmp + 1.0d0
!      enddo
!      phys_evo(1) = tmp/dble(nx)
!      
!return
!end subroutine
!=============================================================
! AnalysisAfterSimu
! Description:
!   Perform post-processing after the time integration has finished.
!   This routine is intended to be called once at the end of the run to produce
!   summary diagnostics and/or derived outputs.
!
! Inputs:
!   dirname  Output directory name where snapshots/logs are stored
!   xv(:)    Cell-center coordinates
!   Q(:,:)   Final primitive variables Q=(rho, v, p)
!   U(:,:)   Final conservative variables U=(rho, mom, E)
!=============================================================
!subroutine AnalysisAfterSimu(time,xf,xv,Q)
!real(8), intent(in)  :: xf(:), xv(:), Q(:,:)
!real(8), intent(in)  :: time
!integer :: i
!real(8) :: error
!
!      error = 0.0d0
!      do i=is,ie
!           error = error + 1.0d0
!      enddo
!      print*, nx, error
!      
!return
!end subroutine
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
character(len=*), intent(in) :: outdir
character(len=256) command
write(command, *) 'if [ ! -d ', trim(outdir), ' ]; then mkdir -p ', trim(outdir), '; fi'
!      write(*, *) trim(command)
      call system(command)
end subroutine makedirs

end program main
