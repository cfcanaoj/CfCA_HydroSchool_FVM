program main
  use omp_lib
use params, only : nxtot, nytot, NVAR, dirname, unitevo, timemax, nevo
implicit none

! time evolution
integer :: ntime = 0    ! counter of the timestep
real(8) :: time = 0.0d0  ! time 
real(8) :: dt   = 0.0d0  ! time width

! definition of arrays 
real(8),dimension(nxtot)::xf,xv
real(8),dimension(nytot)::yf,yv
real(8),dimension(NVAR,nxtot,nytot) :: Uo
real(8),dimension(NVAR,nxtot,nytot) :: U
real(8),dimension(NVAR,nxtot,nytot) :: Q
real(8),dimension(NVAR,nxtot,nytot) :: F
real(8),dimension(NVAR,nxtot,nytot) :: G

! realtime analysis
real(8) :: phys_evo(nevo)

 real(8)::time_begin,time_end
 logical,parameter:: benchmarkmode=.true.
      ! make the directory for output
      call makedirs(trim(dirname))
!$acc data create(dt,xf,xv,yf,yv, Uo, U, Q, F,G)
      write(6,*) "setup grids and initial condition"
      call GenerateGrid(xf, xv, yf, yv)
      call GenerateProblem(xv, yv, Q )
      call Prim2Consv(Q, U)
      call BoundaryCondition(Q)
      call Output( time, .TRUE., xv, yv, Q )

      write(6,*) "Start the simulation"
      open(unitevo,file=trim(dirname)//'/'//'ana.dat', action="write")
! main loop
      ntime = 1
      time_begin = omp_get_wtime()
      mloop: do !ntime=1,ntimemax
         call TimestepControl(xf, yf, Q, dt)
         if( time + dt > timemax ) dt = timemax - time

         call SaveState(U,Uo)
         call NumericalFlux( dt, xf, yf, Q, F, G )
         call UpdateConsv( 0.5d0*dt, xf, yf, F, G, Uo, U )
         call SrcTerms( 0.5d0*dt, dt, Q, U)
         call Consv2Prim( U, Q )
         call BoundaryCondition(Q )
!
         call NumericalFlux( dt, xf, yf, Q, F, G )
         call UpdateConsv( dt, xf, yf, F, G, Uo, U )
         call SrcTerms( dt, dt, Q, U)
         call Consv2Prim( U, Q )
         call BoundaryCondition( Q )

         time=time+dt
         ntime = ntime+1
         
         if(.not. benchmarkmode) call Output( time, .FALSE., xv, yv, Q)

         if(.not. benchmarkmode) print*, "ntime = ",ntime, "time = ",time, dt

         if( mod(ntime,10) .eq. 0 ) then
            if(.not. benchmarkmode) call RealtimeAnalysis(xv,yv,Q,phys_evo)
            if(.not. benchmarkmode) write(unitevo,*) time, phys_evo(1:nevo)
         endif

         if(time >= timemax) exit mloop
      enddo mloop
      time_end = omp_get_wtime()
      print *, "sim time [s]:", time_end-time_begin
      close(unitevo)
      call Output( time, .TRUE.,xv, yv, Q)
!$acc end data
!      write(6,*) "program has been finished"
!contains
end program
!=============================================================
! GenerateGrid
! Description:
!   Generate a 2D uniform grid for a finite-volume scheme.
!   This routine fills:
!     - xf(:), yf(:)  face (cell-boundary) coordinates
!     - xv(:), yv(:)  cell-center coordinates
!   The grid uses global parameters (x1min, x1max, nx, ngh, ...).
!
! Notes:
!   Be careful about array sizes when mixing cell-centered and face-centered
!   quantities. The number of faces is (number of cells + 1).
!=============================================================
subroutine GenerateGrid(xf, xv, yf, yv)
use params, only : nxtot, nytot, ngh, nx, ny, xmax, xmin, ymax, ymin
implicit none
real(8), intent(out) :: xf(nxtot), xv(nxtot)
real(8), intent(out) :: yf(nytot), yv(nytot)
real(8) :: dx,dy
integer::i,j


      dx=(xmax-xmin)/dble(nx)
      do i=1,nxtot
         xf(i) = dx*(i-(ngh+1))+xmin
      enddo
      do i=1,nxtot-1
         xv(i) = 0.5d0*(xf(i+1)+xf(i))
      enddo

      dy=(ymax-ymin)/dble(ny)
      do j=1,nytot
         yf(j) = dy*(j-(ngh+1))+ymin
      enddo
      do j=1,nytot-1
         yv(j) = 0.5d0*(yf(j+1)+yf(j))
      enddo

!$acc update device(xf,xv,yf,yv)
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
subroutine GenerateProblem(xv, yv, Q )
use params, only : IDN, IVX, IVY, IVZ, IPR, IBX, IBY, IBZ, &
                   IMX, IMY, IMZ, IEN, IPS, NVAR, nxtot, nytot, &
                   is, ie, js, je, gam
implicit none
real(8), intent(in ) :: xv(nxtot), yv(nytot)
real(8), intent(out) :: Q(NVAR,nxtot,nytot)
integer::i, j
real(8) :: pi, B0

      pi = dacos(-1.0d0)
      B0 = 10.0d0

      do j=js,je
      do i=is,ie
         Q(IDN,i,j) = 1.0d0
         if( xv(i)**2 + yv(j)**2 < 0.1d0**2 ) then
             Q(IPR,i,j) = 10.0d0
         else 
             Q(IPR,i,j) = 0.1d0
         endif
         Q(IVX,i,j) = 0.0d0
         Q(IVY,i,j) = 0.0d0
         Q(IVZ,i,j) = 0.0d0
         Q(IBX,i,j) = B0/sqrt(2.0d0)
         Q(IBY,i,j) = B0/sqrt(2.0d0)
         Q(IBZ,i,j) = 0.0d0
         Q(IPS,i,j) = 0.0d0
      enddo
      enddo

!$acc update device(Q)
return
end subroutine GenerateProblem
