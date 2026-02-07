module params

real(8), parameter:: timemax=0.15d0 ! simulation end time

! option
integer, parameter :: flag_HDC = 1 ! 1 --> HDC on , 0 --> HDC off
integer, parameter :: flag_flux = 2 ! 1 (HLL), 2 (HLLD)

! coordinate 
integer,parameter::nx=128        ! the number of grids in the simulation box
integer,parameter::ny=128 ! the number of grids in the simulation box
integer,parameter::ngh=2         ! the number of ghost cells
integer,parameter::nxtot=nx+2*ngh+1 ! the total number of grids including ghost cells
integer,parameter::nytot=ny+2*ngh+1 ! the total number of grids including ghost cells
integer,parameter::is=ngh+1         ! the index of the leftmost grid
integer,parameter::js=ngh+1         ! the index of the leftmost grid
integer,parameter::ie=nx+ngh     ! the index of the rightmost grid
integer,parameter::je=ny+ngh     ! the index of the rightmost grid
real(8),parameter::xmin=-0.5d0,xmax=0.5d0
real(8),parameter::ymin=-0.5d0,ymax=0.5d0

real(8),parameter::Ccfl=0.4d0
real(8),parameter::gam=5.0d0/3.0d0 !! adiabatic index

real(8), parameter :: alpha = 0.1d0    ! decay timescale of divergence B

! indices of the conservative variables
integer, parameter :: IDN = 1
integer, parameter :: IMX = 2
integer, parameter :: IMY = 3
integer, parameter :: IMZ = 4
integer, parameter :: IPR = 5
integer, parameter :: IBX = 6
integer, parameter :: IBY = 7
integer, parameter :: IBZ = 8
integer, parameter :: IPS = 9
integer, parameter :: NVAR = 9
integer, parameter :: NFLX = 9

! indices of the primitive variables
integer, parameter :: IVX = 2
integer, parameter :: IVY = 3
integer, parameter :: IVZ = 4
integer, parameter :: IEN = 5

! output 
character(20),parameter::dirname="hdc" ! directory name

! snapshot
integer, parameter :: unitsnap = 17
real(8), parameter:: dtsnap=0.5d-2
logical, parameter :: flag_binary = .false.

! realtime analysis 
integer, parameter :: nevo = 1
integer, parameter :: unitevo =11

end module

program main
!$ use omp_lib
use params, only : nxtot, nytot, NVAR, dirname, unitevo, timemax, nevo
implicit none
include "interfaces_hdc.inc"

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

real(8) :: phys_evo(nevo)

      ! make the directory for output
      call makedirs(trim(dirname))

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
      mloop: do !ntime=1,ntimemax
         dt = TimestepControl(xf, yf, Q)
         if( time + dt > timemax ) dt = timemax - time

         Uo(:,:,:) = U(:,:,:)

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
         call Output( time, .FALSE., xv, yv, Q)

         print*, "ntime = ",ntime, "time = ",time, dt

         if( mod(ntime,10) .eq. 0 ) then
             call RealtimeAnalysis(xv,yv,Q,phys_evo)
             write(unitevo,*) time, phys_evo(1:nevo)
         endif

         if(time >= timemax) exit mloop
      enddo mloop

      close(unitevo)

!      call Output( time, .TRUE.,xv, yv, Q)

!      write(6,*) "program has been finished"
!contains
end program
!-------------------------------------------------------------------
!       Generate coordiantes
!       xf,yf --> cell boundary xf(i) <==> x_{i-1/2}
!       xv,yv --> cell center   xv(i) <==> x_{i}
!-------------------------------------------------------------------
subroutine GenerateGrid(xf, xv, yf, yv)
use params, only : nxtot, nytot, ngh, nx, ny, xmax, xmin, ymax, ymin
implicit none
real(8), intent(out) :: xf(:), xv(:)
real(8), intent(out) :: yf(:), yv(:)
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

return
end subroutine GenerateGrid
!-------------------------------------------------------------------
!       Generate initial condition of the primitive variables
!-------------------------------------------------------------------
subroutine GenerateProblem(xv, yv, Q )
use params, only : IDN, IVX, IVY, IVZ, IPR, IBX, IBY, IBZ, &
                   IMX, IMY, IMZ, IEN, IPS, NVAR, is, ie, js, je, gam
implicit none
real(8), intent(in ) :: xv(:), yv(:)
real(8), intent(out) :: Q(:,:,:)
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

return
end subroutine GenerateProblem

!-------------------------------------------------------------------
!       Boundary Condition of the primitive variables
!-------------------------------------------------------------------
subroutine BoundaryCondition(Q)
use params, only : nxtot, nytot, ngh, is, ie, js, je
implicit none
real(8), intent(inout) :: Q(:,:,:)
integer::i,j

!$omp parallel default(none) &
!$omp shared(Q) &
!$omp private(i,j) 

!$omp do collapse(2) schedule(static)
      do j=1,nytot-1
      do i=1,ngh
          Q(:,is-i,j) = Q(:,ie+1-i,j)
          Q(:,ie+i,j) = Q(:,is+i-1,j)
      enddo
      enddo
!$omp end do

!$omp do collapse(2) schedule(static)
      do j=1,ngh
      do i=1,nxtot-1
          Q(:,i,js-j)  = Q(:,i,je+1-j)
          Q(:,i,je+j)  = Q(:,i,js+j-1)
      enddo
      enddo
!$omp end do

!$omp end parallel

return
end subroutine BoundaryCondition
!-------------------------------------------------------------------
!       Primitive variables ===> Conservative variables
!       Input  : Q
!       Output : U
!-------------------------------------------------------------------
subroutine Prim2Consv(Q, U)
use params, only : IDN, IVX, IVY, IVZ, IPR, IBX, IBY, IBZ, &
                   IMX, IMY, IMZ, IEN, IPS, NVAR, is, ie, js, je, gam
implicit none
real(8), intent(in) :: Q(:,:,:)
real(8), intent(out) :: U(:,:,:)
integer::i,j

!$omp parallel default(none) &
!$omp shared(U,Q) &
!$omp private(i,j) 

!$omp do collapse(2) schedule(static)
      do j=js,je
      do i=is,ie
          U(IDN,i,j) = Q(IDN,i,j)
          U(IMX,i,j) = Q(IDN,i,j)*Q(IVX,i,j)
          U(IMY,i,j) = Q(IDN,i,j)*Q(IVY,i,j)
          U(IMZ,i,j) = Q(IDN,i,j)*Q(IVZ,i,j)
          U(IEN,i,j) = 0.5d0*Q(IDN,i,j)*( Q(IVX,i,j)**2 + Q(IVY,i,j)**2 + Q(IVZ,i,j)**2 ) &
                       + 0.5d0*( Q(IBX,i,j)**2 + Q(IBY,i,j)**2 + Q(IBZ,i,j)**2 ) &
                       + Q(IPR,i,j)/(gam - 1.0d0)
          U(IBX,i,j) = Q(IBX,i,j)
          U(IBY,i,j) = Q(IBY,i,j)
          U(IBZ,i,j) = Q(IBZ,i,j)
          U(IPS,i,j) = Q(IPS,i,j)
      enddo
      enddo
!$omp end do
!$omp end parallel
      
return
end subroutine Prim2Consv
!-------------------------------------------------------------------
!       Conservative variables ===> Primitive variables
!       Input  : U
!       Output : Q
!-------------------------------------------------------------------
subroutine Consv2Prim( U, Q )
use params, only : IDN, IVX, IVY, IVZ, IPR, IBX, IBY, IBZ, &
                   IMX, IMY, IMZ, IEN, IPS, NVAR, is, ie, js, je, gam
implicit none
real(8), intent(in) :: U(:,:,:)
real(8), intent(out) :: Q(:,:,:)
integer::i,j
real(8) :: inv_d;

!$omp parallel default(none) &
!$omp shared(U,Q) &
!$omp private(i,j,inv_d) 

!$omp do collapse(2) schedule(static)
      do j=js,je
      do i=is,ie
           Q(IDN,i,j) = U(IDN,i,j)
           inv_d = 1.0d0/U(IDN,i,j)
           Q(IVX,i,j) = U(IMX,i,j)*inv_d
           Q(IVY,i,j) = U(IMY,i,j)*inv_d
           Q(IVZ,i,j) = U(IMZ,i,j)*inv_d
           Q(IPR,i,j) = ( U(IEN,i,j) &
                        - 0.5d0*(U(IMX,i,j)**2 + U(IMY,i,j)**2 + U(IMZ,i,j)**2)*inv_d  &
                        - 0.5d0*(U(IBX,i,j)**2 + U(IBY,i,j)**2 + U(IBZ,i,j)**2) )*(gam-1.0d0)
           Q(IBX,i,j) = U(IBX,i,j)
           Q(IBY,i,j) = U(IBY,i,j)
           Q(IBZ,i,j) = U(IBZ,i,j)
           Q(IPS,i,j) = U(IPS,i,j)
      enddo
      enddo
!$omp end do
!$omp end parallel

return
end subroutine Consv2Prim
!-------------------------------------------------------------------
!       determine dt 
!-------------------------------------------------------------------
real(8) function TimestepControl(xf, yf, Q)
use params, only : IDN, IVX, IVY, IPR, IBX, IBY, IBZ, is, ie, js, je, Ccfl, gam
implicit none
real(8), intent(in) :: xf(:), yf(:), Q(:,:,:)
real(8)::dtl1
real(8)::dtl2
real(8)::dtmin,cf
integer::i,j

      dtmin=1.0d90

!$omp parallel do default(none) collapse(2) schedule(static) &
!$omp shared(xf,yf,Q) private(dtl1,dtl2,cf) reduction(min:dtmin)
      do j=js,je
      do i=is,ie
         cf = dsqrt( (gam*Q(IPR,i,j) + Q(IBX,i,j)**2 + Q(IBY,i,j)**2 + Q(IBZ,i,j)**2)/Q(IDN,i,j))
         dtl1 =(xf(i+1)-xf(i))/(abs(Q(IVX,i,j)) + cf)
         dtl2 =(yf(j+1)-yf(j))/(abs(Q(IVY,i,j)) + cf)
         dtmin = min(dtl1,dtl2,dtmin)
      enddo
      enddo
!$omp end parallel do

      TimestepControl = Ccfl* dtmin

return
end function TimestepControl

!---------------------------------------------------------------------
!     van Leer monotonicity limiter 
!---------------------------------------------------------------------
subroutine vanLeer(n,dvp,dvm,dv)
use params, only : NVAR 
implicit none
real(8),intent(in)::dvp(NVAR),dvm(NVAR)
integer,intent(in) :: n
real(8),intent(out)::dv(NVAR)
integer :: i

      do i=1,n
         if(dvp(i)*dvm(i) .gt. 0.0d0) then
            dv(i) = 2.0d0*dvp(i)*dvm(i)/(dvp(i)+dvm(i))
         else
            dv(i) = 0.0d0
         endif
      enddo

return
end subroutine vanLeer
!---------------------------------------------------------------------
!     NumericalFlux
!---------------------------------------------------------------------
!     computes the numerical flux at the cell boundary 
!
!     Input: Q: primitive variables at the cell center
!
!     Input: B: magnetic fields
!
!     Output: flux : the numerical flux estimated at the cell boundary
!---------------------------------------------------------------------
subroutine NumericalFlux( dt, xf, yf, Q, F, G )
use params, only : nxtot, nytot, NVAR, NFLX, is, ie, js, je, Ccfl, flag_flux
implicit none
real(8), intent(in) :: dt
real(8), intent(in) :: xf(:), yf(:)
real(8), intent(in) :: Q(:,:,:)
real(8), intent(out) :: F(:,:,:)
real(8), intent(out) :: G(:,:,:)

integer::i,j
real(8),dimension(NFLX,nxtot,nytot):: Ql,Qr
real(8),dimension(NFLX):: flx
real(8) :: dQm(NFLX), dQp(NFLX), dQmon(NFLX)
Real(8) :: ch

      ch = 1.0d0*Ccfl*min( xf(is+1) - xf(is), yf(js+1) - yf(js ) )/dt

!$omp parallel default(none) &
!$omp shared(Q,F,G,Ql,Qr,xf,yf,ch) &
!$omp private(i,j,flx,dQp,dQm,dQmon)

! numerical flux in the x direction
!$omp do collapse(2) schedule(static)
      do j=js,je
      do i=is-1,ie+1
         dQp(1:NVAR) = Q(1:NVAR,i+1,j) - Q(1:NVAR,i  ,j)
         dQm(1:NVAR) = Q(1:NVAR,i  ,j) - Q(1:NVAR,i-1,j)

         call vanLeer(NFLX, dQp, dQm, dQmon)

         ! Ql(i,j) --> W_(i-1/2,j)
         ! Qr(i,j) --> W_(i-1/2,j)
         Ql(1:NVAR,i+1,j) = Q(1:NVAR,i,j) + 0.5d0*dQmon(1:NVAR)
         Qr(1:NVAR,i  ,j) = Q(1:NVAR,i,j) - 0.5d0*dQmon(1:NVAR)

      enddo
      enddo
!$omp end do

! ---- x-direction: Riemann solver ----
!$omp do collapse(2) schedule(static)
  do j=js,je
    do i=is,ie+1
      if (flag_flux == 1) then
        call HLL (1, ch, Ql(:,i,j), Qr(:,i,j), flx)
      else
        call HLLD(1, ch, Ql(:,i,j), Qr(:,i,j), flx)
      end if
      F(:,i,j) = flx(:)
    end do
  end do
!$omp end do

        ! ---- y-direction: reconstruction ----
!$omp do collapse(2) schedule(static)
  do j=js-1,je+1
    do i=is,ie
      dQp(1:NVAR) = Q(1:NVAR,i,j+1) - Q(1:NVAR,i,j  )
      dQm(1:NVAR) = Q(1:NVAR,i,j  ) - Q(1:NVAR,i,j-1)
      call vanLeer(NFLX, dQp, dQm, dQmon)
      Ql(1:NVAR,i,j+1) = Q(1:NVAR,i,j) + 0.5d0*dQmon(1:NVAR)
      Qr(1:NVAR,i,j  ) = Q(1:NVAR,i,j) - 0.5d0*dQmon(1:NVAR)
    end do
  end do
!$omp end do


  ! ---- y-direction: Riemann solver ----
!$omp do collapse(2) schedule(static)
  do j=js,je+1
    do i=is,ie
      if (flag_flux == 1) then
        call HLL (2, ch, Ql(:,i,j), Qr(:,i,j), flx)
      else
        call HLLD(2, ch, Ql(:,i,j), Qr(:,i,j), flx)
      end if
      G(:,i,j) = flx(:)
    end do
  end do
!$omp end do

!$omp end parallel

return
end subroutine Numericalflux

!---------------------------------------------------------------------
!     HLL Riemann Solver
!---------------------------------------------------------------------
!     solve the HLL Riemann solver 
!
!     Input: Ql, Qr: primitive variables containing the perpendicular B fields 
!                    at the left and right states
!            1D array (IDN, IVX, IVY, IVZ, IPR, IBperp1, IBperp2)
!                                 |
!                                 |
!                           Ql    |    Qr
!                                 |
!                                -->
!                                flx
!
!     Input: b1    : magnetic field perpendicular to the initial discontinuity
!
!     Output: flx  : flux estimated at the initial discontinuity
!            index: (IDN, IVX, IVY, IVZ, IPR, IBperp1, IBperp2)
!---------------------------------------------------------------------
subroutine HLL(idir,ch,Ql,Qr,flx)
use params, only : IDN, IVX, IVY, IVZ, IPR, IBX, IBY, IBZ, &
                   IMX, IMY, IMZ, IEN, IPS, NVAR, NFLX, is, ie, js, je, gam, flag_HDC
implicit none
integer, intent(in) :: idir
real(8),intent(in)  :: ch
real(8),intent(in)  :: Ql(NVAR), Qr(NVAR)
real(8),intent(out) :: flx(NFLX)
integer :: IVpara, IVperp1, IVperp2
integer :: IBpara, IBperp1, IBperp2
real(8):: b1
real(8):: Ul(NFLX), Ur(NFLX)
real(8):: Fl(NFLX), Fr(NFLX)
real(8):: cfl,cfr
real(8):: sl, sr
real(8):: pbl, pbr, ptotl, ptotr

      if( idir == 1 ) then
           IVpara  = IVX
           IVperp1 = IVY
           IVperp2 = IVZ
           IBpara  = IBX
           IBperp1 = IBY
           IBperp2 = IBZ
      else if (idir == 2 ) then
           IVpara  = IVY
           IVperp1 = IVZ
           IVperp2 = IVX
           IBpara  = IBY
           IBperp1 = IBZ
           IBperp2 = IBX
      endif
          
          b1 = 0.5d0*( Ql(IBpara) + Qr(IBpara) )
          pbl = 0.5d0*(b1**2 + Ql(IBperp1)**2 + Ql(IBperp2)**2)
          pbr = 0.5d0*(b1**2 + Qr(IBperp1)**2 + Qr(IBperp2)**2)
          ptotl = Ql(IPR) + pbl
          ptotr = Qr(IPR) + pbr

          ! conserved variables in the left and right states
          Ul(IDN) = Ql(IDN)
          Ul(IVpara) = Ql(IDN)*Ql(IVpara)
          Ul(IVperp1) = Ql(IDN)*Ql(IVperp1)
          Ul(IVperp2) = Ql(IDN)*Ql(IVperp2)
          Ul(IEN) = 0.5d0*Ql(IDN)*( Ql(IVpara)**2 + Ql(IVperp1)**2 + Ql(IVperp2)**2) & 
                  + pbl + Ql(IPR)/(gam - 1.0d0)
          Ul(IBperp1) = Ql(IBperp1)
          Ul(IBperp2) = Ql(IBperp2)

          Ur(IDN) = Qr(IDN)
          Ur(IVpara) = Qr(IDN)*Qr(IVpara)
          Ur(IVperp1) = Qr(IDN)*Qr(IVperp1)
          Ur(IVperp2) = Qr(IDN)*Qr(IVperp2)
          Ur(IEN) = 0.5d0*Qr(IDN)*( Qr(IVpara)**2 + Qr(IVperp1)**2 + Qr(IVperp2)**2) & 
                  + pbr + Qr(IPR)/(gam - 1.0d0)
          Ur(IBperp1) = Qr(IBperp1)
          Ur(IBperp2) = Qr(IBperp2)

    !--- Step 3.  Compute L/R fluxes
          Fl(IDN) = Ul(IVpara)
          Fl(IVpara) = Ul(IVpara)*Ql(IVpara) + ptotl - b1**2
          Fl(IVperp1) = Ul(IVperp1)*Ql(IVpara) - b1*Ql(IBperp1)
          Fl(IVperp2) = Ul(IVperp2)*Ql(IVpara) - b1*Ql(IBperp2)
          Fl(IEN) = ( Ul(IEN) + ptotl - b1**2 )*Ql(IVpara) &
                  - b1*( Ql(IBperp1)*Ql(IVperp1) + Ql(IBperp2)*Ql(IVperp2) )
          Fl(IBperp1) = Ql(IBperp1)*Ql(IVpara) - b1*Ql(IVperp1)
          Fl(IBperp2) = Ql(IBperp2)*Ql(IVpara) - b1*Ql(IVperp2)

          Fr(IDN) = Ur(IVpara)
          Fr(IVpara) = Ur(IVpara)*Qr(IVpara) + ptotr - b1**2
          Fr(IVperp1) = Ur(IVperp1)*Qr(IVpara) - b1*Qr(IBperp1)
          Fr(IVperp2) = Ur(IVperp2)*Qr(IVpara) - b1*Qr(IBperp2)
          Fr(IEN) = ( Ur(IEN) + ptotr - b1**2 )*Qr(IVpara) &
                  - b1*( Qr(IBperp1)*Qr(IVperp1) + Qr(IBperp2)*Qr(IVperp2) )
          Fr(IBperp1) = Qr(IBperp1)*Qr(IVpara) - b1*Qr(IVperp1)
          Fr(IBperp2) = Qr(IBperp2)*Qr(IVpara) - b1*Qr(IVperp2)

!          cfl = dsqrt( (gam*Ql(IPR) + Ql(IBperp1)**2 + Ql(IBperp2)**2 + b1**2)/Ql(IDN))
!          cfr = dsqrt( (gam*Qr(IPR) + Qr(IBperp1)**2 + Qr(IBperp2)**2 + b1**2)/Qr(IDN))
          cfl = dsqrt( 0.5d0*( 2.0d0*pbl + gam*Ql(IPR) &
                     + dsqrt( (2.0d0*pbl - gam*Ql(IPR))**2 &
                     + 4.0d0*gam*Ql(IPR)*( Ql(IBperp1)**2 + Ql(IBperp2)**2 ) ) )/Ql(IDN) )
          cfr = dsqrt( 0.5d0*( 2.0d0*pbr + gam*Qr(IPR) &
                    + dsqrt( (2.0d0*pbr - gam*Qr(IPR))**2 &
                     + 4.0d0*gam*Qr(IPR)*( Qr(IBperp1)**2 + Qr(IBperp2)**2 ) ) )/Qr(IDN) )

          sl = min(Ql(IVpara) - cfl,Qr(IVpara) - cfr)
          sr = max(Ql(IVpara) + cfl,Qr(IVpara) + cfr)

          if( sl > 0.0d0 ) then
               flx(:) = Fl(:)
          else if (sr <= 0.0d0 ) then
               flx(:) = Fr(:)
          else 
               flx(:)  = (sr*Fl(:) - sl*Fr(:) + sl*sr*( Ur(:) - Ul(:) ))/(sr - sl)
          endif

          flx(IBpara) = flag_HDC*0.5d0*(Ql(IPS) + Qr(IPS) - ch*(Qr(IBpara) - Ql(IBpara)) )
          flx(IPS) = flag_HDC*0.5d0*(ch*ch*(Ql(IBpara) + Qr(IBpara)) - ch*(Qr(IPS) - Ql(IPS)) )
    

return
end subroutine HLL
!---------------------------------------------------------------------
!     HLLD Riemann Solver
!---------------------------------------------------------------------
!     solve the HLL Riemann solver 
!
!     Input: Ql, Qr: primitive variables containing the perpendicular B fields 
!                    at the left and right states
!            1D array (IDN, IVX, IVY, IVZ, IPR, IBperp1, IBperp2)
!                                 |
!                                 |
!                           Ql    |    Qr
!                                 |
!                                -->
!                                flx
!
!     Input: b1    : magnetic field perpendicular to the initial discontinuity
!
!     Output: flx  : flux estimated at the initial discontinuity
!            index: (IDN, IVX, IVY, IVZ, IPR, IBperp1, IBperp2)
!---------------------------------------------------------------------
subroutine HLLD(idir,ch,Ql,Qr,flx)
use params, only : IDN, IVX, IVY, IVZ, IPR, IBX, IBY, IBZ, &
                   IMX, IMY, IMZ, IEN, IPS, NVAR, NFLX, is, ie, js, je, gam, flag_HDC
implicit none
integer, intent(in) :: idir
real(8),intent(in)  :: ch
real(8),intent(in)  :: Ql(NVAR), Qr(NVAR)
real(8),intent(out) :: flx(NFLX)
integer :: IVpara, IVperp1, IVperp2
integer :: IBpara, IBperp1, IBperp2
real(8):: b1
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
      Fl = 0.0d0
      Fr = 0.0d0
      Ul = 0.0d0
      Ur = 0.0d0
      Ulst = 0.0d0
      Urst = 0.0d0
      Uldst = 0.0d0
      Urdst = 0.0d0

      if( idir == 1 ) then
           IVpara  = IVX
           IVperp1 = IVY
           IVperp2 = IVZ
           IBpara  = IBX
           IBperp1 = IBY
           IBperp2 = IBZ
      else if (idir == 2 ) then
           IVpara  = IVY
           IVperp1 = IVZ
           IVperp2 = IVX
           IBpara  = IBY
           IBperp1 = IBZ
           IBperp2 = IBX
      endif
          b1 = 0.5d0*( Ql(IBpara) + Qr(IBpara) )
          
          pbl = 0.5d0*(b1**2 + Ql(IBperp1)**2 + Ql(IBperp2)**2)
          pbr = 0.5d0*(b1**2 + Qr(IBperp1)**2 + Qr(IBperp2)**2)
          ptotl = Ql(IPR) + pbl
          ptotr = Qr(IPR) + pbr

          cfl = dsqrt( 0.5d0*( 2.0d0*pbl + gam*Ql(IPR) &
                     + dsqrt( (2.0d0*pbl - gam*Ql(IPR))**2 &
                     + 4.0d0*gam*Ql(IPR)*( Ql(IBperp1)**2 + Ql(IBperp2)**2 ) ) )/Ql(IDN) )
          cfr = dsqrt( 0.5d0*( 2.0d0*pbr + gam*Qr(IPR) &
                    + dsqrt( (2.0d0*pbr - gam*Qr(IPR))**2 &
                     + 4.0d0*gam*Qr(IPR)*( Qr(IBperp1)**2 + Qr(IBperp2)**2 ) ) )/Qr(IDN) )
!          cfl = dsqrt( (gam*Ql(IPR) + Ql(IBperp1)**2 + Ql(IBperp2)**2 + b1**2)/Ql(IDN))
!          cfr = dsqrt( (gam*Qr(IPR) + Qr(IBperp1)**2 + Qr(IBperp2)**2 + b1**2)/Qr(IDN))
!
          S0 = min( Ql(IVpara) - cfl, Qr(IVpara) - cfr)
          S4 = max( Ql(IVpara) + cfl, Qr(IVpara) + cfr)

          ! conserved variables in the left and right states
          Ul(IDN) = Ql(IDN)
          Ul(IVpara) = Ql(IDN)*Ql(IVpara)
          Ul(IVperp1) = Ql(IDN)*Ql(IVperp1)
          Ul(IVperp2) = Ql(IDN)*Ql(IVperp2)
          Ul(IEN) = 0.5d0*Ql(IDN)*( Ql(IVpara)**2 + Ql(IVperp1)**2 + Ql(IVperp2)**2) & 
                  + pbl + Ql(IPR)/(gam - 1.0d0)
          Ul(IBperp1) = Ql(IBperp1)
          Ul(IBperp2) = Ql(IBperp2)

          Ur(IDN) = Qr(IDN)
          Ur(IVpara) = Qr(IDN)*Qr(IVpara)
          Ur(IVperp1) = Qr(IDN)*Qr(IVperp1)
          Ur(IVperp2) = Qr(IDN)*Qr(IVperp2)
          Ur(IEN) = 0.5d0*Qr(IDN)*( Qr(IVpara)**2 + Qr(IVperp1)**2 + Qr(IVperp2)**2) & 
                  + pbr + Qr(IPR)/(gam - 1.0d0)
          Ur(IBperp1) = Qr(IBperp1)
          Ur(IBperp2) = Qr(IBperp2)

    !--- Step 3.  Compute L/R fluxes
          Fl(IDN) = Ul(IVpara)
          Fl(IVpara) = Ul(IVpara)*Ql(IVpara) + ptotl - b1**2
          Fl(IVperp1) = Ul(IVperp1)*Ql(IVpara) - b1*Ql(IBperp1)
          Fl(IVperp2) = Ul(IVperp2)*Ql(IVpara) - b1*Ql(IBperp2)
          Fl(IEN) = ( Ul(IEN) + ptotl - b1**2 )*Ql(IVpara) &
                  - b1*( Ql(IBperp1)*Ql(IVperp1) + Ql(IBperp2)*Ql(IVperp2) )
          Fl(IBperp1) = Ql(IBperp1)*Ql(IVpara) - b1*Ql(IVperp1)
          Fl(IBperp2) = Ql(IBperp2)*Ql(IVpara) - b1*Ql(IVperp2)

          Fr(IDN) = Ur(IVpara)
          Fr(IVpara) = Ur(IVpara)*Qr(IVpara) + ptotr - b1**2
          Fr(IVperp1) = Ur(IVperp1)*Qr(IVpara) - b1*Qr(IBperp1)
          Fr(IVperp2) = Ur(IVperp2)*Qr(IVpara) - b1*Qr(IBperp2)
          Fr(IEN) = ( Ur(IEN) + ptotr - b1**2 )*Qr(IVpara) &
                  - b1*( Qr(IBperp1)*Qr(IVperp1) + Qr(IBperp2)*Qr(IVperp2) )
          Fr(IBperp1) = Qr(IBperp1)*Qr(IVpara) - b1*Qr(IVperp1)
          Fr(IBperp2) = Qr(IBperp2)*Qr(IVpara) - b1*Qr(IVperp2)

    !--- Step 4.  Compute middle and Alfven wave speeds
          Cl = S0 - Ql(IVpara)
          Cr = S4 - Qr(IVpara)

          S2 = ( Cr*Ur(IVpara) - Cl*Ul(IVpara) + (ptotl - ptotr) ) &
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

          S1 = S2 - dabs(b1)/sqrtdl
          S3 = S2 + dabs(b1)/sqrtdr

    !--- Step 5.  Compute intermediate states
         ptot_stl = ptotl + Ul(IDN)*Cl*(S2 - Ql(IVpara))
         ptot_str = ptotr + Ur(IDN)*Cr*(S2 - Qr(IVpara))

         ptot_st = 0.5d0*(ptot_stl + ptot_str)

         Ulst(IVpara) = Ulst(IDN)*S2
         if( dabs( Ul(IDN)*Cl*Cml-b1**2) < 1.0d-8*ptot_st ) then
             Ulst(IVperp1) = Ulst(IDN)*Ql(IVperp1)
             Ulst(IVperp2) = Ulst(IDN)*Ql(IVperp2)

             Ulst(IBperp1) = Ul(IBperp1)
             Ulst(IBperp2) = Ul(IBperp2)
         else 
             tmp = b1*( Cl - Cml )/(Ul(IDN)*Cl*Cml - b1**2)
             Ulst(IVperp1) = Ulst(IDN)*( Ql(IVperp1) - Ul(IBperp1)*tmp )
             Ulst(IVperp2) = Ulst(IDN)*( Ql(IVperp2) - Ul(IBperp2)*tmp )

             tmp = (Ul(IDN)*Cl**2 - b1**2)/( Ul(IDN)*Cl*Cml - b1**2)
             Ulst(IBperp1) = Ul(IBperp1)*tmp
             Ulst(IBperp2) = Ul(IBperp2)*tmp
         endif

         v_dot_B_stl = ( Ulst(IVpara)*b1 + Ulst(IVperp1)*Ulst(IBperp1) + Ulst(IVperp2)*Ulst(IBperp2) )*Ulst_d_inv
         Ulst(IEN) = ( Cl*Ul(IEN) - ptotl*Ql(IVpara) + ptot_st*S2 &
                     + b1*( Ql(IVpara)*b1 + Ql(IVperp1)*Ul(IBperp1) + Ql(IVperp2)*Ul(IBperp2) - v_dot_B_stl) )*Cml_inv

         Urst(IVpara) = Urst(IDN)*S2
         if( dabs( Ur(IDN)*Cr*Cmr-b1**2) < 1.0d-8*ptot_st ) then
             Urst(IVperp1) = Urst(IDN)*Qr(IVperp1)
             Urst(IVperp2) = Urst(IDN)*Qr(IVperp2)

             Urst(IBperp1) = Ur(IBperp1)
             Urst(IBperp2) = Ur(IBperp2)
         else 
             tmp = b1*( Cr - Cmr )/(Ur(IDN)*Cr*Cmr - b1**2)
             Urst(IVperp1) = Urst(IDN)*( Qr(IVperp1) - Ur(IBperp1)*tmp )
             Urst(IVperp2) = Urst(IDN)*( Qr(IVperp2) - Ur(IBperp2)*tmp )

             tmp = (Ur(IDN)*Cr**2 - b1**2)/( Ur(IDN)*Cr*Cmr - b1**2)
             Urst(IBperp1) = Ur(IBperp1)*tmp
             Urst(IBperp2) = Ur(IBperp2)*tmp
         endif

         v_dot_B_str = ( Urst(IVpara)*b1 + Urst(IVperp1)*Urst(IBperp1) + Urst(IVperp2)*Urst(IBperp2) )*Urst_d_inv
         Urst(IEN) = ( Cr*Ur(IEN) - ptotr*Qr(IVpara) + ptot_st*S2 &
                     + b1*( Qr(IVpara)*b1 + Qr(IVperp1)*Ur(IBperp1) + Qr(IVperp2)*Ur(IBperp2) - v_dot_B_str) )*Cmr_inv
       
         if( 0.5d0*b1**2 < 1.0d-8*ptot_st )then
             Uldst(:) = Ulst(:)
             Urdst(:) = Urst(:)
         else 
             sum_sqrtd_inv = 1.0d0/(sqrtdl + sqrtdr) 
             if (b1 > 0.0d0 ) then 
                 bxsgn = 1.0d0
             else 
                 bxsgn = -1.0d0
             endif

             Uldst(IDN) = Ulst(IDN)
             Urdst(IDN) = Urst(IDN)

             Uldst(IVpara) = Ulst(IVpara)
             Urdst(IVpara) = Urst(IVpara)

             tmp = sum_sqrtd_inv*(  sqrtdl*(Ulst(IVperp1)*Ulst_d_inv) + sqrtdr*(Urst(IVperp1)*Urst_d_inv) &
                                  + bxsgn*(Urst(IBperp1) - Ulst(IBperp1)) )
             Uldst(IVperp1) = Uldst(IDN)*tmp
             Urdst(IVperp1) = Urdst(IDN)*tmp
!

             tmp = sum_sqrtd_inv*(  sqrtdl*(Ulst(IVperp2)*Ulst_d_inv) + sqrtdr*(Urst(IVperp2)*Urst_d_inv) &
                                  + bxsgn*(Urst(IBperp2) - Ulst(IBperp2)) )
             Uldst(IVperp2) = Uldst(IDN)*tmp
             Urdst(IVperp2) = Urdst(IDN)*tmp

             tmp = sum_sqrtd_inv*(  sqrtdl*Urst(IBperp1) + sqrtdr*Ulst(IBperp1) &
                       + bxsgn*sqrtdl*sqrtdr*( (Urst(IVperp1)*Urst_d_inv) - (Ulst(IVperp1)*Ulst_d_inv) ) )
             Uldst(IBperp1) = tmp
             Urdst(IBperp1) = tmp

             tmp = sum_sqrtd_inv*(  sqrtdl*Urst(IBperp2) + sqrtdr*Ulst(IBperp2) &
                       + bxsgn*sqrtdl*sqrtdr*( (Urst(IVperp2)*Urst_d_inv) - (Ulst(IVperp2)*Ulst_d_inv) ) )
             Uldst(IBperp2) = tmp
             Urdst(IBperp2) = tmp
!
             tmp = S2*b1 + (Uldst(IVperp1)*Uldst(IBperp1) + Uldst(IVperp2)*Uldst(IBperp2))/Uldst(IDN)
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
           flx(IBpara) = 0.0d0

           flx(IBpara) = flag_HDC*0.5d0*(Ql(IPS) + Qr(IPS) - ch*(Qr(IBpara) - Ql(IBpara)) )
           flx(IPS) = flag_HDC*0.5d0*(ch*ch*(Ql(IBpara) + Qr(IBpara)) - ch*(Qr(IPS) - Ql(IPS)) )
    

return
end subroutine HLLD
!-------------------------------------------------------------------
!       Update consevative variables U using numerical flux F
!-------------------------------------------------------------------
subroutine UpdateConsv( dt1, xf, yf, F, G, Uo, U)
use params, only : IDN, IVX, IVY, IVZ, IPR, IBX, IBY, IBZ, &
                   IMX, IMY, IMZ, IEN, IPS, NVAR, is, ie, js, je, gam, alpha, Ccfl
implicit none
real(8), intent(in) :: dt1 
real(8), intent(in)  :: xf(:), yf(:)
real(8), intent(in)  :: F(:,:,:), G(:,:,:)
real(8), intent(in)  :: Uo(:,:,:)
real(8), intent(inout) :: U(:,:,:)
integer::i,j

!$omp parallel default(none) &
!$omp shared(U,Uo,F,G,xf,yf) &
!$omp shared(dt1,dt0) &
!$omp private(i,j)

!$omp do collapse(2) schedule(static)
      do j=js,je
      do i=is,ie
         U(:,i,j) = Uo(:,i,j) + dt1*(- F(:,i+1,j) + F(:,i,j))/(xf(i+1)-xf(i)) &
                              + dt1*(- G(:,i,j+1) + G(:,i,j))/(yf(j+1)-yf(j))
      enddo
      enddo
!$omp end do

!$omp end parallel

return
end subroutine UpdateConsv
!-------------------------------------------------------------------
!       Update consevative variables U using numerical flux F
!-------------------------------------------------------------------
subroutine SrcTerms( dt1, dt0, Q, U )
use params, only : IDN, IVX, IVY, IVZ, IPR, IBX, IBY, IBZ, &
                   IMX, IMY, IMZ, IEN, IPS, NVAR, is, ie, js, je, gam, alpha, Ccfl
implicit none
real(8), intent(in) :: dt1, dt0
real(8), intent(in)  :: Q(:,:,:)
real(8), intent(inout) :: U(:,:,:)
integer :: i,j

      ! Source term
!$omp parallel 
      !$omp do private( i, j )
      do j=js,je
      do i=is,ie
         U(IPS,i,j) = U(IPS,i,j)*dexp(-alpha*Ccfl*dt1/dt0)
      enddo
      enddo
      !$omp end do
!$omp end parallel

end subroutine SrcTerms
!-------------------------------------------------------------------
!       Output snapshot files 
!       Input  : flag, dirname, xf, xv, Q
!
!       flag = .true.  --> output snapshot when calling this subroutine
!       flag = .false. --> output snapshot every dtsnap
!-------------------------------------------------------------------
subroutine Output( time, flag, xv, yv, Q )
use params, only : IDN, IVX, IVY, IVZ, IPR, IBX, IBY, IBZ, IPS, NVAR, & 
                   nx, ny, is, ie, js, je, gam, &
                   flag_binary, dirname, dtsnap, unitsnap
implicit none
real(8), intent(in) :: time! false --> output per dtsnap, true --> force to output
logical, intent(in) :: flag
real(8), intent(in) :: xv(:), yv(:), Q(:,:,:)

integer::i,j
character(100)::filename
real(8), save :: tsnap = - dtsnap
integer, save :: nsnap = 0


    if( .not.flag) then
        if( time + 1.0d-14.lt. tsnap+dtsnap) return
    endif

    write(filename,'(i5.5)') nsnap
    if( flag_binary ) then
        filename = trim(dirname)//"/snap"//trim(filename)//".bin"
        open(unitsnap,file=filename,form='unformatted',access="stream",action="write")
        write(unitsnap) time
        write(unitsnap) nx
        write(unitsnap) ny
        write(unitsnap) 5
        write(unitsnap) NVAR - 5
        write(unitsnap) xv(is:ie)
        write(unitsnap) yv(js:je)
        write(unitsnap) real(Q(1:5,is:ie,js:je)) ! single precision
        write(unitsnap) real(Q(6:NVAR,is:ie,js:je)) ! single precision
        close(unitsnap)
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
      endif

    write(6,*) "output binary file:  ",filename,time

    nsnap=nsnap+1
    tsnap=tsnap + dtsnap

return
end subroutine Output
!
!!-------------------------------------------------------------------
!!       create directory
!!       Input  : the directory to be created
!!-------------------------------------------------------------------
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
!-------------------------------------------------------------------
!       Realtime Analysis
!       Input  : xf, xv
!       Output : phys_evo(nevo)
!-------------------------------------------------------------------
subroutine RealtimeAnalysis(xv,yv,Q,phys_evo)
use params, only : IDN, IVX, IVY, IVZ, IPR, IBX, IBY, IBZ, &
                   IMX, IMY, IMZ, IEN, IPS, NVAR, is, ie, js, je, gam, nevo
implicit none
real(8), intent(in)  :: xv(:), yv(:), Q(:,:,:)
real(8), intent(out) :: phys_evo(:)
integer::i,j
real(8) :: tmp

      
      tmp = 0.0d0
      do j=js,je
      do i=is,ie
          tmp = tmp + Q(IDN,i,j)*Q(IPR,i,j)/(xv(i)+yv(j))
      enddo
      enddo
      phys_evo(1:nevo) = 0.0d0
      
return
end subroutine
