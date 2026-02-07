module params
    real(8),parameter:: timemax=0.15d0 ! simulation end time
    
    ! option
    integer, parameter :: flag_flux = 2 ! 1 (HLL), 2 (HLLD)
    
    ! coordinate 
    integer,parameter::nx=128! the number of grids in the simulation box
    integer,parameter::ny=128 ! the number of grids in the simulation box
    integer,parameter::nz=1          ! the number of grids in the simulation box
    integer,parameter::ngh=2         ! the number of ghost cells
    integer,parameter::nxtot=nx+2*ngh+1 ! the total number of grids including ghost cells
    integer,parameter::nytot=ny+2*ngh+1 ! the total number of grids including ghost cells
    integer,parameter::nztot=2 ! the total number of grids including ghost cells
    integer,parameter::is=ngh+1         ! the index of the leftmost grid
    integer,parameter::js=ngh+1         ! the index of the leftmost grid
    integer,parameter::ks=1         ! the index of the leftmost grid
    integer,parameter::ie=nx+ngh     ! the index of the rightmost grid
    integer,parameter::je=ny+ngh     ! the index of the rightmost grid
    integer,parameter::ke=1 ! the index of the rightmost grid
    real(8),parameter::xmin=-0.5d0,xmax=0.5d0
    real(8),parameter::ymin=-0.5d0,ymax=0.5d0
    real(8),parameter::zmin=0.0d0,zmax=1.0d0
    
    real(8),parameter::Ccfl=0.4d0
    
    
    ! indices of the conservative variables
    integer, parameter :: IDN = 1
    integer, parameter :: IMX = 2
    integer, parameter :: IMY = 3
    integer, parameter :: IMZ = 4
    integer, parameter :: IPR = 5
    integer, parameter :: IBX = 6
    integer, parameter :: IBY = 7
    integer, parameter :: IBZ = 8
    integer, parameter :: NVAR = 5
    integer, parameter :: NFLX = 8
    
    ! indices of the primitive variables
    integer, parameter :: IVX = 2
    integer, parameter :: IVY = 3
    integer, parameter :: IVZ = 4
    integer, parameter :: IEN = 5
    
    real(8),parameter::gam=5.0d0/3.0d0 !! adiabatic index

    ! output 
    character(20),parameter::dirname="ct" ! directory name
    logical, parameter :: flag_binary = .true.
    
    ! snapshot
    integer, parameter :: unitsnap = 17
    real(8), parameter:: dtsnap=0.5d-2
    
    ! realtime analysis 
    integer, parameter :: nevo = 3
    integer, parameter :: unitevo =11
    integer, parameter :: unitbin =13
end module

program main
!$ use omp_lib
use params, only : nxtot, nytot, nztot, NVAR, dirname, unitevo, timemax, nevo
implicit none
include "interfaces_ct.inc"

! time evolution
integer :: ntime = 0    ! counter of the timestep
real(8) :: time = 0.0d0  ! time 
real(8) :: dt   = 0.0d0  ! time width


! definition of arrays 
real(8),dimension(nxtot)::xf,xv
real(8),dimension(nytot)::yf,yv
real(8),dimension(nztot)::zf,zv
real(8),dimension(NVAR,nxtot,nytot,nztot) :: Uo
real(8),dimension(NVAR,nxtot,nytot,nztot) :: U
real(8),dimension(NVAR,nxtot,nytot,nztot) :: Q
real(8),dimension(3,nxtot,nytot,nztot) :: Bso 
real(8),dimension(3,nxtot,nytot,nztot) :: Bs 
real(8),dimension(3,nxtot,nytot,nztot) :: Bc
real(8),dimension(NVAR,nxtot,nytot,nztot) :: F
real(8),dimension(NVAR,nxtot,nytot,nztot) :: G
real(8),dimension(NVAR,nxtot,nytot,nztot) :: H
real(8),dimension(3,nxtot,nytot,nztot) :: E 

real(8) :: phys_evo(nevo)

!logical :: flag_binary = .false.

      ! make the directory for output
      call makedirs(trim(dirname))

    write(6,*) "setup grids and initial condition"
    call GenerateGrid(xf, xv, yf, yv, zf, zv)
    call GenerateProblem(xv, yv, zv, Q, Bs, Bc)
    call Prim2Consv(Q, Bc, U)
    call BoundaryCondition( Q, Bs, Bc )
    call Output( time, .TRUE., dirname, xf, xv, yf, yv, Q, Bc )


  write(6,*) "Start the simulation"
  open(unitevo,file=trim(dirname)//'/'//'ana.dat', action="write")
! main loop
  do 
        dt = TimestepControl(xf, yf, zf, Q, Bc)
        if( time + dt > timemax ) dt = timemax - time

        Uo(:,:,:,:) = U(:,:,:,:)
        Bso(:,:,:,:) = Bs(:,:,:,:)


        call NumericalFlux( 0.5*dt, xf, yf, zf, Q, Bc, Bs, F, G, H, E )
        call UpdateConsv( 0.5d0*dt,  xf, xv, yf, yv, zf, zv, F, G, H, E, Q, U, Bs, U, Bs )
        call Consv2Prim( U, Bs, Q, Bc )
        call BoundaryCondition( Q, Bs, Bc )

        call NumericalFlux( dt, xf, yf, zf, Q, Bc, Bs, F, G, H, E )
        call UpdateConsv( dt, xf, xv, yf, yv, zf, zv, F, G, H, E, Q, Uo, Bso, U, Bs )
        call Consv2Prim( U, Bs, Q, Bc )
        call BoundaryCondition( Q, Bs, Bc )

         time=time+dt
         ntime = ntime + 1

         if( mod(ntime,10) .eq. 0 ) then
             call RealtimeAnalysis(xv,yv,Q,Bc,Bs,phys_evo)
             write(unitevo,*) time, phys_evo(1:nevo)
         endif
         call Output( time, .FALSE., dirname, xf, xv, yf, yv, Q, Bc)
!         call Output( .true., xf, xv, yf, yv, Q, Bc)

         print*, "time = ",time, "dt = ",dt

         if(time >= timemax) exit 
      enddo 
      close(unitevo)
      call Output( time, .TRUE., dirname, xf, xv, yf, yv, Q, Bc)


!      write(6,*) "program has been finished"

end program main
!-------------------------------------------------------------------
!       Generate coordiantes
!       xf,yf,zf --> cell boundary xf(i) <==> x_{i-1/2}
!       xv,yv,zv --> cell center   xv(i) <==> x_{i}
!-------------------------------------------------------------------
subroutine GenerateGrid(xf, xv, yf, yv, zf, zv)
use params, only : nxtot, nytot, nztot, nx, ny, ngh, xmax, xmin, ymax, ymin 
implicit none
real(8), intent(out) :: xf(:), xv(:)
real(8), intent(out) :: yf(:), yv(:)
real(8), intent(out) :: zf(:), zv(:)
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
         yf(j) = dx*(j-(ngh+1))+ymin
    enddo
    do j=1,nytot-1
         yv(j) = 0.5d0*(yf(j+1)+yf(j))
    enddo

return
end subroutine GenerateGrid
!-------------------------------------------------------------------
!       Generate initial condition of the primitive variables
!-------------------------------------------------------------------
subroutine GenerateProblem(xv, yv, zv, Q, Bs, Bc )
use params, only : is, ie, js, je, ks, ke, &
                   IDN, IVX, IVY, IVZ, IPR
implicit none
integer::i, j, k
real(8), intent(in ) :: xv(:), yv(:), zv(:)
real(8), intent(out) :: Q(:,:,:,:)
real(8), intent(out) :: Bs(:,:,:,:)
real(8), intent(out) :: Bc(:,:,:,:)
real(8) :: pi, B0

    B0 = 10d0
    pi = dacos(-1.0d0)

    do k=ks,ke
    do j=js,je
    do i=is,ie
         Q(IDN,i,j,k) = 1.0d0
         if( xv(i)**2 + yv(j)**2 < 0.1d0**2 ) then
             Q(IPR,i,j,k) = 10.0d0
         else 
             Q(IPR,i,j,k) = 0.1d0
         endif
         Q(IVX,i,j,k) = 0.0d0
         Q(IVY,i,j,k) = 0.0d0
         Q(IVZ,i,j,k) = 0.0d0
    enddo
    enddo
    enddo

    do k=ks,ke
    do j=js,je
    do i=is,ie+1
        Bs(1,i,j,k) = B0/sqrt(2.0d0)
    enddo
    enddo
    enddo

    do k=ks,ke
    do j=js,je+1
    do i=is,ie
        Bs(2,i,j,k) = B0/sqrt(2.0d0)
    enddo
    enddo
    enddo

    do k=ks,ke+1
    do j=js,je
    do i=is,ie
        Bs(3,i,j,k) = 0.0d0
    enddo
    enddo
    enddo

    call CellCenterMagneticField(is, ie, js, je, ks, ke, Bs, Bc)


return
end subroutine GenerateProblem
!-------------------------------------------------------------------
!       Boundary Condition of the primitive variables
!-------------------------------------------------------------------
subroutine BoundaryCondition( Q, Bs, Bc )
use params, only : is, ie, js, je, ks, ke, ngh, &
                   IDN, IVX, IVY, IVZ, IPR
implicit none
real(8), intent(inout) :: Q(:,:,:,:)
real(8), intent(inout) :: Bs(:,:,:,:)
real(8), intent(inout) :: Bc(:,:,:,:)
integer::i,j,k

    ! x inner boundary
    do k=ks,ke
    do j=js-ngh,je+ngh
    do i=1,ngh
        Q(:,is-i,j,k)  = Q(:,ie+1-i,j,k)
    enddo
    enddo
    enddo

    do k=ks,ke
    do j=js-ngh,je+ngh
    do i=1,ngh
          Bs(1,is-i,j,k) = Bs(1,ie+1-i,j,k)
    enddo
    enddo
    enddo

    do k=ks,ke
    do j=js-ngh,je+ngh+1
    do i=1,ngh
          Bs(2,is-i,j,k) = Bs(2,ie+1-i,j,k)
    enddo
    enddo
    enddo

    do k=ks,ke+1
    do j=js-ngh,je+ngh
    do i=1,ngh
          Bs(3,is-i,j,k) = Bs(3,ie+1-i,j,k)
    enddo
    enddo
    enddo

    ! x outer boundary
    do k=ks,ke
    do j=js-ngh,je+ngh
    do i=1,ngh
        Q(:,ie+i,j,k) = Q(:,is+i-1,j,k)
    enddo
    enddo
    enddo

    do k=ks,ke
    do j=js-ngh,je+ngh
    do i=1,ngh
        Bs(1,ie+i+1,j,k) = Bs(1,is+i,j,k)
    enddo
    enddo
    enddo

    do k=ks,ke
    do j=js-ngh,je+ngh+1
    do i=1,ngh
        Bs(2,ie+i,j,k) = Bs(2,is+i-1,j,k)
    enddo
    enddo
    enddo

    do k=ks,ke+1
    do j=js-ngh,je+ngh
    do i=1,ngh
        Bs(3,ie+i,j,k) = Bs(3,is+i-1,j,k)
    enddo
    enddo
    enddo

    ! y inner boundary
    do k=ks,ke
    do j=1,ngh
    do i=is-ngh,ie+ngh
        Q(:,i,js-j,k)  = Q(:,i,je+1-j,k)
    enddo
    enddo
    enddo

    do k=ks,ke
    do j=1,ngh
    do i=is-ngh,ie+ngh+1
        Bs(1,i,js-j,k) = Bs(1,i,je+1-j,k)
    enddo
    enddo
    enddo

    do k=ks,ke
    do j=1,ngh
    do i=is-ngh,ie+ngh
        Bs(2,i,js-j,k) = Bs(2,i,je+1-j,k)
    enddo
    enddo
    enddo


    do k=ks,ke+1
    do j=1,ngh
    do i=is-ngh,ie+ngh
        Bs(3,i,js-j,k) = Bs(3,i,je+1-j,k)
    enddo
    enddo
    enddo

    ! y outer boundary
    do k=ks,ke
    do j=1,ngh
    do i=is-ngh,ie+ngh
       Q(:,i,je+j,k)  = Q(:,i,js+j-1,k)
    enddo
    enddo
    enddo

    do k=ks,ke
    do j=1,ngh
    do i=is-ngh,ie+ngh+1
          Bs(1,i,je+j,k) = Bs(1,i,js+j-1,k)
    enddo
    enddo
    enddo

    do k=ks,ke
    do j=1,ngh
    do i=is-ngh,ie+ngh
        Bs(2,i,je+j+1,k) = Bs(2,i,js+j,k)
    enddo
    enddo
    enddo

    do k=ks,ke+1
    do j=1,ngh
    do i=is-ngh,ie+ngh
        Bs(3,i,je+j,k) = Bs(3,i,js+j-1,k)
    enddo
    enddo
    enddo

    ! boundary condition for the cell centered B field
    call CellCenterMagneticField(is-ngh, is-1, js-ngh, je+ngh, ks, ke, Bs, Bc)
    call CellCenterMagneticField(ie+1,   ie+ngh, js-ngh, je+ngh, ks, ke, Bs, Bc)
    call CellCenterMagneticField(is-ngh, ie+ngh, js-ngh, js-1, ks, ke, Bs, Bc)
    call CellCenterMagneticField(is-ngh, ie+ngh, je+1, je+ngh, ks, ke, Bs, Bc)

return
end subroutine BoundaryCondition
!-------------------------------------------------------------------
!       Primitive variables ===> Conservative variables
!       Input  : Q
!       Output : U
!-------------------------------------------------------------------
subroutine Prim2Consv(Q, Bc, U)
use params, only : is, ie, js, je, ks, ke, &
                   IDN, IVX, IVY, IVZ, IPR, &
                   IMX, IMY, IMZ, IEN, gam
implicit none
real(8), intent(in) :: Q(:,:,:,:)
real(8), intent(in) :: Bc(:,:,:,:)
real(8), intent(out) :: U(:,:,:,:)
integer::i,j,k

        do k=ks,ke
       !$omp parallel do private(i)
        do j=js,je
        do i=is,ie
            U(IDN,i,j,k) = Q(IDN,i,j,k)
            U(IMX,i,j,k) = Q(IDN,i,j,k)*Q(IVX,i,j,k)
            U(IMY,i,j,k) = Q(IDN,i,j,k)*Q(IVY,i,j,k)
            U(IMZ,i,j,k) = Q(IDN,i,j,k)*Q(IVZ,i,j,k)
            U(IEN,i,j,k) = 0.5d0*Q(IDN,i,j,k)*( Q(IVX,i,j,k)**2 + Q(IVY,i,j,k)**2 + Q(IVZ,i,j,k)**2 ) &
                       + 0.5d0*( Bc(1,i,j,k)**2 + Bc(2,i,j,k)**2 + Bc(3,i,j,k)**2 ) &
                       + Q(IPR,i,j,k)/(gam - 1.0d0)
        enddo
        enddo
       !$omp end parallel do
        enddo
      
return
end subroutine Prim2Consv
!-------------------------------------------------------------------
!       Conservative variables ===> Primitive variables
!       Input  : U
!       Output : Q
!-------------------------------------------------------------------
subroutine Consv2Prim( U, Bs, Q, Bc )
use params, only : is, ie, js, je, ks, ke, &
                   IDN, IVX, IVY, IVZ, IPR, &
                   IMX, IMY, IMZ, IEN, gam
implicit none
real(8), intent(in) :: U(:,:,:,:),Bs(:,:,:,:)
real(8), intent(out) :: Q(:,:,:,:),Bc(:,:,:,:)
integer::i,j,k
real(8) :: inv_d;

        call CellCenterMagneticField(is, ie, js, je, ks, ke, Bs, Bc)

        do k=ks,ke
       !$omp parallel do private(i, inv_d)
        do j=js,je
        do i=is,ie
            Q(IDN,i,j,k) = U(IDN,i,j,k)
            inv_d = 1.0d0/U(IDN,i,j,k)
            Q(IVX,i,j,k) = U(IMX,i,j,k)*inv_d
            Q(IVY,i,j,k) = U(IMY,i,j,k)*inv_d
            Q(IVZ,i,j,k) = U(IMZ,i,j,k)*inv_d
            Q(IPR,i,j,k) = ( U(IEN,i,j,k) &
                        - 0.5d0*(U(IMX,i,j,k)**2 + U(IMY,i,j,k)**2 + U(IMZ,i,j,k)**2)*inv_d  &
                        - 0.5d0*(Bc(1,i,j,k)**2 + Bc(2,i,j,k)**2 + Bc(3,i,j,k)**2) )*(gam-1.0d0)
        enddo
        enddo
        !$omp end parallel do
        enddo

return
end subroutine Consv2Prim
!-------------------------------------------------------------------
!       Cell Surface B field ===> Cell Center B field
!-------------------------------------------------------------------
subroutine CellCenterMagneticField(ibeg, ifin, jbeg, jfin, kbeg, kfin, Bs, Bc )
use params, only : nxtot, nytot, nztot
implicit none
integer, intent(in) :: ibeg, ifin, jbeg, jfin, kbeg, kfin
real(8), intent(in) :: Bs(3,nxtot,nytot,nztot)
real(8), intent(out) :: Bc(3,nxtot,nytot,nztot)
integer::i,j,k
real(8) :: inv_d;

        do k=kbeg,kfin
       !$omp parallel do private(i, inv_d)
        do j=jbeg,jfin
        do i=ibeg,ifin
            Bc(1,i,j,k) = 0.5d0*( Bs(1,i+1,j,k) + Bs(1,i,j,k) )
            Bc(2,i,j,k) = 0.5d0*( Bs(2,i,j+1,k) + Bs(2,i,j,k) )
            Bc(3,i,j,k) = 0.5d0*( Bs(3,i,j,k+1) + Bs(3,i,j,k) )
        enddo
        enddo
        !$omp end parallel do
        enddo

return
end subroutine CellCenterMagneticField

!-------------------------------------------------------------------
!       determine dt 
!-------------------------------------------------------------------
real(8) function TimestepControl(xf, yf, zf, Q, Bc )
use params, only : is, ie, js, je, ks, ke, &
                   IDN, IVX, IVY, IVZ, IPR, gam, Ccfl
implicit none
real(8), intent(in) :: xf(:), yf(:), zf(:), Q(:,:,:,:), Bc(:,:,:,:)
real(8)::dtl1,dtl2,dtl3
real(8)::dtlocal,dtmin,cf
integer::i,j,k

        dtmin=1.0d90

        do k=ks,ke
     !$omp parallel do private(i,dtl1,dtl2,cf) reduction (min: dtmin)
        do j=js,je
        do i=is,ie
            cf = dsqrt( (gam*Q(IPR,i,j,k) + Bc(1,i,j,k)**2 + Bc(2,i,j,k)**2 + Bc(3,i,j,k)**2)/Q(IDN,i,j,k))
         
            dtl1 =(xf(i+1)-xf(i))/(abs(Q(IVX,i,j,k)) + cf)
            dtl2 =(yf(j+1)-yf(j))/(abs(Q(IVY,i,j,k)) + cf)
!            dtlocal = min(dtl1,dtl2)
!            if(dtlocal .lt. dtmin) dtmin = dtlocal
            dtmin = min(dtl1,dtl2,dtmin)
        enddo
        enddo
      !$omp end parallel do
        enddo

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
real(8) :: sgn
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
subroutine NumericalFlux( dt, xf, yf, zf, Q, Bc, Bs, F, G, H, E)
use params, only : nxtot, nytot, nztot, NFLX, NVAR, is, ie, js, je, ks, ke, &
                   IVX, IVY, IVZ, IBX, IBY, IBZ, flag_flux
implicit none
real(8), intent(in) :: dt,xf(:), yf(:), zf(:)
real(8), intent(in) :: Q(:,:,:,:)
real(8), intent(in) :: Bc(:,:,:,:)
real(8), intent(in) :: Bs(:,:,:,:)
real(8), intent(out) :: F(:,:,:,:)
real(8), intent(out) :: G(:,:,:,:)
real(8), intent(out) :: H(:,:,:,:)
real(8), intent(out) :: E(:,:,:,:)
    
integer::i,j,k
real(8),dimension(NFLX,nxtot,nytot,nztot):: Ql,Qr
real(8),dimension(NFLX):: flx
real(8) :: dQm(NVAR), dQp(NVAR), dQmon(NVAR)
real(8) :: ddmon, dvmon, dpmon
real(8) :: Qltest(NFLX), Qrtest(NFLX);
    
real(8),dimension(nxtot,nytot,nztot) :: e2_xf, e3_xf
real(8),dimension(nxtot,nytot,nztot) :: e1_yf, e3_yf
real(8),dimension(nxtot,nytot,nztot) :: weight1, weight2, weight3
real(8) :: wghtCT

real(8) :: Etmp(nxtot,nytot,nztot) 
real(8) :: de3_l1, de3_r1, de3_l2, de3_r2

    
!$omp parallel 
    ! numerical flux in the x direction
    ! hydro part
        do k=ks,ke
      !$omp do private( i, j, dQp, dQm, dQmon )
        do j=js-1,je+1
        do i=is-1,ie+1
            dQp(1:NVAR) = Q(1:NVAR,i+1,j,k) - Q(1:NVAR,i  ,j,k)
            dQm(1:NVAR) = Q(1:NVAR,i  ,j,k) - Q(1:NVAR,i-1,j,k)
    
            call vanLeer(NVAR, dQp, dQm, dQmon)
    
             ! Ql(i,j,k) --> W_(i-1/2,j,k)
             ! Qr(i,j,k) --> W_(i-1/2,j,k)
            Ql(1:NVAR,i+1,j,k) = Q(1:NVAR,i,j,k) + 0.5d0*dQmon(1:NVAR)
            Qr(1:NVAR,i  ,j,k) = Q(1:NVAR,i,j,k) - 0.5d0*dQmon(1:NVAR)

        enddo
        enddo
      !$omp end do
        enddo
    
        ! B field part
        do k=ks,ke
      !$omp do private( i, dQp, dQm, dQmon )
        do j=js-1,je+1
        do i=is-1,ie+1
            dQp(1:3) = Bc(1:3,i+1,j,k) - Bc(1:3,i  ,j,k)
            dQm(1:3) = Bc(1:3,i  ,j,k) - Bc(1:3,i-1,j,k)
    
            call vanLeer(3, dQp, dQm, dQmon)
    
             ! Ql(i,j,k) --> W_(i-1/2,j,k)
             ! Qr(i,j,k) --> W_(i-1/2,j,k)
            Ql(NVAR+1:NFLX,i+1,j,k) = Bc(1:3,i,j,k) + 0.5d0*dQmon(1:3)
            Qr(NVAR+1:NFLX,i  ,j,k) = Bc(1:3,i,j,k) - 0.5d0*dQmon(1:3)
        enddo
        enddo
      !$omp end do
        enddo
    
        if (flag_flux == 1 ) then
            do k=ks,ke
           !$omp do private( i, flx, wghtCT )
            do j=js-1,je+1
            do i=is,ie+1
                call HLL(1,Ql(:,i,j,k),Qr(:,i,j,k),Bs(1,i,j,k),dt,xf(i+1)-xf(i),flx,wghtCT)
        
                 F(1:NVAR,i,j,k)  = flx(1:NVAR)
                 e3_xf(i,j,k) =  -flx(IBY)
                 e2_xf(i,j,k) =  +flx(IBZ)

!                 if(flx(IDN).ne.flx(IDN)) then
!                     print*,i,j,k,flx(IDN),Ql(i,j,k,IDN),Qr(i,j,k,IDN)
!                 endif
        
                 weight1(i,j,k) = wghtCT
            enddo
            enddo
          !$omp end do
            enddo
        else if (flag_flux == 2 ) then
            do k=ks,ke
           !$omp do private( i, flx, wghtCT )
            do j=js-1,je+1
            do i=is,ie+1
                call HLLD(1,Ql(:,i,j,k),Qr(:,i,j,k),Bs(1,i,j,k),dt,xf(i+1)-xf(i),flx,wghtCT)
        
                 F(1:NVAR,i,j,k)  = flx(1:NVAR)
                 e3_xf(i,j,k) =  -flx(IBY)
                 e2_xf(i,j,k) =  +flx(IBZ)

!                 if(flx(IDN).ne.flx(IDN)) then
!                     print*,i,j,k,flx(IDN),Ql(i,j,k,IDN),Qr(i,j,k,IDN)
!                 endif
        
        
                 weight1(i,j,k) = wghtCT
            enddo
            enddo
          !$omp end do
            enddo
        end if 
    
    
          ! numerical flux in the y direction
          do k=ks,ke
      !$omp do private( i, dQp, dQm, dQmon )
          do j=js-1,je+1
          do i=is-1,ie+1
             dQp(1:NVAR) = Q(1:NVAR,i,j+1,k) - Q(1:NVAR,i,j  ,k)
             dQm(1:NVAR) = Q(1:NVAR,i,j  ,k) - Q(1:NVAR,i,j-1,k)
    
             call vanLeer(NVAR, dQp, dQm, dQmon)
    
             ! Ql(i,j,k) --> W_(i-1/2,j,k)
             ! Qr(i,j,k) --> W_(i-1/2,j,k)
             Ql(1:NVAR,i,j+1,k) = Q(1:NVAR,i,j,k) + 0.5d0*dQmon(1:NVAR)
             Qr(1:NVAR,i,j  ,k) = Q(1:NVAR,i,j,k) - 0.5d0*dQmon(1:NVAR)
          enddo
          enddo
          !$omp end do
          enddo
    
          ! B field part
          do k=ks,ke
      !$omp do private( i, dQp, dQm, dQmon )
          do j=js-1,je+1
          do i=is-1,ie+1
             dQp(1:3) = Bc(1:3,i,j+1,k) - Bc(1:3,i,j  ,k)
             dQm(1:3) = Bc(1:3,i,j  ,k) - Bc(1:3,i,j-1,k)
    
             call vanLeer(3, dQp, dQm, dQmon)
    
             ! Ql(i,j,k) --> W_(i-1/2,j,k)
             ! Qr(i,j,k) --> W_(i-1/2,j,k)
             Ql(NVAR+1:NFLX,i,j+1,k) = Bc(1:3,i,j,k) + 0.5d0*dQmon(1:3)
             Qr(NVAR+1:NFLX,i,j  ,k) = Bc(1:3,i,j,k) - 0.5d0*dQmon(1:3)
          enddo
          enddo
          !$omp end do
          enddo
    
        if( flag_flux == 1 ) then
          do k=ks,ke
         !$omp do private( i, flx, wghtCT )
          do j=js,je+1
          do i=is-1,ie+1
             call HLL(2,Ql(:,i,j,k),Qr(:,i,j,k),Bs(2,i,j,k),dt,yf(j+1) - yf(j), flx,wghtCT)
    
             G(1:NVAR,i,j,k) = flx(1:NVAR)
             e1_yf(i,j,k) =  - flx(IBZ)
             e3_yf(i,j,k) =    flx(IBX)
    
             weight2(i,j,k) = wghtCT
          enddo
          enddo
          !$omp end do
          enddo
        else if (flag_flux == 2 ) then
          do k=ks,ke
         !$omp do private( i, flx, wghtCT )
          do j=js,je+1
          do i=is-1,ie+1
             call HLLD(2,Ql(:,i,j,k),Qr(:,i,j,k),Bs(2,i,j,k),dt,yf(j+1) - yf(j), flx,wghtCT)

    
             G(1:NVAR,i,j,k) = flx(1:NVAR)
             e1_yf(i,j,k) =  - flx(IBZ)
             e3_yf(i,j,k) =    flx(IBX)
    
    
             weight2(i,j,k) = wghtCT
          enddo
          enddo
          !$omp end do
          enddo
        endif

!$omp end parallel
    
!       call ElectricField( Q, Bc, e2_xf, e3_xf, e3_yf, e1_yf, weight1, weight2, E )

!$omp parallel
      do k=ks,ke
      !$omp do private(i)
      do j=js-1, je+1
      do i=is-1, ie+1
           Etmp(i,j,k) = Q(IVY,i,j,k)*Bc(1,i,j,k) - Q(IVX,i,j,k)*Bc(2,i,j,k)
      enddo
      enddo
      !$omp end do 
      enddo

      !$omp do private(i)
      do j=js, je
      do i=is, ie+1
           E(2,i,j,ke+1) = e2_xf(i,j,ks)
           E(2,i,j,ks  ) = e2_xf(i,j,ks)
      enddo
      enddo
      !$omp end do 

      !$omp do private(i)
      do j=js, je+1
      do i=is, ie
           E(1,i,j,ke+1) = e1_yf(i,j,ks)
           E(1,i,j,ks  ) = e1_yf(i,j,ks)
      enddo
      enddo
      !$omp end do 


      do k=ks,ke
      !$omp do private(i, de3_l2, de3_r2, de3_l1, de3_r1 )
      do j=js, je+1
      do i=is, ie+1
          de3_l2 = (1.0d0-weight1(i,j-1,k))*(e3_yf(i  ,j,k) - Etmp(i  ,j-1,k)) + &
                   (      weight1(i,j-1,k))*(e3_yf(i-1,j,k) - Etmp(i-1,j-1,k))

          de3_r2 = (1.0d0-weight1(i,j  ,k))*(e3_yf(i  ,j,k) - Etmp(i  ,j  ,k)) + &
                   (      weight1(i,j  ,k))*(e3_yf(i-1,j,k) - Etmp(i-1,j  ,k))

          de3_l1 = (1.0d0-weight2(i-1,j,k))*(e3_xf(i,j  ,k) - Etmp(i-1,j  ,k)) + &
                   (      weight2(i-1,j,k))*(e3_xf(i,j-1,k) - Etmp(i-1,j-1,k))

          de3_r1 = (1.0d0-weight2(i  ,j,k))*(e3_xf(i,j  ,k) - Etmp(i  ,j  ,k)) + &
                   (      weight2(i  ,j,k))*(e3_xf(i,j-1,k) - Etmp(i  ,j-1,k))

         E(3,i,j,k) = 0.25d0*( de3_l1 + de3_r1 + de3_l2 + de3_r2 + &
                          e3_yf(i-1,j,k) + e3_yf(i,j,k) + e3_xf(i,j-1,k) + e3_xf(i,j,k))
      enddo
      enddo
      !$omp end do 
      enddo

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
subroutine HLL(idir,Ql,Qr,b1,dt,dx,flx,wghtCT)
use params, only : IDN, IVX, IVY, IVZ, IPR, IBX, IBY, IBZ, IEN, NFLX, gam
implicit none
integer, intent(in) :: idir
real(8),intent(in)  :: Ql(NFLX), Qr(NFLX)
real(8),intent(in)  :: dt,dx, b1
real(8),intent(out) :: flx(NFLX), wghtCT
integer :: IVpara, IVperp1, IVperp2
integer :: IBpara, IBperp1, IBperp2
real(8):: Ul(NFLX), Ur(NFLX)
real(8):: Fl(NFLX), Fr(NFLX)
real(8):: Ust(NFLX)
real(8):: Fst(NFLX)
real(8):: cfl,cfr
real(8):: sl, sr
real(8):: pbl, pbr, ptotl, ptotr, v_over_c
integer :: i, n

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

!          sl = min(Ql(IVX),Qr(IVX)) - max(cfl,cfr)
!          sr = max(Ql(IVX),Qr(IVX)) + max(cfl,cfr)
          sl = min(Ql(IVpara) - cfl,Qr(IVpara) - cfr)
          sr = max(Ql(IVpara) + cfl,Qr(IVpara) + cfr)
!          Fst(:)  = (sr*Fl(:) - sl*Fr(:) + sl*sr*( Ur(:) - Ul(:) ))/(sr - sl)
!          Ust(:) = ( sr*Ur(:) - sl*Ul(:) - Fr(:) + Fl(:) )/(sr - sl)

          if( sl > 0.0d0 ) then
               flx(:) = Fl(:)
          else if (sr <= 0.0d0 ) then
               flx(:) = Fr(:)
          else 
               flx(:)  = (sr*Fl(:) - sl*Fr(:) + sl*sr*( Ur(:) - Ul(:) ))/(sr - sl)
          endif

           v_over_c = 1024.0d0*dt*flx(IDN)/( dx*( Ql(IDN) + Qr(IDN) ) )
           wghtCT = 0.5d0 + max( -0.5d0, min(0.5d0, v_over_c) )

!          flx(IPS) = 0.0d0

!         do i=1,NFLX1D 
!         if( flx(i) .ne. flx(i) ) then 
!             print*,(sr*Fl(:) - sl*Fr(:) + sl*sr*( Ur(:) - Ul(:) ))/(sr - sl)
!         stop
!         endif
!         enddo
!
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
subroutine HLLD(idir,Ql,Qr,b1,dt,dx,flx,wghtCT)
use params, only : IDN, IVX, IVY, IVZ, IPR, IBX, IBY, IBZ, IEN, NFLX, gam
implicit none
integer, intent(in) :: idir
real(8),intent(in)  :: dt
real(8),intent(in)  :: Ql(NFLX), Qr(NFLX)
real(8),intent(in)  :: dx
real(8),intent(in) :: b1
real(8),intent(out) :: flx(NFLX)
real(8),intent(out) :: wghtCT
integer :: IVpara, IVperp1, IVperp2
integer :: IBpara, IBperp1, IBperp2, id
real(8):: Ul(NFLX), Ur(NFLX)
real(8):: Ulst(NFLX), Urst(NFLX)
real(8):: Uldst(NFLX), Urdst(NFLX)
real(8):: Fl(NFLX), Fr(NFLX)
real(8):: test(NFLX)
real(8):: cfl,cfr
real(8):: S0, S1, S2, S3, S4
real(8):: pbl, pbr, ptotl, ptotr
real(8) :: sqrtdl, sqrtdr, v_dot_B_stl, v_dot_B_str
real(8) :: Ulst_d_inv, Urst_d_inv, sum_sqrtd_inv, tmp
real(8) :: ptot_stl, ptot_str,ptot_st, Cl, Cr, Cml, Cmr, Cml_inv, Cmr_inv, bxsgn
real(8) :: v_over_c
integer :: i, n
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

!         test = (S4 - S3)*Urst + (S3 - S2)*Urdst + (S2 - S1)*Uldst + (S1 - S0)*Ulst - S4*Ur + S0*Ul + Fr - Fl
!         do id=1,NFLX
!         if( abs(test(id)) > 1.0d-10 ) then
!         print*,test(id)
!         stop
!         endif
!         enddo
         

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

           v_over_c = 1024.0d0*dt*flx(IDN)/( dx*( Ql(IDN) + Qr(IDN) ) )
           wghtCT = 0.5d0 + max( -0.5d0, min(0.5d0, v_over_c) )

return
end subroutine HLLD
!---------------------------------------------------------------------
!     ElectricField
!---------------------------------------------------------------------
!     computes the numerical flux at the cell boundary 
!
!     Input: Q: primitive variables at the cell center
!
!     Input: B: magnetic fields
!
!     Output: flux : the numerical flux estimated at the cell boundary
!---------------------------------------------------------------------
subroutine ElectricField( Q, Bc, e2_xf, e3_xf, e3_yf, e1_yf, weight1, weight2, E )
use params, only : nxtot, nytot, nztot, is, ie, js, je, ks, ke, IVX, IVY, IVZ
implicit none
integer::i,j,k
real(8), intent(in)  :: Q(:,:,:,:), Bc(:,:,:,:)
real(8), intent(in)  :: e2_xf(:,:,:)
real(8), intent(in)  :: e3_xf(:,:,:)
real(8), intent(in)  :: e1_yf(:,:,:)
real(8), intent(in)  :: e3_yf(:,:,:)
real(8), intent(in)  :: weight1(:,:,:)
real(8), intent(in)  :: weight2(:,:,:)
real(8), intent(out) :: E(:,:,:,:)
real(8) :: Etmp(nxtot,nytot,nztot) 
real(8) :: de3_l1, de3_r1, de3_l2, de3_r2

!$omp parallel
      do k=ks,ke
      !$omp do private(i)
      do j=js-1, je+1
      do i=is-1, ie+1
           Etmp(i,j,k) = Q(IVY,i,j,k)*Bc(1,i,j,k) - Q(IVX,i,j,k)*Bc(2,i,j,k)
      enddo
      enddo
      !$omp end do 
      enddo

      !$omp do private(i)
      do j=js, je
      do i=is, ie+1
           E(2,i,j,ke+1) = e2_xf(i,j,ks)
           E(2,i,j,ks  ) = e2_xf(i,j,ks)
      enddo
      enddo
      !$omp end do 

      !$omp do private(i)
      do j=js, je+1
      do i=is, ie
           E(1,i,j,ke+1) = e1_yf(i,j,ks)
           E(1,i,j,ks  ) = e1_yf(i,j,ks)
      enddo
      enddo
      !$omp end do 


      do k=ks,ke
      !$omp do private(i, de3_l2, de3_r2, de3_l1, de3_r1 )
      do j=js, je+1
      do i=is, ie+1
          de3_l2 = (1.0d0-weight1(i,j-1,k))*(e3_yf(i  ,j,k) - Etmp(i  ,j-1,k)) + &
                   (      weight1(i,j-1,k))*(e3_yf(i-1,j,k) - Etmp(i-1,j-1,k))

          de3_r2 = (1.0d0-weight1(i,j  ,k))*(e3_yf(i  ,j,k) - Etmp(i  ,j  ,k)) + &
                   (      weight1(i,j  ,k))*(e3_yf(i-1,j,k) - Etmp(i-1,j  ,k))

          de3_l1 = (1.0d0-weight2(i-1,j,k))*(e3_xf(i,j  ,k) - Etmp(i-1,j  ,k)) + &
                   (      weight2(i-1,j,k))*(e3_xf(i,j-1,k) - Etmp(i-1,j-1,k))

          de3_r1 = (1.0d0-weight2(i  ,j,k))*(e3_xf(i,j  ,k) - Etmp(i  ,j  ,k)) + &
                   (      weight2(i  ,j,k))*(e3_xf(i,j-1,k) - Etmp(i  ,j-1,k))

         E(3,i,j,k) = 0.25d0*( de3_l1 + de3_r1 + de3_l2 + de3_r2 + &
                          e3_yf(i-1,j,k) + e3_yf(i,j,k) + e3_xf(i,j-1,k) + e3_xf(i,j,k))
      enddo
      enddo
      !$omp end do 
      enddo

!$omp end parallel

return
end subroutine ElectricField 
!!=====================================================================
subroutine UpdateConsv( dt1, xf, xv, yf, yv, zf, zv, F, G, H, E, Q, Uo, Bso, U, Bs)
use params, only : is, ie, js, je, ks, ke
implicit none
real(8), intent(in) :: dt1
real(8), intent(in)  :: xf(:), yf(:), zf(:)
real(8), intent(in)  :: xv(:), yv(:), zv(:)
real(8), intent(in)  :: F(:,:,:,:), G(:,:,:,:), H(:,:,:,:)
real(8), intent(in)  :: Uo(:,:,:,:), Q(:,:,:,:), Bso(:,:,:,:)
real(8), intent(out) :: U(:,:,:,:), Bs(:,:,:,:), E(:,:,:,:)

real(8) :: divB 
integer::i,n,j,k

!$omp parallel 

      !$omp do private( i,j,k )
      do k=ks,ke
      do j=js,je
      do i=is,ie
         U(:,i,j,k) = Uo(:,i,j,k) + dt1*(- F(:,i+1,j,k) + F(:,i,j,k))/(xf(i+1)-xf(i)) &
                                  + dt1*(- G(:,i,j+1,k) + G(:,i,j,k))/(yf(j+1)-yf(j)) 
      enddo
      enddo
      enddo
      !$omp end do

      do k=ks,ke
      !$omp do private( i )
      do j=js,je
      do i=is,ie+1
           Bs(1,i,j,k) = Bso(1,i,j,k) &
                       - dt1*(E(3,i,j+1,k) - E(3,i,j,k))/(yv(j+1) - yv(j))
      enddo
      enddo
      !$omp end do
      enddo

      do k=ks,ke
      !$omp do private( i )
      do j=js,je+1
      do i=is,ie
           Bs(2,i,j,k) = Bso(2,i,j,k) &
                       + dt1*(E(3,i+1,j,k) - E(3,i,j,k))/(xv(i+1) - xv(i))
      enddo
      enddo
      !$omp end do
      enddo

      do k=ks,ke+1
      !$omp do private( i )
      do j=js,je
      do i=is,ie
           Bs(3,i,j,k) = Bso(3,i,j,k) &
                       - dt1*(E(2,i+1,j,k) - E(2,i,j,k))/(xf(i+1) - xf(i)) &
                       + dt1*(E(1,i,j+1,k) - E(1,i,j,k))/(yf(j+1) - yf(j))
      enddo
      enddo
      !$omp end do
      enddo


      !$omp end parallel

return
end subroutine UpdateConsv
!-------------------------------------------------------------------
!       Output snapshot files 
!       Input  : flag, dirname, xf, xv, Q
!
!       flag = .true.  --> output snapshot when calling this subroutine
!       flag = .false. --> output snapshot every dtsnap
!-------------------------------------------------------------------
subroutine Output( time, flag, dirname, xf, xv, yf, yv, Q, Bc )
use params, only : nx, ny, nz, is, ie, js, je, ks, ke, NFLX, NVAR, unitsnap, unitbin, &
                  IDN, IVX, IVY, IVZ, IPR, dtsnap, flag_binary
implicit none
real(8), intent(in) :: time
logical, intent(in) :: flag 
character(20), intent(in) :: dirname 
real(8), intent(in) :: xf(:), xv(:), yf(:), yv(:)
real(8), intent(in) :: Q(:,:,:,:), Bc(:,:,:,:)
integer::i,j,k
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
        write(unitsnap) NVAR
        write(unitsnap) NFLX-NVAR
        write(unitsnap) xv(is:ie)
        write(unitsnap) yv(js:je)
        write(unitsnap) real(Q(1:NVAR,is:ie,js:je,ks:ke)) ! single precision
        write(unitsnap) real(Bc(1:3,is:ie,js:je,ks:ke)) ! single precision
        close(unitsnap)
    else 
          filename = trim(dirname)//"/snap"//trim(filename)//".dat"
          open(unitbin,file=filename,form='formatted',action="write")
          write(unitbin,*) "# time = ",time
          write(unitbin,*) "#nx, ny = ", nx, ny
          do k=ks,ke
          do j=js,je
          do i=is,ie
              write(unitbin,*) xv(i), yv(j), Q(IDN,i,j,k), Q(IVX,i,j,k), Q(IVY,i,j,k), Q(IVZ,i,j,k), Q(IPR,i,j,k), &
                  Bc(1,i,j,k), Bc(2,i,j,k), Bc(3,i,j,k)
          enddo
          enddo
          enddo
     endif

     close(unitbin)

     write(6,*) "output:  ",filename,time

      nsnap=nsnap+1
      tsnap=tsnap + dtsnap

return
end subroutine Output
!-------------------------------------------------------------------
!       create directory
!       Input  : the directory to be created
!-------------------------------------------------------------------
subroutine makedirs(outdir)
implicit none
character(len=*), intent(in) :: outdir
character(len=256) command
write(command, *) 'if [ ! -d ', trim(outdir), ' ]; then mkdir -p ', trim(outdir), '; fi'
!      write(*, *) trim(command)
call system(command)
end subroutine makedirs
!-------------------------------------------------------------------
!       Realtime Analysis
!       Input  : xf, xv
!       Output : phys_evo(nevo)
!-------------------------------------------------------------------
subroutine RealtimeAnalysis(xv,yv,Q,Bc,Bs,phys_evo)
use params, only : is, ie, js, je, ks, ke, gam
implicit none
real(8), intent(in) :: xv(:), yv(:), Q(:,:,:,:), Bc(:,:,:,:), Bs(:,:,:,:)
real(8), intent(out) :: phys_evo(:)
integer::i,j,k
real(8) :: dvy, er_divBc, er_divBs

      do k=ks,ke
      do j=js,je
      do i=is,ie
      enddo
      enddo
      enddo
      phys_evo(1) = 0.0d0
      phys_evo(2) = 0.0d0
      phys_evo(3) = 0.0d0
      
return
end subroutine RealtimeAnalysis

