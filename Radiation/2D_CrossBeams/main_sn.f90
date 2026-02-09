!> @file
!! @brief
!===============================================================================
!! 2D Boltzmann transport equation (discrete ordinates S_N in 2D angle), IMEX time stepping
!!   df_m/dt + c*(mux_m * df_m/dx + muy_m * df_m/dy) = rho*kappa*c*(f_LTE - f_m)
!!
!! @details
!! Spatial advection: 5th-order WENO (LF flux splitting), explicit (IMEX)
!! Absorption+emission (LTE relaxation): implicit (local)
!!   df/dt = rho*kappa*c*(f_LTE - f)
!!   => f^{n+1} = (f* + dt*rho*kappa*c*f_LTE) / (1 + dt*rho*kappa*c)
!!
!! @author Tomoya Takiwaki
!! @date 2026-02-06
!! @version 1.0
!===============================================================================

module constants
  implicit none
  real(8), parameter :: pi = acos(-1.d0)
  real(8), parameter :: cl = 1.d0  ! cm/s (use any consistent units)
end module constants

module modelpara
  implicit none
  real(8),parameter:: rho0  = 1.0d0 ! [g/cm^3]
  real(8),parameter:: kap0  = 0.0d0 ! [cm^2/g]
  real(8),parameter:: erad0 = 1.0d-4 ! [erg/cm^3]
  real(8),parameter:: erad1 = 1.0d0 ! [erg/cm^3]
end module modelpara

module weno5
  implicit none
contains
  pure real(8) function weno5_plus(vmm, vm, v0, vp, vpp) result(vh)
    real(8), intent(in) :: vmm, vm, v0, vp, vpp
    real(8) :: b0,b1,b2,a0,a1,a2,w0,w1,w2,eps
    real(8) :: p0,p1,p2,den
    eps = 1.d-40

    p0 = ( 2.d0*vmm - 7.d0*vm + 11.d0*v0)/6.d0
    p1 = (-1.d0*vm  + 5.d0*v0 +  2.d0*vp)/6.d0
    p2 = ( 2.d0*v0  + 5.d0*vp -  1.d0*vpp)/6.d0

    b0 = (13.d0/12.d0)*(vmm-2.d0*vm+v0)**2 + 0.25d0*(vmm-4.d0*vm+3.d0*v0)**2
    b1 = (13.d0/12.d0)*(vm-2.d0*v0+vp)**2  + 0.25d0*(vm-vp)**2
    b2 = (13.d0/12.d0)*(v0-2.d0*vp+vpp)**2 + 0.25d0*(3.d0*v0-4.d0*vp+vpp)**2

    a0 = 0.1d0/( (eps+b0)**2 )
    a1 = 0.6d0/( (eps+b1)**2 )
    a2 = 0.3d0/( (eps+b2)**2 )

    den = a0+a1+a2
    w0 = a0/den; w1 = a1/den; w2 = a2/den
    vh = w0*p0 + w1*p1 + w2*p2
  end function weno5_plus

  pure real(8) function weno5_minus(vmm, vm, v0, vp, vpp) result(vh)
    real(8), intent(in) :: vmm, vm, v0, vp, vpp
    ! right-biased reconstruction via symmetry (reverse stencil)
    vh = weno5_plus(vpp, vp, v0, vm, vmm)
  end function weno5_minus
end module weno5

module gridmod
  implicit none
  integer, parameter :: mgn = 3
  integer, parameter :: izones = 200
  integer, parameter :: jzones = 200
  integer, parameter :: kzones = 1
  integer, parameter :: in = izones + 2*mgn +1
  integer, parameter :: jn = jzones + 2*mgn +1
  integer, parameter :: kn = 1
  integer, parameter :: is = 1 + mgn
  integer, parameter :: ie = izones + mgn
  integer, parameter :: js = 1 + mgn
  integer, parameter :: je = jzones + mgn
  integer, parameter :: ks = 1
  integer, parameter :: ke = 1

  real(8), parameter :: x1min=0.d0, x1max=5.0d0
  real(8), parameter :: x2min=0.d0, x2max=5.0d0
  real(8), parameter :: x3min=0.d0, x3max=0.0d0

  real(8):: dx ! (xmax-xmin)/izones
  real(8):: dy ! (ymax-ymin)/jzones
  real(8):: dz ! z-grid spacing (unused here; kn=1, effectively 2D)
  
  real(8),dimension(in)::x1a,x1b
  real(8),dimension(jn)::x2a,x2b
  real(8),dimension(kn)::x3a,x3b

contains
  subroutine GenerateGrid
    implicit none
    integer::i,j,k
    dx=(x1max-x1min)/izones
    do i=1,in
       x1a(i) = dx*(i-(mgn+1))+x1min
    enddo
    do i=1,in-1
       x1b(i) = 0.5d0*(x1a(i+1)+x1a(i))
    enddo

    dy=(x2max-x2min)/jzones
    do j=1,jn
       x2a(j) = dy*(j-(mgn+1))+x2min
    enddo
    do j=1,jn-1
       x2b(j) = 0.5d0*(x2a(j+1)+x2a(j))
    enddo
    
    return
  end subroutine GenerateGrid
end module gridmod

!> Angular grid and quadrature for 2D discrete ordinates (S_N).
!! Directions are specified by (mux,muy) on the unit circle with weights dth.
module angmod
  use constants, only: pi
  implicit none
  ! 2D distribution function
  ! fdis [erg/cm^3/rad]
  integer, parameter :: mang = 4*4 !  
  real(8), dimension(mang) :: dth
  real(8), dimension(mang+1) :: mua
  real(8), dimension(mang) :: mux, muy
  integer,parameter:: grid_aligned=1, staggered=2 
  integer,parameter:: angsetup = grid_aligned
contains
  subroutine GenerateAngleGrid
    implicit none
    integer :: m
    real(8) :: th
    select case(angsetup)
    case(grid_aligned)
       do m=1,mang+1
          mua(m) = 2.0d0*pi * (dble(m)-1.5d0)/dble(mang) !! Angle grid edges in [0, 2*pi]
       enddo
       do m=1,mang
          th = 2.d0*pi*(dble(m)-1.0d0)/dble(mang)
          mux(m) = cos(th)
          muy(m) = sin(th)
          dth(m) = mua(m+1)-mua(m)
       enddo
    case(staggered)
       do m=1,mang+1
          mua(m) = 2.0d0*pi * (m-1) /mang !! Angle grid edges in [0, 2*pi]
       enddo
       do m=1,mang
          th = 2.d0*pi*(m-0.5d0)/mang
          mux(m) = cos(th)
          muy(m) = sin(th)
          dth(m) = mua(m+1)-mua(m)
       enddo
    end select
  end subroutine  GenerateAngleGrid
end module angmod

!> State variables (distribution function, sources, fluxes, and moments).
module statemod
  use gridmod
  use angmod
  implicit none
  real(8), dimension(mang,in,jn,kn):: fdis, fstar
  !! fdis : angular distribution function (one value per discrete direction)
  !! fstar : the half updated angular distribution function, advection effect is only considered
  real(8), dimension(mang,in,jn,kn):: fx, fy   ! interface fluxes
  !! fx,fy : interface fluxes reconstructed by WENO5 (stored on the cell index)   ! interface fluxes; stored at cell indices
  real(8),dimension(mang,in,jn,kn):: flte
  
  real(8),dimension(in,jn,kn):: d,ei,tempK
  real(8),dimension(in,jn,kn):: kappa
  !! kappa : absorption opacity [cm^2/g] (used via rho*kappa)
  real(8):: eimin,fradmin

  real(8),dimension(    in,jn,kn)::Erad
  !! Erad : angular integral (energy-like density), E = sum_m dth(m) * fdis_m
  integer,parameter:: xdir=1,ydir=2,zdir=3
  real(8),dimension(xdir:zdir,in,jn,kn)::Frad
  integer,parameter:: nxx=1,nyy=2,nzz=3,nxy=4,nyz=5,nzx=6
  real(8),dimension(nxx:nzx,in,jn,kn)::Prad
contains
  subroutine GenerateProblem()
    use constants
    use gridmod
    use modelpara
    implicit none
    integer :: i,j,k,m

    do k=ks,ke
    do j=js-mgn,je+mgn
    do i=is-mgn,ie+mgn
       d(i,j,k) = rho0
       kappa(i,j,k) = kap0
       ei(i,j,k)   = 0.0d0
    enddo
    enddo
    enddo
    fradmin = 1.0d10
    do k=ks,ke
    do j=js,je
    do i=is,ie
        do m=1,mang
           fdis(m,i,j,k) = erad0/dth(m)
           flte(m,i,j,k) = 0.0d0 ! [erg/cm^3]
           fradmin = min (fradmin, fdis(m,i,j,k))
       enddo
    enddo
    enddo
    enddo
 
  end subroutine GenerateProblem
end module statemod


subroutine RadBoundaryCondition
  use constants
  use statemod
  use modelpara
  implicit none
  integer::i,j,k,m

  do k=ks,ke
  do j=js,je
  do i=1,mgn
     do m=1,mang
        if(mux(m) >= 0.0d0) then
           fdis(m     ,is-i,j,k) = erad0/dth(m) ! fixed
        else
           fdis(m     ,is-i,j,k) = fdis(m,is,j,k) ! outflos
        endif
     enddo
  enddo
  enddo
  enddo

  do k=ks,ke
  do j=js,je
  do i=1,mgn
     do m=1,mang
        if(mux(m) > 0.0d0) then
           fdis(m     ,ie+i,j,k) = fdis(m,is,j,k) ! outflos
        else
           fdis(m     ,ie+i,j,k) = erad0/dth(m) ! fixed
        endif
     enddo
  enddo
  enddo
  enddo

  ! j-direction

  select case(angsetup)
  case(grid_aligned)
     do k=ks,ke
     do i=is,ie
     do j=1,mgn
        if(x1b(i) < 1.0d0 )then
           do m=1,mang
              fdis(m,i,js-j,k) = erad0/dth(m)
              if (m == 3 ) fdis(m,i,js-j,k) = erad1/dth(m)
           enddo
        else if(x1b(i) > 4.0d0 )then
           do m=1,mang
              fdis(m,i,js-j,k) = erad0/dth(m)
              if (m == 7 ) fdis(m,i,js-j,k) = erad1/dth(m)
           enddo
        else
           do m=1,mang
              if(muy(m) >= 0.0d0) then
                 fdis(m,i,js-j,k) = erad0/dth(m) ! fixed
              else
                 fdis(m,i,js-j,k) = fdis(m,i,js,k) ! outflow
              endif
           enddo
        endif
     enddo
     enddo
     enddo

  case(staggered)
     do k=ks,ke
     do i=is,ie
     do j=1,mgn
        do m=1,mang
           fdis(m,i,js-j,k) = fdis(m,i,js,k) ! adhoc
        enddo
     enddo
     enddo
     enddo

  end select
  
  do k=ks,ke
  do i=is,ie
  do j=1,mgn
     do m=1,mang
        if(muy(m) > 0.0d0) then
           fdis(m,i,je+j,k) = fdis(m,i,je-j+1,k) ! outflow
        else
           fdis(m,i,je+j,k) = erad0/dth(m) ! fixed
        endif
     enddo
   enddo
   enddo
   enddo

   return
end subroutine RadBoundaryCondition

module fluxmod
  use gridmod
  use angmod
  use statemod, only: fdis, fx, fy
  use constants, only: cl
  use weno5
  implicit none
contains
  !> Compute x-direction numerical fluxes using WENO5 + LF splitting.
  !! The flux array `fx` is defined at cell interfaces (i+1/2) and stored at index i:
  !! `fx(m,i,j,k) = F_x(i+1/2,j,k; direction m)`.
subroutine RadFlux1
    use omp_lib
    integer :: i,j,k,m
    real(8) :: a, alpha
    real(8) :: fp_mm,fp_m,fp_0,fp_p,fp_pp
    real(8) :: fm_mm,fm_m,fm_0,fm_p,fm_pp

    ! fx(m,i,j,k) is flux at i+1/2 (stored at i)

!$omp parallel do default(shared) private(i,j,k,m,a,alpha,fp_mm,fp_m,fp_0,fp_p,fp_pp,fm_mm,fm_m,fm_0,fm_p,fm_pp) collapse(3) schedule(static)
    do k=ks,ke
    do j=js,je
    do i=is-1,ie
       do m=1,mang
          a = cl*mux(m); alpha = abs(a)

          fp_mm = 0.5d0*(a+alpha)*fdis(m,i-2,j,k)
          fp_m  = 0.5d0*(a+alpha)*fdis(m,i-1,j,k)
          fp_0  = 0.5d0*(a+alpha)*fdis(m,i  ,j,k)
          fp_p  = 0.5d0*(a+alpha)*fdis(m,i+1,j,k)
          fp_pp = 0.5d0*(a+alpha)*fdis(m,i+2,j,k)

          ! Right-biased stencil for the negative flux at i+1/2:
          fm_mm = 0.5d0*(a-alpha)*fdis(m,i-1,j,k)
          fm_m  = 0.5d0*(a-alpha)*fdis(m,i  ,j,k)
          fm_0  = 0.5d0*(a-alpha)*fdis(m,i+1,j,k)
          fm_p  = 0.5d0*(a-alpha)*fdis(m,i+2,j,k)
          fm_pp = 0.5d0*(a-alpha)*fdis(m,i+3,j,k)

          fx(m,i,j,k) =   weno5_plus(fp_mm,fp_m,fp_0,fp_p,fp_pp)  &
                       + weno5_minus(fm_mm,fm_m,fm_0,fm_p,fm_pp)
       enddo
    enddo
    enddo
    enddo
!$omp end parallel do

  end subroutine RadFlux1

  !> Compute y-direction numerical fluxes using WENO5 + LF splitting.
  !! The flux array `fy` is defined at cell interfaces (j+1/2) and stored at index j:
  !! `fy(m,i,j,k) = F_y(i,j+1/2,k; direction m)`.
subroutine RadFlux2
    use omp_lib
    integer :: i,j,k,m
    real(8) :: b, beta
    real(8) :: gp_mm,gp_m,gp_0,gp_p,gp_pp
    real(8) :: gm_mm,gm_m,gm_0,gm_p,gm_pp

    ! fy(m,i,j,k) is flux at j+1/2 (stored at j)

!$omp parallel do default(shared) private(i,j,k,m,b,beta,gp_mm,gp_m,gp_0,gp_p,gp_pp,gm_mm,gm_m,gm_0,gm_p,gm_pp) collapse(3) schedule(static)
    do k=ks,ke
    do i=is,ie
    do j=js-1,je
       do m=1,mang
          b = cl*muy(m); beta = abs(b)

          gp_mm = 0.5d0*(b+beta)*fdis(m,i,j-2,k)
          gp_m  = 0.5d0*(b+beta)*fdis(m,i,j-1,k)
          gp_0  = 0.5d0*(b+beta)*fdis(m,i,j  ,k)
          gp_p  = 0.5d0*(b+beta)*fdis(m,i,j+1,k)
          gp_pp = 0.5d0*(b+beta)*fdis(m,i,j+2,k)

          ! Right-biased stencil for the negative flux at j+1/2:
          gm_mm = 0.5d0*(b-beta)*fdis(m,i,j-1,k)
          gm_m  = 0.5d0*(b-beta)*fdis(m,i,j  ,k)
          gm_0  = 0.5d0*(b-beta)*fdis(m,i,j+1,k)
          gm_p  = 0.5d0*(b-beta)*fdis(m,i,j+2,k)
          gm_pp = 0.5d0*(b-beta)*fdis(m,i,j+3,k)

          fy(m,i,j,k) =   weno5_plus(gp_mm,gp_m,gp_0,gp_p,gp_pp)  &
                       + weno5_minus(gm_mm,gm_m,gm_0,gm_p,gm_pp)
       enddo
    enddo
    enddo
    enddo
!$omp end parallel do

  end subroutine RadFlux2
end module fluxmod

module timemod
  use gridmod
  use angmod
  use statemod
  use constants, only: cl
  implicit none
  integer:: ntime
  integer,parameter::ntimemax=100000
  real(8)::time,dt
  real(8),parameter:: cfl=0.1d0
  data time / 0.0d0 /
  real(8),parameter:: timemax=20.0d0
  real(8),parameter:: dtout=timemax/100
contains
  subroutine TimeStepControl
  real(8)::dtsrc
  real(8)::dtl1
  real(8)::dtl2
  real(8)::dtl3
  real(8)::dtlocal
  real(8)::dtmin
  integer::i,j,k
  dtmin=1.0d90
  do k=ks,ke
  do j=js,je
  do i=is,ie
     dtsrc = 1.0d0/(d(i,j,k)*kappa(i,j,k)*cl) 
     dtl1  = (x1a(i+1)-x1a(i))/cl
     dtl2  = (x2a(j+1)-x2a(j))/cl
!         dtl3 =(x1a(i+1)-x1a(i))/(abs(v1(i,j,k)) +cs(i,j,k))
!         dtlocal = min (dtl1,dtl2)
!         if(i.eq.is)write(6,*) "dt",dtsrc,dtl1
     dtlocal = min(dtl1,dtl2)
!     dtlocal = min(dtl1,dtl2,dtsrc)
     if(dtlocal .lt. dtmin) dtmin = dtlocal
  enddo
  enddo
  enddo
  dt = cfl * dtmin
!      write(6,*)"dt",dt
  end subroutine TimeStepControl
end module timemod

module stepmod
  use timemod
  use gridmod
  use angmod
  use statemod
  use constants, only: cl
  implicit none
contains
  !> Explicit transport update (advection only): computes f^{*} from flux divergence.
  !! Uses interface fluxes `fx` and `fy` (WENO5 reconstruction) and advances
  !! `fdis -> fdis - dt*(dFdx + dFdy)`; the source term is applied separately.
subroutine UpdateRad
    use omp_lib
    integer :: i,j,k,m

!$omp parallel do default(shared) private(i,j,k,m) collapse(3) schedule(static)
    do k=ks,ke
    do j=js,je
    do i=is,ie
       do m=1,mang
          fstar(m,i,j,k) = fdis(m,i,j,k) &
               - dt/dx * ( fx(m,i,j,k) - fx(m,i-1,j,k) ) &
               - dt/dy * ( fy(m,i,j,k) - fy(m,i,j-1,k) )
       enddo
    enddo
    enddo
    enddo
!$omp end parallel do

  end subroutine UpdateRad

  !> Local absorption/emission source term (IMEX style).
  !! This routine computes a relaxation toward LTE in each cell and direction.
  !! The stiffness parameter is `alpha = rho*kappa*c*dt`.
  !! IMPORTANT: this matches the "Kirchhoff" form (absorption + emission)
  !! used in many M1 implementations at the moment level.
subroutine UpdateRadWithSource
    use omp_lib
    use modelpara
    use constants
    integer:: i,j,k,m
    real(8):: alpha
!$omp parallel do default(shared) private(i,j,k,m,alpha) collapse(3) schedule(static)
    do k=ks,ke
    do j=js,je
    do i=is,ie
       kappa(i,j,k) = kap0/d(i,j,k)    ! [cm^2/g] absorption opacity per mass alpha = d(i,j,k)*kappa(i,j,k) *cl*dt
       alpha = d(i,j,k)*kappa(i,j,k) *cl*dt
       do m=1,mang
          flte(m,i,j,k) = 0.0d0 ! [erg/cm^3/rad]
          fdis(m,i,j,k) = (fstar(m,i,j,k) + alpha* flte(m,i,j,k))/(1.0d0+alpha)
       enddo
    enddo
    enddo
    enddo
!$omp end parallel do
  end subroutine UpdateRadWithSource
end module stepmod

subroutine Output(flag_force,flag_binary,dirname)
  use gridmod
  use timemod
  use angmod
  use statemod, only: fdis, Erad,Frad
  implicit none
  logical,intent(in):: flag_force
  logical,intent(in):: flag_binary
  character(len=*),intent(in):: dirname
  integer::i,j,k,m
  character(40)::filename
  real(8),save::tout
  data tout / 0.0d0 /
  integer,save::nout
  data nout / 1 /
  integer:: unitbin,unitasc
  integer,parameter:: gs = 1 !! we do not have to write all ghost zone
  integer,parameter:: nrad=mang
  integer,parameter:: nhyd=2
  real(8)::x1out(is-gs:ie+gs,2)
  real(8)::x2out(js-gs:je+gs,2)
  real(8)::thout(1:mang,3)
  real(4)::radout(mang,is-gs:ie+gs,js-gs:je+gs,ks)
  real(4)::hydout(is-gs:ie+gs,js-gs:je+gs,ks,nhyd)

  logical, save:: is_inited
  data is_inited /.false./

  if (.not. is_inited) then
     call makedirs(dirname)
     is_inited =.true.
  endif

  if(.not. flag_force .and. time < tout+dtout) return

  if (flag_binary) then
     write(filename,'(a3,i5.5,a4)')"unf",nout,".dat"
     filename = trim(dirname)//filename
     
     x1out(is:ie,1) = x1b(is:ie)
     x1out(is:ie,2) = x1a(is:ie)

     x2out(js:je,1) = x2b(js:je)
     x2out(js:je,2) = x2a(js:je)
      
     thout(1:mang,1) = dth(1:mang)
     thout(1:mang,2) = mux(1:mang)
     thout(1:mang,3) = muy(1:mang)
     
     radout(1:mang,is:ie,js:je,ks) = fdis(1:mang,is:ie,js:je,ks)
     hydout(is:ie,js:je,ks,1) =    d(  is:ie,js:je,ks)
     hydout(is:ie,js:je,ks,2) =   ei(  is:ie,js:je,ks)
      
     write(filename,'(a4,i5.5,a4)')"snap",nout,".bin"
     filename = trim(dirname)//filename
     open(newunit=unitbin,file=filename,status='replace',form='binary') 
     write(unitbin) time
     write(unitbin) izones
     write(unitbin) jzones
     write(unitbin) mang,nhyd
     write(unitbin) x1out(:,:)
     write(unitbin) x2out(:,:)
     write(unitbin) thout(:,:)
     write(unitbin) radout(:,:,:,:)
     write(unitbin) hydout(:,:,:,:)
     close(unitbin)
  else
     !------------------------------------------------------------
     ! Moment computation (parallel over space; angle loop is serial)
     !------------------------------------------------------------
!$omp parallel do default(shared) private(i,j,k,m) collapse(3) schedule(static)
     do k=ks,ke
     do j=js,je
     do i=is,ie
        Erad(i,j,k) = 0.d0
        Frad(xdir,i,j,k) = 0.d0
        Frad(ydir,i,j,k) = 0.d0
        do m=1,mang
           Erad(i,j,k)      =      Erad(i,j,k) +        fdis(m,i,j,k) * dth(m)
           Frad(xdir,i,j,k) = Frad(xdir,i,j,k) + mux(m)*fdis(m,i,j,k) * dth(m)
           Frad(ydir,i,j,k) = Frad(ydir,i,j,k) + muy(m)*fdis(m,i,j,k) * dth(m)
        enddo
     enddo
     enddo
     enddo
!$omp end parallel do
         
     write(filename,'(a4,i5.5,a4)')"snap",nout,".dat"
     filename = trim(dirname)//filename
     open(newunit=unitasc,file=filename,status='replace',form='formatted',access="stream",action="write") 
     write(unitasc,"(a1,(1x,(A)),(1x,1PE15.4))") "#","time=",time
     write(unitasc,"(a1,(1x,(A)),(1x,i0))") "#","nx=", izones
     write(unitasc,"(a1,(1x,(A)),(1x,i0))") "#","ny=", jzones
     write(unitasc,"(a1,(A))") "#"," x y E Fx Fy"
     k=ks
     do j=js,je
     do i=is,ie
        write(unitasc,"(1x,5(1x,E15.6e3))") x1b(i),x2b(j),Erad(i,j,k),Frad(xdir,i,j,k),Frad(ydir,i,j,k)
     enddo
        write(unitasc,*)
     enddo
     close(unitasc)
         
  endif
  print *, "output:",nout,time
  
  nout=nout+1
  tout=time

  return
end subroutine Output

program main
  use gridmod
  use angmod
  use statemod
  use fluxmod
  use timemod
  use stepmod
  implicit none
  character(20),parameter:: dirname="sn/"
  logical :: flag_binary = .false.
  logical,parameter:: force_on = .true.,force_off = .false.
  write(6,*) "setup grids and fields"
  call GenerateGrid
  call GenerateAngleGrid()
  call GenerateProblem
  call Output(force_on,flag_binary,dirname)
  write(6,*) "entering main loop"
  ! main loop
  write(6,*)"step","time","dt"
  mloop: do ntime=1,ntimemax
     call TimestepControl
     if(mod(ntime,100) .eq. 0 ) then
        print *, "step=",ntime," time=",time," dt=",dt
        call flush(6)
     endif
     call RadBoundaryCondition
     call RadFlux1
     call RadFlux2
     call UpdateRad
     call UpdateRadWithSource
     time=time+dt
     call Output(force_off,flag_binary,dirname)
     if(time > timemax) exit mloop
  enddo mloop
  call Output(force_on,flag_binary,dirname)

  write(6,*) "program has been finished"
  
end program main



subroutine makedirs(outdir)
  implicit none
  character(len=*), intent(in) :: outdir
  character(len=256) command
  write(command, *) 'if [ ! -d ', trim(outdir), ' ]; then mkdir -p ', trim(outdir), '; fi'
  write(*, *) trim(command)
  call system(command)
end subroutine makedirs
