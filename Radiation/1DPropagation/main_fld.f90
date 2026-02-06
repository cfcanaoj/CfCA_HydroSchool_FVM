module units
  implicit none
  real(8),parameter::cl=3.0d10 ! cm/s
end module units

module modelpara
  implicit none
  real(8),parameter:: rho0=0.25d0! [g/cm^3]
  real(8),parameter:: kap0=0.04d0! [cm^2/g]
  real(8),parameter:: erad0=1.0d10! [erg/cm^3]
  real(8),parameter:: erad1=1.0d20! [erg/cm^3]
end module modelpara

module commons
  use units
  implicit none

  integer::ntime
  integer,parameter::ntimemax=1000000
  real(8)::time,dt
  real(8),parameter:: cfl=0.1d0
  data time / 0.0d0 /
  real(8),parameter:: timemax=1.0d-10
  real(8),parameter:: dtout=timemax/100
  
  integer,parameter::izones=100
  integer,parameter::jzones=1
  integer,parameter::mgn=2
  integer,parameter::in=izones+2*mgn+1 &
       &            ,jn=1 &
       &            ,kn=1
  integer,parameter::is=mgn+1 &
       &            ,js=1 &
       &            ,ks=1
  integer,parameter::ie=izones+mgn &
       &            ,je=1 &
       &            ,ke=1
  
  real(8),parameter:: x1min=  0.0d0,x1max=1.0d0
  real(8),dimension(in)::x1a,x1b
  real(8),dimension(jn)::x2a,x2b
  real(8),dimension(kn)::x3a,x3b
  
  real(8),dimension(in,jn,kn):: d
  real(8),dimension(in,jn,kn):: kappa
  real(8),dimension(in,jn,kn):: Elte
      
  real(8),dimension(    in,jn,kn)::Erad
  integer,parameter::xdir=1,ydir=2,zdir=3
  real(8),dimension(xdir:zdir,in,jn,kn)::Frad
  real(8),dimension(6,in,jn,kn)::Prad
      

  real(8),dimension(    in,jn,kn)::LFLD1
  integer,parameter:: nxx=1,nyy=2,nzz=3,nxy=4,nyz=5,nzx=6

  real(8)::denup,dendn
      
end module commons
     
module eosmod
  implicit none
  ! adiabatic
  real(8),parameter::gam=1.4d0 !! adiabatic index
  ! isothermal
  !      real(8)::csiso  !! isothemal sound speed
end module eosmod

module fluxmod
  use commons, only : in,jn,kn
  implicit none
  integer,parameter::nerd=1,nfr1=2,nfr2=3,nfr3=4
  integer,parameter::nrad=4
  real(8),dimension(nrad,in,jn,kn):: radsvc

  integer,parameter::muerd=1,mufru=2,mufrv=3,mufrw=4 &
       &            ,mferd=5,mffru=6,mffrv=7,mffrw=8 &
       &            ,mcsp=9,mcsm=10

  integer,parameter::mradflx=4,mradadd=2
  
  integer,parameter:: merd=1,mfr1=2,mfr2=3,mfr3=4 &
       &            , mfru=mufru,mfrv=mufrv,mfrw=mufrw
  real(8),dimension(mradflx,in,jn,kn):: radnflux1,radnflux2,radnflux3

end module fluxmod
!===============================================================
!  Chebyshev Super-Time-Stepping (STS) substep generator
!  - Input : dt_exp  (the usual explicit stable timestep for diffusion)
!            M       (number of STS substeps, e.g. 10-30)
!            nu      (damping parameter, e.g. 1d-3 to 1d-2; must be > 0)
!  - Output: dt_sub(1:M)  (substep timesteps)
!            dt_tot       (sum of dt_sub; the effective large step)
!
!  Formula (Alexiades et al. type STS):
!     dt_j = dt_exp / [ (1+nu) - (1-nu) * cos( (2j-1) * pi / (2M) ) ]
!===============================================================
module stsmod
  implicit none
  integer, parameter :: dp = kind(1.0d0)
contains

  subroutine MakeSTS(dt_exp, M, nu, dt_sub, dt_tot)
    implicit none
    real(dp), intent(in)  :: dt_exp, nu
    integer,  intent(in)  :: M
    real(dp), intent(out) :: dt_sub(M), dt_tot

    integer :: j
    real(dp) :: pi, theta, denom, nu_eff

    if (M <= 0) then
       dt_tot = 0.0_dp
       return
    end if

    if (dt_exp <= 0.0_dp) then
       dt_sub(:) = 0.0_dp
       dt_tot    = 0.0_dp
       return
    end if

    ! Safety: nu must be > 0 (otherwise some denom can get too small)
    nu_eff = max(nu, 1.0d-12)

    pi = 4.0_dp * atan(1.0_dp)
    dt_tot = 0.0_dp

    do j = 1, M
      theta = (2.0_dp*j - 1.0_dp) * pi / (2.0_dp*M)
      denom = (1.0_dp + nu_eff) - (1.0_dp - nu_eff) * cos(theta)
      dt_sub(j) = dt_exp / denom
      dt_tot    = dt_tot + dt_sub(j)
    end do
  end subroutine MakeSTS

end module stsmod

program main
  use commons
  use stsmod
  implicit none
  character(20),parameter:: dirname="fld/"
  logical :: flag_binary = .false.
  logical,parameter:: force_on = .true.,force_off = .false.
  ! super time step
  logical,parameter:: sts_on = .true.
  integer,parameter:: Msts=12
  real(8),parameter:: nusts=1.0d-3
  real(8),dimension(Msts):: dt_sub
  real(8):: dt_exp,dt_eff
  integer:: nsts
  write(6,*) "setup grids and fields"
  call GenerateGrid
  call GenerateProblem
  call Output(force_on,flag_binary,dirname)
  
  write(6,*) "entering main loop"
  ! main loop
  write(6,*)"step","time","dt"
  mloop: do ntime=1,ntimemax
     if(sts_on) then
        call TimestepControl(dt_exp)
        if(mod(ntime,100) .eq. 0 ) write(6,*)ntime,time,dt,Erad(is,js,ks)
        call MakeSTS(dt_exp, Msts, nusts, dt_sub, dt_eff)
        do nsts=1,Msts
           dt = dt_sub(nsts)
           call RadBoundaryCondition
           call StateRad
           call RadFlux1
           call UpdateRadAdvection
           call UpdateRadSource
           time=time+dt
        enddo
        call Output(force_off,flag_binary,dirname)
        if(time > timemax) exit mloop
     else
        call TimestepControl(dt)
        if(mod(ntime,100) .eq. 0 ) write(6,*)ntime,time,dt,Erad(is,js,ks)
        call RadBoundaryCondition
        call StateRad
        call RadFlux1
        call UpdateRadAdvection
        call UpdateRadSource
        time=time+dt
        call Output(force_off,flag_binary,dirname)
        if(time > timemax) exit mloop
     endif 
  enddo mloop
  call Output(force_on,flag_binary,dirname)
     
  write(6,*) "program has been finished"
  !      call Debug
end program main

subroutine GenerateGrid
  use commons
  implicit none
  real(8)::dx,dy
  integer::i,j,k
  dx=(x1max-x1min)/izones
  do i=1,in
     x1a(i) = dx*(i-(mgn+1))+x1min
  enddo
  do i=1,in-1
     x1b(i) = 0.5d0*(x1a(i+1)+x1a(i))
  enddo

  return
end subroutine GenerateGrid

subroutine GenerateProblem
  use commons
  use eosmod
  use modelpara
  implicit none
  integer::i,j,k
  
  do k=1,kn
  do j=1,jn
  do i=1,in

     d(i,j,k) = rho0   ! [g cm^-3]
     kappa(i,j,k) = kap0   ! [cm^2/g]
     Elte(i,j,k) = erad0 ! [erg/cm^3]

  enddo
  enddo
  enddo
      
  
  Erad(:,:,:) = erad0

  do k=ks,ke
  do j=js,je
  do i=is,ie

     Erad(  i,j,k) = erad0
     Frad(1,i,j,k) = 0.0d0
     Frad(2,i,j,k) = 0.0d0
     Frad(3,i,j,k) = 0.0d0

  enddo
  enddo
  enddo
      
  call RadBoundaryCondition
  
  return
end subroutine GenerateProblem


subroutine RadBoundaryCondition
  use modelpara
  use commons
  implicit none
  integer::i,j,k

  do k=ks,ke
  do j=js,je
  do i=1,mgn
     Erad(    is-i,j,k) = erad1
     Frad(1,  is-i,j,k) = Erad(    is-i,j,k)/sqrt(3.0d0)
     Frad(2:3,  is-i,j,k) = 0.0d0
  enddo
  enddo
  enddo

  do k=ks,ke
  do j=js,je
  do i=1,mgn
     Erad(    ie+i,j,k) = Erad(ie,j,k) 
     Frad(:,  ie+i,j,k) = Frad(:,ie,j,k)
  enddo
  enddo
  enddo

  return
end subroutine RadBoundaryCondition

subroutine TimestepControl(dtcfl)
  use commons
  implicit none
  real(8),intent(inout)::dtcfl
  real(8)::dtsrc
  real(8)::dtl1,dtdif1
  real(8)::dtl2,dtdif2
  real(8)::dtl3,dtdif3
  real(8)::dtlocal
  real(8)::dtmin
  integer::i,j,k
  dtmin=1.0d90
  do k=ks,ke
  do j=js,je
  do i=is,ie
! dt < 1/2 dx^2/D [s] , D = c/(kappa*rho) [cm^2/s],  dt < 1/2 dx^2/c*(kappa rho) 
     dtdif1 = 0.5d0*(x1a(i+1)-x1a(i))**2*d(i,j,k)*kappa(i,j,k)/cl 
     dtl1  = (x1a(i+1)-x1a(i))/(cl)
     !         dtl2 =(x2a(j+1)-x2a(j))/(cl)
     !         dtl3 =(x1a(i+1)-x1a(i))/(abs(v1(i,j,k)) +cs(i,j,k))
     !         dtlocal = min (dtl1,dtl2)
     dtlocal = min(dtl1,dtdif1)
     !         dtlocal = dtl1

     if(dtlocal .lt. dtmin) dtmin = dtlocal
  enddo
  enddo
  enddo

  dtcfl = cfl * dtmin
  !      write(6,*)"dt",dt
  return
end subroutine TimestepControl

subroutine StateRad
  use commons
  use fluxmod
  use eosmod
  implicit none
  integer::i,j,k
  real(8)::fnl,fsq,chi
  do k=ks,ke
  do j=js,je
  do i=1,in-1
     radsvc(nerd,i,j,k) = Erad( i,j,k) 
     radsvc(nfr1,i,j,k) = Frad(1,i,j,k) 
     radsvc(nfr2,i,j,k) = Frad(2,i,j,k) 
     radsvc(nfr3,i,j,k) = Frad(3,i,j,k) 

!         if(i .eq. is) print *, "StateRad",radsvc(nerd,i,j,k),radsvc(nerd,i-1,j,k)
  enddo
  enddo
  enddo

  return
end subroutine StateRad

subroutine RadFlux1
  use commons
  use fluxmod
  use units
  implicit none
  integer::i,j,k
  real(8):: Rfc,lambda
  real(8),parameter::tiny=1.0d-30


  do k=ks,ke
  do j=js,je
  do i=is,ie+1
     Rfc = abs((Erad(i,j,k)-Erad(i-1,j,k))/(x1b(i)-x1b(i-1))) &
    &  * 2.0d0/(Erad(i,j,k)+Erad(i-1,j,k)) &
    &  * 2.0d0/(d(i,j,k)*kappa(i,j,k)+d(i-1,j,k)*kappa(i-1,j,k))
     lambda = (2.0+Rfc)/(6+3*Rfc+Rfc**2)
     LFLD1(i,j,k) = lambda
     radnflux1(merd,i,j,k)= -(Erad(i,j,k)-Erad(i-1,j,k))/(x1b(i)-x1b(i-1)) &
          & /(d(i,j,k)*kappa(i,j,k)+d(i-1,j,k)*kappa(i-1,j,k))*2.0d0*lambda
     if( radnflux1(merd,i,j,k) .ge. 0)then
        radnflux1(merd,i,j,k) = min(radnflux1(merd,i,j,k), Erad(i-1,j,k))
     else
        radnflux1(merd,i,j,k) = max(radnflux1(merd,i,j,k),-Erad(i  ,j,k))
     endif

  enddo
  enddo
  enddo

  return
end subroutine Radflux1


subroutine UpdateRadSource
  use commons
  use fluxmod
  implicit none
  integer::i,j,k
  real(8)::pi
  real(8)::alpha
      
  pi = acos(-1.0d0)

  do k=ks,ke
  do j=js,je
  do i=is,ie
     alpha = d(i,j,k)*kappa(i,j,k)*cl*dt
     Erad(i,j,k) = (Erad(i,j,k) + alpha*Elte(i,j,k))/(1+alpha)
  enddo
  enddo
  enddo

  return
end subroutine UpdateRadSource

subroutine UpdateRadAdvection
  use commons
  use fluxmod
  implicit none
  integer::i,j,k

  do k=ks,ke
  do j=js,je
  do i=is,ie
         
         Erad(i,j,k) = Erad(i,j,k)                    &
     & +dt*cl*(                                       &
     & +(- radnflux1(merd,i+1,j,k)                    &
     &   + radnflux1(merd,i  ,j,k))/(x1a(i+1)-x1a(i)) &
     &      )

         Frad(1,i,j,k) = 0.5d0*(radnflux1(merd,i+1,j,k)+radnflux1(merd,i,j,k))
         Frad(2:3,i,j,k) = 0.0
  enddo
  enddo
  enddo

  return
end subroutine UpdateRadAdvection

subroutine Output(flag_force,flag_binary,dirname)
  use commons
  implicit none
  logical,intent(in):: flag_force
  logical,intent(in):: flag_binary
  character(len=*),intent(in):: dirname
  integer::i,j,k
  character(40)::filename
  real(8),save::tout
  data tout / 0.0d0 /
  integer,save::nout
  data nout / 1 /
  integer:: unitasc, unitbin
  integer,parameter:: gs=1
  integer,parameter:: nvar=2
  real(8)::x1out(is-gs:ie+gs,2)
  real(8)::radout(is-gs:ie+gs,js,ks,nvar)

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

     x1out(is-gs:ie+gs,1) = x1b(is-gs:ie+gs)
     x1out(is-gs:ie+gs,2) = x1a(is-gs:ie+gs)

     radout(is-gs:ie+gs,js,ks,1) = Erad(  is-gs:ie+gs,js,ks)
     radout(is-gs:ie+gs,js,ks,2) = Frad(1,is-gs:ie+gs,js,ks)
  
     write(filename,'(a4,i5.5,a4)')"snap",nout,".bin"
     filename = trim(dirname)//filename
     open(newunit=unitbin,file=filename,status='replace',form='unformatted',access="stream",action="write")
     write(unitbin) time
     write(unitbin) izones+2*gs
     write(unitbin) 2
     write(unitbin) x1out(:,:)
     write(unitbin) radout(:,:,:,:)
     close(unitbin)
  else
     write(filename,'(a4,i5.5,a4)')"snap",nout,".dat"
     filename = trim(dirname)//filename
     open(newunit=unitasc,file=filename,status='replace',form='formatted',access="stream",action="write") 
     write(unitasc,"(a1,(1x,(A)),(1x,1PE15.4))") "#","time=",time
     write(unitasc,"(a1,(1x,(A)),(1x,i0))") "#","nx=", izones+2*gs
     write(unitasc,"(a1,(A))") "#"," x E Fx"
     k=ks
     j=js
     do i=is-gs,ie+gs
        write(unitasc,*) x1b(i),Erad(i,j,k),Frad(xdir,i,j,k)
     enddo
     close(unitasc)
  endif
  write(6,*) "output:",nout,ntime,time,Erad(is,js,ks)/1e20

  nout=nout+1
  tout=time
      
  return
end subroutine Output

subroutine makedirs(outdir)
  implicit none
  character(len=*), intent(in) :: outdir
  character(len=256) command
  write(command, *) 'if [ ! -d ', trim(outdir), ' ]; then mkdir -p ', trim(outdir), '; fi'
  write(*, *) trim(command)
  call system(command)
end subroutine makedirs

subroutine Debug
  use commons
  implicit none
  integer::i,j,k

  do k=ks,ke
  do j=js,je
  do i=is,ie
     write(6,*) x1b(i),Erad(i,j,k),Frad(1,i,j,k),LFLD1(i,j,k)
  enddo
  enddo
  enddo

  return
end subroutine Debug


