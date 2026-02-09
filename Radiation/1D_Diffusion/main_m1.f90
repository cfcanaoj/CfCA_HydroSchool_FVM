module units
  implicit none
  real(8),parameter::cl=3.0d10 ! cm/s
end module units

module modelpara
  implicit none
  real(8),parameter::  rho0 = 1.0d0! [g/cm^3]
  real(8),parameter::  kap0 = 1.0d3! [cm^2/g]
  real(8),parameter:: erad0 = 1.0d0 ! [erg/cm^3]
  real(8),parameter:: erad1 = 1.0d20! [erg/cm^3]
end module modelpara

module commons
  use units
  implicit none

  integer::ntime
  integer,parameter::ntimemax=100000
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
  real(8),dimension(in,jn,kn):: kappaA
  real(8),dimension(in,jn,kn):: kappaS
  real(8),dimension(in,jn,kn):: kappa
  real(8),dimension(in,jn,kn):: Elte

  real(8),dimension(    in,jn,kn)::Erad
  integer,parameter:: xdir=1,ydir=2,zdir=3
  real(8),dimension(xdir:zdir,  in,jn,kn)::Frad
  integer,parameter:: nxx=1,nyy=2,nzz=3,nxy=4,nyz=5,nzx=6
  real(8),dimension(nxx:nzx,in,jn,kn)::Prad

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
  logical,parameter::flagEdd=.false.
  real(8),parameter::fluxfactormax=0.999d0

end module fluxmod
module closure
  implicit none
contains
  real(8) function chi(f)
    implicit none
    real(8),intent(in)::f
    !          chi = (3+4*f**2)/(5+2*sqrt(4-3*f**2))
    chi  = ( 5.0d0 + 6.0d0*f**2 -2.0d0*f**3 + 6.0d0*f**4) /15.0d0
    
  end function chi
  real(8) function chiovf(f)
    implicit none
    real(8),intent(in)::f
    !          chiovf = (6*f*(3+4*f**2))/(sqrt(4-3*f**2)*(5+2*(4-3*f**2))**2) &
    !             & +  (8*f)/(5+2*sqrt(4-3*f**2))
    
    chiovf  = (        12.0d0*f    -6.0d0*f**2 +24.0d0*f**3) /15.0d0
    
  end function chiovf
  real(8) function finvchiovf(f)
    implicit none
    real(8),intent(in)::f
    !          chiovf = (6*f*(3+4*f**2))/(sqrt(4-3*f**2)*(5+2*(4-3*f**2))**2) &
    !             & +  (8*f)/(5+2*sqrt(4-3*f**2))
    
    finvchiovf  = (        12.0d0    -6.0d0*f +24.0d0*f**2) /15.0d0
    
  end function finvchiovf

  subroutine lambdaPlusMinus(f,mu,lambdaPlus,lambdaMins)
    implicit none
    real(8),intent(in)::f,mu
    real(8),intent(out)::lambdaPlus,lambdaMins
    real(8)::lambdaPlusOne,lambdaMinsOne,lambdaPlusZer,lambdaMinsZer
    lambdaPlusOne =  chiovf(f)/2 + 1.0/2*sqrt(chiovf(f)**2+4*(chi(f) - f* chiovf(f)))
    lambdaMinsOne =  chiovf(f)/2 - 1.0/2*sqrt(chiovf(f)**2+4*(chi(f) - f* chiovf(f)))
    lambdaPlusZer =  0.5d0*sqrt(2.0d0*(1-chi(f)) + finvchiovf(f)*(1+2*f**2-3*chi(f)) )
    lambdaMinsZer = -0.5d0*sqrt(2.0d0*(1-chi(f)) + finvchiovf(f)*(1+2*f**2-3*chi(f)) )

    if(mu > 0)then
       lambdaPlus = lambdaPlusZer + abs(mu)*(lambdaPlusOne-lambdaPlusZer)
       lambdaMins = lambdaMinsZer + abs(mu)*(lambdaMinsOne-lambdaMinsZer)
    else
       lambdaPlus = -(lambdaMinsZer + abs(mu)*(lambdaMinsOne-lambdaMinsZer))
       lambdaMins = -(lambdaPlusZer + abs(mu)*(lambdaPlusOne-lambdaPlusZer))
    endif

  end subroutine lambdaPlusMinus
end module closure

program main
  use commons
  implicit none
  character(20),parameter:: dirname="m1/"
  logical :: flag_binary = .false.
  logical,parameter:: force_on = .true.,force_off = .false.
  
  write(6,*) "setup grids and fields"
  call GenerateGrid
  call GenerateProblem
  call Output(force_on,flag_binary,dirname)
  write(6,*) "entering main loop"
  ! main loop
  write(6,*)"step","time","dt"
  mloop: do ntime=1,ntimemax
     call TimestepControl
     if(mod(ntime,100) .eq. 0 ) write(6,*)ntime,time,dt,Erad(is,js,ks)
     call RadBoundaryCondition
     call StateRad
     call RadFlux1
     call UpdateRadAdvection
     call UpdateRadSource
     time=time+dt
     call Output(force_off,flag_binary,dirname)
     if(time > timemax) exit mloop
  enddo mloop
      
  call Output(force_on,flag_binary,dirname)
  
  write(6,*) "program has been finished"
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
      
      do k=ks,ke
      do j=js,je
      do i=1,in
          d(i,j,k) = rho0   ! [g cm^-3]
     kappaA(i,j,k) = 0.0d0   ! [cm^2/g]
     kappaS(i,j,k) = kap0   ! [cm^2/g]
     kappa(i,j,k) = kappaA(i,j,k)+kappaS(i,j,k)
       Elte(i,j,k) = erad0 ! [erg/cm^3]

      enddo
      enddo
      enddo
      
      do k=ks,ke
      do j=js,je
      do i=1,in

       Erad(  i,j,k) = erad1*exp(-(x1b(i)/0.1)**2/2)
       Frad(xdir,i,j,k) = 0.0d0
       Frad(ydir,i,j,k) = 0.0d0
       Frad(zdir,i,j,k) = 0.0d0

      enddo
      enddo
      enddo
      
   call RadBoundaryCondition
 !  write(6,*) "Gene",Erad(  is-1,js,ks),Erad(  is,js,ks)

   return
end subroutine GenerateProblem


subroutine RadBoundaryCondition
  use commons
  use modelpara
  use fluxmod
  implicit none
  integer::i,j,k

  do k=ks,ke
  do j=js,je
  do i=1,mgn
  ! Reflection      
     Erad(     is-i,j,k) = Erad(     is+i-1,j,k)
     Frad(xdir,is-i,j,k) = Frad(xdir,is+i-1,j,k) 
  enddo
  enddo
  enddo 

  do k=ks,ke
  do j=js,je
  do i=1,mgn
  ! outflow      
     Erad(            ie+i,j,k) = Erad(          ie,j,k) 
     Frad(xdir:zdir,  ie+i,j,k) = Frad(xdir:zdir,ie,j,k)
  enddo
  enddo
  enddo

  return
end subroutine RadBoundaryCondition

      subroutine TimestepControl
      use commons
      implicit none
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
         dtl1  = (x1a(i+1)-x1a(i))/(cl)
!         dtl2 =(x2a(j+1)-x2a(j))/(cl)
!         dtl3 =(x1a(i+1)-x1a(i))/(abs(v1(i,j,k)) +cs(i,j,k))
!         dtlocal = min (dtl1,dtl2)
!         if(i.eq.is)write(6,*) "dt",dtsrc,dtl1
         dtlocal = min(dtl1,dtsrc*10.0d0)
         if(dtlocal .lt. dtmin) dtmin = dtlocal
      enddo
      enddo
      enddo

      dt = cfl * dtmin
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
         radsvc(nerd,i,j,k) = Erad(  i,j,k) 
         radsvc(nfr1,i,j,k) = Frad(xdir,i,j,k)/Erad(  i,j,k)  
         radsvc(nfr2,i,j,k) = Frad(ydir,i,j,k)/Erad(  i,j,k) 
         radsvc(nfr3,i,j,k) = Frad(zdir,i,j,k)/Erad(  i,j,k)

!         if(i .eq. is) print *, "StateRad",radsvc(nerd,i,j,k),radsvc(nerd,i-1,j,k)
      enddo
      enddo
      enddo

      return
      end subroutine StateRad

      subroutine minmod(a,b,d)
      use fluxmod, only : nrad
      implicit none
      real(8),dimension(nrad),intent(in)::a,b
      real(8),dimension(nrad),intent(out)::d
      integer:: n

      do n=1,nrad
         d(n) = sign(1.0d0,a(n))*max(0.0d0,min(abs(a(n))                &
     &                                        ,sign(1.0d0,a(n))*b(n)))
      enddo

      return
      end subroutine minmod


      subroutine vanLeer(dvp,dvm,dv)
      use fluxmod, only : nrad
      implicit none
      real(8),dimension(nrad),intent(in)::dvp,dvm
      real(8),dimension(nrad),intent(out)::dv
      integer:: n

      do n=1,nrad
         if(dvp(n)*dvm(n) .gt. 0.0d0)then
            dv(n) =2.0d0*dvp(n)*dvm(n)/(dvp(n)+dvm(n))
         else
            dv(n) = 0.0d0
         endif

      enddo

      return
      end subroutine vanLeer

      subroutine MClimiter(a,b,c,d)
      use fluxmod, only : nrad
      implicit none
      real(8),dimension(nrad),intent(in)::a,b,c
      real(8),dimension(nrad),intent(out)::d
      integer:: n

      do n=1,nrad
         d(n) = sign(1.0d0,a(n))*max(0.0d0,min(abs(a(n))         &
     &                                  ,sign(1.0d0,a(n))*b(n)   &
     &                                  ,sign(1.0d0,a(n))*c(n))) 
      enddo

      return
      end subroutine MClimiter

      subroutine HLLERAD(leftst,rigtst,nflux)
      use fluxmod
      implicit none
      real(8),dimension(2*mradflx+mradadd),intent(in)::leftst,rigtst
      real(8),dimension(mradflx),intent(out)::nflux
      real(8),dimension(mradflx)::ul,ur,fl,fr
      real(8)::cslp,cslm,csrp,csrm
      real(8):: vl, vr
      real(8):: sl, sr

      ul(1:mradflx) = leftst(1:mradflx)
      fl(1:mradflx) = leftst(mradflx+1:2*mradflx)
      ur(1:mradflx) = rigtst(1:mradflx)
      fr(1:mradflx) = rigtst(mradflx+1:2*mradflx)
      cslp = leftst(mcsp)
      cslm = leftst(mcsm)

      csrp = rigtst(mcsp)
      csrm = rigtst(mcsm)

      sl = min(cslm,csrm,0.0)
      sr = max(cslp,csrp,0.0d0)
      
      nflux(:) = (sr*fl(:)-sl*fr(:) +sl*sr*(ur(:)-ul(:)))/(sr-sl)

      return
      end subroutine HLLERAD


subroutine RadFlux1
  use commons, only: is,ie,js,je,ks,ke
  use fluxmod
  use units
  use closure
  implicit none
  integer::i,j,k
  real(8),dimension(nrad):: dsv,dsvm,dsvp
  real(8),dimension(nrad):: Pleftc1,Pleftc2,Plefte
  real(8),dimension(nrad):: Prigtc1,Prigtc2,Prigte
  real(8),dimension(2*mradflx+mradadd):: leftco,rigtco
  real(8),dimension(mradflx):: nflux
  real(8):: fsq,ffc,pxx,pyy,pzz,pxy,pyz,pzx
  real(8):: ff1,ff2,ff3
  real(8):: chil,lp,lm
  real(8),parameter::tiny=1.0d-30

  do k=ks,ke
  do j=js,je
  do i=is,ie+1

     Pleftc1(:) = radsvc(:,i-2,j,k)
     Pleftc2(:) = radsvc(:,i-1,j,k)
     Prigtc1(:) = radsvc(:,i  ,j,k)
     Prigtc2(:) = radsvc(:,i+1,j,k)

! | Pleftc1   | Pleftc2 =>| Prigtc1   | Prigtc2   |
!                     You are here (x1-interface at i)

!====================
! Left state (from i-1 cell)
!====================
     dsvm(:) = Pleftc2(:) - Pleftc1(:)
     dsvp(:) = Prigtc1(:) - Pleftc2(:)
     call vanLeer(dsvp,dsvm,dsv)
     Plefte(:) = Pleftc2(:) + 0.5d0*dsv(:)

!====================
! Right state (from i cell)
!====================
     dsvm(:) = Prigtc1(:) - Pleftc2(:)
     dsvp(:) = Prigtc2(:) - Prigtc1(:)
     call vanLeer(dsvp,dsvm,dsv)
     Prigte(:) = Prigtc1(:) - 0.5d0*dsv(:)

!====================
! Convert primitives -> conservative/flux + characteristic speeds (M1 closure)
!====================
! ----- left
     fsq = (Plefte(nfr1)**2 +Plefte(nfr2)**2 +Plefte(nfr3)**2) ! (F/E/c)^2
     ffc = min(1.0d0,max(tiny,sqrt(fsq)))
     ff1 = min(1.0d0,max(-1.0d0,Plefte(nfr1)/ffc))
     ff2 = min(1.0d0,max(-1.0d0,Plefte(nfr2)/ffc))
     ff3 = min(1.0d0,max(-1.0d0,Plefte(nfr3)/ffc))

     chil = chi(ffc)
     if(flagEdd) chil = 1.0d0/3.0d0
     call lambdaPlusMinus(ffc,ff1,lp,lm)

     leftco(mcsp) = min(lp+tiny, 1.0d0)
     leftco(mcsm) = max(lm-tiny,-1.0d0)

     pxx = 0.5d0*(1.0d0-chil) + ff1**2 *0.5d0*(3.0d0*chil-1.0d0)
     pyy = 0.5d0*(1.0d0-chil) + ff2**2 *0.5d0*(3.0d0*chil-1.0d0)
     pzz = 0.5d0*(1.0d0-chil) + ff3**2 *0.5d0*(3.0d0*chil-1.0d0)
     pxy =                      ff1*ff2*0.5d0*(3.0d0*chil-1.0d0)
     pyz =                      ff2*ff3*0.5d0*(3.0d0*chil-1.0d0)
     pzx =                      ff3*ff1*0.5d0*(3.0d0*chil-1.0d0)

     leftco(muerd) =              Plefte(nerd) ! E
     leftco(mufru) = Plefte(nfr1)*Plefte(nerd) ! F_u
     leftco(mufrv) = Plefte(nfr2)*Plefte(nerd) ! F_v
     leftco(mufrw) = Plefte(nfr3)*Plefte(nerd) ! F_w

     leftco(mferd) = Plefte(nfr1)*Plefte(nerd) ! F
     leftco(mffru) =          pxx*Plefte(nerd) ! P_u
     leftco(mffrv) =          pxy*Plefte(nerd) ! P_v
     leftco(mffrw) =          pzx*Plefte(nerd) ! P_w

! ----- right
     fsq = (Prigte(nfr1)**2 +Prigte(nfr2)**2 +Prigte(nfr3)**2) ! (F/E/c)^2
     ffc = min(1.0d0,max(tiny,sqrt(fsq)))
     ff1 = min(1.0d0,max(-1.0d0,Prigte(nfr1)/ffc))
     ff2 = min(1.0d0,max(-1.0d0,Prigte(nfr2)/ffc))
     ff3 = min(1.0d0,max(-1.0d0,Prigte(nfr3)/ffc))

     chil = chi(ffc)
     if(flagEdd) chil = 1.0d0/3.0d0
     call lambdaPlusMinus(ffc,ff1,lp,lm)

     rigtco(mcsp) = min(lp+tiny, 1.0d0)
     rigtco(mcsm) = max(lm-tiny,-1.0d0)

     pxx = 0.5d0*(1.0d0-chil) + ff1**2 *0.5d0*(3.0d0*chil-1.0d0)
     pyy = 0.5d0*(1.0d0-chil) + ff2**2 *0.5d0*(3.0d0*chil-1.0d0)
     pzz = 0.5d0*(1.0d0-chil) + ff3**2 *0.5d0*(3.0d0*chil-1.0d0)
     pxy =                      ff1*ff2*0.5d0*(3.0d0*chil-1.0d0)
     pyz =                      ff2*ff3*0.5d0*(3.0d0*chil-1.0d0)
     pzx =                      ff3*ff1*0.5d0*(3.0d0*chil-1.0d0)

     rigtco(muerd) =              Prigte(nerd) ! E
     rigtco(mufru) = Prigte(nfr1)*Prigte(nerd) ! F_u
     rigtco(mufrv) = Prigte(nfr2)*Prigte(nerd) ! F_v
     rigtco(mufrw) = Prigte(nfr3)*Prigte(nerd) ! F_w

     rigtco(mferd) = Prigte(nfr1)*Prigte(nerd) ! F
     rigtco(mffru) =          pxx*Prigte(nerd) ! P_u
     rigtco(mffrv) =          pxy*Prigte(nerd) ! P_v
     rigtco(mffrw) =          pzx*Prigte(nerd) ! P_w

!====================
! Riemann solver
!====================
     call HLLERAD(leftco,rigtco,nflux)
     radnflux1(merd,i,j,k) = nflux(merd)
     radnflux1(mfr1,i,j,k) = nflux(mfru)
     radnflux1(mfr2,i,j,k) = nflux(mfrv)
     radnflux1(mfr3,i,j,k) = nflux(mfrw)

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
  real(8)::fnl,flm
  
  pi = acos(-1.0d0)
  
  do k=ks,ke
  do j=js,je
  do i=is,ie
     alpha = d(i,j,k)*kappaA(i,j,k)*cl*dt
     Erad(i,j,k) = (Erad(i,j,k) + alpha*Elte(i,j,k))/(1+alpha)
     alpha = d(i,j,k)*kappa(i,j,k)*cl*dt
     Frad(xdir,i,j,k) = Frad(xdir,i,j,k)/(1.0d0+alpha)
     Frad(ydir,i,j,k) = Frad(ydir,i,j,k)/(1.0d0+alpha)
     Frad(zdir,i,j,k) = Frad(zdir,i,j,k)/(1.0d0+alpha) 
     
     Erad(i,j,k) =  max(Erad(i,j,k),Elte(i,j,k))

     fnl = sqrt(Frad(xdir,i,j,k)**2 +Frad(ydir,i,j,k)**2 +Frad(zdir,i,j,k)**2)
     flm =  fluxfactormax*Erad(i,j,k)
     if(fnl .gt. flm )then
        Frad(xdir,i,j,k) =  Frad(xdir,i,j,k)*flm/fnl
        Frad(ydir,i,j,k) =  Frad(ydir,i,j,k)*flm/fnl
        Frad(zdir,i,j,k) =  Frad(zdir,i,j,k)*flm/fnl
     endif
  enddo
  enddo
  enddo

  return
end subroutine UpdateRadSource

subroutine UpdateRadAdvection
  use commons
  use fluxmod
  implicit none
  real(8)::fnl,flm
  integer::i,j,k

  do k=ks,ke
  do j=js,je
  do i=is,ie
         
         Erad(i,j,k) = Erad(i,j,k)                    &
     & +dt*cl*(                                          &
     & +(- radnflux1(merd,i+1,j,k)                    &
     &   + radnflux1(merd,i  ,j,k))/(x1a(i+1)-x1a(i)) &
     &      )

         Frad(xdir,i,j,k) = Frad(xdir,i,j,k)          &
     & +dt*cl*(                                          &
     & +(- radnflux1(mfr1,i+1,j,k)                    &
     &   + radnflux1(mfr1,i  ,j,k))/(x1a(i+1)-x1a(i)) &
     &      )

         Frad(ydir,i,j,k) = Frad(ydir,i,j,k)          &
     & +dt*cl*(                                          &
     & +(- radnflux1(mfr2,i+1,j,k)                    &
     &   + radnflux1(mfr2,i  ,j,k))/(x1a(i+1)-x1a(i)) &
     &      )

         Frad(zdir,i,j,k) = Frad(zdir,i,j,k)          &
     & +dt*cl*(                                          &
     & +(- radnflux1(mfr3,i+1,j,k)                    &
     &   + radnflux1(mfr3,i  ,j,k))/(x1a(i+1)-x1a(i)) &
     &      )

         Erad(i,j,k) =  max(Erad(i,j,k),Elte(i,j,k))

         fnl = sqrt(Frad(xdir,i,j,k)**2 +Frad(ydir,i,j,k)**2 +Frad(zdir,i,j,k)**2) ! dimensional
         flm =  fluxfactormax*Erad(i,j,k)
         if(fnl .gt. flm )then
            Frad(xdir,i,j,k) =  Frad(xdir,i,j,k)*flm/fnl
            Frad(ydir,i,j,k) =  Frad(ydir,i,j,k)*flm/fnl
            Frad(zdir,i,j,k) =  Frad(zdir,i,j,k)*flm/fnl
         endif
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
         x1out(is-gs:ie+gs,1) = x1b(is-gs:ie+gs)
         x1out(is-gs:ie+gs,2) = x1a(is-gs:ie+gs)

         radout(is-gs:ie+gs,js,ks,1) = Erad(  is-gs:ie+gs,js,ks)
         radout(is-gs:ie+gs,js,ks,2) = Frad(xdir,is-gs:ie+gs,js,ks)
         
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
         write(unitasc,"((a1,1x),((A),1x),(1PE15.4,1x))") "#","time=", time
         write(unitasc,"((a1,1x),((A),1x),(i0,1x))")      "#","  nx=", izones+2*gs
         write(unitasc,"(A)") "# x E Fx "
         k=ks
         j=js
         do i=is-gs,ie+gs
            write(unitasc,"(5(E15.6e3,1x))") x1b(i),Erad(i,j,k),Frad(xdir,i,j,k)
         enddo
         close(unitasc)
      endif
      print *, "output:",nout,time

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
