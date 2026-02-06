module units
  implicit none
  real(8),parameter::cl=3.0d10 ! cm/s
  real(8),parameter::arad=7.56e-15 ! erg/cm^3 deg^-4
end module units

module modelpara
  implicit none
  real(8),parameter::rho0=1.0d0,rho1=1.0d3 ! [g/cm^3]
  real(8),parameter::tempRad=1740 ! [K]
  real(8),parameter::tempMed=290! [K]
  real(8),parameter::Cv=20.79! [erg/cm^3/K]
end module modelpara

module commons
  use units
  implicit none
  
  integer::ntime
  integer,parameter::ntimemax=100000
  real(8)::time,dt
  real(8),parameter:: Coul=0.1d0
  data time / 0.0d0 /
  real(8),parameter:: timemax=1.0d-10
  real(8),parameter:: dtout=timemax/100
  
  integer,parameter::izones = 100*4
  integer,parameter::jzones = 24*4
  integer,parameter::mgn=2
  integer,parameter::in=izones+2*mgn+1 &
 &                  ,jn=jzones+2*mgn+1 &
 &                  ,kn=1
  integer,parameter::is=mgn+1 &
 &                  ,js=mgn+1 &
 &                  ,ks=1
  integer,parameter::ie=izones+mgn &
 &                  ,je=jzones+mgn &
 &                  ,ke=1

  real(8),parameter:: x1min=  0.0d0,x1max=1.0d0
  real(8),parameter:: x2min=  0.0d0,x2max=0.24d0
  real(8),dimension(in)::x1a,x1b
  real(8),dimension(jn)::x2a,x2b
  real(8),dimension(kn)::x3a,x3b

  real(8),dimension(in,jn,kn):: d,ei,tempK
  real(8),dimension(in,jn,kn):: kappa
  real(8),dimension(in,jn,kn):: Elte
  real(8):: eimin,eradmin

  real(8),dimension(    in,jn,kn)::Erad
  integer,parameter:: xdir=1,ydir=2,zdir=3
  real(8),dimension(xdir:zdir,  in,jn,kn)::Frad
  integer,parameter:: nxx=1,nyy=2,nzz=3,nxy=4,nyz=5,nzx=6
  real(8),dimension(nxx:nzx,in,jn,kn)::Prad
      
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
 &                  ,mferd=5,mffru=6,mffrv=7,mffrw=8 &
 &                  ,mcsp=9,mcsm=10

  integer,parameter::mradflx=4,mradadd=2
  
  integer,parameter:: merd=1,mfr1=2,mfr2=3,mfr3=4 &
 &                  , mfru=mufru,mfrv=mufrv,mfrw=mufrw
  real(8),dimension(mradflx,in,jn,kn):: radnflux1,radnflux2,radnflux3
  real(8),dimension(mradflx,in,jn,kn):: srcrad
  logical,parameter::flagEdd=.false.
  real(8),parameter::fluxfactormax=0.9d0

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
  write(6,*) "setup grids and fields"
  call GenerateGrid
  call GenerateProblem
  write(6,*) "entering main loop"
! main loop
  write(6,*)"step","time","dt"
  mloop: do ntime=1,ntimemax
     call TimestepControl
     if(mod(ntime,100) .eq. 0 ) write(6,*)ntime,time,dt
     call RadBoundaryCondition
     call StateRad
     call RadFlux1
     call RadFlux2
     call UpdateRadAdvection
     call UpdateRadSource
     time=time+dt
     call Output
     if(time > timemax) exit mloop
  enddo mloop

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

  dy=(x2max-x2min)/jzones
  do j=1,jn
     x2a(j) = dy*(j-(mgn+1))+x2min
  enddo
  do j=1,jn-1
     x2b(j) = 0.5d0*(x2a(j+1)+x2a(j))
  enddo

  return
end subroutine GenerateGrid

subroutine GenerateProblem
  use commons
  use eosmod
  use modelpara
  implicit none
  integer::i,j,k
  real(8)::delta

  write(6,*) "Erad med[erg/cm^3]",arad*tempMed**4
  write(6,*) "Erad in [erg/cm^3]",arad*tempRad**4
  eimin = rho0*Cv*tempMed
  eradmin = arad*tempMed**4

  do k=ks,ke
  do j=js,je
  do i=is,ie
     tempK(i,j,k) = tempMed
     delta = 10.0d0*(((x1b(i)-0.5)/0.1)**2+((x2b(j)-0.0)/0.06)**2-1.0d0)
          d(i,j,k) = rho0 + (rho1 - rho0)/(1+exp(delta))  ! [g cm^-3]
      kappa(i,j,k) = 0.1 *(tempK(i,j,k)/tempMed)**(-3.5d0) &
    &                    *(    d(i,j,k)/rho0   )**2/d(i,j,k)    ! [cm^2/g]

     ei(i,j,k)   = d(i,j,k)*Cv*tempK(i,j,k)
     Elte(i,j,k) = arad*(tempK(i,j,k))**4 ! [erg/cm^3]
       
  enddo
  enddo
  enddo

  Erad(:,:,:) = arad*(tempMed)**4

  do k=ks,ke
  do j=js,je
  do i=is,ie

!       Erad(  i,j,k) = 1.0d8
     Frad(xdir,i,j,k) = 0.0d0
     Frad(ydir,i,j,k) = 0.0d0
     Frad(zdir,i,j,k) = 0.0d0
       
  enddo
  enddo
  enddo
      
  call RadBoundaryCondition
!      write(6,*) "Gene",Erad(  is-1,js,ks),Erad(  is,js,ks)

  return
end subroutine GenerateProblem


subroutine RadBoundaryCondition
  use commons
  use fluxmod
  use modelpara
  implicit none
  integer::i,j,k

  do k=ks,ke
  do j=js,je
  do i=1,mgn
     Erad(    is-i,j,k) = arad*tempRad**4
     Frad(xdir,  is-i,j,k) = fluxfactormax*Erad(is-i,j,k)
     Frad(ydir:zdir,  is-i,j,k) = 0.0d0
 
  enddo
  enddo
  enddo

  do k=ks,ke
  do j=js,je
  do i=1,mgn
     Erad(    ie+i,j,k) = Erad(  ie,j,k) 
     Frad(:,  ie+i,j,k) = Frad(:,ie,j,k)
  enddo
  enddo
  enddo

! j-direction
  do k=ks,ke
  do i=is,ie
  do j=1,mgn
     Erad(  i,js-j,k) =  Erad(  i,js+j-1,k)
     Frad(xdir,i,js-j,k) =  Frad(xdir,i,js+j-1,k)
     Frad(ydir,i,js-j,k) =  Frad(ydir,i,js+j-1,k)
     Frad(zdir,i,js-j,k) =  Frad(zdir,i,js+j-1,k)
  enddo
  enddo
  enddo

  do k=ks,ke
  do i=is,ie
  do j=1,mgn
     Erad(  i,je+j,k) =  Erad(  i,je-j+1,k)
     Frad(xdir,i,je+j,k) =  Frad(xdir,i,je-j+1,k)
     Frad(ydir,i,je+j,k) =  Frad(ydir,i,je-j+1,k)
     Frad(zdir,i,je+j,k) =  Frad(zdir,i,je-j+1,k)
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

  dt = Coul * dtmin
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
  do j=1,jn-1
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
  real(8)::tiny=1.0d-30
  
  ul(1:mradflx) = leftst(1:mradflx)
  fl(1:mradflx) = leftst(mradflx+1:2*mradflx)
  ur(1:mradflx) = rigtst(1:mradflx)
  fr(1:mradflx) = rigtst(mradflx+1:2*mradflx)
  cslp = leftst(mcsp)
  cslm = leftst(mcsm)
  
  csrp = rigtst(mcsp)
  csrm = rigtst(mcsm)

  sl = min(cslm,csrm,-tiny)
  sr = max(cslp,csrp, tiny)
!  if(sr.eq.sl) write(6,*) "errorin HLL",sl,sr
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


subroutine RadFlux2
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
  do i=is,ie
  do j=js,je+1

     Pleftc1(:) = radsvc(:,i,j-2,k)
     Pleftc2(:) = radsvc(:,i,j-1,k)
     Prigtc1(:) = radsvc(:,i,j  ,k)
     Prigtc2(:) = radsvc(:,i,j+1,k)

! | Pleftc1   | Pleftc2 =>| Prigtc1   | Prigtc2   |
!                     You are here (x2-interface at j)

!====================
! Left state (from j-1 cell)
!====================
     dsvm(:) = Pleftc2(:) - Pleftc1(:)
     dsvp(:) = Prigtc1(:) - Pleftc2(:)
     call vanLeer(dsvp,dsvm,dsv)
     Plefte(:) = Pleftc2(:) + 0.5d0*dsv(:)

!====================
! Right state (from j cell)
!====================
     dsvm(:) = Prigtc1(:) - Pleftc2(:)
     dsvp(:) = Prigtc2(:) - Prigtc1(:)
     call vanLeer(dsvp,dsvm,dsv)
     Prigte(:) = Prigtc1(:) - 0.5d0*dsv(:)

!====================
! Convert primitives -> conservative/flux + characteristic speeds (M1 closure)
! (Note: swap of components follows original RadFlux2)
!====================
! ----- left
     fsq = (Plefte(nfr1)**2 +Plefte(nfr2)**2 +Plefte(nfr3)**2) ! (F/E/c)^2
     ffc = min(1.0d0,max(tiny,sqrt(fsq)))
     ff1 = min(1.0d0,max(-1.0d0,Plefte(nfr1)/ffc))
     ff2 = min(1.0d0,max(-1.0d0,Plefte(nfr2)/ffc))
     ff3 = min(1.0d0,max(-1.0d0,Plefte(nfr3)/ffc))

     chil = chi(ffc)
     if(flagEdd) chil = 1.0d0/3.0d0
     call lambdaPlusMinus(ffc,ff2,lp,lm)

     leftco(mcsp) = min(lp+tiny, 1.0d0)
     leftco(mcsm) = max(lm-tiny,-1.0d0)

     pxx = 0.5d0*(1.0d0-chil) + ff1**2 *0.5d0*(3.0d0*chil-1.0d0)
     pyy = 0.5d0*(1.0d0-chil) + ff2**2 *0.5d0*(3.0d0*chil-1.0d0)
     pzz = 0.5d0*(1.0d0-chil) + ff3**2 *0.5d0*(3.0d0*chil-1.0d0)
     pxy =                      ff1*ff2*0.5d0*(3.0d0*chil-1.0d0)
     pyz =                      ff2*ff3*0.5d0*(3.0d0*chil-1.0d0)
     pzx =                      ff3*ff1*0.5d0*(3.0d0*chil-1.0d0)

     leftco(muerd) =              Plefte(nerd) ! E
     leftco(mufrw) = Plefte(nfr1)*Plefte(nerd) ! F_u  (stored as w)
     leftco(mufru) = Plefte(nfr2)*Plefte(nerd) ! F_v  (stored as u)
     leftco(mufrv) = Plefte(nfr3)*Plefte(nerd) ! F_w  (stored as v)

     leftco(mferd) = Plefte(nfr2)*Plefte(nerd) ! F (normal component)
     leftco(mffrw) =          pxy*Plefte(nerd) ! P_u
     leftco(mffru) =          pyy*Plefte(nerd) ! P_v
     leftco(mffrv) =          pyz*Plefte(nerd) ! P_w

! ----- right
     fsq = (Prigte(nfr1)**2 +Prigte(nfr2)**2 +Prigte(nfr3)**2) ! (F/E/c)^2
     ffc = min(1.0d0,max(tiny,sqrt(fsq)))
     ff1 = min(1.0d0,max(-1.0d0,Prigte(nfr1)/ffc))
     ff2 = min(1.0d0,max(-1.0d0,Prigte(nfr2)/ffc))
     ff3 = min(1.0d0,max(-1.0d0,Prigte(nfr3)/ffc))

     chil = chi(ffc)
     if(flagEdd) chil = 1.0d0/3.0d0
     call lambdaPlusMinus(ffc,ff2,lp,lm)

     rigtco(mcsp) = min(lp+tiny, 1.0d0)
     rigtco(mcsm) = max(lm-tiny,-1.0d0)

     pxx = 0.5d0*(1.0d0-chil) + ff1**2 *0.5d0*(3.0d0*chil-1.0d0)
     pyy = 0.5d0*(1.0d0-chil) + ff2**2 *0.5d0*(3.0d0*chil-1.0d0)
     pzz = 0.5d0*(1.0d0-chil) + ff3**2 *0.5d0*(3.0d0*chil-1.0d0)
     pxy =                      ff1*ff2*0.5d0*(3.0d0*chil-1.0d0)
     pyz =                      ff2*ff3*0.5d0*(3.0d0*chil-1.0d0)
     pzx =                      ff3*ff1*0.5d0*(3.0d0*chil-1.0d0)

     rigtco(muerd) =                    Prigte(nerd) ! E
     rigtco(mufrw) = Prigte(nfr1)*Prigte(nerd) ! F_u  (stored as w)
     rigtco(mufru) = Prigte(nfr2)*Prigte(nerd) ! F_v  (stored as u)
     rigtco(mufrv) = Prigte(nfr3)*Prigte(nerd) ! F_w  (stored as v)

     rigtco(mferd) = Prigte(nfr2)*Prigte(nerd) ! F (normal component)
     rigtco(mffrw) =          pxy*Prigte(nerd) ! P_u
     rigtco(mffru) =          pyy*Prigte(nerd) ! P_v
     rigtco(mffrv) =          pyz*Prigte(nerd) ! P_w

!====================
! Riemann solver
!====================
     call HLLERAD(leftco,rigtco,nflux)
     radnflux2(merd,i,j,k) = nflux(merd)
     radnflux2(mfr1,i,j,k) = nflux(mfrw)
     radnflux2(mfr2,i,j,k) = nflux(mfru)
     radnflux2(mfr3,i,j,k) = nflux(mfrv)

  enddo
  enddo
  enddo
  return
end subroutine Radflux2



subroutine UpdateRadSource
  use commons
  use fluxmod
  use modelpara
  implicit none
  integer::i,j,k
  real(8)::pi
  real(8)::alpha
  real(8)::fnl,flm
  real(8)::etmp
  pi = acos(-1.0d0)

  do k=ks,ke
  do j=js,je
  do i=is,ie   
         kappa(i,j,k) = 0.1 *(tempK(i,j,k)/tempMed)**(-3.5) &
    &                       *(    d(i,j,k)/rho0   )**2/d(i,j,k)    ! [cm^2/g]
         Elte(i,j,k) = arad*(tempK(i,j,k))**4 ! [erg/cm^3]

         alpha = d(i,j,k)*kappa(i,j,k) *cl*dt
         etmp = Erad(i,j,k)
         Erad(i,j,k) = ( Erad(i,j,k)  + alpha* Elte(i,j,k) )/(1.0d0+alpha)
         Frad(xdir,i,j,k) =  Frad(xdir,i,j,k)/(1.0d0+alpha)
         Frad(ydir,i,j,k) =  Frad(ydir,i,j,k)/(1.0d0+alpha)
         Frad(zdir,i,j,k) =  Frad(zdir,i,j,k)/(1.0d0+alpha)           

         Erad(i,j,k) = max(Eradmin,Erad(i,j,k))

         ei(i,j,k) = ei(i,j,k) - (Erad(i,j,k) - etmp)
         ei(i,j,k) = max(eimin,ei(i,j,k))
         ! tempK(i,j,k) = ei(i,j,k)/d(i,j,k)/Cv
         
         fnl = sqrt(Frad(1,i,j,k)**2 +Frad(2,i,j,k)**2 +Frad(3,i,j,k)**2) ! dimensional
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
end subroutine  UpdateRadSource

subroutine UpdateRadAdvection
  use commons
  use fluxmod
  use modelpara
  implicit none
  integer::i,j,k
  !      write(6,*)ntime,Erad(is-1,js,ks),Erad(is,js,ks),Erad(is+1,js,ks)
      
  do k=ks,ke
  do j=js,je
  do i=is,ie
         
         Erad(i,j,k) = Erad(i,j,k)                          &
     & +(                                                   &
     & +(- radnflux1(merd,i+1,j,k)                          &
     &   + radnflux1(merd,i  ,j,k))/(x1a(i+1)-x1a(i))*dt*cl &
     & +(- radnflux2(merd,i,j+1,k)                          &
     &   + radnflux2(merd,i,j  ,k))/(x2a(j+1)-x2a(j))*dt*cl &
     &      )

         Frad(xdir,i,j,k) = Frad(xdir,i,j,k)                      &
     & +(                                                   &
     & +(- radnflux1(mfr1,i+1,j,k)                          &
     &   + radnflux1(mfr1,i  ,j,k))/(x1a(i+1)-x1a(i))*dt*cl &
     & +(- radnflux2(mfr1,i,j+1,k)                          &
     &   + radnflux2(mfr1,i,j  ,k))/(x2a(j+1)-x2a(j))*dt*cl &
     &      )

         Frad(ydir,i,j,k) = Frad(ydir,i,j,k)                      &
     & +(                                                   &
     & +(- radnflux1(mfr2,i+1,j,k)                          &
     &   + radnflux1(mfr2,i  ,j,k))/(x1a(i+1)-x1a(i))*dt*cl &
     & +(- radnflux2(mfr2,i,j+1,k)                          &
     &   + radnflux2(mfr2,i,j  ,k))/(x2a(j+1)-x2a(j))*dt*cl &
     &      )
         
         Frad(zdir,i,j,k) = Frad(zdir,i,j,k)                      &
     & +(                                                   &
     & +(- radnflux1(mfr3,i+1,j,k)                          &
     &   + radnflux1(mfr3,i  ,j,k))/(x1a(i+1)-x1a(i))*dt*cl &
     & +(- radnflux2(mfr3,i,j+1,k)                          &
     &   + radnflux2(mfr3,i,j  ,k))/(x2a(j+1)-x2a(j))*dt*cl &
     &      )
  enddo
  enddo
  enddo
      
  return
end subroutine UpdateRadAdvection

      subroutine Output
      use commons
      implicit none
      integer::i,j,k
      character(20),parameter::dirname="bindata/"
      character(40)::filename
      real(8),save::tout
      data tout / 0.0d0 /
      integer::nout
      data nout / 1 /
      integer,parameter:: unitout=17
      integer,parameter:: unitbin=13
      integer,parameter:: gs=1
      integer,parameter:: nvar=6
      real(8)::x1out(is-gs:ie+gs,2)
      real(8)::x2out(js-gs:je+gs,2)
      real(8)::radout(is-gs:ie+gs,js-gs:je+gs,ks,nvar)

      logical, save:: is_inited
      data is_inited /.false./

      if (.not. is_inited) then
         call makedirs("bindata")
         is_inited =.true.
      endif

      if(time .lt. tout+dtout) return

      write(filename,'(a3,i5.5,a4)')"unf",nout,".dat"
      filename = trim(dirname)//filename

      open(unitout,file=filename,status='replace',form='formatted')
      write(unitout,*) "# ",time,dt
      write(unitout,*) "# ",izones,gs
      write(unitout,*) "# ",jzones,gs
!      write(unitout,*) "# ",denup,dendn
      close(unitout)

      x1out(is-gs:ie+gs,1) = x1b(is-gs:ie+gs)
      x1out(is-gs:ie+gs,2) = x1a(is-gs:ie+gs)

      x2out(js-gs:je+gs,1) = x2b(js-gs:je+gs)
      x2out(js-gs:je+gs,2) = x2a(js-gs:je+gs)

      radout(is-gs:ie+gs,js-gs:je+gs,ks,1) = Erad(  is-gs:ie+gs,js-gs:je+gs,ks)
      radout(is-gs:ie+gs,js-gs:je+gs,ks,2) = Frad(1,is-gs:ie+gs,js-gs:je+gs,ks)
      radout(is-gs:ie+gs,js-gs:je+gs,ks,3) = Frad(2,is-gs:ie+gs,js-gs:je+gs,ks)
      radout(is-gs:ie+gs,js-gs:je+gs,ks,4) = Frad(3,is-gs:ie+gs,js-gs:je+gs,ks)
      radout(is-gs:ie+gs,js-gs:je+gs,ks,5) =    d(  is-gs:ie+gs,js-gs:je+gs,ks)
      radout(is-gs:ie+gs,js-gs:je+gs,ks,6) =   ei(  is-gs:ie+gs,js-gs:je+gs,ks)

      write(filename,'(a3,i5.5,a4)')"bin",nout,".dat"
      filename = trim(dirname)//filename
      open(unitbin,file=filename,status='replace',form='binary') 
      write(unitbin) x1out(:,:)
      write(unitbin) x2out(:,:)
      write(unitbin) radout(:,:,:,:)
      close(unitbin)

      write(6,*) "output:",nout,time

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
