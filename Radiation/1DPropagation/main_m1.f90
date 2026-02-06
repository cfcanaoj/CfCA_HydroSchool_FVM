module units
  implicit none
  real(8),parameter::cl=3.0d10 ! cm/s
end module units

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
  real(8),dimension(mradflx,in,jn,kn):: srcrad

end module fluxmod

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
     call SourceRad
     call RadFlux1
     call UpdateRad
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
      implicit none
      integer::i,j,k
      
      do k=ks,ke
      do j=js,je
      do i=1,in

!          d(i,j,k) = 25.0   ! [g cm^-3]
          d(i,j,k) = 0.25   ! [g cm^-3]
      kappa(i,j,k) = 0.04   ! [cm^2/g]
       Elte(i,j,k) = 1.0d10 ! [erg/cm^3]

      enddo
      enddo
      enddo
      
      do k=ks,ke
      do j=js,je
      do i=1,in

       Erad(  i,j,k) = 1.0d10
       Frad(xdir,i,j,k) = 0.0d0
       Frad(ydir,i,j,k) = 0.0d0
       Frad(zdir,i,j,k) = 0.0d0

      enddo
      enddo
      enddo
      
   call RadBoundaryCondition
   write(6,*) "Gene",Erad(  is-1,js,ks),Erad(  is,js,ks)

   return
end subroutine GenerateProblem


      subroutine RadBoundaryCondition
      use commons
      implicit none
      integer::i,j,k

      do k=ks,ke
      do j=js,je
      do i=1,mgn
          Erad(    is-i,j,k) = 1.0d20
          if(d(is-i,j,k) .gt. 1.0d0) then
             Frad(xdir,  is-i,j,k) = Erad(    is-i,j,k)*cl/sqrt(3.0d0)
          else
             Frad(xdir,  is-i,j,k) = Erad(    is-i,j,k)*cl*0.99
          endif
          Frad(ydir:zdir,  is-i,j,k) = 0.0d0
 
      enddo
      enddo
      enddo

      do k=ks,ke
      do j=js,je
      do i=1,mgn
          Erad(    ie+i,j,k) = Erad(  ie,j,k) 
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
      do j=js,je
      do i=1,in-1
         radsvc(nerd,i,j,k) = Erad(  i,j,k) 
         radsvc(nfr1,i,j,k) = Frad(xdir,i,j,k) 
         radsvc(nfr2,i,j,k) = Frad(ydir,i,j,k) 
         radsvc(nfr3,i,j,k) = Frad(zdir,i,j,k) 

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
      use commons, only: is,ie,in,js,je,jn,ks,ke,kn
      use fluxmod
      use units
      implicit none
      integer::i,j,k
      real(8),dimension(nrad):: dsv,dsvm,dsvp
      real(8),dimension(nrad,in,jn,kn):: leftpr,rigtpr
      real(8),dimension(2*mradflx+mradadd,in,jn,kn):: leftco,rigtco
      real(8),dimension(2*mradflx+mradadd):: leftst,rigtst
      real(8),dimension(mradflx):: nflux
      real(8):: fsq,ffc,chi,pxx,pyy,pzz,pxy,pyz,pzx
      real(8),parameter::tiny=1.0d-30

      do k=ks,ke
      do j=js,je
      do i=is-1,ie+1

           dsvm(:) =                   radsvc(:,i,j,k)-radsvc(:,i-1,j,k)
           dsvp(:) = radsvc(:,i+1,j,k)-radsvc(:,i,j,k)

          call  vanLeer(dsvp,dsvm,dsv)
         leftpr(:,i+1,j,k) = radsvc(:,i,j,k)+0.5d0*dsv(:)
         rigtpr(:,i  ,j,k) = radsvc(:,i,j,k)-0.5d0*dsv(:)
      enddo
      enddo
      enddo

      do k=ks,ke
      do j=js,je
      do i=is,ie+1
! left
         fsq = (leftpr(nfr1,i,j,k)**2 +leftpr(nfr2,i,j,k)**2 +leftpr(nfr3,i,j,k)**2) ! dimensional
         ffc = fsq/cl**2/leftpr(nerd,i,j,k)**2 ! non-dimensional
         ffc = min(1.0,max(0.0,ffc))
         chi = (3.0+4.0*ffc)/(5.0+2.0*sqrt(4.0-3.0*ffc))
!         chi = 1.0d0/3.0 ! Eddington approximation

         leftco(mcsp,i,j,k) = (+1.0d0/sqrt(3.0d0)*(3.0*(1.0-chi)/2.0)+max(-1.0d0,min(1.0d0,leftco(nfr1,i,j,k)/sqrt(fsq+tiny)))*((3.0*chi-1.0)/2))
         leftco(mcsm,i,j,k) = (-1.0d0/sqrt(3.0d0)*(3.0*(1.0-chi)/2.0)+max(-1.0d0,min(1.0d0,leftco(nfr1,i,j,k)/sqrt(fsq+tiny)))*((3.0*chi-1.0)/2))
         leftco(mcsp,i,j,k) = cl*min(leftco(mcsp,i,j,k), 1.0d0)
         leftco(mcsm,i,j,k) = cl*max(leftco(mcsm,i,j,k),-1.0d0)

         pxx = leftpr(nerd,i,j,k)*(1.0/3.0d0*(3.0*(1.0-chi)/2.0) +leftpr(nfr1,i,j,k)**2/(fsq+tiny)*((3.0*chi-1.0)/2.0))
         pyy = leftpr(nerd,i,j,k)*(1.0/3.0d0*(3.0*(1.0-chi)/2.0) +leftpr(nfr2,i,j,k)**2/(fsq+tiny)*((3.0*chi-1.0)/2.0))
         pzz = leftpr(nerd,i,j,k)*(1.0/3.0d0*(3.0*(1.0-chi)/2.0) +leftpr(nfr3,i,j,k)**2/(fsq+tiny)*((3.0*chi-1.0)/2.0))
         pxy = leftpr(nerd,i,j,k)*leftpr(nfr1,i,j,k)*leftpr(nfr2,i,j,k)/(fsq+tiny)*((3.0*chi-1.0)/2)
         pyz = leftpr(nerd,i,j,k)*leftpr(nfr2,i,j,k)*leftpr(nfr3,i,j,k)/(fsq+tiny)*((3.0*chi-1.0)/2)
         pzx = leftpr(nerd,i,j,k)*leftpr(nfr3,i,j,k)*leftpr(nfr1,i,j,k)/(fsq+tiny)*((3.0*chi-1.0)/2)

         leftco(muerd,i,j,k) = leftpr(nerd,i,j,k) ! E
         leftco(mufru,i,j,k) = leftpr(nfr1,i,j,k) ! F_u
         leftco(mufrv,i,j,k) = leftpr(nfr2,i,j,k) ! F_v
         leftco(mufrw,i,j,k) = leftpr(nfr3,i,j,k) ! F_w

         leftco(mferd,i,j,k) = leftpr(nfr1,i,j,k) ! F
         leftco(mffru,i,j,k) = pxx*cl**2 ! P_u
         leftco(mffrv,i,j,k) = pxy*cl**2 ! P_v
         leftco(mffrw,i,j,k) = pzx*cl**2 ! P_w

! right 
         fsq = (rigtpr(nfr1,i,j,k)**2 +rigtpr(nfr2,i,j,k)**2 +rigtpr(nfr3,i,j,k)**2) ! dimensional
         ffc = fsq/cl**2/rigtpr(nerd,i,j,k)**2 ! non-dimensional
         ffc = min(1.0,max(0.0,ffc))
         chi = (3.0+4.0*ffc)/(5.0+2.0*sqrt(4.0-3.0*ffc))
!         chi = 1.0d0/3.0 ! Eddington approximation

         rigtco(mcsp,i,j,k) = (+1.0d0/sqrt(3.0d0)*(3.0*(1.0-chi)/2.0)+max(-1.0d0,min(1.0d0,rigtco(nfr1,i,j,k)/sqrt(fsq+tiny)))*((3.0*chi-1.0)/2))
         rigtco(mcsm,i,j,k) = (-1.0d0/sqrt(3.0d0)*(3.0*(1.0-chi)/2.0)+max(-1.0d0,min(1.0d0,rigtco(nfr1,i,j,k)/sqrt(fsq+tiny)))*((3.0*chi-1.0)/2))
         rigtco(mcsp,i,j,k) = cl*min(rigtco(mcsp,i,j,k), 1.0d0)
         rigtco(mcsm,i,j,k) = cl*max(rigtco(mcsm,i,j,k),-1.0d0)

         pxx = rigtpr(nerd,i,j,k)*(1.0/3.0d0*(3.0*(1.0-chi)/2.0) +rigtpr(nfr1,i,j,k)**2/(fsq+tiny)*((3.0*chi-1.0)/2.0))
         pyy = rigtpr(nerd,i,j,k)*(1.0/3.0d0*(3.0*(1.0-chi)/2.0) +rigtpr(nfr2,i,j,k)**2/(fsq+tiny)*((3.0*chi-1.0)/2.0))
         pzz = rigtpr(nerd,i,j,k)*(1.0/3.0d0*(3.0*(1.0-chi)/2.0) +rigtpr(nfr3,i,j,k)**2/(fsq+tiny)*((3.0*chi-1.0)/2.0))
         pxy = rigtpr(nerd,i,j,k)*rigtpr(nfr1,i,j,k)*rigtpr(nfr2,i,j,k)/(fsq+tiny)*((3.0*chi-1.0)/2.0)
         pyz = rigtpr(nerd,i,j,k)*rigtpr(nfr2,i,j,k)*rigtpr(nfr3,i,j,k)/(fsq+tiny)*((3.0*chi-1.0)/2.0)
         pzx = rigtpr(nerd,i,j,k)*rigtpr(nfr3,i,j,k)*rigtpr(nfr1,i,j,k)/(fsq+tiny)*((3.0*chi-1.0)/2.0)

         rigtco(muerd,i,j,k) = rigtpr(nerd,i,j,k) ! E
         rigtco(mufru,i,j,k) = rigtpr(nfr1,i,j,k) ! F_u
         rigtco(mufrv,i,j,k) = rigtpr(nfr2,i,j,k) ! F_v
         rigtco(mufrw,i,j,k) = rigtpr(nfr3,i,j,k) ! F_w

         rigtco(mferd,i,j,k) = rigtpr(nfr1,i,j,k) ! F
         rigtco(mffru,i,j,k) = pxx*cl**2 ! P_u
         rigtco(mffrv,i,j,k) = pxy*cl**2 ! P_v
         rigtco(mffrw,i,j,k) = pzx*cl**2 ! P_w

!         if(i .eq. is) print *, rigtco(muerd,i,j,k), leftco(muerd,i,j,k)
      enddo
      enddo
      enddo

      do k=ks,ke
      do j=js,je
      do i=is,ie+1
         leftst(:)=leftco(:,i,j,k)
         rigtst(:)=rigtco(:,i,j,k)
 !        if(i .eq. is) print *, rigtst(merd), leftst(merd)
         call HLLERAD(leftst,rigtst,nflux)
         radnflux1(merd,i,j,k)=nflux(merd)
         radnflux1(mfr1,i,j,k)=nflux(mfru)
         radnflux1(mfr2,i,j,k)=nflux(mfrv)
         radnflux1(mfr3,i,j,k)=nflux(mfrw)
  !       if(i .eq. is) print *, radnflux1(merd,i,j,k)
  !       stop
      enddo
      enddo
      enddo

      return
      end subroutine Radflux1

      subroutine SourceRad
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
         srcrad(merd,i,j,k) = (Elte(i,j,k)-Erad(i,j,k))*alpha/(1+alpha)/dt
         srcrad(mfr1,i,j,k) =  -alpha/(1.0d0+alpha)/dt*Frad(xdir,i,j,k)
         srcrad(mfr2,i,j,k) =  -alpha/(1.0d0+alpha)/dt*Frad(ydir,i,j,k)
         srcrad(mfr3,i,j,k)=   -alpha/(1.0d0+alpha)/dt*Frad(zdir,i,j,k)

      enddo
      enddo
      enddo

      return
      end subroutine SourceRad

      subroutine UpdateRad
      use commons
      use fluxmod
      implicit none
      real(8)::fnl,flm
      integer::i,j,k

      do k=ks,ke
      do j=js,je
      do i=is,ie
         
         Erad(i,j,k) = Erad(i,j,k)                    &
     & +dt*(                                          &
     & +(- radnflux1(merd,i+1,j,k)                    &
     &   + radnflux1(merd,i  ,j,k))/(x1a(i+1)-x1a(i)) &
     &   +srcrad(merd,i,j,k)                          &
     &      )

         Frad(xdir,i,j,k) = Frad(xdir,i,j,k)          &
     & +dt*(                                          &
     & +(- radnflux1(mfr1,i+1,j,k)                    &
     &   + radnflux1(mfr1,i  ,j,k))/(x1a(i+1)-x1a(i)) &
     &   + srcrad(mfr1,i,j,k)                         &
     &      )

         Frad(ydir,i,j,k) = Frad(ydir,i,j,k)          &
     & +dt*(                                          &
     & +(- radnflux1(mfr2,i+1,j,k)                    &
     &   + radnflux1(mfr2,i  ,j,k))/(x1a(i+1)-x1a(i)) &
     &   + srcrad(mfr2,i,j,k)                         &
     &      )

         Frad(zdir,i,j,k) = Frad(zdir,i,j,k)          &
     & +dt*(                                          &
     & +(- radnflux1(mfr3,i+1,j,k)                    &
     &   + radnflux1(mfr3,i  ,j,k))/(x1a(i+1)-x1a(i)) &
     &   + srcrad(mfr3,i,j,k)                         &
     &      )

         Erad(i,j,k) =  max(Erad(i,j,k),Elte(i,j,k))

         fnl = sqrt(Frad(xdir,i,j,k)**2 +Frad(ydir,i,j,k)**2 +Frad(zdir,i,j,k)**2) ! dimensional
         flm =  0.99*Erad(i,j,k)*cl
         if(fnl/cl .gt. flm )then
            Frad(xdir,i,j,k) =  Frad(xdir,i,j,k)*flm/fnl
            Frad(ydir,i,j,k) =  Frad(ydir,i,j,k)*flm/fnl
            Frad(zdir,i,j,k) =  Frad(zdir,i,j,k)*flm/fnl
         endif
      enddo
      enddo
      enddo

      return
      end subroutine UpdateRad
    
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
      integer::nout
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
         write(unitbin) izones,gs
         write(unitbin) 4
         write(unitbin) x1out(:,:)
         write(unitbin) radout(:,:,:,:)
         close(unitbin)
      else

         write(filename,'(a4,i5.5,a4)')"snap",nout,".dat"
         filename = trim(dirname)//filename
         open(newunit=unitasc,file=filename,status='replace',form='formatted',access="stream",action="write") 
         write(unitasc,"(a1,(1x,(A)),(1x,1PE15.4))") "#","time=",time
         write(unitasc,"(a1,(1x,(A)),(1x,i0))") "#","nx=", izones+2*gs
         k=ks
         j=js
         do i=is-gs,ie+gs
            write(unitasc,*) x1b(i),Erad(i,j,k),Erad(i,j,k),Frad(xdir,i,j,k)
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
