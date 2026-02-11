module hyddata
  implicit none
  integer:: ntime
  integer,parameter:: ntimemax=20000
  real(8):: time, dt
  data time / 0.0d0 /
  real(8),parameter:: timemax = 1.0d0
  real(8),parameter:: dtout = timemax/10.0d0
  
  integer,parameter:: nx = 30
  integer,parameter:: ny = 30
  integer,parameter:: nz = 60
  integer,parameter:: ngs= 2
  integer,parameter:: nxtot = nx+2*ngs+1 &
       &             ,nytot = ny+2*ngs+1 &
       &             ,nztot = nz+2*ngs+1
  integer,parameter::is=ngs+1 &
       &            ,js=ngs+1 &
       &            ,ks=ngs+1
  integer,parameter::ie=nx+ngs &
       &            ,je=ny+ngs &
       &            ,ke=nz+ngs
  
  real(8),parameter:: xmin=-0.5d0, xmax=0.5d0
  real(8),parameter:: ymin=-0.5d0, ymax=0.5d0
  real(8),parameter:: zmin=-1.0d0, zmax=1.0d0
  real(8),dimension(nxtot)::xf,xv
  real(8),dimension(nytot)::yf,yv
  real(8),dimension(nztot)::zf,zv

  integer, parameter :: IDN = 1
  integer, parameter :: IPR = 2
  integer, parameter :: IVX = 3
  integer, parameter :: IVY = 4
  integer, parameter :: IVZ = 5
  integer, parameter :: NVAR = 5
  real(8),dimension(nxtot,nytot,nztot,NVAR)::Q
end module hyddata

program main
  use hyddata
  implicit none
  
  dt =1.0d-3
  call GenerateGrid
  call GenerateProblem
  mloop: do ntime=1,ntimemax
     write(6,*)ntime,time,dt
     time = time+dtout
     call UpdatePrim
     call Output
     if(time > timemax) exit mloop
  enddo mloop
  
end program main

!-------------------------------------------------------------------
!       generate grid 
!-------------------------------------------------------------------
subroutine GenerateGrid
  use hyddata
  implicit none
  real(8) :: dx,dy,dz
  integer::i,j,k

    dx=(xmax-xmin)/nx
    do i=1,nxtot
         xf(i) = dx*(i-(ngs+1))+xmin
    enddo
    do i=1,nxtot-1
         xv(i) = 0.5d0*(xf(i+1)+xf(i))
    enddo

    dy=(ymax-ymin)/ny
    do j=1,nytot
         yf(j) = dy*(j-(ngs+1))+ymin
    enddo
    do j=1,nytot-1
         yv(j) = 0.5d0*(yf(j+1)+yf(j))
    enddo

    dz=(zmax-zmin)/nz
    do k=1,nztot
         zf(k) = dz*(k-(ngs+1))+zmin
    enddo
    do i=1,nxtot-1
         zv(k) = 0.5d0*(zf(k+1)+zf(k))
    enddo

return
end subroutine GenerateGrid

!-------------------------------------------------------------------
!       generate initial condition 
!-------------------------------------------------------------------
subroutine GenerateProblem
  use hyddata
  implicit none
  integer:: i,j,k
  real(8):: r

  do k=ks,ke
  do j=js,je
  do i=is,ie
     r = sqrt(xv(i)**2 +yv(j)**2 +zv(k)**2) 
     if( r .le. 0.05) then
        Q(i,j,k,IDN) = 2
        Q(i,j,k,IPR) = 4
        Q(i,j,k,IVX) = 0.0d0
        Q(i,j,k,IVY) = 0.0d0
        Q(i,j,k,IVZ) = 0.0d0
     else if( r .le. 0.1) then
        Q(i,j,k,IDN) = 1
        Q(i,j,k,IPR) = 2
        Q(i,j,k,IVX) = 0.0d0
        Q(i,j,k,IVY) = 0.0d0
        Q(i,j,k,IVZ) = 0.0d0
     else
        Q(i,j,k,IDN) = 1.0d-2
        Q(i,j,k,IPR) = 2.0d-2
        Q(i,j,k,IVX) = 0.0d0
        Q(i,j,k,IVY) = 0.0d0
        Q(i,j,k,IVZ) = 0.0d0
     endif

  enddo
  enddo
  enddo
      
  return
end subroutine GenerateProblem

!-------------------------------------------------------------------
!       Boundary Condition 
!-------------------------------------------------------------------
subroutine BoundaryCondition
  use hyddata
  implicit none
  integer:: i,j,k

!-------------------------------------------------------------------
  do k=ks,ke
  do j=js,je
  do i=1,ngs
     Q(is-i,j,k,IDN) = Q(is-1+i,j,k,IDN)
     Q(is-i,j,k,IPR) = Q(is-1+i,j,k,IPR)
     Q(is-i,j,k,IVX) = Q(is-1+i,j,k,IVX)
     Q(is-i,j,k,IVY) = Q(is-1+i,j,k,IVY)
     Q(is-i,j,k,IVZ) = Q(is-1+i,j,k,IVZ)
  enddo
  enddo
  enddo

  do k=ks,ke
  do j=js,je
  do i=1,ngs
     Q(ie+i,j,k,IDN) = Q(ie-i+1,j,k,IDN)
     Q(ie+i,j,k,IPR) = Q(ie-i+1,j,k,IPR)
     Q(ie+i,j,k,IVX) = Q(ie-i+1,j,k,IVX)
     Q(ie+i,j,k,IVY) = Q(ie-i+1,j,k,IVY)
     Q(ie+i,j,k,IVZ) = Q(ie-i+1,j,k,IVZ)
  enddo
  enddo
  enddo
!-------------------------------------------------------------------
  do k=ks,ke
  do i=is,ie
  do j=1,ngs 
     Q(i,js-j,k,IDN) = Q(i,js-1+j,k,IDN)
     Q(i,js-j,k,IPR) = Q(i,js-1+j,k,IPR)
     Q(i,js-j,k,IVX) = Q(i,js-1+j,k,IVX)
     Q(i,js-j,k,IVY) = Q(i,js-1+j,k,IVY)
     Q(i,js-j,k,IVZ) = Q(i,js-1+j,k,IVZ)
  enddo
  enddo
  enddo

  do k=ks,ke
  do i=is,ie
  do j=1,ngs
     Q(i,je+j,k,IDN) = Q(i,je-j+1,k,IDN)
     Q(i,je+j,k,IPR) = Q(i,je-j+1,k,IPR)
     Q(i,je+j,k,IVX) = Q(i,je-j+1,k,IVX)
     Q(i,je+j,k,IVY) = Q(i,je-j+1,k,IVY)
     Q(i,je+j,k,IVZ) = Q(i,je-j+1,k,IVZ)
  enddo
  enddo
  enddo
!-------------------------------------------------------------------
  do j=js,je
  do i=is,ie
  do k=1,ngs 
     Q(i,j,ks-k,IDN) = Q(i,js-1+j,k,IDN)
     Q(i,j,ks-k,IPR) = Q(i,js-1+j,k,IPR)
     Q(i,j,ks-k,IVX) = Q(i,js-1+j,k,IVX)
     Q(i,j,ks-k,IVY) = Q(i,js-1+j,k,IVY)
     Q(i,j,ks-k,IVZ) = Q(i,js-1+j,k,IVZ)
  enddo
  enddo
  enddo

  do j=js,je
  do i=is,ie
  do k=1,ngs
     Q(i,j,ke+k,IDN) = Q(i,j,ke-k+1,IDN)
     Q(i,j,ke+k,IPR) = Q(i,j,ke-k+1,IPR)
     Q(i,j,ke+k,IVX) = Q(i,j,ke-k+1,IVX)
     Q(i,j,ke+k,IVY) = Q(i,j,ke-k+1,IVY)
     Q(i,j,ke+k,IVZ) = Q(i,j,ke-k+1,IVZ)
  enddo
  enddo
  enddo
!-------------------------------------------------------------------

  return
end subroutine BoundaryCondition
!-------------------------------------------------------------------
!       Primitive variables
!-------------------------------------------------------------------
subroutine UpdatePrim
  use hyddata
  implicit none
  
  real(8):: r,xc,yc,zc
  real(8),parameter:: freq=1.0d0
  real(8):: pi
  integer:: i,j,k
  pi  = acos(-1.d0)
  xc = 0.0d0
  yc = 0.0d0
  zc = zmax*sin(2*pi*freq*time)

  do k=ks,ke
  do j=js,je
  do i=is,ie
     r = sqrt((xf(i)-xc)**2 +(yf(j)-yc)**2 +(zf(k)-zc)**2) 
     if( r .le. 0.05) then
        Q(i,j,k,IDN) = 2
        Q(i,j,k,IPR) = 4
        Q(i,j,k,IVX) = 0.0d0
        Q(i,j,k,IVY) = 0.0d0
        Q(i,j,k,IVZ) = 0.0d0
     else if( r .le. 0.1) then
        Q(i,j,k,IDN) = 1
        Q(i,j,k,IPR) = 2
        Q(i,j,k,IVX) = 0.0d0
        Q(i,j,k,IVY) = 0.0d0
        Q(i,j,k,IVZ) = 2*pi*freq*zmax*cos(2*pi*freq*time)
     else
        Q(i,j,k,IDN) = 1.0d-2
        Q(i,j,k,IPR) = 2.0d-2
        Q(i,j,k,IVX) = 0.0d0
        Q(i,j,k,IVY) = 0.0d0
        Q(i,j,k,IVZ) = 0.0d0
     endif
  enddo
  enddo
  enddo
  call BoundaryCondition
  return
end subroutine UpdatePrim

subroutine Output
  use hyddata
  use hdf5
  implicit none 
  integer,parameter::ndim=4                   ! max size of the dimension for output file
  integer(hid_t) :: file_id                   ! file identifier
  integer :: rank                             ! dataset rank
  integer(hsize_t), dimension(ndim) :: dims   ! dataset dimensions
  character(80) :: dsetname                   ! dataset name
  integer :: error
  character*(80) :: fname,fpath
  character*(80) :: fnxmf,fpxmf
  integer unitxmf
  integer,save:: imax,jmax,kmax
  character*(80) :: DIM3D,DIM1D,TIMEEXP
  
  character(40),parameter:: hdfdir = "./hdfdata/" 

  character(3):: headerid

  integer,parameter:: nsfld=2
  integer,parameter:: nvfld=1
  integer,save:: timeindex
  data timeindex /0/

! Variable to esatimate NODE value

  real*4,dimension(:,:,:,:), allocatable,save:: QShdfF
  real*4,dimension(:,:,:,:,:), allocatable,save:: QVhdfF
  real*4,dimension(:),allocatable,save::XhdfF,YhdfF,ZhdfF
  character(len=20,kind=1):: QSname(nsfld),QVname(nvfld)

  integer:: i,j,k,n
  
  logical,save::is_inited
  data is_inited / .false. /

  !      if(.not. mod(ntime,10) ==0 )return

!-------------------------------------------------
! interface part begin
!-------------------------------------------------

! Initialize
  if(.not. is_inited)then
     imax = nxtot-2*ngs-1
     jmax = nytot-2*ngs-1
     kmax = nztot-2*ngs-1

     allocate(XhdfF(imax+1))
     allocate(YhdfF(jmax+1))
     allocate(ZhdfF(kmax+1))
     allocate(QShdfF(  1:imax  ,1:jmax  ,1:kmax  ,nsfld))
     allocate(QVhdfF(3,1:imax  ,1:jmax  ,1:kmax  ,nvfld))
         
     XhdfF(1:imax+1)= real(xf(1:imax+1))
     YhdfF(1:jmax+1)= real(yf(1:jmax+1))
     ZhdfF(1:kmax+1)= real(zf(1:kmax+1))
  
     call makedirs(hdfdir)

     is_inited = .true.
  endif

  headerid="hyd"

  write(TIMEEXP,"(f6.2)") time

  QSname(1) = "density"
  QSname(2) = "pressure"
  QVname(1) = "velocity"

!  print *, Qname(1)
!  print *, Qname(2)
!  print *, Qname(3)
!  print *, Qname(4)
!  print *, Qname(5)

  do k=1,kmax
  do j=1,jmax
  do i=1,imax
     QShdfF(i,j,k,1) = real( Q(ngs+i,ngs+j,ngs+k,1) )
     QShdfF(i,j,k,2) = real( Q(ngs+i,ngs+j,ngs+k,2) )

     QVhdfF(1,i,j,k,1) = real( Q(ngs+i,ngs+j,ngs+k,3) )
     QVhdfF(2,i,j,k,1) = real( Q(ngs+i,ngs+j,ngs+k,4) )
     QVhdfF(3,i,j,k,1) = real( Q(ngs+i,ngs+j,ngs+k,5) )
  enddo
  enddo
  enddo

      
  timeindex = timeindex + 1

!-------------------------------------------------
! interface part end
!-------------------------------------------------

      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Writing HDF5 file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Determine the output filename
  write(fname,"(a3,i5.5,a3)") headerid,timeindex,'.h5'
  fpath =  trim(hdfdir)//fname

  write(fnxmf,"(a3,i5.5,a4)") headerid,timeindex,'.xmf'
  fpxmf =  trim(hdfdir)//fnxmf

! Open HDF5
  CALL h5open_f (error)
! Note:
! H5F_ACC_TRUNC_F: delete the file if that is exist
! H5F_ACC_EXCL_F : fail writing to the file if that  is exist
  CALL h5fcreate_f(fpath, H5F_ACC_TRUNC_F, file_id, error)


      dims(:) = 0
      rank = 1
      dims(1) = imax+1
      dsetname="/x"
      call  hdf5OUT(file_id,rank,dims,dsetname,XhdfF)
      rank=1
      dims(1) = jmax+1
      dsetname="/y"
      call  hdf5OUT(file_id,rank,dims,dsetname,YhdfF)
      rank=1
      dims(1) = kmax+1
      dsetname="/z"
      call  hdf5OUT(file_id,rank,dims,dsetname,ZhdfF)

      rank=3
      dims(1) = imax
      dims(2) = jmax
      dims(3) = kmax
      do n=1,nsfld
         dsetname=QSname(n)
         call  hdf5OUT(file_id,rank,dims,dsetname,QShdfF(1,1,1,n))
      enddo

      rank=4
      dims(1) = 3
      dims(2) = imax
      dims(3) = jmax
      dims(4) = kmax

      do n=1,nvfld
         dsetname=QVname(n)
         call  hdf5OUT(file_id,rank,dims,dsetname,QVhdfF(1,1,1,1,n))
      enddo


  unitxmf=1221
  open(unitxmf,file=fpxmf,status='unknown',form='formatted')

! Header
      write(unitxmf,'(a)')'<?xml version="1.0" ?>'
      write(unitxmf,'(a)')'<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
      write(unitxmf,'(a)')'<Xdmf Version="2.0">'
      write(unitxmf,'(a)')'  <Domain>'
      write(unitxmf,'(a)')'    <Grid Name="mesh" GridType="Uniform">'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Time
      write(unitxmf,'(a)')'      <Time Value="'//trim(TIMEEXP)//'"/>'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Coordinate
      write(DIM3D,"(I0,1x,I0,1x,I0)") kmax+1,jmax+1,imax+1
      write(unitxmf,'(a)')'      <Topology' &
     &               //' TopologyType="3DRectMesh" NumberOfElements="' &
     &               //trim(DIM3D)//'"/>'
      write(unitxmf,'(a)')'      <Geometry GeometryType="VXVYVZ">'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! x coordinate
      write(DIM1D,"(I0)") imax+1
      dsetname = "x"
      write(unitxmf,'(a)')&
     &                '        <DataItem Dimensions="' // trim(DIM1D)&
     &              //'" Name="' // trim(dsetname) // '"'
      write(unitxmf,'(a)') '          ' // &
     &               'NumberType="Float" Precision="4" Format="HDF">'
      write(unitxmf,'(a)')'	  ' // trim(fname) //":/x"
      write(unitxmf,'(a)')'        </DataItem>'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! y coordinate
      write(DIM1D,"(I0)") jmax+1
      dsetname="y"
      write(unitxmf,*)'        <DataItem Dimensions="'//trim(DIM1D)&
     &              //'" Name="' // trim(dsetname) // '"'
      write(unitxmf,*) '          ' // &
     &              'NumberType="Float" Precision="4" Format="HDF">'
      write(unitxmf,'(a)')'	  '//trim(fname)//":/y"
      write(unitxmf,'(a)')'        </DataItem>'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! z coordinate
      write(DIM1D,"(I0)") kmax+1
      dsetname="z"
      write(unitxmf,'(a)')'        <DataItem Dimensions="'//trim(DIM1D) &
     &              //'" Name="'//trim(dsetname) // '"'
      write(unitxmf,'(a)') '          ' // &
     &              'NumberType="Float" Precision="4" Format="HDF">'
      write(unitxmf,'(a)')'	  '//trim(fname)//":/z"
      write(unitxmf,'(a)')'        </DataItem>'
      write(unitxmf,'(a)')'      </Geometry>'

      do n=1,nsfld
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      dsetname=Qsname(n)
      write(unitxmf,'(a)')'      <Attribute Name="'//trim(dsetname) &
     &              //'" AttributeType="Scalar" Center="Cell">'
      write(DIM3D,"(I0,1x,I0,1x,I0)") kmax,jmax,imax
      write(unitxmf,'(a)')'        ' &
     &                    //'<DataItem Dimensions="'//trim(DIM3D)//'"'
      write(unitxmf,'(a)')'          ' // &
     &                'NumberType="Float" Precision="4" Format="HDF">'
      write(unitxmf,'(a)')'	  '//trim(fname)//":/" // QSname(n)
      write(unitxmf,'(a)')'        </DataItem>'
      write(unitxmf,'(a)')'      </Attribute>'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      enddo

      do n=1,nvfld
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      dsetname=Qvname(n)
      write(unitxmf,'(a)')'      <Attribute Name="'//trim(dsetname) &
     &              //'" AttributeType="Vector" Center="Cell">'
      write(DIM3D,"(I0,1x,I0,1x,I0,1x,I0)") kmax,jmax,imax,3
      write(unitxmf,'(a)')'        ' &
     &                    //'<DataItem Dimensions="'//trim(DIM3D)//'"'
      write(unitxmf,'(a)')'          ' // &
     &                'NumberType="Float" Precision="4" Format="HDF">'
      write(unitxmf,'(a)')'	  '//trim(fname)//":/" // Qvname(n)
      write(unitxmf,'(a)')'        </DataItem>'
      write(unitxmf,'(a)')'      </Attribute>'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      enddo

! Close the File
      call  h5Fclose_f(file_id,error)
! Footer
      write(unitxmf,'(a)')'    </Grid>'
      write(unitxmf,'(a)')'  </Domain>'

      write(unitxmf,'(a)')'</Xdmf>'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      close(unitxmf)
      return
      end subroutine Output

      subroutine hdf5OUT(file_id,rank,dims,dsetname,dset)
      use hdf5
      implicit none
      integer(hid_t),intent(in) :: file_id                    ! file identifier
      integer,intent(in)  :: rank                             ! dataset rank
      integer(hsize_t), dimension(rank),intent(in)  :: dims   ! dataset dimensions
      character(80),intent(in):: dsetname                     ! dataset name
      real(4),intent(in) :: dset
      integer(hid_t) :: dset_id       ! dataset identifier
      integer(hid_t) :: dspace_id     ! dataspace identifier
      integer :: error
! Create the data space for the data set
      call h5Screate_simple_f(rank, dims, dspace_id, error)
! Get dset_id for data set
      call h5Dcreate_f(file_id,dsetname,h5t_native_real,dspace_id, &
     &                 dset_id,error)
! Write the data
      call h5Dwrite_f(dset_id, h5t_native_real, dset, dims, error)
! Close the data id
      call h5Dclose_f(dset_id, error)
! Close the data space
      call h5Sclose_f(dspace_id, error)

      return
end subroutine hdf5out

subroutine makedirs(outdir)
  character(len=*), intent(in) :: outdir
  character(len=256) command
  write(command, *) 'if [ ! -d ', trim(outdir), ' ]; then mkdir -p ', trim(outdir), '; fi'
  write(*, *) trim(command)
  call system(command)
end subroutine makedirs
