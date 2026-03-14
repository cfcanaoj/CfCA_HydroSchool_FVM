module params

real(8), parameter:: timemax=0.15d0 ! simulation end time

! option
integer, parameter :: flag_HDC = 1 ! 1 --> HDC on , 0 --> HDC off
integer, parameter :: flag_flux = 2 ! 1 (HLL), 2 (HLLD)

! coordinate 
integer,parameter::nxtotal=1280 ! the number of grids in the simulation box
integer,parameter::nytotal=1280 ! the number of grids in the simulation box

integer,parameter::nxdiv=2 ! the number of grids in the simulation box
integer,parameter::nydiv=1 ! the number of grids in the simulation box


integer,parameter::nx=nxtotal/nxdiv ! the number of grids in the simulation box
integer,parameter::ny=nytotal/nydiv ! the number of grids in the simulation box
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
