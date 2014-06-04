subroutine init

! Expansion of the solar wind into vacuum in 2 or 3 dimensions
! 2june2014 don
!=======================================================================

! GLOBALS
use global
use zone
use mpi

IMPLICIT NONE

! LOCALS
INTEGER :: i, j, jp, m, k, mpierr
REAL :: xmin, xmax, ymin, ymax, zmin, zmax
REAL :: rodt, ridt, xvel, yvel, zvel, width, widthz, widthy
REAL :: dleft, pleft, dright, pright, plane

!--------------------------------------------------------------------------------
! Set up geometry and boundary conditions of grid
!
! Boundary condition flags : nleft, nright
!   = 0  :  reflecting boundary condition
!   = 1  :  inflow/outflow boundary condition
!   = 2  :  fixed inflow boundary condition
!   = 3  :  periodic
! Geometry flag : ngeom                         |  Cartesian:
!   = 0  :  planar                              |    gx = 0, gy = 0, gz = 0
!   = 1  :  cylindrical radial                  |  Cylindrical:
!   = 2  :  spherical   radial             3D= {     gx = 1, gy = 3, gz = 0
!   = 3  :  cylindrical angle                   |
!   = 4  :  spherical polar angle (theta)       |  Spherical:
!   = 5  :  spherical azimu angle (phi)         |    gx = 2, gy = 4, gz = 5

! Define the problem...

! don : Setup 2D cylindrical coordinates (r,z) with inflow/outflow boundary
! conditions everywhere.
ngeomx = 0 ! Cylindrical planar z-axis (from Raph)
ngeomy = 1 ! Cylindrical radial coord (from Raph)
ngeomz = 3 ! Cylindrical angle phi

! don : set all boundaries to inflow/outflow (zero gradient)
nleftx = 1
nrightx= 1
nlefty = 1
nrighty= 1
nleftz = 1
nrightz= 1

!don : The sun is centered at (z,r) = (0,0)
xmin   = -10.0*kmperau ! Min. z
xmax   = 100.0*kmperau ! Max. z
ymin   = 0.0*kmperau ! Min. r 
ymax   = 25.0*kmperau ! Max. r
zmin   = 0.0 ! Min. phi
zmax   = 2.0 ! Max. phi: units of pi, doesn't matter since kmax=1

! If any dimension is angular, multiply coordinates by pi...
if(ngeomy > 2) then
   ymin = ymin * pi
   ymax = ymax * pi
endif
if(ngeomz > 2) then
   zmin = zmin * pi
   zmax = zmax * pi
endif

!======================================================================
! Set up parameters from the problem (Expansion into vacuum)
!don : Insert correct parameters in vh1mods.f90
gam    = 4. / 3.
gamm   = gam - 1.0

Rinject = 0.5*kmperau ! Injection radius
vsphereinject = 4.5*sqrt(RT)
pinject = dinject*RT


!=======================================================================
! set time and cycle counters

time   = 0.0
timep  = 0.0
timem  = 0.0
ncycle = 0
ncycp  = 0
ncycd  = 0
ncycm  = 0
nfile  = 1000

! Set up grid coordinates 

call grid(imax,xmin,xmax,zxa,zxc,zdx)
call grid(jmax,ymin,ymax,zya,zyc,zdy)
call grid(kmax,zmin,zmax,zza,zzc,zdz)

if (ndim == 2) zzc(1) = 0.0

!=======================================================================
! Log parameters of problem in history file

if (mype == 0) then
  write (8,*) 'Expansion of solar wind into a vacuum '
  if (ndim == 3) then
   write (8,"(' Grid dimensions: ',i4,' x ',i4,' x ',i4)") imax,jmax,kmax
  else
   write (8,"(' Grid dimensions: ',i4,' x ',i4)") imax,jmax
  endif
  write (8,*) 
  write (8,*) 'Adiabatic index, gamma = ', gam
  write (8,*) 'Initial Background Pressure is ', pinject
  write (8,*) 'Initial Background Density is ', dinject
  write (8,*) 'Spherical radial velocity is ', vsphereinject
  write (8,*) 
endif

! initialize background grid: I'm actually going to have a dead zone and
! injection region, but those will be initialized on each timestep before
! calling the sweep functions.

do k = 1, ks
 do j = 1, js
  do i = 1, imax
    zro(i,j,k) = dinject
    zpr(i,j,k) = pinject
    !zro(i,mypey*js+j,k) = dinject
    !zpr(i,mypey*js+j,k) = pinject
    !don : set the velocity direction (ux = upolarrad, uy = upolarz)
    ! vsphereinject is the spherical radial velocity, zxa is the polar radius
    ! coordinate array, and zya is the polar z coordinate array. Here I find the
    ! polar radial and z velocities.
    zux(i,j,k) = vsphereinject*zxc(i)/(sqrt(zxc(i)**2 + zyc(mypey*js+j)**2))
    zuy(i,j,k) = vsphereinject*zyc(mypey*js+j)/(sqrt(zxc(i)**2 + zyc(mypey*js+j)**2))
    zuz(i,j,k) = 0.0
    zfl(i,j,k) = 0.0
    !zux(i,mypey*js+j,k) = vsphereinject*zxc(i)/(sqrt(zxc(i)**2 + zyc(mypey*js+j)**2))
    !zuy(i,mypey*js+j,k) = vsphereinject*zyc(mypey*js+j)/(sqrt(zxc(i)**2 + zyc(mypey*js+j)**2))
    !zuz(i,mypey*js+j,k) = 0.0
    !zfl(i,mypey*js+j,k) = 0.0
  enddo
 enddo
enddo

call fillmasks


!########################################################################
! Compute Courant-limited timestep

ridt = 0.

if (ndim == 2) then

 do j = 1, js
  do i = 1, imax
    widthy = zdy(j+mypey*js)
    if(ngeomy > 2) widthy = widthy*zxc(i)
    width  = min(zdx(i),widthy)
    svel = sqrt(gam*zpr(i,j,1)/zro(i,j,1))/width
    xvel = abs(zux(i,j,1)) / zdx(i)
    yvel = abs(zuy(i,j,1)) / widthy
    ridt = max(xvel,yvel,svel,ridt)
  enddo
 enddo

else

 do k = 1, ks
  do j = 1, js
   do i = 1, imax
     widthy = zdy(j+mypey*js)
     widthz = zdz(k+mypez*ks)
     if(ngeomy.gt.2) widthy = widthy*zxc(i)
     if(ngeomz.gt.2) widthz = widthz*zxc(i)
     if(ngeomz.eq.5) widthz = widthz*sin(zyc(j+mypey*js))
     width  = min(zdx(i),widthy,widthz)
     svel = sqrt(gam*zpr(i,j,k)/zro(i,j,k))/width
     xvel = abs(zux(i,j,k)) / zdx(i)
     yvel = abs(zuy(i,j,k)) / widthy
     zvel = abs(zuz(i,j,k)) / widthz
     ridt = max(xvel,yvel,zvel,svel,ridt)
   enddo
  enddo
 enddo

endif

call MPI_ALLREDUCE( ridt, rodt, 1, VH1_DATATYPE, MPI_MAX, MPI_COMM_WORLD, mpierr )
dt = courant / rodt

return
end

!#####################################################################

subroutine grid( nzones, xmin, xmax, xa, xc, dx )

! Create grid to cover physical size from xmin to xmax
! number of physical grid zones is nzones
!
! xa(1) is left boundary location - at xmin
! xa(nzones+1) is right boundary location - at xmax
!----------------------------------------------------------------------
IMPLICIT NONE

! LOCALS
integer :: nzones, n
real, dimension(nzones) :: xa, dx, xc
real :: dxfac, xmin, xmax

!=======================================================================

dxfac = (xmax - xmin) / float(nzones)
do n = 1, nzones
  xa(n) = xmin + (n-1)*dxfac
  dx(n) = dxfac
  xc(n) = xa(n) + 0.5*dx(n)
enddo

return
end

!#####################################################################

subroutine fillmasks

! GLOBALS
use global
use zone
use mpi

IMPLICIT NONE

! LOCALS
INTEGER :: i, jp, j, m, k
REAL :: hdradsph, radsph, raddiff

! Expects Rinject to be defined in vh1mods.f90 so it's global.
! Assumes 2-D cylindrical geometry and finds the grid cells at injection.
! The contents of the grid subroutine ensure VH-1 is using constant grid
! separation dx...

injectimin = imax
injectimax = 1
injectjmin = jmax
injectjmax = 1

hdradsph = 2.0*sqrt(zdx(1)**2 + zdy(1)**2)
do k = 1, ks
 do m = 1, npey
 do jp = 1, js
!  j = 6 + jp + js*(m-1)
  j = mypey*js + jp
  do i = 1, imax
    ! Compute spherical radial coord for cell center
    radsph = sqrt(zxc(i)**2 + zyc(j)**2)
    raddiff = radsph-Rinject
    if ((raddiff <= hdradsph) .AND. (raddiff >= -1.0*hdradsph)) then
        ! This is an injection region cell
        injectmask(i,j) = 1
        if (i < injectimin) then
            injectimin = i
        endif
        if (i > injectimax) then
            injectimax = i
        endif
        if (j < injectjmin) then
            injectjmin = j
        endif
        if (j > injectjmax) then
            injectjmax = j
        endif
    else
        ! This is not an injection region cell
        injectmask(i,j) = 0
    endif
    if (raddiff < -1.0*hdradsph) then
        ! This is a dead zone cell
        deadmask(i,j) = 1
    else
        ! This is not a dead zone cell
        deadmask(i,j) = 0
    endif
  enddo
 enddo
enddo
enddo
end
