module zone 
!=======================================================================
! (formerly zone.h) global (3D) data arrays
!======================================================================= 
 
!don : I configure this for 2D cylindrical (r,z) coordinates
!don : i = r, j = z, k = phi. Use imax=256, jmax=1028
 INTEGER, PARAMETER :: imax = 128, jmax = 512, kmax = 1    ! Memory dimensions
 INTEGER, PARAMETER :: pey = 8, pez = 1      ! number of MPI tasks
!              ####### for 2D:  ^^^  IF kmax=1, MUST HAVE pez=1   #############
 INTEGER, PARAMETER :: nvar = 6              ! number of primitive fluid variables

 INTEGER :: isy, isz, js, ks, Ya2abuff_size, Za2abuff_size
 INTEGER :: npe, npey, npez, mype, mypey, mypez  ! # of pes and local pe number
 INTEGER :: MPI_COMM_ROW, MPI_COMM_COL
 INTEGER :: VH1_DATATYPE
 
 INTEGER :: ngeomx, ngeomy, ngeomz       ! XYZ Geometry flag
 INTEGER :: nleftx, nlefty, nleftz       ! XYZ Lower Boundary Condition
 INTEGER :: nrightx,nrighty,nrightz      ! XYZ Upper Boundary Condition

 REAL, DIMENSION(imax,jmax/pey,kmax/pez) :: zro, zpr, zux, zuy, zuz, zfl
 REAL, DIMENSION(imax) :: zxa, zdx, zxc
 REAL, DIMENSION(jmax) :: zya, zdy, zyc
 REAL, DIMENSION(kmax) :: zza, zdz, zzc

 !don : Declare arrays containing indices for injection region mask
 INTEGER, DIMENSION(imax,jmax) :: injectmask
 !don : store boundaries of the injection region within the mask
 INTEGER :: injectimin
 INTEGER :: injectimax
 INTEGER :: injectjmin
 INTEGER :: injectjmax

 !don : Declare arrays containing indices for dead zone mask
 INTEGER, DIMENSION(imax,jmax) :: deadmask
 
 REAL, DIMENSION(nvar,kmax/pez,jmax/pey,imax) :: send1
 REAL, DIMENSION(nvar,kmax/pez,imax/pey,jmax) :: send2
 REAL, DIMENSION(nvar,jmax/pey,kmax/pez,imax) :: send3
 REAL, DIMENSION(nvar,jmax/pey,imax/pez,kmax) :: send4
 REAL, DIMENSION(nvar,kmax/pez,jmax/pey,imax/pey,pey) :: recv1
 REAL, DIMENSION(nvar,kmax/pez,imax/pey,jmax/pey,pey) :: recv2
 REAL, DIMENSION(nvar,jmax/pey,kmax/pez,imax/pez,pez) :: recv3
 REAL, DIMENSION(nvar,jmax/pey,imax/pez,kmax/pez,pez) :: recv4
 
 EQUIVALENCE ( send1(1,1,1,1),   send2(1,1,1,1),   send3(1,1,1,1),   send4(1,1,1,1) )
 EQUIVALENCE ( recv1(1,1,1,1,1), recv2(1,1,1,1,1), recv3(1,1,1,1,1), recv4(1,1,1,1,1) )
  
end module zone

