subroutine avanzamento

!****************************************************************!
! COMPUTE: time advance of the equations                         ! 
!          [2D version]                                          !
! -------------------------------------------------------------- ! 
! original version:             F. Califano, 2006                !
! MPI parallel version:         M. Faganello/F. Valentini, 2008  !                  
! Anisotropic/FLR-MHD version:  S. S. Cerri, 2011                !
!****************************************************************!

use parameter_mod
use box_mod
use fields_UJ_mod
use fields_DP_mod
use fields_EB_mod
use dom_distr_mod
use parallel_mod

IMPLICIT NONE

INTEGER :: ix, iy, iz, ip,i

REAL(dp), ALLOCATABLE :: zx(:)
REAL(dp), ALLOCATABLE :: RHS_x(:,:,:), RHS_y(:,:,:), Den_intermedio(:,:,:)
REAL(dp), ALLOCATABLE :: pi_para_int(:,:,:), pi_perp_int(:,:,:)
REAL(dp), ALLOCATABLE :: pe_para_int(:,:,:), pe_perp_int(:,:,:)
REAL(dp), ALLOCATABLE :: trt1(:,:,:), trt2(:,:,:), trt3(:,:,:)
REAL(dp), ALLOCATABLE :: pi_para_n(:,:,:), pi_perp_n(:,:,:)
REAL(dp), ALLOCATABLE :: pe_para_n(:,:,:), pe_perp_n(:,:,:)
REAL(dp), ALLOCATABLE :: Tracciante_n(:,:,:)
REAL(dp), ALLOCATABLE :: Den_n(:,:,:)
REAL(dp), ALLOCATABLE :: nU_z_n(:,:,:), nU_x_n(:,:,:), nU_y_n(:,:,:)
REAL(dp), ALLOCATABLE :: Bx_n(:,:,:), By_n(:,:,:), Bz_n(:,:,:)

allocate( zx(nxl) ) 
allocate( RHS_x(nxl,ny,nz) ) !when needed, RHS_x is used instead of RHS_z to save memory 
allocate( RHS_y(nxl,ny,nz) )
allocate( Den_intermedio(nxl,ny,nz) )
allocate( pi_para_int(nxl,ny,nz) )
allocate( pi_perp_int(nxl,ny,nz) )
allocate( pe_para_int(nxl,ny,nz) )
allocate( pe_perp_int(nxl,ny,nz) )
allocate( trt1(nx,nyt,nz) ) 
allocate( trt2(nx,nyt,nz) )
allocate( trt3(nx,nyt,nz) )
allocate( pi_para_n(nxl,nyl,nzl) ) 
allocate( pi_perp_n(nxl,nyl,nzl) )
allocate( pe_para_n(nxl,nyl,nzl) ) 
allocate( pe_perp_n(nxl,nyl,nzl) )
allocate( Tracciante_n(nxl,nyl,nzl) )
allocate( Den_n(nxl,nyl,nzl))
allocate( nU_z_n(nxl,nyl,nzl) )  
allocate( nU_x_n(nxl,nyl,nzl) ) 
allocate( nU_y_n(nxl,nyl,nzl) ) 
allocate( Bx_n(nxl,nyl,nzl) ) 
allocate( By_n(nxl,nyl,nzl) ) 
allocate( Bz_n(nxl,nyl,nzl) ) 


!********************************************************!
! open boundary condition based on MHD characteristics   !  
! (see: Faganello et al., New J. Phys. 11, 063008 (2009) !
!********************************************************!
!
IF (ibc .EQ. 4) THEN
  call BC
ENDIF

!support/auxiliary arrays
pi_para_n = pi_para
pi_perp_n = pi_perp
pe_para_n = pe_para
pe_perp_n = pe_perp
Tracciante_n = Tracciante
Den_n = Den
nU_z_n = nU_z
nU_x_n = nU_x
nU_y_n = nU_y
Bx_n = Bx
By_n = By
Bz_n = Bz


!************************!
!*** Runge-Kutta loop ***!
!************************!

DO i = rk_ord,1,-1  !RK loop: begins


  !******************************!
  !  Compute: gyroviscous tensor  !
  !******************************!
  !
  ! note: update gyroviscous tensor only if you have FLR_ON = 1
  ! (otherwise, they are kept = 0 from initial condition)
  !
  if (flr_on .eq. 1) then
    call flr_i
  endif
  

  !*** LANDAU-FLUID closures: not available in this release
  qi_para = 0.0
  qi_perp = 0.0
  qe_para = 0.0
  qe_perp = 0.0


  !*********************************************************************!
  ! Advancing: ion-pressure equations                                   !
  ! ------------------------------------------------------------------- ! 
  ! (pressures are updated after the momentum equation is advanced too) ! 
  !*********************************************************************!

  call pressure_i(RHS_x,RHS_y) !RHS_x per pi_para e RHS_y per pi_perp

  pi_para_int = pi_para_n + RHS_x * dt / float(i)
  pi_perp_int = pi_perp_n + RHS_y * dt / float(i)

  !--apply filters to ion pressures
  !
  !along x (open)
  call traspdist(pi_para_int,trt1,1)
  call traspdist(pi_perp_int,trt2,1)
  do iz = 1, nz 
    do iy = 1, nyt
      call filtro_open_x(trt1(:,iy,iz))
      call filtro_open_x(trt2(:,iy,iz))
    enddo
  enddo
  call traspdist(pi_para_int,trt1,-1)
  call traspdist(pi_perp_int,trt2,-1)
  !
  !along y (periodic)
  do iz = 1, nz 
    do ix = 1, nxl
      call filtro_per_y(pi_para_int(ix,:,iz))
      call filtro_per_y(pi_perp_int(ix,:,iz))
    enddo
  enddo
  !
  !along z (periodic) 
  ![for 3D version: not released yet]
  IF (nz.GT.1) THEN
    do iy = 1, ny 
      do ix = 1, nxl
        call filtro_per_z(pi_para_int(ix,iy,:))
        call filtro_per_z(pi_perp_int(ix,iy,:))
      enddo
    enddo
  ENDIF


  !*********************************************************************!
  ! Advancing: electron-pressure equations                              !
  ! ------------------------------------------------------------------- ! 
  ! (pressures are updated after the momentum equation is advanced too) ! 
  !*********************************************************************!

  call pressure_e(RHS_x,RHS_y) !RHS_x per pe_para e RHS_y per pe_perp

  pe_para_int = pe_para_n + RHS_x * dt / float(i)
  pe_perp_int = pe_perp_n + RHS_y * dt / float(i)

  !--apply filters to electron pressures
  !
  !along x (open)
  call traspdist(pe_para_int,trt1,1)
  call traspdist(pe_perp_int,trt2,1)
  do iz = 1, nz 
    do iy = 1, nyt
      call filtro_open_x(trt1(:,iy,iz))
      call filtro_open_x(trt2(:,iy,iz))
    enddo
  enddo
  call traspdist(pe_para_int,trt1,-1)
  call traspdist(pe_perp_int,trt2,-1)
  !
  !along y (periodic)
  do iz = 1, nz 
    do ix = 1, nxl
      call filtro_per_y(pe_para_int(ix,:,iz))
      call filtro_per_y(pe_perp_int(ix,:,iz))
    enddo
  enddo
  !
  !along z (periodic) 
  ![for 3D version: not released yet]
  IF (nz.GT.1) THEN
    do iy = 1, ny 
      do ix = 1, nxl
        call filtro_per_z(pe_para_int(ix,iy,:))
        call filtro_per_z(pe_perp_int(ix,iy,:))
      enddo
    enddo
  ENDIF


  !************************************!
  ! Advancing: passive-tracer equation !
  !************************************!

  call trac(RHS_x)

  Tracciante = Tracciante_n + RHS_x * dt / float(i)

  !--apply filters to passive tracer 
  !
  !along x (open)
  call traspdist(Tracciante,trt1,1)
  do iz = 1, nz 
    do iy = 1, nyt
      call filtro_open_x(trt1(:,iy,iz))
    enddo
  enddo
  call traspdist(Tracciante, trt1,-1)
  !
  !along y (periodic)
  do iz = 1, nz 
    do ix = 1, nxl
      call filtro_per_y(Tracciante(ix,:,iz))
    enddo
  enddo
  !
  !along z (periodic) 
  ![for 3D version: not released yet]
  IF (nz.GT.1) THEN
    do iy = 1, ny 
      do ix = 1, nxl
        call filtro_per_z(Tracciante(ix,iy,:))
      enddo
    enddo
  ENDIF


  !******************************************************************!
  ! Advancing: continuity equation                                   !
  ! ---------------------------------------------------------------- ! 
  ! (density is updated after the momentum equation is advanced too) ! 
  !******************************************************************!

  call cont_n_U(RHS_x)

  Den_intermedio = Den_n + RHS_x * dt / float(i)

  !--apply filters to density 
  !
  !along x (open)
  call traspdist(Den_intermedio,trt1,1)
  do iz = 1, nz 
    do iy = 1, nyt
      call filtro_open_x(trt1(:,iy,iz))
    enddo
  enddo
  call traspdist(Den_intermedio, trt1,-1)
  !
  !along y (periodic)
  do iz = 1, nz 
    do ix = 1, nxl
      call filtro_per_y(Den_intermedio(ix,:,iz))
    enddo
  enddo
  !
  !along z (periodic) 
  ![for 3D version: not released yet]
  IF (nz.GT.1) THEN
    do iy = 1, ny 
      do ix = 1, nxl
        call filtro_per_z(Den_intermedio(ix,iy,:))
      enddo
    enddo
  ENDIF


  !****************************************!
  ! Advancing: one-fluid momentum equation !
  !****************************************!
  
  !--momentum equation along z
  call moto_z(RHS_x)  ! using RHS_x rather than RHS_z to save memory 

  nU_z = nU_z_n + RHS_x * dt / float(i)

  !--momentum equation in (x,y) plane
  call moto_xy(RHS_x, RHS_y)

  nU_x = nU_x_n + RHS_x * dt / float(i)
  nU_y = nU_y_n + RHS_y * dt / float(i)

  !--apply filters to one-fluid momentum 
  !
  !along x (open)
  call traspdist(nU_x,trt1,1)
  call traspdist(nU_y,trt2,1)
  call traspdist(nU_z,trt3,1)
  do iz = 1, nz 
    do iy = 1, nyt
      call filtro_open_x(trt1(:,iy,iz))
      call filtro_open_x(trt2(:,iy,iz))       
      call filtro_open_x(trt3(:,iy,iz))
    enddo
  enddo
  call traspdist(nU_x,trt1,-1)
  call traspdist(nU_y,trt2,-1)
  call traspdist(nU_z,trt3,-1)
  !
  !along y (periodic)
  do iz = 1, nz 
    do ix = 1, nxl
      call filtro_per_y(nU_x(ix,:,iz))
      call filtro_per_y(nU_y(ix,:,iz))
      call filtro_per_y(nU_z(ix,:,iz))
    enddo
  enddo
  !
  !along z (periodic) 
  ![for 3D version: not released yet]  
  IF (nz.GT.1) THEN
    do iy = 1, ny 
      do ix = 1, nxl
        call filtro_per_z(nU_x(ix,iy,:))
        call filtro_per_z(nU_y(ix,iy,:))
        call filtro_per_z(nU_z(ix,iy,:))
      enddo
    enddo
  ENDIF


  !*****************************!
  ! Advancing: Faraday equation !
  !*****************************!

  !--Faraday equation along x
  call faraday_x(RHS_x)

  Bx = Bx_n + RHS_x * dt / float(i)

  !--Faraday equation along y
  call faraday_y(RHS_y)

  By = By_n + RHS_y * dt / float(i)

  !--Faraday equation along z
  call faraday_z(RHS_x) ! using RHS_x rather than RHS_z to save memory

  Bz = Bz_n + RHS_x * dt / float(i)

  !--apply filters to B 
  !
  !along x (open)
  call traspdist(Bx,trt1,1)
  call traspdist(By,trt2,1)
  call traspdist(Bz,trt3,1)
  do iz = 1, nz 
    do iy = 1, nyt
      call filtro_open_x(trt1(:,iy,iz))
      call filtro_open_x(trt2(:,iy,iz))
      call filtro_open_x(trt3(:,iy,iz))
    enddo
  enddo
  call traspdist(Bx,trt1,-1)
  call traspdist(By,trt2,-1)
  call traspdist(Bz,trt3,-1)
  !
  !along y (periodic)
  do iz = 1, nz 
    do ix = 1, nxl
      call filtro_per_y(Bx(ix,:,iz))
      call filtro_per_y(By(ix,:,iz))
      call filtro_per_y(Bz(ix,:,iz))
    enddo
  enddo
  !
  !along z (periodic) 
  ![for 3D version: not released yet]  
  IF (nz.GT.1) THEN
    do iy = 1, ny 
      do ix = 1, nxl
        call filtro_per_z(Bx(ix,iy,:))
        call filtro_per_z(By(ix,iy,:))
        call filtro_per_z(Bz(ix,iy,:))
      enddo
    enddo
  ENDIF


  !****************!
  ! Update density !
  !****************!

  Den = Den_intermedio
  Dinv = 1.0d0 / Den


  !**************************************************!
  ! Compute one-fluid velocity U                     !
  ! ------------------------------------------------ !
  ! (from one-fluid momentum, using updated density) ! 
  !**************************************************!

  Ux = nU_x * Dinv
  Uy = nU_y * Dinv
  Uz = nU_z * Dinv


  !***********************************!
  ! Update pressures (with smoothing) ! 
  !***********************************!
 
  !--"effective" isotropic pressures at x-boundaries ( = 1/3 Tr[P] ) 
  pi_para = terzo * ( 2.0d0 * pi_perp_int + pi_para_int )
  pe_para = terzo * ( 2.0d0 * pe_perp_int + pe_para_int )

  !--smoothing to recover the effective pressures at x-boundaries
  !
  !ions
  do ix = 1, nxl
    pi_para(ix,:,:) = lambda(ixlg + ix - 1) * pi_para_int(ix,:,:) + &
                     ( 1.0d0 - lambda(ixlg + ix - 1) ) * pi_para(ix,:,:)
  enddo
  do ix = 1, nxl
    pi_perp(ix,:,:) = lambda(ixlg + ix - 1) * pi_perp_int(ix,:,:) + &
                     ( 1.0d0 - lambda(ixlg + ix - 1) ) * pi_para(ix,:,:)
  enddo
  !
  !electrons
  do ix = 1, nxl
    pe_para(ix,:,:) = lambda(ixlg + ix - 1) * pe_para_int(ix,:,:) + &
                     ( 1.0d0 - lambda(ixlg + ix - 1) ) * pe_para(ix,:,:)
  enddo
  do ix = 1, nxl
    pe_perp(ix,:,:) = lambda(ixlg + ix - 1) * pe_perp_int(ix,:,:) + &
                     ( 1.0d0 - lambda(ixlg + ix - 1) ) * pe_para(ix,:,:)
  enddo

  !--associated temperatures
  Ti_para = Dinv * pi_para
  Ti_perp = Dinv * pi_perp
  Te_para = Dinv * pe_para
  Te_perp = Dinv * pe_perp


  !*************************************!
  ! Compute current density J = curl(B) ! 
  !*************************************!

  call corrente

  !--apply filters to J 
  !
  !along x (open)
  call traspdist(Jx,trt1,1)
  call traspdist(Jy,trt2,1)
  call traspdist(Jz,trt3,1)
  do iz = 1, nz 
    do iy = 1, nyt
      call filtro_open_x(trt1(:,iy,iz))
      call filtro_open_x(trt2(:,iy,iz))
      call filtro_open_x(trt3(:,iy,iz))
    enddo
  enddo
  call traspdist(Jx,trt1,-1)
  call traspdist(Jy,trt2,-1)
  call traspdist(Jz,trt3,-1)
  !
  !along y (periodic)
  do iz = 1, nz 
    do ix = 1, nxl
      call filtro_per_y(Jx(ix,:,iz))
      call filtro_per_y(Jy(ix,:,iz))
      call filtro_per_y(Jz(ix,:,iz))
    enddo
  enddo
  !
  !along z (periodic) 
  ![for 3D version: not released yet]  
  IF (nz.GT.1) THEN
    do iy = 1, ny 
      do ix = 1, nxl
        call filtro_per_z(Jx(ix,iy,:))
        call filtro_per_z(Jy(ix,iy,:))
        call filtro_per_z(Jz(ix,iy,:))
      enddo
    enddo
  ENDIF


  !************************************!
  ! Computing species' fluid velocitiy ! 
  !************************************!

  call u_ei


  !****************************************************!
  ! Compute: Electric field from generalized Ohm's law !
  !****************************************************!

  call ohm

  !--apply filters to E 
  !
  !along x (open)
  call traspdist(Ex,trt1,1)
  call traspdist(Ey,trt2,1)
  call traspdist(Ez,trt3,1)
  do iz = 1, nz 
    do iy = 1, nyt
      call filtro_open_x(trt1(:,iy,iz))
      call filtro_open_x(trt2(:,iy,iz))
      call filtro_open_x(trt3(:,iy,iz))
    enddo
  enddo
  call traspdist(Ex,trt1,-1)
  call traspdist(Ey,trt2,-1)
  call traspdist(Ez,trt3,-1)
  !
  !along y (periodic)
  do iz = 1, nz 
    do ix = 1, nxl
      call filtro_per_y(Ex(ix,:,iz))
      call filtro_per_y(Ey(ix,:,iz))
      call filtro_per_y(Ez(ix,:,iz))
    enddo
  enddo
  !
  !along z (periodic) 
  ![for 3D version: not released yet]  
  IF (nz.GT.1) THEN
    do iy = 1, ny 
      do ix = 1, nxl
        call filtro_per_z(Ex(ix,iy,:))
        call filtro_per_z(Ey(ix,iy,:))
        call filtro_per_z(Ez(ix,iy,:))
      enddo
    enddo
  ENDIF


ENDDO !RK loop: ends


deallocate(Den_intermedio)
deallocate(pi_para_int)
deallocate(pi_perp_int)
deallocate(pe_para_int)
deallocate(pe_perp_int)
deallocate(pi_para_n) 
deallocate(pi_perp_n)
deallocate(pe_para_n) 
deallocate(pe_perp_n)
deallocate(Tracciante_n)
deallocate(Den_n)
deallocate(nU_z_n) 
deallocate(nU_y_n) 
deallocate(nU_x_n)
deallocate(Bx_n) 
deallocate(By_n) 
deallocate(Bz_n)
deallocate(zx) 
deallocate(RHS_x) 
deallocate(RHS_y)
deallocate(trt1)
deallocate(trt2)
deallocate(trt3)

end subroutine

