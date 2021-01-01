MODULE parameter_mod 

  IMPLICIT NONE
  SAVE

  INTEGER, PARAMETER :: dp = selected_real_kind(14,300)

  INTEGER, PARAMETER :: nx=4096, ny=2048, nz=1
  INTEGER, PARAMETER :: nx1=nx-1, ny1=ny-1, nz1=nz-1
  INTEGER, PARAMETER :: nx2=nx-2
  INTEGER, PARAMETER :: nx3=nx-3
  INTEGER, PARAMETER :: nx4=nx-4
  INTEGER, PARAMETER :: n2x=nx/2, n2y=ny/2
! n2z=nz/2 computed in pstartup (so to switch 2D/3D cases)

  INTEGER, PARAMETER :: l_rst=500, l_cpu=501, l_div=502, l_Eng=503 , &
                        l_U=504, l_DPJ=505, l_EB=506, l_Gi=507

  INTEGER, PARAMETER :: nprocx = 20000 
  INTEGER :: l_EB_me, l_U_me, l_DPJ_me, l_Gi_me, l_rst_me, n2z
  REAL(dp)    :: pgreco, ninvy, ninvz

  CHARACTER*4 :: prefix

  CHARACTER(LEN=2), PARAMETER :: space_dim = '2D'

END MODULE parameter_mod

MODULE box_mod

  use parameter_mod, only: dp, nx, ny, nz

  IMPLICIT NONE
  SAVE

  integer :: istart, ioutt, ioutx, nwx, nwy, nwz, ibc
  integer :: rk_ord, nx_lmbd, hall_on, divPe_on, flr_on

  REAL(dp)  :: rapm, de, de2, de12inv, gam, gm1, eta, ampl, alpha
  REAL(dp)  :: dt, ab1, ab2, ab3
  REAL(dp)  :: xl, yl, zl, dx, dy, dz, dxyz, hyl, hzl
  REAL(dp)  :: Tot_inv_ny, Tot_inv_nz, Tot_inv_nxnynz
  REAL(dp)  :: tempo, tmax, tt_last, tt_w, tx_last, tx_w
  REAL(dp)  :: nx_lambda, terzo, grdlmbd

  REAL(dp),DIMENSION (nx) :: x
  REAL(dp),DIMENSION (ny) :: y
  REAL(dp),DIMENSION (nz) :: z
  REAL(dp),DIMENSION (nx) :: lambda

END MODULE box_mod

MODULE fields_UJ_mod

  use parameter_mod, only: dp, nx, ny, nz

  IMPLICIT NONE
  SAVE

  REAL(dp) :: fac1_i, fac2_i, fac3_i
  REAL(dp) :: fac1_e, fac2_e, fac3_e

  REAL(dp),ALLOCATABLE, DIMENSION (:,:,:) :: Uex, Uey, Uez, Uix, Uiy, Uiz
  REAL(dp),ALLOCATABLE, DIMENSION (:,:,:) :: Ux, Uy, Uz, Jx, Jy, Jz
  REAL(dp),ALLOCATABLE, DIMENSION (:,:,:) :: Tracciante
  REAL(dp),ALLOCATABLE, DIMENSION (:,:,:) :: Ti_para, Ti_perp, Te_para, Te_perp

CONTAINS

  SUBROUTINE allocate_fields_UJ( nxl, nyl, nzl )
     INTEGER :: nxl, nyl, nzl
     ALLOCATE( Uex( nxl, nyl, nzl ) )
     ALLOCATE( Uey( nxl, nyl, nzl ) )
     ALLOCATE( Uez( nxl, nyl, nzl ) )
     ALLOCATE( Uix( nxl, nyl, nzl ) )
     ALLOCATE( Uiy( nxl, nyl, nzl ) )
     ALLOCATE( Uiz( nxl, nyl, nzl ) )
     ALLOCATE( Ux( nxl, nyl, nzl ) )
     ALLOCATE( Uy( nxl, nyl, nzl ) )
     ALLOCATE( Uz( nxl, nyl, nzl ) )
     ALLOCATE( Jx( nxl, nyl, nzl ) )
     ALLOCATE( Jy( nxl, nyl, nzl ) )
     ALLOCATE( Jz( nxl, nyl, nzl ) )
     ALLOCATE( Tracciante( nxl, nyl, nzl ) )
     ALLOCATE( Ti_para( nxl, nyl, nzl ) )
     ALLOCATE( Ti_perp( nxl, nyl, nzl ) )
     ALLOCATE( Te_para( nxl, nyl, nzl ) )
     ALLOCATE( Te_perp( nxl, nyl, nzl ) )
     RETURN
  END SUBROUTINE

END MODULE fields_UJ_mod

MODULE fields_DP_mod

  use parameter_mod, only: dp, nx, ny, nz

  IMPLICIT NONE
  SAVE

  REAL(dp),ALLOCATABLE, DIMENSION (:,:,:) :: ax, ay, az, a_inv, c_s, c_s_inv
  REAL(dp),ALLOCATABLE, DIMENSION (:,:,:) :: s, f, f_inv
  REAL(dp),ALLOCATABLE, DIMENSION (:,:,:) :: alpha_1  , alpha_2  , segno   
  REAL(dp),ALLOCATABLE, DIMENSION (:,:,:) :: L_s_m, L_a_m, L_f_m
  REAL(dp),ALLOCATABLE, DIMENSION (:,:,:) :: L_s_p, L_a_p, L_f_p
  REAL(dp),ALLOCATABLE, DIMENSION (:,:,:) :: L_0
  REAL(dp),ALLOCATABLE, DIMENSION (:,:,:) :: Den, Dinv
  REAL(dp),ALLOCATABLE, DIMENSION (:,:,:) :: pe_para, pe_perp, pi_para, pi_perp
  REAL(dp),ALLOCATABLE, DIMENSION (:,:,:) :: nU_x, nU_y, nU_z
  REAL(dp),ALLOCATABLE, DIMENSION (:,:,:) :: qi_para, qi_perp, qe_para, qe_perp
  REAL(dp),ALLOCATABLE, DIMENSION (:,:,:) :: Gi_xx, Gi_xy, Gi_xz, Gi_yz

CONTAINS

  SUBROUTINE allocate_fields_DP( nxl, nyl, nzl )
     INTEGER :: nxl, nyl, nzl
     ALLOCATE( ax( 2, nyl, nzl ) )
     ALLOCATE( ay( 2, nyl, nzl ) )
     ALLOCATE( az( 2, nyl, nzl ) )
     ALLOCATE( a_inv( 2, nyl, nzl ) )
     ALLOCATE( c_s( 2, nyl, nzl ) )
     ALLOCATE( c_s_inv( 2, nyl, nzl ) )
     ALLOCATE( s( 2, nyl, nzl ) )
     ALLOCATE( f( 2, nyl, nzl ) )
     ALLOCATE( f_inv( 2, nyl, nzl ) )
     ALLOCATE( alpha_1( 2, nyl, nzl ) )
     ALLOCATE( alpha_2( 2, nyl, nzl ) )
     ALLOCATE( segno( 2, nyl, nzl ) )
     ALLOCATE( L_s_m( 2, nyl, nzl ) )
     ALLOCATE( L_a_m( 2, nyl, nzl ) )
     ALLOCATE( L_f_m( 2, nyl, nzl ) )
     ALLOCATE( L_s_p( 2, nyl, nzl ) )
     ALLOCATE( L_a_p( 2, nyl, nzl ) )
     ALLOCATE( L_f_p( 2, nyl, nzl ) )
     ALLOCATE( L_0( 2, nyl, nzl ) )
     ALLOCATE( Den( nxl, nyl, nzl ) )
     ALLOCATE( Dinv( nxl, nyl, nzl ) )
     ALLOCATE( pi_para( nxl, nyl, nzl ) )
     ALLOCATE( pi_perp( nxl, nyl, nzl ) )
     ALLOCATE( pe_para( nxl, nyl, nzl ) )
     ALLOCATE( pe_perp( nxl, nyl, nzl ) )
     ALLOCATE( nU_x( nxl, nyl, nzl ) )
     ALLOCATE( nU_y( nxl, nyl, nzl ) )
     ALLOCATE( nU_z( nxl, nyl, nzl ) )
     ALLOCATE( qi_para( nxl, nyl, nzl ) )
     ALLOCATE( qi_perp( nxl, nyl, nzl ) )
     ALLOCATE( qe_para( nxl, nyl, nzl ) )
     ALLOCATE( qe_perp( nxl, nyl, nzl ) )
     ALLOCATE( Gi_xx( nxl, nyl, nzl ) )
     ALLOCATE( Gi_xy( nxl, nyl, nzl ) )
     ALLOCATE( Gi_xz( nxl, nyl, nzl ) )
     ALLOCATE( Gi_yz( nxl, nyl, nzl ) )
     RETURN
  END SUBROUTINE
  

END MODULE fields_DP_mod


MODULE fields_EB_mod

  use parameter_mod, only: dp, nx, ny, nz

  IMPLICIT NONE
  SAVE

   !REAL(dp) :: Bx0, By0, Bz0
   REAL(dp) :: B00, angolo, Aeq, Leq, DeltaDen

  REAL(dp), ALLOCATABLE, DIMENSION (:,:,:) :: Ex, Ey, Ez, Bx, By, Bz

CONTAINS

  SUBROUTINE allocate_fields_EB( nxl, nyl, nzl )
    INTEGER :: nxl, nyl, nzl
    ALLOCATE( Ex( nxl, nyl, nzl ) )
    ALLOCATE( Ey( nxl, nyl, nzl ) )
    ALLOCATE( Ez( nxl, nyl, nzl ) )
    ALLOCATE( Bx( nxl, nyl, nzl ) )
    ALLOCATE( By( nxl, nyl, nzl ) )
    ALLOCATE( Bz( nxl, nyl, nzl ) )
    RETURN
  END SUBROUTINE

END MODULE fields_EB_mod

MODULE deriv_mod

  use parameter_mod, only: dp, nx, ny, nz, nx1, nx2, n2z

  IMPLICIT NONE
  SAVE

! Derivate open (lungo x)

CHARACTER*1 TRANS

integer :: nrhs, ndb, info

integer :: ipv_d(nx), ipv_w(nx), ipv_q(nx)

real(dp) :: a_1_3x3, a_1_3x5, b_1_3x5
real(dp) :: a_1_bc, b_1_bc, c_1_bc, d_1_bc

REAL(dp),DIMENSION (nx)  :: dd_1
REAL(dp),DIMENSION (nx1) :: dm_1, dp_1
REAL(dp),DIMENSION (nx2) :: dw_1

REAL(dp),DIMENSION (nx)  :: wd_1
REAL(dp),DIMENSION (nx1) :: wm_1, wp_1
REAL(dp),DIMENSION (nx2) :: ww_1

REAL(dp),DIMENSION (nx)  :: qd_1
REAL(dp),DIMENSION (nx1) :: qm_1, qp_1
REAL(dp),DIMENSION (nx2) :: qq_1

! Derivate periodiche (lungo y,z)

REAL(dp):: ky_1(ny/2 - 1), ky_4(ny/2 - 1), savey(2 * ny + 15)
REAL(dp):: savez(2 * nz + 15)
REAL(dp), ALLOCATABLE, DIMENSION (:) :: kz_1, kz_4


END MODULE deriv_mod

MODULE filtro_per_mod

  use parameter_mod, only: dp, nx, ny, nz

  IMPLICIT NONE
  SAVE

real(dp) :: omega_x, omega_y, omega_z

!perche'? omega_x def anche qui????? refuso double-periodic??

REAL(dp):: work4x(nx/2), work4y(ny/2)
REAL(dp), ALLOCATABLE, DIMENSION (:) :: work4z

END MODULE filtro_per_mod

MODULE filtro_open_mod

  use parameter_mod, only: dp, nx, nx1, nx2

  IMPLICIT NONE
  SAVE

CHARACTER*1 TARNS

integer :: lrhs, ldb, info

integer :: ipv_f(nx)

real(dp) :: omega_x
real(dp) :: a_filt, b_filt, c_filt, d_filt
real(dp) :: a11_filt, a12_filt, a13_filt, a14_filt, a15_filt
real(dp) :: a21_filt, a22_filt, a23_filt, a24_filt, a25_filt
real(dp) :: a31_filt, a32_filt, a33_filt, a34_filt, a35_filt

REAL(dp),DIMENSION (nx)  :: filt_d
REAL(dp),DIMENSION (nx1) :: filt_m, filt_p
REAL(dp),DIMENSION (nx2) :: filt_2

END MODULE filtro_open_mod

MODULE poisson_open
  use parameter_mod, only: dp, nx, ny, nz, nx1, nx2, n2y, n2z

  IMPLICIT NONE
  SAVE
CHARACTER*1 TRNAS
integer :: irhs, idb, ifno
REAL(dp),DIMENSION (nx2) :: poiss_w

  integer , ALLOCATABLE, DIMENSION (:,:,:) :: ipv_p 
  REAL(dp), ALLOCATABLE, DIMENSION (:,:)   :: pss_rhs_0, pss_rhs_1
  REAL(dp), ALLOCATABLE, DIMENSION (:,:,:) :: poiss_d
  REAL(dp), ALLOCATABLE, DIMENSION (:,:,:) :: poiss_m, poiss_p

CONTAINS

  SUBROUTINE allocate_poisson( n2z )
    INTEGER :: n2z

    ALLOCATE( ipv_p(nx,0:n2y-1,0:n2z-1) )
    ALLOCATE( pss_rhs_0(0:n2y-1,0:n2z-1) )
    ALLOCATE( pss_rhs_1(0:n2y-1,0:n2z-1) )
    ALLOCATE( poiss_d(nx,0:n2y-1,0:n2z-1) )
    ALLOCATE( poiss_m(nx-1,0:n2y-1,0:n2z-1) )
    ALLOCATE( poiss_p(nx-1,0:n2y-1,0:n2z-1) )

    RETURN
  END SUBROUTINE

END MODULE poisson_open


MODULE dom_distr_mod

  use parameter_mod, only: nprocx

  IMPLICIT NONE
  SAVE

  INTEGER :: nxl, nxg, ixlg, nxlm1, nxlm2, nxlm3, nxlm4 &
           , nyl, nyg, iylg, nylm1, nylm2, nylm3, nylm4 &
           , nyt,      iytg &
           , nzl &
           , nxlmax, nytmax

  INTEGER, DIMENSION (nprocx) :: nxlp, nylp, ixlgp, iylgp, iprow, ipcol
  INTEGER, DIMENSION (nprocx) :: nytp, iytgp


END MODULE dom_distr_mod

MODULE parallel_mod
  IMPLICIT NONE
  SAVE
  INTEGER :: mpime, nproc, root, group, nprow
  INTEGER :: npcol, myrow, mycol
  LOGICAL :: parallel_build
END MODULE

