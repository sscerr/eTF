subroutine inizializza 

!--------------------------------------------------- !
!  F. Califano, 2006                                 !
!--------------------------------------------------- !

!**************************************************
!  MPI PARALLEL VERSION: VALENTINI-FAGANELLO 2009
!   3D PARALLEL VERSION: FAGANELLO 2010
!**************************************************

use box_mod
use deriv_mod
use parameter_mod
use poisson_open
use dom_distr_mod 

IMPLICIT NONE

integer  :: iy, iz, ily, ilz

real(dp) :: zz, c1, p_a0, p_a1, p_b0, p_b1, p_a12
real(dp) :: alpha_1_3x3, alpha_1_3x5, alpha_1_bc

! --------------------------------- !
! *** Set up derivate Fourier       !
! --------------------------------- !

     call drffti( ny, savey )
     call drffti( nz, savez )

!---------------------------------------------- !
!   Fourier wave vector for y and z-derivative  !
!---------------------------------------------- !

     do iy = 2, ny1, 2
        ky_1(iy/2) = iy * hyl / (2.0 * ny)
        ky_4(iy/2) = iy * hyl / (2.0)
     enddo 

      ky_4 = ky_4**4 / ny

     do iz = 2, nz1, 2
        kz_1(iz/2) = iz * hzl / (2.0 * nz)
        kz_4(iz/2) = iz * hzl / (2.0)
     enddo 

      kz_4 = kz_4**4 / nz

! --------------------------------- !
! *** Set up derivate Compatte      !
! *** Parametri solver tridiagonale !
! --------------------------------- !

    TRANS = 'N'
    info  = 0
    nrhs  = 1
    ndb   = nx

!    L.H.S. coefficients of implicit I

        alpha_1_3x3 = 1.0d0 / 4.0d0
        alpha_1_3x5 = 1.0d0 / 3.0d0

!    R.H.S. coefficients of implicit I derivative.

        a_1_3x3 = 3.0d0 / (4.0  * dx)
        a_1_3x5 = 7.0d0 / (9.0  * dx)
        b_1_3x5 = 1.0d0 / (36.0 * dx)

!    L.H.S coefficients for FREE boundary conditions, I derivative.

        alpha_1_bc = 3.0d0

!    R.H.S coefficients boundary conditions of I derivative.

        a_1_bc = - (11.0d0 + 2.0d0 * alpha_1_bc) / (6.0d0 * dx)
        b_1_bc = (6.0d0 - alpha_1_bc) / (2.0d0 * dx)
        c_1_bc = (2.0d0 * alpha_1_bc - 3.0d0) / (2.0d0 * dx)
        d_1_bc = (2.0d0 - alpha_1_bc) / (6.0d0 * dx)

!        Differenze finite compatte 3x5
!        R.H.S arrays for LU decomposition, I derivative
!        TRE CASI, d.., w.., q.. a seconda delle cond. contorno

!             Diagonale

        dd_1 = 1.0
        wd_1 = 1.0
        qd_1 = 1.0

!        dm_1 e dp_1, wm_1 e wp_1, qm_1 e qp_1 sono 
!        sotto e sopra diagonale, dim. nx - 1

        dm_1 = alpha_1_3x5
        dp_1 = alpha_1_3x5

        wm_1 = alpha_1_3x5
        wp_1 = alpha_1_3x5

        qm_1 = alpha_1_3x5
        qp_1 = alpha_1_3x5

!    Differenze finite compatte 3x3 in i = 2, nx-1

        dm_1(1)   = alpha_1_3x3
        dp_1(2)   = alpha_1_3x3
        dm_1(nx2) = alpha_1_3x3
        dp_1(nx1) = alpha_1_3x3

        wm_1(1)   = alpha_1_3x3
        wp_1(2)   = alpha_1_3x3
        wm_1(nx2) = alpha_1_3x3
        wp_1(nx1) = alpha_1_3x3

        qm_1(1)   = alpha_1_3x3
        qp_1(2)   = alpha_1_3x3
        qm_1(nx2) = alpha_1_3x3
        qp_1(nx1) = alpha_1_3x3

!     Working arrays

        dw_1 = 0.0
        ww_1 = 0.0
        qq_1 = 0.0

!!! **** CONDIZIONI AL CONTORNO *** !!!

!     Boundary condition calculated with internal points
!        (left and right), I derivative, i=1, i=nx

        dp_1(1)   = alpha_1_bc
        dm_1(nx1) = alpha_1_bc

!     FREE sleep boundary condition (left and right)
!       d/dx = 0 ai bordi

        wp_1(1)   = 0.0
        wm_1(nx1) = 0.0

!     FREE^2 sleep boundary condition (left and right)
!       d/dx = cste ai bordi (f'_1 - f'_2 = 0 e idem a destra)

        qp_1(1)   = - 1.0
        qm_1(nx1) = - 1.0

! DERIVATA SENZA COND. AL CONTORNO (si usano i punti interni al bordo)
!     Fattorizzazione LU derivata I

    CALL DGTTRF(nx, dm_1, dd_1, dp_1, dw_1, ipv_d, info)

      if (info > 0 .or. info < 0) then
      write(*,*) 'Problemi fattorizzazione LU der I, info:', info
      stop
      endif

! DERIVATA CON COND. AL CONTORNO di derivata nulla
!     Fattorizzazione LU derivata I

    CALL DGTTRF(nx, wm_1, wd_1, wp_1, ww_1, ipv_w, info)

      if (info > 0 .or. info < 0) then
      write(*,*) 'Problemi fattorizzazione LU der I, info:', info
      stop
      endif

! DERIVATA CON COND. AL CONTORNO di derivata costante 
!     Fattorizzazione LU derivata I

    CALL DGTTRF(nx, qm_1, qd_1, qp_1, qq_1, ipv_q, info)

      if (info > 0 .or. info < 0) then
      write(*,*) 'Problemi fattorizzazione LU der I, info:', info
      stop
      endif



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  SOLVER DI POISSON direzione con cond. cont. open

! ------------------------------------------------------------ !
! *** Set up coefficienti derivate Compatte per Poisson solver !
! *** Parametri solver tridiagonale                            !
! ------------------------------------------------------------ !

    TRNAS = 'N'
    ifno  = 0
    irhs  = 1
    idb   = nx

do ilz = 0, n2z-1
 do ily = 0, ny/2-1

! ===================================================================== !
! coefficients for the EQUATION:                                        !
!                                                                       !
! E_i + (a1/a0) (E_i-1 + E_i+1) = (b1/a0) (f_i-1 + f_i+1) + (b0/a0) f_i !
!                                                                       !
!         viene da (1 + de^2 nabla^2) E = RHS                           !
!     vedi Twofluids_equation                                           ! 
! ===================================================================== !

 c1 = rapm + ily * ily * hyl * hyl + ilz * ilz * hzl * hzl     ! rapm = 1 / de^2

p_a0 = 10.0 * c1 * dx * dx + 24.0
p_a1 =        c1 * dx * dx - 12.0
p_b0 = 10.0 * dx * dx * rapm
p_b1 =        dx * dx * rapm

! ======================================= !
! R.H.S. coefficients                     !
! ======================================= !

pss_rhs_0(ily,ilz) = p_b0 / p_a0
pss_rhs_1(ily,ilz) = p_b1 / p_a0

! ======================================= !
! L.H.S. coefficients                     !
! ======================================= !

p_a12 = p_a1 / p_a0

  poiss_d(:,ily,ilz) = 1.0
  poiss_m(:,ily,ilz) = p_a12
  poiss_p(:,ily,ilz) = p_a12

  poiss_w = 0.d0

! condizioni al contorno: E = cost a sinistra e a destra per tutti i modi

  poiss_p(1,ily,ilz)   = 0.0d0
  poiss_m(nx1,ily,ilz) = 0.0d0

! condizioni al contorno: dE/dx = 0.0 a sinistra e a destra per il modo ky=0

!  poiss_p(1,il)   = - 1.0d0
!  poiss_m(nx1,il) = - 1.0d0


 CALL DGTTRF(nx, poiss_m(:,ily,ilz), poiss_d(:,ily,ilz), poiss_p(:,ily,ilz), poiss_w, ipv_p(:,ily,ilz), ifno)

      if (ifno > 0 .or. ifno < 0) then
      write(*,*) 'Problemi fattorizzazione LU, info:', info
      stop
      endif

 enddo
enddo

end subroutine
