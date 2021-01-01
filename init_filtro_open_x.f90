! ---------------------------------------------------- !
!      Set up LU diagonalizzazione filtro              !
!             F. Califano, Marzo 2001                  !
! ---------------------------------------------------- !
! A (f_i-1 + f_i+1) + f_i = a_filt f_i +               !
!      b_filt (f_i-1 + f_i+1) + c_filt (f_i-2 + f_i+2) ! 
! ---------------------------------------------------- !

           subroutine init_filtro_open_x

use box_mod
use parameter_mod
use filtro_open_mod

IMPLICIT NONE

! --------------------------------- !
! *** Set up filtro Compatto        !
! *** Parametri solver tridiagonale !
! --------------------------------- !

    TARNS = 'N'
    info = 0
    lrhs  = 1
    ldb   = nx

!    Parametro filtraggio

!    FILTRO IV ordine: R.H.S. coefficients of the tridiag. filter
!    Lele 1992, pag. 40, C.2.4

        !a_filt =   (5.0d0 + 6.0d0 * omega_x) / 8.0d0
        !b_filt =   (1.0d0 + 2.0d0 * omega_x) / 4.0d0
        !c_filt = - (1.0d0 - 2.0d0 * omega_x) / 16.0d0
        !d_filt = 0.0d0

!    FILTRO VI ordine: R.H.S. coefficients of the tridiag. filter
!    Lele 1992, pag. 40, C.2.5

        a_filt = (11.0d0 + 10.0d0 * omega_x) / 16.0d0
        b_filt = (15.0d0 + 34.0d0 * omega_x) / 64.0d0
        c_filt = (- 3.0d0 + 6.0d0 * omega_x) / 32.0d0
        d_filt = (  1.0d0 - 2.0d0 * omega_x) / 64.0d0

!    R.H.S boundary conditions coefficients

        a11_filt =  15.0d0 / 16.0d0
        a12_filt =   1.0d0 / 4.0d0
        a13_filt = - 3.0d0 / 8.0d0
        a14_filt =   1.0d0 / 4.0d0
        a15_filt = - 1.0d0 / 16.0d0

        a21_filt =   1.0d0/ 16.0d0
        a22_filt =   3.0d0 / 4.0d0
        a23_filt =   3.0d0 / 8.0d0
        a24_filt = - 1.0d0 / 4.0d0
        a25_filt =   1.0d0/ 16.0d0

        a31_filt = - 1.0d0/ 16.0d0
        a32_filt =   1.0d0 / 4.0d0
        a33_filt =   5.0d0 / 8.0d0
        a34_filt =   1.0d0 / 4.0d0
        a35_filt = - 1.0d0/ 16.0d0

!    Diagonale filt_d, sotto diagonale filt_m e super diagonale filt_p
!    (filt_2 e' un working array)

        filt_d = 1.0

        filt_p = omega_x
        filt_m = omega_x

        filt_p(1)   = 0.0
        filt_m(nx1) = 0.0

!    Filtro R.H.S boundary conditions coefficients: 
!    f(x) unchanged at i = 1, 2 and the same on the right  OR
!    f(x) Lele boundary conditions at i = 1, 2 and the same on the right

        filt_m(1)   = 0.0
        filt_p(2)   = 0.0
        filt_m(nx2) = 0.0
        filt_p(nx1) = 0.0

!    f(x) unchanged also at i = 3 and the same on the right OR
!    f(x) Lele boundary conditions also at i = 3 and the same on the right

        filt_m(2)   = 0.0
        filt_p(3)   = 0.0
        filt_m(nx3) = 0.0
        filt_p(nx2) = 0.0

!    Working array

        filt_2 = 0.0

!     Fattorizzazione LU

    CALL DGTTRF(nx, filt_m, filt_d, filt_p, filt_2, ipv_f, info)

      if (info > 0 .or. info < 0) then
      write(*,*) 'Problemi fattorizzazione LU, Filtro: info:', info
      stop
      endif

 end subroutine
