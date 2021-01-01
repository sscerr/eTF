! ------------------------------------------------------ !
!             filtro compatto 3x5                        !
!             F. Califano, Marzo 2001                    !
! ------------------------------------------------------ !
! A (f_i-1 + f_i+1) + f_i = a_filt f_i +                 !
!       b_filt (f_i-2 + f_i+2) + c_filt (f_i-1 + f_i+1)  !
! ------------------------------------------------------ !

!**************************************************
!  MPI PARALLEL VERSION: VALENTINI-FAGANELLO 2009
!**************************************************

           subroutine filtro_open_x(QW)

use box_mod
use parameter_mod
use filtro_open_mod

IMPLICIT NONE

  integer :: i
  REAL(dp), dimension (nx) :: QW, TW

!     Soluzione del sistema lineare tridiagonale fattorizzato LU
!     Differenze finite compatte 3x5

     TW = QW

!    Filtro IV ordine, Lele 1992, pag. 40, C.2.4

  !!  do i = 3, nx2
  !!   QW(i) = a_filt * TW(i) + b_filt * (TW(i+1) + TW(i-1)) &
  !!                          + c_filt * (TW(i+2) + TW(i-2))
  !!  enddo

!    Filtro VI ordine, Lele 1992, pag. 40, C.2.5

    do i = 4, nx3
     QW(i) = a_filt * TW(i) + b_filt * (TW(i+1) + TW(i-1)) &
           + c_filt * (TW(i+2) + TW(i-2)) + d_filt * (TW(i+3) + TW(i-3))
    enddo

!   Boundary conditions

     !! QW(1)   = TW(1)
     !! QW(2)   = TW(2)
     !! QW(nx)  = TW(nx)
     !! QW(nx1) = TW(nx1)

!   Boundary conditions VI ordine

     !! QW(3)   = TW(3)
     !! QW(nx2) = TW(nx2)

!   Boundary conditions VI ordine, Lele pag. 41, C.2.11.a,b,c

     QW(1) = a11_filt * TW(1) + a12_filt * TW(2) + a13_filt * TW(3) &
           + a14_filt * TW(4) + a15_filt * TW(5)
     QW(2) = a21_filt * TW(1) + a22_filt * TW(2) + a23_filt * TW(3) &
           + a24_filt * TW(4) + a25_filt * TW(5)
     QW(3) = a31_filt * TW(1) + a32_filt * TW(2) + a33_filt * TW(3) &
           + a34_filt * TW(4) + a35_filt * TW(5)

     QW(nx ) = a11_filt * TW(nx ) + a12_filt * TW(nx1) + a13_filt * TW(nx2) &
             + a14_filt * TW(nx3) + a15_filt * TW(nx4)
     QW(nx1) = a21_filt * TW(nx ) + a22_filt * TW(nx1) + a23_filt * TW(nx2) &
             + a24_filt * TW(nx3) + a25_filt * TW(nx4)
     QW(nx2) = a31_filt * TW(nx ) + a32_filt * TW(nx1) + a33_filt * TW(nx2) &
             + a34_filt * TW(nx3) + a35_filt * TW(nx4)


!   Boundary conditions IV ordine, Lele 1992, pag. 40, C.2.11a,b,c

!    QW(1) = a11_filt * TW(1) + a12_filt * TW(2) + a13_filt * TW(3) &
!                             + a14_filt * TW(4) + a15_filt * TW(5)
!
!    QW(2) = a21_filt * TW(1) + a22_filt * TW(2) + a23_filt * TW(3) &
!                             + a24_filt * TW(4) + a25_filt * TW(5)
!
!    QW(nx1) = a21_filt * TW(nx) + a22_filt * TW(nx1) &
!     + a23_filt * TW(nx2) + a24_filt * TW(nx-3) + a25_filt * TW(nx-4)
!
!    QW(nx)  = a11_filt * TW(nx) + a12_filt * TW(nx1) &
!     + a13_filt * TW(nx2) + a14_filt * TW(nx-3) + a15_filt * TW(nx-4)


 CALL  DGTTRS(TARNS,nx,LRHS,filt_m,filt_d,filt_p,filt_2,ipv_f,QW,ldb,INFO)

   if (info > 0 .or. info < 0) then
    write(*,*) 'Problemi soluzione, info:', info
    stop
   endif

 end subroutine
