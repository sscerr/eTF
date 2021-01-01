! ----------------------------------------- !
!      Set up filtro periodico              !
!             F. Califano, Marzo 2001       !
! ----------------------------------------- !

!**************************************************
!  MPI PARALLEL VERSION: VALENTINI-FAGANELLO 2009
!   3D PARALLEL VERSION: FAGANELLO 2010
!**************************************************

           subroutine init_filtro_per

use box_mod
use parameter_mod
use filtro_per_mod

IMPLICIT NONE

  integer  :: kk
  real(dp) :: tr1, tr2, tr3, tr4, tr5, w0, w11, w22

!   Set up coefficienti filtro periodico in x !!!!!!!!!!!!!!!!refuso!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!   filter parameters

       tr1 = (2.0 + 3.0 * omega_x) / 4.0
       tr2 = (6.0 + 7.0 * omega_x) / 8.0
       tr3 = (6.0 + omega_x)       / 20.0
       tr4 = (2.0 - 3.0 * omega_x) / 40.0
       tr5 = (3.0 - 2.0 * omega_x) / 10.0


 do kk = 1, n2x
  w0 = 2.d0 * pgreco * dfloat(kk) / dfloat(nx)
  w11 = tr1 + tr2*cos(w0) + tr3*cos(2.d0*w0) + tr4*cos(3.d0*w0)
  w22 = 1.0 +2.0 * omega_x*cos(w0) + 2.0 * tr5*cos(2.d0*w0)
  work4x(kk) = w11 / w22
 enddo



!   Set up coefficienti filtro periodico in y

!   filter parameters

       tr1 = (2.0 + 3.0 * omega_y) / 4.0
       tr2 = (6.0 + 7.0 * omega_y) / 8.0
       tr3 = (6.0 + omega_y)       / 20.0
       tr4 = (2.0 - 3.0 * omega_y) / 40.0
       tr5 = (3.0 - 2.0 * omega_y) / 10.0


 do kk = 1, n2y
  w0 = 2.d0 * pgreco * dfloat(kk) / dfloat(ny)
  w11 = tr1 + tr2*cos(w0) + tr3*cos(2.d0*w0) + tr4*cos(3.d0*w0)
  w22 = 1.0 +2.0 * omega_y*cos(w0) + 2.0 * tr5*cos(2.d0*w0)
  work4y(kk) = w11 / w22
 enddo

!   Set up coefficienti filtro periodico in z

!   filter parameters

       tr1 = (2.0 + 3.0 * omega_z) / 4.0
       tr2 = (6.0 + 7.0 * omega_z) / 8.0
       tr3 = (6.0 + omega_z)       / 20.0
       tr4 = (2.0 - 3.0 * omega_z) / 40.0
       tr5 = (3.0 - 2.0 * omega_z) / 10.0


 do kk = 1, n2z
  w0 = 2.d0 * pgreco * dfloat(kk) / dfloat(nz)
  w11 = tr1 + tr2*cos(w0) + tr3*cos(2.d0*w0) + tr4*cos(3.d0*w0)
  w22 = 1.0 +2.0 * omega_z*cos(w0) + 2.0 * tr5*cos(2.d0*w0)
  work4z(kk) = w11 / w22
 enddo


 end subroutine
