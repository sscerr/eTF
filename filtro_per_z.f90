          subroutine filtro_per_z(zk)

! ------------------------------------------------------ !
!             filtro periodico                           !
!             F. Califano, Marzo 2001                    !
! ------------------------------------------------------ !

!**************************************************
!  MPI PARALLEL VERSION: VALENTINI-FAGANELLO 2009
!   3D PARALLEL VERSION: FAGANELLO 2010
!**************************************************

use parameter_mod
use box_mod
use filtro_per_mod
use deriv_mod, only: savez

IMPLICIT NONE

  integer :: m, m1, m2
  REAL(dp), dimension(nz)  :: zk

!  filtro 

! **    z Filtering

        call drfftf(nz, zk, savez)

        m = 1

          do m1 = 2, nz - 1, 2

          m2 = m1 + 1

          zk(m1) = zk(m1) * work4z(m)
          zk(m2) = zk(m2) * work4z(m)

        m = m + 1

!write(*,898) m, work4z(m)

        enddo

!898     format(1x, 1i3, 1x, 1e12.5)
!stop

          zk = zk * ninvz

! **    Filtered function

        call drfftb(nz, zk, savez)
       

       return
       end 
