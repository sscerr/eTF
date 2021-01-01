          subroutine filtro_per_y(yk)

! ------------------------------------------------------ !
!             filtro periodico                           !
!             F. Califano, Marzo 2001                    !
! ------------------------------------------------------ !

!**************************************************
!  MPI PARALLEL VERSION: VALENTINI-FAGANELLO 2009
!**************************************************

use parameter_mod
use box_mod
use filtro_per_mod
use deriv_mod, only: savey

IMPLICIT NONE

  integer :: m, m1, m2
  REAL(dp), dimension(ny)  :: yk

!  filtro 

! **    y Filtering

        call drfftf(ny, yk, savey)

        m = 1

          do m1 = 2, ny - 1, 2

          m2 = m1 + 1

          yk(m1) = yk(m1) * work4y(m)
          yk(m2) = yk(m2) * work4y(m)

        m = m + 1

!write(*,898) m, work4y(m)

        enddo

!898     format(1x, 1i3, 1x, 1e12.5)
!stop

          yk = yk * ninvy

! **    Filtered function

        call drfftb(ny, yk, savey)
       

       return
       end 
