       subroutine dery_1(q0, q1)

!----------------------------------------------------------- !
!               F. Califano, 2006                            !
!  First y-derivative by FFT (periodic boundary conditions)  !
!----------------------------------------------------------- !

!**************************************************
!  MPI PARALLEL VERSION: VALENTINI-FAGANELLO 2009
!**************************************************

use parameter_mod
use box_mod
use deriv_mod

IMPLICIT NONE

integer :: ik, il

REAL(dp),DIMENSION (ny) :: q0, q1, qq

    qq = q0

    call drfftf(ny, qq, savey)

       il = 1
          do ik = 2, ny1, 2
            q1(ik)   = - qq(ik+1) * ky_1(il)
            q1(ik+1) =   qq(ik)   * ky_1(il)
           il = il + 1
          enddo

            q1(1)  = 0.0
            q1(ny) = 0.0

    call drfftb(ny, q1, savey)

end subroutine
