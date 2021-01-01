       subroutine derx_BC(q3, q4)

!----------------------------------------------------------- !
!               F. Califano, 2006                            !
!  First x-derivative by FFT (periodic boundary conditions)  !
!----------------------------------------------------------- !

!**************************************************
!  MPI PARALLEL VERSION: VALENTINI-FAGANELLO 2009
!  3D PARALLEL VERSION: FAGANELLO 2010
!**************************************************

use parameter_mod
use box_mod
use deriv_mod
use dom_distr_mod
use parallel_mod

IMPLICIT NONE

! *************N.B. ny = nyl*****************************

REAL(dp) :: frac
REAL(dp),DIMENSION (nxl,nyl,nzl) :: q3 !!!!!!!!!!!!!! con o senza l vari ???????????????????????????????????
REAL(dp),DIMENSION (2,nyl,nzl) :: q4 !!!!!!!!!!!!!!!!!!! con o senza l vari ????????????????????????????????

!integer :: ix
!do ix = 1, 5
!q3(ix,:) = exp(x(ix))
!enddo      
!do ix = nx4, nx
! q3(ix,:) = exp(x(ix))
!enddo      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Derivata al bordo sx calcolata con i soli punti interni            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


frac = 1.0d0 / ( 12.0d0 * dx )

     q4(1,:,:) = (- 25.0d0 * q3(1,:,:) + 48.0d0 * q3(2,:,:) - 36.0d0  * q3(3,:,:) &
                + 16.0d0 * q3(4,:,:) - 3.0d0 * q3(5,:,:) ) * frac 

     q4(2,:,:) = ( 25.0d0 * q3(nxl,:,:) - 48.0d0 * q3(nxlm1,:,:) + 36.0d0 &
               * q3(nxlm2,:,:) - 16.0d0 * q3(nxlm3,:,:) + &
                       3.0d0 * q3(nxlm4,:,:) ) * frac 


!write(*,*), 'derivata numerica'
!write(*,*), q4(1,1)
!write(*,*), 'valore analitico'
!write(*,*), exp(x(1))

!write(*,*), 'derivata numerica'
!write(*,*), q4(2,1)
!write(*,*), 'valore analitico'
!write(*,*), exp(x(nx))
!stop

end subroutine
