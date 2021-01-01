       subroutine derx_1(q3, q4)

!----------------------------------------------------------- !
!               F. Califano, 2006                            !
!  First x-derivative by FFT (periodic boundary conditions)  !
!----------------------------------------------------------- !

!**************************************************
!  MPI PARALLEL VERSION: VALENTINI-FAGANELLO 2009
!*************************************************


use parameter_mod
use box_mod
use deriv_mod

IMPLICIT NONE

integer :: i

REAL(dp),DIMENSION (nx) :: q3, q4

!do i = 1, nx
!q3(i) = exp(x(i))
!enddo      



!     Soluzione del sistema lineare tridiagonale fattorizzato LU
!     Differenze finite compatte 3x5

    do i = 3, nx2
     q4(i) = a_1_3x5 * (q3(i+1) - q3(i-1)) + b_1_3x5 * (q3(i+2) - q3(i-2))
    enddo

!     Differenze finite compatte 3x3 in i=2, nx-1

    q4(2)   = a_1_3x3 * (q3(3) - q3(1))
    q4(nx1) = a_1_3x3 * (q3(nx) - q3(nx2))

SELECT CASE (ibc) ! *** 1 = free, 2 = fixed ***

 CASE (1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     FREE boundary condition (left and right), I derivative, i=1, i=nx  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     q4(1)  = a_1_bc * q3(1) + b_1_bc * q3(2) &
            + c_1_bc * q3(3) + d_1_bc * q3(4)

     q4(nx) = - a_1_bc * q3(nx)  - b_1_bc * q3(nx1) &
              - c_1_bc * q3(nx2) - d_1_bc * q3(nx3)

    CALL  DGTTRS(TRANS,nx,NRHS,dm_1,dd_1,dp_1,dw_1,ipv_d,q4,ndb,INFO)

      if (info > 0 .or. info < 0) then
      write(*,*) 'Problemi soluzione, info:', info
      stop
      endif

 CASE (2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    Derivata nulla al bordo (left and right), I derivative, i=1, i=nx  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    q4(1)  = 0.0
    q4(nx) = 0.0

    CALL  DGTTRS(TRANS,nx,NRHS,wm_1,wd_1,wp_1,ww_1,ipv_w,q4,ndb,INFO)

      if (info > 0 .or. info < 0) then
      write(*,*) 'Problemi soluzione, info:', info
      stop
      endif

 CASE (3)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    Derivata costante al bordo (left and right), I derivative, i=1, i=nx  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    q4(1)  = 0.0
    q4(nx) = 0.0

    CALL  DGTTRS(TRANS,nx,NRHS,qm_1,qd_1,qp_1,qq_1,ipv_q,q4,ndb,INFO)

      if (info > 0 .or. info < 0) then
      write(*,*) 'Problemi soluzione, info:', info
      stop
      endif

 CASE (4)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    Boundary Condition Basate sulle caratteristiche                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! in avanzamento e' attiva subroutine BC, per ora poniamo 
! la derivata in 1,nx e' calcolata con i punti interni        
! e serve solo per poter calcolare la derivata in 2,nx-1


     q4(1)  = a_1_bc * q3(1) + b_1_bc * q3(2) &
            + c_1_bc * q3(3) + d_1_bc * q3(4)

     q4(nx) = - a_1_bc * q3(nx)  - b_1_bc * q3(nx1) &
              - c_1_bc * q3(nx2) - d_1_bc * q3(nx3)

    CALL  DGTTRS(TRANS,nx,NRHS,dm_1,dd_1,dp_1,dw_1,ipv_d,q4,ndb,INFO)

      if (info > 0 .or. info < 0) then
      write(*,*) 'Problemi soluzione, info:', info
      stop
      endif


  END SELECT

!write(*,*) 'derivata numerica'
!write(*,*) q4(1)
!write(*,*) 'valore analitico'
!write(*,*) exp(x(1))

!write(*,*) 'derivata numerica'
!write(*,*) q4(10)
!write(*,*) 'valore analitico'
!write(*,*) exp(x(10))

!write(*,*) 'derivata numerica'
!write(*,*) q4(nx)
!write(*,*) 'valore analitico'
!write(*,*) exp(x(nx))
!call syncronize()
!stop

end subroutine
