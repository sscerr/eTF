subroutine parametri 

!*** INPUT/SETUP: simulation's parameters / output files ***!


!----------------------------------------------------!
!  ORIGINAL MODEL:    F. Califano, 2006              !
!                     M. Faganello, 2008             !
!                                                    !
!  Faganello et al., New J. Phys. 11, 063008 (2009)  !
!                                                    !
!  ================================================  !
!  EXTENDED MODEL:    S. S. Cerri, 2011              !
!                                                    !
!  Cerri et al., Phys. Plasmas 20, 112112 (2013)     !
!                                                    !
!--------------------------------------------------- !

!**************************************************!
!  MPI PARALLEL VERSION: VALENTINI-FAGANELLO 2009  !
!  3D PARALLEL VERSION: FAGANELLO 2010             !
!  3D ANISOTROPIC/FLR H-MHD/TF: CERRI 2011         !
!**************************************************!

use parameter_mod
use box_mod
use filtro_per_mod
use fields_EB_mod
use parallel_mod, only: npcol, nprow, mycol, myrow, &
                        mpime, nproc, root


IMPLICIT NONE

INTEGER :: ix

CHARACTER*15 filen_U
CHARACTER*16  filen_EB
CHARACTER*16  filen_Gi
CHARACTER*17  filen_DPJ
CHARACTER*1 F_1
CHARACTER*2 F_2
CHARACTER*3 F_3
CHARACTER*4 F_4

! *** Set up parameters

pgreco = dacos(-1.0d0)

ninvy = 1.0d0 / ny
ninvz = 1.0d0 / nz

if( mpime == root ) then

  open(unit=7,status='old',file='2fl.com')

  read (7,*)
  read (7,*) prefix
  read (7,*)
  read (7,*) istart, nwx, nwy, nwz, ibc, hall_on, divPe_on, flr_on
  read (7,*)
  read (7,*) xl, yl, zl, B00, angolo, Aeq, Leq, DeltaDen
  read (7,*)
  read (7,*) dt, rk_ord, tmax, tt_w, tx_w, nx_lambda, grdlmbd
  read (7,*)
  read (7,*) omega_x, omega_y, omega_z, eta, ampl, alpha
  close(7)

  ! finite-electron-inertia effects: not available in this release
  rapm = 9999.0d0
  ! polytropic sub-version: not available in this release
  gam = 1.0d0

end if

nx_lmbd = int(nx_lambda)

call BCAST_CHARACTER( prefix, 4, root )
call BCAST_INTEGER( istart, 1, root )
call BCAST_INTEGER( nwx, 1, root )
call BCAST_INTEGER( nwy, 1, root )
call BCAST_INTEGER( nwz, 1, root )
call BCAST_INTEGER( ibc, 1, root )
call BCAST_INTEGER( hall_on, 1, root )
call BCAST_INTEGER( divPe_on, 1, root )
call BCAST_INTEGER( flr_on, 1, root )
call BCAST_REAL( xl, 1, root )
call BCAST_REAL( yl, 1, root )
call BCAST_REAL( zl, 1, root )
call BCAST_REAL( B00, 1, root )
call BCAST_REAL( angolo, 1, root )
call BCAST_REAL( Aeq, 1, root )
call BCAST_REAL( Leq, 1, root )
call BCAST_REAL( DeltaDen, 1, root )
call BCAST_REAL( dt, 1, root )
call BCAST_INTEGER( rk_ord, 1, root )
call BCAST_REAL( tmax, 1, root )
call BCAST_REAL( tt_w, 1, root )
call BCAST_REAL( tx_w, 1, root )
call BCAST_REAL( rapm, 1, root )
call BCAST_REAL( gam, 1, root )
call BCAST_REAL( omega_x, 1, root )
call BCAST_REAL( omega_y, 1, root )
call BCAST_REAL( omega_z, 1, root )
call BCAST_REAL( eta, 1, root )
call BCAST_REAL( ampl, 1, root )
call BCAST_REAL( alpha, 1, root )
call BCAST_REAL( grdlmbd, 1, root )
call BCAST_INTEGER( nx_lmbd, 1, root )
call BCAST_REAL( nx_lambda, 1, root )

if (mpime==0) then
  write(*,*) ' ********************************* '
  write(*,*) ' Starting Run DUE - FLUIDI:', prefix
  write(*,*) ' ********************************* '

  write(*,*) 'xl, yl, zl, B00, angolo, Aeq, Leq'
  write(*,805) xl, yl, zl, B00, angolo, Aeq, Leq
  write(*,*) ' ------------------- '
  write(*,*) 'dt, rk_ord, tmax, tt_w, tx_w'
  write(*,815) dt, rk_ord, tmax, tt_w, tx_w
  write(*,*) ' ------------------- '
  write(*,*) 'mp/me, gam, eta, ampl '
  write(*,804) rapm, gam, eta, ampl
  write(*,*) ' ------------------- '
  write(*,*) 'n2x, n2y, n2z, nx1, ny1, nz1'
  write(*,706) n2x, n2y, n2z, nx1, ny1, nz1
  write(*,*) ' ------------------- '
endif

! finite-electron-inertia effects: 
! not available in this release
de2 = 0.0
de  = 0.0
de12inv = 1.0


if (mpime==0) then
   write(*,*) ' ************************* '
   write(*,*) '    GENERALIZED OHM LAW    '
   write(*,*) '  -----------------------  '
   write(*,*) '        [ setting ]        '
   write(*,*) '                           '
   if (hall_on .eq. 1) then
     write(*,*) '  Hall term:  ON           '
   else
     write(*,*) '  Hall term:  OFF          '
   endif
   if (divPe_on .eq. 1) then
     write(*,*) '  div[Pe] term:  ON        '
   else
     write(*,*) '  div[Pe] term:  OFF       '
   endif
   if (eta .GT. 0) then
     write(*,*) '  Resistivity:  ON         '
   else
     write(*,*) '  Resistivity:  OFF        '
   endif
   write(*,*) '                           '
   write(*,*) ' ************************* '
   write(*,*) '   FINITE LARMOR RADIUS    '
   write(*,*) '  -----------------------  '
   write(*,*) '        [ setting ]        '
   write(*,*) '                           '
   if (flr_on .eq. 1) then
     write(*,*) '  ion FLR:  ON             '
   else
     write(*,*) '  ion FLR:  OFF            '
   endif
   write(*,*) '                           '
   write(*,*) ' ************************* '
endif


! Set up smoothing function

  terzo = 1.0d0 / 3.0d0
  lambda = 1.0d0
  

IF (eta .GE. 0.0d0) then
  IF (nx_lambda .GT. 0.0d0) then

    do ix = 1, nx_lambda + 1
      lambda(ix) =  0.5d0 * ( 1.0d0 + &
           tanh( grdlmbd * ( ix - 1 - 0.5d0 * nx_lambda) / nx_lambda ))
    enddo
    do ix = nx - nx_lambda, nx
      lambda(ix) =  0.5d0 * ( 1.0d0 - &
           tanh( grdlmbd * ( ix - nx + 0.5d0 * nx_lambda) / nx_lambda ))
    enddo

  ENDIF
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! *** Set up time scheme 

        ab1 = dt * 5.d0 / 12.d0
        ab2 = - dt * 4.d0 / 3.d0
        ab3 = dt * 23.d0 / 12.d0

! *** Set up numerical box (in 2*pi*d_i units) along periodic directions
! *** x: open boundary
! *** y,z: periodic boundaries 

        hyl = 1.0d0 / yl
        hzl = 1.0d0 / zl

        yl = 2.0d0 * pgreco * yl
        zl = 2.0d0 * pgreco * zl

        dx  = xl / (nx-1)   ! direzione NON periodica
        dy  = yl / ny
       if (nz .GT. 1) then
        dz  = zl / nz
       else
        dz  = 1.0d0
       endif
        dxyz = dx * dy * dz  


        Tot_inv_ny   = 1. / (float(ny))
        Tot_inv_nz   = 1. / (float(nz))          
        Tot_inv_nxnynz = 1. / (float(nx) * float(ny) * float(nz))  


      if (mpime==0) then
        write(*,*) 'dx, dy, dz, xl, yl, zl, Nx, Ny, Nz'
        write(*,836) dx, dy, dz, xl, yl, zl, Nx, Ny, Nz
        write(*,*) '  '
       endif

!  Open files for output 

l_EB_me = l_EB + mpime
l_U_me = l_U + mpime
l_DPJ_me = l_DPJ + mpime
l_Gi_me = l_Gi + mpime

if ( mpime.lt.10 ) then
  write( F_1, '(1i1)' ) mpime
  filen_U   = prefix//'_000'//F_1//'_U.dat'
  filen_EB  = prefix//'_000'//F_1//'_EB.dat'
  filen_DPJ = prefix//'_000'//F_1//'_DPJ.dat'
  filen_Gi  = prefix//'_000'//F_1//'_Gi.dat'
endif
if ( (mpime.ge.10) .and. (mpime.le.99) ) then
  write( F_2, '(1i2)' ) mpime
  filen_U   = prefix//'_00'//F_2//'_U.dat'
  filen_EB  = prefix//'_00'//F_2//'_EB.dat'
  filen_DPJ = prefix//'_00'//F_2//'_DPJ.dat'
  filen_Gi  = prefix//'_00'//F_2//'_Gi.dat'
endif
if ( (mpime.gt.99) .and. (mpime.le.999) ) then
  write( F_3, '(1i3)' ) mpime
  filen_U   = prefix//'_0'//F_3//'_U.dat'
  filen_EB  = prefix//'_0'//F_3//'_EB.dat'
  filen_DPJ = prefix//'_0'//F_3//'_DPJ.dat'
  filen_Gi  = prefix//'_0'//F_3//'_Gi.dat'
endif
if ( (mpime.gt.999) .and. (mpime.le.9999) ) then
  write( F_4, '(1i4)' ) mpime
  filen_U   = prefix//'_'//F_4//'_U.dat'
  filen_EB  = prefix//'_'//F_4//'_EB.dat'
  filen_DPJ = prefix//'_'//F_4//'_DPJ.dat'
  filen_Gi  = prefix//'_'//F_4//'_Gi.dat'
endif


SELECT CASE (istart)

! first RUN (from initial condition, at t=0)

  CASE (0)

        tempo   = 0.0 
        tt_last = 0.0
        tx_last = 0.0


!!! OPEN (unit = l_U_me,   STATUS = 'unknown', form='unformatted', FILE = filen_U)
!!! OPEN (unit = l_EB_me,  STATUS = 'unknown', form='unformatted', FILE = filen_EB)
!!! OPEN (unit = l_DPJ_me, STATUS = 'unknown', form='unformatted', FILE = filen_DPJ)
!!! OPEN (unit = l_Gi_me,  STATUS = 'unknown', form='unformatted', FILE = filen_Gi)

 OPEN (unit = l_U_me,   STATUS = 'unknown', FILE = filen_U)
 OPEN (unit = l_EB_me,  STATUS = 'unknown', FILE = filen_EB)
 OPEN (unit = l_DPJ_me, STATUS = 'unknown', FILE = filen_DPJ)
 OPEN (unit = l_Gi_me,  STATUS = 'unknown', FILE = filen_Gi)

        if( mpime == root ) then
        OPEN (unit = l_div, STATUS = 'unknown', FILE = prefix//'_Div.dat')
        OPEN (unit = l_Eng, STATUS = 'unknown', FILE = prefix//'_Energie.dat')
        OPEN (unit = l_cpu, STATUS = 'unknown', FILE = prefix//'_cpu.dat')
        WRITE (l_cpu, 706) nx, ny, nz, nwx, nwy, nwz
        WRITE (l_cpu, 803) xl, yl, zl
        end if

! Restart of an existing RUN

  CASE(+1,+2)

        if( mpime == root ) then

 OPEN (unit = l_cpu, STATUS = 'unknown', POSITION = 'append' &
                             , FILE = prefix//'_cpu.dat')

 OPEN (unit = l_div, STATUS = 'unknown', POSITION = 'append' &
                             , FILE = prefix//'_Div.dat')

 OPEN (unit = l_Eng, STATUS = 'unknown', POSITION = 'append' &
                             , FILE = prefix//'_Energie.dat')
       endif

!!! OPEN (unit = l_U_me,   STATUS = 'unknown', form='unformatted', FILE = filen_U)
!!! OPEN (unit = l_EB_me,  STATUS = 'unknown', form='unformatted', FILE = filen_EB)
!!! OPEN (unit = l_DPJ_me, STATUS = 'unknown', form='unformatted', FILE = filen_DPJ)
!!! OPEN (unit = l_Gi_me,  STATUS = 'unknown', form='unformatted', FILE = filen_Gi)

 OPEN (unit = l_U_me,   STATUS = 'unknown', FILE = filen_U)
 OPEN (unit = l_EB_me,  STATUS = 'unknown', FILE = filen_EB)
 OPEN (unit = l_DPJ_me, STATUS = 'unknown', FILE = filen_DPJ)
 OPEN (unit = l_Gi_me,  STATUS = 'unknown', FILE = filen_Gi)


  END SELECT


702   format(1x, 2(2x, 1i5))
703   format(1x, 3(2x, 1i5))
704   format(1x, 4(2x, 1i5))
706   format(1x, 6(2x, 1i5))
802   format(1x, 2(1x, 1e11.4))
803   format(1x, 3(1x, 1e11.4))
804   format(1x, 4(1x, 1e11.4))
815   format(1x, 1e11.4, 1(1x, 1i5), 3(1x, 1e11.4))
836   format(1x, 6(1x, 1e11.4), 3(1x, 1i4))
805   format(1x, 5(1x, 1e11.4))
806   format(1x, 6(1x, 1e11.4))

end subroutine
