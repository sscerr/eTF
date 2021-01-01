program TWOFLUIDS

!*****************************************************!
!                                                     !
! SOLVER: extended-Two-Fluid (eTF) model              !
!                                                     !
! RELEASE: version 1.0                                !
!                                                     !
! =================================================== !
!                                                     !
!  ______________[ Model's equations ]______________  !
!                                                     !
!  ORIGINAL TF MODEL:    F. Califano, 2006            !
!                        M. Faganello, 2008           !
!                                                     !
!   Faganello et al., New J. Phys. 11, 063008 (2009)  !
!                                                     !
!  - - - - - - - - - - - - - - - - - - - - - - - - -  !
!                                                     !
!  EXTENDED TF MODEL:    S. S. Cerri, 2011            !
!                                                     !
!   Cerri et al., Phys. Plasmas 20, 112112 (2013)     !
!                                                     !
! =================================================== !
!                                                     !
! ______________[ development history ]______________ !
!                                                     !
!  MPI PARALLEL VERSION: VALENTINI-FAGANELLO 2009     !
!  3D PARALLEL VERSION: FAGANELLO 2010                !
!  3D ANISOTROPIC/FLR H-MHD/TF: CERRI 2011            !
!                                                     !
!*****************************************************!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!                                             !!!
!!!  *** NOTICE: features not yet released ***  !!!
!!!                                             !!!
!!!   - 3D version                              !!!
!!!   - finite-m_e contributions                !!!
!!!   - Landau-fluid closures                   !!!
!!!                                             !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! -------------------------------------------------  !
!   nx = number of grid points in the x-direction    !
!   ny = number of grid points in the y-direction    !
!                                                    !
!   xl = box dimension in the x-direction            !
!   yl =  box dimension in the y-direction           !
! -------------------------------------------------- !

use parameter_mod
use box_mod
use fields_UJ_mod
use fields_DP_mod
use fields_EB_mod
use deriv_mod
use filtro_per_mod
use poisson_open
use parallel_mod, only: root, mpime, nprow, npcol, nproc, group, &
                        parallel_build
use dom_distr_mod

implicit none

integer :: it1, it2, it3, it4, it5, iclk
integer :: ix

real(dp) :: timeinit, timestart, timesave, timetot, timeloop

real(dp), allocatable :: RHS_x(:,:,:), RHS_y(:,:,:), AAsmooth(:,:,:)

character(len=6), external :: int_to_char
character*12 filen_DT
character*1 F_1
character*2 F_2
character*3 F_3
character*4 F_4


!*******************************************!
!**  -----------------------------------  **!
!**  SETTING UP ENVIRONMENT & SIMULATION  **!
!**  -----------------------------------  **!
!*******************************************!  

! initialize MPI environment
call parallel_startup( nproc, mpime, root, group )
call parallel_commlib( parallel_build )
call parallel_init()

! allocate auxiliary arrays for time advance,
! boundary conditions, and of physical quantities
allocate( RHS_x( nxl, nyl, nzl ) )  
allocate( RHS_y( nxl, nyl, nzl ) )
if (nz.GT.2) then !note: 3D version not yet released.
  allocate( kz_1( nz/2 - 1 ) )
  allocate( kz_4( nz/2 - 1 ) )
  allocate( work4z( nz/2 ) )
else
  allocate( kz_1(1) )
  allocate( kz_4(1) )
  allocate( work4z(1) )
endif
call allocate_fields_UJ( nxl, nyl, nzl )
call allocate_fields_DP( nxl, nyl, nzl )
call allocate_fields_EB( nxl, nyl, nzl )
call allocate_poisson( n2z )

call system_clock(it1)

! Initialization 
call parametri          ! read parameters 
call griglia            ! grid
call inizializza        ! fourier coeff, etc. 
call init_filtro_per    ! periodic filter
call init_filtro_open_x ! open filter

call system_clock(it2)

 
if (mpime==0) then
  write(*,*) ' --------------------- '
  write(*,*) '   Initial Condition   '
  write(*,*) ' --------------------- ' 
endif

!--set up initial condition
call condinit   

!--output box-integrated quantities
!call outt

!--compute: n*U terms
nU_x = Den * Ux  
nU_y = Den * Uy  
nU_z = Den * Uz


!--------------------------------------------------------!
! smoothing -> ideal (isotropic) MHD at the x-boundaries !
!                                                        ! 
! [compute "effective" isotropic pressures = 1/3 Tr(P)]  !
!--------------------------------------------------------!

allocate( AAsmooth(nxl,ny,nz) )

!ions
AAsmooth = terzo * ( 2.0d0 * pi_perp + pi_para ) 
do ix = 1, nxl
  pi_para(ix,:,:) = lambda(ixlg + ix - 1) * pi_para(ix,:,:) + &
                   ( 1.0d0 - lambda(ixlg + ix - 1) ) * AAsmooth(ix,:,:)
  pi_perp(ix,:,:) = lambda(ixlg + ix - 1) * pi_perp(ix,:,:) + &
                   ( 1.0d0 - lambda(ixlg + ix - 1) ) * AAsmooth(ix,:,:)
enddo
Ti_para = Dinv * pi_para
Ti_perp = Dinv * pi_perp

!electrons
AAsmooth = terzo * ( 2.0d0 * pe_perp + pe_para )
do ix = 1, nxl
  pe_para(ix,:,:) = lambda(ixlg + ix - 1) * pe_para(ix,:,:) + &
                   ( 1.0d0 - lambda(ixlg + ix - 1) ) * AAsmooth(ix,:,:)
  pe_perp(ix,:,:) = lambda(ixlg + ix - 1) * pe_perp(ix,:,:) + &
                   ( 1.0d0 - lambda(ixlg + ix - 1) ) * AAsmooth(ix,:,:)
enddo
Te_para = Dinv * pe_para
Te_perp = Dinv * pe_perp

deallocate(AAsmooth)


if (mpime==0) then
  write(*,*) ' ---------------------------------------- '
  write(*,*) ' ---------------------------------------- '
  write(*,*) '      extended Two-Fluids (eTF) code      '
  write(*,*) ' ........................................ '
  write(*,*) '          [ OPEN x-boundaries ]           '
  write(*,*) '  [ SMOOTH OUT anisotropy @ boundaries ]  '
  write(*,*) ' ---------------------------------------- '
  write(*,*) ' ---------------------------------------- '
endif



!*******************************************!
!**  -----------------------------------  **!
!**  SETTING UP PRELIMINARY TIME ADVANCE  **! 
!**  -----------------------------------  **!
!*******************************************!

SELECT CASE (istart)

  CASE (0) !FIRST RUN (t=0, from initial condition)

    ioutt = 0
    ioutx = 0

    !--initial output at t=0 
    call outx ! grid-based quantities 
    call outt ! box-averaged quantities

    ioutt = 1
    ioutx = 1

    !--First time advance by dt 

    ! set up parametrers for 1st time advance 
    ab1 = 0.0
    ab2 = 0.0
    ab3 = dt 

    ! time-advance routine
    call avanzamento

    ! update time
    tempo = dt

    !--Second time advance (Adams - Basfort 2)

    ! set up parametrers for 2nd time advance 
    ab1 = 0.0
    ab2 = - dt * 1.0 / 2.0
    ab3 = dt * 3.0 / 2.0

    ! time-advance routine
    call avanzamento

    ! update time
    tempo = tempo + dt


  !CASE (+1) !RESTART with same dt as previous run
  !
  ! ** not needed anymore, just go to main loop **
  

  CASE (+2) !RESTART with different dt w.r.t. previous run

    !--First time advance with new dt 

    ! set up parametrers for 1st time advance with new dt 
    ab1 = 0.0
    ab2 = 0.0
    ab3 = dt

    ! time-advance routine
    call avanzamento

    ! update time
    tempo = tempo + dt

    !--Second time advance with new dt (Adams - Basfort 2)

    ! set up parametrers for 2nd time advance 
    ab1 = 0.0
    ab2 = - dt * 1.0 / 2.0
    ab3 = dt * 3.0 / 2.0

    ! time-advance routine
    call avanzamento

    ! update time
    tempo = tempo + dt

END SELECT

call system_clock(it3)


!******************************************************!
!**  ----------------------------------------------  **!
!**  MAIN LOOP for time advance (Adams - Basfort 3)  **! 
!**  ----------------------------------------------  **!
!******************************************************!

! set up parametrers for main time loop 
ab1 = dt * 5.0 / 12.0
ab2 = - dt * 4.0 / 3.0
ab3 = dt * 23.0 / 12.0

! enter main time loop
do while (tempo <= tmax)

  ! time-advance routine
  call avanzamento

  ! update time
  tempo = tempo + dt

  !-----------------!
  !  write outputs  !
  !-----------------!
  if(tempo >= tt_last + tt_w) then
    call outt  
    tt_last = tempo
    ioutt = ioutt + 1
  endif
  if(tempo >= tx_last + tx_w) then
    call outx  
    tx_last = tempo  
    ioutx = ioutx + 1
  endif

enddo

if (mpime==0) then
  write(l_cpu, 702) ioutt, ioutx, tempo
endif

if(mpime==0) then
  close(unit = l_cpu)
  close(unit = l_div)
  close(unit = l_Eng)
endif

close(unit = l_U_me)
close(unit = l_EB_me)
close(unit = l_DPJ_me)
close(unit = l_Gi_me)

deallocate(RHS_x)
deallocate(RHS_y)

call system_clock(it4)

!-------------------------!
! saving data for RESTART !
!-------------------------!

l_rst_me = l_rst + mpime

if ( mpime.lt.10 ) then
  write( F_1, '(1i1)' ) mpime
  filen_DT  = 'DATA_TF_000'//F_1
endif
if ( (mpime.ge.10) .and. (mpime.le.99) ) then
  write( F_2, '(1i2)' ) mpime
  filen_DT  = 'DATA_TF_00'//F_2
endif
if ( (mpime.gt.99) .and. (mpime.le.999) ) then
  write( F_3, '(1i3)' ) mpime
  filen_DT  = 'DATA_TF_0'//F_3
endif
if ( (mpime.gt.999) .and. (mpime.le.9999) ) then
  write( F_4, '(1i4)' ) mpime
  filen_DT  = 'DATA_TF_'//F_4
endif


open(l_rst_me,status='unknown', form='unformatted',file=filen_DT)

write(l_rst_me) ioutt, ioutx
write(l_rst_me) tempo, tt_last, tx_last
write(l_rst_me) Ti_para, Ti_perp, Te_para, Te_perp
write(l_rst_me) Uex, Uey, Uez, Uix, Uiy, Uiz, Ux, Uy, Uz, Jx, Jy, Jz
write(l_rst_me) Den, Dinv, pe_para, pe_perp, pi_para, pi_perp, nU_x, nU_y, nU_z
write(l_rst_me) Ex, Ey, Ez, Bx, By, Bz, Tracciante!, Tracciante_e

close(l_rst_me)

call system_clock(it5,count_rate=iclk)
if(iclk .gt. 0) then
    timestart = dble(it2-it1)/dble(iclk)
    timeinit  = dble(it3-it2)/dble(iclk)
    timeloop  = dble(it4-it3)/dble(iclk)
    timesave  = dble(it5-it4)/dble(iclk)
    timetot   = dble(it5-it1)/dble(iclk)
end if

if(iclk .gt. 0) then
    timestart = dble(it2-it1)/dble(iclk)
    timeinit  = dble(it3-it2)/dble(iclk)
    timeloop  = dble(it4-it3)/dble(iclk)
    timesave  = dble(it5-it4)/dble(iclk)
    timetot   = dble(it5-it1)/dble(iclk)
end if

if (mpime==0) then
  write(*,*) ' '
  write(*,*) ' ------------------------------------ '
  write(*,*) ' Two Fluids program timing: '
  write(*,*) ' '
  write(*,*) '  Startup       Init       Loop      Write      Total'
  write(*,855) timestart, timeinit, timeloop, timesave, timetot
  write(*,*) ' '
  write(*,*) ' END Two Fluids: '
  write(*,*) ' ------------------------------------ '
  write(*,*) ' '
endif

702   format(1x, 2(2x, 1i5), 1x, 1e11.4)
802   format(1x, 2(1x, 1e11.4))
805   format(1x, 5(1x, 1e11.4))
855   format(1x, 5(1x, 1e10.3))


call hangup()
stop 'end.'

end program TWOFLUIDS

