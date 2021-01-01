subroutine condinit 

 ! ------------------------------------ !
 !           Initial Condition          !
 !    [ Kelvi-Helmholtz instability ]   !
 ! ------------------------------------ !

 !**************************************************!
 !  MPI PARALLEL VERSION: VALENTINI-FAGANELLO 2009  !
 !   3D PARALLEL VERSION: FAGANELLO 2010            !
 !**************************************************!

 !**************************************************!
 !  FLR corrected MHD equilibria:   CERRI 2011      !
 !                                                  !
 !  for details see:                                !
 !                                                  !
 !  Cerri et al., Phys. Plasmas 20, 112112 (2013)   !
 !                                                  !
 !**************************************************!  

use parameter_mod
use dom_distr_mod
use box_mod
use fields_UJ_mod
use fields_DP_mod
use fields_EB_mod
use parallel_mod

IMPLICIT NONE

integer :: ix, iy, iz, ir, iky, ikz

REAL(dp) :: ky, kz, zkk, pph, xx, ab, zz
REAL(dp) :: n0, B0
REAL(dp) :: T0i_para, T0i_perp, T0e_para, T0e_perp
REAL(dp) :: gam_i_perp, gam_i_para, gam_e_perp, gam_e_para
REAL(dp) :: gaminv_i_perp, gam_tilde  
REAL(dp) :: beta_i_perp, beta_e_perp
REAL(dp) :: beta_i_perp_tilde, beta_e_perp_tilde
REAL(dp) :: C_0, C_1, C_inf 
REAL(dp) :: Den_inf, Dinv_inf!, DeltaDen
REAL(dp) :: pi_para_inf, pi_perp_inf, pe_para_inf, pe_perp_inf, anis_i_0, anis_e_0
REAL(dp) :: Ti_para_inf, Ti_perp_inf, Te_para_inf, Te_perp_inf

REAL(dp),ALLOCATABLE :: equil(:)  
REAL(dp),ALLOCATABLE :: f_0(:), fnct(:), equilDen(:), equilB(:)
REAL(dp),ALLOCATABLE :: gx1(:,:,:)
REAL(dp),ALLOCATABLE :: gy1(:), phs(:)
REAL(dp),ALLOCATABLE :: zx(:)
REAL(dp),ALLOCATABLE :: g0(:,:,:), gt(:,:,:)

 CHARACTER*13  filen_DT
 CHARACTER*1 F_1
 CHARACTER*2 F_2
 CHARACTER*3 F_3
 CHARACTER*4 F_4


allocate(equil(nxl))
allocate(f_0(nxl))
allocate(fnct(nxl))
allocate(equilB(nxl))
allocate(equilDen(nxl))
allocate(g0(nxl,ny,nz))
allocate(gx1(nxl,ny,nz))
allocate(gy1(ny))
allocate(phs(ny))
allocate(gt(nx,nyt,nz))
allocate(zx(nx))


SELECT CASE (istart)


 CASE (0)   !*** istart=0 ==> first RUN ***!


  !**********************************************!
  !  MHD one-fluid velocity: shear flow along x  !
  !                                              !
  !  Uy(x) = Aeq * tanh[ (x - Lx/2) / Leq ]      !    
  !                                              !
  !  (Leq in d_i units, Aeq in v_A units)        !
  !**********************************************!

   Ux = 0.
   Uz = 0.

   equil = Aeq * dtanh((x(ixlg:ixlg+nxl-1) - xl*0.5d0 ) / Leq)

   if (mpime==0) then
   write(6,*) ' ----------------------------- '
   write(6,*) ' U_0, L_U  = ', Aeq, Leq 
   write(6,*) ' ----------------------------- '
   endif

   do iz=1, nz
    do iy=1, ny  !  nyl == ny
     Uy(:,iy,iz) = equil(:)
    enddo
   enddo


  !------------------------------------------------------------!
  ! 2D incompressible perturbations on U: dU = e_z x grad(phi) !
  !------------------------------------------------------------! 

   ir = 1

   call random_seed(ir)
   
   call random_number(phs)


  ! Sum up modes on phi (2D --> no kz) 

   g0 = 0.

   do iky = 1, MIN(Ny/8,16)
     ky  = float(iky)
     zkk = 1.0 / ky
     pph = (phs(ir) - 0.5) * 2.0 * pgreco
     ir = ir + 1
     do iz = 1, nz 
       do iy = 1, ny
         do ix = 1, nxl
           xx = x(ixlg+ix-1) - xl/2.0
           ab = abs(xx)/Leq
           IF (ab .GE. 9.0)  then 
             g0(ix,iy,iz)= 0.0
           ELSE   
             g0(ix,iy,iz) = g0(ix,iy,iz) + cos(ky * y(iy) * hyl + pph) * zkk *&       
             exp(-1.0 * xx ** 2.0 / Leq ** 2.0)
           ENDIF
         enddo
       enddo
     enddo
   enddo

  ! compute: dU = ez x grad(phi) 
  ! (multiplied by ampl, so that max(|dU|) ~ ampl)

   CALL traspdist( g0, gt, 1 )
 
   do iz = 1, nz 
     do iy = 1, nyt
       call  derx_1( gt(:,iy,iz), zx )
       gt(:,iy,iz) = zx
     enddo
   enddo

   CALL traspdist( gx1, gt, -1 )

   Uy = Uy + ampl * gx1

   do iz = 1, nz 
     do ix = 1, nxl
       call  dery_1( g0(ix,:,iz), gy1 )
       Ux(ix,:,iz) = - ampl * gy1 
     enddo
   enddo

   Uz = 0.0 

deallocate(g0)
deallocate(gx1)
deallocate(gy1)
deallocate(equil)
deallocate(phs)
deallocate(gt)
deallocate(zx)


  !**************************************************!
  ! Compute FLR-corrected MHD equilibrium            !
  ! ------------------------------------------------ !
  !  for details see:                                !
  !                                                  !
  !  Cerri et al., Phys. Plasmas 20, 112112 (2013)   !
  !                                                  !
  !**************************************************!  

   !-- politropic index of the equilibrium
   !   (default: isothermal initial condition)

   gam_i_perp = 1.0d0
   gam_i_para = 1.0d0
  
   gam_e_perp = 1.0d0
   gam_e_para = 1.0d0

   gaminv_i_perp = 1.0d0 / gam_i_perp

   !-- equilibrium quantities "at infinity" 

   !density
   Den_inf = 1.0d0
   Dinv_inf = 1.0d0 / Den_inf
   !ion pressures and temperatures
   anis_i_0 = 0.0d0     ! ion anisotropy at infinity (A = Pperp/Ppara - 1)
   pi_para_inf = 0.5d0  ! (note: beta_inf = 0.5*p_inf/B0 => 0.5d0 -> beta = 1)
   pi_perp_inf = (1.0d0 + anis_i_0)*pi_para_inf
   Ti_para_inf = Dinv_inf * pi_para_inf
   Ti_perp_inf = Dinv_inf * pi_perp_inf
   !electron pressures and temperatures
   anis_e_0 = 0.0d0     ! electron anisotropy at infinity (A = Pperp/Ppara - 1) 
   pe_para_inf = 0.5d0 ! (0.5d0 --> beta = 1)
   pe_perp_inf = (1.0d0 + anis_e_0)*pe_para_inf
   Te_para_inf = Dinv_inf * pe_para_inf
   Te_perp_inf = Dinv_inf * pe_perp_inf
   !perpendicular plasma betas 
   beta_i_perp = (2.0d0 * pi_perp_inf) / B00    ! ions
   beta_e_perp = (2.0d0 * pe_perp_inf) / B00    ! elettrons

   ! useful quantities: "reduced" betas and gamma 
   ! [see: Cerri et al., PoP 20, 112112 (2013)]  
   gam_tilde = (gam_e_perp - gam_i_perp)*gaminv_i_perp
   beta_i_perp_tilde = beta_i_perp / (1.0d0 + beta_i_perp + beta_e_perp) 
   beta_e_perp_tilde = beta_e_perp / (1.0d0 + beta_i_perp + beta_e_perp)

   if (mpime==0) then
   write(6,*) ' ------------------------------------------------------------------ '
   write(6,*) ' gam_i_perp,  gam_i_para,  gam_e_perp,  gam_e_para,  gam_tilde  '
   write(6,*) ' ',gam_i_perp, gam_i_para, gam_e_perp, gam_e_para, gam_tilde
   write(6,*) ' Den_inf,  pi_para_inf,  pi_perp_inf,  pe_para_inf,  pe_perp_inf  '
   write(6,*) ' ',Den_inf, pi_para_inf, pi_perp_inf, pe_para_inf, pe_perp_inf
   write(6,*) ' Ti_para_inf,  Ti_perp_inf,  Te_para_inf,  Te_perp_inf  '
   write(6,*) ' ',Ti_para_inf, Ti_perp_inf, Te_para_inf, Te_perp_inf
   write(6,*) ' beta_i_perp,  beta_e_perp,  beta_i_perp_tilde,  beta_e_perp_tilde  '
   write(6,*) ' ',beta_i_perp, beta_e_perp, beta_i_perp_tilde, beta_e_perp_tilde 
   write(6,*) ' ------------------------------------------------------------------ '
   endif


   if (mpime==0) then
   write(6,*) ' ----------------------------- '
   write(6,*) ' Delta_N  = ', DeltaDen
   write(6,*) ' ----------------------------- '
   endif
   

  ! compute FLR-correction function for equivalent MHD profiles
  ! [see: Cerri et al., PoP 20, 112112 (2013)]  

  if (flr_on .eq. 1) then
    C_0 = (0.5d0/dsqrt(B00))*(Aeq / Leq)*beta_i_perp_tilde
    equilDen(:) = 1.0d0 - 0.5d0*DeltaDen*Dinv_inf*(1.0d0 - dtanh((x(ixlg:ixlg+nxl-1) - xl*0.5d0 ) / Leq) )
    f_0(:) = 1.0d0 / (1.0d0 - C_0*equilDen(:)*(1.0d0 - (dtanh((x(ixlg:ixlg+nxl-1) - xl*0.5d0 ) / Leq) )**2.0)  )
  else 
    f_0(:) = 1.0d0
  endif

  fnct(:) = equilDen(:) * f_0(:) 


  !*** Compute actual equilibrium profiles ***! 

   !density
   do iz=1, nz
     do iy=1, ny  !  nyl == ny
      Den(:,iy,iz) = Den_inf * (fnct(:)**gaminv_i_perp)
     enddo
   enddo

   Dinv = 1.0d0 / Den

   !gyrotropic ion pressures
   do iz = 1, nz
     do iy = 1, ny  ! nyl == ny
       pi_perp(:,iy,iz) = pi_perp_inf * fnct(:) 
     enddo
   enddo 
   do iz=1, nz
     do iy=1, ny  ! nyl == ny
       pi_para(:,iy,iz) = pi_para_inf * (fnct(:)**(gam_i_para*gaminv_i_perp))
     enddo
   enddo

   !gyrotropic electron pressures
   do iz = 1, nz
     do iy = 1, ny  ! nyl == ny
       pe_perp(:,iy,iz) = pe_perp_inf * (fnct(:)**(gam_e_perp*gaminv_i_perp))
     enddo
   enddo
   do iz=1, nz
     do iy=1, ny  ! nyl == ny
       pe_para(:,iy,iz) = pe_para_inf * (fnct(:)**(gam_e_para*gaminv_i_perp))
     enddo
   enddo

  !gyrotropic temperatures: simply T = p/n 
  Ti_perp = Dinv * pi_perp
  Ti_para = Dinv * pi_para
  Te_perp = Dinv * pe_perp
  Te_para = Dinv * pe_para


 !*** equilibrium magnetic field

  !make sure that B is in (y,z) plane
  !(i.e., not perpendicular to the flow)
  Bx = 0.d0 
 
  !FLR-corrected equilibrium profile
  equilB(:) = (1.0d0+beta_i_perp+beta_e_perp)*f_0(:) - (beta_i_perp+beta_e_perp)*fnct(:)
  do iz=1, nz
    do iy=1, ny  !  nyl == ny
      do ix=1, nxl
        Bz(ix,iy,iz) = dsqrt(B00*equilB(ix)) * dcos(angolo)
      enddo
    enddo
  enddo
  do iz=1, nz
    do iy=1, ny  !  nyl == ny
      do ix=1, nxl
        By(ix,iy,iz) = dsqrt(B00*equilB(ix)) * dsin(angolo)
      enddo
    enddo
  enddo

deallocate(f_0)
deallocate(fnct)
deallocate(equilB)
deallocate(equilDen)



 !*****************************************************************!
 ! Compute: passive tracer (defined on the initial velocity shear) !
 !   [ note: this is problem dependent, adjust at will ]           !
 !*****************************************************************!

  do iz = 1, nz
    do ix = 1, nxl
      Tracciante(ix,:,iz) = 0.6 + 0.4 * dtanh(( x(ix + ixlg - 1)-xl*0.5d0) / Leq)
    enddo
  enddo

 !***********************************!
 ! Compute: initial current density  !
 !***********************************!

  call corrente


 !********************************************************!
 ! Compute: initial fluid velocities (ions and electrons) !
 !********************************************************!

  call u_ei

 !*********************************!
 ! Compute: initial electric field !
 !*********************************!

  call ohm


 !**************************!
 ! Compute: inital momentum !
 !**************************!

  nU_x = Den * Ux
  nU_y = Den * Uy
  nU_z = Den * Uz 



 CASE(+1)   !istart = 1 -> Restart from existing RUN using the same timestep dt 

  !***********!
  !  RESTART  !
  ! (same dt) !
  !***********!

  l_rst_me=l_rst+mpime

 if ( mpime.lt.10 ) then
   write( F_1, '(1i1)' ) mpime
   filen_DT  = 'DATA_TF_000'//F_1
 endif

 if ( mpime.ge.10.and.mpime.le.99 ) then
   write( F_2, '(1i2)' ) mpime
   filen_DT  = 'DATA_TF_00'//F_2
 endif

 if ( mpime.gt.99 ) then
   write( F_3, '(1i3)' ) mpime
   filen_DT  = 'DATA_TF_0'//F_3
 endif

 if ( (mpime.gt.999) .and. (mpime.le.9999) ) then
   write( F_4, '(1i4)' ) mpime
   filen_DT  = 'DATA_TF_'//F_4
 endif


 OPEN(l_rst_me,status='unknown', form='unformatted',file=filen_DT)

 READ(l_rst_me) ioutt, ioutx
 READ(l_rst_me) tempo, tt_last, tx_last
 READ(l_rst_me) Ti_para, Ti_perp, Te_para, Te_perp
 READ(l_rst_me) Uex, Uey, Uez, Uix, Uiy, Uiz, Ux, Uy, Uz, Jx, Jy, Jz
 READ(l_rst_me) Den, Dinv, pe_para, pe_perp, pi_para, pi_perp, nU_x, nU_y, nU_z
 READ(l_rst_me) Ex, Ey, Ez, Bx, By, Bz, Tracciante!, Tracciante_e

 CLOSE(l_rst_me)



 CASE(+2)   !istart = 2 -> Restart from existing RUN, but using a different timestep dt 

  !****************!
  !     RESTART    !
  ! (different dt) !
  !****************!

  l_rst_me=l_rst+mpime

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


 OPEN(l_rst_me,status='unknown', form='unformatted',file=filen_DT)

 READ(l_rst_me) ioutt, ioutx
 READ(l_rst_me) tempo, tt_last, tx_last
 READ(l_rst_me) Ti_para, Ti_perp, Te_para, Te_perp
 READ(l_rst_me) Uex, Uey, Uez, Uix, Uiy, Uiz, Ux, Uy, Uz, Jx, Jy, Jz
 READ(l_rst_me) Den, Dinv, pe_para, pe_perp, pi_para, pi_perp, nU_x, nU_y, nU_z
 READ(l_rst_me) Ex, Ey, Ez, Bx, By, Bz, Tracciante

 CLOSE(l_rst_me)


END SELECT


 !*************************************!
 ! Compute: initial gyroviscous tensor !
 !*************************************!

  if (flr_on .eq. 1) then
    write(6,*) ''
    write(6,*) ' IONS FLR: ON! '
    call flr_i
  else !if (flr_on .eq. 0) then
    write(6,*) ''
    write(6,*) ' NO IONS FLR: G_i set to zero! '
    Gi_xx = 0.0d0
    Gi_xy = 0.0d0
    Gi_xz = 0.0d0
    Gi_yz = 0.0d0
  !else
  !  write(6,*) ''
  !  write(6,*) '*****************************************'
  !  write(6,*) '***          !!! WARNING !!!          ***'
  !  write(6,*) '*****************************************'
  !  write(6,*) '*  invalid value for FLR_ON parameter!  *'
  !  write(6,*) '* ------------------------------------- *'
  !  write(6,*) '*    FLR_ON must have 1 or 0 value...   *'
  !  write(6,*) '*****************************************'
  endif


    

 if (mpime==0) then
 write(6,*) ' --------------------------------------------- '
 write(6,*) 'SUBROUTINE CONDINIT, B00, angolo, Lx, Ly, Lz :'
 write(6,805) B00, angolo,  x(nx)-x(1), y(ny) - y(1), z(nz) - z(1)
 write(6,*) ' --------------------------------------------- '
 endif

622  format(1x, 2(2x, 1i4), 2(2x, 1e11.4))
822  format(1x, 2(2x, 1f5.1), 2(2x, 1e11.4))
802  format(1x, 2(2x, 1e11.4))
803  format(1x, 3(2x, 1e11.4))
805  format(1x, 5(2x, 1e11.4))
 
end subroutine

